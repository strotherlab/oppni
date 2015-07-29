function spatial_normalization(InputStruct,reference_file,input_voxelsize,flag_step,DEOBLIQUE)
%==========================================================================
% SPATIAL_NORMALIZATION : registers data to common template
%==========================================================================
%
% SYNTAX:
%
%   spatial_normalization(InputStruct,reference_file,input_voxelsize,flag_step,DEOBLIQUE)
%
% INPUT:
%
%  InputStruct = subject input textfile
%  reference_file = path+name of reference anatomical
%  input_voxelsize = desired voxel dimensions for final resampling (functional data)
%
%  flag_step=0 : runs flag= 1, 2
%  flag_step=1 : STRUCTRUCL to REFERENCE
%  flag_step=2 : fMRI RESULTS to REFERENCE
%  flag_step=3 : fMRI preprocessed data (step1) to REFERENCE
%  it runs when the switch --dospnormfirst used in Run_Pipeline.py 
%
%  DEOBLIQUE = binary flag, corrects for scans acquired at oblique angles
%

global NUMBER_OF_CORES MULTI_RUN_INPUTFILE CODE_PROPERTY
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if (~exist('OCTAVE_VERSION','builtin') && exist('maxNumCompThreads'))
    maxNumCompThreads(NUMBER_OF_CORES);
end

global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('spatial_normalization.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end
if isempty(AFNI_PATH) || isempty(FSL_PATH)
    read_settings;
end
addpath(CODE_PATH)
addpath([CODE_PATH '/NIFTI_tools'])
read_version;

if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

if numel(InputStruct(1).run)>1
    MULTI_RUN_INPUTFILE = true;
end

%% Reading optional input arguments, or giving default assignments to variables

% check if nifti-output option is specified
if nargin<4 || isempty(flag_step)
    flag_step = 0;
else
    if isnumeric(flag_step)
        flag_step = num2str(flag_step);
    end
    flag_step = str2num(flag_step(1)); 
end
% check if data needs to be "de-obliqued" (default = off)
if nargin<5 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
end

% To run on HPCVL
% setenv('PATH',[getenv('PATH') ':' FSL_PATH ':' AFNI_PATH ':' FSL_PATH '/bin/']);
% setenv('FSLDIR',FSL_PATH);
% unix(['source ' FSL_PATH '/etc/fslconf/fsl.sh']);
% setenv('FSLDIR',FSL_PATH);
setenv('FSLOUTPUTTYPE','NIFTI')

for ksub = 1:numel(InputStruct)
    for krun = 1:numel(InputStruct(ksub).run)
        if ~exist(InputStruct(ksub).run(1).STRUCT_File,'file')
            sge_exit(100,sprintf('Spatial normalization failed due to the following error:\n The structural image %s not found',InputStruct(ksub).run(1).STRUCT_File));
        end
    end
end

Nsubject = length(InputStruct); % Count the number of all runs and subjects
for ksub = 1:Nsubject
    mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm']);
end

% Find the transforms
if flag_step==0 || flag_step==1
    
    % check whether voxel size is provided
    % the voxel size is used in the final normalized nifti images
    % the voxel size format should be a string, e.g. voxelsize='0.4 0.4 0.6'
    if nargin>=3
        voxelsize=input_voxelsize;
        if isnumeric(voxelsize)  % the input is voxel size
            v3 = [];
            for i = 1:length(voxelsize)
                v3 = [v3 ' ' num2str(voxelsize(i))];
            end
            voxelsize = v3;
            voxelsize_type = 1; % type = numeric 3d input
        else
            
            if(isempty(strfind(voxelsize,' ')))
                 v2 = cellstr(voxelsize);
            else v2 = regexp(voxelsize,' ','split');
            end
            
            if length(v2)==1
                v2 = repmat(v2,1,3);
            end
            v3 = [];
            for i = 1:length(v2)
                if ~isempty(v2{i})
                    if ~isempty(str2num(v2{i}))
                        v3 = [v3 ' ' v2{i}];
                    end
                end
            end
            voxelsize = v3;
            voxelsize_type = 1; % type = string input, split into cells
        end
        if isempty(voxelsize)   % the input is the master file
            if strcmpi(input_voxelsize,'NONE')
                  voxelsize = [];
                  voxelsize_type = 0; % otherwise --> use original input format
            else
                voxelsize = input_voxelsize;
                voxelsize_type = 2;
            end
        end
    else
        voxelsize = [];
        voxelsize_type = 0; % otherwise --> use original input format
    end
    
    if ~isempty(voxelsize)
        display(sprintf('output voxel sizes: %s',voxelsize));
    else
        display('output voxel sizes: remain intact');
    end
    
%     % skull stripping the reference --> TURNED OFF, AS TEMPLATES ARE USUALLY PRE-STRIPPED
%     [path_temp,name,ext] = fileparts(InputStruct(1).run(1).STRUCT_File);
%     STRUCT_Name = name;
%     [ref_path,ref_name,ref_ext] = fileparts(reference_file);
%     ref_path2 = [InputStruct(1).run(1).Output_nifti_file_path '/spat_norm'];
%     reference_file2 = sprintf('%s/%s_%s_brain%s',ref_path2,ref_name,STRUCT_Name,ref_ext);
%     if ~exist(reference_file2,'file');
%         unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',reference_file2,reference_file)]);
%     end
%     reference_file = reference_file2;
        
    Nsubject = length(InputStruct);
    % go through subjects, create transforms
    for ksub = 1:Nsubject
        
        [path_temp,name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        STRUCT_Name = name;
        
        % this line preven re-registeration of T1 to Reference (faster code)
        strip_struct = sprintf('%s/intermediate_processed/spat_norm/%s_strip.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        % check for stripped T1
        if ~exist(strip_struct,'file')
            
            % check if data is oblique; correct if requested
            if(DEOBLIQUE==1)
                disp('deobliquing!');
                struct_debobl = sprintf('%s/intermediate_processed/spat_norm/%s_deob.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
                unix([AFNI_PATH '3dWarp -oblique2card -prefix ' struct_debobl ' -cubic ' InputStruct(ksub).run(1).STRUCT_File]);
                unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',strip_struct,struct_debobl)]);
            else
                unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',strip_struct,InputStruct(ksub).run(1).STRUCT_File)]);
            end
        end
        % get transformation of T1 to reference volume
        trans_t1_ref = sprintf('%s/intermediate_processed/spat_norm/Transmat_T1toREF_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        if ~exist(trans_t1_ref,'file')
            unix([FSL_PATH sprintf('flirt -in %s/intermediate_processed/spat_norm/%s_strip.nii -ref %s -out %s/intermediate_processed/spat_norm/%s_T1toREF.nii -omat %s/intermediate_processed/spat_norm/Transmat_T1toREF_%s.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,reference_file,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
        % unzip to ensure afni compatibility
        if exist(sprintf('%s/intermediate_processed/spat_norm/%s_T1toREF.nii.gz',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name),'file')
            unix(['gunzip ' sprintf('%s/intermediate_processed/spat_norm/%s_T1toREF.nii.gz',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
        % check if normed, downsampled T1 exists --> and create it as a reference for functional data 
        if ~exist(sprintf('%s/intermediate_processed/spat_norm/%s_T1toREF_downsamp.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name),'file')
            if voxelsize_type==0
                hdr = load_nii_hdr([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(ksub).run(1).subjectprefix(2:end) '_baseproc.nii']);
                dim1 = hdr.dime.dim(2);dim2 = hdr.dime.dim(3);dim3 = hdr.dime.dim(4);
                pixdim1 = hdr.dime.pixdim(2);pixdim2 = hdr.dime.pixdim(3);pixdim3 = hdr.dime.pixdim(4);pixdim4 = hdr.dime.pixdim(5);
                
                % create eye.mat
                eye_file = [InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/eye.mat'];
                File = fopen(eye_file,'w');
                fprintf(File,'1 0 0 0\n');fprintf(File,'0 1 0 0\n');
                fprintf(File,'0 0 1 0\n');fprintf(File,'0 0 0 1\n');
                fclose(File);
                
                % create blank vol
                unix([FSL_PATH sprintf('fslcreatehd %.1f %.1f %.1f 1 %d %d %d %d 0 0 0 16 %s/intermediate_processed/spat_norm/blankvol.nii',dim1,dim2,dim3,pixdim1,pixdim2,pixdim3,pixdim4,InputStruct(ksub).run(1).Output_nifti_file_path)]);
                unix([FSL_PATH 'flirt -in ' InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/' STRUCT_Name '_T1toREF.nii -applyxfm -interp sinc -ref ' InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/blankvol.nii -init ' InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/eye.mat -out ' InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/' STRUCT_Name '_T1toREF_downsamp.nii']);
                
            elseif voxelsize_type==1 % resample to chosen voxel size
                unix([AFNI_PATH sprintf('3dresample -dxyz%s -inset %s/intermediate_processed/spat_norm/%s_T1toREF.nii -prefix %s/intermediate_processed/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
                if exist([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz'],'file')
                    unix(['gunzip -f -d ',InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz']);
                end
            elseif voxelsize_type==2
                unix([AFNI_PATH sprintf('3dresample -master %s -inset %s/intermediate_processed/spat_norm/%s_T1toREF.nii -prefix %s/intermediate_processed/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
                if exist([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz'],'file')
                    unix(['gunzip -f -d ',InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz']);
                end
            end
        end
    end
end

% Apply the estimated transforms to the Nifti files
if flag_step==0 || flag_step==2 || flag_step==3
    
    % go through list of subjects
    for ksub = 1:numel(InputStruct)
        [tmp,STRUCT_Name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        mean_file_name = sprintf('%s/intermediate_processed/spat_norm/mean_%s_baseproc.nii',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
        % check if mean functional volume exists, create if not
        if ~exist(mean_file_name,'file')
            nii      = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_baseproc.nii']);
            nii_mask = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/masks/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_mask.nii']);
            nii.img =  mean(double(nii.img),4);
            nii.img(nii_mask.img==0) = 0;
            nii.hdr.dime.dim([1 5]) = [3 1];
            nii.hdr.dime.datatype = 16;
            nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
            save_untouch_nii(nii,mean_file_name);
        end
        % spatial norm - transform mean epi volume to match stripped T1; create transform matrix
        unix([FSL_PATH sprintf('flirt -in %s -out %s/intermediate_processed/spat_norm/mean_%s_sNorm.nii -ref %s/intermediate_processed/spat_norm/%s_strip.nii -omat %s/intermediate_processed/spat_norm/Transmat_EPItoT1_%s.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6 -interp sinc -sincwidth 7 -sincwindow hanning', ... 
            mean_file_name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
        
        if numel(InputStruct(ksub).run)>1
           aligned_suffix= '_aligned';
        else
            aligned_suffix = [];
        end
        
        % go through list of runs for this subject        
        for krun=1:numel(InputStruct(ksub).run)

            % create net-transform matrix
            if ~exist(sprintf('%s/intermediate_processed/spat_norm/Transmat_EPItoREF_%s.mat',InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix),'file')
            unix([FSL_PATH sprintf('convert_xfm -omat %s/intermediate_processed/spat_norm/Transmat_EPItoREF_%s.mat -concat %s/intermediate_processed/spat_norm/Transmat_T1toREF_%s.mat %s/intermediate_processed/spat_norm/Transmat_EPItoT1_%s.mat', ... 
                InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
            end
            [path_temp,STRUCT_Name] = fileparts(InputStruct(ksub).run(krun).STRUCT_File);
            ref_file                = sprintf('%s/intermediate_processed/spat_norm/%s_T1toREF_downsamp.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,STRUCT_Name);
            transform               = sprintf('%s/intermediate_processed/spat_norm/Transmat_EPItoREF%s.mat',InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).subjectprefix);

            % spatially transform "baseproc" run, and create the mask
            if krun==1  
                input_nifti_file        = [InputStruct(ksub).run(krun).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc' aligned_suffix '.nii'];
                output_nifti_file       = [InputStruct(ksub).run(krun).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc' aligned_suffix '_sNorm.nii'];

                if ~exist(output_nifti_file,'file')
                unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp sinc -ref %s -init %s -out %s',input_nifti_file,ref_file,transform,output_nifti_file)]);
                unix(['gunzip -f -d ' output_nifti_file '.gz']);                
                end
                unix([AFNI_PATH sprintf('3dAutomask -prefix %s/intermediate_processed/spat_norm/%s_mask_sNorm.nii %s',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).subjectprefix(2:end),output_nifti_file)]);
            end
       
%%          

            if flag_step==3  % whether spatial normalization is performed initially
                input_nifti_file_path = strcat(InputStruct(ksub).run(krun).Output_nifti_file_path, '/intermediate_processed/afni_processed');
                input_nifti_file_list = strcat(input_nifti_file_path,'/*',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'*.nii');
                list = dir(input_nifti_file_list); % collect list of all files in specified directory
            else
                input_nifti_file_path = strcat(InputStruct(ksub).run(krun).Output_nifti_file_path,'/optimization_results');
                input_nifti_file_list = strcat(input_nifti_file_path,'/processed/*',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'*.nii');
                  list1 = dir(input_nifti_file_list); % collect list of all files in specified directory
                input_nifti_file_list = strcat(input_nifti_file_path,'/spms/*',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'*.nii');
                  list2 = dir(input_nifti_file_list); % collect list of all files in specified directory
                for(i=1:length(list1)) list1(i).name = ['processed/' list1(i).name]; end
                for(i=1:length(list2)) list2(i).name = ['spms/'      list2(i).name]; end
                list = [list1;list2]; % concatenated
            end
            
            for kpip = 1:length(list)
                if ~isempty(strfind(list(kpip).name,'sNorm'));
                    continue;
                end
                input_nifti_file        = strcat(input_nifti_file_path,'/',list(kpip).name);
                output_nifti_file       = strcat(input_nifti_file_path,'/',list(kpip).name(1:end-4),'_sNorm.nii');

                % if spatially normalized equivalent does not exist, create it!
                if ~exist(output_nifti_file,'file')
                    if ~exist([output_nifti_file,'.gz'],'file')
                    unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp sinc -ref %s -init %s -out %s',input_nifti_file,ref_file,transform,output_nifti_file)]);
                    end
                    
                    unix(['gunzip -f -d ' output_nifti_file '.gz']);

                    if ~exist(output_nifti_file,'file')
                        sge_exit(100,'spatial normalization step has failed, check above for the errors from AFNI and FSL modules');
                    end
                    
                    hdr_orig = load_nii_hdr(input_nifti_file);
                    nii = load_untouch_nii(output_nifti_file);
                    nii.hdr.hist.descrip = hdr_orig.hist.descrip;
                    if ~isempty(strfind(upper(hdr_orig.hist.descrip),'KEEPMEAN')) % If the discription of nifti files is KEEPMEAN
                        nii.img(nii.img<0) = 0;                      % remove small zero values from the output nifti files
                        save_untouch_nii(nii,output_nifti_file);
                    end

                end
            end
            if flag_step==3
                for kpip = 1:length(list)
                    if ~isempty(strfind(list(kpip).name,'sNorm'));
                        continue;
                    end
                    input_nifti_file        = strcat(input_nifti_file_path,'/',list(kpip).name);
                    delete(input_nifti_file);
                    output_nifti_file       = strcat(input_nifti_file_path,'/',list(kpip).name(1:end-4),'_sNorm.nii');
                    movefile(output_nifti_file,input_nifti_file,'f');
                end
            end
%%
        end
    end
end



