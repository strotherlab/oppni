function DTI_spat_norm(InputStruct,reference_file,input_voxelsize,DEOBLIQUE)

% flag=0 : runs flag= 1, 2
% flag=1 : STRUCTRUCL to REFERENCE
% flag=2 : fMRI RESULTS to REFERENCE
% flag=3 : fMRI preprocessed data (step1) to REFERENCE
% it runs when the switch --dospnormfirst used in Run_Pipeline.py 

global NUMBER_OF_CORES
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

if ~isstruct(InputStruct)
    [InputStruct] = Read_Input_DTI(InputStruct);
end

%% Reading optional input arguments, or giving default assignments to variables

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
    if ~exist(InputStruct(ksub).run(1).STRUCT_File,'file')
        sge_exit(100,sprintf('Spatial normalization failed due to the following error:\n The structural image %s not found',InputStruct(ksub).run(1).STRUCT_File));
    end
end

Nsubject = length(InputStruct); % Count the number of all runs and subjects
for ksub = 1:Nsubject
    mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm']);
end

%% Find the transforms

    
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
            voxelsize = input_voxelsize;
            voxelsize_type = 2;
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
        strip_struct = sprintf('%s/spat_norm/%s_strip.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        % check for stripped T1
        if ~exist(strip_struct,'file')
            
            % check if data is oblique; correct if requested
            if(DEOBLIQUE==1)
                disp('deobliquing!');
                struct_debobl = sprintf('%s/spat_norm/%s_deob.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
                unix([AFNI_PATH '3dWarp -oblique2card -prefix ' struct_debobl ' -cubic ' InputStruct(ksub).run(1).STRUCT_File]);
                unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',strip_struct,struct_debobl)]);
            else
                unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',strip_struct,InputStruct(ksub).run(1).STRUCT_File)]);
            end
        end
        % get transformation of T1 to reference volume
        trans_t1_ref = sprintf('%s/spat_norm/Transmat_T1toREF_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        if ~exist(trans_t1_ref,'file')
            unix([FSL_PATH sprintf('flirt -in %s/spat_norm/%s_strip.nii -ref %s -out %s/spat_norm/%s_T1toREF.nii -omat %s/spat_norm/Transmat_T1toREF_%s.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,reference_file,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
        % unzip to ensure afni compatibility
        if exist(sprintf('%s/spat_norm/%s_T1toREF.nii.gz',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name),'file')
            unix(['gunzip ' sprintf('%s/spat_norm/%s_T1toREF.nii.gz',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
        % check if normed, downsampled T1 exists --> and create it as a reference for DTI data 
        if ~exist(sprintf('%s/spat_norm/%s_T1toREF_downsamp.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name),'file')
            if voxelsize_type==0
                hdr = load_nii_hdr([InputStruct(ksub).run(1).Output_nifti_file_path '/' InputStruct(ksub).run(1).Output_nifti_file_prefix '/dti_intermediate/',DTI_TYPE,'_fit_S0.nii.gz']);
                dim1 = hdr.dime.dim(2);dim2 = hdr.dime.dim(3);dim3 = hdr.dime.dim(4);
                pixdim1 = hdr.dime.pixdim(2);pixdim2 = hdr.dime.pixdim(3);pixdim3 = hdr.dime.pixdim(4);pixdim4 = hdr.dime.pixdim(5);
                
                % create eye.mat
                eye_file = [InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm/eye.mat'];
                File = fopen(eye_file,'w');
                fprintf(File,'1 0 0 0\n');fprintf(File,'0 1 0 0\n');
                fprintf(File,'0 0 1 0\n');fprintf(File,'0 0 0 1\n');
                fclose(File);
                
                % create blank vol
                unix([FSL_PATH sprintf('fslcreatehd %.1f %.1f %.1f 1 %d %d %d %d 0 0 0 16 %s/spat_norm/blankvol.nii',dim1,dim2,dim3,pixdim1,pixdim2,pixdim3,pixdim4,InputStruct(ksub).run(1).Output_nifti_file_path)]);
                unix([FSL_PATH 'flirt -in ' InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm/' STRUCT_Name '_T1toREF.nii -applyxfm -interp sinc -ref ' InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm/blankvol.nii -init ' InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm/eye.mat -out ' InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm/' STRUCT_Name '_T1toREF_downsamp.nii']);
                
            elseif voxelsize_type==1 % resample to chosen voxel size
                unix([AFNI_PATH sprintf('3dresample -dxyz%s -inset %s/spat_norm/%s_T1toREF.nii -prefix %s/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
                unix(['gunzip -f -d ',InputStruct(ksub).run(1).Output_nifti_file_path,'/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz']);
            elseif voxelsize_type==2
                unix([AFNI_PATH sprintf('3dresample -master %s -inset %s/spat_norm/%s_T1toREF.nii -prefix %s/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
                unix(['gunzip -f -d ',InputStruct(ksub).run(1).Output_nifti_file_path,'/spat_norm/',STRUCT_Name,'_T1toREF_downsamp.nii.gz']);
            end
        end
    end


%% Apply the estimated transforms to the Diffusion-weighted files

    disp('transforming diffusion-weighted images...');

    % go through list of subjects
    for ksub = 1:numel(InputStruct)
        
        ksub,

        % determine appended form of data
        if( InputStruct(ksub).run(1).MULTI_RUN_HARDI ) DTI_TYPE = 'HARDI';
        else                                           DTI_TYPE = 'DTI';
        end
                
        
        [tmp,STRUCT_Name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        mean_file_name = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix '/',DTI_TYPE,'_fit_S0.nii.gz'];
        % spatial norm - transform mean epi volume to match stripped T1; create transform matrix
        unix([FSL_PATH sprintf('flirt -in %s -out %s/spat_norm/%s_fit_S0_%s_sNorm.nii -ref %s/spat_norm/%s_strip.nii -omat %s/spat_norm/Transmat_dtiT2toT1_%s.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6 -interp sinc -sincwidth 7 -sincwindow hanning', ... 
            mean_file_name,InputStruct(ksub).run(1).Output_nifti_file_path,DTI_TYPE,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
        
        % create net-transform matrix
        if ~exist(sprintf('%s/spat_norm/Transmat_dtiT2toREF_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix),'file')
        unix([FSL_PATH sprintf('convert_xfm -omat %s/spat_norm/Transmat_dtiT2toREF_%s.mat -concat %s/spat_norm/Transmat_T1toREF_%s.mat %s/spat_norm/Transmat_dtiT2toT1_%s.mat', ... 
            InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
        end
        [path_temp,STRUCT_Name] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        ref_file                = sprintf('%s/spat_norm/%s_T1toREF_downsamp.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        transform               = sprintf('%s/spat_norm/Transmat_dtiT2toREF_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
        
        % spatially transform "non-vectorized" volumes
        listv = {'S0','L1','L2','L3','FA','MD','MO'};
        for(i=1:length(listv))
            input_nifti_file        = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix, '/',DTI_TYPE,'_fit_',listv{i},'.nii.gz'];
            output_nifti_file       = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix, '/',DTI_TYPE,'_fit_',listv{i},'_sNorm.nii'];

            if ~exist(output_nifti_file,'file')
            unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp sinc -ref %s -init %s -out %s',input_nifti_file,ref_file,transform,output_nifti_file)]);
            unix(['gunzip -f -d ' output_nifti_file '.gz']);                
            end
        end
        
        % spatially transform "vectorized" volumes
        listv = {'V1','V2','V3'};
        for(i=1:length(listv))
            input_nifti_file        = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix, '/',DTI_TYPE,'_fit_',listv{i},'.nii.gz'];
            output_nifti_file       = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix, '/',DTI_TYPE,'_fit_',listv{i},'_sNorm.nii'];

            if ~exist(output_nifti_file,'file')
            unix([FSL_PATH sprintf('vecreg -i %s -o %s -r %s -t %s',input_nifti_file,output_nifti_file,ref_file,transform)]);
            unix(['gunzip -f -d ' output_nifti_file '.gz']);                
            end
        end     
        
%%
    end
