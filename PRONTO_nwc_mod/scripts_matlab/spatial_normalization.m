function spatial_normalization(InputStruct,reference_file,input_voxelsize,flag_step)

% flag=0 : runs flag= 1, 2
% flag=1 : STRUCTRUCL to REFERENCE
% flag=2 : fMRI RESULTS to REFERENCE
% flag=3 : fMRI preprocessed data (step1) to REFERENCE
% it runs when the switch --dospnormfirst used in Run_Pipeline.py 

global NUMBER_OF_CORES MULTI_RUN_INPUTFILE
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
    CODE_PATH = fileparts(which('Pipeline_PART1_afni_steps.m'));
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
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

if numel(InputStruct(1).run)>1
    MULTI_RUN_INPUTFILE = true;
end

if nargin<4
    flag_step = 0; % do both steps
else
    if isnumeric(flag_step)
        flag_step = num2str(flag_step);
    end
    flag_step = str2num(flag_step(1)); 
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
    mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/spat_norm']);
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
            voxelsize_type = 1;
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
            voxelsize_type = 1;
        end
        if isempty(voxelsize)   % the input is the master file
            voxelsize = input_voxelsize;
            voxelsize_type = 2;
        end
    else
        voxelsize = [];
        voxelsize_type = 0;
    end
    
    if ~isempty(voxelsize)
        display(sprintf('output voxel sizes: %s',voxelsize));
    else
        display('output voxel sizes: remain intact');
    end
    
    
    % skull stripping the reference
    [ref_path,ref_name,ref_ext] = fileparts(reference_file);
    ref_path2 = [InputStruct(1).run(1).Output_nifti_file_path '/spat_norm'];
    reference_file2 = sprintf('%s/%s_%s_brain%s',ref_path2,ref_name,InputStruct(1).run(1).STRUCT_Name,ref_ext);
    if ~exist(reference_file2,'file');
        unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',reference_file2,reference_file)]);
    end
    reference_file = reference_file2;
    
        
    Nsubject = length(InputStruct);
    for ksub = 1:Nsubject
        
        [path_temp,name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        STRUCT_Name = name;
        
        
        % this line preven re-registeration of T1 to Reference (faster
        % code)
        strip_struct = sprintf('%s/spat_norm/%s_strip.nii',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        if ~exist(strip_struct,'file')
            unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s -input %s',strip_struct,InputStruct(ksub).run(1).STRUCT_File)]);
        end
        trans_t1_ref = sprintf('%s/spat_norm/Transmat_T1toREF_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name);
        if ~exist(trans_t1_ref,'file')
            unix([FSL_PATH sprintf('flirt -in %s/spat_norm/%s_strip.nii -ref %s -out %s/spat_norm/%s_T1toREF.nii -omat %s/spat_norm/Transmat_T1toREF_%s.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp sinc -sincwidth 7 -sincwindow hanning',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,reference_file,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
        if voxelsize_type==0
            hdr = load_nii_hdr([InputStruct(ksub).run(1).Output_nifti_file_path '/' InputStruct(ksub).run(1).subjectprefix(2:end) '_baseproc.nii']);
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
            
        elseif voxelsize_type==1
            unix([AFNI_PATH sprintf('3dresample -dxyz%s -inset %s/spat_norm/%s_T1toREF.nii -prefix %s/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        else
            unix([AFNI_PATH sprintf('3dresample -master %s -inset %s/spat_norm/%s_T1toREF.nii -prefix %s/spat_norm/%s_T1toREF_downsamp.nii -rmode Cu',voxelsize,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name)]);
        end
            
    end
end

% Apply the estimated transforms to the Nifti files
if flag_step==0 || flag_step==2 || flag_step==3
    for ksub = 1:numel(InputStruct)
        [tmp,STRUCT_Name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
        mean_file_name = sprintf('%s/matfiles/mean_%s.nii',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
        if ~exist(mean_file_name,'file')
            nii      = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_baseproc.nii']);
            nii_mask = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/masks/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_mask.nii']);
            nii.img =  mean(double(nii.img),4);
            nii.img(nii_mask.img==0) = 0;
            nii.hdr.dime.dim([1 5]) = [3 1];
            nii.hdr.dime.datatype = 16;
            mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/matfiles']);
            save_untouch_nii(nii,mean_file_name);
        end
        unix([FSL_PATH sprintf('flirt -in %s -out %s/spat_norm/mean_%s_sNorm.nii -ref %s/spat_norm/%s_strip.nii -omat %s/spat_norm/Transmat_EPItoT1_%s.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6 -interp sinc -sincwidth 7 -sincwindow hanning', ... 
            mean_file_name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
        
        if numel(InputStruct(ksub).run)>1
           aligned_suffix= '_aligned';
        else
            aligned_suffix = [];
        end
                
        for krun=1:numel(InputStruct(ksub).run)

            unix([FSL_PATH sprintf('convert_xfm -omat %s/spat_norm/Transmat_EPItoREF_%s.mat -concat %s/spat_norm/Transmat_T1toREF_%s.mat %s/spat_norm/Transmat_EPItoT1_%s.mat', ... 
                InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);

            [path_temp,STRUCT_Name] = fileparts(InputStruct(ksub).run(krun).STRUCT_File);
            ref_file                = sprintf('%s/spat_norm/%s_T1toREF_downsamp.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,STRUCT_Name);
            transform               = sprintf('%s/spat_norm/Transmat_EPItoREF%s.mat',InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).subjectprefix);
            
            if krun==1  % Building a mask 
                input_nifti_file        = [InputStruct(ksub).run(krun).Output_nifti_file_path '/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc' aligned_suffix '.nii'];
                output_nifti_file       = [InputStruct(ksub).run(krun).Output_nifti_file_path '/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc' aligned_suffix '_sNorm.nii'];
                unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp sinc -ref %s -init %s -out %s',input_nifti_file,ref_file,transform,output_nifti_file)]);
                unix([AFNI_PATH sprintf('3dAutomask -prefix %s/spat_norm/%s_mask_sNorm.nii %s',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).subjectprefix(2:end),output_nifti_file)]);
            end
            if flag_step==3  % whether spatial normalization is performed initially
                input_nifti_file_path = strcat(InputStruct(ksub).run(krun).Output_nifti_file_path);
                input_nifti_file_list = strcat(input_nifti_file_path,'/*',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'*.nii');
            else
                input_nifti_file_path = strcat(InputStruct(ksub).run(krun).Subject_OutputDirectory,'/niftis_',InputStruct(ksub).run(krun).Output_nifti_file_prefix);
                input_nifti_file_list = strcat(input_nifti_file_path,'/','*_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'*.nii');
            end
            list = dir(input_nifti_file_list);
            for kpip = 1:length(list)
                if ~isempty(strfind(list(kpip).name,'sNorm'));
                    continue;
                end
                input_nifti_file        = strcat(input_nifti_file_path,'/',list(kpip).name);
                output_nifti_file       = strcat(input_nifti_file_path,'/',list(kpip).name(1:end-4),'_sNorm.nii');
                if ~exist(output_nifti_file,'file')
                    unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp sinc -ref %s -init %s -out %s',input_nifti_file,ref_file,transform,output_nifti_file)]);
                    if exist(['gzip -f -d ' output_nifti_file '.gz'],'file')
                        unix(['gzip -f -d ' output_nifti_file '.gz']);
                    end

                    if ~exist(output_nifti_file,'file')
                        sge_exit(100,'spatial normalization step has failed, check above for the errors from AFNI and FSL modules');
                    end
                    
                    hdr_orig = load_nii_hdr(input_nifti_file);
                    nii = load_untouch_nii(output_nifti_file);
                    nii.hdr.hist.descrip = hdr_orig.hist.descrip;
                    if ~isempty(strfind(upper(hdr_orig.hist.descrip),'KEEPMEAN')) % If the discription of nifti files is PRONTO-KEEPMEAN
                        nii.img(nii.img<0) = 0;                      % remove small zero values from the output nifti files
                    else
                        nii.img = bsxfun(@minus,double(nii.img),mean(double(nii.img),4));
                        nii.hdr.dime.datatype = 16;  
                    end
                    save_untouch_nii(nii,output_nifti_file);

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
        end
    end
end



