function ASL_prepare( InputStruct, DEOBLIQUE )

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
if ~isempty(AFNI_PATH) && AFNI_PATH(end)~='/'
	AFNI_PATH = [AFNI_PATH '/'];
end
if ~isempty(FSL_PATH)  && FSL_PATH(end)~='/'
	FSL_PATH = [FSL_PATH '/'];
end

addpath(CODE_PATH)
addpath([CODE_PATH '/NIFTI_tools'])
read_version;

% parse input structure if not already done
if( ~isstruct(InputStruct) )
    [InputStruct] = Read_Input_ASL(InputStruct);
end

disp('now running asl preparation...');

% iterate through all subjects/runs
for ksub = 1:numel(InputStruct)

    inname  = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1},'.nii'];
    outstr  = [InputStruct(ksub).run(1).Output_nifti_file_path,'/asl_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1}];
    
    mkdir_r(outstr)

    % first, smooth file:
    unix([AFNI_PATH '3dmerge -prefix ' outstr '/raw_smo6.nii -doall -1blur_fwhm 6 ' inname]);   
    % generate perfusion estimates
    asl_perfusion_est( [outstr '/raw_smo6.nii'],  1, [], 2, [], -1, [outstr,'/proc'] ); %% v. crude mask approx
    
    typelist={'aCBF','fCBF','BOLD'};
    if( DEOBLIQUE>0 )
        for(i=1:numel(typelist))
            % deobliqueing
            unix([AFNI_PATH '3dWarp -oblique2card -prefix ' outstr '/proc_',typelist{i},'_deob.nii -cubic ' outstr '/proc_',typelist{i},'.nii']);
            unix([AFNI_PATH '3dWarp -oblique2card -prefix ' outstr '/proc_',typelist{i},'_avg_deob.nii -cubic ' outstr '/proc_',typelist{i},'_avg.nii']);
            delete([outstr '/proc_',typelist{i},'.nii']);
            delete([outstr '/proc_',typelist{i},'_avg.nii']);
        end
        % masking the volume
        unix([AFNI_PATH '3dAutomask -prefix ' outstr '/bold_mask.nii ' outstr,'/proc_BOLD_deob.nii']);
    else   
        % masking the volume
        unix([AFNI_PATH '3dAutomask -prefix ' outstr '/bold_mask.nii ' outstr,'/proc_BOLD.nii']);
    end

end

