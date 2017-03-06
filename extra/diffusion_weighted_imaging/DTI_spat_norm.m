function DTI_spat_norm(InputStruct,TEMPLATE_VOL, append_dir)

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
if ~isempty(AFNI_PATH) && AFNI_PATH(end)~='/'
	AFNI_PATH = [AFNI_PATH '/'];
end
if ~isempty(FSL_PATH)  && FSL_PATH(end)~='/'
	FSL_PATH = [FSL_PATH '/'];
end

addpath(CODE_PATH)
addpath([CODE_PATH '/NIFTI_tools'])

if(nargin<2) %% default is not assuming appended files
    append_dir = [];
end

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
setenv('FSLOUTPUTTYPE','NIFTI');

%%% 
Nsubject = length(InputStruct); % Count the number of all runs and subjects
grouptemp_path = [InputStruct(1).run(1).Output_nifti_file_path '/dti_processed/group_template'];
mkdir_r(grouptemp_path);
%%%

%- NB current version approximates tbss protocol, albeit with additional affine template generation step -% 

if( ~isempty(append_dir) )
    
    disp('Running spatial normalization on pre-defined template!')
    
    % handle accidental trailing slash
    if(strcmpi(append_dir(end),'/') || strcmpi(append_dir(end),'\')) append_dir=append_dir(1:end-1); end
    % full path of reference template
    grouptemp_path = [append_dir, '/dti_processed/group_template'];
    
    for ksub = 1:Nsubject

        % path to dti data
        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        % unzip "fitted" FA (other?) for subsequent analysis
        unix(['gunzip ',InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'/*_fit_FA.nii.gz']);
        % spatial norm directory
        mkdir_r([dti_path,'/spat_norm-append']);
        % "stripped", eroded image file for cleaner normalization
        strip_struct = [dti_path, '/spat_norm-append/FA_cln'];
        % collect matrix dimensions of FA image
        anam = dir([dti_path '/*_fit_FA.nii']);
        anam = anam(1).name; %% name of FA file (HARDI or DTI)
        hdr  = load_nii_hdr([dti_path,'/',anam]);
        imdim= hdr.dime.dim(2:4);
        % run masking, erosion etc.
        unix([FSL_PATH 'fslmaths ',[dti_path,'/',anam],' -min 1 -ero -roi 1 ',num2str(imdim(1)-2),' 1 ',num2str(imdim(2)-2),' 1 ',num2str(imdim(3)-2),' 0 1 ',strip_struct]);
        % remask
        unix([FSL_PATH 'fslmaths ',strip_struct,' -bin ',strip_struct,'_remask']); 
        % --exclude for now?--
        %%unix([FSL_PATH 'fslmaths ',strip_struct,'_remask -dilD -dilD -sub 1 -abs -add ',strip_struct,'_remask ',strip_struct,'_remask -odt char']);
    end

    % --> skips the to_T (affine) and to_T2 (intermed non-lin) registration steps 
   
    for ksub = 1:Nsubject
        ksub,
        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        strip_struct = [dti_path, '/spat_norm-append/FA_cln'];
        % second nonlinear transform, to final template (check if exist first)
        if ~exist([strip_struct,'_to_T3.nii'],'file') && ~exist([strip_struct,'_to_T3.nii.gz'],'file')        
            unix([FSL_PATH sprintf('fsl_reg %s %s/template_FA_nl_symm %s_to_T3 -fnirt "--config=FA_2_FMRIB58_1mm.cnf"',strip_struct,grouptemp_path,strip_struct)]);
        end
    end

    otherlist={'FA','MO','MD','ADx','RDx'};
    %% masking and warping other volumes
    for ksub = 1:Nsubject

        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        strip_struct = [dti_path, '/spat_norm-append/FA_cln'];

        %%pull list of other files, apply transform
        for(i=1:length(otherlist))

            innorm = [dti_path,'/DTI_fit_',otherlist{i}];
            outnorm = [dti_path,'/spat_norm-append/DTI_fit_',otherlist{i}];

            if( ~exist( [outnorm,'_to_T3.nii'], 'file' ) && ~exist( [outnorm,'_to_T3.nii.gz'], 'file' ) )
                % second nonlinear transform, to final template (check if exist first)
                unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
            end
        end
    end 
    
else %% ================standard, full registration =======================

    %% cleaning up FA volumes & masks

    tomerge=[];
    for ksub = 1:Nsubject

        % path to dti data
        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        % unzip "fitted" FA (other?) for subsequent analysis
        unix(['gunzip ',InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'/*_fit_FA.nii.gz']);
        % spatial norm directory
        mkdir_r([dti_path,'/spat_norm']);
        % "stripped", eroded image file for cleaner normalization
        strip_struct = [dti_path, '/spat_norm/FA_cln'];
        % collect matrix dimensions of FA image
        anam = dir([dti_path '/*_fit_FA.nii']);
        anam = anam(1).name; %% name of FA file (HARDI or DTI)
        hdr  = load_nii_hdr([dti_path,'/',anam]);
        imdim= hdr.dime.dim(2:4);
        % run masking, erosion etc.
        unix([FSL_PATH 'fslmaths ',[dti_path,'/',anam],' -min 1 -ero -roi 1 ',num2str(imdim(1)-2),' 1 ',num2str(imdim(2)-2),' 1 ',num2str(imdim(3)-2),' 0 1 ',strip_struct]);
        % remask
        unix([FSL_PATH 'fslmaths ',strip_struct,' -bin ',strip_struct,'_remask']); 
        % --exclude for now?--
        %%unix([FSL_PATH 'fslmaths ',strip_struct,'_remask -dilD -dilD -sub 1 -abs -add ',strip_struct,'_remask ',strip_struct,'_remask -odt char']);

        % affine transform (check if exist first)
        if ~exist([strip_struct,'_to_T.nii'],'file') && ~exist([strip_struct,'_to_T.nii.gz'],'file')        
            unix([FSL_PATH sprintf('fsl_reg %s ${FSLDIR}/data/standard/FMRIB58_FA_1mm %s_to_T -a',strip_struct,strip_struct)]);
        end

        tomerge = [tomerge sprintf(' %s_to_T',strip_struct)];    
    end

    %% generating affine group template estimate
    if ~exist([grouptemp_path,'/template_FA_aff_symm.nii'],'file') && ~exist([grouptemp_path,'/template_FA_aff_symm.nii.gz'],'file')        
        unix([FSL_PATH sprintf('fslmerge -t %s/concat_FA_to_T %s',grouptemp_path,tomerge)]);
        unix([FSL_PATH sprintf('fslmaths %s/concat_FA_to_T -Tmean %s/template_FA_aff',grouptemp_path,grouptemp_path)]);
        unix([FSL_PATH sprintf('fslswapdim %s/template_FA_aff -x y z %s/template_FA_aff_flipped',grouptemp_path,grouptemp_path)]);
        unix([FSL_PATH sprintf('fslmaths %s/template_FA_aff -add %s/template_FA_aff_flipped -div 2 %s/template_FA_aff_symm',grouptemp_path,grouptemp_path,grouptemp_path)]);
        % delete concat and flipped templates
        unix(sprintf('rm %s/concat_FA_to_T*',grouptemp_path));
        unix(sprintf('rm %s/template_FA_aff_flipped*',grouptemp_path));
    end

    tomerge = []; %% list of volumes to merge
    for ksub = 1:Nsubject

        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        strip_struct = [dti_path, '/spat_norm/FA_cln'];
        % nonlinear transform (check if exist first)
        if ~exist([strip_struct,'_to_T2.nii'],'file') && ~exist([strip_struct,'_to_T2.nii.gz'],'file')        
           unix([FSL_PATH sprintf('fsl_reg %s %s/template_FA_aff_symm %s_to_T2 -fnirt "--config=FA_2_FMRIB58_1mm.cnf"',strip_struct,grouptemp_path,strip_struct)]);
        end

        tomerge = [tomerge sprintf(' %s_to_T2',strip_struct)];
    end


    %% generating nonlinear group template estimate
    if ~exist([grouptemp_path,'/template_FA_nl_symm.nii'],'file') && ~exist([grouptemp_path,'/template_FA_nl_symm.nii.gz'],'file')        
        unix([FSL_PATH sprintf('fslmerge -t %s/concat_FA_to_T2 %s',grouptemp_path,tomerge)]);
        unix([FSL_PATH sprintf('fslmaths %s/concat_FA_to_T2 -Tmean %s/template_FA_nl',grouptemp_path,grouptemp_path)]);
        unix([FSL_PATH sprintf('fslswapdim %s/template_FA_nl -x y z %s/template_FA_nl_flipped',grouptemp_path,grouptemp_path)]);
        unix([FSL_PATH sprintf('fslmaths %s/template_FA_nl -add %s/template_FA_nl_flipped -div 2 %s/template_FA_nl_symm',grouptemp_path,grouptemp_path,grouptemp_path)]);
        % delete concat and flipped templates
        unix(sprintf('rm %s/concat_FA_to_T2*',grouptemp_path));
        unix(sprintf('rm %s/template_FA_nl_flipped*',grouptemp_path));
    end

    for ksub = 1:Nsubject
        ksub,
        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        strip_struct = [dti_path, '/spat_norm/FA_cln'];
        % second nonlinear transform, to final template (check if exist first)
        if ~exist([strip_struct,'_to_T3.nii'],'file') && ~exist([strip_struct,'_to_T3.nii.gz'],'file')        
            unix([FSL_PATH sprintf('fsl_reg %s %s/template_FA_nl_symm %s_to_T3 -fnirt "--config=FA_2_FMRIB58_1mm.cnf"',strip_struct,grouptemp_path,strip_struct)]);
        end
    end

    otherlist={'FA','MO','MD','ADx','RDx'};
    %% masking and warping other volumes
    for ksub = 1:Nsubject

        dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
        strip_struct = [dti_path, '/spat_norm/FA_cln'];

        %%pull list of other files, apply transform
        for(i=1:length(otherlist))

            innorm = [dti_path,'/DTI_fit_',otherlist{i}];
            outnorm = [dti_path,'/spat_norm/DTI_fit_',otherlist{i}];

            if( ~exist( [outnorm,'_to_T3.nii'], 'file' ) && ~exist( [outnorm,'_to_T3.nii.gz'], 'file' ) )
                % second nonlinear transform, to final template (check if exist first)
                unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
            end
        end
    end

end