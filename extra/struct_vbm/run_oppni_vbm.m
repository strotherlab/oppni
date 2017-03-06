function run_oppni_vbm(InputStruct,DEOBLIQUE,datatype, append_dir)

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
    CODE_PATH = fileparts(which('run_oppni_vbm.m'));
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

if ~isstruct(InputStruct)
    [InputStruct] = Read_Input_struct_vbm(InputStruct);
end

%% Reading optional input arguments, or giving default assignments to variables

% check if data needs to be "de-obliqued" (default = off)
if nargin<2 || isempty(DEOBLIQUE)
    disp('assuming data are non-oblique');
    DEOBLIQUE = 0;
end
% check for input data type / convert to numeric
if nargin<3 || isempty(datatype) 
    disp('assuming T1 input data as default');
    datatype =  1;
elseif( strcmpi(datatype,'template') )
    datatype =  0;
elseif( strcmpi(datatype,'T1') )
    datatype =  1;
elseif( strcmpi(datatype,'T2') )
    datatype =  2;    
elseif( ~isnumeric(dataype) || (sum(datatype==[0 1 2])==0) )
    error('data type not recognized');
end

if nargin<4
    append_dir = [];
end

% To run on HPCVL
% setenv('PATH',[getenv('PATH') ':' FSL_PATH ':' AFNI_PATH ':' FSL_PATH '/bin/']);
% setenv('FSLDIR',FSL_PATH);
% unix(['source ' FSL_PATH '/etc/fslconf/fsl.sh']);
% setenv('FSLDIR',FSL_PATH);
setenv('FSLOUTPUTTYPE','NIFTI')

Nsubject = length(InputStruct); % Count the number of all runs and subjects
grouptemp_path = [InputStruct(1).run(1).Output_nifti_file_path '/struct_vbm/group_template'];
mkdir_r(grouptemp_path);

if( ~isempty(append_dir) )

    % handle accidental trailing slash
    if(strcmpi(append_dir(end),'/') || strcmpi(append_dir(end),'\')) append_dir=append_dir(1:end-1); end
    % full path of reference template
    grouptemp_path = [append_dir, '/struct_vbm/group_template'];
    
    %% STEP1 - brain segmentation and initial transform

    for ksub = 1:Nsubject

        ksub,
        % subdirectories -- each tissue segmentation
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/CSF_warp-append']); %% make directory
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/GM_warp-append']); %% make directory
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/WM_warp-append']); %% make directory

        % untransform path
        vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
        mkdir_r(vbm_path); %% make directory
        strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];

        % check for stripped T1
        if ~exist([strip_struct,'.nii'],'file') && ~exist([strip_struct,'.nii.gz'],'file')
            unix([FSL_PATH sprintf('bet2 %s %s -m', [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}],strip_struct)]);            
            unix(sprintf('gunzip %s.nii.gz',strip_struct));    
        else
            'already stripped',
        end        
        % check for segmentation, and run
        if ~exist([strip_struct,'_GM.nii'],'file') && ~exist([strip_struct,'_GM.nii.gz'],'file')
            unix([FSL_PATH sprintf('fast -R 0.3 -H 0.1 -t %s %s',num2str(datatype),strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_0 %s_CSF',strip_struct,strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_1 %s_GM',strip_struct,strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_2 %s_WM',strip_struct,strip_struct)]);
        else
            'already segmented',
        end
    end

    %%
    %%

    tissueT={'CSF','GM','WM'}; tissueF={'csf','gray','white'}; % <type and fsl label>

    for(ttype = 1:3)

        disp(['running (appended) vbm on ',tissueT{ttype},' tissue']);

        for ksub = 1:Nsubject
            ksub,
            % untransform path
            vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
            mkdir_r(vbm_path); %% make directory
            strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];
            %
            reg_struct = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/',tissueT{ttype},'_warp-append/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1}];
            % second nonlinear transform, to final template (check if exist first)
            if ~exist([reg_struct,'_to_T3_mod.nii'],'file') && ~exist([strip_struct,'_to_T3_mod.nii.gz'],'file')        
                unix([FSL_PATH sprintf('fsl_reg %s_%s %s/template_%s_nl_symm %s_to_T3 -fnirt "--config=GM_2_MNI152GM_2mm.cnf --jout=%s_JAC_T3"',strip_struct,tissueT{ttype},grouptemp_path,tissueT{ttype},reg_struct,reg_struct)]);
                unix([FSL_PATH sprintf('fslmaths %s_to_T3 -mul %s_JAC_T3 %s_to_T3_mod -odt float',reg_struct,reg_struct,reg_struct)]);
            end
        end

    end
    
else

    %% STEP1 - brain segmentation and initial transform

    for ksub = 1:Nsubject

        ksub,
        % subdirectories -- each tissue segmentation
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/CSF_warp']); %% make directory
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/GM_warp']); %% make directory
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/WM_warp']); %% make directory

        % untransform path
        vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
        mkdir_r(vbm_path); %% make directory
        strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];

        % check for stripped T1
        if ~exist([strip_struct,'.nii'],'file') && ~exist([strip_struct,'.nii.gz'],'file')
            unix([FSL_PATH sprintf('bet2 %s %s -m', [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}],strip_struct)]);            
            unix(sprintf('gunzip %s.nii.gz',strip_struct));    
        else
            'already stripped',
        end        
        % check for segmentation, and run
        if ~exist([strip_struct,'_GM.nii'],'file') && ~exist([strip_struct,'_GM.nii.gz'],'file')
            unix([FSL_PATH sprintf('fast -R 0.3 -H 0.1 -t %s %s',num2str(datatype),strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_0 %s_CSF',strip_struct,strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_1 %s_GM',strip_struct,strip_struct)]);
            unix([FSL_PATH sprintf('immv %s_pve_2 %s_WM',strip_struct,strip_struct)]);
        else
            'already segmented',
        end
    end

    %%
    %%

    tissueT={'CSF','GM','WM'}; tissueF={'csf','gray','white'}; % <type and fsl label>

    for(ttype = 1:3)

        disp(['running vbm on ',tissueT{ttype},' tissue']);

        tomerge=[];
        % 1.1 affine transformation to template
        for(ksub = 1:Nsubject)
            ksub,
            % untransform path
            vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
            mkdir_r(vbm_path); %% make directory
            strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];
            %
            reg_struct = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/',tissueT{ttype},'_warp/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1}];
            % affine transform (check if exist first)
            if ~exist([reg_struct,'_to_T.nii'],'file') && ~exist([reg_struct,'_to_T.nii.gz'],'file')        
                unix([FSL_PATH sprintf('fsl_reg %s_%s ${FSLDIR}/data/standard/tissuepriors/avg152T1_%s %s_to_T -a',strip_struct,tissueT{ttype},tissueF{ttype},reg_struct)]);
            end
            tomerge = [tomerge sprintf(' %s_to_T',reg_struct)];
        end
        % 1.2 generating affine group template estimate
        if ~exist([grouptemp_path,'/template_',tissueT{ttype},'_aff_symm.nii'],'file') && ~exist([grouptemp_path,'/template_',tissueT{ttype},'_aff_symm.nii.gz'],'file')        
            unix([FSL_PATH sprintf('fslmerge -t %s/concat_%s_to_T %s',grouptemp_path,tissueT{ttype},tomerge)]);
            unix([FSL_PATH sprintf('fslmaths %s/concat_%s_to_T -Tmean %s/template_%s_aff',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            unix([FSL_PATH sprintf('fslswapdim %s/template_%s_aff -x y z %s/template_%s_aff_flipped',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            unix([FSL_PATH sprintf('fslmaths %s/template_%s_aff -add %s/template_%s_aff_flipped -div 2 %s/template_%s_aff_symm',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            % delete concat and flipped templates
            unix(sprintf('rm %s/concat_%s_to_T*',grouptemp_path,tissueT{ttype}));
            unix(sprintf('rm %s/template_%s_aff_flipped*',grouptemp_path,tissueT{ttype}));
        end

        tomerge = []; %% list of volumes to merge
        % 2.1 fnirting, round-1
        for ksub = 1:Nsubject

            ksub,
            % untransform path
            vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
            mkdir_r(vbm_path); %% make directory
            strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];
            %
            reg_struct = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/',tissueT{ttype},'_warp/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1}];
            % nonlinear transform (check if exist first)
            if ~exist([reg_struct,'_to_T2.nii'],'file') && ~exist([reg_struct,'_to_T2.nii.gz'],'file')        
               unix([FSL_PATH sprintf('fsl_reg %s_%s %s/template_%s_aff_symm %s_to_T2 -fnirt "--config=GM_2_MNI152GM_2mm.cnf"',strip_struct,tissueT{ttype},grouptemp_path,tissueT{ttype},reg_struct)]);
            end

            tomerge = [tomerge sprintf(' %s_to_T2',reg_struct)];
        end
        % 2.2 generating nonlinear group template estimate
        if ~exist([grouptemp_path,'/template_',tissueT{ttype},'_nl_symm.nii'],'file') && ~exist([grouptemp_path,'/template_',tissueT{ttype},'_nl_symm.nii.gz'],'file')        
            unix([FSL_PATH sprintf('fslmerge -t %s/concat_%s_to_T2 %s',grouptemp_path,tissueT{ttype},tomerge)]);
            unix([FSL_PATH sprintf('fslmaths %s/concat_%s_to_T2 -Tmean %s/template_%s_nl',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            unix([FSL_PATH sprintf('fslswapdim %s/template_%s_nl -x y z %s/template_%s_nl_flipped',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            unix([FSL_PATH sprintf('fslmaths %s/template_%s_nl -add %s/template_%s_nl_flipped -div 2 %s/template_%s_nl_symm',grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype},grouptemp_path,tissueT{ttype})]);
            % delete concat and flipped templates
            unix(sprintf('rm %s/concat_%s_to_T2*',grouptemp_path,tissueT{ttype}));
            unix(sprintf('rm %s/template_%s_nl_flipped*',grouptemp_path,tissueT{ttype}));
        end

        for ksub = 1:Nsubject
            ksub,
            % untransform path
            vbm_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/untransformed'];
            mkdir_r(vbm_path); %% make directory
            strip_struct = [vbm_path, '/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1},'_strip'];
            %
            reg_struct = [InputStruct(ksub).run(1).Output_nifti_file_path '/struct_vbm/',tissueT{ttype},'_warp/',InputStruct(ksub).run(1).Output_nifti_file_prefix{1}];
            % second nonlinear transform, to final template (check if exist first)
            if ~exist([reg_struct,'_to_T3_mod.nii'],'file') && ~exist([strip_struct,'_to_T3_mod.nii.gz'],'file')        
                unix([FSL_PATH sprintf('fsl_reg %s_%s %s/template_%s_nl_symm %s_to_T3 -fnirt "--config=GM_2_MNI152GM_2mm.cnf --jout=%s_JAC_T3"',strip_struct,tissueT{ttype},grouptemp_path,tissueT{ttype},reg_struct,reg_struct)]);
                unix([FSL_PATH sprintf('fslmaths %s_to_T3 -mul %s_JAC_T3 %s_to_T3_mod -odt float',reg_struct,reg_struct,reg_struct)]);
            end
        end

    end

end