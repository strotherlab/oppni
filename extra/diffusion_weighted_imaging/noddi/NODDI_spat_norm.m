function NODDI_spat_norm(InputStruct,reference_file)

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
    [InputStruct] = Read_Input_noddi(InputStruct);
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
grouptemp_path = [InputStruct(1).run(1).Output_nifti_file_path '/noddi_coreg/group_template'];
mkdir_r(grouptemp_path);
%%%

%- NB current version approximates tbss protocol, albeit with additional affine template generation step -% 

%% cleaning up ODI volumes & masks

tomerge=[];
for ksub = 1:Nsubject

    % path to dti data
    in_path = [InputStruct(ksub).run(1).Input_nifti_file_path '/',InputStruct(ksub).run(1).Input_nifti_file_prefix];   
    out_path= [InputStruct(ksub).run(1).Output_nifti_file_path '/noddi_coreg/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
    mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/noddi_coreg']);
    
    % "stripped", eroded image file for cleaner normalization
    %
    % collect matrix dimensions of ODI image
    hdr  = load_nii_hdr([in_path, '_odi.nii']);
    imdim= hdr.dime.dim(2:4);
    % run masking, erosion etc.
    unix([FSL_PATH 'fslmaths ',in_path,'_odi -min 1 -ero -roi 1 ',num2str(imdim(1)-2),' 1 ',num2str(imdim(2)-2),' 1 ',num2str(imdim(3)-2),' 0 1 ',out_path,'_odi_cln']);
    % remask
    unix([FSL_PATH 'fslmaths ',out_path,'_odi_cln -bin ',out_path, '_odi_remask']); 
    % taking the complement
    unix([FSL_PATH 'fslmaths ',out_path,'_odi_cln -mul -1 -add 1 -mul ',out_path,'_odi_remask ',out_path,'_odiC']);
    
    % name of eroded complement
    strip_struct = [out_path,'_odiC'];
    
    % affine transform (check if exist first)
    if ~exist([strip_struct,'_to_T.nii'],'file') && ~exist([strip_struct,'_to_T.nii.gz'],'file')        
        unix([FSL_PATH sprintf('fsl_reg %s ${FSLDIR}/data/standard/FMRIB58_FA_1mm %s_to_T -a',strip_struct,strip_struct)]);
    end

    tomerge = [tomerge sprintf(' %s_to_T',strip_struct)];    
end

%% generating affine group template estimate
if ~exist([grouptemp_path,'/template_ODIC_aff_symm.nii'],'file') && ~exist([grouptemp_path,'/template_ODIC_aff_symm.nii.gz'],'file')        
    unix([FSL_PATH sprintf('fslmerge -t %s/concat_ODIC_to_T %s',grouptemp_path,tomerge)]);
    unix([FSL_PATH sprintf('fslmaths %s/concat_ODIC_to_T -Tmean %s/template_ODIC_aff',grouptemp_path,grouptemp_path)]);
    unix([FSL_PATH sprintf('fslswapdim %s/template_ODIC_aff -x y z %s/template_ODIC_aff_flipped',grouptemp_path,grouptemp_path)]);
    unix([FSL_PATH sprintf('fslmaths %s/template_ODIC_aff -add %s/template_ODIC_aff_flipped -div 2 %s/template_ODIC_aff_symm',grouptemp_path,grouptemp_path,grouptemp_path)]);
    % delete concat and flipped templates
    unix(sprintf('rm %s/concat_ODIC_to_T*',grouptemp_path));
    unix(sprintf('rm %s/template_ODIC_aff_flipped*',grouptemp_path));
end

tomerge = []; %% list of volumes to merge
for ksub = 1:Nsubject

    % path to dti data
    in_path = [InputStruct(ksub).run(1).Input_nifti_file_path '/',InputStruct(ksub).run(1).Input_nifti_file_prefix];   
    out_path= [InputStruct(ksub).run(1).Output_nifti_file_path '/noddi_coreg/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
    % name of eroded complement
    strip_struct = [out_path,'_odiC'];
    
    % nonlinear transform (check if exist first)
    if ~exist([strip_struct,'_to_T2.nii'],'file') && ~exist([strip_struct,'_to_T2.nii.gz'],'file')        
       unix([FSL_PATH sprintf('fsl_reg %s %s/template_ODIC_aff_symm %s_to_T2 -fnirt "--config=FA_2_FMRIB58_1mm.cnf"',strip_struct,grouptemp_path,strip_struct)]);
    end

    tomerge = [tomerge sprintf(' %s_to_T2',strip_struct)];
end


%% generating nonlinear group template estimate
if ~exist([grouptemp_path,'/template_ODIC_nl_symm.nii'],'file') && ~exist([grouptemp_path,'/template_ODIC_nl_symm.nii.gz'],'file')        
    unix([FSL_PATH sprintf('fslmerge -t %s/concat_ODIC_to_T2 %s',grouptemp_path,tomerge)]);
    unix([FSL_PATH sprintf('fslmaths %s/concat_ODIC_to_T2 -Tmean %s/template_ODIC_nl',grouptemp_path,grouptemp_path)]);
    unix([FSL_PATH sprintf('fslswapdim %s/template_ODIC_nl -x y z %s/template_ODIC_nl_flipped',grouptemp_path,grouptemp_path)]);
    unix([FSL_PATH sprintf('fslmaths %s/template_ODIC_nl -add %s/template_ODIC_nl_flipped -div 2 %s/template_ODIC_nl_symm',grouptemp_path,grouptemp_path,grouptemp_path)]);
    % delete concat and flipped templates
    unix(sprintf('rm %s/concat_ODIC_to_T2*',grouptemp_path));
    unix(sprintf('rm %s/template_ODIC_nl_flipped*',grouptemp_path));
end

for ksub = 1:Nsubject
    
    % path to dti data
    in_path = [InputStruct(ksub).run(1).Input_nifti_file_path '/',InputStruct(ksub).run(1).Input_nifti_file_prefix];   
    out_path= [InputStruct(ksub).run(1).Output_nifti_file_path '/noddi_coreg/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
    % name of eroded complement
    strip_struct = [out_path,'_odiC'];

    % second nonlinear transform, to final template (check if exist first)
    if ~exist([strip_struct,'_to_T3.nii'],'file') && ~exist([strip_struct,'_to_T3.nii.gz'],'file')        
        unix([FSL_PATH sprintf('fsl_reg %s %s/template_ODIC_nl_symm %s_to_T3 -fnirt "--config=FA_2_FMRIB58_1mm.cnf"',strip_struct,grouptemp_path,strip_struct)]);
    end
end

%%
%%
%%

otherlist={'ficvf','fiso','odi','fmin','kappa'};
%% masking and warping other volumes
for ksub = 1:Nsubject
    
    % path to dti data
    in_path = [InputStruct(ksub).run(1).Input_nifti_file_path '/',InputStruct(ksub).run(1).Input_nifti_file_prefix];   
    out_path= [InputStruct(ksub).run(1).Output_nifti_file_path '/noddi_coreg/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
    % name of eroded complement
    strip_struct = [out_path,'_odiC'];
    
    %%pull list of other files, apply transform
    for(i=1:length(otherlist))

        innorm  = [in_path,'_',otherlist{i}];
        outnorm = [out_path,'_',otherlist{i}];
        % second nonlinear transform, to final template (check if exist first)
        unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
    end
end

