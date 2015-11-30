function Run_ASL_noSGE( subject_inputs, TEMPLATE_VOL, VOXDIMS, DEOBLIQUE )

addpath scripts_matlab;
addpath scripts_matlab/optimization;
addpath scripts_matlab/NIFTI_tools;

if nargin<3 || isempty(TEMPLATE_VOL) || isempty(VOXDIMS)
     normflag = false;
else normflag = true;
end

if nargin<4 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
end

% run preprocessing pipeline
ASL_prepare(subject_inputs,DEOBLIQUE);

% run spatial normalization
if(normflag)
ASL_spat_norm(subject_inputs,TEMPLATE_VOL,VOXDIMS,DEOBLIQUE);
end
