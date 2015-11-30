function Run_DTI_noSGE( subject_inputs, TEMPLATE_VOL, VOXDIMS, DEOBLIQUE )

%%% current script ignores TEMPLATE_VOL, VOXDIMS! 
%%% --> uses FSL template, in their native res.

addpath scripts_matlab;
addpath scripts_matlab/optimization;
addpath scripts_matlab/NIFTI_tools;

% if nargin<3 || isempty(TEMPLATE_VOL) || isempty(VOXDIMS)
%      normflag = false;
% else normflag = true;
% end
normflag = true;

if nargin<5 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
end

% run preprocessing pipeline
%DTI_prepare(subject_inputs);

% run spatial normalization
if(normflag)
DTI_spat_norm(subject_inputs,TEMPLATE_VOL);
end

%%% make it so that default unless otherwise specified (TEMPLATE,VOXDIMS)
%%% ... need explicit instruct to not run spatnorm