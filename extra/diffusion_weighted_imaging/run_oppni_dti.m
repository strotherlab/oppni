function run_oppni_dti( subject_inputs, TEMPLATE_VOL, VOXDIMS, DEOBLIQUE, run_NODDI, list, append_dir )

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

if nargin<4 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
end
if nargin<5 || isempty( run_NODDI )
    run_NODDI=false;
end
if nargin<6
    list=[];
end
if nargin<7
    append_dir = [];
end


if(~run_NODDI)
    
    % run preprocessing pipeline
    DTI_prepare(subject_inputs);

    % run spatial normalization
    if(normflag)
    DTI_spat_norm(subject_inputs,TEMPLATE_VOL,append_dir);
    end
end

if(run_NODDI)
    run_noddi_wmtract(subject_inputs,TEMPLATE_VOL, list,append_dir);
end

%%% make it so that default unless otherwise specified (TEMPLATE,VOXDIMS)
%%% ... need explicit instruct to not run spatnorm