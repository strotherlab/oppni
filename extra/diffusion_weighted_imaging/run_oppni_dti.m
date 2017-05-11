function run_oppni_dti( subject_inputs, EC_type, run_NODDI,corDiffParam, list, append_dir )

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

if nargin<3 || isempty( run_NODDI )
    run_NODDI=false;
end
if( ~isempty( run_NODDI ) && double(run_NODDI)~=0 && (nargin<4 || isempty(corDiffParam)) )
    error('NODDI correction parameter not determined');
end
if nargin<5
    list=[];
end
if nargin<6
    append_dir = [];
end


if(~run_NODDI)
    
    % run preprocessing pipeline
    DTI_prepare(subject_inputs,EC_type);

    % run spatial normalization
    if(normflag)
    DTI_spat_norm(subject_inputs,EC_type,append_dir);
    end
end

if(run_NODDI)
    run_noddi_wmtract(subject_inputs,EC_type,corDiffParam, list,append_dir);
end

%%% make it so that default unless otherwise specified (TEMPLATE,VOXDIMS)
%%% ... need explicit instruct to not run spatnorm
