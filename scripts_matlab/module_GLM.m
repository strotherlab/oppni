function output = module_GLM( datamat, split_info )
%
% =========================================================================
% MODULE_GLM: module that performs General Linear Model analyses
% =========================================================================
%
%   Syntax:
%           output = module_erGLM( datamat, split_info )
%
%

% matrix dimensions
[Nvox Ntime] = size( datamat );

X_design = split_info.design_mat;

% BABAK: the convolve switch has been taken care of in interpret_contrast_list_str.m
% % convolution
% if( split_info.convolve == 0 )
%     % get unconvolved design
%     X_design = split_info.design_mat;
% else
%     % convolve with HRF (must convert into seconds!)
%     X_design = design_to_hrf( split_info.design_mat, split_info.TR_MSEC/1000, [5.0 15.0] );
% end

% split the data matrix
datasplit_1 = datamat(:,1:ceil(Ntime/2));
datasplit_2 = datamat(:,ceil(Ntime/2)+1:end);
% split the hrf-convolved design vector
X_design_1 = X_design(1:ceil(Ntime/2),:);
X_design_2 = X_design(ceil(Ntime/2)+1:end,:);
% run split glm analyses
out1 = GLM_model_fmri( datasplit_1, 0, [], X_design_1, [] );
out2 = GLM_model_fmri( datasplit_2, 0, [], X_design_2, [] );

% reweight the betas
beta1 = bsxfun(@times,out1.Beta_signl,split_info.spat_weight);
beta2 = bsxfun(@times,out2.Beta_signl,split_info.spat_weight);

% reproducibility and SPM:
[ RR, rSPMZ ] = get_rSPM( beta1, beta2, 1 );

% reproducibility
output.metrics.R = RR;
% optimal eigenimage
output.images  = rSPMZ;
% CV score timeseries, on unit-normed eigenimage
output.temp    = datamat'  * bsxfun(@rdivide,rSPMZ,sqrt(sum(rSPMZ.^2)));
