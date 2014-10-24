function output = module_erGNB_group( datamat, split_info )
%
% =========================================================================
% MODULE_erGNB: module that performs event-related gaussian naive bayes analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_erGNB( datamat, split_info )
%
%

% initialized parameters
params.TR     = split_info{1}.TR_MSEC;
params.delay  = split_info{1}.TR_MSEC ./ 2;
params.WIND   = split_info{1}.WIND;
params.Nsplit = 1; %% now, only 1 split per subject!
params.norm   = 1;
% number of subjects
N_subject = length(datamat);
Nvox      = size(datamat{1},1);

% initialize data matrix
windowAvg = zeros( Nvox, split_info{1}.WIND, N_subject );
% load from averaged blocks from cell arrays
for(n=1:N_subject)
    % perform stimulus-locked averaging on each subject
    [windowAvg(:,:,n)] = simple_averaging_for_ER(datamat{n}, split_info{n}.onsetlist, params );
end

% % % count #leftovers, if we do even multiples of splits
% % leftover = rem( size(datamat,2) , params.Nsplit);
% % % trim "overhang" from the data matrix
% % datamat_trim  = datamat(:,1:end-leftover);

%  Time-windowed averaging of data:
% *NB: simple averaging does nearest-neighbour interpolation to nearest TR
%      interval. You can replace with "interp_averaging_for_ER" to do
%      linear time interpolation; the tends to increase reproducibility,
%      but at the apparent expense of decreased prediction!
%
% run optimized analysis
out = ergnb_optimization( windowAvg , split_info.spat_weight);

% keep signed image for representative spatial map
output.images           = out.sens_hrf1;
% get all metrics --> but D( R[signed_1st_hrf],  P[classification] )
output.metrics.R_global =  out.R_global;
output.metrics.R_1hrf   =  out.R_hrf1;
output.metrics.P_class  =  out.P_class;
output.metrics.P_rms    =  out.P_rms;
output.metrics.Dneg     = -sqrt( (1-out.R_hrf1).^2 + (1-out.P_class).^2 );
% special format of temporal output
tempstruct.hrf          = out.HRF1;
tempstruct.sens         = out.sens_global;
tempstruct.pmap         = out.pmap;
output.temp             = tempstruct;



