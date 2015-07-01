function output = module_erGNB( datamat, split_info )
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
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

if( ~isfield(split_info,'WIND')   || isempty(split_info.WIND)   ) 
    disp('erGNB uses default window-size WIND=10');
    split_info.WIND =10; 
end
if( ~isfield(split_info,'Nblock') || isempty(split_info.Nblock) ) 
    disp('erGNB uses default number of blocks Nblock=4');
    split_info.Nblock=4; 
end
if( ~isfield(split_info,'norm')   || isempty(split_info.norm)   ) 
    disp('erGNB uses default normalization turned on (norm=1)');
    split_info.norm  =1; 
end

params.TR     = split_info.TR_MSEC;
params.delay  = split_info.TR_MSEC./2;
params.WIND   = split_info.WIND;
params.Nblock = split_info.Nblock;
params.norm   = 1;

% % % count #leftovers, if we do even multiples of splits
% % leftover = rem( size(datamat,2) , params.Nblock);
% % % trim "overhang" from the data matrix
% % datamat_trim  = datamat(:,1:end-leftover);

%  Time-windowed averaging of data:
% *NB: simple averaging does nearest-neighbour interpolation to nearest TR
%      interval. You can replace with "interp_averaging_for_ER" to do
%      linear time interpolation; the tends to increase reproducibility,
%      but at the apparent expense of decreased prediction!
%
[windowAvg] = simple_averaging_for_ER( datamat, split_info.onsetlist, params );
% run optimized analysis
out = ergnb_optimization( windowAvg , split_info.spat_weight);

% keep signed image for representative spatial map
output.images           = out.sens_hrf1;
% get all metrics --> but D( R[signed_1st_hrf],  P[classification] )
output.metrics.R_global =  out.R_global;
output.metrics.R_1hrf   =  out.R_hrf1;
output.metrics.P_class  =  out.P_class;
output.metrics.P_rms    =  out.P_rms;
output.metrics.dPR      = -sqrt( (1-out.R_hrf1).^2 + (1-out.P_class).^2 );
% special format of temporal output
tempstruct.hrf          = out.HRF1;
tempstruct.sens         = out.sens_global;
tempstruct.pmap         = out.pmap;
output.temp             = tempstruct;



