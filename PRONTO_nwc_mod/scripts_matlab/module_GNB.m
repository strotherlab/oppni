function output = module_GNB( datamat, split_info )
%
% =========================================================================
% MODULE_GNB: module that performs gaussian naive bayes analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_GNB( datamat, split_info )
%
%

% define 4 different task blocks for analysis:
block_cond1_sp1 = datamat(:,split_info.idx_cond1_sp1);
block_cond1_sp2 = datamat(:,split_info.idx_cond1_sp2);
block_cond2_sp1 = datamat(:,split_info.idx_cond2_sp1);
block_cond2_sp2 = datamat(:,split_info.idx_cond2_sp2);

% gnb analysis, under two different splitting structures...

% define split-half task design
design1 = [-ones( size(block_cond1_sp1,2), 1 ); ones( size(block_cond2_sp1,2), 1 )];
design2 = [-ones( size(block_cond1_sp2,2), 1 ); ones( size(block_cond2_sp2,2), 1 )];
% analysis
results = gnb_optimization( [block_cond1_sp1 block_cond2_sp1], [block_cond1_sp2 block_cond2_sp2], design1, design2, split_info.decision_model, split_info.spat_weight );

% Euclid. distance from (P=1,R=1)
DD = sqrt( (1-results.R).^2 + (1-results.P).^2 );
% select PC subspace that minimizes D(P,R)
[vd id]  = min(DD);

% [Record optimal statistics + eigenimages]
%
output.metrics.R    =  results.R(id);
output.metrics.P    =  results.P(id);
output.metrics.Dneg = -vd;
% optimal eigenimage
output.images  = results.eig(:,id);

% CV score timeseries, on unit-normed eigenimage
output.temp    = datamat'  * (output.images ./ sqrt(sum(output.images.^2)));
