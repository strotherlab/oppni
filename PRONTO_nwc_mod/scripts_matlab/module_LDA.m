function output = module_LDA( datamat, split_info )
%
% =========================================================================
% MODULE_LDA: module that performs linear discriminant analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_LDA( datamat, split_info )
%
%

% define 4 different task blocks for analysis:
block_cond1_sp1 = datamat(:,split_info.idx_cond1_sp1);
block_cond1_sp2 = datamat(:,split_info.idx_cond1_sp2);
block_cond2_sp1 = datamat(:,split_info.idx_cond2_sp1);
block_cond2_sp2 = datamat(:,split_info.idx_cond2_sp2);

% define split-half task designs:
design1 = [-ones( size(block_cond1_sp1,2), 1 ); ones( size(block_cond2_sp1,2), 1 )];
design2 = [-ones( size(block_cond1_sp2,2), 1 ); ones( size(block_cond2_sp2,2), 1 )];
% linear discriminant analysis, under single-split structure
results = lda_optimization( [block_cond1_sp1 block_cond2_sp1], [block_cond1_sp2 block_cond2_sp2], design1,design2, split_info.drf );

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

% [CV scores]
%
% CV score timeseries, from reference eigenimage
output.temp.CV_ref = results.CV(:,id);
% CV score timeseries, on unit-normed rSPM eigenimage
output.temp.CV_alt = datamat' * (output.images ./ sqrt(sum(output.images.^2)));

% [Fractional Variance Explained by eigenimage basis]
%
% the scaled projection
svect = output.temp.CV_alt;
% and normed projection
uvect = svect ./ sum(svect.^2);
% get back out the scaling factor relative to normed eig
svar = var( svect ) ./ var( uvect );
% total data variance
tvar = trace( datamat'*datamat );
% fraction
output.temp.CV_alt_varfract = svar ./ tvar;
