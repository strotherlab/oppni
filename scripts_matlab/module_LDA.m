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
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

if( ~isfield(split_info,'drf') || isempty(split_info.drf) )
    disp('LDA uses default data reduction drf=0.5');
    split_info.drf = 0.5;
end

% define 4 different task blocks for analysis:
block_cond1_sp1 = datamat(:,split_info.idx_cond1_sp1);
block_cond1_sp2 = datamat(:,split_info.idx_cond1_sp2);
block_cond2_sp1 = datamat(:,split_info.idx_cond2_sp1);
block_cond2_sp2 = datamat(:,split_info.idx_cond2_sp2);

% define split-half task designs:
design1 = [ ones( size(block_cond1_sp1,2), 1 ); -ones( size(block_cond2_sp1,2), 1 )];
design2 = [ ones( size(block_cond1_sp2,2), 1 ); -ones( size(block_cond2_sp2,2), 1 )];
disp('Checkpoint Kb')
% linear discriminant analysis, under single-split structure
results = lda_optimization( [block_cond1_sp1 block_cond2_sp1], [block_cond1_sp2 block_cond2_sp2], design1,design2, split_info.drf );
disp('Checkpoint Kc')
% Euclid. distance from (P=1,R=1)
DD = sqrt( (1-results.R).^2 + (1-results.P).^2 );
% select PC subspace that minimizes D(P,R)
[vd id]  = min(DD);
disp('Checkpoint L')
% [Record optimal statistics + eigenimages]
%
disp('Andrew Add In Temp');

output.temp = results.temp;
output.metrics.R    =  results.R(id);
output.metrics.P    =  results.P(id);
output.metrics.Acc  =  results.Acc(id); % alt: fractional classif. accuracy
output.metrics.dPR  = -vd;
% optimal eigenimage
output.images  = results.eig(:,id);

% [CV scores]
%
% CV score timeseries, from reference eigenimage
output.temp.CV_ref = results.CV(:,id);
% CV score timeseries, on unit-normed rSPM eigenimage
output.temp.CV_alt = datamat' * (output.images ./ sqrt(sum(output.images.^2)));
disp('Checkpoint M')
% [Fractional Variance Explained by eigenimage basis]
%
% the scaled projection
svect = output.temp.CV_alt;
disp('Checkpoint X')
% and normed projection
uvect = svect ./ sum(svect.^2);
% get back out the scaling factor relative to normed eig
svar = var( svect ) ./ var( uvect );
disp('Checkpoint Y')
% total data variance
tvar = trace( datamat'*datamat );
% fraction
output.temp.CV_alt_varfract = svar ./ tvar;
disp('Checkpoint Z')
