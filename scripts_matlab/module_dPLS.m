function output = module_dPLS( datamat, split_info )
%
% =========================================================================
% MODULE_DPLS: module that performs linear discriminant analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_dPLS( datamat, split_info )
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

% define 4 different task blocks for analysis:
block_cond1_sp1 = datamat(:,split_info.idx_cond1_sp1);
block_cond1_sp2 = datamat(:,split_info.idx_cond1_sp2);
block_cond2_sp1 = datamat(:,split_info.idx_cond2_sp1);
block_cond2_sp2 = datamat(:,split_info.idx_cond2_sp2);

% define split-half task designs:
design1 = [ ones( size(block_cond1_sp1,2), 1 ); -ones( size(block_cond2_sp1,2), 1 )];
design2 = [ ones( size(block_cond1_sp2,2), 1 ); -ones( size(block_cond2_sp2,2), 1 )];
% data splits
xx_sp1 = [block_cond1_sp1 block_cond2_sp1]; xx_sp1=bsxfun(@minus,xx_sp1,mean(xx_sp1,2));
xx_sp2 = [block_cond1_sp2 block_cond2_sp2]; xx_sp2=bsxfun(@minus,xx_sp2,mean(xx_sp2,2));

%%

% sp1
cl_av = [mean(xx_sp1(:,design1<0),2) mean(xx_sp1(:,design1>0),2)]; % class avg.
[u l] = svd( cl_av,'econ' ); % eigen-decomp
u_sp1 = u(:,1) * double( sign((cl_av(:,2)'*u(:,1)) - (cl_av(:,1)'*u(:,1))) ); % sign-flip

% predict(1)
% scores: samples x 1
scores_clav= cl_av'*u_sp1;
scores_sp2 = xx_sp2' * u_sp1;
% pooled variance
sig2 = ( (sum(design1<0)-1)*var(xx_sp1(:,design1<0)' * u_sp1) + (sum(design1>0)-1)*var(xx_sp1(:,design1>0)' * u_sp1) )./(length(design1)-2); %TEMP
% unnormalized probabilities
pp1_nopriors = exp(-((scores_sp2 - scores_clav(1)).^2)./(sig2*2));
pp2_nopriors = exp(-((scores_sp2 - scores_clav(2)).^2)./(sig2*2));
%
pp1_priors   = pp1_nopriors .* (sum(design1<0)/length(design1));
pp2_priors   = pp2_nopriors .* (sum(design1>0)/length(design1));
% normalized
pp1_priors_norm = pp1_priors./(pp1_priors+pp2_priors);
pp2_priors_norm = pp2_priors./(pp1_priors+pp2_priors);
%
pp1_priors_norm(~isfinite(pp1_priors_norm)) = 0.50;
pp2_priors_norm(~isfinite(pp2_priors_norm)) = 0.50;
% probs -- sample x K-size
sum_prob_sp2on1 = (  sum( pp1_priors_norm(design2<0,:) ) + sum( pp2_priors_norm(design2>0,:) )  )';
% simple classification accuracy
sum_correct_sp2on1 = (  sum( pp1_priors_norm(design2<0,:) >0.5 ) + sum( pp2_priors_norm(design2>0,:) >0.5 )  )';

% sp2
cl_av = [mean(xx_sp2(:,design2<0),2) mean(xx_sp2(:,design2>0),2)]; % class avg.
[u l] = svd( cl_av,'econ' ); % eigen-decomp
u_sp2 = u(:,1) * double( sign((cl_av(:,2)'*u(:,1)) - (cl_av(:,1)'*u(:,1))) ); % sign-flip

% predict(2)
% scores: samples x 1
scores_clav= cl_av'*u_sp2;
scores_sp1 = xx_sp1' * u_sp2;
% pooled variance
sig2 = ( (sum(design2<0)-1)*var(xx_sp2(:,design1<0)' * u_sp2) + (sum(design2>0)-1)*var(xx_sp2(:,design2>0)' * u_sp2) )./(length(design2)-2); %TEMP
% unnormalized probabilities
pp1_nopriors = exp(-((scores_sp1 - scores_clav(1)).^2)./(sig2*2));
pp2_nopriors = exp(-((scores_sp1 - scores_clav(2)).^2)./(sig2*2));
%
pp1_priors   = pp1_nopriors .* (sum(design2<0)/length(design2));
pp2_priors   = pp2_nopriors .* (sum(design2>0)/length(design2));
% normalized
pp1_priors_norm = pp1_priors./(pp1_priors+pp2_priors);
pp2_priors_norm = pp2_priors./(pp1_priors+pp2_priors);
%
pp1_priors_norm(~isfinite(pp1_priors_norm)) = 0.50;
pp2_priors_norm(~isfinite(pp2_priors_norm)) = 0.50;
% probs -- sample x K-size
sum_prob_sp1on2 = (  sum( pp1_priors_norm(design1<0,:) ) + sum( pp2_priors_norm(design1>0,:) )  )';
% simple classification accuracy
sum_correct_sp1on2 = (  sum( pp1_priors_norm(design1<0,:) >0.5 ) + sum( pp2_priors_norm(design1>0,:) >0.5 )  )';

% average posterior prob.
results.P = (sum_prob_sp2on1 + sum_prob_sp1on2)./ (length(design1) + length(design2));
% fractional accuracy
results.Acc = (sum_correct_sp2on1 + sum_correct_sp1on2)./ (length(design1) + length(design2));
% reproducibility
[results.R results.eig] = get_rSPM( u_sp1, u_sp2, 1 );

%%

% [Record optimal statistics + eigenimages]
%
output.metrics.R    =  results.R;
output.metrics.P    =  results.P;
output.metrics.Acc  =  results.Acc; % alt: fractional classif. accuracy
output.metrics.dPR  = -sqrt( (1-results.R).^2 + (1-results.P).^2 );
% optimal eigenimage
output.images  = results.eig;

% [CV scores]
%
% CV score timeseries, from reference eigenimage (NOT IN DPLS)
% output.temp.CV_ref = results.CV(:,id);
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
