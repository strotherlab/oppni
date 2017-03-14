function output = module_SVM( datamat, split_info )
%
% =========================================================================
% MODULE_SVM: module that performs linear discriminant analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_SVM( datamat, split_info )
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
xx_sp1 = [block_cond1_sp1 block_cond2_sp1]; xx_sp1=bsxfun(@minus,xx_sp1,mean(xx_sp1,2))';
xx_sp2 = [block_cond1_sp2 block_cond2_sp2]; xx_sp2=bsxfun(@minus,xx_sp2,mean(xx_sp2,2))';

%%
disp('half-scale');
% SVM discriminant analysis, under single-split structure
Cparam = 10.^[-2 0 2];
% 2 rounds of try/catch blocks
for(c=1:length(Cparam))
    
    % sp1
    try
        svm1 = svmtrain(xx_sp1,design1,'kernel_function', 'linear', 'boxconstraint', Cparam(c));%,'autoscale',false);
    catch
        disp('termin. issue (split-1), retrying'); 
        good=0; cnew=Cparam(c);
        while(good<1)
            disp('trying...');
            cnew = cnew * 10^(-0.5); % subtract -0.1 from log-scale
            try 
                svm1 = svmtrain(xx_sp1,design1,'kernel_function', 'linear', 'boxconstraint', cnew);%,'autoscale',false);
                good = 1;
            catch
                disp('failed');
                good = 0;
            end
        end
    end
    map1 = xx_sp1(svm1.SupportVectorIndices,:)' * svm1.Alpha;
    
    % sp2
    try    
        svm2 = svmtrain(xx_sp2,design2,'kernel_function', 'linear', 'boxconstraint', Cparam(c));%,'autoscale',false);
    catch
        disp('termin. issue (split-2), retrying'); 
        good=0; cnew=Cparam(c);
        while(good<1)
            disp('trying...');
            cnew = cnew * 10^(-0.5); % subtract -0.1 from log-scale
            try 
                svm2 = svmtrain(xx_sp2,design2,'kernel_function', 'linear', 'boxconstraint', cnew);%,'autoscale',false);
                good = 1;
            catch
                disp('failed');                
                good = 0;
            end
        end
    end
    map2 = xx_sp2(svm2.SupportVectorIndices,:)' * svm2.Alpha;
    
    % prediction
    acc(1,1) = sum( svmclassify( svm1, xx_sp2 ) == design2 )./length(design2);
    acc(2,1) = sum( svmclassify( svm2, xx_sp1 ) == design1 )./length(design1);
    results.P(c,1) = mean(acc);
    
    % reproducibility
    [results.R(c,1) results.eig(:,c)] = get_rSPM( map1, map2, 1 );
end

%%

% Euclid. distance from (P=1,R=1)
DD = sqrt( (1-results.R).^2 + (1-results.P).^2 );
% select PC subspace that minimizes D(P,R)
[vd id]  = min(DD);

% [Record optimal statistics + eigenimages]
%
output.metrics.R    =  results.R(id);
output.metrics.P    =  results.P(id);
output.metrics.dPR  = -vd;
% optimal eigenimage
output.images  = results.eig(:,id);

% [CV scores]
%
% CV score timeseries, from reference eigenimage (NOT IN SVM)
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
