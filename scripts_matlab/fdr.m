function [pcritSet threshMat] = fdr( dataMat, datType, qval, cv, dof )
% 
% False-Discovery Rate correction for multiple tests
% for matrix of some test statistic, get back binary matrix where
% 1=signif. at FDR thresh / 0=nonsignif
% 
%   [pcrit threshMat] = fdr( dataMat, datType, qval, cv, dof )
%
%   where datType -> 'p': p-value
%                    't': t-score -- requires a dof value
%                    'z': z-score
%
%         qval    -> level of FDR control (expected fdr rate)
% 
%         cv      -> constant of 0:[cv=1]  1:[cv= sum(1/i)]
%                    (1) applies under any joint distribution of pval
%                    (0) requires relative test indep & Gaussian noise with
%                        nonnegative correlation across tests
%         dof     -> degrees of freedom, use to assess significance of t-stat
%
%   * get back testing matrix threshMat & critical pvalue pcrit
%   * testing is currently 2-tail only!
%   * dof value is only relevant if you are using tstats
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

[Ntest Nk] = size(dataMat);

if    ( datType == 'p' )  probMat = dataMat;
elseif( datType == 't' )  probMat = 1-tcdf( abs(dataMat), dof );
elseif( datType == 'z' )  probMat = 1-normcdf( abs(dataMat) );    
end

threshMat = zeros( Ntest, Nk );

for( K=1:Nk )

    % (1) order pvalues smallest->largest
    pvect = sort( probMat(isfinite( probMat(:,K) ), K), 'ascend' );
    Ntest2= length(pvect);
    % (2) find highest index meeting limit criterion
    if(cv == 0) c_V = 1;
    else        c_V = log(Ntest2) + 0.5772 + 1/Ntest2; % approx soln.
    end
    % index vector
    indxd = (1:Ntest2)';
    % get highest index under adaptive threshold
    r = sum( pvect./indxd <= qval/(Ntest2*c_V) );

    if( r > 0 )
        
        % limiting p-value
        pcrit = pvect(r);
        % threshold matrix values
        if    ( datType == 'p' )
            threshMat(:,K) = double(dataMat(:,K)  <= pcrit);
        elseif( datType == 't' )
            vcrit = tinv(1-pcrit,dof);
            threshMat(:,K) = double(abs(dataMat(:,K)) >= vcrit);
        elseif( datType == 'z' )
            vcrit = norminv(1-pcrit);
            threshMat(:,K) = double(abs(dataMat(:,K)) >= vcrit);
        end
        
        pcritSet(K,1) = pcrit;
    else
        threshMat(:,K) = zeros(Ntest,1);
        pcritSet(K,1)  = NaN;
    end
    
end
