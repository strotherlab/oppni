function [ output ] = rSVD_splithalf( X, Npc, Niters, threshold )
%
% RSVD_SPLITHALF:  this script takes a matrix of image vectors and
% performs Singular Value Decomposition (SVD) in a split-half resampling
% framework. This allows you to identify spatially reproducible components,
% and obtain Z-scored component maps.
%
% Syntax:
%         [ rSPM_match, R_match ] = rSVD_splithalf( X, Niters, threshold )
%
% Input:
%              X  : matrix, of (Nvoxel x Nsubject) column vectors
%          Niters : number of split-half resamples, used to estimate
%                   components (at least Niters=50 recommended)
%       threshold : applies False-Discovery Rate significance threshold 
%                  (threshold=0.05 is typical). This will threshold rSPM maps, 
%                   so that nonsignificant voxels have weight=0, and discards 
%                   rSPM components with <0.01% significant voxels
%                  *If you don't want this, just specify "threshold=[]" (eg make it empty)
%
% Output: an "output" structure with the following fields:
%
%        rSPM_match : set of Z-scored eigenimage maps (Nvoxel x Nsubject)
%      subj_weights : subject loadings on each rSPM pattern (normalized to unit variance)
%         Var_match : median percent variance accounted for by each component
%           R_match : associated median reproducibility values
%           R_distr : full distribution of reproducibility values for each PC
%          R_signif : significance of reproducibility value, measured as
%                     fraction of resamples with R>0. You will probably need 
%                     (Niters > 500) to get a stable, meaningful
%                     significance estimate here. This is generally more
%                     conservative than the FDR cutoff specified in "threshold".
%
% * NB: we match SVD components across using split-halves using a "greedy search" 
%       implementation of the Procrustes matching algorithm (a la PLSNPAIRS)
% * NB: we scale the eigenimages by their singular values, so that the
%       procrustes matching accounts for variance magnitude.
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
% ------------------------------------------------------------------------%% ------------------------------------------------------------------------%
%

 % get dimensions, and number of subjects in 1 split-half
 [Nvox Nsubj] = size(X);  Nsp = floor(Nsubj/2);
 % SVD on full data matrix for reference
 [Vref Sref] = svd( X'*X );
 % these are the reference eigenimages, we use these to match split-half components
 Uref = X * Vref;% variance scaled eigenimage
 % reduce to top "Nsp" components
 if( Npc > Nsp ) Npc = Nsp; end
 %
 Uref = Uref(:,1:Npc);
 
 % initialize variables
 rSPM_match = zeros( Nvox, Npc );
 repr_set   = zeros( Npc, Niters );
 svar       = zeros( Npc, Niters );
 
for( i=1:Niters )

   disp(i);
    
   % randomize subject-list & split in half
   prmlist = randperm(Nsubj);
   list1   = prmlist(1:Nsp);
   list2   = prmlist(Nsp+1:end);
   % create split-half data matrices
   X1 = X(:,list1);
   X2 = X(:,list2);
   
   % run split 1 SVD
   [V1 S1] = svd( X1'*X1 ); s11 = diag(S1)./trace(S1);
   % these are the reference eigenimages
   U1 = X1 * V1(:,1:Npc); % variance scaled eigenimage
   % run split 2 SVD
   [V2 S2] = svd( X2'*X2 ); s22 = diag(S2)./trace(S2);
   % these are the reference eigenimages
   U2 = X2 * V2(:,1:Npc); % variance scaled eigenimage
   
   % Procrustes matching transforms of spatial components, across splits
   % currently matching on "correlation":

%    [ Out1 ] = mini_procrust_ex( Uref, U1, 'rss' );
%    [ Out2 ] = mini_procrust_ex( Uref, U2, 'rss' );
%    % apply the spatial transforms
%    U1 = U1(:,Out1.index) * diag(Out1.flip);
%    U2 = U2(:,Out2.index) * diag(Out2.flip);
   
   for(k=1:Npc)
      U1(:,k) = U1(:,k) .* sign( corr( U1(:,k), Uref(:,k) ) );
      U2(:,k) = U2(:,k) .* sign( corr( U2(:,k), Uref(:,k) ) );       
   end
%    
%    % get mean %variance too
%    svar(:,i) = ( s11(Out1.index) + s22(Out2.index) )./2;
%    
   % estimate reproducible SPM for each component dimension:
   for(k=1:Npc)
       [ rep, rSPM ]   = get_rSPM_ex( U1(:,k), U2(:,k), 1 );
       rSPM_match(:,k) =  rSPM_match(:,k) + rSPM;
       repr_set(k,i)   =  rep;
   end   
   
   
end

% normalize Z-scores
rSPM_match = rSPM_match ./ Niters;
% median reproducibility
R_match = median(repr_set,2);
R_distr = repr_set';
% fraction of R estimates greater than zero
R_signif= sum( double(repr_set>0), 2 )./Niters;
% project data onto SPMs to get subject scores
subj_weights = X' * rSPM_match;
% normalize subject scores to have unit variance
subj_weights = bsxfun(@rdivide,subj_weights,sqrt(sum(subj_weights.^2)));
% get median %variance of each component
Var_match    = median( svar, 2 );

if( ~isempty( threshold ) )
    
    [pcrit thresh] = fdr_ex( rSPM_match, threshold, 0 );
    % apply thresholding
    rSPM_match     = rSPM_match .* thresh;
    % components with >.1% significant voxels
    signif_comps   = sum(thresh) > (0.001 * Nvox);

    % trim off non-significant components
    rSPM_match     = rSPM_match(:,signif_comps);
    R_match        = R_match( signif_comps );
    R_distr        = R_distr( :, signif_comps );
    R_signif       = R_signif( signif_comps );
    
    disp(strcat('found_',num2str( sum(signif_comps) ),'_significant PCs.'));
else
    disp('no thresholding...');
end

% record to output structure
output.rSPM_match = rSPM_match;
output.R_match    = R_match;
output.R_distr    = R_distr;
output.R_signif   = R_signif;
output.Var_match  = Var_match;
output.subj_weights = subj_weights;

%% ===================================================================== %%
% code snippets below are called in this function;
% they also exist independently in my code, minus the "_ex" bit

function [ Out ] = mini_procrust_ex( refVects, subVects, type )
%
% VERY simple version of procrustes - matches subVects to most appropriate
% refVect, in order to minimize global SumSquares difference / correlation criteria
%

% get dimensions from subspace-vectors
nVct     = size( subVects,2);
subV_idx = zeros(nVct,1);

if( strcmp( type , 'rss' ) )

    % get reference vector ordering, by decreasing variance
    ordRef  = sortrows( [ (1:nVct)' std( refVects )'], -2 );
    ordRef  = ordRef(:,1);

    % step through the reference vectors
    for( ir=1:nVct )
        % replicate out the reference
        tmpRef   = repmat( refVects(:,ir), [1 nVct] );
        % get the sum-of-squares difference from each reference vector (flipping to match by sign)
        SS(ir,:) = min( [sum((subVects - tmpRef).^2)', sum((-subVects - tmpRef).^2)'], [], 2 );
    end

    % we have sum-of-sqr difference matrix SS = ( nref x nsub )
    
    % step through reference vectors again (now by amt of explained var.)
    for(j=1:nVct)    
        % find the sub-vector index minimizing deviation from Ref
        [vs is] = min( SS(ordRef(j),:) );
        % for this Ref-vector, get index of best match SubVect
        subV_idx( ordRef(j) ) = is;
        % then "blank out" this option for all subsequent RefVects
        SS( :, is ) = max(SS(SS~=0)) + 1;
    end

    % reordered to match their appropriate RefVects
    subVects_reord = subVects(:,subV_idx);
    % now recontstruct what the sign was via index
    [vvv iii] = min([sum((-subVects_reord - refVects).^2)', sum((subVects_reord - refVects).^2)'], [], 2);
    % convert to actual sign
    flip= sign(iii-1.5);
    % output:
    % 
    Out.index  = subV_idx(:);
    Out.flip   = flip(:);

elseif( strcmp( type , 'corr' ) )
    
    ordRef  = (1:nVct)';

    % full correlations [ref x sub]
    CC = abs(corrcoef( [refVects, subVects] ));    
    CC = CC(1:nVct, nVct+1:end);
    
    remainRef = (1:nVct)'; ordRef = [];
    remainSub = (1:nVct)'; ordSub = [];
    
    CCtmp = CC;
    
    for( i=1:nVct )
        
        % get max correlation of ea. ref
        [vMax iMax] = max( CCtmp,[], 2 );
        % find Ref with highest match
        [vOpt iRef] = max( vMax     );
        % also get "sub" index:
              iSub  = iMax(iRef);
        
        ordRef = [ordRef remainRef( iRef )];
        remainRef( iRef ) = [];
        %
        ordSub = [ordSub remainSub( iSub )];
        remainSub( iSub ) = [];
        %
        CCtmp(iRef,:) = [];
        CCtmp(:,iSub) = [];
    end
    
    CCnew = corrcoef( [refVects(:,ordRef) subVects(:,ordSub)] );
    flip  = sign( diag( CCnew(1:nVct,nVct+1:end) ) );
    
    resort = sortrows([ordRef(:) ordSub(:) flip(:)], 1);

    Out.index = resort(:,2);
    Out.flip   = resort(:,3);    
end
%%
function [ rep, rSPM ] = get_rSPM_ex( vect1, vect2, keepMean )
%
%    [ rep, rSPM ] = get_rSPM( vect1, vect2, keepMean )
%
%  * This script takes in 2 vectors, returns reproducibility (rep)
%    and reproducible SPM (rSPM)
%
%  * keepMean: reinsert the mean offset present in vectors
%

cc  = corrcoef (vect1, vect2);
rep = cc (1, 2);

%(1) getting the mean offsets (normed by SD)
normedMean1 = mean(vect1)./std(vect1);
normedMean2 = mean(vect2)./std(vect2);
%    and rotating means into signal/noise axes
sigMean = (normedMean1 + normedMean2)/sqrt(2);
%noiMean = (normedMean1 - normedMean2)/sqrt(2);
%(2) getting  signal/noise axis projections of (zscored) betamaps
sigProj = ( zscore(vect1) + zscore(vect2) ) / sqrt(2);
noiProj = ( zscore(vect1) - zscore(vect2) ) / sqrt(2);
% noise-axis SD
noiStd = std(noiProj);
%(3) norming by noise SD:
%     ...getting the (re-normed) mean offsets
sigMean = sigMean./noiStd;
%noiMean = noiMean./noiStd; 
%  getting the normed signal/noise projection maps
sigProj = sigProj ./ noiStd;
%noiProj = noiProj ./ noiStd;

% Produce the rSPM:
if    ( keepMean == 1 )   rSPM = sigProj + sigMean;
elseif( keepMean == 0 )   rSPM = sigProj;
end
%%
function [pcritSet threshMat] = fdr_ex( dataMat, qval, cv )
% 
% trimmed-down FDR code
%

[Ntest Nk] = size(dataMat);
probMat    = 1-normcdf( abs(dataMat) );    
threshMat  = zeros( Ntest, Nk );

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
        % threshold matrix values based on prob.        
        threshMat(:,K) = double(probMat(:,K)  <= pcrit);
        % critical p-values        
        pcritSet(K,1)  = pcrit;
    else
        threshMat(:,K) = zeros(Ntest,1);
        pcritSet(K,1)  = NaN;
    end
end
