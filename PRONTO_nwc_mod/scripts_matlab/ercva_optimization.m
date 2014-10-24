function [out] = ercva_optimization( DATA, List, opt_crit )
%
% ERCVA_OPTIMIZATION:  hybrid event-related Canonical Variates Analysis (CVA).
% This is a multivariate analysis model that takes in a set of time-windowed
% data blocks, then performs CVA analysis to identify event-related components.
%
% Syntax:
%            out = ercva_optimization( DATA, List, opt_crit )
%
% Input:
%            DATA : 3D matrix of time-windowed data blocks, of size (voxel x timelag x splits)
%            List : 2D integer vector, e.g. List=[#full-data PCs, #split-half PCs]
%        opt_crit : string, specifying Procrustes matching criterion, 
%                   'temp'=temporal(CV timeseries), 'spat'=spatial(eigenimages)
%
% Output:  this model estimates a component subspace, for each split-half PC dimensionality k=1...K
%          where number of components C = min([k, (timelag-1)])
%
%          out.P_global : K x 1 vector of "global" prediction values estimated over all CV dimensions
%                         (Bayes posterior probability), one per PC subspace size k
%          out.P_per_cv : K x 1 cell array. Each entry is a vector of length C, giving individual component 
%                         prediction (normalized RMS accuracy)
%          out.R        : K x 1 cell array. Each entry is a vector of length C, giving individual component 
%                         SPM reproducibilities
%          out.Tset     : K x 1 cell array. Each entry is a (timelag x C) matrix of average HRF timecourse 
%                         basis vectors expressed by each of the components
%          out.SPMs     : K x 1 cell array. Each entry is a (voxel x C) matrix of Z-scored eigenimage SPMs,
%                         one per CV dimension
%

PCA1 = List(1); % Number of PCs to keep after avg; 0=skip
PCA2 = List(2); % range 2-PCAsp for split-halves
% data dimensions
[Nvox Nlag Nfull] = size(DATA);
% mean-center the individual data windows
DATA = DATA - repmat( mean(DATA,2), [1 Nlag 1] );
for(n=1:Nfull) DATA(:,:,n) = DATA(:,:,n)./mean(std(DATA(:,:,n),0,2)); end
% DATA = DATA ./ repmat( std(DATA,0,2), [1 Nlag 1] );
DATA(~isfinite(DATA)) = eps;
% number of window samples per split
Nhalf = ceil(Nfull/2);

if( Nfull < 4 ) error('At least 4 data splits are required!'); end
%%
%--------------------------------------------------
% Splitting Structure Design: (random for Nsplit>6)
switch Nfull
               
    case 4    
          TopIter = 3;
          check   = [1 2, 3 4; 1 3, 2 4; 1 4, 2 3];
    case 5    
          TopIter = 10;
          check   = [1 2 3, 4 5; 1 2 4, 3 5; 1 2 5, 3 4; 1 3 4, 2 5; 1 3 5, 2 4; 
                     1 4 5, 2 3; 2 3 4, 1 5; 2 3 5, 1 4; 2 4 5, 1 3; 3 4 5, 1 2 ];
    case 6    
          TopIter = 10;
          check   = [1 2 3, 4 5 6; 1 2 4, 3 5 6; 1 2 5, 3 4 6; 1 2 6, 3 4 5; 1 3 4, 2 5 6; 
                     1 3 5, 2 4 6; 1 3 6, 2 4 5; 1 4 5, 2 3 6; 1 4 6, 2 3 5; 1 5 6, 2 3 4];       
    otherwise
          TopIter = 30;
end
%--------------------------------------------------

% get fulldata PCA
DATA2D = reshape( DATA, Nvox,[],1 );
[Vfull Sfull temp] = svd( DATA2D'*DATA2D );
% ----------
% reduced #fulldata pcs, according to DRF
 Vfull = Vfull(:,1:PCA1);
 Sfull = sqrt(Sfull(1:PCA1,1:PCA1));
 Ufull = DATA2D * Vfull *inv(Sfull);
 clear Sfull Vfull;
% ----------
QSET = zeros( PCA1, Nlag, Nfull );
% ----------
for(i=1:Nfull)
    % get pc-space representation of each window
    QSET(:,:,i) = Ufull'*DATA(:,:,i);
end

% initialize design matrixes
DESIGN_ful  = repmat( eye( Nlag ), [Nfull       1] ); %% fulldata
DESIGN_spl1 = repmat( eye( Nlag ), [Nhalf       1] ); %% split#1
DESIGN_spl2 = repmat( eye( Nlag ), [Nfull-Nhalf 1] ); %% split#2
%
PClist = 1:PCA2;
% initializing result structures, in cell array format
SPMset = cell(PCA2,1); %% SPMs
TCRset = cell(PCA2,1); %% HRF timecourses
REPset = cell(PCA2,1); %% reproducibility values
RMSset = cell(PCA2,1); %% rms temporal prediction
% initialize global CV prediction
Pprob  = zeros( PCA2, TopIter );

% initialize full-data representations
IMG0   = cell(PCA2,1);
TCR0avg= cell(PCA2,1);

% pca computations on fulldata:
QSET_ful         = reshape( QSET, PCA1,[],1 );  %%%  
QSET_ful         = QSET_ful - repmat( mean(QSET_ful,2), [1 size(QSET_ful,2)] );
[Uful Sful temp] = svd( QSET_ful ); %%%
Q2SET_ful        = Uful'*QSET_ful; %%%
    
for( w=1:TopIter ) %% for each resampling split...
    
    % define splitting structure
    if( Nfull <= 6 ) list = check(w,:);
    else             list = randperm(Nfull);
    end
    
    QSET_trn = reshape( QSET(:,:,list(1:Nhalf)),     PCA1,[],1 );
    QSET_tst = reshape( QSET(:,:,list(Nhalf+1:end)), PCA1,[],1 );
    %
    QSET_trn = QSET_trn - repmat( mean(QSET_trn,2), [1 size(QSET_trn,2)] );
    QSET_tst = QSET_tst - repmat( mean(QSET_tst,2), [1 size(QSET_tst,2)] );
    %
    [Utrn Strn temp] = svd( QSET_trn );
    [Utst Stst temp] = svd( QSET_tst );
    
    % split-half PCA space representations
    Q2SET_trn = Utrn'*QSET_trn;
    Q2SET_tst = Utst'*QSET_tst;
    % project other split onto PCA space, for prediction
    Q2SET_trnOtst = Utst'*QSET_trn;
    Q2SET_tstOtrn = Utrn'*QSET_tst;
    
    for(q=PClist) %% for each PC subspace dimensionality...
        
        % max. number of CV dimensions
        Ncap = min( [q Nlag-1] );
        % preinitialize
        TCR1avg = zeros( Nlag, Ncap );
        TCR2avg = zeros( Nlag, Ncap );
        
        % initializing full-data results - only on first resampling iter.
        if( w==1 )
            
            % full data transformation
            warning off;
            [lin_discr_ful,b]   = canoncorr(Q2SET_ful(1:q,:)',DESIGN_ful); 
            warning on;
            IMG0{q} = Ufull * Uful (:, 1:q) * lin_discr_ful(:,1:Ncap);
            TCR0    = Q2SET_ful(1:q,:)' * lin_discr_ful(:,1:Ncap);
            %
            for(i=1:Ncap)
                TCR0avg{q}(:,i) = mean(reshape( TCR0(:,i), Nlag, Nfull ),2);
            end
            
            %% -- initialize cell contents at iteration #1 -- %%
            SPMset{q} = zeros( Nvox, Ncap );
            TCRset{q} = zeros( Nlag, Ncap );
            REPset{q} = zeros( Ncap, TopIter );
            RMSset{q} = zeros( Ncap, TopIter );
        end
        
        %% CVA Analysis:
        warning off;
        % ---
        % Running CVA on split#1
        [lin_discr_trn,b]   = canoncorr(Q2SET_trn(1:q,:)',DESIGN_spl1); 
         lin_discr_trn = lin_discr_trn(:,1:Ncap);
        IMG1 = Ufull * Utrn (:, 1:q) * lin_discr_trn;
        TCR1 = Q2SET_trn(1:q,:)' * lin_discr_trn;
        % Running CVA on split#2
        [lin_discr_tst,b]   = canoncorr(Q2SET_tst(1:q,:)',DESIGN_spl2); 
         lin_discr_tst = lin_discr_tst(:,1:Ncap);
        IMG2 = Ufull * Utst (:, 1:q) * lin_discr_tst;
        TCR2 = Q2SET_tst(1:q,:)' * lin_discr_tst;
        % ---
        warning on;
        
        % compute the mean CV timeseries
        for(i=1:Ncap)
            TCR1avg(:,i) = mean(reshape( TCR1(:,i), Nlag, Nhalf       ),2);
            TCR2avg(:,i) = mean(reshape( TCR2(:,i), Nlag, Nfull-Nhalf ),2);
        end

        % Matching CV dimensions between splits, uses a restricted procrustes transformation
        %
        if( strcmp( opt_crit, 'temp' ) ) %option 1: match on timeseries

            [ Out ]      = mini_procrust_ex( TCR0avg{q}, TCR1avg, 'rss'  );
            IMG1         = IMG1(:,Out.index) * diag( Out.flip );
            TCR1         = TCR1(:,Out.index) * diag( Out.flip );
            TCR1avg      = TCR1avg(:,Out.index) * diag( Out.flip );
            lin_discr_trn= lin_discr_trn(:,Out.index) * diag(Out.flip);
            %
            [ Out ]      = mini_procrust_ex( TCR0avg{q}, TCR2avg, 'rss'  );
            IMG2         = IMG2(:,Out.index) * diag( Out.flip );
            TCR2         = TCR2(:,Out.index) * diag( Out.flip );
            TCR2avg      = TCR2avg(:,Out.index) * diag( Out.flip );
            lin_discr_tst= lin_discr_tst(:,Out.index) * diag(Out.flip);
            
        elseif( strcmp( opt_crit, 'spat' ) ) % optiont 2: match on spatial brain patterns
            
            [ Out ]      = mini_procrust_ex( IMG0{q}, IMG1, 'corr' ); 
            IMG1         = IMG1(:,Out.index) * diag( Out.flip );
            TCR1avg      = TCR1avg(:,Out.index) * diag( Out.flip );
            lin_discr_trn= lin_discr_trn(:,Out.index) * diag(Out.flip);
            %
            [ Out ]      = mini_procrust_ex( IMG0{q}, IMG2, 'corr' );
            IMG2         = IMG2(:,Out.index) * diag( Out.flip );
            TCR2avg      = TCR2avg(:,Out.index) * diag( Out.flip );
            lin_discr_tst= lin_discr_tst(:,Out.index) * diag(Out.flip);
        end

%% -----------------------------------------------------------------------

        % REPRODUCIBILITY
        
        for(i=1:Ncap)
            %
            [rep SPM] = get_rSPM( IMG1(:,i), IMG2(:,i), 1 );
            %
            REPset{q}(i,w) = rep;
            SPMset{q}(:,i) = SPMset{q}(:,i) + SPM;
        end
        
        % averaged CV timecourse, across splits and resamples
        TCRset{q} = TCRset{q} + (TCR1avg + TCR2avg)./2;
                
%% -----------------------------------------------------------------------

        % PREDICTION

        %% CV space transforms, preparing for Prediction

        % training data projections
        CV_trnOtrn = TCR1; % sampl x dim
        CV_tstOtst = TCR2; % sampl x dim
        % test data projections
        CV_trnOtst = Q2SET_trnOtst(1:q,:)' * lin_discr_tst; % sampl x dim
        CV_tstOtrn = Q2SET_tstOtrn(1:q,:)' * lin_discr_trn; % sampl x dim
        % training class averages
        CLavg_trn = TCR1avg; % class x dim (row-vects)
        CLavg_tst = TCR2avg;        
        
        % (FIRST PREDICTION METRIC): Global prediction across all CVs, on classification probability

        % initialize true-class probability matrix
        truprob = zeros( size( CV_trnOtst,1 )+size( CV_tstOtrn,1 ), 1 );
        kq = 0;

        % (1) Predict Test from Train data
        for(i=1:size( CV_tstOtrn,1 )) % for each (multi-dimensional) data point

            kq=kq+1; %%increment

            % rows of DIFF are difference-vector from point-i to each class' average
            %DIFF = ( repmat( CV_tstOtrn(i,:), [Nlag 1] ) - CLavg_trn );
            DIFF = bsxfun(@minus,CV_tstOtrn(i,:),CLavg_trn);
            % EXPO vector: get squared distance from class-average to point-i, convert to Gaussian prob
            EXPO = exp( -0.5* diag( DIFF*DIFF' ) );
            % catch for all prob ~0; renorm total prob=1
            if( sum(EXPO) == 0 ) EXPO = (1./Nlag) * ones(Nlag,1);
            else                 EXPO = EXPO ./ sum(EXPO);
            end
            % record true class' probability
            truprob(kq) = EXPO( DESIGN_spl2(i,:) >0 );
        end
        % (2) Predict Train from Test data
        for(i=1:size( CV_trnOtst,1 )) % for each (multi-dimensional) data point

            kq=kq+1; %%increment

            % rows of DIFF are difference-vector from point-i to each class' average
            %DIFF = ( repmat( CV_trnOtst(i,:), [Nlag 1] ) - CLavg_tst ); 
            DIFF = bsxfun(@minus,CV_trnOtst(i,:),CLavg_tst);
            % EXPO vector: get squared distance from class-average to point-i, convert to Gaussian prob
            EXPO = exp( -0.5* diag( DIFF*DIFF' ) );
            % catch for all prob ~0; renorm total prob=1
            if( sum(EXPO) == 0 ) EXPO = (1./Nlag) * ones(Nlag,1);
            else                 EXPO = EXPO ./ sum(EXPO);
            end
            % record true class' probability
            truprob(kq) = EXPO( DESIGN_spl1(i,:) >0 );
        end

        % record mean posterior-probability
        Pprob(q,w) = mean( truprob(:) );

        % (SECOND PREDICTION METRIC): Per-component prediction, on RMS fit of estimated HRF 

        % renormalize to constrain RMS distribution
        Z_trnOtst = zscore( CV_trnOtst );
        Z_tstOtrn = zscore( CV_tstOtrn );
        Z_trnOtrn = zscore( CV_trnOtrn );
        Z_tstOtst = zscore( CV_tstOtst );
    
        % iterate through each CV dimension
        for(i=1:Ncap)

            % (1) Predict Test from Training
            %
            % training mean, across samples
            Z0 = mean( reshape( Z_trnOtrn(:,i), Nlag, [] ), 2);
            % test-data reshaped into (lags x splits) matrix
            Zx = reshape( Z_tstOtrn(:,i), Nlag, [] );
            % RMS distance
            DRMS1 = sqrt( (Zx - repmat(Z0, [1 size(Zx,2)])).^2 );

            % (2) Predict Training from Test
            %
            % test mean, across samples
            Z0 = mean( reshape( Z_tstOtst(:,i), Nlag, [] ), 2);
            % train-data reshaped into (lags x splits) matrix
            Zx = reshape( Z_trnOtst(:,i), Nlag, [] );
            % RMS distance
            DRMS2 = sqrt( (Zx - repmat(Z0, [1 size(Zx,2)])).^2 );

            % record RMS error in predictions
            RMSset{q}(i,w) = mean( [DRMS1(:); DRMS2(:)] );
        end
    end
end

% Computing average across resamples, for different metrics
for(q=PClist)
    % divide by # resamples
    SPMset{q} = SPMset{q}./w;
    TCRset{q} = TCRset{q}./w;
    REPset{q} = median( REPset{q}, 2 );
    
    % renormalize RMS for output:
    tmp_rms   = median( RMSset{q}, 2 );
    % convert to prediction (higher=better), normalized to range from 0-1
    RMSset{q} = (sqrt(2) - tmp_rms)./sqrt(2);
end

% Record results to output structure
out.R        = REPset;
out.P_global = median(Pprob,2); 
out.P_per_cv = RMSset;
% ---------------
out.Tset = TCRset;
out.SPMs = SPMset;

%%
function [ Out ] = mini_procrust_ex( refVects, subVects, type )
%
% VERY simple version of procrustes - matches subVects to most appropriate
% refVect, in order to minimize global SumSquares difference criteria
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
    CC = abs( corr( refVects, subVects ) );    
    
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
        
        ordSub = [ordSub remainSub( iSub )];
        remainSub( iSub ) = [];
        
        CCtmp(iRef,:) = [];
        CCtmp(:,iSub) = [];        
    end
    
    CCnew = corr( refVects(:,ordRef), subVects(:,ordSub) );
    flip  = sign( diag( CCnew ) );
    
    resort = sortrows([ordRef(:) ordSub(:) flip(:)], 1);

    Out.index = resort(:,2);
    Out.flip   = resort(:,3);    
end
