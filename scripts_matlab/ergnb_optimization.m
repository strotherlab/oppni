function result_set = ergnb_optimization( data_all, spatial_prior )
% 
%  Performs event related GNB-NPAIRS analysis
%
%  res = gnb_optimization( data_trn, data_tst, design_trn, design_tst, decision_model, spatial_prior )
%

%% ============= (I) DEFINE LABEL INPUTS AND CLASS MATRICES ================

% dimensions
[Nvox Nclass Nfull] = size( data_all );
% mean-center the individual data windows
data_all = data_all - repmat( mean(data_all,2), [1 Nclass 1] );
% number of window samples per split
Nhalf = ceil(Nfull/2);

if( Nfull < 4 ) error('At least 4 data splits are required!'); end
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

% initializations
RR    = zeros(TopIter,1);
rSPMZ = zeros( Nvox, TopIter );
pSPM  = zeros( Nvox, TopIter );
PP    = zeros(TopIter,1);

% ============= (2) COMPUTE STATISTICS + SENSITIVITY MAPS ON SPLITS ================
for( w=1:TopIter ) %% for each resampling split...
    
    if( Nfull <= 6 ) list = check(w,:);
    else             list = randperm(Nfull);
    end
    
    % define splits
    DSET_trn = data_all(:,:,list(1:Nhalf))    ;
    DSET_tst = data_all(:,:,list(Nhalf+1:end));

    % normalize to "global" standard deviation of 1 for prediction estimates
    % first split
    tmp = sqrt( sum(sum( DSET_trn.^2, 3 ),2) ./ ( size(DSET_trn,2)*size(DSET_trn,3) - 1 ) );
    tmp = repmat( tmp, [1 Nclass] );
    %
    for(n=1:size(DSET_trn,3))  DSET_trn(:,:,n) = DSET_trn (:,:,n)./ tmp;  end
    % second split
    tmp = sqrt( sum(sum( DSET_tst.^2, 3 ),2) ./ ( size(DSET_tst,2)*size(DSET_tst,3) - 1 ) );
    tmp = repmat( tmp, [1 Nclass] );
    %
    for(n=1:size(DSET_tst,3))  DSET_tst(:,:,n) = DSET_tst (:,:,n)./ tmp;  end
    % correct for ill-posed voxels
    DSET_trn(~isfinite(DSET_trn))=eps;
    DSET_tst(~isfinite(DSET_tst))=eps;
    
    % within-class variance, pooled

    % initialized
    Wmat_trn = zeros(Nvox,1);
    Wmat_tst = zeros(Nvox,1);
    %
    for(c=1:Nclass) % adding variance from each class
        Wmat_trn = Wmat_trn + var(  permute( DSET_trn(:,c,:), [1 3 2] ), 0,2  );
        Wmat_tst = Wmat_tst + var(  permute( DSET_tst(:,c,:), [1 3 2] ), 0,2  );
    end
    % averaging
    Wmat_trn = Wmat_trn./Nclass;
    Wmat_tst = Wmat_tst./Nclass; 
         
    % between-class variance

    % class means
    Mmat_trn = mean(DSET_trn,3);
    Mmat_tst = mean(DSET_tst,3);
    % variance across means
    Bmat_trn = var( Mmat_trn,  0,2 );
    Bmat_tst = var( Mmat_tst,  0,2 );
    
    % --------------------- sensitivity mapping
    
    s_trn = Bmat_trn ./ Wmat_trn;
    s_tst = Bmat_tst ./ Wmat_tst;

    if( ~isempty( spatial_prior ) )
        %
        s_trn = s_trn.* spatial_prior;
        s_tst = s_tst.* spatial_prior;
    end
    %
    % reproducible spm
    [ RR(w,1) rSPMZ(:,w) ] = get_rSPM( s_trn, s_tst, 1 );
    
    % ---------------------- predictions, global
    
    CORR_sum = 0;

    for(n=1:size(DSET_tst,3) )
    for(c=1:Nclass) 

        % testing volume
        x = DSET_tst(:,c,n);    
        % get prob for each training class
        for(t=1:Nclass)
            %
            logPP(:,t) = -((x - Mmat_trn(:,t)).^2) ./ (2.*Wmat_trn);
        end
        % correct for illposed
        logPP( ~isfinite(sum(logPP,2)),: ) = 0;

        % identify the most probable class
        [vc ic]  = max( sum(logPP,1) );
        % increment if trueclass is most likely
        CORR_sum = CORR_sum + double(ic==c);
    end
    end

    for(n=1:size(DSET_trn,3) )
    for(c=1:Nclass) 

        % testing volume
        x = DSET_trn(:,c,n);    
        % get prob for each training class
        for(t=1:Nclass)
            %
            logPP(:,t) = -((x - Mmat_tst(:,t)).^2) ./ (2.*Wmat_tst);
        end
        % correct for illposed
        logPP( ~isfinite(sum(logPP,2)),: ) = 0;

        % identify the most probable class
        [vc ic]  = max( sum(logPP,1) );
        % increment if trueclass is most likely
        CORR_sum = CORR_sum + double(ic==c);
    end
    end

    % classifier accuracy:
    PP(w,1) = CORR_sum ./ (Nfull*Nclass);   
    
    % ---------------------- predictions, localer
    
    CORR_sum =             0; % global prediction
    PMAP_sum = zeros(Nvox,1); % predictive SPM

    for(n=1:size(DSET_tst,3) )
        %
        x = DSET_tst(:,:,n);
        d = (x-Mmat_trn).^2;
        %
        CORR_sum = CORR_sum + sum(  sqrt(mean(d,1))  );
        PMAP_sum = PMAP_sum + sqrt(mean(d,2));
    end
    for(n=1:size(DSET_trn,3) )
        %
        x = DSET_trn(:,:,n);
        d = (x-Mmat_tst).^2;
        %
        CORR_sum = CORR_sum + sum(  sqrt(mean(d,1))  );
        PMAP_sum = PMAP_sum + sqrt(mean(d,2));
    end
        
    % classifier accuracy:
    TT(w,1)   = CORR_sum ./ (Nfull*Nclass);
    pSPM(:,w) = PMAP_sum ./ (Nfull*Nclass);
end

% re-center TT and pSPM to make it a positive predictive term
TT   = (sqrt(2) -   TT)./sqrt(2);
pSPM = (sqrt(2) - pSPM)./sqrt(2);

% now record results for output:
result_set.R_global    =      mean(RR);
result_set.P_class     =      mean(PP);
result_set.P_rms       =      mean(TT);
% spatial maps: global sensitivity, prediction
result_set.sens_global = mean(rSPMZ,2);
result_set.pmap        = mean(pSPM, 2);

% -------- single component representation starts -------- %

    % average time-window, for normalized vox-variance
    DSET_all = data_all ./ repmat( std(data_all,0,2), [1 size(data_all,2) 1] );
    DSET_all = mean(DSET_all,3);
    DSET_all(~isfinite(DSET_all))=eps;
    % rescale by activation sensitivity
    DSET_all = DSET_all .* repmat( result_set.sens_global, [1 size(DSET_all,2)] );
    % svd on the map
    [u s v] = svd(DSET_all'*DSET_all); v1comp = v(:,1);
    %
    clear Rnorm DSET_all;

    % re-initializations
    RR    = zeros(TopIter,1);
    rSPMZ = zeros( Nvox, TopIter );
    % now testing reproducibility of spatial representation
    for( w=1:TopIter ) %% for each resampling split...

        if( Nfull <= 6 ) list = check(w,:);
        else             list = randperm(Nfull);
        end

        % define splits
        DSET_trn = data_all(:,:,list(1:Nhalf))    ;
        DSET_tst = data_all(:,:,list(Nhalf+1:end));

        DSET_trn = DSET_trn ./ repmat( std(DSET_trn,0,2), [1 size(DSET_trn,2) 1] );
        DSET_tst = DSET_tst ./ repmat( std(DSET_tst,0,2), [1 size(DSET_tst,2) 1] );
        %
        DSET_trn = mean(DSET_trn ,3);     DSET_trn(~isfinite(DSET_trn))=eps;
        DSET_tst = mean(DSET_tst ,3);     DSET_tst(~isfinite(DSET_tst))=eps;
        %
        DSET_trn = DSET_trn .* repmat( result_set.sens_global, [1 size(DSET_trn,2)] );
        DSET_tst = DSET_tst .* repmat( result_set.sens_global, [1 size(DSET_tst,2)] );
        
        sa = DSET_trn * v1comp;
        sb = DSET_tst * v1comp;
        
        % reproducible spm
        [ RR(w,1) rSPMZ(:,w) ] = get_rSPM( sa, sb, 1 );
    end

    % record in outputs
    result_set.sens_hrf1 = mean(rSPMZ,2);
    result_set.R_hrf1    = mean(RR);
    result_set.HRF1      = v1comp;

% -------- single component representation done -------- %
