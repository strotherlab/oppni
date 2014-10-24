function result_set = lda_optimization_group ( data_class1, data_class2, drf, Resampling_Index )
%
% This is LDA for multiple data splits, when class structure is (possibly) different between splits
% * accounts for differences in model priors / design-matrix structure
%
% result_set = lda_optimization_group( data_class1, data_class2, drf, N_iters )
%

%% 1. concatenate
    N_iters = size(Resampling_Index,1);
    % number of voxels
    Nvox    = size(data_class1{1},1);
    % number of subjects (splits)
    Nsubj1  = length( data_class1 );
    Nsubj2  = length( data_class2 );
    % splits must be balanced
    if( Nsubj1 ~= Nsubj2 ) error('splits not balanced!');
    else                   Nsubj = Nsubj1;
    end
    
    % create fulldata matrices; mean-center 'em
    data_class1_full = [];
    data_class2_full = [];
    % removing subject means
    for(n = 1:Nsubj)
        %
        avg_temp = mean([data_class1{n} data_class2{n}],2);
        %
        data_class1{n} = data_class1{n} - repmat( avg_temp, [1 size(data_class1{n},2)] );
        data_class2{n} = data_class2{n} - repmat( avg_temp, [1 size(data_class2{n},2)] );        
        %
        data_class1_full = [data_class1_full data_class1{n}];
        data_class2_full = [data_class2_full data_class2{n}];
    end
    %
    N_class1_full = size(data_class1_full,2);
    N_class2_full = size(data_class2_full,2);
    %
    data_full = [data_class1_full data_class2_full];
    %
    design_full = [-ones( N_class1_full, 1 ); ones( N_class2_full, 1 )];
    
%% 2. initial data reduction factor + generate fulldata matrix

    % initial SVD -> run on full dataset
    [v s temp] = svd( data_full'*data_full ); s = sqrt(s);
    % now reduce dimensionality before further steps
    
    % first round feature selection
    drfPCs     = floor (size (data_full, 2) * drf);
    % catch any issues with DRF specification:
    %
    % if zero PCs, re-set to PC=1
    if    ( drfPCs < 1 )                              drfPCs = 1;
	% if more than 95% of PCs specified, adjust --> corrects for possible rank-defic
    elseif( drfPCs > floor(size(data_full, 2)*0.95) ) drfPCs = floor(size(data_full, 2)*0.95);
    end
    % define available PC range for splitting  
    K_max = min( [ floor((N_class1_full+N_class2_full)/2) round(drfPCs/2)] ) - 1;

    % now reduce PCA data
    v = v(:,1:drfPCs);
    s = s(1:drfPCs,1:drfPCs);
    % get image bases + PC-space representation
    img_bases      = data_full * v * inv(s);
    Z_full         = s*v'; % [pcs x time] matrix
    % SVD#2 on full data set (used for reference)
    Z_full     = Z_full - repmat (mean (Z_full, 2), [1 size(Z_full, 2)]);
    % record splits for later
    Z_full_class1 = Z_full(:, 1:N_class1_full );
    Z_full_class2 = Z_full(:, N_class1_full+1:end);
    %
    [u_full, s, v]   = svd (Z_full, 0);
    Z_reproj_full = s*v';
    
    % REFERENCE
    % reference set: full data (same # components!!)
  [lin_discr] = LDA_train_mini (Z_reproj_full, design_full, K_max);
    % --------------------------
    dx_full = u_full(:,1:K_max) * lin_discr;  % discriminant in PC space
    % aligning image based on CV scores:    
    CV_avg_sc1 = mean(Z_full_class1' * dx_full); % mean cv, class 1
    CV_avg_sc2 = mean(Z_full_class2' * dx_full); % mean cv, class 2
    CV_dif = CV_avg_sc2 - CV_avg_sc1;
    % flipem so that sign matches direction
    CV_dif = CV_dif.*sign( CV_dif );
    
    % SPLITDATA get the reduced features
    for(n = 1:Nsubj)
        %
        Z_class1{n} = img_bases' *data_class1{n};
        Z_class2{n} = img_bases' *data_class2{n};

        %% recording CV scores for output
        res_cv{n}.class1 = Z_class1{n}' * dx_full;
        res_cv{n}.class2 = Z_class2{n}' * dx_full;
    end
    
    % initializing data matrices...
    res_r      = zeros(K_max,N_iters);
    res_p      = zeros(K_max,N_iters);
    res_rSPMZ  = zeros(Nvox,K_max);
    
for( iter = 1:N_iters )
    
    %----%
    list1= Resampling_Index(iter,:);
    list2= setdiff(1:Nsubj,list1);

%% 3. define split1/2 data halves after initial feature reduction

    Z_sp1_class1 = [];
    Z_sp1_class2 = [];
    %
    for( i=1:length(list1) )
        %
        Z_sp1_class1 = [Z_sp1_class1 Z_class1{list1(i)}];
        Z_sp1_class2 = [Z_sp1_class2 Z_class2{list1(i)}];        
    end
    %
    Z_sp1 = [Z_sp1_class1 Z_sp1_class2];

    Z_sp2_class1 = [];
    Z_sp2_class2 = [];
    %
    for( i=1:length(list2) )
        %
        Z_sp2_class1 = [Z_sp2_class1 Z_class1{list2(i)}];
        Z_sp2_class2 = [Z_sp2_class2 Z_class2{list2(i)}];
    end
    %
    Z_sp2 = [Z_sp2_class1 Z_sp2_class2];

    % define split-half task design
    design_sp1 = [-ones( size(Z_sp1_class1,2), 1 ); ones( size(Z_sp1_class2,2), 1 )];
    design_sp2 = [-ones( size(Z_sp2_class1,2), 1 ); ones( size(Z_sp2_class2,2), 1 )];
    % centering PCs
    %
    avg_temp     = mean(Z_sp1,2);
    Z_sp1        = Z_sp1       -repmat( avg_temp, [1        size(Z_sp1,2)]);
    Z_sp1_class1 = Z_sp1_class1-repmat( avg_temp, [1 size(Z_sp1_class1,2)]);
    Z_sp1_class2 = Z_sp1_class2-repmat( avg_temp, [1 size(Z_sp1_class2,2)]);
    %
    avg_temp     = mean(Z_sp2,2);
    Z_sp2        = Z_sp2       -repmat( avg_temp, [1        size(Z_sp2,2)]);
    Z_sp2_class1 = Z_sp2_class1-repmat( avg_temp, [1 size(Z_sp2_class1,2)]);
    Z_sp2_class2 = Z_sp2_class2-repmat( avg_temp, [1 size(Z_sp2_class2,2)]);
    
    % number 'o' scans
    n_sp1_cl1 = size( Z_sp1_class1, 2 );
    n_sp1_cl2 = size( Z_sp1_class2, 2 );
    n_sp2_cl1 = size( Z_sp2_class1, 2 );
    n_sp2_cl2 = size( Z_sp2_class2, 2 );
    %
    N_sp1     = n_sp1_cl1 + n_sp1_cl2;
    N_sp2     = n_sp2_cl1 + n_sp2_cl2;
    
    % SVD#2 on split matrices
    [u_sp1, s, v] = svd (Z_sp1, 0); Z_reproj_sp1 = s*v';
    [u_sp2, s, v] = svd (Z_sp2, 0); Z_reproj_sp2 = s*v';
    
%%  -- 4 Compute linear discriminants
    
    % SPLIT #1
    % calculating SPM for the 1st set
  [lin_discr] = LDA_train_mini (Z_reproj_sp1, design_sp1, K_max);
    % --------------------------
     dx_sp1 = u_sp1(:,1:K_max) * lin_discr; % discriminant direction vector
    map_sp1 = img_bases * dx_sp1; 
    %
    CV_avg_sc1 = mean(Z_sp1_class1' * dx_sp1); % mean cv, class 1
    CV_avg_sc2 = mean(Z_sp1_class2' * dx_sp1); % mean cv, class 2
    CV_dif_sp1 = CV_avg_sc2 - CV_avg_sc1;
    % flipem so that sign matches direction
    CV_flip = sign( CV_dif_sp1 .* CV_dif );
    %
    map_sp1=map_sp1*diag( CV_flip );
        
    % SPLIT #2
    % calculating SPM for the 2nd set
  [lin_discr] = LDA_train_mini (Z_reproj_sp2, design_sp2, K_max);
    % --------------------------
     dx_sp2 = u_sp2(:,1:K_max) * lin_discr; % discriminant direction vector
    map_sp2 = img_bases * dx_sp2;     
    %
    CV_avg_sc1 = mean(Z_sp2_class1' * dx_sp2); % mean cv, class 1
    CV_avg_sc2 = mean(Z_sp2_class2' * dx_sp2); % mean cv, class 2
    CV_dif_sp2 = CV_avg_sc2 - CV_avg_sc1;
    % flipem so that sign matches direction
    CV_flip = sign( CV_dif_sp2 .* CV_dif );
    %
    map_sp2=map_sp2*diag( CV_flip );
    
%%  -- 3.3 Reproducibility and rSPM estimation
    
    clear spm;
    %
    for(k=1:K_max)
        [ res_r(k,iter), spm(:,k) ] = get_rSPM( map_sp1(:,k), map_sp2(:,k),1 );
    end
    %
    res_rSPMZ = res_rSPMZ + spm;
    
%%  -- 3.4 Prediction

    warning off;
    
    % (I) PREDICTION sp2 on sp1
    % mean CVscores, sp1 --- 1xK-size
    CV_sp1_avg_sc1 = mean(Z_sp1_class1' * dx_sp1); % mean cv, class 1
    CV_sp1_avg_sc2 = mean(Z_sp1_class2' * dx_sp1); % mean cv, class 2
        
    % scores: samples x K-size
    scores_sp2 = Z_sp2' * dx_sp1 ;    
    % unnormalized probabilities
    pp1_nopriors = exp(-((scores_sp2 - repmat(CV_sp1_avg_sc1, [size(scores_sp2,1) 1])).^2)./2);
    pp2_nopriors = exp(-((scores_sp2 - repmat(CV_sp1_avg_sc2, [size(scores_sp2,1) 1])).^2)./2);
    %
    pp1_priors   = pp1_nopriors .* (n_sp1_cl1/N_sp1);
    pp2_priors   = pp2_nopriors .* (n_sp1_cl2/N_sp1);
    % normalized
    pp1_priors_norm = pp1_priors./(pp1_priors+pp2_priors);
    pp2_priors_norm = pp2_priors./(pp1_priors+pp2_priors);
    %
    pp1_priors_norm(~isfinite(pp1_priors_norm)) = 0.50;
    pp2_priors_norm(~isfinite(pp2_priors_norm)) = 0.50;
    % probs -- sample x K-size
    sum_prob_sp2on1 = (  sum( pp1_priors_norm(design_sp2<0,:) ) + sum( pp2_priors_norm(design_sp2>0,:) )  )';

    % (I) PREDICTION sp1 on sp2
    % mean CVscores, sp2 --- 1xK-size
    CV_sp2_avg_sc1 = mean(Z_sp2_class1' * dx_sp2); % mean cv, class 1
    CV_sp2_avg_sc2 = mean(Z_sp2_class2' * dx_sp2); % mean cv, class 2
    % scores: samples x K-size
    scores_sp1 = Z_sp1' * dx_sp2 ;
    % unnormalized probabilities
    pp1_nopriors = exp(-((scores_sp1 - repmat(CV_sp2_avg_sc1, [size(scores_sp1,1) 1])).^2)./2);
    pp2_nopriors = exp(-((scores_sp1 - repmat(CV_sp2_avg_sc2, [size(scores_sp1,1) 1])).^2)./2);
    %
    pp1_priors   = pp1_nopriors .* (n_sp2_cl1/N_sp2);
    pp2_priors   = pp2_nopriors .* (n_sp2_cl2/N_sp2);
    % normalized
    pp1_priors_norm = pp1_priors./(pp1_priors+pp2_priors);
    pp2_priors_norm = pp2_priors./(pp1_priors+pp2_priors);
    %
    pp1_priors_norm(~isfinite(pp1_priors_norm)) = 0.50;
    pp2_priors_norm(~isfinite(pp2_priors_norm)) = 0.50;
    % probs -- sample x K-size
    sum_prob_sp1on2 = (  sum( pp1_priors_norm(design_sp1<0,:) ) + sum( pp2_priors_norm(design_sp1>0,:) )  )';

    warning on;
    
    res_p(:,iter) = (sum_prob_sp2on1 + sum_prob_sp1on2)./ (N_sp1 + N_sp2);
end
    
% divide by iterations
res_rSPMZ = res_rSPMZ./N_iters;

% now record results for output:
res_r = median( res_r, 2);
res_p = median( res_p, 2);
%
result_set.R    = res_r;
result_set.P    = res_p;
result_set.CV   = res_cv;
result_set.eig  = res_rSPMZ;

function [lin_discr] = LDA_train_mini (data, T_class, Range)
% performs compact linear discriminant:

warning off;

N = size (data, 2);
coord1 = data (:, T_class > 0);
coord0 = data (:, T_class < 0);
coord_norm = data - repmat (mean (data, 2), [1 size(data,2)]);
coord1_norm = coord1 - repmat (mean (coord1, 2), [1 size(coord1,2)]);
coord0_norm = coord0 - repmat (mean (coord0, 2), [1 size(coord0,2)]);

lin_discr = zeros( Range, Range );

for(k=1:Range)
    % -----
    W_ssp = coord1_norm(1:k,:)*coord1_norm(1:k,:)' + coord0_norm(1:k,:)*coord0_norm(1:k,:)';
    T_ssp = coord_norm(1:k,:)*coord_norm(1:k,:)';
    B_ssp = T_ssp - W_ssp;
    if rank(W_ssp)==length(W_ssp)
        [e, l, temp] = svd (inv (W_ssp) * B_ssp, 0);
        ee = e (:, 1);
        % -----
        lin_discr(1:k,k) = ee / (sqrt (ee' * W_ssp * ee / (N - 2)));
    else
        lin_discr(1:k,k) = 0;
    end
end

warning on;