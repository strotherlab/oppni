function result_set = gnb_optimization( data_trn, data_tst, design_trn, design_tst, decision_model, spatial_prior )
% 
%  Performs general 2-class GNB-NPAIRS analysis
%
%  res = gnb_optimization( data_trn, data_tst, design_trn, design_tst, decision_model, spatial_prior )
%

%% ============= (I) DEFINE LABEL INPUTS AND CLASS MATRICES ================

% truncate "empty" scans
keeptrn = find( design_trn ~= 0 );
keeptst = find( design_tst ~= 0 );
% now discard truncated (in design + data matrices)
design_trn = design_trn(keeptrn); data_trn = data_trn(:,keeptrn);
design_tst = design_tst(keeptst); data_tst = data_tst(:,keeptst);
% then check for zero variance voxels to drop out!!
allFinite = ( std(data_trn,0,2) > 0 ) .* ( std(data_tst,0,2) > 0 );

% define fMRI matrix subsets for NPAIRS train/test:
trnMat1 = data_trn(allFinite>0, design_trn == -1 );  % training (base) vols =first half of (base) indexed
trnMat2 = data_trn(allFinite>0, design_trn ==  1 );  % training (task) vols =first half of (task) indexed
tstMat1 = data_tst(allFinite>0, design_tst == -1 );  % test (base) vols =last half of (base) indexed
tstMat2 = data_tst(allFinite>0, design_tst ==  1 );  % test (task) vols =last half of (task) indexed

% ============= (2) COMPUTE STATISTICS + SENSITIVITY MAPS ON SPLITS ================

if( strcmp( decision_model, 'linear' ) )
    
    % * computed pooled variance statistic between the two classes

    % compute sample means and variances for train/test
    avg1_trn = mean(trnMat1,    2);   avg1_tst = mean(tstMat1,    2);
    avg2_trn = mean(trnMat2,    2);   avg2_tst = mean(tstMat2,    2);
    std1_trn = std (trnMat1, 0, 2);   std1_tst = std (tstMat1, 0, 2);
    std2_trn = std (trnMat2, 0, 2);   std2_tst = std (tstMat2, 0, 2);  
    % compute number of samples (image-vectors) for each group
    n1_trn = size (trnMat1, 2);       n1_tst = size (tstMat1,  2);
    n2_trn = size (trnMat2, 2);       n2_tst = size (tstMat2, 2);      

    % (R) SPATIAL REPRODUCIBILITY =========================================
    
    % calculate sensitivity map on the training set
    std_pooled_trn = sqrt ((n1_trn-1)*std1_trn.^2 + (n2_trn-1)*std2_trn.^2) / sqrt (n1_trn + n2_trn - 2);
    t_trn          = (avg2_trn - avg1_trn) ./ (std_pooled_trn .^2);
    % calculate sensitivity map on the test set
    std_pooled_tst = sqrt ((n1_tst-1)*std1_tst.^2 + (n2_tst-1)*std2_tst.^2) / sqrt (n1_tst + n2_tst - 2);
    t_tst          = (avg2_tst - avg1_tst) ./ (std_pooled_tst.^2);

    if( ~isempty( spatial_prior ) )
        %
        t_trn = t_trn.* spatial_prior(allFinite>0);
        t_tst = t_tst.* spatial_prior(allFinite>0);        
    end
    %
    % reproducible spm
    [ RR, rSPMZ ] = get_rSPM( t_trn, t_tst, 1 );
    eig = allFinite; eig(eig>0) = rSPMZ;

    % (P) PREDICTION ACCURACY =========================================    
    
    % preparatory: get the log-pp at each voxel / class ... done for full data matrix
    CORR_sum = 0;
    %
    for(w=1:n1_tst)
        logPP_t = - ((tstMat1(:,w) - avg1_trn).^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = - ((tstMat1(:,w) - avg2_trn).^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n2_tst)
        logPP_t = - ((tstMat2(:,w) - avg2_trn).^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = - ((tstMat2(:,w) - avg1_trn).^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n1_trn)
        logPP_t = - ((trnMat1(:,w) - avg1_tst).^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = - ((trnMat1(:,w) - avg2_tst).^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n2_trn)
        logPP_t = - ((trnMat2(:,w) - avg2_tst).^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = - ((trnMat2(:,w) - avg1_tst).^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    
    % classifier accuracy:
    PP = CORR_sum ./ (n1_trn+n2_trn+n1_tst+n2_tst);   
    
elseif( strcmp( decision_model, 'nonlinear' ) )
    
    % * based on the class-specific variance estimates

    % compute sample means and variances for train/test
    avg1_trn = mean(trnMat1,    2);   avg1_tst = mean(tstMat1,    2);
    avg2_trn = mean(trnMat2,    2);   avg2_tst = mean(tstMat2,    2);
    std1_trn = std (trnMat1, 0, 2);   std1_tst = std (tstMat1, 0, 2);
    std2_trn = std (trnMat2, 0, 2);   std2_tst = std (tstMat2, 0, 2);  
    % compute number of samples (image-vectors) for each group
    n1_trn = size (trnMat1, 2);       n1_tst = size (tstMat1,  2);
    n2_trn = size (trnMat2, 2);       n2_tst = size (tstMat2, 2);      

    % (R) SPATIAL REPRODUCIBILITY =========================================
    
    % calculate sensitivity map on the training set
    d_trn = 0;
    for(w=1:n1_trn) d_trn = d_trn + ( (trnMat1(:,w)-avg2_trn)./std2_trn.^2 - (trnMat1(:,w)-avg1_trn)./std1_trn.^2  ); end
    for(w=1:n2_trn) d_trn = d_trn + ( (trnMat2(:,w)-avg2_trn)./std2_trn.^2 - (trnMat2(:,w)-avg1_trn)./std1_trn.^2  );  end
    t_trn = d_trn ./ (n1_trn+n2_trn);
    % calculate sensitivity map on the test set
    d_tst = 0;
    for(w=1:n1_tst) d_tst = d_tst + ( (tstMat1(:,w)-avg2_tst)./std2_tst.^2 - (tstMat1(:,w)-avg1_tst)./std1_tst.^2  ); end
    for(w=1:n2_tst) d_tst = d_tst + ( (tstMat2(:,w)-avg2_tst)./std2_tst.^2 - (tstMat2(:,w)-avg1_tst)./std1_tst.^2  );  end
    t_tst = d_tst ./ (n1_tst+n2_tst);

    if( ~isempty( spatial_prior ) )
        %
        t_trn = t_trn.* spatial_prior(allFinite>0);
        t_tst = t_tst.* spatial_prior(allFinite>0);        
    end
    % reproducible spm
    [ RR, rSPMZ ] = get_rSPM( t_trn, t_tst, 1 );
    eig = allFinite; eig(eig>0) = rSPMZ;
    
    % (P) PREDICTION ACCURACY =========================================    
    
    % preparatory: get the log-pp at each voxel / class ... done for full data matrix
    CORR_sum = 0;
    %
    for(w=1:n1_tst)
        logPP_t = -log(std1_trn) - ((tstMat1(:,w) - avg1_trn).^2)./(2*std1_trn.^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = -log(std2_trn) - ((tstMat1(:,w) - avg2_trn).^2)./(2*std2_trn.^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n2_tst)
        logPP_t = -log(std2_trn) - ((tstMat2(:,w) - avg2_trn).^2)./(2*std2_trn.^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = -log(std1_trn) - ((tstMat2(:,w) - avg1_trn).^2)./(2*std1_trn.^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n1_trn)
        logPP_t = -log(std1_tst) - ((trnMat1(:,w) - avg1_tst).^2)./(2*std1_tst.^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = -log(std2_tst) - ((trnMat1(:,w) - avg2_tst).^2)./(2*std2_tst.^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end
    %
    for(w=1:n2_trn)
        logPP_t = -log(std2_tst) - ((trnMat2(:,w) - avg2_tst).^2)./(2*std2_tst.^2); logPP_t(~isfinite(logPP_t)) =0;
        logPP_f = -log(std1_tst) - ((trnMat2(:,w) - avg1_tst).^2)./(2*std1_tst.^2); logPP_f(~isfinite(logPP_f)) =0;
        % add to "correct" sum: weight=(1 for correct / 0.5 for undecided)
        CORR_sum  = CORR_sum + double( sum(logPP_t) > sum(logPP_f) );
    end

    % classifier accuracy:
    PP = CORR_sum ./ (n1_trn+n2_trn+n1_tst+n2_tst);   
end

% now record results for output:
result_set.R     = RR;
result_set.P     = PP;
result_set.eig   = eig;
