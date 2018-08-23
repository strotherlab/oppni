function result_set = lda_optimization ( data_sp1, data_sp2, design_sp1, design_sp2, drf )

% This is LDA for 2 data splits, when class structure is different between splits
% * accounts for differences in model priors / design-matrix structure
%
% result_set = lda_optimization2( data_sp1, data_sp2, design_sp1, design_sp2, drf )
%
% design vectors = single-string vectors:
%     -1= condition 1 // +1= condition2
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

%% 1. drop censored scans + concatenate

Nvox    = size(data_sp1,1);
% initial selection of non-transition scans
keep_sp1   = find( design_sp1 ~= 0 );
keep_sp2   = find( design_sp2 ~= 0 );
% combine data + ideal-file matrices
data_sp1   =  data_sp1(:,keep_sp1);
data_sp2   =  data_sp2(:,keep_sp2);
%
design_sp1 =  design_sp1(keep_sp1);
design_sp2 =  design_sp2(keep_sp2);
%
N_sp1      =    length(design_sp1);
n_sp1_cl1  = sum( design_sp1 < 0 );
n_sp1_cl2  = sum( design_sp1 > 0 );
N_sp2      =    length(design_sp2);
n_sp2_cl1  = sum( design_sp2 > 0 );
n_sp2_cl2  = sum( design_sp2 < 0 );

% concatenated data matrix + design
data_full   = [data_sp1       data_sp2];
design_full = [design_sp1;  design_sp2];


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
K_max = min( [N_sp1 N_sp2 round(drfPCs/2)] ) - 1;

% now reduce PCA data
v = v(:,1:drfPCs);
s = s(1:drfPCs,1:drfPCs);
% get image bases + PC-space representation
img_bases      = data_full * v * inv(s);
Z_full         = s*v'; % [pcs x time] matrix

% SVD on full data set (used for reference)
Z_full         = bsxfun(@minus,Z_full,mean (Z_full, 2));
[u_full, s, v] = svd (Z_full, 0);
Z_reproj_full  = s*v';


%% 3. define split1/2 data halves after initial feature reduction

% re-split data once feature selection is complete
Z_sp1 = Z_full(:,      1:N_sp1);
Z_sp2 = Z_full(:,N_sp1+1:end   );
% centering PCs
Z_sp1  = bsxfun(@minus,Z_sp1,mean (Z_sp1, 2));
Z_sp2  = bsxfun(@minus,Z_sp2,mean (Z_sp2, 2));
% get class-splits for laters
Z_full_class1 = Z_full(:,design_full<0); 
Z_full_class2 = Z_full(:,design_full>0);
Z_sp1_class1  =  Z_sp1(:,design_sp1 <0); 
Z_sp1_class2  =  Z_sp1(:,design_sp1 >0);
Z_sp2_class1  =  Z_sp2(:,design_sp2 <0); 
Z_sp2_class2  =  Z_sp2(:,design_sp2 >0);

% svd on split matrices
[u_sp1, s, v] = svd (Z_sp1, 0); Z_reproj_sp1 = s*v';
[u_sp2, s, v] = svd (Z_sp2, 0); Z_reproj_sp2 = s*v';

%%  -- 4 Compute linear discriminants

% REFERENCE
% reference set: full data (same # components!!)
[lin_discr,K_max] = LDA_train_mini (Z_reproj_full, design_full, K_max);
[lin_discr1,K_max] = LDA_train_mini (Z_reproj_sp1, design_sp1, K_max);
[lin_discr2,K_max] = LDA_train_mini (Z_reproj_sp2, design_sp2, K_max);

lin_discr = lin_discr(1:K_max,1:K_max);
lin_discr1 = lin_discr1(1:K_max,1:K_max);
lin_discr2 = lin_discr2(1:K_max,1:K_max);

res_r      = zeros(K_max,1);
res_rSPMZ  = zeros(Nvox,K_max);

% --------------------------
dx_full = u_full(:,1:K_max) * lin_discr;  % discriminant in PC space
% aligning image based on CV scores:
CV_avg_sc1 = mean(Z_full_class1' * dx_full); % mean cv, class 1
CV_avg_sc2 = mean(Z_full_class2' * dx_full); % mean cv, class 2
CV_dif = CV_avg_sc2 - CV_avg_sc1;

% flipem so that sign matches direction
CV_dif = CV_dif.*sign( CV_dif );
% recording CV scores::
res_cv = Z_full' * dx_full;


% SPLIT #1
% calculating SPM for the 1st set
% --------------------------
dx_sp1 = u_sp1(:,1:K_max) * lin_discr1; % discriminant direction vector
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
% --------------------------
dx_sp2 = u_sp2(:,1:K_max) * lin_discr2; % discriminant direction vector
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

for(k=1:K_max)
    [ res_r(k,1), res_rSPMZ(:,k) ] = get_rSPM( map_sp1(:,k), map_sp2(:,k),1 );
end

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

%Added Aug 2018, to stop rounding determining result 
pp1_nopriors(find(pp1_nopriors <= 1*10^-16)) = 1*10^-16;
pp2_nopriors(find(pp2_nopriors <= 1*10^-16)) = 1*10^-16;

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
% simple classification accuracy
sum_correct_sp2on1 = (  sum( pp1_priors_norm(design_sp2<0,:) >0.5 ) + sum( pp2_priors_norm(design_sp2>0,:) >0.5 )  )';

if( min(pp1_nopriors(:))<1E-6 )
   fprintf('OPPNI LDA debug: prediction of split2 on split1\n posterior prob. of class 1 (unnormalized) is <10-6\n Value is: %f\n', ( min(pp1_nopriors(:))));
end
if( min(pp2_nopriors(:))<1E-6 )
   fprintf('OPPNI LDA debug: prediction of split2 on split1\n posterior prob. of class 2 (unnormalized) is <10-6\n Value is: %f\n', ( min(pp2_nopriors(:))));
end

% (I) PREDICTION sp1 on sp2
% mean CVscores, sp2 --- 1xK-size
CV_sp2_avg_sc1 = mean(Z_sp2_class1' * dx_sp2); % mean cv, class 1
CV_sp2_avg_sc2 = mean(Z_sp2_class2' * dx_sp2); % mean cv, class 2
% scores: samples x K-size
scores_sp1 = Z_sp1' * dx_sp2 ;
% unnormalized probabilities
pp1_nopriors = exp(-((scores_sp1 - repmat(CV_sp2_avg_sc1, [size(scores_sp1,1) 1])).^2)./2);
pp2_nopriors = exp(-((scores_sp1 - repmat(CV_sp2_avg_sc2, [size(scores_sp1,1) 1])).^2)./2);

%Added Aug 2018, to stop rounding determining result 
pp1_nopriors(find(pp1_nopriors <= 1*10^-16)) = 1*10^-16;
pp2_nopriors(find(pp2_nopriors <= 1*10^-16)) = 1*10^-16;

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
% simple classification accuracy
sum_correct_sp1on2 = (  sum( pp1_priors_norm(design_sp1<0,:) >0.5 ) + sum( pp2_priors_norm(design_sp1>0,:) >0.5 )  )';

if( min(pp1_nopriors(:))<1E-6 )
   fprintf('OPPNI LDA debug: prediction of split1 on split2\n posterior prob. of class 1 (unnormalized) is <10-6\n Value is: %f\n', ( min(pp1_nopriors(:))));
end
if( min(pp2_nopriors(:))<1E-6 )
   fprintf('OPPNI LDA debug: prediction of split1 on split2\n posterior prob. of class 2 (unnormalized) is <10-6\n Value is: %f\n', ( min(pp2_nopriors(:))));
end

warning on;

% average posterior prob.
res_p = (sum_prob_sp2on1 + sum_prob_sp1on2)./ (N_sp1 + N_sp2);
% fractional accuracy
res_acc = (sum_correct_sp2on1 + sum_correct_sp1on2)./ (N_sp1 + N_sp2);

% drop PC#1 for single-subject
res_r=res_r(2:end);
res_p=res_p(2:end);
res_acc=res_acc(2:end);
res_rSPMZ =res_rSPMZ(:,2:end);

% now record results for output:
result_set.R    = res_r;
result_set.P    = res_p;
result_set.Acc  = res_acc;
result_set.CV   = res_cv;
result_set.eig  = res_rSPMZ;

% quick optimization on D metric:
D      = sqrt( (1-res_r).^2 + (1-res_p).^2 );
[dv i] = min( D );
%
result_set.OPT.RPD = [res_r(i) res_p(i) dv];
result_set.OPT.CV  = res_cv(:,i);
result_set.OPT.eig = res_rSPMZ(:,i);


function [lin_discr,Kmax_out] = LDA_train_mini (data, T_class, Range)
% performs compact linear discriminant:

warning off;

N = size (data, 2);
coord1 = data (:, T_class > 0);
coord0 = data (:, T_class < 0);
coord_norm  = bsxfun(@minus,data,  mean (data, 2));
coord1_norm = bsxfun(@minus,coord1,mean (coord1, 2));
coord0_norm = bsxfun(@minus,coord0,mean (coord0, 2));

lin_discr = zeros( Range, Range );
Kmax_out = Range;
for(k=1:Range)
    % -----
    W_ssp = coord1_norm(1:k,:)*coord1_norm(1:k,:)' + coord0_norm(1:k,:)*coord0_norm(1:k,:)';
    T_ssp = coord_norm(1:k,:)*coord_norm(1:k,:)';
    B_ssp = T_ssp - W_ssp;
    if rank(W_ssp)==length(W_ssp)  
        
        %Octave Adjustment
        % [e, l, temp] = svd (inv (W_ssp) * B_ssp, 0);
        [e, l, temp] = svd (W_ssp \ B_ssp);
 
        % Add a double check to confirm that the svd was done correctly
        if any(any(((W_ssp \ B_ssp) - (e*l*temp')) >= 0.001))
            error('SVD done incorrectly in this version of MATLAB or OCTAVE is off by >1, check lda_optimization.');
        end
      
        ee = e (:, 1);
       
        % -----
        lin_discr(1:k,k) = ee / (sqrt (ee' * W_ssp * ee / (N - 2)));
     else
        Kmax_out = k-1;
        break;
    end
end

warning on;
