function out = module_ANOVA1( datamat, signif )

des =[]; kk = length(datamat);
Xmat=[];
for(i=1:kk) 
    des = [des; i*ones(size(datamat{i},2),1)]; 
    Xmat=[Xmat datamat{i}];
    nn(i) = size(datamat{i},2); 
end
% grand mean
GM = mean(Xmat,2);
% total sum of squares
SST= sum( bsxfun(@minus, Xmat,GM).^2, 2);
SSW= 0;
for(i=1:kk) 
    % treatment means
    TM(:,i) = mean(Xmat(:,des==i),2);
    % add to within-class variance
    SSW = SSW + sum(bsxfun(@minus,Xmat(:,des==i),TM(:,i)).^2,2);
end
% between-class variance
SSB= (bsxfun(@minus, TM,GM).^2)*nn(:);
% degrees of freedom
df_b = kk-1;
df_w = sum(nn) - kk;
% mean squared errors
MSB= SSB./df_b;
MSW= SSW./df_w;
% f statistic and likelihood
out.Fstat = MSB./MSW;
out.p_f = 1-fcdf(out.Fstat,(kk-1),sum(nn) - kk);
% effect size
out.omega2 = (SSB - (kk-1)*MSW)./(SST + MSW);

% thresholding
sig_omni = out.p_f <= signif;

% pairwise tests
pairs=[];
for(i=1:kk-1)
for(j=i+1:kk)
    pairs=[pairs; [i j]];
end
end
ncmp = size(pairs,1);

disp('Post hoc...');
% post-hoc testing
for(i=1:ncmp)
    F_prot = ((TM(:,pairs(i,1))-TM(:,pairs(i,2))).^2)  ./ (MSW .*(1/nn(pairs(i,1)) + 1/nn(pairs(i,2)))); %% protected ftest
    out.p_fpost(:,i) = (1-fcdf(F_prot,1,sum(nn) - kk)) ./ ncmp;    
    out.t_post(:,i)  = (TM(:,pairs(i,1))-TM(:,pairs(i,2))) ./ sqrt(MSW .*(1/nn(pairs(i,1)) + 1/nn(pairs(i,2)))); %% protected tstat
end
% zero out posthoc values
out.t_post = out.t_post.* double(out.p_fpost <= signif);
% drop voxels that did not reach omnibus significance
out.p_fpost(~sig_omni,:)=NaN;
out.t_post(~sig_omni,:) =NaN;
out.Fstat(~sig_omni) = NaN;
out.pairs = pairs;
