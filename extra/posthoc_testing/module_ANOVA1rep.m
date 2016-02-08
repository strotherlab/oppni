function out = module_ANOVA1rep( datamat, signif )

%number of groups, voxels timepoints
G = length(datamat); 
[Nvox, ~, T] = size(datamat{1});
%subjects per group 
S_g = zeros(G, 1);
for g=1:G
    S_g(g) = size(datamat{g},2);    
end
S = sum(S_g);

% overall mean (vox x 1)
y = 0;
for g=1:G
    y = y + sum(sum(datamat{g},2),3);
end
y = y ./ (S * T);

% allocate means
y_g = zeros(Nvox,G);
y_t = zeros(Nvox,T);
y_gt = zeros(Nvox,G,T);
y_gs = cell(G,1);
for g=1:G
    y_gs{g} = zeros(Nvox,S_g(g));
end

% group means
for g=1:G
    y_g(:,g) = sum(sum(datamat{g},2),3) / (S_g(g) * T);
end
% follow-up means
for t=1:T
    y_t(:,t) = 0;
    for g=1:G
        y_t(:,t) = y_t(:,t) + sum(datamat{g}(:,:,t),2);
    end
    y_t(:,t) = y_t(:,t) ./ S;
end
% group g and time t mean
for g=1:G
    for t=1:T
        y_gt(:,g,t) = sum(datamat{g}(:,:,t) ./ S_g(g),2);
    end
end
% subject s'th of group g mean
for g=1:G
    for s=1:S_g(g)
        y_gs{g}(:,s) = sum(datamat{g}(:,s,:),3) ./ T;
    end
end

% calculate the sum of squares
ssG  = 0;
ssSG = 0;
ssT  = 0;
ssGT = 0;
ssR  = 0;

for g=1:G
    for s=1:S_g(g)
        for t=1:T
            ssG  = ssG  + (y_g(:,g) - y).^2;
            ssSG = ssSG + (y_gs{g}(:,s) - y_g(:,g)).^2;
            ssT  = ssT  + (y_t(:,t) - y).^2;
            ssGT = ssGT + (y_gt(:,g,t) - y_g(:,g) - y_t(:,t) + y).^2;
            ssR  = ssR  + (datamat{g}(:,s,t) - y_gt(:,g,t) - y_gs{g}(:,s) + y_g(:,g)).^2;
        end
    end
end

% calculate means
if G > 1
    msG  = ssG  ./ (G-1);
    msGT = ssGT ./ ((G-1)*(T-1));
end
msSG = ssSG ./ (S-G);
msT  = ssT  ./ (T-1);
msR  = ssR  ./ ((S-G)*(T-1));

% calculate the F-statistics
if G > 1
    FG  = msG  ./ msSG;
    FGT = msGT ./ msR;
end
FT  = msT  ./ msR;
FSG = msSG ./ msR;

if G > 1 %% if multiple groups
    
    pG  = 1 - fcdf(FG, G-1, S-G);
    pT  = 1 - fcdf(FT, T-1, (S-G)*(T-1));
    pGT = 1 - fcdf(FGT, (G-1)*(T-1), (S-G)*(T-1));
    pSG = 1 - fcdf(FSG, S-G, (S-G)*(T-1));

    out.p_f   = [pT, pG, pGT];
    out.Fstat = [FT, FG, FGT];
else
    pT  = 1 - fcdf(FT, T-1, (S-G)*(T-1));
    pSG = 1 - fcdf(FSG, S-G, (S-G)*(T-1));

    out.p_f   = [pT];
    out.Fstat = [FT];
end
