function out = ANOVA1rep( datamat, design )
%
% . 1-way repeated measures ANOVA --> 1 experimental manipulation, multiple levels, 
%                                     + multiple samples over time,

% signif = 0.05;

if( size(design,2) ==1 )
    error('needs both subject and design factors');
elseif( size(design,2) ==2 )
    design = [ones(size(design,1),1) design];
elseif( size(design,2)  >2 )
    error('can only have 1 repeated factor (+1 random factor)');
end

% structured: (group)/(time)/(subj)
ix = unique( design(:,1) );%unique groups
jx = unique( design(:,2) );%unique timepoints
for( i=1:numel(ix)) % per group
   
    tmpmat = datamat( :, design(:,1)==ix(i) ); %all data belonging to this group
    tmpdes = design( design(:,1)==ix(i), 2:end); %design for all data in group (excl. group)
    
    for( j =1:numel(jx) ) % per timepoint
       
        tmpmat2 = tmpmat(:, tmpdes(:,1)==jx(j) ); %all data belonging to this timepoint
        tmpdes2 = tmpdes( tmpdes(:,1)==jx(j), end); %subject list belonging to this timepoint
        ord     = sortrows( [tmpdes2, (1:length(tmpdes2))'],1 ); % ordering numerically
        %
        tmpcell{i}(:,:,j) = tmpmat2(:,ord(:,2)); % store reordered data in cell matrix
    end
end
datamat = tmpcell; clear tmpcell;

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
  
    % group effect
    out.fstat_G    = FG;
    out.fstat_G_p  = pG;
    % time effect
    out.fstat_T    = FT;
    out.fstat_T_p  = pT;
    % interaction effect
    out.fstat_GT   = FGT;
    out.fstat_GT_p = pGT;
    
else
    pT  = 1 - fcdf(FT, T-1, (S-G)*(T-1));
    pSG = 1 - fcdf(FSG, S-G, (S-G)*(T-1));

    % time effect
    out.fstat_T    = FT;
    out.fstat_T_p  = pT;
end


out.testname = 'anova1rep';
