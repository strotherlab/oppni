function out = module_FRIEDMAN( datamat, alpha )

% wherein [ respBlock = 'treatments' x 'observations' ]
% e.g. each column is a block of observed treatments.

for(c=1:length(datamat)) XX(:,:,c) = datamat{c}; end
% tests x obs x treatments

[tests N s] = size( XX ); % (tests x obs x treats)

for(n=1:N)
    n,
    tmp = tiedrank(permute( XX(:,n,:), [3 1 2] )); %%treat x tests x obs
    XXrnk(:,:,n) = tmp;
end
XXrnk = permute(XXrnk,[2 3 1]);
meanRank = permute(mean( XXrnk, 2 ),[1 3 2] ); %% tests x treat -->avg across obs.
grandMean = (s+1)/2;                       % compute grand mean on ranks:

% get friedman statistic:
% sum-of-squares difference between all treatment means and grand mean,
% with factor giving normal approximation for large samples:
Q = 12*N/( s*(s+1) ) * sum( (meanRank - grandMean).^2, 2 );

% probability of difference >= Q approximated by the Chi-Square 
% cumulative distribution (for df = s-1):
prob   = 1 - chi2cdf( Q, s-1 );

% computing critical difference:
% (1) set critical inv-tscore
thresh = tinv(1-alpha/2, (N-1)*(s-1));
% (2) compute C-value, based on design
C = N*s*(s+1)*(s-1)/12;
% (3) now compute critical difference in ranks for sig. difference
sigdiff = (thresh/N) * sqrt( ( 2*N*C/( (N-1)*(s-1)) ) * (1 - Q/(N*(s-1))) );

%%

out.prob = prob;
out.sigdiff = sigdiff;
