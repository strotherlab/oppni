function [ prob, sigdiff ] = friedman_test( respBlock, varargin )
% performs the Friedman test statistic on multiple-treatment blocks:
% 
%         [prob sigdiff] =  friedman_test( respBlock, (alpha) ) 
%
% wherein [ respBlock = 'treatments' x 'observations' ]
% e.g. each column is a block of observed treatments.
% 
% * Note that observations should be independant
% * Note also: test statistic given by Chi-Square approximation;
%   ideal performance is given for large samples and/or many treatments
% * sigdiff = critical difference at given alpha, for test
% * if alpha not specified, default is 0.05

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%


if( nargin == 2 )
    alpha = varargin{1};
else
    alpha = 0.05;
end

N = length( respBlock(1,:) ); % no. cols = no. observations
s = length( respBlock(:,1) ); % no. rows = no. treatments

% take observation blocks and rank treatments in each:
rankBlock = tiedrank( respBlock ); % rank treatment for each observed set
meanRank  = mean( rankBlock, 2 );          % compute mean rank, per treatment
grandMean = (s+1)/2;                       % compute grand mean on ranks:

% get friedman statistic:
% sum-of-squares difference between all treatment means and grand mean,
% with factor giving normal approximation for large samples:
Q = 12*N/( s*(s+1) ) * sum( (meanRank - grandMean).^2 );

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
