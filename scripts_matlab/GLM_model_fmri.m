function [ output ] = GLM_model_fmri( dataVol, detrend_order, Xnoise, Xsignl, econflag, keepmean )


% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

% option to exclude analysis results -- defaut off
if( nargin < 5 ) econflag = []; end
% option to re-add mean during signal reconstruction -- default off
if( nargin < 6 ) keepmean =  0; end

[Nmeas Nsamp] = size(dataVol);

% computing the Legendre polynomial bases spanning this space:
Xdet   = get_legendre( detrend_order, Nsamp );
Xnoise = [Xdet Xnoise];
%
Xall   = [Xnoise Xsignl];
%
Nnoise   = size(Xnoise,2);
Nsignl   = size(Xsignl,2);
idxNoise = 1:Nnoise;
idxSignl = Nnoise+1:Nnoise+Nsignl;

% now perform regression, and get the matrix
BetaWeights = dataVol * Xall /( Xall'*Xall );
vol_estim   = BetaWeights * Xall';
noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
%
% filtered data matrices
if( keepmean == 0 )
    output.vol_resid   = dataVol - vol_estim;
    output.vol_denoi   = dataVol - noi_estim;
elseif( keepmean == 1 ) %% re-add mean (zero-order legendre poly.)
    output.vol_resid   = dataVol - vol_estim + ( BetaWeights(:,idxNoise(1))* Xall(:,idxNoise(1))' );
    output.vol_denoi   = dataVol - noi_estim + ( BetaWeights(:,idxNoise(1))* Xall(:,idxNoise(1))' );
end    

% for computational efficiency:
% check for 'econ' flag before doing Beta/Tstat estimation
if( ~strcmp(econflag,'econ') )

    % now, estimate t-statistics on signal
    residvar    = var(dataVol - vol_estim,0,2);
    % catch instances of zero-variance
    residvar(var(dataVol,0,2)==0) = 1E-10;
    %
    Xinvdiag    = diag(inv(Xall'*Xall));
    TmapWeights = zeros( Nmeas, Nnoise+Nsignl );
    %
    for( i=1:(Nnoise+Nsignl) )
        %
        TmapWeights(:,i) = BetaWeights(:,i) ./ sqrt( Xinvdiag(i) * residvar );
        %
    end

    % beta and t-stat maps
    output.Beta_signl  = BetaWeights(:,idxSignl);
    output.Tmap_signl  = TmapWeights(:,idxSignl);
    %
    output.Beta_noise  = BetaWeights(:,idxNoise);
    output.Tmap_noise  = TmapWeights(:,idxNoise);    
end

%%
function DL = random_polynomial( ord, N )


DL = [ones(N,1)];
for i = 0:ord
    ll = round((N-1)/2);
    A = 1./([1:ll].^(2*0.8-1));
    A = [0 A A(end:-1:1)];
    P = 2*pi*rand(1,ll)-pi;
    P = [0 P -P(end:-1:1)];
    P = exp(-sqrt(-1)*P);
    XX = (ifft([A.*P]));
    
    DL = [DL XX'];
end
%%

function DL = get_legendre( ord, N )
%
x  = linspace(-1,1,N)';
%
DL = [];
%
if (ord>=0)
    DL = [];
    for i = 0:ord
        X = generate_legendre(i,N);
        DL = [DL X];
    end
end


% if( ord>=0 ) d0 = ones(N,1);                                       DL=[DL d0]; end
% if( ord>=1 ) d1 = x;                                               DL=[DL d1]; end
% if( ord>=2 ) d2 = 0.5    *(3*x.^2   - 1);                          DL=[DL d2]; end
% if( ord>=3 ) d3 = 0.5    *(5*x.^3   - 3*x);                        DL=[DL d3]; end
% if( ord>=4 ) d4 = 0.125  *(35*x.^4  - 30*x.^2  + 3);               DL=[DL d4]; end
% if( ord>=5 ) d5 = 0.125  *(63*x.^5  - 70*x.^3  + 15*x);            DL=[DL d5]; end
% if( ord>=6 ) d6 = 0.0625 *(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);    DL=[DL d6]; end
% if( ord>=7 ) d7 = 0.0625 *(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x); DL=[DL d7]; end
%

function X = generate_legendre(ord,N)

x  = linspace(-1,1,N)';
pk = LegendrePoly(ord);
xp = bsxfun(@power,x,[ord:-1:0]);    
X = bsxfun(@times,xp,pk');    
X = sum(X,2);




