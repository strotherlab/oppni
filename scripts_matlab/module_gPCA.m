function output = module_gPCA( datamat, split_info )
%
% =========================================================================
% MODULE_GSVD: measures generalization SVD on data
% =========================================================================
%
%   Syntax:
%           output = module_gPCA( datamat, split_info )
%
%
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

% matrix dimensions
[Nvox Ntime] = size(datamat);
% split into two equal halves
datamat_sp1 = datamat(:,1:ceil(Ntime/2));
datamat_sp2 = datamat(:,ceil(Ntime/2)+1:end);
% mean-center splits, to do correlation later
datamat_sp1 = bsxfun(@minus, datamat_sp1, mean(datamat_sp1,2) );
datamat_sp2 = bsxfun(@minus, datamat_sp2, mean(datamat_sp2,2) );

if( ~isfield(split_info,'num_PCs') || isempty(split_info.num_PCs) )
    disp('gPCA uses default 50% of total PCs');
    split_info.num_PCs = floor( 0.5 * Ntime );
end

% svd on splits
[u1 l1 v1]  = svd( datamat_sp1,'econ');
[u2 l2 v2]  = svd( datamat_sp2,'econ');
% total variance of splits
Tvar1 = trace(l1.^2);
Tvar2 = trace(l2.^2);
% fulldata reference
[u0 l0 v0] = svd( [datamat_sp1 datamat_sp2],'econ');
%
u0 = u0(:,1:split_info.num_PCs);
l0 = l0(1:split_info.num_PCs,1:split_info.num_PCs);
v0 = v0(:,1:split_info.num_PCs);

%%

for(kk = split_info.num_PCs)

    % IF generalization error
    %p_sse(1) = evalSSE_MOD(datamat_sp2,u1(:,1:kk),diag(l1(1:kk,1:kk)));
    %p_sse(2) = evalSSE_MOD(datamat_sp1,u2(:,1:kk),diag(l2(1:kk,1:kk)));
    % Gaussian generalization error
    p_gge(1) =  PGauss_MOD(datamat_sp2,u1(:,1:kk),diag(l1(1:kk,1:kk)).^2, kk, Tvar1);
    p_gge(2) =  PGauss_MOD(datamat_sp1,u2(:,1:kk),diag(l2(1:kk,1:kk)).^2, kk, Tvar2);
    % spatial reproducibility
    rv = RV_coef( u1(:,1:kk), u2(:,1:kk) ); 
end

output.metrics.R = rv;
output.metrics.P = -mean(p_gge);
output.temp = v0 * l0;

U1set = u1(:,1:split_info.num_PCs);
U2set = u2(:,1:split_info.num_PCs);
oo    = mini_procrust_ex( u0, U1set, 'corr' );
U1set = U1set(:,oo.index) * diag(oo.flip);
oo    = mini_procrust_ex( u0, U2set, 'corr' );
U2set = U2set(:,oo.index) * diag(oo.flip);

for(kk = 1:split_info.num_PCs)
    [ output.metrics.R_percomp(kk,1), output.images(:,kk) ] = get_rSPM( U1set(:,kk), U2set(:,kk), 1 );
end


%%
function SSE=evalSSE_MOD(X,A,c_diag)
% General evaluation of the indirect fitting objective
% SSE=\sum_k ||X{k}*X{k}'-A diag(C_k) H'*H diag(C_k) A'||_F^2

c_diag = c_diag(:); %%MODICATE

AtA=A'*A;
HtH=eye(length(c_diag));
SST=sum(sum((X'*X).^2));    
HCkAtX= (diag(c_diag))*A'*X; %(H*diag(c_diag))*A'*X;
SSXM=sum(sum(HCkAtX.^2));            
Q=AtA*diag(c_diag)*HtH*diag(c_diag);
RR2=sum(sum(Q.*Q'));

SSE=SST-2*SSXM+RR2;

%%
function PP=PGauss_MOD(X,Us,L_diag,K,Tvar)
%
% Gaussian likelihood

% basic params
Ndat   =  size(X,2); % #datapoints
Nvox   =  size(X,1); %number of voxels

% renorm var
L_diag = L_diag ./ Ndat;
% variance
sigma2 = ( (Tvar./Ndat) - sum(L_diag)) ./ ( Nvox-K ); % residual variance
% test variance
proj   = Us'*X;
varj   = sum( proj.^2, 2 );
% -log posterior prob.
PP =        0.5*log(2*pi) ...
          + 0.5*sum(log(L_diag)) ...
          + 0.5*( Nvox-K )*log(sigma2) ...
          + ( sum(sum(X.^2)) - sum( varj .* (L_diag-sigma2)./L_diag ) )./(2*sigma2*Ndat);

      
%% --------------------------------------------------
function rv = RV_coef( X, Y )
%
% rv = RV_coef( X, Y )
%

X = bsxfun(@minus,X,mean(X));
Y = bsxfun(@minus,Y,mean(Y));

SIG_xx = X'*X;
SIG_yy = Y'*Y;
SIG_xy = X'*Y;
SIG_yx = Y'*X;

covv = trace( SIG_xy * SIG_yx );
vavx = trace( SIG_xx * SIG_xx );
vavy = trace( SIG_yy * SIG_yy );

rv = covv ./ sqrt( vavx * vavy );

%% --------------------------------------------------
function [ Out ] = mini_procrust_ex( refVects, subVects, type )
%
% VERY simple version of procrustes - matches subVects to most appropriate
% refVect, in order to minimize global SumSquares difference criteria
%

% get dimensions from subspace-vectors
nVct     = size( subVects,2);
subV_idx = zeros(nVct,1);

if( strcmp( type , 'rss' ) )

    % get reference vector ordering, by decreasing variance
    ordRef  = sortrows( [ (1:nVct)' std( refVects )'], -2 );
    ordRef  = ordRef(:,1);

    % step through the reference vectors
    for( ir=1:nVct )
        % replicate out the reference
        tmpRef   = repmat( refVects(:,ir), [1 nVct] );
        % get the sum-of-squares difference from each reference vector (flipping to match by sign)
        SS(ir,:) = min( [sum((subVects - tmpRef).^2)', sum((-subVects - tmpRef).^2)'], [], 2 );
    end

    % we have sum-of-sqr difference matrix SS = ( nref x nsub )
    
    % step through reference vectors again (now by amt of explained var.)
    for(j=1:nVct)    
        % find the sub-vector index minimizing deviation from Ref
        [vs is] = min( SS(ordRef(j),:) );
        % for this Ref-vector, get index of best match SubVect
        subV_idx( ordRef(j) ) = is;
        % then "blank out" this option for all subsequent RefVects
        SS( :, is ) = max(SS(SS~=0)) + 1;
    end

    % reordered to match their appropriate RefVects
    subVects_reord = subVects(:,subV_idx);
    % now recontstruct what the sign was via index
    [vvv iii] = min([sum((-subVects_reord - refVects).^2)', sum((subVects_reord - refVects).^2)'], [], 2);
    % convert to actual sign
    flip= sign(iii-1.5);
    % output:
    % 
    Out.index  = subV_idx(:);
    Out.flip   = flip(:);

elseif( strcmp( type , 'corr' ) )
    
    ordRef  = (1:nVct)';

    % full correlations [ref x sub]
    CC = abs( corr( refVects, subVects ) );    
    
    remainRef = (1:nVct)'; ordRef = [];
    remainSub = (1:nVct)'; ordSub = [];
    
    CCtmp = CC;
    
    for( i=1:nVct )
        
        % get max correlation of ea. ref
        [vMax iMax] = max( CCtmp,[], 2 );
        % find Ref with highest match
        [vOpt iRef] = max( vMax     );
        % also get "sub" index:
              iSub  = iMax(iRef);
        
        ordRef = [ordRef remainRef( iRef )];
        remainRef( iRef ) = [];
        
        ordSub = [ordSub remainSub( iSub )];
        remainSub( iSub ) = [];
        
        CCtmp(iRef,:) = [];
        CCtmp(:,iSub) = [];        
    end
    
    CCnew = corr( refVects(:,ordRef), subVects(:,ordSub) );
    flip  = sign( diag( CCnew ) );
    
    resort = sortrows([ordRef(:) ordSub(:) flip(:)], 1);

    Out.index = resort(:,2);
    Out.flip   = resort(:,3);    
end
