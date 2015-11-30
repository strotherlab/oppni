function output = module_gconn( datamat, split_info )
%
% =========================================================================
% MODULE_GCONN: global functional connectivity maps
% =========================================================================
%
%   Syntax:
%           output = module_gconn( datamat, split_info )
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

% dimensions
[Nvox Ntime] =  size(datamat); 
% split into two equal halves
datamat_sp1 = datamat(:,1:ceil(Ntime/2));
datamat_sp2 = datamat(:,ceil(Ntime/2)+1:end);
% mean-center splits
datamat_sp1 = bsxfun(@minus,datamat_sp1,mean(datamat_sp1,2));
datamat_sp2 = bsxfun(@minus,datamat_sp2,mean(datamat_sp2,2));

% set default type if unspecified
if( ~isfield(split_info,'conn_metric') || isempty(split_info.conn_metric) )
    disp('gconn default: computed on correlation');
    split_info.conn_metric='corr';
end
% variance normalize if done on correlation/partial correlation
if( strcmp(split_info.conn_metric,'corr') || strcmp(split_info.conn_metric,'pcorr') )
    datamat_sp1 = bsxfun(@rdivide,datamat_sp1,std(datamat_sp1,0,2));
    datamat_sp2 = bsxfun(@rdivide,datamat_sp2,std(datamat_sp2,0,2));
end

if( strcmp(split_info.conn_metric,'cov') || strcmp(split_info.conn_metric,'corr') )
    % eigendecomposition
    [u l v] = svd(datamat_sp1,'econ');
    % global estime
    gconn1 = ((u*(l.^2)) * mean(u,1)')./ ( size(datamat_sp1,1) - 1 );

    % eigendecomposition
    [u l v] = svd(datamat_sp2,'econ');
    % global estime
    gconn2 = ((u*(l.^2)) * mean(u,1)')./ ( size(datamat_sp2,1) - 1 );
else
    % computing partial correlation
    
    % eigendecomposition
    [u l v] = svd(datamat_sp1,'econ');
    % per voxel precision
    sig2 = (u.^2) * diag(pinv(l.^2));
    % rescale by sqrt-precision
    us = bsxfun(@rdivide,u,sqrt(sig2));
    % global estime
    gconn1= -(us*pinv(l.^2))*mean(us)' + 2./Nvox;

    % eigendecomposition
    [u l v] = svd(datamat_sp2,'econ');
    % per voxel precision
    sig2 = (u.^2) * diag(pinv(l.^2));
    % rescale by sqrt-precision
    us = bsxfun(@rdivide,u,sqrt(sig2));
    % global estime
    gconn2= -(us*pinv(l.^2))*mean(us)' + 2./Nvox; %% adjustment for p(i,i)=-1  
end
% catch for bad-valued voxels
gconn1(~isfinite(gconn1))=eps;
gconn2(~isfinite(gconn2))=eps;

% reproducibility estimates on "most independent" splits
[output.metrics.R rSPMZ] = get_rSPM( gconn1, gconn2, 1 );

if    ( strcmp( split_info.spm, 'zscore' ) ) output.images  = rSPMZ;
elseif( strcmp( split_info.spm, 'raw'   ) )  output.images  = (gconn1+gconn2)./2;
else  error( 'invalid output type for gconn (should be zscore or raw).' );
end

output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));

    
   
    