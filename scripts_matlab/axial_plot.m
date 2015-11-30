function axial_plot( data, mask, Nsq, bound, transp )
%
% AXIAL_PLOT: simple script to plot brain images as axial slices
% 
% Syntax:
%              axial_plot( data, mask, Nsq, bound )
%
% where: 
%              data  = image, as single vector
%              mask  = 3d binary image volume
%              Nsq   = plots Nsq x Nsq image slices
%              bound = 2d vector specifying [min max] colourmap bounds
%              transp= transpose image slices
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: Sept. 1 2014
% ------------------------------------------------------------------------%

% mask dimensions
[Nx Ny Nz] = size(mask);
% reshape volume into 3D
tmp=mask; tmp(tmp>0)= double( data );
% get index on z-axis
s       = round( (Nz-1)/Nsq^2);
indlist = 2:s:Nz;
% adjust if it over-runs volume size
if( length(indlist) > Nsq^2 ) 
    indlist = indlist(1:Nsq^2 );
end

% define colourmap range
if( length(bound)==2 )
    pct = bound;
elseif(bound==1)
    refdist = abs( data(isfinite(data)) );
    pct     = prctile( refdist(refdist>0), 95 );
    pct     = [-pct pct];
elseif(bound==2)
    refdist = abs( data(isfinite(data)) );
    pct     = prctile( refdist, [5 95] );
end

if( isempty(transp) || transp==0 )

    % define full "image-panel"
    fullpanel = zeros( Nsq*Nx, Nsq*Ny );

    for(j=1:length(indlist))
        coldx =   rem(j,Nsq); %zero-indexed
        rowdx = floor(j/Nsq);
        fullpanel( rowdx*Nx+1:(rowdx+1)*Nx, coldx*Ny+1:(coldx+1)*Ny ) = tmp(:,:,indlist(j));
    end

else
    
    % define full "image-panel"
    fullpanel = zeros( Nsq*Ny, Nsq*Nx );

    for(j=1:length(indlist))
        coldx =   rem(j,Nsq); %zero-indexed
        rowdx = floor(j/Nsq);
        fullpanel( rowdx*Ny+1:(rowdx+1)*Ny, coldx*Nx+1:(coldx+1)*Nx ) = flipud(tmp(:,:,indlist(j))');
    end    
end
% initialize plot
figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1 0.1 0.50 0.80]);
imagesc( fullpanel, pct );
set(gca,'xtick',[],'ytick',[]);

