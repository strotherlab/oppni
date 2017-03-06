function component_auto_filter( volname, outvolname, maskname, COMP_imag_path, COMP_temp_path )
%
% =========================================================================
% COMPONENT_AUTO_FILTER:  script for component based filtering of motion
% "spikes". Looks for temporal outliers and outliers in spatial patterns
% =========================================================================
%
% Syntax:
%           component_auto_filter( volname, outvolname, maskname, icapath )
% Input:
%           volname     = string, giving path/name of fMRI data (NIFTI or ANALYZE format)
%           outvolname  = string, giving path/name of output, interpolated fMRI data
%           maskname    = string, giving path/name of mask
%           icapath     = string, giving path of ICA data.
%
% Output:
%          -creates a new fMRI data volume, with name "outvolname", where
%           outlier components have been regressed out
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

if( exist('OCTAVE_VERSION','builtin') )
    % load stats packages
    pkg load statistics;
    pkg load optim; 
end
 
% load fMRI NIFTI-format data into MatLab
VV = load_untouch_nii(volname);
MM = load_untouch_nii(maskname); mask=double(MM.img);
VV.img = double(VV.img); %% format as double, for compatibility
[nx ny nz nt] = size(VV.img);

% load ICA data
Tmodes = load(COMP_temp_path);
% unit-norm temporal mode vectors
Tmodes = bsxfun(@rdivide,Tmodes,sqrt(sum(Tmodes.^2)));
II = load_untouch_nii(COMP_imag_path);
icamat = nifti_to_mat( II, MM );

fname = [COMP_temp_path,'_outlier_list.txt'];
oo = fopen(fname,'w');

%% (1) temporal spiking criteria

% run through timepoints (1...Ntime)
for(k=1:size(icamat,2))
for(t=1:nt)
    % estimate 15-TR time window
    wind=(t-7):(t+7);
    if(wind( 1 )< 1) wind(wind<1)  = []; end
    if(wind(end)>nt) wind(wind>nt) = []; end
                     wind(wind==t) = [];
	% get distance estimate at each timepoint, relative to median coordinate
    Tmode_dist(t,1) = (Tmodes(t,k) - median(Tmodes(wind,k))).^2;
end
% largest spike for this component
Tmode_max(k,1) = max(Tmode_dist);
end
% identify components that are significantly more "spiky" than the rest
par_ab  = gamfit( Tmode_max(:) );
outThr  = gaminv( 0.95, par_ab(1), par_ab(2) );
temp_spike_index   = find( Tmode_max > outThr );

if( isempty(temp_spike_index) )
      fprintf(oo,['Temporal spikes in Components#: (none) \n']);
else  fprintf(oo,['Temporal spikes in Components#: ', num2str(temp_spike_index(:)'),'\n']);
end
%% (2) spatial gradient criteria

% spatial gradient computed
[fx fy fz] = gradient( mean(VV.img,4) );
FXYZ = [fx(mask>0) fy(mask>0) fz(mask>0)];
% searching components
for(i=1:size(icamat,2)) 
    % measure multivariate corr with gradient axes
    [a b RHO(i,1)] = canoncorr( icamat(:,i).^2, FXYZ.^2 ); 
end
% distribution on edge-correlate
par_ab  = gamfit( RHO );
outThr  = gaminv( 0.95, par_ab(1), par_ab(2) );
spat_spike_index = find( RHO > outThr );

if( isempty(spat_spike_index) )
     fprintf(oo,['Spatial-pattern motion in Components#: (none) \n']);
else fprintf(oo,['Spatial-pattern motion in Components#: ', num2str(spat_spike_index(:)'),'\n']);
end
    
%% (3) filtering

all_spike_index = unique( [temp_spike_index(:); spat_spike_index(:)] );
T_spikes        = Tmodes(:,all_spike_index);

% detrend each slice
for( z=1:nz )
    %
    % take slice from 4D fMRI volumes; convert to vox x time matrix
    slcmat = reshape( VV.img(:,:,z,:), [],nt );
    slcmat = ols_regress_ex( slcmat, T_spikes, 1 );
end

fprintf(oo,['removed ',num2str(length(all_spike_index)),' components in total\n']);
fclose(oo);

% save the output results
save_untouch_nii(VV,outvolname);   

%%
function [ detrVol ] = ols_regress_ex( dataVol, regVcts, keepmean )
% 
%  OLS regression of timeseries from data matrix
%

% matrix dimensions
[Nmeas Nsamp] = size(dataVol);
% regressors + mean timecourse
X         = [ones(Nsamp,1) regVcts];
% beta map estimation
BetaFits  = inv( X' * X ) * X' * dataVol'; 

if( keepmean == 0 )
    %
    % OLS reconstruction of data, based on regressors
    detr_est  = ( X * BetaFits )';         
    % residual data
    detrVol   = dataVol - detr_est;    
else
    %
    % OLS reconstruction of data, WITHOUT mean-weight
    detr_est  = ( X(:,2:end) * BetaFits(2:end,:) )';         
    % residual data
    detrVol   = dataVol - detr_est;        
end
