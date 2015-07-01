function output = module_erGLM( datamat, split_info )
%
% =========================================================================
% MODULE_ERGLM: module that performs General Linear Modelling of a single
% event-related stimulus type
% =========================================================================
%
%   Syntax:
%           output = module_erGLM( datamat, split_info )
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
[Nvox Ntime] = size( datamat );

% building design vector: we want to subsample at 100 ms (faster than human RT)
Nsubs  = max([1 round(split_info.TR_MSEC/100)]); % (#samples per TR) catch in case TR<100ms
design = zeros( Ntime*Nsubs, 1 );     % initialize design matrix
% index allocates onsets to appropriate design-points
didx   = unique(round( split_info.onsetlist./(split_info.TR_MSEC/Nsubs) ));
% catch + adjust overruns, setvalue=1 on design vector
didx(didx==0)          = 1;
didx(didx>Ntime*Nsubs) = Ntime*Nsubs;
design( didx )         = 1;
% convolve with HRF (must convert into seconds!)
HRFdesign = design_to_hrf( design, (split_info.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
% now, subsample back to get HRF at actual fMRI sampling rate
HRFdesign = HRFdesign( round(Nsubs/2): Nsubs : end );      

% split the data matrix
datasplit_1 = datamat(:,1:round(Ntime/2));
datasplit_2 = datamat(:,round(Ntime/2)+1:end);
% split the hrf-convolved design vector
HRFdes_1 = HRFdesign(1:round(Ntime/2));
HRFdes_2 = HRFdesign(round(Ntime/2)+1:end);
% run split glm analyses
out1 = GLM_model_fmri( datasplit_1, 0, [], HRFdes_1, [] );
out2 = GLM_model_fmri( datasplit_2, 0, [], HRFdes_2, [] );

% reweight the betas
beta1 = out1.Beta_signl .* split_info.spat_weight;
beta2 = out2.Beta_signl .* split_info.spat_weight;

% reproducibility and SPM:
[ RR, rSPMZ ] = get_rSPM( beta1, beta2, 1 );

% reproducibility
output.metrics.R = RR;
% optimal eigenimage
output.images  = rSPMZ;
% CV score timeseries, on unit-normed eigenimage
output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
