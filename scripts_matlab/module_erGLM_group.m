function output = module_erGLM_group( datamat, split_info, Resampling_Index )
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

N_subject = numel(datamat);

for k = 1:N_subject
    [Nvox,Ntime] = size( datamat{k});
    
    % building design vector: we want to subsample at 100 ms (faster than human RT)
    Nsubs  = max([1 round(split_info{k}.TR_MSEC/100)]); % (#samples per TR) catch in case TR<100ms
    design = zeros( Ntime*Nsubs, 1 );     % initialize design matrix
    % index allocates onsets to appropriate design-points
    didx   = unique(round( split_info{k}.onsetlist./(split_info{k}.TR_MSEC/Nsubs) ));
    % catch + adjust overruns, setvalue=1 on design vector
    didx(didx==0)          = 1;
    didx(didx>Ntime*Nsubs) = Ntime*Nsubs;
    design( didx )         = 1;
    % convolve with HRF (must convert into seconds!)
    HRFdesign = design_to_hrf( design, (split_info{k}.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
    % now, subsample back to get HRF at actual fMRI sampling rate
    X_design{k} = HRFdesign( round(Nsubs/2): Nsubs : end );
end
N_resample = size(Resampling_Index,1);

for i = 1:N_resample
    % split the data matrix

    index = randperm(N_subject);
    set1 = Resampling_Index(k,:);
    set2 = setdiff(1:N_subject,set1);
    
    datasplit1 = [];
    X_design1  = [];
    for k = 1:length(set1)
        datasplit1 = [datasplit1 datamat{set1(k)}];
        X_design1  = [X_design1;X_design{set1(k)}];
    end
    
    datasplit2 = [];
    X_design2  = [];
    for k = 1:length(set2)
        datasplit2 = [datasplit2 datamat{set2(k)}];
        X_design2  = [X_design2;X_design{set2(k)}];
    end
    
    % run split glm analyses

    out1 = GLM_model_fmri( datasplit1, 0, [], X_design1, [] );
    out2 = GLM_model_fmri( datasplit2, 0, [], X_design2, [] );
    
    % reweight the betas

    beta1 = bsxfun(@times,out1.Beta_signl,split_info.spat_weight);
    beta2 = bsxfun(@times,out2.Beta_signl,split_info.spat_weight);
    [ RR_temp, rSPMZ_temp ] = get_rSPM( beta1, beta2, 1 );
    rSPMZ = rSPMZ + rSPMZ_temp/N_resample;
    RR_stack(i,:) = RR_temp;
end

RR = median(RR_stack,1); 

% reproducibility
output.metrics.R = RR;
% optimal eigenimage
output.images  = rSPMZ;
% CV score timeseries, on unit-normed eigenimage
output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
