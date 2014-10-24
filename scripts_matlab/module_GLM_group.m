function output = module_GLM_group( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_GLM: module that performs General Linear Model analyses
% =========================================================================
%
%   Syntax:
%           output = module_erGLM( datamat, split_info )
%
%

% convolution
N_subject = numel(datamat);

for k = 1:N_subject
    if( split_info{k}.convolve == 0 )
        % get unconvolved design
        X_design{k} = split_info{k}.design_mat;
    else
        % convolve with HRF (must convert into seconds!)
        X_design{k} = design_to_hrf( split_info{k}.design_mat, split_info{k}.TR_MSEC/1000, [5.0 15.0] );
    end
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
    
% reproducibility and SPM:


% reproducibility
output.metrics.R = RR;
% optimal eigenimage
output.images    = rSPMZ;
% CV score timeseries, on unit-normed eigenimage
output.temp      = datamat'  * (rSPMZ ./ repmat( sqrt(sum(rSPMZ.^2)), [1 size(rSPMZ,2)] ) );
