function output = falff( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_FALFF: fractional amplitude of low-frequency fluctuations
% =========================================================================
%
%   Syntax:
%           output = module_falff( datamat, split_info )
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

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name  = 'falff';
output.attributes.design_type = 'nocontrast';
output.attributes.model_type  = 'univariate';
output.attributes.num_comp    = 'one_component';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%

if( ~isfield(split_info{1},'spm') || isempty(split_info{1}.spm) )
    disp('falff uses default raw maps (instead of z-scored)');
    split_info{1}.spm = 'raw';
end
if( ~isfield(split_info{1},'norm')  || isempty(split_info{1}.norm) )
    disp('falff uses default normalized variability');
    split_info{1}.norm = 1; 
end
%-------------------------------------------------------------------------%

%% INDIVIDUAL SUBJECT ANALYSIS
if( ~iscell(datamat) || length(datamat)==1 )
    
    % in cases where it is a single cell, revert to single-subject anaylsis
    if( iscell(datamat) && length(datamat)==1 ) datamat = datamat{1}; end
    split_info = split_info{1}; %% take entry from single cell

    % dimensions
    [Nvox Ntime] =  size(datamat); 
    Nhalf        = floor(Ntime/2);
    TR           = round(split_info.TR_MSEC./1000);
    % parameters for spectral estimation
    Fny    = 0.5 * (1/TR);                % nyquist frequency
    NFFT   = 2^nextpow2( Nhalf );         % Next power of 2 from length of time-axis
    f      = Fny*linspace(0,1,NFFT/2+1);  % fourier data corresponds to these frequency points

    kq=0;
    for(t=1:10:Nhalf)

        kq=kq+1; %%increment
        % get sum of spectral power values, for f > dataInfo.FreqCut
        powMat  = abs( fft( datamat(:,t:t+Nhalf) , NFFT ,2) ) / Nhalf; 
        powMat  = powMat(:,1:NFFT/2+1);
        
        % fractional power
        if( split_info.norm >0 )
             falff(:,kq) = sum( powMat(:, (f>0.015)&(f<0.08)),2) ./ sum(powMat(:, (f>0.0001)),2);
        else falff(:,kq) = sum( powMat(:, (f>0.015)&(f<0.08)),2);
        end            
    end
    % catch for bad-valued voxels
    falff(~isfinite(falff))=eps;

    % reproducibility estimates on "most independent" splits
    [output.metrics.R rSPMZ] = get_rSPM( falff(:,1), falff(:,end), 1 );

    if    ( strcmp( split_info.spm, 'zscore' ) ) output.images  = rSPMZ;
    elseif( strcmp( split_info.spm, 'raw'   ) )  output.images  =  mean(falff,2);
    else  error( 'invalid output type for falff (should be zscore or raw).' );
    end

    output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
    
%% GROUP LEVEL ANALYSIS
else
    
    N_resample = size(Resampling_Index,1);
    N_subject = numel(datamat);
    
    % compute per-subject connectivity patterns
    for(i=1:N_subject)

        % dimensions
        [Nvox Ntime] =  size(datamat{i}); 
        Nhalf        = floor(Ntime/2);
        TR           = round(split_info{i}.TR_MSEC./1000);
        % parameters for spectral estimation
        Fny    = 0.5 * (1/TR);                % nyquist frequency
        NFFT   = 2^nextpow2( Nhalf );         % Next power of 2 from length of time-axis
        f      = Fny*linspace(0,1,NFFT/2+1);  % fourier data corresponds to these frequency points

        kq=0; clear falff;
        for(t=1:10:Nhalf)

            kq=kq+1; %%increment
            % get sum of spectral power values, for f > dataInfo.FreqCut
            powMat  = abs( fft( datamat{i}(:,t:t+Nhalf) , NFFT ,2) ) / Nhalf; 
            powMat  = powMat(:,1:NFFT/2+1);
            % fractional power
            if( split_info{1}.norm >0 )
                falff(:,kq) = sum( powMat(:, (f>0.015)&(f<0.08)),2) ./ sum(powMat(:, (f>0.0001)),2);
            else
                falff(:,kq) = sum( powMat(:, (f>0.015)&(f<0.08)),2);
            end
        end
        falff_set(:,i) = mean(falff,2);
    end
    % remove non-finite computations
    falff_set(~isfinite(falff_set))=eps;

    rSPMZ=0; %% initialize rSPMz

    for k = 1:N_resample

        % split the data matrix
        set1 = Resampling_Index(k,:);
        set2 = setdiff(1:N_subject,set1);
    
        % split-half averages
        falff_map1 = mean( falff_set(:,set1),2);
        falff_map2 = mean( falff_set(:,set2),2);

        %% map manipulations for results

        % reweight the correlations, controlling vasc variance
        falff_map1 = falff_map1 .* split_info{1}.spat_weight;
        falff_map2 = falff_map2 .* split_info{1}.spat_weight;
        % average correlation map

        % reproducibility and SPM:
        [ R_allvox, rSPMZ_temp ] = get_rSPM( falff_map1, falff_map2, 1 );
        % record optima --> exclude seed ROI, although the impact is generally small
        rSPMZ = rSPMZ + rSPMZ_temp/N_resample;
        R(k) = R_allvox;
    end

    output.metrics.R =  median(R);
    % optimal eigenimage
    if    ( strcmp( split_info{1}.spm, 'zscore' ) ) output.images  = rSPMZ;
    elseif( strcmp( split_info{1}.spm, 'raw'   ) ) output.images  =  mean(falff_set,2);
    else  error( 'invalid output type for falff (should be zscore or raw).' );
    end
    
    for(k=1:N_subject)
    % seed timeseries, on unit-normed eigenimage
    output.temp{k}    = datamat{k}'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
    end
end
