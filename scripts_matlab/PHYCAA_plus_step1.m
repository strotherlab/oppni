function output = PHYCAA_plus_step1( dataMat, dataInfo )
%
%==========================================================================
%  PHYCAA_PLUS_STEP1: this data-driven procedure takes a set of fMRI data 
%  matrices and estimates a map, that down-weights voxels with probable non
%  -neuronal tissue content (e.g. vasculature, ventricles). This map is 
%  used as a spatial prior for the PHYCAA_PLUS_STEP2 regression algorithm.
%==========================================================================
%
% SYNTAX:
%
%   output   = PHYCAA_plus_step1( dataMat, dataInfo )
%
% INPUT:
%
%   dataMat  = cell array of fMRI data splits (e.g. independent runs) used 
%              to estimate non-neuronal tissue. each cell entry should be 
%              an fMRI data matrix of size (voxels x time), although the 
%              number of timepoints can vary between splits. We recommend 
%              at least 2 splits (cell entries) to get a stabilized 
%              estimate of non-neuronal tissue.
%   dataInfo = structure containing the following fields:
%
%              dataInfo.TR           : rate of data acquisition (in sec.)
%              dataInfo.FreqCut      : high-frequency cutoff in Hz, for which f >FreqCut is primarily
%                                      physiological noise. DEFAULT is a recommended FreqCut=0.10 (Hz)
%              dataInfo.thresh_method: defined the cutoff to mask out "non-neuronal" tissue (voxels 
%                                      weighted =0) based on deviation from linearity (RSD) values:
%                                      'nothreshold'= no voxels masked out; threshold given by the max
%                                                     (RSD) value, and all voxels smoothly vary 
%                                                     between [0,1] - conservative approach
%                                      'noprior'    = declares top 5% of voxels non-neuronal (weight=0)
%                                                     This is the average threshold observed across 
%                                                     subjects and tasks in (Churchill & Strother, 2013)
%                                      'prior'      = maximizes overlap with a spatial prior, supplied 
%                                                     by 'dataInfo.priorMask'. Must be a binary vector.
%                                      DEFAULT (thresh_method=[]) is 'noprior' which gives robust results
%              dataInfo.priorMask    : OPTIONAL argument, a binary (voxels x 1) vector mask, mapping 
%                                      probable non-neuronal tissues (1=non-neuronal, 0=neuronal tissue) 
%                                      e.g. a CSF atlas. Used to theshold the non-neuronal weighting map
%                                      (see threshold_method above)
%              dataInfo.out_format   : format of output, noise corrected data (in outputs.dataMat_denoised)
%                                         0 = no output (dataMat_denoised=[]
%                                         1 = output downweighted data matrices
%                                      If not specified, DEFAULT value is out_format=1.
% OUTPUT:
%
%   output = structure containing the following fields:
%
%            output.NN_weight       : (voxel x 1) map of voxel weights, of range [0,1]. High weights 
%                                     (of ~1) are given to neuronal tissue, and weights -->0 for 
%                                     increasing non-neuronal content. This map is used to re-weight 
%                                     voxels and decrease variance contributed by non-neuronal tissue
%            output.NN_mask         : (voxel x 1) binary mask, indicating regions with more non-neuronal
%                                     than neuronal tissue - useful for comparing vascular territories
%            output.NN_rep          : average reproducibility of NN_wt estimates across the
%                                     data splits (measured via Pearson correlation)
%            output.dataMat_weighted: cell array of fMRI data splits (e.g. independent runs) after applying 
%                                     non-neuronal downweighting (empty if dataInfo.out_format==0)
%
% ------------------------------------------------------------------------%
%   Copyright 2013 Baycrest Centre for Geriatric Care
%
%   This file is part of the PHYCAA+ program. PHYCAA+ is free software: you 
%   can redistribute it and/or modify it under the terms of the GNU Lesser 
%   General Public License as published by the Free Software Foundation, 
%   either version 3 of the License, or (at your option) any later version.
% 
%   PHYCAA+ is distributed in the hope that it will be useful, but WITHOUT 
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
%   for more details.
% 
%   You should have received a copy of the GNU General Lesser Public 
%   License along with PHYCAA+. If not, see <http://www.gnu.org/licenses/>.
% 
%   This code was developed by Nathan Churchill Ph.D., University of Toronto,
%   during his doctoral thesis work. Email: nchurchill@research.baycrest.org
%
%   Any use of this code should cite the following publication:
%      Churchill & Strother (2013). "PHYCAA+: An Optimized, Adaptive Procedure for 
%      Measuring and Controlling Physiological Noise in BOLD fMRI". NeuroImage 82: 306-325
%
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%% ------------------------------------------------------------------------%

disp('Running PHYCAA+, Step-1: spatial weighting of non-neuronal tissues' );

% TR - inter-scan interval must be specified
if( ~isfield(dataInfo, 'TR' ) || isempty(dataInfo.TR) ) 
    %
    error('need to specify "dataInfo.TR" value (in sec.)');
end
% frequency cutoff (default 0.10 Hz)
if( ~isfield(dataInfo, 'FreqCut' ) || isempty(dataInfo.FreqCut) ) 
    %
    disp('    (default FreqCut=0.10 Hz)');
    %
    dataInfo.FreqCut       = 0.10;      
end
% thresholding method (default 'noprior')
if( ~isfield(dataInfo, 'thresh_method' ) || isempty(dataInfo.thresh_method)) 
    %
    disp('    (default "noprior" setting; threshold top 5%)');
    %
    dataInfo.thresh_method = 'noprior'; 
end
% if no choice on output, default is not to output processed data (only non-neuronal map)
if( ~isfield(dataInfo, 'out_format') || isempty(dataInfo.out_format)  ) 
    %
    disp('    (default out_format=1)');
    %
    dataInfo.out_format = 1; 
end
% If renormalization chosen (removes variance), implement ... otherwise retain var
if( ~isfield(dataInfo, 'renorm') || isempty(dataInfo.renorm)  ) 
    %
    disp('    (default renorm=0)');
    %
    dataInfo.renorm = 0; 
end

% store unaltered data, if chosen to output later
if( dataInfo.out_format > 0 ) dataMat_ref = dataMat; end

% get data dimensions
Nruns = length( dataMat );
for(r=1:Nruns) 
    % dimensions of each split
    [Nvox Ntime(r,1)] = size( dataMat{r} ); 
    % subtract temporal means
    dataMat{r} = bsxfun(@minus,dataMat{r},mean(dataMat{r},2));
end

% parameters for threshold estimation
NumList = (1:Nvox)'; 
Nqart   = round(Nvox/4);
% initialize matrix values
WTset     = zeros(Nvox,Nruns);
RSDset    = zeros(Nvox,Nruns);
percentNT = zeros(Nruns,1);
percentNN = zeros(Nruns,1);

for( r=1:Nruns ) %% estimate on each split

    % parameters for spectral estimation
    Fny    = 0.5 * (1/dataInfo.TR);                % nyquist frequency
    NFFT   = 2^nextpow2( Ntime(r) );           % Next power of 2 from length of time-axis
    f      = Fny*linspace(0,1,NFFT/2+1);  % fourier data corresponds to these frequency points
    numlow = sum( f <= dataInfo.FreqCut );           % count the number of frequency points below threshold

% =========================================================================
%   get ordered spectral estimates    

    % get sum of spectral power values, for f > dataInfo.FreqCut
    powMat  = abs( fft( dataMat{r} , NFFT ,2) ) / Ntime(r); 
    powMat  = powMat(:,1:NFFT/2+1);
    % index for non-ill-posed voxels
    ikeep = find( sum(powMat.^2,2) > eps );
    % if any ill-posed voxels, drop from the index
    if(length(ikeep)<Nvox)
        disp(['you have ',num2str(Nvox-length(ikeep)),' "empty" voxel timeseries; correcting, but check your masks!']);
        powMat = powMat(ikeep,:);
    end
    % parameters for threshold estimation
    NumList = (1:size(powMat,1))'; 
    Nqart   = round(Nvox/4);
    
    % take variance (possibly renormalized)
    if    ( dataInfo.renorm == 0 ) powHgh  = sum( sqrt(powMat(:,numlow+1:end)), 2 );
    elseif( dataInfo.renorm == 1 ) powHgh  = sum( sqrt(powMat(:,numlow+1:end)), 2 ) ./ sum( sqrt(powMat(:,3:end)), 2 );
    else  error('if renorm is specified, must be =0,1');
    end

    % order high-frequency power values, smallest to largest
    PowSort = sort(powHgh);
    PowIndx = sortrows([powHgh NumList],1); 
    PowIndx = PowIndx(:,2);
    % Get linear fit for central 50% of voxel values (linear pieces)
    P            = polyfit(NumList(Nqart:end-Nqart),PowSort(Nqart:end-Nqart),1);
    PowLQ        = polyval(P,NumList);
    RSD          = zeros( length( powHgh ),1 );
    RSD_sorted   = PowSort-PowLQ;
    RSD(PowIndx) = RSD_sorted;
    
% =========================================================================
%   identifying the "neuronal-tissue" threshold    

    % identify threshold for deviation >linearity, corresponding to nonneuronal  tissue
    stdlin = std( RSD_sorted(Nqart:end-Nqart) );
    NCut   = norminv(0.99,0,stdlin);
    % percentile threshold, for significant deviation from linearity
    MinThr = 100 * sum( RSD < NCut ) ./ length( RSD );
    % threshold RSD value
    minRSD =  prctile( RSD, MinThr );    
    
% =========================================================================
%   identifying the "non-neuronal tissue" threshold
%
%   --- method of thresholding depends on user choice ---

    if    ( strcmp( dataInfo.thresh_method, 'nothreshold' ) )
        %
        % no voxels declared non-neuronal (most conservative)
        MaxThr = 100;
        %
    elseif( strcmp( dataInfo.thresh_method, 'noprior'     ) )
        %
        % top 5% of voxels declared non-neuronal
        MaxThr = 95.0;
        %
    elseif( strcmp( dataInfo.thresh_method, 'prior' )       )
        %
        if( ~isfield( dataInfo, 'priorMask' ) )
            %
            error( 'Need to specify binary spatial prior in "dataInfo.priorMask"' );
        else
            NN_msk = double( dataInfo.priorMask > 0.5 );
            % drop out "ill-posed" voxels from mask
            if(length(ikeep)<Nvox)
               NN_msk = NN_msk(ikeep,1);
            end
        end
        %
        % find threshold that maximizes overlap with the binary "prior" mask
        k=0;
        % list of percentile thresholds to test
        PCTLIST = MinThr(1):0.20:99.50;
        % initialize list of overlap values
        Overlap = zeros( length(PCTLIST), 1 );
        % test different thresholds
        for( qc=PCTLIST )
            %
            k=k+1;
            % map of supra-threshold voxels (non-neuronal regions)
            Z=double(RSD>prctile(RSD,qc));        
            % Jaccard overlap of (thresholded RSD, binary prior mask)
            Overlap(k,1) = ( sum( Z.*NN_msk ) )./ sum( double((Z+NN_msk)>0) ); 
        end

        % convolve with simple running average smoother (width of 1%), to avoid local minima in curve
        Overlap = conv( Overlap, ones(1,5)./5 ); 
        Overlap = Overlap(3:end-2);
        % point of maximized overlap
        [omax imax]   = max( Overlap ); 
        % percentile of maximized overlap
        MaxThr = PCTLIST(imax);
        % catch for extreme low values
        if( MaxThr < 90.0 ) 
            disp( 'Warning: non-neuronal threshold is overly conservative at_',num2str(MaxThr),'%...readjusting to_90%');
            MaxThr = 90.0; 
        end
    else
        disp( 'Error: need to specify thresholding method.');
        return;
    end

    % threshold RSD value
    maxRSD =  prctile( RSD, MaxThr );
    
% =========================================================================
%   converting into map of spatial weights, of values [0, 1]
    
    % convert RSD values into [0,1] scale
    % where 0= greater than maxRSD (non-neuronal threshold)
    %       1= less than minRSD (neuronal threshold)
    HI_dev =   RSD  - minRSD; 
    HI_dev(HI_dev<0) = 0; 
    HI_wt = HI_dev ./ maxRSD; 
    HI_wt(HI_wt>1)=1; 
    HI_wt = 1-HI_wt;

    % if missing voxels, adjust:
    if(length(ikeep)<Nvox)
        tmp = 0.5*ones(Nvox,1); % default setting of 0.5=equal likelihood of either type, all vox
        tmp(ikeep) = HI_wt; % for retained voxels, replace with actual values
        HI_wt = tmp; % swap in for original HI_wt
    end
    
    % recording split information:
    WTset(:,r)   = HI_wt;
    percentNT(r) = MinThr;
    percentNN(r) = 100 - MaxThr;
end

% ====== Output-1: physiological mask parameters ====== %

    output.NN_weight = mean( WTset, 2 ); % average non-neuronal weighting map
    output.NN_mask   = prod( double(WTset>=0.5), 2 ); % intersection of mask estimates
    % reproducibility across splits
    CC = triu( corr(WTset), 1 );
    output.NN_rep  = mean( CC(CC~=0) );

% ====== Output-2: reweighted fMRI data, if selected ====== %

    if    (dataInfo.out_format==0) % do not produce de-noised data
        %
        weightMat = [];

    elseif(dataInfo.out_format==1) % non-neuronal downweighting

        % initialize
        weightMat = cell(Nruns,1);
        % rescale based on physio. noise
        for(r=1:Nruns)  weightMat{r} = bsxfun(@times,dataMat_ref{r},output.NN_weight);  end
    end

    % down-weighted data matrices
    output.dataMat_weighted = weightMat;
