function output = PHYCAA_plus_step2( dataMat, dataInfo )
% 
%==========================================================================
%  PHYCAA_PLUS_STEP2: regresses physiological timecourses from fmri data, 
%  after applying PHYCAA_PLUS_STEP1. This algorithm uses multivariate 
%  Canonical Autocorrelations Analysis (CAA) to identify autocorrelated 
%  timeseries with spatially reproducible patterns in non-neuronal tissue
%  (macro-vasculature, ventricles, sinuses), and regresses them out.
%==========================================================================
%
% SYNTAX:
%
%   output   = PHYCAA_plus_step2( dataMat, dataInfo )
%
% INPUT:
%
%   dataMat  = cell array of fMRI data splits (e.g. independent runs) used 
%              to estimate physiological noise components. each cell entry 
%              should be an fMRI data matrix of size (voxels x time), 
%              although the number of timepoints can vary between splits.
%              This model requires at least 2 splits (cell entries) in 
%              order to estimate reproducibility.
%   dataInfo = structure containing the following fields:
%
%              dataInfo.physio_map: a (voxels x 1) brain image vector, mapping out probable sources of 
%                                   physiological noise (non-neuronal tissues). This should be the 
%                                   'NN_weight' map produced by PHYCAA_plus_step1, but you can input 
%                                   any inverse probabilistic map (spatial map with values of range 
%                                   [0,1], where values->0 denote increased non-neuronal tissue 
%                                   likelihood), to target other sources of noise (e.g. white matter).
%              dataInfo.task_SPMs : (voxels x K) matrix of K brain activation map(s) related to BOLD 
%                                   signal(s) of interest. Activation maps are formatted as a set of 
%                                   column vectors. These are regressed from 'dataMat' in order to 
%                                   estimate physiological noise in the residual data; protects against 
%                                   over-regression of BOLD signal. If you don't want to estimate the 
%                                   residual subspace, make this entry empty (dataInfo.task_SPMs = []).
%              dataInfo.out_format: format of output, noise corrected data (in outputs.dataMat_denoised)
%                                    0= no output (dataMat_denoised=[])
%                                   -1= regression only (don't down-weight non-neuronal tissues).
%                                    1= regression + non-neuronal downweighting 
%                                   If not specified, DEFAULT value is out_format=1.
%              dataInfo.comp_crit : parameter determines how conservative physiological component 
%                                   selection is. Value can range (0 <= comp_crit < 1), where a larger 
%                                   comp_crit indicates more conservative selection, i.e. variance must 
%                                   be more concentrated in non-neuronal tissue. If not specified, 
%                                   DEFAULT value is comp_crit=0, which generally gives robust results.
%              dataInfo.keepmean  : PHYCAA+ subtracts voxel means of each run, before component estimation.
%                                     0= voxel means are discarded 
%                                     1= voxel means re-added after noise regression
%                                   If not specified, DEFAULT value is keepmean=0.
%              dataInfo.max_PC_dim: cutoff on the maximum explored PC dimensionality, used to identify 
%                                   a physio. noise subspace. The model can explore a range of 
%                                   [3PCs-->50% of data dimensionality]. If your data is large (e.g. 100s of timepoints), 
%                                   you may want to set a lower threshold to decrease compute time. 
%                                   If not specified, DEFAULT value is max_PC_dim=50.
%
% OUTPUT: 
%
%   output = structure containing the following fields:
%
%            output.Physio_SPM      : (voxels x 1) vector, Z-scored map of reproducibible 
%                                     physiological component variance
%            output.Physio_rep      : reproducibility of physiological noise map
%            output.PC_CV           : 2D vector of optimal [#PCs; #CVs] for physiological estimation
%            output.Physio_Tset     : cell array of physiological noise components. Each cell contains 
%                                     (time x components) matrix of physiological timecourses for that split
%            output.dataMat_denoised: cell array of fMRI data splits (e.g. independent runs) after  
%                                     physiological noise correction (empty if dataInfo.out_format==0)
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
% ------------------------------------------------------------------------%%

disp('Running PHYCAA+, Step-2: component regression' );

%% (1) Pre-specifying model parameters

% threshold for significant correlation in CVs
signif_cutoff  = 0.05;
% time-lag in autocorrelation model (offSet>1 does not appears to improve estimation)
offSet         =    1;

% if user has not included a non-neuronal reference, terminate
if( ~isfield(dataInfo, 'physio_map') || isempty(dataInfo.physio_map)  )  
    %
    error('need to specify non-neuronal tissue map in "dataInfo.physio_map"');
end
% if no choice on output, default is not to output processed data (only physiological maps + matfiles)
if( ~isfield(dataInfo, 'out_format') || isempty(dataInfo.out_format)  ) 
    %
    disp('    (default out_format=1)');
    %
    dataInfo.out_format = 1; 
end
% cat if out_format is not an acceptable value
if( sum( dataInfo.out_format == [-1 0 1] ) ==0 )
    %
    error('need to set "dataInfo.out_format" to {-1, 0 or 1}');
end
% catch if comp_crit not specified
if( ~isfield(dataInfo, 'comp_crit') || isempty(dataInfo.comp_crit)  ) 
    %
    disp('    (default comp_crit=0)');
    %
    dataInfo.comp_crit = 0.00; 
elseif( (dataInfo.comp_crit < 0) || (dataInfo.comp_crit >= 1.0) )
    %
    % catch for bad comp_crit values
    error('"dataInfo.comp_crit" value must be non-negative and less than 1!'); 
end
% if user has not specified to restore mean, this is NOT done
if( ~isfield(dataInfo, 'keepmean') || isempty(dataInfo.keepmean)  )  
    %
    disp('    (default keepmean=0)');
    %
    dataInfo.keepmean=0;  
end

% get data dimensions
Nruns = length( dataMat );
% store data without mean removed, if this is selected
if( dataInfo.keepmean > 0 ) dataMat_keepmean = dataMat; end
Ntime = zeros(Nruns,1);
for(r=1:Nruns) 
    % dimensions of each split
    [Nvox,Ntime(r,1)] = size( dataMat{r} ); 
    % subtract temporal means
    dataMat{r} = bsxfun(@minus,dataMat{r},mean(dataMat{r},2));
    %dataMat{r} = dataMat{r} - repmat( mean(dataMat{r},2), [1 size(dataMat{r},2)] );      
end

% shortest time length, limits range of PCs
Ntime_min = min(Ntime);

% if user has not specified maximal PC range, default cut of 50
if( ~isfield(dataInfo, 'max_PC_dim') || isempty(dataInfo.max_PC_dim)  )  
    % reset PC range to smallest of 50 pcs / 50% data dimensionality
    dataInfo.max_PC_dim= min( [50 round(Ntime_min/2 - 1)] );  
    %
elseif( dataInfo.max_PC_dim < 3 ) 
    % if value is too small, reset to minimum of 3
    dataInfo.max_PC_dim=3;  
else
    % reset PC range to smallest of specified #pcs / 50% data dimensionality
    dataInfo.max_PC_dim= min( [dataInfo.max_PC_dim round(Ntime_min/2 - 1)] );  
end

%% (2) Regressing out task effect

task_course = cell(Nruns,1);
regMethod   = 0;

% if task spm(s) entry is given, and non-empty, estimate signal timecourse
if(  isfield( dataInfo,'task_SPMs' ) && ~isempty(dataInfo.task_SPMs)  ) 
    %
    % then flag indicates residual subspace estimation
    regMethod   = 1;
    % preallocate signal timecourse matrix
    task_course = cell(Nruns,1);
    % estimates the signal timecourse(s) for each data split
    for(r=1:Nruns)
        % downweight the non-neuronal tissue regions (minimize bias on signal estimation)
        dataNrm        = bsxfun(@times, dataMat{r}, dataInfo.physio_map );
        % project SPMs onto fMRI data, to get signal timecourses
        task_course{r} = dataNrm'* zscore( dataInfo.task_SPMs );
    end
end

% Initialize task residual matrix
residualMat = dataMat;
% initialize PCA-space representation
Q_c = cell(Nruns,1);
% regress signal, if specified
if( regMethod>0 )
    %
    disp( '    (task SPM included - performing residual subspace estimation)')
    %
    for(r=1:Nruns) 
	    % Validate data matrix - LMP
	    if (any(isinf(dataMat{r}) == 1) == 1)
	       disp("WARNING dataMat contains Inf elements")
	    end

	    if (any(isnan(dataMat{r}) == 1) == 1)
	       disp("WARNING dataMat contains NaN elements")
	    end
      
        % regress out estimated signal effect
        residualMat{r} = ols_regress_ex( dataMat{r}, task_course{r}, 0 );
        % using SVD to reduce dimensionality of residual data
        [temp,S_c,V_c] = svd( residualMat{r},'econ' );  
        % PC-space representation
        Q_c{r} = V_c * sqrt( S_c );
    end
else
    disp( '    (no task SPM included)')
    %
    for(r=1:Nruns)            
        % using SVD to reduce dimensionality of residual data
        [temp,S_c,V_c] = svd( residualMat{r},'econ' );
        % PC-space representation
        Q_c{r} = V_c * sqrt( S_c );
    end
end

%% (3) Iterative estimation of autocorrelated timeseries

% select PC subsets:
% Maximum should be < 50% of time-points (samples) to ensure stable CAA solution    
pcsets = 3:dataInfo.max_PC_dim;
% initialize reproducibility measure + counting measure
opt_REP  =-1;
opt_PCs  = 0;
opt_CVs  = 0;
rspm_opt = zeros(Nvox,1);
% initialize cell matrix of physiological timecourses
Vset_phy = cell( Nruns, 1 );
Tset_phy = cell( Nruns, 1 );
opt_Tset = cell( Nruns, 1 );
%
% count # PC subspaces explored
iter_s    =   0;
compcount = zeros( Nruns,1 );

% now iterate through each PC subspace dimensionality

Ind1 = dataInfo.physio_map < ( 0.50 - (dataInfo.comp_crit/2) );
Ind2 = dataInfo.physio_map > ( 0.50 + (dataInfo.comp_crit/2) );
for(pcs=pcsets)

    iter_s=iter_s+1;

    for(r=1:Nruns) 
        %
        % estimate temporal autocorrelation maximized "sources"
        Q1 = Q_c{r}( 1:(end-offSet) , 1:pcs ); % un-offset
        Q2 = Q_c{r}( (offSet+1):end , 1:pcs ); % offsetted timeseries
        % canonical correlations on time-lagged data
        [A,B,R,U,V,stats] = canoncorr(Q1,Q2); 

        % getting stable "average" autocorrelated timeseries
        a=[U(1,:)]; b=[U(2:end,:) + V(1:end-1,:)./2]; c=[V(end,:)]; 
        % normalizing the timeseries to unit variance
        tset = unitnorm_ex([a;b;c],1);

        % selection based on temporal criteria:
        % threshold keeps only significantly autocorrelated components
        Tset = tset(:, stats.pChisq(:) < signif_cutoff );

        if( ~isempty( Tset ) )

            % spatial mapping of components:
            % 1. beta maps corresponding to each component (simple because Tset is unit-normed)
            Bset       = residualMat{r} * Tset;
            % 2. total variance at each voxel
            % For speed up- I changed these lines (Babak)
%             tic;dat_var1   = repmat( sum(residualMat{r}.^2,2), [1 size(Bset,2)] );
%             % catch for zero-variance
%             dat_var1(dat_var1==0) = eps;
%             % 3. variance explained by each component, per voxel
%             est_var1   = Bset.^2;
%             % 4. maps: fraction of variance explained at each voxel, for every component
%             Vset       = est_var1 ./ dat_var1;toc
            
            est_var1   = Bset.^2;
            dat_var1   = sum(residualMat{r}.^2,2);
            dat_var1(dat_var1==0) = eps;
            Vset       = bsxfun(@rdivide,est_var1,dat_var1);

            % selection based on spatial criteria:
            % median variance in non-neuronal vs. neuronal tissues 
            % comp_crit -> dictates the threshold for non-neuronal / neuronal tissue
            vmed_nonneur=median( Vset(Ind1 , : ) );
            vmed_neuronl=median( Vset(Ind2 , : ) );
            % index of components where ratio of median variances >=1.0 (mostly noise)
            phyIdx = find( (vmed_nonneur./vmed_neuronl) >= 1.0 );
            % keep spatial variance maps + timeseries of "physiological" components
            Vset_phy{r}  = Vset(:, phyIdx );
            Tset_phy{r}  = Tset(:, phyIdx );
            compcount(r) = length( phyIdx );
        else
            compcount(r) = 0;
        end

    end

    % if we can't find physio components for one split, set reproducibility to zero
    if( sum( compcount==0 ) > 0 )
          phy_REP = -1;
    else        
        % for every possible pair of splits, compute reproducibility and
        % Z-scored map of physiological regions

        rspm_phy = zeros(Nvox,1);
        phy_REP  = 0;
        kij      = 0;
        %
        for(i=1:Nruns-1)
        for(j=i+1:Nruns)
            %
            kij = kij+1;
            %
            [new_REP, rspm_new] = get_rSPM_ex( sum(Vset_phy{i},2),sum(Vset_phy{j},2), 1); 
            % 
            phy_REP  = phy_REP  + new_REP;
            rspm_phy = rspm_phy + rspm_new;
        end
        end
        % now divide by #pairwise comparisons, to get mean
        phy_REP  = phy_REP/kij;
        rspm_phy = rspm_phy./kij;

        % update if this PC space gives greater reproducibility than the last-best
        if( phy_REP > opt_REP )
            %
            opt_REP  = phy_REP;   % update component with greatest Corr
            rspm_opt = rspm_phy;   % update rSPM of optimal components
            opt_PCs  = pcs;       % record this PC-subspace size
            opt_CVs  = [size(Vset_phy{1},2) size(Vset_phy{2},2)] ;       % record this PC-subspace size
            %
            opt_Tset = Tset_phy;
        end
    end
end

if    (dataInfo.out_format == 0) % do not produce de-noised data
    %
    denoisedMat = [];

elseif(dataInfo.out_format ==-1) % regress out physiology
    %
    % initialize
    denoisedMat = cell(Nruns,1);
    % regress physiological components:
    for(r=1:Nruns)
        %
        if( dataInfo.keepmean > 0 )
              %
              denoisedMat{r} = ols_regress_ex( dataMat_keepmean{r}, opt_Tset{r}, 1 );
        else  denoisedMat{r} = ols_regress_ex( dataMat{r},          opt_Tset{r}, 0 );
        end
    end

elseif(dataInfo.out_format == 1) % regress out physiology + downweighting
    %
    % initialize
    denoisedMat = cell(Nruns,1);
    % regress physiological components:
    for(r=1:Nruns)
        %
        if( dataInfo.keepmean > 0 )
              %
              denoisedMat{r} = ols_regress_ex( dataMat_keepmean{r}, opt_Tset{r}, 1 );
        else  denoisedMat{r} = ols_regress_ex( dataMat{r},          opt_Tset{r}, 0 );
        end
        denoisedMat{r} = bsxfun(@times, denoisedMat{r}, dataInfo.physio_map); 
        %denoisedMat{r} = denoisedMat{r} .* repmat( dataInfo.physio_map, [1 Ntime(r)] );
    end
end

% ============= OUTPUTS ================== %
output.Physio_rep       = opt_REP;     % reproducibility of physiological noise
output.Physio_SPM       = rspm_opt;    % Z-scored physiologica map
output.PC_CV            = [mean(opt_PCs) mean(opt_CVs)]; % optimal PC#/CV#
output.Physio_Tset      = opt_Tset;    % cell matrix of physiological timecourses
output.dataMat_denoised = denoisedMat; % data with physiological effects removed


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

%%
function [ Xnorm ] = unitnorm_ex( X, flag )
%
% quick normalization procedure
%

% subtract mean if specified
if(flag>0)
   X = bsxfun(@minus,X,mean(X));
end
% normalize to unit variane
norm  = sqrt(sum( X.^2 ));
Xnorm = bsxfun(@rdivide,X,norm);

%%
function [ rep, rSPM ] = get_rSPM_ex( vect1, vect2, keepMean )
%
% get reproducible, Z-scored activation map
%

rep = corr(vect1, vect2);

%(1) getting the mean offsets (normed by SD)
normedMean1 = mean(vect1)./std(vect1);
normedMean2 = mean(vect2)./std(vect2);
%    and rotating means into signal/noise axes
sigMean = (normedMean1 + normedMean2)/sqrt(2);
%noiMean = (normedMean1 - normedMean2)/sqrt(2);
%(2) getting  signal/noise axis projections of (zscored) betamaps
zvect1 = zscore(vect1);
zvect2 = zscore(vect2);

sigProj = ( zvect1 + zvect2 ) / sqrt(2);
noiProj = ( zvect1 - zvect2 ) / sqrt(2);
% noise-axis SD
noiStd = std(noiProj);
%(3) norming by noise SD:
%     ...getting the (re-normed) mean offsets
sigMean = sigMean./noiStd;
%noiMean = noiMean./noiStd; 
%  getting the normed signal/noise projection maps
sigProj = sigProj ./ noiStd;
%noiProj = noiProj ./ noiStd;

% Produce the rSPM:
if    ( keepMean == 1 )   rSPM = sigProj + sigMean;
elseif( keepMean == 0 )   rSPM = sigProj;
end
