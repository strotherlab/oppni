function output = run_analyses_wrapper( datamat, split_info, analysis_model )
%
% =========================================================================
% RUN_ANALYSES_WRAPPER: general "wrapper" script for fMRI analysis, to be
% used in combination with Pipeline_STEP2_MM. Allows users to specify different 
% analysis models, and can substitute them as "modules".
% This is the wrapper for single-subject analysis. For group-level
% analyses, see "group_analyses_wrapper"!
% =========================================================================
%
% Syntax:   output = run_analyses_wrapper( datamat, split_info, analysis_model )
%  Input:
%
%          datamat : (voxels x time) fmri data matrix
%       split_info : structure with fields required to run individual analysis modules.
%                    The required fields for current modules are listed below
%   analysis_model : string specifying chosen analysis model. Options also listed below.
%
% Output: 
%           output : structure with following fields
%
%                      output.metrics  : structure containing performance metrics
%                      output.images   : (voxels x components) matrix of brain maps
%                      output.temp     : (time x components) matrix of timeseries
%                      output.modeltype: string specifying either
%                                          'onecomp'   = single component image
%                                          'multicomp' = many component images
%
% =========================================================================
%      analysis_model: 'LDA'   (linear discriminant)
%    * predictive multivariate analysis for block design, with 2 task conditions *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.cond1_sp1 = vector of volumes indices for task condition 1, for data split 1
%                  split_info.cond2_sp1 = vector of volumes indices for task condition 2, for data split 1
%                  split_info.cond1_sp2 = vector of volumes indices for task condition 1, for data split 2
%                  split_info.cond2_sp2 = vector of volumes indices for task condition 2, for data split 2
%                  split_info.drf       = scalar value of range (0,1), indicating the fraction of full-date
%                                         PCA subspace to keep during PCA-LDA analysis
%
%                  split_info.type      = 'block'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'GNB'  (gaussian naive bayes)
%    * predictive univariate analysis for block design, with 2 task conditions *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.cond1_sp1 = vector of volumes indices for task condition 1, for data split 1
%                  split_info.cond2_sp1 = vector of volumes indices for task condition 2, for data split 1
%                  split_info.cond1_sp2 = vector of volumes indices for task condition 1, for data split 2
%                  split_info.cond2_sp2 = vector of volumes indices for task condition 2, for data split 2
%
%                  split_info.decision_model = string specifying type of decision boundary. either:
%                                                'linear'   : pooled covariance model
%                                                'nonlinear': class-specific covariances
%
%                  split_info.type      = 'block'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'GLM'  (standard General Linear Model)
%    * basic univariate analysis for general stimuli *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.design_mat = design column-matrix, with dimensions (T timepoints  x  K stimulus types)
%                  split_info.convolve   = binary value, for whether design matrix should be convolved 
%                                          with a standard SPMG1 HRF. 0= do not convolve / 1=perform convolution
%
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'erCVA'  (event-related Canonical Variates Analysis)
%    * predictive multivariate analysis for event-related design, with 1 task type onset *
% -------------------------------------------------------------------------
% 'split_info' fields:
% 
%                  split_info.onsetlist = vector list of integers representing stimulus onset times
%                                         in milliseconds, zero-relative to the first TR volume
%                  split_info.Nblock    = number of equal sized splits to break the data into, to do 
%                                         time-locked averaging. Must be at least 2, with even numbers >=4 
%                                         recommended to obtain robust covariance estimates
%                  split_info.WIND      = window size to average on, in TR (usually in range 6-10 TR)
%                  split_info.drf       = scalar value of range (0,1), indicating the fraction of full-date
%                                         PCA subspace to keep during PCA-LDA analysis
%                  split_info.subspace  = string specifying either:
%                                         'onecomp'   = only optimize on CV#1 
%                                         'multicomp' = full multidimensional subspace
%                                         (DEFAULT = 'onecomp')
%
%                  split_info.type      = 'event'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'erGNB'  (event-related gaussian naive bayes)
%    * predictive univariate analysis for event-related design, with 1 task type onset *
% -------------------------------------------------------------------------
% 'split_info' fields:
% 
%                  split_info.onsetlist = vector list of integers representing stimulus onset times
%                                         in milliseconds, zero-relative to the first TR volume
%                  split_info.Nblock    = number of equal sized splits to break the data into, to do 
%                                         time-locked averaging. Must be at least 2, with even numbers >=4 
%                                         recommended to obtain robust covariance estimates
%                  split_info.WIND      = window size to average on, in TR (usually in range 6-10 TR)
%
%                  split_info.type      = 'event'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'erGLM'  (event-related General Linear Model with HRF estimation)
%    * basic univariate analysis for event-related task design; models a single HRF *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.onsetlist = vector list of integers representing stimulus onset times
%                                         in milliseconds, zero-relative to the first TR volume
%
%                  split_info.type      = 'event'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'SCONN'   (seed-based connectivity analysis)
%    * voxelwise correlations with seed Region of Interest (ROI) *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.seed_name = string giving location of a seed ROI volume. This must
%                                         be a 3D binary volume, of same dimensions as input fMRI data.
%                                               Voxels in ROI=1 / non-ROI voxels = 0
%                  split_info.spm       = string specifying format of output SPM. Options include: 
%                                         'raw'  : map of voxelwise seed correlations  
%                                         'zscore': Z-scored map of reproducible correlation values 
%
%                  split_info.type      = 'nocontrast'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'gPCA'   (generalization of Principal Components)
%    * measures overall reliability of K-dimensional covariance subspace *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.num_PCs   = integer defining size of Principal Components subspace,
%                                         where we retain the first "1 to num_PCs" 
%                  split_info.type      = 'nocontrast'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%

% MODIFIER: performs spatial re-weighting for multivariate models

% "secret" extra function analysis models --> testing in progress
%
% =========================================================================
%      analysis_model: 'FALFF'   (fractional amplitude of low-frequency fluctuations)
%    * regional activity in bold power band *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.spm       = string specifying format of output SPM. Options include: 
%                                         'raw'  : map of fractional energy
%                                         'zscore': Z-scored map of reproducible correlation values 
%
%                  split_info.type      = 'nocontrast'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'GCONN'   (global connectivity analysis)
%    * average correlation of each voxel with all other brain regions *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.conn_metric = string giving type of connectivity metric. Options include: 
%                                         'cov'  : covariance
%                                         'corr' : correlation
%                                         'pcorr': partial correlation
%
%                  split_info.spm       = string specifying format of output SPM. Options include: 
%                                         'raw'  : map of voxelwise seed correlations  
%                                         'zscore': Z-scored map of reproducible correlation values 
%
%                  split_info.type      = 'nocontrast'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
% =========================================================================
%      analysis_model: 'HURST'   (regional estimate of Hurst exponent)
%    * quantifies "long-memory" process in BOLD signal *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.spm       = string specifying format of output SPM. Options include: 
%                                         'raw'  : map of voxelwise seed correlations  
%                                         'zscore': Z-scored map of reproducible correlation values 
%
%                  split_info.type      = 'nocontrast'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
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

% list of univariate models to exclude
univar_list = {'gnb', 'glm', 'erglm', 'ergnb', 'sconn', 'falff', 'gconn', 'hurst' };

% check that model is NOT on the exclusion list before reweighting
if( sum(strcmpi( analysis_model, univar_list ))==0 )
    %
    % vascular down-weighting to current, preprocessed data
    datamat = bsxfun(@times,datamat,split_info.spat_weight);
end

switch lower( analysis_model )
    %==========================================================================
    case 'lda'
        output = module_LDA( datamat, split_info );
        % --
        output.modeltype = 'one_component';
    case 'gnb'
        output = module_GNB( datamat, split_info );
        % --
        output.modeltype = 'one_component';
    case 'glm'
        output = module_GLM( datamat, split_info );
        % --
        if    ( size(split_info.design_mat,2)==1 ) output.modeltype = 'one_component';
        elseif( size(split_info.design_mat,2) >1 ) output.modeltype = 'multi_component';
        end
    case 'ercva'
        output = module_erCVA( datamat, split_info );
        % --
        if    ( ~isfield( split_info, 'subspace' ) || strcmp(split_info.subspace, 'onecomp')   ) output.modeltype = 'one_component';
        elseif(                                       strcmp(split_info.subspace, 'multicomp') ) output.modeltype = 'multi_component';
        end
    case 'ergnb'
        output = module_erGNB( datamat, split_info );
        % --
        output.modeltype = 'one_component';
    case 'erglm'
        output = module_erGLM( datamat, split_info );
        % --
        output.modeltype = 'one_component';
    case 'sconn'
        output = module_SCONN( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --
    case 'gpca'
        output = module_gPCA( datamat, split_info );
        % --
        output.modeltype = 'multi_component';
        % --
    case 'multi-lda'
        output = module_MultiClassLDA( datamat, split_info );
        % --
        output.modeltype = 'multi_component';
        % --
    case 'svm'
        output = module_SVM( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --
    case 'dpls'
        output = module_dPLS( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --
    case 'falff'
        output = module_falff( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --        
    case 'gconn'
        output = module_gconn( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --        
    case 'hurst'
        output = module_hurst( datamat, split_info );
        % --
        output.modeltype = 'one_component';
        % --        
        %==========================================================================        
    otherwise
        error('Model not on the list! Please specify another.');
end


