function output = group_analyses_wrapper( datamat, split_info, analysis_model )
%
% =========================================================================
% GROUP_ANALYSES_WRAPPER: general "wrapper" script for fMRI analysis, to be
% used in combination with Pipeline_STEP2_MM_GROUP.m
% Allows users to specify different analysis models, and can substitute
% them as "modules".
% This is the wrapper for multi-subject analysis. For single-subject
% analyses, see "run_analyses_wrapper"!
% =========================================================================
%
% Syntax:   output = group_analyses_wrapper( datamat, split_info, analysis_model )
%  Input:
%
%          datamat : cell array, where each entry corresponds to a subject.
%                    Each entry is a (voxels x time) fmri data matrix
%                    you can have different #timepoints for different subjects
%       split_info : cell array, where each entry corresponds to a subject.
%                    Each entry is a structure with fields required to run analysis modules.
%                    The required fields for current modules are listed below
%   analysis_model : string specifying chosen analysis model. Options also listed below.
%
% Output: 
%           output : structure with following fields
%
%                      output.metrics  : structure containing performance metrics
%                      output.images   : (voxels x components) matrix of brain maps
%                      output.temp     : cell array, where each entry corresponds to a subject.
%                                        each entry is (time x components) matrix of subject timeseries
%                                        expressin the component brain maps
%                      output.modeltype: string specifying either
%                                          'onecomp'   = single component image
%                                          'multicomp' = many component images
%
% =========================================================================
%      analysis_model: 'LDA'   (linear discriminant)
%    * multivariate analysis for block design, with 2 task conditions *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.task1 = vector of volumes indice for task condition 1
%                  split_info.task2 = vector of volumes indice for task condition 2
%                  split_info.drf   = scalar value of range (0,1), indicating the fraction of full-date
%                                     PCA subspace to keep during PCA-LDA analysis
%             split_info.N_resample = scalar indicating the number of resampling splits
%                                     to perform (and average over)
%                  split_info.type  = 'block'  (this needs to be a fixed parameter!)
%
% =========================================================================
%   analaysis_model: 'gCCA' (generalized canonical correlation analysis)
%      Fully exploratory, unsupervised multivariate analysis technique
%       
%         split_info.drf   = scalar value of range (0,1), indicating the fraction of signle subject-date
%                                     PCA subspace to keep during PCA-LDA analysis
%         split_info.N_resample = scalar indicating the number of resampling splits
%         split_info.Spatial_downsampling = scalar indicating spatial down
%                                       sampling rate 
%                                     (if not provided is 1 which the algorithm includes all voxels in calculation of MI)
%         split_info.k  = k in the k-nearest neighbor MI estimation
%                                      (default 27)
% -------------------------------------------------------------------------

%
% =========================================================================
% GROUP_ANALYSES_WRAPPER: general "wrapper" script for fMRI analysis, to be
% used in combination with Pipeline_STEP2_MM_GROUP.m
% Allows users to specify different analysis models, and can substitute
% them as "modules".
% This is the wrapper for multi-subject analysis. For single-subject
% analyses, see "run_analyses_wrapper"!
% =========================================================================
%
% Syntax:   output = group_analyses_wrapper( datamat, split_info, analysis_model )
%  Input:
%
%          datamat : cell array, where each entry corresponds to a subject.
%                    Each entry is a (voxels x time) fmri data matrix
%                    you can have different #timepoints for different subjects
%       split_info : cell array, where each entry corresponds to a subject.
%                    Each entry is a structure with fields required to run analysis modules.
%                    The required fields for current modules are listed below
%   analysis_model : string specifying chosen analysis model. Options also listed below.
%
% Output: 
%           output : structure with following fields
%
%                      output.metrics  : structure containing performance metrics
%                      output.images   : (voxels x components) matrix of brain maps
%                      output.temp     : cell array, where each entry corresponds to a subject.
%                                        each entry is (time x components) matrix of subject timeseries
%                                        expressin the component brain maps
%                      output.modeltype: string specifying either
%                                          'onecomp'   = single component image
%                                          'multicomp' = many component images
%
% =========================================================================
%      analysis_model: 'LDA'   (linear discriminant)
%    * multivariate analysis for block design, with 2 task conditions *
% -------------------------------------------------------------------------
% 'split_info' fields:
%
%                  split_info.idx_cond1 = vector of volumes indices for task condition 1
%                  split_info.idx_cond2 = vector of volumes indices for task condition 2
%                  split_info.drf   = scalar value of range (0,1), indicating the fraction of full-date
%                                     PCA subspace to keep during PCA-LDA analysis
%             split_info.N_resample = scalar indicating the number of resampling splits
%                                     to perform (and average over)
%
%                  split_info.type      = 'block'  (this needs to be a fixed parameter!)
%                  split_info.TR_MSEC   = integer specifying TR (acquisition rate) in milliseconds
%
%  ** {drf, N_resample, type} fields are only required in cell split_info{1}
%
% =========================================================================
%      analysis_model: 'erCVA'  (event-related Canonical Variates Analysis)
%    * predictive multivariate analysis for event-related design, with 1 task type onset *
% -------------------------------------------------------------------------
% 'split_info' fields:
% 
%                  split_info.onsetlist = vector list of integers representing stimulus onset times
%                                         in milliseconds, zero-relative to the first TR volume
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
%  ** {Nblock, WIND, drf, subspace, type} fields are only required in cell split_info{1}
%

%==========================================================================

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

% MODIFIER: performs spatial re-weighting for multivariate models

% list of univariate models to exclude
univar_list = {'gnb', 'glm', 'erglm', 'ergnb', 'sconn' };

% check that model is NOT on the exclusion list before reweighting
if( sum(strcmpi( analysis_model, univar_list ))==0 )
    
    for( is=1:length(datamat) )
        % vascular down-weighting applied to current, preprocessed data
        datamat{is} = bsxfun(@times,datamat{is},split_info{1}.spat_weight);  % denoised volume x weight
    end
end

N_subject = numel(datamat);
if ~isfield(split_info{1},'N_resample')
    N_resample = 10;
else
    N_resample = split_info{1}.N_resample;
end
[Resampling_Index] = generate_split_half_list(N_subject,N_resample);
N_resample = size(Resampling_Index,1);
 

switch lower( analysis_model )
    %==========================================================================
    case 'lda'
        output = module_LDA_group( datamat, split_info, Resampling_Index );
        % --
        output.modeltype = 'one_component';
    case 'gnb'
        output = module_GNB_group( datamat, split_info, Resampling_Index );
        % --
        output.modeltype = 'one_component';
    case 'glm'
        output = module_GLM_group( datamat, split_info, Resampling_Index );
        % --
        if    ( size(split_info.design_mat,2)==1 ) output.modeltype = 'one_component';
        elseif( size(split_info.design_mat,2) >1 ) output.modeltype = 'multi_component';
        end
    case 'ercva'
        output = module_erCVA_group( datamat, split_info );
        % --
        if    ( ~isfield( split_info, 'subspace' ) || strcmp(split_info.subspace, 'onecomp')   ) output.modeltype = 'one_component';
        elseif(                                       strcmp(split_info.subspace, 'multicomp') ) output.modeltype = 'multi_component';
        end
    case 'ergnb'
        output = module_erGNB_group( datamat, split_info );
        % --
        output.modeltype = 'one_component';
    case 'erglm'
        output = module_erGLM_group( datamat, split_info, Resampling_Index );
        % --
        output.modeltype = 'one_component';
    case 'sconn'
        output = module_SCONN_group( datamat, split_info, Resampling_Index );
        % --
        output.modeltype = 'one_component';
        % --
    case 'multi-lda'
        
        output = module_MultiClassLDA_group( datamat, split_info, Resampling_Index);
        % --
        output.modeltype = 'multi_component';
        % --
        %==========================================================================
    otherwise
                error('Model not on the list! Please specify another.');
end


