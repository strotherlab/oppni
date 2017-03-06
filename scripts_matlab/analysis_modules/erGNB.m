function output = erGNB( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_erGNB: module that performs event-related gaussian naive bayes analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_erGNB( datamat, split_info )
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
output.attributes.model_name  = 'erGNB';
output.attributes.design_type = 'event';
output.attributes.model_type  = 'univariate';
output.attributes.num_comp    = 'one_component';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
if( ~isfield(split_info{1},'WIND')   || isempty(split_info{1}.WIND)   ) 
    disp('erGNB uses default window-size WIND=10');
    split_info{1}.WIND =10; 
end
if( ~isfield(split_info{1},'Nblock') || isempty(split_info{1}.Nblock) ) 
    disp('erGNB uses default number of blocks Nblock=4');
    split_info{1}.Nblock=4; 
end
if( ~isfield(split_info{1},'norm')   || isempty(split_info{1}.norm)   ) 
    disp('erGNB uses default normalization turned on (norm=1)');
    split_info{1}.norm  =1; 
end
%-------------------------------------------------------------------------%

params.TR     = split_info.TR_MSEC;
params.delay  = split_info.TR_MSEC./2;
params.WIND   = split_info.WIND;
params.norm   = 1;

%% INDIVIDUAL SUBJECT ANALYSIS
if( ~iscell(datamat) || length(datamat)==1 )
    
    % in cases where it is a single cell, revert to single-subject anaylsis
    if( iscell(datamat) && length(datamat)==1 ) datamat = datamat{1}; end
    split_info = split_info{1}; %% take entry from single cell

    % % % count #leftovers, if we do even multiples of splits
    % % leftover = rem( size(datamat,2) , params.Nblock);
    % % % trim "overhang" from the data matrix
    % % datamat_trim  = datamat(:,1:end-leftover);

    %  Time-windowed averaging of data:
    % *NB: simple averaging does nearest-neighbour interpolation to nearest TR
    %      interval. You can replace with "interp_averaging_for_ER" to do
    %      linear time interpolation; the tends to increase reproducibility,
    %      but at the apparent expense of decreased prediction!
    %

    % blocks pre-specified
    params.Nblock = split_info.Nblock;

    [windowAvg] = simple_averaging_for_ER( datamat, split_info.onsetlist, params );
    % run optimized analysis
    out = ergnb_optimization( windowAvg , split_info.spat_weight);

    % keep signed image for representative spatial map
    output.images           = out.sens_hrf1;
    % get all metrics --> but D( R[signed_1st_hrf],  P[classification] )
    output.metrics.R_global =  out.R_global;
    output.metrics.R_1hrf   =  out.R_hrf1;
    output.metrics.P_class  =  out.P_class;
    output.metrics.P_rms    =  out.P_rms;
    output.metrics.dPR      = -sqrt( (1-out.R_hrf1).^2 + (1-out.P_class).^2 );
    % special format of temporal output
    tempstruct.hrf          = out.HRF1;
    tempstruct.sens         = out.sens_global;
    tempstruct.pmap         = out.pmap;
    output.temp             = tempstruct;

%% GROUP LEVEL ANALYSIS
else

    % number of subjects
    N_subject = length(datamat);
    Nvox      = size(datamat{1},1);

    if( N_subject >=4 ) %% general case -- runs = resampling units

        % Enforce 1 split per subject for group analysis
        params.Nblock = 1; 

        % initialize data matrix
        windowAvg = zeros( Nvox, split_info{1}.WIND, N_subject );
        % load from averaged blocks from cell arrays
        for(n=1:N_subject)
            % perform stimulus-locked averaging on each subject
            [windowAvg(:,:,n)] = simple_averaging_for_ER(datamat{n}, split_info{n}.onsetlist, params );
        end
        
    else %% cases where only 2-3 runs ... need enough for variance estimation within each split (2 samples per class)
        
        % Enforce 2 splits per subject for group analysis
        params.Nblock = 2; 

        % initialize data matrix
        windowAvg = zeros( Nvox, split_info{1}.WIND, 2*N_subject );
        % load from averaged blocks from cell arrays
        for( n=1:N_subject )
            % perform stimulus-locked averaging on each subject
            [windowAvg(:,:,2*n-1:2*n)] = simple_averaging_for_ER( datamat{n}, split_info{n}.onsetlist, params );
        end
        % load from averaged blocks from cell arrays
        for(n=1:N_subject)
            % perform stimulus-locked averaging on each subject
            [windowAvg(:,:,n)] = simple_averaging_for_ER(datamat{n}, split_info{n}.onsetlist, params );
        end
    end

    % % % count #leftovers, if we do even multiples of splits
    % % leftover = rem( size(datamat,2) , params.Nblock);
    % % % trim "overhang" from the data matrix
    % % datamat_trim  = datamat(:,1:end-leftover);

    %  Time-windowed averaging of data:
    % *NB: simple averaging does nearest-neighbour interpolation to nearest TR
    %      interval. You can replace with "interp_averaging_for_ER" to do
    %      linear time interpolation; the tends to increase reproducibility,
    %      but at the apparent expense of decreased prediction!
    %
    % run optimized analysis
    out = ergnb_optimization( windowAvg , split_info{1}.spat_weight);

    % keep signed image for representative spatial map
    output.images           = out.sens_hrf1;
    % get all metrics --> but D( R[signed_1st_hrf],  P[classification] )
    output.metrics.R_global =  out.R_global;
    output.metrics.R_1hrf   =  out.R_hrf1;
    output.metrics.P_class  =  out.P_class;
    output.metrics.P_rms    =  out.P_rms;
    output.metrics.dPR      = -sqrt( (1-out.R_hrf1).^2 + (1-out.P_class).^2 );
    % special format of temporal output
    tempstruct.hrf          = out.HRF1;
    tempstruct.sens         = out.sens_global;
    tempstruct.pmap         = out.pmap;
    output.temp             = tempstruct;
end
