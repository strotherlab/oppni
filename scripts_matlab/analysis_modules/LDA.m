function output = LDA( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_LDA: module that performs linear discriminant analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_LDA( datamat, split_info )
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
output.attributes.model_name  = 'LDA';
output.attributes.design_type = 'block';
output.attributes.model_type  = 'multivariate';
output.attributes.num_comp    = 'one_component';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
    if( ~isfield(split_info{1},'drf') || isempty(split_info{1}.drf) )
        disp('LDA uses default data reduction drf=0.5');
        split_info{1}.drf = 0.5;
    end
%-------------------------------------------------------------------------%
    
%% INDIVIDUAL SUBJECT ANALYSIS
if( ~iscell(datamat) || length(datamat)==1 )
    
    % in cases where it is a single cell, revert to single-subject anaylsis
    if( iscell(datamat) && length(datamat)==1 ) datamat = datamat{1}; end
    split_info = split_info{1}; %% take entry from single cell
    
    % define 4 different task blocks for analysis:
    block_cond1_sp1 = datamat(:,split_info.idx_cond1_sp1);
    block_cond1_sp2 = datamat(:,split_info.idx_cond1_sp2);
    block_cond2_sp1 = datamat(:,split_info.idx_cond2_sp1);
    block_cond2_sp2 = datamat(:,split_info.idx_cond2_sp2);

    % define split-half task designs:
    design1 = [ ones( size(block_cond1_sp1,2), 1 ); -ones( size(block_cond2_sp1,2), 1 )];
    design2 = [ ones( size(block_cond1_sp2,2), 1 ); -ones( size(block_cond2_sp2,2), 1 )];
    % linear discriminant analysis, under single-split structure
    results = lda_optimization( [block_cond1_sp1 block_cond2_sp1], [block_cond1_sp2 block_cond2_sp2], design1,design2, split_info.drf );

    % Euclid. distance from (P=1,R=1)
    DD = sqrt( (1-results.R).^2 + (1-results.P).^2 );
    % select PC subspace that minimizes D(P,R)
    [vd id]  = min(DD);

    % [Record optimal statistics + eigenimages]
    %
    output.metrics.R    =  results.R(id);
    output.metrics.P    =  results.P(id);
    output.metrics.Acc  =  results.Acc(id); % alt: fractional classif. accuracy
    output.metrics.dPR  = -vd;
    % optimal eigenimage
    output.images  = results.eig(:,id);

    % [CV scores]
    %
    % CV score timeseries, from reference eigenimage
    output.temp.CV_ref = results.CV(:,id);
    % CV score timeseries, on unit-normed rSPM eigenimage
    output.temp.CV_alt = datamat' * (output.images ./ sqrt(sum(output.images.^2)));

    % [Fractional Variance Explained by eigenimage basis]
    %
    % the scaled projection
    svect = output.temp.CV_alt;
    % and normed projection
    uvect = svect ./ sum(svect.^2);
    % get back out the scaling factor relative to normed eig
    svar = var( svect ) ./ var( uvect );
    % total data variance
    tvar = trace( datamat'*datamat );
    % fraction
    output.temp.CV_alt_varfract = svar ./ tvar;

%% GROUP LEVEL ANALYSIS
else
    
    % number of subjects
    N_subject = length(datamat);
    % split into cell arrays
    for( n=1:N_subject )
        %
        data_class1{n} = datamat{n}(:, split_info{n}.idx_cond1);
        data_class2{n} = datamat{n}(:, split_info{n}.idx_cond2);    
    end
    % run lda
    results = lda_optimization_group ( data_class1, data_class2, split_info{1}.drf, Resampling_Index );

    % optimization stats
    DD       = sqrt( (1-results.R).^2 + (1-results.P).^2 );
    % select PC subspace that minimizes D(P,R)
    [vd id]  = min(DD);

    % [Record optimal statistics + eigenimages]
    %
    output.metrics.R    =  results.R(id);
    output.metrics.P    =  results.P(id);
    output.metrics.dPR  = -vd;
    % optimal eigenimage
    output.images       = results.eig(:,id);

    % [CV scores]
    %
    % CV score timeseries, from reference eigenimage

    % CV score timeseries, on unit-normed rSPM eigenimage
    for(is=1:N_subject) 
        %
        output.temp.CV_alt{is} = datamat{is}' * (output.images ./ sqrt(sum(output.images.^2)));    
        %
        %--- now get fractional explained variance ---
        %
        % the scaled projection
        svect = output.temp.CV_alt{is};
        % and normed projection
        uvect = svect ./ sum(svect.^2);
        % get back out the scaling factor relative to normed eig
        svar = var( svect ) ./ var( uvect );
        % total data variance
        tvar = trace( datamat{is}'* datamat{is} );
        % fraction
        output.temp.CV_alt_varfract{is} = svar ./ tvar;
    end
end