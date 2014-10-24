function output = module_LDA_group( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_LDA_GROUP: module that performs group-level linear discriminant analysis 
% in split-half NPAIRS framework
% =========================================================================
%
%   Syntax:
%           output = module_LDA_group( datamat, split_info )
%
%

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
output.metrics.Dneg = -vd;
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
