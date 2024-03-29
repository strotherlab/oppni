function output = module_GNB_group( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_GNB: module that performs gaussian naive bayes analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_GNB( datamat, split_info )
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

if( ~isfield(split_info{1},'decision_model') || isempty(split_info{1}.decision_model) )
    disp('GNB uses default linear decision model');
    split_info{1}.decision_model = 'linear';
end

N_resample = size(Resampling_Index,1);
N_subject = length(datamat);
eig = 0;
% split into cell arrays
for( k=1:N_resample )
    %
    
    index = randperm(N_subject);
    set1 = Resampling_Index(k,:);
    set2 = setdiff(1:N_subject,set1);

    data_class1 = [];
    design1     = [];
    for i = 1:length(set1)
        data_class1 = [data_class1 datamat{set1(i)}(:,split_info{set1(i)}.idx_cond1) datamat{set1(i)}(:,split_info{set1(i)}.idx_cond2)];
        design1     = [design1; ones(length(split_info{set1(i)}.idx_cond1),1); -ones(length(split_info{set1(i)}.idx_cond2),1)];
    end
    data_class2 = [];
    design2     = [];
    for i = 1:length(set2)
        data_class2   = [data_class2 datamat{set2(i)}(:,split_info{set2(i)}.idx_cond1) datamat{set2(i)}(:,split_info{set2(i)}.idx_cond2)];
        design2     = [design2; ones(length(split_info{set2(i)}.idx_cond1),1); -ones(length(split_info{set2(i)}.idx_cond2),1)];
    end
    results_temp = gnb_optimization( data_class1, data_class2, design1, design2, split_info{1}.decision_model, split_info{1}.spat_weight );
    eig = eig + results_temp.eig;
    P(k,:) = results_temp.P;
    R(k,:) = results_temp.R;
end
eig = eig /  N_resample;
results.eig = eig;
results.P = median(P,1);
results.R = median(R,1);

    % gnb analysis, under two different splitting structures...

    % define split-half task design

    % analysis

% Euclid. distance from (P=1,R=1)
DD = sqrt( (1-results.R).^2 + (1-results.P).^2 );
% select PC subspace that minimizes D(P,R)
[vd id]  = min(DD);

% [Record optimal statistics + eigenimages]
%
output.metrics.R    =  results.R(id);
output.metrics.P    =  results.P(id);
output.metrics.dPR  = -vd;
% optimal eigenimage
output.images  = results.eig(:,id);

% CV score timeseries, on unit-normed eigenimage
output.temp    = datamat'  * (output.images ./ sqrt(sum(output.images.^2)));
