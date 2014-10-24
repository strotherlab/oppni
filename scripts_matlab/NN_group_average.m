function nn_weight_avg = NN_group_average( nn_mat )
%
%==========================================================================
%  NN_GROUP_AVERAGE: performs robust averaging across multiple subjects' 
%  non-neuronal (NN) weighting maps - derived in PHYCAA+ step-1. Used to 
%  get a single, stable map of vasculature/CSF, for multi-subject analyses
%==========================================================================
%
% SYNTAX:
%
%   nn_weight_avg = robust_avg_NNweight( nn_mat )
%
% INPUT:
%
%   nn_mat = (voxels x subjects) matrix, consisting of individuals' NN 
%            weighting maps (output.NN_weight from "PHYCAA_plus_step1.m")
%
% OUTPUT:
%
%   nn_weight_avg = (voxels x 1) vector, giving consensus NN weighting 
%                   values across individual subjects' maps
%
% ------------------------------------------------------------------------
% Notes:  this is a simple procedure that does the following steps:
%
%      (1) identifies and removes outlier subject maps, by fitting a Gamma 
%          distribution on correlation distance between subject NN-maps
%      (2) takes a probabilistically-weighted mean of the NN-maps, based on
%          Gamma-likelihood of the correlation distance
%      (3) rescales the weighted-mean NN map, do have the same fraction of
%          zero-weighted voxels as the individual subject ones 
%      (4) outputs the averaged, rescaled NN-map
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
% version history: 2013/07/21
% ------------------------------------------------------------------------%

disp('Robust averaging of non-neuronal maps...');

% get dimensions of data
[Nvox Ncell] = size( nn_mat    );

% correlation matrix on maps
cc = corr( nn_mat );

% drop diagonal entries (cross-corr only)
for(w=1:Ncell)
    %
    cctmp    = cc(w,:);
    cctmp(w) = [];
    cc2(w,:) = cctmp; 
    %
end

% get the mean cross-correlation distance
cdist_avg = 1 - mean( cc2,2 );
% Gamma fit on correlation distance values
% (distribution on set of strictly-positive values)
par_ab    = gamfit( cdist_avg );
% get Gamma probability of each point + outlier threshold
probGam   = 1-gamcdf( cdist_avg, par_ab(1), par_ab(2) );

disp('Outlier subjects at p=.05 (atypical non-neuronal maps):');
disp(find( probGam(:)' <= 0.05 ) );

% trim outliers
nn_mat    = nn_mat(:, probGam > 0.05 );
probGam   = probGam( probGam > 0.05 );
% net similarity weighting (sum to 1)
probNorm  = probGam ./ sum(probGam);

% quick estimate of average similarity between NN maps
cc_new = triu( corr( nn_mat ), 1);
cc_avg = mean(cc_new(cc_new~=0));
cc_std =  std(cc_new(cc_new~=0));

disp(['Average correlation between non-neuronal maps: (', num2str(cc_avg,2),') +/- (', num2str(cc_std,2),')']);

% get weighted average
nn_weight_avg  = nn_mat * probNorm;
% fraction of voxels ~ zero
zeropct = 100* median( sum( nn_mat < 0.001 )./ size(nn_mat,1) );
% threshold on values, to set =0
mincut  = prctile( nn_weight_avg,  zeropct );
% recenter/rescale to correct for this
nn_weight_avg = (nn_weight_avg - mincut)./(1-mincut);
% zero out negative values
nn_weight_avg(nn_weight_avg<0)=0;
