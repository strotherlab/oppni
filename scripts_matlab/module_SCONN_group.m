function output = module_SCONN_group( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_GNB: module that performs gaussian naive bayes analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_SCONN( datamat, split_info )
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


% load in seed volume
SS      = load_untouch_nii( split_info.seed_name  ); 
seedvol = double(SS.img > eps);

% get dilated seed for masking
dilvol = seedvol;
for(nz=1:size(seedvol,3)) % dilate each slice
   %
   slc = seedvol(:,:,nz);

   %LMP modified for Octave
   %slc = imdilate(slc,strel('disk',1));
   slc = imdilate(slc,strel('disk',1,0));      
   dilvol(:,:,nz) = slc;
end

% vectorize both
seedvect = seedvol( split_info.mask_vol>0 );
dilvect  =  dilvol( split_info.mask_vol>0 );

N_resample = size(Resampling_Index,1);
N_subject = numel(datamat);
for k = 1:N_resample
    
    index = randperm(N_subject);
    set1 = Resampling_Index(k,:);
    set2 = setdiff(1:N_subject,set1);
    
    datamat_sp1 = [];
    for k = 1:length(set1)
        datamat_sp1 = [datamat_sp1 datamat{set1(k)}];
    end
    
    datamat_sp2 = [];
    for k = 1:length(set2)
        datamat_sp2 = [datamat_sp2 datamat{set2(k)}];
    end
    
    % mean-center splits, to do correlation later
    datamat_sp1 = bsxfun(@minus,datamat_sp1,mean(datamat_sp1,2));
    datamat_sp2 = bsxfun(@minus,datamat_sp2,mean(datamat_sp2,2));
    
    % get seed timeseries
    tseed_1  = mean( datamat_sp1(seedvect>0,:), 1 );
    tseed_2  = mean( datamat_sp2(seedvect>0,:), 1 );
    % re-mean center them
    tseed_1=tseed_1-mean(tseed_1);
    tseed_2=tseed_2-mean(tseed_2);
        
    % corr-map #1
    corr_map1 = sum( bsxfun(@times,datamat_sp1,tseed_1), 2) ./ ( bsxfun(@times,std(datamat_sp1,0,2),std(tseed_1,0,2)) .* ( ceil(Ntime/2) - 1) );
    corr_map1(~isfinite(corr_map1))=0;
    % corr-map #2
    corr_map2 = sum( bsxfun(@times,datamat_sp2,tseed_2), 2) ./ ( bsxfun(@times,std(datamat_sp2,0,2),std(tseed_2,0,2)) .* ( ceil(Ntime/2) - 1) );
    corr_map2(~isfinite(corr_map2))=0;
    
    
    %% map manipulations for results
    
    % reweight the correlations, controlling vasc variance
    corr_map1 = corr_map1 .* split_info.spat_weight;
    corr_map2 = corr_map2 .* split_info.spat_weight;
    % average correlation map
    CORR = (corr_map1 + corr_map2)./2;
    
    % reproducibility and SPM:
    [ R_allvox, rSPMZ_temp ] = get_rSPM( corr_map1, corr_map2, 1 );
    % reproducibility without central seed-ROI
    R_noseed = corr( corr_map1( dilvect==0 ), corr_map2( dilvect==0 ) );
    
    % record optima --> exclude seed ROI, although the impact is generally small
    rSPMZ = rSPMZ + rSPMZ_temp/N_resample;
    R(k) = R_noseed;
end

output.metrics.R =  median(R);
% optimal eigenimage
if( ~isfield(split_info,'spm') )
    disp('SCONN default SPM type is z-scored');
    split_info.spm = 'zscore';
end
if    ( strcmp( split_info.spm, 'zscore' ) ) output.images  = rSPMZ;
elseif( strcmp( split_info.spm, 'raw'   ) ) output.images  =  CORR;
else  error( 'invalid output type for seed correlations (should be zscore or corr).' );
end
% seed timeseries, on unit-normed eigenimage
output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
