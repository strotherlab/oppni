function imed_out = min_displace_brick( volname, maskname, outname )
%
% MIN_DISPLACE_BRICK:  This script reads in fMRI 4D dataset and brain mask file.
% It uses these data to find the timepoint with smallest overall displacement,
% measured in Principal Component space. This is used to choose a "reference"
% volume for motion correction, which all other volumes will be matched against. 
%
% Syntax:
%           min_displace_brick( volname, maskname, outname )
% Input:
%           volname  = string, giving path/name of fMRI data (NIFTI or ANALYZE format)
%           maskname = string, giving path/name of binary brain mask (excluding non-brain tissue)
%           outname  = string, giving path/name of output textfile that states minimum-displacement volume
% Output:
%          -creates a text file with name "<outname>_0ref_motbrick.txt"
%           with a single integer. This denotes the timepoint in the input data 
%           (specified by volname) with minimum displacement
%          *note that output is zero-indexed, eg a value of '0' corresponds
%           to timepoint #1 in the fMRI data, etc.
%

% add the path to use code for converting NIFTI files into MatLab
addpath scripts_matlab;
addpath scripts_matlab/NIFTI_tools;

% load fMRI data volumes (in NIFTI format)
VV     = load_untouch_nii(volname);
MM     = load_untouch_nii(maskname);
% convert to 4D fMRI data volume (VV.img) into 2D matrix (voxels x time) for computational purposes
rawepi = convert_nii_to_mat( VV,MM ); 
% subtract the mean from each voxel timeseries, so that this doesn't
% dominate the Principal Component Analysis (PCA)
epimat = rawepi - repmat( mean(rawepi,2), [1 size(rawepi,2)] );
% get fMRI data matrix dimensions
[Nvox Ntime] = size(rawepi);
% singular value decomposition -- gives the Principal Component vectors (V)
% and associated (squared) eigenvalues (S2)
[V S2 temp] = svd( epimat'*epimat ); 
% transform data into PC space (each timepoint = PC-space column vector)
Qdat     = (V * sqrt(S2))';
% compute the median PC-space coordinate (robust measure of the center of the dataset)
Qmed     = median( Qdat,2 );
% measure distnace of each data-point from Qmed
Dist     = sqrt(sum((Qdat - repmat( Qmed, [1 Ntime] )).^2));
% find datapoint closest to Qmed; this is the "most central" datapoint
[v imed] = min( Dist );
% output the index of this datapoint to a text file
dlmwrite( strcat(outname,'_0ref_motbrick.txt'), [imed-1] );
imed_out = imed-1;
%%
function dataMat = convert_nii_to_mat( niiVol, niiMask )
%
% take niiVol and niiMask (nifti) format, and convert
% to matlab vector/matrix structure:
%
% dataMat = nifti_to_mat( niiVol, niiMask )
%
vol = double(niiVol.img);
msk = double(niiMask.img);

dataMat = zeros( sum(msk(:)>0), size(vol,4) );

for(t=1:size(vol,4))
    tmp=vol(:,:,:,t);
    dataMat(:,t) = tmp(msk>0);
end
