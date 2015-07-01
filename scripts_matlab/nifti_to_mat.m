function dataMat = nifti_to_mat( niiVol, niiMask )

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

% take niiVol and niiMask (nifti) format, and convert
% to matlab vector/matrix structure:
%
% dataMat = nifti_to_mat( niiVol, niiMask )
%
msk     = double(niiMask.img);
dataMat = zeros( sum(msk(:)>0), size(niiVol.img,4) );

for(t=1:size(niiVol.img,4))
    tmp=double(niiVol.img(:,:,:,t));
    dataMat(:,t) = tmp(msk>0);
end