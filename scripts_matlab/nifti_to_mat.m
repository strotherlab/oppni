function dataMat = nifti_to_mat( niiVol, niiMask )
%
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