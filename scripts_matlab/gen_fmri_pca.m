function gen_fmri_pca( volname, outdir, maskname, normflag )

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

mkdir(outdir);

% load fMRI NIFTI-format data into MatLab
VV = load_untouch_nii(volname);
MM = load_untouch_nii(maskname); mask=double(MM.img);
VV.img = double(VV.img); %% format as double, for compatibility
[nx ny nz nt] = size(VV.img);

volmat = nifti_to_mat( VV,MM );

if( normflag==1 ) volmat = bsxfun( @minus, volmat, mean(volmat,2) ); end
if( normflag==2 ) volmat = bsxfun( @rdivide, volmat, std(volmat,0,2) ); end

[u l v] = svd( volmat, 'econ' );

varfract = diag(l.^2)./trace(l.^2);

Npc = sum( cumsum(varfract) < 0.95 )+1;

u=u(:,1:Npc);
v=v(:,1:Npc);

%%
tcomp = v;
save([outdir '/temporal_PCs'],'tcomp','-ascii');
%%
TMPVOL = zeros( [size(mask) Npc] );
for(n=1:Npc)
   tmp=mask; tmp(tmp>0) = u(:,n);
   TMPVOL(:,:,:,n) = tmp;
end

nii=VV;
nii.img = TMPVOL;
nii.hdr.dime.datatype = 16;
nii.hdr.hist = VV.hdr.hist;
nii.hdr.dime.dim(5) = Npc;
save_untouch_nii(nii,[outdir '/spatial_PCs.nii']);