function posthoc_analysis( output, test_type, maskfile, files1, files2 )
%
% test_type = 'ttest,1sample'
%             'ttest,2sample'
%             'ttest,paired'
%             'boots,1sample'
%             'boots,2sample'
%             'boots,paired'


if( nargin<5 ) data2=[]; end

ix   = strfind( test_type,',' );
test = test_type(1:ix-1);
type = test_type(ix+1:end);

% load brain mask
M=load_untouch_nii(maskfile);

% load all files into matrices, group 1
for(n=1:length(files1))
    
   V=load_untouch_nii(files1{n});
   tmp = nifti_to_mat( V,M );
   mat1(:,:,n) = tmp;    
end
% load all files into matrices, group 2
if(~isempty(files2))
for(n=1:length(files1))
    
   V=load_untouch_nii(files2{n});
   tmp = nifti_to_mat( V,M );
   mat2(:,:,n) = tmp;    
end
end

mat1 = permute(mat1,[1 3 2]);

if( ~isempty(files2) )

mat2 = permute(mat2,[1 3 2]);

    for(i=1:size(mat1,3))
    if( strcmp(test,'ttest') )
       out = module_ttest( {mat1(:,:,i), mat2(:,:,i)}, type );
       Pmap(:,i) = out.pmap;
       Zmap(:,i) = out.tmap;
    end

    if( strcmp(test,'boots') )
       out = module_bootstrap( {mat1(:,:,i), mat2(:,:,i)}, type );
       Pmap(:,i) = out.pmap;
       Zmap(:,i) = out.bsrmap;
    end
    end

else
    
    for(i=1:size(mat1,3))
    if( strcmp(test,'ttest') )
       out = module_ttest( {mat1(:,:,i)}, type );
       Pmap(:,i) = out.pmap;
       Zmap(:,i) = out.tmap;
    end

    if( strcmp(test,'boots') )
       out = module_bootstrap( {mat1(:,:,i)}, type );
       Pmap(:,i) = out.pmap;
       Zmap(:,i) = out.bsrmap;
    end
    end   
end

[pcrit, th] = fdr( Pmap,'p',0.05,0 );

disp('number of significant voxels');
sum(th),


[opath,oname,oext] = fileparts(output);

% now save files
nii = M;
nii.img = zeros([size(M.img) size(Zmap,2)]);
Z = zeros(size(M.img));
for nvol = 1:size(Zmap,2)
    Z(M.img~=0) = Zmap(:,nvol);
    nii.img(:,:,:,nvol) = Z;
end
nii.hdr.dime.datatype = 16;
nii.hdr.dime.dim([1 5]) = [4 nvol];
if(~isempty(opath))
save_untouch_nii(nii,[opath,'/',oname,'_unthresh.nii']);
else
save_untouch_nii(nii,[oname,'_unthresh.nii']); 
end

% now save files
nii = M;
nii.img = zeros([size(M.img) size(Zmap,2)]);
Z = zeros(size(M.img));
for nvol = 1:size(Zmap,2)
    Z(M.img~=0) = Zmap(:,nvol) .* th(:,nvol);
    nii.img(:,:,:,nvol) = Z;
end
nii.hdr.dime.datatype = 16;
nii.hdr.dime.dim([1 5]) = [4 nvol];
if(~isempty(opath))
save_untouch_nii(nii,[opath,'/',oname,'_FDR.05.nii']);
else
save_untouch_nii(nii,[oname,'_FDR.05.nii']); 
end

