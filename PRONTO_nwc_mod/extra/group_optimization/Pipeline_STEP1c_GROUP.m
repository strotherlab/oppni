function Pipeline_STEP1c_GROUP( inputfile, newmaskname )
%
%==========================================================================
% PIPELINE_STEP1C_GROUP: third piece of group-level preprocessing pipeline,
% generates a consensus group-level brain mask, and outputs some diagnostic
% information about the quality of between-subject alignment.
%==========================================================================
%
% SYNTAX:    
%
%   Pipeline_STEP1c_GROUP( inputfile, newmaskname )
%            
% INPUT:
%
%   inputfile  = string specifying "input" textfile (path+name), containing
%                input subject information
%   newmaskname= string specifying name/path of the new group-level 
%                consensus mask being produced as output
%
% OUTPUT:  
%
%         > 3D binary volume, giving consensus group mask, 'newmaskname'
%         > set of plots, examining consistency of subjects' fMRI data 
%           after spatial normalization
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: 2013/08/15
% ------------------------------------------------------------------------%

if( exist('OCTAVE_VERSION','builtin') )
    % load stats packages
    pkg load statistics;
    pkg load optim; 
end
 
addpath scripts_matlab;
addpath scripts_matlab/NIFTI_tools;
   
% opens the inputfile (includes subject/dataset names that preprocessing is performed on...
fid   = fopen(inputfile);
tline = fgetl(fid);
ksub  = 0;
while ischar(tline) % for every input line in textfile...

    ksub=ksub+1; % index #subjects(/datasets)
       
    % parse output directory
    ifile = strfind( tline, 'OUT=' ); ifile = ifile+4;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>ifile);
    fullline = tline(ifile:ips(1));
    isepr    = strfind( fullline, '/' );
    isepr    = isepr(end);
    prefix   = fullline(isepr+1:end);
    outdir   = fullline(1:isepr-1);
    
    xbase_string  = strcat(outdir,'/spat_norm/',prefix,'_mask_sNorm.nii');
    MX     = load_untouch_nii(  xbase_string );
    %
    maskSet(:,:,:,ksub) = double(MX.img);
    
    if(ksub==1) [Nx Ny Nz] = size(MX.img); end
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    
    % then get next line (subject)...
    tline = fgetl(fid);
end

fclose(fid);

% 
N_subject = ksub;
%
for(w=1:N_subject) 
    vol          = maskSet(:,:,:,w);
    maskvct(:,w) = vol(:);
end
%
J = jaccard_ovl( maskvct, 0 );
for(w=1:N_subject) 
    %
    Jtmp    = J(w,:);
    Jtmp(w) = [];
    J2(w,:) = Jtmp;
end

% median overlap
jmed      = median(J2,2);
% convert into distance measure
jmed_dist = 1 - jmed;
% Gamma fit on correlation distance values
par_ab    = gamfit( jmed_dist );
% get Gamma probability of each point + outlier threshold
probGam   = 1 - gamcdf( jmed_dist, par_ab(1), par_ab(2) );
% threshold of min. jaccard overlap
jmed_thr  = 1 - gaminv( 0.95, par_ab(1), par_ab(2) );
%
keepindex = find( probGam > 0.05 );
N_adjust  = length( keepindex );
maskFract = sum(maskSet(:,:,:,keepindex),4)./N_adjust;
%

% ----------- (Fig.0) Display Temporal Variance information ----------- %

% Only plot output if operating in matlab
if( exist('OCTAVE_VERSION','builtin') )
    %
    disp(['figures not plotted in Octave']);
else
    h0=figure;%('visible','off');
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.10 0.15 0.60 0.65]);

    stdvol = maskFract;
    slcx = permute( stdvol(round(Nx/2),:,:), [3 2 1] );
    slcy = permute( stdvol(:,round(Ny/2),:), [3 1 2] );
    slcz = stdvol(:,:,round(Nz/2));

    subplot(2,3,1); imagesc( flipud(slcx), [0 1] );
    subplot(2,3,2); imagesc( flipud(slcy), [0 1] );
    title('Fraction of subjects included in mask (=brain tissue)');
    subplot(2,3,3); imagesc( slcz, [0 1] );
    %
    subplot(2,1,2); boxplot( J2' ); hold on;
    plot( [0.5, N_subject+0.5], [1 1].*(jmed_thr), '--r');
    xlabel('subjects');
    ylabel('Jaccard overlap');
    title('average overlap of individual subjects mask with all others (with 95% CIs)');
end

% get threshold
consensus_mask = double( maskFract > 0.50 );
%
nii     = MX; % copy nifti struct
nii.img = consensus_mask; % replace volume
nii.hdr.dime.dim(5) = 1; % adjust for #timepoints
%
save_untouch_nii(nii,newmaskname);  

%%
   
% opens the inputfile (includes subject/dataset names that preprocessing is performed on...
fid   = fopen(inputfile);
tline = fgetl(fid);
ksub  = 0;
while ischar(tline) % for every input line in textfile...

    ksub=ksub+1; % index #subjects(/datasets)
       
    % parse output directory
    ifile = strfind( tline, 'OUT=' ); ifile = ifile+4;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>ifile);
    fullline = tline(ifile:ips(1));
    isepr    = strfind( fullline, '/' );
    isepr    = isepr(end);
    prefix   = fullline(isepr+1:end);
    outdir   = fullline(1:isepr-1);
    
    xbase_string  = strcat(outdir,'/',prefix,'_baseproc_sNorm.nii');
    VX     = load_untouch_nii(  xbase_string );
    vxmat  = nifti_to_mat(VX,nii);
    [Nvox Ntime] = size(vxmat);
    %
    meanset(:,ksub) = mean(vxmat,  2);
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    
    % then get next line (subject)...
    tline = fgetl(fid);
end

fclose(fid);

% put all subjects on same overall scale (global norming)
meanset = meanset ./ repmat( sqrt(sum(meanset.^2)), [size(meanset,1) 1] );

cc_mean = corrcoef( meanset );

for(w=1:N_subject) 
    %
    cc_temp       = cc_mean(w,:);
    cc_temp(w)    = [];
    cc_mean2(w,:) = cc_temp; 
end

% median overlap
cmed      = median(cc_mean2,2);
% convert into distance measure
cmed_dist = 1 - cmed;
% Gamma fit on correlation distance values
par_ab    = gamfit( cmed_dist );
% get Gamma probability of each point + outlier threshold
probGam   = 1 - gamcdf( cmed_dist, par_ab(1), par_ab(2) );
% threshold of min. jaccard overlap
cmed_thr  = 1 - gaminv( 0.95, par_ab(1), par_ab(2) );

% Only plot output if operating in matlab
if( exist('OCTAVE_VERSION','builtin') )
    %
    disp(['figures not plotted in Octave']);
else

    h1=figure;%('visible','off');
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.20 0.05 0.50 0.85]);

    stdvol = consensus_mask;
    stdvol(stdvol>0) = mean(meanset,2);
    slcx = permute( stdvol(round(Nx/2),:,:), [3 2 1] );
    slcy = permute( stdvol(:,round(Ny/2),:), [3 1 2] );
    slcz = stdvol(:,:,round(Nz/2));

    subplot(3,3,1); imagesc( flipud(slcx), prctile( mean(meanset,2), [1 99] ) );
    subplot(3,3,2); imagesc( flipud(slcy), prctile( mean(meanset,2), [1 99] ) );
    title('Average of mean BOLD images, across subjects');
    subplot(3,3,3); imagesc( slcz, prctile( mean(meanset,2), [1 99] ) );

    stdvol = consensus_mask;
    stdvol(stdvol>0) = std(meanset,0,2)./mean(meanset,2);
    slcx = permute( stdvol(round(Nx/2),:,:), [3 2 1] );
    slcy = permute( stdvol(:,round(Ny/2),:), [3 1 2] );
    slcz = stdvol(:,:,round(Nz/2));

    subplot(3,3,4); imagesc( flipud(slcx), prctile( std(meanset,0,2)./mean(meanset,2), [1 99] ) );
    subplot(3,3,5); imagesc( flipud(slcy), prctile( std(meanset,0,2)./mean(meanset,2), [1 99] ) );
    title('StDev of mean BOLD images, across subjects');
    subplot(3,3,6); imagesc( slcz, prctile( std(meanset,0,2)./mean(meanset,2), [1 99] ) );
    %
    subplot(3,1,3); boxplot( cc_mean2' ); hold on;
    plot( [0.5, N_subject+0.5], [1 1].*(cmed_thr), '--r');
    xlabel('subjects');
    ylabel('Correlation');
    title('average correlation of individual subjects mean volume with all others (with 95% CIs)');
end

% record some outputs
mask_ovl     = jmed;
meanvol_corr = cmed;
%
outlier_mask    = ( jmed < jmed_thr );
outlier_meanvol = ( cmed < cmed_thr );

% save results to matfile
if( exist('OCTAVE_VERSION','builtin') )
     % matlab-compatible
     save([newmaskname(1:end-4),'_outputs.mat'],'maskFract','mask_ovl','meanvol_corr','outlier_mask','outlier_meanvol', '-mat7-binary');
else save([newmaskname(1:end-4),'_outputs.mat'],'maskFract','mask_ovl','meanvol_corr','outlier_mask','outlier_meanvol');
end
