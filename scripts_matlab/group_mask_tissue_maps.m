function group_mask_tissue_maps( inputfile, newmaskname, reference_file, WARP_TYPE )
%
%==========================================================================
% GROUP_MASK_TISSUE_MAPS: second part of spatial normalization process,
% generating a consensus brain mask and statistics on quality of transforms
%==========================================================================
%
% SYNTAX:    
%
%   group_mask_tissue_maps( inputfile, newmaskname )
%            
% INPUT:
%
%   inputfile  = string specifying "input" textfile (path+name), containing
%                input subject information
%   newmaskname= string specifying name/path of the new group-level 
%                consensus mask being produced as output
%                * If running OPPNI, leave it empty *
%                  Default newmaskname=[] puts it in first output folder specified in "inputfile" 
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
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%
if( exist('OCTAVE_VERSION','builtin') )
    % load stats packages
    pkg load statistics;
    pkg load optim; 
end
   
global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('group_mask_tissue_maps.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end
if ~isdeployed
    addpath(CODE_PATH);
    addpath([CODE_PATH 'NIFTI_tools']);
    addpath([CODE_PATH 'toolbox'])
end

if nargin < 2
     newmaskname = [];
end

[~, reference_prefix, ~]=fileparts( reference_file );

%% Load individual masks

[InputStruct, MULTI_RUN_INPUTFILE] = Read_Input_File(inputfile);

% check if multi-run -- aligned appended
if MULTI_RUN_INPUTFILE
   aligned_suffix= '_aligned';
else
    aligned_suffix = [];
end
% group mask directory
if( isempty(newmaskname) )
    %%% create directory + new mask files, in first output dir.
    newpath = [InputStruct(1).run(1).Output_nifti_file_path, '/GroupMasks',reference_prefix,'-',WARP_TYPE];
    mkdir_r(newpath);
    newmaskname = [newpath '/group'];                
else
    [apath,aprefix,aext] = fileparts(newmaskname);
    apath = [apath,'_',reference_prefix,'-',WARP_TYPE];
    mkdir_r(apath);
    %%% make sure mask-name terminates with .nii string
    newmaskname = [apath,'/',aprefix];
end

for(is=1:numel(InputStruct))
    % based on run-1 of each subject, load the mask
    xbase_string = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/',reference_prefix,'-',WARP_TYPE,'/', InputStruct(is).run(1).Output_nifti_file_prefix '_mask_sNorm.nii'];
    MX           = load_untouch_nii(  xbase_string );
    % save into 4D matrix
    maskSet(:,:,:,is) = double(MX.img);
    maskVct(:,is)     = double(MX.img(:)); %% vectorizing to compare
end

%% Estimate consensus mask

J = jaccard_ovl( maskVct, 0 );
for(w=1:size(maskVct,2)) 
    % dropping self-overlap for each mask
    Jtmp    = J(w,:);
    Jtmp(w) = [];
    J2(w,:) = Jtmp;
end
% median overlap, converted to distance
jmed_dist = 1-median(J2,2);
[thr_j]   = quick_gamma_test( jmed_dist );
% index within-threshold maps
cens_j  = (jmed_dist < thr_j);
% get fraction of subjects included in mask, after discarding outliers
maskFract = sum(maskSet(:,:,:,cens_j),4)./sum(cens_j>0);

% get threshold --> more than 50% of subjects must include
consensus_mask = double( maskFract > 0.50 );
maskFract_vect = maskFract(consensus_mask>0); %% for output
%
nii     = MX; % copy nifti struct
nii.img = consensus_mask; % replace volume
nii.hdr.dime.dim(5) = 1; % adjust for #timepoints
%
save_untouch_nii(nii,[newmaskname '_consensus_mask.nii']);  

%% Load transformed data

for(is=1:numel(InputStruct))

    xbase_string  = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm/',reference_prefix,'-',WARP_TYPE,'/afni_processed/' InputStruct(is).run(1).Output_nifti_file_prefix '_baseproc' aligned_suffix '_sNorm.nii'];
    VX     = load_untouch_nii(  xbase_string );
    vxmat  = nifti_to_mat(VX,nii);

    % collecting mean volume
    meanset(:,is) = mean(vxmat,2);

    % formatting for tissue estimation in new normalized space....
    out = GLM_model_fmri( vxmat, 1,[],[], 'econ' ); % just noise
    volcell{1}    = out.vol_denoi(:,1:ceil(size(vxmat,2)/2));
    volcell{2}    = out.vol_denoi(:,ceil(size(vxmat,2)/2)+1:end);

    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------

    % re-computation of tissue priors
    
    % INFO: for physiological (vascular) downweighting maps --> need TRs
    [split_info_struct] = Parse_Split_Info(InputStruct(is).run(1).split_info_file);
    %
    dataInfo.TR            = split_info_struct.TR_MSEC/1000;
    dataInfo.FreqCut       = 0.10;
    dataInfo.thresh_method = 'noprior';
    dataInfo.out_format    = 0; %% don't make output volumes
    % estimate vascular map
    outwt = PHYCAA_plus_step1( volcell, dataInfo );
    %white matter tissue mask
    outwm = WM_weight( volcell, dataInfo );

    % GROUP: get priors into 2d matrices
    NN_weight_set(:,is) = outwt.NN_weight;
    WM_weight_set(:,is) = outwm.WM_weight;    
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
end

%% Summary statistics on transformed data

% ------- volume stats:

% MEAN: correlations
C = corr( meanset );
for(w=1:numel(InputStruct)) 
    %
    Ctmp    = C(w,:);
    Ctmp(w) = [];
    C2(w,:) = Ctmp;
end
% median correlation, converted to distance
cmed_mean       = 1-median(C2,2);
[thr_mean] = quick_gamma_test( cmed_mean );

% NN: correlations
C = corr( NN_weight_set );
for(w=1:numel(InputStruct)) 
    %
    Ctmp    = C(w,:);
    Ctmp(w) = [];
    C2(w,:) = Ctmp;
end
% median correlation, converted to distance
cmed_nn          = 1-median(C2,2);
[thr_nn] = quick_gamma_test( cmed_nn );

% WM: correlations
C = corr( WM_weight_set );
for(w=1:numel(InputStruct)) 
    %
    Ctmp    = C(w,:);
    Ctmp(w) = [];
    C2(w,:) = Ctmp;
end
% median correlation, converted to distance
cmed_wm          = 1-median(C2,2);
[thr_wm] = quick_gamma_test( cmed_wm );

% -------- robust averaging

MEAN_avg      = mean( meanset,2 );
NN_weight_avg = NN_group_average( NN_weight_set );
WM_weight_avg = NN_group_average( WM_weight_set );

% -------- recording summary results

% Saving to volumes: concatenate
tmp1=double(consensus_mask); tmp1(tmp1>0) = MEAN_avg;
tmp2=double(consensus_mask); tmp2(tmp2>0) = 1-NN_weight_avg;
tmp3=double(consensus_mask); tmp3(tmp3>0) = 1-WM_weight_avg;
%
nii=VX; % copy nifti struct
nii.img = cat(3,tmp1,tmp2,tmp3); % replace volume
nii.hdr.dime.datatype = 16;
nii.hdr.hist = VX.hdr.hist;
nii.hdr.dime.dim(5) = 3;
%
save_untouch_nii(nii,[newmaskname, '_mean_NN_WM.nii']);  

% subject statistics on similarity of maps
volume_stats.MEAN_distance = cmed_mean;
volume_stats.MEAN_outlier_thresh = thr_mean;
volume_stats.NN_distance = cmed_nn;
volume_stats.NN_outlier_thresh = thr_nn;
volume_stats.WM_distance = cmed_wm;
volume_stats.WM_outlier_thresh = thr_wm;
volume_stats.MASK_distance = jmed_dist;
volume_stats.MASK_outlier_thresh = thr_j;

% mean and variability of maps (without robustifying average)
brain_volumes.MEAN_avg   = mean(meanset,2);
brain_volumes.MEAN_std   = std(meanset,0,2);
brain_volumes.NN_avg     = mean(NN_weight_set,2);
brain_volumes.NN_std     = std(NN_weight_set,0,2);
brain_volumes.WM_avg     = mean(WM_weight_set,2);
brain_volumes.WM_std     = std(WM_weight_set,0,2);
brain_volumes.MASK_fract = maskFract_vect;

% save results to matfile
     % matlab-compatible
save([newmaskname,'_spat_norm_qc.mat'],'brain_volumes','volume_stats', '-v7');


%%
function [thresh] = quick_gamma_test( dist_vector )

% vectorize
dist_vector = dist_vector(:);
%
ratiotest = range(dist_vector)/min(dist_vector);

if( exist('OCTAVE_VERSION','builtin') && (ratiotest < (1/5)) )
    % in case differences are v. small (+run in octave), default out:
    thresh = max(dist_vector) + eps;
else
    % Gamma fit on correlation distance values
    par_ab   = gamfit( dist_vector );
    % threshold of max dist.
    thresh   = gaminv( 0.95, par_ab(1), par_ab(2) );
end

