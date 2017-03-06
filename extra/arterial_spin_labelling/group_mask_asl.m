function group_mask_tissue_maps( InputStruct, newmaskname )
%

if( exist('OCTAVE_VERSION','builtin') )
    % load stats packages
    pkg load statistics;
    pkg load optim; 
end
   
global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('group_mask_asl.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end
if ~isdeployed
    addpath(CODE_PATH);
    addpath([CODE_PATH 'NIFTI_tools']);
    addpath([CODE_PATH 'toolbox'])
end

%% Load individual masks
% parse input structure if not already done

[InputStruct] = Read_Input_ASL(InputStruct);

% group mask directory
if( isempty(newmaskname) )
    %%% create directory + new mask files, in first output dir.
    newpath = [InputStruct(1).run(1).Output_nifti_file_path, '/asl_processed/GroupMasks'];
    mkdir_r(newpath);
    newmaskname = [newpath '/group'];                
else
    [apath,aprefix,aext] = fileparts(newmaskname);
    mkdir_r(apath);
    %%% make sure mask-name terminates with .nii string
    newmaskname = [apath,'/',aprefix];
end

% % for(s=1:23)
% %     s,
% %     unix(['rm /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Clist(s)),'C/bold_mask_sNorm.nii']);
% %     unix(['rm /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Tlist(s)),'T/bold_mask_sNorm.nii']);
% %     
% %     unix(['3dAutomask -prefix /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Clist(s)),'C/bold_mask_sNorm.nii -clfrac 0.6 /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Clist(s)),'C/proc_BOLD_avg_sNorm.nii']);
% %     unix(['3dAutomask -prefix /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Tlist(s)),'T/bold_mask_sNorm.nii -clfrac 0.6 /Users/nchurchill/Documents/MATLAB/TBI/proc_acute/asl_processed/tbi_',num2str(Tlist(s)),'T/proc_BOLD_avg_sNorm.nii']);
% % end

for(is=1:numel(InputStruct))
    outstr  = [InputStruct(is).run(1).Output_nifti_file_path,'/asl_processed/',InputStruct(is).run(1).Output_nifti_file_prefix{1}];
    xbase_string  = [outstr,'/proc_BOLD_avg_sNorm.nii'];
    ybase_string  = [outstr,'/bold_mask_sNorm.nii'];
    unix(['3dAutomask -prefix ',ybase_string,' -clfrac 0.6 ',xbase_string]);
    MX           = load_untouch_nii(  ybase_string );
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
    outstr  = [InputStruct(is).run(1).Output_nifti_file_path,'/asl_processed/',InputStruct(is).run(1).Output_nifti_file_prefix{1}];
    xbase_string  = [outstr,'/proc_TCBF_avg_sNorm.nii'];
    VX            = load_untouch_nii(  xbase_string );
    meanset(:,is)  = nifti_to_mat(VX,nii);
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

% -------- robust averaging

MEAN_avg      = mean( meanset,2 );

% -------- recording summary results

% Saving to volumes: concatenate
tmp1=double(consensus_mask); tmp1(tmp1>0) = MEAN_avg;
%
nii=VX; % copy nifti struct
nii.img = tmp1; % replace volume
nii.hdr.dime.datatype = 16;
nii.hdr.hist = VX.hdr.hist;
nii.hdr.dime.dim(5) = 1;
%
save_untouch_nii(nii,[newmaskname, '_mean_TCBF.nii']);  

% subject statistics on similarity of maps
volume_stats.MEAN_distance = cmed_mean;
volume_stats.MEAN_outlier_thresh = thr_mean;
volume_stats.MASK_distance = jmed_dist;
volume_stats.MASK_outlier_thresh = thr_j;

% mean and variability of maps (without robustifying average)
brain_volumes.MEAN_avg   = mean(meanset,2);
brain_volumes.MEAN_std   = std(meanset,0,2);
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

