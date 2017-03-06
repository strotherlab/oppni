close all;
clear;

addpath(genpath('niftimatlib-1.2'));
addpath(genpath('NODDI_example_dataset'));
addpath(genpath('NODDI_toolbox_v0.9'));

addpath(genpath('/Users/Schweizerlab/Documents/Nathan_WorkSpace/OPPNI'));
% CreateROI('NODDI_example_dataset/NODDI_DWI.hdr', 'NODDI_example_dataset/roi_mask.hdr', 'NEW_EXAMPLE/NODDI_roi.mat');
% 
% protocol = FSL2Protocol('NODDI_example_dataset/NODDI_protocol.bval', 'NODDI_example_dataset/NODDI_protocol.bvec');
% noddi    = MakeModel('WatsonSHStickTortIsoV_B0');
% batch_fitting_single('NEW_EXAMPLE/NODDI_roi.mat', protocol, noddi, 'NEW_EXAMPLE/FittedParams.mat');
% SaveParamsAsNIfTI('NEW_EXAMPLE/FittedParams.mat', 'NEW_EXAMPLE/NODDI_roi.mat', 'brain_mask.hdr', 'example')


% subjlist=[5 23]
% delete('noddi_dat/FittedParams.mat');
% 
% for(s=1:length(subjlist))
% 
% gunzip(['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii.gz']);
% gunzip(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.nii.gz']);
% 
% CreateROI(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.nii'], ['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat']);
% protocol = FSL2Protocol(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.bval'], ['tbi_',num2str(subjlist(s)),'C/HARDI_ec.bvec']);
% noddi    = MakeModel('WatsonSHStickTortIsoV_B0');
% batch_fitting_single(['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat'], protocol, noddi, 'noddi_dat/FittedParams.mat');
% SaveParamsAsNIfTI('noddi_dat/FittedParams.mat', ['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat'], ['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'C']);
% delete('noddi_dat/FittedParams.mat');
% end


% subjlist=[10] %8
% % delete('noddi_dat/FittedParams.mat');
% 
% for(s=1:length(subjlist))
% 
%     if( exist(['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii.gz'],'file') )
% gunzip(['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii.gz']);
% gunzip(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.nii.gz']);
%     end
% 
% CreateROI(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.nii'], ['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat']);
% protocol = FSL2Protocol(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.bval'], ['tbi_',num2str(subjlist(s)),'T/HARDI_ec.bvec']);
% noddi    = MakeModel('WatsonSHStickTortIsoV_B0');
% batch_fitting_single(['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat'], protocol, noddi, 'noddi_dat/FittedParams.mat');
% SaveParamsAsNIfTI('noddi_dat/FittedParams.mat', ['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat'], ['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'T']);
% delete('noddi_dat/FittedParams.mat');
% end
% 
% subjlist=[2:7 9] %8
% % delete('noddi_dat/FittedParams.mat');
% 
% for(s=1:length(subjlist))
% 
%     if( exist(['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii.gz'],'file') )
% gunzip(['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii.gz']);
% gunzip(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.nii.gz']);
%     end
% 
% CreateROI(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.nii'], ['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat']);
% protocol = FSL2Protocol(['tbi_',num2str(subjlist(s)),'T/HARDI_ec.bval'], ['tbi_',num2str(subjlist(s)),'T/HARDI_ec.bvec']);
% noddi    = MakeModel('WatsonSHStickTortIsoV_B0');
% batch_fitting_single(['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat'], protocol, noddi, 'noddi_dat/FittedParams.mat');
% SaveParamsAsNIfTI('noddi_dat/FittedParams.mat', ['noddi_dat/tbi_',num2str(subjlist(s)),'T_roi.mat'], ['tbi_',num2str(subjlist(s)),'T/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'T']);
% delete('noddi_dat/FittedParams.mat');
% end


subjlist=[36:38 40]
%delete('noddi_dat/FittedParams.mat');

for(s=1:length(subjlist))

gunzip(['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii.gz']);
gunzip(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.nii.gz']);

M=load_untouch_nii(['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii']);
M.img(1:2:end,:,:)=0;
M.img(:,1:2:end,:)=0;
M.img(:,:,1:2:end)=0;
save_untouch_nii(M,['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask_decime.nii']);

CreateROI(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.nii'], ['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat']);
protocol = FSL2Protocol(['tbi_',num2str(subjlist(s)),'C/HARDI_ec.bval'], ['tbi_',num2str(subjlist(s)),'C/HARDI_ec.bvec']);
noddi    = MakeModel('WatsonSHStickTortIsoV_B0');
batch_fitting_single(['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat'], protocol, noddi, 'noddi_dat/FittedParams.mat');
SaveParamsAsNIfTI('noddi_dat/FittedParams.mat', ['noddi_dat/tbi_',num2str(subjlist(s)),'C_roi.mat'], ['tbi_',num2str(subjlist(s)),'C/HARDI_bet_mask_decime.nii'], ['noddi_dat/tbi_',num2str(subjlist(s)),'C']);
delete('noddi_dat/FittedParams.mat');
end

