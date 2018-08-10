%Tester
% CHECK ALOMOST EVERYTHING, REQUIES MANUAL FILE INPUT.
% CAN ONLY CHECK DESIGNATED OUTPUT FILES.

%% Intermediate Processed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagnostics
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/processing_in_rec_yng1t_LDA_task_A-baseline/intermediate_processed/diagnostic/session1_ID2382_run3_mc+smo_QC_output.mat')
MCSMO_M = output;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/processing_in_rec_yng1t_LDA_task_A-baseline/intermediate_processed/diagnostic/session1_ID2382_run3_smo_QC_output.mat')
SMO_M = output;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/processing_in_rec_yng1t_LDA_task_A-baseline/Octave_mc.mat')
MCSMO_O = output;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/processing_in_rec_yng1t_LDA_task_A-baseline/Octave_smo.mat')
SMO_O = output;
clear output;

% Split_info
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/processing_in_rec_yng1t_LDA_task_A-baseline/intermediate_processed/split_info/session1_ID2382_run3.mat')
SPLIT_M = split_info;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1ot/processing_in_rec_yng1ot_LDA_task_A-baseline/intermediate_processed/split_info/session1_ID2382_run3.mat')
SPLIT_O = split_info;
clear split_info;


%% Optimization Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spat Norm optimazation summery
%LOADING PROBLEM
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/optimization_results/matfiles/optimization_summary.mat')
OS_M1 = CODE_PROPERTY;
OS_M2 = METRIC_opt;
OS_M3 = pipeline_sets;
OS_M4 = TEMP_opt;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/optimization_results/matfiles/optimization_summary.mat')
OS_O1 = CODE_PROPERTY;
OS_O2 = METRIC_opt;
OS_O3 = pipeline_sets;
OS_O4 = TEMP_opt;
clear CODE_PROPERTY METRIC_opt pipeline_sets TEMP_opt;


%% Intermediate Metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Res0 folder with split into
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/intermediate_metrics/res0_params/params_session1_ID2382_run3.mat')
NNw_M = split_info_set{1,1}.NN_weight_avg;
WNw_M = split_info_set{1,1}.WM_weight_avg;
NNm_M = split_info_set{1,1}.NN_mask_avg;
WMm_M = split_info_set{1,1}.WM_mask_avg;
FXYZa_M = split_info_set{1,1}.FXYZ_avg;
SW_M = split_info_set{1,1}.spat_weight;
MV_M = split_info_set{1,1}.mask_vol;
Xn_M = Xnoise{1,1};
Xs_M = Xsignal{1,1};
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/intermediate_metrics/res0_params/params_session1_ID2382_run3.mat')
NNw_O = split_info_set{1,1}.NN_weight_avg;
WNw_O = split_info_set{1,1}.WM_weight_avg;
NNm_O = split_info_set{1,1}.NN_mask_avg;
WMm_O = split_info_set{1,1}.WM_mask_avg;
FXYZa_O = split_info_set{1,1}.FXYZ_avg;
SW_O = split_info_set{1,1}.spat_weight;
MV_O = split_info_set{1,1}.mask_vol;
Xn_O = Xnoise{1,1};
Xs_O = Xsignal{1,1};

%res1 spm images
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/intermediate_metrics/res1_spms/spms_session1_ID2382_run3.mat')
ImgSet_M = IMAGE_set;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/intermediate_metrics/res1_spms/spms_session1_ID2382_run3.mat')
ImgSet_O = IMAGE_set;

%res2 temp
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/intermediate_metrics/res2_temp/temp_session1_ID2382_run3.mat')
TCV_M = TEMP_set;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/intermediate_metrics/res2_temp/temp_session1_ID2382_run3.mat')
TCV_O = TEMP_set;
clear TEMP_set;

%res 3 stats
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/intermediate_metrics/res3_stats/stats_session1_ID2382_run3.mat')
MsAP_M = METRIC_set;        %ERROR STARTS HERE
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/intermediate_metrics/res3_stats/stats_session1_ID2382_run3.mat')
MsAP_O = METRIC_set;

%% QC1 Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/QC1_results/output_qc1.mat')
QC1_M = output_qc1;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/QC1_results/output_qc1.mat')
QC1_O = output_qc1;
clear output_qc1;



%% QC2 Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/QC2_results/group_pca.mat')
QCG2_M = group_pca;
% load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/QC2_results/group_pca.mat')
% QCG2_O = group_pca;
load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1/processing_in_rec_yng1_LDA_task_A-baseline/QC2_results/output_qc2.mat')
QC2_M = output_qc2;
%load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1o/processing_in_rec_yng1o_LDA_task_A-baseline/QC2_results/output_qc2.mat')
%QC2_O = output_qc2;
clear group_pca output_qc2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Intermediate Processed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split Info
SPLIT_design = SPLIT_M.design_mat - SPLIT_O.design_mat;                         if any(SPLIT_design); disp('SPLIT design diff');SPLIT_design_m= max(max(abs(SPLIT_design))); end
SPLIT_baseline = SPLIT_M.baseline - SPLIT_O.baseline;                           if any(SPLIT_baseline); disp('SPLIT baseline diff');SPLIT_basline_m = max(max(abs(SPLIT_baseline))); end
SPLIT_sp1_1 = SPLIT_M.single.idx_cond(1).sp1 - SPLIT_O.single.idx_cond(1).sp1;  if any(SPLIT_sp1_1); disp('SPLIT sp1_1 diff');SPLIT_sp1_1_m = max(max(abs(SPLIT_sp1_1))); end
SPLIT_sp1_2 = SPLIT_M.single.idx_cond(2).sp1 - SPLIT_O.single.idx_cond(2).sp1;  if any(SPLIT_sp1_2); disp('SPLIT sp1_2 diff');SPLIT_sp1_2_m = max(max(abs(SPLIT_sp1_2))); end
SPLIT_sp2_1 = SPLIT_M.single.idx_cond(1).sp2 - SPLIT_O.single.idx_cond(1).sp2;  if any(SPLIT_sp2_1); disp('SPLIT sp2_1 diff');SPLIT_sp2_1_m = max(max(abs(SPLIT_sp2_1))); end
SPLIT_sp2_2 = SPLIT_M.single.idx_cond(2).sp2 - SPLIT_O.single.idx_cond(2).sp2;  if any(SPLIT_sp2_2); disp('SPLIT sp2_2 diff');SPLIT_sp2_2_m = max(max(abs(SPLIT_sp2_2))); end
SPLIT_idx_1 = SPLIT_M.group.idx_cond(1).sp - SPLIT_O.group.idx_cond(1).sp;                          if any(SPLIT_idx_1); disp('SPLIT sidx_1.sp diff');SPLIT_idx_1_m = max(max(abs(SPLIT_idx_1))); end
SPLIT_idx_2 = SPLIT_M.group.idx_cond(2).sp - SPLIT_O.group.idx_cond(2).sp;                          if any(SPLIT_idx_2); disp('SPLIT sidx_2.sp diff');SPLIT_idx_2_m = max(max(abs(SPLIT_idx_2))); end
SPLIT_ubidx_1 = SPLIT_M.group.unbalanced_idx_cond(1).sp - SPLIT_O.group.unbalanced_idx_cond(1).sp;  if any(SPLIT_ubidx_1); disp('SPLIT ubidx_1.sp diff');SPLIT_ubidx_1_m = max(max(abs(SPLIT_ubidx_1))); end
SPLIT_ubidx_2 = SPLIT_M.group.unbalanced_idx_cond(2).sp - SPLIT_O.group.unbalanced_idx_cond(2).sp;  if any(SPLIT_ubidx_2); disp('SPLIT ubidx_2.sp diff');SPLIT_ubidx_2_m = max(max(abs(SPLIT_ubidx_2))); end



% Diagnostics file 1
A1 = MCSMO_M.censor_mot - MCSMO_O.censor_mot ;          if any(A1); disp('A1 diff');A1m = max(max(abs(A1))); end
B1 = MCSMO_M.censor_vol - MCSMO_O.censor_vol ;          if any(B1); disp('B1 diff');B1m = max(max(abs(B1))); end
C1 = MCSMO_M.censor_slc - MCSMO_O.censor_slc ;          if any(C1); disp('C1 diff');C1m = max(max(abs(C1))); end
D1 = MCSMO_M.censor_volmot - MCSMO_O.censor_volmot ;    if any(D1); disp('D1 diff');D1m = max(max(abs(D1))); end
E1 = MCSMO_M.censor_slcmot - MCSMO_O.censor_slcmot ;    if any(E1); disp('E1 diff');E1m = max(max(abs(E1))); end
F1 = MCSMO_M.eigfract_fmri - MCSMO_O.eigfract_fmri;     if any(F1); disp('F1 diff');F1m = max(max(abs(F1))); end
G1 = MCSMO_M.eigfract_mot - MCSMO_O.eigfract_mot;       if any(G1); disp('G1 diff');G1m = max(max(abs(G1))); end
H1 = MCSMO_M.eigimages_fmri - MCSMO_O.eigimages_fmri;   if any(H1); disp('H1 diff');H1m = max(max(abs(H1))); end
I1 = MCSMO_M.eigvect_fmri - MCSMO_O.eigvect_fmri;       if any(I1); disp('I1 diff');I1m = max(max(abs(I1))); end
J1 = MCSMO_M.eigweights_mot - MCSMO_O.eigweights_mot;   if any(J1); disp('J1 diff');J1m = max(max(abs(J1))); end
K1 = MCSMO_M.eigvect_mot - MCSMO_O.eigvect_mot;         if any(K1); disp('K1 diff');K1m = max(max(abs(K1))); end
L1 = MCSMO_M.global_signal - MCSMO_O.global_signal;     if any(L1); disp('L1 diff');L1m = max(max(abs(L1))); end

% Diagnostics file 2
A2 = SMO_M.censor_mot - SMO_O.censor_mot ;          if any(A2); disp('A2 diff');A2m = max(max(abs(A2))); end
B2 = SMO_M.censor_vol - SMO_O.censor_vol ;          if any(B2); disp('B2 diff');B2m = max(max(abs(B2))); end
C2 = SMO_M.censor_slc - SMO_O.censor_slc ;          if any(C2); disp('C2 diff');C2m = max(max(abs(C2))); end
D2 = SMO_M.censor_volmot - SMO_O.censor_volmot ;    if any(D2); disp('D2 diff');D2m = max(max(abs(D2))); end
E2 = SMO_M.censor_slcmot - SMO_O.censor_slcmot ;    if any(E2); disp('E2 diff');E2m = max(max(abs(E2))); end
F2 = SMO_M.eigfract_fmri - SMO_O.eigfract_fmri;     if any(F2); disp('F2 diff');F2m = max(max(abs(F2))); end
G2 = SMO_M.eigfract_mot - SMO_O.eigfract_mot;       if any(G2); disp('G2 diff');G2m = max(max(abs(G2))); end
H2 = SMO_M.eigimages_fmri - SMO_O.eigimages_fmri;   if any(H2); disp('H2 diff');H2m = max(max(abs(H2))); end
I2 = SMO_M.eigvect_fmri - SMO_O.eigvect_fmri;       if any(I2); disp('I2 diff');I2m = max(max(abs(I2))); end
J2 = SMO_M.eigweights_mot - SMO_O.eigweights_mot;   if any(J2); disp('J2 diff');J2m = max(max(abs(J2))); end
K2 = SMO_M.eigvect_mot - SMO_O.eigvect_mot;         if any(K2); disp('K2 diff');K2m = max(max(abs(K2))); end
L2 = SMO_M.global_signal - SMO_O.global_signal;     if any(L2); disp('L2 diff');L2m = max(max(abs(L2))); end

%% Optimization Resaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATfiles

% Code property
OSc_1 = isequal(OS_M1,OS_O1);                if ~OSc_1; disp('OSc_1 diff');OS_1_m = OSc_1; end

% Metric opt
%con
OSc_2_R = OS_M2.con.R - OS_O2.con.R ;        if any(OSc_2_R); disp('OSc_2 R diff');OSc_2_R_m = max(max(abs(OSc_2_R))); end
OSc_2_P = OS_M2.con.P - OS_M2.con.P;         if any(OSc_2_P); disp('OSc_2 P diff');OSc_2_P_m = max(max(abs(OSc_2_P))); end
OSc_2_Acc = OS_M2.con.Acc - OS_M2.con.Acc;   if any(OSc_2_Acc); disp('OSc_2 Acc diff');OSc_2_Acc_m = max(max(abs(OSc_2_Acc))); end
OSc_2_dPR = OS_M2.con.dPR - OS_M2.con.dPR ;  if any(OSc_2_dPR); disp('OSc dPr diff');OSc_2_dPR_m = max(max(abs(OSc_2_dPR))); end
%fix
OSf_2_R = OS_M2.fix.R - OS_O2.fix.R ;        if any(OSf_2_R); disp('OSf_2 R diff');OSf_2_R_m = max(max(abs(OSf_2_R))); end
OSf_2_P = OS_M2.fix.P - OS_M2.fix.P;         if any(OSf_2_P); disp('OSf_2 P diff');OSf_2_P_m = max(max(abs(OSf_2_P))); end
OSf_2_Acc = OS_M2.fix.Acc - OS_M2.fix.Acc;   if any(OSf_2_Acc); disp('OSf_2 Acc diff');OSf_2_Acc_m = max(max(abs(OSf_2_Acc))); end
OSf_2_dPR = OS_M2.fix.dPR - OS_M2.fix.dPR ;  if any(OSf_2_dPR); disp('OSf_2 dPR diff');OSf_2_dPR_m = max(max(abs(OSf_2_dPR))); end
%min
OSmin_2_R = OS_M2.min.R - OS_O2.min.R ;        if any(OSmin_2_R); disp('OSmin_2 R diff');OSmin_2_R_m = max(max(abs(OSmin_2_R))); end
OSmin_2_P = OS_M2.min.P - OS_M2.min.P;         if any(OSmin_2_P); disp('OSmin_2 P diff');OSmin_2_P_m = max(max(abs(OSmin_2_P))); end
OSmin_2_Acc = OS_M2.min.Acc - OS_M2.min.Acc;   if any(OSmin_2_Acc); disp('OSmin_2 Acc diff');OSmin_2_Acc_m = max(max(abs(OSmin_2_Acc))); end
OSmin_2_dPR = OS_M2.min.dPR - OS_M2.min.dPR ;  if any(OSmin_2_dPR); disp('OSmin_2 dPR diff');OSmin_2_dPR_m = max(max(abs(OSmin_2_dPR))); end
%max
OSmax_2_R = OS_M2.max.R - OS_O2.max.R ;        if any(OSmax_2_R); disp('OSmax_2 R diff');OSmax_2_R_m = max(max(abs(OSmax_2_R))); end
OSmax_2_P = OS_M2.max.P - OS_M2.max.P;         if any(OSmax_2_P); disp('OSmax_2 P diff');OSmax_2_P_m = max(max(abs(OSmax_2_P))); end
OSmax_2_Acc = OS_M2.max.Acc - OS_M2.max.Acc;   if any(OSmax_2_Acc); disp('OSmax_2 Acc diff');OSmax_2_Acc_m = max(max(abs(OSmax_2_Acc))); end
OSmax_2_dPR = OS_M2.max.dPR - OS_M2.max.dPR ;  if any(OSmax_2_dPR); disp('OSmax_2 dPR diff');OSmax_2_dPR_m = max(max(abs(OSmax_2_dPR))); end
%ind
OSind_2_R = OS_M2.ind.R - OS_O2.ind.R ;        if any(OSind_2_R); disp('OSind_2 R diff');OSind_2_R_m = max(max(abs(OSind_2_R))); end
OSind_2_P = OS_M2.ind.P - OS_M2.ind.P;         if any(OSind_2_P); disp('OSind_2 P diff');OSind_2_P_m = max(max(abs(OSind_2_P))); end
OSind_2_Acc = OS_M2.ind.Acc - OS_M2.ind.Acc;   if any(OSind_2_Acc); disp('OSind_2 Acc diff');OSind_2_Acc_m = max(max(abs(OSind_2_Acc))); end
OSind_2_dPR = OS_M2.ind.dPR - OS_M2.ind.dPR ;  if any(OSind_2_dPR); disp('OSind_2 dPR diff');OSind_2_dPR_m = max(max(abs(OSind_2_dPR))); end


% Pipeline Sets
if ~(size(OS_M3.Signif_Fix) == size(OS_O3.Signif_Fix)); disp('OS_3 signif_fix diff'); OS_3_SigF_m = false ; end 
OS_3_fix = OS_M3.fix - OS_M3.fix ;  if any(OS_3_fix); disp('OS_3_fix diff');OS_3_fix_m = max(max(abs(OS_3_fix))); end
OS_3_min = OS_M3.min - OS_M3.min ;  if any(OS_3_min); disp('OS_3_min diff');OS_3_min_m = max(max(abs(OS_3_min))); end
OS_3_max = OS_M3.max - OS_M3.max ;  if any(OS_3_max); disp('OS_3_max diff');OS_3_max_m = max(max(abs(OS_3_max))); end
OS_3_con = OS_M3.con - OS_M3.con ;  if any(OS_3_con); disp('OS_3_con diff');OS_3_con_m = max(max(abs(OS_3_con))); end
OS_3_ind = OS_M3.ind - OS_M3.ind ;  if any(OS_3_ind); disp('OS_3_ind diff');OS_3_ind_m = max(max(abs(OS_3_find))); end
if ~all(OS_M3.fix_index == OS_O3.fix_index); disp('OS_3 fix index diff'); OS_3_FixIndx_m = false ; end 
if ~all(OS_M3.min_index == OS_O3.min_index); disp('OS_3 min index diff'); OS_3_MinIndx_m = false ; end 
if ~all(OS_M3.max_index == OS_O3.max_index); disp('OS_3 max index diff'); OS_3_MaxIndx_m = false ; end 
if ~all(OS_M3.ind_index == OS_O3.ind_index); disp('OS_3 ind index diff'); OS_3_IndIndx_m = false ; end 
if ~all(OS_M3.con_index == OS_O3.con_index); disp('OS_3 con index diff'); OS_3_ConIndx_m = false ; end 

% Temp Opp
for i = 1: length(OS_M4) %number of files
    %con
    OSc_4ref(i,:) = OS_M4{1,i}.con.CV_ref - OS_O4{1,i}.con.CV_ref;       
    OSc_4alt(i,:) = OS_M4{1,i}.con.CV_alt - OS_O4{1,i}.con.CV_alt;       
    OSc_4altvar(i) = OS_M4{1,i}.con.CV_alt_varfract - OS_M4{1,i}.con.CV_alt_varfract; 
    %fix
    OSf_4ref(i,:) = OS_M4{1,i}.fix.CV_ref - OS_O4{1,i}.fix.CV_ref;      
    OSf_4alt(i,:) = OS_M4{1,i}.fix.CV_alt - OS_O4{1,i}.fix.CV_alt;       
    OSf_4altvar(i) = OS_M4{1,i}.fix.CV_alt_varfract - OS_M4{1,i}.fix.CV_alt_varfract; 
    %ind
    OSi_4ref(i,:) = OS_M4{1,i}.ind.CV_ref - OS_O4{1,i}.ind.CV_ref;     
    OSi_4alt(i,:) = OS_M4{1,i}.ind.CV_alt - OS_O4{1,i}.ind.CV_alt;      
    OSi_4altvar(i) = OS_M4{1,i}.ind.CV_alt_varfract - OS_M4{1,i}.ind.CV_alt_varfract; 
end
%con
if any(OSc_4ref); disp(['OSc_4ref diff' ]); OSc_4ref_m = max(max(abs( OSc_4ref))); end
if any(OSc_4alt); disp(['OSc_4alt diff' ]); OSc_4alt_m = max(max(abs( OSc_4alt))); end
if any(OSc_4altvar); disp(['OSc_4altvar diff']); OSc_4altvar_m = max(max(abs( OSc_4altvar))); end
%fix
if any(OSf_4ref); disp(['OSf_4ref diff' ]); OSf_4ref_m = max(max(abs( OSf_4ref))); end
if any(OSf_4alt); disp(['OSf_4alt diff' ]); OSf_4alt_m = max(max(abs( OSf_4alt))); end
if any(OSf_4altvar); disp(['OSf_4altvar diff' ]); OSf_4altvar_m = max(max(abs( OSf_4altvar))); end
%ind
if any(OSi_4ref); disp(['OSi_4ref diff' ]); OSi_4ref_m = max(max(abs( OSi_4ref))); end
if any(OSi_4alt); disp(['OSi_4alt diff' ]); OSi_4alt_m = max(max(abs( OSi_4alt))); end
if any(OSi_4altvar); disp(['OSi_4altvar diff']); OSi_4altvar_m = max(max(abs( OSi_4altvar))); end


%% Intermediate Metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%More split info from res0
NNw = NNw_M - NNw_O; if any(NNw); disp('NNw  diff');NNw_m= max(max(abs(NNw))); end
WNw = WNw_M - WNw_O; if any(WNw); disp('WNw  diff');WNw_m= max(max(abs(WNw))); end
NNm = NNm_M - NNm_O; if any(NNm); disp('NNm  diff');NNm_m= max(max(abs(NNm))); end
WMm = WMm_M - WMm_O; if any(WMm); disp('WMm  diff');WMm_m= max(max(abs(WMm))); end
FXYZa = FXYZa_M - FXYZa_O; if any(FXYZa); disp('FXYZa  diff');FXYZa_m= max(max(abs(NNw))); end
SW = SW_M - SW_O; if any(SW); disp('SW  diff');SW_m= max(max(abs(SW))); end
MV = MV_M - MV_O; if any(MV); disp('MV  diff');MV_m= max(max(abs(MV))); end
Xn = Xn_M - Xn_O; if any(Xn); disp('Xn  diff');Xn_m= max(max(abs(Xn))); end
Xs = Xs_M - Xs_O; if any(Xs); disp('Xs  diff');Xs_m= max(max(abs(Xs))); end

%res1 spm images 1024(pipeset)
for i = 1:size(ImgSet_M)
    ImgSet(i,:) = ImgSet_M{i,1} - ImgSet_O{i,1};
end
if any(any(ImgSet)); disp(['ImgSet diff']); ImgSet_m = max(max(abs( ImgSet))); end

%res2 temp 1024(pipeset)
for i = 1:size(TCV_M)
    TCVref(i,:) = TCV_M{i,1}.CV_ref - TCV_O{i,1}.CV_ref; 
    TCValt(i,:) = TCV_M{i,1}.CV_alt - TCV_O{i,1}.CV_alt;  
    TCValtvar(i) = TCV_M{i,1}.CV_alt_varfract - TCV_O{i,1}.CV_alt_varfract;
end
if any(any(TCVref)); disp(['TCVref diff']); TCVref_m = max(max(abs( TCVref))); end
if any(any(TCValt)); disp(['TCValt diff']); TCValt_m = max(max(abs( TCValt))); end
if any(TCValtvar); disp(['TCValtvar diff']); TCValtvar_m = max(max(abs( TCValtvar))); end


%res3 stats 1024(pipeset)
for i = 1:size(MsAP_M)
    MsAP_MOT(i)  = MsAP_M{i,1}.artifact_prior.MOT_corr - MsAP_O{i,1}.artifact_prior.MOT_corr;
    MsAP_WM(i)  = MsAP_M{i,1}.artifact_prior.WM_zscor - MsAP_O{i,1}.artifact_prior.WM_zscor;
    MsAP_GS(i)  = MsAP_M{i,1}.artifact_prior.GS_fract - MsAP_O{i,1}.artifact_prior.GS_fract;
    
    Ms_R(i)  = MsAP_M{i,1}.R - MsAP_O{i,1}.R;
    Ms_P(i)  = MsAP_M{i,1}.P - MsAP_O{i,1}.P;                %ERROR MUST START HERE
    Ms_Acc(i)  = MsAP_M{i,1}.Acc - MsAP_O{i,1}.Acc;
    Ms_dPR(i)  = MsAP_M{i,1}.dPR - MsAP_O{i,1}.dPR;
end
if any(any(MsAP_MOT));  disp(['MsAP_MOT diff']); MsAP_MOT_m = max(max(abs( MsAP_MOT))); end
if any(any(MsAP_WM));   disp(['MsAP_WM diff']);  MsAP_WM_m = max(max(abs( MsAP_WM))); end
if any(MsAP_GS);        disp(['MsAP_GS diff']);  MsAP_GS_m = max(max(abs( MsAP_GS))); end
if any(any(Ms_R));      disp(['Ms_R diff']); Ms_R_m = max(max(abs( Ms_R))); end
if any(any(Ms_P));      disp(['Ms_P diff']);  Ms_P_m = max(max(abs( Ms_P))); end                %THIS FINDS THE ERROR
if any(Ms_Acc);         disp(['Ms_Acc diff']);  Ms_Acc_m = max(max(abs( Ms_Acc))); end
if any(Ms_dPR);         disp(['Ms_dPR diff']);  Ms_dPR_m = max(max(abs( Ms_dPR))); end


%% QC1 Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QC1_ms_ad = QC1_M.motion_stats.avg_disp - QC1_O.motion_stats.avg_disp;          if any(QC1_ms_ad); disp('QC1_ms_ad diff');QC1_ms_ad_m = max(max(abs(QC1_ms_ad))); end
QC1_ms_ns = QC1_M.motion_stats.num_spikes - QC1_O.motion_stats.num_spikes;      if any(QC1_ms_ns); disp('QC1_ms_ns diff');QC1_ms_ns_m = max(max(abs(QC1_ms_ns))); end
QC1_ms_st = QC1_M.motion_stats.task_corr - QC1_O.motion_stats.task_corr;        if any(QC1_ms_st); disp('QC1_ms_st diff');QC1_ms_st_m = max(max(abs(QC1_ms_st))); end
QC1_sa_mc = QC1_M.spm_artifact.motion_corr - QC1_O.spm_artifact.motion_corr;    if any(QC1_sa_mc); disp('QC1_sa_mc diff');QC1_sa_mc_m = max(max(abs(QC1_sa_mc))); end
QC1_sa_pf = QC1_M.spm_artifact.pos_fraction - QC1_O.spm_artifact.pos_fraction;  if any(QC1_sa_pf); disp('QC1_sa_pf diff');QC1_sa_pf_m = max(max(abs(QC1_sa_pf))); end
QC1_sa_wz = QC1_M.spm_artifact.wm_zscore - QC1_O.spm_artifact.wm_zscore;        if any(QC1_sa_wz); disp('QC1_sa_wz diff');QC1_sa_wz_m = max(max(abs(QC1_sa_wz))); end
QC1_m_asc = QC1_M.metrics.avg_spm_corr - QC1_O.metrics.avg_spm_corr;            if any(QC1_m_asc); disp('QC1_m_asc diff');QC1_m_asc_m = max(max(abs(QC1_m_asc))); end

for i = 1:size(QC1_M.metrics.perform_metrics)
    QC1_m_pm = QC1_M.metrics.perform_metrics{i,1} - QC1_O.metrics.perform_metrics{i,1};
end
if any(QC1_m_pm); disp(['QC1_m_pm diff']);  QC1_m_pm_m = max(max(abs(QC1_m_pm))); end


%% QC2 Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the PDFs, why octave is not exporting

%Group PCA file
% QC2_con_SPM = QCG2_M.con.rSPMZs - QCG2_O.con.rSPMZs;    if any(QC2_con_SPM); disp('QC2_con_SPM diff');QC2_con_SPM_m = max(max(abs(QC2_con_SPM))); end
% QC2_con_rep = QCG2_M.con.rep - QCG2_O.con.rep;          if any(QC2_con_rep); disp('QC2_con_rep diff');QC2_con_rep_m = max(max(abs(QC2_con_rep))); end
% QC2_con_sc = QCG2_M.con.scores -QCG2_O.con.scores;      if any(QC2_con_sc); disp('QC2_con_sc diff');QC2_con_sc_m = max(max(abs(QC2_con_sc))); end
% QC2_con_var = QCG2_M.con.var - QCG2_O.con.var;          if any(QC2_con_var); disp('QC2_con_var diff');QC2_con_var_m = max(max(abs(QC2_con_var))); end
% QC2_fix_SPM  = QCG2_M.fix.rSPMZs -QCG2_O.fix.rSPMZs;    if any(QC2_fix_SPM); disp('QC2_fix_SPM diff');QC2_fix_SPM_m = max(max(abs(QC2_fix_SPM))); end
% QC2_fix_rep = QCG2_M.fix.rep -QCG2_O.fix.rep;           if any(QC2_fix_rep); disp('QC2_fix_rep diff');QC2_fix_rep_m = max(max(abs(QQC2_fix_rep))); end
% QC2_fix_sc = QCG2_M.fix.scores - QCG2_O.fix.scores;     if any(QC2_fix_sc); disp('QC2_fix_sc diff');QC2_fix_sc_m = max(max(abs(QC2_fix_sc))); end
% QC2_fix_var = QCG2_M.fix.var - QCG2_O.fix.var;          if any(QC2_fix_var); disp('QC2_fix_var diff');QC2_fix_var_m = max(max(abs(QC2_fix_var))); end
% QC2_ind_SPM  = QCG2_M.ind.rSPMZs - QCG2_O.ind.rSPMZs;   if any(QC2_ind_SPM); disp('QC2_ind_SPM diff');QC2_ind_SPM_m = max(max(abs(QC2_ind_SPM))); end
% QC2_ind_rep = QCG2_M.ind.rep - QCG2_O.ind.rep;          if any(QC2_ind_rep); disp('QC2_ind_rep diff');QC2_ind_rep_m = max(max(abs(QC2_ind_rep))); end
% QC2_ind_sc = QCG2_M.ind.scores - QCG2_O.ind.scores;     if any(QC2_ind_sc); disp('QC2_ind_sc diff');QC2_ind_sc_m = max(max(abs(QC2_ind_sc))); end
% QC2_ind_var = QCG2_M.ind.var -QCG2_O.ind.var;           if any(QC2_ind_var); disp('QC2_ind_var diff');QC2_ind_var_m = max(max(abs(QC2_ind_var))); end
% 
% %QC2 group file
% %1
% QC2_spmo_con = QC2_M.spm_ovl.con - QC2_O.spm_ovl.con;    if any(QC2_spmo_con); disp('QC2_spmo_con diff');QC2_spmo_con_m = max(max(abs(QC2_spmo_con))); end
% QC2_spmo_fix = QC2_M.spm_ovl.fix - QC2_O.spm_ovl.fix;    if any(QC2_spmo_fix); disp('QC2_spmo_fix diff');QC2_spmo_fix_m = max(max(abs(QC2_spmo_fix))); end
% QC2_spmo_ind = QC2_M.spm_ovl.ind - QC2_O.spm_ovl.ind;    if any(QC2_spmo_ind); disp('QC2_spmo_ind diff');QC2_spmo_ind_m = max(max(abs(QC2_spmo_ind))); end
% %2
% QC2_spmc_con = QC2_M.spm_corr.con - QC2_O.spm_corr.con;    if any(QC2_spmc_con); disp('QC2_spmc_con diff');QC2_spmc_con_m = max(max(abs(QC2_spmoc_con))); end
% QC2_spmc_fix = QC2_M.spm_corr.fix - QC2_O.spm_corr.fix;    if any(QC2_spmc_fix); disp('QC2_spmc_fix diff');QC2_spmc_fix_m = max(max(abs(QC2_spmoc_fix))); end
% QC2_spmc_ind = QC2_M.spm_corr.ind - QC2_O.spm_corr.ind;    if any(QC2_spmc_ind); disp('QC2_spmc_ind diff');QC2_spmc_ind_m = max(max(abs(QC2_spmoc_ind))); end
% %3
% QC2_ps_SF = QC2_M.pipeline_sets.Signif_Fix - QC2_M.pipeline_sets.Signif_Fix;    if any(QC2_spmc_con); disp('QC2_ps_SF diff');QC2_ps_SF_m = max(max(abs(QC2_ps_SF))); end
% QC2_ps_fix = QC2_M.pipeline_sets.fix - QC2_M.pipeline_sets.fix;    if any(QC2_ps_fix); disp('QC2_ps_fix diff');QC2_ps_fix_m = max(max(abs(QC2_ps_fix))); end
% QC2_ps_min = QC2_M.pipeline_sets.min - QC2_M.pipeline_sets.min;    if any(QC2_ps_min); disp('QC2_ps_min diff');QC2_ps_min_m = max(max(abs(QC2_ps_min))); end
% QC2_ps_max = QC2_M.pipeline_sets.max - QC2_M.pipeline_sets.max;    if any(QC2_ps_max); disp('QC2_ps_max diff');QC2_ps_max_m = max(max(abs(QC2_ps_max))); end
% QC2_ps_con = QC2_M.pipeline_sets.con - QC2_M.pipeline_sets.con;    if any(QC2_ps_con); disp('QC2_ps_con diff');QC2_ps_con_m = max(max(abs(QC2_ps_con))); end
% QC2_ps_ind = QC2_M.pipeline_sets.ind - QC2_M.pipeline_sets.ind;    if any(QC2_ps_max); disp('QC2_ps_ind diff');QC2_ps_ind_m = max(max(abs(QC2_ps_ind))); end
% QC2_ps_ind_index = QC2_M.pipeline_sets.ind_index - QC2_M.pipeline_sets.ind_index;    if any(QC2_ps_ind_index); disp('QC2_ps_ind_index diff');QC2_ps_ind_index_m = max(max(abs(QC2_ps_ind_index))); end
% 
% %4 and 5
% for i = 1:size(QC2_M.SPM_opt_norm)
%     %4
%     QC2_spmopt_con(i,:) = QC2_M.SPM_opt_norm{1,i}.con - QC2_O.SPM_opt_norm{1,i}.con;
%     QC2_spmopt_fix(i,:) = QC2_M.SPM_opt_norm{1,i}.fix - QC2_O.SPM_opt_norm{1,i}.fix;
%     QC2_spmopt_ind(i,:) = QC2_M.SPM_opt_norm{1,i}.ind - QC2_O.SPM_opt_norm{1,i}.ind;
%     %5
%     QC2_tempopt_con(i,:) = QC2_M.TEMP_opt{1,i}.con - QC2_O.TEMP_opt{1,i}.con;
%     QC2_tempopt_fix(i,:) = QC2_M.TEMP_opt{1,i}.fix - QC2_O.TEMP_opt{1,i}.fix;
%     QC2_temppt_ind(i,:) = QC2_M.TEMP_opt{1,i}.ind - QC2_O.TEMP_opt{1,i}.ind;
% end
% %4
% if any(QC2_spmopt_con); disp('QC2_spmopt_con diff');  QC2_spmopt_con_m = max(max(abs(QC2_spmopt_con))); end
% if any(QC2_spmopt_fix); disp('QC2_spmopt_fix diff');  QC2_spmopt_fix_m = max(max(abs(QC2_spmopt_fix))); end
% if any(QC2_spmopt_ind); disp('QC2_spmopt_ind diff');  QC2_spmopt_ind_m = max(max(abs(QC2_spmopt_ind))); end
% %5
% if any(QC2_tempopt_con); disp('QC2_tempopt_con diff');  QC2_tempopt_con_m = max(max(abs(QC2_tempopt_con))); end
% if any(QC2_tempopt_fix); disp('QC2_tempopt_fix diff');  QC2_tempopt_fix_m = max(max(abs(QC2_tempopt_fix))); end
% if any(QC2_tempopt_ind); disp('QC2_tempopt_ind diff');  QC2_tempopt_ind_m = max(max(abs(QC2_tempopt_ind))); end











