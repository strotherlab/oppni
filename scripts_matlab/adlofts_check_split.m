%CHECKS ALL THE INPUTS FOR THE "ANALYSIS WRAPPER FUNTION"

%% THIS WAS DONE WITH MATLABS INPUTS AT WRAPPER LEVEL, the result is the
% same. --> this would indicate that the enterence values are already
% different when entering the wrapper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allPM = load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/allMATLAB.mat');
% allPO = load('/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t/allOCTAVE.mat');
% D = allPM.allP - allPO.allP;
% max(D);

%|
%|
%V
%This running  on 352/all files should show that the inputs are different.

%% CHOOSE THE FILE TO USE - LOOP THROUGH ALL?
for i = 1:371 %371
    
    if (mod(i,2)) %odd
        MP = load(['FindPat_' num2str(i) '.mat']);
        OP = load(['FindPat_' num2str(i) '_Octave.mat']);
    else
        MP = load(['FindPat_' num2str(i) 'a.mat']);
        OP = load(['FindPat_' num2str(i) 'a_Octave.mat']);    
    end
    
    %^^ above are the same inputs
    
    % Optional Test 
    output_temp_M = run_analyses_wrapper( MP.vol_filt, MP.split_info, MP.analysis_model );
    
    %%CANT _ THIS MUST BE DONE IN OCTAVE!!!
    output_temp_O = run_analyses_wrapper( OP.vol_filt, OP.split_info, OP.analysis_model );
    
    % P test
    Pdiff = output_temp_M.metrics.P - output_temp_O.metrics.P
    if(max(max(Pdiff)) >= 0.0001) ; disp('found PDIFF!');end;
    

%%Not Vol_Filt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pdiff = MP.vol_filt - OP.vol_filt;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

%%Struct Elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pdiff = MP.split_info.NN_weight_avg - OP.split_info.NN_weight_avg;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.WM_weight_avg - OP.split_info.WM_weight_avg;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.NN_mask_avg - OP.split_info.NN_mask_avg;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.WM_mask_avg - OP.split_info.WM_mask_avg;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.FXYZ_avg - OP.split_info.FXYZ_avg;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.spat_weight - OP.split_info.spat_weight;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.idx_cond1_sp1 - OP.split_info.idx_cond1_sp1;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.idx_cond1_sp2 - OP.split_info.idx_cond1_sp2;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.idx_cond2_sp1 - OP.split_info.idx_cond2_sp1;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.idx_cond2_sp2 - OP.split_info.idx_cond2_sp2;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.mask_vol - OP.split_info.mask_vol;
if(max(max(max(Pdiff))) >= 0.0001) ; disp('found something');end;

Pdiff = MP.split_info.baseline - OP.split_info.baseline;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

%group
Pdiff = MP.split_info.group.idx_cond(1).sp - OP.split_info.group.idx_cond(1).sp;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.group.idx_cond(2).sp - OP.split_info.group.idx_cond(2).sp;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.group.unbalanced_idx_cond(1).sp - OP.split_info.group.unbalanced_idx_cond(1).sp;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.group.unbalanced_idx_cond(2).sp - OP.split_info.group.unbalanced_idx_cond(2).sp;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

%single
Pdiff = MP.split_info.single.idx_cond(1).sp1 - OP.split_info.single.idx_cond(1).sp1;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.single.idx_cond(2).sp1 - OP.split_info.single.idx_cond(2).sp1;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.single.idx_cond(1).sp2 - OP.split_info.single.idx_cond(1).sp2;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.single.idx_cond(2).sp2 - OP.split_info.single.idx_cond(2).sp2;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;

%cond
Pdiff = MP.split_info.cond(1).onsetlist - OP.split_info.cond(1).onsetlist;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.cond(1).blklength - OP.split_info.cond(1).blklength;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.cond(2).onsetlist - OP.split_info.cond(2).onsetlist;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;
Pdiff = MP.split_info.cond(2).blklength - OP.split_info.cond(2).blklength;
if(max(max(Pdiff)) >= 0.0001) ; disp('found something');end;


end