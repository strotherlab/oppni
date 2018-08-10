   addpath(genpath('/global/home/hpc4253/OG_oppni/oppni-0.7.3.1_06JUL2017'));
%Check the fuinction and the Acc,P,R outputs
for i = 1:1024

    filename = ['_session2_ID6230_run3Outtemp_indx_' num2str(i) '.mat'];
    load(filename);
    
    output_temp_M = run_analyses_wrapper(vol_filt,split_info,'LDA');
    
    savename = ['Indx_' num2str(i) '_matlab.mat'];
    save(savename,'output_temp_M');

end


%Testing Temp
for i = 1000:1024
    disp(i)
    filename = ['Indx_' num2str(i) '_octave.mat'];
    load(filename);
    filename = ['Indx_' num2str(i) '_matlab.mat'];
    load(filename);
    
A = output_temp_M.temp.pp1_2on1 - output_temp_O.temp.pp1_2on1;
B = output_temp_M.temp.pp2_2on1 - output_temp_O.temp.pp2_2on1;
C = output_temp_M.temp.sc1_2on1 - output_temp_O.temp.sc1_2on1;
D = output_temp_M.temp.sc1_2on1 - output_temp_O.temp.sc1_2on1;

E = output_temp_M.temp.pp1_1on2 - output_temp_O.temp.pp1_1on2;
F = output_temp_M.temp.pp2_1on2 - output_temp_O.temp.pp2_1on2;
G = output_temp_M.temp.sc1_1on2 - output_temp_O.temp.sc1_1on2;
H = output_temp_M.temp.sc1_1on2 - output_temp_O.temp.sc1_1on2;

R = output_temp_M.metrics.R - output_temp_O.metrics.R;
P = output_temp_M.metrics.P - output_temp_O.metrics.P;
Acc = output_temp_M.metrics.Acc - output_temp_O.metrics.Acc;

diffA(i) = max(max(A));
diffB(i) = max(max(B));
diffC(i) = max(max(C));
diffD(i) = max(max(D));
diffE(i) = max(max(E));
diffF(i) = max(max(F));
diffG(i) = max(max(G));
diffH(i) = max(max(H));
index(i) = i;



diffR(i) = max(max(R));

diffP(i) = max(max(P));

diffAcc(i) = max(max(Acc));

end








