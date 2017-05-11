function posthoc_analysis_fmri_func( inputfile, newmaskname, ppl_opt, contr_num, analysis1, analysis2, suppDat, design, pathstuff )
%
%
% Syntax:
%            posthoc_analysis( inputfile, newmaskname, ppl_opt, contr_num, analysis, design, threshold )
%
% Inputs:
%
%              inputfile: string specifying path+name of file used to run full pipeline optimization
%            newmaskname: string specifying binary brain tissue mask. Leave empty (newmaskname=[])
%                         if you want to use default group_consensus_mask.nii produced by OPPNI 
%                ppl_opt: string specifying which set of optimized pipelines to run analysis on
%                         options: 'CON', 'FIX' or 'IND'
%              contr_num: number specifying which contrast, in cases of multi-run optimization 
%                         e.g. if you run 'taskA-baseline,taskB-baseline,taskA-taskB', contr_num=2 or contr_num='2' will analyze taskB-baseline contrast
%                              if you want to analyze differences between contrasts, do list, e.g. [1 2] or '1,2'  
%               analysis: defines your analysis method. Options:
%                           'Ttest'    : standard t-test (1-group,2-group unpaired, 2-group paired)
%                           'Bootstrap': non-parametric bootstrapped analysis of differences (1-way,2-way unpaired, 2-way paired)
%                           'Splithalf': non-parametric split-half Z-scored analysis of differences (1-way,2-way unpaired, 2-way paired)
%                           'ANOVA1'   : 1-way analysis of variance (1 random factor; main effect) 
%                           'ANOVA2'   : 2-way analysis of variance (2 random factors; main effects + interaction) 
%                           'ANOVA1rep': 1-way repeated measures (1 random factor, 1 repeated measure; main effects + interactions) 
%                  design: design matrix specifying factors/covariates of interest. the nth row show correspond to the nth
%                         line in the 'inputfile'  
%
% Outputs: Generates a "Posthoc analyses" folder in your output directory. 
%          Stores brain SPMs produced by your analysis model of choice. 
%
% -------------------------------------------------------------------------------------------
% EXAMPLES: running post-hoc testing
%
%  .consider a dataset where you analyzed 20 subjects, with three task contrasts
%   ('taskA-baseline,taskB-baseline,taskA-taskB'), with inputfile 'inputs.txt'
%
%
%  1. Testing for consistent group-level activation (1-sample testing).
%     The following command is submitted:
%     
%         posthoc_analysis( 'inputs.txt', [], 'IND', 1, 'Ttest', [] );
%
%     this will analyze individually optimized pipelines (ppl_opt='IND') for the 
%     'taskA-baseline' contrast (contr_num=1). You are running a t-test based on
%     the fact that (a) contr_num gives a single contrast and (b) design=[]
%     means there are no other factors of interest. This will default to
%     1-way analysis. 
%     NB1: substitute 'Bootstrap' for 'Ttest' to run the equivalent non-parametric estimates
%
%  2. Testing for consistent between-group differences, 2 groups (2-sample unpaired).
%     Say that subjects 1-10 are controls, and subjects 11-20 are patients
%     Define the following design matrix:
%                
%         design = [zeros(1,10), ones(1,10)]';
%
%     The following command is submitted:
%     
%         posthoc_analysis( 'inputs.txt', [], 'IND', 1, 'Ttest', design );
%
%     Same as (ex.1), but now you have specified a design factor of interest via design.
%     NB1: you can again swap 'Bootstrap' for 'Ttest' to run the equivalent non-parametric estimates
%     NB2: t-test and bootstrap only work for 2-group comparisons. They will fail if design codes >2 groups
%          in this case, you have to run a 1-way ANOVA. For example, consider that 
%          patients 11-15 are distinct from patients 16-20. Then you would
%          code the design matrix:
%           
%              design = [zeros(1,10), 1*ones(1,5), 2*ones(1,5)]';
%
%          and submit the command
%
%              posthoc_analysis( 'inputs.txt', [], 'IND', 1, 'ANOVA1', design, [0 20] );
%
%  3. Testing for consistent between-group differences, 2 matched groups (2-sample paired).
%     From (ex.2), say that patients and controls are now individually matched, e.g. 
%     subj-1 (control) matches subj-11 (patient), subj-2 (control) matches subj-12 (patient) ... 
%     Now we define the following design matrix:
%                
%         design = [zeros(1,10), ones(1,10); 1:10, 1:10]';
%
%     The following command is submitted:
%     
%         posthoc_analysis( 'inputs.txt', [], 'IND', 1, 'Ttest', design );
%
%     Same as (ex.2), but now you have specified an additional subject-specific factor for column 2 of design.
%     NB1: you can again swap 'Bootstrap' for 'Ttest' to run the equivalent non-parametric estimates
%     NB2: as in (ex.2), if you have >2 groups, you will need to run a 1-way (*repeated-measures*) ANOVA. 
%          in this case, you would code the design matrix:
%           
%              design = [zeros(1,10), 1*ones(1,5), 2*ones(1,5); 1:10, 1:10 ]';
%
%          and submit the command
%
%              posthoc_analysis( 'inputs.txt', [], 'IND', 1, 'ANOVA1rep', design, [0 20] );
%
%  4. Testing for consistent between-contrast differences (2-sample paired)
%     If you want to directly compare contrasts produced by multi-run analysis, there is
%     special syntax. As in (ex.1) we treat all 20 subjects as a single group
%     The following command is submitted:
%     
%         posthoc_analysis( 'inputs.txt', [], 'IND', '1,2', 'Ttest', [] );
%
%     This will run 2-group paired analysis contrasting SPMs from 'taskA-baseline' vs. 'taskB-baseline' 
%     NB1: again, 'Ttest' can be exchanged for 'Bootstrap' for non-parametric analyses 
%     NB2: if you want to compare multiple contrasts, you could code e.g. contr_num='1,2,3'
%          and run ANOVA1rep.
%
%


% get location for result files
[InputStruct, MULTI_RUN_INPUTFILE] = Read_Input_File(inputfile);
% output path for results (same as group mask)
PH_folder = [InputStruct(1).run(1).Output_nifti_file_path,'/Posthoc_analyses'];
mkdir (PH_folder);
mkdir ([PH_folder,'/tempfiles']);


% name for group masks --> default should be empty
if( isempty(newmaskname) )
    %%% create directory + new mask files, in first output dir.
    newpath = [InputStruct(1).run(1).Output_nifti_file_path, '/GroupMasks'];
    mkdir_r(newpath);
    newmaskname = [newpath '/group'];                
else
    [apath,aprefix,aext] = fileparts(newmaskname);
    mkdir_r(apath);
    %%% make sure mask-name terminates with .nii string
    newmaskname = [apath,'/',aprefix,'.nii'];
end

% % pipeline optimization method - convert to numeric index
% if    ( ischar(ppl_opt) && strcmpi(ppl_opt,'PIPE1') ) ppl_opt=1;
% elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'CON') )   ppl_opt=1;
% elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'FIX') )   ppl_opt=2;
% elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'IND') )   ppl_opt=3;
% end

% list of contrasts to analyze
%
if(isnumeric(contr_num)) %% if a number, convert to string
    costr = num2str( contr_num );
else %% if a string, format appropriately
    
    costr = contr_num;
    
    % check if comma-separated, indicating multiple contrasts to analyze
    if( ~isempty(strfind(contr_num,',')) )
        cpair = regexp(contr_num,',','split');
        contr_num = str2num(cell2mat(cpair'));
    else
        contr_num = str2num(contr_num);
    end 
end
% end result is numeric vector listing contrasts of interest

MM = load_untouch_nii([newmaskname]);
mask = double(MM.img);

if(ischar(design)   ) design = dlmread(design); end
if( numel(design)==1) design = [];              end

%% upload relevant [pipeline files]

% check if files already stored:
filename = [PH_folder,'/tempfiles/tempvol',num2str(numel(InputStruct)),'.mat'];
if(~exist(filename,'file'))
    disp('loading + converting functional runs...');
    for(is=1:numel(InputStruct))

        is,

        % load split_info data for further analysis
        Subject_OutputDirIntermed = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_metrics'];
        subjectprefix = InputStruct(is).run(1).subjectprefix;
        load([Subject_OutputDirIntermed,'/',pathstuff{1},'/res0_params/params' subjectprefix '.mat']);

        outdir = InputStruct(is).run(1).Output_nifti_file_path;
        prefix = InputStruct(is).run(1).Output_nifti_file_prefix;

        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------
        VV       = load_untouch_nii(strcat(outdir,'/optimization_results/',pathstuff{1},'/images_',pathstuff{2},'/Proc_',prefix,'_',ppl_opt,'_sNorm.nii'));  
        volmat   = nifti_to_mat( VV,MM ); %% load volumes
        filename = [PH_folder,'/tempfiles/tempvol',num2str(is),'.mat'];
        save(filename,'volmat');
        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------
    end
end

analysis1 = regexp(analysis1,'-','split');
Ksel      = str2num(analysis1{2});


disp('now reloading + preparatory analyses...');

if(strcmpi(analysis1{1},'PCA') || strcmpi(analysis1{1},'CLUST'))
    
    mkdir([ PH_folder,'/tempfiles/PCA']);
    
    %%% CHECK: group-level PCA, component estimation
    if(~exist([ PH_folder,'/tempfiles/PCA/U_group.mat'],'file'))
        
        % load functional volumes
        for(is=1:numel(InputStruct))
            filename = [PH_folder,'/tempfiles/tempvol',num2str(is),'.mat'];
            load(filename);
            [u l v]  = svd( volmat,'econ' );
            uset{is} = u(:,1:40);
            ntime(is)= size(volmat,2);
        end
        kmin = floor(0.1*min(ntime));
        
        % rv coefficient testing for pca-space outliers
        disp('checking for PCA outliers...');
        rvmat = ones( numel(InputStruct));
        for(i=1:numel(InputStruct)-1)
            for(j=i+1:numel(InputStruct))
                [i j],
                rvmat(i,j) = RV_coef( uset{i}(:,1:kmin), uset{j}(:,1:kmin), [] );
                rvmat(j,i) = rvmat(i,j);
            end
        end
        g    = zscore(  (sum(rvmat)-1)./(size(rvmat,1)-1)  );
        disp(['outliers dropped: ',num2str( sum( g < -1.65 ) )])
        uset = uset( g > -1.65 );
        
        % rv coefficient dimensionality estimation
        disp('dimensionality selection...');
        k=0;
        for(i=1:numel(uset)-1)
            for(j=i+1:numel(uset))
                k=k+1;
                [i j],
                for(t=1:40)
                   rvset(t,k) = RV_coef( uset{i}(:,1:t), uset{j}(:,1:t), [] ); 
                end
            end
        end
        figure,plot( mean(rvset,2),'o-k'); ylabel('RV coefficient'); xlabel('dimensionality');
        [vx kmax] = max( mean(rvset,2) );
        if(kmax<kmin) 
            disp('dimensionality very low! Resetting to minimum threshold of 10% input dimensionality');
            kmax=kmin; 
        end
        
        % estimation of concat-PCA components
        ucat   = [];
        for(is=1:numel(uset))
            ucat = [ucat uset{is}(:,1:kmax)];
        end
        [Ugrp Lgrp Vgrp] = svd( ucat,'econ' );
        Ugrp = Ugrp(:,1:kmax);
        
        % estimation of timeseries weights, etc
        for(is=1:numel(InputStruct))
            filename = [PH_folder,'/tempfiles/tempvol',num2str(is),'.mat'];
            load(filename);
            tser{is} = zscore(volmat') * Ugrp;
            tvar(is,:) = sum(tser{is}.^2);
            regmaps{is} = zscore(volmat')'*zscore(tser{is}) ./ (size(tser{is},1)-1);
        end
        
        save([PH_folder,'/tempfiles/PCA/U_group.mat'],'Ugrp','kmax','tvar','rvset');
        save([PH_folder,'/tempfiles/PCA/u_subject.mat'],'kmax','tser','regmaps');
        
        [ clust_out ] = kmeans_quick( Ugrp, kmax, 100 );
        for(i=1:kmax) BINMAT(:,i) = double(clust_out.LAB==i); end
        lix = sortrows([(1:kmax)' sum(BINMAT)'],-2);
        BINMAT = BINMAT(:,lix(:,1));
        LABMAT = BINMAT * (1:kmax)';
        
        % estimation of timeseries weights, etc
        for(is=1:numel(InputStruct))
            filename = [PH_folder,'/tempfiles/tempvol',num2str(is),'.mat'];
            load(filename);
            tserK{is} = zscore(volmat') * bsxfun(@rdivide,BINMAT,sum(BINMAT));
            regmapsK{is} = zscore(volmat')'*zscore(tserK{is}) ./ (size(tserK{is},1)-1);
        end
        
        save([PH_folder,'/tempfiles/PCA/K_group.mat'],'BINMAT','LABMAT');       
        save([PH_folder,'/tempfiles/PCA/k_subject.mat'],'tserK','regmapsK');       
    end
    
    
    if( str2num(analysis1{2})==0 )

        load([PH_folder,'/tempfiles/PCA/U_group.mat']);
        load([PH_folder,'/tempfiles/PCA/K_group.mat']);       
        
        if( strcmpi(analysis1{1},'PCA') )
        components_disp( Ugrp, mask )
        TMPVOL = zeros([size(mask), kmax]);
        for(k=1:kmax) tmp=mask;tmp(tmp>0)=Ugrp(:,k)./std(Ugrp(:,k)); TMPVOL(:,:,:,k)=tmp; end
        % predefining nifti fields
        nii                   = MM;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist          = MM.hdr.hist;
        nii.hdr.dime.dim(5)   = kmax;
        nii.img               = TMPVOL;
        mkdir([PH_folder,'/func']);
        save_untouch_nii(nii,strcat(PH_folder,'/func/PCA.nii'));
        elseif( strcmpi(analysis1{1},'CLUST') )
        axial_plot( LABMAT,mask,6,[0 kmax],1); colorbar;
        % predefining nifti fields
        nii                   = MM;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist          = MM.hdr.hist;
        nii.hdr.dime.dim(5)   = 1;
        tmp = mask; tmp(tmp>0) = LABMAT;
        nii.img               = tmp;
        mkdir([PH_folder,'/func']);
        save_untouch_nii(nii,strcat(PH_folder,'/func/CLUST.nii'));
        end
    
        spm_set=[];
    else
        load([PH_folder,'/tempfiles/PCA/u_subject.mat']);
        load([PH_folder,'/tempfiles/PCA/k_subject.mat']);
        
        Ksel = str2num(analysis1{2});
        
        for(is=1:numel(InputStruct))
            if( strcmpi(analysis1{1},'PCA') )
            spm_set(:,is) = regmaps{is}(:,Ksel);
            elseif( strcmpi(analysis1{1},'CLUST') )
            spm_set(:,is) = regmapsK{is}(:,Ksel);
            end
        end
    end
end

disp('now');

%%

if(~isempty(spm_set))
%% now: run analysis

analysis_model.filepath = which(analysis2);
% run analaysis --> ensure code is in path, call it from string
addpath( analysis_model.filepath );
pfun = str2func( analysis2 );

% run the analysis
output = pfun( spm_set, design );  

% pull list of all outputs from analysis
e = fieldnames( output );
% drop the "testname" field from the list
ix = cellfun(@(s) ~isempty(strfind('testname', s)), e);
e(ix) = [];
% split into test statistics and associated p-value stats
ix = cellfun(@(s) ~isempty(strfind('_p', s(end-1:end) )), e);
e_stat  = e(ix==0); 
e_pval0 = e(ix >0); 
%%% check enough paired stats/pvalues
if( numel(e_stat) ~= numel(e_pval0) )
    disp('warning: pstats not available for all test statistics!');
end

% predefining nifti fields
nii                   = MM;
nii.hdr.dime.datatype = 16;
nii.hdr.hist          = MM.hdr.hist;
nii.hdr.dime.dim(5)   = 4;
nii.img               = zeros( [size(mask) 4] );

% predefined cluster size (based on FWHM=6 smoothing, standard task data, 2x2x2 mm res.)
MinClustSize = 73;

for(i=1:length(e_stat))
       
    map = output.(e_stat{i});
    pmap= output.([e_stat{i},'_p']);
    
    % (1) unthresholded
    tmp=mask; tmp(tmp>0)=map; 
    nii.img(:,:,:,1) = tmp;
    % (2) thresholded p<0.005 uncorr
    tmp=mask; tmp(tmp>0)=map .* double(pmap<=0.005); 
    nii.img(:,:,:,2) = tmp;
    if( sum(abs(tmp(:)))>0)
    % (3) cluster-size thresholded
    tmp=clust_up(tmp,MinClustSize);
    nii.img(:,:,:,3) = tmp;
    end
    % (4) fdr threshold
    [pcrit thr0] = fdr( pmap, 'p', 0.05, 0 );
    tmp=mask; tmp(tmp>0)=map .* thr0; 
    nii.img(:,:,:,4) = tmp;
    % save the nifti file
    mkdir([PH_folder,'/func']);
    save_untouch_nii(nii,strcat(PH_folder,'/func/',output.testname,'_',e{i},'_',analysis1{1},'.',analysis1{2},'.nii'));
end

end
%%
function components_disp( datamat, mask )

N   = size(datamat,2);
Nsq = 2;

figure; 
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1 0.1 0.50 0.80]);

for(n=1:min([9 N]))
    
    data=datamat(:,n)./std(datamat(:,n)); % variance normalized

    % mask dimensions
    [Nx Ny Nz] = size(mask);
    % reshape volume into 3D
    tmp=mask; tmp(tmp>0)= double( data );
    % get index on z-axis
    indlist = floor(Nz/7).*[2:5];

    fullpanel = [tmp(:,:,indlist(1)), tmp(:,:,indlist(2)); tmp(:,:,indlist(3)), tmp(:,:,indlist(4))]';


    % initialize plot
    subplot(3,3,n);
    imagesc( fullpanel, [-2.8 2.8] );
    set(gca,'xtick',[],'ytick',[]);
    title(['PC#',num2str(n)]);

end

function rv = RV_coef( X, Y, colt )
%
% rv = RV_coef( X, Y )
%

X = X - repmat( mean(X), [size(X,1) 1] );
Y = Y - repmat( mean(Y), [size(Y,1) 1] );

SIG_xx = X'*X;
SIG_yy = Y'*Y;
SIG_xy = X'*Y;
SIG_yx = Y'*X;

covv = trace( SIG_xy * SIG_yx );
vavx = trace( SIG_xx * SIG_xx );
vavy = trace( SIG_yy * SIG_yy );

rv = covv ./ sqrt( vavx * vavy );
