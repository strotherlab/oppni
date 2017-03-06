function posthoc_analysis_fmri_spms( inputfile, newmaskname, ppl_opt, contr_num, analysis, design )
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
    newmaskname = [apath,'/',aprefix];
end

% pipeline optimization method - convert to numeric index
if    ( ischar(ppl_opt) && strcmpi(ppl_opt,'PIPE1') ) ppl_opt=1;
elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'CON') )   ppl_opt=1;
elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'FIX') )   ppl_opt=2;
elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'IND') )   ppl_opt=3;
end

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

MM = load_untouch_nii([newmaskname '_consensus_mask.nii']);
mask = double(MM.img);

if(ischar(design)   ) design = dlmread(design); end
if( numel(design)==1) design = [];              end

%% upload relevant SPMs
for(is=1:numel(InputStruct))

    is,
    
    % load split_info data for further analysis
    Subject_OutputDirIntermed = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_metrics'];
    subjectprefix = InputStruct(is).run(1).subjectprefix;
    load([Subject_OutputDirIntermed '/res0_params/params' subjectprefix '.mat']);

    outdir = InputStruct(is).run(1).Output_nifti_file_path;
    prefix = InputStruct(is).run(1).Output_nifti_file_prefix;

    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------

    VV=load_untouch_nii(strcat(outdir,'/optimization_results/spms/rSPM_',prefix,'_CON_FIX_IND_sNorm.nii'));  
    
    if(is==1)
       if( rem(size(VV.img,4),3) >0 ) error('does not have equal conditions for CON, FIX, IND'); end
       Ncont = round( size(VV.img,4)/3 );
       %-- initialize
       spm_set = zeros( sum(mask(:)), numel(InputStruct), length(contr_num) ); % vox x subj x #contrasts
    end
    
    tmp = nifti_to_mat( VV,MM ); %% load volumes CON( contr ), FIX( contr ), IND( contr )
    
    % coolect spms, (vox x subj), contrasts stacked on 3rd dimension
    if( length(contr_num)==1 )
        spm_set(:,is) = tmp(:, (ppl_opt-1)*Ncont + contr_num );
    else
        for(i=1:length(contr_num))
            spm_set(:,is,i) = tmp(:, (ppl_opt-1)*Ncont + contr_num(i) );
        end
    end

    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
end

disp('reshaping...');

% .if >1 condition being analyzed reconcatenate spms into 2d matrix
% .creates design matrix of [condition, subj] for paired analysis
if( length(contr_num)>1 ) 
    
    if( length(contr_num)>2 && (strcmpi(analysis,'Ttest') || strcmpi(analysis,'Bootstrap')))
        error([analysis, ' can only compare between 2 contrasts max.']);
    end
    
    spm_set = reshape(spm_set, size(spm_set,1),[] ); 
    
    design2=design; %% start with original input matrix
    
    % now add contrast as design factor (plus per-run effect, paired measure)
    for( i=1:length(contr_num))
    design2  = [design2; [i*ones(numel(InputStruct),1), (1:numel(InputStruct))'] ];
    end
    
    design = design2; %% swap out the original "inmat" now
    clear design2;
end

%% now: run analysis

analysis_model.filepath = which(analysis);
% run analaysis --> ensure code is in path, call it from string
addpath( analysis_model.filepath );
pfun = str2func( analysis );

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
nii                   = VV;
nii.hdr.dime.datatype = 16;
nii.hdr.hist          = VV.hdr.hist;
nii.hdr.dime.dim(5)   = 4;
nii.hdr.hist.descrip  = CODE_PROPERTY.NII_HEADER;
nii.img               = zeros( [size(mask) 4] );

% predefined cluster size (based on FWHM=6 smoothing, standard task data, 2x2x2 mm res.)
MinClustSize = 95;

for(i=1:length(e_stat))
       
    map = output.(e_stat{i});
    pmap= output.([e_stat{i},'_p']);
    
    % (1) unthresholded
    tmp=mask; tmp(tmp>0)=map; 
    nii.img(:,:,:,1) = tmp;
    % (2) thresholded p<0.005 uncorr
    tmp=mask; tmp(tmp>0)=map .* double(pmap<=0.005); 
    nii.img(:,:,:,2) = tmp;
    % (3) cluster-size thresholded
    tmp=clust_up(tmp,MinClustSize);
    nii.img(:,:,:,3) = tmp;
    % (4) fdr threshold
    [pcrit thr0] = fdr( pmap, 'p', 0.05, 0 );
    tmp=mask; tmp(tmp>0)=map .* thr0; 
    nii.img(:,:,:,4) = tmp;
    
    % store thresholded map
    nii.img = TMPVOL;
    % save the nifti file
    save_untouch_nii(nii,strcat(PH_folder,'/',output.testname,'_',e{i},'_Contrast.',costr,'.nii'));
end
