function run_oppni_posthoc( datatype, inputfile, newmaskname, img_selector, analysis, design, contrast, rmeas, threshold, outlier_reject )


if(nargin<10) 
    outlier_reject=0; 
elseif(~isnumeric(outlier_reject))
    outlier_reject=str2num(outlier_reject);
end

% convert to string
if(isnumeric(threshold)) threshold = num2str(threshold); end

%% 1. loading and formatting design matrix inputs

% load, convert design matrix
if(ischar(design)   ) design = dlmread(design); end
if( numel(design)==1) design = [];              end
% check if enough subjects to analyze
if(size(design,1)>1 && size(design,1)<8) error('<8 subjects. Cannot do reliable analysis'); end

% check for repeated measures, reformat
if(~isempty(rmeas))
    
    rm=1;
    % strip off "run" and "subject" labels
    runvect = design(:,end-1);
    subvect = design(:,end);
    rmstr   = rmeas;

    if(size(design,2)==1) % check if both run + subject label cols
         error('for repeated measures, needs at least 2 columns: run list, subject list');
    elseif(size(design,2)==2) % special case where repeated-measures only
         design  = [];            
    else design  = design(:,1:end-2); % strip off rm columns
    end
else
    rm=0;
    runvect = [];
    subvect = [];
    rmstr   = '1run';
end

%dimensions
Nlev = size(design,2);
%check for constant columns, label binary (or continuous)
for(i=1:Nlev) 
    tmp= sort( unique(design(:,i)),'ascend') ;
    if(numel(tmp)==1) error(['design column #',num2str(i),' has only a single value!']); end
    if( numel(tmp)==2 && tmp(1)==0 && tmp(2)==1 ) 
         binlev(i,1) = 1; %binary
    else binlev(i,1) = 0; %non-binary
         if(numel(tmp)==3 && tmp(1)==0 && tmp(2)==1 && ~isfinite(tmp(3))) error('no missing values allowed in binary design'); end
% %          tmp=design(:,i); 
% %          if( sum( ~isfinite(tmp)>0 ) ) disp(['total of ',num2str(sum(~isfinite(tmp))),' missing values detected in col#',num2str(i),'. Performing median imputation...']); end
% %          tmp(~isfinite(tmp)) = median( tmp(isfinite(tmp)) );
% %          design(:,i) = tmp;
    end
end

if(isempty(contrast) || isempty(design) ) 
    costr         = '1way'; %contrast string
    design_newmat =[];
else
    %taking design+contrast, creating new design matrix for analysis
    costr         = contrast; %full contrast string
    contrast      = regexp(contrast,',','split'); % split into separate conditions
    design_newmat = zeros( size(design,1), numel(contrast) );
    
    for(i=1:numel(contrast)) %go through list of contrasts
        i,
        contrtemp = contrast{i};
        
        allcond   = regexp(contrtemp,'-|+','split'); %get all condition indices
        %convert to numeric indices for checking purposes
        clear allcix;
        for(j=1:numel(allcond)) allcix(j,1) = str2num(allcond{j}); end
        %check if levels are in bounds
        if( max(allcix)>Nlev ) error(['variable index #',num2str(max(allcix)),' exceeds total number of columns']); end
        %check if contrast mixes binary/continuous variables
        if( numel(unique( binlev(allcix) ))>1 ) 
            error('cannot mix binary and continuous values in contrast'); 
        else
            contr_type(i,1) = unique(binlev(allcix));
        end
        
        if(isempty(strfind(contrtemp,'-'))) %nonsubtractive
            
            contrtemp_sub=regexp(contrtemp,'+','split');
            for(j=1:numel(contrtemp_sub))
                design_newmat(:,i) = design_newmat(:,i) + design(:,str2num(contrtemp_sub{j}));
            end
            
        else % subtractive difference
            
            if( numel(strfind(contrtemp,'-'))>1 ) error('multiple - in contrast!'); end
            contrtemp = regexp(contrtemp,'-','split');
            %
            contrtemp_sub=regexp(contrtemp{1},'+','split');
            newreg=0;
            for(j=1:numel(contrtemp_sub))
                newreg=newreg + design(:,str2num(contrtemp_sub{j}));
            end
            if( contr_type(i,1)==1 ) design_newmat(:,i) = design_newmat(:,i) + double(newreg==j);
            else                     design_newmat(:,i) = design_newmat(:,i) + newreg;
            end
            contrtemp_sub=regexp(contrtemp{2},'+','split');
            newreg=0;
            for(j=1:numel(contrtemp_sub))
                newreg=newreg + design(:,str2num(contrtemp_sub{j}));
            end
            if( contr_type(i,1)==1 ) design_newmat(:,i) = design_newmat(:,i) - double(newreg==j);
            else                     design_newmat(:,i) = design_newmat(:,i) - newreg;
            end            
        end
    end
end

disp('checking for missing values...');
for(i=1:size(design_newmat,2))
    if( contr_type(i)~=1 && sum(~isfinite(design_newmat(:,i)))>0 )
        disp(['median-values imputated on col#',num2str(i),'.'])
        tmp = design_newmat(:,i);
        tmp(~isfinite(tmp))= median(tmp(isfinite(tmp)));
        design_newmat(:,i)=tmp;
    end
end

%% 2. loading and formatting neuroimaging data

if(     strcmpi(datatype,'fMRI') ) [InputStruct,~] = Read_Input_File(inputfile);
                                   subdir          = [];
elseif( strcmpi(datatype,'VBM')  ) [InputStruct]   = Read_Input_struct_vbm(inputfile);
                                   subdir          = '/struct_vbm';
elseif( strcmpi(datatype,'DTI')  ) [InputStruct]   = Read_Input_DTI(inputfile);
                                   subdir          = '/dti_processed';
elseif( strcmpi(datatype,'ASL')  ) [InputStruct]   = Read_Input_ASL(inputfile);
                                   subdir          = '/asl_processed';
else                               error('unrecognized datatype');
end

if(isempty(InputStruct)) error('cannot read inputs -- is your inputfile correctly defined??'); end
if(~isempty(design) && size(design,1)~=numel(InputStruct)) error('design size doesnt match the number of subjects in the inputfile'); end
% output path for results (same as group mask)
PH_folder = [InputStruct(1).run(1).Output_nifti_file_path,subdir,'/Posthoc_analyses'];
mkdir (PH_folder);
Nsubject = numel(InputStruct);

%================================================================== fMRI =%
if( strcmpi(datatype,'fMRI') )

    img_selectorSP = regexp(img_selector,',','split');
    contr_typ =img_selectorSP{1};
    anl_typ   =img_selectorSP{2};
    warp_typ  =img_selectorSP{3};
    ppl_opt   =img_selectorSP{4};
    contr_num =str2num(img_selectorSP{5});

    % pipeline optimization method - convert to numeric index
    if    ( ischar(ppl_opt) && strcmpi(ppl_opt,'PIPE1') ) ppl_opt=1;
    elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'CON') )   ppl_opt=1;
    elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'FIX') )   ppl_opt=2;
    elseif( ischar(ppl_opt) && strcmpi(ppl_opt,'IND') )   ppl_opt=3;
    end

        % name for group masks --> default should be empty
    if( isempty(newmaskname) )
        newmaskname = [InputStruct(1).run(1).Output_nifti_file_path, '/GroupMasks',warp_typ,'/group_consensus_mask.nii'];           
        M = load_nii(newmaskname);
    else
        [apath,aprefix,aext] = fileparts(newmaskname);
        mkdir_r(apath);
        %%% make sure mask-name terminates with .nii string
        newmaskname = [apath,'/',aprefix,'.nii'];
        M = load_nii(newmaskname);
    end

    % upload relevant SPMs
    for(is=1:numel(InputStruct))

        % load split_info data for further analysis
        Subject_OutputDirIntermed = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_metrics/'];
        subjectprefix = InputStruct(is).run(1).subjectprefix;
        load([Subject_OutputDirIntermed, contr_typ,'-',anl_typ,'/res0_params/params' subjectprefix '.mat']);

        outdir = InputStruct(is).run(1).Output_nifti_file_path;
        prefix = InputStruct(is).run(1).Output_nifti_file_prefix;

        V=load_nii(strcat(outdir,'/optimization_results/',contr_typ,'-',anl_typ,'/images_',warp_typ,'/rSPM_',prefix,'_CON_FIX_IND_sNorm.nii'));  

        if(is==1)
           if( rem(size(V.img,4),3) >0 ) error('does not have equal conditions for CON, FIX, IND'); end
           Ncont = round( size(V.img,4)/3 );
        end
        
        tmp = nifti_to_mat( V,M ); %% load volumes CON( contr ), FIX( contr ), IND( contr )
        vimg_smooth(:,is) = tmp(:, (ppl_opt-1)*Ncont + contr_num );
    end
    
%=================================================================== VBM =%    
elseif( strcmpi(datatype,'VBM') )

    % name for group masks --> default should be empty
    if( isempty(newmaskname) )
        grouptemp_path = [InputStruct(1).run(1).Output_nifti_file_path '/struct_vbm/group_template'];
        T = load_nii(sprintf('%s/template_%s_nl_symm.nii',grouptemp_path,img_selector));
        M = T;
        M.img = double(T.img>0.2);
    else
        [apath,aprefix,aext] = fileparts(newmaskname);
        mkdir_r(apath);
        %%% make sure mask-name terminates with .nii string
        newmaskname = [apath,'/',aprefix];
        M = load_nii(newmaskname);
    end

    disp('loading subjects...');
    clear vimg;
    for(is=1:Nsubject)
        reg_struct = [InputStruct(is).run(1).Output_nifti_file_path '/struct_vbm/',img_selector,'_warp/',InputStruct(is).run(1).Output_nifti_file_prefix{1}];
        V=load_nii([reg_struct,'_to_T3_mod.nii']);
        vimg(:,is) = nifti_to_mat(V,M);
    end
    disp('...done!');

    disp('now smoothing images...');
    for(is=1:Nsubject)

        tmp=M.img; tmp(tmp>0)=vimg(:,is);
        tmp=smooth3(tmp,'gaussian',[7 7 7],1.5);
        vimg_smooth(:,is) = tmp(M.img>0);
    end
    disp('...done!');
    
%=================================================================== DTI =%
elseif( strcmpi(datatype,'DTI') )
    
    % name for group masks --> default should be empty
    if( isempty(newmaskname) )
        T = load_nii(sprintf('%s/template_FA_nl_symm',grouptemp_path));
        M = T;
        M.img = double(T.img>0.2);
    else
        [apath,aprefix,aext] = fileparts(newmaskname);
        mkdir_r(apath);
        %%% make sure mask-name terminates with .nii string
        newmaskname = [apath,'/',aprefix];
        M = load_nii(newmaskname);
    end

    disp('loading subjects...');
    clear vimg;
    for(is=1:Nsubject)
        reg_struct = [InputStruct(is).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(is).run(1).Output_nifti_file_prefix,'/spat_norm/DTI_fit_',img_selector];
        V=load_nii([reg_struct,'_to_T3_mod.nii']);
        vimg(:,is) = nifti_to_mat(V,M);
    end
    disp('...done!');

    disp('now smoothing images...');
    for(is=1:Nsubject)

        tmp=M.img; tmp(tmp>0)=vimg(:,is);
        tmp=smooth3(tmp,'gaussian',[7 7 7],1.5);
        vimg_smooth(:,is) = tmp(M.img>0);
    end
    disp('...done!');

%=================================================================== ASL =%    
elseif( strcmpi(datatype,'ASL') )
    
    % name for group masks --> default should be empty
    if( isempty(newmaskname) )
        newmaskname = [InputStruct(1).run(1).Output_nifti_file_path, '/asl_processed/GroupMasks/group_consensus_mask.nii'];           
        M = load_nii(newmaskname);
    else
        [apath,aprefix,aext] = fileparts(newmaskname);
        mkdir_r(apath);
        %%% make sure mask-name terminates with .nii string
        newmaskname = [apath,'/',aprefix];
        M = load_nii(newmaskname);
    end

    disp('loading subjects...');
    clear vimg;
    for(is=1:Nsubject)
        reg_struct = [InputStruct(is).run(1).Output_nifti_file_path,'/asl_processed/',InputStruct(is).run(1).Output_nifti_file_prefix{1},'/proc_',img_selector];
        V=load_nii([reg_struct,'_sNorm.nii']);
        if( size(V.img,4)>1 ) error('posthoc ASL analysis can only be run on "mean" volumes'); end
        vimg_smooth(:,is) = nifti_to_mat(V,M);
    end
    disp('...done!');
end

mask =double(M.img);
    
%% 3. running analyses

% figure, imagesc( corr(vimg_smooth), [0 1] );

% if repeated measures, run on paired differences
if(rm>0)
    
    if(isempty(strfind(rmeas,'-'))) %nonsubtractive
        RUNV          = str2num(rmeas);
        %---
        % just takes design values, images from that run
        vimg_smooth   = vimg_smooth(:,runvect==RUNV);
        design_newmat = design_newmat(runvect==RUNV,:);
    else
        rmeas=regexp(rmeas,'-','split');
        RUNV    = [str2num(rmeas{1}), str2num(rmeas{2})];
        %---
        % go through list of unique subjects, get diff in images / design matrices
        subun=unique(subvect);     
        for(i=1:length(subun))
           ix1=find( subvect==subun(i) & runvect==RUNV(1) );
           ix2=find( subvect==subun(i) & runvect==RUNV(2) );
           vimg_smooth_diff(:,i) = vimg_smooth(:,ix1)-vimg_smooth(:,ix2);
           design_newmat_s1(i,:) = design_newmat(ix1,:);
           design_newmat_s2(i,:) = design_newmat(ix2,:);       
        end        
        % difference image update
        vimg_smooth = vimg_smooth_diff; clear vimg_smooth_diff;
        % design difference ...
        for(i=1:size(design_newmat,2))
           if( mean(abs( design_newmat_s1(:,i)-design_newmat_s2(:,i) )>eps ) ==0 )
               % no difference - fixed measures
               design_newmat_diff(:,i) = design_newmat_s1(:,i);
           else
               if( ( mean(abs( design_newmat_s1(:,i)-design_newmat_s2(:,i) )>eps ) <0.10 ) )
                   warning(['repeated-measures contrast #',num2str(i),' has <10% changed values...make sure there are no errors']);
               end
               % difference difference - changing measure
               design_newmat_diff(:,i) = design_newmat_s2(:,i) - design_newmat_s1(:,i);
           end
        end
        design_newmat = design_newmat_diff; clear design_newmat_diff;        
    end
end

% subvect = list of subject IDs (if paired-measures)
% binlev  = list of which design factors are binary (contr_type = same list for contrasts)
% costr   = string of contrasts
% design  = raw design matrix
% design_newmat = matrix of chosen contrasts

%%% call the model
if(isempty(strfind(analysis,'_phm')))
    analysis = [analysis,'_phm'];
end
analysis_model.filepath = which(analysis);
% run analaysis --> ensure code is in path, call it from string
addpath( analysis_model.filepath );
pfun = str2func( analysis );

if(strcmpi(analysis,'Ttest_phm') || strcmpi(analysis,'Bootstrap_phm') || strcmpi(analysis,'Splithalf_phm'))

    if(~isempty(design_newmat))
        if(size(design_newmat,2)>1) error('cannot run t-test or bootstrap on multiple design contrasts. Use a GLM!'); end
        if( contr_type ~=1 ) error('contrasts must be on binary class labels');
        end
        X = vimg_smooth(:, design_newmat==-1 | design_newmat==1 );
        y = design_newmat( design_newmat==-1 | design_newmat==1 );
    else
        X = vimg_smooth;
        y = [];
    end
    
    if( outlier_reject>0 )
        C2 = corr( X );
        % median correlation, converted to distance
        cmed_nn          = 1-median(C2,2);
        [thr_nn] = quick_gamma_test( cmed_nn );
        X = X(:, cmed_nn<thr_nn );
        y = y(cmed_nn<thr_nn,:);
        nout = sum(cmed_nn>thr_nn);
        
        if(nout==0) disp('no outliers found!');
        else disp(['found ',num2str(nout),' outliers. Subject lines #:']);
             disp(find(cmed_nn>=thr_nn));
        end
    else
        nout =0;
    end
    
    % run the analysis
    output = pfun( X, y );    
else
    X = vimg_smooth;
    y = design_newmat;
    
    if(strcmpi(analysis,'GLM'))
       ycheck = [ones(size(y,1),1) y];
       ycheck = bsxfun(@rdivide,ycheck,sqrt(sum(ycheck.^2)));
       cnum   = cond( ycheck );
       if    ( cnum > 1000 )   error(['GLM design matrix is multi-collinear! condition number=',num2str(cnum),'.']);
       elseif( cnum > 100  ) warning(['GLM design matrix is approaching multi-collinear! condition number=',num2str(cnum),'.']);
       else                     disp(['GLM design matrix is ok! condition number=',num2str(cnum),'.']);
       end
    end
    
    if( outlier_reject>0 )
        C2 = corr( X );
        % median correlation, converted to distance
        cmed_nn          = 1-median(C2,2);
        [thr_nn] = quick_gamma_test( cmed_nn );
        X = X(:, cmed_nn<thr_nn );
        y = y(cmed_nn<thr_nn,:);
        nout = sum(cmed_nn>thr_nn);
    else
        nout = 0;
    end

    output = pfun( X, y );        
end

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
nii                   = V;
nii.hdr.dime.datatype = 16;
nii.hdr.hist          = V.hdr.hist;
nii.hdr.dime.dim(5)   = size( output.(e_stat{1}),2 );
nii.img               = zeros( [size(mask) size( output.(e_stat{1}),2 )] );

% predefined cluster size (based on fixed smoothing params, 2x2x2 mm res.)
MinClustSize = 70;

for(i=1:numel(e_stat))
       
    map = output.(e_stat{i});
    pmap= output.([e_stat{i},'_p']);
    
    for(j=1:size(map,2))
        
        if    ( strcmpi(threshold,'0') || strcmpi(threshold,'no_thresh') )% (1) unthresholded
            tmp=mask; tmp(tmp>0)=map(:,j); 
            nii.img(:,:,:,j) = tmp;
        elseif( strcmpi(threshold,'1') || strcmpi(threshold,'uncorrected') ) % (2) thresholded p<0.005 uncorr
            tmp=mask; tmp(tmp>0)=map(:,j) .* double(pmap(:,j)<=0.005); 
            nii.img(:,:,:,j) = tmp;
        elseif( strcmpi(threshold,'2') || strcmpi(threshold,'cluster') ) % (3) cluster-size thresholded
            tmp=mask; tmp(tmp>0)=map(:,j) .* double(pmap(:,j)<=0.005); 
            tmp=clust_up(tmp,MinClustSize);
            nii.img(:,:,:,j) = tmp;
        elseif( strcmpi(threshold,'3') || strcmpi(threshold,'fdr') ) % (4) fdr threshold
            [pcrit thr0] = fdr( pmap(:,j), 'p', 0.05, 0 );
            tmp=mask; tmp(tmp>0)=map .* thr0; 
            nii.img(:,:,:,j) = tmp;
        else
            error('unrecognized thresholding option!');
        end
    end
    % save the nifti file
    if(outlier_reject>0) outlcat = ['OutlierCorrect_n',num2str(nout)];
    else                 outlcat = [];
    end
    strcat(PH_folder,'/',output.testname,'_',e{i},'_Tissue.',img_selector,'_Contrast.',costr,'_Rmeas.',rmstr,'_Thresh.',threshold,outlcat,'.nii'),
    save_nii(nii,strcat(PH_folder,'/',output.testname,'_',e{i},'_Tissue.',img_selector,'_Contrast.',costr,'_Rmeas.',rmstr,'_Thresh.',threshold,outlcat,'.nii'));
end

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

