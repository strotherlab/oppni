function Pipeline_PART1(InputStruct, input_pipeset, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst)
%
%==========================================================================
% PIPELINE_PART1 : second piece of the optimized preprocessing
% pipeline, specifically for group-level (multi-subject) analysis models
% Runs specified set of pipeline steps for all subjects, and performs
% analysis with chosen model. Multi-model: user can choose from  a
% list of different analysis models (see below).
%==========================================================================
%
% SYNTAX:
%
%   Pipeline_PART1( InputStruct, input_pipeset, analysis_model, niiout, contrast_list_str )
%
% INPUT:
%
%   InputStruct    = string specifying "input" textfile (path/name),
%                   containing subject information
%   input_pipeset  = string specifying "pipeline" textfile (path/name),
%                   listing all preprocessing pipelines to test
%   analysis_model= string specifying choice of pipeline analysis model.
%                   Choices include:
%
%                     'LDA'  : linear discriminant (2-class block design)
%
%                   specify model parameters (task onsets, split structure etc)
%                   as fields in the "split_info" structure; check out
%                   'help run_analyses_wrapper' for details on each model
%   niiout        = binary value specifying output format,
%                   0=matfile only, 1=nifti file outputs too
%  contrast_list_str(optional) = string that specifies the desired task contrasts, e.g. for the
%                                the 1st task vs the 2nd task and the 2nd task vs the 3rd task, you may use contrast_list_str='1-2,2-3'
%                                                                                   
%
% OUTPUT:
%
% Output:
%
%   three matfiles (+niftis, if input niiout=0):
%
%   <path>/results1_group_spms_<subjectprefix>.mat, contains the following:
%
%     IMAGE_set       = activation maps computed for each pipeline. This is
%                       a (pipelines x 1) cell array, where each entry is a
%                       (voxels x components) matrix of images
%     prior_brain_maps= structure with following fields:
%
%                       wm_mask  : (voxels x 1) binary vector of brain regions
%                                  with predominantly predominantly white matter
%                       wm_weight: (voxels x 1) map down-weighting voxels by
%                                  likelihood of white matter content
%                       nn_mask  : (voxels x 1) binary vector of estimated regions
%                                  with predominantly non-neuronal tissue
%                       nn_weight: (voxels x 1) map down-weighting voxels
%                                  by likelihood of non-neuronal tissue
%                       mot_deriv: (voxels x 3) matrix, maps of spatial derivatives
%                                  on X,Y,Z axes, used to estimate head motion;
%
%   <path>/results2_group_temp_<subjectprefix>.mat, contains the following:
%
%     TEMP_set        = BOLD timeseries associated with the activations in
%                       IMAGE_set. This is a (pipelines x 1) cell array, where
%                       each entry is a (time x components) timeseries matrix
%
%   <path>/results3_group_stats_<subjectprefix>.mat, contains the following:
%
%     METRIC_set      = cell array with one entry per pipeline, providing
%                       performance metrics for a given analysis model. See
%                       'help run_analyses_wrapper' for details on metrics
%                       e.g. the LDA model produces metrics of prediction (P),
%                       reproducibility (R) and negative of distance (Dneg)
%                       from (P=1,R=1). Thus for pipeline q, we have entries:
%
%                       METRIC_set{q}.P, METRIC_set{q}.R, METRIC_set{q}.Dneg
%
%     artifact_priors = structure with following fields:
%
%                       motionCorr: (pipelines x 1) cell array. Each cell contains a
%                                   vector of correlation for IMAGE_set maps with brain-edge
%                                   derivatives from "prior_brain_maps.mot_deriv" indicating
%                                   motion artifact (usually significant at rho>0.30)
%                       wmfract   : (pipelines x 1) cell array. Each cell contains a vector
%                                   corresponding to IMAGE_set maps, where values are the
%                                   fraction of WM voxels with positive activation, indicating
%                                   potential white matter signal bias
%
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% version history: 1.00 2014/10/24
% ------------------------------------------------------------------------%


% Check Parameters
%%
global NUMBER_OF_CORES
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if ( ~exist('OCTAVE_VERSION','builtin') && exist('maxNumCompThreads') )
    maxNumCompThreads(NUMBER_OF_CORES);
end
%setenv('LD_LIBRARY_PATH','/usr/local/ge2011.11/lib/linux-x64:/opt/lib.exported:/opt/lib.exported')
% In some environment, like HPCVL, AFNI does not run properly since matlab changes library paths, the above line make it run.
%%

global CODE_PATH AFNI_PATH FSL_PATH MULTI_RUN_INPUTFILE
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART1.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
        addpath(CODE_PATH);
        addpath([CODE_PATH '/NIFTI_tools'])
    end
end
if isempty(AFNI_PATH) || isempty(FSL_PATH)
    read_settings;
end

version = '1.00 2014_10_24';

if ischar(niiout)
    niiout = str2double(niiout);
end
if nargin<6
    contrast_list_str = 'NONE';
end
if isempty(contrast_list_str)
        contrast_list_str = 'NONE';
end
if nargin<7
    dospnormfirst = false;
end
    
%% Read Inputfiles
if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

[pipeset_half, detSet, mprSet, tskSet, phySet, gsSet, Nhalf, Nfull] = get_pipe_list(input_pipeset);
pipeset_full = zeros( Nfull, 9 );

if numel(InputStruct(1).run)>1
    MULTI_RUN_INPUTFILE = true;
    aligned_suffix = '_aligned';
else
    aligned_suffix = '';
end

%% Interpret the contrast, and build the split-half data if neccessary
    
check_input_file_integrity(InputStruct,max(pipeset_half(:,3)),max(tskSet)); % check whether the analysis model and split info are matched.

%%
Pipeline_PART1_afni_steps(InputStruct,pipeset_half,dospnormfirst);

%%

InputStruct = interpret_contrast_list_str(InputStruct,modelparam,analysis_model,contrast_list_str);             % generate contrast list for each subject and run

spatial_normalization_noise_roi(InputStruct); % Transform user defined 

%%
% save generate split_info files

for ksub = 1:numel(InputStruct)
    for krun = 1:numel(InputStruct(ksub).run)
        mkdir_r([InputStruct(ksub).run(krun).Output_nifti_file_path '/split_info']);
        split_info = InputStruct(ksub).run(krun).split_info;
        save([InputStruct(ksub).run(krun).Output_nifti_file_path '/split_info/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '.mat'],'split_info','-v7');
    end
end
clear split_info

%%
for ksub = 1:numel(InputStruct)
    
    subjectmask = InputStruct(ksub).run(1).subjectmask;
    Subject_OutputDirectory = [InputStruct(ksub).run(1).Subject_OutputDirectory];
    subjectprefix = InputStruct(ksub).run(1).subjectprefix;
    mkdir_r([Subject_OutputDirectory]);
    mkdir_r([Subject_OutputDirectory '/regressors' subjectprefix]);
    kcount = 0;
    
    clear Output_nifti_file_path Output_nifti_file_prefix split_info_file
    N_run = numel(InputStruct(ksub).run);
    
    for krun = 1:numel(InputStruct(ksub).run)
        kcount = kcount + 1;
        Output_nifti_file_path{kcount}    = InputStruct(ksub).run(krun).Output_nifti_file_path;
        Output_nifti_file_prefix{kcount}  = InputStruct(ksub).run(krun).Output_nifti_file_prefix;
        Noise_ROI{kcount}                 = InputStruct(ksub).run(krun).Noise_ROI;
    end
    
    
    %% ==== Step 2.2(c): load brain mask, standard fmri data === %%%
    %
    % load in brain mask
    MM   = load_untouch_nii( subjectmask );
    mask = double(MM.img);
    
    %%%%% I. First Iteration through subjects -- load all prep. information,
    %%%%%    including MPEs, HRF, tissue maps etc
    
    krun=0;
    if N_run==1
         aligned_suffix = '';
    else aligned_suffix = '_aligned';
    end
    
    clear FXYZ
    for krun = 1:N_run
        
        %% Step 2.2: preparatory steps before pipeline testing
        
        %%% ==== 2.2(a) Read in output strings + task strings ==== %%%
        
        outdir   = Output_nifti_file_path{krun};
        
        % load in "split_info" structure with information about task onset and
        % analysis model parameters
        
        % load file
        split_info = InputStruct(ksub).run(krun).split_info;
        split_info_set{krun} = split_info;
        Contrast_List = InputStruct(ksub).run(krun).split_info.Contrast_List;
        
        %% GROUP: load an array of split_info files (one per subject)
        
        %% create output directory for matfiles
        warning off;
        mkdir_r( strcat(outdir,'/matfiles' ) );
        mkdir_r(strcat(Subject_OutputDirectory,'/niftis',subjectprefix));
        warning on;
        
        %%% ==== Step 2.2(b): load and prepare motion MPEs === %%%
        %
        % load MPEs and generate PCA subset (>=85% variance),
        % currently not recorded to output (record=0)
        mpe_instring  = strcat(outdir,'/mpe/',Output_nifti_file_prefix{krun},'_mpe'    );
        mpe_outstring = strcat(outdir,'/mpe/',Output_nifti_file_prefix{krun},'_mpe_PCs');
        %% GROUP: PCA of head motion parameters -- load into cell array
        motPCs{krun} = motion_to_pcs( mpe_instring, mpe_outstring, 0.85,0 );
        
        %
        % load "baseproc" (basic preprocessing) dataset for
        % (a) estimating characteristic head motion patterns
        % (b) estimating vascular and white-matter maps
        %
        xbase_string  = strcat(outdir,'/',Output_nifti_file_prefix{krun},'_baseproc',aligned_suffix,'.nii');
        VX            = load_untouch_nii(  xbase_string );
        vxmat{krun}   = nifti_to_mat(VX,MM);
        [Nvox Ntime]  = size(vxmat{krun});
        
        %%% ==== Step 2.2(d): simple task modelling === %%%
        %
        % used if pipeline step TASK=1, and there an overt task design to model.
        %
        % design_mat field was created using interpret_contrast_list_str.m
        % function in the new version, 
        if ~isfield(split_info,'design_mat') % In the new version of split_info the code generates design_mat from onsets.
            
            if (strcmp(split_info.type,'block') || strcmp(split_info.type,'multitask-block'))
                
                % build task-design vector from input information
                % this is a vector of signed values; -1=task condition1, 1=task condition2
                design = zeros(Ntime,size(Contrast_List,1));
                for contrast_counter = 1:size(Contrast_List,1)
%                     design(split_info.group.unbalanced_idx_cond(contrast_counter,1).sp, contrast_counter) = -1;
%                     design(split_info.group.unbalanced_idx_cond(contrast_counter,2).sp, contrast_counter) =  1;
                    design([ split_info.single.idx_cond(contrast_counter,1).sp1 split_info.single.idx_cond(contrast_counter,1).sp2], contrast_counter) = -1;
                    design([ split_info.single.idx_cond(contrast_counter,2).sp1 split_info.single.idx_cond(contrast_counter,2).sp2], contrast_counter) =  1;

                end
                
                
                % WARNING: MUST check whether design is full-rank
                
                %% GROUP: smooth with standard HRF function, into cell array
                HRFdesign_temp = design_to_hrf( design, (split_info.TR_MSEC/1000), [5.0 15.0] );
                HRFdesign{krun} = HRFdesign_temp(1:Ntime,:);
            elseif (strcmp( split_info.type,'event'))
                
                
                % building design vector: we want to subsample at 100 ms (faster than human RT)
                
                Nsubs  = max([1 round(split_info.TR_MSEC/100)]); % (#samples per TR) catch in case TR<100ms
                design = zeros( Ntime*Nsubs, size(split_info.cond,1));     % initialize design matrix
                % index allocates onsets to appropriate design-points
                for contrast_counter = 1:size(split_info.cond,1)
                    didx = unique(round( split_info.cond(contrast_counter).onsetlist./(split_info.TR_MSEC/Nsubs) ));
                    % catch + adjust overruns, setvalue=1 on design vector
                    didx(didx==0)          = 1;
                    didx(didx>Ntime*Nsubs) = Ntime*Nsubs;
                    design( didx , contrast_counter)         = 1;
                end
                % convolve with HRF (must convert into seconds!)
                
                HRFdesign_temp = design_to_hrf( design, (split_info.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
                % now, subsample back to get HRF at actual fMRI sampling rate
                HRFdesign{krun} = HRFdesign_temp( round(Nsubs/2): Nsubs : end );
                % -------------------- %
                
            else
                % GROUP: otherwise leave empty for each subject cell array
                HRFdesign{krun} = [];
            end
        else
            HRFdesign{krun} = split_info.design_mat;
        end
        
        %%% ==== Step 2.2(e): head motion spatial derivative maps ==== %%%
        %
        % used to detect large head motion artifact in brain maps
        % getting spatial derivative maps of brain, for motion artifact detection
        avg3d{krun}         =  mask;
        avg3d{krun}(mask>0) = mean(vxmat{krun},2);
        [fx fy fz]     = gradient( avg3d{krun} );
        % get the vector of derivatives on each axis
        %% GROUP: load into 3d matrix, across subjects
        FXYZ(:,:,krun) = [fx(mask>0) fy(mask>0) fz(mask>0)];
        outdir_stack{krun} = outdir;
    end
    
    
    %%% ==== Step 2.2(f): data-driven tissue segmentation/mapping
    
    avg3d_avg = 0;
    for krun = 1:N_run
        % INFO: for physiological (vascular) downweighting maps
        % remove mean/linear signal
        out   = GLM_model_fmri( vxmat{krun}, 1,[],[], 'econ' ); % just noise
        % cell formatting
        vxmat{krun} = out.vol_denoi;
        avg3d_avg = avg3d_avg + avg3d{krun};
    end
    avg3d_avg = avg3d_avg/N_run;
    
    dataInfo.TR            = (split_info_set{1}.TR_MSEC./1000);
    dataInfo.FreqCut       = 0.10;
    dataInfo.thresh_method = 'noprior';
    dataInfo.out_format    = 0;
    
    % estimate vascular map
    if N_run==1
        len2 = ceil(size(vxmat{1},2)/2);
        volcel{1} = vxmat{1}(:,1:len2);
        volcel{2} = vxmat{1}(:,len2+1:end);
        outwt      = PHYCAA_plus_step1( volcel, dataInfo );
        outwm      = WM_weight( volcel, dataInfo );
    else
        outwt      = PHYCAA_plus_step1( vxmat, dataInfo );
        outwm      = WM_weight( vxmat, dataInfo );
    end
    clear vxmat volcel
    
    % GROUP: get priors into 2d matrices
    NN_weight_avg = outwt.NN_weight;
    WM_weight_avg = outwm.WM_weight;
    NN_mask_avg   = outwt.NN_mask;
    WM_mask_avg   = outwm.WM_mask;
    if size(FXYZ,3)>1
         FXYZ_avg      = mean( FXYZ, 3 );
    else FXYZ_avg = FXYZ;
    end
    Xsignal = HRFdesign;
    Xnoise  = motPCs;
    
    for krun = 1:N_run
        split_info_set{krun}.spat_weight = NN_weight_avg;
        split_info_set{krun}.mask_vol    = mask;
    end
    
    save([Subject_OutputDirectory '/parameters' subjectprefix '.mat'],'Xsignal','Xnoise','FXYZ_avg','NN_mask_avg','WM_mask_avg','NN_weight_avg','WM_weight_avg','modelparam','-v7');
    save_untouch_nii(MM,[Subject_OutputDirectory '/masks' subjectprefix '.nii']);
    AVG = MM;
    AVG.img = avg3d_avg;
    AVG.hdr.dime.datatype=16;
    save_untouch_nii(AVG,[Subject_OutputDirectory '/mean' subjectprefix '.nii']);
    
    % initialize cell array for activation maps, one cell per pipeline
    IMAGE_set_0  = cell( Nfull, 1 );
    TEMP_set_0   = cell( Nfull, 1 );
    METRIC_set_0 = cell( Nfull, 1 );
    % if phycaa+ is being performed, define additional cell array for resulting images
    if( sum( phySet ) > 0 )
        IMAGE_set_y  = cell( Nfull, 1 );
        TEMP_set_y   = cell( Nfull, 1 );
        METRIC_set_y = cell( Nfull, 1 );
    end
    
    EmptyCell= cell(N_run,1);
    for(is=1:N_run) EmptyCell{is}=[]; end
    
    %% Step 2.3: performing pipeline testing
    %% run through pipeline options, generate output
    %  iterate through already-processed pipelines,
    %  load, run further processing and analyze...
    kall = 0;
    for i=1:Nhalf
        
        PipeHalfList = strcat( 'm',num2str( pipeset_half(i,1) ), ...
            'c',num2str( pipeset_half(i,2) ), ...
            'p',num2str( pipeset_half(i,3) ), ...
            't',num2str( pipeset_half(i,4) ), ...
            's',num2str( pipeset_half(i,5) )      );
        
        display([num2str(krun), ' --- ' PipeHalfList, ' doing:', num2str(i),'/',num2str(Nhalf)])
        
        %%%%% II. Second iteration level. Load all subjects for a given
        %%%%%     "pre-made" pipeline -- eg every pipeline made in Step-1
        
        for krun = 1:N_run
            
            Contrast_List = InputStruct(ksub).run(krun).split_info.Contrast_List;
            outdir   =  Output_nifti_file_path{krun};
            
            %%% ==== Step 2.3(a): load fmri volume
            %
            % specify volume and mask file names
            volname  = strcat( outdir,'/',Output_nifti_file_prefix{krun},'_',PipeHalfList,aligned_suffix,'.nii' );
            % load the nifti files
            VV = load_untouch_nii(  volname );
            % convert nifti volume into matfile
            volmat{krun}  = nifti_to_mat(VV,MM);
            
            % ------------------------------------------------------------------------------------------------
            
            % load the nifti files
            if ~isempty(Noise_ROI{krun})
                [tmp,roi_name,ext] = fileparts(Noise_ROI{krun});
                roi_full_name = [outdir '/noise_roi/' roi_name '.nii'];
                VV = load_untouch_nii(roi_full_name);
                noise_roi{krun} = nifti_to_mat(VV,MM);
            else
                noise_roi{krun}= [];
            end
            % 
        end
        
        %%%%% III. Third iteration level. run through each pipeline combination that
        %%%%%      is performed on pre-loaded data, analyze + save results
        
        %% run through additional proccessing choices
        for DET = detSet
            for MPR = mprSet
                for TASK= tskSet
                    for GS  = gsSet
                        kall = kall + 1;
                        if N_run>1
                            [IMAGE_set_0{kall},TEMP_set_0{kall},METRIC_set_0{kall},IMAGE_set_y{kall},TEMP_set_y{kall},METRIC_set_y{kall},pipeset_full(kall,:),modeltype] = apply_regression_step_group(volmat,PipeHalfList,pipeset_half(i,:),DET,MPR,TASK,GS,phySet, Xsignal, Xnoise, noise_roi, NN_weight_avg, split_info_set, analysis_model,Subject_OutputDirectory,subjectprefix,Contrast_List,VV);
                        else
                            [IMAGE_set_0{kall},TEMP_set_0{kall},METRIC_set_0{kall},IMAGE_set_y{kall},TEMP_set_y{kall},METRIC_set_y{kall},pipeset_full(kall,:),modeltype] = apply_regression_step(volmat,PipeHalfList,pipeset_half(i,:),DET,MPR,TASK,GS,phySet, Xsignal{1}, Xnoise{1}, noise_roi{1}, NN_weight_avg, split_info_set, analysis_model,Subject_OutputDirectory,subjectprefix,Contrast_List,VV);
                        end
                    end
                end
            end
        end
    end
    
    if( strcmp(analysis_model,'NONE') )
        
        % structure to save maps of motion / non-neuronal / white matter
        prior_brain_maps.wm_mask  =   WM_mask_avg;
        prior_brain_maps.wm_weight=   WM_weight_avg;
        prior_brain_maps.nn_mask  =   NN_mask_avg;
        prior_brain_maps.nn_weight=   NN_weight_avg;
        prior_brain_maps.mot_deriv=   FXYZ_avg;

        % save output matfiles
        %
        suffix = '';
        if  ~exist('OCTAVE_VERSION','builtin')
             save(strcat(Subject_OutputDirectory,'/results0_noanalysis', subjectprefix,suffix,'.mat'),'prior_brain_maps','modelparam', 'pipechars', 'pipenames', 'pipeset','modeltype','version');
        else save(strcat(Subject_OutputDirectory,'/results0_noanalysis', subjectprefix,suffix,'.mat'),'prior_brain_maps','modelparam', 'pipechars', 'pipenames', 'pipeset','modeltype','version', '-mat7-binary'););
        end

        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------
        
        disp('done preprocessing.');
    else
    
    %% Step 2.4: consolidate results for output
    
    % combine images/metrics for results with or without PHYCAA+
    if( length(phySet) > 1 )
        IMAGE_set  = [IMAGE_set_0; IMAGE_set_y];
        TEMP_set   = [TEMP_set_0; TEMP_set_y];
        METRIC_set = [METRIC_set_0; METRIC_set_y];
        pipeset    = [ pipeset_full zeros(Nfull,1); pipeset_full ones(Nfull,1) ];
        
    elseif( phySet == 0 )
        
        IMAGE_set  = [IMAGE_set_0];
        TEMP_set   = [TEMP_set_0];
        METRIC_set = [METRIC_set_0];
        pipeset    = [ pipeset_full zeros(Nfull,1) ];
        
    elseif( phySet == 1 )
        
        IMAGE_set  = [IMAGE_set_y];
        TEMP_set   = [TEMP_set_y];
        METRIC_set = [METRIC_set_y];
        pipeset    = [ pipeset_full ones(Nfull,1) ];
    end
    
    % pipeline information
    pipechars = ['m' 'c' 'p' 't' 's' 'd' 'r' 'x' 'g' 'y'];
    pipenames = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK',    'GSPC1', 'PHYPLUS'};
    
    % initializing cell arrays:
    rho_IMAGE    = cell(Nfull, 1);
    wmfrac_IMAGE = cell(Nfull, 1);
    % testing for spatial motion artifact and white matter bias
    %
    for(n=1:(Nfull*length(phySet)))   % for each pipeline entry...
        
        %initialize:
        rho = zeros( size(IMAGE_set{n},2), 1 );
        %
        for(p=1:size(IMAGE_set{n},2)) % and each image from set...
            % ------ test for motion (correlation with spatial derivatives
            [Ao,Bo,rho(p,1)] = canoncorr( IMAGE_set{n}(:,p), FXYZ_avg );
            %
        end
        
        % correlation with spatial derivatives (motion effects)
        METRIC_set{n}.artifact_prior.MOT_corr = rho;
        if N_run==1
            % average white matter z-score
            METRIC_set{n}.artifact_prior.WM_zscor = mean( IMAGE_set{n}(outwm.WM_mask > 0, : ) );
        else
            % GROUP: take mask as >50% subjects with declared white matter tissue
            METRIC_set{n}.artifact_prior.WM_zscor = sum( IMAGE_set{n}( WM_mask_avg > 0.5, : ) > 0 )./ sum( WM_mask_avg > 0.5 )';
        end
        % fraction of voxels >0 (global signal effects)
        METRIC_set{n}.artifact_prior.GS_fract = sum( IMAGE_set{n}>0 )./Nvox;
    end
    
    % structure to save maps of motion / non-neuronal / white matter
    prior_brain_maps.wm_mask  =   WM_mask_avg;
    prior_brain_maps.wm_weight=   WM_weight_avg;
    prior_brain_maps.nn_mask  =   NN_mask_avg;
    prior_brain_maps.nn_weight=   NN_weight_avg;
    prior_brain_maps.mot_deriv=   FXYZ_avg;
    
    % save output matfiles
    %
    suffix = '';
    if  ~exist('OCTAVE_VERSION','builtin')
        save(strcat(Subject_OutputDirectory,'/results1_spms', subjectprefix,suffix,'.mat'),'IMAGE_set','prior_brain_maps','modelparam');
        save(strcat(Subject_OutputDirectory,'/results2_temp', subjectprefix,suffix,'.mat'),'TEMP_set','modelparam');
        save(strcat(Subject_OutputDirectory,'/results3_stats',subjectprefix,suffix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','analysis_model','version','modelparam');
    else
        save(strcat(Subject_OutputDirectory,'/results1_spms', subjectprefix,suffix,'.mat'),'IMAGE_set','prior_brain_maps','modelparam','-mat7-binary');
        save(strcat(Subject_OutputDirectory,'/results2_temp', subjectprefix,suffix,'.mat'),'TEMP_set','modelparam', '-mat7-binary');
        save(strcat(Subject_OutputDirectory,'/results3_stats',subjectprefix,suffix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','analysis_model','version','modelparam', '-mat7-binary');
    end
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    
    if    (niiout>0)
        
        %---------- now, brain maps
        
        if( strcmp( modeltype, 'one_component') )
            
            disp('Only 1 image per pipeline. Concatenating all NIFTIS into single 4D matrix');
            
            TMPVOL = zeros( [size(mask), Nfull] );
            
            for(n=1:Nfull )
                tmp=mask;tmp(tmp>0)=IMAGE_set{n};
                TMPVOL(:,:,:,n) = tmp;
            end
            
            nii=VV;
            nii.img = TMPVOL;
            nii.hdr.dime.datatype = 16;
            nii.hdr.hist = VV.hdr.hist;
            nii.hdr.dime.dim(5) = Nfull;
            save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/Images',subjectprefix,suffix,'_pipelines_all.nii'));
        else
            
            disp('Multiple images per pipeline. Producing 4D volume for each pipeline');
            
            
            for(n=1:Nfull )
                
                TMPVOL = zeros( [size(mask), size(IMAGE_set{n},2)] );
                
                for(p=1:size(IMAGE_set{n},2) )
                    tmp=mask;tmp(tmp>0)=IMAGE_set{n}(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                
                nii=VV;
                nii.img = TMPVOL;
                nii.hdr.dime.datatype = 16;
                nii.hdr.hist = VV.hdr.hist;
                nii.hdr.dime.dim(5) = size(IMAGE_set{n},2);
                pipeline_name = generate_pipeline_name(pipechars,pipeset(n,:));
                save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/Images',subjectprefix,'_pipeline_',num2str(n),'_',num2str(p),'vols.nii'));
            end
        end
    end
    
    end
end

%%
%%
function name = generate_pipeline_name(pipechars,pipeset)
for i = 1:length(pipechars)
    name(i*2-1) = pipechars(i);
    name(i*2) = pipeset(i);
end

function [IMAGE_set_0,TEMP_set_0,METRIC_set_0,IMAGE_set_y,TEMP_set_y,METRIC_set_y,pipeset_full,modeltype] = apply_regression_step_group(volmat,PipeHalfList,pipeset_half,DET,MPR,TASK,GS,phySet, Xsignal, Xnoise, noise_roi, NN_weight_avg, split_info_set, analysis_model,Subject_OutputDirectory,subjectprefix,Contrast_List,VV)

N_run = length(volmat);
EmptyCell = cell(1,N_run);
% build pipeline prefix name -- and define current noise matrix
nomen=[PipeHalfList 'd' num2str(DET)];
% save list for labelling purposes later
pipeset_full = [pipeset_half DET MPR TASK GS];

%%% ==== Step 2.3(b): regression-based preprocessing ==== %%%
%
% add motion parameter timecourse as noise regressor
if( MPR==1 ) nomen=[nomen 'r1']; Xnoi_curr = Xnoise;
else         nomen=[nomen 'r0']; Xnoi_curr = EmptyCell;
end
% add signal timecourse as additional regressor
if(TASK==1 ) nomen=[nomen 'x1']; Xsig_curr = Xsignal;
else         nomen=[nomen 'x0']; Xsig_curr = EmptyCell;
end
%
for run_counter = 1:N_run
    Regressors{run_counter}.MP       = Xnoi_curr{run_counter};
    Regressors{run_counter}.Signal   = Xsig_curr{run_counter};
    Regressors{run_counter}.NOISEROI      = [];    
    Regressors{run_counter}.DET      = DET;
    Regressors{run_counter}.GSPC1    =   [];
    Regressors{run_counter}.PHYPLUS    = [];
end
% General Linear Model regression:
out_vol_denoi  = apply_glm( volmat, Regressors);

%%% ==== Step 2.3(c): global signal removed ==== %%%
%
% if GlobalSignal regression is "on", estimate and regress from data
nomen=[nomen 'g' num2str(GS)];

if fix(GS/2)==1
    for( is=1:length(volmat) )
        ln = find(noise_roi>0);
        weight = noise_roi(ln)/sum(noise_roi(ln));
        Regressors{is}.NOISEROI = [Regressors{is}.NOISEROI  mean(bsxfun(@times,out_vol_denoi(ln,:),weight))'];
    end
    GS = mod(GS,2);
end

out_vol_denoi  = apply_glm( volmat, Regressors);

if( GS == 1 )
    
    for( is=1:length(volmat) )
        % get PC components on vascular-masked data
        volmat_temp = bsxfun(@times,out_vol_denoi{is},NN_weight_avg);
        [vx sx temp]     = svd( volmat_temp'*volmat_temp );
        % regress out "global" PC component
        out                = GLM_model_fmri( out_vol_denoi{is}, 0, vx(:,1), Xsig_curr{is}, 'econ' );
        out_vol_denoi{is}  = out.vol_denoi;
        Regressors{is}.GSPC1 = vx(:,1);
    end
    %
end


%%% ==== Step 2.3(d): non-neuronal voxel downweighting ==== %%%
%
for( is=1:N_run)
    % vascular down-weighting applied to current, preprocessed data
    volmat_current{is} = bsxfun(@times,out_vol_denoi{is},NN_weight_avg);  % denoised volume x weight
end


%% ANALYSIS I
%%
%%% ==== Step 2.3(e): run analysis with multiple contrasts==== %%%

if( ~strcmp(analysis_model,'NONE') )

    for contrast_counter = 1:size(Contrast_List,1)
        if (strcmpi(split_info_set{1}.type,'block') || strcmpi(split_info_set{1}.type,'multitask-block'))
            for k = 1:length(split_info_set)
                split_info_set{k}.idx_cond1 = split_info_set{k}.group.idx_cond(contrast_counter,1).sp;
                split_info_set{k}.idx_cond2 = split_info_set{k}.group.idx_cond(contrast_counter,2).sp;
            end
        end
        if strcmpi(split_info_set{1}.type,'event')
            for k = 1:length(split_info_set)
                split_info_set{k}.onsetlist    = split_info_set{k}.cond(max(Contrast_List(contrast_counter,:))).onsetlist;
            end
        end
        output_temp = group_analyses_wrapper( volmat_current, split_info_set, analysis_model );
        if contrast_counter>1
            output.images = [output.images output_temp.images];
            output.temp   = [output.temp output_temp.temp];
            names = fieldnames(output_temp.metrics);
            for fld_counter = 1:length(names)
                output.metrics.(names{fld_counter}) = [output.metrics.(names{fld_counter}) output_temp.metrics.(names{fld_counter})];
            end
        else
            output = output_temp;
        end
    end
    %
    % determine if single-component or multicomponent
    modeltype = output.modeltype;
    if contrast_counter>1
        modeltype = 'multi_component';
        output.modeltype = modeltype;
    end
    % record optimal eigenimages / metrics
    IMAGE_set_0  = output.images;
    TEMP_set_0   = output.temp;
    METRIC_set_0 = output.metrics;
else
    modeltype = [];
    for( is=1:N_run)
        %% save files as full 4D .nii volumes (vol_concat)
        TMPVOL = zeros( [size(split_info_set{is}.mask_vol), size(volmat_current{is},2)] );

        for(p=1:size(volmat_current{is},2) )
            tmp=split_info.mask_vol;tmp(tmp>0)=volmat_current{is}(:,p);
            TMPVOL(:,:,:,p) = tmp;
        end

        nii=VV;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = VV.hdr.hist;
        nii.hdr.dime.dim(5) = size(vol_concat,2);
        save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/processed',subjectprefix,'_run',num2str(is),'_',nomen, 'y0.nii'));

        %% declare empty datasets
        IMAGE_set_0  = [];
        TEMP_set_0   = [];
        METRIC_set_0 = [];
    end
end

num_phyca_reg = 0;
save([Subject_OutputDirectory '/regressors' subjectprefix  '/' nomen 'y0.mat'],'Regressors','-v7');

%% PHYCAA+ option
%%
if( ~isempty(find( phySet == 1 )) ) % perform if PHYCAA+ is being tested
    
    %%% ==== Step 2.3(f): PHYCAA+ physiological regression ==== %%%
    %
    %==============================================
    taskInfo.physio_map = NN_weight_avg;   % include vascular prior
    taskInfo.task_SPMs  = IMAGE_set_0; % include reference SPM (Edited by babak)
    taskInfo.comp_crit  = 0;                 % less conservative threshold for noise components
    taskInfo.out_format = 1;                 % outputs downweight+regressed dataset
    %==============================================
    
    
    for( is=1:N_run )
        % rename the detrend/mpr'd/etc. matrices
        matt{is} = out_vol_denoi{is};
    end
    
    % run phycaa+
    Q=PHYCAA_plus_step2( matt, taskInfo );
    % get "denoised" + downweighted matrices out
    volmat_regress = Q.dataMat_denoised;
    for is = 1:N_run
        Regressors{is}.PHYPLUS = Q.Physio_Tset{is};
    end
    
    %% ANALYSIS II (if phycaa+ turned on)
    %%
    %%% ==== Step 2.3(g): run analysis ==== %%%
    
if( ~strcmp(analysis_model,'NONE') )

    for contrast_counter = 1:size(Contrast_List,1)
        if isfield(split_info_set{1},'group')
            for k = 1:length(split_info_set)
                split_info_set{k}.idx_cond1 = split_info_set{k}.group.idx_cond(contrast_counter,1).sp;
                split_info_set{k}.idx_cond2 = split_info_set{k}.group.idx_cond(contrast_counter,2).sp;
            end
        end
        output_temp = group_analyses_wrapper( volmat_regress, split_info_set, analysis_model );
        if contrast_counter>1
            output.images = [output.images output_temp.images];
            output.temp   = [output.temp output_temp.temp];
            names = fieldnames(output_temp.metrics);
            for fld_counter = 1:length(names)
                output.metrics.(names{fld_counter}) = [output.metrics.(names{fld_counter}) output_temp.metrics.(names{fld_counter})];
            end
        else
            output = output_temp;
        end
    end

    % record optimal eigenimages / metrics
    IMAGE_set_y  = output.images;
    TEMP_set_y   = output.temp;
    METRIC_set_y = output.metrics;
    
else
    modeltype = [];
    for( is=1:N_run)
        %% save files as full 4D .nii volumes (vol_concat)
        TMPVOL = zeros( [size(split_info_set{is}.mask_vol), size(volmat_current{is},2)] );

        for(p=1:size(volmat_current{is},2) )
            tmp=split_info.mask_vol;tmp(tmp>0)=volmat_current{is}(:,p);
            TMPVOL(:,:,:,p) = tmp;
        end

        nii=VV;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = VV.hdr.hist;
        nii.hdr.dime.dim(5) = size(vol_concat,2);
        save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/processed',subjectprefix,'_run',num2str(is),'_',nomen, 'y1.nii'));

        %% declare empty datasets
        IMAGE_set_y  = [];
        TEMP_set_y   = [];
        METRIC_set_y = [];
    end
end
    save([Subject_OutputDirectory '/regressors' subjectprefix '/' nomen 'y1.mat'],'Regressors','-v7');
else
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end

function [IMAGE_set_0,TEMP_set_0,METRIC_set_0,IMAGE_set_y,TEMP_set_y,METRIC_set_y,pipeset_full,modeltype] = apply_regression_step(volmat,PipeHalfList,pipeset_half,DET,MPR,TASK,GS,phySet,Xsignal, Xnoise, noise_roi, NN_weight_avg, split_info_set, analysis_model,Subject_OutputDirectory, subjectprefix, Contrast_List,VV)

% build pipeline prefix name -- and define current noise matrix

volmat = volmat{1};
Ntime = size(volmat,2);
split_info = split_info_set{1};

% build pipeline prefix name -- and define current noise matrix
nomen=[PipeHalfList 'd' num2str(DET)];
% save list for labelling purposes later
pipeset_full = [pipeset_half DET MPR TASK GS];

%%% ==== Step 2.3(b): regression-based preprocessing ==== %%%
%
% add motion parameter timecourse as noise regressor; separated into split-halves
if( MPR==1 ) 
    nomen=[nomen 'r1']; 
    Xnoi_sp1 = Xnoise(1:ceil(Ntime/2),:);
    Xnoi_sp2 = Xnoise(ceil(Ntime/2)+1:end,:);
else         
    nomen=[nomen 'r0'];
    Xnoi_sp1 = [];
    Xnoi_sp2 = [];
end
% add signal timecourse as additional regressor; separated into split-halves
if(TASK==1 ) 
    nomen=[nomen 'x1']; 
    Xsig_sp1 = Xsignal(1:ceil(Ntime/2)    ,:);
    Xsig_sp2 = Xsignal(ceil(Ntime/2)+1:end,:);
else         nomen=[nomen 'x0']; Xsig_sp1 = [];
    Xsig_sp2 = [];
end

if ~isempty(Xsig_sp1)
    sr1 = [Xsig_sp1;zeros(size(Xsig_sp2,1),size(Xsig_sp1,2))];
    sr2 = [zeros(size(Xsig_sp1,1),size(Xsig_sp2,2));Xsig_sp2];
    Regressors.Signal       = [sr1 sr2];
else
    Regressors.Signal = [];
end

if ~isempty(Xnoi_sp1)
    nr1 = [Xnoi_sp1;zeros(size(Xnoi_sp2,1),size(Xnoi_sp1,2))];
    nr2 = [zeros(size(Xnoi_sp1,1),size(Xnoi_sp2,2));Xnoi_sp2];
    Regressors.MP       = [nr1 nr2];
else
    Regressors.MP       = [];
end

% General Linear Model regression on splits:
Regressors.DET     = DET;
Regressors.NOISEROI= [];
Regressors.GSPC1   = [];
Regressors.PHYPLUS = [];


%%% ==== Step 2.3(c): global signal removed ==== %%%
%
% if Global Signal regression is "on", estimate and regress from data
nomen=[nomen 'g' num2str(GS)];

if fix(GS/2)==1
    [out_sp]   = apply_glm(volmat,Regressors);
    ln = find(noise_roi>0);
    weight = noise_roi(ln)/sum(noise_roi(ln));
    NOISEREG =  mean(bsxfun(@times,out_sp(ln,:),weight))';
    ind_sp1 = 1:ceil(size(out_sp,2)/2);
    ind_sp2 = ceil(size(out_sp,2)/2)+1:size(out_sp,2);    
    nr1 = [NOISEREG(ind_sp1,1);zeros(length(ind_sp2),1)];
    nr2 = [zeros(length(ind_sp1),1);NOISEREG(ind_sp2,1)];
    Regressors.NOISEROI = [nr1 nr2];
    GS = mod(GS,2);
end
if( GS == 1 )
    [out_sp]   = apply_glm(volmat,Regressors);
    out_sp1 = out_sp(:,1:ceil(size(out_sp,2)/2));
    out_sp2 = out_sp(:,ceil(size(out_sp,2)/2)+1:end);
    % (SPLIT-1)
    % get PC components on vascular-masked data
    volmat_temp  = bsxfun(@times,out_sp1,NN_weight_avg);
    [vx sx temp] = svd( volmat_temp'*volmat_temp );
    nr1 = [vx(:,1);zeros(size(out_sp2,2),1)];
    
    % (SPLIT-2)
    % get PC components on vascular-masked data
    volmat_temp  = bsxfun(@times,out_sp2,NN_weight_avg);
    [vx sx temp] = svd( volmat_temp'*volmat_temp );    
    nr2 = [zeros(size(out_sp1,2),1);vx(:,1)];
    Regressors.GSPC1 = [nr1 nr2];
end


%% ANALYSIS I
%%
% reconcatenate splits
vol_concat = apply_glm(volmat,Regressors);
out_sp1 = vol_concat(:,1:ceil(size(vol_concat,2)/2));
out_sp2 = vol_concat(:,ceil(size(vol_concat,2)/2)+1:end);
%%% ==== Step 2.3(e): run analysis with multiple contrasts==== %%%

if( ~strcmp(analysis_model,'NONE') )

    for contrast_counter = 1:size(Contrast_List,1)
        if strcmpi(split_info.type,'block') || strcmpi(split_info.type,'multitask-block')
            split_info.idx_cond1_sp1 = split_info.single.idx_cond(contrast_counter,1).sp1;
            split_info.idx_cond1_sp2 = split_info.single.idx_cond(contrast_counter,1).sp2;
            split_info.idx_cond2_sp1 = split_info.single.idx_cond(contrast_counter,2).sp1;
            split_info.idx_cond2_sp2 = split_info.single.idx_cond(contrast_counter,2).sp2;
        end
        if strcmpi(split_info.type,'event')
            split_info.onsetlist    = split_info.cond(max(Contrast_List(contrast_counter,:))).onsetlist;
        end
        output_temp = run_analyses_wrapper( vol_concat, split_info, analysis_model );
        if contrast_counter>1
            output.images = [output.images output_temp.images];
            output.temp   = [output.temp output_temp.temp];
            names = fieldnames(output_temp.metrics);
            for fld_counter = 1:length(names)
                output.metrics.(names{fld_counter}) = [output.metrics.(names{fld_counter}) output_temp.metrics.(names{fld_counter})];
            end
        else
            output = output_temp;
        end
    end

    % determine if single-component or multicomponent
    modeltype = output.modeltype;
    if contrast_counter>1
        modeltype = 'multi_component';
    end
    % record optimal eigenimages / metrics
    IMAGE_set_0  = output.images;
    TEMP_set_0   = output.temp;
    METRIC_set_0 = output.metrics;

else
    %% save files as full 4D .nii volumes (vol_concat)
    modeltype = [];
    TMPVOL = zeros( [size(split_info.mask_vol), size(vol_concat,2)] );

    for(p=1:size(vol_concat,2) )
        tmp=split_info.mask_vol;tmp(tmp>0)=vol_concat(:,p);
        TMPVOL(:,:,:,p) = tmp;
    end
    
    nii=VV;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = VV.hdr.hist;
    nii.hdr.dime.dim(5) = size(vol_concat,2);
    save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/processed',subjectprefix,'_',nomen, 'y0.nii'));
    
    %% declare empty datasets
    IMAGE_set_0  = [];
    TEMP_set_0   = [];
    METRIC_set_0 = [];
end
    
save([Subject_OutputDirectory '/regressors' subjectprefix '/' nomen 'y0.mat'],'Regressors','-v7');

%% PHYCAA+ option
%%
if( ~isempty(find( phySet == 1 )) ) % perform if PHYCAA+ is being tested
    
    % format matrices for phycaa+ --> re-used cell array structure
    volcell{1} = out_sp1;
    volcell{2} = out_sp2;
    
    %%% ==== Step 2.3(f): PHYCAA+ physiological regression ==== %%%
    %
    %==============================================
    taskInfo.physio_map = NN_weight_avg;   % include vascular prior
    taskInfo.task_SPMs  = IMAGE_set_0;     % include reference SPM
    taskInfo.comp_crit  = 0;                 % less conservative threshold for noise components
    taskInfo.out_format =-1;                 % outputs regressed ONLY dataset
    taskInfo.keepmean   = 0;                 % don't keep the mean
    %==============================================
    
    % run phycaa+
    Q=PHYCAA_plus_step2( volcell, taskInfo );
    % reconcatenate "denoised" matrices
    vol_concat = [Q.dataMat_denoised{1} Q.dataMat_denoised{2}];
    
    phycaa_reg1 = [Q.Physio_Tset{1};zeros(size(volcell{2},2),size(Q.Physio_Tset{1},2))];
    phycaa_reg2 = [zeros(size(volcell{1},2),size(Q.Physio_Tset{2},2));Q.Physio_Tset{2}];
    
    Regressors.PHYPLUS = [phycaa_reg1 phycaa_reg2];
    
    %% ANALYSIS II (if phycaa+ turned on)
    %%
    %%% ==== Step 2.3(g): run analysis ==== %%%
    
if( ~strcmp(analysis_model,'NONE') )
    
    for contrast_counter = 1:size(Contrast_List,1)
        if isfield(split_info,'single')
            split_info.idx_cond1_sp1 = split_info.single.idx_cond(contrast_counter,1).sp1;
            split_info.idx_cond1_sp2 = split_info.single.idx_cond(contrast_counter,1).sp2;
            split_info.idx_cond2_sp1 = split_info.single.idx_cond(contrast_counter,2).sp1;
            split_info.idx_cond2_sp2 = split_info.single.idx_cond(contrast_counter,2).sp2;
        end
        output_temp = run_analyses_wrapper( vol_concat, split_info, analysis_model );

        if contrast_counter>1
            output.images = [output.images output_temp.images];
            output.temp   = [output.temp output_temp.temp];
            names = fieldnames(output_temp.metrics);
            for fld_counter = 1:length(names)
                output.metrics.(names{fld_counter}) = [output.metrics.(names{fld_counter}) output_temp.metrics.(names{fld_counter})];
            end
        else
            output = output_temp;
        end
    end
        
    % record optimal eigenimages / metrics
    IMAGE_set_y  = output.images;
    TEMP_set_y   = output.temp;
    METRIC_set_y = output.metrics;
    
else
    %% save files as full 4D .nii volumes (vol_concat)
    modeltype = [];
    TMPVOL = zeros( [size(split_info.mask_vol), size(vol_concat,2)] );

    for(p=1:size(vol_concat,2) )
        tmp=split_info.mask_vol;tmp(tmp>0)=vol_concat(:,p);
        TMPVOL(:,:,:,p) = tmp;
    end
    
    nii=VV;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = VV.hdr.hist;
    nii.hdr.dime.dim(5) = size(vol_concat,2);
    save_untouch_nii(nii,strcat(Subject_OutputDirectory,'/niftis',subjectprefix,'/processed',subjectprefix,'_',nomen, 'y1.nii'));
    
    %% declare empty datasets
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end
    
    save([Subject_OutputDirectory '/regressors' subjectprefix '/' nomen 'y1.mat'],'Regressors','-v7');
else
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end
