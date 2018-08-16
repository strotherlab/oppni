function Pipeline_PART1(InputStruct, input_pipeset, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE, TPATTERN, TOFWHM, KEEPMEAN)
%
%==========================================================================
% PIPELINE_PART1 : main script used for running pipelines and obtaining
% performance metrics
%==========================================================================
%
% SYNTAX:
%
%   Pipeline_PART1(InputStruct, input_pipeset, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE, TPATTERN)
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
%  dospnormfirst = binary flag, determines if spatial normalization is done immediately after AFNI steps 
%                  instead of after optimization (Required for group-level analysis
%  DEOBLIQUE     = binary flag, determines if we correct for oblique imaging axes
%  TPATTERN      = multi-input flag, used to determine slice-timing, if this info is not in NIFTI headers
%                  0 or [] = read slice-timing from header
%                  1 or 'altplus' = default interleaved ascending
%                  'altminus' = interleaved descending
%                  'seqplus', 'seqminus' = sequential ascending, descending
%
% OUTPUT:
%
% Output:
%
%   three sets matfiles (+niftis, if input niiout=0):
%
%   <path>/intermediate_metrics/res1_spms/spms_<subjectprefix>.mat, contains the following:
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
%   <path>/intermediate_metrics/res2_temp/temp_<subjectprefix>.mat, contains the following:
%
%     TEMP_set        = BOLD timeseries associated with the activations in
%                       IMAGE_set. This is a (pipelines x 1) cell array, where
%                       each entry is a (time x components) timeseries matrix
%
%   <path>/intermediate_metrics/res3_stats/stats_<subjectprefix>.mat, contains the following:
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
% CODE_PROPERTY = '$Revision: 165 $';
% CODE_DATE    = '$Date: 2014-12-04 18:33:31 -0500 (Thu, 04 Dec 2014) $';
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

global CODE_PATH AFNI_PATH FSL_PATH MULTI_RUN_INPUTFILE CODE_PROPERTY
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART1.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
        if ~isdeployed
            addpath(CODE_PATH);
            addpath([CODE_PATH '/toolbox'])
            addpath([CODE_PATH '/NIFTI_tools'])
        end
        
    end
end
if isempty(AFNI_PATH) || isempty(FSL_PATH)
    read_settings;
end
if ~isempty(AFNI_PATH) && AFNI_PATH(end)~='/'
	AFNI_PATH = [AFNI_PATH '/'];
end
if ~isempty(FSL_PATH)  && FSL_PATH(end)~='/'
	FSL_PATH = [FSL_PATH '/'];
end

read_version;

%% Reading optional input arguments, or giving default assignments to variables

% check if nifti-output option is specified 
if nargin<5 || isempty(niiout)
    niiout=0;
elseif ischar(niiout)
    niiout = str2double(niiout);
end
if( niiout >0 ) 
    disp('WARNING: activation maps from all preprocessing pipelines are being created.')
    disp('This is not recommended for large numbers of pipelines, as you will quickly run out of disk space!');
end
% check if analysis model is specified
if isempty(analysis_model)
    analysis_model = 'NONE';
end
if( strcmp(analysis_model,'NONE') )
    disp('WARNING: no analysis model selected. You will only get preprocessed data in your output');
    disp('This is not recommended for large numbers of pipelines, as you will quickly run out of disk space!');
end
% check if task contrast is specified
if nargin<6 || isempty(contrast_list_str)
    contrast_list_str = 'NONE';
end
% check if spatial normalization is done first
if nargin<7 || isempty(dospnormfirst)
    dospnormfirst = false;
elseif( ischar(dospnormfirst) )
    dospnormfirst = str2double(dospnormfirst);   
end
% check if data needs to be "de-obliqued" (default = off)
if nargin<8 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
elseif strcmpi(DEOBLIQUE,'None')
    DEOBLIQUE = 0;
end
% check if slice-timing pattern is defined
if nargin<9 
    TPATTERN = [];
    fprintf('ERROR: User must specify the slice-timing pattern for fMRI data (TPATTERN) as part of inputs.\n');
    fprintf('TPATTERN options include:\n\t alt+z (or altplus)  = alternating in the plus direction\n\t alt+z2              = alternating, starting at slice #1 instead of #0\n\t alt-z (or altminus) = alternating in the minus direction\n\t alt-z2              = alternating, starting at slice #nz-2 instead of #nz-1\n\t seq+z (or seqplus)  = sequential in the plus direction\n\t seq-z (or seqminus) = sequential in the minus direction\n\n');
    error('Terminating OPPNI');
else
    tpatlist={'alt+z','altplus','alt+z2','alt-z','altminus','alt-z2','seq+z','seqplus','seq-z','seqminus','auto_hdr'};
%     compar=0; 
%     if( ischar(TPATTERN) )
%     for(i=1:numel(tpatlist)) 
%         if(strcmp(TPATTERN,tpatlist{i})) compar=1; end; 
%     end
%     end
%     if( compar == 0 )
%         fprintf('ERROR: Invalid slice-timing pattern (TPATTERN).\n');
%         fprintf('TPATTERN options include:\n\t alt+z (or altplus)  = alternating in the plus direction\n\t alt+z2              = alternating, starting at slice #1 instead of #0\n\t alt-z (or altminus) = alternating in the minus direction\n\t alt-z2              = alternating, starting at slice #nz-2 instead of #nz-1\n\t seq+z (or seqplus)  = sequential in the plus direction\n\t seq-z (or seqminus) = sequential in the minus direction\n\n');
%         error('Terminating OPPNI');
%     end        

    if ~ischar(TPATTERN) || ~ismember(TPATTERN, tpatlist)
        fprintf('ERROR: Invalid slice-timing pattern (TPATTERN).\n');
        fprintf('TPATTERN options include:\n\t alt+z (or altplus)  = alternating in the plus direction\n\t alt+z2              = alternating, starting at slice #1 instead of #0\n\t alt-z (or altminus) = alternating in the minus direction\n\t alt-z2              = alternating, starting at slice #nz-2 instead of #nz-1\n\t seq+z (or seqplus)  = sequential in the plus direction\n\t seq-z (or seqminus) = sequential in the minus direction\n\n');
        error('Terminating OPPNI');
    end
end
if nargin<10 || isempty(TOFWHM) || ~exist('TOFWHM','var')
    TOFWHM = 0;
end

if ischar(TOFWHM)
    TOFWHM = str2double(TOFWHM);
end

if ~ismember(TOFWHM,[0 1])
    TOFWHM = 0;
end

if nargin<11 || isempty(KEEPMEAN)
    KEEPMEAN = 0;
else
    if ischar(KEEPMEAN)
        KEEPMEAN = str2num(KEEPMEAN);
    end    
end

% %%========== TEMPORARY: uncomment this option to turn on BlurToFWHM
% 
% TOFWHM=1; %% --> this option force on BlurToFWHM smoothing
% 
% %%========== TEMPORARY: uncomment this option to turn on BlurToFWHM

%% Read Inputfiles
if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

%% acquire a list of all pipeline choices
[pipeset_half, detSet, mprSet, tskSet, phySet, gsSet, lpSet, Nhalf, Nfull] = get_pipe_list(input_pipeset);
pipeset_full = zeros( Nfull, 10 );

if numel(InputStruct(1).run)>1
    MULTI_RUN_INPUTFILE = true;
    aligned_suffix = '_aligned';
else
    aligned_suffix = '';
end
    
%% check whether all specified files exist (func, struct, physio, split_info)
check_input_file_integrity(InputStruct,max(pipeset_half(:,3)),max(tskSet)); 
%% now, defining contrasts for analysis
InputStruct = interpret_contrast_list_str(InputStruct,modelparam,analysis_model,contrast_list_str);             % generate contrast list for each subject and run

%% run all AFNI-based preprocessing steps
Pipeline_PART1_afni_steps(InputStruct, pipeset_half, dospnormfirst,DEOBLIQUE,TPATTERN,TOFWHM );

spatial_normalization_noise_roi(InputStruct); % Transform user defined 

%% save generated split_info files
for ksub = 1:numel(InputStruct)
    for krun = 1:numel(InputStruct(ksub).run)
        mkdir_r([InputStruct(ksub).run(krun).Output_nifti_file_path '/intermediate_processed/split_info']);
        split_info = InputStruct(ksub).run(krun).split_info;
        save([InputStruct(ksub).run(krun).Output_nifti_file_path '/intermediate_processed/split_info/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '.mat'],'split_info','CODE_PROPERTY','-v7');
    end
end
clear split_info


%%
for ksub = 1:numel(InputStruct)
    
    subjectmask = InputStruct(ksub).run(1).subjectmask;
    Subject_OutputDirIntermed = [InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_metrics'];
    Subject_OutputDirOptimize = [InputStruct(ksub).run(1).Output_nifti_file_path '/optimization_results'];
    subjectprefix = InputStruct(ksub).run(1).subjectprefix;
    mkdir_r([Subject_OutputDirIntermed]);
    mkdir_r([Subject_OutputDirIntermed '/regressors/reg' subjectprefix]);
    kcount = 0;
    
    %%%% Check if output files already exist %%%%
    suffix='';
    chk0 = exist( [Subject_OutputDirIntermed,'/res0_params/params' subjectprefix '.mat'], 'file' );
    chk1 = exist( [Subject_OutputDirIntermed,'/res1_spms/spms', subjectprefix,suffix,'.mat'], 'file' );
    chk2 = exist( [Subject_OutputDirIntermed,'/res2_temp/temp', subjectprefix,suffix,'.mat'], 'file' );
    chk3 = exist( [Subject_OutputDirIntermed,'/res3_stats/stats',subjectprefix,suffix,'.mat'], 'file' );
    
    if( chk0 && chk1 && chk2 && chk3 ) %% skip processing for this subject
        
        disp(['skipped regression processing for ',subjectprefix,'.' ]);
        disp(['--> output files already exist!'])
    else
        disp(['regression processing for ',subjectprefix,'.' ]);

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

            warning off;
            %% a. create output directories for metrics
            mkdir_r(strcat(Subject_OutputDirIntermed,'/res0_params'));
            if( ~strcmpi(analysis_model,'NONE') )
            mkdir_r(strcat(Subject_OutputDirIntermed,'/res1_spms'));
            mkdir_r(strcat(Subject_OutputDirIntermed,'/res2_temp'));
            mkdir_r(strcat(Subject_OutputDirIntermed,'/res3_stats'));
            end
            %% b. create output directories for optimization results
            mkdir_r(strcat(Subject_OutputDirOptimize,'/processed'));
            if( ~strcmpi(analysis_model,'NONE') )
            mkdir_r(strcat(Subject_OutputDirOptimize,'/spms'));
            mkdir_r(strcat(Subject_OutputDirOptimize,'/matfiles'));
            end
            warning on;

            %%% ==== Step 2.2(b): load and prepare motion MPEs === %%%
            %
            % load MPEs and generate PCA subset (>=85% variance),
            % currently not recorded to output (record=0)
            mpe_instring  = strcat(outdir,'/intermediate_processed/mpe/',Output_nifti_file_prefix{krun},'_mpe'    );
            mpe_outstring = strcat(outdir,'/intermediate_processed/mpe/',Output_nifti_file_prefix{krun},'_mpe_PCs');
            %% GROUP: PCA of head motion parameters -- load into cell array
            motPCs{krun} = motion_to_pcs( mpe_instring, mpe_outstring, 0.85,0 );

            %
            % load "baseproc" (basic preprocessing) dataset for
            % (a) estimating characteristic head motion patterns
            % (b) estimating vascular and white-matter maps
            %
            xbase_string  = strcat(outdir,'/intermediate_processed/afni_processed/',Output_nifti_file_prefix{krun},'_baseproc',aligned_suffix,'.nii');
            VX            = load_untouch_nii(  xbase_string );
            vxmat{krun}   = nifti_to_mat(VX,MM);
            [Nvox Ntime]  = size(vxmat{krun});

            %%% ==== Step 2.2(d): simple task modelling === %%%

            HRFdesign{krun} = split_info.design_mat;

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
        clear vxmat volcel;

        % GROUP: get priors into 2d matrices
        split_info_set{1}.NN_weight_avg = outwt.NN_weight;
        split_info_set{1}.WM_weight_avg = outwm.WM_weight;
        split_info_set{1}.NN_mask_avg   = outwt.NN_mask;
        split_info_set{1}.WM_mask_avg   = outwm.WM_mask;
        split_info_set{1}.FXYZ_avg      = mean( FXYZ, 3 );        
        
        Xsignal = HRFdesign;
        Xnoise  = motPCs;

        for krun = 1:N_run
            % check if argument for vascular masking exists and is turned off
            if( isfield(split_info_set{1},'VASC_MASK') && split_info_set{1}.VASC_MASK==0 )
                 split_info_set{1}.spat_weight = ones(size(split_info_set{1}.NN_weight_avg)); %%replace w/ unit weights
            else split_info_set{1}.spat_weight = split_info_set{1}.NN_weight_avg;
            end
            split_info_set{1}.mask_vol    = mask;
        end

        save([Subject_OutputDirIntermed '/res0_params/params' subjectprefix '.mat'],'split_info_set','Xsignal','Xnoise','modelparam','CODE_PROPERTY','-v7');

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
                volname  = strcat( outdir,'/intermediate_processed/afni_processed/',Output_nifti_file_prefix{krun},'_',PipeHalfList,aligned_suffix,'.nii' );
                % load the nifti files
                VV = load_untouch_nii(  volname );
                % convert nifti volume into matfile
                volmat{krun}  = nifti_to_mat(VV,MM);

                % ------------------------------------------------------------------------------------------------

                % load the nifti files
                if ~isempty(Noise_ROI{krun})
                    [tmp,roi_name,ext] = fileparts(Noise_ROI{krun});
                    roi_full_name = [outdir '/intermediate_processed/noise_roi/' roi_name '.nii'];
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
                    for TASK = tskSet
                        for GS  = gsSet
                            for LP = lpSet

                                kall = kall + 1;
                                % full list of preprocessing choices for this pipeline step
                                pipeset_full(kall,:) = [pipeset_half(i,:) DET MPR TASK GS LP];

                                if N_run>1
                                  [IMAGE_set_0{kall},TEMP_set_0{kall},METRIC_set_0{kall},IMAGE_set_y{kall},TEMP_set_y{kall},METRIC_set_y{kall},modeltype] = apply_regression_step_group(volmat,PipeHalfList,DET,MPR,TASK,GS,LP,phySet, Xsignal, Xnoise, noise_roi, split_info_set, analysis_model, InputStruct(ksub).run(1).Output_nifti_file_path,subjectprefix,Contrast_List,VV,KEEPMEAN);
                                else
				                  [IMAGE_set_0{kall},TEMP_set_0{kall},METRIC_set_0{kall},IMAGE_set_y{kall},TEMP_set_y{kall},METRIC_set_y{kall},modeltype,aa] = apply_regression_step(volmat,PipeHalfList,DET,MPR,TASK,GS,LP,phySet, Xsignal{1}, Xnoise{1}, noise_roi{1}, split_info_set, analysis_model, InputStruct(ksub).run(1).Output_nifti_file_path,subjectprefix,Contrast_List,VV,KEEPMEAN);
                                end
                            end
                        end
                    end
                end
            end
        end

        % pipeline information
        pipechars = ['m' 'c' 'p' 't' 's' 'd' 'r' 'x' 'g' 'l' 'y'];
        pipenames = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK', 'GSPC1', 'LOWPASS', 'PHYPLUS'};    
        % matrix of pipeline designations
        if( length(phySet) > 1 ) pipeset = [ pipeset_full zeros(Nfull,1); pipeset_full ones(Nfull,1) ];
        elseif( phySet == 0 )    pipeset = [ pipeset_full zeros(Nfull,1) ];
        elseif( phySet == 1 )    pipeset = [ pipeset_full ones(Nfull,1) ];
        end

        if( ~strcmpi(analysis_model,'NONE') )

            %% Step 2.4: consolidate results for output

            % combine images/metrics for results with or without PHYCAA+
            if( length(phySet) > 1 )
                IMAGE_set  = [IMAGE_set_0; IMAGE_set_y];
                TEMP_set   = [TEMP_set_0; TEMP_set_y];
                METRIC_set = [METRIC_set_0; METRIC_set_y];

            elseif( phySet == 0 )

                IMAGE_set  = [IMAGE_set_0];
                TEMP_set   = [TEMP_set_0];
                METRIC_set = [METRIC_set_0];

            elseif( phySet == 1 )

                IMAGE_set  = [IMAGE_set_y];
                TEMP_set   = [TEMP_set_y];
                METRIC_set = [METRIC_set_y];
            end

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
                    [Ao,Bo,rho(p,1)] = canoncorr( IMAGE_set{n}(:,p), split_info_set{1}.FXYZ_avg );
                    %
                end
                % correlation with spatial derivatives (motion effects)
                METRIC_set{n}.artifact_prior.MOT_corr = rho;
                % GROUP: take mask as >50% subjects with declared white matter tissue
                METRIC_set{n}.artifact_prior.WM_zscor = sum( IMAGE_set{n}( split_info_set{1}.WM_mask_avg > 0.5, : ) > 0 )./ sum( split_info_set{1}.WM_mask_avg > 0.5 )';
                % fraction of voxels >0 (global signal effects)
                METRIC_set{n}.artifact_prior.GS_fract = sum( IMAGE_set{n}>0 )./Nvox;
            end

            % save output matfiles
            %
            suffix = '';
            if  ~exist('OCTAVE_VERSION','builtin')
                save(strcat(Subject_OutputDirIntermed,'/res1_spms/spms', subjectprefix,suffix,'.mat'),'IMAGE_set','modelparam','CODE_PROPERTY');
                save(strcat(Subject_OutputDirIntermed,'/res2_temp/temp', subjectprefix,suffix,'.mat'),'TEMP_set','modelparam','CODE_PROPERTY');
                save(strcat(Subject_OutputDirIntermed,'/res3_stats/stats',subjectprefix,suffix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','analysis_model','CODE_PROPERTY','modelparam','CODE_PROPERTY');
            else
                save(strcat(Subject_OutputDirIntermed,'/res1_spms/spms', subjectprefix,suffix,'.mat'),'IMAGE_set','modelparam','CODE_PROPERTY','-mat7-binary');
                save(strcat(Subject_OutputDirIntermed,'/res2_temp/temp', subjectprefix,suffix,'.mat'),'TEMP_set','modelparam','CODE_PROPERTY', '-mat7-binary');
                save(strcat(Subject_OutputDirIntermed,'/res3_stats/stats',subjectprefix,suffix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','analysis_model','CODE_PROPERTY','modelparam','CODE_PROPERTY', '-mat7-binary');
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
                    nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
                    save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/spms/rSPM_',subjectprefix,suffix,'_pipelines_all.nii'));
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
                        nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
                        save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/spms/rSPM_',subjectprefix,'_pipeline_',num2str(n),'_',num2str(p),'vols.nii'));
                    end
                end
            end

        end
    end
end

disp('OPPNI PART 1 CODE COMPLETE');

%%
function name = generate_pipeline_name(pipechars,pipeset)

for i = 1:length(pipechars)
    name(i*2-1) = pipechars(i);
    name(i*2) = pipeset(i);
end


function [IMAGE_set_0,TEMP_set_0,METRIC_set_0,IMAGE_set_y,TEMP_set_y,METRIC_set_y,modeltype] = apply_regression_step_group(volmat,PipeHalfList,DET,MPR,TASK,GS,LP,phySet, Xsignal, Xnoise, noise_roi, split_info_set, analysis_model,OutputDirPrefix,subjectprefix,Contrast_List,VV,KEEPMEAN)

global CODE_PROPERTY

Subject_OutputDirIntermed = [OutputDirPrefix '/intermediate_metrics'];
Subject_OutputDirOptimize = [OutputDirPrefix '/optimization_results'];

N_run = length(volmat);
EmptyCell = cell(1,N_run);
% build pipeline prefix name -- and define current noise matrix
if(DET==-1)
nomen=[PipeHalfList 'dA'];
else
nomen=[PipeHalfList 'd' num2str(DET)];  
end

for(n=1:N_run)
    if(KEEPMEAN>0) mean_volmat{n} = mean(volmat{n},2);
    else           mean_volmat{n} = 0;
    end
end

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
    Regressors{run_counter}.NOISEROI = []; 
    if( DET==-1 )
        Regressors{run_counter}.DET = 1+floor( (split_info_set{1}.TR_MSEC./1000) * (Ntime/2) ./ 150 );    
    else
        Regressors{run_counter}.DET      = DET;
    end    
    Regressors{run_counter}.GSPC1    = [];
    Regressors{run_counter}.PHYPLUS  = [];
    Regressors{run_counter}.TR = (split_info_set{1}.TR_MSEC./1000);
end
% General Linear Model regression:
out_vol_denoi  = apply_glm( volmat, Regressors);

%%% ==== Step 2.3(c): global signal removed ==== %%%
%
% if GlobalSignal regression is "on", estimate and regress from data
nomen=[nomen 'g' num2str(GS)];

% customreg included (GS=2,3)
if fix(GS/2)==1
    for( is=1:length(volmat) )
        ln = find(noise_roi>0);
        weight = noise_roi(ln)/sum(noise_roi(ln));
        Regressors{is}.NOISEROI = [Regressors{is}.NOISEROI  mean(bsxfun(@times,out_vol_denoi(ln,:),weight))'];
    end
    GS = mod(GS,2);
end

out_vol_denoi  = apply_glm( volmat, Regressors);
% global estimation included (GS=0,1)
if( GS == 1 )
    
    for( is=1:length(volmat) )
        % get PC components on vascular-masked data
        volmat_temp = bsxfun(@times,out_vol_denoi{is},split_info_set{1}.spat_weight);
        [vx sx temp]     = svd( volmat_temp'*volmat_temp );
        % regress out "global" PC component
        out                = GLM_model_fmri( out_vol_denoi{is}, 0, vx(:,1), Xsig_curr{is}, 'econ' );
        out_vol_denoi{is}  = out.vol_denoi;
        Regressors{is}.GSPC1 = vx(:,1);
    end
    %
end

% low-pass filtering
if( LP == 1 )
    %
    nomen=[nomen 'l1'];
    %filters above 0.10 Hz
    for( is=1:length(volmat) )
        [ out_vol_filt{is} ] = quick_lopass( out_vol_denoi{is}, (split_info_set{1}.TR_MSEC./1000) );
    end
else
    out_vol_filt = out_vol_denoi;
    nomen=[nomen 'l0'];
end

%% ANALYSIS I
%%
%%% ==== Step 2.3(e): run analysis with multiple contrasts==== %%%

if( ~strcmpi(analysis_model,'NONE') )

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
        output_temp = group_analyses_wrapper( out_vol_filt, split_info_set, analysis_model );
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
        TMPVOL = zeros( [size(split_info_set{is}.mask_vol), size(out_vol_filt{is},2)] );
        % vascular weighting
        for(p=1:size(out_vol_filt{is},2) )
            tmp=split_info.mask_vol;tmp(tmp>0)= (out_vol_filt{is}(:,p) .* split_info_set{1}.spat_weight) + mean_volmat{is};
            TMPVOL(:,:,:,p) = tmp;
        end
        
        if(KEEPMEAN>0)
            nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' subjectprefix,'_run',num2str(is) '-'  [nomen,'y0']   ];
        else
            nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' subjectprefix,'_run',num2str(is) '-' [nomen,'y0']   ];
        end        

        nii=VV;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = VV.hdr.hist;
        nii.hdr.dime.dim(5) = size(vol_concat,2);
        nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
        save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/processed/Proc',subjectprefix,'_run',num2str(is),'_',nomen, 'y0.nii'));

        %% declare empty datasets
        IMAGE_set_0  = [];
        TEMP_set_0   = [];
        METRIC_set_0 = [];
    end
end

save([Subject_OutputDirIntermed '/regressors/reg' subjectprefix  '/' nomen 'y0.mat'],'Regressors','CODE_PROPERTY','-v7');
METRIC_set_0.cond_struc = design_cond(volmat,Regressors);

%% PHYCAA+ option
%%
if( ~isempty(find( phySet == 1 )) ) % perform if PHYCAA+ is being tested
    
    %%% ==== Step 2.3(f): PHYCAA+ physiological regression ==== %%%
    %
    %==============================================
    taskInfo.physio_map = split_info_set{1}.NN_weight_avg;   % include vascular prior
    taskInfo.task_SPMs  = IMAGE_set_0; % include reference SPM (Edited by babak)
    taskInfo.comp_crit  = 0;                 % less conservative threshold for noise components
    taskInfo.out_format =-1;                 % outputs regressed ONLY dataset
    %==============================================
    
    % run phycaa+
    Q=PHYCAA_plus_step2( out_vol_denoi, taskInfo );
    % get "denoised" + downweighted matrices out
    volmat_regress = Q.dataMat_denoised;
    for is = 1:N_run
        Regressors{is}.PHYPLUS = Q.Physio_Tset{is};
    end
    
    %% ANALYSIS II (if phycaa+ turned on)

    % low-pass filtering
    if( LP == 1 )
        %
        %filters above 0.10 Hz
        for( is=1:length(volmat) )
            [ out_vol_filt{is} ] = quick_lopass( volmat_regress{is}, (split_info_set{1}.TR_MSEC./1000) );
        end
    else
        out_vol_filt = volmat_regress;
    end
    
if( ~strcmpi(analysis_model,'NONE') )

    for contrast_counter = 1:size(Contrast_List,1)
        if isfield(split_info_set{1},'group')
            for k = 1:length(split_info_set)
                split_info_set{k}.idx_cond1 = split_info_set{k}.group.idx_cond(contrast_counter,1).sp;
                split_info_set{k}.idx_cond2 = split_info_set{k}.group.idx_cond(contrast_counter,2).sp;
            end
        end
        output_temp = group_analyses_wrapper( out_vol_filt, split_info_set, analysis_model );
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
        TMPVOL = zeros( [size(split_info_set{is}.mask_vol), size(out_vol_filt{is},2)] );

        for(p=1:size(out_vol_filt{is},2) )
            tmp=split_info.mask_vol;tmp(tmp>0)= (out_vol_filt{is}(:,p) .* split_info_set{1}.spat_weight) + mean_volmat{is};
            TMPVOL(:,:,:,p) = tmp;
        end
        
        if(KEEPMEAN>0)
            nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' subjectprefix,'_run',num2str(is) '-'  [nomen,'y1']   ];
        else
            nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' subjectprefix,'_run',num2str(is) '-' [nomen,'y1']   ];
        end         

        nii=VV;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = VV.hdr.hist;
        nii.hdr.dime.dim(5) = size(vol_concat,2);
        nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
        save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/processed/Proc',subjectprefix,'_run',num2str(is),'_',nomen, 'y1.nii'));

        %% declare empty datasets
        IMAGE_set_y  = [];
        TEMP_set_y   = [];
        METRIC_set_y = [];
    end
end
    save([Subject_OutputDirIntermed '/regressors/reg' subjectprefix '/' nomen 'y1.mat'],'Regressors','CODE_PROPERTY','-v7');
    METRIC_set_y.cond_struc = design_cond(volmat,Regressors);
else
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end


function [IMAGE_set_0,TEMP_set_0,METRIC_set_0,IMAGE_set_y,TEMP_set_y,METRIC_set_y,modeltype] = apply_regression_step(volmat,PipeHalfList,DET,MPR,TASK,GS,LP,phySet,Xsignal, Xnoise, noise_roi, split_info_set, analysis_model,OutputDirPrefix, subjectprefix, Contrast_List,VV,KEEPMEAN)

% build pipeline prefix name -- and define current noise matrix
global CODE_PROPERTY

Subject_OutputDirIntermed = [OutputDirPrefix '/intermediate_metrics'];
Subject_OutputDirOptimize = [OutputDirPrefix '/optimization_results'];

volmat = volmat{1};
if(KEEPMEAN>0) mean_volmat = mean(volmat,2);
else           mean_volmat = 0;
end
Ntime = size(volmat,2);
split_info = split_info_set{1};

% build pipeline prefix name -- and define current noise matrix
if(DET==-1)
nomen=[PipeHalfList 'dA'];
else
nomen=[PipeHalfList 'd' num2str(DET)];  
end

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
    Xsig_sp1(:,std(Xsig_sp1)<(1e-4/Ntime))=[];
    Xsig_sp2(:,std(Xsig_sp2)<(1e-4/Ntime))=[]; 
else nomen=[nomen 'x0']; Xsig_sp1 = [];
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
if( DET==-1 )
    Regressors.DET = 1+floor( (split_info_set{1}.TR_MSEC./1000) * (Ntime/2) ./ 150 );    
else
    Regressors.DET     = DET;
end
Regressors.NOISEROI= [];
Regressors.GSPC1   = [];
Regressors.PHYPLUS = [];

%%% ==== Step 2.3(c): global signal removed ==== %%%
%
% if Global Signal regression is "on", estimate and regress from data
nomen=[nomen 'g' num2str(GS)];
% regression of customreg (GS=2,3)
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
% regression of global effect (GS=1,3)
if( GS == 1 )
    [out_sp]   = apply_glm(volmat,Regressors);
    out_sp1 = out_sp(:,1:ceil(size(out_sp,2)/2));
    out_sp2 = out_sp(:,ceil(size(out_sp,2)/2)+1:end);
    % (SPLIT-1)
    % get PC components on vascular-masked data
    volmat_temp  = bsxfun(@times,out_sp1,split_info.spat_weight);
    [vx sx temp] = svd( volmat_temp'*volmat_temp );
    nr1 = [vx(:,1);zeros(size(out_sp2,2),1)];
    
    % (SPLIT-2)
    % get PC components on vascular-masked data
    volmat_temp  = bsxfun(@times,out_sp2,split_info.spat_weight);
    [vx sx temp] = svd( volmat_temp'*volmat_temp );    
    nr2 = [zeros(size(out_sp1,2),1);vx(:,1)];
    Regressors.GSPC1 = [nr1 nr2];
end

%% ANALYSIS I
%%
% reconcatenate splits
vol_concat = apply_glm(volmat,Regressors);
% low-pass filtering
if( LP == 1 )
    %
    nomen=[nomen 'l1'];
    %filters above 0.10 Hz
    [ vol_filt ] = quick_lopass( vol_concat, (split_info_set{1}.TR_MSEC./1000) );
else
    vol_filt = vol_concat;
    nomen=[nomen 'l0'];
end

%%% ==== Step 2.3(e): run analysis with multiple contrasts==== %%%
% save('indiscriminate_dump.mat')

if( ~strcmpi(analysis_model,'NONE') )

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
       
        output_temp = run_analyses_wrapper( vol_filt, split_info, analysis_model );
        
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
    %% save files as full 4D .nii volumes (vol_filt)
    modeltype = [];
    TMPVOL = zeros( [size(split_info.mask_vol), size(vol_concat,2)] );

    for(p=1:size(vol_concat,2) )
        tmp=split_info.mask_vol;tmp(tmp>0)= (vol_filt(:,p).*split_info.spat_weight) + mean_volmat;
        TMPVOL(:,:,:,p) = tmp;
    end
    
    if(KEEPMEAN>0)
        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' subjectprefix '-'  [nomen,'y0']   ];
    else
        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' subjectprefix '-' [nomen,'y0']   ];
    end
    
    nii=VV;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = VV.hdr.hist;
    nii.hdr.dime.dim(5) = size(vol_filt,2);
    nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
    save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/processed/Proc',subjectprefix,'_',nomen, 'y0.nii'));
    
    %% declare empty datasets
    IMAGE_set_0  = [];
    TEMP_set_0   = [];
    METRIC_set_0 = [];
end
    
save([Subject_OutputDirIntermed '/regressors/reg' subjectprefix '/' nomen 'y0.mat'],'Regressors','CODE_PROPERTY','-v7');
METRIC_set_0.cond_struc = design_cond(volmat,Regressors);

%% PHYCAA+ option
%%
if( ~isempty(find( phySet == 1 )) ) % perform if PHYCAA+ is being tested
    
    % format matrices for phycaa+ --> re-used cell array structure
    volcell{1} = vol_concat(:,1:ceil(size(vol_concat,2)/2));
    volcell{2} = vol_concat(:,ceil(size(vol_concat,2)/2)+1:end);
        
    %%% ==== Step 2.3(f): PHYCAA+ physiological regression ==== %%%
    %
    %==============================================
    taskInfo.physio_map = split_info.NN_weight_avg;   % include vascular prior
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
    if( LP == 1 )
        %
        %filters above 0.10 Hz
        [ vol_filt ] = quick_lopass( vol_concat, (split_info_set{1}.TR_MSEC./1000) );
    else
        vol_filt = vol_concat;
    end    

if( ~strcmpi(analysis_model,'NONE') )
    
    for contrast_counter = 1:size(Contrast_List,1)
        if isfield(split_info,'single')
            split_info.idx_cond1_sp1 = split_info.single.idx_cond(contrast_counter,1).sp1;
            split_info.idx_cond1_sp2 = split_info.single.idx_cond(contrast_counter,1).sp2;
            split_info.idx_cond2_sp1 = split_info.single.idx_cond(contrast_counter,2).sp1;
            split_info.idx_cond2_sp2 = split_info.single.idx_cond(contrast_counter,2).sp2;
        end
                       
        output_temp = run_analyses_wrapper( vol_filt, split_info, analysis_model );
        
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
    %% save files as full 4D .nii volumes (vol_filt)
    modeltype = [];
    TMPVOL = zeros( [size(split_info.mask_vol), size(vol_filt,2)] );

    for(p=1:size(vol_concat,2) )
        tmp=split_info.mask_vol;tmp(tmp>0)= (vol_filt(:,p).*split_info.spat_weight) + mean_volmat;
        TMPVOL(:,:,:,p) = tmp;
    end
    
    if(KEEPMEAN>0)
        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' subjectprefix '-'  [nomen,'y1']   ];
    else
        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' subjectprefix '-' [nomen,'y1']   ];
    end    
    
    nii=VV;
    nii.img = TMPVOL;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = VV.hdr.hist;
    nii.hdr.dime.dim(5) = size(vol_filt,2);
    nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
    save_untouch_nii(nii,strcat(Subject_OutputDirOptimize,'/processed/Proc',subjectprefix,'_',nomen, 'y1.nii'));
    
    %% declare empty datasets
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end
    
    save([Subject_OutputDirIntermed '/regressors/reg' subjectprefix '/' nomen 'y1.mat'],'Regressors','CODE_PROPERTY','-v7');
    METRIC_set_y.cond_struc = design_cond(volmat,Regressors);
else
    IMAGE_set_y  = [];
    TEMP_set_y   = [];
    METRIC_set_y = [];
end

