function Pipeline_STEP2_MM_GROUP( inputfile, pipelinefile, groupmask, groupprefix, analysis_model, niiout )
%
%==========================================================================
% PIPELINE_STEP2_MM_GROUP : second piece of the optimized preprocessing 
% pipeline, specifically for group-level (multi-subject) analysis models
% Runs specified set of pipeline steps for all subjects, and performs 
% analysis with chosen model. MM = multi-model, user can choose from  a 
% list of different analysis models (see below).
%==========================================================================
%
% SYNTAX:    
%
%   Pipeline_STEP2_MM_GROUP( inputfile, pipelinefile, groupmask, groupprefix, analysis_model, niiout )
%            
% INPUT:
%
%   inputfile     = string specifying "input" textfile (path/name), 
%                   containing subject information
%   pipelinefile  = string specifying "pipeline" textfile (path/name), 
%                   listing all preprocessing pipelines to test
%   groupmask     = string specifying (path+name) of group brain mask
%   groupprefix   = string specifying prefix to insert in result output
%   analysis_model= string specifying choice of pipeline analysis model. 
%                   Choices include:
%
%                     'LDA'  : linear discriminant (2-class block design)
%                     'erCVA': event-related Canonical Variates Analysis
%
%                   specify model parameters (task onsets, split structure etc) 
%                   as fields in the "split_info" structure; check out 
%                   'help run_analyses_wrapper' for details on each model
%   niiout        = binary value specifying output format, 
%                   0=matfile only, 1=nifti file outputs too
%
% OUTPUT:
%
% Output:  
%
%   three matfiles (+niftis, if input niiout=0):
%            
%   <path>/results1_group_spms_<groupprefix>.mat, contains the following:
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
%   <path>/results2_group_temp_<groupprefix>.mat, contains the following:
%
%     TEMP_set        = BOLD timeseries associated with the activations in 
%                       IMAGE_set. This is a (pipelines x 1) cell array, where 
%                       each entry is a (time x components) timeseries matrix 
%
%   <path>/results3_group_stats_<groupprefix>.mat, contains the following:
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
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: 2013/10/31
% ------------------------------------------------------------------------%

addpath scripts_matlab;
addpath scripts_matlab/NIFTI_tools;
version = '2013_11_26';

%% ==== Step 2.1: generate list of pipelines from Step1 === %%%
%
% load the "pipelinefile" textfile, get information on pipelines run in Step1 
% .expressed in "pipeset_half": each row = pipeline performed in Step1, expressed as binary/integer values
%
% Flags {detSet,mprSet,tskSet,phySet,gsSet} dictate choices on the
% remaining pipeline steps, performed in this code
[pipeset_half, detSet, mprSet, tskSet, phySet, gsSet, Nhalf, Nfull] = get_pipe_list( pipelinefile );

% pre-initialize names for full pipelines (similar to pipeset_half)
% each row=full pipeline performed in this code, expressed as list of binary/integer values
 pipeset_full = zeros( Nfull, 9 );
 
%%% ==== Step 2.2(c): load brain mask, standard fmri data === %%%
%
% load in brain mask
MM   = load_untouch_nii( groupmask ); 
mask = double(MM.img);
  
%%%%% I. First Iteration through subjects -- load all prep. information,
%%%%%    including MPEs, HRF, tissue maps etc

% opens the inputfile (includes subject/dataset names that preprocessing is performed on...
fid = fopen(inputfile);
% read in first line
tline = fgetl(fid);
ksub=0;
while ischar(tline) % for every input line in textfile...

    ksub=ksub+1; % index #subjects(/datasets)
       
%% Step 2.2: preparatory steps before pipeline testing
    
    %%% ==== 2.2(a) Read in output strings + task strings ==== %%%

    % parse output directory
    ifile = strfind( tline, 'OUT=' ); ifile = ifile+4;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>ifile);
    fullline = tline(ifile:ips(1));
    isepr    = strfind( fullline, '/' );
    isepr    = isepr(end);
    prefix   = fullline(isepr+1:end);
    outdir   = fullline(1:isepr-1);
    % load in "split_info" structure with information about task onset and
    % analysis model parameters
    itask = strfind( tline, 'TASK=' ); itask = itask+5;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>itask);
    fullline = tline(itask:ips(1));    
    % load file
    load( fullline );
    
    %% GROUP: load an array of split_info files (one per subject)
    split_info_set{ksub} = split_info;
    
    %% create output directory for matfiles
    warning off;
    mkdir( strcat(outdir,'/matfiles' ) );
    warning on;
    
    %%% ==== Step 2.2(b): load and prepare motion MPEs === %%%
    %
    % load MPEs and generate PCA subset (>=85% variance),
    % currently not recorded to output (record=0)
    mpe_instring  = strcat(outdir,'/mpe/',prefix,'_mpe'    );
    mpe_outstring = strcat(outdir,'/mpe/',prefix,'_mpe_PCs');
    %% GROUP: PCA of head motion parameters -- load into cell array
    motPCs{ksub} = motion_to_pcs( mpe_instring, mpe_outstring, 0.85,0 );

    %
    % load "baseproc" (basic preprocessing) dataset for 
    % (a) estimating characteristic head motion patterns
    % (b) estimating vascular and white-matter maps
    %
    xbase_string  = strcat(outdir,'/',prefix,'_baseproc_sNorm.nii');
    VX     = load_untouch_nii(  xbase_string );
    vxmat  = nifti_to_mat(VX,MM);
    [Nvox Ntime] = size(vxmat);
    
    %%% ==== Step 2.2(d): simple task modelling === %%%
    %
    % used if pipeline step TASK=1, and there an overt task design to model. 
    %
    if   ( strcmp( split_info.type, 'block' ) )

        % build task-design vector from input information
        % this is a vector of signed values; -1=task condition1, 1=task condition2
        design = zeros(Ntime,1);
        design( split_info.idx_cond1(:) ) = -1;
        design( split_info.idx_cond2(:) ) =  1; 
        %% GROUP: smooth with standard HRF function, into cell array
        HRFdesign{ksub} = design_to_hrf( design, (split_info.TR_MSEC/1000), [5.0 15.0] );        
        
    elseif( strcmp( split_info.type, 'event' ) )
        
        [Nx Ny Nz] = size( VX.img(:,:,:,1) );
        % First: build task-design vector
        % design vector, init at finest subsampling rate TE ~ TR/Nz
        design = zeros(Ntime*Nz,1);
        % allocate onsets to appropriate design-points
        didx   = unique(round( split_info.onsetlist./(split_info.TR_MSEC/Nz) ));
        design(didx( (didx<Ntime*30) & (didx>0) )) = 1;
        % convolve with HRF (must convert into seconds!)
        HRFdesign{ksub} = design_to_hrf( design, (split_info.TR_MSEC/Nz)/1000, [5.0 15.0] );
        % now, subsample back to get HRF at actual fMRI sampling rate
        %% GROUP: load into cell array
        HRFdesign{ksub} = HRFdesign{ksub}( round(Nz/2): Nz : end );        
    else
        % GROUP: otherwise leave empty for each subject cell array
        HRFdesign{ksub} = [];
    end
    
    %%% ==== Step 2.2(e): head motion spatial derivative maps ==== %%%
    %
    % used to detect large head motion artifact in brain maps
    % getting spatial derivative maps of brain, for motion artifact detection
     avg3d         =  mask; 
     avg3d(mask>0) = mean(vxmat,2);
    [fx fy fz]     = gradient( avg3d );
    % get the vector of derivatives on each axis
    %% GROUP: load into 3d matrix, across subjects
	FXYZ(:,:,ksub) = [fx(mask>0) fy(mask>0) fz(mask>0)];
    
    %%% ==== Step 2.2(f): data-driven tissue segmentation/mapping
    %
    % INFO: for physiological (vascular) downweighting maps
    dataInfo.TR            = (split_info.TR_MSEC./1000);
    dataInfo.FreqCut       = 0.10;
    dataInfo.thresh_method = 'noprior';
    dataInfo.out_format    = 0;
    % remove mean/linear signal
    out   = GLM_model_fmri( vxmat, 1,[],[], 'econ' ); % just noise
    % cell formatting
    volcell{1} = out.vol_denoi(:,1:round(Ntime/2));
    volcell{2} = out.vol_denoi(:,round(Ntime/2)+1:end);
    % estimate vascular map
    outwt      = PHYCAA_plus_step1( volcell, dataInfo );
    %white matter tissue mask (NEW!!!)
    outwm      = WM_weight( volcell, dataInfo );
    
    % GROUP: get priors into 2d matrices
    NN_weight_set(:,ksub) = outwt.NN_weight;
    WM_weight_set(:,ksub) = outwm.WM_weight;
    %
    NN_mask_set(:,ksub)   = outwt.NN_mask;
    WM_mask_set(:,ksub)   = outwm.WM_mask;
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    
    % then get next line (subject)...
    tline = fgetl(fid);
end

fclose(fid);

kall=0;       %% indexing for each pipeline step
    
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

%%% GROUP averaging of the tissue templates
N_subject = ksub; %%total # subjects
% for prior weighting
NN_weight_avg = NN_group_average( NN_weight_set );
WM_weight_avg = NN_group_average( WM_weight_set );
% mask, just for interest
NN_mask_avg = mean( NN_mask_set,2 );
WM_mask_avg = mean( WM_mask_set,2 );
% for motion artifact
FXYZ_avg      = mean( FXYZ, 3 );

% NOW: copy a few structures into split_info fields
% ------------------------------------- %
% include vascular spatial weighting as split_info parameter, to integrate with analysis models
split_info_set{1}.spat_weight = outwt.NN_weight;
% include mask volume, e.g. for models requiring spatial characterization
split_info_set{1}.mask_vol    = mask;
% ------------------------------------- %
    
% define signal+noise matrices
Xsignal = HRFdesign;
Xnoise  = motPCs;
%
EmptyCell= cell(N_subject,1);
for(is=1:N_subject) EmptyCell{is}=[]; end

%% Step 2.3: performing pipeline testing
%% run through pipeline options, generate output
%  iterate through already-processed pipelines,
%  load, run further processing and analyze...
for(i=1:Nhalf) 

    PipeHalfList = strcat( 'm',num2str( pipeset_half(i,1) ), ...
                           'c',num2str( pipeset_half(i,2) ), ...
                           'p',num2str( pipeset_half(i,3) ), ...
                           't',num2str( pipeset_half(i,4) ), ...
                           's',num2str( pipeset_half(i,5) )      );

    [num2str(ksub), ' --- ' PipeHalfList, ' doing:', num2str(i),'/',num2str(Nhalf)],

    %%%%% II. Second iteration level. Load all subjects for a given
    %%%%%     "pre-made" pipeline -- eg every pipeline made in Step-1
 
    % opens the inputfile (includes subject/dataset names that preprocessing is performed on...
    fid = fopen(inputfile);
    % read in first line
    tline = fgetl(fid);
    ksub=0;
    while ischar(tline) % for every input line in textfile...

        ksub=ksub+1; % index #subjects(/datasets)

        %%% ==== 2.2(a) Read in output strings + task strings ==== %%%

        % parse output directory
        ifile = strfind( tline, 'OUT=' ); ifile = ifile+4;
        ips   = [strfind( tline, ' ' )-1 length(tline)];
        ips   = ips(ips>ifile);
        fullline = tline(ifile:ips(1));
        isepr    = strfind( fullline, '/' );
        isepr    = isepr(end);
        prefix   = fullline(isepr+1:end);
        outdir   = fullline(1:isepr-1);
        
        %%% ==== Step 2.3(a): load fmri volume
        %
        % specify volume and mask file names
        volname  = strcat( outdir,'/',prefix,'_',PipeHalfList,'_sNorm.nii' );
        % load the nifti files
        VV = load_untouch_nii(  volname );
        % convert nifti volume into matfile
        volmat{ksub}  = nifti_to_mat(VV,MM);

        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------

        % then get next line (subject)...
        tline = fgetl(fid);
    end

    fclose(fid);
        
    %%%%% III. Third iteration level. run through each pipeline combination that
    %%%%%      is performed on pre-loaded data, analyze + save results
    
    %% run through additional proccessing choices
    for( DET = detSet )
    for( MPR = mprSet )
    for( TASK= tskSet )
    for( GS  = gsSet  )

        % build pipeline prefix name -- and define current noise matrix
        nomen=[PipeHalfList 'd' num2str(DET)];

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
        for( is=1:N_subject )
            % General Linear Model regression:
            out  = GLM_model_fmri( volmat{is}, DET, Xnoi_curr{is}, Xsig_curr{is}, 'econ' );
            out_vol_denoi{is} = out.vol_denoi;
        end
        
        %%% ==== Step 2.3(c): global signal removed ==== %%%
        %
        % if GlobalSignal regression is "on", estimate and regress from data
        if( GS == 1 )
            
             for( is=1:N_subject )
                 % get PC components on vascular-masked data
                 volmat_temp = (out_vol_denoi{is}) .* repmat(NN_weight_avg, [1 size(out_vol_denoi{is},2)]);
                 [vx sx temp]     = svd( volmat_temp'*volmat_temp );
                 % regress out "global" PC component
                 out                = GLM_model_fmri( out_vol_denoi{is}, 0, vx(:,1), Xsig_curr{is}, 'econ' );
                 out_vol_denoi{is}  = out.vol_denoi;
             end
             %
             nomen=[nomen 'g1'];
        else nomen=[nomen 'g0'];
        end            

        kall=kall+1; %% increment

        % save list for labelling purposes later
        pipeset_full(kall,:) = [pipeset_half(i,:) DET MPR TASK GS];

%% ANALYSIS I
%%
        %%% ==== Step 2.3(e): run analysis ==== %%%
        output = group_analyses_wrapper( out_vol_denoi, split_info_set, analysis_model );

        % determine if single-component or multicomponent
        if(kall==1) modeltype = output.modeltype; end

        % record optimal eigenimages / metrics
        IMAGE_set_0{kall}  = output.images;
        TEMP_set_0{kall}   = output.temp;
        METRIC_set_0{kall} = output.metrics;

%% PHYCAA+ option
%%
        if( ~isempty(find( phySet == 1 )) ) % perform if PHYCAA+ is being tested

        %%% ==== Step 2.3(f): PHYCAA+ physiological regression ==== %%%
        %
        %==============================================
        taskInfo.physio_map = NN_weight_avg;   % include vascular prior
        taskInfo.task_SPMs  = IMAGE_set_0{kall}; % include reference SPM
        taskInfo.comp_crit  = 0;                 % less conservative threshold for noise components
        taskInfo.out_format =-1;                 % outputs ONLY regressed dataset
        taskInfo.keepmean   = 0;                 % don't keep the mean
        %==============================================
        % run phycaa+
        Q=PHYCAA_plus_step2( out_vol_denoi, taskInfo ); 
        % update "denoised" matrices
        out_vol_denoi = Q.dataMat_denoised;

%% ANALYSIS II (if phycaa+ turned on)
%%
        %%% ==== Step 2.3(g): run analysis ==== %%%
        output = group_analyses_wrapper( out_vol_denoi, split_info_set, analysis_model );
        % record optimal eigenimages / metrics
        IMAGE_set_y{kall}  = output.images;
        TEMP_set_y{kall}   = output.temp;
        METRIC_set_y{kall} = output.metrics;

        end
    end
    end
    end
    end
end

%% Step 2.4: consolidate results for output
    
    % combine images/metrics for results with or without PHYCAA+
    if( length(phySet) > 1 )
        %
        IMAGE_set  = [IMAGE_set_0; IMAGE_set_y];
        TEMP_set   = [TEMP_set_0; TEMP_set_y];
        METRIC_set = [METRIC_set_0; METRIC_set_y];
        %
        pipeset = [ pipeset_full zeros(Nfull,1); pipeset_full ones(Nfull,1) ];                
        
    elseif( phySet == 0 )    

        IMAGE_set  = [IMAGE_set_0];        
        TEMP_set   = [TEMP_set_0];        
        METRIC_set = [METRIC_set_0];
        %
        pipeset = [ pipeset_full zeros(Nfull,1) ];        

    elseif( phySet == 1 )
        
        IMAGE_set  = [IMAGE_set_y];
        TEMP_set   = [TEMP_set_y];        
        METRIC_set = [METRIC_set_y];
        %
        pipeset = [ pipeset_full zeros(Nfull,1) ];        
    end
    
    % pipeline information
    pipechars = ['m' 'c' 'p' 't' 's' 'd' 'r' 'x' 'g' 'y'];
    pipenames = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK',    'GSPC1', 'PHYPLUS'};

    % testing for spatial motion artifact and white matter bias
    %
    for(n=1:(Nfull*length(phySet)))   % for each pipeline entry...
        
        %initialize:
        rho = zeros( size(IMAGE_set{n},2), 1 );
        %
        for(p=1:size(IMAGE_set{n},2)) % and each image from set...
            %
            % ------ test for motion (correlation with spatial derivatives
            [Ao,Bo,rho(p,1)] = canoncorr( IMAGE_set{n}(:,p), FXYZ_avg );
            %
        end
        
        % correlation with spatial derivatives (motion effects)
        METRIC_set{n}.artifact_prior.MOT_corr = rho; 
        % average white matter z-score
        METRIC_set{n}.artifact_prior.WM_zscor = mean( IMAGE_set{n}(outwm.WM_mask > 0, : ) );
        % fraction of voxels >0 (global signal effects)
        METRIC_set{n}.artifact_prior.GS_fract = sum( IMAGE_set{n}>0 )./Nvox;
    end
                
    % structure to save maps of motion / non-neuronal / white matter
    prior_brain_maps.wm_mask  =   WM_mask_avg;
    prior_brain_maps.wm_weight= WM_weight_avg;
    prior_brain_maps.nn_mask  =   NN_mask_avg;
    prior_brain_maps.nn_weight= NN_weight_avg;
    prior_brain_maps.mot_deriv= FXYZ_avg;

    % save output matfiles
    %
    if( exist('OCTAVE_VERSION','builtin') )
         % matlab-compatible
        save(strcat(outdir,'/matfiles/results1_group_spms_', groupprefix,'.mat'),'IMAGE_set','prior_brain_maps', '-mat7-binary');
        save(strcat(outdir,'/matfiles/results2_group_temp_', groupprefix,'.mat'),'TEMP_set', '-mat7-binary');
        save(strcat(outdir,'/matfiles/results3_group_stats_',groupprefix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','version', '-mat7-binary');
    else 
        save(strcat(outdir,'/matfiles/results1_group_spms_', groupprefix,'.mat'),'IMAGE_set','prior_brain_maps');
        save(strcat(outdir,'/matfiles/results2_group_temp_', groupprefix,'.mat'),'TEMP_set');
        save(strcat(outdir,'/matfiles/results3_group_stats_',groupprefix,'.mat'),'METRIC_set', 'pipechars', 'pipenames', 'pipeset','modeltype','version');
    end

    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------        
    
    % if requested, saving eigenimages in nifti data format
    %
    if    (niiout>0)

        %---------- now, brain maps
        
        if( strcmp( modeltype, 'one_component') )
        
            disp('Only 1 image per pipeline. Concatenating all NIFTIS into single 4D matrix');
            
            TMPVOL = zeros( [size(mask), Nfull] );

            for(n=1:Nfull )
                tmp=mask;tmp(tmp>0)=IMAGE_set{n};
                TMPVOL(:,:,:,n) = tmp;
            end            
            %
            nii     = VV; % copy nifti struct
            nii.img = TMPVOL; % replace volume
            nii.hdr.dime.dim(5) = size(TMPVOL,4); % adjust for #timepoints
            %
            save_untouch_nii(nii,strcat(outdir,'/matfiles/Images_',groupprefix,'_pipelines_all.nii'));  
        else
            
            disp('Multiple images per pipeline. Producing 4D volume for each pipeline');

            % create output directory
            mkdir(strcat(outdir,'/matfiles/niftis_',prefix));
            
            for(n=1:Nfull )
            
                TMPVOL = zeros( [size(mask), size(IMAGE_set{n},2)] );

                for(p=1:size(IMAGE_set{n},2) )
                    tmp=mask;tmp(tmp>0)=IMAGE_set{n}(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                %
                nii     = VV; % copy nifti struct
                nii.img = TMPVOL; % replace volume
                nii.hdr.dime.dim(5) = size(TMPVOL,4); % adjust for #timepoints
                %
                save_untouch_nii(nii,strcat(outdir,'/matfiles/niftis_',groupprefix,'/Images_',groupprefix,'_pipeline_',num2str(n),'_',num2str(p),'vols.nii'))
            end
        end
    end
    
    
%%
%%

function [pipeset_half, detSet, mprSet, tskSet, phySet, gsSet, Nhalf, Nfull] = get_pipe_list( filename )

% reads in the inputfile
fid = fopen(filename);
newline = fgetl(fid);
pipelinelist=[];
while ischar(newline)
    %
    pipelinelist = [pipelinelist ' ' newline];
    newline      = fgetl(fid);
end
fclose(fid);

proclist = {'MOTCOR='; 'CENSOR='; 'RETROICOR='; 'TIMECOR=';'SMOOTH='; 'DETREND=';'MOTREG=';'TASK=';'PHYPLUS='; 'GSPC1='};
ileft  = strfind( pipelinelist, '[' );
iright = strfind( pipelinelist, ']' );
for(s=1:10)
    iStep = strfind( pipelinelist, proclist{s} );
    bleft       = ileft (ileft >iStep);   bleft= bleft(1);
    bright      = iright(iright>iStep);  bright=bright(1); 
    substr{s,1} = pipelinelist(bleft:bright);
end

% PIPE-1:Motion
pipeset_old=[];
pipeset_new=[]; K=1;
% 
if( ~isempty( strfind(substr{K},'0') ) )  pipeset_new = [pipeset_new; 0];   end
if( ~isempty( strfind(substr{K},'1') ) )  pipeset_new = [pipeset_new; 1];   end

% PIPE-2:Censor
pipeset_old = pipeset_new;
pipeset_new = []; K=2;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
   tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end

% PIPE-3:Retroicor
pipeset_old = pipeset_new;
pipeset_new = []; K=3;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
   tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end

% PIPE-4:Timecor
pipeset_old = pipeset_new;
pipeset_new = []; K=4;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
   tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
   pipeset_new=[pipeset_new; tmpset];
end

% PIPE-5:Smooth
pipeset_old = pipeset_new;
pipeset_new = []; K=5;
fulidx  = [1  strfind(substr{K},',')  length(substr{K})];
numscal = length(fulidx)-1;
for(i=1:numscal)
    scaltok = substr{K}( fulidx(i)+1:fulidx(i+1)-1 );
    tmpset  =[ pipeset_old, str2num(scaltok)*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end


pipeset_half = pipeset_new; % everything that was already done!!
Nhalf = size( pipeset_half, 1 );

% -------------------------------------------------------------------

% PIPE-6:Detrend
detSet = [];
K=6;
fulidx  = [1  strfind(substr{K},',')  length(substr{K})];
numscal = length(fulidx)-1;
for(i=1:numscal)
    scaltok = substr{K}( fulidx(i)+1:fulidx(i+1)-1 );
    detSet = [detSet str2num(scaltok)];
end
% PIPE-7:Motreg
mprSet = [];
K=7;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   mprSet     =[ mprSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
   mprSet     =[ mprSet, 1];
end

% PIPE-8:Taskreg
tskSet = [];
K=8;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   tskSet     =[ tskSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
   tskSet     =[ tskSet, 1];
end

% PIPE-9:phycaa+
phySet = [];
K=9;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   phySet     =[ phySet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
   phySet     =[ phySet, 1];
end

% PIPE-10:globalsig
gsSet = [];
K=10;
% 
if( ~isempty( strfind(substr{K},'0') ) )
   gsSet     =[ gsSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
   gsSet     =[ gsSet, 1];
end


Nfull = Nhalf * length(detSet) * length(mprSet) * length(tskSet) * length(gsSet);
