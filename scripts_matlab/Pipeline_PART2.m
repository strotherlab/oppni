function Pipeline_PART2( InputStruct, optimize_metric, mot_gs_control, process_out, keepmean)
%
%==========================================================================
% PIPELINE_PART2 : this step identifies optimal pipelines, and produces
% optimization results, plus figures. Note that script can only optimize
% single-componen analyses for now.
%==========================================================================
%
% SYNTAX:
%
%   Pipeline_PART2( InputStruct, optimize_metric, mot_gs_control, process_out, keepmean)
%
% INPUT:
%
%   inputfile      = string specifying "input" textfile (path/name),
%                    containing subject information
%   optimize_metric= string specifying the metric used to select optimal pipelines.
%                    See "run_analysis_wrapper" for metrics of each analysis
%                    model. Pipelines chosen to MAXIMIZE metric of interest
%   mot_gs_control = 2D binary vector, indicating whether we control for
%                    (1) motion artifact and (2) white matter bias, using
%                    spatial priors:
%                                     [0 0] = no artifact control
%                                     [1 0] = control motion artifact
%                                     [0 1] = control white matter bias
%                                     [1 1] = control both
%                    *recommended to correct for at least motion i.e. [1 0]
%   process_out    = binary value specifying whether optimally preprocessed
%                    4D nifti data is output,
%                       0=no, 1=yes (includes 3 fixed pipelines; 1 individually optimized)
%   keepmean       = binary flag, determines if voxel means are re-added to
%                    optimally processed data
%
% OUTPUT:
%
%   (1) Set of Matlab/Octave data, in .mat file named:
%       [outputdirectory,'/optimization_results/matfiles/optimization_summary.mat']
%
%       Includes three cell arrays, (1) METRIC_opt  (2) SPM_opt  (3) TEMP_opt
%       containing optimal pipeline metrics, brain maps and timecourses,
%       respectively. Each cell entry has 4 sub-fields of "con","fix","ind", "min", "max"
%       for the 4 different pipelines optimization methods. Cell entries
%       correspond to individual subjects, as listed in 'inputfile'. For
%       example, the kth subject has optimal "fixed" pipeline results of:
%
%         SPM_opt{k}.fix    = the (voxels x 1) activation map
%         TEMP_opt{k}.fix   = associated (time x 1) BOLD timecourse vector
%         METRIC_opt{k}.fix = structure containing scalar performance metrics,
%                             e.g. for an 'LDA' analysis model:
%
%                             METRIC_opt{k}.fix.R  METRIC_opt{k}.fix.P  METRIC_opt{k}.fix.dPR
%
%       the four different optimization pipelines are:
%
%         (con) = conservative pipeline (all standard steps ON)
%         (fix) = highest-ranked fixed pipeline
%         (ind) = individually optimized subject pipelines
%         (max) = optimal fixed pipeline with most extensive preprocessing
%         (min) = optimal fixed pipeline with least preprocessing
%
%       Also Includes "pipeline_sets" structure with the following fields:
%
%         pipeline_sets.(con/fix/min/max): (pipelines x 10) design matrix, specifying
%                                           pipeline steps for each optimization method
%         pipeline_sets.ind          : (subject x pipeline step) matrix of design
%                                      vectors for each indivially optimized subject
%         pipeline_sets.Signif_Fix   : (pipelines x pipelin step) design matrix
%                                      of all optimal fixed pipeline choices
%
%         pipeline_sets.(con/fix/min/max)_index: index number of each optimal fixed pipeline
%         pipeline_sets.ind_index          : vector of indices denoting each subject's optimal pipeline
%
%   (2) if (niiout==1), produces 4D nifti file of three optimized pipeline SPMs,
%       one per subject. Syntax:  [outdir,'/matfiles/spms',prefix,'_optimized_CON_FIX_IND.nii']
%
%       *note that MIN and MAX pipelines are omitted here
%       *for each volume, timepoints correspond to pipelines:  1=(con) 2=(fix) 3=(ind)
%
%   (3) if (process_out==1), produces the preprocessed 4D fMRI data, under
%       the four optimal preprocessing pipelines. Syntax is:
%
%         [outdir,'/matfiles/niftis_',prefix,'/Preprocessed_',prefix,'_CON.nii']
%         [outdir,'/matfiles/niftis_',prefix,'/Preprocessed_',prefix,'_FIX.nii']
%         [outdir,'/matfiles/niftis_',prefix,'/Preprocessed_',prefix,'_IND.nii']
%

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 167 $';
% CODE_DATE    = '$Date: 2014-12-15 17:31:32 -0500 (Mon, 15 Dec 2014) $';
% ------------------------------------------------------------------------%

global NUMBER_OF_CORES CODE_PROPERTY
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if (~exist('OCTAVE_VERSION','builtin') && exist('maxNumCompThreads'))
    maxNumCompThreads(NUMBER_OF_CORES);
end

global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART2.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end
addpath(CODE_PATH);
addpath([CODE_PATH 'NIFTI_tools']);

if isempty(AFNI_PATH) || isempty(FSL_PATH)
    read_settings;
end

read_version;

if nargin<6
    reference_file = '';
end

%%
if ischar(process_out)
    process_out = str2double(process_out);
end
if ischar(mot_gs_control)
    mot_gs_control = str2double(mot_gs_control);
    if mot_gs_control==0
        mot_gs_control = [0,0];
    end
    if mot_gs_control==1
        mot_gs_control = [0,1];
    end
    if mot_gs_control==10
        mot_gs_control = [1,0];
    end
    if mot_gs_control==11
        mot_gs_control = [1,1];
    end
end
if nargin<5
    keepmean = 0;
else
    if ischar(keepmean)
        keepmean = str2num(keepmean);
    end
end

output_notes{1} = CODE_PROPERTY.NII_HEADER;
output_notes{2} = ['optimization metric: ' optimize_metric];

% labels of different pipeline steps, for plotting
el_list{1} = 'Motcor';
el_list{2} = 'Censor';
el_list{3} = 'Retroi';
el_list{4} = 'Tslice';
el_list{5} = 'Smooth';
el_list{6} = 'Detrend';
el_list{7} = 'Motreg';
el_list{8} = 'Taskreg';
el_list{9} = 'GSPC1';
el_list{10}= 'Lowpass';
el_list{11}= 'Phyplus';

%% Read Inputfiles
if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end
if size(InputStruct,2)~=1
    MULTI_RUN_INPUTFILE = true;
    aligned_suffix = '_aligned';
else
    aligned_suffix = '';
end
Nsubject = numel(InputStruct);

%%
for ksub = 1:Nsubject
    Nrun(ksub) = numel(InputStruct(ksub).run);
end

%% check for multiple pipelines - load first subject
load(strcat(InputStruct(1).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats',InputStruct(ksub).run(1).subjectprefix,'.mat'));
metric_names = fieldnames( METRIC_set{1} );

if(length(METRIC_set)>1) %% if more than one pipeline found, we do optimization...
    
    disp('Now selecting optimal pipelines...');
    
    %% [I] FIRST RUN-THROUGH, LOADING PERFORMANCE METRICS

    % opens the inputfile (includes subject/dataset names that preprocessing is performed on...
    for ksub=1:Nsubject

        %% ----- load STATS only ----- %%
        load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats',InputStruct(ksub).run(1).subjectprefix,'.mat'));

        for i = 1:length(METRIC_set)
            for j = 1:length(metric_names)
                if ~strcmp(metric_names{j},'artifact_prior')
                    METRIC_set{i}.(metric_names{j}) = mean(METRIC_set{i}.(metric_names{j}));
                end
            end

            METRIC_set{i}.artifact_prior.MOT_corr = mean(METRIC_set{i}.artifact_prior.MOT_corr);
            METRIC_set{i}.artifact_prior.WM_zscor = mean(METRIC_set{i}.artifact_prior.WM_zscor);
            METRIC_set{i}.artifact_prior.GS_fract = mean(METRIC_set{i}.artifact_prior.GS_fract);
        end
        % initial check on parameters, done after 1st subject is loaded
        if( ksub==1 )

            % must be single-component analysis, otherwise terminate
            if( ~strcmp( modeltype, 'one_component') )
                warning('Multi-component analysis currently not supported!');
            end
            %
            % get the list of available metrics
            metric_names = fieldnames( METRIC_set{1} );
            % drop out the "artifact priors" entry
            metric_names(strcmp(metric_names,'artifact_prior'))=[];

            N_metric     = length(metric_names);
            %
            % optimization metric must be present on list
            if( sum(strcmp( metric_names, optimize_metric )) == 0 )
                error(['The metric ', optimize_metric, ' is not an output in this dataset. Choose another!']);
            end
            % index of metric used for optimization
            idx_optMet   = find( strcmp( metric_names, optimize_metric ) );
            % initialize cell array of different performance metrics
            allMetrics   = cell(N_metric,1);
        end

        % SUBJECT STATISTICS
        for(ip=1:length( METRIC_set ) ) % iterate through pipelines
            for(im=1:N_metric)         % and through metrics
                allMetrics{im}( ip, ksub ) = METRIC_set{ip}.(metric_names{im});
            end
            % degree of motion-related edge artifact
            spat_mot(ip, ksub) = METRIC_set{ip}.artifact_prior.MOT_corr;
            % amount of white matter signal bias
            spat_wm(ip, ksub) = METRIC_set{ip}.artifact_prior.WM_zscor;

            spat_gsf(ip, ksub) = METRIC_set{ip}.artifact_prior.GS_fract;
        end
    end


    %% [II] Model Optimization, using performance metrics
    %

    % useful info: TR and number of timepoints
    [split_info] = Parse_Split_Info(InputStruct(1).run(1).split_info_file);
    TR_MSEC = split_info.TR_MSEC;
    Ntime = get_numvols([InputStruct(1).run(1).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(1).run(1).Output_nifti_file_prefix '_baseproc.nii']);

    %=================================================================================================
    %=================================================================================================

    % CONSERVATIVE PIPELINE:

    % smoothing (special consideration): closest to 6mm FWHM
    smolist  = unique(pipeset(:,5));
    smodif   = abs(smolist-6.00);
    [v i]= min( smodif ); v_smo = smolist(i);
    % detrend (special consideration): closest to heuristic recommended (for split length of Ntime/2)
    k_ord    = 1+floor( (TR_MSEC/1000) * (Ntime/2) ./ 150 );
    detlist  = unique(pipeset(:,6));
    detdif   = abs(detlist-k_ord);
    [v i]= min( detdif ); v_det = detlist(i);

    % ideal conservative pipeline [MC=1 CENS=1 RET=1 TSLC=1, (smooth~6mm) (detrend~heuristic) MPR=1, TASK=0,GS=0,PHY+=0,LP=0] 
    true_conserv = [1 1 1 1, v_smo v_det, 1 0 0 0 0];
    % find closest possible match
    [v ibase] = max( sum( pipeset == repmat(true_conserv,[size(pipeset,1) 1]), 2 ) );
    % notify if not perfect fit
    if( v<size(pipeset,2) ) disp('WARNING: full conservative pipeline not on list. We can only approximate it!'); end

    %=================================================================================================
    %=================================================================================================

    % FIXED PIPELINES:

    % ranking (higher rank=better), take the median across subjects
    % NB: optimization performed on chosen metric
    MED_rank     = median( tiedrank( allMetrics{idx_optMet} ) , 2);
    % test for significant ordering
    [prob cdist] =  friedman_test( allMetrics{idx_optMet} ) ;
    % recording rank significance
    output_notes{3} = ['fixed ranking significance at p=',num2str(prob)];

    % index the set of significant pipelines
    optfix0   = find( MED_rank >= (max(MED_rank)-cdist) );
    % get the pipeline names
    optpipes0 = pipeset( optfix0, : );
    % index optimal fix pipelines with both MC and PHY+ included
    optkeeps = intersect( find(optpipes0(:,1) ==1), find(optpipes0(:,end) ==1) );

    % if MC and PHYCAA+ are significantly detrimental, warn & re-evaluate:
    if( isempty(optkeeps) )
        %%
        disp( 'UNEXPECTED PIPELINE TRENDS: pipelines with both MC and PHY+ not selected' );

        if( sum(pipeset(:,1))>0 && sum(pipeset(:,end))>0 ) % if both exist

            disp('Readjusting FIX choice...');

            % get highest ranked MC&PHY+ pipeline
            gomax = max( MED_rank( pipeset(:,1) & pipeset(:,end) ) );
            % index the set of significant pipelines
            optfix0   = find( MED_rank >= gomax );
            % get the pipeline names
            optpipes0 = pipeset( optfix0, : );
            % index optimal fix pipelines with both MC and PHY+ included
            optkeeps = intersect( find(optpipes0(:,1) ==1), find(optpipes0(:,end) ==1) );
            % added flag to notify of unexpected trends
            output_notes{4} = ['WARNING: adjusted FIX threshold to include MC & PHY+'];

        elseif( sum(pipeset(:,1))>0 || sum(pipeset(:,end))>0 )  % just MC/PHY+

            % get highest ranked MC or PHY+ pipeline
            gomax = max( MED_rank( pipeset(:,1) | pipeset(:,end) ) );
            % index the set of significant pipelines
            optfix0   = find( MED_rank >= gomax );
            % get the pipeline names
            optpipes0 = pipeset( optfix0, : );
            % index optimal fix pipelines with both MC and PHY+ included
            optkeeps = union( find(optpipes0(:,1) ==1), find(optpipes0(:,end) ==1) );
            % added flag to notify of unexpected trends
            output_notes{4} = ['WARNING: FIX options may not include MC/PHY+'];        
        else
            % just restate the index list
            optkeeps = 1:size(optpipes0,1);
            % added flag to notify of unexpected trends
            output_notes{4} = ['WARNING: FIX options do not include MC/PHY+'];        
        end
    end

    % sub-select pipeline indices with MC/PHY+
    optfix   = optfix0( optkeeps );
    optpipes = optpipes0( optkeeps,: );

    % selection 1: optimally ranked pipeline
    [v ix] = max( MED_rank(optfix) ); irank= optfix(ix);
    % count net #pipeline steps applied in each case
    stuffdone= sum(optpipes,2);
    % selection 2: minimal processing in optimal set
    [v ix] = min(stuffdone);  iless= optfix(ix);
    % selection 3: most extensive processing in optimal set
    [v ix] = max(stuffdone);  imore= optfix(ix);

    %=================================================================================================
    %=================================================================================================

    % INDIVIDUAL PIPELINES:

    if( (mot_gs_control(1)>0) && (mot_gs_control(2)==0) ) %% 1. Only control for motion artifact

        disp(['IND, correcting for motion artifact']);
        % record to output
        output_notes{5} = ['IND: correct for motion only'];

        % best set with no steps held fixed
        [vind,iind0] = max( allMetrics{idx_optMet},[],1 );
        % best set with MC fixed on
        [vind,iindM] = max( allMetrics{idx_optMet} - 1000*repmat( (1-pipeset(:,1)), [1 Nsubject] ),[],1 );

        % spatial correlation with edge artifact for each subject / pair of pipelines
        for(is=1:Nsubject) cormot(is,:) = [ spat_mot( iind0(is),is ) spat_mot( iindM(is),is ) ]; end
        % test: cormot value > sigthr indicates sinificant increase in motion artifact
        sigthr = mean(cormot(:,2)) + 2.3*std(cormot(:,2));
        % now dictate optimal subject pipelines accordingly

        for( is=1:Nsubject )
            if( cormot(is,1) > sigthr ) iind(is,1) = iindM(is);
            else                        iind(is,1) = iind0(is);
            end
        end
        if( exist('OCTAVE_VERSION','builtin') )
            %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE

            disp(['figures not plotted in Octave']);
        elseif (usejava('desktop') && usejava('jvm'))

            % Check whether desktop and JAVA are available in MATLAB
            %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE
            figure;
            subplot(1,2,1);
            fence_plot( cormot, 0, 0.001 );
            hold on; plot( [0.75 2.25], sigthr*[1 1], '--k' ); xlim([0.75 2.25]);
            pairflag = (cormot(:,1)==cormot(:,2));
            plot( ones(sum(pairflag),1) , cormot(pairflag,1), 'ok', 2*ones(sum(pairflag),1) , cormot(pairflag,2), 'ok','markerfacecolor','g', 'markersize',4);
            text(2.05, sigthr*1.05, 'p=.01');
            title('Spatial motion artifact');
            set(gca,'Xtick',1:2,'XTickLabel',{'no MC', 'with MC'});
            ylabel('edge artifact correlation');
        end
    elseif( (mot_gs_control(1)==0) && (mot_gs_control(2)>0) )    %% 2. Only control for global signal bias

        disp(['IND, correcting for GS bias']);
        % record to output
        output_notes{5} = ['IND: correct for GS bias only'];

        % best set with no steps held fixed
        [vind iind0] = max( allMetrics{idx_optMet},[],1 );
        % best set with GSPC1 fixed on
        [vind iindG] = max( allMetrics{idx_optMet} - 1000*repmat( (1-pipeset(:,end-1)), [1 Nsubject] ) ,[],1);

        % spatial correlation with edge artifact for each subject / pair of pipelines
        for(is=1:Nsubject)
            cor_gs(is,:) = abs([ spat_gsf( iind0(is),is ) spat_gsf( iindG(is),is ) ] - 0.5) + 0.5;
        end
        % test: cormot value > sigthr indicates sinificant increase in motion artifact
        sigthr = mean(cor_gs(:,2)) + 2.3*std(cor_gs(:,2));
        % now dictate optimal subject pipelines accordingly
        for( is=1:Nsubject )
            if( cor_gs(is,1) > sigthr ) iind(is,1) = iindG(is);
            else                        iind(is,1) = iind0(is);
            end
        end

        if( exist('OCTAVE_VERSION','builtin') )
            %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE

            disp(['figures not plotted in Octave']);
        elseif (usejava('desktop') && usejava('jvm'))

            % Check whether desktop and JAVA are available in MATLAB
            %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE
            figure;
            subplot(1,2,2);
            fence_plot( cor_gs, 0, 0.001 );
            hold on; plot( [0.75 2.25], sigthr*[1 1], '--k' ); xlim([0.75 2.25]);
            pairflag = (cor_gs(:,1)==cor_gs(:,2));
            plot( ones(sum(pairflag),1) , cor_gs(pairflag,1), 'ok', 2*ones(sum(pairflag),1) , cor_gs(pairflag,2), 'ok','markerfacecolor','g', 'markersize',4);
            text(2.05, sigthr*1.05, 'p=.01');
            title('Global signal bias');
            set(gca,'Xtick',1:2,'XTickLabel',{'no GSPC1', 'with GSPC1'});
            ylabel('fraction of voxels with same sign');
        end

    elseif( (mot_gs_control(1)>0) && (mot_gs_control(2)>0) )  %% 3. Control for both motion artifact and global signal bias

        disp(['IND, correcting for motion artifact AND GS bias']);
        % record to output
        output_notes{5} = ['IND: correct for motion and GS bias'];

        % best set with no steps held fixed
        [vind iind0] = max( allMetrics{idx_optMet},[],1 );
        % best set with MC fixed on
        [vind iindM] = max( allMetrics{idx_optMet} - 1000*repmat( (1-pipeset(:,1)), [1 Nsubject] ),[],1 );

        % spatial correlation with edge artifact for each subject / pair of pipelines
        for(is=1:Nsubject) cormot(is,:) = [ spat_mot( iind0(is),is ) spat_mot( iindM(is),is ) ]; end
        % test: cormot value > sigthr indicates sinificant increase in motion artifact
        sigthr = mean(cormot(:,2)) + 2.3*std(cormot(:,2)); saa = sigthr;

        % null out any significant motion data
        Mset_star = allMetrics{idx_optMet} - double(spat_mot > sigthr) *1000;

        %=================

        % best set with no steps held fixed
        [vind iind0] = max( Mset_star,[],1 );
        % best set with GSPC1 fixed on
        [vind iindG] = max( Mset_star - 1000*repmat( (1-pipeset(:,end-1)), [1 Nsubject] ),[],1 );

        % spatial correlation with edge artifact for each subject / pair of pipelines
        for(is=1:Nsubject) cor_gs(is,:) = abs([ spat_gsf( iind0(is),is ) spat_gsf( iindG(is),is ) ] - 0.5) + 0.5; end
        % test: cormot value > sigthr indicates sinificant increase in motion artifact
        sigthr = mean(cor_gs(:,2)) + 2.3*std(cor_gs(:,2)); sbb = sigthr;
        % now dictate optimal subject pipelines accordingly
        for( is=1:Nsubject )
            if( cor_gs(is,1) > sigthr ) iind(is,1) = iindG(is);
            else                        iind(is,1) = iind0(is);
            end
        end

        % FIGURE PLOTTING ONLY FOR MATLAB
        if( exist('OCTAVE_VERSION','builtin') )
            %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE

            disp(['figures not plotted in Octave']);
        elseif (usejava('desktop') && usejava('jvm'))
            figure;

            subplot(1,2,1);
            fence_plot( cormot, 0, 0.001 );
            hold on; plot( [0.75 2.25], saa*[1 1], '--k' ); xlim([0.75 2.25]);
            pairflag = (cormot(:,1)==cormot(:,2));
            plot( ones(sum(pairflag),1) , cormot(pairflag,1), 'ok',2*ones(sum(pairflag),1) , cormot(pairflag,2), 'ok','markerfacecolor','g', 'markersize',4);
            text(2.05, saa*1.05, 'p=.01');
            title('Spatial motion artifact (stage1)');
            set(gca,'Xtick',1:2,'XTickLabel',{'no MC', 'with MC'});
            ylabel('edge artifact correlation');
            subplot(1,2,2);
            fence_plot( cor_gs, 0, 0.001 );
            hold on; plot( [0.75 2.25], sbb*[1 1], '--k' ); xlim([0.75 2.25]);
            pairflag = (cor_gs(:,1)==cor_gs(:,2));
            plot( ones(sum(pairflag),1) , cor_gs(pairflag,1), 'ok',2*ones(sum(pairflag),1) , cor_gs(pairflag,2), 'ok','markerfacecolor','g', 'markersize',4);
            text(2.05, sbb*1.05, 'p=.01');
            title('Global signal bias (stage2)');
            set(gca,'Xtick',1:2,'XTickLabel',{'no GSPC1', 'with GSPC1'});
            ylabel('fraction of GS voxels with same sign');
        end

    else %% 4. No control of artifact

        disp('IND, no apriori artifact control - BE CAREFUL WHEN INTERPRETING!');
        % record to output
        output_notes{5} = ['IND: no apriori artifact control'];

        % best set with no steps held fixed
        [vind iind] = max( allMetrics{idx_optMet},[],1 );
    end


    % SUMMARIZE METRICS
    if Nsubject==1
        irank = iind(1);
        iless = iind(1);
        imore = iind(1);
    end

    for(is=1:Nsubject)

        for(im=1:N_metric)         % and through metrics
            %
            METRIC_opt.con.(metric_names{im})(is,1) = allMetrics{im}(ibase   ,is);
            METRIC_opt.fix.(metric_names{im})(is,1) = allMetrics{im}(irank   ,is);
            METRIC_opt.min.(metric_names{im})(is,1) = allMetrics{im}(iless   ,is);
            METRIC_opt.max.(metric_names{im})(is,1) = allMetrics{im}(imore   ,is);
            METRIC_opt.ind.(metric_names{im})(is,1) = allMetrics{im}(iind(is),is);
        end
    end

    % FIGURE PLOTTING ONLY FOR MATLAB
    if( exist('OCTAVE_VERSION','builtin') )
        %%% FIGURE PLOTTING TURNED OFF FOR OCTAVE

        disp(['figures not plotted in Octave']);
    elseif (usejava('desktop') && usejava('jvm'))
        %plotting performance metrics
        figure;
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.10 0.15 0.70 0.60]);

        for(im=1:N_metric) % iterate through metrics
            %
            subplot(1,N_metric,im); fence_plot(  [METRIC_opt.con.(metric_names{im}) ...
                METRIC_opt.fix.(metric_names{im}) ...
                METRIC_opt.ind.(metric_names{im})], 0, 0.001 );
            %
            set(gca,'Xtick',1:3,'XTickLabel',{'CON', 'FIX', 'IND'});
            ylabel(metric_names{im}); xlim([0.5 3.5]);
        end
    end


    % list of significant fixed pipelines
    pipeline_sets.Signif_Fix = optpipes0;
    % record 3 optimal fixed choices
    pipeline_sets.fix = pipeset(irank  ,:);
    pipeline_sets.min = pipeset(iless  ,:);
    pipeline_sets.max = pipeset(imore  ,:);
    pipeline_sets.con = pipeset(ibase  ,:);

    % list of individual optimized subject pipeline choices
    for(is=1:Nsubject)
        pipeline_sets.ind(is,:) = pipeset(iind(is),:);
    end
    % get the specific pipeline indices
    pipeline_sets.fix_index = irank;
    pipeline_sets.min_index = iless;
    pipeline_sets.max_index = imore;
    pipeline_sets.ind_index = iind;
    pipeline_sets.con_index = ibase;

    pipeline_sets.optimize_metric = optimize_metric;

    %  save pipeline names / representative characters
    pipeline_sets.pipechars = pipechars;
    pipeline_sets.pipenames = pipenames;

    %% [IV] NOW GENERATE PREPROCESSED DATA

    if( process_out > 0 ) %% only if option turned on

        optType   = {'CON','FIX','IND','MIN','MAX'};

        for ksub = 1:Nsubject

            % Read the subject directory
            MM = load_untouch_nii( InputStruct(ksub).run(1).subjectmask );
            load([InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res0_params/params' InputStruct(ksub).run(1).subjectprefix '.mat']);        
            load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res1_spms/spms', InputStruct(ksub).run(1).subjectprefix ,'.mat')); %IMAGE_set
            load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res2_temp/temp', InputStruct(ksub).run(1).subjectprefix,'.mat')); %TEMP_set

            % load optimal data (spms+timeseries) into cell arrays
            SPM_opt{ksub}.con = IMAGE_set{ibase     };
            SPM_opt{ksub}.fix = IMAGE_set{irank     };
            SPM_opt{ksub}.ind = IMAGE_set{iind(ksub)};

            %---
            TEMP_opt{ksub}.con = TEMP_set{ibase      };
            TEMP_opt{ksub}.fix = TEMP_set{irank      };
            TEMP_opt{ksub}.ind = TEMP_set{iind(ksub) };

            SPMs = [SPM_opt{ksub}.con SPM_opt{ksub}.fix SPM_opt{ksub}.ind];
            nii = MM;
            nii.img = zeros([size(MM.img) size(SPMs,2)]);
            Z = zeros(size(MM.img));
            for nvol = 1:size(SPMs,2)
                Z(MM.img~=0) = SPMs(:,nvol);
                nii.img(:,:,:,nvol) = Z;
            end
            nii.hdr.dime.datatype = 16;
            nii.hdr.dime.dim([1 5]) = [4 nvol];
            nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
            save_untouch_nii(nii,strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/optimization_results/spms/rSPM',InputStruct(ksub).run(1).subjectprefix,'_CON_FIX_IND.nii'));

            %% select additional preprocessing choices (Step2) for current pipeline
            pipe_temp = [pipeline_sets.con;pipeline_sets.fix;pipeline_sets.ind(ksub,:);pipeline_sets.min; pipeline_sets.max;];

            %% Save optimum parameters for each subject in the pipeline
            SV.pipeline_sets = pipeline_sets;
            SV.pipeline_sets.ind = pipeline_sets.ind(ksub,:);
            SV.pipeline_sets.ind = pipeline_sets.ind(ksub,:);
            SV.pipeline_sets.ind_for_the_group = pipeline_sets.ind;
            %SV.METRIC_opt = METRIC_opt;
            for metric_counter = 1:length(metric_names)
                if ~strcmp(metric_names{metric_counter},'artifact_prior')
                    SV.METRIC_opt.con.(metric_names{metric_counter})=METRIC_opt.con.(metric_names{metric_counter})(ksub,1);
                    SV.METRIC_opt.fix.(metric_names{metric_counter})=METRIC_opt.fix.(metric_names{metric_counter})(ksub,1);
                    SV.METRIC_opt.ind.(metric_names{metric_counter})=METRIC_opt.ind.(metric_names{metric_counter})(ksub,1);
                    SV.METRIC_opt.min.(metric_names{metric_counter})=METRIC_opt.min.(metric_names{metric_counter})(ksub,1);
                    SV.METRIC_opt.max.(metric_names{metric_counter})=METRIC_opt.max.(metric_names{metric_counter})(ksub,1);
                end
            end

            SV.SPM_opt  = SPM_opt{ksub};
            SV.TEMP_opt = TEMP_opt{ksub};
            SV.nii_mask = MM;
            SV.ksub     = ksub;
            save([InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/matfiles/opt_' optimize_metric InputStruct(ksub).run(1).subjectprefix '.mat'],'-struct','SV','-v7');
            save([InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/matfiles/opt_' optimize_metric InputStruct(ksub).run(1).subjectprefix '.mat'],'CODE_PROPERTY','-append');

            for ik = 1:size(pipe_temp,1)-2   % Only con, fix and ind files are generated
                Invol_name = ['m',num2str(pipe_temp(ik,1)),'c',num2str(pipe_temp(ik,2)),'p',num2str(pipe_temp(ik,3)),'t',num2str(pipe_temp(ik,4)),'s',num2str(pipe_temp(ik,5))];
                nomem = sprintf('m%dc%dp%dt%ds%dd%dr%dx%dg%dl%dy%d',pipe_temp(ik,:));
                NR = load(sprintf('%s/intermediate_metrics/regressors/reg%s/m%dc%dp%dt%ds%dd%dr%dx%dg%dl%dy%d.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).subjectprefix,pipe_temp(ik,:)));

                for krun = 1:Nrun(ksub)

                    aligned_suffix_alt = '_aligned';
                    volname_alt    = strcat( InputStruct(ksub).run(krun).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',Invol_name,aligned_suffix_alt,'.nii' );
                    if exist(volname_alt,'file')
                        aligned_suffix= '_aligned';
                    else
                        aligned_suffix = '';
                    end

                    volname    = strcat( InputStruct(ksub).run(krun).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',Invol_name,aligned_suffix,'.nii' );
                    VV         = load_untouch_nii(  volname );
                    % convert nifti volume into matfile
                    volmat  = nifti_to_mat(VV,MM);
                    % store voxelwise means
                    mean_volmat = mean(volmat,2);

                    if iscell(NR.Regressors) % Multi-Run design
                        volmat_temp{1} = volmat;
                        volmat_temp = apply_glm(volmat_temp,NR.Regressors(krun));
                        volmat = volmat_temp{1};
                    else
                        volmat = apply_glm(volmat,NR.Regressors);
                    end

                    if( pipe_temp(ik,end-1) == 1 ) % check for lowpass filtering option
                        %
                        %filters above 0.10 Hz
                        [ volmat ] = quick_lopass( volmat, (TR_MSEC./1000) );
                    end    

                    volmat  = bsxfun(@times,volmat,NN_weight_avg);

                    nii = MM;
                    if keepmean
                        volmat  = bsxfun(@plus,volmat,mean_volmat);
                        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' optType{ik} '-'   nomem   ];
                    else
                        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' optType{ik} '-' nomem   ];
                    end    

                    nii.img = zeros([size(MM.img) size(volmat,2)]);
                    Z = zeros(size(MM.img));
                    for nvol = 1:size(volmat,2)
                        Z(MM.img~=0) = volmat(:,nvol);
                        nii.img(:,:,:,nvol) = Z;
                    end
                    nii.hdr.dime.datatype = 16;
                    nii.hdr.dime.dim([1 5]) = [4 nvol];
                    output_nii_path        = strcat(InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed');
                    output_nii_file        = strcat(output_nii_path,'/Proc_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',optType{ik},'.nii');
                    mkdir_r(output_nii_path)
                    %normalized_output_file =  strcat(InputStruct(ksub).run(1).Output_nifti_file_path,'/niftis_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'/',dataoutfix,'_opt_',optimize_metric,'_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',optType{ik},'_sNorm.nii');
                    save_untouch_nii(nii,output_nii_file);
                end
            end
        end
        save(strcat(InputStruct(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat'),'SPM_opt','TEMP_opt', 'METRIC_opt','pipeline_sets','CODE_PROPERTY','-v7');
    end

else %% If only 1 pipeline being tested, this becomes the default output
    

    pipeline_sets.pipe1 = pipeset; %
    %  save pipeline names / representative characters
    pipeline_sets.pipechars = pipechars;
    pipeline_sets.pipenames = pipenames;

    %% [IV] NOW GENERATE PREPROCESSED DATA

    if( process_out > 0 ) %% only if option turned on

        for ksub = 1:Nsubject

            % Read the subject directory
            MM = load_untouch_nii( InputStruct(ksub).run(1).subjectmask );
            load([InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res0_params/params' InputStruct(ksub).run(1).subjectprefix '.mat']);        
            load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res1_spms/spms', InputStruct(ksub).run(1).subjectprefix ,'.mat')); %IMAGE_set
            load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res2_temp/temp', InputStruct(ksub).run(1).subjectprefix,'.mat')); %TEMP_set
            load(strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats',InputStruct(ksub).run(1).subjectprefix,'.mat')); %METRIC_set        

            if( (length(IMAGE_set)>1) || (length(TEMP_set)>1) || (size(pipeset,1)>1) ) error('number of images/timeseries/pipelines does not match metrics'); end
            
            % load optimal data (spms+timeseries) into cell arrays
            SPM_opt{ksub}.pipe1     = IMAGE_set{1};
            TEMP_opt{ksub}.pipe1    = TEMP_set{1};
            METRIC_opt{ksub}.pipe1 = METRIC_set{1};


            SPMs = [SPM_opt{ksub}.pipe1];
            nii = MM;
            nii.img = zeros([size(MM.img) size(SPMs,2)]);
            Z = zeros(size(MM.img));
            for nvol = 1:size(SPMs,2)
                Z(MM.img~=0) = SPMs(:,nvol);
                nii.img(:,:,:,nvol) = Z;
            end
            nii.hdr.dime.datatype = 16;
            nii.hdr.dime.dim([1 5]) = [4 nvol];
            nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
            save_untouch_nii(nii,strcat(InputStruct(ksub).run(1).Output_nifti_file_path, '/optimization_results/spms/rSPM',InputStruct(ksub).run(1).subjectprefix,'_PIPE1.nii'));

            %% select additional preprocessing choices (Step2) for current pipeline
            pipe_temp = [pipeline_sets.pipe1];
            optType{1} = 'PIPE1';

            SV.SPM_opt  = SPM_opt{ksub};
            SV.TEMP_opt = TEMP_opt{ksub};
            SV.METRIC_opt = METRIC_opt;
            
            SV.nii_mask = MM;
            SV.ksub     = ksub;
            save([InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/matfiles/opt_' optimize_metric InputStruct(ksub).run(1).subjectprefix '.mat'],'-struct','SV','-v7');
            save([InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/matfiles/opt_' optimize_metric InputStruct(ksub).run(1).subjectprefix '.mat'],'CODE_PROPERTY','-append');

            for ik = 1   % Only pipe1 files are generated
                Invol_name = ['m',num2str(pipe_temp(ik,1)),'c',num2str(pipe_temp(ik,2)),'p',num2str(pipe_temp(ik,3)),'t',num2str(pipe_temp(ik,4)),'s',num2str(pipe_temp(ik,5))];
                nomem = sprintf('m%dc%dp%dt%ds%dd%dr%dx%dg%dl%dy%d',pipe_temp(ik,:));
                NR = load(sprintf('%s/intermediate_metrics/regressors/reg%s/m%dc%dp%dt%ds%dd%dr%dx%dg%dl%dy%d.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).subjectprefix,pipe_temp(ik,:)));

                for krun = 1:Nrun(ksub)

                    aligned_suffix_alt = '_aligned';
                    volname_alt    = strcat( InputStruct(ksub).run(krun).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',Invol_name,aligned_suffix_alt,'.nii' );
                    if exist(volname_alt,'file')
                        aligned_suffix= '_aligned';
                    else
                        aligned_suffix = '';
                    end

                    volname    = strcat( InputStruct(ksub).run(krun).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',Invol_name,aligned_suffix,'.nii' );
                    VV         = load_untouch_nii(  volname );
                    % convert nifti volume into matfile
                    volmat  = nifti_to_mat(VV,MM);
                    % store voxelwise means
                    mean_volmat = mean(volmat,2);

                    if iscell(NR.Regressors) % Multi-Run design
                        volmat_temp{1} = volmat;
                        volmat_temp = apply_glm(volmat_temp,NR.Regressors(krun));
                        volmat = volmat_temp{1};
                    else
                        volmat = apply_glm(volmat,NR.Regressors);
                    end

                    if( pipe_temp(ik,end-1) == 1 ) % check for lowpass filtering option
                        %
                        %filters above 0.10 Hz
                        [ volmat ] = quick_lopass( volmat, (TR_MSEC./1000) );
                    end    

                    volmat  = bsxfun(@times,volmat,NN_weight_avg);

                    nii = MM;
                    if keepmean
                        volmat  = bsxfun(@plus,volmat,mean_volmat);
                        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER  ' PREPROCESSING: KEEPMEAN-' optType{ik} '-'   nomem   ];
                    else
                        nii.hdr.hist.descrip = [CODE_PROPERTY.NII_HEADER ' PREPROCESSING: REMOVEMEAN-' optType{ik} '-' nomem   ];
                    end    

                    nii.img = zeros([size(MM.img) size(volmat,2)]);
                    Z = zeros(size(MM.img));
                    for nvol = 1:size(volmat,2)
                        Z(MM.img~=0) = volmat(:,nvol);
                        nii.img(:,:,:,nvol) = Z;
                    end
                    nii.hdr.dime.datatype = 16;
                    nii.hdr.dime.dim([1 5]) = [4 nvol];
                    output_nii_path        = strcat(InputStruct(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed');
                    output_nii_file        = strcat(output_nii_path,'/Proc_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',optType{ik},'.nii');
                    mkdir_r(output_nii_path)
                    %normalized_output_file =  strcat(InputStruct(ksub).run(1).Output_nifti_file_path,'/niftis_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'/',dataoutfix,'_opt_',optimize_metric,'_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_',optType{ik},'_sNorm.nii');
                    save_untouch_nii(nii,output_nii_file);
                end
            end
        end
        save(strcat(InputStruct(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat'),'SPM_opt','TEMP_opt', 'METRIC_opt','pipeline_sets','CODE_PROPERTY','-v7');
    else
        disp('You ran 1 pipeline, and chose not to create outputs? Why are you running this step again?');    
    end
           
end
    
function x = get_numvols(file)

hdr = load_nii_hdr(file);
x = hdr.dime.dim(5);