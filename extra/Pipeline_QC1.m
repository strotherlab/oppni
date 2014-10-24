function Pipeline_QC1( inputfile, optimize_matfile_name, prefix_result, save_prefix )
%
%==========================================================================
% PIPELINE_QC1: Quality Control plotting tool. Shows impact of running
% pipelines on subject data, including SPMs (brain maps) and performance metrics.
%==========================================================================
%
% SYNTAX:    
%
%   Pipeline_QC1( inputfile, save_prefix )
%            
% INPUT:
%
%   inputfile             = string specifying "input" textfile (path/name), 
%                           containing subject information
%   optimize_matfile_name = string specifying FULL path+name of summary optimization 
%                           file created during pipeline step-3
%                           NOTE: this is stictly optional! If you haven't run optimization, 
%                           you can still do QC. Just set it as empty:
%                               optimize_matfile_name = []
%   prefix_result         = string specifying the exact prefix used when
%                           runnin the pipeline 
%   save_prefix           = string specifying prefix, for saved plots 
%
% OUTPUT:  set of results saved to "QC1_results" local directory
%
% 1. series of plots with prefix <save_prefix>:
%
%    ['FIG1_motion_statistics.png']: plots estimates of head motion, #motion spikes, and task correlation
%    ['FIG2_spm_artifact.png']: plots measure of artifact in SPM, including motion, global signal and white matter
%    ['FIG3_pipeline_similarity_by_dataset.png']: plots distribution of performance metrics and SPM similarity for pipelines
%    ['FIG4_effects_of_pipeline_steps.png']: shows relative effects of individual pipeline steps on SPMs and metrics
%
% 2. summary results of QC analysis, saved in .mat files:
%
%   [save_prefix,'_output_qc1.mat']: 
%                             output_qc1.motion_stats = head motion statistics (shown in fig.1)
%                             output_qc1.spm_artifact = measures of artifact in SPM (shown in fig.2)
%                             output_qc1.metrics      = performance metrics and average SPM correlation for each pipeline
%                             output_qc1.pipelines    = contains list of
%                             all tested pipelines
%  
% 3. set of diagnostic figures showing optimally preprocessed datasets.
%    Only available if "optimize_matfile_name" input is provided!
% 
%    ['QC1_results/diagnostics/STD/',prefix]: plots for standard (convservative) pipeline
%    ['QC1_results/diagnostics/FIX/',prefix]: plots for optimal fixed pipeline
%    ['QC1_results/diagnostics/IND/',prefix]: plots for individually optimized pipelines 
%

% add paths
% addpath MatFiles;
% addpath MatFiles/NIFTI_tools;
% mkdir QC1_results;
%%%%%%%%%%%%%%%%%%%%%%%%%addpath MatFiles/NIFTI_tools;
[pathstr] = which('Pipeline_QC1.m');
pathstr = fileparts(pathstr);
pathstr = fileparts(pathstr);
addpath([pathstr '/scripts_matlab']);
addpath([pathstr '/scripts_matlab/NIFTI_tools']);

%%%%%%%%%%%%%%%%%%%%%%%%%
QC1_results = ['QC1_results_',save_prefix];
mkdir (QC1_results);

% load pipeline optimization data --> only if this is specified
if( isempty(optimize_matfile_name) )
    optflag = 0;
else
    load(strcat(optimize_matfile_name));
    optflag = 1;
end

% opens the inputfile (includes subject/dataset names that preprocessing is performed on...
fid = fopen(inputfile);
% read in first line
tline = fgetl(fid);
ksub=0;
while ischar(tline) % for every input line in textfile...

    % index #subjects(/datasets)
    ksub=ksub+1,
       
    %% Read in subject paths

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
    
%%  % ------------------------------------------------------------------------------------------------
    % Loading information on motion effects in data
    
    % load spiking information 
    qc_instring  = strcat(outdir,'/diagnostic/',prefix,'_smo_QC_output.mat'    ); % without motion correction
    load(qc_instring); o1 = output;
    qc_instring  = strcat(outdir,'/diagnostic/',prefix,'_mc+smo_QC_output.mat'    ); % with motion correction
    load(qc_instring); o2 = output;
    % record spike information
    Ntime = length( o1.censor_mot );
    spikeset(ksub,:) = [sum(o1.censor_mot==0) sum(o1.censor_vol==0) sum(o2.censor_vol==0) sum(o1.censor_volmot==0) sum(o2.censor_volmot==0)];
    
    % load motion parameter estimates (MPEs)
    mpe_instring  = strcat(outdir,'/mpe/',prefix,'_mpe'    );
    X = load(mpe_instring);
    motor(ksub,:) = mean( abs(X) ); % average displacement values

%% produces full diagnostic outputs on optimally preprocessed pipelines
if( optflag>0 )

    if(ksub==1) 
        mkdir (QC1_results,'/diagnostics/STD');
        mkdir (QC1_results,'/diagnostics/FIX');
        mkdir (QC1_results,'/diagnostics/IND');
    end
    
    MM=(strcat(outdir,'/masks/',prefix,'_mask.nii'));  
    % STD pipeline
    VV=(strcat(outdir,'/matfiles/niftis_',prefix,'/',prefix_result,'_opt_',pipeline_sets.optimize_metric,'_',prefix,'_STD.nii'));  
    
    diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_results,'/diagnostics/STD/',prefix] );%Saman
    % FIX pipeline
    VV=(strcat(outdir,'/matfiles/niftis_',prefix,'/',prefix_result,'_opt_',pipeline_sets.optimize_metric,'_',prefix,'_FIX.nii'));%Saman
    diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_results,'/diagnostics/FIX/',prefix] );
    % IND pipeline
    VV=(strcat(outdir,'/matfiles/niftis_',prefix,'/',prefix_result,'_opt_',pipeline_sets.optimize_metric,'_',prefix,'_IND.nii'));%Saman
    diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_results,'/diagnostics/IND/',prefix] );
end    
%%

    %% Acquire task-design information, to measure task-coupled motion
    %
    if    ( isfield(split_info,'type') && strcmp( split_info.type, 'block' ) )

        % build task-design vector from input information
        % this is a vector of signed values; -1=task condition1, 1=task condition2
        design = zeros(Ntime,1);
        design( [split_info.idx_cond1_sp1(:); split_info.idx_cond1_sp2(:)] ) = -1;
        design( [split_info.idx_cond2_sp1(:); split_info.idx_cond2_sp2(:)] ) =  1; 
        % smooth with standard HRF function
        HRFdesign = design_to_hrf( design, (split_info.TR_MSEC/1000), [5.0 15.0] );        
        
    elseif( isfield(split_info,'type') && strcmp( split_info.type, 'event' ) )
        
        % building design vector: we want to subsample at 100 ms (faster than human RT)
        Nsubs  = max([1 round(split_info.TR_MSEC/100)]); % (#samples per TR) catch in case TR<100ms
        design = zeros( Ntime*Nsubs, 1 );     % initialize design matrix
        % index allocates onsets to appropriate design-points
        didx = unique(round( split_info.onsetlist./(split_info.TR_MSEC/Nsubs) ));
        % catch + adjust overruns, setvalue=1 on design vector
        didx(didx==0)          = 1;
        didx(didx>Ntime*Nsubs) = Ntime*Nsubs;
        design( didx )         = 1;
        % convolve with HRF (must convert into seconds!)
        HRFdesign = design_to_hrf( design, (split_info.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
        % now, subsample back to get HRF at actual fMRI sampling rate
        HRFdesign = HRFdesign( round(Nsubs/2): Nsubs : end );
    else
        % otherwise leave empty
        HRFdesign = [];
    end
    % record task-MPE correlation
    if( ~isempty(HRFdesign) ) tcorrset(ksub,:) = abs( corr(HRFdesign, o1.eigvect_mot(:,1)) );
    else                      tcorrset(ksub,:) = 0;
    end
           
%%  % ------------------------------------------------------------------------------------------------
    % Loading information on SPM (brain map) and performance metrics

    % load subject SPM data
    mat_instring  = strcat(outdir,'/matfiles/results1_spms_',prefix,'.mat'    );
    load(mat_instring);

    % compute cross-correlation between pipeline SPMs
    eigo=zeros( size(IMAGE_set{1},1), length(IMAGE_set) );
    for(i=1:length(IMAGE_set))
       eigo(:,i)=zscore(IMAGE_set{i});
    end
    cormat = eigo'*eigo;
    cormat = cormat ./ ( size(IMAGE_set{1},1) - 1 );
    % record median + (95%CI) on SPM correlations
    corpct(ksub,:) = prctile(cormat(cormat<1),[2.5 50 97.5]);
    % for recording later
    corref(:,ksub) = (sum(cormat, 2)-1)./size(cormat, 2);

    if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
        % temporary indexing:
        ixtemp = [pipeline_sets.std_index ... 
                  pipeline_sets.fix_index ...
                  pipeline_sets.ind_index(ksub)];
        %
        cor_opt(ksub,:) = median( cormat(ixtemp,:), 2 );
    end
    
    mat_instring  = strcat(outdir,'/matfiles/results3_stats_',prefix,'.mat'    );
    load(mat_instring);
    
    % get the list of available performance metrics
    if( ksub==1 )
        %
        % get the list of available metrics
        metric_names = fieldnames( METRIC_set{1} );
        % drop out the "artifact priors" entry
        metric_names(strcmp(metric_names,'artifact_prior'))=[];
        %
        N_metric     = length(metric_names);
        % initialize cell array of different performance metrics
        allMetrics   = cell(N_metric,1);
    end
    %%%
    for(ip=1:length( METRIC_set ) ) % iterate through pipelines
        for(im=1:N_metric)         % and through metrics
           allMetrics{im}( ip, ksub ) = METRIC_set{ip}.(metric_names{im});
        end
        % degree of motion-related edge artifact
        spat_mot( ip, ksub ) = METRIC_set{ip}.artifact_prior.MOT_corr;
        % amount of white matter signal bias
        spat_wmz( ip, ksub ) = METRIC_set{ip}.artifact_prior.WM_zscor;
        % amount of global signal bias
        spat_gsf( ip, ksub ) = METRIC_set{ip}.artifact_prior.GS_fract;
        
        if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
            % temporary indexing:
            ixtemp = [pipeline_sets.std_index ... 
                      pipeline_sets.fix_index ...
                      pipeline_sets.ind_index(ksub)];
            % get artifact measures for optimal pipelines
            spat_mot_optSFI(ksub,:) = [METRIC_set{ixtemp(1)}.artifact_prior.MOT_corr,...
                                       METRIC_set{ixtemp(2)}.artifact_prior.MOT_corr,...
                                       METRIC_set{ixtemp(3)}.artifact_prior.MOT_corr];
            spat_wmz_optSFI(ksub,:) = [METRIC_set{ixtemp(1)}.artifact_prior.WM_zscor,...
                                       METRIC_set{ixtemp(2)}.artifact_prior.WM_zscor,...
                                       METRIC_set{ixtemp(3)}.artifact_prior.WM_zscor];
            spat_gsf_optSFI(ksub,:) = [METRIC_set{ixtemp(1)}.artifact_prior.GS_fract,...
                                       METRIC_set{ixtemp(2)}.artifact_prior.GS_fract,...
                                       METRIC_set{ixtemp(3)}.artifact_prior.GS_fract];
        end

    end

    % compute euclidean distance matrix between pipelines, in combined metrics
    for(im=1:N_metric) metricmat(:,im) = allMetrics{im}(:,ksub); end
    metricmat = zscore(metricmat);
    orcmat=zeros( size(metricmat,1));
    for(i=1:size(metricmat,1)-1)
        for(j=i+1:size(metricmat,1))
            orcmat(i,j)= sqrt( sum((metricmat(i,:)-metricmat(j,:)).^2) );
        end
    end
    
    % measure relative influence of each pipeline step on (spm correlation; metric similarity)
    % we record average matrix value, from all entries comparing step=ON vs step=OFF
    for(i=1:size(pipeset,2))
        
        % take design-column for a given step
        xx1  = pipeset(:,i);
        % check to see if step is varied at all     
        if( length(unique(xx1)) > 1 )
            %
            % if step is varied, measure effect on metrics + spm correlation
            xx1  = (xx1 - min(xx1))./(max(xx1)-min(xx1));
            xx2  = 1-xx1;
            xx12 = xx1*xx2' + xx2*xx1';
            %
            CDIF(i,ksub) = mean( 1-cormat(xx12>0) );
            MDIF(i,ksub) = mean( orcmat(xx12>0) );
        else
            % otherwise just set =0 
            CDIF(i,ksub) = 0;
            MDIF(i,ksub) = 0;            
        end
    end

    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    % SUB-SECTION: PLOTTING OUTLIER TESTING
    
    % ------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------
    
    % then get next line (subject)...
    tline = fgetl(fid);
end

%% Plot#1: motion artifact measures

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1 0.1 0.80 0.80]);
subplot(3,1,1);
    imagesc( motor', prctile( motor(:), [1 99]) ); colorbar;
    title('Average rigid-body displacement');
    %xlabel('dataset #');
    set(gca, 'YTick', 1:6 );
    set(gca,'YTickLabel',{'ROLL (deg)';'PITCH (deg)';'YAW (deg)';'IS (z-axis; mm)';'RL (x-axis; mm)';'AP (y-axis; mm)'});
subplot(3,1,2);
    imagesc( spikeset', [0 ceil(0.1*Ntime)] ); colorbar;
    title(['Number of motion spikes per run (p<0.01), for T=',num2str(Ntime),' time points']);
    %xlabel('dataset #');
    set(gca, 'YTick', 1:5 );
    set(gca,'YTickLabel',{'MPE';'(NO MC): fmri';'(MC): fmri'; '(NO MC): both';'(MC): both'});
subplot(3,1,3);
    bar( tcorrset );
    title('Correlation between MPE and task design');
    xlabel('dataset #');
    ylabel('correlation');
    ylim([0 0.5]);
    xlim([0.5 length(tcorrset)+0.5]);

savestring = [QC1_results,'/',save_prefix,'_FIG1_motion_statistics.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);

%% Plot#2: measures of artifact on pipeline SPMs

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1 0.1 0.80 0.80]);

% choose to plot distribution medians
plot_med=1;

% first plot
subplot(1,4,1);
if( length(unique( pipeset(:,1) ))>1 )
    plot( 0,0,'r','linewidth',2); plot(0,0,'b','linewidth',2);
    prto1 = prctile(spat_mot( pipeset(:,1)==0,:),[5 50 95] )';
    prto2 = prctile(spat_mot( pipeset(:,1)==1,:),[5 50 95] )';
    my_error_plot( [prto1(:,2) prto2(:,2)], [prto1(:,1) prto2(:,1)], [prto1(:,3) prto2(:,3)],'rb',plot_med );
    plot( [0.5 size(spat_mot,1)+0.5], [0.0 0.0], ':k' );
    title('SPM correlation with motion artifact');
    xlabel('dataset #'); ylim([0 1]); legend('MOTCOR ON','MOTCOR OFF');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)
else
    prto1 = prctile(spat_mot,[5 50 95] )';
    my_error_plot( [prto1(:,2) ], [prto1(:,1) ], [prto1(:,3) ],'r',plot_med );
    plot( [0.5 size(spat_mot,1)+0.5], [0.0 0.0], ':k' );
    title('SPM correlation with motion artifact');
    xlabel('dataset #'); ylim([0 1]);
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
end
if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
    % plot for STD/FIX/IND
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,1), 'xk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,2), 'vk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,3), 'ok', 'markersize', 5, 'color', 'k', 'linewidth',2, 'markerfacecolor', 'k' );
end

% second plot
subplot(1,4,2);
if( length(unique( pipeset(:,7) ))>1 )
    plot( 0,0,'r','linewidth',2); plot(0,0,'b','linewidth',2);
    prto1 = prctile(spat_mot( pipeset(:,7)==0,:),[5 50 95] )';
    prto2 = prctile(spat_mot( pipeset(:,7)==1,:),[5 50 95] )';
    my_error_plot( [prto1(:,2) prto2(:,2)], [prto1(:,1) prto2(:,1)], [prto1(:,3) prto2(:,3)],'rb',plot_med );
    plot( [0.5 size(spat_mot,1)+0.5], [0.0 0.0], ':k' );
    title('SPM correlation with motion artifact');
    xlabel('dataset #'); ylim([0 1]); legend('MOTREG ON','MOTREG OFF');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)
else
    prto1 = prctile(spat_mot,[5 50 95] )';
    my_error_plot( [prto1(:,2) ], [prto1(:,1) ], [prto1(:,3) ],'b',plot_med );
    plot( [0.5 size(spat_mot,1)+0.5], [0.0 0.0], ':k' );
    title('SPM correlation with motion artifact');
    xlabel('dataset #'); ylim([0 1]);
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
end
if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
    % plot for STD/FIX/IND
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,1), 'xk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,2), 'vk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_mot_optSFI,1))', spat_mot_optSFI(:,3), 'ok', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
end

% third plot
subplot(1,4,3);
if( length(unique( pipeset(:,9) ))>1 )
    plot( 0,0,'r','linewidth',2); plot(0,0,'b','linewidth',2);
    prto1 = prctile(spat_wmz( pipeset(:,9)==0,:),[5 50 95] )';
    prto2 = prctile(spat_wmz( pipeset(:,9)==1,:),[5 50 95] )';
    my_error_plot( [prto1(:,2) prto2(:,2)], [prto1(:,1) prto2(:,1)], [prto1(:,3) prto2(:,3)],'rb',plot_med );
    plot( [0.5 size(spat_wmz,1)+0.5], [0.0 0.0], ':k' );
    title('Average white matter Z-score');
    xlabel('dataset #'); ylim([-3 3]); legend('GSPC1 ON','GSPC1 OFF');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)
else
    prto1 = prctile(spat_wmz,[5 50 95] )';
    my_error_plot( [prto1(:,2) ], [prto1(:,1) ], [prto1(:,3) ],'b',plot_med );
    title('Average white matter Z-score');
    xlabel('dataset #'); ylim([-3 3]);
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
end
if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
    % plot for STD/FIX/IND
    plot( (1:size(spat_wmz_optSFI,1))', spat_wmz_optSFI(:,1), 'xk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_wmz_optSFI,1))', spat_wmz_optSFI(:,2), 'vk', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_wmz_optSFI,1))', spat_wmz_optSFI(:,3), 'ok', 'markersize', 5, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
end

% fourth plot
subplot(1,4,4);
if( length(unique( pipeset(:,9) ))>1 )
    plot( 0,0,'r','linewidth',2); plot(0,0,'b','linewidth',2);
    prto1 = prctile(spat_gsf( pipeset(:,9)==0,:),[5 50 95] )';
    prto2 = prctile(spat_gsf( pipeset(:,9)==1,:),[5 50 95] )';
    my_error_plot( [prto1(:,2) prto2(:,2)], [prto1(:,1) prto2(:,1)], [prto1(:,3) prto2(:,3)],'rb',plot_med );
    plot( [0.5 size(spat_gsf,1)+0.5], [0.5 0.5], ':k' );
    title('Fraction of SPM voxels >0');
    xlabel('dataset #'); ylim([0 1]); legend('GSPC1 ON','GSPC1 OFF');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)
else
    prto1 = prctile(spat_gsf,[5 50 95] )';
    my_error_plot( [prto1(:,2) ], [prto1(:,1) ], [prto1(:,3) ],'b',plot_med );
    title('Fraction of SPM voxels >0');
    xlabel('dataset #'); ylim([0 1]);
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
end
if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
    % plot for STD/FIX/IND
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,1), 'xk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,2), 'vk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,3), 'ok', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
end

savestring = [QC1_results,'/',save_prefix,'_FIG2_spm_artifact.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);

%% Plot#3: summary metric and SPM similarity measures

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1 0.1 0.80 0.80]);

% first plot
subplot(1,N_metric+1, 1 );
my_error_plot( corpct(:,2), corpct(:,1), corpct(:,3),'b',1 );

if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
    
    % plot for STD/FIX/IND
    plot( (1:size(cor_opt,1))', cor_opt(:,1), 'xk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(cor_opt,1))', cor_opt(:,2), 'vk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(cor_opt,1))', cor_opt(:,3), 'ok', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );

    plot( [0.5 size(cor_opt,1)+0.5], [0.0 0.0], ':k' );
end
title('correlation of pipeline SPMs');
xlabel('dataset #');
%%%%
set(gca,'FontSize',11)
set(findall(gcf,'type','text'),'FontSize',11)

% the rest of plots
for(n=1:N_metric)
    
    subplot(1,N_metric+1,n+1);  
    prto = prctile(allMetrics{n},[5 50 95] )';
    my_error_plot( prto(:,2), prto(:,1), prto(:,3),'r',1 );

    if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
        
        for(i=1:size(allMetrics{n},2))
            % temporary indexing:
            ixtemp = [pipeline_sets.std_index ... 
                      pipeline_sets.fix_index ...
                      pipeline_sets.ind_index(i)];
            prt_opt(i,:) = allMetrics{n}(ixtemp,i);
        end        
        
        % plot for STD/FIX/IND
        plot( (1:size(prt_opt,1))', prt_opt(:,1), 'xk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
        plot( (1:size(prt_opt,1))', prt_opt(:,2), 'vk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
        plot( (1:size(prt_opt,1))', prt_opt(:,3), 'ok', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    end
    
    title(['metric:  ', metric_names{n}]);
    xlabel('dataset #');
    %%%%
    set(gca,'FontSize',11);
    set(findall(gcf,'type','text'),'FontSize',11);
    % for prediction, reproducibility -- rescale plot limits to [0,1]
    if( strcmp(metric_names{n},'R') || strcmp(metric_names{n},'P'))  ylim([0.0 1.0]);    end
end

savestring = [QC1_results,'/',save_prefix,'_FIG3_pipeline_similarity_by_dataset.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);

%% Plot#4: influence of different pipeline steps on results

    figure;
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.1 0.1 0.80 0.80]);
    % list of pipeline steps
    el_list = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK',    'GSPC1', 'PHYPLUS'};
    % number of steps with varied parameters
    totest   = std(pipeset) >0;
    
    % measure effects on spatial patterns
    CDIF(~isfinite(CDIF))=0;
	RNK      = mean( tiedrank(CDIF(totest>0,:)) , 2 );
   [prob sigdiff] =  friedman_test( CDIF(totest>0,:) );
    RNK2     = double(totest); RNK2(totest>0)=RNK;

	subplot(1,2,1);barh( RNK2 );
    title('Pipeline steps, ranked by effect on SPM correlation');
    set(gca, 'YTick', 1:10 );
    set(gca,'YTickLabel',el_list);
    hold on; 
    plot( max(RNK).*[1 1], [0 11], ':r', (max(RNK)-sigdiff).*[1 1], [0 11], ':r' );
    xlim([0 sum(totest)+0.5]); ylim([0 11]); xlabel('average rank');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
    
    % measure effects on performance metrics (combined)
    MDIF(~isfinite(MDIF))=0;
	RNK      = mean( tiedrank(MDIF(totest>0,:)) , 2 );
   [prob sigdiff] =  friedman_test( MDIF(totest>0,:) );
    RNK2     = double(totest); RNK2(totest>0)=RNK;

	subplot(1,2,2); barh( RNK2 );
    title('Pipeline steps, ranked by effect on perform. metrics (combined)');
    set(gca, 'YTick', 1:10 );
    set(gca,'YTickLabel',el_list);
    hold on; 
    plot( max(RNK).*[1 1], [0 11], ':r', (max(RNK)-sigdiff).*[1 1], [0 11], ':r' );
    xlim([0 sum(totest)+0.5]); ylim([0 11]); xlabel('average rank');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    

savestring = [QC1_results,'/',save_prefix,'_FIG4_effects_of_pipeline_steps.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);
    
fclose(fid);

%% Recording outputs to save:
%
% motion parameter results
output_qc1.motion_stats.avg_disp    = motor;    %displacement
output_qc1.motion_stats.num_spikes  = spikeset; %spikes
output_qc1.motion_stats.task_corr   = tcorrset; %task-coupling of motion
% spatial artifact measures
output_qc1.spm_artifact.motion_corr = spat_mot;
output_qc1.spm_artifact.pos_fraction= spat_gsf;
output_qc1.spm_artifact.wm_zscore   = spat_wmz;
% other metrics
output_qc1.metrics.perform_metrics  = allMetrics;
output_qc1.metrics.avg_spm_corr     = corref;
% pipeline listings
output_pc1.pipelines.pipechars      = pipechars;
output_pc1.pipelines.pipenames      = pipenames;
output_pc1.pipelines.pipeset        = pipeset;
% save results --> must be in Matlab format, so no Octave checking required
savestring = [QC1_results,'/',save_prefix,'_output_qc1.mat'];
save(savestring,'output_qc1');

%%
function my_error_plot( X, L, H, colos, plot_med )

[samp dim] = size(X);
% offset for each entry
shift = linspace(0, 0.1, dim);
shift = shift - mean(shift);

hold on;

for(dd=1:dim)
    
    % errorbars
    for(j=1:size(L,1))
       plot( j*[1 1]+shift(dd), [H(j,dd)  L(j,dd)], 'color', colos(dd), 'linewidth',2 );
    end
    
    % dotplot -- optional
    if( plot_med > 0 )
    plot( (1:samp)+shift(dd), X(:,dd), 'o', 'markersize', 3, 'color',  colos(dd), 'linewidth',1, 'markerfacecolor',colos(dd) );
    end
end

xlim([0.5 samp+0.5]);
