function Pipeline_QC1( inputfile )
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
%    ['QC1_results/diagnostics/CON/',prefix]: plots for standard (convservative) pipeline
%    ['QC1_results/diagnostics/FIX/',prefix]: plots for optimal fixed pipeline
%    ['QC1_results/diagnostics/IND/',prefix]: plots for individually optimized pipelines 
%

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

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

if(~exist(inputfile,'file'))
    error('cannot locate input file - make sure path & name is correct. If it still fails and you are using OSX, replace all "\" with "/".');
end

[InputStruct, MULTI_RUN_INPUTFILE] = Read_Input_File(inputfile);

% --> Detect Optimize Results <--
summaryname = [InputStruct(1).run(1).Output_nifti_file_path, '/optimization_results/matfiles/optimization_summary.mat'];
if( isempty( ls(summaryname) ) )
      optflag=0;
else  optflag=1;
      load(summaryname);
end
% --> create outputs folder
QC1_folder= [InputStruct(1).run(1).Output_nifti_file_path, '/QC1_results'];
mkdir (QC1_folder);

% loading sample of optimization metrics
load([InputStruct(1).run(1).Output_nifti_file_path,'/intermediate_metrics/res3_stats/stats_',InputStruct(1).run(1).Output_nifti_file_prefix,'.mat']);
% get the list of available metrics
metric_names = fieldnames( METRIC_set{1} );
% drop out the "artifact priors" entry
metric_names(strcmp(metric_names,'artifact_prior'))=[];
%
N_metric     = length(metric_names);
% initialize cell array of different performance metrics
allMetrics   = cell(N_metric,1);

spikeset=[];
motor   =[];

for(is=1:numel(InputStruct))

    % load split_info data for further analysis
    Subject_OutputDirIntermed = [InputStruct(is).run(1).Output_nifti_file_path '/intermediate_metrics'];
    subjectprefix = InputStruct(is).run(1).subjectprefix;
    load([Subject_OutputDirIntermed '/res0_params/params' subjectprefix '.mat']);

    outdir = InputStruct(is).run(1).Output_nifti_file_path;
    prefix = InputStruct(is).run(1).Output_nifti_file_prefix;
        
    %% Loading information on motion effects in data
    for( ir=1:numel(InputStruct(is).run) )

        % load spiking information 
        qc_instring  = strcat(outdir,'/intermediate_processed/diagnostic/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_smo_QC_output.mat'    ); % without motion correction
        load(qc_instring); o1 = output;
        qc_instring  = strcat(outdir,'/intermediate_processed/diagnostic/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_mc+smo_QC_output.mat'    ); % with motion correction
        load(qc_instring); o2 = output;
        % record spike information
        spikeset = [spikeset; [sum(o1.censor_mot==0) sum(o1.censor_vol==0) sum(o2.censor_vol==0) sum(o1.censor_volmot==0) sum(o2.censor_volmot==0)]];

        % load motion parameter estimates (MPEs)
        mpe_instring  = strcat(outdir,'/intermediate_processed/mpe/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_mpe'    );
        X = load(mpe_instring);
        motor = [motor; mean( abs(X) )]; % average displacement values
    
        %% reproduces full diagnostic outputs, on spatially normed optimally preprocessed pipelines
        if( optflag>0 )
            
            if(is==1 && ir==1) 
                mkdir_r ([QC1_folder,'/diagnostics_CON']);
                mkdir_r ([QC1_folder,'/diagnostics_FIX']);
                mkdir_r ([QC1_folder,'/diagnostics_IND']);
            end
            %
            MM=(strcat(outdir,'/intermediate_processed/masks/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_mask.nii'));  
            
            if( ~exist( [QC1_folder,'/diagnostics_CON/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_QC_output.mat'], 'file' ) )
            % CON pipeline
            VV=(strcat(outdir,'/optimization_results/processed/Proc_',InputStruct(is).run(ir).Output_nifti_file_prefix,'_CON.nii'));  
            diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_folder,'/diagnostics_CON/',InputStruct(is).run(ir).Output_nifti_file_prefix] );%Saman
            end
            if( ~exist( [QC1_folder,'/diagnostics_FIX/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_QC_output.mat'], 'file' ) )
            % FIX pipeline
            VV=(strcat(outdir,'/optimization_results/processed/Proc_',InputStruct(is).run(ir).Output_nifti_file_prefix,'_FIX.nii'));  
            diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_folder,'/diagnostics_FIX/',InputStruct(is).run(ir).Output_nifti_file_prefix] );
            end
            if( ~exist( [QC1_folder,'/diagnostics_IND/',InputStruct(is).run(ir).Output_nifti_file_prefix,'_QC_output.mat'], 'file' ) )
            % IND pipeline
            VV=(strcat(outdir,'/optimization_results/processed/Proc_',InputStruct(is).run(ir).Output_nifti_file_prefix,'_IND.nii'));  
            diagnostic_fmri_pca( VV, MM, mpe_instring, [QC1_folder,'/diagnostics_IND/',InputStruct(is).run(ir).Output_nifti_file_prefix] );
            end
        end
    
        %% record task-MPE correlations --> is principal motion vector in span of contrast design?
        if( ~isempty(split_info_set{ir}.design_mat) ) 
            [a,b, tcorrset(is,ir), u,v] = canoncorr( split_info_set{ir}.design_mat, o1.eigvect_mot(:,1) );
        else      tcorrset(is,ir)       = 0;
        end
    end
    
    %% Loading information on SPMs (brain maps) and performance metrics

    % load subject SPM data
    mat_instring  = strcat(outdir,'/intermediate_metrics/res1_spms/spms_',prefix,'.mat'    );
    load(mat_instring);

    % compute cross-correlation between pipeline SPM sets
    % currently assumes SPMs are matched between pipelines!
    kq=0;
    corrmat = zeros( length(IMAGE_set) ); %% fix diagonals=0
    for(i=1:length(IMAGE_set)-1)
        A = zscore(IMAGE_set{i});
        for(j=i+1:length(IMAGE_set)) 
            kq = kq+1; %% index counter
            B  = zscore(IMAGE_set{j});
            cc = diag( A'*B )./ (size(IMAGE_set{1},1) - 1);
            %-----
            corrmat(i,j) = mean(cc);
            corrmat(j,i) = mean(cc);
        end
    end
    % record median + (95%CI) on SPM correlations
    corpct(is,:) = prctile( corrmat(abs(corrmat)>eps),[2.5 50 97.5]);
    % for recording later
    corref(:,is) = sum(corrmat,2)./(size(corrmat,2)-1);

    if(optflag>0) %---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
        % temporary indexing:
        ixtemp = [pipeline_sets.con_index ... 
                  pipeline_sets.fix_index ...
                  pipeline_sets.ind_index(is)];
        %
        cor_opt(is,:) = median( corrmat(ixtemp,:), 2 );
    end

    mat_instring  = strcat(outdir,'/intermediate_metrics/res3_stats/stats_',prefix,'.mat'    );
    load(mat_instring);
    %%%
    for(ip=1:length( METRIC_set ) ) % iterate through pipelines
        for(im=1:N_metric)         % and through metrics
           allMetrics{im}( ip, is ) = mean( METRIC_set{ip}.(metric_names{im}) );
        end
        % degree of motion-related edge artifact
        spat_mot( ip, is ) = mean( METRIC_set{ip}.artifact_prior.MOT_corr );
        % amount of white matter signal bias
        spat_wmz( ip, is ) = mean( METRIC_set{ip}.artifact_prior.WM_zscor );
        % amount of global signal bias
        spat_gsf( ip, is ) = mean( METRIC_set{ip}.artifact_prior.GS_fract );

        if(optflag>0) %%%---------------- ONLY TAKE VALUES IF OPTIMIZED RESULTS AVAILABLE
            % temporary indexing:
            ixtemp = [pipeline_sets.con_index ...
                      pipeline_sets.fix_index ...
                      pipeline_sets.ind_index(is)];
            % get artifact measures for optimal pipelines
            spat_mot_optSFI(is,:) = [mean( METRIC_set{ixtemp(1)}.artifact_prior.MOT_corr ),...
                                     mean( METRIC_set{ixtemp(2)}.artifact_prior.MOT_corr ),...
                                     mean( METRIC_set{ixtemp(3)}.artifact_prior.MOT_corr )];
            spat_wmz_optSFI(is,:) = [mean( METRIC_set{ixtemp(1)}.artifact_prior.WM_zscor ),...
                                     mean( METRIC_set{ixtemp(2)}.artifact_prior.WM_zscor ),...
                                     mean( METRIC_set{ixtemp(3)}.artifact_prior.WM_zscor )];
            spat_gsf_optSFI(is,:) = [mean( METRIC_set{ixtemp(1)}.artifact_prior.GS_fract ),...
                                     mean( METRIC_set{ixtemp(2)}.artifact_prior.GS_fract ),...
                                     mean( METRIC_set{ixtemp(3)}.artifact_prior.GS_fract )];
        end

    end

    % compute euclidean distance matrix between pipelines, in combined metrics
    for(im=1:N_metric) metricmat(:,im) = allMetrics{im}(:,is); end
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
            CDIF(i,is) = mean( 1-corrmat(xx12>0) );
            MDIF(i,is) = mean( orcmat(xx12>0) );
        else
            % otherwise just set =0 
            CDIF(i,is) = 0;
            MDIF(i,is) = 0;            
        end
    end    
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
    imagesc( spikeset'); colorbar; %, [0 ceil(0.1*Ntime)] ); colorbar;
    title(['Number of motion spikes per run (p<0.01)']);%, for T=',num2str(Ntime),' time points']);
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

savestring = [QC1_folder,'/FIG1_motion_statistics.png'];
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
    % plot for CON/FIX/IND
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
    % plot for CON/FIX/IND
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
    % plot for CON/FIX/IND
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
    % plot for CON/FIX/IND
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,1), 'xk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,2), 'vk', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
    plot( (1:size(spat_gsf_optSFI,1))', spat_gsf_optSFI(:,3), 'ok', 'markersize', 4, 'color', 'k', 'linewidth',1, 'markerfacecolor', 'k' );
end

savestring = [QC1_folder,'/FIG2_spm_artifact.png'];
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
    
    % plot for CON/FIX/IND
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
            ixtemp = [pipeline_sets.con_index ... 
                      pipeline_sets.fix_index ...
                      pipeline_sets.ind_index(i)];
            prt_opt(i,:) = allMetrics{n}(ixtemp,i);
        end        
        
        % plot for CON/FIX/IND
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

savestring = [QC1_folder,'/FIG3_pipeline_similarity_by_dataset.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);

%% Plot#4: influence of different pipeline steps on results

    figure;
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.1 0.1 0.80 0.80]);
    % list of pipeline steps
    el_list = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK','GSPC1','LOWPASS','PHYPLUS'};
    % number of steps with varied parameters
    totest   = std(pipeset) >0;
    
    % measure effects on spatial patterns
    CDIF(~isfinite(CDIF))=0;
	RNK      = mean( tiedrank(CDIF(totest>0,:)) , 2 );
   [prob sigdiff] =  friedman_test( CDIF(totest>0,:) );
    RNK2     = double(totest); RNK2(totest>0)=RNK;

	subplot(1,2,1);barh( RNK2 );
    title('Pipeline steps, ranked by effect on SPM correlation');
    set(gca, 'YTick', 1:11 );
    set(gca,'YTickLabel',el_list);
    hold on; 
    plot( max(RNK).*[1 1], [0 12], ':r', (max(RNK)-sigdiff).*[1 1], [0 12], ':r' );
    xlim([0 sum(totest)+0.5]); ylim([0 12]); xlabel('average rank');
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
    set(gca, 'YTick', 1:11 );
    set(gca,'YTickLabel',el_list);
    hold on; 
    plot( max(RNK).*[1 1], [0 12], ':r', (max(RNK)-sigdiff).*[1 1], [0 12], ':r' );
    xlim([0 sum(totest)+0.5]); ylim([0 12]); xlabel('average rank');
    %%%%
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    

savestring = [QC1_folder,'/FIG4_effects_of_pipeline_steps.png'];
img = getframe(gcf);
imwrite(img.cdata, [savestring]);
    
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
savestring = [QC1_folder,'/output_qc1.mat'];
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
