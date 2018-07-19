function Pipeline_QC2( inputfile, newmaskname, Npcs )
%
%==========================================================================
% PIPELINE_QC2: for examining the results of optimized subject data.
%==========================================================================
%
% SYNTAX:    
%
%   Pipeline_QC2( inputfile, newmaskname, Npcs )
%            
% INPUT:
%
%   inputfile             = string specifying "input" textfile (path/name), 
%                           containing subject information
%   newmaskname           = string specifying name/path of the new group-level 
%                           consensus mask
%                           ** IF running OPPNI, set newmaskname=[], and it
%                           will automatically detect new masks
%   Npcs                  = # of principal components to produce from subject SPMs
%                           for group-level PCA (rSVD analyses).
%                           All selected PCs are plotted AND saved as output.
%                           Produces Npcs results for CON, FIX and IND pipelines
%                           ** Setting Npcs=[] gives a default Npcs=2
%
% OUTPUT:  set of results saved to "QC2_results" local directory
%
% 1. series of plots with prefix <save_prefix>:
%
%    ['FIG1_optimized_pipeline_steps.png']: shows which pipeline steps were turned on/off for each subject 
%    ['FIG2_optimized_performance_metrics.png']: shows subject performance metrics, for CON vs FIX vs IND pipelines 
%    ['FIG3_spatial_norm_statistics.png']: quality control stats for evaluating spatial normalization 
%    ['FIG4_inter_subject_SPM_similarity.png']:  shows  between-subject SPM correlation and overlap, for CON vs FIX vs IND pipelines 
%    ['FIG5_group_pca_CON_PCs_<n1>-<n2>.png']: results of group PCA, for CON pipeline
%    ['FIG6_group_pca_FIX_PCs_<n1>-<n2>.png']: results of group PCA, for FIX pipeline
%    ['FIG7_group_pca_IND_PCs_<n1>-<n2>.png']: results of group PCA, for IND pipeline
%
% 2. summary results of QC analysis, saved in .mat files:
%
%   [save_prefix,'_output_qc2.mat']: 
%                             output_qc2.METRIC_opt    = structure containing scalar performance metrics;
%                             output_qc2.SPM_opt_norm  = optimized subject SPMs;
%                             output_qc2.TEMP_opt      = the (time x 1) BOLD timecourse vector associated with SPMs;
%                             output_qc2.pipeline_sets = contains information about optimal pipelines
%                             output_qc2.spm_corr      = cross-correlation matrices between subject SPMs
%                             output_qc2.spm_ovl       = pairwise jaccard overlap matrices between subject SPMs
%
%    [save_prefix,'_group_pca.mat']: summary results of group PCA, including...
%                              group_pca.(con/fix/ind).rSPMZs = z-scored PCA eigenimages 
%                              group_pca.(con/fix/ind).rep    = vector of reproducibility values for each PC  
%                              group_pca.(con/fix/ind).scores = matrix of subject PC scores (column vectors) 
%
% 3. eigenimages from group PCA, saved as .nii files:
%
%    [save_prefix,'_CON_group_PCs1-',(Npcs),'.nii']
%    [save_prefix,'_FIX_group_PCs1-',(Npcs),'.nii']
%    [save_prefix,'_IND_group_PCs1-',(Npcs),'.nii']
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

% % add paths
% addpath_oppni MatFiles;
% addpath_oppni MatFiles/NIFTI_tools;
% mkdir QC2_results;
%%%%%%%%%%%%%%%%%%%%%%%%%addpath_oppni MatFiles/NIFTI_tools;

[pathstr] = which('Pipeline_QC2.m');
pathstr = fileparts(pathstr);
pathstr = fileparts(pathstr);
addpath_oppni([pathstr '/scripts_matlab']);
addpath_oppni([pathstr '/scripts_matlab/NIFTI_tools']);

% quick open of input file to get directory
fid = fopen(inputfile);
% read in first line
tline = fgetl(fid);
ksub=0;
while ischar(tline) % for every input line in textfile...

    % index #subjects(/datasets)
    ksub=ksub+1,
       
    if(ksub==1)
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
            
    % --> find optimized results
    summaryname = [outdir, '/optimization_results/matfiles/optimization_summary.mat'];
    load(summaryname);
    % --> create outputs folder
    QC2_folder= [outdir, '/QC2_results'];
    mkdir (QC2_folder);

    % name for group masks
    if( isempty(newmaskname) )
        %%% create directory + new mask files, in first output dir.
        newpath = [outdir, '/GroupMasks'];
        mkdir_r(newpath);
        newmaskname = [newpath '/group'];                
    else
        [apath,aprefix,aext] = fileparts(newmaskname);
        mkdir_r(apath);
        %%% make sure mask-name terminates with .nii string
        newmaskname = [apath,'/',aprefix];
    end
    end
    % then get next line (subject)...
    tline = fgetl(fid);
end
fclose(fid);

if( ksub < 4 )
    disp(['Only ',num2str(ksub),' subjects. Insufficient for group-level brain maps - skipping this step']);
    Npcs = [];
elseif( Npcs > floor(ksub/2) )
    disp(['Too many PCs for sample size. Requested ',num2str(Npcs),', but maximum is ',num2str(floor(ksub/2)),'. Resetting...']);
    Npcs = floor(ksub/2);
end


%% 1: plotting pipeline lists

el_list = {'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH', 'DETREND', 'MOTREG', 'TASK', 'GSPC1', 'LOWPASS', 'PHYPLUS'};

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.10 0.15 0.80 0.60]);
bnd = [0 max( pipeline_sets.ind(:) )];
subplot(1,4,1); imagesc( pipeline_sets.con',bnd );
    set(gca, 'YTick', 1:10 );
    set(gca,'YTickLabel',el_list);
    title('Conservative');
    xlabel('all datasets'); set(gca, 'XTick', []);
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
subplot(1,4,2); imagesc( pipeline_sets.fix',bnd );
    title('Fixed');
    xlabel('all datasets'); 
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
subplot(1,2,2); imagesc( pipeline_sets.ind',bnd ); colorbar;
    title('Individually-optimized');
    xlabel('dataset #');
    set(gca, 'YTick', []);
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    

% savestring = [QC2_folder,'/FIG1_optimized_pipeline_steps.png'];
% img = getframe(gcf);
% imwrite(img.cdata, [savestring]);
print_to_file(gcf, QC2_folder,'FIG1_optimized_pipeline_steps');

%% 2: optimization metrics
    
% get list of metrics + number thereof
% plot all metrics, per optimization scheme

% get the list of available metrics
metric_names = fieldnames( METRIC_opt.con );        
N_metric     = length(metric_names);

%plotting performance metrics
figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.10 0.15 0.70 0.60]);

for(im=1:N_metric) % iterate through metrics
    %
    subplot(1,N_metric,im); 
    fence_plot(  [METRIC_opt.con.(metric_names{im}) METRIC_opt.fix.(metric_names{im}) METRIC_opt.ind.(metric_names{im})], 0, 0.001 );
    %
    set(gca,'Xtick',1:3,'XTickLabel',{'CON', 'FIX', 'IND'});
    title(['metric:  ',metric_names{im}]); xlim([0.5 3.5]);
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
    % for prediction, reproducibility -- rescale plot limits to [0,1]
    if( strcmp(metric_names{im},'R') || strcmp(metric_names{im},'P'))  ylim([0.0 1.0]);    end
end

% savestring = [QC2_folder,'/FIG2_optimized_performance_metrics.png'];
% img = getframe(gcf);
% imwrite(img.cdata, [savestring]);
print_to_file(gcf, QC2_folder,'FIG2_optimized_performance_metrics');

load([newmaskname,'_spat_norm_qc.mat']);

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.10 0.10 0.70 0.80]);

subplot(2,2,1); plot( volume_stats.MASK_distance, 'o-k', 'linewidth',2,'markerfacecolor','k');
    hold on; plot( [0.5 length(volume_stats.MASK_distance)+0.5], volume_stats.MASK_outlier_thresh.*[1 1], ':r' );
    title('Mask Volumes');
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11) 
	ylabel('overlap distance');
    xlim( [0.5 length(volume_stats.MASK_distance)+0.5] ); ylim([min(volume_stats.MASK_distance)-0.05 max(volume_stats.MASK_distance)+0.05]);
    text(1.0,volume_stats.MASK_outlier_thresh+0.01,'p=.05 outlier threshold','color','r');
subplot(2,2,2); plot( volume_stats.MEAN_distance, 'o-k', 'linewidth',2,'markerfacecolor','k');
    hold on; plot( [0.5 length(volume_stats.MEAN_distance)+0.5], volume_stats.MEAN_outlier_thresh.*[1 1], ':r' );
    title('Mean Volumes');
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
	ylabel('corr. distance');
    xlim( [0.5 length(volume_stats.MEAN_distance)+0.5] ); ylim([min(volume_stats.MEAN_distance)-0.05 max(volume_stats.MEAN_distance)+0.05]);
    text(1.0,volume_stats.MEAN_outlier_thresh+0.01,'p=.05 outlier threshold','color','r');
subplot(2,2,3); plot( volume_stats.NN_distance, 'o-k', 'linewidth',2,'markerfacecolor','k');
    hold on; plot( [0.5 length(volume_stats.NN_distance)+0.5], volume_stats.NN_outlier_thresh.*[1 1], ':r' );
    title('Non-Neuronal Tissue');
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)
    xlabel('dataset #'); ylabel('corr. distance');
    xlim( [0.5 length(volume_stats.NN_distance)+0.5] ); ylim([min(volume_stats.NN_distance)-0.05 max(volume_stats.NN_distance)+0.05]);
    text(1.0,volume_stats.NN_outlier_thresh+0.01,'p=.05 outlier threshold','color','r');
subplot(2,2,4); plot( volume_stats.WM_distance, 'o-k', 'linewidth',2,'markerfacecolor','k');
    hold on; plot( [0.5 length(volume_stats.WM_distance)+0.5], volume_stats.WM_outlier_thresh.*[1 1], ':r' );
    title('White Matter');
    set(gca,'FontSize',11)
    set(findall(gcf,'type','text'),'FontSize',11)    
    xlabel('dataset #'); ylabel('corr. distance');
    xlim( [0.5 length(volume_stats.WM_distance)+0.5] ); ylim([min(volume_stats.WM_distance)-0.05 max(volume_stats.WM_distance)+0.05]);
    text(1.0,volume_stats.WM_outlier_thresh+0.01,'p=.05 outlier threshold','color','r');

% savestring = [QC2_folder,'/FIG3_spatial_norm_statistics.png'];
% img = getframe(gcf);
% imwrite(img.cdata, [savestring]);
print_to_file(gcf, QC2_folder,'FIG3_spatial_norm_statistics');

% 3:
% load mask; tissue stats
% load spms
% plot overlap+corr
% generate rSVD

MM = load_untouch_nii([newmaskname '_consensus_mask.nii']);
mask = double(MM.img);

    % opens the inputfile (includes subject/dataset names that preprocessing is performed on...
    fid = fopen(inputfile);
    % read in first line
    tline = fgetl(fid);
    ksub=0;
    while ischar(tline) % for every input line in textfile...

        ksub=ksub+1, % index #subjects(/datasets)

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
        [split_info] = Parse_Split_Info(fullline);

        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------

        VV=load_untouch_nii(strcat(outdir,'/optimization_results/spms/rSPM_',prefix,'_CON_FIX_IND_sNorm.nii'));  
        
        spm_set(:,:,ksub) = nifti_to_mat( VV,MM );
        [p thr_set(:,:,ksub)] = fdr( spm_set(:,:,ksub), 'z', 0.05, 0 );
        
        % ------------------------------------------------------------------------------------------------
        % ------------------------------------------------------------------------------------------------

        % storing as outputs
        SPM_opt_norm{ksub}.con = spm_set(:,1,ksub);
        SPM_opt_norm{ksub}.fix = spm_set(:,2,ksub);
        SPM_opt_norm{ksub}.ind = spm_set(:,3,ksub);
        
        % then get next line (subject)...
        tline = fgetl(fid);
    end
    fclose(fid);

    spm_set = permute( spm_set,[1 3 2] );
    thr_set = permute( thr_set,[1 3 2] );

    for(j=1:3)
        
        ctmp = corr( spm_set(:,:,j) );
        jtmp = jaccard_ovl( thr_set(:,:,j),0 );
        %
        cc(:,j) = (sum(ctmp,2)-1)./size(ctmp,2);
        jj(:,j) = (sum(jtmp,2)-1)./size(jtmp,2);
        %
        if    (j==1) spm_corr.con = ctmp; spm_ovl.con = jtmp;
        elseif(j==2) spm_corr.fix = ctmp; spm_ovl.fix = jtmp;
        elseif(j==3) spm_corr.ind = ctmp; spm_ovl.ind = jtmp;
        end
    end    
    
    figure; 
    subplot(1,2,1); fence_plot( cc, 0, 0.001 );
        set(gca,'Xtick',1:3,'XTickLabel',{'CON', 'FIX', 'IND'});
        title('between-subject SPM correlation'); 
        xlim([0.5 3.5]); ylim([0 max(cc(:))+0.1]);
        set(gca,'FontSize',11)
        set(findall(gcf,'type','text'),'FontSize',11)    

    subplot(1,2,2); fence_plot( jj, 0, 0.001 );    
        set(gca,'Xtick',1:3,'XTickLabel',{'CON', 'FIX', 'IND'});
        title('between-subject SPM overlap (FDR=.05)'); 
        xlim([0.5 3.5]); ylim([0 max(jj(:))+0.1]);
        set(gca,'FontSize',11)
        set(findall(gcf,'type','text'),'FontSize',11) 
        
% savestring = [QC2_folder,'/FIG4_inter_subject_SPM_similarity.png'];
% img = getframe(gcf);
% imwrite(img.cdata, [savestring]);
print_to_file(gcf, QC2_folder,'FIG4_inter_subject_SPM_similarity');

if(~isempty(Npcs))

%%---------------- CON

    [ out ] = rSVD_splithalf( spm_set(:,:,1), Npcs, 50, [] );
    
    %% -------- making plots -------- %%
    for(KPC=1:ceil(Npcs/2))
        figure;
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.03 0.1 0.90 0.80]);
        quind = round( linspace( 1, size(mask,3), 10) );
        ns_pc = [2*(KPC-1)+1 2*KPC]; %% 2 pcs being plotted in this fig.
        for(k=1:2)
            %
            if( (rem(Npcs,2)>0) && (KPC == ceil(Npcs/2)) && (k==2) )
                %
                % skip step --> we've overrun Npcs (odd number), so second panel = empty
            else
                %
                signos = sign( mean(out.subj_weights(:,ns_pc(k))) );
                tmp=mask;tmp(tmp>0)=out.rSPM_match(:,ns_pc(k));
                slc1=[]; for(w=1:4) slc1=[slc1 tmp(:,:,quind(w+1))']; end 
                slc2=[]; for(w=1:4) slc2=[slc2 tmp(:,:,quind(w+5))']; end 
                slcall = [slc1;slc2];
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11)    
                subplot(2,2,1+(k-1));%('position',[0.04 0.42 0.45 0.45]);    
                imagesc( slcall .* signos, [-3.5 3.5] ); colorbar;
                set(gca, 'YTick', [], 'XTick', [] );
                title(['CON:  z-scored eigenimage, PC#',num2str(ns_pc(k))]);
                subplot(2,2,3+(k-1));%('position',[0.04 0.09 0.45 0.3]);    
                plot( out.subj_weights(:,ns_pc(k)).* signos, 'o-k', 'linewidth',2,'markerfacecolor','k');
                xlabel('dataset #'); ylabel('PC scores');
                title(['CON:  subject loadings, PC#',num2str(ns_pc(k))]);
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11) 
            end
        end
        %savestring = [QC2_folder,'/FIG5_group_pca_CON_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)),'.png'];
        %img = getframe(gcf);
        %imwrite(img.cdata, [savestring]);
        out_fig_name = [ 'FIG5_group_pca_CON_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)) ];
        print_to_file(gcf, QC2_folder,out_fig_name);
    end
    %% -------- making plots -------- %%
    
    % store results
    group_pca.con.rSPMZs = out.rSPM_match;
    group_pca.con.rep    = out.R_match;    
    group_pca.con.scores = out.subj_weights;
    group_pca.con.var    = out.Var_match;
    
%%---------------- FIX

    [ out ] = rSVD_splithalf( spm_set(:,:,2), Npcs, 50, [] );
    
    %% -------- making plots -------- %%
    for(KPC=1:ceil(Npcs/2))
        figure;
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.03 0.1 0.90 0.80]);
        quind = round( linspace( 1, size(mask,3), 10) );
        ns_pc = [2*(KPC-1)+1 2*KPC]; %% 2 pcs being plotted in this fig.
        for(k=1:2)
            %
            if( (rem(Npcs,2)>0) && (KPC == ceil(Npcs/2)) && (k==2) )
                %
                % skip step --> we've overrun Npcs (odd number), so second panel = empty
            else
                %
                signos = sign( mean(out.subj_weights(:,ns_pc(k))) );
                tmp=mask;tmp(tmp>0)=out.rSPM_match(:,ns_pc(k));
                slc1=[]; for(w=1:4) slc1=[slc1 tmp(:,:,quind(w+1))']; end 
                slc2=[]; for(w=1:4) slc2=[slc2 tmp(:,:,quind(w+5))']; end 
                slcall = [slc1;slc2];
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11)    
                subplot(2,2,1+(k-1));%('position',[0.04 0.42 0.45 0.45]);    
                imagesc( slcall .* signos, [-3.5 3.5] ); colorbar;
                set(gca, 'YTick', [], 'XTick', [] );
                title(['FIX:  z-scored eigenimage, PC#',num2str(ns_pc(k))]);
                subplot(2,2,3+(k-1));%('position',[0.04 0.09 0.45 0.3]);    
                plot( out.subj_weights(:,ns_pc(k)).* signos, 'o-k', 'linewidth',2,'markerfacecolor','k');
                xlabel('dataset #'); ylabel('PC scores');
                title(['FIX:  subject loadings, PC#',num2str(ns_pc(k))]);
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11) 
            end
        end
        %         savestring = [QC2_folder,'/FIG6_group_pca_FIX_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)),'.png'];
        %         img = getframe(gcf);
        %         imwrite(img.cdata, [savestring]);
        out_fig_name = [ 'FIG6_group_pca_FIX_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)) ];
        print_to_file(gcf, QC2_folder,out_fig_name);
    end
    %% -------- making plots -------- %%
    
    % store results
    group_pca.fix.rSPMZs = out.rSPM_match;
    group_pca.fix.rep    = out.R_match;    
    group_pca.fix.scores = out.subj_weights;
    group_pca.fix.var    = out.Var_match;    
    
%%---------------- IND

    [ out ] = rSVD_splithalf( spm_set(:,:,3), Npcs, 50, [] );
    
    %% -------- making plots -------- %%
    for(KPC=1:ceil(Npcs/2))
        figure;
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.03 0.1 0.90 0.80]);
        quind = round( linspace( 1, size(mask,3), 10) );
        ns_pc = [2*(KPC-1)+1 2*KPC]; %% 2 pcs being plotted in this fig.
        for(k=1:2)
            %
            if( (rem(Npcs,2)>0) && (KPC == ceil(Npcs/2)) && (k==2) )
                %
                % skip step --> we've overrun Npcs (odd number), so second panel = empty
            else
                %
                signos = sign( mean(out.subj_weights(:,ns_pc(k))) );
                tmp=mask;tmp(tmp>0)=out.rSPM_match(:,ns_pc(k));
                slc1=[]; for(w=1:4) slc1=[slc1 tmp(:,:,quind(w+1))']; end 
                slc2=[]; for(w=1:4) slc2=[slc2 tmp(:,:,quind(w+5))']; end 
                slcall = [slc1;slc2];
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11)    
                subplot(2,2,1+(k-1));%('position',[0.04 0.42 0.45 0.45]);    
                imagesc( slcall .* signos, [-3.5 3.5] ); colorbar;
                set(gca, 'YTick', [], 'XTick', [] );
                title(['IND:  z-scored eigenimage, PC#',num2str(ns_pc(k))]);
                subplot(2,2,3+(k-1));%('position',[0.04 0.09 0.45 0.3]);    
                plot( out.subj_weights(:,ns_pc(k)).* signos, 'o-k', 'linewidth',2,'markerfacecolor','k');
                xlabel('dataset #'); ylabel('PC scores');
                title(['IND:  subject loadings, PC#',num2str(ns_pc(k))]);
                    set(gca,'FontSize',11)
                    set(findall(gcf,'type','text'),'FontSize',11) 
            end
        end
        %savestring = [QC2_folder,'/FIG7_group_pca_IND_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)),'.png'];
        %img = getframe(gcf);
        %imwrite(img.cdata, [savestring]);
        out_fig_name = [ 'FIG7_group_pca_IND_PCs',num2str(ns_pc(1)),'-',num2str(ns_pc(2)) ];
        print_to_file(gcf, QC2_folder, out_fig_name);
        
    end
    
    %% -------- making plots -------- %%
    
    % store results
    group_pca.ind.rSPMZs = out.rSPM_match;
    group_pca.ind.rep    = out.R_match;    
    group_pca.ind.scores = out.subj_weights;
    group_pca.ind.var    = out.Var_match;

end
    
%% Recording outputs to save:

% re-store metrics
output_qc2.METRIC_opt    = METRIC_opt;
output_qc2.SPM_opt_norm  = SPM_opt_norm;
output_qc2.TEMP_opt      = TEMP_opt;
output_qc2.pipeline_sets = pipeline_sets;
% inter-subject similarities
output_qc2.spm_corr      = spm_corr;
output_qc2.spm_ovl       = spm_ovl;
% save results --> must be in Matlab format, so no Octave checking required
savestring = [QC2_folder,'/output_qc2.mat'];
save(savestring,'output_qc2');

if(~isempty(Npcs))
    % save results --> must be in Matlab format, so no Octave checking required
    savestring = [QC2_folder,'/group_pca.mat'];
    save(savestring,'group_pca');

    % also: make into nifti volumes
    %
    % (CON)
    TMPVOL = zeros( [size(mask), Npcs] );
    for(pc=1:Npcs)  
        tmp=mask; 
        tmp(tmp>0)= group_pca.con.rSPMZs(:,pc);  
        TMPVOL(:,:,:,pc) = tmp;
    end
    nii     = VV; % copy nifti struct
    nii.img = TMPVOL; % replace volume
    nii.hdr.dime.dim(5) = Npcs; % adjust for #timepoints
    %
    savestring = [QC2_folder,'/CON_group_PCs1-',num2str(Npcs),'.nii'];
    save_untouch_nii(nii,savestring);

    % (FIX)
    TMPVOL = zeros( [size(mask), Npcs] );
    for(pc=1:Npcs)  
        tmp=mask; 
        tmp(tmp>0)= group_pca.fix.rSPMZs(:,pc);  
        TMPVOL(:,:,:,pc) = tmp;
    end
    nii     = VV; % copy nifti struct
    nii.img = TMPVOL; % replace volume
    nii.hdr.dime.dim(5) = Npcs; % adjust for #timepoints
    %
    savestring = [QC2_folder,'/FIX_group_PCs1-',num2str(Npcs),'.nii'];
    save_untouch_nii(nii,savestring);

    % (IND)
    TMPVOL = zeros( [size(mask), Npcs] );
    for(pc=1:Npcs)  
        tmp=mask; 
        tmp(tmp>0)= group_pca.ind.rSPMZs(:,pc);  
        TMPVOL(:,:,:,pc) = tmp;
    end
    nii     = VV; % copy nifti struct
    nii.img = TMPVOL; % replace volume
    nii.hdr.dime.dim(5) = Npcs; % adjust for #timepoints
    %
    %
    savestring = [QC2_folder,'/IND_group_PCs1-',num2str(Npcs),'.nii'];
    save_untouch_nii(nii,savestring);

end