function out = OPPNI_version_check_v3( InputFile1, InputFile2, out_name, distatis_flag )
%
%
% Syntax:
%             out = OPPNI_version_check( InputFile1, InputFile2, out_name, distatis_flag );
%
% Input:
%             InputFile1, InputFile2: strings, giving path+name of the input text-files
%                                     used to generate OPPNI outputs being compared
%             out_name              : name of output textfile summarizing comparison
%                                     between OPPNI outputs
%             distatis_flag         : (optional) argument specifying if DISTATIS should be run
%                                     to compare optimized SPMs.
%
% Output:     textfile "out_name", which summarizes comparisons between
%             the two OPPNI results
%

if(nargin<4 || isempty(distatis_flag)) distatis_flag=0; end

%% Read Inputfiles

[InputStruct1,MULTI_RUN_INPUTFILE] = Read_Input_File(InputFile1);
[InputStruct2,MULTI_RUN_INPUTFILE] = Read_Input_File(InputFile2);

% first check -- number of lines match?
Nsubject1 = numel(InputStruct1);
Nsubject2 = numel(InputStruct2);

% initialize output file
fin = fopen(out_name,'wt');
fprintf(fin,'Comparing OPPNI version outputs from:\n    Output-1: %s\n    Output-2: %s\n\n',InputFile1,InputFile2);
fprintf(fin,'============== 1. CONSISTENCY CHECKS ==============\n\n');

if( Nsubject1 ~= Nsubject2 ) %% check if number of subjects matches
    %
    disp('Early termination. Number of subjects does not agree between versions');
    fprintf(fin,'Number of data runs: %u (Output-1) %u (Output-2)..................[X] \n\tINCONSISTENT NUMBER OF LINES. ENDING EARLY.\n',Nsubject1,Nsubject2);
    clear Nsubject1 Nsubject2;
else
    fprintf(fin,'Number of data runs: %u (Output-1) %u (Output-2)..................[OK]\n',Nsubject1,Nsubject2);
    Nsubject= Nsubject1;
    clear Nsubject1 Nsubject2;
    
    %%%%%%%%
    CCmat=[];
    
    optcheck1 = exist( strcat(InputStruct1(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat'), 'file' );
    optcheck2 = exist( strcat(InputStruct2(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat'), 'file' );
    
    if( optcheck1 && optcheck2 )
        
        %% PART-1: checking optimization results
        
        % 1.1 loading and consistency checking on data
        
        % loading optimization summary file
        dout1 = load([InputStruct1(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat']);
        dout2 = load([InputStruct2(1).run(1).Output_nifti_file_path,'/optimization_results/matfiles/optimization_summary.mat']);
        
        % consistency checking: SPM dims
        Nspm1 = size( dout1.SPM_opt{1}.con,2 );
        Nspm2 = size( dout2.SPM_opt{1}.con,2 );
        
        if( Nspm1 ~= Nspm2 )
            disp('Early termination. Inconsistent number of contrast SPMs');
            fprintf(fin,'Number of contrast SPMs: %u (Output-1) %u (Output-2)................[X] \n\tINCONSISTENT NUMBER OF SPMS. ENDING EARLY.\n',Nspm1,Nspm2);
            clear Nspm1 Nspm2;
        else
            fprintf(fin,'Number of contrast SPMs: %u (Output-1) %u (Output-2)................[OK]\n',Nspm1,Nspm2);
            Nspm = Nspm1; clear Nspm1 Nspm2;
            
            % consistency checking: available metrics
            metric_names1 = fieldnames( dout1.METRIC_opt.con );
            metric_names2 = fieldnames( dout2.METRIC_opt.con );
            metric_names  = intersect( metric_names1, metric_names2 );
            
            if( isempty(metric_names) )
                disp('Early termination. No consistent metrics between code sets');
                fprintf(fin,'All metrics consistent?...........................................[X] \n\t NO COMMON METRIC NAMES BETWEEN OUTPUTS. ENDING EARLY.\n');
            else
                
                % non-fatal consistency checks
                
                % Number of matching metrics
                diff = setdiff( metric_names1,metric_names2 ); %% in 1 but not 2
                if(~isempty(diff))
                    out.flags.missing_metr_1 =diff;
                else out.flags.missing_metr_1 =[];
                end
                diff = setdiff( metric_names2,metric_names1 ); %% in 2 but not 1
                if(~isempty(diff))
                    out.flags.missing_metr_2 =diff;
                else out.flags.missing_metr_2 =[];
                end
                if( isempty(out.flags.missing_metr_1) &&  isempty(out.flags.missing_metr_2) )
                    fprintf(fin,'All metrics consistent?...........................................[OK]\n');
                else      fprintf(fin,'All metrics consistent?...........................................[!!]\n');
                    if(~isempty(out.flags.missing_metr_1))
                        fprintf(fin,'\tOutput-1 metrics (%s) not found in Output-2, omitted from comparisons.\n',strjoin(out.flags.missing_metr_1(:)',','));
                    end
                    if(~isempty(out.flags.missing_metr_2))
                        fprintf(fin,'\tOutput-2 metrics (%s) not found in Output-1, omitted from comparison.\n',strjoin(out.flags.missing_metr_2(:)',','));
                    end
                end
                clear metric_names1 metric_names2;
                
                % Voxel consistency testing
                Nvox1=zeros(Nsubject,1);
                Nvox2=zeros(Nsubject,1);
                for(i=1:Nsubject)
                    Nvox1(i,1) = size( dout1.SPM_opt{i}.con, 1 );
                    Nvox2(i,1) = size( dout2.SPM_opt{i}.con, 1 );
                end
                if( sum(Nvox1 == Nvox2)<Nsubject )
                    tmpval = abs(Nvox1-Nvox2); tmpval=tmpval(tmpval>0);
                    fprintf(fin,'Same number of brain voxels?......................................[!!]\n');
                    fprintf(fin,'\tMismatched voxels for %u/%u subjects - have the brain masks changed?\n\tMedian [min, max] in runs showing change: %u [%u, %u]\n\tSPM comparisons will not be performed\n',...
                        median(tmpval),min(tmpval),max(tmpval),sum(Nvox1 == Nvox2),Nsubject);
                    out.flags.diffvox=1;
                else
                    fprintf(fin,'Same number of brain voxels?......................................[OK]\n');
                    out.flags.diffvox=0;
                end
                
                % Tested pipeline set consistency
                pipe_names1 = dout1.pipeline_sets.pipenames;
                pipe_names2 = dout2.pipeline_sets.pipenames;
                pipe_names  = intersect( pipe_names1, pipe_names2,'stable' );
                % for later examination of pipes
                pipe_idx12  = zeros( length(pipe_names),2 );
                for(i=1:size(pipe_idx12,1))
                    pipe_idx12(i, :) = [ strmatch(pipe_names{i},pipe_names1) strmatch(pipe_names{i},pipe_names2) ];
                end
                % checking for specific discrepancies
                diff = setdiff( pipe_names1,pipe_names2 ); %% in 1 but not 2
                if(~isempty(diff))
                    out.flags.missing_pipe_1 =diff;
                else out.flags.missing_pipe_1 =[];
                end
                diff = setdiff( pipe_names2,pipe_names1 ); %% in 2 but not 1
                if(~isempty(diff))
                    out.flags.missing_pipe_2 =diff;
                else out.flags.missing_pipe_2 =[];
                end
                if( isempty(out.flags.missing_pipe_1) &&  isempty(out.flags.missing_pipe_2) )
                    fprintf(fin,'Number of pipeline steps consistent?..............................[OK]\n');
                else      fprintf(fin,'Number of pipeline steps consistent?..............................[!!]\n');
                    if(~isempty(diff_struct.checks.missing_metr_1))
                        fprintf(fin,'\tOutput-1 steps (%s) not found in Output-2\n\tFurther testing only examines consistent steps.\n',strjoin(out.flags.missing_pipe_1(:)',','));
                    end
                    if(~isempty(out.flags.missing_metr_2))
                        fprintf(fin,'\tOutput-2 steps (%s) not found in Output-1\n\tFurther testing only examines consistent steps.\n',strjoin(out.flags.missing_pipe_2(:)',','));
                    end
                end
                clear pipe_names1 pipe_names2;
                
                if( ~strcmp(dout1.pipeline_sets.optimize_metric,dout2.pipeline_sets.optimize_metric) )
                    fprintf(fin,'Same optimization metric?.........................................[!!]\n');
                    fprintf(fin,'\tOutput-1 optimizes on (%s) but Output-2 optimizes on (%s).\n',dout1.pipeline_sets.optimize_metric,dout2.pipeline_sets.optimize_metric);
                else
                    fprintf(fin,'Same optimization metric?.........................................[OK]\n');
                end
                
                fprintf(fin,'\n============== 2. OPTIMIZATION OUTPUTS ==============\n\n');
                
                % consistency checking: number of pipeline steps
                
                % diff on pipelines
                diffs_con = (dout1.pipeline_sets.con(pipe_idx12(:,1))-dout2.pipeline_sets.con(pipe_idx12(:,2)) );
                diffs_fix = (dout1.pipeline_sets.fix(pipe_idx12(:,1))-dout2.pipeline_sets.fix(pipe_idx12(:,2)) );
                diffs_ind = (dout1.pipeline_sets.ind(:,pipe_idx12(:,1))-dout2.pipeline_sets.ind(:,pipe_idx12(:,2)) );
                diffs_ppl = [max(abs(diffs_con)) max(abs(diffs_fix)) max(max(abs(diffs_ind)))];
                
                if( max(diffs_ppl) > 1E-6 ) %% any major differences in pipes
                    
                    fprintf(fin,'Comparing optimal pipelines.......................................[!!]\n');
                    fprintf(fin,'\tDifferences in: ');
                    if(max(abs(diffs_con))     >1E-6) fprintf(fin,'CON, '); end
                    if(max(abs(diffs_fix))     >1E-6) fprintf(fin,'FIX, '); end
                    if(max(max(abs(diffs_ind)))>1E-6) olix = find(max(abs(diffs_ind),[],2) > 1E-6 );
                        fprintf(fin,'IND (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                    end
                else
                    fprintf(fin,'Comparing optimal pipelines.......................................[OK]\n');
                end
                
                % diff on metrics
                ppltype={'con','fix','ind'};
                for i = 1:length(metric_names)
                    for j = 1:length(ppltype)
                        val1 = dout1.METRIC_opt.(ppltype{j}).(metric_names{i});
                        val2 = dout2.METRIC_opt.(ppltype{j}).(metric_names{i});
                        diffs_metric(:,i,j) = abs(val1 - val2);
                    end
                end
                
                if( max(diffs_metric(:)) > 1E-6 ) %% any major differences in pipes
                    
                    fprintf(fin,'Comparing optimized metrics.......................................[!!]\n');
                    fprintf(fin,'\tDifferences in: \n');
                    if(max(max(diffs_metric(:,:,1))) >1E-6)
                        olix=find( max(diffs_metric(:,:,1),[],2) > 1E-6);
                        fprintf(fin,'\tCON (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        cato=[];
                        for(i=1:length(metric_names)) 
                            cato = [cato,sprintf( '(%s)=%0.3f [%0.3f, %0.3f], ', ...
                                metric_names{i},median(diffs_metric(olix,i,1)), ...
                                min(diffs_metric(olix,i,1)),max(diffs_metric(olix,i,1)) )]; 
                        end
                        fprintf(fin,'\t\tMedian [min, max] abs. difference, in runs showing changes: %s\n',cato);
                        
                    end
                    if(max(max(diffs_metric(:,:,2))) >1E-6)
                        olix=find( max(diffs_metric(:,:,2),[],2) > 1E-6);
                        fprintf(fin,'\tFIX (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        cato=[];
                        for(i=1:length(metric_names)) cato = [cato,sprintf( '(%s)=%0.3f [%0.3f, %0.3f], ',metric_names{i},median(diffs_metric(olix,i,2)),min(diffs_metric(olix,i,2)),max(diffs_metric(olix,i,2)) )]; end
                        fprintf(fin,'\t\tMedian [min, max] abs. difference, in runs showing changes: %s\n',cato);
                    end
                    if(max(max(diffs_metric(:,:,3))) >1E-6)
                        olix=find( max(diffs_metric(:,:,3),[],2) > 1E-6);
                        fprintf(fin,'\tIND (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        cato=[];
                        for(i=1:length(metric_names)) cato = [cato,sprintf( '(%s)=%0.3f [%0.3f, %0.3f], ',metric_names{i},median(diffs_metric(olix,i,3)),min(diffs_metric(olix,i,3)),max(diffs_metric(olix,i,3)) )]; end
                        fprintf(fin,'\t\tMedian [min, max] abs. difference, in runs showing changes: %s\n',cato);
                    end
                else
                    fprintf(fin,'Comparing optimized metrics.......................................[OK]\n');
                end
                
                % diff on spms
                if( Nvox1 == Nvox2 )
                    
                    for i = 1:Nsubject
                        for j = 1:length(ppltype)
                            
                            val1 = dout1.SPM_opt{i}.(ppltype{j});
                            val2 = dout2.SPM_opt{i}.(ppltype{j});
                            diff_spms(i,:,j) = max( abs(val1 - val2) ); %./abs(val1 + val2) );
                            corr_spms(i,:,j) = diag(corr(val1,val2));
                        end
                    end
                    
                    if( max(diff_spms(:)) > 1E-6 ) %% any major differences in pipes
                        
                        fprintf(fin,'Comparing optimized SPMs..........................................[!!!]\n');
                        fprintf(fin,'\tDifferences in: \n');
                        if(max(max(diff_spms(:,:,1))) >1E-6)
                            olix=find( max(diff_spms(:,:,1),[],2) > 1E-6);
                            fprintf(fin,'\tCON (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                            cato=[];
                            for(i=1:size(diff_spms,2)) cato = [cato,sprintf( '(contrast-%u)=%0.3f [%0.3f, %0.3f], ',i,median(corr_spms(olix,i,1)),min(corr_spms(olix,i,1)),max(corr_spms(olix,i,1)) )]; end
                            fprintf(fin,'\t\tMedian [min, max] SPM corr., in runs showing changes: %s\n',cato);
                        end
                        if(max(max(diff_spms(:,:,2))) >1E-6)
                            olix=find( max(diff_spms(:,:,2),[],2) > 1E-6);
                            fprintf(fin,'\tFIX (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                            cato=[];
                            for(i=1:size(diff_spms,2)) cato = [cato,sprintf( '(contrast-%u)=%0.3f [%0.3f, %0.3f], ',i,median(corr_spms(olix,i,2)),min(corr_spms(olix,i,2)),max(corr_spms(olix,i,2)) )]; end
                            fprintf(fin,'\t\tMedian [min, max] SPM corr., in runs showing changes: %s\n',cato);
                        end
                        if(max(max(diff_spms(:,:,3))) >1E-6)
                            olix=find( max(diff_spms(:,:,3),[],2) > 1E-6);
                            fprintf(fin,'\tIND (%u of %u runs): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                            cato=[];
                            for(i=1:size(diff_spms,2)) cato = [cato,sprintf( '(contrast-%u)=%0.3f [%0.3f, %0.3f], ',i,median(corr_spms(olix,i,3)),min(corr_spms(olix,i,3)),max(corr_spms(olix,i,3)) )]; end
                            fprintf(fin,'\t\tMedian [min, max] SPM corr., in runs showing changes: %s\n',cato);
                        end
                    else
                        fprintf(fin,'Comparing optimized SPMs..........................................[OK]\n');
                    end
                    
                    %%%%%%%%%%% OPTIONAL DISTATIS PREP. %%%%%%%%%%%%%%%%
                    
                    % On spms
                    CCmat=[];
                    if( distatis_flag>0 ) %% any major differences in metrics
                        if( max(diff_spms(:)) > 1E-6 )
                            
                            disp('Setting up DISTATIS...');
                            %
                            CCmat = zeros( 6,6, Nsubject );
                            for i = 1:Nsubject
                                vall12set = zeros( Nvox1(i), 2*Nspm );
                                for j = 1:length(ppltype)
                                    vall12set(:, j ) = dout1.SPM_opt{i}.(ppltype{j})(:);
                                    vall12set(:,j+3) = dout2.SPM_opt{i}.(ppltype{j})(:);
                                end
                                CCmat(:,:,i) = corr( vall12set );
                            end
                        else
                            disp('Skipping DISTATIS - no significant differences in SPMs');
                        end
                    end
                    %%%%%%%%%%% OPTIONAL DISTATIS PREP %%%%%%%%%%%%%%%%
                    
                end
                
                fprintf(fin,'\n============== 3. INTERMEDIATE OUTPUTS (SPMS) ==============\n\n');
                
                % first check if mismatched on metrics
                dout1 = load([InputStruct1(1).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats' InputStruct1(1).run(1).subjectprefix '.mat']);
                dout2 = load([InputStruct2(1).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats' InputStruct2(1).run(1).subjectprefix '.mat']);
                
                % checking to see if pipeline steps being tested are consistent
                if( size(dout1.pipeset,1) ~= size(dout2.pipeset,1) )
                    fprintf(fin,'Same set of tested steps?.........................................[X]\n');
                    fprintf(fin,'\tDifferent total #pipelines: %u (Output-1) and %u (Output-2) - check your input files\n',size(dout1.pipeset,1),size(dout2.pipeset,1));
                elseif( abs( tiedrank(pipe_idx12(:,1))-tiedrank(pipe_idx12(:,2)) ) > 1E-6 )
                    fprintf(fin,'Same set of tested steps?.........................................[X]\n');
                    fprintf(fin,'\tpermuted pipeline ordering, cannot do comparison - check your input files\n');
                elseif( max(abs(dout1.pipeset - dout2.pipeset)) > 1E-6 )
                    fprintf(fin,'Same set of tested steps?.........................................[X]\n');
                    fprintf(fin,'\Same total #pipelines, but not using common set of steps - check your input files\n');
                else
                    fprintf(fin,'Same set of tested steps?.........................................[OK]\n');
                    
                    % -- common pipeline list --
                    pipeset = dout1.pipeset;
                    clear diffx;
                    for( ksub=1:Nsubject )
                        
                        disp(['runs_,',num2str(ksub),'_of_',num2str(Nsubject)]),
                        
                        dout1 = load([InputStruct1(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats' InputStruct1(ksub).run(1).subjectprefix '.mat']);
                        dout2 = load([InputStruct2(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res3_stats/stats' InputStruct2(ksub).run(1).subjectprefix '.mat']);
                        
                        for i = 1:length(dout1.METRIC_set) % pipelines
                            for j = 1:length(metric_names) % metrix
                                diffx(j,i,ksub) = mean(abs( dout1.METRIC_set{i}.(metric_names{j}) - dout2.METRIC_set{i}.(metric_names{j}) ));
                            end
                        end
                    end
                    
                    diffx = permute( max(diffx,[],1), [3 2 1]);
                    
                    if( max(diffx(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x metric)
                        
                        fprintf(fin,'Comparing all pipeline metrics....................................[!!]\n');
                        fprintf(fin,'\tDifferences in: \n');
                        
                        olix=find( max(max(diffx,[],3),[],2) > 1E-6);
                        fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        olix=find( max(max(diffx,[],3),[],1) > 1E-6);
                        fprintf(fin,'\tPipelines (%u of %u)\n',length(olix),size(pipeset,1));
                        %
                        nonzero_diffs = diffx(diffx>1E-6);
                        fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(nonzero_diffs),min(nonzero_diffs),max(nonzero_diffs));
                        
                        fprintf('\t\tQuantiles [0:0.1:1] of abs. difference, in runs showing changes: \n\t\t  %s\n',strcat(num2str(quantile(nonzero_diffs(:), 0:0.1:1))));
                        
                        %============= more detailed determination
                        diffx_ppl = max(max( diffx,[],1 ),[],3 );
                        for(t=1:size(pipeset,2))
                            uniq=unique(pipeset(:,t));
                            if( length(uniq)>1 )
                                clear vll;
                                for(l=1:length(uniq))
                                    vll(l,:) = min(diffx_ppl(pipeset(:,t)==uniq(l)));
                                end
                                effecte(t,:) = (min(vll)==0) & (max(vll)>1E-6);
                            else
                                effecte(t,1) = false;
                            end
                        end
                        cato=[];
                        for(i=1:length(effecte)) if(effecte(i)) cato=[cato pipe_names{i},',']; end; end
                        fprintf(fin,'\tSteps that alter output: %s\n',cato);
                        fprintf(fin,'\ti.e., step OFF sometimes gives identical outputs,\n\t ON always causes difference; or vice-versa)\n');
                    else
                        fprintf(fin,'Comparing all pipeline metrics....................................[OK]\n');
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % diff on spms
                    if( Nvox1 == Nvox2 )
                        
                        clear diffx corrx;
                        for( ksub=1:Nsubject )
                            
                            disp(['runs_,',num2str(ksub),'_of_',num2str(Nsubject)]),
                            
                            dout1 = load([InputStruct1(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res1_spms/spms' InputStruct1(ksub).run(1).subjectprefix '.mat']);
                            dout2 = load([InputStruct2(ksub).run(1).Output_nifti_file_path, '/intermediate_metrics/res1_spms/spms' InputStruct2(ksub).run(1).subjectprefix '.mat']);
                            
                            for i = 1:length(dout1.IMAGE_set) % contr x pipelines x subj
                                diffx(:,i,ksub) = mean(abs(dout1.IMAGE_set{i}-dout2.IMAGE_set{i}));
                            end
                        end
                        
                        diffx = permute( max(diffx,[],1), [3 2 1]);
                        
                        if( max(diffx(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x contr)
                            
                            fprintf(fin,'Comparing all pipeline SPMs.......................................[!!]\n');
                            fprintf(fin,'\tDifferences in: \n');
                            
                            olix=find( max(max(diffx,[],3),[],2) > 1E-6);
                            fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                            olix=find( max(max(diffx,[],3),[],1) > 1E-6);
                            fprintf(fin,'\tPipelines (%u of %u)\n',length(olix),size(pipeset,1));
                            %
                            fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(diffx(diffx>1E-6)),min(diffx(diffx>1E-6)),max(diffx(diffx>1E-6)));
                            
                            %============= more detailed determination
                            diffx_ppl = max(max( diffx,[],1 ),[],3 );
                            for(t=1:size(pipeset,2))
                                uniq=unique(pipeset(:,t));
                                if( length(uniq)>1 )
                                    clear vll;
                                    for(l=1:length(uniq))
                                        vll(l,:) = min(diffx_ppl(pipeset(:,t)==uniq(l)));
                                    end
                                    effecte(t,:) = (min(vll)==0) & (max(vll)>1E-6);
                                else
                                    effecte(t,1) = false;
                                end
                            end
                            cato=[];
                            for(i=1:length(effecte)) if(effecte(i)) cato=[cato pipe_names{i},',']; end; end
                            fprintf(fin,'\tSteps that alter output: %s\n',cato);
                            fprintf(fin,'\ti.e., step OFF sometimes gives identical outputs,\n\t ON always causes difference; or vice-versa)\n');
                        else
                            fprintf(fin,'Comparing all pipeline SPMs.......................................[OK]\n');
                        end
                        
                        
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                end
                
                fprintf(fin,'\n============== 4. INTERMEDIATE OUTPUTS (PROCESSED) ==============\n\n');
                
                clear Nvox1 Nvox2;
                %per subjects: count number of outputted proc-data (v1/v2)
                for( ksub=1:Nsubject )
                    
                    M1=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask);
                    M2=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask);
                    
                    Nvox1(ksub,1) = sum(M1.img(:));
                    Nvox2(ksub,1) = sum(M2.img(:));
                end
                
                % %                 if( sum(Nvox1 == Nvox2)<Nsubject )
                % %                     fprintf(fin,'Same number of brain voxels?......................................[!!]\n');
                % %                     fprintf(fin,'\tMismatched voxels for %u of %u subjects - have the brain masks changed?\n\tProc comparisons will not be performed\n',sum(Nvox1 == Nvox2),Nsubject);
                % %                     out.flags.diffvox=1;
                % %                 else
                % %                     fprintf(fin,'Same number of brain voxels?......................................[OK]\n');
                % %                     out.flags.diffvox=0;
                % %                 end
                
                if( prod(Nvox1==Nvox2) ) %% if all matching
                    
                    for(ksub=1:Nsubject)
                        
                        disp(['--runs_,',num2str(ksub),'_of_',num2str(Nsubject)]),
                        
                        M1=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask); mask1=double(M1.img);
                        M2=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask); mask2=double(M2.img);
                        
                        %%% 1
                        e=dir([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/Proc',InputStruct1(ksub).run(1).subjectprefix,'*']);
                        clear eu1 es1; k1=0;k2=0;
                        for(i=1:numel(e))
                            if(isempty(strfind(e(i).name,'_sNorm')))
                                k1=k1+1; eu1{k1}=e(i).name;
                            else  k2=k2+1; es1{k2}=e(i).name;
                            end
                        end
                        %%% 2
                        e=dir([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/Proc',InputStruct2(ksub).run(1).subjectprefix,'*']);
                        clear eu2 es2; k1=0;k2=0;
                        for(i=1:numel(e))
                            if(isempty(strfind(e(i).name,'_sNorm')))
                                k1=k1+1; eu2{k1}=e(i).name;
                            else  k2=k2+1; es2{k2}=e(i).name;
                            end
                        end
                        
                        if( prod(strcmpi(eu1,eu2)) )
                            %unnormalized files are same ... proceed to compare contents
                            for(i=1:numel(eu1))
                                V1=load_untouch_nii([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',eu1{i}]);
                                V2=load_untouch_nii([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',eu2{i}]);
                                tmpval = abs(V1.img(:)-V2.img(:));
                                voldif_u(ksub,i) = median( tmpval(tmpval>eps) );
                            end
                        else
                            %file mismatch
                            erstr1=[]; for(i=1:numel(eu1)) erstr1=[erstr1,', ',eu1{i}]; end
                            erstr2=[]; for(i=1:numel(eu2)) erstr2=[erstr2,', ',eu2{i}]; end
                            error('cannot match optimally processed, unnormalized datasets for %s.\nFiles for pipeline 1:%s\nFiles for pipeline 2:%s\n\n',InputStruct1(ksub).run(1).subjectprefix,errstr1(2:end),errstr2(2:end));
                        end
                        %
                        if( prod(strcmpi(es1,es2)) )
                            %unnormalized files are same ... proceed to compare contents
                            for(i=1:numel(es1))
                                V1=load_untouch_nii([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',es1{i}]);
                                V2=load_untouch_nii([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',es2{i}]);
                                tmpval = abs(V1.img(:)-V2.img(:));
                                voldif_s(ksub,i) = median( tmpval(tmpval>eps) );
                            end
                        else
                            %file mismatch
                            erstr1=[]; for(i=1:numel(es1)) erstr1=[erstr1,', ',es1{i}]; end
                            erstr2=[]; for(i=1:numel(es2)) erstr2=[erstr2,', ',es2{i}]; end
                            error('cannot match optimally processed, spatially normalized datasets for %s.\nFiles for pipeline 1:%s\nFiles for pipeline 2:%s\n\n',InputStruct1(ksub).run(1).subjectprefix,errstr1(2:end),errstr2(2:end));
                        end
                        
                        
                    end
                    %%% now looking at outputs...
                    
                    if( max(voldif_u(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x contr)
                        
                        fprintf(fin,'Comparing all pipeline Procs(unnorm)..............................[!!]\n');
                        fprintf(fin,'\tDifferences in: \n');
                        olix=find( max(voldif_u,[],2) > 1E-6);
                        fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        olix=find( max(voldif_u,[],1) > 1E-6);
                        if(numel(olix)>10)
                            fprintf(fin,'\tPipelines (%u of %u): too many to list here!\n',length(olix),numel(eu1));
                        else
                            escat=[];
                            for(i=1:numel(eu1))
                                if(sum(i==olix)>0)
                                    etmp=eu1{i};
                                    ixct=strfind(etmp,  InputStruct1(ksub).run(1).subjectprefix);
                                    etmp(ixct:ixct+numel(  InputStruct1(ksub).run(1).subjectprefix )-1)=[];
                                    escat=[escat,', ',etmp];
                                end
                            end
                            fprintf(fin,'\tPipelines (%u of %u):%s\n',length(olix),numel(eu1),escat(2:end));
                        end
                        fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(voldif_u(voldif_u>1E-6)),min(voldif_u(voldif_u>1E-6)),max(voldif_u(voldif_u>1E-6)));
                    else
                        fprintf(fin,'Comparing all pipeline Procs(unnorm)..............................[OK]\n');
                    end
                    if( max(voldif_s(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x contr)
                        
                        fprintf(fin,'Comparing all pipeline Procs(spnorm)..............................[!!]\n');
                        fprintf(fin,'\tDifferences in: \n');
                        olix=find( max(voldif_s,[],2) > 1E-6);
                        fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                        olix=find( max(voldif_s,[],1) > 1E-6);
                        if(numel(olix)>10)
                            fprintf(fin,'\tPipelines (%u of %u): too many to list here!\n',length(olix),numel(es1));
                        else
                            escat=[];
                            for(i=1:numel(es1))
                                if(sum(i==olix)>0)
                                    etmp=es1{i};
                                    ixct=strfind(etmp,  InputStruct1(ksub).run(1).subjectprefix);
                                    etmp(ixct:ixct+numel(  InputStruct1(ksub).run(1).subjectprefix )-1)=[];
                                    escat=[escat,', ',etmp];
                                end
                            end
                            fprintf(fin,'\tPipelines (%u of %u):%s\n',length(olix),numel(es1),escat(2:end));
                        end
                        fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(voldif_s(voldif_s>1E-6)),min(voldif_s(voldif_s>1E-6)),max(voldif_s(voldif_s>1E-6)));
                    else
                        fprintf(fin,'Comparing all pipeline Procs(spnorm)..............................[OK]\n');
                    end
                    
                end
            end
        end
    else
        
        %%% NOTED THAT MULTIPLE STEPS ARE BEING SKIPPED FOR NO-ANALYSIS
        %%% PIPELINE...
        
        clear Nvox1 Nvox2;
        %per subjects: count number of outputted proc-data (v1/v2)
        for( ksub=1:Nsubject )
            
            M1=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask);
            M2=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask);
            
            Nvox1(ksub,1) = sum(M1.img(:));
            Nvox2(ksub,1) = sum(M2.img(:));
        end
        
        if( sum(Nvox1 == Nvox2)<Nsubject )
            tmpval = abs(Nvox1-Nvox2); tmpval=tmpval(tmpval>0);
            fprintf(fin,'Same number of brain voxels?......................................[!!]\n');
            fprintf(fin,'\tMismatched voxels for %u/%u subjects - have the brain masks changed?\n\tMedian [min, max] in runs showing change: %u [%u, %u]\n\tSPM comparisons will not be performed\n',...
                median(tmpval),min(tmpval),max(tmpval),sum(Nvox1 == Nvox2),Nsubject);
            out.flags.diffvox=1;
        else
            fprintf(fin,'Same number of brain voxels?......................................[OK]\n');
            out.flags.diffvox=0;
        end
        
        fprintf(fin,'\n============== 2. OPTIMIZATION OUTPUTS ==============\n\n');
        fprintf(fin,'\n\nNo optimization metric outputs - no analysis model!\n\n');
        fprintf(fin,'\n============== 3. INTERMEDIATE OUTPUTS (SPMS) ==============\n\n');
        fprintf(fin,'\n\nNo optimization SPM outputs - no analysis model!\n\n');
        fprintf(fin,'\n============== 4. INTERMEDIATE OUTPUTS (PROCESSED) ==============\n\n');
        
        if( prod(Nvox1==Nvox2) ) %% if all matching
            
            for(ksub=1:Nsubject)
                
                
                M1=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask); mask1=double(M1.img);
                M2=load_untouch_nii(InputStruct1(ksub).run(1).subjectmask); mask2=double(M2.img);
                
                %%% 1
                e=dir([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/Proc',InputStruct1(ksub).run(1).subjectprefix,'*']);
                clear eu1 es1; k1=0;k2=0;
                for(i=1:numel(e))
                    if(isempty(strfind(e(i).name,'_sNorm')))
                        k1=k1+1; eu1{k1}=e(i).name;
                    else  k2=k2+1; es1{k2}=e(i).name;
                    end
                end
                %%% 2
                e=dir([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/Proc',InputStruct2(ksub).run(1).subjectprefix,'*']);
                clear eu2 es2; k1=0;k2=0;
                for(i=1:numel(e))
                    if(isempty(strfind(e(i).name,'_sNorm')))
                        k1=k1+1; eu2{k1}=e(i).name;
                    else  k2=k2+1; es2{k2}=e(i).name;
                    end
                end
                
                if( prod(strcmpi(eu1,eu2)) )
                    %unnormalized files are same ... proceed to compare contents
                    for(i=1:numel(eu1))
                        V1=load_untouch_nii([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',eu1{i}]);
                        V2=load_untouch_nii([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',eu2{i}]);
                        tmpval = abs(V1.img(:)-V2.img(:));
                        voldif_u(ksub,i) = median( tmpval(tmpval>eps) );
                    end
                else
                    %file mismatch
                    erstr1=[]; for(i=1:numel(eu1)) erstr1=[erstr1,', ',eu1{i}]; end
                    erstr2=[]; for(i=1:numel(eu2)) erstr2=[erstr2,', ',eu2{i}]; end
                    error('cannot match optimally processed, unnormalized datasets for %s.\nFiles for pipeline 1:%s\nFiles for pipeline 2:%s\n\n',InputStruct1(ksub).run(1).subjectprefix,errstr1(2:end),errstr2(2:end));
                end
                %
                if( prod(strcmpi(es1,es2)) )
                    %unnormalized files are same ... proceed to compare contents
                    for(i=1:numel(es1))
                        V1=load_untouch_nii([InputStruct1(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',es1{i}]);
                        V2=load_untouch_nii([InputStruct2(ksub).run(1).Output_nifti_file_path,'/optimization_results/processed/',es2{i}]);
                        tmpval = abs(V1.img(:)-V2.img(:));
                        voldif_s(ksub,i) = median( tmpval(tmpval>eps) );
                    end
                else
                    %file mismatch
                    erstr1=[]; for(i=1:numel(es1)) erstr1=[erstr1,', ',es1{i}]; end
                    erstr2=[]; for(i=1:numel(es2)) erstr2=[erstr2,', ',es2{i}]; end
                    error('cannot match optimally processed, spatially normalized datasets for %s.\nFiles for pipeline 1:%s\nFiles for pipeline 2:%s\n\n',InputStruct1(ksub).run(1).subjectprefix,errstr1(2:end),errstr2(2:end));
                end
                
                
            end
            %%% now looking at outputs...
            
            if( max(voldif_u(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x contr)
                
                fprintf(fin,'Comparing all pipeline Procs(unnorm)..............................[!!]\n');
                fprintf(fin,'\tDifferences in: \n');
                olix=find( max(voldif_u,[],2) > 1E-6);
                fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                olix=find( max(voldif_u,[],1) > 1E-6);
                if(numel(olix)>10)
                    fprintf(fin,'\tPipelines (%u of %u): too many to list here!\n',length(olix),numel(eu1));
                else
                    escat=[];
                    for(i=1:numel(eu1))
                        if(sum(i==olix)>0)
                            etmp=eu1{i};
                            ixct=strfind(etmp,  InputStruct1(ksub).run(1).subjectprefix);
                            etmp(ixct:ixct+numel(  InputStruct1(ksub).run(1).subjectprefix )-1)=[];
                            escat=[escat,', ',etmp];
                        end
                    end
                    fprintf(fin,'\tPipelines (%u of %u):%s\n',length(olix),numel(eu1),escat(2:end));
                end
                fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(voldif_u(voldif_u>1E-6)),min(voldif_u(voldif_u>1E-6)),max(voldif_u(voldif_u>1E-6)));
            else
                fprintf(fin,'Comparing all pipeline Procs(unnorm)..............................[OK]\n');
            end
            
            if( max(voldif_s(:)) > 1E-6 ) %% any major differences in pipes (subj x ppln x contr)
                
                fprintf(fin,'Comparing all pipeline Procs(spnorm)..............................[!!]\n');
                fprintf(fin,'\tDifferences in: \n');
                olix=find( max(voldif_s,[],2) > 1E-6);
                fprintf(fin,'\tRuns (%u of %u): %s\n',length(olix),Nsubject,num2str(olix(:)','%u,'));
                olix=find( max(voldif_s,[],1) > 1E-6);
                if(numel(olix)>10)
                    fprintf(fin,'\tPipelines (%u of %u): too many to list here!\n',length(olix),numel(es1));
                else
                    escat=[];
                    for(i=1:numel(es1))
                        if(sum(i==olix)>0)
                            etmp=es1{i};
                            ixct=strfind(etmp,  InputStruct1(ksub).run(1).subjectprefix);
                            etmp(ixct:ixct+numel(  InputStruct1(ksub).run(1).subjectprefix )-1)=[];
                            escat=[escat,', ',etmp];
                        end
                    end
                    fprintf(fin,'\tPipelines (%u of %u):%s\n',length(olix),numel(es1),escat(2:end));
                end
                fprintf(fin,'\tMedian [min, max] abs. difference, in runs showing changes: %0.3f [%0.3f, %0.3f]\n',median(voldif_s(voldif_s>1E-6)),min(voldif_s(voldif_s>1E-6)),max(voldif_s(voldif_s>1E-6)));
            else
                fprintf(fin,'Comparing all pipeline Procs(spnorm)..............................[OK]\n');
            end
        end
        
    end
    %%%%%%%%%%%%%%
end

if( distatis_flag>0 )
    fprintf(fin,'\n\nChecking average cross-correlations between optimized SPMs...\n');
    if( isempty(CCmat) )
        fprintf(fin,'No measurable differences (all correlations ~1)\n');
    else
        fprintf(fin,'\n          CON(1) FIX(1) IND(1) CON(2) FIX(2) IND(2)\n');
        fprintf(fin,'       --------------------------------------------\n');
        colist={'CON(1)','FIX(1)','IND(1)','CON(2)','FIX(2)','IND(2)'};
        for(l=1:6)
            cvct = mean( CCmat(l,:,:),3 );
            fprintf(fin,' %s|   %s|\n',colist{l},num2str(cvct,'%0.2f   '));
            if(l<6)
                fprintf(fin,'       |                                          |\n');
            end
        end
        fprintf(fin,'       --------------------------------------------\n');
        fprintf(fin,'\n');
    end
end

fclose(fin);

%%%====================================================%%%
if( ~isempty(CCmat) )
    % running distatis to test for differences
    out = DISTATIS_quick( CCmat, 'dist', 'krbkrb', '---:::', 500 );
end
%%%====================================================%%%

%%
function x = get_numvols(file)

[p,f,e] = fileparts(file);
if(isempty(strfind(e,'.gz'))) %if not a zip file, read the header direct
    hdr = load_nii_hdr(file);
else %otherwise need to inflate and load .nii
    v = load_untouch_nii(file);
    hdr=v.hdr; clear v;
end

x = hdr.dime.dim(5);

function out = DISTATIS_quick( CCmat, type, colorset, typestr, Niter )
%
% =========================================================================
% DISTATIS_QUICK: rapid 3-way multidimensional scaling code, with
% bootstrapped confidence intervals
% =========================================================================
%
% Syntax:
%         out = DISTATIS_quick( CCmat, type, colorset, typestr, Niter )
%
% Input:
%         DISTATIS requires (KxK) matrices of distance/similarity measures
%         (e.g. a correlation matix), between K different conditions.
%         It requires a (KxK) matrix per subject, with at least 5 input
%         subjects, to obtain a stable solution
%
%          CCmat    : 3D matrix of "stacked" 2D distance/similarity matrices.
%                     For N subjects, CCmat has dimensions (K x K x N)
%          type     : string specifying what type of measure CCmat contains.
%                       'sim' = similarity matrix (e.g. correlations)
%                       'dist'= distance matrix   (e.g. euclidean distance)
%
%          colorset : a (K x 3) matrix, where the kth row denotes the RGB colour
%                     when plotting condition k (k=1...K).
%                     If colorset=[], colour values are randomly assigned
%          typestr  : vector of characters, where the kth element determines
%                     the line formatting for the ellipse of condition k (k=1...K)
%                     If typestr=[], all ellipses have format '-'
%          Niters   : number of bootstrap resamples, used to estimate
%                     components (at least Niters=50 recommended)
%
% Output:
%          Produces a multidimensional scaling plot of the K 95% CI ellipses
%          of each condition (obtained via Boostrap resampling). This plot
%          shows the first 2PCs of greatest variance in DISTATIS space.
%
%          Non-overlapping ellipses represent significantly different
%          groups; overlapping ellipses are not distinguishable in this
%          space (although they may be significantly different along some
%          other PC dimension).
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: Jan 15 2014
% ------------------------------------------------------------------------%
%

out     = [];
CI_prob = 0.95; % probability bounds
PCmax   = 2;    % max. number of PCdims to record in resampling

% input matrix dimensions
[K1 K2 N] = size(CCmat);
% catches for bad data formatting
if(K1==K2) K=K1;
else       error('similarity matrices must be square!');
end
if(N<5)    error('not enough subjects (dim3) to do Bootstrap resampling');
end

if( isempty(colorset) ) colorset = rand(K,3); end
if( isempty(typestr)  ) typestr  = repmat( '-', 1,K ); end
if( isempty(Niter)    ) Niter    = 1000; end

%% PREPARING DATA MATRICES

% centring matrix
E = eye(K) - ones(K)./K;
% initialize scp matrix
SCP = zeros( size(CCmat) );
%
if    ( strcmp(type, 'sim' ) ) for(n=1:N) SCP(:,:,n) = 0.5*E*CCmat(:,:,n)*E; end
elseif( strcmp(type, 'dist') ) for(n=1:N) SCP(:,:,n) =-0.5*E*CCmat(:,:,n)*E; end
else                           error('invalid datatype');
end
% normalize each matrix by first eigenvalue
for(n=1:N)
    [v l] = svd(SCP(:,:,n));
    SCP(:,:,n) = SCP(:,:,n) ./ l(1,1);
end

%% POPULATION COMPROMISE

% initialize RV matrix
RVmat = ones(N);
% populate RV matrix
for(i=1:N-1)
    for(j=i+1:N)
        rv = trace( SCP(:,:,i)*SCP(:,:,j) )./ sqrt( trace( SCP(:,:,i)*SCP(:,:,i) ) * trace( SCP(:,:,j)*SCP(:,:,j) ) );
        %
        RVmat(i,j) = rv;
        RVmat(j,i) = rv;
    end
end
% decomp on RV matrix
[p theta] = svd( RVmat );
% comproise weights
alfa = p(:,1) ./ sum(p(:,1));
% compromise SCP
S_plus = zeros(K);
for(n=1:N) S_plus = S_plus + SCP(:,:,n) .* alfa(n); end
%
[V L] = svd( S_plus );
% percentvar
PercLoad = diag(L(1:2,1:2))./trace(L);
% factor loadings
F_plus = V(:,1:PCmax) * sqrt(    L(1:PCmax,1:PCmax) );
RP     = V(:,1:PCmax) * sqrt(inv(L(1:PCmax,1:PCmax)));

%% BOOTSTRAPPED COMPROMISE

% initialize bootstrapped projection score matrix
F_boot = zeros( K, PCmax, Niter );

for( bb=1:Niter )
    
    disp(['boostrap iter:_',num2str(bb),'/',num2str(Niter)]);
    
    % bootstrapped SCP matrix
    list     = ceil(N*rand(N,1));
    SCP_boot = SCP(:,:,list);
    
    % initialize RV matrix
    RVmat = ones(N);
    % populate RV matrix
    for(i=1:N-1)
        for(j=i+1:N)
            rv = trace( SCP_boot(:,:,i)*SCP_boot(:,:,j) )./ sqrt( trace( SCP_boot(:,:,i)*SCP_boot(:,:,i) ) * trace( SCP_boot(:,:,j)*SCP_boot(:,:,j) ) );
            %
            RVmat(i,j) = rv;
            RVmat(j,i) = rv;
        end
    end
    % decomp on RV matrix
    [p theta] = svd( RVmat );
    % comproise weights
    alfa = p(:,1) ./ sum(p(:,1));
    % compromise SCP
    S_plus = zeros(K);
    for(n=1:N) S_plus = S_plus + SCP_boot(:,:,n) .* alfa(n); end
    
    % bootstrapped projection scores
    F_boot(:,:,bb) = S_plus * RP;
end

% dimensions: (boot x 2) x Kgroups
F_boot = permute( F_boot, [3 2 1] );
% set nbound for CIs
Nidx = round( CI_prob*Niter );

%% FIGURE PLOTTING

figure, hold on;
% plotting limits
pc1lim  = max(max(abs(F_boot(:,1,:)),[],3),[],1);
pc2lim  = max(max(abs(F_boot(:,2,:)),[],3),[],1);
pc12lim = max([pc1lim pc2lim]);

% go through and estimate CI parameters
for(k=1:K)
    
    % bootstrapped
    F_temp = F_boot(:,:,k);
    
    % xy coordinates
    x=F_temp(:,1);
    y=F_temp(:,2);
    % get mean
    xo = mean(x);
    yo = mean(y);
    % center coords
    x = x-xo;
    y = y-yo;
    
    % pca decomposition
    [u l v] = svd( [x y] );
    % PCspace coordinates, centered
    q  = u * l;
    % mahalanobis distance
    Md = sqrt(sum( q.^2 ./ repmat( var(q), [size(q,1) 1] ), 2));
    % list by increasing MD value
    index  = sortrows( [(1:Niter)' Md], 2 );
    ibound = index(Nidx,1);
    % point on CI boundary in pc-space
    qx_bnd = q(ibound,1);
    qy_bnd = q(ibound,2);
    
    % NB: 'a' = major (x) axis, 'b' = minor (y) axis
    % fractional scaling ratio 'c' of b/a
    c = l(2,2)./l(1,1);
    % scaling on 'a'
    a = sqrt( qx_bnd.^2 + (qy_bnd.^2)./(c.^2) );
    % scaling on 'b'
    b = a.*c;
    
    % angle of major axis (relative to x)
    phi = atan( v(2,1) ./ v(1,1) );
    % degrees in rad.
    t=0:0.01:2*pi;
    % trace of ellipse in (x,y) coordinates
    xe = xo + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    ye = yo + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    
    plot( xo,yo,'ok', 'markersize',4, 'markerfacecolor', colorset(k,:) ); hold on;
    plot( xe,ye, typestr(k), 'color', colorset(k,:), 'linewidth', 2);
end

plot( 0.9*sqrt(L(1,1))*[-1 1], [0 0], '-k', [0 0], 0.9*sqrt(L(2,2))*[-1 1], '-k' );

text(0.9*pc12lim, 0.1*pc12lim ,['var:',num2str(round(100*PercLoad(1))),'%']);
text(0.1*pc12lim, 0.9*pc12lim ,['var:',num2str(round(100*PercLoad(2))),'%']);

%
xlim([-1.1*pc12lim 1.1*pc12lim]);
ylim([-1.1*pc12lim 1.1*pc12lim]);
