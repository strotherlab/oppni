function Pipeline_PART1_afni_steps(InputStruct,input_pipeset_half,dospnormfirst,DEOBLIQUE,TPATTERN,TOFWHM)

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 163 $';
% CODE_DATE    = '$Date: 2014-12-03 17:30:16 -0500 (Wed, 03 Dec 2014) $';
% ------------------------------------------------------------------------%

global NUMBER_OF_CORES MULTI_RUN_INPUTFILE CODE_PROPERTY
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if ( ~exist('OCTAVE_VERSION','builtin') && exist('maxNumCompThreads') )
    maxNumCompThreads(NUMBER_OF_CORES);
end

global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART1_afni_steps.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
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

addpath(CODE_PATH)
addpath([CODE_PATH '/NIFTI_tools'])

%%
if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

if ~isnumeric(input_pipeset_half)
    [input_pipeset_half] = get_pipe_list(input_pipeset_half);
end
if numel(InputStruct(1).run)>1
    MULTI_RUN_INPUTFILE = true;
end
if isempty(InputStruct) % Exit if Input file is empty or not found
    return;
end

if nargin<3
    dospnormfirst = false;
end
if nargin<4
    DEOBLIQUE = 0;
end
if nargin<5
    TPATTERN = [];
end
if nargin<6
    TOFWHM   =  0;
end
setenv('FSLOUTPUTTYPE','NIFTI')
read_version;

%% read input file
N_Subject = numel(InputStruct);

for subject_counter = 1:N_Subject
    N_run = length(InputStruct(subject_counter).run);
    for run_counter = 1:N_run
        % Check whether any stage preprocessing has been performed
        pipeline_index = [];
        for pipe_counter = 1:size(input_pipeset_half,1)
            final_preprocessed_filename = sprintf('%s/intermediate_processed/afni_processed/%s_m%dc%dp%dt%ds%d.nii',InputStruct(subject_counter).run(run_counter).Output_nifti_file_path,InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
            if ~exist(final_preprocessed_filename,'file')
                pipeline_index = [pipeline_index;pipe_counter];
            end
        end
        pipeset_half = input_pipeset_half(pipeline_index,:);
        if isempty(pipeset_half)
            continue;
        end
        %

        M_select = unique(pipeset_half(:,1));
        C_select = unique(pipeset_half(:,2));
        R_select = unique(pipeset_half(:,3));
        T_select = unique(pipeset_half(:,4));
        smoset   = unique(pipeset_half(:,5));
        
        OUTstr = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
        
        OUTstr_sub1 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/intermediate_processed'];
        OUTstr_sub2 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
                
        display(InputStruct(subject_counter).run(run_counter).Output_nifti_file_path)
        mkdir_r([OUTstr_sub1 '/afni_processed']);
        mkdir_r([OUTstr_sub1 '/masks']);
        mkdir_r([OUTstr_sub1 '/diagnostic']);
        mkdir_r([OUTstr_sub1 '/mpe']);
        
        num_vol = get_numvols([InputStruct(subject_counter).run(run_counter).Input_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Input_nifti_file_prefix]);
        
        newSTART = InputStruct(subject_counter).run(run_counter).DROP_first;
        newEND = num_vol - InputStruct(subject_counter).run(run_counter).DROP_last - 1;
        
        delete([OUTstr_sub1 '/afni_processed/temp*.nii']);
        
        if(DEOBLIQUE==1)
            disp('deobliquing...');
            unix([AFNI_PATH '3dWarp -oblique2card -prefix ' OUTstr '_deob.nii -cubic ' InputStruct(subject_counter).run(run_counter).Input_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Input_nifti_file_prefix]);
            unix([AFNI_PATH '3dTcat -prefix ' OUTstr '_drop.nii ''' OUTstr '_deob.nii[' num2str(newSTART) '..' num2str(newEND) ']''']);
        else
            unix([AFNI_PATH '3dTcat -prefix ' OUTstr '_drop.nii ''' InputStruct(subject_counter).run(run_counter).Input_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Input_nifti_file_prefix '[' num2str(newSTART) '..' num2str(newEND) ']''']);
        end
                
        unix([AFNI_PATH '3dmerge -prefix ' OUTstr '_tempsmo.nii -doall -1blur_fwhm 6 ' OUTstr '_drop.nii']);
        unix([AFNI_PATH '3dAutomask -prefix ' OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii ' OUTstr '_tempsmo.nii']);
        maskpath_m0   =  [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii'];
        outputpath_m0 =  [OUTstr_sub1 '/diagnostic/' OUTstr_sub2 '_smo'];
        bricknum = min_displace_brick([OUTstr '_tempsmo.nii'],maskpath_m0,outputpath_m0);
        
        unix([AFNI_PATH '3dvolreg -prefix ' OUTstr '_drop+mc.nii -1Dfile ' OUTstr_sub1 '/mpe/' OUTstr_sub2 '_mpe -maxdisp1D ' OUTstr_sub1 '/mpe/' OUTstr_sub2 '_maxdisp -base '  num2str(bricknum)  ' '  OUTstr  '_drop.nii']);
        unix([AFNI_PATH '3dmerge -prefix ' OUTstr '_tempsmo_mc.nii -doall -1blur_fwhm 6 ' OUTstr '_drop+mc.nii']);
        unix([AFNI_PATH '3dAutomask -prefix ' OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii ' OUTstr '_tempsmo_mc.nii']);
        
        mpepath = [OUTstr_sub1 '/mpe/' OUTstr_sub2 '_mpe'];
        display('[~] Running Diagnostic 1 .. no MC .. in Matlab... ');
        
        diagnostic_fmri_pca([OUTstr '_tempsmo.nii'],maskpath_m0,mpepath,outputpath_m0);
        
        maskpath_m1   = [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii'];
        outputpath_m1 = [OUTstr_sub1 '/diagnostic/' OUTstr_sub2 '_mc+smo'];
        
        diagnostic_fmri_pca([OUTstr '_tempsmo_mc.nii'],maskpath_m1,mpepath,outputpath_m1);
        copyfile([OUTstr '_tempsmo_mc.nii'],[OUTstr '_baseproc.nii']);
                
        kk = 0;
        if ( ~isempty(find(M_select==0)) )
            if ( ~isempty(find(C_select==0)) ) % no censor
                stringval = [OUTstr '_m0c0'];
                copyfile([OUTstr '_drop.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ( ~isempty(find(C_select==1)) ) % full-volume despike
                stringval = [OUTstr '_m0c1'];
                censorpath  = [outputpath_m0  '_QC_output.mat'];
                interpolate_fmri([OUTstr '_drop.nii'],[stringval '.nii'],censorpath,'volume+motion');
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ( ~isempty(find(C_select==2)) ) % full-volume + pca cleanup
                
                % check if despiked exists
                tempstringval = [OUTstr '_tmp_cens'];
                % .make if not
                if ~exist([OUTstr '_m0c1.nii'],'file')
                    tempstringval = [OUTstr '_tmp_cens'];
                    censorpath  = [outputpath_m0  '_QC_output.mat'];
                    interpolate_fmri([OUTstr '_drop.nii'],[tempstringval '.nii'],censorpath,'volume+motion');
                % .otherwise copy over
                else
                    copyfile([OUTstr '_m0c1.nii'],[tempstringval '.nii']);
                end
                % smooth with standard 6mm kernel
                unix([AFNI_PATH '3dmerge -prefix ' tempstringval '+smo.nii -doall -1blur_fwhm 6 ' tempstringval '.nii']);
                % run pca
                disp('Running PCA decomposition for more aggressive spike correction...');
                gen_fmri_pca([tempstringval '+smo.nii'],[OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_nomc'],maskpath_m0,1);
                % 
                stringval = [OUTstr '_m0c2'];
                comppath_imag = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_nomc/spatial_PCs.nii'];
                comppath_temp = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_nomc/temporal_PCs'];
                component_auto_filter([tempstringval '+smo.nii'],[stringval,'.nii'], maskpath_m0, comppath_imag, comppath_temp);
                delete([tempstringval '.nii']);
                delete([tempstringval '+smo.nii']);
                
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ( ~isempty(find(C_select==3)) ) % full-volume + ica cleanup
                
                % check if despiked exists
                tempstringval = [OUTstr '_tmp_cens'];
                % .make if not
                if ~exist([OUTstr '_m0c1.nii'],'file')
                    tempstringval = [OUTstr '_tmp_cens'];
                    censorpath  = [outputpath_m0  '_QC_output.mat'];
                    interpolate_fmri([OUTstr '_drop.nii'],[tempstringval '.nii'],censorpath,'volume+motion');
                % .otherwise copy over
                else
                    copyfile([OUTstr '_m0c1.nii'],[tempstringval '.nii']);
                end
                % smooth with standard 6mm kernel
                unix([AFNI_PATH '3dmerge -prefix ' tempstringval '+smo.nii -doall -1blur_fwhm 6 ' tempstringval '.nii']);
                % run melodic
                disp('Running ICA decomposition for more aggressive spike correction...(this may take a while)...');
                if ~exist([OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_nomc/melodic_IC.nii'],'file')
                    unix([ FSL_PATH 'melodic -i ',[tempstringval '+smo.nii'],' -m ',maskpath_m0,' -o ',[OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_nomc'] ]);
                    unix(['gunzip ' [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_nomc/melodic_IC.nii.gz']]);
                end
                % 
                stringval = [OUTstr '_m0c3'];
                comppath_imag = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_nomc/melodic_IC.nii'];
                comppath_temp = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_nomc/melodic_Tmodes'];
                % run filtering on unsmoothed, interpolated data
                component_auto_filter([tempstringval '.nii'],[stringval,'.nii'], maskpath_m0, comppath_imag, comppath_temp);
                delete([tempstringval '.nii']);
                delete([tempstringval '+smo.nii']);
                
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
        end
        if ( ~isempty(find(M_select==1)) )
            if (  ~isempty(find(C_select==0))  )
                stringval = [OUTstr '_m1c0'];
                copyfile([OUTstr '_drop+mc.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if (  ~isempty(find(C_select==1))  )
                stringval = [OUTstr '_m1c1'];
                censorpath  = [outputpath_m1  '_QC_output.mat'];
                interpolate_fmri([OUTstr '_drop+mc.nii'],[stringval '.nii'],censorpath,'volume+motion');
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ( ~isempty(find(C_select==2)) ) % full-volume + pca cleanup
                
                % check if despiked exists
                tempstringval = [OUTstr '_tmp_mc+cens'];
                % .make if not
                if ~exist([OUTstr '_m1c1.nii'],'file')
                    tempstringval = [OUTstr '_tmp_cens'];
                    censorpath  = [outputpath_m1  '_QC_output.mat'];
                    interpolate_fmri([OUTstr '_drop+mc.nii'],[tempstringval '.nii'],censorpath,'volume+motion');
                % .otherwise copy over
                else
                    copyfile([OUTstr '_m1c1.nii'],[tempstringval '.nii']);
                end
                % smooth with standard 6mm kernel
                unix([AFNI_PATH '3dmerge -prefix ' tempstringval '+smo.nii -doall -1blur_fwhm 6 ' tempstringval '.nii']);
                % run pca               
                disp('Running PCA decomposition for more aggressive spike correction...');
                gen_fmri_pca([tempstringval '+smo.nii'],[OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_mc'],maskpath_m1,1);
                % 
                stringval = [OUTstr '_m1c2'];
                comppath_imag = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_mc/spatial_PCs.nii'];
                comppath_temp = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_pca_mc/temporal_PCs'];
                component_auto_filter([tempstringval '+smo.nii'],[stringval,'.nii'], maskpath_m1, comppath_imag, comppath_temp);
                delete([tempstringval '.nii']);
                delete([tempstringval '+smo.nii']);
                
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ( ~isempty(find(C_select==3)) ) % full-volume + ica cleanup
                
                % check if despiked exists
                tempstringval = [OUTstr '_tmp_mc+cens'];
                % .make if not
                if ~exist([OUTstr '_m1c1.nii'],'file')
                    tempstringval = [OUTstr '_tmp_cens'];
                    censorpath  = [outputpath_m1  '_QC_output.mat'];
                    interpolate_fmri([OUTstr '_drop+mc.nii'],[tempstringval '.nii'],censorpath,'volume+motion');
                % .otherwise copy over
                else
                    copyfile([OUTstr '_m1c1.nii'],[tempstringval '.nii']);
                end
                % smooth with standard 6mm kernel
                unix([AFNI_PATH '3dmerge -prefix ' tempstringval '+smo.nii -doall -1blur_fwhm 6 ' tempstringval '.nii']);
                % run melodic
                disp('Running ICA decomposition for more aggressive spike correction...(this may take a while)...');
                if ~exist([OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_mc/melodic_IC.nii'],'file')
                    unix([ FSL_PATH 'melodic -i ',[tempstringval '+smo.nii'],'-m ',maskpath_m1,' -o ',[OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_mc'] ]);
                    unix(['gunzip ' [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_mc/melodic_IC.nii.gz']]);
                end
                % 
                stringval = [OUTstr '_m0c3'];
                comppath_imag = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_mc/melodic_IC.nii'];
                comppath_temp = [OUTstr_sub1 '/diagnostic/',OUTstr_sub2,'_ica_mc/melodic_Tmodes'];
                % run filtering on unsmoothed, interpolated data
                component_auto_filter([tempstringval '.nii'],[stringval,'.nii'], maskpath_m1, comppath_imag, comppath_temp);
                delete([tempstringval '.nii']);
                delete([tempstringval '+smo.nii']);
                
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
        end
        kk = 0;
        if( ~isempty(find(R_select==0)) )
            for nn=1:length(pplList1)
                stringval = [pplList1{nn} 'p0'];
                copyfile([pplList1{nn} '.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList2{kk} = stringval;
            end
        end
        if( ~isempty(find(R_select==1)) )
            for nn=1:length(pplList1)
                stringval = [pplList1{nn} 'p1'];
                unix([AFNI_PATH '3dretroicor -prefix ' stringval '.nii -resp ' InputStruct(subject_counter).run(run_counter).PHYstr  '.resp.1D -card ' InputStruct(subject_counter).run(run_counter).PHYstr  '.puls.1D ' pplList1{nn} '.nii']);
                kk = kk + 1;
                pplList2{kk} = stringval;
            end
        end
        kk=0;
        if( ~isempty(find(T_select==0)) )
            for nn=1:length(pplList2)
                stringval = [pplList2{nn} 't0'];
                copyfile([pplList2{nn} '.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList3{kk} = stringval;
            end
        end
        if( ~isempty(find(T_select==1)) )
            for nn=1:length(pplList2)
                stringval = [pplList2{nn} 't1'];
                
                if( isempty(TPATTERN) || strcmpi(TPATTERN,'auto_hdr') ) % allow users to input slice-timing pattern if not available
                        unix([AFNI_PATH '3dTshift -prefix ' stringval '.nii '  pplList2{nn} '.nii']);
                else    unix([AFNI_PATH '3dTshift -prefix ' stringval '.nii -tpattern ',TPATTERN,' ', pplList2{nn} '.nii']);                    
                end
                kk = kk + 1;
                pplList3{kk} = stringval;
            end
        end
        kk = 0;
        for nn=1:length(smoset)
            for mm = 1:length(pplList3)
                stringval = [pplList3{mm} 's' num2str(smoset(nn))];
                if smoset(nn)==0
                    copyfile([pplList3{mm} '.nii'],[stringval '.nii']);
                else
                    
                    if( TOFWHM== 0 )
                    
                        unix([AFNI_PATH '3dmerge -prefix ' stringval '.nii -doall -1blur_fwhm ' num2str(smoset(nn)) ' ' pplList3{mm} '.nii']);
                    
                    elseif( TOFWHM== 1 )
                    
                        VV = load_untouch_nii([pplList3{mm} '.nii']);
                        VV.img = double(VV.img); %% format as double, for compatibility

                        Time_series = reshape(VV.img,size(VV.img,1)*size(VV.img,2)*size(VV.img,3), size(VV.img,4));
                        split_info = InputStruct(subject_counter).run(run_counter).split_info;
                        Design_matrix = split_info.design_mat; %% replaces BA's HRFdesign
                        trends = get_legendre(1+floor( (split_info.TR_MSEC/1000) * (size(VV.img,4)/2) ./ 150 ),size(VV.img,4));
                        Null_space  = [Design_matrix trends];Pn = eye(size(VV.img,4))-Null_space*inv(Null_space'*Null_space)*Null_space';
                        Time_series = Time_series*Pn;
                        Time_series = reshape(Time_series,size(VV.img,1),size(VV.img,2),size(VV.img,3), size(VV.img,4));
                        VV.img = Time_series;
                        VV.hdr.dime.datatype = 16;
                        VV.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
                        save_untouch_nii(VV,[pplList3{mm} '_blurmaster.nii']);
                        if isempty(strfind(pplList3{mm},'_m0'))

                            mmmask_name = [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii'];
                        else
                            mmmask_name = [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii'];
                        end

                        unix([AFNI_PATH '3dBlurToFWHM -blurmaster ' [pplList3{mm} '_blurmaster.nii'] ' -mask ' mmmask_name ' -input ' pplList3{mm} '.nii -prefix ' stringval '.nii -FWHM ' num2str(smoset(nn))]);
                        %unix([AFNI_PATH '3dmerge -prefix ' stringval '.nii -doall -1blur_fwhm ' num2str(smoset(nn)) ' ' pplList3{mm} '.nii']);
                        delete([pplList3{mm} '_blurmaster.nii']);                    
                    else
                        error('unrecognized smoothing TOFWHM option');
                    end
                end
                kk = kk + 1;
                pplList4{kk} = stringval;
            end
        end
        
        delete([OUTstr '_drop*.nii']);
        for kk=1:length(pplList4)
            display(pplList4{kk});
        end
        pplList_delete = [pplList1 pplList2 pplList3];
        for kk=1:length(pplList_delete)
            delete([pplList_delete{kk} '.nii']);
        end
        % remove "no-mc" mask as it is no longer needed 
        % --> COMMENTED OUT FOR COMPATIBILITY WITH SGE QUEUEING
        delete([OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii']);
    end
end


for subject_counter = 1:numel(InputStruct)
    for run_counter = 1:numel(InputStruct(subject_counter).run)
        OUTstr = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
        
        if dospnormfirst
            
            In_temp = InputStruct(subject_counter);
            In_temp.run = In_temp.run(run_counter);
            spatial_normalization(In_temp,[],[],3);
            % generate masks once again for the spatially normalized data
            
            OUTstr_sub1 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/intermediate_processed'];
            OUTstr_sub2 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
            
            % removing previously build mask
            delete([OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii']);
            delete([OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii']);
            
            % generate new masks in the common space
            unix([AFNI_PATH '3dAutomask -prefix ' OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii ' OUTstr '_tempsmo.nii']);
            unix([AFNI_PATH '3dAutomask -prefix ' OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii ' OUTstr '_tempsmo_mc.nii']);
        end
        
        % removing extra files
        delete([OUTstr '_tempsmo*.nii']);
    end
end


if MULTI_RUN_INPUTFILE && dospnormfirst   % If the spatial normalization has performed, then there is no need to realign runs 
    % change name to aligned
    for ksub=1:numel(InputStruct)
        for krun=1:numel(InputStruct(ksub).run)
            final_preprocessed_filename = sprintf('%s/intermediate_processed/afni_processed/%s_baseproc.nii',InputStruct(ksub).run(1).Output_nifti_file_path, InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            final_preprocessed_filename_aligned = sprintf('%s/intermediate_processed/afni_processed/%s_baseproc_aligned.nii',InputStruct(ksub).run(1).Output_nifti_file_path, InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            movefile(final_preprocessed_filename,final_preprocessed_filename_aligned,'f');
            for pipe_counter=1:size(input_pipeset_half,1)
                final_preprocessed_filename = sprintf('%s/intermediate_processed/afni_processed/%s_m%dc%dp%dt%ds%d.nii',InputStruct(ksub).run(krun).Output_nifti_file_path, InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                final_preprocessed_filename_aligned = sprintf('%s/intermediate_processed/afni_processed/%s_m%dc%dp%dt%ds%d_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_path, InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                movefile(final_preprocessed_filename,final_preprocessed_filename_aligned,'f');
            end
        end
        %
    end
end
     
if MULTI_RUN_INPUTFILE && ~dospnormfirst  % If the spatial normalization has not performed, then realignment is performed.
    
    % create unmasked means as a reference, in order to do alignments
    for ksub = 1:numel(InputStruct)
        mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/align_multirun']);
        for krun = 1:numel(InputStruct(ksub).run)

            if exist([InputStruct(ksub).run(krun).Output_nifti_file_path '/intermediate_processed/afni_processed/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc.nii' ],'file')

                mean_file_name = [InputStruct(ksub).run(krun).Output_nifti_file_path,'/intermediate_processed/align_multirun/mean_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_unmasked.nii' ];

                if ~exist(mean_file_name,'file')
                    unix(sprintf('%sfslmaths %s/intermediate_processed/afni_processed/%s_baseproc.nii -Tmean %s',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix,mean_file_name));
                    unix(['gunzip -f -d ' mean_file_name '.gz']);
                end
            end
        end
    end
    
    for ksub = 1:numel(InputStruct)
        if numel(InputStruct(ksub).run)==1
            copyfile(sprintf('%s/intermediate_processed/align_multirun/mean_%s_unmasked.nii',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix),sprintf('%s/intermediate_processed/align_multirun/reg_mean_%s_unmasked.nii',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix));
            continue;
        end
        for krun = 1:numel(InputStruct(ksub).run) 
            if isempty(InputStruct(ksub).run(krun).Output_nifti_file_path)
                continue;
            end
            display(sprintf('Run 3dvolreg subject=%d,run=%d',ksub,krun))
            unix(sprintf('%sflirt -out %s/intermediate_processed/align_multirun/reg_mean_%s_unmasked.nii -omat %s/intermediate_processed/align_multirun/%03d_reg.mat -ref %s/intermediate_processed/align_multirun/mean_%s_unmasked.nii -in %s/intermediate_processed/align_multirun/mean_%s_unmasked.nii -dof 6',FSL_PATH,InputStruct(ksub).run(ksub).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Output_nifti_file_path,krun,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix));
            
            for pipe_counter=1:size(input_pipeset_half,1)
                final_preprocessed_filename = sprintf('%s_m%dc%dp%dt%ds%d.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                if ~exist(sprintf('%s/intermediate_processed/afni_processed/%s_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4)),'file')
                    display(sprintf('Run flirt %d',pipe_counter))
                    unix(sprintf('%sflirt -out %s/intermediate_processed/afni_processed/%s_aligned.nii -applyxfm -init %s/intermediate_processed/align_multirun/%03d_reg.mat -in %s/intermediate_processed/afni_processed/%s -ref %s/intermediate_processed/afni_processed/%s',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4),InputStruct(ksub).run(krun).Output_nifti_file_path,krun,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename));
                end
            end
            final_preprocessed_filename = sprintf('%s_baseproc.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            if ~exist(sprintf('%s/intermediate_processed/afni_processed/%s_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4)),'file')
                display(sprintf('Run 3drotate %d',pipe_counter))
                unix(sprintf('%sflirt -out %s/intermediate_processed/afni_processed/%s_aligned.nii -applyxfm -init %s/intermediate_processed/spat_norm/%03d_reg.mat -in %s/intermediate_processed/afni_processed/%s -ref %s/intermediate_processed/afni_processed/%s',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4),InputStruct(ksub).run(krun).Output_nifti_file_path,krun,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename));
            end
        end
    end
end



function x = get_numvols(file)

[p,f,e] = fileparts(file);
if(isempty(strfind(e,'.gz'))) %if not a zip file, read the header direct
    hdr = load_nii_hdr(file);
else %otherwise need to inflate and load .nii
    v = load_untouch_nii(file);
    hdr=v.hdr; clear v;
end

x = hdr.dime.dim(5);