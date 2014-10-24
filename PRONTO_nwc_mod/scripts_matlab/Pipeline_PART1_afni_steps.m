function Pipeline_PART1_afni_steps(InputStruct,input_pipeset_half,dospnormfirst)


global NUMBER_OF_CORES MULTI_RUN_INPUTFILE
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
setenv('FSLOUTPUTTYPE','NIFTI')


%% DEBUGING PURPOSE

% file = fopen(temp_InputStruct,'r');
% tline = fgetl(file);
% fclose(file);
% i1 = strfind(tline,'OUT=');
% i2 = strfind(tline,'STRUCT=');
% test_dir =  ['/scratch/strother_lab/babak/pipetestdir/' sprintf('%ld',(int64(round(now*1e9))))];
% mkdir_r(test_dir);
% tline = [tline(1:i1-1) ' OUT=' test_dir '/' InputStruct(1).Output_nifti_file_prefix ' ' tline(i2:end)];
% test_file = [sprintf('%ld',(int64(round(now*1e9)))) '.txt'];
% file = fopen(test_file,'w');
% fprintf(file,'%s',tline);
% fclose(file);
% unix(['./Pipeline_STEP1_test ' test_file ' ' temp_input_pipeset_half ' octave']);



%% read input file
N_Subject = numel(InputStruct);

for subject_counter = 1:N_Subject
    N_run = length(InputStruct(subject_counter).run);
    for run_counter = 1:N_run
        % Check whether any stage preprocessing has been performed
        pipeline_index = [];
        for pipe_counter = 1:size(input_pipeset_half,1)
            final_preprocessed_filename = sprintf('%s/%s_m%dc%dp%dt%ds%d.nii',InputStruct(subject_counter).run(run_counter).Output_nifti_file_path,InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
            if ~exist(final_preprocessed_filename,'file')
                pipeline_index = [pipeline_index;pipe_counter];
            end
        end
        pipeset_half = input_pipeset_half(pipeline_index,:);
        if isempty(pipeset_half)
            continue;
        end
        %
        
        uu = unique(pipeset_half(:,1));
        if length(uu)==2
            M_select = 2;
        else
            M_select = uu;
        end
        uu = unique(pipeset_half(:,2));
        if length(uu)==2
            C_select = 2;
        else
            C_select = uu;
        end
        uu = unique(pipeset_half(:,3));
        if length(uu)==2
            R_select = 2;
        else
            R_select = uu;
        end
        uu = unique(pipeset_half(:,4));
        if length(uu)==2
            T_select = 2;
        else
            T_select = uu;
        end
        smoset = unique(pipeset_half(:,5));
        
        OUTstr = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
        
        OUTstr_sub1 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path];
        OUTstr_sub2 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
        
        
        display(InputStruct(subject_counter).run(run_counter).Output_nifti_file_path)
        mkdir_r(InputStruct(subject_counter).run(run_counter).Output_nifti_file_path);
        mkdir_r([InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/masks']);
        mkdir_r([InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/diagnostic']);
        
        num_vol = get_numvols([InputStruct(subject_counter).run(run_counter).Input_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Input_nifti_file_prefix]);
        
        newSTART = InputStruct(subject_counter).run(run_counter).DROP_first;
        newEND = num_vol - InputStruct(subject_counter).run(run_counter).DROP_last - 1;
        
        
        
        delete([InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/temp*.nii']);
        unix([AFNI_PATH '3dTcat -prefix ' OUTstr '_drop.nii ''' InputStruct(subject_counter).run(run_counter).Input_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Input_nifti_file_prefix '[' num2str(newSTART) '..' num2str(newEND) ']''']);
        unix([AFNI_PATH '3dmerge -prefix ' OUTstr '_tempsmo.nii -doall -1blur_fwhm 6 ' OUTstr '_drop.nii']);
        unix([AFNI_PATH '3dAutomask -prefix ' OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii ' OUTstr '_tempsmo.nii']);
        maskpath_m0   =  [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask_nomc.nii'];
        outputpath_m0 =  [OUTstr_sub1 '/diagnostic/' OUTstr_sub2 '_smo'];
        bricknum = min_displace_brick([InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix '_tempsmo.nii'],maskpath_m0,outputpath_m0);
        
        %outputpath_m0 = [outputpath_m0 '_0ref_motbrick.txt'];
        mkdir_r([InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/mpe']);
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
        if ((M_select==0 || M_select==2))
            if ((C_select==0) || (C_select==2))
                stringval = [OUTstr '_m0c0'];
                copyfile([OUTstr '_drop.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            if ((C_select==1) || (C_select==2))
                stringval = [OUTstr '_m0c1'];
                maskpath_m1 = [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii'];
                censorpath  = [outputpath_m0  '_QC_output.mat'];
                interpolate_fmri([OUTstr '_drop.nii'],[stringval '.nii'],censorpath,'volume+motion');
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
        end
        if ((M_select==1 || M_select==2))
            if ((C_select==0) || (C_select==2))
                stringval = [OUTstr '_m1c0'];
                copyfile([OUTstr '_drop+mc.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
            
            if ((C_select==1) || (C_select==2))
                stringval = [OUTstr '_m1c1'];
                maskpath_m1 = [OUTstr_sub1 '/masks/' OUTstr_sub2 '_mask.nii'];
                censorpath  = [outputpath_m1  '_QC_output.mat'];
                interpolate_fmri([OUTstr '_drop+mc.nii'],[stringval '.nii'],censorpath,'volume+motion');
                kk = kk + 1;
                pplList1{kk} = stringval;
            end
        end
        kk = 0;
        if( (R_select==0) || (R_select==2) )
            for nn=1:length(pplList1)
                stringval = [pplList1{nn} 'p0'];
                copyfile([pplList1{nn} '.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList2{kk} = stringval;
            end
        end
        if( (R_select==1) || (R_select==2) )
            for nn=1:length(pplList1)
                stringval = [pplList1{nn} 'p1'];
                unix([AFNI_PATH '3dretroicor -prefix ' stringval '.nii -resp ' InputStruct(subject_counter).run(run_counter).PHYstr  '.resp.1D -card ' InputStruct(subject_counter).run(run_counter).PHYstr  '.puls.1D ' pplList1{nn} '.nii']);
                kk = kk + 1;
                pplList2{kk} = stringval;
            end
        end
        kk=0;
        if( (T_select==0) || (T_select==2) )
            for nn=1:length(pplList2)
                stringval = [pplList2{nn} 't0'];
                copyfile([pplList2{nn} '.nii'],[stringval '.nii']);
                kk = kk + 1;
                pplList3{kk} = stringval;
            end
        end
        if( (T_select==1) || (T_select==2) )
            for nn=1:length(pplList2)
                stringval = [pplList2{nn} 't1'];
                unix([AFNI_PATH '3dTshift -prefix ' stringval '.nii '  pplList2{nn} '.nii']);
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
                    unix([AFNI_PATH '3dmerge -prefix ' stringval '.nii -doall -1blur_fwhm ' num2str(smoset(nn)) ' ' pplList3{mm} '.nii']);
                end
                kk = kk + 1;
                pplList4{kk} = stringval;
            end
        end
        
        %% DEBUGING PURPOSE
        %     clear x
        %     list = dir([test_dir '/*.nii']);
        %     length(list)
        %     for i = 1:length(list)
        %         c = load_untouch_nii([test_dir '/' list(i).name]);
        %         m = load_untouch_nii([InputStruct(1).Output_nifti_file_path '/' list(i).name]);
        %         x{i,2} = sum((abs(double(m.img(:))-double(c.img(:))))~=0);
        %         x{i,1} = list(i).name;
        %     end
        %     x
        %     save(test_file,'x','-mat7-binary');
        %%
        delete([OUTstr '_drop*.nii']);
        for kk=1:length(pplList4)
            display(pplList4{kk});
        end
        pplList_delete = [pplList1 pplList2 pplList3];
        for kk=1:length(pplList_delete)
            delete([pplList_delete{kk} '.nii']);
        end
    end
end

for ksub = 1:numel(InputStruct)
    mkdir_r([InputStruct(ksub).run(1).Subject_OutputDirectory '/realignment']);
    for krun = 1:numel(InputStruct(ksub).run)
        if exist([InputStruct(ksub).run(krun).Output_nifti_file_path '/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc.nii' ],'file')
            if ~exist([InputStruct(ksub).run(krun).Subject_OutputDirectory,'/realignment/mean_',InputStruct(ksub).run(krun).Output_nifti_file_prefix,'_baseproc.nii'],'file')
                unix(sprintf('%sfslmaths %s/%s_baseproc.nii -Tmean %s/realignment/mean_%s_baseproc.nii',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,InputStruct(ksub).run(krun).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Subject_OutputDirectory,InputStruct(ksub).run(krun).Output_nifti_file_prefix));
            end
        end
    end
end

for subject_counter = 1:numel(InputStruct)
    for run_counter = 1:numel(InputStruct(subject_counter).run)
        OUTstr = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path '/' InputStruct(subject_counter).run(run_counter).Output_nifti_file_prefix];
        
        if dospnormfirst
            
            In_temp = InputStruct(subject_counter);
            In_temp.run = In_temp.run(run_counter);
            spatial_normalization(In_temp,[],[],3);
            % generate masks once again for the spatially normalized data
            
            OUTstr_sub1 = [InputStruct(subject_counter).run(run_counter).Output_nifti_file_path];
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
            final_preprocessed_filename = sprintf('%s_baseproc.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            final_preprocessed_filename_aligned = sprintf('%s_baseproc_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            movefile(final_preprocessed_filename,final_preprocessed_filename_aligned,'f');
            for pipe_counter=1:size(input_pipeset_half,1)
                final_preprocessed_filename = sprintf('%s_m%dc%dp%dt%ds%d.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                final_preprocessed_filename_aligned = sprintf('%s_m%dc%dp%dt%ds%d_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                movefile(final_preprocessed_filename,final_preprocessed_filename_aligned,'f');
            end
        end
        %
    end
end
     
if MULTI_RUN_INPUTFILE && ~dospnormfirst  % If the spatial normalization has not performed, then realignment is performed.
    for ksub = 1:numel(InputStruct)
        if numel(InputStruct(ksub).run)==1
            copyfile(sprintf('%s/realignment/mean_%s_baseproc.nii',InputStruct(ksub).run(1).Subject_OutputDirectory,InputStruct(ksub).run(1).Output_nifti_file_prefix),sprintf('%s/realignment/reg_mean_%s_baseproc.nii',InputStruct(ksub).run(1).Subject_OutputDirectory,InputStruct(ksub).run(1).Output_nifti_file_prefix));
            continue;
        end
        for krun = 1:numel(InputStruct(ksub).run) 
            if isempty(InputStruct(ksub).run(krun).Output_nifti_file_path)
                continue;
            end
            display(sprintf('Run 3dvolreg subject=%d,run=%d',ksub,krun))
            unix(sprintf('%sflirt -out %s/realignment/reg_mean_%s_baseproc.nii -omat %s/realignment/%03d_reg.mat -ref %s/realignment/mean_%s_baseproc.nii -in %s/realignment/mean_%s_baseproc.nii -dof 6',FSL_PATH,InputStruct(ksub).run(ksub).Subject_OutputDirectory,InputStruct(ksub).run(krun).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Subject_OutputDirectory,krun,InputStruct(ksub).run(krun).Subject_OutputDirectory,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(krun).Subject_OutputDirectory,InputStruct(ksub).run(krun).Output_nifti_file_prefix));
            
            for pipe_counter=1:size(input_pipeset_half,1)
                final_preprocessed_filename = sprintf('%s_m%dc%dp%dt%ds%d.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix,input_pipeset_half(pipe_counter,:));
                if ~exist(sprintf('%s/%s_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4)),'file')
                    display(sprintf('Run flirt %d',pipe_counter))
                    unix(sprintf('%sflirt -out %s/%s_aligned.nii -applyxfm -init %s/realignment/%03d_reg.mat -in %s/%s -ref %s/%s',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4),InputStruct(ksub).run(krun).Subject_OutputDirectory,krun,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename));
                end
            end
            final_preprocessed_filename = sprintf('%s_baseproc.nii',InputStruct(ksub).run(krun).Output_nifti_file_prefix);
            if ~exist(sprintf('%s/%s_aligned.nii',InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4)),'file')
                display(sprintf('Run 3drotate %d',pipe_counter))
                unix(sprintf('%sflirt -out %s/%s_aligned.nii -applyxfm -init %s/realignment/%03d_reg.mat -in %s/%s -ref %s/%s',FSL_PATH,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename(1:end-4),InputStruct(ksub).run(krun).Subject_OutputDirectory,krun,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename,InputStruct(ksub).run(krun).Output_nifti_file_path,final_preprocessed_filename));
            end
        end
    end
end



function x = get_numvols(file)

hdr = load_nii_hdr(file);
x = hdr.dime.dim(5);