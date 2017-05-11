function run_noddi_wmtract(InputStruct,EC_type,corrDiffParam,list,append_dir)


global NUMBER_OF_CORES
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
    CODE_PATH = fileparts(which('spatial_normalization.m'));
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

if ~isstruct(InputStruct)
    [InputStruct] = Read_Input_DTI(InputStruct);
end

if( nargin<4 ) list=[];
end

if(nargin<5) %% default is not assuming appended files
    append_dir = [];
end

if(strcmpi(EC_type,'eddy_correct')) EC_type='ec';
elseif( ~strcmpi(EC_type,'eddy') ) error('unrecognized correction method.');
end

% To run on HPCVL
% setenv('PATH',[getenv('PATH') ':' FSL_PATH ':' AFNI_PATH ':' FSL_PATH '/bin/']);
% setenv('FSLDIR',FSL_PATH);
% unix(['source ' FSL_PATH '/etc/fslconf/fsl.sh']);
% setenv('FSLDIR',FSL_PATH);
setenv('FSLOUTPUTTYPE','NIFTI');

%%% 
Nsubject = length(InputStruct); % Count the number of all runs and subjects
grouptemp_path = [InputStruct(1).run(1).Output_nifti_file_path '/dti_processed/group_template'];

if( ~isempty(append_dir))
    

    % handle accidental trailing slash
    if(strcmpi(append_dir(end),'/') || strcmpi(append_dir(end),'\')) append_dir=append_dir(1:end-1); end
    % full path of reference template
    grouptemp_path = [append_dir, '/dti_processed/group_template'];

    for ksub = 1:Nsubject

        if( isempty(list) || sum(list==ksub)>0 )

            disp(['running (appended) subj.#',num2str(ksub),' ID=',InputStruct(1).run(1).Output_nifti_file_path]);

            % path to dti data
            dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
            % noddi output directory
            mkdir_r([dti_path,'/noddi_out']);

            % 1.generating required transforms

            if( ~exist([dti_path,'/noddi_out/template_FA_native_thr0.20.nii'],'file') )

                % generate the inverse tranform (T3 to native space)
                unix([FSL_PATH [ 'invwarp --ref=',dti_path, '/spat_norm-append/FA_cln.nii --out=',dti_path,'/noddi_out/FA_T3_to_cln_invwarp --warp=',dti_path, '/spat_norm-append/FA_cln_to_T3_warp']]);
                % apply to the template
                unix([FSL_PATH ['applywarp --ref=',dti_path, '/spat_norm-append/FA_cln.nii --in=',grouptemp_path,'/template_FA_nl_symm.nii',' --warp=',dti_path,'/noddi_out/FA_T3_to_cln_invwarp --out=',dti_path,'/noddi_out/template_FA_native.nii']]);
                % maksing the transformed template
                unix([FSL_PATH ['fslmaths ',dti_path,'/noddi_out/template_FA_native.nii -thr 0.20 ',dti_path,'/noddi_out/template_FA_native_thr0.20.nii']]);
            else
                disp('transform exists.');
            end

            % 2.noddi steps now:
            unix(['gunzip ',dti_path,'/noddi_out/DTI_Multi_',EC_type,'.nii.gz']);
            
            if(corrDiffParam==0)
                noddi_inprefix = ['DTI_Multi_',EC_type];
            elseif(corrDiffParam==1)
                noddi_inprefix = ['DTI_Multi_',EC_type,'_adjust']; 
                % create adjusted file...
                if(~exist([dti_path,'/noddi_out/DTI_Multi_',EC_type,'_adjust.nii'],'file'))
                    V=load_untouch_nii([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.nii']);
                    M=load_untouch_nii([dti_path,'/noddi_out/template_FA_native_thr0.20.nii']); maskvol=double(M.img>0);
                    N_runs = length(InputStruct(ksub).run(1).Input_nifti_file_prefix);
                    clear numvol;
                    for(i=1:N_runs)
                        hdr=load_nii_hdr([InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{i},'.nii']);
                        numvol(i,1) = hdr.dime.dim(5);
                    end
                    ncumul     = cumsum(numvol);
                    numix(:,1) = [0; ncumul(1:end-1)]+1;
                    numix(:,2) = ncumul(1:end);
                    vimg_temp = V.img; 
                    for(i=1:N_runs)
                        vimg_temp(:,:,:,numix(i,1):numix(i,2)) = bsxfun(@times, bsxfun(@rdivide,  vimg_temp(:,:,:,numix(i,1):numix(i,2)), vimg_temp(:,:,:,numix(i,1))  ), maskvol);
                    end
                    V.img = vimg_temp;
                    save_untouch_nii(V,[dti_path,'/noddi_out/',noddi_inprefix,'.nii']);
                end
            end

            if( ~exist([dti_path,'/noddi_out/NODDI_fit_odi.nii'],'file'))

                % creating .mat ROI set
                CreateROI( [dti_path,'/noddi_out/',noddi_inprefix,'.nii'], [dti_path,'/noddi_out/template_FA_native_thr0.20.nii'], [dti_path,'/noddi_out/DTI_Multi_ROI.mat'] );
                % defining the NODDI protocol
                if(strcmpi(EC_type,'ec'))       protocol = FSL2Protocol([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.bval'], [dti_path,'/noddi_out/DTI_Multi_ec.bvec']);
                elseif(strcmpi(EC_type,'eddy')) protocol = FSL2Protocol([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.bval'], [dti_path,'/noddi_out/DTI_Multi_eddy.nii.eddy_rotated_bvecs']);
                end

                % make model
                noddi = MakeModel('WatsonSHStickTortIsoV_B0');

                % batch fitting, all brain voxels
                batch_fitting_single([dti_path,'/noddi_out/DTI_Multi_ROI'], protocol, noddi, [dti_path,'/noddi_out/paramFit.mat']);
                % store as nifti file
                SaveParamsAsNIfTI([dti_path,'/noddi_out/paramFit.mat'],[dti_path,'/noddi_out/DTI_Multi_ROI.mat'],[dti_path,'/noddi_out/template_FA_native_thr0.20.nii'],[dti_path,'/noddi_out/NODDI_fit']);

            else
                disp('noddi params already exist.');
            end

            if( ~exist([dti_path,'/spat_norm-append/NODDI_fit_odi_to_T3.nii'],'file') && ~exist([dti_path,'/spat_norm-append/NODDI_fit_odi_to_T3.nii.gz'],'file') )

                % transform back to template
                strip_struct = [dti_path, '/spat_norm-append/FA_cln'];
                % list of params to transform
                otherlist = {'ficvf','fiso','fmin','kappa','odi'};

                for(i=1:length(otherlist))
                    innorm = [dti_path,'/noddi_out/NODDI_fit_',otherlist{i}];
                    outnorm = [dti_path,'/spat_norm-append/NODDI_fit_',otherlist{i}];
                    unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
                end
            else
                disp('transformed noddi params already exist.');
            end

            %% additional transform of mask for distortion correction
            if( ~exist([dti_path,'/spat_norm-append/template_FA_native_thr0.20_to_T3.nii'],'file') && ~exist([dti_path,'/spat_norm-append/template_FA_native_thr0.20_to_T3.nii.gz'],'file') )

                % transform back to template
                strip_struct = [dti_path, '/spat_norm-append/FA_cln'];
                innorm = [dti_path,'/noddi_out/template_FA_native_thr0.20'];
                outnorm = [dti_path,'/spat_norm-append/template_FA_native_thr0.20'];
                unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
            else
                disp('transformed noddi mask already exist.');
            end
        end
    end

else %% ================standard, full registration =======================
    
    for ksub = 1:Nsubject

        if( isempty(list) || sum(list==ksub)>0 )

            disp(['running subj.#',num2str(ksub),' ID=',InputStruct(1).run(1).Output_nifti_file_path]);

            % path to dti data
            dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
            % noddi output directory
            mkdir_r([dti_path,'/noddi_out']);

            % 1.generating required transforms

            if( ~exist([dti_path,'/noddi_out/template_FA_native_thr0.20.nii'],'file') )

                % generate the inverse tranform (T3 to native space)
                unix([FSL_PATH [ 'invwarp --ref=',dti_path, '/spat_norm/FA_cln.nii --out=',dti_path,'/noddi_out/FA_T3_to_cln_invwarp --warp=',dti_path, '/spat_norm/FA_cln_to_T3_warp']]);
                % apply to the template
                unix([FSL_PATH ['applywarp --ref=',dti_path, '/spat_norm/FA_cln.nii --in=',grouptemp_path,'/template_FA_nl_symm.nii',' --warp=',dti_path,'/noddi_out/FA_T3_to_cln_invwarp --out=',dti_path,'/noddi_out/template_FA_native.nii']]);
                % maksing the transformed template
                unix([FSL_PATH ['fslmaths ',dti_path,'/noddi_out/template_FA_native.nii -thr 0.20 ',dti_path,'/noddi_out/template_FA_native_thr0.20.nii']]);
            else
                disp('transform exists.');
            end

            % 2.noddi steps now:
            unix(['gunzip ',dti_path,'/noddi_out/DTI_Multi_',EC_type,'.nii.gz']);

            if(corrDiffParam==0)
                noddi_inprefix = ['DTI_Multi_',EC_type];
            elseif(corrDiffParam==1)
                noddi_inprefix = ['DTI_Multi_',EC_type,'_adjust']; 
                % create adjusted file...
                if(~exist([dti_path,'/noddi_out/DTI_Multi_',EC_type,'_adjust.nii'],'file'))
                    V=load_untouch_nii([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.nii']);
                    M=load_untouch_nii([dti_path,'/noddi_out/template_FA_native_thr0.20.nii']); maskvol=double(M.img>0);
                    N_runs = length(InputStruct(ksub).run(1).Input_nifti_file_prefix);
                    clear numvol;
                    for(i=1:N_runs)
                        hdr=load_nii_hdr([InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{i},'.nii']);
                        numvol(i,1) = hdr.dime.dim(5);
                    end
                    ncumul     = cumsum(numvol);
                    numix(:,1) = [0; ncumul(1:end-1)]+1;
                    numix(:,2) = ncumul(1:end);
                    vimg_temp = V.img; 
                    for(i=1:N_runs)
                        vimg_temp(:,:,:,numix(i,1):numix(i,2)) = bsxfun(@times, bsxfun(@rdivide,  vimg_temp(:,:,:,numix(i,1):numix(i,2)), vimg_temp(:,:,:,numix(i,1))  ), maskvol);
                    end
                    V.img = vimg_temp;
                    save_untouch_nii(V,[dti_path,'/noddi_out/',noddi_inprefix,'.nii']);
                end
            end
            
            if( ~exist([dti_path,'/noddi_out/NODDI_fit_odi.nii'],'file'))

                % creating .mat ROI set
                CreateROI( [dti_path,'/noddi_out/',noddi_inprefix,'.nii'], [dti_path,'/noddi_out/template_FA_native_thr0.20.nii'], [dti_path,'/noddi_out/DTI_Multi_ROI.mat'] );
                % defining the NODDI protocol
                if(strcmpi(EC_type,'ec'))       protocol = FSL2Protocol([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.bval'], [dti_path,'/noddi_out/DTI_Multi_ec.bvec']);
                elseif(strcmpi(EC_type,'eddy')) protocol = FSL2Protocol([dti_path,'/noddi_out/DTI_Multi_',EC_type,'.bval'], [dti_path,'/noddi_out/DTI_Multi_eddy.nii.eddy_rotated_bvecs']);
                end
                % make model
                noddi = MakeModel('WatsonSHStickTortIsoV_B0');

                % batch fitting, all brain voxels
                batch_fitting_single([dti_path,'/noddi_out/DTI_Multi_ROI'], protocol, noddi, [dti_path,'/noddi_out/paramFit.mat']);
                % store as nifti file
                SaveParamsAsNIfTI([dti_path,'/noddi_out/paramFit.mat'],[dti_path,'/noddi_out/DTI_Multi_ROI.mat'],[dti_path,'/noddi_out/template_FA_native_thr0.20.nii'],[dti_path,'/noddi_out/NODDI_fit']);

            else
                disp('noddi params already exist.');
            end

            if( ~exist([dti_path,'/spat_norm/NODDI_fit_odi_to_T3.nii'],'file') && ~exist([dti_path,'/spat_norm/NODDI_fit_odi_to_T3.nii.gz'],'file') )

                % transform back to template
                strip_struct = [dti_path, '/spat_norm/FA_cln'];
                % list of params to transform
                otherlist = {'ficvf','fiso','fmin','kappa','odi'};

                for(i=1:length(otherlist))
                    innorm = [dti_path,'/noddi_out/NODDI_fit_',otherlist{i}];
                    outnorm = [dti_path,'/spat_norm/NODDI_fit_',otherlist{i}];
                    unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
                end
            else
                disp('transformed noddi params already exist.');
            end

            %% additional transform of mask for distortion correction
            if( ~exist([dti_path,'/spat_norm/template_FA_native_thr0.20_to_T3.nii'],'file') && ~exist([dti_path,'/spat_norm/template_FA_native_thr0.20_to_T3.nii.gz'],'file') )

                % transform back to template
                strip_struct = [dti_path, '/spat_norm/FA_cln'];
                innorm = [dti_path,'/noddi_out/template_FA_native_thr0.20'];
                outnorm = [dti_path,'/spat_norm/template_FA_native_thr0.20'];
                unix([FSL_PATH sprintf('applywarp --ref=%s/template_FA_nl_symm --in=%s --warp=%s_to_T3_warp.nii --out=%s_to_T3',grouptemp_path,innorm,strip_struct,outnorm)]);
            else
                disp('transformed noddi mask already exist.');
            end
        end
    end
end