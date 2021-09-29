function spatial_normalization_noise_roi(InputStruct)

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 163 $';
% CODE_DATE    = '$Date: 2014-12-03 17:30:16 -0500 (Wed, 03 Dec 2014) $';
% ------------------------------------------------------------------------%


global NUMBER_OF_CORES CODE_PROPERTY
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if ( ~exist('OCTAVE_VERSION','builtin') && exist('maxNumCompThreads') )
    maxNumCompThreads(NUMBER_OF_CORES);
end

global CODE_PATH AFNI_PATH FSL_PATH MULTI_RUN_INPUTFILE
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('spatial_normalization_noise_roi.m'));
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

if ~isdeployed
    addpath(CODE_PATH)
    addpath([CODE_PATH '/NIFTI_tools'])
end
read_version;

if ~isstruct(InputStruct)
    %DEBUG
    disp(sprintf("LMP-DEBUG: RE-Read_Input_File"));

    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

setenv('FSLOUTPUTTYPE','NIFTI')

for ksub = 1:numel(InputStruct)
    if ~isempty(InputStruct(ksub).run(1).Noise_ROI)
        input_nifti_file  = InputStruct(ksub).run(1).Noise_ROI;

        %DEBUG
        disp(sprintf("LMP-DEBUG: Noise_ROI file = %s",input_nifti_file));

        [path_temp,input_nifti_name]  = fileparts(input_nifti_file);input_nifti_name = [input_nifti_name '.nii'];
        output_nifti_file = [InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/noise_roi/' input_nifti_name ];
        if ~exist(output_nifti_file,'file')
            [tmp,STRUCT_Name,ext] = fileparts(InputStruct(ksub).run(1).STRUCT_File);
            mean_file_name = sprintf('%s/intermediate_processed/spat_norm/mean_%s_baseproc.nii',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
            if ~exist(mean_file_name,'file') % generate mean volume if not exist, then mask it.
                nii      = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_baseproc.nii']);
                nii_mask = load_untouch_nii([InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/masks/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_mask.nii']);
                nii.img =  mean(double(nii.img),4);
                nii.img(nii_mask.img==0) = 0;
                nii.hdr.dime.dim([1 5]) = [3 1];
                nii.hdr.dime.datatype = 16;
                mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm']);
                nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
                save_untouch_nii(nii,mean_file_name);
            end
            
            mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/spat_norm']);
            mkdir_r([InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/noise_roi']);
            out_trans = sprintf('%s/intermediate_processed/spat_norm/Transmat_T1toEPI_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
            if ~exist(out_trans,'file')
                unix([AFNI_PATH sprintf('3dSkullStrip -prefix %s/intermediate_processed/spat_norm/%s_strip.nii -input %s',InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).STRUCT_File)]);
                
                unix([FSL_PATH sprintf('flirt -in %s -ref %s/intermediate_processed/spat_norm/%s_strip.nii -omat %s/intermediate_processed/spat_norm/Transmat_EPItoT1_%s.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6 -interp sinc -sincwidth 7 -sincwindow hanning', ...
                    mean_file_name,InputStruct(ksub).run(1).Output_nifti_file_path,STRUCT_Name,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
                
                unix([FSL_PATH sprintf('convert_xfm -inverse  %s/intermediate_processed/spat_norm/Transmat_EPItoT1_%s.mat -omat %s/intermediate_processed/spat_norm/Transmat_T1toEPI_%s.mat',...
                    InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix,InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix)]);
            end
            
            input_nifti_file  = InputStruct(ksub).run(1).Noise_ROI;
            [path_temp,input_nifti_name,ext]  = fileparts(input_nifti_file);
            if strcmpi(ext,'.gz')
                 gunzip(input_nifti_file);
                [path_temp,input_nifti_name,ext]  = fileparts(input_nifti_name);
            end
            
            input_nifti_name = [input_nifti_name ext];
            output_nifti_file = [InputStruct(ksub).run(1).Output_nifti_file_path '/intermediate_processed/noise_roi/' input_nifti_name ];
            
            nii = load_untouch_nii(input_nifti_file);
            nii.img = smooth(double(nii.img));
            nii.hdr.hist.descrip = CODE_PROPERTY.NII_HEADER;
            save_untouch_nii(nii,[output_nifti_file(1:end-4) '_smoothed.nii']);
            
            input_nifti_file = [output_nifti_file(1:end-4) '_smoothed.nii'];
            reference_file    = [InputStruct(ksub).run(1).Output_nifti_file_path,'/intermediate_processed/afni_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix,'_baseproc.nii'];
            transform         = sprintf('%s/intermediate_processed/spat_norm/Transmat_T1toEPI_%s.mat',InputStruct(ksub).run(1).Output_nifti_file_path,InputStruct(ksub).run(1).Output_nifti_file_prefix);
            unix([FSL_PATH sprintf('flirt -in %s -applyxfm -interp trilinear -ref %s -init %s -out %s',input_nifti_file,reference_file,transform,output_nifti_file)]);
            delete(input_nifti_file);
        end
    end
end



function b = smooth(a)

Ind = find(a~=0);
[I,J,K] = ind2sub(size(a),Ind);
for i = 1:length(I)
    a(I(i)-1:I(i)+1,J(i)-1:J(i)+1,K(i)-1:K(i)+1) = a(I(i)-1:I(i)+1,J(i)-1:J(i)+1,K(i)-1:K(i)+1) + 2;
    a(I(i)-2:I(i)+2,J(i)-2:J(i)+2,K(i)-2:K(i)+2) = a(I(i)-2:I(i)+2,J(i)-2:J(i)+2,K(i)-2:K(i)+2) + 1;
end
b = a;
    

