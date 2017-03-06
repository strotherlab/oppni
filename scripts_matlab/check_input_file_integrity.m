function check_input_file_integrity(InputStruct, contrast_list_str, retroicor_flag, task_flag)
%
% Syntax:
%         check_input_file_integrity(InputStruct, contrast_list_str, retroicor_flag, task_flag)
%
% .runs a series of checks to ensure arguments (e.g. contrasts, physio correction, task regression)
%  are compatible with input files
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

% if "input" is still just a string, look for + load the file as a struct
if ~isstruct(InputStruct)
    [InputStruct,MULTI_RUN_INPUTFILE] = Read_Input_File(InputStruct);
end

% check whether the file exists
for ksub = 1:numel(InputStruct)
    for krun=1:numel(InputStruct(ksub).run)
        in_nii  = [InputStruct(ksub).run(krun).Input_nifti_filename];
        str_nii = InputStruct(ksub).run(krun).STRUCT_File;
        phystr1 = [InputStruct(ksub).run(krun).PHYstr '.resp.1D'];
        phystr2 = [InputStruct(ksub).run(krun).PHYstr '.puls.1D'];
        taskstr = InputStruct(ksub).run(krun).split_info_file;
        
        % does input file exist
        if ~exist(in_nii,'file')
            display(sprintf('ERROR: The input file %s does not exist',in_nii));
            sge_exit(100);
        end
        % does T1 structural file exist -need for spatial normalization
        if ~exist(str_nii,'file')
            display(sprintf('WARNING: The structrul input file %s does not exist, The spatial normalization can not be run!',str_nii));
        end
        % do properly formatted physio. data exist -need to run RETROICOR
        if retroicor_flag~=0
            if ~exist(phystr1,'file')
                display(sprintf('ERROR: The physiological input file %s does not exist, RETROICOR can not be run!',phystr1));
                sge_exit(100);
            end
            if ~exist(phystr2,'file')
                display(sprintf('ERROR: The physiological input file %s does not exist, RETROICOR can not be run!',phystr2));
                sge_exit(100);
            end
        end
        % is a task contrast specified -need to run TASK regression step
        if task_flag ~=0
           if isempty(contrast_list_str) || strcmpi(contrast_list_str,'NONE') 
                display(sprintf('ERROR: No contrast analysis, so TASK pipeline step cannot be run!'));
                sge_exit(100);            
           end
        end
        % does split_info file exist -need for all future analysis
        if ~exist(taskstr,'file')
            display(sprintf('ERROR: The split_info file %s does not exist',taskstr));
            sge_exit(100);            
        end
    end
end