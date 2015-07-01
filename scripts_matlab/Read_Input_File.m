function [InputStruct, MULTI_RUN_INPUTFILE] = Read_Input_File(inputfile)

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

fid     = fopen(inputfile);
if fid==-1
    InputStruct = [];
    MULTI_RUN_INPUTFILE = [];
    return;
end
% read in first line
tline               = fgetl(fid);
if ~ischar(tline)
    InputStruct = [];
    MULTI_RUN_INPUTFILE = [];
    return;
end

ksub     = 0;
MULTI_RUN_INPUTFILE = false;

%% read input file
while ischar(tline) 
        ksub = ksub  +1;
        [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix,split_info_file,STRUCT_File,PHYstr,Noise_ROI,DROP_first,DROP_last] = Parse_Input_File(tline);
        for krun = 1:length(Input_nifti_file_path)
            InputStruct(ksub).run(krun).Input_nifti_file_path             = Input_nifti_file_path{krun};
            InputStruct(ksub).run(krun).Input_nifti_file_prefix           = Input_nifti_file_prefix{krun};
            InputStruct(ksub).run(krun).Output_nifti_file_path            = Output_nifti_file_path{krun};
            InputStruct(ksub).run(krun).Output_nifti_file_prefix          = Output_nifti_file_prefix{krun};
            InputStruct(ksub).run(krun).split_info_file                   = split_info_file{krun};
            InputStruct(ksub).run(krun).STRUCT_File                       = STRUCT_File{krun};
            InputStruct(ksub).run(krun).PHYstr                            = PHYstr{krun};
            InputStruct(ksub).run(krun).DROP_first                        = DROP_first{krun};
            InputStruct(ksub).run(krun).DROP_last                         = DROP_last{krun};
            InputStruct(ksub).run(krun).Noise_ROI                         = Noise_ROI{krun};
            InputStruct(ksub).run(krun).subjectprefix                     = ['_' Output_nifti_file_prefix{1}];
            InputStruct(ksub).run(krun).subjectmask                       = [Output_nifti_file_path{1}  '/intermediate_processed/masks/' Output_nifti_file_prefix{1} '_mask.nii'];
            
            % load to get task-type for integrity check
            [split_info] = Parse_Split_Info(split_info_file{krun});
            InputStruct(ksub).run(krun).TaskType                          = split_info.type;
        end
        if length(Input_nifti_file_path)>1
            MULTI_RUN_INPUTFILE=true;
        end
 
    tline = fgetl(fid);
    if isempty(tline)
        tline = fgetl(fid);
    end
end
fclose(fid); 


% Covert subjects or runs

%InputStruct = Apply_Optimization_Scheme(InputStruct,'subject');
% for furthur development

% 
outdir = InputStruct(1).run(1).Output_nifti_file_path;
l1 = length(outdir);
for is = 2:length(InputStruct)
    l2 = length(InputStruct(is).run(1).Output_nifti_file_path);
    l1 = length(outdir);
    l = min(l1,l2);
    ind = find(outdir(1:l)~=InputStruct(is).run(1).Output_nifti_file_path(1:l));
    if isempty(ind) 
        ind(1) = l+1; 
    end
    outdir = outdir(1:ind(1)-1);
end
Group_OutputDirectory = outdir;


for ksub = 1:size(InputStruct)
    for krun = 1:length(InputStruct(ksub).run)
        InputStruct(ksub).run(krun).Group_OutputDirectory = Group_OutputDirectory;
    end
end

% function OutputStruct = Apply_Optimization_Scheme(InputStruct,optimization_scheme)
% 
% count = 0;
% for ksub = 1:length(InputStruct)
%     for krun =1:numel(InputStruct(ksub).run)
%         count = count + 1;
%         Input_nifti_file_path{count}    = InputStruct(ksub).run(krun).Input_nifti_file_path;
%         Input_nifti_file_prefix{count}  = InputStruct(ksub).run(krun).Input_nifti_file_prefix;
%         Output_nifti_file_path{count}   = InputStruct(ksub).run(krun).Output_nifti_file_path;
%         Output_nifti_file_prefix{count} = InputStruct(ksub).run(krun).Output_nifti_file_prefix;
%         split_info_file{count}          = InputStruct(ksub).run(krun).split_info_file;
%     
%         STRUCT_File{count}  = InputStruct(ksub).run(krun).STRUCT_File;
%         
%         PHYstr{count}      = InputStruct(ksub).run(krun).PHYstr;
%         DROP_first{count}   = InputStruct(ksub).run(krun).DROP_first;
%         DROP_last{count}    = InputStruct(ksub).run(krun).DROP_last;
%         Noise_ROI{count}    = InputStruct(ksub).run(krun).Noise_ROI;
%         Subject_OutputDirectory{count} = InputStruct(ksub).run(krun).Subject_OutputDirectory;
%         subjectprefix{count}  = InputStruct(ksub).run(krun).subjectprefix;
%         subjectmask{count}    = InputStruct(ksub).run(krun).subjectmask;
%     end
% end
% Num_runs = length(Input_nifti_file_path);
% % those subjects that have same 
% 
% if strcmpi(optimization_scheme,'run')
%     for count = 1:Num_runs
%         OutputStruct(count).run(1).Input_nifti_file_path = Input_nifti_file_path{count};
%         OutputStruct(count).run(1).Input_nifti_file_prefix = Input_nifti_file_prefix{count};
%         OutputStruct(count).run(1).Output_nifti_file_path = Output_nifti_file_path{count};
%         OutputStruct(count).run(1).Output_nifti_file_prefix = Output_nifti_file_prefix{count};
%         OutputStruct(count).run(1).split_info_file = split_info_file{count};
%         OutputStruct(count).run(1).STRUCT_File     =  STRUCT_File{count};
%         OutputStruct(count).run(1).PHYstr   = PHYstr{count};
%         OutputStruct(count).run(1).DROP_first     =  DROP_first{count};
%         OutputStruct(count).run(1).DROP_last   = DROP_last{count};
%         OutputStruct(count).run(1).Subject_OutputDirectory   = [Output_nifti_file_path{count} '/matfiles'];
%         OutputStruct(count).run(1).subjectprefix     =  ['_' Output_nifti_file_prefix{count}];
%         OutputStruct(count).run(1).subjectmask   = [Output_nifti_file_path{count}  '/masks/' Output_nifti_file_prefix{count} '_mask.nii'];      
%         OutputStruct(ksub).run(krun).Noise_ROI    = Noise_ROI{count};
%     end
% else
%     
%     SO = [STRUCT_File;Output_nifti_file_path];
%     
%     u = my_unique(SO);
%     
%     
%     N_Subject = size(u,2);
%     N_run     = zeros(1,N_Subject);
%     for count = 1:Num_runs
%         ksub = find(strcmp(STRUCT_File{count},u(1,:)) & strcmp(Output_nifti_file_path{count},u(2,:)));
%         N_run(ksub) = N_run(ksub) + 1;
%         krun = N_run(ksub);
%         OutputStruct(ksub).run(krun).Input_nifti_file_path = Input_nifti_file_path{count};
%         OutputStruct(ksub).run(krun).Input_nifti_file_prefix = Input_nifti_file_prefix{count};
%         OutputStruct(ksub).run(krun).Output_nifti_file_path = Output_nifti_file_path{count};
%         OutputStruct(ksub).run(krun).Output_nifti_file_prefix = Output_nifti_file_prefix{count};
%         OutputStruct(ksub).run(krun).split_info_file = split_info_file{count};
%         OutputStruct(ksub).run(krun).STRUCT_File     =  STRUCT_File{count};
%         OutputStruct(ksub).run(krun).PHYstr   = PHYstr{count};
%         OutputStruct(ksub).run(krun).DROP_first     =  DROP_first{count};
%         OutputStruct(ksub).run(krun).DROP_last   = DROP_last{count};
%         OutputStruct(ksub).run(krun).Subject_OutputDirectory   = [InputStruct(ksub).run(1).Output_nifti_file_path '/matfiles'];
%         OutputStruct(ksub).run(krun).subjectprefix     =  ['_' InputStruct(ksub).run(1).Output_nifti_file_prefix];
%         OutputStruct(ksub).run(krun).subjectmask   = [OutputStruct(ksub).run(1).Output_nifti_file_path '/masks/' OutputStruct(ksub).run(1).Output_nifti_file_prefix '_mask.nii'];    
%         OutputStruct(ksub).run(krun).Noise_ROI    = Noise_ROI{count};
%     end
% end
% 
% 
% function u = my_unique(SO)
% sz = 0;
% u = {[];[]};
% for count = 1:size(SO,2)
%     if ~any(strcmp(SO{1,count},u(1,:)) & strcmp(SO{2,count},u(2,:)))
%         sz  =sz + 1;
%         u{1,sz} = SO{1,count};
%         u{2,sz} = SO{2,count};
%     end
% end
