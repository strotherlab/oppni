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

% initialize, number of subjects / multi-run?
ksub     = 0;
MULTI_RUN_INPUTFILE = false;

%% read each line of input file, parse into relevant structure
while ischar(tline) 

    ksub = ksub  +1;
    % for each run in line, store file names / locations

%%------------------- PARSE INPUT FILE LINES INTO FIELDS ----------------%%
%%-----------------------------------------------------------------------%%

        % --> the OUT field is the only one with both path and prefix
        fieldlist = {'IN','OUT','STRUCT','TASK','PHYSIO','CUSTOMREG'}; %% list of fields - IN should always be first!!!
        tline     = regexprep(tline, '\t', ' ');
        ispaces   = [strfind( tline, ' ' )-1 length(tline)];           %% index all spaces (-1), and eol
        temp_path_struc = [];                                          %% initialize structure for storing paths
        
        for(i=1:length(fieldlist)) %% go through each field

            if(~isempty(strfind( upper(tline),[fieldlist{i},'='] ))) %% check if exists
            
                % find field, index just past '=' sign // first space after this field
                ifile = strfind( upper(tline),[fieldlist{i},'='] ) + (numel(fieldlist{i})+1);
                ips   = ispaces(ispaces>ifile);
                % take substring 
                filestring_temp = tline(ifile:ips(1));
                % split csvs --> convert to cell array
                if(isempty(strfind(filestring_temp,','))) 
                     filestring_temp = {filestring_temp};
                else filestring_temp = regexp(filestring_temp,',','split'); 
                end

                % go through runs and store
                for(j=1:length(filestring_temp))
                    %
                    if( strcmpi(fieldlist{i},'OUT') ) %--------------------
                        
                        % store path and prefix separately
                        [temp_path_struc.([fieldlist{i},'_path']){j},temp_path_struc.([fieldlist{i},'_prefix']){j},ext] = fileparts(filestring_temp{j});
                        % if only entry #1 has a path, assume other files have same
                        if( j>1 && isempty(temp_path_struc.([fieldlist{i},'_path']){j}) ) 
                            temp_path_struc.([fieldlist{i},'_path']){j} = temp_path_struc.([fieldlist{i},'_path']){1}; 
                        end                        
                        
                    else %-------------------------------------------------
                    
                        temp_path_struc.([fieldlist{i},'_filename']){j} = filestring_temp{j};
                        % store path and prefix temporarily                        
                        [temppath{j},tempprefix1,ext1] = fileparts(filestring_temp{j});
                        % if multi-run and only entry #1 has a path, assume other files have same
                        if( j>1 && isempty(temppath{j}) ) 
                            temp_path_struc.([fieldlist{i},'_filename']){j} = [temppath{1}, '/', temp_path_struc.([fieldlist{i},'_filename']){j}]; 
                        end
                    end
                end

                if(i==1) %% initialize number of runs
                    N_run = length(filestring_temp); 
                elseif(length(filestring_temp) ~= N_run ) %% check all subsequent fields
                    display(sprintf('Error, number of runs in %s= does not match IN=, please check the line %s',fieldlist{i},tline));
                    sge_exit(100);            
                end
            else
                % if field not specified, this is only OK for PHYSIO and CUSTOMREG
                if(  strcmpi( fieldlist{i}, 'PHYSIO' ) || strcmpi( fieldlist{i}, 'CUSTOMREG' )  )
                    display(sprintf('Warning: field %s not specified, make sure this is correct for line %s',fieldlist{i},tline));
                    temp_path_struc.([fieldlist{i},'_filename']) = [];
                else
                    display(sprintf('Error mandatory field %s not specified, please check the line %s',fieldlist{i},tline));
                    sge_exit(100);                    
                end
            end
        end
        
        %% number of volumes to drop from start/end of run
        isdrop = strfind(  upper(tline), 'DROP=[' );
        if ~isempty(isdrop)
            isdrop = isdrop+6;
            xline=tline(isdrop:end);
            ips = strfind(xline,',');
            DROP_first_temp = str2num(xline(1:ips(1)-1));
            ips2 = strfind(xline,']');
            DROP_last_temp  = str2num(xline(ips(1)+1:ips2(1)-1));
        else
            disp('Warning DROP field not specified -- assumes you want to keep all run volumes');
            DROP_first_temp = 0;
            DROP_last_temp = 0;
        end

%%------------------- DONE PARSING LINES, NOW STORE INFO ----------------%%
%%-----------------------------------------------------------------------%%
        
        % for each run, format and store file names / locations
        for krun = 1:N_run

            InputStruct(ksub).run(krun).Output_nifti_file_prefix  = temp_path_struc.OUT_prefix{krun};
            InputStruct(ksub).run(krun).Output_nifti_file_path    = temp_path_struc.OUT_path{krun};            
            InputStruct(ksub).run(krun).Input_nifti_filename      = temp_path_struc.IN_filename{krun};
            InputStruct(ksub).run(krun).STRUCT_File               = temp_path_struc.STRUCT_filename{krun};
            InputStruct(ksub).run(krun).split_info_file           = temp_path_struc.TASK_filename{krun};
            if ~isempty(temp_path_struc.PHYSIO_filename)
                  InputStruct(ksub).run(krun).PHYstr              =  temp_path_struc.PHYSIO_filename{krun};
            else  InputStruct(ksub).run(krun).PHYstr              = [];
            end
            if ~isempty(temp_path_struc.CUSTOMREG_filename)
                  InputStruct(ksub).run(krun).Noise_ROI           = temp_path_struc.CUSTOMREG_filename{krun};
            else  InputStruct(ksub).run(krun).Noise_ROI           = [];
            end            
            InputStruct(ksub).run(krun).DROP_first                = DROP_first_temp;
            InputStruct(ksub).run(krun).DROP_last                 = DROP_last_temp;
        end
        % extra fields --> for single output, produced from multiple runs
        InputStruct(ksub).run(1).subjectprefix  = ['_' temp_path_struc.OUT_prefix{1}];
        InputStruct(ksub).run(1).subjectmask    = [temp_path_struc.OUT_path{1}  '/intermediate_processed/masks/' temp_path_struc.OUT_prefix{1} '_mask.nii'];
        
        % store information: multiple runs?
        if krun>1
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
