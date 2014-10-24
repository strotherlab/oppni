function extract_subjects(InputStruct,OutputStruct_FileName)


global NUMBER_OF_CORES
NUMBER_OF_CORES = str2double(getenv('PIPELINE_NUMBER_OF_CORES'));
if isnan(NUMBER_OF_CORES)
    NUMBER_OF_CORES = 1;
end
display(sprintf('The number of cores used by the code=%d',NUMBER_OF_CORES));
if ~exist('OCTAVE_VERSION','builtin')
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

if ~isstruct(InputStruct)
    [InputStruct,NEW_INPUTFILE_VERSION] = Read_Input_File(InputStruct,'subject*run');
end
if size(InputStruct,2)~=1
    NEW_INPUTFILE_VERSION = true;
end

count = 0;
for k = 1:numel(InputStruct)
    if ~isempty(InputStruct(k).STRUCT_File)
        count = count  +1;
        TempStruct(count) = InputStruct(k);
    end
end

count = 0;
for k = 1:numel(TempStruct)
    Input_STRUCT_Name{count}     = TempStruct(k).STRUCT_File;    
end

[uInput_STRUCT_Name, uInd] = unique(Input_STRUCT_Name);
InputStruct = TempStruct(uInd);

[pathstr,name,ext] = fileparts(OutputStruct_FileName);
mkdir_r(pathstr);
save(OutputStruct_FileName,'InputStruct','-v7');


