function [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix,STRUCT_File] = Parse_Input_DOALL(tline)

global CODE_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Parse_Input_ASL.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
        addpath(CODE_PATH);
        addpath([CODE_PATH '/NIFTI_tools'])
    end
end

% parse output directory
ifile = strfind( upper(tline), 'IN=' ); ifile = ifile+3;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
Input_nifti_file_temp = tline(ifile:ips(1));
Input_nifti_file_temp = strrep(Input_nifti_file_temp,'.nii','');
[Input_nifti_file_path,Input_nifti_file_prefix,ext] = fileparts(Input_nifti_file_temp);
% split into runs
if(isempty(strfind(Input_nifti_file_prefix,',')))
     Input_nifti_file_prefix = cellstr(Input_nifti_file_prefix);
else Input_nifti_file_prefix = regexp(Input_nifti_file_prefix,',','split');
end

% parse output directory
ifile = strfind( upper(tline), 'OUT=' ); ifile = ifile+4;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
[Output_nifti_file_path,Output_nifti_file_prefix,ext] = fileparts(tline(ifile:ips(1)));
% split into runs
if(isempty(strfind(Output_nifti_file_prefix,',')))
     Output_nifti_file_prefix = cellstr(Output_nifti_file_prefix);
else Output_nifti_file_prefix = regexp(Output_nifti_file_prefix,',','split');
end

if( length(Output_nifti_file_prefix) ~= length(Input_nifti_file_prefix) )
    error('number of input prefixes should match number of output prefixes');
end

% parse structural file
istruct = strfind(  upper(tline), 'STRUCT=' ); 
if(isempty(istruct))
    STRUCT_File=[];
else
    istruct = istruct+7;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>istruct);
    STRUCT_File = tline(istruct:ips(1));
end