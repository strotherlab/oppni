function [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix] = Parse_Input_noddi(tline)

global CODE_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Parse_Input_struct_vbm.m'));
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
Input_nifti_file_prefix = cellstr(Input_nifti_file_prefix);

% parse output directory
ifile = strfind( upper(tline), 'OUT=' ); ifile = ifile+4;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
[Output_nifti_file_path,Output_nifti_file_prefix,ext] = fileparts(tline(ifile:ips(1)));
Output_nifti_file_prefix = cellstr(Output_nifti_file_prefix);

