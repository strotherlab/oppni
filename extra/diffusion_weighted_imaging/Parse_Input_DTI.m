function [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix,STRUCT_File] = Parse_Input_DTI(tline)

global CODE_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Parse_Input_DTI.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
        addpath(CODE_PATH);
        addpath([CODE_PATH '/NIFTI_tools'])
    end
end

% parse output directory
ifile = strfind( upper(tline), 'IN=' ); 
ifile = ifile+3;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
Input_nifti_file_temp = tline(ifile:ips(1));
Input_nifti_file_temp = strrep(Input_nifti_file_temp,'.nii','');
[Input_nifti_file_path,Input_nifti_file_prefix,ext] = fileparts(Input_nifti_file_temp);
% split in DWI runs
if(isempty(strfind(Input_nifti_file_prefix,',')))
     Input_nifti_file_prefix = cellstr(Input_nifti_file_prefix);
else Input_nifti_file_prefix = regexp(Input_nifti_file_prefix,',','split');
end

if length(Input_nifti_file_prefix)>2
    display(sprintf('Error you can only specify 1 or 2 diffusion-weighted input files, please check the line %s',tline));
    sge_exit(100);
end

% parse output directory
ifile = strfind( upper(tline), 'OUT=' ); ifile = ifile+4;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
[Output_nifti_file_path,Output_nifti_file_prefix,ext] = fileparts(tline(ifile:ips(1)));

if(~isempty(strfind(Output_nifti_file_prefix,',')))
    display(sprintf('Error you can only specify 1 output destination, please check the line %s',tline));
    sge_exit(100);
end

% parse structural file
istruct = strfind(  upper(tline), 'STRUCT=' ); istruct = istruct+7;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>istruct);
STRUCT_File = tline(istruct:ips(1));

    
    
    
