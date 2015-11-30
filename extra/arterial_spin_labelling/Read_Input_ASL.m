function [InputStruct] = Read_Input_ASL(inputfile)

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
    return;
end
% read in first line
tline               = fgetl(fid);
if ~ischar(tline)
    InputStruct = [];
    return;
end

ksub     = 0;

%% read input file
while ischar(tline) 

    ksub = ksub  +1;
    
    [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix,STRUCT_File] = Parse_Input_ASL(tline);
    InputStruct(ksub).run(1).Input_nifti_file_path             = Input_nifti_file_path;
    InputStruct(ksub).run(1).Input_nifti_file_prefix           = Input_nifti_file_prefix;
    InputStruct(ksub).run(1).Output_nifti_file_path            = Output_nifti_file_path;
    InputStruct(ksub).run(1).Output_nifti_file_prefix          = Output_nifti_file_prefix;
    InputStruct(ksub).run(1).STRUCT_File                       = STRUCT_File;
 
    tline = fgetl(fid);
    if isempty(tline)
        tline = fgetl(fid);
    end
end
fclose(fid); 

