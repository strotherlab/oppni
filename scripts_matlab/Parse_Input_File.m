function [Input_nifti_file_path,Input_nifti_file_prefix,Output_nifti_file_path,Output_nifti_file_prefix,split_info_file,STRUCT_File,PHYstr,NOISE_ROI,DROP_first,DROP_last] = Parse_Input_File(tline)

global CODE_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Parse_Input_File.m'));
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
%%Input_nifti_file_temp = strrep(Input_nifti_file_temp,'.nii','');
[Input_nifti_file_path_temp,Input_nifti_file_prefix,infile_ext] = fileparts(Input_nifti_file_temp);

if(isempty(strfind(Input_nifti_file_prefix,',')))
     Input_nifti_file_prefix = cellstr(Input_nifti_file_prefix);
else Input_nifti_file_prefix = regexp(Input_nifti_file_prefix,',','split');
end

N_run = length(Input_nifti_file_prefix);


% parse output directory
ifile = strfind( upper(tline), 'OUT=' ); ifile = ifile+4;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>ifile);
[Output_nifti_file_path_temp,Output_nifti_file_prefix,ext] = fileparts(tline(ifile:ips(1)));

if(isempty(strfind(Output_nifti_file_prefix,',')))
     Output_nifti_file_prefix = cellstr(Output_nifti_file_prefix);
else Output_nifti_file_prefix = regexp(Output_nifti_file_prefix,',','split');
end

if length(Output_nifti_file_prefix)~=N_run
    display(sprintf('Error the number of output prefixes does not match with the number of inputs, please check the line %s',tline));
    sge_exit(100);
end


% load in "split_info" structure with information about task onset and
% analysis model parameters
itask = strfind(  upper(tline), 'TASK=' );
if ~isempty(itask)
    itask = itask+5;
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>itask);
    split_info_file_temp = tline(itask:ips(1));
    % file parts
    [split_info_file_path_temp,split_info_file_prefix,split_info_file_ext] = fileparts(split_info_file_temp);
    
    if(isempty(strfind(split_info_file_prefix,',')))
        split_info_file_prefix = cellstr(split_info_file_prefix);
    else split_info_file_prefix = regexp(split_info_file_prefix,',','split');
    end
    
    if length(split_info_file_prefix)~=N_run
        display(sprintf('Error the number of split info files does not match with the number of inputs, please check the line %s',tline));
        sge_exit(100);
    end
else
    split_info_file_path_temp = [];
    split_info_file_prefix = cell(1,N_run);
end

istruct = strfind(  upper(tline), 'STRUCT=' ); istruct = istruct+7;
ips   = [strfind( tline, ' ' )-1 length(tline)];
ips   = ips(ips>istruct);
STRUCT_File_temp = tline(istruct:ips(1));


isphysio = strfind(  upper(tline), 'PHYSIO=' ); isphysio = isphysio+7;
if ~isempty(isphysio)
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>isphysio);
    [PHYstr_path_temp,PHYstr_prefix,ext] = fileparts(tline(isphysio:ips(1)));
    
    if(isempty(strfind(PHYstr_prefix,',')))
         PHYstr_prefix = cellstr(PHYstr_prefix);
    else PHYstr_prefix = regexp(PHYstr_prefix,',','split');
    end
    
    if length(PHYstr_prefix)~=N_run
        display(sprintf('Error the number of Physio data files does not match with the number of inputs, please check the line %s',tline));
        sge_exit(100);
    end
   
else
    PHYstr_prefix = cell(N_run,1);
    PHYstr_path_temp = [];
end

isdrop = strfind(  upper(tline), 'DROP=[' );
if ~isempty(isdrop)
    isdrop = isdrop+6;
    xline=tline(isdrop:end);
    ips = strfind(xline,',');
    DROP_first_temp = str2num(xline(1:ips(1)-1));
    ips2 = strfind(xline,']');
    DROP_last_temp  = str2num(xline(ips(1)+1:ips2(1)-1));
else
    DROP_first_temp = 0;
    DROP_last_temp = 0;
end

% parse output directory
ifile = strfind(  upper(tline), 'CUSTOMREG=' ); ifile = ifile+10;
if ~isempty(ifile)
    ips   = [strfind( tline, ' ' )-1 length(tline)];
    ips   = ips(ips>ifile);
    noise_roi_file = tline(ifile:ips(1));
    %LMP 13/09/2021
    [noise_roi_path_temp,noise_roi_file_prefix,noise_roi_ext] = fileparts(noise_roi_file);
    %[noise_roi_path_temp,noise_roi_file_prefix,ext] = fileparts(noise_roi_file);
    
    if(isempty(strfind(noise_roi_file_prefix,',')))
         noise_roi_file_prefix = cellstr(noise_roi_file_prefix);
    else 
         noise_roi_file_prefix = regexp(noise_roi_file_prefix,',','split');
    end
    
    if length(noise_roi_file_prefix)~=N_run
        display(sprintf('Error the number of NOISE ROI files does not match with the number of inputs, please check the line %s',tline));
        sge_exit(100);
    end
else
    noise_roi_file_prefix = cell(N_run,1);
    noise_roi_path_temp   = [];
end

for k = 1:N_run
    Input_nifti_file_prefix{k} = [Input_nifti_file_prefix{k}, infile_ext]; %% add in ext
    Input_nifti_file_path{k} = Input_nifti_file_path_temp;
    Output_nifti_file_path{k} = Output_nifti_file_path_temp;
    STRUCT_File{k} = STRUCT_File_temp;
    if ~isempty(PHYstr_path_temp)
        PHYstr{k} = [PHYstr_path_temp '/' PHYstr_prefix{k}];
    else
        PHYstr{k} = [];
    end
    DROP_first{k} = DROP_first_temp;
    DROP_last{k}  = DROP_last_temp;
    if ~isempty(split_info_file_prefix{k})
        % re-generate split_info path+filename        
        if(isempty(split_info_file_path_temp))
            split_info_file{k} = [split_info_file_prefix{k} split_info_file_ext];
        else  
            split_info_file{k} = [split_info_file_path_temp '/' split_info_file_prefix{k} split_info_file_ext];
        end
    else
        error('split_info file missing for one of the runs');
    end

    if ~isempty(noise_roi_path_temp)
        %LMP 13/09/2021
        NOISE_ROI{k} = [noise_roi_path_temp '/' noise_roi_file_prefix{k} noise_roi_ext];
        %NOISE_ROI{k} = [noise_roi_path_temp '/' noise_roi_file_prefix{k} '.nii'];
    else
        NOISE_ROI{k} = noise_roi_file_prefix{k};
    end
    disp(sprintf("LMP-DEBUG: NOISE_ROI{k} = %s", NOISE_ROI{k}));

end
    
    
    
