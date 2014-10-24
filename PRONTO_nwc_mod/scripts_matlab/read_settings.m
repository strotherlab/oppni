function read_settings

global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART1_afni_steps.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end
if exist([cd '/SETTINGS.txt'],'file')
    File = fopen([cd '/SETTINGS.txt'],'r');
else
    File = fopen([CODE_PATH '/SETTINGS.txt'],'r');
end
if File~=-1
    tline = fgetl(File);
    
    while(ischar(tline))
        
        index_afni = strfind(tline,'AFNI_PATH');
        if ~isempty(index_afni)
            ine = strfind(tline,'=');
            AFNI_PATH = tline(ine+1:end);
        end
        
        index_fsl  = strfind(tline,'FSL_PATH');
        if ~isempty(index_fsl)
            ine = strfind(tline,'=');
            FSL_PATH = tline(ine+1:end);
        end
        tline = fgetl(File);
    end
end