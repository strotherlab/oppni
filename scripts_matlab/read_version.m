function read_version

global CODE_PATH CODE_PROPERTY
if isempty(CODE_PROPERTY)
    if isempty(CODE_PATH)
        CODE_PATH = fileparts(which('read_version.m'));
        if CODE_PATH(end)~='/'
            CODE_PATH = [CODE_PATH '/'];
            addpath(CODE_PATH);
            addpath([CODE_PATH '/NIFTI_tools'])
        end
    end
    cp = CODE_PATH;
    if cp(end)=='/' || cp(end)=='\'
        cp(end)=[];
    end
    [pathstr,name] = fileparts(cp);
    try
        fid = fopen([pathstr '/_documentation/UPDATE_HISTORY.txt'],'r');
        for i = 1:4
            tline               = fgetl(fid);
            eval(tline);
        end
        CODE_PROPERTY.VERSION         = VERSION;
        CODE_PROPERTY.DATE            = DATE;
        CODE_PROPERTY.REVISION        = REVISION;
        CODE_PROPERTY.DESCRIPTION     = DESCRIPTION;
        str = ['PRONTO ' sprintf('%.2f',CODE_PROPERTY.VERSION) 'rev' CODE_PROPERTY.REVISION(12:end-2) '-Date:' CODE_PROPERTY.DATE(7:end-2)];
        CODE_PROPERTY.NII_HEADER = str;
    catch
        
        CODE_PROPERTY.VERSION         = 'unknown';
        CODE_PROPERTY.DATE            = 'unknown';
        CODE_PROPERTY.REVISION        = 'unknown';
        CODE_PROPERTY.DESCRIPTION     = 'unknown';
        CODE_PROPERTY.NII_HEADER = 'PRONTO unknwon version';
    end
end
