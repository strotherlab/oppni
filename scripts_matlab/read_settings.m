function read_settings

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

global CODE_PATH AFNI_PATH FSL_PATH
if isempty(CODE_PATH)
    CODE_PATH = fileparts(which('Pipeline_PART1_afni_steps.m'));
    if CODE_PATH(end)~='/'
        CODE_PATH = [CODE_PATH '/'];
    end
end

% if exist([cd '/SETTINGS.txt'],'file')
%     File = fopen([cd '/SETTINGS.txt'],'r');
% else
%     File = fopen([CODE_PATH '/SETTINGS.txt'],'r');
% end
% if File~=-1
%     tline = fgetl(File);
%     
%     while(ischar(tline))
%         
%         index_afni = strfind(tline,'AFNI_PATH');
%         if ~isempty(index_afni)
%             ine = strfind(tline,'=');
%             AFNI_PATH = tline(ine+1:end);
%             if ~isempty(AFNI_PATH)
%                 if AFNI_PATH(end)~='/'
%                     AFNI_PATH = [AFNI_PATH '/'];
%                 end
%             end
%         end
%         
%         index_fsl  = strfind(tline,'FSL_PATH');
%         if ~isempty(index_fsl)
%             ine = strfind(tline,'=');
%             FSL_PATH = tline(ine+1:end);
%             if ~isempty(FSL_PATH)
%                 if FSL_PATH(end)~='/'
%                     FSL_PATH = [FSL_PATH '/'];
%                 end
%             end
%         end
%         tline = fgetl(File);
%     end
%     fclose(File);
% end

AFNI_PATH = getenv('AFNI_PATH');
FSL_PATH  = getenv('FSL_PATH');
