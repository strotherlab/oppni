function mkdir_r(pathstr)
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

ind = strfind(pathstr,'/');
for i = 1:length(ind)
    current_path = pathstr(1:ind(i)-1);
    if ~isempty(current_path)
        if ~exist(current_path,'dir')
            mkdir(current_path);
        end
    end
end
if ~exist(pathstr,'dir')
    mkdir(pathstr);
end

