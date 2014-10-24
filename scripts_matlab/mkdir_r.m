function mkdir_r(pathstr)
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

