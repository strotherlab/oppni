function addpath_oppni(folder_path)
% path management helper

if ~isdeployed
    addpath(folder_path);
end
