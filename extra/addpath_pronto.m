function addpath_pronto(folder_path)

if ~isdeployed
    addpath(folder_path);
end
