function pronto(proc,varargin)

if strcmpi(proc,'PART1')
    Pipeline_PART1(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7},varargin{8},varargin{9});
end
if strcmpi(proc,'PART2')
    Pipeline_PART2(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
end

if strcmpi(proc,'SPNORM')  
    spatial_normalization(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
end

if strcmpi(proc,'GMASK') 
    if length(varargin)==2
        group_mask_tissue_maps(varargin{1},varargin{2});
    else
        group_mask_tissue_maps(varargin{1},'');
    end
end
