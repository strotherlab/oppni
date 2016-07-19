function pronto(proc,varargin)

% for ii = 1 : length(varargin)
%     if ~ischar(varargin{ii})
%         str = num2str(varargin{ii});
%     else
%         str = varargin{ii};
%     end
%     fprintf('var arg %d : %s\n',ii, str)
% end

if strcmpi(proc,'PART1')
    if nargin < 11
        error('Insufficient number of arguments for Part 1 - must supply:\n InputStruct,input_pipeset, analysis_model, modelparam, niiout,     contrast_list_str, dospnormfirst, DEOBLIQUE,  TPATTERN,   TOFWHM');
    end
    Pipeline_PART1(varargin{1},varargin{2},  varargin{3},   varargin{4},  varargin{5},varargin{6},      varargin{7},    varargin{8},varargin{9},varargin{10});
    % Pipeline_PART1(InputStruct,input_pipeset, analysis_model, modelparam, niiout,     contrast_list_str, dospnormfirst, DEOBLIQUE,  TPATTERN,   TOFWHM)
elseif strcmpi(proc,'PART2')
    if nargin < 7
         error('Insufficient number of arguments for Part 2 - must supply:\n InputStruct, optimize_metric, mot_gs_control, process_out, keepmean,   whichpipes');
    end
	% Pipeline_PART2(InputStruct, optimize_metric, mot_gs_control, process_out, keepmean,   whichpipes)
    Pipeline_PART2(varargin{1},  varargin{2},		varargin{3},	varargin{4},varargin{5},varargin{6});
    
elseif strcmpi(proc,'SPNORM')
    if nargin < 6
         error('Insufficient number of arguments for SPNORM - must supply:\n InputStruct,reference_file,input_voxelsize,flag_step,DEOBLIQUE');
    end
    % spatial_normalization(InputStruct,reference_file,input_voxelsize,flag_step,DEOBLIQUE)
    spatial_normalization(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
    
elseif strcmpi(proc,'GMASK') 
    if length(varargin)==2
        group_mask_tissue_maps(varargin{1},varargin{2});
    else
        group_mask_tissue_maps(varargin{1},'');
    end
    
elseif strcmpi(proc,'QC1') || strcmpi(proc,'QC2')
    % first argument coming from compiled code is an integer indicating step
    QC_wrapper(varargin{1}, varargin{2},varargin{3}, varargin{4});

% elseif strcmpi(proc,'QC0')
%     QC_wrapper(0, varargin{1},varargin{2},varargin{3}); 

else
    error('Unrecognized part name: must be one of PART1, PART2, SPNORM, GMASK, QC1 and QC2.');
end


