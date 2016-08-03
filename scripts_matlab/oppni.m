function oppni(proc,varargin)

% varargins:
%
% 1 IN path
% 2 OUT path
% 3 STRUCT path
% 4 PHYSIO path
% 5 DROP value1
% 6 DROP value2
% 7 TASK JSON (BIDS)
% 8 EVENTS TSV
% 9 ATLAS PATH

% setting defaults, part1
task_type         = 'BLOCK';
analysis_model    = 'LDA';
modelparam        =[];
contrast_list_str =[];
dospnormfirst     =0;
DEOBLIQUE         =0;
TPATTERN          ='auto_hdr';
TOFWHM            =0;
niiout            =0;
% setting defaults, part2
optimize_metric   ='dPR';
mot_gs_control    =[1 0];
process_out       =1;
keepmean          =0;
whichpipes        ='ALL';
% setting defaults, spatnorm
input_voxelsize   =[];
flag_step         =0;

if( nargin >9 )
    reference_file=varargin{9};
else
    reference_file=[];
end

[outpath,prefix,~] = fileparts( varargin{2} );
mkdir_r(outpath);
% populate pipeline file
newpipefile = fullfile( outpath, '/pipeline_file.txt');
if ~exist(newpipefile,'file')
    make_pipeline_file( newpipefile );
end
% create task file
newtaskfile = fullfile( outpath, strcat(prefix,'_task_file.txt'));
bids_to_oppni_task(varargin{7}, varargin{8}, task_type, newtaskfile);
% create input file
newinputfile = fullfile( outpath, strcat(prefix,'_input_file.txt') );
make_input_file( newinputfile, varargin(1:6) );

if strcmpi(proc,'PART1')

    Pipeline_PART1(newinputfile,newpipefile, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE,  TPATTERN,   TOFWHM);

elseif strcmpi(proc,'PART2')

	Pipeline_PART2(newinputfile, optimize_metric, mot_gs_control, process_out, keepmean,   whichpipes)
    
elseif strcmpi(proc,'SPNORM')

    spatial_normalization(newinputfile,reference_file,input_voxelsize,flag_step,DEOBLIQUE)
    
elseif strcmpi(proc,'GMASK') 

    group_mask_tissue_maps(newinputfile,'');
    
elseif strcmpi(proc,'QC1') || strcmpi(proc,'QC2')

    group_mask_tissue_maps(newinputfile,'');
    QC_wrapper(flag_step,newinputfile, [], 2);
else
    error('Unrecognized part name: must be one of PART1, PART2, SPNORM, GMASK, QC1 and QC2.');
end
