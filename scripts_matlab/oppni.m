function oppni(proc,varargin)

% 1 IN path
% 2 OUT path
% 3 STRUCT path
% 4 PHYSIO path
% 5 DROP value1
% 6 DROP value2
% 7 TASK JSON (BIDS)
% 8 EVENTS TSV

% setting defaults
task_type = 'BLOCK';
analysis_model = 'LDA';
modelparam=[];
contrast_list_str=[];
dospnormfirst=0;
DEOBLIQUE=0;
TPATTERN='auto_hdr';
TOFWHM=0;
niiout=0;
optimize_metric='dPR';
mot_gs_control=[1 0];
process_out=1;
keepmean=0;
whichpipes='ALL';
input_voxelsize=[];
flag_step=0;
reference_file=[];

% populate pipeline file
newpipefile = [varargin{2},'_pipeline_file.txt'];
make_pipeline_file( newpipefile );
% create task file
newtaskfile = [varargin{2},'_task_file.txt'];
bids_to_oppni_task(varargin{7}, varargin{8}, task_type, newtaskfile);
% create input file
newinputfile = [varargin{2},'_input_file.txt'];
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

    QC_wrapper(flag_step,newinputfile, [], 2);

else
    error('Unrecognized part name: must be one of PART1, PART2, SPNORM, GMASK, QC1 and QC2.');
end
