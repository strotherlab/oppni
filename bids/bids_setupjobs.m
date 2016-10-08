function bids_setupjobs(proc,outpath,varargin)
%
% BIDS_SETUPJOBS: takes in formatted list of bids files and prepares
% requisite oppni files to running the pipelines
%
%  Syntax:
%
%     bids_setupjobs(proc,outpath,varargin)
%

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
% 10 contrast
% 11 task design
% 12 analysis model

%% now setting parameter defaults

% setting defaults, part1
task_type         = varargin{11};%'BLOCK';
analysis_model    = varargin{12};%'LDA';
modelparam        =[];
contrast_list_str = varargin{10};%[];
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
reference_file=varargin{9};
input_voxelsize   =[2 2 2]; %% temp
flag_step         =0;

mkdir_r(outpath);
mkdir( fullfile( outpath, 'task_files') );
mkdir( fullfile( outpath, 'input_files') );

% populate pipeline file
newpipefile = fullfile( outpath, 'pipeline_combinations.txt');
if ~exist(newpipefile,'file')
    make_pipeline_file( newpipefile );
end

nsub=numel( varargin{2}    );
nrun=numel( varargin{2}{1} );

% create task file
for(i=1:nsub)
for(j=1:nrun)
    [path,prefix,ext] = fileparts(varargin{2}{i}{j});
    newtaskfile{i}{j} = fullfile( outpath, 'task_files', strcat(prefix,'.txt') ); %% list of taskfile names
end
end
bids_to_oppni_task(varargin{7}, varargin{8}, task_type, newtaskfile); %% populate full set of taskfiles

% create "all subjects" input file
if(nsub>1)
    newinputfile = fullfile( outpath, 'input_files', 'all_sub.txt' ); %% change to either "all_sub" or "single-sub"
    make_input_file( newinputfile, varargin(1:6), newtaskfile ); %% modify to generate multi-line input file
elseif(nsub==1)
    % make "single-subject" input files
    ix=strfind(varargin{2}{1}{1},'_'); ix=ix( ix>length(outpath)+2 );
    prefix = varargin{2}{1}{1}(length(outpath)+2:ix(1)-1);
    newinputfile = fullfile( outpath, 'input_files', strcat(prefix,'.txt') ); %% change to either "all_sub" or "single-sub"
    make_input_file( newinputfile, varargin(1:6), newtaskfile ); %% modify to generate multi-line input file
end


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
