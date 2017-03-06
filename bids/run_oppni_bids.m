function run_oppni_bids(proc,varargin)

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

%% some early value checks
if( nargin < 9 )
    error('too few input arguments (1=IN,2=OUT,3=STRUCT,4=PHYSIO,5=DROP1,6=DROP2,7=TASK,8=EVENTS,9=ATLAS)');
end
if( isempty(varargin{1}) || ~exist(varargin{1},'file') ) error(strcat('input file: ',varargin{1},' does not exist')); end
if( isempty(varargin{3}) || ~exist(varargin{3},'file') ) error(strcat('structural file: ',varargin{3},' does not exist')); end
if( isempty(varargin{4}) )
    warning(strcat('physio file not specified for: ',varargin{1},'. Cannot do RETROICOR'));
    varargin{4}='None';
elseif(~exist(varargin{4},'file') && (~exist([varargin{4},'.resp.1D'],'file') && ~exist([varargin{4},'.card.1D'],'file'))  ) 
    warning(strcat('physio file: ',varargin{4},' does not exist. Cannot do RETROICOR'));
    varargin{4}='None';
end
if( isempty(varargin{5}) || (isnumeric(varargin{5}) && varargin{5}==0) ) 
    varargin{5}=0;
    warning(strcat('No scans dropped from start of: ',varargin{1},'.Is this OK?'));
elseif( all(ismember(varargin{5},'0123456789 ')) )
    varargin{5} = str2num(varargin{5});
elseif( ~isnumeric(varargin{5}))
    error(strcat('Invalid DROP(1) value: ',varargin{5}))
end
if( isempty(varargin{6}) || (isnumeric(varargin{6}) && varargin{6}==0) ) 
    varargin{6}=0;
    warning(strcat('No scans dropped from end of: ',varargin{1},'.Is this OK?'));
elseif( all(ismember(varargin{6},'0123456789 ')) )
    varargin{6} = str2num(varargin{6});
elseif( ~isnumeric(varargin{6}))
    error(strcat('Invalid DROP(2) value: ',varargin{6}))
end
if( isempty(varargin{7}) || ~exist(varargin{7},'file') ) error(strcat('Task .json file: ',varargin{7},' does not exist')); end
if( isempty(varargin{8}) || ~exist(varargin{8},'file') ) error(strcat('Event .tsv file: ',varargin{8},' does not exist')); end

%% now setting parameter defaults

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
mkdir( fullfile( outpath, 'task_files') );
mkdir( fullfile( outpath, 'input_files') );

% populate pipeline file
newpipefile = fullfile( outpath, 'pipeline_combinations.txt');
if ~exist(newpipefile,'file')
    make_pipeline_file( newpipefile );
end
% create task file
newtaskfile = fullfile( outpath, 'task_files', strcat(prefix,'.txt') );
bids_to_oppni_task(varargin{7}, varargin{8}, task_type, newtaskfile);
% create input file
newinputfile = fullfile( outpath, 'input_files', strcat(prefix,'.txt') );
make_input_file( newinputfile, varargin(1:6), newtaskfile );

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
