function run_oppni_bids( in_dir, out_dir, level, varargin )
%
% RUN_OPPNI_BIDS: wrapper script to execute oppni jobs in BIDS framework.
%
%  Syntax:
%
%     run_oppni_bids( in_dir, out_dir, level, varargin )
%

% number of supplementary arguments
numvars = round(numel(varargin)/2);
% ensure there is an even number of supplemental arguments
if( mod(numel(varargin),2) ~=0 ) error('supplementary arguments must come in name/value pairs'); end
% ensure there arent too many supplementayl args (7 max)
if(numvars>7) error('too many arguments. full list: {--run_name, --task_design, --participant, --contrast, --analysis_model, --ndrop, --atlasfile}'); end
% ensure there are minimum number (at least 2 required
if(numvars<2) error('not enough supplementary arguments. Need to at least specify {--run_name, --task_design}'); end

% start with empty supplemental fields
run_name       =[];
task_design    =[];
participant    =[];
contrast       =[];
analysis_model =[];
ndrop          =[];
atlasfile      =[];

for(i=1:numvars) %% go through each supplemental field
    
   valpairs{1} = varargin{ (i-1)*2+1 }; 
   valpairs{2} = varargin{ (i-1)*2+2 };
   
   if( isempty(strfind(valpairs{1},'--'))) error(['out of order. expected ', valpairs{1} ' to be an argument name']);  end
   if(~isempty(strfind(valpairs{2},'--'))) error(['out of order. expected ', valpairs{2} ' to be an argument value']); end
   
   % check mandatory fields
   if( strcmpi(valpairs{1},'--run_name')       ) run_name       = valpairs{2}; end
   if( strcmpi(valpairs{1},'--task_design')    ) task_design    = valpairs{2}; end
   % check optional fields
   if( strcmpi(valpairs{1},'--participant')    ) participant    = valpairs{2}; end
   if( strcmpi(valpairs{1},'--contrast')       ) contrast       = valpairs{2}; end
   if( strcmpi(valpairs{1},'--analysis_model') ) analysis_model = valpairs{2}; end
   if( strcmpi(valpairs{1},'--ndrop')          ) ndrop          = valpairs{2}; end
   if( strcmpi(valpairs{1},'--atlasfile')      ) atlasfile      = valpairs{2}; end
end

%% checks on settings

if(isempty(run_name)   ) error('need to specify a --run_name argument'); end
if(isempty(task_design)) error('need to specify a --task_design argument'); end

if(isempty(participant))
    disp('WARNING: all subjects are submitted as single job');
end
if(isempty(contrast))
    disp('WARNING: no contrast is specified. By default, will run all task conditions vs. baseline');
end
if(isempty(analysis_model))
    disp('WARNING: no analysis specified. Using default model selection');
    if    (strcmpi(task_design,'block'))      analysis_model='LDA';
    elseif(strcmpi(task_design,'event'))      analysis_model='erCVA';
    elseif(strcmpi(task_design,'nocontrast')) analysis_model='falff';
    else error('the --task_design value is not valid. Needs to be either "block", "event" or "nocontrast"');
    end
end
if(isempty(ndrop)) 
    disp('WARNING: the --ndrop field is empty. The default is setting ndrop=0 (no non-equilibrium scans dropped)');
    ndrop=0; 
end
if(isempty(atlasfile))
    disp('WARNING: the --atlasfile field is empty. No spatial normalization is possible');
end

%% run files
bids_parsejobs( in_dir, out_dir, level, participant, run_name, contrast, task_design, analysis_model, ndrop, atlasfile );


