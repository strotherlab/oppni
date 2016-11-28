function run_oppni_local( subject_inputs, pipelines, analysis_model, TR_MSEC, contrast_list, TEMPLATE_VOL, VOXDIMS, DEOBLIQUE, TPATTERN, TOFWHM, opt_metric, WARP_TYPE )
%
% RUN_OPPNI_LOCAL: script for running pipelines directly in matlab, for
% instance where SGE system is not available.
%
%  Syntax:
%
%     run_oppni_local( subject_inputs, pipelines, analysis_model, TR_MSEC, contrast_list, TEMPLATE_VOL, VOXDIMS, DEOBLIQUE, TPATTERN )
%
%

CODE_VERSION = '$Revision: 158 $';
CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';

addpath scripts_matlab;
addpath scripts_matlab/optimization;
addpath scripts_matlab/NIFTI_tools;

if nargin<5
    contrast_list=[];
end

if nargin<6 || isempty(TEMPLATE_VOL) || isempty(VOXDIMS)
     normflag = false;
else normflag = true;
end
if nargin<8 || isempty(DEOBLIQUE)
    DEOBLIQUE = 0;
end
if nargin<9 
    TPATTERN = [];
end
if nargin<10
    TOFWHM = 0; 
end
if nargin<11
    def_flag   = 1;
    opt_metric = 'dPR';
else
    def_flag   = 0;
end

%quick check to see if split_info specified, if not add it
[subject_inputs make_splitfile] = Read_Inputs(subject_inputs,analysis_model);
% now create the split_info file
if( make_splitfile>0 )
    
   split_info.TR_MSEC = TR_MSEC;
   split_info.type    = 'nocontrast';
   save('split_info_AUTO.mat','split_info');
end

%if pipelines=[], make "conservative"
if( isempty(pipelines) )
    error('Need to specify pipeline list');
end

% run preprocessing pipeline
Pipeline_PART1(subject_inputs,pipelines,analysis_model,[],0,contrast_list,false,DEOBLIQUE,TPATTERN);
% run optimization, if analysis model included
if( ~strcmp(analysis_model ,'NONE') )
    if(def_flag>0) disp('WARNING: pipeline optimization uses default metric of dPR'); end
    Pipeline_PART2(subject_inputs,analysis_model,contrast_list,opt_metric,[1 0],1,0);
end

% run spatial normalization
if(normflag)
    spatial_normalization(subject_inputs,analysis_model,contrast_list,TEMPLATE_VOL,VOXDIMS,0,DEOBLIQUE,WARP_TYPE);
    group_mask_tissue_maps( subject_inputs, [], TEMPLATE_VOL, WARP_TYPE );
end

%%
function [inputfile_update make_splitfile] = Read_Inputs(inputfile,analysis_model)

fid     = fopen(inputfile);
if fid==-1
    InputStruct = [];
    return;
end
% read in first line
tline               = fgetl(fid);
if ~ischar(tline)
    InputStruct = [];
    return;
end

% open new file to copy-over
[opath,oprefix,ext] = fileparts(inputfile);
inputfile_new = [opath, oprefix, '_new', ext];
fin = fopen(inputfile_new,'wt');

ksub = 0;
tflag= 0;
%% read subject_inputs
while ischar(tline) 
    
        ksub = ksub  +1;
        % check existence of TASK flag
        ifile = strfind( upper(tline), 'TASK=' ); 
        
        if( isempty(ifile) )
           if( strcmp(analysis_model ,'NONE') )
               
                tflag=1;
                % create temporary local split_info name
                ifile = strfind( upper(tline), 'IN=' ); 
                ifile = ifile+3;
                ips   = [strfind( tline, ' ' )-1 length(tline)];
                ips   = ips(ips>ifile);
                Input_nifti_file_temp = tline(ifile:ips(1));
                Input_nifti_file_temp = strrep(Input_nifti_file_temp,'.nii','');
                [Input_path,Input_prefix,ext] = fileparts(Input_nifti_file_temp);

                path=pwd;
                tline_new = [tline, ' TASK=',path,'/split_info_AUTO.mat'];
                fprintf(fin,'%s\n',tline_new);
           else
               error('cannot specify an analysis model without defining split-info structure (TASK in subject_inputs)');
           end
        end

        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
end
fclose(fid); 
fclose(fin); 

if( tflag==0 )
   delete(inputfile_new);
   inputfile_update = inputfile;
   make_splitfile   = 0;
else
   inputfile_update = inputfile_new;
   make_splitfile   = 1;
end
