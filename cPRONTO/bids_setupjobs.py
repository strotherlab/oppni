# -*- coding: utf-8 -*-

"""
BIDS_SETUPJOBS: takes in formatted list of bids files and prepares
requisite oppni input files for running the pipelines
"""
__author__      = "L. Mark Prati, Nathan, Pradeep"
__copyright__   = "Copyright 2019, The OPPNI Project"
__credits__     = ["Stepen Strother"]
__maintainer__  = "Mark Prati"
__email__       = "mprati@research.baycrest.org"
__status__      = "Development"
__license__     = ""
__version__     = "0.9"

import pprint
from os import path
from os import makedirs
from make_input_file import make_input_file
from bids_to_oppni_task import bids_to_oppni_task 
        
#def bids_setupjobs(proc, outpath, varargin):
def bids_setupjobs(proc, outpath, fmri_in_list, fmri_out_list, struct_list, physio_list, drop1, drop2, jsonfile_list, tsvfile_list, reference_file, task_type):

# BIDS_SETUPJOBS: takes in formatted list of bids files and prepares
# requisite oppni input files for running the pipelines
#
#  Syntax:
#
#     bids_setupjobs(proc,outpath,varargin)
#

# varargins:
#
# 1 IN path
# 2 OUT path
# 3 STRUCT path
# 4 PHYSIO path
# 5 DROP value1
# 6 DROP value2
# 7 TASK JSON (BIDS)
# 8 EVENTS TSV
# 9 ATLAS PATH
# 10 contrast
# 11 task design
# 12 analysis model

## now setting parameter defaults

    # setting defaults, part1
    #analysis_model    = varargin[12]    #'LDA';
    #task_type         = varargin[11]    #'BLOCK';
    #contrast_list_str = varargin[10]    #[];
    modelparam        = [];
    dospnormfirst     = 0
    DEOBLIQUE         = 0
    TPATTERN          = 'auto_hdr';
    TOFWHM            = 0
    niiout            = 0

    # setting defaults, part2
    optimize_metric   = 'dPR'
    mot_gs_control    = [1,0]
    process_out       = 1
    keepmean          = 0
    whichpipes        = 'ALL'

    # setting defaults, spatnorm
    #reference_file    = varargin[9]
    input_voxelsize   = [2, 2, 2] ## temp
    flag_step         = 0
    
    makedirs(outpath, exist_ok=True)
    makedirs(path.join(outpath, 'task_files') , exist_ok=True)
    makedirs(path.join(outpath, 'input_files'), exist_ok=True)

    # pipeline file to be provided by oppni wrapper.
    # populate pipeline file
    #newpipefile = path.join( outpath, 'pipeline_combinations.txt')
    #if (path.isfile(newpipefile) == False):
    #    make_pipeline_file( newpipefile )
    print('create task file')
    pprint.pprint(fmri_out_list)
    
    
    #nsub = len(fmri_out_list)
    #nrun = len(fmri_out_list[0])
    
    # create task file 
    # need to check indexing on args matlab vs python
    print('\nBIDS_SETUPJOBS OPPNI OUT = ')
    pprint.pprint(fmri_out_list)
    
    newtaskfile = {}
    for sub in fmri_out_list.keys():
         newtaskfile[sub] = {}
         for tsk in fmri_out_list[sub].keys():
            newtaskfile[sub][tsk] = {}
            for ses in fmri_out_list[sub][tsk].keys():
                newtaskfile[sub][tsk][ses] = []
                for i in range(0 , len(fmri_out_list[sub][tsk][ses])):
                    #[path,prefix,ext] = fileparts(varargin[2][i][j]);
                    print("fmri_out_list[{}][{}][{}] = {}\n".format(sub,tsk,ses,fmri_out_list[sub][tsk][ses][i]))            
                    (opath,fn) = path.split(fmri_out_list[sub][tsk][ses][i])
                    prefix = fn.split('.')          
                    newtaskfile[sub][tsk][ses].append(path.join( outpath, 'task_files', prefix[0] + '.txt' ))  ## list of taskfile names

    print('NEW task file list')
    pprint.pprint(newtaskfile)
    
    bids_to_oppni_task(jsonfile_list, tsvfile_list, task_type, newtaskfile); ## populate full set of taskfiles
    
    # create "all subjects" input file or individual input file
    nsub = len(fmri_out_list)    
    if (nsub > 1):
        newinputfile = path.join( outpath, 'input_files', 'input_all_sub.txt' ) ## change to either "all_sub" or "single-sub"
    elif (nsub == 1):
        sub = list(fmri_out_list.keys())[0]
        tsk = list(fmri_out_list[sub].keys())[0]
        ses = list(fmri_out_list[sub][tsk].keys())[0][2:] #skip first 2 characters in ses
        prefix = sub + "_" + ses + '_' + tsk
        newinputfile = path.join( outpath, 'input_files', 'input_' + prefix +'.txt' ) ## change to either "all_sub" or "single-sub"
    pass
    make_input_file( newinputfile, fmri_in_list, fmri_out_list, struct_list, physio_list, drop1, drop2, newtaskfile ) ## modify to generate multi-line input file

    return newinputfile

    """
    proc = proc.upper()
    if proc == 'PART1'
        Pipeline_PART1(newinputfile,newpipefile, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE,  TPATTERN,   TOFWHM);
    elif proc == 'PART2'
    	Pipeline_PART2(newinputfile, optimize_metric, mot_gs_control, process_out, keepmean,   whichpipes)
    elif proc == 'SPNORM'    
        spatial_normalization(newinputfile,reference_file,input_voxelsize,flag_step,DEOBLIQUE)
    elif strcmpi(proc,'GMASK') 
        group_mask_tissue_maps(newinputfile,'');
    elif strcmpi(proc,'QC1') || strcmpi(proc,'QC2')
        group_mask_tissue_maps(newinputfile,'');
        QC_wrapper(flag_step,newinputfile, [], 2);
    else
        error('Unrecognized part name: must be one of PART1, PART2, SPNORM, GMASK, QC1 and QC2.');
    end
    """
    