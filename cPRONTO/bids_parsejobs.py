#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    BIDS_PARSEJOBS: script takes in bids input arguments then constructs a set of files in preparation of setting up oppni analyses
"""
__author__      = "L. Mark Prati, Nathan, Pradeep"
__copyright__   = "Copyright 2019, The OPPNI Project"
__credits__     = ["Stepen Strother"]
__maintainer__  = "Mark Prati"
__email__       = "mprati@research.baycrest.org"
__status__      = "Development"
__license__     = ""
__version__     = "0.9"


TEST_MODE = True

import json
import pprint
from os import path
from bids_setupjobs import bids_setupjobs
from bids import BIDSLayout



def bids_parsejobs(bids_dir, input_dir, output_dir, level, participant_label, task_name, task_design, drop1, drop2, atlasfile ):
    '''
    BIDS_PARSEJOBS: script takes in bids input arguments, then constructs a set of files, in preparation of setting up oppni analyses

    Syntax:
    
    bids_parsejobs( bids_dir, input_dir, output_dir, level, participant_label, task_name, task_design, drop1, drop2, atlasfile )
    '''
    # check if we need to crate an input file of group processing
    if level.startswith("group"):
        if (level == "group"):
            groupName = "input_group_subs.txt"
        else:
            groupName = "input_group_{}_subs.txt".format(level[5:])
                
        newinputfile = path.join( input_dir, 'input_files', groupName )
        '''
        #LMP All this is now handled by main OPPNI wrapper 
        elif (level == 'group'):
            #call bids_setupjobs (PART2,spnorm,groupmask)
            bids_setupjobs( 'PART2', output_dir, fmri_in_list, fmri_out_list, struct_list, {}, drop1, 0, json_list, tsv_list, [atlasfile], task_design);
            if (atlasfile):
                bids_setupjobs( 'SPNORM',output_dir, fmri_in_list, fmri_out_list, struct_list, [], drop1, 0, json_list, tsv_list, [atlasfile], task_design);
                bids_setupjobs( 'GMASK', output_dir, fmri_in_list, fmri_out_list, struct_list, [], drop1, 0, json_list, tsv_list, [atlasfile], task_design);
            pass
        '''         
    else:                      
        fmri_in_list = {}
        fmri_out_list = {}
        struct_list = {}
        tsv_list = {}
        json_list = {}
        
        # Get lists of subjects, sessions, tasks
        layout = BIDSLayout(bids_dir,validate=True)
        sublist = layout.get_subjects()
        seslist = layout.get_sessions()
        if task_name:
            tasklist = task_name.split(",")
        else:
            tasklist = layout.get_tasks()
        
        #BIDS spec subject at top level        
        for subj in sublist:
            if level.startswith("participant"):        
                if (participant_label and (subj != participant_label)):  #if individual:            
                    continue
                
            if subj not in fmri_in_list: 
               fmri_in_list[subj] = {}
            if subj not in fmri_out_list: 
               fmri_out_list[subj] = {}
            if subj not in struct_list:
               struct_list[subj] = {}
            if subj not in tsv_list:
               tsv_list[subj] = {}
            
            for tsk in tasklist:
                runlist = []
                runfiles = layout.get(subject=subj, task=tsk, extension='nii.gz', return_type='file')
                for f in runfiles:
                    #get runname : substring between 'run-' and '_bold' in the filename
                    i1 = f.find('run-')
                    i2 = f.find('_bold')
                    #print("i1={}".format(i1))
                    runnumber = f[int(i1 + 4):int(i2)]
                    if runnumber not in runlist:                 
                        runlist.append(runnumber)            
                
                if tsk not in fmri_in_list[subj]:
                   fmri_in_list[subj][tsk] = {}
                if tsk not in fmri_out_list[subj]:
                   fmri_out_list[subj][tsk] = {}
                if tsk not in struct_list[subj]:
                   struct_list[subj][tsk] = {}
                if tsk not in tsv_list[subj]:
                   tsv_list[subj][tsk] = {}
                if tsk not in json_list:
                   json_list[tsk] = {}
                
                if (seslist):         
                    for sess in seslist:
                        if sess not in fmri_in_list[subj][tsk]:
                            fmri_in_list[subj][tsk][sess] = []
                        if sess not in fmri_out_list[subj][tsk]:
                            fmri_out_list[subj][tsk][sess] = []
                        if sess not in struct_list[subj][tsk]:
                            struct_list[subj][tsk][sess] = []
                        if sess not in tsv_list[subj][tsk]:
                            tsv_list[subj][tsk][sess] = []
                        if sess not in json_list[tsk]:
                            json_list[tsk][sess] = {}
                            
                        json_list[tsk][sess] = layout.get(task=tsk, session=sess, suffix='bold', extension='json',return_type='filename') 
                        for runnumber in runlist:
                            print("tsk = {} : subj = {} : sess = {} : run = {}".format(tsk,subj,sess,runnumber)) 
                            niifile_list = layout.get(subject=subj, session=sess, run=int(runnumber), task=tsk, extension='nii.gz', return_type='file')
                            #skip if no data
                            if niifile_list:
                                fmri_in_list[subj][tsk][sess].append(niifile_list[0])
                                fmri_out_list[subj][tsk][sess].append(path.join(output_dir, subj + '_' + sess + '_task-' + tsk + 'run-' + runnumber))
    
                                structfile_list = layout.get(subject=subj, session=sess, suffix='T1w', extension='nii.gz', return_type='file')
                                if structfile_list:
                                    struct_list[subj][tsk][sess].append(structfile_list[0])
                                
                                tsvfile_list = layout.get(subject=subj, session=sess, run=int(runnumber), task=tsk, suffix='events', extension='tsv', return_type='file')    
                                if tsvfile_list:    
                                    tsv_list[subj][tsk][sess].append(tsvfile_list[0])
    
                else:
                    #no sessions so set session = '01' for lists.
                    fmri_in_list[subj][tsk]['01'] = []
                    fmri_out_list[subj][tsk]['01'] = []
                    struct_list[subj][tsk]['01'] = []
                    tsv_list[subj][tsk]['01'] = []
                    json_list[tsk]['01'] = {}                    
                    json_list[tsk]['01'] = layout.get(task=tsk, suffix='bold', extension='json',return_type='filename') #should only be one json file
                    for runnumber in runlist:
                        niifile_list = layout.get(subject=subj, session=sess, run=int(runnumber), task=tsk, extension='nii.gz', return_type='file')
                        if niifile_list:                                             
                            fmri_in_list[subj][tsk]['01'].append(layout.get(subject=subj, run=int(runnumber), task=tsk, extension='nii.gz', return_type='file')[0])
                            fmri_out_list[subj][tsk]['01'].append(path.join(output_dir, subj + '_task-' + tsk + 'run-' + runnumber))
    
                            structfile_list = layout.get(subject=subj, suffix='T1w', extension='nii.gz', return_type='file')
                            if structfile_list:  
                                struct_list[subj][tsk]['01'].append(structfile_list[0])
    
                            tsvfile_list = layout.get(subject=subj, run=int(runnumber), task=tsk, suffix='events', extension='tsv', return_type='file')    
                            if tsvfile_list:                                
                                tsv_list[subj][tsk]['01'].append(tsvfile_list[0])
        
        if TEST_MODE:
            print('\nOPPNI IN = ')
            pprint.pprint(fmri_in_list)
            print('\nOPPNI OUT = ')
            pprint.pprint(fmri_out_list)
            print('\nOPPNI STRUCT = ')
            pprint.pprint(struct_list)
            print('\nOPPNI TASK = ')
            pprint.pprint(tsv_list)
            print('\nJSONFILE = ')
            pprint.pprint(json_list)
        
        #call bids_setupjobs to generate processed data (PART1)
        newinputfile = bids_setupjobs(input_dir, fmri_in_list, fmri_out_list, struct_list, {}, drop1, 0, json_list, tsv_list, [atlasfile], task_design);
    
    return newinputfile
