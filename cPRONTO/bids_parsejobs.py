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

#
TEST_MODE = True

import json
import pprint
from bids_setupjobs import bids_setupjobs
from bids import BIDSLayout



def bids_parsejobs(bids_dir, output_dir, level, participant, run_name, task_design, drop1, drop2, atlasfile ):
    '''
    BIDS_PARSEJOBS: script takes in bids input arguments, then constructs a set of files, in preparation of setting up oppni analyses

    Syntax:
    
    bids_parsejobs( bids_dir, output_dir, level, participant, run_name, task_design, drop1, drop2, atlasfile )
    '''
        
    layout = BIDSLayout(bids_dir,validate=True)
    # Get lists of subjects, sessions
    sublist = layout.get_subjects()
    seslist = layout.get_sessions()
    tasklist = layout.get_tasks()
    fmri_in_list = {}
    fmri_out_list = {}
    struct_list = {}
    tsvfile_list = {}
    jsonfile_list = {}
    '''
    #build the run list
    #runfiles = layout.get(task='rest', extension='nii.gz', return_type='file')
    runfiles = layout.get(extension='nii.gz', return_type='file')
    
    for f in runfiles:
        #print("file = {}".format(f))
        #get runname : substring between 'run-' and '_bold' in the filename
        i1 = f.find('run-')
        i2 = f.find('_bold')
        #print("i1={}".format(i1))
        runnumber = f[int(i1 + 4):int(i2)]
        if runnumber not in runlist:                 
        	runlist.append(runnumber)
    '''        
    #BIDS spec subject at top level        
    for subj in sublist:
        if (level == "participant"):        
            if (participant and (subj != participant)):  #if individual:            
                continue
            
        if subj not in fmri_in_list: 
           fmri_in_list[subj] = {}
        if subj not in fmri_out_list: 
           fmri_out_list[subj] = {}
        if subj not in struct_list:
           struct_list[subj] = {}
        if subj not in tsvfile_list:
           tsvfile_list[subj] = {}
        
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
            if tsk not in tsvfile_list[subj]:
               tsvfile_list[subj][tsk] = {}
            if tsk not in jsonfile_list:
               jsonfile_list[tsk] = {}
            
            if (seslist):         
                for sess in seslist:
                    if sess not in fmri_in_list[subj][tsk]:
                       fmri_in_list[subj][tsk][sess] = []
                    if sess not in fmri_out_list[subj][tsk]:
                       fmri_out_list[subj][tsk][sess] = []
                    if sess not in struct_list[subj][tsk]:
                       struct_list[subj][tsk][sess] = []
                    if sess not in tsvfile_list[subj][tsk]:
                       tsvfile_list[subj][tsk][sess] = []
                    if sess not in jsonfile_list[tsk]:
                        jsonfile_list[tsk][sess] = {}
                        
                    jsonfile_list[tsk][sess] = layout.get(task=tsk, session=sess, suffix='bold', extension='json',return_type='filename') 
                    for runnumber in runlist:
                        print("tsk = {} : subj = {} : sess = {} : run = {}".format(tsk,subj,sess,runnumber)) 
                        fmri_in_list[subj][tsk][sess].append(layout.get(subject=subj, session=sess, run=int(runnumber), task=tsk, extension='nii.gz', return_type='file')[0])
                        #fmri_out_list[subj][tsk][sess].append(output_dir + '/' + subj + '_' + sess + '_task-' + tsk + 'run-' + runnumber)
                        fmri_out_list[subj][tsk][sess].append('/bidsResult' + '/' + subj + '_' + sess + '_task-' + tsk + 'run-' + runnumber)
                        struct_list[subj][tsk][sess].append(layout.get(subject=subj, session=sess, suffix='T1w', extension='nii.gz', return_type='file')[0])
                        tsvfile_list[subj][tsk][sess].append(layout.get(subject=subj, session=sess, run=int(runnumber), task=tsk, suffix='events', extension='tsv', return_type='file')[0])

            else:
                #no sessions so set session = '01' for lists.
                fmri_in_list[subj][tsk]['01'] = []
                fmri_out_list[subj][tsk]['01'] = []
                struct_list[subj][tsk]['01'] = []
                tsvfile_list[subj][tsk]['01'] = []
                jsonfile_list[tsk]['01'] = {}                    
                jsonfile_list[tsk]['01'] = layout.get(task=tsk, suffix='bold', extension='json',return_type='filename') #should only be one json file
                for runnumber in runlist:                                                   
                    fmri_in_list[subj][tsk]['01'].append(layout.get(subject=subj, run=int(runnumber), task=tsk, extension='nii.gz', return_type='file')[0])
                    #fmri_out_list[subj][tsk]['01'].append(output_dir + '/' + subj + '_task-' + tsk + 'run-' + runnumber)
                    fmri_out_list[subj][tsk]['01'].append('/bidsResult' + '/' + subj + '_task-' + tsk + 'run-' + runnumber)
                    struct_list[subj][tsk]['01'].append(layout.get(subject=subj, suffix='T1w', extension='nii.gz', return_type='file')[0])
                    tsvfile_list[subj][tsk]['01'].append(layout.get(subject=subj, run=int(runnumber), task=tsk, suffix='events', extension='tsv', return_type='file')[0])
    
    if TEST_MODE:
        print('\nOPPNI IN = ')
        pprint.pprint(fmri_in_list)
        print('\nOPPNI OUT = ')
        pprint.pprint(fmri_out_list)
        print('\nOPPNI STRUCT = ')
        pprint.pprint(struct_list)
        print('\nOPPNI TASK = ')
        pprint.pprint(tsvfile_list)
        print('\nJSONFILE = ')
        pprint.pprint(jsonfile_list)
        
    if (level == "participant"):
        """            
        if(participant):  #if individual:            
            fmri_in_list  = fmri_in_list[participant];
            fmri_out_list = fmri_out_list[participant];
            tsvfile_list  = tsvfile_list[participant];
            struct_list   = struct_list[participant];
        """ 
        #call bids_setupjobs to generate processed data (PART1)
        newinputfile = bids_setupjobs( 'PART1', output_dir, fmri_in_list, fmri_out_list, struct_list, {}, drop1, 0, jsonfile_list, tsvfile_list, [atlasfile], task_design);

    #LMP this is handled by wrapper and should not be called from within oppni wrapper
    elif (level == 'group'):
        #call bids_setupjobs (PART2,spnorm,groupmask)
        bids_setupjobs( 'PART2', output_dir, fmri_in_list, fmri_out_list, struct_list, {}, drop1, 0, jsonfile_list, tsvfile_list, [atlasfile], task_design);
        if (atlasfile):
            bids_setupjobs( 'SPNORM',output_dir, fmri_in_list, fmri_out_list, struct_list, [], drop1, 0, jsonfile_list, tsvfile_list, [atlasfile], task_design);
            bids_setupjobs( 'GMASK', output_dir, fmri_in_list, fmri_out_list, struct_list, [], drop1, 0, jsonfile_list, tsvfile_list, [atlasfile], task_design);
        pass
    pass

    return newinputfile


if __name__ == '__main__':
    
    if TEST_MODE:
        print("\nPython Test routine to build OPPNI input-file from BIDS structured data\n")
        bids_data_dir = "/global/home/hpc4346/BIDS_data_test/output_OND01_BYC_N"
        output_dir = "/global/home/hpc4346/BIDS_data_test/test_output"
    
    bids_parsejobs(bids_data_dir, output_dir, 'participant', 'OND01BYC1008', 'rest', "", 0, 0, "")
    
