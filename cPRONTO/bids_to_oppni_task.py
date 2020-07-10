#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    BIDS_TO_OPPNI_TASK: script to reformat bids json and tsv event files and constructs a set of OPPMI task files
"""
__author__      = "L. Mark Prati, Nathan, Pradeep"
__copyright__   = "Copyright 2019, The OPPNI Project"
__credits__     = ["Stephen Strother"]
__maintainer__  = "Mark Prati"
__email__       = "mprati@research.baycrest.org"
__status__      = "Development"
__license__     = ""
__version__     = "0.9"

import pprint
import json
import sys

def is_number(s):
    if s.isdigit():
        return True
    
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def bids_to_oppni_task(jsonfilelist, tsvfilelist, task_type, newtaskfilelist):
    '''
    Accept BIDS JSON and tsv files to produce OPPNI taskfile
    '''    
    #nsub = len(newtaskfilelist)
    #nrun = len(newtaskfilelist[1])
    
    for sub in newtaskfilelist.keys():
        for tsk in newtaskfilelist[sub].keys():
            for ses in newtaskfilelist[sub][tsk].keys():
                for i in range(0,len(newtaskfilelist[sub][tsk][ses])):
                    with open(jsonfilelist[tsk][ses][0]) as json_file:  
                        jdata = json.load(json_file)    
                        # check for 'RepetitionTime' in JSON data
                        try:
                            RepetitionTime = int(jdata['RepetitionTime'])
                        except:
                            print("ERROR: missing RepetitionTime field in JSON task file- {}".format(jsonfilelist[tsk][ses][0]))
                            exit()
                        json_file.close()
                    
                    if newtaskfilelist[sub][tsk][ses][i]:
                        with open(newtaskfilelist[sub][tsk][ses][i],'w+') as fout:
                            fout.write("UNIT=[msec]\n")    
                            fout.write("TR_MSEC=[{}]\n".format(str(1000*RepetitionTime))) 
                            fout.write("TYPE=[{}]\n".format(task_type)) 
                               
                            #2. task conditions, onsets, durations
                            if len(tsvfilelist[sub][tsk][ses]) > i:
                                tsvfile = tsvfilelist[sub][tsk][ses][i] 
                                with open(tsvfile) as fid:
                                    C_text = fid.readlines()
                                    # remove whitespace characters `\n` at the end of each line
                                    C_text = [x.strip('\n') for x in C_text] 
                                    fid.close();
                                    
                                headr = C_text[0]
                                body  = list(C_text) #duplicate content
                                #print('\BODY = ')
                                #pprint.pprint(body)
                                
                                colnames = headr.split('\t')
                                if ( colnames[0] != 'onset'): 
                                    print('ERROR: onset should be 1st column - {}'.format(tsvfile))
                                    exit()
                                    
                                if ( colnames[1] != 'duration' ): 
                                    print('ERROR: duration should be 2nd column - {}'.format(tsvfile))
                                    exit()
                                                                
                                ncol = len(colnames);
                               
                                nrow = len(body);                                                    
                                if nrow < 2:
                                    print('WARNING - tsv file contain no usable data : {}'.format(tsvfile))
                                    continue
                                                        
                                body.pop(0)             #remove header
                                body.pop()              #remove null line
                                nrow = len(body);       #reset length                                             
                                
                                #print('\POP-BODY = ')
                                #pprint.pprint(body)
                                
                                # find the trial type column
                                try:
                                    ixtrial = colnames.index('trial_type')
                                except ValueError:
                                    print('ERROR: tsv file missing trial_type specification : {}'.format(tsvfile))
                                    sys.exit()
                                    break
                                
                                #create a 2D array of values
                                bodycell= []
                                for k in range(0,nrow):
                                    bodycell.append(body[k].split('\t',))
                                    
                                #print('\ntrilcol is {} bodycell = '.format(ixtrial))
                                #pprint.pprint(bodycell)
                                
                                #transpose the cells
                                t= list(zip(*bodycell))
                                
                                #print('\ntransposed bodycell = ')
                                #pprint.pprint(t)
                                
                                #get take the trial_type row
                                d = t[ixtrial]    #list of trails             #matlab d=bodycell(:,ixtrial);
                                tasklist = list(set(d))       #unique list
                                #print('\ntasklist = ')
                                #print(tasklist)
                                
                                # collect all onsets, durations affiliated with task                
                                for k in range(0, len(tasklist)):
                                    onslist=[]
                                    durlist=[]
                                    for l in range (0, len(bodycell)):
                                        if d[l] == tasklist[k]:
                                            #print("\nbodycell[{}][1] = {}".format(l,bodycell[l][1]))
                                            if (is_number(bodycell[l][0]) and is_number(bodycell[l][1])):
                                                    onslist.append(bodycell[l][0])
                                                    durlist.append(bodycell[l][1])
                                    
                                    # print to file                                       
                                    if (onslist and durlist):        
                                        fout.write('NAME=[{}]\n'.format(tasklist[k].replace('-','_')))    
                                        fout.write('ONSETS={}\n'.format(onslist).replace("'","").replace(" ",""))    
                                        fout.write('DURATIONS={}\n'.format(durlist).replace("'","").replace(" ",""))  
                                    pass
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    pass