# -*- coding: utf-8 -*-

"""
Description here
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

def make_input_file( filename, fmri_in_list, fmri_out_list, struct_list, physio_list, drop1, drop2, taskfile_list ):

# 1 IN         fmri_in_list[sub][tsk][ses]
# 2 OUT        fmri_out_list[sub][tsk][ses]
# 3 STRUCT path
# 4 PHYSIO path
# 5 DROP value1
# 6 DROP value2
    """
    print('\nmake_input_file:  fmri_in_list = ')
    pprint.pprint(fmri_in_list)
    print('\nmake_input_file:  fmri_out_list = ')
    pprint.pprint(fmri_out_list)
    print('\nmake_input_file:  struct_list = ')
    pprint.pprint(struct_list)
    """
    with open(filename,'wt') as fout:    
        #nsub = len(taskfile_list)
        #nrun = len(taskfile_list[1])
        for sub in taskfile_list.keys():        
            for tsk in taskfile_list[sub].keys():
                for ses in taskfile_list[sub][tsk].keys():
                    if (fmri_in_list[sub][tsk][ses]):
                        fout.write("IN={} ".format(fmri_in_list[sub][tsk][ses][0]))
                        fout.write("OUT={} ".format(fmri_out_list[sub][tsk][ses][0]))
                        if len(struct_list[sub][tsk][ses]):
                            fout.write("STRUCT={} ".format(struct_list[sub][tsk][ses][0]))
                        if (physio_list):
                            if len(physio_list[sub][tsk][ses]): 
                                fout.write("PHYSIO={} ".format(physio_list[sub][tsk][ses][0]))
                        fout.write("DROP=[{},{}] TASK={}\n".format(drop1, drop2, taskfile_list[sub][tsk][ses][0]))         
        fout.close()
