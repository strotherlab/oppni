
from collections import namedtuple
import os
import re
from oppni import get_out_dir_line
from os.path import join as pjoin
from os.path import dirname as pdirname
from os.path import basename as pbasename
from os.path import isdir as pisdir

# Set up
#####################################
# Looping Conditions
EXE_METHODS           = ['matlab','octave']
CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'erCVA', 'erGNB', 'erGLM', 'SCONN']
SLICE_TIMING_PATTERNS = ['alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z', 'auto_hdr']
DEOBLIQUE             = ['--DEOBLIQUE', '']

# Consistent Conditions
metric                  = "dPR"
optim_schemes           = "ALL"

# Manual Paths
# TODO: Change these to reflect the use system location (local vs HPC)
oppni_ver      = '/home/username/Documents/Octave_Testing/oppni_octave_git/cPRONTO/oppni.py'
base_work_dir = '/home/username/Documents/Octave_Testing'
input_file_og  = pjoin(base_work_dir, 'example_input_file.txt')
pipeline_file  = '/home/username/Documents/Octave_Testing/pipeline_file.txt'
reference_file = '/home/username/Documents/Octave_Testing/MNI_orient_Nathan.nii'

# Optional syntax(--call and value)
# TODO: Change these to reflect the use system location (local vs HPC)
# Mandatory Conditions
# TODO: Change this as needed
contrast                = "task_A-baseline"
# Optional
volume                  = ' -v "3.125 3.125 5.0"'
convolve                = ""
cluster                 = ""
vasc_mask               = ""
blur                    = ""
control_motion_artifact = ""
control_wm_bias         = ""
output_nii_also         = ""
run_locally             = "--run_locally"
dry_run                 = "--dry_run"

# Sub options for analysis models
#TODO: Add sub options in for each option into the loop
CODES_SUBSPACES  = ['onecomp','multicomp']     #erCVA
CODES_SPM_TYPES  = ['zscore', 'corr']            #SCONN
NBLOCKS          = ['None','4','10']             #erGNB and erCVA
WIND             = ['None','6','10']             #erGNB and erCVA
DRF              = ['None','0','1']              #erCVA
DECISION_MODELS  = ['None','linear','nonlinear'] #GNB
SPMS             = ['None','corr','zscore']      #SCONN
SUBSPACES        = ['None','onecomp','multicomp']#erCVA

# Initialize
line_text = []
out_text = []

def main():
    # Read and Save the Original Input File to Work From
    with open(input_file_og, 'r') as ipf_og:
        for idx, line in enumerate(ipf_og.readlines()):
            line_text.append(line.strip())
            out_text.append(get_out_dir_line(line))  # Method from oppni.py
    base_out_dir = pdirname(out_text[0])

    # Main Loop : Apply Conditions
    for model in CODES_ANALYSIS_MODELS:
        for pattern in SLICE_TIMING_PATTERNS:
            for dflag in DEOBLIQUE:

                # Generate Suffix
                suffix = "Opt_" + str(model) + '_' + str(pattern) + '_' + str(dflag)

                # 2 Execution Methods
                for env in EXE_METHODS:

                    # Create New Folders
                    base_cond_dir = pjoin(base_out_dir, suffix)
                    base_cond_env_dir = pjoin(base_cond_dir, str(env))
                    if not pisdir(base_cond_dir):
                        os.mkdir(base_cond_dir)
                    if not pisdir(base_cond_env_dir):
                        os.mkdir(base_cond_env_dir)

                    # Make and Write a Modified Input File
                    input_file_mod = pjoin(base_cond_env_dir, 'input_file.txt')
                    with open(input_file_mod, 'w+') as ipf_mod:
                        for idx in range(len(line_text)):
                            mod_out_text = pjoin(pdirname(out_text[idx]), suffix, str(env),pbasename(out_text[idx]))
                            mod_line_text = line_text[idx].replace(out_text[idx], mod_out_text)
                            ipf_mod.write(mod_line_text + "\n")

                    # Make and Write a Modified Command File
                    commandfile = pjoin(pjoin(base_cond_env_dir, "cmd_" + suffix + "_" + str(env) + '.sh'))
                    with open(commandfile,"w+") as cmd_fid:
                         submit_string = ('python ' + oppni_ver +
                                          ' -i ' + input_file_mod +
                                          ' -c ' + pipeline_file +
                                          ' -a ' + str(model) +
                                          ' -r ' + reference_file +
                                          ' '    + volume +
                                          ' --TPATTERN ' + str(pattern) +
                                          ' '    + convolve +
                                          ' -m ' + metric +
                                          ' -e ' + str(env) +
                                          ' '    + str(dflag) +
                                          ' '    + cluster +
                                          ' --contrast ' + contrast +
                                          ' '    + vasc_mask +
                                          ' '    + blur +
                                          ' '    + control_motion_artifact +
                                          ' '    + control_wm_bias +
                                          ' '    + output_nii_also +
                                          ' '    + dry_run +
                                          ' '    + run_locally)
                         cmd_fid.write(submit_string)

                    # Run OPPNI in Terminal
                    # Make Sure that System Settings are Correct
                    os.system(submit_string) # TURN THIS ON AND OFF

if __name__ == '__main__':
    main()
