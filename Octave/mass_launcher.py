#
# Mass Launch OPPNI against a data set cycling throuh a multiple parameter sets (see Looping Conditions) 
#
# Command line args to overright defaults added by Mark Prati (mprati@research.baycrest.org) 2019/01/08
# for instruction type:
# ]python mass_launcher.py -h
#

from collections import namedtuple
import os, re, getopt, sys, json, ast

from os.path import join as pjoin, exists as pexists
from os.path import dirname as pdirname
from os.path import basename as pbasename
from os.path import isdir as pisdir

# Default Parameter Set up
#####################################
# Looping Conditions
EXE_METHODS           = ['matlab','octave']
#CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'erCVA', 'erGNB', 'erGLM', 'SCONN']
CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'SCONN']
SLICE_TIMING_PATTERNS = ['alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z', 'auto_hdr']
DEOBLIQUE             = ['--DEOBLIQUE', '']

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
CONVOLVE_RANGE   = [0, 1]

# Consistent Conditions
metric                  = "dPR"
optim_schemes           = "ALL"

# Optional syntax(--call and value)
# TODO: Change these to reflect the use system location (local vs HPC)
# Mandatory Conditions
# TODO: Change this as needed
contrast                = "task_A-baseline"
# Optional
volume                  = ' -v "3.125 3.125 5.0"'
convolve                = "0"
cluster                 = ""
vasc_mask               = ""
blur                    = ""
control_motion_artifact = ""
control_wm_bias         = ""
output_nii_also         = ""
run_locally             = "--run_locally"
dry_run                 = "--dry_run"

###########################################


# Initialize
line_text = []
out_text = []

class color:
   BLUE = '\033[34m'
   MAGNETA = '\033[35m'
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   LIGHTBLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

# defining regexes that may be useful in various functions
# regex to extract the relevant parts of the input file
# Method from oppni.py
reOut = re.compile(r"OUT=([\w\./+_-]+)[\s]*")
def get_out_dir_line(line):
    """Returns the value of OUT=section in the input line."""

    return os.path.abspath(reOut.search(line).group(1))

def main(oppni_ver,input_file_og,pipeline_file,reference_file):
    # Read and Save the Original Input File to Work From
    with open(input_file_og, 'r') as ipf_og:
        for idx, line in enumerate(ipf_og.readlines()):
            line_text.append(line.strip())
            out_text.append(get_out_dir_line(line))  
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

                    os.makedirs(base_cond_dir, exist_ok=True)
                    os.makedirs(base_cond_env_dir, exist_ok=True)

                    # Make and Write a Modified Input File
                    input_file_mod = pjoin(base_cond_env_dir, 'input_file.txt')
                    with open(input_file_mod, 'w+') as ipf_mod:
                        for idx in range(len(line_text)):
                            # prepending it with OUT= to restrict the sub to only OUT, and not elsewhere such as TASK=
                            mod_out_text = 'OUT=' + pjoin(pdirname(out_text[idx]), suffix, str(env),pbasename(out_text[idx]))

                            out_text[idx] = 'OUT={}'.format(base_out_dir)
                            #DEBUG
                            #if input("\nout_text = [{}] Continue? (Y)".format(out_text[idx])) != "Y":
                            #    sys.exit()
 
                            mod_line_text = line_text[idx].replace(out_text[idx], mod_out_text)
                            #DEBUG
                            #if input("\nmod_line_text = [{}] Continue? (Y)".format(mod_line_text)) != "Y":
                            #    sys.exit()

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
                                          ' --convolve ' + convolve +
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
                    exit_code = os.system(submit_string) # TURN THIS ON AND OFF
                    if exit_code != 0:
                        print("\nExit Code {} was returned for command [{}]".format(exit_code,submit_string))
                        if (input("\nContinue with next command? (Y/N):").upper() == 'N'):
                            sys.exit();			

if __name__ == '__main__':

    # Manual Paths
    # Use command line args to reflect the use system location (local vs HPC)
    oppni_path     = './oppni/cPRONTO/oppni.py'  # path to oppni
    base_work_dir  = './Results'                 # Base working directory
    input_file     = 'input_file.txt'
    pipeline_file  = 'pipeline_file.txt'
    reference_file = 'reference_file.nii' #'MNI_orient_Nathan.nii'

    #TODO - accept sub option as command line args
    ANALYSIS_SUB_OPTIONS = {'Subspaces':CODES_SUBSPACES, 'SPM Types':CODES_SPM_TYPES, 'nblocks':NBLOCKS, 'wind':WIND, 'drf':DRF, 'Decision models':DECISION_MODELS, 'SPMS':SPMS, 'Subspaces':SUBSPACES, 'Convolve range':CONVOLVE_RANGE}

    #process command line args
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:o:i:p:r:e:m:d:l:",["base=","oppni=","input=","pipeline=","reference=","ex=","models=","dry=","local="])
    except getopt.GetoptError:
        print ('Options [ -b <Base Directory> ][ -o <OPPNI Directory> ][ -p <OPPNI pipeline file path> ][ -r <OPPNI reference file path> ][-m <analysis models array>][-d <dry run [Y/N]][-l <run locally [Y/N]]' )
        sys.exit(2)
        
    for opt, arg in opts:
        arg = arg.strip()
        if opt == '-h':
            print ('\n',color.UNDERLINE,'Command line syntax:',color.END,'\n')
            print ('mass_launcher.py [ -b <Base Working Directory> ][ -o <OPPNI path> ][ -i <OPPNI input file path> ][ -p <OPPNI pipeline file path> ][ -r <OPPNI reference file path> ][-m <analysis models array>][-m <exececution methods array>][-d <dry run [Y/N]]\n' )
            print (color.BLUE,'Default paths:',color.END)
            print ('\n<Base Working Directory>    : ', base_work_dir)
            print ('<OPPNI path>                : ', oppni_path)
            print ('<OPPNI input file path>     : ', input_file)
            print ('<OPPNI pipeline file path>  : ', pipeline_file)
            print ('<OPPNI reference file path> : ', reference_file)
            print ('use -b "" to eliminate default working directory')
            print ('use -b <Base Working Directory>" to alter base working directory')
            print ('file paths = <Base Working Directory>/./<file path>\n')
            print ('Options requiring an array of values must be entered as in the following example:\n')
            print (color.UNDERLINE,color.MAGNETA,'Example:',color.END)
            print (color.MAGNETA,'\nmass_launcher.py -m "[\'None\', \'LDA\', \'GNB\', \'GLM\']" -d N\n',color.END) 
            print (color.BLUE,'\nDefault options:',color.END) 
            print ('\n',color.UNDERLINE,'Execution Methods',color.END,'\n',EXE_METHODS,'\n')
            print (color.UNDERLINE,'Analysis Models',color.END,'\n',CODES_ANALYSIS_MODELS,'\n')
            print (color.UNDERLINE,'Analysis Sub-Options',color.END,'\n',ANALYSIS_SUB_OPTIONS,'\n')
            print (color.UNDERLINE,'Slice Timing Patterns',color.END,'\n',SLICE_TIMING_PATTERNS,'\n')

                        
            sys.exit()
        elif opt in ("-b", "--base"):
            base_work_dir  = arg
            print ('Base Working Directory set to: ', base_work_dir)
        elif opt in ("-o", "--oppni"):
            oppni_path  = arg
        elif opt in ("-i", "--input"):
            input_file  = arg
        elif opt in ("-p", "--pipeline"):
            pipeline_file  = arg
        elif opt in ("-r", "--reference"):
            reference_file  = arg
        elif opt in ("-e", "--ex"):
            EXE_METHODS  = ast.literal_eval(arg)
        elif opt in ("-m", "--models"):
            CODES_ANALYSIS_MODELS  = ast.literal_eval(arg)
        elif opt in ("-d", "--dry"):
            dry  = arg
            if dry == "N":
                dry_run = "";
        elif opt in ("-l", "--local"):
            dry  = arg
            if local == "N":
                run_locally = "";

    input_file_path     = pjoin(base_work_dir, input_file)
    pipeline_file_path  = pjoin(base_work_dir, pipeline_file)
    reference_file_path = pjoin(base_work_dir, reference_file)


    print (color.PURPLE,'\n\nExecuting OPPNI mass_launcher with the following parameters:\n')
    print ('\n<Base Working Directory>    : ', base_work_dir)
    print ('<OPPNI path>                : ', oppni_path)
    print ('<OPPNI input file path>     : ', input_file_path)
    print ('<OPPNI pipeline file path>  : ', pipeline_file_path)
    print ('<OPPNI reference file path> : ', reference_file_path, color.END)
    print ('\n',color.UNDERLINE,'Execution Methods',color.END,'\n',EXE_METHODS,'\n')
    print (color.UNDERLINE,'Analysis Models',color.END,'\n',CODES_ANALYSIS_MODELS,'\n')
    #print (color.UNDERLINE,'Analysis Sub-Options',color.END,'\n',json.dumps(ANALYSIS_SUB_OPTIONS, indent=4),'\n')
    print (color.UNDERLINE,'Analysis Sub-Options',color.END,'\n',ANALYSIS_SUB_OPTIONS,'\n')
    print (color.UNDERLINE,'Slice Timing Patterns',color.END,'\n',SLICE_TIMING_PATTERNS,'\n')


    # validate paths to files
    if os.path.isfile(oppni_path) == False:
        print(color.RED,'Invalid File Path Reference : ',oppni_path,color.END)
        sys.exit();
    if os.path.isfile(input_file_path) == False:
        print(color.RED,'Invalid File Path Reference : ',input_file_path,color.END)
        sys.exit();
    if os.path.isfile(pipeline_file_path) == False:
        print(color.RED,'Invalid File Path Reference : ',pipeline_file_path,color.END)
        sys.exit();
    if os.path.isfile(reference_file_path) == False:
        print(color.RED,'Invalid File Path Reference : ',reference_file_path,color.END)
        sys.exit();

    if (input("\nContinue with the paths and options shown? (Y/N):").upper() != 'Y'):
        sys.exit(); 

    main(oppni_path, input_file_path, pipeline_file_path,reference_file_path)

