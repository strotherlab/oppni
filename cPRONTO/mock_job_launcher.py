
from collections import namedtuple
import os

# Set up of the key varibles
##############################################################################
#Looping
CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'erCVA', 'erGNB', 'erGLM', 'SCONN']

CODES_OPTIM_SCHEMES = [ 'CON', 'FIX', 'IND' ]

CODES_OPTIM_SCHEMES_OPTIONS = [ 'CON', 'FIX', 'IND', 'ALL' ]

SLICE_TIMING_PATTERNS = ['alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z', 'auto_hdr']


##Sub options for analysis models
CODES_SUBSPACES = [ 'onecomp','multicomp' ]  #erCVA
CODES_SPM_TYPES = ['zscore', 'corr']         #SCONN

##############################################################################
#Test Once
CODES_PRONTO_STEPS = [ 'PART1', 'QC1', 'PART2', 'QC2', 'GMASK', 'SPNORM' ]



##############################################################################
#Consistent

CODES_PREPROCESSING_STEPS = ["MOTCOR", "CENSOR", "RETROICOR", "TIMECOR",
             "SMOOTH", "DETREND", "MOTREG", "TASK",
             "GSPC1", "PHYPLUS", "CUSTOMREG", "LOWPASS"]

CODES_INPUT_SECTIONS = ['IN', 'OUT', 'DROP', 'STRUCT', 'PHYSIO', 'TASK', 'CUSTOMREG']

CODES_METRIC_LIST = ["dPR", "P", "R"]



TASK_MANDATORY_FIELDS = ['UNIT', 'TR_MSEC', 'TYPE']  #?????????????????????????????????????

WARP_TYPES = [ "affine", "nonlinear" ]   #????????????????????????????????????????????????

def main():
##Main Loop

     #Create String

     #Submit string to file
     oppni_ver = '/global/hpc4253/oppni-0.7.3.1_06JUL2017/cPRONTO/oppni.py'


     #Input file base
     input_file_base = '/global/home/hpc4253/YNG-REC-test90_inOctave.txt'
     #Load and Change the suffix
     #There maybe a case where the file name is not changed based on the options selected



     #Make and Load File to fill
     filename = '/global/home/hpc4253' + cmd_for_
     submit_string = 'python ' + oppni_ver + ' -i ' + input_file + ' -c '
     f = open()
     f.write(submit_string)

##EXAMPLE
python /global/home/hpc4253/oppni-0.7.3.1_06JUL2017/cPRONTO/oppni.py
-i /global/home/hpc4253/YNG-REC-test90_inOctave.txt
-c /global/home/hpc3194/frontenac/pipeline_file.txt -
a LDA --contrast "task_A-baseline"
-r /global/home/hpc3194/atlases/MNI_orient_Nathan.nii -
v "3.125 3.125 5.0"  --TPATTERN auto_hdr -m dPR -e octave  --cluster FRONTENAC




if __name__ == '__main__':
    main()
