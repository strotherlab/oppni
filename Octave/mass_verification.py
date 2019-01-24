#'/home/mprati/Results/RecYng_matlab_local4/processing_input_RecYng_matlab_local4_LDA_task_A-baseline/QC2_results/input_file.txt'
#script to run varification betwen OPPNI matlab results and octave results for a dataset 
#
from collections import namedtuple
import os, getopt, sys
import subprocess
from os.path import basename as pbasename
from os.path import dirname as pdirname
from os.path import splitext as psplit
from os.path import join as pjoin

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

#version_check_location = base_work_dir = '/home/adlofts/Documents/Octave_Testing/test_rec_yng1t'
#checkscript_ver = '/home/adlofts/Documents/Octave_Testing/oppni_octave_git/extra/OPPNI_version_check_v3.m'
#matlab_path = '/home/adlofts/Documents/Octave_Testing/oppni_octave_git/'
#version_check_location = "/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t"

checkscript_ver = './oppni_git/extra/OPPNI_version_check_v3.m'
matlab_path = './oppni_git/'
version_check_location = "./Results"

matlab_input_file_name = "input_file.txt"
octave_input_file_name = "input_file.txt"

def main():

    # Set Up

    checkscript,ext = psplit(pbasename(checkscript_ver))
    os.chdir(version_check_location)
    summary_file = pjoin(version_check_location,"Mass_Verification_Summary.txt")

    # Main Caller
    # Crawl through folders starting at base searching for matlab and octave results to verify
    matlabfile = ""
    octavefile = ""
    for option in os.listdir(version_check_location):
        if os.path.isdir(pjoin(version_check_location,option)):
            for environment in os.listdir(pjoin(version_check_location,option)):
                if os.path.isdir(pjoin(version_check_location,option,environment)):
                    print("Checking {} for OPPNI input files".format(pjoin(version_check_location,option,environment)))
#                    for processed in os.listdir(pjoin(version_check_location, option, environment)):
#                        if os.path.isdir(pjoin(version_check_location,option,environment,processed)):
#                            input_txt_file_path = pjoin(version_check_location, option, environment, processed) 
                    input_txt_file_path = pjoin(version_check_location, option, environment) 
                    # Gets input files ready
                    if "matlab" in environment :
                         matlabfile = pjoin(input_txt_file_path, matlab_input_file_name)
                    if "octave" in environment :
                         octavefile = pjoin(input_txt_file_path, octave_input_file_name)
#                        pass
                    pass
                pass

                if matlabfile == "":
                    print("Unable to locate OPPNI Matlab input.txt: ",input_txt_file_path)
                    continue;
                if octavefile == "":
                    print("Unable to locate OPPNI Octave input.txt: ",input_txt_file_path)
                    continue;
            pass
            if matlabfile == "" or octavefile == "" :
                print("Unable to locate file to varify, moving to next directory path")
                continue

            # Get call ready
            print("matlabfile = ",matlabfile)
            print("octavefile = ",octavefile)  
          
            savefile = pjoin(version_check_location,option,"MvsO_results.txt")
            print("Results will be saved to: ", savefile)
             

            call_command = "addpath(genpath('" + matlab_path + "')); " + checkscript + "('" + matlabfile + "','" + octavefile + "','" + savefile + "'); exit;"

            print("validation command line: ",call_command)

            #Call matlab from command terminal to run Nathans Version Checker
            subprocess.call(["matlab","-nodesktop","-nosplash","-r",call_command])

            # Load file and check for all [OK]
            totalok = 0
            with open(savefile, 'r') as result:
                for idx, line in enumerate(result.readlines()):
                    if "[OK]" in line:
                        totalok += 1
            summary_string = "Verified " + str(totalok) + "/14 [OK] for: " + pbasename(pdirname(savefile))

            # Save results into summary
            with open(summary_file, 'w+') as summary:
                 summary.write(summary_string)


if __name__ == '__main__':

    #process command line args - LMP
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hb:v:m:",["base=","script=","mat="])
    except getopt.GetoptError:
        print ('Options [ -b <Verification Check Base Directory > ][ -v <verification script path> ][ -m <matlab path> ]\n' )
        sys.exit(2)
        
    for opt, arg in opts:
        arg = arg.strip()
        if opt == '-h':
            print ('\n',color.UNDERLINE,'Command line syntax:',color.END,'\n')
            print ('mass_verification.py [ -b <Verification Check Base Directory> ][ -v <verification script path> ][ -m <matlab path> ]\n' )
            print (color.BLUE,'Default paths:',color.END)
            print ('\n<Base Directory>            : ', version_check_location)
            print ('<verification script path>  : ', checkscript_ver)
            print ('<matlab path>               : ', matlab_path)
            print ('use -b "" to eliminate default base directory')
            print ('use -b <Base Directory>" to alter base working directory')
            print ('file paths = <Base Working Directory>/./<file path>\n')                       
            sys.exit()
        elif opt in ("-b", "--base"):
            version_check_location  = arg
        elif opt in ("-v", "--script"):
            checkscript_ver  = arg
        elif opt in ("-m", "--mat"):
            matlab_path  = arg
 
    print (color.PURPLE,'\n\nExecuting OPPNI mass_verification with the following parameters:\n')
    print ('\n<Base Directory>            : ', version_check_location)
    print ('<verification script path>  : ', checkscript_ver)
    print ('<matlab path>               : ', matlab_path)

    if (input("\nContinue with the paths and options shown? (Y/N):").upper() != 'Y'):
        print (color.END)
        sys.exit(); 
 
    main()
