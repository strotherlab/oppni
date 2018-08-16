
from collections import namedtuple
import os
import subprocess
from os.path import basename as pbasename
from os.path import dirname as pdirname
from os.path import splitext as psplit
from os.path import join as pjoin

#version_check_location = base_work_dir = '/home/adlofts/Documents/Octave_Testing/test_rec_yng1t'
checkscript_ver = '/home/adlofts/Documents/Octave_Testing/oppni_octave_git/extra/OPPNI_version_check_v3.m'
matlab_path = '/home/adlofts/Documents/Octave_Testing/oppni_octave_git/'

version_check_location = "/home/adlofts/Documents/Octave_Testing/trials_oppni/test_rec_yng1t"

def main():

    # Set Up

    checkscript,ext = psplit(pbasename(checkscript_ver))
    os.chdir(version_check_location)
    summary_file = pjoin(version_check_location,"Mass_Verification_Summary.txt")

    # Main Caller
    # Cycles through folders
    for option in os.listdir(version_check_location):
        if os.path.isdir(pjoin(version_check_location,option)):
            for environment in os.listdir(pjoin(version_check_location,option)):
                if os.path.isdir(pjoin(version_check_location,option,environment)):
                    for processed in os.listdir(pjoin(version_check_location, option, environment)):
                        if os.path.isdir(pjoin(version_check_location,option,environment,processed)):
                            # Gets input files ready
                            if "matlab" in environment :
                                matlabfile = pjoin(version_check_location, option, environment, processed, "input.txt")
                            if "octave" in environment :
                                octavefile = pjoin(version_check_location, option, environment, processed, "input.txt")
            # Get call ready
            savefile = pjoin(version_check_location,option,"MvsO_results.txt")
            print(savefile)
            call_command = "addpath(genpath('" + matlab_path + "')); " + checkscript + "('" + matlabfile + "','" + octavefile + "','" + savefile + "'); exit;"

            # Call matlab from command terminal to run Nathans Version Checker
            #subprocess.call(["matlab","-nodesktop","-nosplash","-r",call_command])

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
    main()
