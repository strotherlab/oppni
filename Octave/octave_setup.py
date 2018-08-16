
from collections import namedtuple
import os
import subprocess
from os.path import basename as pbasename
from os.path import dirname as pdirname
from os.path import splitext as psplit
from os.path import join as pjoin

#version_check_location = base_work_dir = '/home/adlofts/Documents/Octave_Testing/test_rec_yng1t'

def main():
    """To be called when using octave for the first time"""

    # Load the CAC module
    print "Loading octave module..."
    try:
        subprocess.call("module load octave")
        subprocess.call("export OPPNI_PATH = /global/home/hpc3194/oppni-0.7.3.1_06JUL2017")
        print "HPC Environment Found..."
    except:
        print "Local Environment Found... set OPPNI_PATH yourself"

    # Install the packages from within octave, in order
    # If you have SUDO, and/or you dont have lib-dev Octave
    try:
         subprocess.call("apt-get install -y octave-io octave-control octave-struct octave-statistics octave-signal octave-optim")
         print "Installed Octave Packages with apt-get"
    except:
        # If you dont have SUDO, but you do have lib-dev Octave
        subprocess.call(["octave-cli","--eval","disp('Searching Octave forge...')"])
        subprocess.call(["octave-cli","--eval","pkg install -forge io"])
        subprocess.call(["octave-cli","--eval","pkg install -forge control"])
        subprocess.call(["octave-cli","--eval","pkg install -forge struct"])
        subprocess.call(["octave-cli","--eval","pkg install -forge statistics"])
        subprocess.call(["octave-cli","--eval","pkg install -forge signal"])
        subprocess.call(["octave-cli","--eval","pkg install -forge optim"])
        print "Installed Octave Packages using Octave Forge and no SUDO"

    print "Creating an .octaverc..."

    octaverc_string = """display('Octave OPPNI setting up...') 
                    addpath(genpath('/global/home/hpc3194/oppni-0.7.3.1_06JUL2017/scripts_matlab')) 
                    addpath(genpath('/global/home/hpc3194/oppni-0.7.3.1_06JUL2017/extra')) 
                    % Modify the --traditional settings 
                    try 
                        OCTAVE_VERSION; 
                            pkg load io 
                            pkg load control 
                            pkg load struct 
                            pkg load statistics 
                            pkg load signal 
                            pkg load optim 
                            disp('Octave Packages Loaded'); 
                            save_default_options("-mat7-binary"); 
                            disp('Changing default save to matlab version 7 -mat7-binary'); 
                            disable_range(false); 
                            disp('Disable_range was turned to false to correct median BUG'); 
                            confirm_recursive_rmdir(false); 
                            disp('Allowing rm -r commands without confirmation'); 
                            if exist ('do_braindead_shortcircuit_evaluation','builtin') 
                                 do_braindead_shortcircuit_evaluation(true); 
                                 disp('Enabling Matlab short cicuit evaluations'); 
                            end 
                    catch 
                        disp('No Octave version Found or Packages not installed'); 
                    end"""

    # Create .octaverc
    with open('~/.octaverc', 'w+') as makeoctaverc:
        makeoctaverc.write(octaverc_string)

    print "Moved .octaverc to ~/ /home/user"


if __name__ == '__main__':
    main()
