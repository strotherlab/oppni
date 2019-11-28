
from collections import namedtuple
import os,sys,subprocess
from os.path import basename as pbasename
from os.path import dirname as pdirname
from os.path import splitext as psplit
from os.path import join as pjoin

def main():
    """To be called when using octave for the first time"""

    print("Setting Up Octave Packages - Loading octave module...")
    try:
        subprocess.call("module load octave/4.4.1")
    except:
        print("WARNING: Unable to load octave")
        sys.exit()

    try:
        # Install the packages from within octave, in order
        subprocess.call(["octave-cli","--eval","disp('Searching Octave forge...')"])
        subprocess.call(["octave-cli","--eval","pkg install -forge io"])
        subprocess.call(["octave-cli","--eval","pkg install -forge control"])
        subprocess.call(["octave-cli","--eval","pkg install -forge struct"])
        subprocess.call(["octave-cli","--eval","pkg install -forge statistics"])
        subprocess.call(["octave-cli","--eval","pkg install -forge signal"])
        subprocess.call(["octave-cli","--eval","pkg install -forge optim"])
        print("Installed Octave Packages using Octave Forge")
    except:
        print("WARNING: Unable to load octave Packages")
        sys.exit()
        
    print("Creating an .octaverc...")
    #This needs to be modified to reflect the "Release source paths for scripts and extras"
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
    opath = os.path.expanduser('~/.octaverc')
    with open(opath, 'w+') as makeoctaverc:
        makeoctaverc.write(octaverc_string)

    print("Moved .octaverc to ~/ /home/user")
    print("Octave Package setup is complete")


if __name__ == '__main__':
    main()
