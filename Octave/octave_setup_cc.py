#
# Mark Prati
# Modified to set up an Octave environment for running OPPNI on Compute Canada (Cedar)
#
import os
import subprocess

def main():
    """To be called when using octave for the first time"""

    # Load the CC module
    print ("Loading octave version 4.4.1 ...")
    try:
        subprocess.call("module load nixpkgs/16.09  gcc/7.3.0  octave/4.4.1")
        subprocess.call("export OPPNI_PATH=/home/mprati/workspace/oppni/cPronto")
        print ("HPC Environment Found...")
    except:
        print ("Local Environment Found... set OPPNI_PATH yourself")

    # Install the packages from within octave, in order
    subprocess.call(["octave-cli","--eval","disp('Searching Octave forge...')"])
    subprocess.call(["octave-cli","--eval","pkg install -forge io"])
    subprocess.call(["octave-cli","--eval","pkg install -forge control"])
    subprocess.call(["octave-cli","--eval","pkg install -forge struct"])
    subprocess.call(["octave-cli","--eval","pkg install -forge statistics"])
    subprocess.call(["octave-cli","--eval","pkg install -forge signal"])
    subprocess.call(["octave-cli","--eval","pkg install -forge optim"])
    print ("Installed Octave Packages using Octave Forge")
    print ("Creating an .octaverc...")

    octaverc_string = """display('Octave OPPNI setting up...') 
                    addpath(genpath('/home/mprati/workspace/oppni/scripts_matlab')) 
                    addpath(genpath('/home/mprati/workspace/oppni/extra')) 
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

    print ("Moved .octaverc to ~/ /home/user")


if __name__ == '__main__':
    main()
