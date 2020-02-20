#!/usr/bin/env python

#
#    OPPNI environment setup script.
#    Adjust the following files to allow OPPNI to run in a variety of environments:
#
#    .bash_profile
#    .octaverc
#    matlab/startup.m
#
#    This script assumes the following required software paths are accessible from your terminal:
#    oppni,afni,fsl,matlab,octave. 
#
#    If any required software can not be located by this script it will be reported as an error and the script will terminate
#    Once all prerequisite software packages have been installed, this script should be re-run.
#
#    Modified to support HPC environments using Octave
#    Also see "install_octave_packages.sh" bash script which is called from this script.
#    L. Mark Prati mprati@research.baycrest.org
#
# Updates:
# 24-07-19 - added octave image package
# 20-02-20 - tidy up..

import os, sys, subprocess, argparse, shutil
from time import localtime, strftime
from distutils.spawn import find_executable

__author__      = "L. Mark Prati, Pradeep Reddy Raamana"
__copyright__   = "Copyright 2019, The OPPNI Project"
__credits__     = ["Stephen Strother"]
__maintainer__  = "Mark Prati"
__email__       = "mprati@research.baycrest.org"
__status__      = "Development"
__license__     = ""
__version__     = "1.0.0"

identifier_matlab_native_mode = 'USING-MATLAB-CODE'.upper()

# default cvmfs computecanada stack versions of software available on CC and CAC clusters
CVMFS_AFNI =    '/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/afni/20180404' 
CVMFS_FSL =     '/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2018.3/fsl/5.0.11/fsl'
CVMFS_OCTAVE =  '/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc7.3/octave/4.4.1/bin'
CVMFS_MATLAB =  ''

def parse_args_check():
    "Parse inputs and different flags"

    parser = argparse.ArgumentParser(prog="setup_oppni")

    parser.add_argument("-o", "--oppni", action="store", dest="oppni_path",
                        default=None,
                        help="Absolute path to OPPNI installation location e.g. /home/user/oppni.")

    parser.add_argument("-a", "--afni", action="store", dest="afni_path",
                        default=None,
                        help="Absolute path to AFNI installation location e.g. /opt/afni.")

    parser.add_argument("-f", "--fsl", action="store", dest="fsl_path",
                        default=None,
                        help="Absolute path to FSL installation location e.g. /opt/fsl. "
                             "Do not include the bin/ directory at the end i.e. specify /opt/fsl, but not /opt/fsl/bin")
    parser.add_argument("--oct", action="store", dest="octave_path",
                        default=None,
                        help="Absolute path to octave installation location e.g. /opt/octave "
                             "This is needed only when you would like to run OPPNI using octave.")                        
    parser.add_argument("-m", "--mcr", action="store", dest="mcr_path",
                        default=None,
                        help="Absolute path to Matlab Component Runtime (MCR) installation location e.g. /opt/mcr_v80. "
                             "This is needed only when you would like to run OPPNI in the compiled mode. "
                             "If you would like to use OPPNI in its native form (which requires a Matlab license on the cluster), "
                             "please supply {} as the value for this argument.".format(identifier_matlab_native_mode))

    parser.add_argument("--assume_rotman_env", action="store_true", dest="assume_rotman_env",
                        default=False,
                        help="Preconfigured setup for the users located at the Rotman Research Institute making use of Rotman SGE.")
    
    parser.add_argument("-e", "--env", action="store", choices=['ROTMAN','CAC','CC', 'FRONTENAC', 'CEDAR', 'GRAHAM'], dest="env",
                        default='CUSTOM',
                        help="Select a preconfigured Strother-Lab setup for either Rotman: SGE, CAC: FRONTENAC, CC: CEDAR or GRAHAM choose from {'CAC,'CC','ROTMAN'}"
                             "Strother-Lab uses specific versions of supporting software such as AFNI and FSL")


    if len(sys.argv) < 2:
        print('Too few arguments!')
        parser.print_help()
        parser.exit(1)

    # parsing
    try:
        options = parser.parse_args()
    except:
        parser.print_help()
        parser.exit(1)
        
    #for backward compatibility
    if options.assume_rotman_env == True:
        print ('--assume_rotman_env overrides --env')
        options.env = 'ROTMAN'
        
    if options.env == 'CUSTOM':
        if (options.oppni_path is None) or (options.afni_path is None) or (options.fsl_path is None) or (options.mcr_path is None and options.octave_path is None):
            raise ValueError('Some of the the required paths are not specified.')
    elif options.env.upper() == 'ROTMAN' : #historical 
        print ('Preconfigured setup for Rotman is chosen. Ignoring the other paths provided, if any.')
        options.oppni_path  = '/opt/oppni'
        options.afni_path   = '/opt/afni'
        options.fsl_path    = '/opt/fsl/bin'        
        options.fsl_dir     = '/opt/fsl'
        options.mcr_path    = '/opt/mcr/v80'
        options.octave_path = '/opt/octave'
    elif options.env.upper() in [ 'CC', 'CEDAR', 'GRAHAM'] :    
        print ('Preconfigured Strother lab setup for CC CEDAR / GRAHAM chosen. Ignoring the other paths provided, if any.')
        options.oppni_path  = '/home/raamana/software/oppni'
        options.afni_path   = '/home/raamana/software/afni'
        options.fsl_path    = '/home/raamana/software/fsl/bin'
        options.fsl_dir     = '/home/raamana/software/fsl'
        options.mcr_path    = '/home/raamana/software/mcr/v80'
        options.octave_path = CVMFS_OCTAVE
    elif options.env.upper() in ['CAC', 'FRONTENAC'] :
        print ('Preconfigured Strother lab setup for CAC FRONTENAC chosen. Ignoring the other paths provided, if any.')
        options.oppni_path  = '/global/home/hpc3194/software/oppni'
        options.afni_path   = '/global/home/hpc3194/software/afni'
        options.fsl_path    = '/global/home/hpc3194/software/fsl/bin'
        options.fsl_dir     = '/global/home/hpc3194/software/fsl'
        options.mcr_path    = '/global/home/hpc3194/software/mcr/v80'
        options.octave_path = CVMFS_OCTAVE        
    else:
        print ('No Preconfigured Strother lab setup was chosen. The paths you provided (if any) will be applied.')

    #NOTE: when looking for "oppni.py" it would be found in sub-dir cPRONTO under the install directory
    # therefore backup one path level for correct oppni_path. 
    if not os.path.isdir(options.oppni_path):
        executible = shutil.which("oppni.py")
        if (executible):
            print("OPPNI was located here: {}".format(os.path.dirname(os.path.dirname(executible))))
            if input("Use this path? (Yes/No): ").upper() in ['Y','YES']:
                options.oppni_path = os.path.dirname(os.path.dirname(executible))
            else:
                raise ValueError("Specified OPPNI path {} doesn't exist!".format(options.oppni_path))
        else:
            #check if OPPNI is available in this setup path
            print("\nOPPNI was not found on the default path: {}\n".format(options.oppni_path)) 
            setup_dir = os.path.dirname(os.path.abspath(__file__))
            if os.path.isfile(os.path.join(setup_dir,"oppni.py")):
                print("However! OPPNI was located here: {}\n".format(setup_dir))
                if input("\nUse this path? (Yes/No): ").upper() in ['Y','YES']:
                    options.oppni_path =  os.path.dirname(setup_dir)
                else:
                    print("ERROR: A valid OPPNI path is required")
                    sys.exit()
            else:
                print("ERROR: A valid OPPNI path is required")
                sys.exit()
                    
    if not os.path.isdir(options.afni_path):
        #try to find afni on an exisiting path
        executible = shutil.which("afni")
        if (executible):
            print("AFNI was located here: {}".format(executible))
            if input("Use this path? (Yes/No): ").upper() in ['Y','YES']:
                options.afni_path = os.path.dirname(executible)
            else:
                print("\nSpecified AFNI path {} doesn't exist!".format(options.afni_path))
                print("WARNING: AFNI must be installed for OPPNI to run\n")
        elif options.env.upper() in ['CC','CAC','FRONTENAC','GRAHAM','CEDAR']:
            #check for cvmfs stack version
            if os.path.isdir(CVMFS_AFNI):
                print("Using AFNI located here: {}".format(CVMFS_AFNI))
                options.afni_path = CVMFS_AFNI
        else:    
            print("\nSpecified AFNI path {} doesn't exist!".format(options.afni_path))
            print("WARNING: AFNI must be installed for OPPNI to run\n")
            
    if not os.path.isdir(options.fsl_dir):
        executible = shutil.which("fsl")
        if (executible):
            print("FSL was located here: {}".format(executible))
            if input("Use this path? (Yes/No): ").upper() in ['Y','YES']:
                options.fsl_path = os.path.dirname(executible)
                options.fsl_dir = os.path.abspath(os.path.join(options.fsl_path, '..'))
            else:
                print("\nSpecified FSL path {} doesn't exist!".format(options.fsl_path))
                print("WARNING: FSL may be needed for some OPPNI functions\n")
        elif options.env.upper() in ['CC','CAC','FRONTENAC','GRAHAM','CEDAR']:
            #check for cvmfs stack version
            if os.path.isdir(CVMFS_FSL):
                print("Using FSL located here: {}".format(CVMFS_FSL))
                options.fsl_path = CVMFS_FSL
        else:    
            print("\nSpecified FSL path {} doesn't exist!".format(options.fsl_path))
            print("WARNING: FSL may be needed for some OPPNI functions\n")
        
    if not os.path.isdir(options.octave_path):
        executible = shutil.which("octave")
        if (executible):
            print("Octave was located here: {}".format(executible))
            if input("Use this path? (Yes/No: )").upper() in ['Y','YES']:
                options.octave_path = os.path.dirname(executible)
            else:
                print("Specified Octave path {} doesn't exist. Octave option will be unavailable!".format(options.octave_path))
                options.octave_path = None
        elif options.env.upper() in ['CC','CAC','FRONTENAC','GRAHAM','CEDAR']:
            #check for cvmfs stack version
            if os.path.isdir(CVMFS_OCTAVE):
                print("Using Octave located here: {}".format(CVMFS_OCTAVE))
                options.octave_path = CVMFS_OCTAVE

        else:    
            print("Specified Octave path {} doesn't exist. Octave option will be unavailable!".format(options.octave_path))
            options.octave_path = None

        
    if options.env.upper() != 'ROTMAN' or options.mcr_path.upper() == identifier_matlab_native_mode:
        print ('Ignoring the MCR setting according to user request')
    else:
        if not os.path.isdir(options.mcr_path):
            executible = shutil.which("mcr")
            if (executible):
                print("MCR was located here: {}".format(executible))
                if input("Use this path? (Yes/No): ").upper() in ['Y','YES']:
                    options.mcr_path = os.path.dirname(executible)
                else:
                    raise ValueError("Specified MCR path {} doesn't exist!".format(options.mcr_path))
            else:    
                raise ValueError("Specified MCR path {} doesn't exist!".format(options.mcr_path))

    return options


def validate_env_var(var):
    "Making sure changes will be reflected when these files are sourced."

    pass
    # assert os.getenv(var) is not None, "Path {} is not defined. Rerun the setup correctly.".format(var)


def validate_user_env(variables):
    "Making sure changes will be reflected when these files are sourced."

    for var in variables:
        validate_env_var(var)


def run_system_checks(options):
    "Runs a separate check script for each software installation, when available."

    pass


def add_matlab_paths(ms, oppni_path):
    #Adds the oppni matlab codes to the user matlab path.

    long_code_str="""
if ~isdeployed
    oppni_path = '{}';
    addpath(genpath(fullfile(oppni_path,'scripts_matlab')));
    addpath(genpath(fullfile(oppni_path,'extra')));
end
    """.format(oppni_path)

    ms.write('\n{}\n'.format(long_code_str))

def add_octave_paths(oc, oppni_path):
    #create and save octave startup commands in octaverc file
    octaverc_string = """
display('Octave OPPNI setting up...') 
addpath(genpath('""" + oppni_path + """/scripts_matlab')) 
addpath(genpath('""" + oppni_path + """/extra')) 
% Modify the --traditional settings 
try 
    OCTAVE_VERSION; 
    pkg load io 
    pkg load control 
    pkg load struct 
    pkg load statistics 
    pkg load signal 
    pkg load optim
    pkg load image 
    disp('Octave Packages Loaded'); 
    save_default_options("-mat7-binary"); 
    disp('Changing default save to matlab version 7 -mat7-binary'); 
    disable_range(false); 
    disp('Disable_range was turned to false to correct median BUG'); 
    confirm_recursive_rmdir(false); 
    disp('Allowing rm -r commands without confirmation'); 
    if exist ('do_braindead_shortcircuit_evaluation','builtin') 
        do_braindead_shortcircuit_evaluation(true); 
        disp('Enabling Matlab short circuit evaluations'); 
    end 
catch 
    disp('No Octave version Found or Packages not installed'); 
end
    """

    oc.write('\n{}\n'.format(octaverc_string))

def add_path_user_env(bp, path, var_name):
    "Adds a given path to the user environment and exports it ENV variable also"

    bp.write('\nexport {}="{}"'.format(var_name, path))
    bp.write('\nexport PATH=${{{}}}:$PATH'.format(var_name))

def add_user_env(fh, var_value, var_name):
    "exports ENV variable"

    fh.write('\nexport {}="{}"'.format(var_name, var_value))

def backup_file(file_path):
    "Make a time-stamped backup, if the file exists."

    if os.path.exists(file_path):
        time_stamp = strftime('%Y%m%d-T%H%M%S', localtime())
        shutil.copy(file_path, file_path + '_backup_{}'.format(time_stamp))


def make_needed_files_dirs():
    "Create ~/.bash_profile, and ~/matlab/startup.m if necessary."

    home_dir = os.getenv('HOME')

    #backup and create user bash_profile
    bash_profile = os.path.join(home_dir,'.bash_profile')
    backup_file(bash_profile)
    remove_section_from_textfile(bash_profile,"# Start of section added by setup_oppni:","# End of section added by setup_oppni:")
    if os.path.exists(bash_profile):
        bp = open(bash_profile,'a')
    else:
        bp = open(bash_profile, 'w')
    
    #backup and create matlab startup file
    matlab_path = os.path.join(home_dir,'matlab')
    if not os.path.isdir(matlab_path):
        os.mkdir(matlab_path)
    matlab_startup = os.path.join(matlab_path,'startup.m')
    backup_file(matlab_startup)
    if os.path.exists(matlab_startup):
        ms = open(matlab_startup, 'a')
    else:
        ms = open(matlab_startup, 'w')
        
    #backup and create octave startup file    
    octave_startup = os.path.join(home_dir,'.octaverc')
    backup_file(octave_startup)
    if os.path.exists(octave_startup):
        oc = open(octave_startup, 'a')
    else:
        oc = open(octave_startup, 'w')

    return bp, ms, oc


def remove_section_from_textfile(file_path, start_text, end_text):
    '''Remove a section of text from a file begining with the first line starting the "start_text" up to and including
    the line begining with "end_text. Note: only one section witll be removed per call
    NOTE: this will replace the oringinal file - to preserve the original file use "backup_file" before calling this function. 
    LMP''' 
      
    s_found = False
    e_found = False   
        
    # get a iterator over the lines in the file:
    with open(file_path, 'rt') as fin:
    # check that both start_text and end_text can be found
        for line in fin:
            if s_found == False and line.startswith(start_text):
                s_found = True
            
            if s_found == True and e_found == False and line.startswith(end_text):
                e_found = True
                break;    
        
        #only remove if a bound section is found
        if s_found and e_found:
            s_found = e_found = False   #reset flags
            fin.seek(0) #back to the beginning
            with open(file_path + ".new",'w') as fout:
                for line in fin:
                    #write output until the section start is found
                    if s_found == False and line.startswith(start_text) == False:                        
                        fout.write(line)
                    else:
                        s_found = True        

                    if s_found == True and e_found == False and line.startswith(end_text):
                        e_found = True
                        continue
                    
                    # write out the remainder of the file
                    if e_found:
                        fout.write(line)
                        
            fin.close()
            fout.close()
            try:
                shutil.copy(file_path + ".new", file_path)
            except:
                return False
                
            os.remove(file_path + ".new")
                
        else:
            fin.close()
            return False
               
    return True



def setup_paths():
    #Modifying or creating the necessary startup files or bash profiles to add different software to user environment.

    options = parse_args_check()
    if check_current_setup() == False:
        return False
    
    bp, ms, oc = make_needed_files_dirs()
    
    time_stamp = strftime('%Y%m%d-T%H%M%S', localtime())
    bp.write("\n# Start of section added by setup_oppni: {}\n".format(time_stamp))
    add_path_user_env(bp, options.oppni_path, 'OPPNI_PATH')
    add_path_user_env(bp, options.afni_path , 'AFNI_PATH')
    #provide path for optional load libaraies 
    bp.write('\nexport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${{{}}}'.format(options.afni_path + '/lib'))

    add_path_user_env(bp, options.octave_path, 'OCTAVE_PATH')
        
    add_path_user_env(bp, options.fsl_path      , 'FSL_PATH')
    add_user_env(bp, 'NIFTI_GZ', 'FSLOUTPUTTYPE')
    add_user_env(bp, options.fsl_dir, 'FSLDIR')
    


    # adding only what is necessary
    if options.mcr_path.upper() == identifier_matlab_native_mode:
        add_matlab_paths(ms, options.oppni_path)
    else:
        add_path_user_env(bp, options.mcr_path, 'MCR_PATH')
        validate_user_env(['MCR_PATH',])

    # checking if SGE is properly setup
    if options.env.upper() == 'ROTMAN': 
        if os.getenv('SGE_ROOT') is None:
            add_path_user_env(bp, '/usr/local/ge2011.11', 'SGE_ROOT')
            add_path_user_env(bp, '/usr/local/ge2011.11/bin/linux-x64', 'SGE_ROOT_BIN')

        if find_executable('qsub') is None:
            raise SystemError('SGE paths for user are not setup properly!')
            
    # checking if SLURM is properly setup for HPC clusters    
    if options.env.upper() in ['CAC','CC','FRONTENAC','GRAHAM','CEDAR']:
        #need to load octave module on HPC clusters
        bp.write('\nmodule load nixpkgs/16.09 gcc/7.3.0 octave/4.4.1')
        
        #load pyton 3.6
        bp.write('\nmodule load python/3.6\n')
        
        # checking if SLURM is properly setup for HPC clusters    
        if find_executable('sbatch') is None:
            #raise SystemError('SLURM (sbatch) paths for user are not setup properly!')
            pass

    # checking if SGE is properly setup    
    if options.env.upper() in ['ROTMAN']:
        if find_executable('qsub') is None:
            raise SystemError('SGE (qsub) paths for user are not setup properly!')
            pass
        
    # double check
    validate_user_env(['AFNI_PATH', 'FSL_PATH', 'OPPNI_PATH', 'OCTAVE_PATH','SGE_ROOT'])
    for execble in ['afni', 'fsl', 'oppni']:
        if find_executable(execble) is None:
            #raise SystemError('Executable {} can not be found. \n\t{} is not setup '.format(execble, execble.upper()))
            pass
        
    if find_executable('matlab') is None and find_executable('octave') is None:
        #raise SystemError('Neither MATLAB or OCTAVE executables can be found. OPPNI required one of these to be installed')
        pass
    
    #create and save octave startup commands in octaverc file
    add_octave_paths(oc, options.oppni_path)
    
    # properly closing them
    bp.write("\n# End of section added by setup_oppni: {}\n".format(time_stamp))
    
    bp.close()
    ms.close()
    oc.close()
    
    #finally install octave support pacages
    if options.octave_path:
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        subprocess.call(os.path.join(setup_dir,'install_octave_packages.sh'))
    else:
        print ('\e[31m\nWARNING: Skipping - install of octave packages\e[0m\n')
    return True


def check_current_setup():
    #Check for a pre-existing OPPNI setup, return True if existing setup is to be maintained
    oppni_path = os.getenv("OPPNI_PATH")
    afni_path = os.getenv("AFNI_PATH")
    fsl_path = os.getenv("FSL_PATH")
    octave_path = os.getenv("OCTAVE_PATH")
    
    if oppni_path or afni_path or fsl_path or octave_path:
        print("\nThe following OPPNI setup PATH environment variable exist:\n")
        if oppni_path:
            print("   OPPNI_PATH = {}".format(oppni_path))
        if afni_path:
            print("   AFNI_PATH = {}".format(afni_path))
        if fsl_path:
            print("   FSL_PATH = {}".format(fsl_path))
        if octave_path:
            print("   OCTAVE_PATH = {}".format(octave_path))
            
        #check for existence of some expected files:    
        home_dir = os.getenv('HOME')
        bash_profile = os.path.join(home_dir,'.bash_profile')
        if os.path.exists(bash_profile):
            print("\n   bash profile: {}".format(bash_profile))

        octaverc = os.path.join(home_dir,'.octaverc')    
        if os.path.exists(octaverc):
            print("   octaverc: {}".format(octaverc))
            
        return input("\nDo you want to overwrite the existing environment [Yes]/No: ").upper() in ['Y','Yes']
    else:
        print("\nNo previous setup found - continuing with OPPNI setup...")
        return True
        

if __name__ == '__main__':
    if setup_paths():
        print('OPPNI environment setup is almost complete!')
        print('Log out of your terminal and then log back in to activate the updated ".bash_profile"')
    else:
        print("\nYour existing environment has not been altered!\n")
exit(0)
            
        