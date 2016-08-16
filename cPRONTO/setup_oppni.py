#!/usr/bin/env python

import os, sys, argparse, shutil
from time import localtime, strftime
from distutils.spawn import find_executable

identifier_matlab_native_mode = 'USING-MATLAB-CODE'.upper()

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

    parser.add_argument("-m", "--mcr", action="store", dest="mcr_path",
                        default=None,
                        help="Absolute path to Matlab Component Runtime (MCR) installation location e.g. /opt/mcr_v80. "
                             "This is needed only when you would like to run OPPNI in the compiled mode. "
                             "If you would like to use OPPNI in its native form (which requires a Matlab license on the cluster), "
                             "please supply {} as the value for this argument.".format(identifier_matlab_native_mode))

    parser.add_argument("--assume_rotman_env", action="store_true", dest="assume_rotman_env",
                        default=False,
                        help="Preconfigured setup for the users located at the Rotman Research Institute making use of Rotman SGE.")

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

    if options.assume_rotman_env is False:
        if options.oppni_path is None or options.afni_path is None or options.fsl_path is None or options.mcr_path is None:
            raise ValueError('Some of the the required paths are not specified.')
    else:
        print 'Preconfigured setup for Rotman is chosen. Ignoring the other paths provided, if any.'
        options.oppni_path  = '/home/praamana/oppni/cPRONTO' # '/opt/oppni'
        options.afni_path   = '/opt/afni'
        options.fsl_path    = '/opt/fsl'
        options.mcr_path    = '/data1/strother_lab/praamana/software/mcr_install/v80'

    if not os.path.isdir(options.oppni_path):
        raise ValueError("Specified OPPNI path {} doesn't exist!".format(options.oppni_path))

    if not os.path.isdir(options.afni_path):
        raise ValueError("Specified AFNI path {} doesn't exist!".format(options.afni_path))

    if not os.path.isdir(options.fsl_path):
        raise ValueError("Specified FSL path {} doesn't exist!".format(options.fsl_path))

    if options.mcr_path.upper() == identifier_matlab_native_mode:
        print 'Ignoring the MCR setting according to user request'
    else:
        if not os.path.isdir(options.mcr_path):
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
    "Adds the oppni matlab codes to the user matlab path."

    long_code_str="""
if ~isdeployed
    oppni_path = '{}';
    addpath(genpath(fullfile(oppni_path,'scripts_matlab')));
    addpath(genpath(fullfile(oppni_path,'extra')));
end
    """.format(oppni_path)

    ms.write('\n{}\n'.format(long_code_str))


def add_path_user_env(bp, path, var_name):
    "Adds a given path to the user environment and exports it ENV variable also"

    bp.write('\nexport {}="{}"'.format(var_name, path))
    bp.write('\nexport PATH=${{{}}}:$PATH'.format(var_name))


def backup_file(file_path):
    "Make a time-stamped backup, if the file exists."

    if os.path.exists(file_path):
        time_stamp = strftime('%Y%m%d-T%H%M%S', localtime())
        shutil.copy(file_path, file_path + '_backup_{}'.format(time_stamp))


def make_needed_files_dirs():
    "Create ~/.bash_profile and ~/matlab/startup.m if necessary."

    home_dir = os.getenv('HOME')

    bash_profile = os.path.join(home_dir,'.bash_profile_test')
    backup_file(bash_profile)
    if os.path.exists(bash_profile):
        bp = open(bash_profile,'a')
    else:
        bp = open(bash_profile, 'w')

    matlab_path = os.path.join(home_dir,'matlab')
    if not os.path.isdir(matlab_path):
        os.mkdir(matlab_path)
    matlab_startup = os.path.join(matlab_path,'startup_test.m')
    backup_file(matlab_startup)
    if os.path.exists(matlab_startup):
        ms = open(matlab_startup, 'a')
    else:
        ms = open(matlab_startup, 'w')

    return bp, ms


def setup_paths():
    "Modifying or creating the necessary login or bash profiles to add different software to user environment."

    options = parse_args_check()

    bp, ms = make_needed_files_dirs()

    fsl_path_bin = os.path.join(options.fsl_path, 'bin')
    add_path_user_env(bp, options.oppni_path, 'OPPNI_PATH')
    add_path_user_env(bp, options.afni_path , 'AFNI_PATH')
    add_path_user_env(bp, fsl_path_bin      , 'FSL_PATH')

    # adding only what is necessary
    if options.mcr_path.upper() == identifier_matlab_native_mode:
        add_matlab_paths(ms, options.oppni_path)
    else:
        add_path_user_env(bp, options.mcr_path, 'MCR_PATH')
        validate_user_env(['MCR_PATH',])

    # checking if SGE is properly setup
    if options.assume_rotman_env is True and os.getenv('SGE_ROOT') is None:
        add_path_user_env(bp, '/usr/local/ge2011.11', 'SGE_ROOT')
        add_path_user_env(bp, '/usr/local/ge2011.11/bin/linux-x64', 'SGE_ROOT_BIN')

    if find_executable('qsub') is None:
        raise SystemError('SGE paths for user are not setup properly!')

    # double check
    validate_user_env(['AFNI_PATH', 'FSL_PATH', 'OPPNI_PATH', 'SGE_ROOT'])
    for execble in ['afni', 'fsl', 'oppni', 'matlab', 'qsub']:
        if find_executable(execble) is None:
            raise SystemError('Executable {} can not be found. \n\t{} is not setup properly'.format(execble, execble.upper()))

    # properly closing them
    bp.close()
    ms.close()


if __name__ == '__main__':
    setup_paths()