#!/usr/bin/env python
# OPPNI Tool: for fMRI pReprocessing and OptimizatioN Toolkit
# Author: Pradeep Reddy Raamana <praamana@research.baycrest.org>
# Version 0.6 (May 2016)
from __future__ import print_function
import argparse
import json
import os
import pickle
import random
import re
import stat
import subprocess
import sys
import tempfile
import math
import time
import traceback
import warnings
from collections import OrderedDict
from copy import copy
from distutils.spawn import find_executable
#from shutil import which -adlofts
from multiprocessing import Pool
from shutil import rmtree
from time import localtime, strftime
from datetime import timedelta


# Testing making my own which to allow ver < 3.3 - adlofts
import platform
if platform.python_version() < '3.3':
    print('Using older version of Python ...  defining which function')
    print(platform.python_version())

    def which(file):
        for path in os.environ["PATH"].split(os.pathsep):
            if os.path.exists(os.path.join(path, file)):
                return os.path.join(path, file)

        #return None
else:
    print('Using new version of Python ... using shutil import')
    print(platform.python_version())
    from shutil import which 


# OPPNI related
import cfg_front as cfg_pronto
import proc_status_front as check_proc_status

file_name_hpc_config = 'hpc_config.json'
file_name_job_ids_by_group = 'pronto_job_ids_by_step_prefix.json'
file_name_prev_options = 'pronto-options.pkl'

# descriptive variables
NOT_DONE = False
DONE = True

# noinspection PyGlobalUndefined
global hpc
hpc = {'type': 'SGE',
       'shell': '/bin/bash',
       'prefix': '#$',
       'spec': None,
       'header': '',
       'dry_run': False,
       'hold_jobid_list': [],
       'job_ids_grouped': {}}

# defining regexes that may be useful in various functions
# regex to extract the relevant parts of the input file
reIn = re.compile(r"IN=([\w\./+_-]+)[\s]*")
reOut = re.compile(r"OUT=([\w\./+_-]+)[\s]*")
reDrop = re.compile(r"DROP=\[(\d+),(\d+)\][\s]*")
reTask = re.compile(r"TASK=([\w\./+_-]+)[\s]*")
# to parse the task files
reName = re.compile(r"NAME=\[([\w\./+_-]+)\][\s]*")
reCondNames = re.compile(r"[^+-]+")

rePhysio = re.compile(r"PHYSIO=([\w\./+_-]+)[\s]*")
reStruct = re.compile(r"STRUCT=([\w\./+_-]+)[\s]*")
reCustReg = re.compile(r"CUSTOMREG=([\w\./+_-]+)[\s]*")

# for the pipeline file
rePip = re.compile('([0-9A-Z\s]+)=.+', re.IGNORECASE)
rePip2 = re.compile(r'([0-9A-Z\s]+)=\[([aA\d,]*)\][\s]*')


def get_out_dir_line(line):
    return os.path.abspath(reOut.search(line).group(1))


def get_out_dir_first_line(input_file):
    with open(input_file) as fID:
        # read one line
        first_line = fID.readline()
        # and move on

    return os.path.dirname(get_out_dir_line(first_line))


def make_time_stamp():
    # # with the minute
    # return  strftime('%Y%m%d-T%H%M',localtime())

    # just by the hour
    return strftime('%Y%m%d-T%H', localtime())


def validate_pipeline_file(pipeline_file):
    if not os.path.isfile(pipeline_file):
        raise IOError('Pipeline file specified doesn\'t exist')

    print('Validating the pipeline specs ...')
    steps_list_file = []
    with open(pipeline_file) as pip_f:
        for line in pip_f.readlines():
            steps_list_file.append(line.strip())

    steps_spec = rePip2.findall(' '.join(steps_list_file).upper())
    # steps_no_spec = map( lambda str1: str1.strip(' \n'), steps_no_spec)
    steps_dict = {}
    for step in steps_spec:
        if not step[0] in cfg_pronto.CODES_PREPROCESSING_STEPS:
            print('Error in pipeline file: %s' % pipeline_file)
            raise TypeError('Unrecognized pipeline step: %s' % step[0])

        step_values = step[1].replace(',', '')
        steps_dict[step[0]] = [ int(val) if val.isdigit() else val for val in step_values]

    print('  Done.')

    return steps_dict


def validate_task_file(task_path, cond_names_in_contrast=None):
    """Basic validation of task spec file."""

    with open(task_path, 'r') as tf:
        task_spec = tf.read().splitlines()
        # task_spec = [ line.strip('\n ') for line in task_spec ]
    task_spec = '\n'.join(task_spec)

    for field in cfg_pronto.TASK_MANDATORY_FIELDS:
        if field + '=' not in task_spec:
            raise TypeError('{} is not defined in task file'.format(field))

    if not_unspecified(cond_names_in_contrast):
        cond_names_in_file = reName.findall(task_spec)
        for name in cond_names_in_contrast:
            if name not in cond_names_in_file:
                raise ValueError("Condition {} in contrast is not defined in task file:"
                                 "\n {} \n Defined: {}".format(name, task_path, cond_names_in_file))

    return True


def not_unspecified( var ):
    """ Checks for null values of a give variable! """

    return var not in [ 'None', None, '' ]


def version_strings_differ(ver1, ver2):
    "Method to control the level of version match."

    # removing new lines and spaces to make comparison easy
    ver1 = ver1.lower().strip().decode('utf-8')
    ver2 = ver2.lower().strip().decode('utf-8')

    print(ver1)
    print(ver2)
    versions_differ = not ver1.startswith(ver2)

    return versions_differ


def validate_software_version(version_cmd_list, version_to_match, software_name, verbose = False):
    "Throw a warning if the user software version differs from the one tested."

    version_str = subprocess.check_output(version_cmd_list)
    if version_strings_differ(version_str, version_to_match):
        warnings.warn('\nYour {} version differs from the version tested by developers.'
                      '\nYours \n{} \nTested:\n{}'
                      '\nThis might cause differences in results.'.format(software_name, version_str, version_to_match))

    if verbose:
        print(version_str.strip().splitlines()[0])


def validate_env_var(var, verbose = False):

    upath = os.getenv(var)
    if upath is None:
        raise ValueError("Path {} is not defined. Fix your environment and rerun.".format(var))

    if verbose:
        print('{}: {}'.format(var, upath))

def validate_user_env(opt, verbose = False):
    """Validates set up of user's shell environment variables."""
    if not hpc['dry_run']:
        for var in ['AFNI_PATH', 'FSL_PATH']:
            validate_env_var(var, verbose)

        if opt.environment.lower() in ['compiled']:
            validate_env_var('MCR_PATH', verbose)

        validate_software_version(['afni',  '-ver' ]      , cfg_pronto.AFNI_VERSION_TESTED   , 'AFNI'   , verbose)
        validate_software_version(['flirt', '-version']   , cfg_pronto.FLIRT_VERSION_TESTED  , 'FLIRT'  , verbose)
        validate_software_version(['melodic', '--version'], cfg_pronto.MELODIC_VERSION_TESTED, 'MELODIC', verbose)


def validate_input_file(input_file, options=None, new_input_file=None, cond_names_in_contrast=None):
    """Key function to ensure input file is valid, and creates a copy of the input file in the output folders.
        Also handles the reorganization of output files depending on options chosen."""

    if (new_input_file is None) or (options is None) or (options.use_prev_processing_for_QC):
        # in case of resubmission, or when applying QC on an existing processing from older versions of OPPNI,
        # this should not append additional layer
        new_file = tempfile.TemporaryFile(mode='w')
        # in case of resubmission, this should not append additional layer
        cur_suffix = None
    else:
        new_file = open(new_input_file, 'w')
        cur_suffix = options.suffix

    unique_subjects = OrderedDict()
    invalid_lines = list()
    line_count = 0
    dupl_prefix_count = 0
    with open(input_file, 'r') as ipf:
        for line in ipf.readlines():
            line_count += 1
            subject, new_line = validate_input_line(line, cur_suffix, cond_names_in_contrast)
            if subject is False or subject is None:
                print('Error in line number {}.'.format(line_count))
                invalid_lines.append(line_count)
                continue
            else:
                subject['line'] = new_line
                # if the key doesnt exist, dict returns None
                if unique_subjects.get(subject['prefix']) is not None:
                    print("Potential duplicate prefix in line {}: {}".format(line_count, subject['prefix']))
                    print(" \t Previously processed line contained this prefix.")
                    dupl_prefix_count += 1
                unique_subjects[subject['prefix']] = subject
                new_file.write(new_line)

            # options = None is when this helper script is called from outside of OPPNI
            if options is not None and options.contrast_specified:
                assert subject['task'] is not None, \
                    'Contrast specified, but not a task file in line number {}.'.format(line_count)
            if options is not None and options.reference_specified:
                assert subject['struct'] is not None, \
                    'Reference atlas is specified, but not a structural scan in line number {}.'.format(line_count)
            # if one of the physiological methods are requested
            if options is not None and options.physio_correction_requested:
                assert subject['physio'] is not None, \
                    'RETROICOR is requested but physiological files not specified! Line number {}.'.format(
                        line_count)
            if options is not None and options.custom_mask_requested:
                assert subject['mask'] is not None, \
                    'CUSOMREG is specified but not a binary mask! Line number {}.'.format(line_count)

    new_file.close()

    # if the optional fields are specified, making sure they are specified for all the subjects
    optional_fields = ['task', 'struct', 'physio']
    for optf in optional_fields:
        bool_all_subjects = [sub[optf] is not None for ix, sub in unique_subjects.items()]
        num_runs_specified = sum([1 for ii in bool_all_subjects if ii is True])
        if not (num_runs_specified == 0 or num_runs_specified == len(bool_all_subjects)):
            print('Optional fields can either be specified for ALL the subjects, or NONE at all. ')
            raise ValueError(
                ' ---> {0} specified only for {1} out of {2} subjects'.format(optf, num_runs_specified,
                                                                              len(bool_all_subjects)))

    if len(invalid_lines) > 0:
        raise SyntaxError('Found {} lines erroneous: {}'.format(len(invalid_lines), ' '.join(map(str,invalid_lines))))
    else:
        print('Input file valid: Parsed {0} lines without error.'.format(line_count))
        if dupl_prefix_count > 0:
            print('\t excluding {} duplicate prefixes: {} lines'.format(dupl_prefix_count, len(unique_subjects)))

    return unique_subjects


def validate_input_line(ip_line, suffix='', cond_names_in_contrast=None):
    """Method where the real validation of the input line happens!"""

    line = ip_line.strip()
    LINE = line.upper()

    # commas are not allowed unless it's multi-run analysis
    num_commas = line.count(',')
    if num_commas > 1:
        # split by space
        line_sections = line.split(' ')
        for sec in line_sections:
            if "IN=" in sec.upper():
                in_sec = sec
            elif "OUT=" in sec.upper():
                out_sec = sec
            elif "DROP=" in sec.upper():
                if sec.count(',') > 1:
                    raise ValueError('Only 1 comma allowed in the DROP= section ')
            elif sec.count(',') > 0:
                raise ValueError('Commas are not allowed in the line '
                                 'except in IN= and OUT= sections when doing multi-run analaysis which is not supported yet.')
        if in_sec.count(',') != out_sec.count(','):
            raise ValueError('Number of commans in IN= section does not match those in OUT= section.'
                             'Ensure you provide equal number of output prefixes (comma separated) '
                             'in OUT= section as # runs in IN= section.')


    # defining an empty subject or run
    subject = {
        'nii': None,
        'prefix': None,
        'out': None,
        'drop_beg': None,
        'drop_end': None,
        'task': None,  # optional from here onwards
        'physio': None,
        'struct': None,
        'mask': None
    }

    in_flag = True
    out_flag = True
    drop_flag = True

    # IN part
    if "IN=" not in LINE:
        print("IN= section not defined.")
        return (False, "")
    else:
        nii = reIn.search(line).group(1)
        if not os.path.isfile(nii):
            print("Input file not found: " + nii)
            return (False, "")
        else:
            subject['nii'] = nii

    # OUT part
    if "OUT=" not in LINE:
        print("OUT= section not defined.")
        return (False, "")
    else:
        out = reOut.search(line).group(1)
        base_out_dir = os.path.dirname(out)
        # adding another directory level based on input file name and analysis model
        if suffix not in [None, '']:
            subject['out'] = os.path.join(base_out_dir, suffix)
        else:
            # in case of resubmission, don't alter the previous setup
            subject['out'] = base_out_dir

        if not os.path.exists(subject['out']):
            os.makedirs(subject['out'])

        # prepending it with OUT= to restrict the sub to only OUT, and not elsewhere such as TASK=
        prev_dir = 'OUT={}'.format(base_out_dir)
        curr_dir = 'OUT={}'.format(subject['out'])
        new_line = re.sub(prev_dir, curr_dir, ip_line)

        subject['prefix'] = os.path.basename(out)
        # Parse_Input_File.m is not robust with parsing e.g. an extra / at the end will mess up everything
        # so checking to make sure prefix is not empty
        assert (subject['prefix'] not in [None, '']), \
            'subject prefix can not be empty!'

    # DROP part
    if "DROP=" not in LINE:
        print("DROP= section not defined.")
        return (False, "")
    else:
        idx = reDrop.search(line).groups()
        if not len(idx) == 2:
            print("Atleast two indices must be specified to drop from the start and end.")
            return (False, "")
        elif (not all([ii.isdigit() for ii in idx])):
            print("DROP indices must be numeric!")
            return (False, "")
        elif not all([int(ii) >= 0 for ii in idx]):
            print("All the DROP indices must be non-negative.")
        else:
            subject['drop_beg'] = idx[0]
            subject['drop_end'] = idx[1]

    # the following parts are not mandatory, errors will be raised later on, when inconsistencies are found.
    # TASK part
    if "TASK=" in line:
        task = reTask.search(line).group(1)
        if not os.path.isfile(task):
            print("Task file " + task + " not found.")
            return None, None
        else:
            validate_task_file(task, cond_names_in_contrast)
            subject['task'] = task

    # PHYSIO part
    if "PHYSIO=" in LINE:
        physio = rePhysio.search(line).group(1)
        if not os.path.isfile(physio + '.puls.1D') or not os.path.isfile(physio + '.resp.1D'):
            print("PHYSIO files (puls and/or resp) at " + physio + " not found.")
            return None, None
        else:
            subject['physio'] = physio

    # STRUCT part
    if "STRUCT=" in LINE:
        struct = reStruct.search(line).group(1)
        if not os.path.isfile(struct):
            print("STRUCT file " + struct + " not found.")
            return None, None
        else:
            subject['struct'] = struct

    # PHYSIO part
    if "CUSTOMREG=" in LINE:
        mask = reCustReg.search(line).group(1)
        if not os.path.isfile(mask):
            print("Binary mask " + mask + " not found.")
            return None, None
        else:
            subject['mask'] = mask

    return subject, new_line


def parse_args_check():
    """Parser setup and assessment of different input flags."""
    parser = argparse.ArgumentParser(prog="oppni")

    parser.add_argument("-s", "--status", action="store", dest="status_update_in",
                        default=None,
                        help="Performs a status update on previously submitted processing. "
                             " Supply the output directory where the previous processing has been stored in.")

    parser.add_argument("--validate", action="store", dest="val_input_file_path",
                        default=None,
                        help="Performs a basic validation of an input file (existence of files and consistency!).")

    parser.add_argument("-p", "--part", action="store", dest="part", type=int,
                        default=0, choices=[0, 1, 2, 3, 4],
                        # #  "one", "two", "preproc", "optim", "spnorm", "all", "qc1", "qc2", "qc"]
                        help="select pipeline optimization step,"
                             "\n\t 0: All steps [default]"
                             "\n\t 1: Preprocessing ans statistics estimation step, "
                             "\n\t 2: Optimization step, "
                             "\n\t 3: Spatial normalization,"
                             "\n\t 4: quality control. ")
    parser.add_argument("-i", "--input_data", action="store", dest="input_data_orig",
                        help="File containing the input and output data paths", metavar="input spec file")
    parser.add_argument("-c", "--pipeline", action="store", dest="pipeline_file", metavar="pipeline combination file",
                        help="select the preprocessing steps")

    # parser.add_argument("-o","--output_dir_common", action="store", dest="out_dir_common",
    #                     help="Output folder to store all the processing and results. "
    #                          "This is convenient way to specify once, instead of having to repeat it on every line in the input file. "
    #                          "If you specify this, you don't have to include OUT= in the input file. "
    #                          "If you do, this will overwrite what is included in input file.")

    parser.add_argument("-a", "--analysis", action="store", dest="analysis",
                        default="None",
                        choices=cfg_pronto.CODES_ANALYSIS_MODELS,
                        help="Choose an analysis model :" + ",".join(cfg_pronto.CODES_ANALYSIS_MODELS))
    parser.add_argument("-m", "--metric", action="store", dest="metric",
                        default="dPR",
                        choices=cfg_pronto.CODES_METRIC_LIST,
                        help="Optimization metric")

    parser.add_argument("--os", "--opt_scheme", action="store", dest="opt_scheme",
                        default="ALL",
                        choices=cfg_pronto.CODES_OPTIM_SCHEMES_OPTIONS,
                        help="Optimization scheme to decide on which pipelines to be produced.")

    # parser.add_argument("--autodetect", action="store_true", dest="autodetect",
    #                    help="Automatically detect subjects and optimize for each subject independently, "
    #                         "the lines in input files that have same structrul image (STRUCT) and "
    #                         "output directory (OUT) are considered a subject")

    parser.add_argument("-r", "--reference", action="store", dest="reference",
                        default=None,
                        help="anatomical reference to be used in the spatial normalization step, i.e. -p,--part=3")
    parser.add_argument("--dospnormfirst",
                        action="store_true", dest="dospnormfirst", default=False,
                        help=argparse.SUPPRESS)
                        # help="First normalize the data to a reference (specified by switch -r), then perform the preprocessing optimization.")

    # TODO option to make the contrasts separable.
    parser.add_argument("--contrast", action="store", dest="contrast_list_str",
                        default="None",
                        help="desired task contrast in form of task-baseline, using names as defined in the task file. "
                             "Multi-contrast analysis is disabled for now.")

    # # TODO multi-contrast : doesn't work right now.
    # parser.add_argument("--contrast", action="store", dest="contrast_list_str",
    #                     default="None",
    #                     help="desired task contrast; This argument is necessary when more than two contrasts are "
    #                          "defined in the split info files .. "
    #                          "Syntax: --contrast \"CON1 vs CON2,CON2 vs CON3\". Output will be a consensus ouput "
    #                          "from all the contrasts.")

    parser.add_argument("--vasc_mask", action="store", dest="vasc_mask",
                        default="1", choices = ("0", "1"),
                        help="Toggles estimation of subject-specific vascular mask that would be excluded prior to analysis (0: disble, 1: enable). Recommended.")

    parser.add_argument("-k", "--keepmean", action="store", dest="keepmean",
                        default="0", choices = ("0", "1"),
                        help="(optional) determine whether the ouput nifti files contain the mean scan "
                             "(Default keepmean=0, i.e. remove the mean)")

    parser.add_argument("-v", "--voxelsize", action="store", dest="voxelsize",
                        help="(optional) determine the output voxel size of nifti file")

    parser.add_argument("--convolve", action="store", dest="convolve",
                        default="None",
                        help="VALUE=Binary value, for whether design matrix should be convolved with a standard SPMG1 HRF. "
                             "0 = do not convolve and 1 = perform convolution", metavar="VALUE")
    parser.add_argument("--decision_model", action="store", dest="decision_model",
                        default="None",
                        help="MODEL=string specifying type of decision boundary. Either: linear for a pooled covariance model or "
                             "nonlinear for class-specific covariances", metavar="MODEL")
    parser.add_argument("--drf", action="store", dest="drf",
                        default="None",
                        help="FRACTION=Scalar value of range (0,1), indicating the fraction of full-date PCA subspace to keep during PCA-LDA analysis. A drf of 0.3 is recommended as it has been found to be optimal in previous studies.",
                        metavar="FRACTION")
    parser.add_argument("--Nblock", action="store", dest="Nblock",
                        default="None",
                        help="NUMBER= number of equal sized splits to break the data into, in order to perform time-locked averaging. "
                             "Must be at least 2, with even numbers >=4, recommended to obtain robust covariance estimates",
                        metavar="NUMBER")
    parser.add_argument("--WIND", action="store", dest="WIND",
                        default="None",
                        help="SIZE = window size to average on, in TR (usually in range 6-10 TR)", metavar="SIZE")
    parser.add_argument("--num_PCs", action="store", dest="num_PCs",
                        default="None",
                        help="NUMBER = total number of principal components to retain", metavar="NUMBER")
    parser.add_argument("--subspace", action="store", dest="subspace",
                        default="None", choices=cfg_pronto.CODES_SUBSPACES,
                        help="COMP = string specifying either: 'onecomp'   = only optimize on CV#1 or "
                             "'multicomp' = optimize on full multidimensional subspace", metavar="SIZE")
    parser.add_argument("--spm", action="store", dest="spm",
                        default="None", choices=cfg_pronto.CODES_SPM_TYPES,
                        help="FORMAT =string specifying format of output SPM. "
                             "Options include corr (map of voxelwise seed correlations) or "
                             "zscore (Z-scored map of reproducible correlation values)", metavar="FORMAT")
    parser.add_argument("--N_resample", action="store", dest="N_resample",
                        default="10",
                        help=argparse.SUPPRESS)
                        # help="Specify the number of resamples for multi-run analysis")
    parser.add_argument("--TR_MSEC", action="store", dest="TR_MSEC",
                        default="None",
                        help="Specify TR in msec for all entries in the input file, overides the TR_MSEC in the TASK files")
    parser.add_argument("--DEOBLIQUE", action="store_true", dest="DEOBLIQUE",
                        default="0",
                        help="Correct for oblique scans (DEOBLIQUE) to improve spatial normalization")
    parser.add_argument("--TPATTERN", action="store", dest="TPATTERN",
                        default="",
                        help="Use if data contain no slice-timing information stored in the NIFTI headers (TPATTERN)")
    parser.add_argument("--BlurToFWHM", action="store", dest="BlurToFWHM",
                        default="0", choices=["0", "1"],
                        help="This option will enable adaptive spatial smoothing to to equalize smoothing across multiple sites.")

    # artifact control
    parser.add_argument("--control_motion_artifact", action="store", dest="ctrl_motion_artifact",
                        default='yes', choices=['yes', 'no'],
                        help="Control for motion artifact.")
    parser.add_argument("--control_wm_bias", action="store", dest="ctrl_wm_bias",
                        default='no', choices=['yes', 'no'],
                        help="Control for white matter bias using spatial priors.")

    parser.add_argument("--output_all_pipelines", action="store_true", dest="output_all_pipelines",
                        default="1",
                        help="Whether to output spms for all the optimal pipelines.")

    parser.add_argument("--output_nii_also", action="store", dest="output_nii_also",
                        default="0", choices=['1', '0'],
                        help="Whether to output all pipeline spms in Nifti format. \n WARNING: Be advised the space requirements will be orders of magnitude higher. Default off")

    # HPC related
    parser.add_argument("-e", "--environment", action="store", dest="environment",
                        default="compiled",
                        choices=('matlab','octave','compiled'),
                        help="(optional) determine which software to use to run the code: matlab or compiled(default)")

    parser.add_argument("--cluster", action="store", dest="hpc_type",
                        default=None, choices=('FRONTENAC', 'BRAINCODE', 'CAC', 'SCINET', 'SHARCNET', 'CBRAIN'),
                        help="Please specify the type of cluster you're running the code on.")

    parser.add_argument("--memory", action="store", dest="memory",
                        default="4",
                        help="(optional) determine the minimum amount RAM needed for the job, e.g. --memory 8 (in gigabytes)!")

    parser.add_argument("--walltime", action="store", dest="walltime",
                        default="30:00:00",
                        help="(optional) specify total run time needed for each job, e.g. --walltime 30:00:00 (in hours:minutes:seconds format)!")

    parser.add_argument("-n", "--numcores", action="store", dest="numcores",
                        default=1,
                        help=argparse.SUPPRESS)
    # help="DEPRECATED. number of threads used for the processing of a single run. (restricted to 1 now)")

    parser.add_argument("-q", "--queue", action="store", dest="queue",
                        default=None,
                        help="(optional) SGE queue name. Default is bigmem_16.q, but it is recommended to specify this explicitly.")

    parser.add_argument("-pe", "--parallel_env", action="store", dest="parallel_env",
                        default=None,
                        help="(optional) Name of the parallel environment under which the multi-core jobs gets executed. This must be specified explicitly when numcores > 1.")

    parser.add_argument("--run_locally", action="store_true", dest="run_locally",
                        default=False,
                        help="Run the pipeline on this computer without using SGE. This has not been fully tested yet, and is not recommended."
                             "Specify the number of cores using -n (or --numcores) to allow the program to run in pararallel")
    parser.add_argument("--noSGE", action="store_true", dest="run_locally",
                        default=False,
                        help=argparse.SUPPRESS)
                        # help="Same as --run_locally. Retained for backward compatibility.")
    parser.add_argument("--numprocess", action="store", dest="numprocess",
                        default=1,
                        help=argparse.SUPPRESS)
                        # help="When running the pipeline wihout using SGE, this switch specifies the number of simultaneous processes")

    parser.add_argument("--force_rerun", action="store_true", dest="force_rerun",
                        default=False,
                        help="DANGER: Cleans up existing processing and forces a rerun. Be absolutely sure this is what you want to do.")
    parser.add_argument("--dry_run", action="store_true", dest="dry_run",
                        default=False,
                        help="Generates job files only, but not run/sbumit them. Helps in debugging the queue/HPC options.")

    parser.add_argument("--use_prev_processing_for_QC", action="store_true", dest="use_prev_processing_for_QC",
                        default=False,
                        help=argparse.SUPPRESS)
                        # help="This option enables you to run QC jobs on existing processing generated with older versions of OPPNI.")

    parser.add_argument("--print_options_in", "--po", action="store", dest="print_options_path",
                        default=None,
                        help="Prints the options used in the previous processing of this folder.")

    parser.add_argument("--validate_user_env", action="store_true", dest="val_user_env",
                        default=False,
                        help=argparse.SUPPRESS) # "Performs a basic validation of user environment and version checks.")

    if len(sys.argv) < 2:
        print('Too few arguments!')
        parser.print_help()
        parser.exit(1)

    # parsing
    try:
        options = parser.parse_args()
    except:
        parser.exit(1)

    # updating the status to the user if requested.
    if options.status_update_in is not None and options.input_data_orig is None:
        cur_garage = os.path.abspath(options.status_update_in)
        update_status_and_exit(cur_garage)
    elif options.val_input_file_path is not None and options.input_data_orig is None:
        # performing a basic validation
        try:
            _ = validate_input_file(options.val_input_file_path)
            print(" validation succesful.")
        except:
            traceback.print_exc()
            print(" validation failed.")
        sys.exit(0)
    elif options.print_options_path is not None and options.input_data_orig is None:
        print_options(options.print_options_path)
        sys.exit(0)
    elif options.val_user_env and options.input_data_orig is None:
        try:
            validate_user_env(options, verbose=True)
            print('Versions (major) match those tested.')
            sys.exit(0)
        except UserWarning:
            warnings.warn('The versions of some of your software do not match the tested versions.')
            sys.exit(1)
        except:
            raise
    elif len(sys.argv) < 3:
        print('Invalid single arg! \n'
              'Use only one of \n --status \n--validate \n--print_options_in \n'
              'Or supply a full logical set of arguments for processing.\n'
              'For help, try oppni -h, or simply oppni')
        parser.exit(1)

    global hpc

    # sanity checks
    # on HPC inputs and obtaining the cfg of hpc
    options.numcores = int(options.numcores)
    hpc['type'] = options.hpc_type
    if options.run_locally is False:

        if options.numcores > 1:
            setattr(options, 'numcores', int(1))
            warnings.warn(
                '--numcores is specified. This flag is deprecated, and is restricted to 1 in favor of single-core jobs.')

        hpc['type'] = find_hpc_type(options.hpc_type, options.run_locally)
        hpc['spec'], hpc['header'], hpc['prefix'], hpc['shell'] = get_hpc_spec(hpc['type'], options)
    else:
        if not hpc['type'] in (None, 'LOCAL'):
            raise ValueError('Conflicting options specified: specify either of --run_locally or --cluster CLUSTERTYPE.')

    if hpc['type'] in (None, 'LOCAL'):
        if options.run_locally == False:
            print("Sun grid engine (SGE) has not been detected!")
            print("Use --run_locally switch if you want to run OPPNI without HPC cluster on your computer locally.")
            exit(1)
        else:
            print("Running jobs to the current node")
            print("The code will wait until the jobs are finished.")
            # even though it will be run locally, we will generate the job file
            # which will be executed in a subshell
            hpc['type'] = 'SGE'
    else:
        print("Submitting jobs to Sun Grid Engine (SGE)")

    if options.dry_run:
        hpc['dry_run'] = True

    # validating the pipeline file
    steps_dict = validate_pipeline_file(options.pipeline_file)
    setattr(options, 'pipeline_steps', steps_dict)

    # # enforcing the conditions of proper usage

    # if slice-timing correction is turned on
    if 1 in steps_dict['TIMECOR']:
        # make sure EPI acquisition pattern is specified
        if options.TPATTERN.lower() not in cfg_pronto.SLICE_TIMING_PATTERNS or options.TPATTERN in [None, ""]:
            raise ValueError( "Slice-timing correction is turned on - EPI acquisition pattern has not been specified"
                              " or is invalid. \n Must specify one of {}".format(cfg_pronto.SLICE_TIMING_PATTERNS))

    contrast_specified = not_unspecified(options.contrast_list_str)
    reference_specified = not_unspecified(options.reference)
    physio_correction_requested = (1 in options.pipeline_steps['RETROICOR'])
    custom_mask_requested = 1 in options.pipeline_steps['CUSTOMREG']
    vasc_mask_requested = "1" == options.vasc_mask

    if 1 not in options.pipeline_steps['PHYPLUS'] and vasc_mask_requested:
        warnings.warn("PHYCAA+ is not enabled, which is needed for the requested estimation of subject-specific vascular mask!")

    # saving useful flags for sanity checks on completeness of inputs later
    setattr(options, 'contrast_specified', contrast_specified)
    setattr(options, 'reference_specified', reference_specified)
    setattr(options, 'physio_correction_requested', physio_correction_requested)
    setattr(options, 'vasc_mask_requested', vasc_mask_requested)
    setattr(options, 'custom_mask_requested', custom_mask_requested)

    cur_garage, time_stamp, proc_out_dir, suffix = organize_output_folders(options)
    setattr(options, 'out_dir_common', cur_garage)
    setattr(options, 'suffix', suffix)

    if options.dospnormfirst and not options.reference_specified:
        raise ValueError('Spatial normalization requested, but a reference atlas is not specified.')

    # assert a single contrast to ensure QC doesnt fail either
    options.contrast_list_str = options.contrast_list_str.strip()
    # assert '-' in options.contrast_list_str, "Minus not found in the contrast string. Syntax: conditionA-conditionB"
    if not_unspecified(options.contrast_list_str):
        cond_names_in_contrast = reCondNames.findall(options.contrast_list_str)  # options.contrast_list_str.split('-')
    else:
        cond_names_in_contrast = None

    assert ',' not in options.contrast_list_str, "Multiple contrasts are specified with a comma! Only 1 contrast is allowed for now."

    # ensuring the validity of input file, given the options and pipeline file
    # and compiling a list of unique subjects into a new file with output folder modified
    new_input_file = os.path.join(cur_garage, 'input_file.txt')
    # validation also will create the new input file in the output folder
    unique_subjects = validate_input_file(options.input_data_orig, options, new_input_file, cond_names_in_contrast)
    if options.use_prev_processing_for_QC:
        new_input_file = options.input_data_orig

    # making sure user environment is properly setup before even submitting jobs
    # validate_user_env(options)

    ## -------------- Check the Input parameters  --------------

    # assert N > 3 to ensure QC/split-half mechanism doesnt fail
    assert len(unique_subjects) > 3, 'Too few (N<=3) runs in the input file. Rerun with N>3 runs.'

    if hasattr(options, 'reference') and options.reference is not None:
        reference = os.path.abspath(options.reference)
        if not os.path.isfile(reference):
            raise IOError('Reference file supplied does not exist!')
    else:
        options.reference = ""

    if options.analysis is None or options.analysis == "None":
        print("WARNING: without an analysis model (specified by switch -a), no optimization will be performed")
        print("  OPPNI will only generate the preprocessed data")
        options.contrast_list_str = "None"

    if hasattr(options, 'DEOBLIQUE') and (options.DEOBLIQUE == 1 or options.DEOBLIQUE is True):
        options.DEOBLIQUE = "1"
    else:
        options.DEOBLIQUE = "0"

    ## --------------  Checking the switches --------------
    analysis = options.analysis
    if (analysis.upper() == "LDA") and (options.drf == "None"):
        print("WARNING (Deprecated usage): --drf switch not defined for LDA model. OPPNI will check TASK files for parameter(s)")

    if (analysis.upper() == "ERCVA") and (options.drf == "None"):
        print("WARNING (Deprecated usage): --drf switch not defined for erCVA model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "ERCVA") and (options.Nblock == "None"):
        print("WARNING (Deprecated usage): --Nblock switch not defined for erCVA model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "ERCVA") and (options.WIND == "None"):
        print("WARNING (Deprecated usage): --WIND switch not defined for erCVA model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "ERCVA") and (options.subspace == "None"):
        print("WARNING (Deprecated usage): --subspace switch not defined for erCVA model. OPPNI will check TASK files for parameter(s)")

    if (analysis.upper() == "GNB") and (options.decision_model == "None"):
        print("WARNING (Deprecated usage): --decision_model switch not defined for GNB model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "ERGNB") and (options.Nblock == "None"):
        print("WARNING (Deprecated usage): --Nblock switch not defined for erGNB model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "ERGNB") and (options.WIND == "None"):
        print("WARNING (Deprecated usage): --WIND switch not defined for erGNB model. OPPNI will check TASK files for parameter(s)")

    if (analysis.upper() == "SCONN") and (options.spm == "None"):
        print("WARNING (Deprecated usage): --spm switch has to be used with the SCONN model. OPPNI will check TASK files for parameter(s)")

    if (analysis.upper() == "GLM") and (options.convolve == "None"):
        print("WARNING (Old style usage): --convolve switch has to be used with the GLM model. OPPNI will check TASK files for parameter(s)")
    if (analysis.upper() == "GPCA") and (options.num_PCs == "None"):
        print("WARNING (Deprecated usage): --num_PCs switch not defined for gPCA model. OPPNI will check TASK files for parameter(s)")

    if not (options.convolve in ["1", "0", "None"]):
        print("WARNING (Deprecated usage): --convolve has to be 0 or 1")
    if not (options.decision_model.lower() in ["linear", "nonlinear", "none"]):
        print("WARNING (Deprecated usage): --decision_model has to be linear or nonlinear")
    if not (options.subspace.lower() in ["onecomp", "multicomp", "none"]):
        print("WARNING (Deprecated usage): --subspace has to be onecomp or multicomp")
    if not (options.spm.lower() in ["corr", "zscore", "none"]):
        print("WARNING (Deprecated usage): --spm has to be corr or zscore")

    options.model_param_list_str = "keepmean " + options.keepmean \
                                   + " vasc_mask " + options.vasc_mask \
                                   + " convolve " + options.convolve \
                                   + " decision_model " + options.decision_model \
                                   + " drf " + options.drf \
                                   + " Nblock " + options.Nblock \
                                   + " WIND " + options.WIND \
                                   + " subspace " + options.subspace \
                                   + " spm " + options.spm \
                                   + " N_resample " + options.N_resample \
                                   + " TR_MSEC " + options.TR_MSEC \
                                   + " num_PCs " + options.num_PCs

    # removing unnecessary (leading, trailing or duplicate) white space
    options.model_param_list_str = options.model_param_list_str.strip()
    options.model_param_list_str = re.sub('\s+', ' ', options.model_param_list_str)
    ## -------------- END checking the parameters --------------


    os.environ["PIPELINE_NUMBER_OF_CORES"] = str(options.numcores)
    os.environ["FSLOUTPUTTYPE"] = "NIFTI"

    print("Chosen options: ")
    print(options)

    saved_cfg_path = os.path.join(cur_garage, file_name_prev_options)
    with open(saved_cfg_path, 'wb') as cfg:
        pickle.dump([unique_subjects, options, new_input_file, cur_garage], cfg, protocol=2)

    return unique_subjects, options, new_input_file, cur_garage, time_stamp, proc_out_dir


def organize_output_folders(options):
    """Organization of the folders and cleanup if needed."""

    parent_dir_wrapper = os.path.dirname(os.path.abspath(sys.argv[0]))

    proc_out_dir = get_out_dir_first_line(options.input_data_orig)
    if not os.path.exists(proc_out_dir):
        os.makedirs(proc_out_dir)

    time_stamp = make_time_stamp()

    if options.use_prev_processing_for_QC:
        print('Using the existing processing...')
        cur_garage = proc_out_dir
        suffix = ''
    else:
        # TODO need a way to pass the suffix to matlab core so it can reuse preprocessing for multiple analysis models and contrasts.
        suffix = os.path.splitext(os.path.basename(options.input_data_orig))[0]
        if options.analysis is not "None":
            suffix = suffix + '_' + options.analysis

        if options.contrast_specified:
            # task1-fixation,task2-task1 will be changed to task1-fixation_task2-task1
            suffix = suffix + '_' + options.contrast_list_str.replace(',', '_')

        # suffix = "garage_{0}_{1}".format(time_stamp,suffix)
        suffix = "processing_{0}".format(suffix)

        cur_garage = os.path.join(proc_out_dir, suffix)
        if os.path.exists(cur_garage):
            if options.force_rerun:
                user_confirmation = input("Are you sure you want to delete previous results? (y/[N])")
                if user_confirmation.lower() in ['y', 'yes', 'ye']:
                    print('Removing any existing preprocessing, as requested!')
                    rmtree(cur_garage)
                    os.mkdir(cur_garage)
                else:
                    print('Leaving the existing preprocessing as is!')
        else:
            os.mkdir(cur_garage)

    print("Output processing folder: " + cur_garage)

    return cur_garage, time_stamp, proc_out_dir, suffix


def find_hpc_type(user_supplied_type=None, run_locally=False):
    """Helper to determine (in a hacky way) the type of HPC environment is running on. """

    if user_supplied_type is None and run_locally is False:
        if find_executable('qsub') is None:
            raise SystemError('Requested processing on SGE, but command qsub is not found!')

        sge_root = os.getenv("SGE_ROOT")
        if sge_root is not None and find_executable('qconf') is not None:
            h_type = 'SGE'
        elif find_executable('qmgr') is not None and \
                (find_executable('qconf') is None and find_executable('sbatch') is None):
            h_type = 'TORQUE'
        elif find_executable('sbatch') is not None and find_executable('srun') is not None:
            h_type = 'SLURM'
        else:
            raise ValueError('Can not determine the type of cluster you are one.\n Please specify it with --cluster')
    else:
        h_type = user_supplied_type.upper()

    return h_type


def set_defaults_hpc(options, input_memory, input_queue, input_numcores, input_parallel_env,
        input_walltime='30:00:00'):
    """Assigns known defaults to HPC parameters"""

    if options.memory is None:
        memory = input_memory
    else:
        memory = options.memory

    if options.queue is None:
        queue = input_queue
    else:
        queue = options.queue

    if options.numcores is None or options.numcores < 1:
        numcores = input_numcores
    else:
        numcores = options.numcores

    if options.numcores >= 1:
        if options.parallel_env is None:
            parallel_env = input_parallel_env
        else:
            parallel_env = options.parallel_env
    else:
        parallel_env = None

    if options.walltime is None:
        walltime = input_walltime
    else:
        walltime = options.walltime

    return memory, queue, numcores, parallel_env, walltime


def get_hpc_spec(h_type=None, options=None):
    """Helper to provide the HPC directives for different HPC environments, such as SGE and Torque.
    SLURM is not fully supported yet."""

    if h_type is None:
        h_type = find_hpc_type().upper()

    shell = '/bin/bash'

    # assigning defaults to make it easy for the end user
    # TODO need to tease out the lists of names for different HPC environments into cfg_oppni.py
    if options is not None:
        if h_type in ('ROTMAN', 'ROTMAN-SGE', 'SGE'):
            memory, queue, numcores, parallel_env, walltime = set_defaults_hpc(options, 2, 'all.q', 1, 'npairs')
        elif h_type in ('CAC', 'HPCVL', 'QUEENSU'):
            memory, queue, numcores, parallel_env, walltime = set_defaults_hpc(options, 2, 'abaqus.q', 1, 'shm.pe')
        elif h_type in ('BRAINCODE-SGE', 'BRAINCODE', 'BCODE'):
            memory, queue, numcores, parallel_env, walltime = set_defaults_hpc(options, 2, 'common.q', 1, '')
        elif h_type in ('SCINET', 'PBS', 'TORQUE'):
            warnings.warn('HPC {} has not been tested fully. Use at your own risk!'.format(h_type))
            memory, queue, numcores, parallel_env, walltime = set_defaults_hpc(options, 2, 'batch', 1, '')
        elif h_type in ('FRONTENAC', 'SLURM'):
            memory, queue, numcores, parallel_env, walltime = set_defaults_hpc(options, 2, 'standard', 1, '',
                                                                               input_walltime='30:00:00')
    else:
        memory = '2'
        numcores = 1
        queue = 'all.q'
        parallel_env = None

    spec = {'memory': None, 'numcores': None, 'queue': None}

    if h_type in ('ROTMAN', 'ROTMAN-SGE', 'SGE'):
        prefix = '#$'
        spec['memory'] = ('-l mf=', memory + 'G')
        spec['numcores'] = ('-pe {} '.format(parallel_env), numcores)
        spec['queue'] = ('-q ', queue)
    elif h_type in ('CAC', 'HPCVL', 'QUEENSU'):
        prefix = '#$'
        spec['memory'] = ('-l mf=', memory + 'G')
        spec['numcores'] = ('-pe shm.pe ', numcores)
        spec['queue'] = ('-q ', queue)
        spec['export_user_env'] = ('-V', '')
        spec['workdir'] = '-wd'
        spec['jobname'] = '-N'
        spec['submit_cmd'] = 'qsub'
    elif h_type in ('BRAINCODE-SGE', 'BRAINCODE', 'BCODE'):
        prefix = '#$'
        spec['memory'] = ('-l mf=', memory + 'G')
        spec['numcores'] = ('-pe {} '.format(parallel_env), numcores)
        spec['queue'] = ('-q ', queue)
        spec['export_user_env'] = ('-V', '')
        spec['workdir'] = '-wd'
        spec['jobname'] = '-N'
        spec['submit_cmd'] = 'qsub'
    elif h_type in ('SCINET', 'PBS', 'TORQUE'):
        prefix = '#PBS'
        spec['memory'] = ('-l mem=', memory)
        spec['numcores'] = ('-l ppn=', numcores)
        spec['queue'] = ('-q ', queue)
    elif h_type in ('FRONTENAC', 'SLURM'):
        prefix = '#SBATCH'
        spec['memory'] = ('--mem=', int(memory)*1024) #
        spec['numcores'] = ('-c ', numcores)
        spec['queue'] = ('-p ', queue)
        spec['walltime'] = ('-t ', walltime)
        spec['export_user_env'] = ('--export=', 'ALL')
        spec['workdir'] = '--workdir'
        spec['jobname'] = '--job-name'
        spec['jobname_frontenac'] = '--output=%x_%j.log'
        spec['submit_cmd'] = 'sbatch'
        # slurm does not allow any shell specification
        shell=None
    else:
        raise ValueError('HPC type {} unrecognized or not implemented.'.format(h_type))

    header = list()
    for key in ['export_user_env', 'queue', 'memory', 'numcores', 'walltime']:
        # avoiding unnecessary specifications
        val = spec[key]
        if (key == 'numcores' and int(numcores) == 1) or (key == 'queue' and queue is None) or val is None:
            continue
        else:
            header.append('{0} {1}{2}'.format(prefix, val[0], val[1]))

    if h_type in ('FRONTENAC', 'SLURM'):
        # controlling the output job name
        header.append('{0} --output=%x_%j.log'.format(prefix))

    # not joining them for later use
    # header = '\n'.join(header)

    return spec, header, prefix, shell


def print_options(proc_path):
    """prints the various options specified by the user for previous processing."""
    proc_path = os.path.abspath(proc_path)
    assert os.path.exists(proc_path), "Specified processing folder does not exist!"
    opt_file = os.path.join(proc_path, file_name_prev_options)
    with open(opt_file, 'rb') as of:
        all_subjects, prev_options, prev_input_file, _ = pickle.load(of)
        attr_width = 29
        step_width = 10
        for attr in dir(prev_options):
            if not attr.startswith('_'):
                attr_val = prev_options.__getattribute__(attr)
                if attr == 'pipeline_steps':
                    print("{:>{}} : ".format(attr, attr_width))
                    for step, step_choices in attr_val.items():
                        print(" {:>{}}  {:>{}} : {}".format(' ', attr_width, step, step_width, step_choices))
                    print(" ")
                else:
                    print("{:>{}} : {}".format(attr, attr_width, attr_val))


def estimate_processing_times(input_file_all, options, all_subjects):
    """Helper to make a very rough estimate of processing time for entire workflow."""

    ipfile_modtime = os.path.getmtime(input_file_all)
    opt_summary = os.path.join(options.out_dir_common, 'optimization_results', 'matfiles', 'optimization_summary.mat')
    opt_summary_modtime = os.path.getmtime(opt_summary)

    proc_time_est = opt_summary_modtime - ipfile_modtime
    tdelta = timedelta(seconds=proc_time_est)

    if tdelta.total_seconds() > 0:

        num_days = tdelta.days
        rem_secs = tdelta.total_seconds()-num_days*24*3600

        num_minutes = math.floor(rem_secs/60)
        rem_secs = tdelta.total_seconds()-num_minutes*60

        readable_tdelta = ''
        if num_days > 0:
            readable_tdelta = '{} {:n} days'.format(readable_tdelta, num_days)

        if num_minutes > 1:
            readable_tdelta = '{} {:n} minutes'.format(readable_tdelta, num_minutes)

        if rem_secs > 0:
            readable_tdelta = '{} {:n} seconds'.format(readable_tdelta, rem_secs)

        print('Estimated processing time for whole workflow: \n\t {}.'.format(readable_tdelta))

    else:
        print('Estimated processing is 0, something seems to be wrong!')




def update_status_and_exit(out_dir):
    """Utility to update the user on the status of processing, and resubmit failed jobs for processing if necessary."""

    assert os.path.exists(out_dir), "Processing folder to traverse doesn't exist or not readable!"

    # print('Queue status update requested ... ')
    # update_Q_status(out_dir)
    print('\nNow checking the outputs on disk ...')
    try:
        prev_proc_status, prev_options, prev_input_file_all, \
        failed_sub_file, failed_spnorm_file, all_subjects = update_proc_status(out_dir)

        try:
            estimate_processing_times(prev_input_file_all, prev_options, all_subjects)
        except:
            print('Processing time could not be estimated.')

        if not prev_proc_status.all_done:
            print('Previous processing is incomplete.')
            if failed_sub_file is not None and failed_spnorm_file is not None:
                user_confirmation = input("Would you like to resubmit jobs for failed subjects/runs?  y / [N] : ")
                if user_confirmation.lower() in ['y', 'yes', 'ye']:
                    print(' Yes. \n\nAttempting resubmission ... \n')
                    try:
                        reprocess_failed_subjects(prev_proc_status, prev_options, failed_sub_file, failed_spnorm_file,
                                                  prev_input_file_all, all_subjects, out_dir)
                        print('Resubmission done. Check the status later.')
                    except:
                        print('Resubmission failed. Try again later.')
                        raise
                else:
                    # presenting the user with chosen option.
                    print(' No.')
            else:
                print('Unable to create input files for failed subjects/runs - make sure you have write permissions.')
    except:
        print('Unable to update status or resubmit failed jobs.')
        raise

    # to avoid entering a recursive state
    sys.exit(0)


def update_proc_status(out_dir):
    """Helper to the original update_status routine."""
    opt_file = os.path.join(out_dir, file_name_prev_options)

    with open(opt_file, 'rb') as of:
        all_subjects, options, new_input_file, _ = pickle.load(of)
        print(new_input_file)
        proc_status, failed_sub_file, failed_spnorm_file = check_proc_status.run(
            [new_input_file, options.pipeline_file, '--skip_validation'], options)
    return proc_status, options, new_input_file, failed_sub_file, failed_spnorm_file, all_subjects


# def update_Q_status(out_dir):
#     """Helper to track the status of the submitted jobs on the queue."""
#     if find_executable('qstat') is None:
#         print('\t qstat is not found - perhaps you are not on a headnode for the cluster. \n'
#               '\t Unable to query queue statuses. Run it on headnode to get the complete report.')
#         return
#
#     scheduler = drmaa.Session()
#     scheduler.initialize()
#
#     job_id_file = os.path.join(out_dir, file_name_job_ids_by_group)
#     try:
#         with open(job_id_file, 'rb') as jobs_list:
#             jobs_list_grouped = json.load(jobs_list)
#
#             for step, id_dict in jobs_list_grouped.items():
#                 print('  --> {}'.format(step))
#                 num_jobs_rqh, num_jobs_err = q_status(scheduler, id_dict)
#
#     except IOError as ioe:
#         print('Trouble reading the job IDs saved from last session. Check queue status manually. \n {}'.format(ioe))
#
#     # exiting the session
#     scheduler.exit()
#
#
# def q_status(scheduler, job_ids):
#     """
#      Display the status of each job, and count the number of jobs in error/Okay state
#
#     :param scheduler: drmma scheduler object
#     :param job_ids: dict of job IDs to track
#     :return: number of jobs running or in queue, and failed.
#
#     """
#
#     # decode_status = {drmaa.JobState.UNDETERMINED: 'status cannot be determined',
#     #                  drmaa.JobState.QUEUED_ACTIVE: 'queued and active',
#     #                  drmaa.JobState.SYSTEM_ON_HOLD: 'queued and in system hold',
#     #                  drmaa.JobState.USER_ON_HOLD: 'queued and in user hold',
#     #                  drmaa.JobState.USER_SYSTEM_ON_HOLD: 'queued and in user and system hold',
#     #                  drmaa.JobState.RUNNING: 'running',
#     #                  drmaa.JobState.SYSTEM_SUSPENDED: 'suspended by system',
#     #                  drmaa.JobState.USER_SUSPENDED: 'suspended by user',
#     #                  drmaa.JobState.DONE: 'finished normally',
#     #                  drmaa.JobState.FAILED: 'finished, but failed'}
#
#     num_jobs_rqh = 0  # num. jobs running or in queue
#     num_jobs_err = 0  # num. jobs in error state
#
#     for prefix, job_id in job_ids.items():
#         try:
#             job_state = scheduler.jobStatus(str(job_id))
#         except drmaa.InvalidJobException:
#             # when jobid has been flushed out of the scheduler memory!
#             job_state = drmaa.JobState.UNDETERMINED
#         else:
#             print("error quering the job status"
#             raise
#
#         print('\t {}: {} '.format(prefix, job_state))
#
#         # if job_state != 'unknown':
#         #     if 'e' in job_state:
#         #         num_jobs_err += 1
#         #
#         #     if 'r' in job_state or 'q' in job_state:
#         #         num_jobs_rqh += 1
#
#     # if num_jobs_rqh > 0:
#     #     print(' {} jobs are still either running, held or waiting to be processing.'.format(num_jobs_rqh))
#     #
#     # if num_jobs_err > 0:
#     #     print('{} jobs are stuck in error state. Manage them manually'.format(num_jobs_err))
#
#     return num_jobs_rqh, num_jobs_err

def make_job_file_and_1linecmd(file_path):
    """Generic job file generator, and returns one line version too."""

    hpc_directives = list()
    job_name = os.path.splitext(os.path.basename(file_path))[0]
    if not hpc['type'].upper() == "LOCAL":
        hpc_directives.append('#!/bin/bash')
        if hpc['shell'] is not None:
            hpc_directives.append('{0} -S {1}'.format(hpc['prefix'], hpc['shell']))
        hpc_directives.extend(hpc['header'])
        hpc_directives.append('{0} {1} {2}'.format(hpc['prefix'], hpc['spec']['jobname'], job_name))
        hpc_directives.append('{0} {1} {2}'.format(hpc['prefix'], hpc['spec']['workdir'], os.path.dirname(file_path)))
    else:
        # for jobs to run locally, no hpc directives are needed.
        hpc_directives.append('#!/bin/bash')

        #add adlofts Jun 26 2018
        #must add a directoty change so that the file can be found using MATLAB m file loading scheme
        #if environment.lower() in ('matlab'):
        print(os.path.dirname(file_path))
        print('FOUND MATLAB LOCAL')
        hpc_directives.append('cd {0} '.format(os.path.dirname(file_path)))

        # hpc_directives.append('{0} -S {1}'.format(hpc['prefix'], hpc['shell']))

    with open(file_path, 'w') as jID:
        # one directive per line
        jID.write('\n'.join(hpc_directives))

    st = os.stat(file_path)
    os.chmod(file_path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    return hpc_directives


def make_job_file(file_path):
    """Generic job file generator"""

    hpc_directives = list()
    job_name = os.path.splitext(os.path.basename(file_path))[0]
    if not hpc['type'].upper() == "LOCAL":
        if hpc.get('spec') is None:
            hpc['spec'], hpc['header'], hpc['prefix'], hpc['shell'] = get_hpc_spec(hpc['type'])

        hpc_directives.extend(hpc['header'])
        hpc_directives.append('{0} -N {1}'.format(hpc['prefix'], job_name))
        hpc_directives.append('{0} -wd {1}'.format(hpc['prefix'], os.path.dirname(file_path)))
    else:
        # for jobs to run locally, no directives are needed.
        hpc_directives.append('{0} -S '.format(hpc['shell']))

    with open(file_path, 'w') as jID:
        jID.write('\n'.join(hpc_directives))

    st = os.stat(file_path)
    os.chmod(file_path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def local_exec(script_path):
    """
    Runs a job script locally using subprocess, optionally returning stdout
    """

    # logger = logging.getLogger(script_path + '.log')

    # make it executable

    #adlofts add jun 26 2018
    script_path = list(script_path)
    script_path = script_path[0]
    st = os.stat(script_path)
    #st = os.stat(script_path)            #replaced with lines above
    print('Check the log file for progress:')
    print(script_path + '.log')

    os.chmod(script_path, st.st_mode | stat.S_IXGRP)

    proc = subprocess.Popen(script_path, shell=True, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # communicate waits for the subprocess to finish
    std_output, _ = proc.communicate()
    # # print(outputs and logs
    # logger.info('\n%s\n', std_output)
    print(std_output)

    return -random.randrange(1000)  # proc.returncode


def reprocess_failed_subjects(prev_proc_status, prev_options, failed_sub_file, failed_spnorm_file,
                              prev_input_file_all, all_subjects, garage):
    """Utility to resubmit the jobs to process the failed parts of the processing. """
    global hpc

    try:
        hpc_cfg_file = os.path.join(garage, file_name_hpc_config)
        with open(hpc_cfg_file, 'r') as hpc_f:
            hpc = json.load(hpc_f)
    except:
        print('hpc config was not saved previously or proerly! Unable to resubmit.')
        raise

    hpc['dry_run'] = False

    if prev_proc_status.preprocessing is NOT_DONE:
        failed_sub_p1 = validate_input_file(failed_sub_file, prev_options, None)
        # running the failed subjects throught part 1
        print('Resubmitting preprocessing jobs .. ')
        status_p1, jobs_p1 = run_preprocessing(failed_sub_p1, prev_options, failed_sub_file, garage)

    if prev_proc_status.optimization is NOT_DONE:
        # but optimization will be done entire dataset
        print('Resubmitting stats & optim. jobs .. ')
        status_p2, jobs_p2 = run_optimization(all_subjects, prev_options, prev_input_file_all, garage)

    if prev_proc_status.spnorm is NOT_DONE:
        failed_sub_spn = validate_input_file(failed_spnorm_file, prev_options, None)

        print('Resubmitting spatial normalization jobs .. ')
        sp_norm_step = 0  # both steps 1 and 2
        status_spn, jobs_spn = process_spatial_norm(failed_sub_spn, prev_options, failed_spnorm_file, sp_norm_step,
                                                    garage)

    if prev_proc_status.gmask is NOT_DONE:
        print('Resubmitting the group mask generation job .. ')
        status_gm , jobs_gm  = process_group_mask_generation(all_subjects, prev_options, prev_input_file_all, garage)

    if prev_proc_status.QC1 is NOT_DONE:
        # rerunning QC
        print('Resubmitting jobs for QC 1 .. ')
        status_qc1, jobs_qc1 = run_qc_part_one(all_subjects, prev_options, prev_input_file_all, garage)

    if prev_proc_status.QC2 is NOT_DONE:
        print('Resubmitting jobs for QC 2 .. ')
        status_qc2, jobs_qc2 = run_qc_part_two(all_subjects, prev_options, prev_input_file_all, garage)

    # saving the job ids to facilitate a status update in future
    job_id_file = os.path.join(garage, file_name_job_ids_by_group)
    if os.path.isfile(job_id_file):
        os.remove(job_id_file)
    with open(job_id_file, 'w') as jlist:
        json.dump(hpc['job_ids_grouped'], jlist)


def run_preprocessing(subjects, opt, input_file, garage):
    """
    Generates a job script to run the preprocessing for all combinations of pipeline steps requested.
    """

    if opt.dospnormfirst:
        str_dospnormfirst = '1'
    else:
        str_dospnormfirst = '0'

    # matlab: Pipeline_PART1(InputStruct, input_pipeset, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE, TPATTERN, TOFWHM)
    # input file will be prepended in the process module
    arg_list = [opt.pipeline_file, opt.analysis, opt.model_param_list_str, opt.output_nii_also,
                opt.contrast_list_str, str_dospnormfirst, opt.DEOBLIQUE, opt.TPATTERN, opt.BlurToFWHM, opt.keepmean]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'PART1', 'Pipeline_PART1', arg_list, garage, None)

    return proc_status, job_id_list


def run_qc_part_one(subjects, opt, input_file, garage):
    # maskname is set to be empty
    arg_list = [1, input_file, 'None', opt.num_PCs]
    dependencies = ['PART1', 'PART2', 'GMASK']
    proc_status, job_id_list = process_module_generic(subjects, opt, 'QC1', 'QC_wrapper', arg_list, garage,
                                                      dependencies)

    return proc_status, job_id_list


def run_qc_part_two(subjects, opt, input_file, garage):
    arg_list = [2, input_file, 'None', opt.num_PCs]
    dependencies = ['PART1', 'PART2', 'GMASK']
    proc_status, job_id_list = process_module_generic(subjects, opt, 'QC2', 'QC_wrapper', arg_list, garage,
                                                      dependencies)

    return proc_status, job_id_list


def run_optimization(subjects, opt, input_file, garage):
    """
    Generates a job script to run the optimization on the existing processing.
    """

    job_id_list = None

    # ideal:
    # arg_list = [ input_file, opt.metric, opt.ctrl_motion_artifact, opt.ctrl_wm_bias,
    #               opt.output_all_pipelines, opt.keepmean ]

    # current
    if opt.ctrl_motion_artifact == 'yes':
        mc_str = '1'
    else:
        mc_str = '0'
    if opt.ctrl_wm_bias == 'yes':
        wm_str = '1'
    else:
        wm_str = '0'
    mot_gs_control = mc_str + wm_str

    arg_list = [input_file, opt.metric, mot_gs_control, opt.output_all_pipelines, opt.keepmean, opt.opt_scheme]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'PART2', 'Pipeline_PART2', arg_list, garage,
                                                      'PART1')

    return proc_status, job_id_list


def process_spatial_norm(subjects, opt, input_file, sp_norm_step, garage):
    """
    Generates a job script to register all the MRI's of given subjects to a reference.

    :returns: status of processing and a list of job IDs (or process IDs if running locally).
    """

    if opt.dospnormfirst:
        dependencies = ['PART1']
    else:
        dependencies = ['PART1', 'PART2']

    arg_list = [opt.reference, opt.voxelsize, sp_norm_step, opt.DEOBLIQUE]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'SPNORM', 'spatial_normalization', arg_list,
                                                      garage, dependencies)

    return proc_status, job_id_list


def process_group_mask_generation(subjects, opt, input_file, garage):
    """Module to submit jobs for the generation of a group mask."""
    arg_list = [input_file]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'GMASK', 'group_mask_tissue_maps', arg_list,
                                                      garage, 'SPNORM')
    return proc_status, job_id_list


def construct_full_cmd(environment, step_id, step_cmd_matlab, arg_list, prefix=None, job_dir=None):
    """Helper to construct the necessary complete commands for various parts of the OPPNI processing."""
    if environment.lower() in ('matlab'):
        oppni_matlab_path = os.getenv('OPPNI_PATH_MATLAB_ORIG')
        single_quoted = lambda s: r"'{}'".format(s)

        cmd_options = ', '.join(map(single_quoted, arg_list))

        # make an m-file script
        # hyphen/dash can be treated as an operator
        mfile_name = prefix.replace('-', '_')
        mfile_path = os.path.join(job_dir, mfile_name + '.m')
        with open(mfile_path, 'w') as mfile:
            mfile.write('\n')
            mfile.write("addpath(genpath('{}'));".format(oppni_matlab_path))
            mfile.write("try, {0}({1}); catch ME, exc_report = getReport(ME); display('--->>> reporting exception details ..'); display(exc_report); display(' <<--- Done.'); exit(1); end; exit;".format(step_cmd_matlab,cmd_options))
            mfile.write('\n')

        # # problem: -nojvm omitted to allow future comptibility
        # cmd = 'time -p matlab -nojvm -nodesktop -nosplash'
        # full_cmd = "{0} -r \"try, {1}({2}); catch ME, display(ME.message); exit(1); end; exit;\" ".format(
        #     cmd, step_cmd_matlab, cmd_options)

        # to simplify the trouble with quotes and escape sequences etc
        # -nojvm flag to matlab is very important so it doesnt crash
        # time -p command is removed as it was failing on HPCVL
        setup_cmd = ''
        # if hpc['type'] == 'CAC':
        #    setup_cmd = r'use matlab'
        full_cmd = setup_cmd + "\n" + r"matlab -nodesktop -nosplash -r {0}".format(mfile_name)
    elif environment.lower() in ('octave'):
         
        # Same as Matlab but change out the full_cmd call and replace the getReport
         single_quoted = lambda s: r"'{}'".format(s)
         cmd_options = ', '.join(map(single_quoted, arg_list))
         oppni_octave_path = os.getenv('OPPNI_PATH') # The original path now defaults to octave use

         # make an m-file script
         mfile_name = prefix.replace('-', '_')
         mfile_path = os.path.join(job_dir, mfile_name + '.m')
         with open(mfile_path, 'w') as mfile:
            mfile.write('\n')
            mfile.write("addpath(genpath('{}'));".format(oppni_octave_path))
            mfile.write("try, {0}({1}); catch ; exc_report = lasterror; display('reporting exception details ..'); display(exc_report.message); display(exc_report.identifier); display(exc_report.stack(1).file); display(exc_report.stack(1).line); display(' <<--- Done.'); end; ".format(step_cmd_matlab,cmd_options))
            mfile.write('\n')

         setup_cmd = ''
         full_cmd = setup_cmd + "\n" + r"octave -W --traditional -q {0}".format(mfile_path)


    elif environment.lower() in ('standalone', 'compiled'):
        # first single quote is bash to protect from expansion or globbing
        #   second double quote is for matlab to interpret
        # strong_quoted = lambda s: r""" "'{}'" """.format(s)
        strong_quoted = lambda s: r""" "{}" """.format(s)

        # compiled_exe_path = find_executable('OPPNI')
        # time -p command is removed as it was failing on HPCVL
        compiled_exe_path = 'run_cPRONTO.sh ' + os.getenv('MCR_PATH')

        cmd = '{0} {1} '.format(compiled_exe_path, step_id)
        # enclosing the args with quotes only when needed
        quoted_args = list([])
        for arg in arg_list:
            if isinstance(arg, str) and set('[~!@#$%^&*()+{}":;-\' \t]+$').intersection(
                    arg):  # primitive: ' ' not in arg:
                # when special characters are present, enclose it in quotes
                quoted_args.append(strong_quoted(arg))
            else:
                quoted_args.append(str(arg))
        cmd_options = ' '.join(quoted_args)
        full_cmd = cmd + cmd_options

    else:
        raise ValueError('Unrecognized environment: {}'.format(environment))

    return full_cmd


def process_module_generic(subjects, opt, step_id, step_cmd_matlab, arg_list, garage, depends_on_step):
    """
    Generates a job script (per subject, or per dataset) to register all the MRI's of given subjects to a reference.
    :param step_id: identifier of the step being processed such as SPNORM, PREPROCESS, OPTIM
    :returns: status of processing and a list of job IDs (or process IDs if running locally).
    """
    global hpc

    input_dir = os.path.join(garage, 'input_files')
    job_dir = os.path.join(garage, 'job_files')
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    jobs_dict = {}
    if step_id.upper() in ['QC1', 'PART2', 'QC2', 'GMASK']:
        # these steps operate on the dataset as a whole
        # input list supplied from their individual functions
        arg_list_subset = copy(arg_list)
        prefix = '{0}_all_subjects'.format(step_id.lower())
        # each item will be a tuple (job_path, job_str)
        jobs_dict['all_subjects'] = make_single_job(opt.environment, step_id, step_cmd_matlab, prefix, arg_list_subset,
                                                    job_dir)

    else:
        for idx, subject in enumerate(iter(subjects.values())):
            # subject-wise processing
            # TODO input file doesnt change with step, try refactoring this to have only one input file per run/subject
            prefix = '{1}_s{0:0>3}_{2}'.format(idx + 1, step_id.lower(), subject['prefix'])
            subset_input_file = os.path.join(input_dir, prefix + '.input.txt')
            with open(subset_input_file, 'w') as sif:
                sif.write(subject['line'])

            # adding the input file as the first arg
            arg_list_subset = [subset_input_file] + arg_list
            # each item will be a tuple (job_path, job_str)
            jobs_dict[subject['prefix']] = make_single_job(opt.environment, step_id, step_cmd_matlab, prefix,
                                                           arg_list_subset, job_dir)

    jobs_status, job_id_list = run_jobs(jobs_dict, opt.run_locally, int(opt.numcores), depends_on_step)
    # storing the job ids by group to facilitate a status update in future
    hpc['job_ids_grouped'][step_id] = job_id_list

    return jobs_status, job_id_list


def make_single_job(environment, step_id, step_cmd_matlab, prefix, arg_list_subset, job_dir):
    """
    Helper to generate a standalone job file including the HPC directives and processing commands.
    """
    full_cmd = construct_full_cmd(environment, step_id, step_cmd_matlab, arg_list_subset, prefix, job_dir)

    out_job_filename = prefix + '.job'
    job_path = os.path.join(job_dir, out_job_filename)
    # create header with resource specs
    hpc_dir_1 = make_job_file_and_1linecmd(job_path)

    # add the commands and its arguments
    hpc_dir_2 = list()
    hpc_dir_2.append('')  # to get the newline concat working
    # hpc_directives.append('cd {0} '.format(job_dir))
    hpc_dir_2.append(r"{0}".format(full_cmd))

    with open(job_path, 'a') as jID:
        # remove any single quotes
        if environment.lower() in ('matlab', 'octave'):
            quote_del = [strg.replace(r"'", '') for strg in hpc_dir_2]
        else:
            quote_del = hpc_dir_2
        jID.write('\n'.join(quote_del))

    # complete single line command for qsub
    qsub_opt = hpc_dir_1 + hpc_dir_2
    # ensuring unnecessary prefix #$ #PBS is removed
    qsub_opt = [strg.replace(hpc['prefix'], '') for strg in qsub_opt]

    return job_path, qsub_opt


def run_jobs(job_paths, run_locally, num_procs, depends_on_step):
    """ Tool to either submit jobs or run them locally, as directed by the user."""

    job_id_list = {}
    if run_locally:

        if not hpc['dry_run']:
            # ret_code = local_exec(job_path, log_fname, print_to_screen=True)
            # if ret_code < 0 or ret_code > 0:
            #     raise OSError('Error executing the {} job locally.'.format(step_id))

            # list of all script paths

            print('JobValues')
            print(job_paths)
            print(job_paths.values())

            #add by adlofts Jun 26 2018
            jobs_list = list(job_paths.values())
            #pool = Pool(len(jobs_list))
            pool = Pool()

            #removed
            #jobs_list = job_paths.values()
            # pool = Pool()

            print('JobList')
            print(len(jobs_list))
            for idx_beg in range(0, len(jobs_list), num_procs):

                #pool = Pool(num_procs + 1)   #added

                print('Index')
                print(idx_beg)
                print('Added')
                print(num_procs)
                print('Iteration')
                print(jobs_list[idx_beg:idx_beg + num_procs])
                # running only on a specified num cores, one subject at a time
                pool.map(local_exec, jobs_list[idx_beg:idx_beg + num_procs])
                #print('Finished pool map')
                #pool.close()
                #print('Finished close')
                #pool.join()
                #print('Finished join')
            print('Finished pool map')
            pool.close()
            print('Finished close')
            pool.join()
            print('Finished join')
        else:
            for prefix, job_path in job_paths.items():
                job_id_list[prefix] = make_dry_run(job_path)
    else:
        # num_sub_per_job = 1 # ability to specify >1 subjects per job
        # job_id_list = map(submit_queue, jobs_list)
        txt_out = list()
        for prefix, (job_path, job_details) in job_paths.items():
            job_details_str = ' '.join(job_details)
            # job_id_list[prefix] = submit_queue(job_details_str, depends_on_step)
            job_id_list[prefix] = submit_queue(job_path, depends_on_step)
            txt_out.append('{} (job id: {})'.format(prefix, job_id_list[prefix]))

        try:
            tty_height, tty_width = subprocess.check_output(['stty', 'size']).split()
        except:
            tty_height, tty_width = 120, 80
        max_width = max(map(len,txt_out))
        num_sets  = int(math.floor( int(tty_width) / (max_width+4)))
        for idx in range(0, len(txt_out), num_sets):
            print('\t' + '\t'.join(txt_out[idx:min(idx+num_sets,len(txt_out))]))
            time.sleep(0.08)

    print('')

    return True, job_id_list


# noinspection PyUnusedLocal
def make_dry_run(cmd_str):
    """ Simple function to perform a dry run (indicated by a -ve job id)."""

    # print(cmd_str)
    job_id = random.randrange(1000)

    return -job_id


def submit_queue(job, depends_on_steps):
    """ Helper to submit jobs to the queue, taking care of the inter-dependencies. Returns the job ID."""
    global hpc

    qsub_path = which(hpc['spec']['submit_cmd'])  # find_executable('qsub')

    # encoding dependencies
    if depends_on_steps is not None:
        # making it a list when only one step is specified
        if not isinstance(depends_on_steps, list):
            depends_on_steps = [depends_on_steps, ]

        job_ids = list()
        for step in depends_on_steps:
            if step in hpc['job_ids_grouped']:
                # this condition can be False when rerunning from a existing processing
                # some steps may have been completed already - which means they are not in queue now,
                # in which case this current step doesnt need to wait
                job_ids.extend(hpc['job_ids_grouped'][step].values())

        if job_ids:  # if not empty
            job_id_list_str = map(str, job_ids)
        else:
            job_id_list_str = ''
    else:
        job_id_list_str = ''

    if hpc['type'].upper() in ('ROTMAN', 'BRAINCODE', 'ROTMAN-SGE', 'SGE'):
        # qsub_cmd  = qsub_path + ' -terse '
        qsub_cmd = qsub_path
        terse = '-terse'
        hold_spec = '-hold_jid ' + ",".join(job_id_list_str)
    elif hpc['type'].upper() in ('CAC', 'HPCVL', 'QUEENSU'):
        # qsub_cmd  = qsub_path + ' -terse '
        qsub_cmd = qsub_path
        terse = '-terse'
        hold_spec = '-hold_jid ' + ",".join(job_id_list_str)
    elif hpc['type'].upper() in ('FRONTENAC', 'SLURM'):
        qsub_cmd = qsub_path
        terse = '--parsable'
        hold_spec = '--dependency=afterok:' + ":".join(job_id_list_str)
    elif hpc['type'].upper() in ('SCINET', 'PBS', 'TORQUE'):
        qsub_cmd = qsub_path
        terse = ''
        hold_spec = '-W depend=' + ":".join(job_id_list_str)
    else:
        raise ValueError('HPC type {} unrecognized or not implemented.'.format(hpc['type']))

    if job_id_list_str not in [None, '', []]:
        hold_spec_str = hold_spec
    else:
        # some jobs don't need to wait for others!
        hold_spec_str = ''

    arg_list = [qsub_cmd, terse, hold_spec_str, job]
    # removing empty args (terse and hold list can be empty sometimes)
    arg_list = [arg for arg in arg_list if arg]
    full_cmd = ' '.join(arg_list)
    full_cmd = re.sub('\s+', ' ', full_cmd).strip()

    # splitting on space is required, otherwise subprocess throws an error,
    # as the '-opt value' in a single string gets interepreted as the name of an option
    arg_list = full_cmd.split(' ')

    if not hpc['dry_run']:
        job_id = subprocess.check_output(arg_list)
        job_id = job_id.strip().decode('utf-8')
    else:
        job_id = make_dry_run(full_cmd)

    return job_id


def submit_jobs():
    """
    Gateway to OPPNI preprocessing and optimization tool.
        Coordinates the status checks and runs/reruns what is required to perform the requested processing.
    """

    global hpc
    # check args
    unique_subjects, options, input_file, cur_garage, time_stamp, proc_out_dir = parse_args_check()

    if options.status_update_in is not None and unique_subjects in [[], None]:
        print('Status update done.')
        return

    # we dont need to check status if its a dry run meant to generate job scripts
    if not hpc['dry_run']:
        # check the status of processing
        # so processing can be done only for the unfinished or failed subjects
        # notice the inputs are combined as a list
        print('\n Running a status check on previous processing ...')
        proc_status, rem_input_file, rem_spnorm_file = check_proc_status.run(
            [input_file, options.pipeline_file, '--skip_validation', '--not_verbose'], options)
        for step in cfg_pronto.STEPS_PROCESSING_STATUS:
            if step is not 'all_done':
                if not getattr(proc_status, step):
                    bool_str = 'not done.'
                else:
                    bool_str = '    done.'
                print('{:>15} is {}'.format(step, bool_str))
        print(' ')

        if proc_status.all_done:
            print("All of preprocessing, optimization and QC seem to be finished already.")
            print(" if you'd like to force preprocessing, please rename/remove/move the existing outputs and rerun.")
            return
    else:
        print('This is just a dry run - generating jobs for all steps regardless of their processing status.')
        proc_status = cfg_pronto.initialize_proc_status()
        rem_input_file = None
        rem_spnorm_file = None

    if options.run_locally:
        hpc['type'] = "LOCAL"
    else:
        hpc['type'] = find_hpc_type(options.hpc_type)

    if options.part is 0:
        run_part_one = True
        run_part_two = True
        run_sp_norm = True
        run_gmask = True
        run_qc1 = True
        run_qc2 = True
    elif options.part is 1:
        run_part_one = True
        run_part_two = False
        run_sp_norm = False
        run_gmask = False
        run_qc1 = False
        run_qc2 = False
    elif options.part is 2:
        run_part_one = False
        run_part_two = True
        run_sp_norm = False
        run_gmask = False
        run_qc1 = False
        run_qc2 = False
    elif options.part is 3:
        run_part_one = False
        run_part_two = False
        # run spatial norm only when the reference is specified and exists
        if options.reference_specified:
            run_sp_norm = True
            run_gmask = True
        else:
            print('A reference atlas is not specified - skipping spatial normalization..')
            run_sp_norm = False
            run_gmask = False
        run_qc1 = False
        run_qc2 = False
    elif options.part is 4:
        run_part_one = False
        run_part_two = False
        run_sp_norm = False
        run_gmask = False
        run_qc1 = True
        run_qc2 = True

    if options.part in [1, 2, 4]:
        run_sp_norm = False
        if options.reference_specified:
            print('A reference atlas is specified although SPNORM is not requested - ignoring the atlas.')

    # submitting jobs for preprocessing for all combinations of pipelines
    if run_part_one and proc_status.preprocessing is False:
        # running part 1 only on subjects with incomplete processing
        print('Preprocessing:')
        status_p1, job_ids_pOne = run_preprocessing(unique_subjects, options, rem_input_file, cur_garage)
        if options.run_locally is True and (status_p1 is False or status_p1 is None):
            raise Exception('Preprocessing step failed.')

    spnorm_completed = False
    spnorm_step1_completed = False
    if options.dospnormfirst:
        sp_norm_step = 1
        print('spatial normalization (step 1) BEFORE optimization:')
        status_sp, job_ids_spn = process_spatial_norm(unique_subjects, options, rem_spnorm_file, sp_norm_step,
                                                      cur_garage)
        if options.run_locally is True and (status_sp is False or status_sp is None):
            raise Exception('Spatial normalization - step 1 failed.')
        else:
            spnorm_step1_completed = True

    # submitting jobs for optimization
    if options.analysis in [ "None", None, '' ]:
        print("WARNING: analysis model is NOT specified (specified by switch -a)")
        print("\tNO optimization will be performed, OPPNI will generate ONLY the preprocessed data.\n")
    elif run_part_two and proc_status.optimization is False:
        # optimization is done for ALL the subjects in the input file,
        # even though part 1 may have been rerun just for failed/unfinished subjects
        print('stats and optimization :')
        status_p2, job_ids_pTwo = run_optimization(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_p2 is False or status_p2 is None):
            raise Exception('Optimization failed.')
    else:
        raise ValueError('Unexpected or invalid state of flags in the wrapper.')

    # finishing up the spatial normalization
    if run_sp_norm:
        if spnorm_step1_completed is True:
            sp_norm_step = 2
        else:
            sp_norm_step = 0

        print('spatial normalization (step 2): ')
        status_sp = process_spatial_norm(unique_subjects, options, rem_spnorm_file, sp_norm_step, cur_garage)
        if options.run_locally is True and (status_sp is False or status_sp is None):
            raise Exception('Spatial normalization - steps 2 and later failed.')

    # group mask
    if proc_status.gmask is False and run_gmask is True:
        print('group mask generation: Submitting jobs ..')
        status_gm = process_group_mask_generation(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_gm is False or status_gm is None):
            raise Exception('Group mask generation failed.')

    # generating QC1 if not done already
    if proc_status.QC1 is False and run_qc1 is True:
        print('QC 1 :')
        status_qc1, job_ids_qc1 = run_qc_part_one(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_qc1 is False or status_qc1 is None):
            raise Exception('QC part 1 failed.')

    # generating QC2 if not done already
    if proc_status.QC2 is False and run_qc2 is True:
        print('QC 2 :')
        status_qc2, job_ids_qc2 = run_qc_part_two(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_qc2 is False or status_qc2 is None):
            raise Exception('QC part 2 failed.')

    # saving the job ids and hpc cfg to facilitate a status update in future
    save_hpc_cfg_and_jod_ids(cur_garage)


def save_hpc_cfg_and_jod_ids(cur_garage):
    """Saves the job ids and hpc cfg to facilitate a status update in future"""

    # saving the job ids
    job_id_file = os.path.join(cur_garage, file_name_job_ids_by_group)
    if os.path.isfile(job_id_file):
        os.remove(job_id_file)
    with open(job_id_file, 'w') as jlist:
        json.dump(hpc['job_ids_grouped'], jlist, indent=2)

    # saving the config
    cfg_file = os.path.join(cur_garage, file_name_hpc_config)
    if os.path.isfile(cfg_file):
        os.remove(cfg_file)
    with open(cfg_file, 'w') as hcf:
        json.dump(hpc, hcf, indent=2)


if __name__ == '__main__':
    submit_jobs()
