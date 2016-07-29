#!/usr/bin/env python

import logging, pwd, time
import os, sys, subprocess, stat, re, random, argparse
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool
from collections import namedtuple, OrderedDict
from recordclass import recordclass
from distutils.spawn import find_executable
from time import localtime, strftime
from time import sleep
from os import listdir
from copy import copy
from os.path import isdir, join
from shutil import rmtree
import shutil
import pickle
import tempfile

# PRONTO related
import cfg_front as cfg_pronto
import proc_status_front as check_proc_status
import sgeparse

file_name_hpc_config = 'hpc_config.pkl'
file_name_job_ids_by_group = 'pronto_job_ids_by_step_prefix.pkl'
file_name_prev_options = 'pronto-options.pkl'

# descriptive variables
NOT_DONE = False
DONE = True

global hpc
hpc = {'type': 'SGE',
       'shell' : '/bin/bash',
       'prefix': '#$',
       'spec': None,
       'header': '',
       'dry_run': False,
       'hold_jobid_list': [],
       'job_ids_grouped': {}}

# defining regexes that may be useful in various functions
# regex to extract the relevant parts of the input file
reIn = re.compile(r"IN=([\w\./_-]+)[\s]*")
reOut = re.compile(r"OUT=([\w\./_-]+)[\s]*")
reDrop = re.compile(r"DROP=\[(\d+),(\d+)\][\s]*")
reTask = re.compile(r"TASK=([\w\./_-]+)[\s]*")

rePhysio = re.compile(r"PHYSIO=([\w\./_-]+)[\s]*")
reStruct = re.compile(r"STRUCT=([\w\./_-]+)[\s]*")
reCustReg = re.compile(r"CUSTOMREG=([\w\./_-]+)[\s]*")

# for the pipeline file
rePip = re.compile('([0-9A-Z\s]+)=.+', re.IGNORECASE)
rePip2 = re.compile(r'([0-9A-Z\s]+)=\[([\d,]*)\][\s]*')


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
            steps_list_file.append(line.rstrip(' '))

    steps_spec = rePip2.findall(' '.join(steps_list_file).upper())
    # steps_no_spec = map( lambda str1: str1.strip(' \n'), steps_no_spec)
    steps_dict = {}
    for step in steps_spec:
        if not step[0] in cfg_pronto.CODES_PREPROCESSING_STEPS:
            print('Error in pipeline file: %s' % pipeline_file)
            raise TypeError('Unrecognized pipeline step: %s' % step[0])
        steps_dict[step[0]] = map(int, step[1].replace(',', ''))

    print('  Done.')

    return steps_dict


def validate_task_file(task_path):
    "Basic validation of task spec file."

    with open(task_path, 'r') as tf:
        task_spec = tf.read().splitlines()
        # task_spec = [ line.strip('\n ') for line in task_spec ]
    task_spec = '\n'.join(task_spec)

    for field in cfg_pronto.TASK_MANDATORY_FIELDS:
        if field + '=' not in task_spec:
            raise TypeError('{} is not defined in task file'.format(field))

    return True


def validate_env_var(var):
    assert os.getenv(var) is not None, "Path {} is not defined. Fix your environment and rerun.".format(var)


def validate_user_env(opt):
    if not hpc['dry_run']:
        for var in ['AFNI_PATH', 'FSL_PATH']:
            validate_env_var(var)

        if opt.environment.lower() in ['compiled']:
            validate_env_var('MCR_PATH')


def validate_input_file(input_file, options=None, new_input_file=None):
    if new_input_file is None or options is None:
        new_file = tempfile.TemporaryFile()
        # in case of resubmission, this should not append additional layer
        cur_suffix = None
    else:
        new_file = open(new_input_file, 'w')
        cur_suffix = options.suffix

    unique_subjects = OrderedDict()

    line_count = 0
    dupl_count = 0
    with open(input_file, 'r') as ipf:
        for line in ipf.readlines():
            subject, new_line = validate_input_line(line, cur_suffix)
            if subject is False or subject is None:
                raise ValueError('Error in line number {}.'.format(line_count))
            else:
                line_count += 1
                subject['line'] = new_line
                # if the key doesnt exist, dict returns None
                if unique_subjects.get(subject['prefix']) is not None:
                    print "Potential duplicate prefix in line {}: {} \t Previously processed line contained this prefix.".format(line_count, subject['prefix'])
                    dupl_count += 1
                unique_subjects[subject['prefix']] = subject
                new_file.write(new_line)

            # options = None is when this helper script is called from outside of fRONT
            if options is not None and options.contrast_specified:
                assert subject['task'] is not None, \
                    'Contrast specified, but not a task file in line number {}.'.format(line_count)
            if options is not None and options.reference_specified:
                assert subject['struct'] is not None, \
                    'Reference atlas is specified, but not a structural scan in line number {}.'.format(line_count)
            # if one of the physiological methods are requested
            if options is not None and options.physio_correction_requested:
                assert subject['physio'] is not None, \
                    'RETROICOR and/or PHYPLUS are specified but not the physiological files! Line number {}.'.format(
                        line_count)
            if options is not None and options.custom_mask_requested:
                assert subject['mask'] is not None, \
                    'CUSOMREG is specified but not a binary mask! Line number {}.'.format(line_count)

    new_file.close()

    # if the optional fields are specified, making sure they are specified for all the subjects
    optional_fields = ['task', 'struct', 'physio']
    for optf in optional_fields:
        bool_all_subjects = [sub['task'] is not None for ix, sub in unique_subjects.items()]
        num_true = sum([1 for ii in bool_all_subjects if ii is True])
        if not (num_true == 0 or num_true == len(bool_all_subjects)):
            print('Optional fields can either be specified for ALL the subjects, or NONE at all. ')
            raise ValueError(
                ' ---> {0} specified only for {1} out of {2} subjects'.format(optf, num_true, len(bool_all_subjects)))

    print('Parsed {0} lines without error.'.format(line_count))
    if dupl_count > 0:
        print('\t excluding {} duplicate prefixes: {} lines'.format(dupl_count, len(unique_subjects)))
    return unique_subjects


def validate_input_line(ip_line, suffix=''):
    line = ip_line.strip()
    LINE = line.upper()

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
        print "IN= section not defined."
        return (False, "")
    else:
        nii = reIn.search(line).group(1)
        if not os.path.isfile(nii):
            print "Input file not found: " + nii
            return (False, "")
        else:
            subject['nii'] = nii

    # OUT part
    if ("OUT=" not in LINE):
        print "OUT= section not defined."
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
        new_line = re.sub(base_out_dir, subject['out'], ip_line)

        subject['prefix'] = os.path.basename(out)
        # Parse_Input_File.m is not robust with parsing e.g. an extra / at the end will mess up everything
        # so checking to make sure prefix is not empty
        assert (subject['prefix'] not in [None, '']), \
            'subject prefix can not be empty!'

    # DROP part
    if "DROP=" not in LINE:
        print "DROP= section not defined."
        return (False, "")
    else:
        idx = reDrop.search(line).groups()
        if (not len(idx) is 2):
            print "Atleast two indices must be specified to drop from the start and end."
            return (False, "")
        elif (not all([ii.isdigit() for ii in idx])):
            print "DROP indices must be numeric!"
            return (False, "")
        elif (not all([ii >= 0 for ii in idx])):
            print("All the DROP indices must be non-negative.")
        else:
            subject['drop_beg'] = idx[0]
            subject['drop_end'] = idx[1]

    # the following parts are not mandatory, errors will be raised later on, when inconsistencies are found.
    # TASK part
    if ("TASK=" in line):
        task = reTask.search(line).group(1)
        if not os.path.isfile(task):
            print "Task file " + task + " not found."
        else:
            validate_task_file(task)
            subject['task'] = task

    # PHYSIO part
    if ("PHYSIO=" in LINE):
        physio = rePhysio.search(line).group(1)
        if not os.path.isfile(physio + '.puls.1D') or not os.path.isfile(physio + '.resp.1D'):
            print "PHYSIO files (puls and/or resp) at " + physio + " not found."
        else:
            subject['physio'] = physio

    # STRUCT part
    if ("STRUCT=" in LINE):
        struct = reStruct.search(line).group(1)
        if not os.path.isfile(struct):
            print "STRUCT file " + struct + " not found."
        else:
            subject['struct'] = struct

    # PHYSIO part
    if ("CUSTOMREG=" in LINE):
        mask = reCustReg.search(line).group(1)
        if not os.path.isfile(mask):
            print "Binary mask " + mask + " not found."
        else:
            subject['mask'] = mask

    return subject, new_line


def parse_args_check():
    # parse different input flags
    parser = argparse.ArgumentParser(prog="pronto")

    parser.add_argument("-s", "--status", action="store", dest="status_update_in",
                        default=None,
                        help="Performs a status update on previously submitted processing. "
                             " Supply the output directory where the previous processing has been stored in.")

    parser.add_argument("--validate", action="store", dest="val_input_file_path",
                        default=None,
                        help="Performs a basic validation of an input file (existence of files and consistency!).")

    parser.add_argument("-p", "--part", action="store", dest="part", type=int,
                        default=0, choices=[0, 1, 2, 3],
                        help="select pipeline optimization step, 1: Estimation step, 2: Optimization step, "
                             "3: Spatial normalization, 0: All three steps (default 0)")
    parser.add_argument("-i", "--input_data", action="store", dest="input_data_orig",
                        help="FILE.txt contains the input and output data paths", metavar="input spec file")
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
                        help="Pptimization metric")

    # parser.add_argument("--autodetect", action="store_true", dest="autodetect",
    #                    help="Automatically detect subjects and optimize for each subject independently, "
    #                         "the lines in input files that have same structrul image (STRUCT) and "
    #                         "output directory (OUT) are considered a subject")

    parser.add_argument("-r", "--reference", action="store", dest="reference",
                        default=None,
                        help="anatomical reference to be used in the spatial normalization step, i.e. -p,--part=3")
    parser.add_argument("--dospnormfirst",
                        action="store_true", dest="dospnormfirst", default=False,
                        help="First normalize the data to a reference (specified by switch -r), then perform the preprocessing optimization.")
    parser.add_argument("--contrast", action="store", dest="contrast_list_str",
                        default="None",
                        help="desired task contrast; This argument is necessary when more than two contrasts are "
                             "defined in the split info files .. "
                             "Syntax: --Contrast \"CON1 vs CON2,CON2 vs CON3\". Output will be a consensum ouput "
                             "from all the contrasts.")
    # TODO option to make the contrasts separable.

    parser.add_argument("-k", "--keepmean", action="store", dest="keepmean",
                        default="1",
                        help="(optional) determine whether the ouput nifti files contain the mean scan "
                             "(Default keepmean=1, i.e. not remove the mean)")
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
                        help="FRACTION=Scalar value of range (0,1), indicating the fraction of full-date PCA subspace to keep during PCA-LDA analysis. "
                             "A drf of 0.3 is recommended as it has been found to be optimal in previous studies.",
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
                        help="Specify the number of resamples for multi-run analysis")
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
                        help="Whether to output spms for all the optimal pipelines.")

    # HPC related
    parser.add_argument("-e", "--environment", action="store", dest="environment",
                        default="compiled",
                        choices=('matlab', 'compiled'),
                        help="(optional) determine which software to use to run the code: matlab or compiled(default)")
    parser.add_argument("--memory", action="store", dest="memory",
                        default="4",
                        help="(optional) determine the minimum amount RAM needed for the job, e.g. --memory 8G (G==gigabytes)!")
    parser.add_argument("-n", "--numcores", action="store", dest="numcores",
                        default=1,
                        help="(optional) number of threads used for the process (not allowed for some SGE systems)")
    parser.add_argument("-q", "--queue", action="store", dest="queue",
                        default="abaqus.q",
                        help="(optional) SGE queue name, default is bigmem_16.q")

    parser.add_argument("--cluster", action="store", dest="hpc_type",
                        default='LOCAL', choices=('LOCAL', 'SGE', 'SUNGRID', 'TORQUE', 'PBS'),
                        help="Please specify the type of cluster you're running the code on. SGE/SUNGRID mean the same as TORQUE/PBS.")
    parser.add_argument("--run_locally", action="store_true", dest="run_locally",
                        default=False,
                        help="Run the pipeline on this computer without using SGE. "
                             "Specify the number of cores using -n (or --numcores) to allow the program to run in pararallel")
    parser.add_argument("--noSGE", action="store_true", dest="run_locally",
                        default=False,
                        help="Same as --run_locally. Retained for backward compatibility.")
    parser.add_argument("--numprocess", action="store", dest="numprocess",
                        default=1,
                        help="When running the pipeline wihout using SGE, this switch specifies the number of simultaneous processes")

    parser.add_argument("--force_rerun", action="store_true", dest="force_rerun",
                        default=False,
                        help="DANGER: Cleans up existing processing and forces a rerun. Be absolutely sure this is what you want to do.")
    parser.add_argument("--dry_run", action="store_true", dest="dry_run",
                        default=False,
                        help="Generates job files only, but not run/sbumit them. Helps in debugging the queue/HPC options.")

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

    # updating the status to the user if requested.
    if options.status_update_in is not None:
        cur_garage = os.path.abspath(options.status_update_in)
        update_status_and_exit(cur_garage)
    elif options.val_input_file_path is not None:
        # performing a basic validation
        try:
            unq_sub = validate_input_file(options.val_input_file_path)
            print " validation succesful."
        except:
            print " validation failed."
        sys.exit(0)

    global hpc

    # sanity checks on HPC inputs and obtaining the cfg of hpc
    hpc['type'] = options.hpc_type
    if options.run_locally is False:
        hpc['type'] = find_hpc_type(options.hpc_type, options.run_locally)
        hpc['spec'], hpc['header'], hpc['prefix'] = get_hpc_spec(hpc['type'], options)
    else:
        if not hpc['type'] == 'LOCAL':
            raise ValueError('Conflicting options specified: specify either of --run_locally or --cluster CLUSTERTYPE.')

    if hpc['type'] in (None, 'LOCAL'):
        if options.run_locally == False:
            print "Sun grid engine (SGE) has not been detected!"
            print "Use --run_locally switch if you want to run PRONTO without HPC cluster on your computer locally."
            exit(1)
        else:
            print "Running jobs to the current node"
            print "The code will wait until the jobs are finished."
            # even though it will be run locally, we will generate the job file
            # which will be executed in a subshell
            hpc['type'] = 'SGE'
    else:
        print "Submitting jobs to Sun Grid Engine (SGE)"

    if options.dry_run:
        hpc['dry_run'] = True

    # validating the pipeline file
    steps_dict = validate_pipeline_file(options.pipeline_file)
    setattr(options, 'pipeline_steps', steps_dict)

    contrast_specified = options.contrast_list_str != 'None'
    reference_specified = options.reference != None
    physio_correction_requested = (1 in options.pipeline_steps['RETROICOR']) or (1 in options.pipeline_steps['PHYPLUS'])
    custom_mask_requested = 1 in options.pipeline_steps['CUSTOMREG']

    # saving useful flags for sanity checks on completeness of inputs later
    setattr(options, 'contrast_specified', contrast_specified)
    setattr(options, 'reference_specified', reference_specified)
    setattr(options, 'physio_correction_requested', physio_correction_requested)
    setattr(options, 'custom_mask_requested', custom_mask_requested)

    cur_garage, time_stamp, proc_out_dir, suffix = organize_output_folders(options)
    setattr(options, 'out_dir_common', cur_garage)
    setattr(options, 'suffix', suffix)

    if options.dospnormfirst and not options.reference_specified:
        raise ValueError('Spatial normalization requested, but a reference atlas is not specified.')

    # ensuring the validity of input file, given the options and pipeline file
    # and compiling a list of unique subjects into a new file with output folder modified
    new_input_file = os.path.join(cur_garage, 'input_file.txt')
    unique_subjects = validate_input_file(options.input_data_orig, options, new_input_file)

    # making sure user environment is properly setup before even submitting jobs
    validate_user_env(options)

    ## -------------- Check the Input parameters  --------------


    if hasattr(options, 'reference') and options.reference is not None:
        reference = os.path.abspath(options.reference)
        if not os.path.isfile(reference):
            raise IOError('Reference file supplied does not exist!')
    else:
        options.reference = ""

    if options.analysis is None or options.analysis == "None":
        print "WARNING: without an analysis model (specified by switch -a), no optimization will be performed"
        print "  PRONTO will only generate the preprocessed data"
        options.contrast_list_str = "None"

    if hasattr(options, 'DEOBLIQUE') and options.DEOBLIQUE:
        options.DEOBLIQUE = "1"
    else:
        options.DEOBLIQUE = "0"

    ## --------------  Checking the switches --------------
    analysis = options.analysis
    if (analysis.upper() == "LDA") and (options.drf == "None"):
        print "WARNING (Deprecated usage): --drf switch not defined for LDA model. PRONTO will check TASK files for parameter(s)"

    if (analysis.upper() == "ERCVA") and (options.drf == "None"):
        print "WARNING (Deprecated usage): --drf switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "ERCVA") and (options.Nblock == "None"):
        print "WARNING (Deprecated usage): --Nblock switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "ERCVA") and (options.WIND == "None"):
        print "WARNING (Deprecated usage): --WIND switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "ERCVA") and (options.subspace == "None"):
        print "WARNING (Deprecated usage): --subspace switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"

    if (analysis.upper() == "GNB") and (options.decision_model == "None"):
        print "WARNING (Deprecated usage): --decision_model switch not defined for GNB model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "ERGNB") and (options.Nblock == "None"):
        print "WARNING (Deprecated usage): --Nblock switch not defined for erGNB model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "ERGNB") and (options.WIND == "None"):
        print "WARNING (Deprecated usage): --WIND switch not defined for erGNB model. PRONTO will check TASK files for parameter(s)"

    if (analysis.upper() == "SCONN") and (options.spm == "None"):
        print "WARNING (Deprecated usage): --spm switch has to be used with the SCONN model. PRONTO will check TASK files for parameter(s)"

    if (analysis.upper() == "GLM") and (options.convolve == "None"):
        print "WARNING (Old style usage): --convolve switch has to be used with the GLM model. PRONTO will check TASK files for parameter(s)"
    if (analysis.upper() == "GPCA") and (options.num_PCs == "None"):
        print "WARNING (Deprecated usage): --num_PCs switch not defined for gPCA model. PRONTO will check TASK files for parameter(s)"

    if not (options.convolve in ["1", "0", "None"]):
        print "WARNING (Deprecated usage): --convolve has to be 0 or 1"
    if not (options.decision_model.lower() in ["linear", "nonlinear", "none"]):
        print "WARNING (Deprecated usage): --decision_model has to be linear or nonlinear"
    if not (options.subspace.lower() in ["onecomp", "multicomp", "none"]):
        print "WARNING (Deprecated usage): --subspace has to be onecomp or multicomp"
    if not (options.spm.lower() in ["corr", "zscore", "none"]):
        print "WARNING (Deprecated usage): --spm has to be corr or zscore"

    options.model_param_list_str = "keepmean " + options.keepmean \
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
    print options

    saved_cfg_path = os.path.join(cur_garage, file_name_prev_options)
    with open(saved_cfg_path, 'wb') as cfg:
        pickle.dump([unique_subjects, options, new_input_file, cur_garage], cfg)

    return unique_subjects, options, new_input_file, cur_garage, time_stamp, proc_out_dir


def organize_output_folders(options):
    parent_dir_wrapper = os.path.dirname(os.path.abspath(sys.argv[0]))

    proc_out_dir = get_out_dir_first_line(options.input_data_orig)
    if not os.path.exists(proc_out_dir):
        os.makedirs(proc_out_dir)

    # TODO need a way to pass the suffix to matlab core so it can reuse preprocessing for multiple analysis models and contrasts.
    suffix = os.path.splitext(os.path.basename(options.input_data_orig))[0]
    if options.analysis is not "None":
        suffix = suffix + '_' + options.analysis

    if options.contrast_specified:
        # task1-fixation,task2-task1 will be changed to task1-fixation_task2-task1
        suffix = suffix + '_' + options.contrast_list_str.replace(',', '_')

    time_stamp = make_time_stamp()
    # suffix = "garage_{0}_{1}".format(time_stamp,suffix)
    suffix = "processing_{0}".format(suffix)

    cur_garage = os.path.join(proc_out_dir, suffix)
    if os.path.exists(cur_garage):
        if options.force_rerun:
            user_confirmation = raw_input("Are you sure you want to delete previous results? (y/[N])")
            if user_confirmation.lower() in ['y', 'yes', 'ye']:
                print('Removing any existing preprocessing, as requested!')
                rmtree(cur_garage)
                os.mkdir(cur_garage)
            else:
                print('Leaving the existing preprocessing as is!')
    else:
        os.mkdir(cur_garage)

    print "Output processing folder: " + cur_garage

    return cur_garage, time_stamp, proc_out_dir, suffix


def find_hpc_type(user_supplied_type, run_locally=False):
    if user_supplied_type is None and run_locally is False:
        if find_executable('qsub') is None:
            raise SystemError('Requested processing on SGE, but command qsub is not found!')

        sge_root = os.getenv("SGE_ROOT")
        if sge_root is not None and find_executable('qconf') is not None:
            type = 'SGE'
        elif find_executable('qmgr') is not None and \
                (find_executable('qconf') is None and find_executable('sbatch') is None):
            type = 'TORQUE'
        elif find_executable('sbatch') is not None and find_executable('srun') is not None:
            type = 'SLURM'
        else:
            raise ValueError('Can not determine the type of cluster you are one.\n Please specify it with --cluster')
    else:
        type = user_supplied_type.upper()

    return type


def get_hpc_spec(type=None, options=None):

    if type is None:
        type = find_hpc_type()

    if options is not None:
        memory = options.memory
        numcores = options.numcores
        queue = options.queue
    else:
        memory = '2'
        numcores = 1
        queue = 'all.q'

    spec = {'memory': None, 'numcores': None, 'queue': None}

    if type.upper() in ('ROTMAN-SGE', 'SGE', 'ROTMAN'):
        prefix = '#$'
        spec['memory'] = ('-l h_vmem=', memory + 'g')
        spec['numcores'] = ('-pe npairs ', numcores)
        spec['queue'] = ('-q ', queue)
    elif type.upper() in ('HPCVL', 'QUEENSU'):
        prefix = '#$'
        spec['memory'] = ('-l mf=', memory + 'G')
        spec['numcores'] = ('-pe shm.pe ', numcores)
        spec['queue'] = ('-q ', queue)
    elif type.upper() in ('SCINET', 'PBS', 'TORQUE'):
        prefix = '#PBS'
        spec['memory'] = ('-l mem=', memory)
        spec['numcores'] = ('-l ppn=', numcores)
        spec['queue'] = ('-q ', queue)
    elif type.upper() in ('SLURM'):
        prefix = '#SBATCH'
        spec['memory'] = ('--mem=', memory)
        spec['numcores'] = ('--n ', numcores)
        spec['queue'] = ('-p ', queue)
    else:
        raise ValueError('HPC type {} unrecognized or not implemented.'.format(type))

    header = list()
    header.append('{0} -V'.format(prefix))
    header.append('{0} -b y'.format(prefix))
    header.append('{0} -j y'.format(prefix))
    for key, val in spec.items():
        if key == 'numcores' and int(numcores) == 1:
            continue
        else:
            header.append('{0} {1}{2}'.format(prefix, val[0], val[1]))

    # not joining them for later use
    # header = '\n'.join(header)

    return spec, header, prefix


def update_status_and_exit(out_dir):
    assert os.path.exists(out_dir), \
        "The processing directory for which the status update has been requested doesn't exist!"

    # print('Queue status update requested ... ')
    # update_Q_status(out_dir)
    print('\nNow checking the outputs on disk ...')
    try:
        prev_proc_status, prev_options, prev_input_file_all, failed_sub_file, all_subjects = update_proc_status(out_dir)

        if not prev_proc_status.all_done:
            print('Previous processing is incomplete.')
            user_confirmation = raw_input("Would you like to resubmit jobs for failed subjects/runs? (y/[N])")
            if user_confirmation.lower() in ['y', 'yes', 'ye']:
                print('Attempting resubmission ... ')
                try:
                    reprocess_failed_subjects(prev_proc_status, prev_options, failed_sub_file,
                                              prev_input_file_all, all_subjects, out_dir)
                    print('Resubmission done. Check the status later.')
                except:
                    print('Resubmission failed. Try again later.')
                    raise
    except:
        print('Unable to resubmit.')
        raise

    # to avoid entering a recursive state
    sys.exit(0)


def update_proc_status(out_dir):
    opt_file = os.path.join(out_dir, file_name_prev_options)

    with open(opt_file, 'rb') as of:
        all_subjects, options, new_input_file, _ = pickle.load(of)
        proc_status, failed_sub_file = check_proc_status.run(
            [new_input_file, options.pipeline_file, '--skip_validation'])
    return proc_status, options, new_input_file, failed_sub_file, all_subjects


def update_Q_status(out_dir):
    if find_executable('qstat') is None:
        print('\t qstat is not found - perhaps you are not on a headnode for the cluster. \n'
              '\t Unable to query queue statuses. Run it on headnode to get the complete report.')
        return

    scheduler = drmaa.Session()
    scheduler.initialize()

    job_id_file = os.path.join(out_dir, file_name_job_ids_by_group)
    try:
        with open(job_id_file, 'rb') as jobs_list:
            jobs_list_grouped = pickle.load(jobs_list)

            for step, id_dict in jobs_list_grouped.items():
                print('  --> {}'.format(step))
                num_jobs_rqh, num_jobs_err = q_status(scheduler, id_dict)

    except IOError as ioe:
        print('Trouble reading the job IDs saved from last session. Check queue status manually. \n {}'.format(ioe))

    # exiting the session
    scheduler.exit()


def q_status(scheduler, job_ids):
    "Display the status of each job, and count the number of jobs in error/Okay state"

    # decode_status = {drmaa.JobState.UNDETERMINED: 'status cannot be determined',
    #                  drmaa.JobState.QUEUED_ACTIVE: 'queued and active',
    #                  drmaa.JobState.SYSTEM_ON_HOLD: 'queued and in system hold',
    #                  drmaa.JobState.USER_ON_HOLD: 'queued and in user hold',
    #                  drmaa.JobState.USER_SYSTEM_ON_HOLD: 'queued and in user and system hold',
    #                  drmaa.JobState.RUNNING: 'running',
    #                  drmaa.JobState.SYSTEM_SUSPENDED: 'suspended by system',
    #                  drmaa.JobState.USER_SUSPENDED: 'suspended by user',
    #                  drmaa.JobState.DONE: 'finished normally',
    #                  drmaa.JobState.FAILED: 'finished, but failed'}

    num_jobs_rqh = 0  # num. jobs running or in queue
    num_jobs_err = 0  # num. jobs in error state

    for prefix, job_id in job_ids.items():
        try:
            job_state = scheduler.jobStatus(str(job_id))
        except drmaa.InvalidJobException:
            # when jobid has been flushed out of the scheduler memory!
            job_state = drmaa.JobState.UNDETERMINED
        else:
            print "error quering the job status"
            raise

        # job = sgeparse.job_status(job_id)
        # if 'state' in job:
        #     job_state = job['state'].lower()
        # else:
        #     job_state = 'running or unknown'

        print('\t {}: {} '.format(prefix, job_state))

        # if job_state != 'unknown':
        #     if 'e' in job_state:
        #         num_jobs_err += 1
        #
        #     if 'r' in job_state or 'q' in job_state:
        #         num_jobs_rqh += 1

    # if num_jobs_rqh > 0:
    #     print(' {} jobs are still either running, held or waiting to be processing.'.format(num_jobs_rqh))
    #
    # if num_jobs_err > 0:
    #     print('{} jobs are stuck in error state. Manage them manually'.format(num_jobs_err))

    return num_jobs_rqh, num_jobs_err

def make_job_file_and_1linecmd(file_path):
    "Generic job file generator, and returns one line version too."

    hpc_directives = list()
    job_name = os.path.splitext(os.path.basename(file_path))[0]
    if not hpc['type'].upper() == "LOCAL":
        # hpc_directives.append('{0} -S '.format(hpc['shell']))
        hpc_directives.extend(hpc['header'])
        hpc_directives.append('{0} -N {1}'.format(hpc['prefix'], job_name))
        hpc_directives.append('{0} -wd {1}'.format(hpc['prefix'], os.path.dirname(file_path)))
    else:
        # for jobs to run locally, no hpc directives are needed.
        hpc_directives.append('{0} -S '.format(hpc['shell']))

    with open(file_path, 'w') as jID:
        # one directive per line
        jID.write('\n'.join(hpc_directives))

    st = os.stat(file_path)
    os.chmod(file_path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    return hpc_directives


def make_job_file(file_path):
    "Generic job file generator"

    hpc_directives = list()
    job_name = os.path.splitext(os.path.basename(file_path))[0]
    if not hpc['type'].upper() == "LOCAL":
        if hpc.get('spec') is None:
            hpc['spec'], hpc['header'], hpc['prefix'] = get_hpc_spec(hpc['type'])

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
    st = os.stat(script_path)
    os.chmod(script_path, st.st_mode | stat.S_IXGRP)

    proc = subprocess.Popen(script_path, shell=True, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # communicate waits for the subprocess to finish
    std_output, _ = proc.communicate()
    # # print outputs and logs
    # logger.info('\n%s\n', std_output)
    print std_output

    return -random.randrange(1000)  # proc.returncode


def reprocess_failed_subjects(prev_proc_status, prev_options, failed_sub_file, prev_input_file_all, all_subjects,
                              garage):
    global hpc

    try:
        hpc_cfg_file = os.path.join(garage, file_name_hpc_config)
        with open(hpc_cfg_file, 'rb') as hpc_f:
            hpc = pickle.load(hpc_f)
    except:
        print('hpc config was not saved previously or proerly!')
        raise

    hpc['dry_run'] = False

    if prev_proc_status.preprocessing is NOT_DONE:
        failed_subjects = validate_input_file(failed_sub_file, prev_options, None)
        # running the failed subjects throught part 1
        print('Resubmitting preprocessing jobs .. ')
        status_p1, jobs_p1 = run_preprocessing(failed_subjects, prev_options, failed_sub_file, garage)

    if prev_proc_status.optimization is NOT_DONE:
        # but optimization will be done entire dataset
        print('Resubmitting stats & optim. jobs .. ')
        status_p2, jobs_p2 = run_optimization(all_subjects, prev_options, prev_input_file_all, garage)

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
    with open(job_id_file, 'wb') as jlist:
        pickle.dump(hpc['job_ids_grouped'], jlist)


def run_preprocessing(subjects, opt, input_file, garage):
    """
    Generates a job script to run the preprocessing for all combinations of pipeline steps requested.
    """

    job_id_list = None

    if opt.dospnormfirst:
        str_dospnormfirst = '1'
    else:
        str_dospnormfirst = '0'

    # matlab: Pipeline_PART1(InputStruct, input_pipeset, analysis_model, modelparam, niiout, contrast_list_str, dospnormfirst, DEOBLIQUE, TPATTERN, TOFWHM)
    # input file will be prepended in the process module
    arg_list = [opt.pipeline_file, opt.analysis, opt.model_param_list_str, opt.output_nii_also,
                opt.contrast_list_str, str_dospnormfirst, opt.DEOBLIQUE, opt.TPATTERN, opt.BlurToFWHM]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'PART1', 'Pipeline_PART1', arg_list, garage, None)

    return proc_status, job_id_list


def run_qc_part_one(subjects, opt, input_file, garage):
    # maskname is set to be empty
    arg_list = [1, input_file, 'None', opt.num_PCs]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'QC1', 'QC_wrapper', arg_list, garage, 'PART1')

    return proc_status, job_id_list


def run_qc_part_two(subjects, opt, input_file, garage):
    arg_list = [2, input_file, '', opt.num_PCs]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'QC2', 'QC_wrapper', arg_list, garage, 'PART2')

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

    arg_list = [input_file, opt.metric, mot_gs_control, opt.output_all_pipelines, opt.keepmean]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'PART2', 'Pipeline_PART2', arg_list, garage,
                                                      'PART1')

    return proc_status, job_id_list


def process_spatial_norm(subjects, opt, input_file, sp_norm_step, garage):
    """
    Generates a job script to register all the MRI's of given subjects to a reference.

    :returns: status of processing and a list of job IDs (or process IDs if running locally).
    """

    arg_list = [opt.reference, opt.voxelsize, sp_norm_step, opt.DEOBLIQUE]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'SPNORM', 'spatial_normalization', arg_list,
                                                      garage, None)

    return proc_status, job_id_list


def process_group_mask_generation(subjects, opt, input_file, garage):
    arg_list = [input_file]
    proc_status, job_id_list = process_module_generic(subjects, opt, 'GMASK', 'group_mask_tissue_maps', arg_list,
                                                      garage, 'SPNORM')

    return proc_status, job_id_list


def construct_full_cmd(environment, step_id, step_cmd_matlab, arg_list, prefix = None, job_dir = None):
    if environment.lower() in ('matlab', 'octave'):
        single_quoted = lambda s: r"'{}'".format(s)

        cmd_options = ', '.join(map(single_quoted, arg_list))

        # make an m-file script
        # hyphen/dash can be treated as an operator
        mfile_name = prefix.replace('-','_')
        mfile_path = os.path.join(job_dir, mfile_name + '.m')
        with open(mfile_path, 'w') as mfile:
            mfile.write('\n')
            mfile.write("try, {0}({1}); catch ME, display(ME.message); exit(1); end; exit;".format(step_cmd_matlab, cmd_options))
            mfile.write('\n')

        # # -nojvm omitted to allow future comptibility
        # cmd = 'time -p matlab -nodesktop -nosplash'
        # full_cmd = "{0} -r \"try, {1}({2}); catch ME, display(ME.message); exit(1); end; exit;\" ".format(
        #     cmd, step_cmd_matlab, cmd_options)

        # to simplify the trouble with quotes and escape sequences etc
        full_cmd = r"time -p matlab -nodesktop -nosplash -r {0} ".format(mfile_name)

    elif environment.lower() in ('standalone', 'compiled'):
        # first single quote is bash to protect from expansion or globbing
        #   second double quote is for matlab to interpret
        # strong_quoted = lambda s: r""" "'{}'" """.format(s)
        strong_quoted = lambda s: r""" '{}' """.format(s)

        # compiled_exe_path = find_executable('PRONTO')
        compiled_exe_path = 'run_cPRONTO.sh ' + os.getenv('MCR_PATH')

        cmd = 'time -p {0} {1} '.format(compiled_exe_path, step_id)
        cmd_options = ' '.join(map(strong_quoted, arg_list))
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
        prefix = '{0}-all-subjects'.format(step_id.lower())
        # each item will be a tuple (job_path, job_str)
        jobs_dict['all-subjects'] = make_single_job(opt.environment, step_id, step_cmd_matlab, prefix, arg_list_subset, job_dir)

    else:
        for idx, subject in enumerate(subjects.itervalues()):
            # subject-wise processing
            # TODO input file doesnt change with step, try refactoring this to have only one input file per run/subject
            prefix = '{1}_s{0:0>3}_{2}'.format(idx + 1, step_id.lower(), subject['prefix'])
            subset_input_file = os.path.join(input_dir, prefix + '.input.txt')
            with open(subset_input_file, 'w') as sif:
                sif.write(subject['line'])

            # adding the input file as the first arg
            arg_list_subset = [subset_input_file] + arg_list
            # each item will be a tuple (job_path, job_str)
            jobs_dict[subject['prefix']] = make_single_job(opt.environment, step_id, step_cmd_matlab, prefix, arg_list_subset, job_dir)

    jobs_status, job_id_list = run_jobs(jobs_dict, opt.run_locally, int(opt.numcores), depends_on_step)
    # storing the job ids by group to facilitate a status update in future
    hpc['job_ids_grouped'][step_id] = job_id_list

    return jobs_status, job_id_list


def make_single_job(environment, step_id, step_cmd_matlab, prefix, arg_list_subset, job_dir):
    # making a standalone job file
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
        quote_del = [ str.replace(r"'",'') for str in hpc_dir_2 ]
        jID.write('\n'.join(quote_del))

    # complete single line command for qsub
    qsub_opt = hpc_dir_1 + hpc_dir_2
    # ensuring unnecessary prefix #$ #PBS is removed
    qsub_opt = [ str.replace(hpc['prefix'],'') for str in qsub_opt ]

    return job_path, qsub_opt


def run_jobs(job_paths, run_locally, num_procs, depends_on_step):
    job_id_list = {}
    if run_locally:

        if not hpc['dry_run']:
            # ret_code = local_exec(job_path, log_fname, print_to_screen=True)
            # if ret_code < 0 or ret_code > 0:
            #     raise OSError('Error executing the {} job locally.'.format(step_id))

            # list of all script paths
            jobs_list = job_paths.values()

            pool = Pool()
            for idx_beg in range(0, len(jobs_list), num_procs):
                # running only on a specified num cores, one subject at a time
                pool.map(local_exec, jobs_list[idx_beg:idx_beg + num_procs])
                pool.close()
                pool.join()
        else:
            for prefix, job_path in job_paths.items():
                job_id_list[prefix] = make_dry_run(job_path)
    else:
        # num_sub_per_job = 1 # ability to specify >1 subjects per job
        # job_id_list = map(submit_queue, jobs_list)
        for prefix, (job_path, job_details) in job_paths.items():
            job_details_str = ' '.join(job_details)
            job_id_list[prefix] = submit_queue(job_details_str, depends_on_step)
            print('\t {} \t job id: {}'.format(prefix, job_id_list[prefix]))

    return True, job_id_list


def make_dry_run(cmd_str):
    # """ Simple function to perform a dry run (indicated by a -ve job id)."""

    print(cmd_str)
    job_id = random.randrange(1000)

    return -job_id


def submit_queue(job, depends_on_step):
    global hpc

    qsub_path = 'qsub'  # find_executable('qsub')

    # encoding dependencies
    if depends_on_step is not None and depends_on_step in hpc['job_ids_grouped']:
        # second condition can be False when rerunning from a existing processing
        # some steps may have been completed already - which means they are not in queue now,
        # in which this current step doesnt need to wait
        job_id_list_str = map(str, hpc['job_ids_grouped'][depends_on_step].values())
    else:
        job_id_list_str = ''

    if hpc['type'].upper() in ('ROTMAN-SGE', 'ROTMAN', 'SGE'):
        # qsub_cmd  = qsub_path + ' -terse '
        qsub_cmd = qsub_path
        terse = '-terse'
        hold_spec = '-hold_jid ' + ",".join(job_id_list_str)
    elif hpc['type'].upper() in ('HPCVL', 'QUEENSU'):
        # qsub_cmd  = qsub_path + ' -terse '
        qsub_cmd = qsub_path
        terse = '-terse'
        hold_spec = '-hold_jid ' + ",".join(job_id_list_str)
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
    full_cmd = re.sub( '\s+', ' ', full_cmd).strip()

    # splitting on space is required, otherwise subprocess throws an error,
    # as the '-opt value' in a single string gets interepreted as the name of an option
    arg_list = full_cmd.split(' ')

    if not hpc['dry_run']:
        job_id = subprocess.check_output(arg_list)
        job_id = job_id.strip()
    else:
        job_id = make_dry_run(full_cmd)

    return job_id


def submit_jobs():
    """
    Gateway to PRONTO preprocessing and optimization tool.
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
        is_done, rem_input_file = check_proc_status.run(
            [input_file, options.pipeline_file, '--skip_validation', '--not_verbose'])
        for step in cfg_pronto.STEPS_PROCESSING_STATUS:
            if step is not 'all_done':
                if not getattr(is_done, step):
                    bool_str = 'not done.'
                else:
                    bool_str = '    done.'
                print('{:>15} is {}'.format(step, bool_str))
        print ' '

        if is_done.preprocessing and is_done.optimization and is_done.QC1 and is_done.QC2:
            print "Both preprocessing and optimization seems to be finished already, along with QC."
            print " if you'd like to force preprocessing, please rename/remove/move the existing outputs and rerun."
            return
    else:
        print('This is just a dry run - generating jobs for all steps regardless of their processing status.')
        is_done = cfg_pronto.initialize_proc_status()
        rem_input_file = None

    if options.run_locally:
        hpc['type'] = "LOCAL"
    else:
        hpc['type'] = find_hpc_type(options.hpc_type)

    if options.part is 0:
        run_part_one = True
        run_part_two = True
    elif options.part is 1:
        run_part_one = True
        run_part_two = False
    elif options.part is 2:
        run_part_one = False
        run_part_two = True

    # run spatial norm only when the reference is specified and exists
    if options.reference_specified:
        run_sp_norm = True
    else:
        run_sp_norm = False

    spnorm_completed = False
    spnorm_step1_completed = False
    if options.dospnormfirst:
        sp_norm_step = 1
        print('spatial normalization (step 1) BEFORE preprocessing:')
        status_sp, job_ids_spn = process_spatial_norm(unique_subjects, options, input_file, sp_norm_step, cur_garage)
        if options.run_locally is True and (status_sp is False or status_sp is None):
            raise Exception('Spatial normalization - step 1 failed.')
        else:
            spnorm_step1_completed = True

    # submitting jobs for preprocessing for all combinations of pipelines
    if run_part_one and is_done.preprocessing is False:
        # running part 1 only on subjects with incomplete processing
        print('Preprocessing:')
        status_p1, job_ids_pOne = run_preprocessing(unique_subjects, options, rem_input_file, cur_garage)
        if options.run_locally is True and (status_p1 is False or status_p1 is None):
            raise Exception('Preprocessing step failed.')

    # generating QC1 if not done already
    if is_done.QC1 is False:
        print('QC 1 :')
        status_qc1, job_ids_qc1 = run_qc_part_one(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_qc1 is False or status_qc1 is None):
            raise Exception('QC part 1 failed.')

    # submitting jobs for optimization
    if run_part_two and options.analysis != "None" and is_done.optimization is False:
        # optimization is done for ALL the subjects in the input file,
        # even though part 1 may have been rerun just for failed/unfinished subjects
        print('stats and optimization :')
        status_p2, job_ids_pTwo = run_optimization(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_p2 is False or status_p2 is None):
            raise Exception('Optimization failed.')

    # generating QC2 if not done already
    if is_done.QC2 is False:
        print('QC 2 :')
        status_qc2, job_ids_qc2 = run_qc_part_two(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_qc1 is False or status_qc1 is None):
            raise Exception('QC part 1 failed.')

    # finishing up the spatial normalization
    if run_sp_norm:
        if spnorm_step1_completed is True:
            sp_norm_step = 2
        else:
            sp_norm_step = 0

        print('spatial normalization (step 2): ')
        status_sp = process_spatial_norm(unique_subjects, options, input_file, sp_norm_step, cur_garage)
        if options.run_locally is True and (status_sp is False or status_sp is None):
            raise Exception('Spatial normalization - steps 2 and later failed.')

        print('group mask generation: Submitting jobs ..')
        status_gm = process_group_mask_generation(unique_subjects, options, input_file, cur_garage)
        if options.run_locally is True and (status_gm is False or status_gm is None):
            raise Exception('Group mask generation failed.')

    # saving the job ids and hpc cfg to facilitate a status update in future
    save_hpc_cfg_and_jod_ids(cur_garage)


def save_hpc_cfg_and_jod_ids(cur_garage):
    "Saves the job ids and hpc cfg to facilitate a status update in future"

    # saving the job ids
    job_id_file = os.path.join(cur_garage, file_name_job_ids_by_group)
    if os.path.isfile(job_id_file):
        os.remove(job_id_file)
    with open(job_id_file, 'wb') as jlist:
        pickle.dump(hpc['job_ids_grouped'], jlist)

    # saving the config
    cfg_file = os.path.join(cur_garage, file_name_hpc_config)
    if os.path.isfile(cfg_file):
        os.remove(cfg_file)
    with open(cfg_file, 'wb') as hcf:
        pickle.dump(hpc, hcf)


if __name__ == '__main__':
    submit_jobs()
