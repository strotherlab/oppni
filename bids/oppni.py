#!/usr/bin/env python

import argparse
import os, sys
from glob import glob

import subprocess

# import grabbids
# import cfg_front as cfg_oppni
# import proc_status_front as check_proc_status

global mcr_path, afni_path, fsl_path

def validate_env_var(var):
    assert os.getenv(var) is not None, "Path {} is not defined. Fix your environment and rerun.".format(var)


def validate_user_env():
    """Ensuring the user has proper environment."""
    global mcr_path, afni_path, fsl_path

    for var in ['AFNI_PATH', 'FSL_PATH', 'MCR_PATH']:
        validate_env_var(var)

    fsl_path  = os.getenv('FSLDIR')
    mcr_path  = os.getenv('MCR_PATH')
    afni_path = os.getenv('AFNI_PATH')


def parse_args_check():
    "Parsing the cmd line args and making basic checks."

    parser = argparse.ArgumentParser(description='OPPNI tool')
    parser.add_argument('bids_dir', help='The directory with the input dataset '
                                         'formatted according to the BIDS standard.')
    parser.add_argument('output_dir', help='The directory where the output files '
                                           'should be stored. If you are running group level analysis '
                                           'this folder should be prepopulated with the results of the'
                                           'participant level analysis.')
    parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                                               'Multiple participant level analyses can be run independently '
                                               '(in parallel) using the same output_dir.',
                        choices=['participant', 'group','participant1', 'group1', 'participant2', 'group2'])

    parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                                                    'corresponds to sub-<participant_label> from the BIDS spec '
                                                    '(so it does not include "sub-"). If this parameter is not '
                                                    'provided all subjects should be analyzed. Multiple '
                                                    'participants can be specified with a space separated list.')

    args = parser.parse_args()

    return args


def run_part_one(bids_dir, subject_label, task_name, output_dir):
    """Runs unit-level computation on a given list of subjects/runs."""

    anat_file = 'None'
    physio_file = 'None'
    drop_beg = 0
    drop_end = 0

    # find all bold
    epi_dir = os.path.join(bids_dir, "sub-%s" % subject_label, "func")
    func_pattern = "sub-{}_*task-{}*_bold.nii*".format(subject_label, task_name)
    epi_list = glob(os.path.join(epi_dir, func_pattern))

    physio_pattern = "sub-{}_*task-{}*_physio*".format(subject_label, task_name)

    anat_dir = os.path.join(bids_dir, "sub-%s" % subject_label, "anat")
    anat_pattern = "sub-{}_T1w.nii*".format(subject_label)
    anat_file_list = glob(os.path.join(anat_dir, anat_pattern))

    task_json = os.path.join(bids_dir,'task-{}_bold.json'.format(task_name))
    events_tsv = os.path.join(epi_dir, "sub-{}_task-{}_events.tsv".format(subject_label, task_name))

    # cmd = "run_oppni.sh %s PART1 %s %s %s %s %s %d %d %s %s" % (mcr_path, subject_label, func_file[0],
    #                                                          output_dir, anat_file_list[0],
    #                                                          physio_file, drop_beg, drop_end, task_json,
    #                                                          events_tsv)

    output_dir_sub = os.path.join(output_dir,"sub-%s" % subject_label)
    cmd = "run_oppni.sh %s PART1 %s %s %s %s %d %d %s %s" % (mcr_path, epi_list[0],
                                                                output_dir_sub, anat_file_list[0],
                                                                physio_file, drop_beg, drop_end, task_json,
                                                                events_tsv)

    print(cmd)
    try:
        txt_out = subprocess.check_output(cmd, shell=True)
        print txt_out
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise


def run_oppni():
    """Main function handling different levels of the OPPNI processing and analysis."""

    global mcr_path, afni_path, fsl_path

    args = parse_args_check()

    validate_user_env()

    task_group = 'linebisection'

    # status_check_dir = os.path.join(args.output_dir, 'status_checks')
    # if not os.path.exists(status_check_dir):
    #     os.mkdir(status_check_dir)

    # ref_atlas = os.path.join(fsl_path,'data','standard','MNI152_T1_1mm.nii.gz')

    # running participant level
    if args.analysis_level in [ "participant", "participant1"]:

        assert args.participant_label not in [None, ''], "Subject label must be non-empty."
        run_part_one(args.bids_dir, args.participant_label, task_group, args.output_dir )

        # # TODO generate an error when the processing is not successful
        # input_file = os.path.join(args.output_dir, 'input_files',args.participant_label+'.txt')
        # pipln_comb = os.path.join(args.output_dir, 'pipeline_combinations.txt')
        # if not ( os.path.exists(input_file) and os.path.exists(pipln_comb) ):
        #     raise IOError('Input file and/or pipeline combinations do not exist.')
        # proc_status = check_proc_status.run(input_file, pipln_comb)
        # if not proc_status.preprocessing:
        #     raise IOError('OPPNI preprocessing for {} is not complete.'.format(args.participant_label))


    # running group level
    elif args.analysis_level in ["group", "group1"]:

        # for all subjects
        raise NotImplemented


if __name__ == '__main__':
    run_oppni()