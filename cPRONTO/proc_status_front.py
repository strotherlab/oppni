#!/usr/bin/env python

import fnmatch
import glob
import os
import re
import sys
from argparse import ArgumentParser
from cStringIO import StringIO

import oppni
import cfg_front as cfg_pronto

tick_mark = u'\u2713'.encode('utf-8')
crossed = u'\u2718'.encode('utf-8')

def all_files_exist_in(files):
    for fp in files:
        if not os.path.isfile(fp):
            return False
    return True


def count_variations_pipeline_step(pipelineFile, descr):
    """
    Count the number of configurations specified for a given pipeline step.

    :param pipelineFile: pipeline config file.
    :param descr: string specifying the exact step, such as 'MOTCOR', 'CENSOR', 'RETROICOR', 'TIMECOR', 'SMOOTH'
    :return:
    :rtype: int
    """
    reStep = re.compile(descr + "=\[([\d,]+)\]", re.DOTALL)
    digits = re.compile('(\d)+')
    with open(pipelineFile, 'r') as pipID:
        pipText = pipID.read()
        stepLine = reStep.search(pipText).group(1)
        return len(digits.findall(stepLine))


def is_done_spnorm(sub_prefix, out_dir):
    """
    Checks for the completeness of processing, according to the documentation provided in Pipeline_Part1.m

    :param sub_prefix:
    :param out_dir:
    :return:
    :rtype: bool, str
    """

    # Notes from Nathan:
    #  I would consider two important "checkpoints" in the registration process.
    # The first is generation of "intermediate_processed/spat_norm/Transmat_EPItoREF_[ prefix ].mat", which means you have a complete affine transform from functional to template space.
    #
    # The second is the generation of the transformed, optimized pipeline outputs
    #
    # optimization_results/spms/rSPM_[ prefix ]_CON_FIX_IND_sNorm.nii
    # optimization_results/processed/Proc_[ prefix ]_CON_sNorm.nii
    # optimization_results/processed/Proc_[ prefix ]_FIX_sNorm.nii
    # optimization_results/processed/Proc_[ prefix ]_IND_sNorm.nii
    #
    # These are the essential final outputs needed to go forward with group masking, QC, etc.

    intProcFolder = os.path.join(out_dir, 'intermediate_processed')
    optim_dir = os.path.join(out_dir, 'optimization_results')

    must_exist_list = list()
    must_exist_list.append(os.path.join(intProcFolder, 'spat_norm', 'Transmat_EPItoREF_'+ sub_prefix + '.mat'))
    must_exist_list.append(os.path.join(optim_dir, 'spms', 'rSPM_' + sub_prefix + '_CON_FIX_IND_sNorm.nii'))
    for scheme in cfg_pronto.CODES_OPTIM_SCHEMES:
        must_exist_list.append(os.path.join(optim_dir, 'processed', 'Proc_' + sub_prefix + '_' + scheme + '_sNorm.nii'))

    exist_bool = map(os.path.exists, must_exist_list)
    if not all(exist_bool):
        return False, " spat. norm : Incomplete " + crossed
    else:
        return True , " spat. norm : Done.      " + tick_mark


def is_done_part1_afni(subPrefix, outFolder, numPipelineSteps):
    """
    Checks for the completeness of processing, according to the documentation provided in Pipeline_Part1.m

    :param subPrefix:
    :param outFolder:
    :param numPipelineSteps:
    :return:
    :rtype: bool
    """
    intProcFolder = os.path.join(outFolder, 'intermediate_processed')
    afniFolder = os.path.join(intProcFolder, 'afni_processed')

    file1 = os.path.join(intProcFolder, 'diagnostic', subPrefix + '_mc+smo_QC_output.mat')
    file2 = os.path.join(intProcFolder, 'mpe', subPrefix + '_mpe')
    file3 = os.path.join(intProcFolder, 'mpe', subPrefix + '_maxdisp')
    file4 = os.path.join(intProcFolder, 'afni_processed', subPrefix + '_baseproc.nii')
    mustExistList = (file1, file2, file3, file4)

    # everything starts with MOTCOR, like session1_ID2382_run3_m0c0p0t0s6.nii
    fileList = glob.glob(os.path.join(afniFolder, subPrefix + '_m*.nii'))
    fileList = set(fileList) - set(fnmatch.filter(fileList, '*baseproc*'))
    if len(fileList) != numPipelineSteps:
        # print "unequal AFNI steps - expected : ", numPipelineSteps, " actual: ", len(fileList)
        PipelineStepsComplete = False
    else:
        PipelineStepsComplete = True

    if (not PipelineStepsComplete) or (not all_files_exist_in(mustExistList)):
        return False, " preproc : Incomplete " + crossed
    else:
        return True , " preproc : Done.      " + tick_mark


def is_done_part1_stats(subPrefix, outFolder):
    """
    Checks for the completeness of processing in the fixed part of Part 1,
        according to the documentation provided in Pipeline_Part1.m

    :param subPrefix:
    :param outFolder:
    :return:
    :rtype: bool
    """
    intMetricFolder = os.path.join(outFolder, 'intermediate_metrics')
    param_file1 = os.path.join(intMetricFolder, 'res0_params', 'params_' + subPrefix + '.mat')
    metricFile1 = os.path.join(intMetricFolder, 'res1_spms'  , 'spms_'   + subPrefix + '.mat')
    metricFile2 = os.path.join(intMetricFolder, 'res2_temp'  , 'temp_'   + subPrefix + '.mat')
    metricFile3 = os.path.join(intMetricFolder, 'res3_stats' , 'stats_'  + subPrefix + '.mat')

    mustExistList = (param_file1, metricFile1, metricFile2, metricFile3)
    if not all_files_exist_in(mustExistList):
        return False, " Metrics : Incomplete " + crossed
    else:
        return True , " Metrics : Done.      " + tick_mark


def is_done_group_mask_gen(out_dir):
    """
    Checks for the successful generation of group mask

    :param out_dir:
    :return:
    :rtype: bool
    """
    gmask_dir = os.path.join(out_dir, 'GroupMasks')
    file_names = ['group_consensus_mask.nii', 'group_spat_norm_qc.mat', 'group_mean_NN_WM.nii']

    must_exist_list = list()
    for file in file_names:
        must_exist_list.append(os.path.join(gmask_dir,file))

    if not all_files_exist_in(must_exist_list):
        print('GMASK: failed.')
        return False
    else:
        print('GMASK:   done.')
        return True


def is_done_QC_part_one(out_dir):
    """
    Checks for the completeness of Quality Control (Part 1)

    :param out_dir:
    :return:
    :rtype: bool
    """

    qc1_dir = os.path.join(out_dir, 'QC1_results')
    qc1_res = os.path.join(qc1_dir, 'output_qc1.mat')

    mustExistList = []
    for png in cfg_pronto.outputs_qc_part_one:
        mustExistList.append(os.path.join(qc1_dir, png))

    if not all_files_exist_in(mustExistList):
        print('QC1: Some/all of the visualizations have NOT been exported as images.')
    else:
        print('QC1: figures done.')

    if os.path.isfile(qc1_res) and os.path.getsize(qc1_res) > 0:
        #print('QC1: Necessary diagnostic data to produce QC1 visualizations exist.')
        return True
    else:
        print('QC1: diagnostic data (path below) doesn\'t exist or is empty.\n\t {}'.format(qc1_res))
        return False


def is_done_QC_part_two(out_dir):
    """
    Checks for the completeness of Quality Control (Part 2)

    :param out_dir:
    :return:
    :rtype: bool
    """

    qc2_dir = os.path.join(out_dir, 'QC2_results')
    qc2_res = os.path.join(qc2_dir, 'output_qc2.mat')

    mustExistList = []
    for png in cfg_pronto.outputs_qc_part_two:
        mustExistList.append(os.path.join(qc2_dir, png))

    if not all_files_exist_in(mustExistList):
        print('QC2: Some/all of the visualizations have NOT been exported as images.')
    else:
        print('QC2: figures done.')

    if os.path.isfile(qc2_res) and os.path.getsize(qc2_res) > 0:
        #print('QC2: Necessary diagnostic data to produce visualizations seems to exist.')
        return True
    else:
        print('QC2: diagnostic data (path below) doesn\'t exist or is empty.\n\t {}'.format(qc2_res))
        return False


def is_done_part_two_opt_summary(out_dir):
    """
    Checks for the completeness of processing in optimization (Part 2),
        according to the documentation provided in Pipeline_Part2.m.
        This is rather a loose check - only looking for final optimization_summary.mat.

    :param out_dir:
    :rtype: bool
    """

    optim_dir = os.path.join(out_dir, 'optimization_results')
    file1 = os.path.join(optim_dir, 'matfiles', 'optimization_summary.mat')
    mustExistList = (file1,)

    if not all_files_exist_in(mustExistList):
        #print "P2 : Optim. summary : Incomplete "
        return False
    else:
        #print "P2 : Optim. summary : Finished.  "
        return True


def flush_out_text(old_stdout, my_stdout, not_verbose):

    sys.stdout = old_stdout
    if not not_verbose:
        output_text = my_stdout.getvalue()
        print(output_text)

    my_stdout.close()


def parse_args_check(input_args):
    """
    Arg parser!
    """

    parser = ArgumentParser(prog='proc_status_oppni.py')
    parser.add_argument("InputFile", help="Input file specifying subject's data and output folder.")
    parser.add_argument("PipelineFile", help="Pipeline file specifying choices for various preprocessing steps.")

    parser.add_argument("-s", "--skip_validation", dest="skip_validation",
                        action="store_true", default=False,
                        help="skips validation of input and pipeline files.")
    parser.add_argument("-q", "--not_verbose", dest="not_verbose",
                        action="store_true", default=False,
                        help="Removes the text output to screen.")

    try:
        args = parser.parse_args(input_args)
    except:
        parser.print_help()
        sys.exit(0)

    inputFile = os.path.abspath(args.InputFile)
    assert os.path.exists(inputFile), "Input file doesn't exist!"

    pipFile = os.path.abspath(args.PipelineFile)
    assert os.path.exists(pipFile), "Pipeline file doesn't exist!"

    return inputFile, pipFile, args.not_verbose


def run(input_args):
    """
    Checks for the completeness of processing in both initial preprocessing stage (Part 1),
        as well as in the optimization stage (Part 2), and builds the list of failed subjects for reprocessing.
        Usage: <code> InputFile PipelineFile

    """

    old_stdout = sys.stdout
    sys.stdout = my_stdout = StringIO()

    inputFile, pipFile, not_verbose = parse_args_check(input_args)

    proc_status = cfg_pronto.initialize_proc_status()

    # estimating the number of pipeline steps needed for
    numPipelineSteps = 1
    # "MOTCOR", "CENSOR", "RETROICOR", "TIMECOR", "SMOOTH"
    # only the combinations for the first 5 steps need to be multiplied
    for flag in cfg_pronto.CODES_PREPROCESSING_STEPS[0:5]:
        numPipelineSteps *= count_variations_pipeline_step(pipFile, flag)

    outFolderStatus = os.path.dirname(inputFile)
    name_suffix = os.path.splitext(os.path.basename(inputFile))[0]

    try:
        resubmit_part1_file = os.path.join(outFolderStatus, name_suffix + '-resubmit-part1.txt')
        resub_part1 = open(resubmit_part1_file, 'w')

        resubmit_spnorm_file = os.path.join(outFolderStatus, name_suffix + '-resubmit-spnorm.txt')
        resub_spnorm = open(resubmit_spnorm_file, 'w')
        writable = True
    except:
        print('Unable to write list of failed subjects - will skip resubmission. ')
        writable = False
        resubmit_part1_file = resubmit_spnorm_file = None
        resub_part1 = resub_spnorm = None


    num_subjects = 0
    failed_count_preproc = 0
    failed_count_stats = 0
    failed_count_spnorm = 0

    # # regex to extract the relevant parts of the input file
    reOut = re.compile(r"OUT=([\w\./+_-]+)[\s]*")

    try:

        # TODO reimplement this using front.validate_input_file to parse and iterative over subjects in input file
        common_out_dir = ' '
        with open(inputFile, 'r') as fID:
            for inputLine in fID.readlines():
                num_subjects += 1

                out_results = reOut.search(inputLine)
                if out_results is not None:
                    outFolderSpec = out_results.group(1)
                else:
                    print 'Either OUT= not specifed or contains an invalid path that can not be parsed.'
                    print 'Only Alphanumeric, underscore (_), hyphen (-) and plus (+) characters are allowed. skipping this line {}'.format(num_subjects)
                    continue

                # sometimes the output prefixes are specified with .nii extention
                # stripping it off (like PRONTO does internally)
                outFolderSpec, tmp_ext = os.path.splitext(outFolderSpec)
                subjectPrefix = os.path.basename(outFolderSpec)
                subjectPrefix, tmp_ext = os.path.splitext(subjectPrefix)

                out_dir = os.path.dirname(outFolderSpec)
                if num_subjects == 1:
                    # pronto saves the optimization results in the output folder specified for the first subject
                    common_out_dir = out_dir

                part1_preproc_done, msg1 = is_done_part1_afni(subjectPrefix, out_dir, numPipelineSteps)
                part1_stats_done  , msg2 = is_done_part1_stats(subjectPrefix, out_dir)

                spnorm_done, msg3 = is_done_spnorm(subjectPrefix, out_dir)

                print('{:>15}:  {} \t {} \t {}'.format(subjectPrefix, msg1, msg2, msg3))

                if not part1_preproc_done or not part1_stats_done:
                    failed_count_preproc += 1
                    failed_count_stats += 1
                    if writable:
                        resub_part1.write(inputLine)

                if not spnorm_done:
                    failed_count_spnorm += 1
                    if writable:
                        resub_spnorm.write(inputLine)

        # print out summary
        print "\nSummary: \n# subjects: ", num_subjects
        # part 1
        if failed_count_preproc == 0:
            proc_status.preprocessing = True
            proc_status.stats = True
            print "P1 : all finished."
            if resubmit_part1_file is not None:
                os.remove(resubmit_part1_file)
            if not writable:
                resubmit_part1_file = None
            else:
                resubmit_part1_file = 'writable'
        else:
            proc_status.preprocessing = False
            proc_status.stats = False
            proc_status.rem_input_file = resubmit_part1_file
            print "P1 : incomplete \t  # subjects/runs failed: {} / {} ({:.0f}%)".format(
                failed_count_preproc, num_subjects, (100 * failed_count_preproc / num_subjects))
            print "\t resubmit list : ", resubmit_part1_file

        # part 2
        proc_status.optimization = (failed_count_stats == 0) and is_done_part_two_opt_summary(common_out_dir)
        if not proc_status.optimization:
            print "P2 : Incomplete \n\t # sbujects whose stats need to be computed: {}".format(failed_count_stats)
        else:
            print "P2 : all finished.  "

        # spnorm
        if failed_count_spnorm == 0:
            proc_status.spnorm = True
            print "SPNORM : all finished."
            if resubmit_spnorm_file is not None:
                os.remove(resubmit_spnorm_file)
            if not writable:
                resubmit_spnorm_file = None
            else:
                resubmit_spnorm_file = 'writable'
        else:
            proc_status.spnorm = False
            proc_status.rem_spnorm_file = resubmit_spnorm_file
            print "SPNORM : incomplete \t  # subjects/runs failed: {} / {} ({:.0f}%)".format(
                failed_count_spnorm, num_subjects, (100 * failed_count_spnorm / num_subjects))
            print "\t resubmit list : ", resubmit_spnorm_file

        # GMASK
        proc_status.gmask = is_done_group_mask_gen(common_out_dir)

        # QC
        proc_status.QC1 = is_done_QC_part_one(common_out_dir)
        proc_status.QC2 = is_done_QC_part_two(common_out_dir)

        if writable:
            resub_spnorm.close()
            resub_part1.close()

    except Exception as e:
        print "following exception occurred: \n{}".format(e)
        raise

    finally:
        # making sure output of the partial execution will be shown to the user
        flush_out_text(old_stdout, my_stdout, not_verbose)

    # summary indicator of all flags
    proc_status.all_done = True
    for flag in cfg_pronto.STEPS_PROCESSING_STATUS:
        if getattr(proc_status,flag) is False and flag is not 'all_done':
            proc_status.all_done = False
            break

    return proc_status, resubmit_part1_file, resubmit_spnorm_file


if __name__ == '__main__':
    # passing the command line positional args as such, without the script name
    run(sys.argv[1:])
