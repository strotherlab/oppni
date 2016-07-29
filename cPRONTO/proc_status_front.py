#!/usr/bin/env python

import os, sys, re, glob, fnmatch
from argparse import ArgumentParser
from collections import namedtuple
from cStringIO import StringIO

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
        print('QC1: Some/all of the visualizations have not been exported as images.')

    if os.path.isfile(qc1_res) and os.path.getsize(qc1_res) > 0:
        print('QC1: Necessary diagnostic data to produce QC1 visualizations exist.')
        return True
    else:
        print('QC1: diagnostic data doesn\'t exist or is empty.')
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
        print('QC2: Some/all of the visualizations have not been exported as images.')

    if os.path.isfile(qc2_res) and os.path.getsize(qc2_res) > 0:
        print('QC2: Necessary diagnostic data to produce visualizations seems to exist.')
        return True
    else:
        print('QC2: diagnostic data doesn\'t exist or is empty.')
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
        print "P2 : Optim. summary : Incomplete "
        return False
    else:
        print "P2 : Optim. summary : Finished.  "
        return True


def parse_args_check(input_args):
    """
    Arg parser!
    """

    parser = ArgumentParser(prog='checkProcessingStatusPronto.py')
    parser.add_argument("InputFile", help="Input file specifying subject's data and output folder.")
    parser.add_argument("PipelineFile", help="Pipeline file specifying choices for various preprocessing steps.")

    parser.add_argument("-s", "--skip_validation", dest="skip_validation",
                        action="store_true",
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
    for step in cfg_pronto.CODES_PREPROCESSING_STEPS[0:5]:
        numPipelineSteps *= count_variations_pipeline_step(pipFile, step)

    outFolderStatus = os.path.dirname(inputFile)
    resubmitListPart1File = os.path.join(outFolderStatus, 'InputFile-Resubmit-Part1.prontoStatusCheck')
    resubmitListPart1 = open(resubmitListPart1File, 'w')

    totalNumSubjects = 0
    failed_count_preproc = 0
    failed_count_stats = 0

    # # regex to extract the relevant parts of the input file
    # reIn = re.compile("IN=([\w\./_-]+)\s")
    # reTask = re.compile("TASK=([\w\./_-]+)\s")
    reOut = re.compile("OUT=([\w\./_-]+)\s")

    # TODO reimplement this using front.validate_input_file to parse and iterative over subjects in input file
    commonOutFolder = ' '
    with open(inputFile, 'r') as fID:
        for inputLine in fID.readlines():
            totalNumSubjects += 1
            # inImages  = reIn.search(inputLine).group(1)
            # taskSpec  = reTask.search(inputLine).group(1)
            outFolderSpec = reOut.search(inputLine).group(1)
            # sometimes the output prefixes are specified with .nii extention
            # stripping it off (like PRONTO does internally) 
            outFolderSpec, tmp_ext = os.path.splitext(outFolderSpec)
            subjectPrefix = os.path.basename(outFolderSpec)
            subjectPrefix, tmp_ext = os.path.splitext(subjectPrefix)

            outFolder = os.path.dirname(outFolderSpec)
            if totalNumSubjects == 1:
                # pronto saves the optimization results in the output folder specified for the first subject
                commonOutFolder = outFolder

            part1_preproc_done, msg1 = is_done_part1_afni(subjectPrefix, outFolder, numPipelineSteps)
            part1_stats_done  , msg2 = is_done_part1_stats(subjectPrefix, outFolder)
            both_parts_done = (part1_preproc_done and part1_stats_done)

            print('{}:  {} \t {}'.format(subjectPrefix, msg1, msg2))

            if not part1_preproc_done:
                failed_count_preproc += 1
                resubmitListPart1.write(inputLine)

            if not part1_stats_done:
                failed_count_stats += 1

    # print out summary
    print "\nSummary: \nTotal subjects: ", totalNumSubjects
    if failed_count_preproc == 0:
        proc_status.preprocessing = True
        print "P1 : All finished."
        os.remove(resubmitListPart1File)
        resubmitListPart1File = None

        proc_status.QC1 = is_done_QC_part_one(commonOutFolder)
        proc_status.QC2 = is_done_QC_part_two(commonOutFolder)

        proc_status.optimization = failed_count_stats==0 and is_done_part_two_opt_summary(commonOutFolder)
        if not proc_status.optimization:
            print " Resubmit PRONTO jobs with -p 2 flag, with the same input file."
    else:
        proc_status.preprocessing = False
        proc_status.rem_input_file = resubmitListPart1File
        print "P1: Incomplete \t  # Subjects/Runs Failed: {} / {} ({:.0f}%)".format(failed_count_preproc, totalNumSubjects, (100 * failed_count_preproc / totalNumSubjects))
        print "Resubmit list : ", resubmitListPart1File
        # resubmit jobs is not easily possible, due to the presence of
        # number of other arguments to PRONTO, that are not supplied to this code.
        # e.g. subprocess.call("./Run_Pipelines.py -i {0} -c {1}".format(resubmitListPart1File,pipFile), shell=True)
        #   would be insufficient!

    sys.stdout = old_stdout
    if not not_verbose:
        output_text = my_stdout.getvalue()
        print(output_text)

    my_stdout.close()

    # summary indicator of all flags
    proc_status.all_done = True
    for flag in cfg_pronto.STEPS_PROCESSING_STATUS:
        if getattr(proc_status,flag) is not False and flag is not 'all_done':
            proc_status.all_done = False
            break

    return proc_status, resubmitListPart1File


if __name__ == '__main__':
    # passing the command line positional args as such, without the script name
    run(sys.argv[1:])
