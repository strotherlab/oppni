
from collections import namedtuple

# the set of 5 steps must correspond to the AFNI processing done in matlab
CODES_PREPROCESSING_STEPS = ["MOTCOR", "CENSOR", "RETROICOR", "TIMECOR",
             "SMOOTH", "DETREND", "MOTREG", "TASK",
             "GSPC1", "PHYPLUS", "CUSTOMREG", "LOWPASS"]

TASK_MANDATORY_FIELDS = ['UNIT', 'TR_MSEC', 'TYPE']

CODES_PRONTO_STEPS = [ 'PART1', 'QC1', 'PART2', 'QC2', 'GMASK', 'SPNORM' ]

CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'erCVA', 'erGNB', 'erGLM', 'SCONN']

CODES_METRIC_LIST = ["dPR", "P", "R"]

CODES_SUBSPACES = [ 'onecomp','multicomp' ]

CODES_SPM_TYPES = ['zscore', 'corr']

STEPS_PROCESSING_STATUS = ['preprocessing', 'optimization', 'QC1', 'QC2', 'all_done' ]

outputs_qc_part_one = [ 'FIG1_motion_statistics.png',
                        'FIG2_spm_artifact.png',
                        'FIG3_pipeline_similarity_by_dataset.png',
                        'FIG4_effects_of_pipeline_steps.png' ]

outputs_qc_part_two = [ 'FIG1_optimized_pipeline_steps.png',
                        'FIG2_optimized_performance_metrics.png',
                        'FIG3_spatial_norm_statistics.png',
                        'FIG4_inter_subject_SPM_similarity.png' ]


def initialize_proc_status():

    # creating a namedtuple to store the indicators of status
    flag_list = ' '.join(STEPS_PROCESSING_STATUS) + ' rem_input_file'
    proc_status = namedtuple("status_indicator", flag_list)
    # initializing them to False
    for flag in proc_status._fields:
        setattr(proc_status, flag, False)

    return proc_status