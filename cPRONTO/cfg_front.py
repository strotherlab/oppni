
from collections import namedtuple

# the set of 5 steps must correspond to the AFNI processing done in matlab
CODES_PREPROCESSING_STEPS = ["MOTCOR", "CENSOR", "RETROICOR", "TIMECOR",
             "SMOOTH", "DETREND", "MOTREG", "TASK",
             "GSPC1", "PHYPLUS", "CUSTOMREG", "LOWPASS"]

CODES_INPUT_SECTIONS = ['IN', 'OUT', 'DROP', 'STRUCT', 'PHYSIO', 'TASK', 'CUSTOMREG']

TASK_MANDATORY_FIELDS = ['UNIT', 'TR_MSEC', 'TYPE']

CODES_PRONTO_STEPS = [ 'PART1', 'QC1', 'PART2', 'QC2', 'GMASK', 'SPNORM' ]

CODES_ANALYSIS_MODELS = ['None', 'LDA', 'GNB', 'GLM', 'erCVA', 'erGNB', 'erGLM', 'SCONN']

CODES_METRIC_LIST = ["dPR", "P", "R"]

CODES_OPTIM_SCHEMES = [ 'CON', 'FIX', 'IND' ]
CODES_OPTIM_SCHEMES_OPTIONS = [ 'CON', 'FIX', 'IND', 'ALL' ] # ['CON', 'FIX', 'IND']

CODES_SUBSPACES = [ 'onecomp','multicomp' ]

WARP_TYPES = [ "affine", "nonlinear" ]

CODES_SPM_TYPES = ['zscore', 'corr']

SLICE_TIMING_PATTERNS = ['alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z', 'auto_hdr']

STEPS_PROCESSING_STATUS_WITH_QC = ['preprocessing', 'stats', 'optimization', 'spnorm', 'gmask', 'QC1', 'QC2', 'all_done' ]
STEPS_PROCESSING_STATUS         = ['preprocessing', 'stats', 'optimization', 'spnorm', 'gmask', 'QC1', 'QC2', 'all_done' ]

qc_out_format = 'pdf' # 'png
outputs_qc_part_one = [ 'FIG1_motion_statistics.' + qc_out_format,
                        'FIG2_spm_artifact.' + qc_out_format,
                        'FIG3_pipeline_similarity_by_dataset.' + qc_out_format,
                        'FIG4_effects_of_pipeline_steps.' + qc_out_format ]

outputs_qc_part_two = [ 'FIG1_optimized_pipeline_steps.' + qc_out_format,
                        'FIG2_optimized_performance_metrics.' + qc_out_format,
                        'FIG3_spatial_norm_statistics.' + qc_out_format,
                        'FIG4_inter_subject_SPM_similarity.' + qc_out_format ]

# version check

AFNI_VERSION_TESTED    = 'Version AFNI_2011_12_21_1014\n[[Precompiled binary linux_openmp_64: Aug 31 2012]]\n'
FLIRT_VERSION_TESTED   = 'FLIRT version 6.0\n'
MELODIC_VERSION_TESTED = '\nMELODIC Version 3.14\n\n'

def initialize_proc_status():
    """Method to initialize a namedtuple to store the indicators of processing status"""

    flag_list = ' '.join(STEPS_PROCESSING_STATUS_WITH_QC) + ' rem_input_file rem_spnorm_file'
    proc_status = namedtuple("status_indicator", flag_list)
    # initializing them to False
    for flag in proc_status._fields:
        setattr(proc_status, flag, False)

    return proc_status