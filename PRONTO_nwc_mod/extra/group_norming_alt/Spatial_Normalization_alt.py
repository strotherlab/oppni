#!/usr/bin/python

#$ -cwd -j y
#$ -S /usr/bin/python

# SPATIAL_NORMALIZATION: This script is used to submit Spatial Normalization 
# jobs, as part of the preprocessing pipeline. Consists of 2 parts.
#
# PART-1: generates transformations matrices, and spatially normalized data,
#         for optimal pipelines. Syntax:
#         ./Spatial_Normalization.py -p 1 -i {input_file.txt} -r {reference_volume.nii} -m {optimization metric}
#
# PART-2: obtains group brain mask and tissue maps. Syntax:
#         ./Spatial_Normalization.py -p 2 -i {input_file.txt} -o {output_prefix}
#

import os
from optparse import OptionParser
from datetime import datetime
from time import mktime

def time_stamp():
    t = datetime.now()
    t = mktime(t.timetuple()) + 1e-6 * t.microsecond
    return "{0:.7f}".format(t)

# task_id = int(os.environ['SGE_TASK_ID'])
# job_id = os.environ['JOB_ID']

# parse different input flags
parser = OptionParser()
parser.add_option("-p",action="store", dest="part", type="int")
parser.add_option("-i", action="store", dest="inputdata")
parser.add_option("-r",action="store", dest="reference")
parser.add_option("-c", action="store", dest="pipeline")
parser.add_option("-o", action="store", dest="out_prefix")

(options, args) = parser.parse_args()


# for pipeline steps 1,2
if options.part == 1:

    print 'Spatial normalization onto a pre-defined template'

    # submit into queue
    cmd = "./scripts_other/Prep_Transform.out {0} {1}".format(options.inputdata, options.reference)

    #execute command
    os.system(cmd)

    # create temporary directory for input files, if not yet done
    temp_dir = "submission_temp"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    # read each line of input file, turn into "temp" file, and submit to queue
    with open(options.inputdata) as f:

        lines = f.readlines()

        for line in lines:

            temp_file = "{1}/3_{0}".format(time_stamp(), temp_dir)
            with open(temp_file, 'w') as p:
                p.write(line)

            # submit into queue
            cmd = "qsub -q bigmem.q -cwd -j y -b y ./scripts_other/Run_Transform_alt.out {0} {1} {2}".format(temp_file, options.reference, options.pipeline)

            #execute command
            os.system(cmd)

# for pipeline step 3 (post-processing optimization)
elif options.part == 2:

    # declare pipeline step-3;  NB replaced "-q all.q" with "-q bigmem.q" to call all nodes
    cmd = "qsub -q bigmem.q -cwd -j y -b y ./scripts_other/mask_tissue_wrapper.py -i {0} -o {1}".format(options.inputdata, options.out_prefix)

    #execute command
    os.system(cmd)


