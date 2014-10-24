#!/usr/bin/python
#$ -S /usr/bin/python

# MASK_TISSUE_WRAPPER: This script is used to submit spatial normalization,
# part-2 job into queue (from Spatial_Normalization.py). Syntax:
#
#  ./scripts_other/mask_tissue_wrapper.py -i {input_file.txt} -o {output prefix}
# 

import os
from optparse import OptionParser
import distutils
import distutils.spawn

parser = OptionParser()
parser.add_option("-i", action="store", dest="inputdata")
parser.add_option("-o", action="store", dest="out_prefix")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("-e", action="store", dest="environment")

(options, args) = parser.parse_args()

if hasattr(options,'numcores'):
    numcores = options.numcores
else:
    numcores = '1'
if numcores==None:
    numcores = '1'
if hasattr(options,'environment'):
    environment = options.environment
else:
    environment = 'octave'
if environment==None:
    environment = 'octave'

environment = environment.replace('\\n','\x0A')

os.environ["PIPELINE_NUMBER_OF_CORES"] = options.numcores
os.environ["OMP_NUM_THREADS"] = options.numcores


if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"group_mask_tissue_maps('{1}','{2}');\"".format(environment,options.inputdata, options.out_prefix)
if (environment.find("matlab")>=0):
    cmd = "time -p {0} -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');group_mask_tissue_maps('{2}','{3}');exit\"".format(environment,numcores, options.inputdata, options.out_prefix)

os.system(cmd)
