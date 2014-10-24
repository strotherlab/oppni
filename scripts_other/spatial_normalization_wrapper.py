#!/usr/bin/python
#$ -S /usr/bin/python

# OPTIMIZATION_WRAPPER: this script submits pipeline optimization jobs onto
# queue (part-2 of Run_Pipelines.py). This is done after generating all
# preprocessed outputs. Syntax:
#
# ./scripts_other/optimization_wrapper.py -i {input_files.txt} -m {optimization metric} -o {output prefix}
#

import os
from optparse import OptionParser
import distutils
import distutils.spawn

parser = OptionParser()
parser.add_option("-i", action="store", dest="inputdata")
parser.add_option("-r", action="store", dest="reference")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("-v", action="store", dest="voxelsize")
parser.add_option("-s", action="store", dest="step")
parser.add_option("-e", action="store", dest="environment")

(options, args) = parser.parse_args()

if hasattr(options,'numcores'):
    numcores = options.numcores
else:
    numcores = '1'
if numcores==None:
    numcores = '1'

if hasattr(options,'step'):
    step = options.step
else:
    step = '0'
if numcores==None:
    step = '0'
if hasattr(options,'environment'):
    environment = options.environment
else:
    environment = 'octave'
if environment==None:
    environment = 'octave'

environment = environment.replace('\\n','\x0A')
print options.voxelsize

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["OMP_NUM_THREADS"] = numcores

if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"spatial_normalization('{1}','{2}','{3}','{4}');\"".format(environment,options.inputdata, options.reference, options.voxelsize, step)
if (environment.find("matlab")>=0):
    cmd = "time -p {0}  -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');spatial_normalization('{2}','{3}','{4}','{5}');exit\"".format(environment,numcores,options.inputdata, options.reference, options.voxelsize, step)

os.system(cmd)
