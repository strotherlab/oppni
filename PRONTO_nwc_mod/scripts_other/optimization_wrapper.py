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
parser.add_option("-m", action="store", dest="metric")
parser.add_option("-o", action="store", dest="out_prefix")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("-k", action="store", dest="keepmean")
parser.add_option("-e", action="store", dest="environment")

(options, args) = parser.parse_args()

if hasattr(options,'numcores'):
    numcores = options.numcores
else:
    numcores = '1'
if numcores==None:
    numcores = '1'
if hasattr(options,'keepmean'):
    keepmean = options.keepmean
else: 
    keepmean = 0
if keepmean==None:
    keepmean = 0
if hasattr(options,'environment'):
    environment = options.environment
else:
    environment = 'octave'
if environment==None:
    environment = 'octave'

environment = environment.replace('\\n','\x0A')

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["OMP_NUM_THREADS"] = numcores

if (environment.find("octave")>=0):    
    cmd = "time -p {0} -q -p scripts_matlab --eval \"Pipeline_PART2('{1}','{2}',[1 0],'{3}', 1, {4});\"".format(environment,options.inputdata, options.metric, options.out_prefix, keepmean)
if (environment.find("matlab")>=0):
    cmd = "time -p {0}  -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');Pipeline_PART2('{2}','{3}',[1 0],'{4}', 1, {5});exit\"".format(environment,numcores, options.inputdata, options.metric, options.out_prefix, keepmean)
print cmd
os.system(cmd)
