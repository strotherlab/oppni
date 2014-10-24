#!/usr/bin/python
#$ -S /usr/bin/python

# PIPELINE_WRAPPER: this script submits preprocessing pipeline jobs onto
# queue (part-1 of Run_Pipelines.py). Syntax:
#
# ./scripts_other/pipeline_wrapper.py -i {input_files.txt} -c {pipeline_list.txt} -a {analysis model}
#

import os
from optparse import OptionParser

import distutils
import distutils.spawn

parser = OptionParser()
parser.add_option("-i", action="store", dest="inputdata")
parser.add_option("-c", action="store", dest="pipeline")
parser.add_option("-a", action="store", dest="analysis")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("--contrast", action="store", dest="contrast")
parser.add_option("--dospnormfirst", action="store", dest="dospnormfirst")
parser.add_option("-e", action="store", dest="environment")
parser.add_option("--modelparam", action="store", dest="modelparam")

(options, args) = parser.parse_args()
print options.environment
if hasattr(options,'numcores'):
    numcores = options.numcores
else:
    numcores = '1'
if numcores==None:
    numcores = '1'
if hasattr(options,'contrast'):
    contrast = options.contrast
if hasattr(options,'dospnormfirst'):
    dospnormfirst = options.dospnormfirst
if dospnormfirst==None:
    dospnormfirst = 0
if hasattr(options,'environment'):
    environment = options.environment
else:
    environment = 'octave'
if environment==None:
    environment = 'octave'
if hasattr(options,'modelparam'):
    modelparam = options.modelparam
else:
    modelparam = None

environment = environment.replace('\\n','\x0A')
print environment

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["OMP_NUM_THREADS"] = numcores

if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"Pipeline_PART1('{1}','{2}','{3}','{4}',0,'{5}',{6});\"".format(environment,options.inputdata, options.pipeline, options.analysis,modelparam,contrast,dospnormfirst)
if (environment.find("matlab")>=0):
    cmd = "time -p {0} -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');Pipeline_PART1('{2}','{3}','{4}','{5}',0,'{6}',{7});exit\"".format(environment,numcores, options.inputdata, options.pipeline, options.analysis,modelparam,contrast,dospnormfirst)
os.system(cmd)
