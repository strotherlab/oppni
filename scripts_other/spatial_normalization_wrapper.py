#!/usr/bin/python
#$ -S /usr/bin/python

# OPTIMIZATION_WRAPPER: this script submits pipeline optimization jobs onto
# queue (part-2 of Run_Pipelines.py). This is done after generating all
# preprocessed outputs. Syntax:
#
# ./scripts_other/optimization_wrapper.py -i {input_files.txt} -m {optimization metric} -o {output prefix}
#

import os, sys
from optparse import OptionParser
import distutils
import distutils.spawn

def load_settings():

    codepathname  = os.path.dirname(sys.argv[0])
    codefull_path = os.path.abspath(codepathname)
    setting_name = codefull_path+"/../scripts_matlab/SETTINGS.txt"

    with open(setting_name) as f:
        lines = f.readlines()
        for line in lines:
            temp=line.split("=")
            if len(temp)==2:
                var1=temp[0].rstrip()
                val1=temp[1].rstrip()
                var1=var1.rstrip("\n")
                val1=val1.rstrip("\n")
                if len(val1)>0:
                    oldvar1 = ""
                    try:
                        oldvar1 = os.environ[var1]
                        os.environ[var1] = oldvar1+":"+val1
                    except KeyError:
                        os.environ[var1] = val1



CODE_VERSION = "$Revision: 158 $"
CODE_DATE    = "$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $"

os.system("echo %s >current_job_nosge.txt" % os.getpid())

parser = OptionParser()
parser.add_option("-i", action="store", dest="inputdata")
parser.add_option("-r", action="store", dest="reference")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("-v", action="store", dest="voxelsize")
parser.add_option("-s", action="store", dest="step")
parser.add_option("-e", action="store", dest="environment")
parser.add_option("--DEOBLIQUE", action="store", dest="DEOBLIQUE")

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

if hasattr(options,'DEOBLIQUE'):
    DEOBLIQUE = options.DEOBLIQUE
else:
    DEOBLIQUE = "0"

environment = environment.replace('\\n','\x0A')
print options.voxelsize

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["OMP_NUM_THREADS"] = numcores
os.environ["FSLOUTPUTTYPE"] = "NIFTI"
load_settings()

if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"spatial_normalization('{1}','{2}','{3}','{4}','{5}');\"".format(environment,options.inputdata, options.reference, options.voxelsize, step,DEOBLIQUE)
if (environment.find("matlab")>=0):
    cmd = "time -p {0}  -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');try, spatial_normalization('{2}','{3}','{4}','{5}','{6}');catch,display(lasterr); sge_exit(100);end, exit\"".format(environment,numcores,options.inputdata, options.reference, options.voxelsize, step,DEOBLIQUE)

os.system(cmd)


