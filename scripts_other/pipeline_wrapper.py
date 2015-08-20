#!/usr/bin/python
#$ -S /usr/bin/python

# PIPELINE_WRAPPER: this script submits preprocessing pipeline jobs onto
# queue (part-1 of Run_Pipelines.py). Syntax:
#
# ./scripts_other/pipeline_wrapper.py -i {input_files.txt} -c {pipeline_list.txt} -a {analysis model}
#

import os, sys
from optparse import OptionParser

import distutils
import distutils.spawn
import subprocess


def load_settings():

    codepathname  = os.path.dirname(sys.argv[0])
    codefull_path = os.path.abspath(codepathname)
    setting_name = codefull_path+"/../config/SETTINGS.txt"

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
parser.add_option("-c", action="store", dest="pipeline")
parser.add_option("-a", action="store", dest="analysis")
parser.add_option("-n", action="store", dest="numcores")
parser.add_option("--contrast", action="store", dest="contrast")
parser.add_option("--dospnormfirst", action="store", dest="dospnormfirst")
parser.add_option("-e", action="store", dest="environment")
parser.add_option("--modelparam", action="store", dest="modelparam")
parser.add_option("--TPATTERN", action="store", dest="TPATTERN")
parser.add_option("--DEOBLIQUE", action="store", dest="DEOBLIQUE")

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

if hasattr(options,'DEOBLIQUE'):
    DEOBLIQUE = options.DEOBLIQUE
else:
    DEOBLIQUE = "0"
if hasattr(options,'TPATTERN'):
    TPATTERN = options.TPATTERN
else:
    TPATTERN = ""
if TPATTERN==None:
    TPATTERN = ""

environment = environment.replace('\\n','\x0A')

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["OMP_NUM_THREADS"] = numcores
os.environ["FSLOUTPUTTYPE"] = "NIFTI"
task_id = os.environ.get("SGE_TASK_ID")
if task_id:
	if task_id.isdigit():
		os.environ["SGE_TASK_ID"] = "%04d" % int(task_id)

load_settings()

if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"Pipeline_PART1('{1}','{2}','{3}','{4}',0,'{5}',{6},'{7}','{8}');\"".format(environment,options.inputdata, options.pipeline, options.analysis,modelparam,contrast,dospnormfirst,DEOBLIQUE,TPATTERN)
if (environment.find("matlab")>=0):
    cmd = "time -p {0} -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');try,Pipeline_PART1('{2}','{3}','{4}','{5}',0,'{6}',{7},'{8}','{9}');catch,display(lasterr),sge_exit(100);end,exit\"".format(environment,numcores, options.inputdata, options.pipeline, options.analysis,modelparam,contrast,dospnormfirst,DEOBLIQUE,TPATTERN)
if (environment.find("standalone")>=0):
    mcrpath = os.environ["MCR_PATH"]
    MCRJRE  =  mcrpath+"/sys/java/jre/glnxa64/jre/lib/amd64"
    os.environ["LD_LIBRARY_PATH"] = ".:"+mcrpath+"/runtime/glnxa64:"+mcrpath+"/bin/glnxa64:"+mcrpath+"/sys/os/glnxa64:"+mcrpath+"/native_threads:"+mcrpath+"/server:"+mcrpath+"/client:"+mcrpath
    os.environ["XAPPLRESDIR"]=mcrpath+"/X11/app-defaults";
    cmd = "time -p ./scripts_matlab/compiled/PRONTO PART1 {1} '{2}' '{3}' '{4}' 0 '{5}' '{6}' '{7}' '{8}'".format(environment,options.inputdata, options.pipeline, options.analysis,modelparam,contrast,dospnormfirst,DEOBLIQUE,TPATTERN)
    

print cmd
os.system("echo "+cmd)
os.system(cmd)


