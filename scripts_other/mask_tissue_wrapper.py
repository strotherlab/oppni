#!/usr/bin/python
#$ -S /usr/bin/python

# MASK_TISSUE_WRAPPER: This script is used to submit spatial normalization,
# part-2 job into queue (from Spatial_Normalization.py). Syntax:
#
#  ./scripts_other/mask_tissue_wrapper.py -i {input_file.txt} -o {output prefix}
# 

import os,sys
from optparse import OptionParser
import distutils
import distutils.spawn

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


CODE_VERSION = "$Revision: 171 $"
CODE_DATE    = "$Date: 2014-12-17 10:53:34 -0500 (Wed, 17 Dec 2014) $"

os.system("echo %s >current_job_nosge.txt" % os.getpid())

parser = OptionParser()
parser.add_option("-i", action="store", dest="inputdata")
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
os.environ["FSLOUTPUTTYPE"] = "NIFTI"

load_settings()

if (environment.find("octave")>=0):
    cmd = "time -p {0} -q -p scripts_matlab --eval \"group_mask_tissue_maps('{1}',[]);\"".format(environment,options.inputdata)
if (environment.find("matlab")>=0):
    cmd = "time -p {0} -nodesktop -nojvm -nosplash -r \"maxNumCompThreads({1});addpath('./scripts_matlab');try, group_mask_tissue_maps('{2}',[]);catch,display(lasterr);sge_exit(100);end, exit\"".format(environment,numcores, options.inputdata)
if (environment.find("standalone")>=0):
    mcrpath = os.environ["MCR_PATH"]
    MCRJRE  =  mcrpath+"/sys/java/jre/glnxa64/jre/lib/amd64"
    os.environ["LD_LIBRARY_PATH"] = ".:"+mcrpath+"/runtime/glnxa64:"+mcrpath+"/bin/glnxa64:"+mcrpath+"/sys/os/glnxa64:"+mcrpath+"/native_threads:"+mcrpath+"/server:"+mcrpath+"/client:"+mcrpath
    os.environ["XAPPLRESDIR"]=mcrpath+"/X11/app-defaults";
    cmd = "time -p ./scripts_matlab/compiled/PRONTO GMASK {0}".format(options.inputdata)
print cmd
os.system(cmd)
