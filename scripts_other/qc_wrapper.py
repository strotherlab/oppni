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
import os.path

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
parser.add_option("-p", action="store", dest="part")
parser.add_option("--newmaskname", action="store", dest="newmaskname")
parser.add_option("--Npcs", action="store", dest="Npcs")
parser.add_option("-e", action="store", dest="environment")

(options, args) = parser.parse_args()
environment = options.environment
environment = environment.replace('\\n','\x0A')

load_settings()

codepathname  = os.path.dirname(sys.argv[0])
codefull_path = os.path.abspath(codepathname)
codefull_path = codefull_path + "/../extra/"
if (environment.find("octave")>=0):
    cmd = "{0} -q -p {5} --eval \"QC_wrapper('{1}','{2}','{3}','{4}');\"".format(environment,options.part,options.inputdata, options.newmaskname,options.Npcs,codefull_path)
if (environment.find("matlab")>=0):
    cmd = "{0} -nodesktop -nosplash -r \"addpath('{5}');QC_wrapper('{1}','{2}','{3}','{4}');exit\"".format(environment,options.part,options.inputdata, options.newmaskname,options.Npcs,codefull_path)
if (environment.find("standalone")>=0):
    mcrpath = os.environ["MCR_PATH"]
    MCRJRE  =  mcrpath+"/sys/java/jre/glnxa64/jre/lib/amd64"
    os.environ["LD_LIBRARY_PATH"] = ".:"+mcrpath+"/runtime/glnxa64:"+mcrpath+"/bin/glnxa64:"+mcrpath+"/sys/os/glnxa64:"+mcrpath+"/native_threads:"+mcrpath+"/server:"+mcrpath+"/client:"+mcrpath
    os.environ["XAPPLRESDIR"]=mcrpath+"/X11/app-defaults";
    cmd = "{5}/compiled/QC {1} '{2}' '{3}' '{4}'".format(environment,options.part,options.inputdata, options.newmaskname,options.Npcs,codefull_path)
    

print cmd
os.system("echo "+cmd)
os.system(cmd)

