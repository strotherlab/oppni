#!/usr/bin/python

#$ -cwd -j y
#$ -S /usr/bin/python
import os.path
import sys
codepathname  = os.path.dirname(sys.argv[0])
codefull_path = os.path.abspath(codepathname)
dirname   = codefull_path+"/scripts_gui/"
filename  = codefull_path+"/scripts_gui/pronto_gui"
if os.path.isfile(filename):
    os.system(filename)
else:
    arch  = codefull_path+"/scripts_gui/pronto_gui.tar.gz"
    os.system("tar -zxf %s -C %s" % (arch,dirname))
    os.system(filename)

