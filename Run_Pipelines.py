#!/usr/bin/python

#$ -cwd -j y
#$ -S /usr/bin/python

# RUN_PIPELINES: This script is used to submit fMRI pipeline optimization
# jobs. Consists of 3 parts.
#
# PART-1: runs all preprocessing pipelines & chosen analysis model. Produces
#         SPMs and pipeline performance metrics. Syntax:
#         ./Run_Pipelines.py -p 1 -i {input_file.txt} -c {pipeline_list.txt} -a {analysis model}
#
# PART-2: obtains group brain mask and tissue maps. Syntax:
#         ./Run_Pipelines.py -p 2 -i {input_file.txt} -m {optimization metric} -o {output_prefix}


import os, sys
import subprocess
from optparse import OptionParser
from datetime import datetime
from time import mktime
from time import sleep
from os import listdir
from os.path import isdir, join
from shutil import rmtree

def time_stamp():
    t = datetime.now()
    t = mktime(t.timetuple()) + 1e-6 * t.microsecond
    t = int(t*1e7)
    return "{0}".format(t)

def uniq(input):
    output = []
    list_num = []
    position = []    
    count  = 0
    for x in input:
        count = count + 1
        if x not in output:
            output.append(x)
            list_num.append(count)
            position.append(len(output)-1)
        else:
            inds = [i for (i, val) in enumerate(output) if val==x]
            position.append(inds[0])

    return list_num, position

def merge_lines(lines):
    temp  = lines[0].rstrip()
    stemp = temp.split(' ');
    physi = ''
    struct = ''
    customreg=''
    task=''
    drop='[0,0]'
    for s in stemp:
        if ("STRUCT=" in s.upper()):
            struct = s
        if ("OUT=" in s.upper()):
            outdir = s
        if ("IN=" in s.upper()):
            indir  = s
        if ("PHYSIO=" in s.upper()):
            physi = s
        if ("DROP=" in s.upper()):
            drop  = s
        if ("TASK=" in s.upper()):
            task= s
        if ("CUSTOMREG=" in s.upper()):
            customreg = s
    for line in lines[1:]:
        line = line.rstrip()
        stemp = line.split(' ')
        for s in stemp:
            if ("OUT=" in s.upper()):
                outdir = outdir+","+os.path.split(s)[1]
            if ("IN=" in s.upper()):
                indir  = indir+","+os.path.split(s)[1]
            if ("PHYSIO=" in s.upper()):
                physi = physi + ","+os.path.split(s)[1]
            if ("TASK=" in s.upper()):
                task= task+","+os.path.split(s)[1]
    mline = indir + " "+outdir+" "+physi+" "+ struct+ " "+ task + " " + drop + " " + customreg
    return mline    

def check_pipeline_file(pipeline):
    steps = ["MOTCOR","CENSOR","RETROICOR","TIMECOR","SMOOTH","DETREND","MOTREG","TASK","GSPC1","PHYPLUS","CUSTOMREG","LOWPASS"]
    with open(pipeline) as f:
        lines = f.readlines()
        for line in lines:
            temp=line.split("=")
            temp=temp[0].rstrip()
            k     = [i for i in range(0,12) if (steps[i]==temp)]  
            if not k:   
                print "Unknown preprocessing step in line: "+line.rstrip()
                print "Check file: "+pipeline
                exit(1)

def check_input_lines(line):
    line = line.rstrip()
    uline = line.split(' ')
    for temp in uline:
        t = temp.split("=")
        if ("IN=" in temp.upper()):
            t = temp.split("=")
            if not os.path.isfile(t[1]):
                print "Input file "+t[1]+" not found"
                exit(1)
        if ("TASK=" in temp.upper()):
            if not os.path.isfile(t[1]):
                print "Input file "+t[1]+" not found"
                exit(1)
def check_pid(pid):        
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True        

def wait_until_jobs_done(job_id_list_nosge,used_cores,numcores):
    if used_cores>=int(numcores):
        print used_cores
        print numcores
        while(True):
            for job_id in job_id_list_nosge:
                t = check_pid(job_id)
                if t==False:
                    used_cores = used_cores - 1
                    job_id_list_nosge.remove(job_id)
                    return used_cores
    else:
        return used_cores

def pair_list(l):
    o = []
    print l
    for i in range(0,len(l)-1,2):
        o.append(l[i]+" "+l[i+1])
    return o
    
def load_parrallel_settings():

    codepathname  = os.path.dirname(sys.argv[0])
    codefull_path = os.path.abspath(codepathname)
    setting_name  = codefull_path+"/scripts_matlab/SETTINGS.txt"

    with open(setting_name) as f:
        lines = f.readlines()
        for line in lines:
            temp=line.split("=")
            if len(temp)==2:
                var1=temp[0].rstrip()
                val1=temp[1].rstrip()
                var1=var1.rstrip("\n")
                val1=val1.rstrip("\n")
                if var1.upper()=="PARALLEL_ENVIRONMENT":
                    return val1
    return ""


CODE_VERSION = "$Revision: 173 $"
CODE_DATE    = "$Date: 2014-12-17 11:37:36 -0500 (Wed, 17 Dec 2014) $"
    
# parse different input flags
parser = OptionParser()
parser.add_option("-p","--part",action="store", dest="part", type="int",help="select pipeline optimization step, 1: Estimation step, 2: Optimization step, 3: Spatial normalization, 0: All three steps (default 0)")
parser.add_option("-i","--inputdata",action="store", dest="inputdata",help="FILE.txt contains the input and output data paths",metavar="FILE.txt")
parser.add_option("-c","--pipeline", action="store", dest="pipeline",help="select the preprocessing steps")
parser.add_option("-a","--analysis", action="store", dest="analysis",help="determine analysis model (LDA, GNB, GLM, erCVA, erGNB, erGLM, SCONN)")
parser.add_option("-m", "--metric",action="store", dest="metric",help="optimization metric")

parser.add_option("--autodetect", action="store_true", dest="autodetect",help="Automatically detect subjects and optimize for each subject independently, the lines in input files that have same structrul image (STRUCT) and output directory (OUT) are considered a subject")
parser.add_option("--dospnormfirst", action="store_true", dest="dospnormfirst",help="First normalize the data to a reference (specified by switch -r), then perform the preprocessing optimization.")
parser.add_option("--contrast",action="store", dest="contrast",help="desired task contrast, it is necessary when more than two contrasts are defined in the split info files, syntax: --Contrast \"CON1 vs CON2,CON2 vs CON3\"")
parser.add_option("-r", "--reference",action="store", dest="reference",help="anatomical reference to be used in the spatial normalization step, i.e. -p,--part=3")
parser.add_option("-n", "--numcores",action="store", dest="numcores",help="(optional) number of threads used for the process (not allowed for some SGE systems)")
parser.add_option("-q", "--queue", action="store", dest="queue",help="(optional) SGE queue name, default is bigmem_16.q")
parser.add_option("-k", "--keepmean",action="store", dest="keepmean",help="(optional) determine whether the ouput nifti files contain the mean scan (Default keepmean=0, i.e. remove the mean)")
parser.add_option("-v", "--voxelsize",action="store", dest="voxelsize",help="(optional) determine the output voxel size of nifti file")
parser.add_option("-e", "--environment",action="store", dest="environment",help="(optional) determine which software to use to run the code: matlab or octave(default)")
parser.add_option("--memory",action="store", dest="environment",help="(optional) determine minimum amount RAM needed for the job!")

parser.add_option("--convolve",action="store", dest="convolve",help="VALUE=Binary value, for whether design matrix should be convolved with a standard SPMG1 HRF.  0 = do not convolve and 1 = perform convolution",metavar="VALUE")
parser.add_option("--decision_model",action="store", dest="decision_model",help="MODEL=string specifying type of decision boundary. Either: linear for a pooled covariance model or nonlinear for class-specific covariances",metavar="MODEL")
parser.add_option("--drf",action="store", dest="drf",help="FRACTION=Scalar value of range (0,1), indicating the fraction of full-date PCA subspace to keep during PCA-LDA analysis. A drf of 0.3 is recommended as it has been found to be optimal in previous studies.",metavar="FRACTION")
parser.add_option("--Nblock",action="store", dest="Nblock",help="NUMBER= number of equal sized splits to break the data into, in order to perform time-locked averaging. Must be at least 2, with even numbers >=4, recommended to obtain robust covariance estimates",metavar="NUMBER")
parser.add_option("--WIND",action="store", dest="WIND",help="SIZE = window size to average on, in TR (usually in range 6-10 TR)",metavar="SIZE")
parser.add_option("--num_PCs",action="store", dest="num_PCs",help="NUMBER = total number of principal components to retain",metavar="NUMBER")
parser.add_option("--subspace",action="store", dest="subspace",help="COMP = string specifying either: 'onecomp'   = only optimize on CV#1 or 'multicomp' = optimize on full multidimensional subspace",metavar="SIZE")
parser.add_option("--spm",action="store", dest="spm",help="FORMAT =string specifying format of output SPM.  Options include corr (map of voxelwise seed correlations) or zscore (Z-scored map of reproducible correlation values)",metavar="FORMAT")
parser.add_option("--N_resample",action="store", dest="N_resample",help="Specify the number of resamples for multi-run analysis")
parser.add_option("--TR_MSEC",action="store", dest="TR_MSEC",help="Specify TR in msec for all entries in the input file, overides the TR_MSEC in the TASK files")
parser.add_option("--noSGE",action="store_true", dest="noSGE",help="Run the pipeline wihout using SGE, It is good to specify the number of cores using -n (or --numcores) to allow the program to run in pararallel")
parser.add_option("--numprocess",action="store", dest="numprocess",help="When running the pipeline wihout using SGE, this switch specifies the number of simultaneous processes")
parser.add_option("--DEOBLIQUE",action="store_true", dest="DEOBLIQUE",help="Correct for oblique scans (DEOBLIQUE) to improve spatial normalization")
parser.add_option("--TPATTERN",action="store", dest="TPATTERN",help="Use if data contain no slice-timing information stored in the NIFTI headers (TPATTERN)")

(options, args) = parser.parse_args()

if hasattr(options,'noSGE'):
    noSGE = options.noSGE
else:
    noSGE  = False
if noSGE==None:
    noSGE = False

sge_flag  = os.getenv("SGE_ROOT")
if sge_flag==None and not noSGE==True:
    print "Sun grid engine (SGE) has not been detected!"
    print "Use --noSGE switch if you want to run The PReprocessing OptimizatioN TOolkit (PRONTO) without SGE, this will run the code only on the current node!"
    exit(1)
if not sge_flag==None:
    print "Submitting jobs to Sun Grid Engine (SGE)"

if noSGE==True:
    print "Submitting jobs to a single node"
    print "The code will wait until the jobs are finished."
    

tt_stm = time_stamp()
codepathname  = os.path.dirname(sys.argv[0])
codefull_path = os.path.abspath(codepathname)

input_files_temp = os.getcwd()+"/Submission_temp/input_files_temp_{0}".format(tt_stm)

job_id_file = "jobID_{0}.txt".format(tt_stm)

if not os.path.exists("./Submission_temp"):
    os.mkdir("./Submission_temp")
if os.path.exists(input_files_temp):
    rmtree(input_files_temp)
os.mkdir(input_files_temp)

# Get a copy of input file
F2 = open("{0}/default.txt".format(input_files_temp),"w");
current_subject = 0
with open(options.inputdata) as f:
    lines = f.readlines()
    for line in lines:
        F2.write(line)
F2.close()

list_of_unique_subjects = range(1,current_subject+1)

# otherwise format for input file is used
statinfo = os.stat("{0}/default.txt".format(input_files_temp));
if statinfo.st_size==0:
    os.remove("{0}/default.txt".format(input_files_temp)) 
else:
    struct = []
    outdir = []
    line_stack = [];
    with open("{0}/default.txt".format(input_files_temp)) as f:
        lines = f.readlines()
        for line in lines:
            current_subject = current_subject  + 1
            check_input_lines(line)
            F2 = open("{0}/%04d.txt".format(input_files_temp) % current_subject,"w")
            F2.write(line)
            F2.close()
            # we have to know how many unique structrul MRI exists in the maps so
            uline = line.split(' ')
            for temp in uline:
                if ("STRUCT=" in temp.upper()):
                    struct.append(temp) 
                if ("OUT=" in temp.upper()):
                    outdir.append(os.path.split(temp)[0])
            
    os.remove("{0}/default.txt".format(input_files_temp)) 
    # find unique structures
    stdir = zip(struct,outdir)
    list_of_unique_subjects, position = uniq(stdir)

autodetect = options.autodetect

input_file_for_part2 = options.inputdata
# Detect runs
if autodetect == True:
    print "{0} subjects were detected".format(len(list_of_unique_subjects))
    if os.path.exists(input_files_temp):
        rmtree(input_files_temp)
    os.mkdir(input_files_temp)
    ulines = [lines[i-1] for i in list_of_unique_subjects]
    current_subject = 0
    input_file_for_part2 = "{0}/subject_autodetect.txt".format(input_files_temp)
    F3 = open(input_file_for_part2,"w")
    for i in range(0,len(ulines)):
        current_subject = current_subject  + 1
        inds = [j for (j, val) in enumerate(position) if val==i]
        F2 = open("{0}/%04d.txt".format(input_files_temp) % current_subject,"w")
        mline = [];    
        for j in inds:
            mline.append(lines[j])
        mline = merge_lines(mline)
        F2.write(mline)
        F2.close()
        F3.write(mline+"\n")
    F3.close()
    list_of_unique_subjects = range(1,current_subject+1)
    
########################## Check the Input parameters

if hasattr(options,'part'):
    part = options.part
else:
    part = 0
if part==None:
    part = 0
if hasattr(options,'queue'):
    queue_name = options.queue
else:
    queue_name = 'bigmem_16.q'
if queue_name==None:
    queue_name = 'bigmem_16.q'
if hasattr(options,'numprocess'):
    numprocess = options.numprocess
else:
    numprocess = '1'
if numprocess==None:
    numprocess = '1'

if hasattr(options,'numcores'):
    numcores = options.numcores
else:
    numcores = '1'
if numcores==None:
    numcores = '1'
PARALLEL_ENVIRONMENT = load_parrallel_settings()
if PARALLEL_ENVIRONMENT!="":
    PARALLEL_ENVIRONMENT = "-pe "+PARALLEL_ENVIRONMENT + " " + numcores
if hasattr(options,'reference'):
    reference = options.reference
else:
    reference = ''
if hasattr(options,'keepmean'):
    keepmean = options.keepmean
else: 
    keepmean = "0"
if keepmean==None:
    keepmean = "0"
if hasattr(options,'voxelsize'):
    voxelsize = options.voxelsize
if hasattr(options,'contrast'):
    contrast = options.contrast
if hasattr(options,'TR_MSEC'):
    TR_MSEC = options.TR_MSEC
if TR_MSEC==None:
    TR_MSEC = "None"
if hasattr(options,'dospnormfirst'):
    dospnormfirst = options.dospnormfirst
if dospnormfirst==None:
    dospnormfirst = False
if hasattr(options,'environment'):
    environment = options.environment
else:
    environment = "octave"
if environment==None:
    environment = "octave"
if hasattr(options,'analysis'):
    analysis = options.analysis
else: 
    analysis = None
if analysis==None:
    print "WARNING: without an analysis model (specified by switch -a), no optimization will be performed"
    print "PRONTO only generates the preprocessed data"
    contrast="nocontrast"
    analysis="None"
if hasattr(options,'metric'):
    metric = options.metric
else: 
    metric = ""
if metric==None:
    print "WARNING: optimization metric not defined (specified by switch -m), will use default dPR criterion"
    metric="dPR"
if hasattr(options,'convolve'):
    convolve = options.convolve
else: 
    convolve = None
if convolve==None:
    convolve = "None"
if hasattr(options,'decision_model'):
    decision_model = options.decision_model
else: 
    decision_model = "None"
if decision_model==None:
    decision_model = "None"
if hasattr(options,'drf'):
    drf = options.drf
else: 
    drf = "None"
if drf==None:
    drf = "None"
if hasattr(options,'Nblock'):
    Nblock = options.Nblock
else: 
    Nblock = None
if Nblock==None:
    Nblock = "None"
if hasattr(options,'WIND'):
    WIND = options.WIND
else: 
    WIND = None
if WIND==None:
    WIND = "None"
if hasattr(options,'num_PCs'):
    num_PCs = options.num_PCs
else: 
    num_PCs = None
if num_PCs==None:
    num_PCs = "None"
if hasattr(options,'subspace'):
    subspace = options.subspace
else: 
    subspace = None
if subspace==None:
    subspace = "None"
if hasattr(options,'spm'):
    spm = options.spm
else: 
    spm = None
if spm==None:
    spm = "None"
if hasattr(options,'N_resample'):
    N_resample = options.N_resample
else:
    N_resample  = '10'
if N_resample==None:
    N_resample = "10"	 	
if hasattr(options,'DEOBLIQUE'):
    DEOBLIQUE = options.DEOBLIQUE
else:
    DEOBLIQUE  = False
if DEOBLIQUE==True:
    DEOBLIQUE = "1"
else:
    DEOBLIQUE = "0"
if hasattr(options,'TPATTERN'):
    TPATTERN = options.TPATTERN
else:
    TPATTERN  = "None"
if TPATTERN==None:
    TPATTERN = "None"

if hasattr(options,'memory'):
    memory = "-l mem="+options.memory
else:
    memory  = ""
if memory==None:
    memory  = ""

###############  Checking the switches
if (analysis.upper()=="LDA") and (drf=="None"):
    print "WARNING (Deprecated usage): --drf switch not defined for LDA model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="GNB") and (decision_model=="None"):
    print "WARNING (Deprecated usage): --decision_model switch not defined for GNB model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERCVA") and (drf=="None"):
    print "WARNING (Deprecated usage): --drf switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERCVA") and (Nblock=="None"):
    print "WARNING (Deprecated usage): --Nblock switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERCVA") and (WIND=="None"):
    print "WARNING (Deprecated usage): --WIND switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERCVA") and (subspace=="None"):
    print "WARNING (Deprecated usage): --subspace switch not defined for erCVA model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERGNB") and (Nblock=="None"):
    print "WARNING (Deprecated usage): --Nblock switch not defined for erGNB model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="ERGNB") and (WIND=="None"):
    print "WARNING (Deprecated usage): --WIND switch not defined for erGNB model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="SCONN") and (spm=="None"):
    print "WARNING (Deprecated usage): --spm switch has to be used with the SCONN model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="GLM") and (convolve=="None"):
    print "WARNING (Old style usage): --convolve switch has to be used with the GLM model. PRONTO will check TASK files for parameter(s)"
if (analysis.upper()=="GPCA") and (num_PCs=="None"):
    print "WARNING (Deprecated usage): --num_PCs switch not defined for gPCA model. PRONTO will check TASK files for parameter(s)"

if not (convolve in ["1","0","None"]):
    print "WARNING (Deprecated usage): --convolve has to be 0 or 1"
if not (decision_model.lower() in ["linear","nonlinear","none"]):
    print "WARNING (Deprecated usage): --decision_model has to be linear or nonlinear"
if not (subspace.lower() in ["onecomp","multicomp","none"]):
    print "WARNING (Deprecated usage): --subspace has to be onecomp or multicomp"
if not (spm.lower() in ["corr","zscore","none"]):
    print "WARNING (Deprecated usage): --spm has to be corr or zscore"

print N_resample
model_parameter = "keepmean "+keepmean+" convolve "+convolve+" decision_model "+decision_model+" drf "+drf+" Nblock "+Nblock+" WIND "+WIND+" subspace "+subspace+" spm "+spm+" N_resample "+N_resample+" TR_MSEC "+TR_MSEC+" num_PCs "+num_PCs

#################### END checking the parameters

os.environ["PIPELINE_NUMBER_OF_CORES"] = numcores
os.environ["FSLOUTPUTTYPE"] = "NIFTI"

if part==0:
    job_id_str = "-hold_jid "
else:
    job_id_str  = ""

job_id_str3_0 = ""
job_id_list_nosge = []
if dospnormfirst:
    if reference==None:
        print "Switch -r has to be used for part 3"
    exit(0)
    
    job_id_str3_0 = "-hold_jid "
    used_cores = 0

    # Initial spatial normalization: estimate transform
    for job_counter in list_of_unique_subjects:
        temp_file = input_files_temp + '/' + "%04d.txt" % job_counter
        # submit into queue
        job_name = "spest%04d" % job_counter
        cmd = "qsub -q {0} {10} {11} -N {1} -wd {2} -j y -b y 'python {2}/scripts_other/spatial_normalization_wrapper.py  -n {3} -i {4} -r {5} -v \"{6}\" -s 1 -e \"{7}\" --DEOBLIQUE {8}' >{9}".format(queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file,PARALLEL_ENVIRONMENT,memory)
        
    
        if noSGE == True:
            id_no      = "noSGE"
            job_id_file = job_name+".nosge.txt"
            cmd = "python {2}/scripts_other/spatial_normalization_wrapper.py  -n {3} -i {4} -r {5} -v \"{6}\" -s 1 -e \"{7}\" --DEOBLIQUE {8} >{9}".format(queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file)
            if not job_counter==list_of_unique_subjects[-1]:
                cmd = cmd + " &"
            if os.path.isfile("current_job_nosge.txt"):
                os.remove("current_job_nosge.txt")
            os.system(cmd)
            while(not os.path.isfile("current_job_nosge.txt")):
                sleep(0.1)
            procfile = open("current_job_nosge.txt","r");
            procid_str = procfile.read();
            pid_no = [int(s) for s in procid_str.split() if s.isdigit()]
            used_cores = used_cores + 1
            print pid_no[0]
            job_id_list_nosge.append(pid_no[0])
            print "Your initial spatial normalization (transform estimation) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
            used_cores = wait_until_jobs_done(job_id_list_nosge,used_cores,numprocess)
            
        else:
            os.system(cmd)
            jobfile = open(job_id_file,"r");
            jobid_str = jobfile.read();
            id_no = [int(s) for s in jobid_str.split() if s.isdigit()]
            job_id_str3_0 = job_id_str3_0 + str(id_no[0]) + ","
            print "Your initial spatial normalization (transform estimation) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
        #execute command

    if noSGE == False:
        job_id_str3_0 = job_id_str3_0[:-1]
        os.remove(job_id_file)   

job_id_list_nosge = []
if part==0 or part==1:
    # check whether the pipeline file contain correct steps
    check_pipeline_file(options.pipeline)       
    print "Submitting the estimation jobs:"
    used_cores = 0
    for job_counter in range(1,current_subject+1):
        temp_file = input_files_temp + '/' + "%04d.txt" % job_counter
        # submit into queue
        job_name = "pipopt%04d" % job_counter
        cmd = "qsub {0} -q {1} {15} {16} -N {2} -wd {3} -j y -b y 'python {3}/scripts_other/pipeline_wrapper.py -n {4} -i {5} -c {6} -a {7} --contrast \"{8}\" --dospnormfirst {9} -e \"{10}\" --modelparam \"{11}\" --DEOBLIQUE {12} --TPATTERN {13}' >{14}".format(job_id_str3_0,queue_name,job_name,codefull_path,numcores,temp_file, options.pipeline, options.analysis,contrast,dospnormfirst+0,environment,model_parameter,DEOBLIQUE,TPATTERN,job_id_file,PARALLEL_ENVIRONMENT,memory)
        if noSGE == True:
            id_no = "noSGE"
            job_id_file = job_name + ".nosge.txt"
            cmd = "python {3}/scripts_other/pipeline_wrapper.py -n {4} -i {5} -c {6} -a {7} --contrast \"{8}\" --dospnormfirst {9} -e \"{10}\" --modelparam \"{11}\" --DEOBLIQUE {12} --TPATTERN {13} >{14}".format(job_id_str3_0,queue_name,job_name,codefull_path,numcores,temp_file, options.pipeline, options.analysis,contrast,dospnormfirst+0,environment,model_parameter,DEOBLIQUE,TPATTERN,job_id_file)
            if job_counter<current_subject:
                cmd = cmd + " &"
            if os.path.isfile("current_job_nosge.txt"):
                os.remove("current_job_nosge.txt")
            os.system(cmd)
            while(not os.path.isfile("current_job_nosge.txt")):
                sleep(0.1)
            procfile = open("current_job_nosge.txt","r");
            procid_str = procfile.read();
            pid_no = [int(s) for s in procid_str.split() if s.isdigit()]
            used_cores = used_cores + 1
            print pid_no[0]
            job_id_list_nosge.append(pid_no[0])
            print "     Your pipeline optimization job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
            used_cores = wait_until_jobs_done(job_id_list_nosge,used_cores,numprocess)
        else:
            os.system(cmd)
        
        
        if noSGE == False:
            jobfile = open(job_id_file,"r");
            jobid_str = jobfile.read();
            id_no = [int(s) for s in jobid_str.split() if s.isdigit()]
            job_id_str = job_id_str + str(id_no[0]) + ","

            print "     Your pipeline optimization job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
            #execute command
    if noSGE == False:
        job_id_str = job_id_str[:-1]
        os.remove(job_id_file)

job_id_list_nosge = []
if ((part==0 or part==2) and analysis.upper()!="NONE"):
    # for pipeline step 3 (post-processing optimization)
    print "Submitting the post-processing optimization jobs:"
    if part==0:
        print "(The post-processing optimization jobs will be on hold until the estimation jobs are finished)"
    # declare pipeline step-3;  NB replaced "-q all.q" with "-q bigmem.q" to call all nodes
    
    cmd = "qsub {0} -q {1} {9} {10} -N pipopt_final -wd {2} -j y -b y 'python {2}/scripts_other/optimization_wrapper.py  -n {3} -i {4} -m {5} -k {6} -e \"{7}\"' >{8}".format(job_id_str,queue_name,codefull_path,numcores,input_file_for_part2, metric,keepmean,environment,job_id_file,PARALLEL_ENVIRONMENT,memory)
    if noSGE == True:
        id_no = "noSGE"
        job_id_file = "pipopt_final.nosge.txt"
        print "     Your pipeline post-processing optimization job #{0} has been submitted".format(id_no)
        cmd = "python {2}/scripts_other/optimization_wrapper.py  -n {3} -i {4} -m {5} -k {6} -e \"{7}\" >{8}".format(job_id_str,queue_name,codefull_path,numcores,input_file_for_part2, metric, keepmean,environment,job_id_file)

    os.system(cmd)
    if noSGE == False:
        jobfile = open(job_id_file,"r");
        jobid_str = jobfile.read();
        id_no = [int(s) for s in jobid_str.split() if s.isdigit()]
        job_id_str = "-hold_jid "+str(id_no[0])
    
    print "     Your pipeline post-processing optimization job #{0} has been submitted".format(id_no)
    
    #execute command

job_id_list_nosge = []    
if (part==0 or part==3) and not dospnormfirst:
    job_id_str3_1 = "-hold_jid "
    job_id_str3_2 = "-hold_jid "

    if reference==None:
        print "Switch -r has to be used for part 3 (spatial normalization)"
        exit(0)
    print "Submitting the spatial normalization jobs:"
    # extract runs with similar structrul data
    if part==0:
            print "(The spatial normalization jobs will be on hold until the post processing jobs are finished)"
    
     # for pipeline step 3 (spatial normalization: apply transform)
    used_cores = 0
    for job_counter in list_of_unique_subjects:
        temp_file = input_files_temp + '/' + "%04d.txt" % job_counter
        # submit into queue
        job_name = "spest%04d" % job_counter
        
        cmd = "qsub {0} -q {1} {11} {12} -N {2} -wd {3} -j y -b y 'python {3}/scripts_other/spatial_normalization_wrapper.py  -n {4} -i {5} -r {6} -v \"{7}\" -s 1 -e \"{8}\" --DEOBLIQUE {9}' >{10}".format(job_id_str,queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file,PARALLEL_ENVIRONMENT,memory)
        if noSGE == True:
            id_no = "noSGE"
            job_id_file = job_name +".nosge.txt"
            cmd = "python {3}/scripts_other/spatial_normalization_wrapper.py  -n {4} -i {5} -r {6} -v \"{7}\" -s 1 -e \"{8}\" --DEOBLIQUE {9} >{10}".format(job_id_str,queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file)      
            if not job_counter==list_of_unique_subjects[-1]:
                cmd = cmd + " &"
            if os.path.isfile("current_job_nosge.txt"):
                os.remove("current_job_nosge.txt")
            os.system(cmd)
            while(not os.path.isfile("current_job_nosge.txt")):
                sleep(0.1)
            procfile = open("current_job_nosge.txt","r");
            procid_str = procfile.read();
            pid_no = [int(s) for s in procid_str.split() if s.isdigit()]
            used_cores = used_cores + 1
            print pid_no[0]
            job_id_list_nosge.append(pid_no[0])
            print "     Your spatial normalization (transform estimation) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
            used_cores = wait_until_jobs_done(job_id_list_nosge,used_cores,numprocess)
        else:
            os.system(cmd)
            jobfile = open(job_id_file,"r");
            jobid_str = jobfile.read();
            id_no = [int(s) for s in jobid_str.split() if s.isdigit()]
            job_id_str3_1 = job_id_str3_1 + str(id_no[0]) + ","
            print "     Your spatial normalization (transform estimation) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
    if noSGE == False:  
        job_id_str3_1 = job_id_str3_1[:-1]
        os.remove(job_id_file)   
    
    used_cores = 0
    # for pipeline step 3 (spatial normalization: apply transform)
    for job_counter in range(1,current_subject+1):
        temp_file = input_files_temp + '/' + "%04d.txt" % job_counter
        # submit into queue
        job_name = "spnorm%04d" % job_counter
        cmd = "qsub {0} -q {1} {11} {12} -N {2} -wd {3} -j y -b y 'python {3}/scripts_other/spatial_normalization_wrapper.py  -n {4} -i {5} -r {6} -v \"{7}\" -s 2 -e \"{8}\" --DEOBLIQUE {9}' >{10}".format(job_id_str3_1,queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file,PARALLEL_ENVIRONMENT,memory)
        if noSGE == True:
            id_no = "noSGE"
            job_id_file = job_name + ".nosge.txt"
            cmd = "python {3}/scripts_other/spatial_normalization_wrapper.py  -n {4} -i {5} -r {6} -v \"{7}\" -s 2 -e \"{8}\" --DEOBLIQUE {9} >{10}".format(job_id_str3_1,queue_name,job_name,codefull_path,numcores,temp_file,reference,voxelsize,environment,DEOBLIQUE,job_id_file)        
            if job_counter<current_subject:
                cmd = cmd + " &"
            if os.path.isfile("current_job_nosge.txt"):
                os.remove("current_job_nosge.txt")
            os.system(cmd)
            while(not os.path.isfile("current_job_nosge.txt")):
                sleep(0.1)
            procfile = open("current_job_nosge.txt","r");
            procid_str = procfile.read();
            pid_no = [int(s) for s in procid_str.split() if s.isdigit()]
            used_cores = used_cores + 1
            print pid_no[0]
            job_id_list_nosge.append(pid_no[0])
            print "     Your spatial normalization (apply transform) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter
            used_cores = wait_until_jobs_done(job_id_list_nosge,used_cores,numprocess)            
        else:
            os.system(cmd)
            jobfile = open(job_id_file,"r");
            jobid_str = jobfile.read();
            id_no = [int(s) for s in jobid_str.split() if s.isdigit()]
            job_id_str3_2 = job_id_str3_2 + str(id_no[0]) + ","
            print "     Your spatial normalization (apply transform) job #{0} (Subject %04d) has been submitted".format(id_no) % job_counter

        #execute command
    if noSGE == False:     
        job_id_str3_2 = job_id_str3_2[:-1]
        os.remove(job_id_file)   
    
    print "Submitting the group mask generation job:"
    cmd = "qsub {0} -q {1} {6} {7} -N pipopt_maskgen -wd {2} -j y -b y 'python {2}/scripts_other/mask_tissue_wrapper.py  -n {3} -i {4} -e \"{5}\"' ".format(job_id_str3_2,queue_name,codefull_path,numcores,options.inputdata, environment,PARALLEL_ENVIRONMENT,memory)  
    if noSGE == True:
        cmd = "python {2}/scripts_other/mask_tissue_wrapper.py  -n {3} -i {4} -e \"{5}\" > pipopt_maskgen.nosge.txt".format(job_id_str3_2,queue_name,codefull_path,numcores,options.inputdata, environment)  

    os.system(cmd)



