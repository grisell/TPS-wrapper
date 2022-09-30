'''
Created on May 7, 2010

@author: wolf
'''

import os


class qsubsystem(object):
    def __init__(self):
        self.wolf = "test"
    
    def writeKernelQsubFile(self,basedir,kernel,pythonpath,mode,dirstring,reverse=False,method="tis",qsubs="Standard",walltime="100:00:00",queuename="serial"):
        if qsubs=="Standard":        
            if reverse:
                executable = "GTISkernelR.py"
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev"
                pe = "pbs.rev.err\n" 
                po = "pbs.rev.out\n" 
                
            else:
                if method == "tis":
                    executable = "GTISkernel.py"
                else:
                    executable = "GTPSkernel.py"
                
                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = "k"+method+"-%03d\n" % kernel
                pe = "pbs.%03d.err\n" % kernel
                po = "pbs.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("#!/bin/bash\n")
            wf.write("#$ -N " + qn)
            wf.write("#$ -S /bin/bash\n")
            wf.write("#$ -o "+po)
            wf.write("#$ -e "+pe)
            wf.write("#$ -r n\n")
            wf.write("#$ -V\n")
            wf.write("#$ -q +" + queuename+ "\n")
            wf.write("#$ -pe mpi 8\n")
            wf.write("#$ -cwd\n")
            wf.write("cd "+ basedir+"\n")
            
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")
            
            wf.close()
        elif qsubs=="PBS":
            if reverse:
                executable = "GTISkernelR.py"
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev"
                pe = "pbs.rev.err\n" 
                po = "pbs.rev.out\n" 
                
            else:
                if method == "tis":
                    executable = "GTISkernel.py"
                else:
                    executable = "GTPSkernel.py"
                
                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = "k"+method+"-%03d\n" % kernel
                pe = "pbs.%03d.err\n" % kernel
                po = "pbs.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("#PBS -S /bin/bash\n")
            wf.write("#PBS -lnodes=1:ppn=8:"+queuename+"\n")
            wf.write("#PBS -lwalltime="+walltime+"\n")
            wf.write("#PBS -o "+po+"\n")
            wf.write("#PBS -e "+pe+"\n")
            wf.write("cd "+ basedir+"\n")
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")
            
            wf.close()
        elif qsubs=="lisa":
            if reverse:
                executable = "GTISkernelR.py"
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev"
                pe = "pbs.rev.err\n" 
                po = "pbs.rev.out\n" 
                
            else:
                if method == "tis":
                    executable = "GTISkernel.py"
                else:
                    executable = "GTPSkernel.py"
                
                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = "k"+method+"-%03d\n" % kernel
                pe = "pbs.%03d.err\n" % kernel
                po = "pbs.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("#PBS -S /bin/bash\n")
            wf.write("#PBS -lnodes=1:ppn=8:"+queuename+"\n")
            wf.write("#PBS -lwalltime="+walltime+"\n")
            wf.write("#PBS -o "+po+"\n")
            wf.write("#PBS -e "+pe+"\n")
            wf.write("cd "+ basedir+"\n")
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")
            
            wf.close()
             
        elif qsubs=="LL":
            if reverse:
                executable = "GTISkernelR.py"
                filename = os.path.join(basedir,"kernelll.reverse.sub")
                sk = ""
                qn = "kernel-rev"
                pe = "ll.rev.err\n" 
                po = "ll.rev.out\n" 
        
            else:
                if method == "tis":
                    executable = "GTISkernel.py"
                else:
                    executable = "GTPSkernel.py"

                filename = os.path.join(basedir,"kernelll.%03d.sub" % kernel)
                sk = str(kernel)
                qn = "k"+method+"-%03d\n" % kernel
                pe = "ll.%03d.err\n" % kernel
                po = "ll.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("# Loadleveler directives start with # @\n")
            wf.write("#\n")
            wf.write("# Request one node with 64 processes.\n")
            wf.write("# @ node = 1\n") 
            wf.write("# @ tasks_per_node = 64\n")
            wf.write("#\n")
            wf.write("# Loadleveler can send email, but for this job, we ask Loadleveler not \n")
            wf.write("#     to send any email:\n")
            wf.write("#\n")
            wf.write("# @ notification = never\n")
            wf.write("#\n")
            wf.write("# Define the standard input, output and error for the job:\n")
            wf.write("#\n") 
            wf.write("# @ input = /dev/null\n")
            wf.write("# @ output = "+po)
            wf.write("# @ error = "+pe)
            wf.write("#\n") 
            wf.write("# @ wall_clock_limit = "+walltime+"\n")
            wf.write("#\n")
            wf.write("# Tell Loadleveler that this is a parallel job. Loadleveler will take\n")
            wf.write("# care of defining the proper environment for the job\n")
            wf.write("# \n")
            wf.write("# @ job_type = parallel\n")
            wf.write("#\n")
            wf.write("# The following line ensures that the communication between nodes will\n")
            wf.write("# use the infiniband hardware:\n")
            wf.write("#\n")
            wf.write("# @ network.MPI = sn_all,not_shared,US\n")
            wf.write("# \n")
            wf.write("# This ends the instructions for Loadleveler:\n") 
            wf.write("#\n")
            wf.write("# @ queue\n")
            wf.write("#\n")
            wf.write("# Here the shell script starts. \n")
            wf.write("# go to the working directory: \n")
            wf.write("# \n")
            wf.write("cd "+ basedir+"\n")
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")
            
            wf.close()

        return filename
