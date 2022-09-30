'''
Created on May 7, 2010

@author: wolf
'''

import os
import sys
import options as opt

#create a class object
optobj = opt.tpsoptions()


class qsubsystem(object):
    def __init__(self):
        self.wolf = "test"
    
    def writeKernelQsubFile(self,basedir,kernel,pythonpath,mode,dirstring,reverse=False,method="tis",qsubs="Standard",walltime="100:00:00",queuename="serial",tispath="",wrappath=""):
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
        elif qsubs=="local": 
            if reverse:
                executable = os.path.join(tispath,"GTISkernelR.py")
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev\n"
                pe = "std.rev.err\n" 
                po = "std.rev.out\n" 
                
            else:
                if method == "tis":
                    executable = os.path.join(tispath,"GTISkernel.py")
                else:
                    executable = os.path.join(tispath,"GTPSkernel.py")
                
                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = "k"+method+"-%03d\n" % kernel
                pe = "std.%03d.err" % kernel
                po = "std.%03d.out" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("#!/bin/bash\n")
	    wf.write("source $HOME/.bashrc\n")
	    wf.write("export PYTHONPATH=$PYTHONPATH:"+wrappath+"\n")
	    wf.write("hostname\n")

            wf.write("cd "+ basedir+"\n")
	    wf.write("echo \" \"\n")
	    wf.write("echo \"directory\"\n")
	    wf.write("pwd\n")
            
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + " > " + po + " 2>" + pe + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + " > " + po + " 2>" + pe + "\n")
            
            wf.close()
        elif qsubs=="vulcan":        
            if reverse:
                executable = os.path.join(tispath,"GTISkernelR.py")
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev\n"
                pe = "sge.rev.err\n"
                po = "sge.rev.out\n"

            else:
                if method == "tis":
                    executable = os.path.join(tispath,"GTISkernel.py")
                else:
                    executable = os.path.join(tispath,"GTPSkernel.py")

                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = (("%s"+"-%03d\n") % (optobj.runoptions["jobname"],kernel))
                pe = "sge.%03d.err\n" % kernel
                po = "sge.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            wf.write("#!/bin/bash\n")
            wf.write("#$ -N " + qn)
            wf.write("#$ -S /bin/bash\n")
            wf.write("#$ -o "+po)
            wf.write("#$ -e "+pe)
            wf.write("#$ -r n\n")

            specificq = True
            
            if ((queuename=='any.q') or (queuename=='any')):
                specificq = False

            if specificq:
                wf.write("#$ -q " + queuename + "\n")

            if mode == "tis":
                wf.write("#$ -l h_rt=15:50:00\n")

            wf.write("#$ -pe smp 1\n")
            wf.write("#$ -P ams.p\n")
            wf.write("#$ -cwd\n\n")

            wf.write("source $HOME/.bashrc\n")
            wf.write("export PYTHONPATH=$PYTHONPATH:"+wrappath+"\n")
#            wf.write("module unload 64/intel/12.0.0.08\n")
            wf.write("module load intel/11.1046\n")
            wf.write("module load mpi/intel/openmpi/1.3.3\n")               #DS for using lmp-plumed-joshua-cubic-bussi-nvt-npt executable (bussinpt barostat)
            wf.write("hostname\n\n")
            
            wf.write("cd "+ basedir+"\n")
            
            wf.write("echo \" \"\n")
            wf.write("echo \"directory\"\n")
            wf.write("pwd\n\n")
            
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")

            wf.close()  

        elif qsubs=="vulcan-parallel":        
            if reverse:
                executable = os.path.join(tispath,"GTISkernelR.py")
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev\n"
                pe = "sge.rev.err\n"
                po = "sge.rev.out\n"

            else:
                if method == "tis":
                    executable = os.path.join(tispath,"GTISkernel.py")
                else:
                    executable = os.path.join(tispath,"GTPSkernel.py")

                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = (("%s"+"-%03d\n") % (optobj.runoptions["jobname"],kernel))
                pe = "sge.%03d.err\n" % kernel
                po = "sge.%03d.out\n" % kernel
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            specificq = True
            #srm beter queue handlinf for vulcan parallel
            dum = queuename.split('-')
            if dum[0]=='serial':
                queue = 'serial.q'
                nslots = int(dum[1])
            elif dum[0]=='shorttime':
                queue = 'shorttime.q'
                nslots = int(dum[1])
            elif dum[0]=='edu':
                queue = 'edu.q'
                nslots = int(dum[1])

            #SRM add any queue support
            
            elif dum[0]=='any':
                queue = 'serial.q'
                specificq = False
                nslots = int(dum[1])

            else:
                print 'Wrong queue name for vulcan parallel!'
                print 'Exiting program...!'
                sys.exit(1)
            if nslots>8:
                print 'number of cores requested greater than eight'
                print 'Exiting program...!'
                sys.exit(1)


            wf.write("#!/bin/bash\n")
            wf.write("#$ -N " + qn)
            wf.write("#$ -cwd\n\n")
            wf.write("#$ -S /bin/bash\n")
            wf.write("#$ -o "+po)
            wf.write("#$ -e "+pe)
            wf.write("#$ -r n\n")
            if specificq:
                wf.write("#$ -q " + queue + "\n")
            wf.write("#$ -l h_rt=15:50:00\n")
            wf.write("#$ -pe smp " + str(nslots) + "\n")
            wf.write("#$ -P ams.p\n")
            

            wf.write("source $HOME/.bashrc\n")
            wf.write("export PYTHONPATH=$PYTHONPATH:"+wrappath+"\n")
#            wf.write("module unload 64/intel/12.0.0.08\n")
            wf.write("module load intel/016.0.047\n")
            wf.write("module load mpi/intel/openmpi/1.6.4\n")               #DS for using lmp-plumed-joshua-cubic-bussi-nvt-npt executable (bussinpt barostat)
            wf.write("hostname\n\n")
            
            wf.write("cd "+ basedir+"\n")
            
            wf.write("echo \" \"\n")
            wf.write("echo \"directory\"\n")
            wf.write("pwd\n\n")
            
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")

            wf.close()  

        elif qsubs=="slurm":        
            if reverse:
                #wont work!!!
                executable = os.path.join(tispath,"GTISkernelR.py")
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev\n"
                pe = "sge.rev.err\n"
                po = "sge.rev.out\n"

            else:
                if method == "tis":
                    executable = os.path.join(tispath,"GTISkernel.py")
                else:
                    executable = os.path.join(tispath,"GTPSkernel.py")

                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = (("%s"+"-%03d\n") % (optobj.runoptions["jobname"],kernel))
                pe = "sge.%03d.err\n" % kernel
                po = "sge.%03d.out\n" % kernel
            
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            
            
            if queuename=='shorttime':
                queue = 'shorttime'
                nslots = 1
            else:
                print 'Wrong queue name for slurm!'
                print 'Exiting program...!'
                sys.exit(1)


            wf.write("#!/bin/bash\n")
            wf.write("#SBATCH --job-name=" + qn)
            wf.write("#SBATCH --time=23:59:00\n")
            wf.write("#SBATCH --partition " + queue + "\n")
            wf.write("#SBATCH --ntasks=" + str(nslots) + "\n")
            wf.write("#SBATCH --mem-per-cpu=3GB\n")
            #wf.write("#SBATCH --ntasks-per-node=" + str(16) + "\n")
            wf.write("#SBATCH --hint=nomultithread\n")
            #wf.write("#SBATCH --exclude=zgh040\n")
            wf.write("#SBATCH -o "+po)
            wf.write("#SBATCH -e "+pe)
            
            wf.write("source $HOME/.bashrc\n")
            wf.write("export PYTHONPATH=$PYTHONPATH:"+wrappath+"\n")
            wf.write("module purge\n")
            #wf.write("module load mpi/openmpi/intel/1.10.4\n")
            #wf.write("module load mpi/openmpi/gcc/1.10.4\n")               #DS for using lmp-plumed-joshua-cubic-bussi-nvt-npt executable (bussinpt barostat)
            wf.write("module load gcc/7.3.0\n")
            wf.write("module load mpi/openmpi/gcc/1.10.7-shorttime\n")
            wf.write("hostname\n\n")
            
            wf.write("cd "+ basedir+"\n")
            
            wf.write("echo \" \"\n")
            wf.write("echo \"directory\"\n")
            wf.write("pwd\n\n")
	    wf.write("uss=$(whoami)\n")
            wf.write("find /dev/shm/ -user $uss -type f -mmin +30 -delete\n")
            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")

            wf.close()


        elif qsubs=="slurm-parallel":        
            if reverse:
                #wont work!!!
                executable = os.path.join(tispath,"GTISkernelR.py")
                filename = os.path.join(basedir,"kernelqsub.reverse.sub")
                sk = ""
                qn = "kernel-rev\n"
                pe = "sge.rev.err\n"
                po = "sge.rev.out\n"

            else:
                if method == "tis":
                    executable = os.path.join(tispath,"GTISkernel.py")
                else:
                    executable = os.path.join(tispath,"GTPSkernel.py")

                filename = os.path.join(basedir,"kernelqsub.%03d.sub" % kernel)
                sk = str(kernel)
                qn = (("%s"+"-%03d\n") % (optobj.runoptions["jobname"],kernel))
                pe = "sge.%03d.err\n" % kernel
                po = "sge.%03d.out\n" % kernel
            
            wf = open(filename,"w")
            pythoncommand = os.path.join(pythonpath,"python")
            
            dum = queuename.split('-')
            if dum[0]=='shorttime':
                queue = 'shorttime'
                nslots = int(dum[1])
            else:
                print 'Wrong queue name for vulcan parallel!'
                print 'Exiting program...!'
                sys.exit(1)
            if nslots>55:
                print 'number of cores requested greater than eight'
                print 'Exiting program...!'
                sys.exit(1)


            wf.write("#!/bin/bash\n")
            wf.write("#SBATCH --job-name=" + qn)
            wf.write("#SBATCH --time=23:59:00\n")
            wf.write("#SBATCH --partition " + queue + "\n")
            wf.write("#SBATCH --ntasks=" + str(nslots) + "\n")
            wf.write("#SBATCH --mem-per-cpu=3GB\n")
            #wf.write("#SBATCH --ntasks-per-node=" + str(16) + "\n")
            wf.write("#SBATCH --hint=nomultithread\n")
            #wf.write("#SBATCH --exclude=zgh040\n")
	    wf.write("#SBATCH -o "+po)
            wf.write("#SBATCH -e "+pe)
            
            wf.write("source $HOME/.bashrc\n")
            wf.write("export PYTHONPATH=$PYTHONPATH:"+wrappath+"\n")
            wf.write("module purge\n")
            #wf.write("module load mpi/openmpi/intel/1.10.4\n")
            #wf.write("module load mpi/openmpi/gcc/1.10.4\n")               #DS for using lmp-plumed-joshua-cubic-bussi-nvt-npt executable (bussinpt barostat)
            wf.write("module load gcc/7.3.0\n")
            wf.write("module load mpi/openmpi/gcc/1.10.7-shorttime\n")
            wf.write("hostname\n\n")
            
            wf.write("cd "+ basedir+"\n")
            
            wf.write("echo \" \"\n")
            wf.write("echo \"directory\"\n")
            wf.write("pwd\n\n")
            wf.write("uss=$(whoami)\n")
            wf.write("find /dev/shm/ -user $uss -type f -mmin +30 -delete\n")

            if mode == "initial":
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=initial " + "kernel="+sk + "\n")
            else:
                wf.write(pythoncommand + " " + executable + " basedir="+basedir + " mode=" + method + " " +     "kernel="+sk + " " + "dirnumber="+dirstring + "\n")

            wf.close()  

        return filename
