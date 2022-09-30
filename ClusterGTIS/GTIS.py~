'''
Created on May 3, 2010

@author: Wolfgang Lechner
@status: beta 
'''
import pygromacstps.tis as gtis
#python stuff
import sys
import os
import random
import time
import shutil


def shootingKernel(dirnumber,mtispath=""):
    for i in range(sim.kernels.nkernels):
        pythonpath= sim.paths[0].options.paths["pythonpath"]
        mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
        mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
        queuename = sim.paths[0].options.runoptions["queuename"]
	mwrappath = sim.paths[0].options.paths["wrapperpath"]
        
        filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"tis",str(dirnumber),reverse=False,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)
        if mqsubs == "LL":
            os.system("llsubmit " + filename)
	elif mqsubs == "local":
	    os.system("chmod u+x " + " " + filename)
            pe = "std.head.%03d.err\n" % i
            po = "std.head.%03d.out\n" % i
            os.system(filename + " >" + po + " " + " 2>" + pe)
        elif mqsubs == "local-slurm":
	    os.system("chmod u+x " + " " + filename)
            pe = "std.head.%03d.err\n" % i
            po = "std.head.%03d.out\n" % i
            os.system(filename + " >" + po + " " + " 2>" + pe)
	elif mqsubs == "slurm" or mqsubs == "slurm-parallel" :
            os.system("sbatch " + " " + filename)
        else:
            os.system("qsub " + filename)
        
def shootingKernelReverse(rdirnumber,mtispath=""):
    pythonpath= sim.paths[0].options.paths["pythonpath"]
    mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
    mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
    queuename = sim.paths[0].options.runoptions["queuename"]
    mwrappath = sim.paths[0].options.paths["wrapperpath"]
    
    filename = sim.qsubsystem.writeKernelQsubFile(basedir, 0,pythonpath,"tis",str(rdirnumber),reverse=True,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)
    if mqsubs == "LL":
        os.system("llsubmit " + filename)
    elif mqsubs == "local":
	os.system("chmod u+x " + " " + filename)
        pe = "std.head.%03d.err\n" % i
        po = "std.head.%03d.out\n" % i
        os.system(filename + " >" + po + " " + " 2>" + pe)
    elif mqsubs == "local-slurm":
	os.system("chmod u+x " + " " + filename)
        pe = "std.head.%03d.err\n" % i
        po = "std.head.%03d.out\n" % i
        os.system(filename + " >" + po + " " + " 2>" + pe)
    elif mqsubs == "slurm" or mqsubs == "slurm-parallel":
            os.system("sbatch " + filename)
    else:
        os.system("qsub " + filename)
    

def allFinished(dirstring):
    allF = True
    for i in sim.kernels.kernelPaths:
        lapath = sim.paths[i].nfsladir
        newladir = os.path.join(lapath,dirstring)
        filename = os.path.join(newladir,"finished")
        if not(os.path.exists(filename)):
            allF = False
    return allF

def allMoveDir(dirstring,newdirstring):                                        #DS for changing the directories name
    for i in sim.kernels.kernelPaths:
        lapath = sim.paths[i].nfsladir
        removedir = os.path.join(lapath,dirstring)
        movedir = os.path.join(lapath,newdirstring)
        shutil.rmtree(removedir)                                               #DS delete the old directory
        shutil.move(movedir,removedir)                                         #DS move the new directory to the old one
        fileold = os.path.join(removedir,"trajectory.00."+newdirstring+".dat") #DS after moving the directory we must rename the files acording to the directory name
        fileoldtrial = os.path.join(removedir,"trajectory.00."+newdirstring+".trial") 
        filenew = os.path.join(removedir,"trajectory.00."+dirstring+".dat")
        filenewtrial = os.path.join(removedir,"trajectory.00."+dirstring+".trial")
        shutil.move(fileold,filenew)
        shutil.move(fileoldtrial,filenewtrial)


def allFinishedReverse(dirstring):
    allF = True
    for i in range(len(sim.reversePaths)):
        lapath = sim.reversePaths[i].nfsladir
        newladir = os.path.join(lapath,dirstring)
        filename = os.path.join(newladir,"finished")
        if not(os.path.exists(filename)):
            allF = False
    return allF

def exchangeKernel(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
    sim.log.log.info(dirstring + " exk: starting trialReplicaExchange ")  
    exchange = sim.trialReplicaExchange(dirstring,newdirstring)
    sim.log.log.info(dirstring + " exk: reevaluating trajectories ")  
    for i in sim.kernels.kernelPaths:
        sim.getFullTrajectory(i, newdirstring)
        sim.paths[i].lastAcceptedFullTrajectory =[]
        for opi in range(len(sim.orderparameters.op)):
            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
    sim.log.log.info(dirstring + " exk: printing trajectories")  
    sim.outputAllFullTrajectories(newdirstring)
    return exchange


def exchangeKernelFast(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
#    sim.log.log.info(dirstring + " exk: starting trialReplicaExchange ")  
    exchange = sim.trialReplicaExchangeFast(dirstring,newdirstring)
#    sim.log.log.info(dirstring + " exk: reevaluating trajectories ")  
#    for i in sim.kernels.kernelPaths:
#        sim.getFullTrajectory(i, newdirstring)
#        sim.paths[i].lastAcceptedFullTrajectory =[]
#        for opi in range(len(sim.orderparameters.op)):
#            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
#    sim.log.log.info(dirstring + " exk: printing trajectories")  
#    sim.outputAllFullTrajectories(newdirstring)
    return exchange


def bfexchangeKernel(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
    
    sim.log.log.info(dirstring + " bfk: starting trialReplicaExchange ")  
    exchange = sim.trialBackwardForwardExchange(dirstring, newdirstring)
    
    sim.log.log.info(dirstring + " bfk: reevaluating trajectories ")  
    for i in sim.kernels.kernelPaths:
        sim.getFullTrajectory(i, newdirstring)
        sim.paths[i].lastAcceptedFullTrajectory =[]
        for opi in range(len(sim.orderparameters.op)):
            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
    sim.log.log.info(dirstring + " bfk: printing trajectories")  
    sim.outputAllFullTrajectories(newdirstring)
    return exchange


def bfexchangeKernelFast(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
    
#    sim.log.log.info(dirstring + " bfk: starting trialReplicaExchange ")  
    exchange = sim.trialBackwardForwardExchangeFast(dirstring, newdirstring)
    
#    sim.log.log.info(dirstring + " bfk: reevaluating trajectories ")  
#    for i in sim.kernels.kernelPaths:
#        sim.getFullTrajectory(i, newdirstring)
#        sim.paths[i].lastAcceptedFullTrajectory =[]
#        for opi in range(len(sim.orderparameters.op)):
#            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
#    sim.log.log.info(dirstring + " bfk: printing trajectories")  
#    sim.outputAllFullTrajectories(newdirstring)
    return exchange
 

if __name__ == '__main__':

    basedir = os.getcwd()
    dirnumber = 0
    rdirnumber = 0
    runs = 1
    nrshoots = 1
    mode = "tis"
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "basedir":
            basedir = argument[1]
        elif argument[0] == "mode":
            mode = argument[1]
        elif argument[0] == "startwith":
            dirnumber = int(float(argument[1]))
        elif argument[0] == "rstartwith":
            rdirnumber = int(float(argument[1]))
        elif argument[0] == "runs":
            runs = int(float(argument[1]))
        elif argument[0] == "nrshoots":
            nrshoots = int(float(argument[1])) 
        else:
            print "Arguments should be one of those: basedir, mode, startwith, rstartwith, runs"
            
    if basedir!="":
        if mode == 'initial': 
            sim = gtis.gromacstis(basedir,"initial",kernel="head")
            for i in range(sim.kernels.nkernels):   	#JR orig
#            for i in range(1):   			#JR for initializing reverse path
                pythonpath= sim.paths[0].options.paths["pythonpath"]
                mtispath= sim.paths[0].options.paths["tispath"]
                mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
        	queuename = sim.paths[0].options.runoptions["queuename"]
                mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
		mwrappath = sim.paths[0].options.paths["wrapperpath"]
    
                filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"initial",0,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)	#JR original
               # filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"initial",0,reverse=True,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)	#JR for initializing reverse path

                if mqsubs == "LL":
                    os.system("llsubmit " + " " + filename)
		elif mqsubs == "local":
		    os.system("chmod u+x " + " " + filename)
		    pe = "std.head.%03d.err\n" % i
                    po = "std.head.%03d.out\n" % i
		    os.system(filename + " >" + po + " " + " 2>" + pe)
                elif mqsubs == "local-slurm":
		    os.system("chmod u+x " + " " + filename)
		    pe = "std.head.%03d.err\n" % i
                    po = "std.head.%03d.out\n" % i
		    os.system(filename + " >" + po + " " + " 2>" + pe)
		elif mqsubs == "slurm" or mqsubs == "slurm-parallel":
                    os.system("sbatch " + " " + filename)
                else:
                    os.system("qsub " + " " + filename)
                    
        if mode == 'tis':
            
            dirstring = "%07d" % (dirnumber,)
            newdirstring = "%07d" % (dirnumber+1,)
            
            rdirstring = "%07d" % (rdirnumber,)
            rnewdirstring = "%07d" % (rdirnumber+1,)
            
            os.chdir(basedir)
            sim = gtis.gromacstis(basedir,"tis",kernel="head")
            
#JR start: comment out initial writing of trajectory files
            for i in sim.kernels.kernelPaths:
                        sim.getFullTrajectory(i, dirstring)
                        sim.paths[i].lastAcceptedFullTrajectory = sim.paths[i].fullTrajectory[:]
            sim.outputAllFullTrajectories(dirstring)
#JR end
#	    sys.exit(1)
            mtispath= sim.paths[0].options.paths["tispath"]
            os.chdir(basedir)
            shootingKernel(dirnumber,mtispath)  #JR commented out for testing
	   
            
            #shootingKernelReverse(rdirnumber,mtispath)		#JR stable state shooting, was commented out
#	    sys.exit(1)					#JR just stop here
            
            while True:
                time.sleep(20)
                if allFinished(newdirstring):
                    break
#JR comment out reverse shooting
              #  if allFinishedReverse(rdirstring):
              #      rdirnumber += 1
              #      rdirstring = "%07d" % (rdirnumber,)
              #      rnewdirstring = "%07d" % (rdirnumber+1,)
              #      shootingKernelReverse(rdirnumber,mtispath)
# JR end comment out reverse shooting
#                    
            dirnumber += 1
            
            for run in range(runs):
                os.chdir(basedir)
                dirstring = "%07d" % (dirnumber,)
                newdirstring = "%07d" % (dirnumber+1,)
		for i in range(nrshoots):                     	#DS for printing the trajectories only n-th time
                	os.chdir(basedir)
              		rn = random.randint(0,11)		#JR original, choose between different MC moves
#                	rn = 0					#JR just set this to 0 for testing
#                	if rn == 0:
                	if rn < 5:
		    		mtispath= sim.paths[0].options.paths["tispath"]		#JR set path to python executable
                    		#sim.log.log.info(dirstring + " sh: submitting shoot")
            	    		shootingKernel(dirnumber,mtispath)			#JR set path to python executable
#                    		shootingKernel(dirnumber)		#JR original
                    		ki = 0
                    		while True:
                    			#sim.log.log.info(dirstring + " sh: waiting " + str(ki))
                        		time.sleep(20)
                        		ki += 1
                        		if allFinished(newdirstring):
                            			break
# JR comment out reverse shooting
                       	#	if allFinishedReverse(rdirstring):
                        #    			rdirnumber += 1
                        #    			rdirstring = "%07d" % (rdirnumber,)
                        #    			rnewdirstring = "%07d" % (rdirnumber+1,)
                        #    			if rdirnumber < 300:
                        #        			shootingKernelReverse(rdirnumber,mtispath)
# JR end comment out reverse shooting
                    
                    		sim.log.log.info(dirstring + " sh")
                	elif rn < 11:
                    		#sim.log.log.info(dirstring + " ex: submitting ")  
#                    		exchange = exchangeKernel(dirnumber)
                    		exchange = exchangeKernelFast(dirnumber)
                    		sim.log.log.info(dirstring + " ex " + str(exchange))  
                	elif rn == 11:
                    		#sim.log.log.info(dirstring + " bf: submitting ")  
#                    		exchange = bfexchfangeKernel(dirnumber)
                    		exchange = bfexchangeKernelFast(dirnumber)
                    		sim.log.log.info(dirstring + " bf " + str(exchange)) 
			if i < (nrshoots-1):   
                        	allMoveDir(dirstring,newdirstring)            	#DS create the paths to the old and new directories
	
                dirnumber += 1							#DS the number of diretory increases only after going out from the loop

        if mode=="test":
            sim = gtis.gromacstis(basedir,"initial",kernel="head")
            
            
            
