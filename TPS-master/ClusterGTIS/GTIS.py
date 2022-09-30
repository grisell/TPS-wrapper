'''
Created on May 3, 2010

@author: Wolfgang Lechner
@status: beta 
'''
import sys
import os
import random
import time
import shutil
import logging
import subprocess as sub


#SRM:import the module for selecting lammps/gromacs
import pygromacstps.findmd as md

#SRM:find the module path
module_path,custom_path = md.find_paths()

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False              


#SRM:check if module path is okay and add it to list of paths
if os.path.exists(module_path):
    sys.path.append(module_path)
else:
    logger.error("lammps/gromacs specific modules not found in the code directory")
    raise SystemExit()

#SRM:check if custom path is okay and add it to list of paths
if os.path.exists(custom_path):
    sys.path.append(custom_path)
else:
    logger.error("lammps/gromacs specific order parameter files not found in the code directory")
    raise SystemExit()


#SRM:import rest of the things
import pygromacstps.tis as gtis
import pygromacstps.pathstats as pstat

#rest of the code starts
def shootingKernel(dirnumber,mtispath=""):
    for i in range(sim.kernels.nkernels):
	pythonpath= sim.paths[0].options.paths["pythonpath"]
	mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
	mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
	queuename = sim.paths[0].options.runoptions["queuename"]
	mwrappath = sim.paths[0].options.paths["wrapperpath"]
	
	filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"tis",str(dirnumber),reverse=False,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)
	print mqsubs
	if mqsubs == "LL":
	    os.system("llsubmit " + filename)
	elif mqsubs == "local":
	    rcmd = []
	    os.system("chmod u+x " + " " + filename)
	    pe = "std.head.%03d.err\n" % i
	    po = "std.head.%03d.out\n" % i
	    #os.system(filename + " >" + po + " " + " 2>" + pe)
	    #srm: parallelize calculation in local
	    rcmd.append(filename)
	    rcmd.append(">")
	    rcmd.append(po)
	    rcmd.append("2>")
            rcmd.append(pe)
	    proc = sub.Popen(rcmd, stdin=None,stdout=None,stderr=None,close_fds=True)
	    proc.poll()
	    #1234
	elif mqsubs == "vulcan" or mqsubs=="vulcan-parallel":
	    os.system("qsub " + filename)
        elif mqsubs == "slurm" or mqsubs == "slurm-parallel" :
            os.system("sbatch " + filename)
	
def shootingKernelReverse(rdirnumber,mtispath=""):
    pythonpath= sim.paths[0].options.paths["pythonpath"]
    mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
    mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
    queuename = sim.paths[0].options.runoptions["queuename"]
    mwrappath = sim.paths[0].options.paths["wrapperpath"]
    
    filename = sim.qsubsystem.writeKernelQsubFile(basedir, 0,pythonpath,"tis",str(rdirnumber),reverse=True,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)
    print mqsubs
    if mqsubs == "LL":
	os.system("llsubmit " + filename)
    elif mqsubs == "local":
	os.system("chmod u+x " + " " + filename)
	pe = "std.head.%03d.err\n" % i
	po = "std.head.%03d.out\n" % i
	os.system(filename + " >" + po + " " + " 2>" + pe)
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
    exchange = sim.trialReplicaExchangeFast(dirstring,newdirstring)
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
    exchange = sim.trialBackwardForwardExchangeFast(dirstring, newdirstring)
    return exchange



if __name__ == '__main__':

    basedir = os.getcwd()
    dirnumber = 0
    rdirnumber = 0
    runs = 1
    nrshoots = 1
    mode = "tis"

    runningfile = os.path.join(basedir,'.isrunning')
    if os.path.exists(runningfile):
        logger.info("Seems to be that a code is running here or it did not exit properly. Please check!!")
    else:
        os.system("touch .isrunning")
    
    #SRM:to check the input options
    option_checker = [0,0,0,0,0]
    
    for arg in sys.argv[1:]:
	argument = arg.split("=")
	if argument[0] == "basedir":
	    basedir = argument[1]
	 
	elif argument[0] == "mode":
	    #SRM:check for mode
	    if ((argument[1] != "initial") and (argument[1] != "tis") and (argument[1] != "tps")):
		logger.error("mode should be tis/tps/initial. User entered %s. Program exiting",argument[1])
		raise SystemExit()
	    option_checker[0] = 1
	    mode = argument[1]
	elif argument[0] == "startwith":            
	    #SRM:check for dirnumber
	    try:
		float(argument[1])
	    except ValueError:
		logger.error("startwith should be a number. User entered %s. Program exiting",argument[1])
		raise SystemExit()
	    option_checker[1] = 1
	    dirnumber = int(float(argument[1]))
	elif argument[0] == "rstartwith":
	    #SRM:check for dirnumber
	    try:
		float(argument[1])
	    except ValueError:
		logger.error("rstartwith should be a number. User entered %s. Program exiting",argument[1])
		raise SystemExit()
	    option_checker[2] = 1
	    rdirnumber = int(float(argument[1]))           
	elif argument[0] == "runs":
	    #SRM:check for runs
	    try:
		float(argument[1])
	    except ValueError:
		logger.error("runs should be a number. User entered %s. Program exiting",argument[1])
		raise SystemExit()
	    option_checker[3] = 1
	    runs = int(float(argument[1]))
	elif argument[0] == "nrshoots":
	    #SRM:check for nrshoots
	    try:
		float(argument[1])
	    except ValueError:
		logger.error("nrshoots should be a number. User entered %s. Program exiting",argument[1])
		raise SystemExit()
	    option_checker[4] = 1
	    nrshoots = int(float(argument[1]))
	else:
	    logger.error("Arguments should be one of those: basedir, mode, startwith, rstartwith, runs")
	    raise SystemExit()

    #SRM:code block to print the values out
    #SRM:to be setup to a loop once the problem with the list is identified
    
    if (option_checker[0]==0):
	logger.info('Value set for mode = %s , default value',mode)
    else: 
	logger.info('Value set for mode = %s , user input value',mode)
    if (option_checker[1]==0):
	logger.info('Value set for startwith = %s , default value',dirnumber)
    else: 
	logger.info('Value set for startwith = %s , user input value',dirnumber)
    if (option_checker[2]==0):
	logger.info('Value set for rstartwith = %s , default value',rdirnumber)
    else: 
	logger.info('Value set for rstartwith = %s , user input value',rdirnumber)
    if (option_checker[3]==0):
	logger.info('Value set for runs = %s , default value',runs)
    else: 
	logger.info('Value set for runs = %s , user input value',runs)
    if (option_checker[4]==0):
	logger.info('Value set for nrshoots = %s , default value',nrshoots)
    else: 
	logger.info('Value set for nrshoots = %s , user input value',nrshoots)        
    

    if basedir!="":
	if mode == 'initial': 
	    #add one more argument and pass the full path of the mddriver modules
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
		#filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"initial",0,reverse=True,qsubs=mqsubs,walltime=mwalltime,queuename=queuename,tispath=mtispath,wrappath=mwrappath)	#JR for initializing reverse path

		if mqsubs == "LL":
		    os.system("llsubmit " + " " + filename)
		elif mqsubs == "local":
		    rcmd = []
		    os.system("chmod u+x " + " " + filename)
		    pe = "std.head.%03d.err\n" % i
		    po = "std.head.%03d.out\n" % i
		    rcmd.append(filename)
		    rcmd.append(">")
		    rcmd.append(po)
		    rcmd.append("2>")
		    rcmd.append(pe)
		    proc = sub.Popen(rcmd, stdin=None,stdout=None,stderr=None,close_fds=True)
		    proc.poll()
		    #os.system(filename + " >" + po + " " + " 2>" + pe)
		elif mqsubs == "vulcan" or mqsubs=="vulcan-parallel":
		    os.system("qsub " + " " + filename)
                #elif mqsubs == "slurm":
                    #os.system("sbatch " + " " + filename)
        	elif mqsubs == "slurm" or mqsubs == "slurm-parallel" :
            	    os.system("sbatch " + filename)
		    
	if mode == 'tis':
	    
	    dirstring = "%07d" % (dirnumber,)
	    newdirstring = "%07d" % (dirnumber+1,)
	    
	    rdirstring = "%07d" % (rdirnumber,)
	    rnewdirstring = "%07d" % (rdirnumber+1,)
	    
	    os.chdir(basedir)
	    sim = gtis.gromacstis(basedir,"tis",kernel="head")
            stats = pstat.pathstats()

	    #srm code to check if reversepaths are there
	    if ((len(sim.interfaces.reversepaths[0])!=0) or (len(sim.interfaces.reversepaths[1])!=0)):
		reverserun=True
	    else:
		reverserun=False

	    if (sim.interfaces.ninter>1):
		shoot = sim.paths[0].options.runoptions["shoot"]
		swap = sim.paths[0].options.runoptions['swap']
		bfswap = sim.paths[0].options.runoptions['bfswap']
		logger.info(("Percentage of shooting moves set : %d")% shoot)
		logger.info(("Percentage of swapping moves set : %d")% swap)
		logger.info(("Percentage of bf swapping moves set : %d")% bfswap)
		if ((shoot+swap+bfswap!=100)):
			logger.error("Sum of shoot,swap and bfswap moves is not hundred! Please check the values.")
			logger.error("Only integer values between 0 to 100 are accepted for shoot,swap and bfswap.")
			raise SystemExit()

	    else:
	    	shoot = 100
	    	swap = 0
	    	bfswap = 0
	    	logger.info("Only one interface is found. Switching to only shotting moves")
 
     
#JR start: comment out initial writing of trajectory files
	    #srm option to keep or remove initial evaluation of order parameter trajectory
	    if sim.paths[0].options.runoptions["initialopeval"] == 'True':
	    	for i in sim.kernels.kernelPaths:
			sim.getFullTrajectory(i, dirstring)
			sim.paths[i].lastAcceptedFullTrajectory = sim.paths[i].fullTrajectory[:]
	    	sim.outputAllFullTrajectories(dirstring)
#JR end
#	    sys.exit(1)
	    mtispath= sim.paths[0].options.paths["tispath"]
	    os.chdir(basedir)
	    shootingKernel(dirnumber,mtispath)  #JR commented out for testing
	   
	    if reverserun==True:
		shootingKernelReverse(rdirnumber,mtispath)		#JR stable state shooting, was commented out
#	    sys.exit(1)					#JR just stop here
	    
	    while True:
		time.sleep(20)
		if allFinished(newdirstring):
		    break
		if reverserun==True:
		    if allFinishedReverse(rdirstring):
			rdirnumber += 1
			rdirstring = "%07d" % (rdirnumber,)
			rnewdirstring = "%07d" % (rdirnumber+1,)
			shootingKernelReverse(rdirnumber,mtispath)
# JR end comment out reverse shooting
#                    
            statfrequency = sim.paths[0].options.runoptions["statfrequency"]
	    dirnumber += 1
	    nps = 0

	    for run in range(runs):
		os.chdir(basedir)
		dirstring = "%07d" % (dirnumber,)
		newdirstring = "%07d" % (dirnumber+1,)
		for i in range(nrshoots):                     	#DS for printing the trajectories only n-th time
			nps+=1
                        os.chdir(basedir)
 #             		rn = random.randint(0,11)		#JR original, choose between different MC moves
			rn = random.randint(0,100)
#                	rn = 0					#JR just set this to 0 for testing
#                	if rn == 0:
			if (rn <= shoot):
                                #print "I shot"
                                #print rn
                                #print shoot
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
				if reverserun==True:
				    if allFinishedReverse(rdirstring):
					rdirnumber += 1
					rdirstring = "%07d" % (rdirnumber,)
					rnewdirstring = "%07d" % (rdirnumber+1,)
					if rdirnumber < 300:
					    shootingKernelReverse(rdirnumber,mtispath)
# JR end comment out reverse shooting
		    
				sim.log.log.info(dirstring + " sh")
			#elif rn < 11:
			elif rn <= (shoot+swap):
                                #print "I swapped"
                                #print rn
                                #print shoot+swap
				#sim.log.log.info(dirstring + " ex: submitting ")  
				exchange = exchangeKernelFast(dirnumber)
				sim.log.log.info(dirstring + " ex " + str(exchange))  
			#elif rn == 11:
			elif rn <= (shoot+swap+bfswap):
                                #print "I bfswapped"
                                #print rn
                                #print shoot+swap+bfswap
				#sim.log.log.info(dirstring + " bf: submitting ")
				exchange = bfexchangeKernelFast(dirnumber)
				sim.log.log.info(dirstring + " bf " + str(exchange)) 
			if i < (nrshoots-1):   
				allMoveDir(dirstring,newdirstring)            	#DS create the paths to the old and new directories
	               
                        if (nps%statfrequency==0):
                                stats.PrintReport()
                                
		dirnumber += 1

            os.system("rm .isrunning")

	if mode=="test":
	    sim = gtis.gromacstis(basedir,"initial",kernel="head")
	    
	    
	    
