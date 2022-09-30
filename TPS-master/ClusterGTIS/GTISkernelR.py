'''
Created on July 27, 2010

@author: Wolfgang Lechner 
'''

import sys
import logging
import os

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


if __name__ == '__main__':
    
    basedir = ""
    dirnumber = 0
    mode = "tis"
    kernel = 0
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "basedir":
            basedir = argument[1]
        elif argument[0] == "mode":
            mode = argument[1]
        elif argument[0] == "dirnumber":
            dirnumber = int(float(argument[1]))
#	elif argument[0] == "kernel":
#            kernel = int(float(argument[1]))
        else:
            print "argument should be one of those: basedir,mode, startwith,rstartwith,runs"
#    mode="initial"		#JR artifically reset mode for testing 
    if basedir !="":
        if mode == 'initial':
            dirstring = "%07d" % (dirnumber,)
#            newdirstring = "%07d" % (dirnumber + 1,)                #DS the initial last accepted directory must be also 0000000
            newdirstring = "%07d" % (dirnumber ,)
#            sim = gtis.gromacstis(basedir,"tis",kernel="reverse")   #JR original
            sim = gtis.gromacstis(basedir,mode,kernel="reverse")
            k=0
            while(k<10):
		print 'reverse initial: ',k
            	sim.preperationReverse()
            	sim.shootingQueueReverse()
		sim.checkTisPathsAcceptedReverse()
            	sim.finalizeReverse(dirstring,newdirstring)

		successfull = []
		for i in range(len(sim.reversePaths)):
			successfull.append(sim.reversePaths[i].tisaccepted)
		if all(successfull):
			break
		k+=1

            	for i in range(len(sim.reversePaths)):
                	sim.getFullTrajectoryReverse(i, newdirstring)
                	sim.reversePaths[i].lastAcceptedFullTrajectory = sim.reversePaths[i].fullTrajectory[:]
            	sim.outputAllFullTrajectoriesReverse(newdirstring)
            
#            	sim.deleteScratchFiles()

            
        elif mode == 'tis':
            sim = gtis.gromacstis(basedir,"tis",kernel="reverse")
            dirstring = "%07d" % (dirnumber,)
            newdirstring = "%07d" % (dirnumber + 1,)
            
            sim.preperationReverse()
#            sim.lastReverseToGro(dirstring)
            sim.lastReverseToConf(dirstring)
            
            for i in range(len(sim.reversePaths)):
                sim.getFullTrajectoryReverse(i, dirstring)
                sim.reversePaths[i].lastAcceptedFullTrajectory = sim.reversePaths[i].fullTrajectory[:]
            sim.outputAllFullTrajectoriesReverse(dirstring)
            
#            sim.copyStateABGroFiles(dirstring)
            sim.copyStateABFiles(dirstring)
            sim.pickConfigurationReverse(dirstring)
            sim.shootingQueueReverse()
            sim.finalizeReverse(dirstring,newdirstring)
            sim.writeFinishedFilesReverse(dirstring)
            sim.deleteScratchFiles()
            
            


            
    
    
