'''
Created on July 27, 2010

@author: Wolfgang Lechner 
'''
import pygromacstps.tis as gtis
import sys

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
            sim.copyStateABConfFiles(dirstring)
            sim.pickConfigurationReverse(dirstring)
            sim.shootingQueueReverse()
            sim.finalizeReverse(dirstring,newdirstring)
            sim.writeFinishedFilesReverse(dirstring)
            sim.deleteScratchFiles()
            
            


            
    
    
