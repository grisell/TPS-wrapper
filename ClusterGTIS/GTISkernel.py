'''
Created on July 27, 2010

@author: Wolfgang Lechner 
'''
import pygromacstps.tis as gtis
import sys
from random import shuffle,randint
import os


def replaceFirstInterface():
    for i in sim.kernels.kernelPaths:
        for j in range(len(sim.reversePaths)):
            if sim.paths[i].interface == sim.reversePaths[j].interface:
                filelist = sim.filesystem.getFileList(sim.reversePaths[j].nfsconfigstore, "path*")
                if len(filelist) > 0:
                    sim.log.log.info("swapping with stable state") 
                    shuffle(filelist)
                    sim.filesystem.copyFiles(filelist[0], os.path.join(sim.paths[i].workdir,"conf.dump"))
                    sim.filesystem.copyFiles(filelist[0], os.path.join(sim.paths[i+1].workdir,"conf.dump"))    #DS copy the chosen slice to the backward part too
                    sim.reverseBackwardLammpsFile(i+1)                                                         #DS reverse velocities for backward shooting





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
        elif argument[0] == "kernel":
            kernel = int(float(argument[1]))
        else:
            print "argument should be one of those: basedir,mode, startwith,rstartwith,runs"
    
    if basedir!="":
        if mode == 'initial':
            os.chdir(basedir)
            sim = gtis.gromacstis(basedir,"initial",kernel)
            
            successfull = [False for x in sim.kernels.kernelPaths]

            k = 0
            while(True):   #JR: orig -> run as long as initial path for all interfaces found!
#            while(k<1):    #JR
                os.chdir(basedir)
                sim.preperationFromStart()
#                sim.preperationReverse() 	#JR: comment in for dynamic shooting
                sim.shootingInitialLammpsFiles()
                
                if sim.kernels.dynamicshooting == 1:
                    sim.shootingQueueDynamic()
                else:
                    sim.shootingQueue()
                
                sim.checkAllTisPathsAccepted()
                sim.finalizeInitial()
                successfull = []
                for i in sim.kernels.kernelPaths:
                    successfull.append(sim.paths[i].tisaccepted)
                if all(successfull):
                    break
                k+=1
                sim.outputAllFullTrajectories("%07d" % (0,))
                sim.deleteScratchFiles()    #JR
        
            
        elif mode == 'tis':
            dirstring = "%07d" % (dirnumber,)
            newdirstring = "%07d" % (dirnumber + 1,)
            
            sim = gtis.gromacstis(basedir,"tis",kernel)
            print sim.kernels.kernelPathsAll
            sim.preperationTIS()
#            sim.preperationReverse()	#JR: comment in for dynamic shooting
            sim.readLastAcceptedTrajectories(dirstring)
            sim.lastAcceptedToConf(dirstring)
            sim.pickConfigurationsTIS(dirstring)
#	    JR: swap slice with stable state!
            #if randint(0,2) == 0:
            #    if range(len(sim.reversePaths)) > 0:
            #        replaceFirstInterface()
#	    JR: end comment out
#            sim.log.log.info(dirstring + " shk: do shooting")
            if sim.kernels.dynamicshooting == 1:
                sim.shootingQueueDynamic()
            else:
                sim.shootingQueue()
            
#            sim.log.log.info(dirstring + " shk: checking acceptance")
            sim.checkAllTisPathsAccepted()
            sim.finalizeTIS(dirstring,newdirstring)
#            sim.log.log.info(dirstring + " shk: printing trajectories")
            sim.outputAllFullTrajectories(newdirstring)
#            sim.deleteOldTrrFiles(dirstring)		#JR: check if there are also files to be deleted from lammps
            sim.deleteScratchFiles()
            sim.writeFinishedFiles(newdirstring)
        elif mode=="test":
            sim = gtis.gromacstis(basedir,"tis",kernel)
            print sim.kernels.kernelPathsAll
            
    
