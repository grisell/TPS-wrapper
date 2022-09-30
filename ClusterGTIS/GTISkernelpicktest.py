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
                    shuffle(filelist)
                    sim.filesystem.copyFiles(filelist[0], os.path.join(sim.paths[i].workdir,"conf.gro"))
                    

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
   
    if basedir != "":
        if mode == 'initial':
            sim = gtis.gromacstis(basedir,"initial",kernel)
            successfull = [False for x in sim.kernels.kernelPaths]
            k = 0
            while(True):
                sim.preperationFromStart()
                sim.shootingInitialGroFiles()
                sim.shootingQueue()
                sim.checkAllTisPathsAccepted()
                sim.finalizeInitial()
                
                for i in range(len(sim.kernels.kernelPaths)):
                    if sim.paths[i].tisaccepted:
                        successfull[i] = True
                print k, successfull
                if all(successfull):
                    break
                k+=1
                sim.deleteScratchFiles()     
        elif mode == 'tis':

            dirstring = "%07d" % (dirnumber,)
            newdirstring = "%07d" % (dirnumber + 1,)
            
            sim = gtis.gromacstis(basedir,"tis",kernel)
            sim.preperationTIS()
            sim.readLastAcceptedTrajectories(dirstring)
            sim.lastAcceptedToGro(dirstring)
            sim.pickConfigurationsTIS(dirstring)
            
            if randint(0,10) == 0:
                if range(len(sim.reversePaths)) > 0:
                    replaceFirstInterface()
            
            sim.shootingQueue()
            sim.checkAllTisPathsAccepted()
            sim.finalizeTIS(dirstring,newdirstring)
            sim.outputAllFullTrajectories(newdirstring)
            sim.deleteOldTrrFiles(dirstring)
            sim.deleteScratchFiles()
            sim.writeFinishedFiles(newdirstring)
            
    