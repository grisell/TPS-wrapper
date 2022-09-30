'''
Created on Aug 4, 2010

@author: wolf
'''
import pygromacstps.wrappers as wrappers
import pygromacstps.parser as parser
import pygromacstps.pathdata as pathdata
import pygromacstps.filesystem as filesystem
import pygromacstps.gtpslogging as gtpslogging
import pygromacstps.interfaces as interfaces
import pygromacstps.helpers as helpers
import pygromacstps.kernels as kernels
import pygromacstps.qsubsystem as qsubsystem
import pygromacstps.crossinghistogram as crossinghistogram

import multiprocessing
import os
import sys

class crossingHistoData(object):
    """
    Main TIS class which stores the paths, interfaces, stable states, and  paths. The class 
    also uses the gromacswrapper class, the filesystem class and helper.py. The logger is 
    self.log.
    """
    def __init__(self,basedir=".",mode="initial"):
        self.basedir = basedir
        self.mode = mode
        self.cores = multiprocessing.cpu_count()
        self.wrapper = wrappers.gromacswrapper()
        self.distparser = parser.gdistparser()
        self.filesystem = filesystem.filesystem()
        self.interfaces = interfaces.interfaces()
        self.helper = helpers.helpers()
        self.kernels = kernels.kernels("head")
        self.qsubsystem = qsubsystem.qsubsystem()
        
        
        self.log = gtpslogging.log("info",basedir,"analysis")
        self.interfaces.readInterfaces(os.path.join(basedir,"options","interfaces.txt"))
        self.kernels.readKernelOptions(os.path.join(basedir,"options","kerneloptions.txt"))
        
        
        self.paths = []
        self.crossinghistos = []
        self.freeenergyhistos = []
        self.npaths = 0
        for i in range(self.interfaces.ninterfaces[0]):
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(forward=True))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(forward=True))
            
            self.npaths += 1
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(forward=False))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(forward=False))
            self.npaths += 1
        
        self.kernels.generateKernelLists(self.npaths)

    def readFullTrajectory(self,filename,pathnumber):
        traj = []
        self.paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectory = []
        self.paths[pathnumber+1].lastAcceptedFullTrajectory = []
        
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            if int(float(raw[2])) == 0:
                self.paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
            else:
                self.paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
        
        self.paths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.paths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        
        
if __name__=="__main__":
    basedir = os.path.join(os.getcwd(),"..")
    cdata = crossingHistoData(basedir,"tis")
    
    if len(sys.argv) == 3:
        start = int(float(sys.argv[1]))
        stop =  int(float(sys.argv[2]))
    else:
        print "two arguments required: start and stop"
    
    
    for dirnumber in range(start,stop):
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
            cdata.crossinghistos[i/2].addRangeToHisto(cdata.paths[i].getMaxTrajectoryValue())
            above = cdata.paths[i].getPointsBeyondInterface()
            for point in above:
                cdata.freeenergyhistos[i/2].addPointToHisto(point)
    
    for i in cdata.kernels.kernelPaths:
        filename = os.path.join(cdata.paths[i].nfsladir,"crossinghisto."+str(i)+".txt")
        cdata.crossinghistos[i/2].outputCrossingHisto(filename)
        
        filename = os.path.join(cdata.paths[i].nfsladir,"freeenergyhisto."+str(i)+".txt")
        cdata.freeenergyhistos[i/2].outputCrossingHisto(filename)
        
    
    
    