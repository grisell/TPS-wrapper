'''
Created on Aug 4, 2010

@author: wolf
'''

import multiprocessing
import os
import sys
import time
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.analysis.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

#SRM:import the module for selecting lammps/gromacs
import pygromacstps.findmd as md

#SRM:find the module path
module_path,custom_path = md.find_paths()

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

import wrappers as wrappers
#import pygromacstps.parser as parser
import pygromacstps.pathdata as pathdata
import pygromacstps.filesystem as filesystem
import pygromacstps.interfaces as interfaces
import helpers as helpers
import pygromacstps.kernels as kernels
import pygromacstps.qsubsystem as qsubsystem
import pygromacstps.crossinghistogram as crossinghistogram
import pygromacstps.orderparameters as orderparameters


HISTOSIZE=100
HISTOMAX=500.0
HISTOMIN=0.0

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
#        self.wrapper = wrappers.gromacswrapper()
        self.wrapper = wrappers.wrapper()
#        self.distparser = parser.gdistparser()
        self.filesystem = filesystem.filesystem()
        self.interfaces = interfaces.interfaces()
        self.helper = helpers.helpers()
        self.kernels = kernels.kernels("head")
        self.qsubsystem = qsubsystem.qsubsystem()
        self.orderparameters = orderparameters.orderparameters()                   #DS
        self.interfaces.readInterfaces(os.path.join(basedir,"options","interfaces.txt"))
        logger.info("Interfaces read")
        self.kernels.readKernelOptions(os.path.join(basedir,"options","kerneloptions.txt")) 

        #read the stables states from a file
        self.orderparameters.readOP(os.path.join(basedir,"options","orderparameters.txt"))  #DS
        logger.info("Read OP : " + str(self.orderparameters.op[0]))
      
        self.paths = []
        self.crossinghistos = []
        self.freeenergyhistos = []
        self.boundaryhistos = []      #DS
        self.npaths = 0
        for i in range(self.interfaces.ninterfaces[0]):
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=True))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=True))
            self.boundaryhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=True))      #DS init the boudaryhistogram AB
            
            self.npaths += 1
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=False))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=False))
            self.boundaryhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=False))     #DS init the boudaryhistogram BA
            self.npaths += 1
        
        self.kernels.generateKernelLists(self.npaths)

    def readFullTrajectory(self,filename,pathnumber):
        traj = []
        self.paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.paths[pathnumber].fullTrajectoryblength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.paths[pathnumber].fullTrajectoryflength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectory = []
        self.paths[pathnumber].fullTrajectory = []
        self.paths[pathnumber+1].lastAcceptedFullTrajectory = []
        self.paths[pathnumber+1].fullTrajectory = []
        
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            if int(float(raw[2])) == 0:
                self.paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                self.paths[pathnumber].fullTrajectoryblength += 1
            else:
                self.paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
                self.paths[pathnumber].fullTrajectoryflength += 1
        
        self.paths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.paths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.paths[pathnumber].fullTrajectory.append(traj[:])
        self.paths[pathnumber+1].fullTrajectory.append(traj[:])
        
        
if __name__=="__main__":
    basedir = os.path.join(os.getcwd())
    cdata = crossingHistoData(basedir,"tis")
    
    if len(sys.argv) == 3:
        try:
            int(float(sys.argv[1]))
        except ValueError:
            logger.error("start should be a number. User entered %s. Program exiting",argument[1])
            raise SystemExit()
        start = int(float(sys.argv[1]))
        try:
            int(float(sys.argv[2]))
        except ValueError:
            logger.error("stop should be a number. User entered %s. Program exiting",argument[1])
            raise SystemExit()
        stop =  int(float(sys.argv[2]))
    else:
        print "two arguments required: start and stop"
    
    
    for dirnumber in range(start,stop):
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
#            cdata.crossinghistos[i/2].addRangeToHisto(cdata.paths[i].getMaxTrajectoryValue())
            cdata.crossinghistos[i/2].addRangeToHistoFromInterface(cdata.paths[i].getMaxTrajectoryValue(),cdata.paths[i].interface)
            above = cdata.paths[i].getPointsBeyondInterface()   
            for point in above:
                cdata.freeenergyhistos[i/2].addPointToHisto(point)
#--------------------------------------------------------------------------------------------------------#DS 
  
            ab = cdata.paths[i].checkAcceptedTISAB(cdata.orderparameters.op)
            if ab:
                if cdata.paths[i].forward == True:
                    if cdata.paths[i+2].forward == True :                
                        interfacef = cdata.paths[i].getPointsBeyondInterfaceforward(cdata.paths[i+2].interface)
                        for point in interfacef:  
                            cdata.boundaryhistos[i/2].addPointToHisto(point)             
                    else:
                        above = cdata.paths[i].getPointsBeyondInterface()   
                        for point in above:
                            cdata.boundaryhistos[i/2].addPointToHisto(point)

                if cdata.paths[i].forward == False:
                    if cdata.paths[i-2].forward == False :                
                        interfaceb = cdata.paths[i].getPointsBeyondInterfacebackward(cdata.paths[i-2].interface)
                        for point in interfaceb:  
                            cdata.boundaryhistos[i/2].addPointToHisto(point)             
                    else:
                        above = cdata.paths[i].getPointsBeyondInterface()   
                        for point in above:
                            cdata.boundaryhistos[i/2].addPointToHisto(point) 

#---------------------------------------------------------------------------------------------------------#DS
           
    for i in cdata.kernels.kernelPaths:
        filename = os.path.join(cdata.paths[i].nfsladir,"crossinghisto."+str(i)+".txt")
        cdata.crossinghistos[i/2].outputCrossingHisto(filename)
        
        filename = os.path.join(cdata.paths[i].nfsladir,"freeenergyhisto."+str(i)+".txt")
        cdata.freeenergyhistos[i/2].outputCrossingHisto(filename)

#---------------------------------------------------------------------------------------------------------#DS
        filename = os.path.join(cdata.paths[i].nfsladir,"boundaryhisto."+str(i)+".txt")     
        cdata.boundaryhistos[i/2].outputCrossingHisto(filename)
                  
