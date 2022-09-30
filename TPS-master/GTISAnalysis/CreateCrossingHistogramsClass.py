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

class crossingHistoData(object):
    """
    Main TIS class which stores the paths, interfaces, stable states, and  paths. The class 
    also uses the gromacswrapper class, the filesystem class and helper.py. The logger is 
    self.log.
    """
    def __init__(self,histosize,histomax,histomin,start,stop,basedir=".",mode="initial"):
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

        #new variables
        self.HISTOSIZE = histosize
        self.HISTOMAX = histomax
        self.HISTOMIN = histomin
        self.start = start
        self.stop = stop

        for i in range(self.interfaces.ninterfaces[0]):
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=True))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=True))
            self.boundaryhistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=True))      #DS init the boudaryhistogram AB
            
            self.npaths += 1
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.crossinghistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=False))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=False))
            self.boundaryhistos.append(crossinghistogram.generichistogram(self.HISTOSIZE,self.HISTOMIN,self.HISTOMAX,forward=False))     #DS init the boudaryhistogram BA
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
        
        
    def run_main(self):
    
        for dirnumber in range(self.start,self.stop):
            dirstring = "%07d" % dirnumber
            for i in self.kernels.kernelPaths:
                filename= os.path.join(self.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
                self.readFullTrajectory(filename, i)
#               self.crossinghistos[i/2].addRangeToHisto(self.paths[i].getMaxTrajectoryValue())
                self.crossinghistos[i/2].addRangeToHistoFromInterface(self.paths[i].getMaxTrajectoryValue(),self.paths[i].interface)
                above = self.paths[i].getPointsBeyondInterface()   
                for point in above:
                    self.freeenergyhistos[i/2].addPointToHisto(point)
#--------------------------------------------------------------------------------------------------------#DS 
  
                ab = self.paths[i].checkAcceptedTISAB(self.orderparameters.op)
                if ab:
                    if self.paths[i].forward == True:
                        if self.paths[i+2].forward == True :                
                            interfacef = self.paths[i].getPointsBeyondInterfaceforward(self.paths[i+2].interface)
                            for point in interfacef:  
                                self.boundaryhistos[i/2].addPointToHisto(point)             
                        else:
                            above = self.paths[i].getPointsBeyondInterface()   
                            for point in above:
                                self.boundaryhistos[i/2].addPointToHisto(point)

                    if self.paths[i].forward == False:
                        if self.paths[i-2].forward == False :                
                            interfaceb = self.paths[i].getPointsBeyondInterfacebackward(self.paths[i-2].interface)
                            for point in interfaceb:  
                                self.boundaryhistos[i/2].addPointToHisto(point)             
                        else:
                            above = self.paths[i].getPointsBeyondInterface()   
                            for point in above:
                                self.boundaryhistos[i/2].addPointToHisto(point) 

#---------------------------------------------------------------------------------------------------------#DS
           
        for i in self.kernels.kernelPaths:
            filename = os.path.join(self.paths[i].nfsladir,"crossinghisto."+str(i)+".txt")
            self.crossinghistos[i/2].outputCrossingHisto(filename)
        
            filename = os.path.join(self.paths[i].nfsladir,"freeenergyhisto."+str(i)+".txt")
            self.freeenergyhistos[i/2].outputCrossingHisto(filename)

#---------------------------------------------------------------------------------------------------------#DS
            filename = os.path.join(self.paths[i].nfsladir,"boundaryhisto."+str(i)+".txt")     
            self.boundaryhistos[i/2].outputCrossingHisto(filename)
                  
