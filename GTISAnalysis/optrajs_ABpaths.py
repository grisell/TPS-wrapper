'''
Created on Jan 30, 2015

@author: grisell
'''
import pygromacstps.wrappers as wrappers
import pygromacstps.pathdata as pathdata
import pygromacstps.pathdatareverse as pathdatareverse
import pygromacstps.filesystem as filesystem
import pygromacstps.gtpslogging as gtpslogging
import pygromacstps.interfaces as interfaces
import pygromacstps.helpers_lammps as helpers
import pygromacstps.kernels as kernels
import pygromacstps.qsubsystem as qsubsystem
import pygromacstps.crossinghistogram as crossinghistogram
import pygromacstps.orderparameters as orderparameters

import multiprocessing
import os
import sys
import time
from math import log,fabs, exp
import subprocess as sub


TRAJECTORYFILE="trajectory.00."
SPTRAJECTORYFILE="corens3traj."
#TRAJECTORYFILE = "newnstrajectory."
#SSTRAJFILE="rnsnewtrajectory."


class rpeData(object):
    """
    Main TIS class which stores the paths, interfaces, stable states, and  paths. The class 
    also uses the gromacswrapper class, the filesystem class and helper.py. The logger is 
    self.log.
    """
    def __init__(self,basedir=".",mode="initial"):
        self.basedir = basedir
        self.mode = mode
        self.cores = multiprocessing.cpu_count()
        self.wrapper = wrappers.lammpswrapper()
        self.filesystem = filesystem.filesystem()
        self.interfaces = interfaces.interfaces()
        self.helper = helpers.helpers()
        self.kernels = kernels.kernels("head")
        self.qsubsystem = qsubsystem.qsubsystem()
        self.orderparameters = orderparameters.orderparameters()                   
        
        
        self.log = gtpslogging.log("info",basedir,"analysis")
        self.interfaces.readInterfaces(os.path.join(basedir,"options","interfaces.txt"))
        self.kernels.readKernelOptions(os.path.join(basedir,"options","kerneloptions.txt")) 

        #read the stables states from a file
        self.orderparameters.readOP(os.path.join(basedir,"options","orderparameters.txt"))  
        self.log.log.debug("Read OP : " + str(self.orderparameters.op[0]))
      
        self.paths = []
        self.newpaths = []
        self.npaths = 0
        #append histograms for forward ensemble
        for i in range(self.interfaces.ninterfaces[0]):
            
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.npaths += 1
        #append histograms for backward ensemble
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.newpaths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.newpaths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.npaths += 1
        
        #append histograms for stable state paths (reverse paths)
        self.reversePaths = []
        self.newreversePaths = []

        if len(self.interfaces.reversepaths[0]):
            for rp in self.interfaces.reversepaths[0]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
                self.newreversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
        if len(self.interfaces.reversepaths[1]):
            for rp in self.interfaces.reversepaths[1]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
                self.newreversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
        
        #List of interfaces and (directory) paths 
        self.kernels.generateKernelLists(self.npaths)
        
        #number of interfaces forward and backward 
        self.fninter = self.interfaces.ninterfaces[0]
        self.bninter = self.interfaces.ninterfaces[1]
        

    def readFullTrajectory(self,filename,pathnumber): #read full trajectories for forward and backward ensembles
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
        

    def readreverseFullTrajectory(self,filename,pathnumber): #read full trajectory for stable state trajectories
        traj = []
        self.reversePaths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.reversePaths[pathnumber].fullTrajectoryblength = 0
        self.reversePaths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.reversePaths[pathnumber].fullTrajectoryflength = 0
        self.reversePaths[pathnumber].lastAcceptedFullTrajectory = []
        self.reversePaths[pathnumber].fullTrajectory = []
        #self.reversePaths[pathnumber+1].lastAcceptedFullTrajectory = []
        #self.reversePaths[pathnumber+1].fullTrajectory = []
        
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            self.reversePaths[pathnumber].lastAcceptedFullTrajectoryflength += 1
            self.reversePaths[pathnumber].fullTrajectoryflength += 1
        
        self.reversePaths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        #self.reversePaths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.reversePaths[pathnumber].fullTrajectory.append(traj[:])
        #self.reversePaths[pathnumber+1].fullTrajectory.append(traj[:])
    

#***********************************************create reweighted path ensemble**********************************************************************************
if __name__=="__main__":
    basedir = os.path.join(os.getcwd(),"..")
    cdata = rpeData(basedir,"tis")   
    rstart=0
    rstop=0


#**Number of trajectories to be read for the interfaces and the stable states**
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "start": #starting trajectory to be read from each interface
            start = int(float(argument[1])) 
        elif argument[0] == "stop": #end trajectory to be read from each interface
            stop = int(float(argument[1]))
        elif argument[0] == "rstart": #starting trajectory tobe read from the stable states
            rstart = int(float(argument[1]))
        elif argument[0] == "rstop": #end trajectory to be read from the stable states
            rstop = int(float(argument[1]))
        else:
            print "Arguments required should be one of those: start, stop, rstart, rstop"
                                                 
    if not os.path.exists("AB_trajs"):
        proc = sub.call(["mkdir","AB_trajs"])
    if not os.path.exists("BA_trajs"):
        proc2 = sub.call(["mkdir","BA_trajs"])
    

    for dirnumber in range(start,stop): 
        dirstring = "%07d" % dirnumber
    
        for i in cdata.kernels.kernelPaths:   
            #read trajectories             
            filename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,TRAJECTORYFILE+dirstring+".dat")
            spfilename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,SPTRAJECTORYFILE+dirstring+".dat")
            cdata.readFullTrajectory(filename, i) 
                         
            typepath = cdata.paths[i].checkAcceptedTISpathtype(cdata.orderparameters.op) #assing type to path: AA, AB, BA, BB
            

            if typepath[1]: # AB paths
                if os.path.exists(spfilename): 
                    proc3 = sub.call(["cp",spfilename,"./AB_trajs/path_%d.%07d"%(i/2, dirnumber)])
                                            
            if typepath[3]: #BA 
                if os.path.exists(spfilename): 
                    proc4 = sub.call(["cp",spfilename,"./BA_trajs/path_%d.%07d"%(i/2, dirnumber)])
                
                #of = open(wfilename,"w")
                #of.write("%.30f\n" % (weight))
                #of.close()

          

            
            #for point in cdata.newpaths[i].lastAcceptedFullTrajectory[0]:    
               
                #x=float(point[1])
                #index=float(point[2])
    
                #if (index-oldindex)!= 0:
                  #  pass
                #else: 
                  
                                              
                #oldindex= index
                
                #print filename
            
                

    

   
