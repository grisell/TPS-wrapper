'''
Created on Jan 30, 2015

@author: grisell
'''
import pygromacstps.wrappers as wrappers
#import pygromacstps.parser as parser
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

HISTOSIZE=2000
CUTOFF=0.01
LARGENUMBER=1000000


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
#        self.wrapper = wrappers.gromacswrapper()
        self.wrapper = wrappers.lammpswrapper()
#        self.distparser = parser.gdistparser()
        self.filesystem = filesystem.filesystem()
        self.interfaces = interfaces.interfaces()
        self.helper = helpers.helpers()
        self.kernels = kernels.kernels("head")
        self.qsubsystem = qsubsystem.qsubsystem()
        self.orderparameters = orderparameters.orderparameters()                   #DS
        
        
        self.log = gtpslogging.log("info",basedir,"analysis")
        self.interfaces.readInterfaces(os.path.join(basedir,"options","interfaces.txt"))
        self.kernels.readKernelOptions(os.path.join(basedir,"options","kerneloptions.txt")) 

        #read the stables states from a file
        self.orderparameters.readOP(os.path.join(basedir,"options","orderparameters.txt"))  #DS
        self.log.log.debug("Read OP : " + str(self.orderparameters.op[0]))
      
        self.paths = []
        self.freeenergyhistos = []
        self.npaths = 0
        for i in range(self.interfaces.ninterfaces[0]):
            
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(histosize=2000,histomin=0.0,histomax=2048.0,forward=True))
            self.npaths += 1
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.freeenergyhistos.append(crossinghistogram.generichistogram(histosize=2000,histomin=0.0,histomax=2048.0,forward=False))

            self.npaths += 1
        
        self.reversePaths = []
        
        if len(self.interfaces.reversepaths[0]):
            for rp in self.interfaces.reversepaths[0]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
        if len(self.interfaces.reversepaths[1]):
            for rp in self.interfaces.reversepaths[1]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
        
        self.kernels.generateKernelLists(self.npaths)


        self.fninter = self.interfaces.ninterfaces[0]
        self.bninter = self.interfaces.ninterfaces[1]
        

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

    def readreverseFullTrajectory(self,filename,pathnumber):
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

    def readpathweights(self,ninter,forward):
        PAnorm= [0]*ninter
        weight= [0]*ninter
        filename = forward and os.path.join(os.getcwd(),"selfconsistent","fpathweights.txt") or os.path.join(os.getcwd(),"selfconsistent","bpathweights.txt")       
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[2]) #GDL:Get path weights 
            z = float(line.split()[3]) #GDL: Normalization factor 
            PAnorm[x]=z
            weight[x] = y
            c+= 1
        return (PAnorm, weight)

    def readcacb(self):
        filename = os.path.join(os.getcwd(),"selfconsistent","cacb.dat")       
        c = 0
        for line in open( filename , 'r' ):
            ca = float(line.split()[3]) #GDL:Get path weights 
            cb = float(line.split()[4]) #GDL: Normalization factor     
        return (ca, cb)

############################################################################################################################################################################
if __name__=="__main__":
    basedir = os.path.join(os.getcwd(),"..")
    cdata = rpeData(basedir,"tis")   
    rstart=0
    rstop=0
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "start":
            start = int(float(argument[1])) 
        elif argument[0] == "stop":
            stop = int(float(argument[1]))
        elif argument[0] == "rstart":
            rstart = int(float(argument[1]))
        elif argument[0] == "rstop":
            rstop = int(float(argument[1]))
        else:
            print "Arguments required should be one of those: start, stop, rstart, rstop"

    
    fweight = cdata.readpathweights(cdata.fninter,forward=True)
    bweight = cdata.readpathweights(cdata.bninter,forward=False)
    cvalue=cdata.readcacb() # Note: cA=cvalue[0], cB=cvalue[1]
   
    for dirnumber in range(start,stop):
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
            
            
            typepath = cdata.paths[i].checkAcceptedTISpathtype(cdata.orderparameters.op)
            wfilename=os.path.join(cdata.paths[i].nfsladir,dirstring,"weight.dat")
            
            if typepath[0]: #Assign weights to paths typepath[0]: AA paths
                
                if not os.path.exists(wfilename):
                    maxlambda=cdata.paths[i].getMaxTrajectoryValue()
                    nwi=-1  
                 
                    for k in reversed(range(0,cdata.fninter)):
                        intervalue = cdata.interfaces.interfaces[0][k]
                        if maxlambda > intervalue:
                            nwi=k
                            break
                    weight=cvalue[0]*fweight[1][nwi]
                    of = open(wfilename,"w")
                    of.write("%.18f\n" % (weight))
                    of.close()

                else:                       
                    for line in open( wfilename , 'r' ):
                        weight = float(line.split()[0]) 
                

            if typepath[1]: #Assign weights to paths typepath[1]: AB paths
                if not os.path.exists(wfilename):
                    nwi=cdata.fninter-1
                    weight=cvalue[0]*fweight[1][nwi]

                    of = open(wfilename,"w")
                    of.write("%.18f\n" % (weight))
                    of.close()
                else:
                    for line in open( wfilename , 'r' ):
                        weight = float(line.split()[0])

            if typepath[2]: #Assign weights to paths typepath[2]: BB paths
                if not os.path.exists(wfilename):
                    minlambda=cdata.paths[i].getMaxTrajectoryValue()
                    nwi=-1  
                    for k in range(0,cdata.bninter):
                        intervalue = cdata.interfaces.interfaces[1][k]
                        if minlambda < intervalue:
                            nwi=k
                            break
                    weight=cvalue[1]*bweight[1][nwi]

                    of = open(wfilename,"w")
                    of.write("%.18f\n" % (weight))
                    of.close()
                else:
                    for line in open( wfilename , 'r' ):
                        weight = float(line.split()[0])
                        
            if typepath[3]: #Assign weights to paths typepath[3]: BA paths
                if not os.path.exists(wfilename):
                    nwi=0
                    weight=cvalue[1]*bweight[1][nwi]
                
                    of = open(wfilename,"w")
                    of.write("%.18f\n" % (weight))
                    of.close()
                else:
                    for line in open( wfilename , 'r' ):
                        weight = float(line.split()[0])

            #Create weighted histograms from the trajectories
            for point in cdata.paths[i].lastAcceptedFullTrajectory[0]:
                x=float(point[1])
                cdata.freeenergyhistos[0].addweightedPointToHisto(x,weight)
                    

    #Stable states histograms:
    countA=0
    countB=0
    oldposA=LARGENUMBER
    oldposB=-LARGENUMBER
    numinterpaths=fabs(start-stop)
    
    #compute the actual number of complete stable state paths to re-scale the weight
    for dirnumber in range(rstart,rstop):
        dirstring = "%07d" % dirnumber
        reversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"trajectoryZ.00."+dirstring+".acc")
        cdata.readreverseFullTrajectory(reversefilename, 0)
        
        for point in cdata.reversePaths[0].lastAcceptedFullTrajectory[0]:
            if cdata.reversePaths[0].state == 0:
                newposA=float(point[1])
                if cdata.reversePaths[0].ifcrossstablestate(newposA,oldposA):
                    countA +=1
                oldposA=newposA
            else:
                newposB=float(point[1])
                if cdata.reversePaths[0].ifcrossstablestate(newposB,oldposB):
                    countB +=1
                oldposB=newposB

    #assign weights to the stable states and compute the weighted histograms of the order parameter:            
    for dirnumber in range(rstart,rstop):
        dirstring = "%07d" % dirnumber
        reversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"trajectoryZ.00."+dirstring+".acc")
        wreversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"weight.dat")
        cdata.readreverseFullTrajectory(reversefilename, 0)
        scalenumber=cdata.reversePaths[0].getscalenumberstablestates(countA,countB,numinterpaths)
        
        if cdata.reversePaths[0].state == 0:
            if not os.path.exists(wreversefilename):
                weight=cvalue[0]*fweight[1][0]
                of = open(wreversefilename,"w")
                of.write("%.18f\n" % (weight))
                of.close()
            else:
                for line in open(wreversefilename, 'r' ):
                    weight = float(line.split()[0])
        else:
            if not os.path.exists(wreversefilename):
                weight= cvalue[1]*bweight[1][cdata.bninter-1]
                of = open(wreversefilename,"w")
                of.write("%.18f\n" % (weight))
                of.close()
            else:
                for line in open(wreversefilename, 'r' ):
                    weight = float(line.split()[0])
        
        weight = scalenumber*weight

        for point in cdata.reversePaths[0].lastAcceptedFullTrajectory[0]:
            x=float(point[1])
            #if cdata.reversePaths[0].state == 0 and x < cdata.reversePaths[0].interface:
            cdata.freeenergyhistos[0].addweightedPointToHisto(x,weight)
            #elif cdata.reversePaths[0].state == 1 and x > cdata.reversePaths[0].interface:
            #    cdata.freeenergyhistos[0].addweightedPointToHisto(x,weight)

    ######Print output
    filename = os.path.join(os.getcwd(),"freeenergy.dat")
    histofilename = os.path.join(os.getcwd(),"freeenergyhisto.dat")
    cdata.freeenergyhistos[0].outputfreeenergyHisto(filename)
    cdata.freeenergyhistos[0].outputCrossingHisto(histofilename)
