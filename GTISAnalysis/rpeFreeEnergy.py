'''
Created on Feb 30, 2015

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
import numpy as np

HISTOSIZE=3000
HISTOMIN=0.0
HISTOMAX=8788.0
#CUTOFF=0.001
LARGENUMBER=10000000
temperature=1370.0

SSTRAJFILE="trajectoryZ.00."
TRAJECTORYFILE="trajectory.00."
#TRAJECTORYFILE = "newnstrajectory."
#SSTRAJFILE="rnsnewtrajectory."
p0f=[0]*HISTOSIZE
p0b=[0]*HISTOSIZE

class rpeData(object):
    """
    Main TIS class which stores the paths, interfaces, stable states, and  paths. The class 
    also uses the gromacswrapper/lammps class, the filesystem class and helper.py. 
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
        self.freeenergyhistos = []
        self.committorhistos = []
        self.npaths = 0

        self.histo = np.zeros(HISTOSIZE)
        self.fehisto = np.zeros(HISTOSIZE)
        self.comhisto = np.zeros(HISTOSIZE)
        self.histosize=HISTOSIZE
        self.histomin=HISTOMIN
        self.histomax=HISTOMAX

         
        """
          append histograms for forward ensemble
        """
        for i in range(self.interfaces.ninterfaces[0]):
            
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.npaths += 1

        """
          append histograms for backward ensemble
        """    
        
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.newpaths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.newpaths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.npaths += 1
        
        """
            append histograms for stable state paths (reverse paths)
        """
        self.reversePaths = []
        self.newreversePaths = []

        if len(self.interfaces.reversepaths[0]):
            for rp in self.interfaces.reversepaths[0]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
                self.newreversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
                self.freeenergyhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=True))
        if len(self.interfaces.reversepaths[1]):
            for rp in self.interfaces.reversepaths[1]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
                self.newreversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
                self.freeenergyhistos.append(crossinghistogram.generichistogram(HISTOSIZE,HISTOMIN,HISTOMAX,forward=True))
        
        """
            List of interfaces and (directory) paths 
        """      
        self.kernels.generateKernelLists(self.npaths)
        
        """
            append histograms for forward ensemble
        """

        self.fninter = self.interfaces.ninterfaces[0]
        self.bninter = self.interfaces.ninterfaces[1]
        

    def readFullTrajectory(self,filename,pathnumber): 
        """
           read full trajectories for forward and backward ensembles 
        """
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
        
    def readnewFullTrajectory(self,filename,pathnumber): 
        """
           read full orde parameter trajectories for forward and backward ensembles 
        """
        traj = []
        self.newpaths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.newpaths[pathnumber].fullTrajectoryblength = 0
        self.newpaths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.newpaths[pathnumber].fullTrajectoryflength = 0
        self.newpaths[pathnumber].lastAcceptedFullTrajectory = []
        self.newpaths[pathnumber].fullTrajectory = []
        self.newpaths[pathnumber+1].lastAcceptedFullTrajectory = []
        self.newpaths[pathnumber+1].fullTrajectory = []
        
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            if int(float(raw[2])) == 0:
                self.newpaths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                self.newpaths[pathnumber].fullTrajectoryblength += 1
            else:
                self.newpaths[pathnumber].lastAcceptedFullTrajectoryflength += 1
                self.newpaths[pathnumber].fullTrajectoryflength += 1
        
        self.newpaths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.newpaths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.newpaths[pathnumber].fullTrajectory.append(traj[:])
        self.newpaths[pathnumber+1].fullTrajectory.append(traj[:])

    

    def readreverseFullTrajectory(self,filename,pathnumber): 
        """
          read full trajectory for stable state trajectories
        """
        traj = []
        self.reversePaths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.reversePaths[pathnumber].fullTrajectoryblength = 0
        self.reversePaths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.reversePaths[pathnumber].fullTrajectoryflength = 0
        self.reversePaths[pathnumber].lastAcceptedFullTrajectory = []
        self.reversePaths[pathnumber].fullTrajectory = []
        
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            self.reversePaths[pathnumber].lastAcceptedFullTrajectoryflength += 1
            self.reversePaths[pathnumber].fullTrajectoryflength += 1
        
        self.reversePaths[pathnumber].lastAcceptedFullTrajectory.append(traj[:]) 
        self.reversePaths[pathnumber].fullTrajectory.append(traj[:])
        
    
    def readreversenewFullTrajectory(self,filename,pathnumber): 
        """
          read full order parameter trajectory for stable state trajectories
        """
        traj = []
        self.newreversePaths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.newreversePaths[pathnumber].fullTrajectoryblength = 0
        self.newreversePaths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.newreversePaths[pathnumber].fullTrajectoryflength = 0
        self.newreversePaths[pathnumber].lastAcceptedFullTrajectory = []
        self.newreversePaths[pathnumber].fullTrajectory = []
         
        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
            self.newreversePaths[pathnumber].lastAcceptedFullTrajectoryflength += 1
            self.newreversePaths[pathnumber].fullTrajectoryflength += 1
        
        self.newreversePaths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])       
        self.newreversePaths[pathnumber].fullTrajectory.append(traj[:])
        

    def readpathweights(self,ninter,forward): 
        """
          read path weights obtained previously in selfconsistent algorithm 
        """
        PAnorm= [0]*ninter
        weight= [0]*ninter
        filename = forward and os.path.join(os.getcwd(),"selfconsistent","fpathweights.txt") or os.path.join(os.getcwd(),"selfconsistent","bpathweights.txt")       
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[2]) #Get path weights 
            z = float(line.split()[3]) # Normalization factor 
            PAnorm[x]=z
            weight[x] = y
            c+= 1
        return (PAnorm, weight)

    def readcacb(self): 
        """
          Read scaling factors ca and cb from the previous selfconsistent calculation
        """
        filename = os.path.join(os.getcwd(),"selfconsistent_fe","cacb.dat")       
        c = 0
        for line in open( filename , 'r' ):
            ca = float(line.split()[3]) 
            cb = float(line.split()[4]) 
	#print ca, cb 
        return (ca, cb)
        
    def getBox(self,value): 
        """
        Convert the value to the histogrambox
        """
        return (float(self.histosize)*(float(value) - float(self.histomin)))/(float(self.histomax-self.histomin))
    
    def getValue(self,hbox): 
        """
        Convert the histogram index to the according value
        """
        return (float(hbox)*float(self.histomax-self.histomin) )/float(self.histosize) + float(self.histomin)
    def addweightedPointToHisto(self,value,weight):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self.getBox(value)))
        self.fehisto[hbox] += weight
    
    def addweightedPointTocommittorHisto(self,value,weight):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self.getBox(value)))
        self.comhisto[hbox] += weight

    def outputcommitorgHisto(self,filename):
        """
        Print the histogram
        """
        of = open(filename,"w")
        for i in range(self.histosize):
            if self.fehisto[i] != 0.0:                                                                               
                of.write("%.18f %.18f\n" % (self.getValue(i), self.comhisto[i]/self.fehisto[i]))  
        of.close()
    
    def outputfreeenergyHisto(self,filename,unit):
        """
        Print the histogram
        """
        max=-1.0
        for i in range(self.histosize):
            if  self.fehisto[i] >  max:
                max = self.fehisto[i]

        of = open(filename,"w")
        for i in range(self.histosize):
            if self.fehisto[i]>0.0 :                                                                               
                of.write("%.18f  %.18f  %.18f %.18f\n" % (self.getValue(i),-log((self.fehisto[i])/max)*unit, self.fehisto[i], self.fehisto[i]/max))  
        of.close()

        
    def checknumpaths(self,start,stop,i):
        """
          compute the actual number of complete stable state paths to re-scale the weights
        """       
        pos=[0]*LARGENUMBER
        c=0
        cross_state=0
        count=0
        for dirnumber in range(start,stop):
            #count = 0
            dirstring = "%07d" % dirnumber
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            self.readFullTrajectory(filename, i) 
            for point in self.paths[i].lastAcceptedFullTrajectory[0]:       
                pos[c]=float(point[1])
                
                if c > 0:
                    if self.paths[i].ifcrossinterface(pos[c],pos[c-1]):
                        cross_state=1
                        
                    if self.paths[i].ifinstablestate(pos[c],cdata.orderparameters.op):
                        if cross_state==1:
                            count +=1
                        cross_state = 0
                c += 1
            
        return count
        
#***********************************************create reweighted path ensemble**********************************************************************************
if __name__=="__main__":
    basedir = os.path.join(os.getcwd(),"..")
    cdata = rpeData(basedir,"tis")   
    rstart=0
    rstop=0

    """
        Select number of trajectories to be read for the interfaces and the stable states  
    """
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "start": #starting trajectory - interfaces
            start = int(float(argument[1])) 
        elif argument[0] == "stop": #end trajectory - interfaces
            stop = int(float(argument[1]))
        elif argument[0] == "rstart": #starting trajectory - stable states
            rstart = int(float(argument[1]))
        elif argument[0] == "rstop": #end trajectory -  stable states
            rstop = int(float(argument[1]))
        else:
            print "Arguments required should be one of those: start, stop, rstart, rstop"

    
    fweight = cdata.readpathweights(cdata.fninter,forward=True) #forward path weights
    bweight = cdata.readpathweights(cdata.bninter,forward=False) #backward path weights
    cvalue=cdata.readcacb() # Note: cA=cvalue[0], cB=cvalue[1]  (scaling factors)
   
    """
       Compute scaling factor of the weights (counting real number of paths)
    """
    numinterpaths=fabs(start-stop)
    numinter=len(cdata.kernels.kernelPaths)
    integral=[0.0]*numinter
    numpaths=[0]*numinter

    for i in cdata.kernels.kernelPaths:  
            numpaths[i/2]= cdata.checknumpaths(start,stop,i)
            #print numpaths[i/2]
    """
       Calculate the normalization factors for the histograms per interface
    """
    
    for dirnumber in range(start,stop): 
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
            flagshoot=0
            for point in cdata.paths[i].lastAcceptedFullTrajectory[0]:
                if (point[2]==flagshoot):
                    integral[i/2] += 1
                flagshoot= point[2]
    """
    for dirnumber in range(start,stop): 
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
            above = cdata.paths[i].getPointsBeyondInterface()
            for point in above:
                integral[i/2] += 1
    """         
    """
      Read order parameter trajectories, assign path weights and compute free energy and committor histograms
    """
    for dirnumber in range(start,stop): 
        dirstring = "%07d" % dirnumber
    
        for i in cdata.kernels.kernelPaths:   
            """
               read trajectories 
            """
            
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            newfilename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,TRAJECTORYFILE+dirstring+".dat")
            cdata.readFullTrajectory(filename, i) 
            cdata.readnewFullTrajectory(newfilename, i)
             
            """
               assing type to path: AA, AB, BA, BB
            """
            typepath = cdata.paths[i].checkAcceptedTISpathtype(cdata.orderparameters.op) 
            #print typepath
            """
               Scale the for the real number of paths 
            """
           #scalefactor=1.0
            scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[i/2],numinterpaths)
            #print scalefactor
            
            """
               Calculate/ assign weights to the paths
            """
            wfilename=os.path.join(cdata.paths[i].nfsladir,dirstring,"weight.dat")
            #if os.path.exists(wfilename):
            #     for line in open( wfilename , 'r' ):
            #            weight = float(line.split()[0])
            #else:
            
            
            if typepath[0]: # AA paths                
                    
                maxlambda=cdata.paths[i].getMaxTrajectoryValue()
                nwi=-1
  
                for k in reversed(range(0,cdata.fninter)):
                    intervalue = cdata.interfaces.interfaces[0][k]
                    if maxlambda > intervalue:
                        nwi=k
                        break

                    
                weight=cvalue[0]*fweight[1][nwi]
                scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                weight *= scalefactor
                    
                of = open(wfilename,"w")
                of.write("%.30f\n" % (weight))
                of.close()
                
                
            if typepath[1]: #AB paths
                    
                nwi=cdata.fninter-1
                weight=cvalue[0]*fweight[1][nwi]
                scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                weight *= scalefactor
                    
                of = open(wfilename,"w")
                of.write("%.30f\n" % (weight))
                of.close()
                    
            if typepath[2]: #BB paths
                     
                minlambda=cdata.paths[i].getMaxTrajectoryValue()
                nwi=-1  
                for k in range(0,cdata.bninter):
                    intervalue = cdata.interfaces.interfaces[1][k]
                    if minlambda < intervalue:
                        nwi=k
                        break
                weight=cvalue[1]*bweight[1][nwi]
                scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                weight *= scalefactor
                of = open(wfilename,"w")
                of.write("%.30f\n" % (weight))
                of.close()
                                            
            if typepath[3]: #BA paths
                nwi=0
                weight=cvalue[1]*bweight[1][nwi]
                scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                weight *= scalefactor
                of = open(wfilename,"w")
                of.write("%.30f\n" % (weight))
                of.close()
            """
            Rescale the weights with the normalization factors per interface (see Eq. 6 in JCP 133, 174109 (2010) )
            """   
            
            weight=weight/integral[i/2]

            """
            Create weighted free energy and committor histograms
            """  
            above = cdata.paths[i].getPointsBeyondInterface() 
            
            """
            for point in above:
                
                    if point > cdata.interfaces.interfaces[0][0]:
                        cdata.addweightedPointToHisto(point,weight)                 
                        if  cdata.paths[i].checifendinstateB(cdata.orderparameters.op):
                            cdata.addweightedPointTocommittorHisto(point,weight)
                
            """
            flagshoot=0
            for point in cdata.paths[i].lastAcceptedFullTrajectory[0]:
            
                if (point[2]==flagshoot):
                    x=float(point[1])
                    
                    if x > cdata.interfaces.interfaces[0][0]:
                        cdata.addweightedPointToHisto(x,weight)                 
                        if  cdata.paths[i].checifendinstateB(cdata.orderparameters.op):
                            cdata.addweightedPointTocommittorHisto(x,weight)
                flagshoot= point[2]
            
    """
       include stable state histograms if rstart and rstop != 0
    """

#For stable states ensembles:
    countA=0
    countB=0
    oldposA=LARGENUMBER
    oldposB=-LARGENUMBER
    
    """
       compute the actual number of complete stable state paths to re-scale the weight
    """

    for dirnumber in range(rstart,rstop):
        dirstring = "%07d" % dirnumber
	#reversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"rnsnewtrajectory."+dirstring+".acc")
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
    print countA, countB
    
    """
       assign weights to the stable states and compute the weighted histograms: 
    """
    ssintegral=0 
    for dirnumber in range(rstart,rstop):
        dirstring = "%07d" % dirnumber
	reversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"trajectoryZ.00."+dirstring+".acc") 
        cdata.readreversenewFullTrajectory(reversefilename, -1)
        for point in cdata.reversePaths[0].lastAcceptedFullTrajectory[0]:
            x=float(point[1])
            if len(cdata.interfaces.reversepaths[0]):
                if cdata.reversePaths[0].state == 0: 
                    if x < cdata.reversePaths[0].interface:
                        ssintegral += 1

    for dirnumber in range(rstart,rstop):
        dirstring = "%07d" % dirnumber
        newreversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,SSTRAJFILE+dirstring+".acc")
	reversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"trajectoryZ.00."+dirstring+".acc")
        wreversefilename=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"weight.dat")
        cdata.readreversenewFullTrajectory(reversefilename, -1)
        scalenumber=cdata.reversePaths[0].getscalenumberstablestates(countA,countB,numinterpaths)
        #if os.path.exists(wreversefilename):
        #    for line in open(wreversefilename, 'r' ):
        #            weight = float(line.split()[0])
        #else:    
        if cdata.reversePaths[0].state == 0:
                weight=cvalue[0]*fweight[1][0]
                of = open(wreversefilename,"w")
                of.write("%.18f\n" % (weight))
                of.close()
        else:
                weight= cvalue[1]*bweight[1][cdata.bninter-1]
                of = open(wreversefilename,"w")
                of.write("%.18f\n" % (weight))
                of.close()
        
        weight = scalenumber*weight/ssintegral
        
        for point in cdata.newreversePaths[0].lastAcceptedFullTrajectory[0]:
            x=float(point[1])
            
            if len(cdata.interfaces.reversepaths[0]):
                if cdata.reversePaths[0].state == 0: 
                    if x < cdata.reversePaths[0].interface:
                        cdata.addweightedPointToHisto(x,weight)    
                        #cdata.freeenergyhistos[0].addweightedPointToHisto(x,weight)
            elif  len(cdata.interfaces.reversepaths[1]):
                if cdata.reversePaths[0].state == 1: 
                    if x > cdata.reversePaths[0].interface:
                        cdata.addweightedPointToHisto(x,weight)    
                        #cdata.freeenergyhistos[0].addweightedPointToHisto(x,weight)

                
    """
          print output (freeenergy.dat and commitor.dat)
    """
                
    eV=float(0.00013806488*temperature)/float(1.60217657)
    filename = os.path.join(os.getcwd(),"freeenergy.dat") 
    committorfilename = os.path.join(os.getcwd(),"committor.dat") 
    cdata.outputfreeenergyHisto(filename,eV)
    cdata.outputcommitorgHisto(committorfilename)
    

   
