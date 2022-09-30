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

STATEA=0.0
STATEB=900.0

HISTOSIZE=3000
HISTOSIZECV2=100

HISTOMIN=0.0
HISTOMAX=8788.0

HISTOMINCV2=0.0
HISTOMAXCV2=1.2

CV1COLUMN=1
CV2COLUMN=5
CUTOFF=0.001
LARGENUMBER=1000000
temperature=1370.0


SSTRAJFILE="trajectoryZ.00."
TRAJECTORYFILE="trajectory.00."
TRAJECTORYFILECV1="trajectory.00."
TRAJECTORYFILECV2="spns2traj."
TRAJECTORYFILECV3="ns420traj."
#TRAJECTORYFILECV4="clstrQ6traj."

#SSTRAJFILE="rnsnewtrajectory."


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
        self.cv2paths = []
        self.ns4paths= []
        self.cv4paths = []
        self.freeenergyhistos = []
        self.committorhistos = []
        self.npaths = 0

        self.histo = np.zeros(HISTOSIZE)
        self.fehisto = np.zeros(shape=(HISTOSIZE,HISTOSIZECV2))
        #self.fehisto = np.zeros(HISTOSIZE)
        #self.comhisto = np.zeros(HISTOSIZE)
        self.comhisto = np.zeros(shape=(HISTOSIZE,HISTOSIZECV2))  
        self.histosize=HISTOSIZE
        self.histosizecv2=HISTOSIZECV2
        self.histomin=HISTOMIN
        self.histomax=HISTOMAX
        self.histomincv2=HISTOMINCV2
        self.histomaxcv2=HISTOMAXCV2
         
        """
          append histograms for forward ensemble
        """
        for i in range(self.interfaces.ninterfaces[0]):
            
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.newpaths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.cv2paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.cv2paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.ns4paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.ns4paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.cv4paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.cv4paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
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
            self.cv2paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.cv2paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.ns4paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.ns4paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.cv4paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.cv4paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
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

    def readcv2FullTrajectory(self,filename,pathnumber):
        """                                                                                                                                           read full orde parameter trajectories for forward and backward ensembles
        """
        traj = []
        self.cv2paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.cv2paths[pathnumber].fullTrajectoryblength = 0
        self.cv2paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.cv2paths[pathnumber].fullTrajectoryflength = 0
        self.cv2paths[pathnumber].lastAcceptedFullTrajectory = []
        self.cv2paths[pathnumber].fullTrajectory = []
        self.cv2paths[pathnumber+1].lastAcceptedFullTrajectory = []
        self.cv2paths[pathnumber+1].fullTrajectory = []

        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),float(raw[2]), float(raw[3]),float(raw[4]),float(raw[5]),float(raw[6]),float(raw[7]),float(raw[8]),float(raw[9]),float(raw[10]), float(raw[11]), float(raw[12]),  float(raw[13]), int(float(raw[14]))])
            if int(float(raw[14])) == 0:
                self.cv2paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                self.cv2paths[pathnumber].fullTrajectoryblength += 1
            else:
                self.cv2paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
                self.cv2paths[pathnumber].fullTrajectoryflength += 1

        self.cv2paths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.cv2paths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.cv2paths[pathnumber].fullTrajectory.append(traj[:])
        self.cv2paths[pathnumber+1].fullTrajectory.append(traj[:])


    def readcv4FullTrajectory(self,filename,pathnumber):
        """                                                                                                                                           read full orde parameter trajectories for forward and backward ensembles
        """
        traj = []
        self.cv4paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.cv4paths[pathnumber].fullTrajectoryblength = 0
        self.cv4paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.cv4paths[pathnumber].fullTrajectoryflength = 0
        self.cv4paths[pathnumber].lastAcceptedFullTrajectory = []
        self.cv4paths[pathnumber].fullTrajectory = []
        self.cv4paths[pathnumber+1].lastAcceptedFullTrajectory = []
        self.cv4paths[pathnumber+1].fullTrajectory = []

        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),float(raw[2]), float(raw[3]),float(raw[4]),float(raw[5]),float(raw[6]),float(raw[7])])
            if int(float(raw[7])) == 0:
                self.cv4paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                self.cv4paths[pathnumber].fullTrajectoryblength += 1
            else:
                self.cv4paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
                self.cv4paths[pathnumber].fullTrajectoryflength += 1

        self.cv4paths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.cv4paths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.cv4paths[pathnumber].fullTrajectory.append(traj[:])
        self.cv4paths[pathnumber+1].fullTrajectory.append(traj[:])

    def readns4FullTrajectory(self,filename,pathnumber):
        """                                                                                                                                           read full orde parameter trajectories for forward and backward ensembles
        """
        traj = []
        self.ns4paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.ns4paths[pathnumber].fullTrajectoryblength = 0
        self.ns4paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.ns4paths[pathnumber].fullTrajectoryflength = 0
        self.ns4paths[pathnumber].lastAcceptedFullTrajectory = []
        self.ns4paths[pathnumber].fullTrajectory = []
        self.ns4paths[pathnumber+1].lastAcceptedFullTrajectory = []
        self.ns4paths[pathnumber+1].fullTrajectory = []

        for line in open(os.path.join(filename),"r"):
            raw = line.split()
            traj.append([int(float(raw[0])),float(raw[1]),float(raw[2]), float(raw[3]),float(raw[4]),float(raw[5]),float(raw[6]),float(raw[7]),float(raw[8])])
            if int(float(raw[8])) == 0:
                self.ns4paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                self.ns4paths[pathnumber].fullTrajectoryblength += 1
            else:
                self.ns4paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
                self.ns4paths[pathnumber].fullTrajectoryflength += 1

        self.ns4paths[pathnumber].lastAcceptedFullTrajectory.append(traj[:])
        self.ns4paths[pathnumber+1].lastAcceptedFullTrajectory.append(traj[:])
        self.ns4paths[pathnumber].fullTrajectory.append(traj[:])
        self.ns4paths[pathnumber+1].fullTrajectory.append(traj[:])

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
            return (ca, cb)
        
    #def getBox(self,value): 
        """
        Convert the value to the histogrambox
        """
    #    return (float(self.histosize)*(float(value) - float(self.histomin)))/(float(self.histomax-self.histomin))
    
    def getBox(self,value,histomin,histomax,histosize):
        """                                                                                                    
        Convert the value to the histogrambox                                                                  
        """
        return (float(histosize)*(float(value) - float(histomin)))/(float(histomax-histomin))
    
    #def getValue(self,hbox): 
        """
        Convert the histogram index to the according value
        """
    #    return (float(hbox)*float(self.histomax-self.histomin) )/float(self.histosize) + float(self.histomin)

    def getValue(self,hbox,histomin,histomax,histosize):
        """                                                                                                                                        Convert the histogram index to the according value                                                                                  
        """
        return (float(hbox)*float(histomax-histomin) )/float(histosize) + float(histomin)

    def addweightedPointToHisto(self,value,weight,histomin,histomax):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self.getBox(value,histomin,histomax)))
        self.fehisto[hbox] += weight

    def addweightedPointToHisto2D(self,value,valuecv2, weight,histomin,histomax, histomincv2,histomaxcv2):
        """                                                                                                                                   
        Add a single point to the histogram                                                                                                           """
        hbox = int(round(self.getBox(value,histomin,histomax,self.histosize)))
        hboxcv2 = int(round(self.getBox(valuecv2,histomincv2,histomaxcv2,self.histosizecv2)))
        
        self.fehisto[hbox][hboxcv2] += weight

    def addweightedPointTocommittorHisto(self,value,weight,histomin,histomax):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self.getBox(value,histomin,histomax)))
       
        self.comhisto[hbox] += weight

    def addweightedPointTocommittorHisto2D(self,value,valuecv2, weight,histomin,histomax, histomincv2,histomaxcv2):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self.getBox(value,histomin,histomax,self.histosize)))
        hboxcv2 = int(round(self.getBox(valuecv2,histomincv2,histomaxcv2,self.histosizecv2)))
        self.comhisto[hbox][hboxcv2] += weight

    def outputcommitorgHisto2D(self,filename):
        """
        Print the histogram
        """
        of = open(filename,"w")
        for i in range(self.histosize):
            for j in range(self.histosizecv2):
                
                if self.fehisto[i][j] != 0.0:                                                                               
                    of.write("%.18f %.18f %.18f\n" % (self.getValue(i,self.histomin,self.histomax,self.histosize),self.getValue(j,self.histomincv2,self.histomaxcv2,self.histosizecv2), self.comhisto[i][j]/self.fehisto[i][j]))  
        of.close()
    
    #def outputfreeenergyHisto(self,filename,unit):
        """
        Print the histogram
        """
    #    max=-1.0
    #    for i in range(self.histosize):
    #        if  self.fehisto[i] >  max:
    #            max = self.fehisto[i]

#        of = open(filename,"w")
#        for i in range(self.histosize):
#            if self.fehisto[i]>0.0 :                                                                               
#                of.write("%.18f  %.18f  %.18f %.18f\n" % (self.getValue(i,self.histomin,self.histomax),-log((self.fehisto[i])/max)*unit, self.fehisto[i], self.fehisto[i]/max))  
#        of.close()

    def outputfreeenergyHisto2D(self,filename,unit):
        """                                                                                                                                    
        Print the histogram                                                                                                                    
        """
        max=-1.0
        for i in range(self.histosize):
            for j in range(self.histosizecv2): 
                if  self.fehisto[i][j] >  max:
                    max = self.fehisto[i][j]

        of = open(filename,"w")
        
        for i in range(self.histosize):
            for j in range(self.histosizecv2):
                if self.fehisto[i][j]>0.0 :
                    of.write("%.18f %.18f %.18f  %.18f %.18f\n" % (self.getValue(i,self.histomin,self.histomax,self.histosize),self.getValue(j,self.histomincv2,self.histomaxcv2,self.histosizecv2),-log((self.fehisto[i][j])/max)*unit, self.fehisto[i][j], self.fehisto[i][j]/max))
		#else:
		    #of.write("%.18f %.18f %.18f  %.18f %.18f\n" % (self.getValue(i,self.histomin,self.histomax,self.histosize),self.getValue(j,self.histomincv2,self.histomaxcv2,self.histosizecv2),self.fehisto[i][j], self.fehisto[i][j], self.fehisto[i][j]))
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
    for dirnumber in range(start,stop):
        dirstring = "%07d" % dirnumber
        for i in cdata.kernels.kernelPaths:
            filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
            cdata.readFullTrajectory(filename, i)
            #above = cdata.paths[i].getPointsBeyondInterface()
	    flagshoot=0
            for point in cdata.paths[i].lastAcceptedFullTrajectory[0]:
		if (point[2]==flagshoot):
                    integral[i/2] += 1
                flagshoot= point[2]
                
 
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
            newfilename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,TRAJECTORYFILECV1+dirstring+".dat")
            ns4filename= os.path.join(cdata.ns4paths[i].nfsladir,dirstring,TRAJECTORYFILECV3+dirstring+".dat")
            cv2filename= os.path.join(cdata.cv2paths[i].nfsladir,dirstring,TRAJECTORYFILECV2+dirstring+".dat")
            #cv4filename= os.path.join(cdata.cv4paths[i].nfsladir,dirstring,TRAJECTORYFILECV4+dirstring+".dat")
            cdata.readFullTrajectory(filename, i) 
            cdata.readnewFullTrajectory(newfilename, i)
            cdata.readcv2FullTrajectory(cv2filename, i)
            #cdata.readns4FullTrajectory(ns4filename, i)
           # cdata.readcv4FullTrajectory(cv4filename, i)
             
            """
               assing type to path: AA, AB, BA, BB
            """
            typepath = cdata.paths[i].checkAcceptedTISpathtype(cdata.orderparameters.op) 
            #print typepath
            """
               Scale the for the real number of paths 
            """
            scalefactor=1.0
           # scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[i/2],numinterpaths)
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

            #above = cdata.paths[i].getPointsBeyondInterface() 
            #for point in above:
            ss1=STATEA
            ss2=STATEB
            flagshoot=0
            
                
            """
	    for (point, point2, point3) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0], cdata.ns4paths[i].lastAcceptedFullTrajectory[0]):

                  #if point > cdata.interfaces.interfaces[0][0]:
                if (point[2]==flagshoot):
                    if cdata.paths[i].forward:
                        if point[1]  >  cdata.paths[i].interface:
                              #if (point[1] >=ss1 and point[1]<= ss2):
			    #scalepoint=8788*point2[CV2COLUMN]
			    #cdata.addweightedPointToHisto2D(point[CV1COLUMN],scalepoint,weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                            #cdata.addweightedPointToHisto2D(point3[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2) 
                            cdata.addweightedPointToHisto2D(point[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)            
                    else:
                        if point[1]  <  cdata.paths[i].interface:
                            cdata.addweightedPointToHisto2D(point[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2) 
                            #cdata.addweightedPointToHisto2D(point3[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                flagshoot= point[2]  
                 
            

            for (point, point2) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0]):

                  #if point > cdata.interfaces.interfaces[0][0]:
                if (point[2]==flagshoot):
                    
                              #if (point[1] >=ss1 and point[1]<= ss2):
                    cdata.addweightedPointToHisto2D(point[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)            
                    
                flagshoot= point[2]
            
	    """

	    #for (point, point2, point3,point4) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0], cdata.ns4paths[i].lastAcceptedFullTrajectory[0],cdata.cv4paths[i].lastAcceptedFullTrajectory[0]):
            for (point, point2) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0]):
            #for (point, point2, point3) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0], cdata.ns4paths[i].lastAcceptedFullTrajectory[0]):
                  #if point > cdata.interfaces.interfaces[0][0]:
                if (point[2]==flagshoot):
                   #if point3[CV2COLUMN]<=point[CV1COLUMN]: 
                            if (point[1] >=ss1 and point[1]<= ss2):
                            #scalepoint=8788*point2[CV2COLUMN]
                            #cdata.addweightedPointToHisto2D(point[CV1COLUMN],scalepoint,weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                            	cdata.addweightedPointToHisto2D(point[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2) 
                            	#cdata.addweightedPointToHisto2D(point[CV1COLUMN],point2[CV2COLUMN]*point[1],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                                #cdata.addweightedPointToHisto2D(point[CV1COLUMN],point3[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                                if  cdata.paths[i].checifendinstateB(cdata.orderparameters.op):
                                    #cdata.addweightedPointTocommittorHisto2D(point[CV1COLUMN],point2[CV2COLUMN]*point[1],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                                    cdata.addweightedPointTocommittorHisto2D(point[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
                                
                            #cdata.addweightedPointToHisto2D(point3[CV1COLUMN],point2[CV2COLUMN],weight,cdata.histomin,cdata.histomax,cdata.histomincv2,cdata.histomaxcv2)
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
    filename = os.path.join(os.getcwd(),"freeenergy2D.dat") 
    committorfilename = os.path.join(os.getcwd(),"committor2D.dat") 
    cdata.outputfreeenergyHisto2D(filename,eV)
    cdata.outputcommitorgHisto2D(committorfilename)
    

   
