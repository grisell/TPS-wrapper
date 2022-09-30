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
import random
from math import log,fabs, exp, tanh, sqrt
import numpy as np

NALPHA=6
SEED=3
NCOLV=2
NTOT=200   #Monte Carlo cycles 
REL_DEVIATION=0.2
DISPX= 70.0
DISPY=0.08
iaverage=1387

NVAR=(NCOLV+1)*(NALPHA-2)+2


FILETYPE=2  #type 0= 3 columns; type 1 = 14 columns; type 2= 8 columns
CV2COLUMN=1
LARGENUMBER=1000000
temperature=1370.0
beta=1.0

SSTRAJFILE="trajectoryZ.00."
TRAJECTORYFILE="trajectory.00."
#TRAJECTORYFILECV2="spns2traj."
TRAJECTORYFILECV2="ns420traj."
#TRAJECTORYFILE = "newnstrajectory."
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
        self.freeenergyhistos = []
        self.committorhistos = []
        self.npaths = 0

       
        
        #self.histo = np.zeros(HISTOSIZE)
        #self.fehisto = np.zeros(shape=(HISTOSIZE,HISTOSIZECV2))
        #self.fehisto = np.zeros(HISTOSIZE)
        #self.comhisto = np.zeros(HISTOSIZE)
        #self.histosize=HISTOSIZE
        #self.histosizecv2=HISTOSIZECV2
        #self.histomin=HISTOMIN
        #self.histomax=HISTOMAX
        #self.histomincv2=HISTOMINCV2
        #self.histomaxcv2=HISTOMAXCV2
         

        self.PBhisto = np.zeros(NALPHA)
        self.PAhisto = np.zeros(NALPHA)
        
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

    def getaverageshootingpoints(self, iaverage=100): 
        """
          read path weights obtained previously in selfconsistent algorithm 
        """
        filename = os.path.join(os.getcwd(),"MLE","freeenergy2D.dat")
        numpoints = sum(1 for line in open(filename))
        avpoints=numpoints/iaverage
        cv1= [0]*avpoints
        cv2= [0]*avpoints
    
        c = 0
        s=0
        x=0.0
        y=0.0
        j=0
        for line in open( filename , 'r' ):            
            if c % iaverage != 0:
                x += float(line.split()[0])
                y += float(line.split()[1])
                s += 1
            else:
                if s > 0 :
                    x=x/s
                    y=y/s
                    cv1[j]=x
                    cv2[j]=y
                    #print c,s,x,y
                    j+=1

                s=0
                x=0.0
                y=0.0
            c +=1

        return (cv1, cv2, avpoints)

    def readString(self,nstrpoints,ncvs): 
        """
        read string points 
        """
        stringvar=np.zeros(shape=(nstrpoints,ncvs))
        filename = os.path.join(os.getcwd(),"MLE","string2D.txt")       
        c = 0
        for line in open( filename , 'r' ):
            for i in range(ncvs):
                
                x = float(line.split()[i]) #CV point of the string
                stringvar[c][i]=x
            c+= 1
        assert c == nstrpoints, "MLEString.readString - File length does not match nstrpoints. nstrpoints: {0:d} filelength: {1:d}".format(nrstpoints, c)
        return (stringvar)
        
    def voronoi(self,nstrpoints, stringvar,cv1,cv2, weight, interface,pointtoB):
        
        rmin=100000.0
        close_pt = -1
        x=cv1
        y=cv2
        
        r2=0.0
        for i in range(nstrpoints):
            dx=stringvar[i][0]-x
            dy=stringvar[i][1]-y
            r2= dx*dx + dy*dy
            if(r2 < rmin):
                rmin= r2
                close_pt=i
        
        if(pointtoB==1):
            self.PBhisto[close_pt] += weight
        else:
            if(pointtoB==0):
                self.PAhisto[close_pt] += weight
    
        return(close_pt)

    def mindist_avshootingpts(self, numavpoints, avstringvar,cv1,cv2):
        
        rmin=1000000.0
        close_pt = -1
        x=cv1
        y=cv2
        
        
        for i in range(numavpoints):
            #print avstringvar[0][i], avstringvar[1][i]
            dx=avstringvar[0][i]-x
            dy=avstringvar[1][i]-y
            r2= dx*dx + dy*dy
            if(r2 < rmin):
                rmin= r2
                close_pt=i
        
        rmin=sqrt(r2)
                    
        return(rmin)
    def check_string_consistency(self, nstrpoints, stringvar, deviation):
        
        close_pt1=close_pt2=0
        
        for j in range(1, nstrpoints-1):
            xj=stringvar[j][0]
            yj=stringvar[j][1]
            
            rmin1=1000000
            close_pt2=close_pt1=-1
            
            for i in range(j):
                dx=stringvar[i][0]-xj
                dy=stringvar[i][1]-yj
                r2= dx*dx + dy*dy
                if (r2 < rmin1):
                    rmin1 = r2
                    close_pt1=i

            rmin2=1000000
            for i in range(j+1,nstrpoints):
                dx=stringvar[i][0]-xj
                dy=stringvar[i][1]-yj
                r2 = dx*dx + dy*dy
                if (r2 < rmin2):
                    rmin2 = r2
                    close_pt2=i
            
            okay=( close_pt1 == j-1 ) & (close_pt2 == j+1)
            if not(okay):
                #print "String not okay"
                #print "point", j
                #print "closest interface 1:", close_pt1,"rmin1", rmin1 
                #print "closest interface 2:", close_pt2,"rmin2", rmin2
                return False 
                
            relative_dev=fabs((rmin1-rmin2)/(rmin1+rmin2))
            if(relative_dev > deviation):
                #print "String not okay"
                #print "point", j
                #print "closest interface 1:", close_pt1,"rmin1:", rmin1, "closest interface 2:", close_pt2,"rmin2:", rmin2, "relative_dev:", relative_dev
                return False
            
        return True    
                
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


    def outputpBhistogram(self,filename,nstrpoints,slope):
        of = open(filename,"w")
        for i in range(nstrpoints):
            y = slope*(-0.5 +  (i +0.5) /nstrpoints)
            of.write("%d %.18f %.18f\n" % (i, self.PBhisto[i],0.5*(1.0+tanh(y))))

    def outputstring2D(self,filename,nstrpoints,stringvar,stringrc):
        of = open(filename,"w")
        for i in range(nstrpoints):
            of.write("%d %.18f %.18f %.18f\n" % (i, stringvar[i][0],stringvar[i][1],stringrc[i]))
                       
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
    rewrite=0

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
    #cvalue=cdata.readcacb() # Note: cA=cvalue[0], cB=cvalue[1]  (scaling factors)
    


    """
       Compute scaling factor of the weights (counting real number of paths)
    """
    numinterpaths=fabs(start-stop)
    numinter=len(cdata.kernels.kernelPaths)
    #integral=[0.0]*numinter
    numpaths=[0]*numinter
    
    for i in cdata.kernels.kernelPaths:  
        numpaths[i/2]= cdata.checknumpaths(start,stop,i)
    

    """
       Maximum likelihood estimation along the string variable using voronoi mapping.
    """
    nstrpoints = NALPHA
    ncvs=NCOLV
    lnL=-100000
    lnLrc=-100000
    dispx= DISPX
    dispy=DISPY
    ntot=NTOT
    nlnLrc=0
    sumlnLrc=0.0
    avlnLrc=0.0
    deviation=REL_DEVIATION
    string_old = np.zeros(shape=(nstrpoints,ncvs))
    avstring = np.zeros(shape=(nstrpoints,ncvs))
    reject=0
    accept=0

    nlnL=0
    sumlnL=0.0
    avlnL=0.0
    rejectrc=0
    acceptrc=0
    """
    Reading average shooting points

    """
    avpoints=cdata.getaverageshootingpoints(iaverage)
    print avpoints[0], avpoints[1]    
    
    """
       Getting string points
    """
    
    string =cdata.readString(nstrpoints,ncvs)

    lnLfilename="avlnL_2D.dat"
    lnLrcfilename="avlnLstring_2D.dat"
    oflnL = open(lnLfilename,"w")
    oflnLrc = open(lnLrcfilename,"w")
    Ndata=0.0

    stringrc = np.zeros(shape=(nstrpoints))
    avstringrc = np.zeros(shape=(nstrpoints))

    navstringsteps=0
    
    for m in range(ntot):
        #dispx=DISPX
        #dispy=DISPY
        #if (m%30==0.0):
        #    lnLrc=-100000
            
        lnLrc_old=lnLrc
            
        kpoint=random.randint(1,nstrpoints-2)
        #if (kpoint==0 or kpoint==nstrpoints-1):
        #    dispx= 0.0
        #print kpoint
        """
        Select string point to modify 
        """

        string_old_x=string[kpoint][0]
        string_old_y=string[kpoint][1]
        
        for i in range(nstrpoints):
            string_old[i][0]=string[i][0]
            string_old[i][1]=string[i][1]
            #print string_old[i][0], string_old[i][1]
            
        """
        Check string consistency 
        """

        if not(cdata.check_string_consistency(nstrpoints, string, deviation)):
            print "Error: string is not consistent"
        
        """
        Random displacement of string point 
        """

        dx=(random.random()-0.5)*dispx
        dy=(random.random()-0.5)*dispy
        string[kpoint][0]= string_old_x +  dx
        string[kpoint][1]= string_old_y +  dy
        
         
        while (cdata.check_string_consistency(nstrpoints, string, deviation)==False):
            dx=(random.random()-0.5)*dispx
            dy=(random.random()-0.5)*dispy
            string[kpoint][0]= string_old_x +  dx
            string[kpoint][1]= string_old_y +  dy
        
        if (cdata.check_string_consistency(nstrpoints, string, deviation)==True):
            print "now string is ok"

        print "moving string point:", kpoint, "x=", string[kpoint][0], "y=", string[kpoint][1]
        #for i in range(nstrpoints):
        #    print  string[i][0],  string[i][1]
        
    #"""
    #   Calculate the normalization factors for the histograms per interface
    #"""
    #for dirnumber in range(start,stop): 
    #    dirstring = "%07d" % dirnumber
    #    for i in cdata.kernels.kernelPaths:
    #        filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
    #        cdata.readFullTrajectory(filename, i) 
    #        above = cdata.paths[i].getPointsBeyondInterface()
    #        for point in above:
    #            integral[i/2] += 1
     
    
        """
        Read order parameter trajectories, assign path weights and compute free energy and committor histograms
        """
        cdata.PBhisto = np.zeros(NALPHA)
        cdata.PAhisto = np.zeros(NALPHA)
        
        for dirnumber in range(start,stop): 
            dirstring = "%07d" % dirnumber
    
            for i in cdata.kernels.kernelPaths:   
                """
                read trajectories 
                """
            
                filename= os.path.join(cdata.paths[i].nfsladir,dirstring,"trajectory.00."+dirstring+".dat")
                cdata.readFullTrajectory(filename, i)
                
                if (FILETYPE==0):
                    newfilename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,TRAJECTORYFILECV2+dirstring+".dat")
                    cdata.readnewFullTrajectory(newfilename, i)
                else:
                    if (FILETYPE==1):
                        newfilename= os.path.join(cdata.cv2paths[i].nfsladir,dirstring,TRAJECTORYFILECV2+dirstring+".dat")
                        cdata.readcv2FullTrajectory(newfilename, i)
                    if (FILETYPE==2):
                        newfilename= os.path.join(cdata.ns4paths[i].nfsladir,dirstring,TRAJECTORYFILECV2+dirstring+".dat")
                        cdata.readns4FullTrajectory(newfilename, i) 
                #newfilename= os.path.join(cdata.newpaths[i].nfsladir,dirstring,TRAJECTORYFILE+dirstring+".dat")
                #cv2filename= os.path.join(cdata.cv2paths[i].nfsladir,dirstring,TRAJECTORYFILECV2+dirstring+".dat")

                #cdata.readnewFullTrajectory(newfilename, i)
                #cdata.readcv2FullTrajectory(cv2filename, i)
             
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
                wfilename=os.path.join(cdata.paths[i].nfsladir,dirstring,"weight_rep.dat")
                
                if rewrite==0:
                    if os.path.exists(wfilename):
                        #print "Reading path weights from files"
                        for line in open( wfilename , 'r' ):
                            weight = float(line.split()[0])
                        if typepath[0]:
                            pointtoB= 0
                        if typepath[1]:
                            pointtoB= 1
                        if typepath[2]:
                            pointtoB= 1
                        if typepath[3]:
                            pointtoB= 0 
                if rewrite==1:
                    
                  
                  if typepath[0]: # AA paths                
                    
                     maxlambda=cdata.paths[i].getMaxTrajectoryValue()
                     nwi=-1
  
                     for k in reversed(range(0,cdata.fninter)):
                        intervalue = cdata.interfaces.interfaces[0][k]
                        if maxlambda > intervalue:
                            nwi=k
                            break

                    
                     #weight=cvalue[0]*fweight[1][nwi]
                     weight=fweight[1][nwi]
                     pointtoB= 0
                     scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                     weight *= scalefactor
                    
                     of = open(wfilename,"w")
                     of.write("%.30f\n" % (weight))
                     of.close()
                
                
                  if typepath[1]: #AB paths
                    
                     nwi=cdata.fninter-1
                     #weight=cvalue[0]*fweight[1][nwi]
                     weight=fweight[1][nwi]
                     pointtoB= 1
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
                     #weight=cvalue[1]*bweight[1][nwi]
                     weight=bweight[1][nwi]
                     scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                     weight *= scalefactor
                     pointtoB= 1
                     of = open(wfilename,"w")
                     of.write("%.30f\n" % (weight))
                     of.close()
                                            
                  if typepath[3]: #BA paths
                     nwi=0
                     #weight=cvalue[1]*bweight[1][nwi]
                     weight=bweight[1][nwi]
                     pointtoB= 0
                     scalefactor=cdata.paths[i].getscalenumbersinterfaces(numpaths[nwi],numinterpaths)
                     weight *= scalefactor
                     of = open(wfilename,"w")
                     of.write("%.30f\n" % (weight))
                     of.close()

                """
                Rescale the weights with the normalization factors per interface (see Eq. 6 in JCP 133, 174109 (2010) )
                """   
            
                #weight=weight/integral[i/2]

               
                ss1=20.0
                ss2=750.0
                flagshoot=0

                if (FILETYPE==0):
                    print "do you want to project the FE on ns vs ns?"
                    
                    break

                if (FILETYPE==1):
                    for (tispoint, point2) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.cv2paths[i].lastAcceptedFullTrajectory[0]):
                
                        if (tispoint[2]==flagshoot):
                            
                            if (tispoint[1]>=ss1 and tispoint[1]<=ss2):
                        
                                cdata.voronoi(nstrpoints, string,tispoint[1],point2[CV2COLUMN], weight,i,pointtoB)
                            
                                if m==0:
                                    Ndata += weight
                            
                        flagshoot= tispoint[2]
                        
                if (FILETYPE==2):
                    for (tispoint, point2) in zip(cdata.paths[i].lastAcceptedFullTrajectory[0], cdata.ns4paths[i].lastAcceptedFullTrajectory[0]):
                
                        if (tispoint[2]==flagshoot):

                            if (tispoint[1]>=ss1 and tispoint[1]<=ss2):
                                x=tispoint[1]/750.0
                                y=point2[CV2COLUMN]/450.0
                                #cdata.voronoi(nstrpoints, string,tispoint[1],point2[CV2COLUMN], weight,i,pointtoB)
                                cdata.voronoi(nstrpoints, string,x,y, weight,i,pointtoB)
                                if m==0:
                                    Ndata += weight
                            
                            
                        flagshoot= tispoint[2]

                    
                        
       
        if m==0:
            BIC= 0.5*NVAR*log(Ndata)
            print "BIC=", BIC
            
        """
        Calculating minium distance to average of the shooting points
        """
        
        #avpoints=cdata.getaverageshootingpoints(iaverage)
        numavpoints=avpoints[2]
        
        dmin_old=cdata.mindist_avshootingpts(numavpoints, avpoints,string_old_x,string_old_y)
        dmin_new=cdata.mindist_avshootingpts(numavpoints, avpoints,string[kpoint][0],string[kpoint][1])
        daux=-beta*(dmin_new-dmin_old)
        if (daux > 500.0):
            daux=500.0
        print dmin_old, dmin_new
                                    

        #lnL=0.0
        slope=nstrpoints

        
        #for i in range(nstrpoints):
        #    print  cdata.PBhisto[i], cdata.PAhisto[i], "Image No.", i

        #for i in range(nstrpoints):
        #    x=slope*(-0.5 +  (i +0.5)/nstrpoints)
        #    #print x, cdata.PBhisto[i], cdata.PAhisto[i]
        #    lnL += cdata.PBhisto[i]*log(0.5*(1.0+tanh(x))) + cdata.PAhisto[i]* log(0.5*(1.0-tanh(x)))


        rej=0
        
        #stringrc = np.zeros(shape=(nstrpoints))
        #avstringrc = np.zeros(shape=(nstrpoints))
        rcdisp=0.5
        if m==0:
            for i in range(nstrpoints):
                stringrc[i]=slope*(-0.5 + (i +0.5)/nstrpoints)

        ncycle1= 50;
        ncycle2 =40000;
        avlnL=0.0
        nlnL=0
        sumlnL = 0.0
        for n in range(ncycle2):
            lnL_old=lnL
            
            rcpoint=random.randint(0,nstrpoints-1)
            
            rcold=stringrc[rcpoint]
            count=0
            
            stringrc[rcpoint]=rcold+(random.uniform(0.0,1.0)-0.5)*rcdisp
            
            if(rcpoint > 0 and rcpoint< nstrpoints-1):
                
                while(((stringrc[rcpoint] < stringrc[rcpoint-1]) or (stringrc[rcpoint] > stringrc[rcpoint+1])) and (count < 100)):
                    stringrc[rcpoint]=rcold+(random.uniform(0.0,1.0)-0.5)*rcdisp
                    count +=1
                if (count>=100):
                    stringrc[rcpoint]=rcold
                    rej +=1
            
            if (rcpoint == 0):
                count=0
                while(stringrc[rcpoint] < -10.0 or stringrc[rcpoint] > stringrc[rcpoint+1] and (count < 100)):
                    stringrc[rcpoint]=rcold+(random.uniform(0.0,1.0)-0.5)*rcdisp
                    count +=1
                
                if (count>=100):
                    stringrc[rcpoint]=rcold
                    rej +=1
            
            if (rcpoint == nstrpoints-1):
                count=0
                while(stringrc[rcpoint] > 10.0 or stringrc[rcpoint] < stringrc[rcpoint-1] and (count < 100)):
                    stringrc[rcpoint]=rcold+(random.uniform(0.0,1.0)-0.5)*rcdisp
                    count +=1
                
                if (count>=100):
                    stringrc[rcpoint]=rcold
                    rej +=1
            
            
            lnL=0.0
            for i in range(nstrpoints):
                x=stringrc[i]
                lnL += cdata.PBhisto[i]*log(0.5*(1.0+tanh(x))) + cdata.PAhisto[i]* log(0.5*(1.0-tanh(x)))



            print beta, lnL, lnL_old
            aux= beta*(lnL-lnL_old)

            """
            Check exponential maximum value
            """
            if (aux > 500.0):
                aux=500.0
                
           
            rn= random.uniform(0.0,1.0)
            print "rn=", rn, "aux=",aux, "exp(aux)=",exp(aux)
            
            if (rn > exp(aux)):
                print "value of likelihood, lnL=", lnL, "aux=", aux
                #string[kpoint][0]= string_old_x
                #string[kpoint][1]= string_old_y
                lnL=lnL_old
                rejectrc +=1
                print "REJECTING RC MOVE: now checking old value"
        
                #for i in range(nstrpoints):
                #    string[i][0]= string_old[i][0]
                    #   string[i][1]= string_old[i][1]
                print "restored lnL =", lnL, "oldlnL=", lnL_old
            else:
	        acceptrc+=1
                print "RC move accepted  ln L =", lnL, "oldlnL = ", lnL_old
	        #evstringfile=os.path.join(os.getcwd(),"stringevol.dat")
    	        #cdata.outputstring2D(evstringfile,nstrpoints,string)
                #for i in range(nstrpoints):
                #    print stringrc[i]

            nlnL +=1
            sumlnL += lnL
            avlnL= sumlnL/nlnL
            if (n%500 == 0) :
                oflnL.write("%d %.18f  %.18f\n" % (nlnL, avlnL,lnL))
            

        
	"""
	optimizing lnL
        """
        """
	A=-7.0
        B=0.25
        Aoptsteps=50
        Boptsteps=50
        lnLmax=-1000000
        for j in range(Aoptsteps):
            lnL=0.0
            B=0.25
            for k in range(Boptsteps):
                #print j,k
                for i in range(nstrpoints):
                    x=slope*(-0.5 +  (i +0.5)/nstrpoints)
                    lnL += cdata.PBhisto[i]*log(0.5*(1.0+tanh(B*(x+A)))) + cdata.PAhisto[i]* log(0.5*(1.0-tanh(B*(x+A))))
                if (lnL > lnLmax):
                    lnLmax=lnL 
                    #print lnLmax
                    Amax= A
                    Bmax= B
                B += 0.017
            A += 0.28
        print "lnLmax=",lnLmax, "A=", Amax, "B=", Bmax
        lnL=lnLmax
        #print beta, lnL, lnL_old
        aux= beta*(lnL-lnL_old)
        """

        """
	optimizing lnL
        
	A=-7.0
        B=1.0
        Aoptsteps=50
        Boptsteps=50
        lnLmax=-1000000
        for j in range(Aoptsteps):
            lnL=0.0
            #B=0.25
            #for l in range(Boptsteps):
            for i in range(nstrpoints):
                x=slope*(-0.5 +  (i +0.5)/nstrpoints)
                lnL += cdata.PBhisto[i]*log(0.5*(1.0+tanh(B*(x+A)))) + cdata.PAhisto[i]* log(0.5*(1.0-tanh(B*(x+A))))
            if (lnL > lnLmax):
                lnLmax=lnL 
                    #print lnLmax
                Amax= A
                    #Bmax= B
                #B += 0.017
            A += 0.28

        print "lnLmax=",lnLmax, "A=", Amax
        #print "lnLmax=",lnLmax, "A=", Amax, "B=", Bmax
        lnL=lnLmax
        aux= beta*(lnL-lnL_old)
        """

        lnLrc=lnL
        aux2= beta*(lnLrc-lnLrc_old)


        
        """
        Check exponential maximum value
        """
        if (aux2 > 500.0):
            aux2=500.0
        
        
        rn2 = random.uniform(0.0,1.0)
        
        print "rn2=", rn2, "aux2=",aux2, "exp(aux2)=",exp(aux2), "daux",daux,"exp(daux)=",exp(daux)
        if (rn > exp(aux2)):
        #if (rn > exp(aux2) or rn > exp(daux)):
            print "value of likelihood, lnLrc=", lnLrc, "aux2=", aux2
            string[kpoint][0]= string_old_x
            string[kpoint][1]= string_old_y
            lnLrc=lnLrc_old
            reject +=1
            print "REJECTING STRING MOVE: now checking old value"
        
            for i in range(nstrpoints):
                string[i][0]= string_old[i][0]
                string[i][1]= string_old[i][1]
            print "restored lnL =", lnL, "oldlnL=", lnL_old
        else:
	    accept+=1
            print "string move accepted  ln Lrc =", lnLrc, "oldlnL = ", lnLrc_old
	    #evstringfile=os.path.join(os.getcwd(),"stringevol.dat")
    	    #cdata.outputstring2D(evstringfile,nstrpoints,string)
            
        nlnLrc +=1
        sumlnLrc += lnLrc
        avlnLrc=sumlnLrc/nlnLrc
        oflnLrc.write("%d %.18f %.18f\n" % (nlnLrc, avlnLrc, lnLrc))

        if m > 100:
            
            for i in range(nstrpoints):
                avstring[i][0] += string[i][0]
                avstring[i][1] += string[i][1]
                avstringrc[i] += stringrc[i]
            navstringsteps +=1



    """
    Final results
    """
    for i in range(nstrpoints):
         avstring[i][0]=avstring[i][0]/navstringsteps
         avstring[i][1]=avstring[i][1]/navstringsteps
         avstringrc[i]= avstringrc[i]//navstringsteps

         
    Cbic= avlnLrc - 0.5*NVAR*log(Ndata) 
    print "Accepted string moves=", accept, "rejected string moves=", reject
    print "Accepted rc moves=", acceptrc, "rejected rc moves=", rejectrc
    print "C_BIC=", Cbic

    filename = os.path.join(os.getcwd(),"pBhisto_2D.dat")
    cdata.outputpBhistogram(filename,nstrpoints,slope)
    stringfile=os.path.join(os.getcwd(),"finalstring_2D.dat")
    cdata.outputstring2D(stringfile,nstrpoints,string,stringrc)
    avstringfile=os.path.join(os.getcwd(),"averagestring_2D.dat")
    cdata.outputstring2D(avstringfile,nstrpoints,avstring,avstringrc)
    

        


   
