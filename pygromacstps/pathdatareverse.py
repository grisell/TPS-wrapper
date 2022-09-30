'''
Created on May 3, 2010

@author: wolf
'''

import options
import os
#import crossinghistogram

LARGENUMBER = 1000000000000.0

class pathdatareverse(object):
    def __init__(self,number,basedir,mode,forward,forwardpart,interface = -1,state=0):
        self.number = number
        self.interface = interface
        self.options = options.lammpstpsoptions(basedir,mode)
        
        #self.crossingHisto = crossinghistogram.crossinghisto(1000, -1.5, 1.5,forward)
        self.tisaccepted = False
        self.fullTrajectory = []
        self.fullTrajectoryblength = 0
        self.fullTrajectoryflength = 0
        
        self.lastAcceptedFullTrajectory = []
        self.lastAcceptedFullTrajectoryblength = 0
        self.lastAcceptedFullTrajectoryflength = 0
        
        # shooting time is the shooting position in the last accepted trajectory
        self.shootingTime = 0
       
        #try to use md.in for input data for lammps	
        self.options.readOptions(os.path.join(basedir,"options","runoptions.txt"), self.options.mdpoptions)
        self.options.readOptions(os.path.join(basedir,"options","runoptions.txt"), self.options.runoptions)
        
        self.forward = forward
        self.forwardpart = forwardpart
        self.finishedState = -1
        self.state = state
        
        self.status = 0
        
        if state == 0:
            namepart = "stateA" 
        else:
            namepart = "stateB"
        
        
        self.srcatchbase = os.path.join(self.options.paths["scratchpath"],"GTIS")
        if self.forward:
            self.nfsbaseworkdir = os.path.join(self.options.paths["initialpath"],namepart)
            self.nfsladir = os.path.join(self.options.paths["initialpath"],"la/"+namepart)
            self.nfsalldir = os.path.join(self.options.paths["initialpath"],"all/"+namepart)
            
            self.baseworkdir = os.path.join(self.options.paths["scratchpath"],"GTIS/"+namepart)
            self.ladir = os.path.join(self.options.paths["scratchpath"],"GTIS/la/"+namepart)
            self.alldir = os.path.join(self.options.paths["scratchpath"],"GTIS/all/"+namepart)
            
        else:
            print "Stable State is always forward"
        
        if self.forwardpart:
            self.workdir = os.path.join(self.baseworkdir,"forward")
        else:
            print "Stable State simulation has only forward part"
        self.nfsconfigstore = os.path.join(self.nfsladir,"startconfs")
        self.latemp = os.path.join(self.workdir,"latemp")
    
#    def checkFinishedTIS(self,position,states):	#JR: orig, but states not needed anymore
    def checkFinishedTIS(self,position):
        finished = False
        if self.status == 1:
            if self.state == 0:
                if position > self.interface:
                    finished = True
            else:
                if position < self.interface:
                    finished = True
        else:
            if self.state == 0:
                if position < self.interface:
                    self.status = 1
            else:
                if position > self.interface:
                    self.status = 1
        return finished
        
    
    def checkAcceptedTIS(self,states):
        accepted = False
        if self.state == 0:
            accepted = self.fullTrajectory[-1][1] > self.interface 
        else:
            accepted = self.fullTrajectory[-1][1] < self.interface 
        return accepted
    
    def getMaxTrajectoryValue(self):
        if self.forward:
            maxvalue = -LARGENUMBER
            for tpoint in self.lastAcceptedFullTrajectory:
                if tpoint[1] > maxvalue:
                    maxvalue = tpoint[1]
        else:
            maxvalue = LARGENUMBER
            for tpoint in self.lastAcceptedFullTrajectory:
                if tpoint[1] < maxvalue:
                    maxvalue = tpoint[1]
        return maxvalue
    
    def ifcrossstablestate(self,newpos,oldpos):
        if self.state == 0:
            if newpos >= self.interface  and  oldpos < self.interface:
                return True
        else:
            if newpos <= self.interface  and  oldpos > self.interface:
                return True
        return False
        
    def getscalenumberstablestates(self,numpathsA,numpathsB, interpaths):
        if self.state == 0:
            if numpathsA > 0:
                scalenumber= interpaths/numpathsA
            else:
                scalenumber=1.0
        else:
            if numpathsB > 0:
                scalenumber= interpaths/numpathsB
            else:
                scalenumber=1.0
        return scalenumber
