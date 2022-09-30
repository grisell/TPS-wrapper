'''
Created on May 3, 2010

@author: wolf
'''

import options
import os
#import crossinghistogram

LARGENUMBER = 1000000000000.0

distopi = 0

class pathdata(object):
    def __init__(self,number,basedir,mode,forward,forwardpart,interface = -1):
        self.number = number
        self.interface = interface
        self.options = options.lammpstpsoptions(basedir,mode)

        
        #self.crossingHisto = crossinghistogram.crossinghisto(1000, -2.0, 2.0,forward)
        self.tisaccepted = False
        self.tpsaccepted = False
        self.fullTrajectory = []
        self.fullTrajectoryblength = 0
        self.fullTrajectoryflength = 0
        
        self.lastAcceptedFullTrajectory = []
        self.lastAcceptedFullTrajectoryblength = 0
        self.lastAcceptedFullTrajectoryflength = 0
        
        self.lastAcceptedDirString = 0
        
        # shooting time is the shooting position in the last accepted trajectory
        self.shootingTime = 0
        
        #this is not really used with lammps -> read md options from md.in when running 	
        self.options.readOptions(os.path.join(basedir,"options","runoptions.txt"), self.options.mdpoptions)
	#but this should remain here
        self.options.readOptions(os.path.join(basedir,"options","runoptions.txt"), self.options.runoptions)
        
        self.forward = forward
        self.forwardpart = forwardpart
        self.finishedState = -1
        self.srcatchbase = os.path.join(self.options.paths["scratchpath"],"GTIS")
        if self.forward:
            self.nfsbaseworkdir = os.path.join(self.options.paths["initialpath"],"f%02d" % (number,))
            self.nfsladir = os.path.join(self.options.paths["initialpath"],"la/f%02d" % (number,))
            self.nfsalldir = os.path.join(self.options.paths["initialpath"],"all/f%02d" % (number,))
            
            self.baseworkdir = os.path.join(self.options.paths["scratchpath"],"GTIS/f%02d" % (number,))
            self.ladir = os.path.join(self.options.paths["scratchpath"],"GTIS/la/f%02d" % (number,))
            self.alldir = os.path.join(self.options.paths["scratchpath"],"GTIS/all/f%02d" % (number,))
            
        else:
            self.nfsbaseworkdir = os.path.join(self.options.paths["initialpath"],"b%02d" % (number,))
            self.nfsladir = os.path.join(self.options.paths["initialpath"],"la/b%02d" % (number,))
            self.nfsalldir = os.path.join(self.options.paths["initialpath"],"all/b%02d" % (number,))
            
            self.baseworkdir = os.path.join(self.options.paths["scratchpath"],"GTIS/b%02d" % (number,))
            self.ladir = os.path.join(self.options.paths["scratchpath"],"GTIS/la/b%02d" % (number,))
            self.alldir = os.path.join(self.options.paths["scratchpath"],"GTIS/all/b%02d" % (number,))
        
        
        
        if self.forwardpart:
            self.workdir = os.path.join(self.baseworkdir,"forward")
        else:
            self.workdir = os.path.join(self.baseworkdir,"backward")
        self.latemp = os.path.join(self.workdir,"latemp")
        
        
    def checkFinished(self,position,states):
        finished = position > states[1] or position < states[0]
        self.finishedState = -1
        if self.forwardpart:
            if position > states[1]:
                self.finishedState = 4
        else:
            if position < states[0]:
                self.finishedState = 1
        return finished
    
    def checkFinishedTIS(self,position,states):
        finished = position > states[1] or position < states[0]
        self.finishedState = -1
        if position > states[1]:
            self.finishedState = 4
        if position < states[0]:
            self.finishedState = 1
        return finished
    
    def checkFinishedTPS(self,position,states):
        finished = position > states[1] or position < states[0]
        self.finishedState = -1
        if position > states[1]:
            self.finishedState = 4
        if position < states[0]:
            self.finishedState = 1
        return finished
    
    def checkAcceptedTIS(self,states,log):
        accepted = False
        
        if self.forward:
            if self.fullTrajectory[distopi][0][1] < states[0].A[1]:
                if self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
			accepted = True
		elif self.fullTrajectory[distopi][-1][1] < states[0].A[1]:
                	for position in self.fullTrajectory[0]:
                    		log.log.debug("forwardinner " + str(position[1]) + " " + str(self.interface) + " " + str(self.fullTrajectory[0][0][1]) + " " + str(states[0].A[1]))
                    		if position[1] > self.interface:
                        		accepted = True
        else:
            
            if self.fullTrajectory[distopi][0][1] > states[0].B[0]:
                if self.fullTrajectory[distopi][-1][1] < states[0].A[1]:
			accepted = True
		elif self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
                	for position in self.fullTrajectory[0]:
                    		log.log.debug("backwardinne " + str(position[1]) + " " + str(self.interface) + " " + str(self.fullTrajectory[0][-1][1]) + " " + str(states[0].B[0]))
                    		if position[1] < self.interface:
                        		accepted = True
        log.log.debug("final " + str(accepted))
                    
        return accepted
    
    def checkAcceptedTPS(self,states):
        return self.fullTrajectory[distopi][0][1] < states[0].A[1] and self.fullTrajectory[distopi][-1][1] > states[0].B[0]
    
    
    def getMaxTrajectoryValue(self):
        if self.forward:
            maxvalue = -LARGENUMBER
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] > maxvalue:
                    maxvalue = tpoint[1]
        else:
            maxvalue = LARGENUMBER
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] < maxvalue:
                    maxvalue = tpoint[1]
        return maxvalue
    
    def getPointsBeyondInterface(self):
        points = []
        if self.forward:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] > self.interface:
                    points.append(tpoint[1])
        else:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] < self.interface:
                    points.append(tpoint[1])
        return points


    def getPointsBelowInterface(self):                                                 #DS for FE histogram of the stable states
        points = []
        if self.forward:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] < self.interface:
                    points.append(tpoint[1])
        else:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if tpoint[1] > self.interface:
                    points.append(tpoint[1])
        return points


    def getPointsBeyondInterfaceforward(self,max_interface):                             #DS for boundary histograms
        points = []
        if self.forward == True:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if  tpoint[1] < max_interface and tpoint[1] > self.interface:
                    points.append(tpoint[1])
        return points

    def getPointsBeyondInterfacebackward(self,min_interface):                             #DS for boundary histograms
        points = []
        if self.forward == False:
            for tpoint in self.lastAcceptedFullTrajectory[distopi]:
                if  tpoint[1] > min_interface and tpoint[1] < self.interface:
                    points.append(tpoint[1])
        return points

    
    
    def checkAcceptedTISAB(self,states):
        ab = False
        if self.forward:
            if self.fullTrajectory[distopi][0][1] < states[0].A[1] and self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
                ab = True
        else:
            if self.fullTrajectory[distopi][0][1] > states[0].B[0] and self.fullTrajectory[distopi][-1][1] < states[0].A[1]:
                ab = True
        
        return ab
    
    
    
    
    def checkAcceptedTISpathtype(self,states):
        aa = False
        ab = False
        bb = False
        ba = False

        if self.forward:
            if self.fullTrajectory[distopi][0][1] < states[0].A[1] and self.fullTrajectory[distopi][-1][1] < states[0].A[1]:
                aa = True
            if self.fullTrajectory[distopi][0][1] < states[0].A[1] and self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
                ab = True
                
        else:
            if self.fullTrajectory[distopi][0][1] > states[0].B[0] and self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
                bb = True
            if self.fullTrajectory[distopi][0][1] > states[0].B[0] and self.fullTrajectory[distopi][-1][1] < states[0].A[1]:
                ba = True
        return (aa,ab,bb,ba)

    def ifcrossinterface(self,newpos,oldpos):
        if self.forward == True:
            if newpos >= self.interface  and  oldpos < self.interface:
                return True
        else :
            if newpos <= self.interface  and  oldpos > self.interface:
                return True
        return False

    def ifinstablestate(self, pos, states):
        if pos <= states[0].A[1] or pos >= states[0].B[0]:
            return True
        return False
    def getscalenumbersinterfaces(self,numpaths, interpaths):
        if numpaths > 0:
            scalenumber= interpaths/numpaths
        else:
            scalenumber=1.0
        
        return scalenumber

    def checifendinstateB(self,states):
        if self.fullTrajectory[distopi][-1][1] > states[0].B[0]:
            return True
        return False
