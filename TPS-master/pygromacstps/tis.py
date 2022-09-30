
'''
Created on May 3, 2010

@author: Wolfgang Lechner
'''

import parser
import pathdata
import pathdatareverse
import filesystem
import gtpslogging
import interfaces
import orderparameters
import helpers 
import kernels
import qsubsystem
import wrappers 
import random
import multiprocessing
import os
import subprocess
from math import fabs
from math import exp
from decimal import Decimal
import time
import sys
import logging
import options as opt
import pathstats as pstat

#create a class object
optobj = opt.tpsoptions()
stat = pstat.pathstats()

#SRM:setting up logger
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False

#different logger for logging shoot points
shootlog = logging.getLogger('log2')
handler = logging.FileHandler("gtpslog.shoot.txt")
formatter = logging.Formatter('%(message)s')
handler.setFormatter(formatter)
shootlog.addHandler(handler)
shootlog.setLevel(logging.DEBUG)
shootlog.propagate = False

distopi = 0

class gromacstis(object):
    """
    Main TIS class which stores the paths, interfaces, stable states, and  paths. The class 
    also uses the gromacswrapper class, the filesystem class and helper.py. The logger is 
    self.log.
    """
    def __init__(self,basedir=".",mode="initial",kernel=0):
        self.basedir = basedir
        self.mode = mode
#        self.cores = multiprocessing.cpu_count()  #JR: this cannot be used if run with a queueing system!
        self.cores = 2
        self.wrapper = wrappers.wrapper()
        self.distparser = parser.gdistparser()
        self.filesystem = filesystem.filesystem()
        self.interfaces = interfaces.interfaces()
        self.orderparameters = orderparameters.orderparameters()
        self.helper = helpers.helpers()
        self.kernels = kernels.kernels(kernel)
        self.qsubsystem = qsubsystem.qsubsystem()
   

        
        #Initialize the logger
        if kernel=="head":
            self.log = gtpslogging.log("debug",basedir,kernel)
        elif kernel=="reverse":
            self.log = gtpslogging.log("debug",basedir,kernel)   
        else:
            self.log = gtpslogging.log("info",basedir,kernel)
        
        self.log.log.debug("logfile created")
        self.log.log.debug(str(self.cores) + " CPUs detected")
        
        self.interfaces.readInterfaces(os.path.join(basedir,"options","interfaces.txt"))
        self.log.log.debug("Read Interfaces Forward : " + str(self.interfaces.interfaces[0]))
        self.log.log.debug("Read Interfaces Backward : " + str(self.interfaces.interfaces[1]))
        
        self.kernels.readKernelOptions(os.path.join(basedir,"options","kerneloptions.txt"))
        
        #read the stables states from a file
        self.orderparameters.readOP(os.path.join(basedir,"options","orderparameters.txt"))
        self.log.log.debug("Read OP : " + str(self.orderparameters.op[0]))

        
        """
        This array holds all the information about the paths. Each trajectory consists of a 
        forward-part and a backward-part. The trajectory can either be a forward trajectory
        or a backward trajectory. Forward trajectroies start in A.
        """
        
        self.paths = []
        self.npaths = 0
        
        for i in range(self.interfaces.ninterfaces[0]):
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=True,interface=self.interfaces.interfaces[0][i]))
            self.paths.append(pathdata.pathdata(i,basedir,mode,forward=True,forwardpart=False,interface=self.interfaces.interfaces[0][i]))
            self.npaths += 1
        for i in range(self.interfaces.ninterfaces[1]):
            n = i + self.interfaces.ninterfaces[0]
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=True,interface=self.interfaces.interfaces[1][i]))
            self.paths.append(pathdata.pathdata(n,basedir,mode,forward=False,forwardpart=False,interface=self.interfaces.interfaces[1][i]))
            self.npaths += 1
        
        #SRM: Check directory paths,scratch, and other runoptions
        #---------------------------------------
        #SRM: check for kerneloptions
        if (self.interfaces.ninter < self.kernels.nkernels):
            self.kernels.nkernels = self.interfaces.ninter
            logger.info("No of kernels more than no of interfaces. The value is reset to no of interfaces")

        self.reversePaths = []
        for rp in self.interfaces.reversepaths[0]:
            self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))    
        for rp in self.interfaces.reversepaths[1]:
            self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))
        
        self.kernels.generateKernelLists(self.npaths)
        
        self.interfacedir = int(float(self.paths[0].options.runoptions["interfacecoordinate"]))  
        print self.interfacedir                  
        
    """ 
    ************************************************************************
    
        Initial Procedure
        
        Set of functions used for the generation of initial trajectories
    
    ************************************************************************
    """
    
    
    def preperationFromStart(self,copyfiles = True):
        """
        Creating the directory structure and copying standardfiles there, given that copyfile = True
        """
        for i in self.kernels.kernelPathsAll:
            initdir = self.paths[i].options.paths["initialpath"]
            workdir = self.paths[i].workdir
            dirlist = [self.paths[i].srcatchbase,self.paths[i].baseworkdir,self.paths[i].workdir,os.path.join(workdir,'paths')]
            #dirlist = [self.paths[i].srcatchbase,self.paths[i].workdir,os.path.join(workdir,'paths')]
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            #dirlist = [os.path.join(initdir,'la'),self.paths[i].nfsbaseworkdir,self.paths[i].nfsladir,os.path.join(self.paths[i].nfsladir,"%07d" % (0))]
            dirlist = [os.path.join(initdir,'la'),self.paths[i].nfsladir,os.path.join(self.paths[i].nfsladir,"%07d" % (0))]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                #for myfile in self.paths[i].options.standardfiles:
                for myfile in self.paths[i].options.standardfiles.itervalues():
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",myfile), workdir)

	    if self.paths[i].options.runoptions["mddriver"] in ('lammps','lammps_parallel'):
	    	self.wrapper.modifyMDseed(workdir,self.paths[i].options)
	



    def shootingInitialFiles(self):
        """
        Perform a shooting move on a gro file.
        """
        for i in self.kernels.kernelPaths:
            conffile = os.path.join(self.paths[i].workdir,self.paths[i].options.dumpfiles["conffile"])
            bakfile = os.path.join(self.paths[i].workdir,self.paths[i].options.dumpfiles["bakfile"])
            bconffile = os.path.join(self.paths[i].workdir,self.paths[i].options.dumpfiles["bconffile"])
            self.helper.shootVelocities(self.paths[i].options,conffile, bconffile)  
            self.filesystem.moveFile(conffile, bakfile)
            self.filesystem.moveFile(bconffile, conffile)
            dest = self.paths[i+1].workdir
            self.filesystem.copyFiles(conffile, dest)
            self.reverseBackwardFile(i+1)

            
                 
                         
    def finalizeInitial(self):
        
        def _copyFiles(fromdir,destdir):
	#   self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
            for value in optobj.copyfiles1.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,value), destdir)
        
        dirnumber = "%07d" % 0
        for i in self.kernels.kernelPaths:
            if self.paths[i].tisaccepted:
                
                baseladir = os.path.join(self.paths[i].options.paths["initialpath"],"la")
                lapath = self.paths[i].nfsladir
                newladir = os.path.join(lapath,dirnumber)
                newladirf = os.path.join(newladir,"forward")
                newladirb = os.path.join(newladir,"backward")
                dirlist = [baseladir,lapath,newladir,newladirf,newladirb]
                self.filesystem.createDirList(dirlist, self.log.log)
                
                _copyFiles(os.path.join(self.paths[i].workdir), newladirf)
                _copyFiles(os.path.join(self.paths[i+1].workdir), newladirb)
        
        for i in range(len(self.reversePaths)):
            baseladir = os.path.join(self.reversePaths[i].options.paths["initialpath"],"la")
            lapath = self.reversePaths[i].nfsladir
            newladir = os.path.join(lapath,dirnumber)
            newladirf = os.path.join(newladir,"forward")
            dirlist = [baseladir,lapath,newladir,newladirf]
            self.filesystem.createDirList(dirlist, self.log.log)
            _copyFiles(os.path.join(self.reversePaths[i].workdir), newladirf)
    
    
    def deleteScratchFiles(self):
        for i in self.kernels.kernelPaths:
            scratch = os.path.join(self.paths[i].baseworkdir)
            os.system("rm -r " + scratch)
        
    """ 
    ************************************************************************
    
        TIS Procedure
        
        Set of functions that perform a TIS shooting
    
    ************************************************************************
    """

    def preperationTIS(self,copyfiles = True):
        """
        Creating the directory structure and copying standardfiles there, given that copyfile = True
        """
                       
        for i in self.kernels.kernelPathsAll:
            initdir = self.paths[i].options.paths["initialpath"]
            workdir = self.paths[i].workdir
            dirlist = [self.paths[i].srcatchbase,self.paths[i].baseworkdir,self.paths[i].workdir,os.path.join(workdir,'paths'),os.path.join(workdir,"latemp")]
            #dirlist = [self.paths[i].srcatchbase,self.paths[i].workdir,os.path.join(workdir,'paths'),os.path.join(workdir,"latemp")]
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            #dirlist = [os.path.join(initdir,'all'),os.path.join(initdir,'la'),self.paths[i].nfsbaseworkdir,self.paths[i].nfsladir]
            dirlist = [os.path.join(initdir,'all'),os.path.join(initdir,'la'),self.paths[i].nfsladir]
            
	    self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                #for myfile in self.paths[i].options.standardfiles:
                for myfile in self.paths[i].options.standardfiles.itervalues():
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",myfile), workdir)
           
	    if self.paths[i].options.runoptions["mddriver"] in ('lammps','lammps_parallel'):
		self.wrapper.modifyMDseed(workdir,self.paths[i].options)
  
    

    """ 
    ************************************************************************
    
        Shooting procedure
    
    ************************************************************************
    """
    
    def lastAcceptedToConf(self,dirnumber):
        """
        The last accepted trajectory is transformed into pathXX.dump files 
        """
        for i in self.kernels.kernelPathsAll:
            
            
            mpath = self.paths[i]
	    # JR: this means if forwardpart = True -> assign "forward", else assign "backward"
            fbdir = mpath.forwardpart and "forward" or "backward"
            nfsworkdir = os.path.join(mpath.nfsladir,dirnumber,fbdir)
            latemp = mpath.latemp
            for file in optobj.copyfiles4.itervalues():
                self.filesystem.copyFiles(os.path.join(nfsworkdir,file) , latemp)
            self.wrapper.TrajToConf(mpath.options,latemp)

    def shootMove0(self,dirnumber,pathnumber,sum):
        """
        Normal shooting move. Pick any random slice
        """
        rn = sum/2
        while True:
                rn = random.randint(1,sum-2)
                if self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] > self.orderparameters.op[distopi].A[1] and self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] < self.orderparameters.op[distopi].B[0]:
                    break
                c += 1
                if c > 100:
                    print "ERROR picking shootingpoint"
                    sys.stdout.flush()
                    break
        return rn                               


    def shootMove1(self,dirnumber,pathnumber,sum):
        """
        Pick slice closest to the interface
        """
        mindist = 1000000.0
        rn = sum/2
        for kr in range(1,sum-1):
            act = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][kr][1]
            distance = fabs(act - self.paths[pathnumber].interface)
            accept = (bool(self.paths[pathnumber].forward) == bool(act > self.paths[pathnumber].interface))
            self.log.log.debug("ptest " + " "+ str(accept)  + " "  + str(act)  + " "  + str(self.paths[pathnumber].interface) )  # JR: was commented out
            if distance < mindist and accept:
                mindist = distance
                rn = kr

        return rn


    def shootMove2(self,dirnumber,pathnumber,sum):
        """
        Pick shooting point from adjoint interfaces
        """

        #now loop everything over full traj
        #if a bad interface and shootinterface != -1, shootfrom the shootinterface
        trajfor = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==1]
        trajbak = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==0]
        #remove the first slice from forward traj
        trajfor = trajfor[1:]
        #put it back
        fulltraj = trajbak + trajfor
        
        if (self.paths[pathnumber].badinterface) and (self.paths[pathnumber].shootinterface!=-1):
                if self.paths[pathnumber].forward:
                        interfaces = self.interfaces.interfaces[0]
                else:
                        interfaces = self.interfaces.interfaces[1]
                index = self.paths[pathnumber].shootinterface

                indexbefore = index-1
                indexafter = index+1
                #correct the interfaces if they go out of bounds
                if indexbefore < 0:
                        indexbefore = index
                if indexafter > len(interfaces)-1:
                        indexafter = index
                op_low = interfaces[indexbefore]
                op_high = interfaces[indexafter]            
        else:


                
                if self.paths[pathnumber].forward:
                        interfaces = self.interfaces.interfaces[0]
                        index = interfaces.index(self.paths[pathnumber].interface)
                        #get the interface before and after
                        indexbefore = index-1
                        indexafter = index+1
                        #correct the interfaces if they go out of bounds
                        if indexbefore < 0:
                                indexbefore = index
                        if indexafter > len(interfaces)-1:
                                indexafter = index
                        op_low = interfaces[indexbefore]
                        op_high = interfaces[indexafter]
                        #if its a backward interface
                else:
                        interfaces = self.interfaces.interfaces[1]
                        #get index of current interface
                        index = interfaces.index(self.paths[pathnumber].interface)
                        #get the interface before and after
                        indexbefore = index-1
                        indexafter = index+1
                        #correct the interfaces if they go out of bounds
                        #for baclward, if its the first one, the lower limit is stable state A
                        if indexbefore < 0:
                                indexbefore = index
                        if indexafter > len(interfaces)-1:
                                indexafter = index
                        op_low = interfaces[indexbefore]
                        op_high = interfaces[indexafter]
        rn = sum/2
        possibleshootpoints = []
        possibleshootvalues = []
        #loop over all values, add them to possible shooting points
        for kr in range(1,sum-2):
            act = fulltraj[kr][1]
            if (act>op_low) and (act<op_high):
                possibleshootpoints.append(kr)
                possibleshootvalues.append(act)
            #now pick one point from this list
        if (len(possibleshootpoints)==0):
            logger.error("Unable to pick a shooting point. Check your interface setup. For now, switching to shootfrominterface=1")
            rn = self.shootMove1(dirnumber,pathnumber,sum)
        else:
            random_point = random.randint(0,len(possibleshootpoints)-1)
            rn = possibleshootpoints[random_point]

        return rn

    def shootMove3(self,dirnumber,pathnumber,sum):                
        """
        Shoot from beyond the interface
        """
        trajfor = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==1]
        trajbak = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==0]
        #remove the first slice from forward traj
        trajfor = trajfor[1:]
        #put it back
        fulltraj = trajbak + trajfor
        #now loop everything over full traj
        #now loop everything over full traj
        #shooting beyond interface option
        #first check if its first or last interface. In this case just switch to shootingmove 2
        if ((self.paths[pathnumber].interface==self.orderparameters.op[distopi].A[1]) or (self.paths[pathnumber].interface==self.orderparameters.op[distopi].B[0])): 
            rn = self.shootMove2(dirnumber,pathnumber,sum)
        else:
            #if its a forward interface
            if self.paths[pathnumber].forward:
                op_low = self.paths[pathnumber].interface
                op_high = self.orderparameters.op[distopi].B[0]
            #if its a backward interface
            else:
                op_low = self.orderparameters.op[distopi].A[1]
                op_high = self.paths[pathnumber].interface
        
            rn = sum/2
            possibleshootpoints = []
            possibleshootvalues = []
            #loop over all values, add them to possible shooting points
            for kr in range(1,sum-2):
                act = fulltraj[kr][1]
                if (act>op_low) and (act<op_high):
                    possibleshootpoints.append(kr)
                    possibleshootvalues.append(act)
            if len(possibleshootpoints)==0:
                logger.error("Unable to pick a shooting point. Check your interface setup. For now, switching to shootfrominterface=1")
                rn = self.shootMove1(dirnumber,pathnumber,sum)
            else:
                random_point = random.randint(0,len(possibleshootpoints)-1)
                rn = possibleshootpoints[random_point]
        return rn

    def shootMove4(self,dirnumber,pathnumber,sum):                
        """
        Shoot from before the interface
        """
        trajfor = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==1]
        trajbak = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==0]
        #remove the first slice from forward traj
        trajfor = trajfor[1:]
        #put it back
        fulltraj = trajbak + trajfor
        #now loop everything over full traj
        #now loop everything over full traj
        #shooting beyond interface option
        #if its a forward interface
        if ((self.paths[pathnumber].interface==self.orderparameters.op[distopi].A[1]) or (self.paths[pathnumber].interface==self.orderparameters.op[distopi].B[0])): 
            rn = self.shootMove2(dirnumber,pathnumber,sum)
        else:
            if self.paths[pathnumber].forward:
                op_low = self.orderparameters.op[distopi].A[1]
                op_high = self.paths[pathnumber].interface
            #if its a backward interface
            else:
                op_low = self.paths[pathnumber].interface
                op_high = self.orderparameters.op[distopi].B[0]
        
            rn = sum/2
            possibleshootpoints = []
            possibleshootvalues = []
            #loop over all values, add them to possible shooting points
            for kr in range(1,sum-2):
                act = fulltraj[kr][1]
                if (act>op_low) and (act<op_high):
                    possibleshootpoints.append(kr)
                    possibleshootvalues.append(act)
            if len(possibleshootpoints)==0:
                logger.error("Unable to pick a shooting point. Check your interface setup. For now, switching to shootfrominterface=1")
                rn = self.shootMove1(dirnumber,pathnumber,sum)
            else:
                random_point = random.randint(0,len(possibleshootpoints)-1)
                rn = possibleshootpoints[random_point]
        return rn

    def shootMove5(self,dirnumber,pathnumber,sum):
        """
        Biased shooting with a biasing function. From Jurasczek
        """
        #fout = open("shootest",'w')
        sigma = optobj.runoptions["sigma"]
        minprob = optobj.runoptions["minprob"]
        s2 = 2.0*sigma*sigma
        if (self.paths[pathnumber].shootinterface!=-1):
                biasval = self.paths[pathnumber].shootinterface
        else:
                biasval = self.paths[pathnumber].interface
        #logger.info("biasval")
        #logger.info(biasval)
        trajfor = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==1]
        trajbak = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==0]
        #remove the first slice from forward traj
        trajfor = trajfor[1:]
        #put it back
        fulltraj = trajbak + trajfor
        sliceprob = []
        netprob = 0.00
        krns = []
        shootvalues = []
        if self.paths[pathnumber].forward:
                a = -1
        else:
                a = 1
        #probabilities are calculated and assigned
        for kr in range(1,sum-2):
                act = fulltraj[kr][1]
                if ((act>self.orderparameters.op[distopi].A[1]) and (act<self.orderparameters.op[distopi].B[0])):
                        pp = exp(a*(act - biasval))/s2
                else:
                        pp =0.00
                sliceprob.append(pp)
                krns.append(kr)
                shootvalues.append(act)
                netprob+=pp
        if netprob!=0:
                for i in range(len(sliceprob)):
                        sliceprob[i]/=netprob
        #logger.info("probs")
        #logger.info(sliceprob)
        #filter out really low probs and scale the rest
        netprob = 0.0
        for i in range(len(sliceprob)):
                if sliceprob[i]<minprob:
                        sliceprob[i]=0
                netprob+=sliceprob[i]
        #round off to two decimal places and multiply by 100
        #logger.info("netprob")
        #logger.info(netprob)
        if netprob>0:
                for i in range(len(sliceprob)):
                        sliceprob[i]/=netprob
                        dummy = round(sliceprob[i],2)
                        sliceprob[i] = dummy*100
        #logger.info("corrected probs")
        #logger.info(sliceprob)

        rand_no = random.randint(0,100)
        #logger.info("randint")
        #logger.info(rand_no)

        #now select the slice
        netprob = 0.0
        selected_slice = -1
        #for i in range(len(sliceprob)):
        #        if rand_no <= netprob:
        #                selected_slice = i
        #                break;
        #        netprob+=sliceprob[i]
        for kr in range(1,sum-2):
                if rand_no <= netprob:
                        selected_slice = kr
                        break;
                netprob+=sliceprob[kr-1]


        #logger.info("selected_slice")
        #logger.info(selected_slice)

        if selected_slice!= -1:
                rn = selected_slice
        else:
                rn = self.shootMove1(dirnumber,pathnumber,sum)
        #logger.info("corrected rn")
        #logger.info(rn)
        #logger.info("point value")
        #logger.info(shootvalues[kr-1])
        return rn

    def shootMove6(self,dirnumber,pathnumber,sum):
        """
        The new adaptive shooting move.
        """
        statfreq = optobj.runoptions["statfrequency"]
        minacct = optobj.runoptions["minacceptance"]
        stat.PathInfo()
        totpaths = stat.tot
        specialshoot=False
        #if its not the required frequency
        #pathdata has to be mmodified?
        #if its a forward slice
        if self.paths[pathnumber].forward:
            interfaces = self.interfaces.interfaces[0]
            index = interfaces.index(self.paths[pathnumber].interface)
            #received the index of the path.
            #now check if the acceptance is good or not
            intf_accp = stat.intf_accp[index]
            if intf_accp < minacct:
                self.paths[pathnumber].badinterface = True

            #checked if its a bad interface
            if (self.paths[pathnumber].badinterface) and (totpaths%statfreq==0):
                #this is a correction step
                #now find the reason for this bad acceptance
                #possibilities- f interface and BB paths
                logger.info("bad interface")
                logger.info(index)
                bbps = stat.intf_bb[index]
                aps = stat.intf_aa[index] + stat.intf_ab[index]
                #if too many bb paths - beyond transition state
                if (bbps>aps):
                        #find an interface close to transition state
                        maxabpaths = 0
                        goodintf = -1
                        for i in range(len(stat.intf_no)):
                                if (stat.intf_ab[i]>maxabpaths) and (stat.intf_accp>minacct):
                                        maxabpaths = stat.intf_ab[i]
                                        goodintf = stat.intf_no[i]
                        if goodintf==-1:
                                logger.info("No good interface at all!")
                        else:
                                self.paths[pathnumber].shootinterface = goodintf
                                specialshoot=True
                        logger.info("replacement interface")
                        logger.info(goodintf)

        else:
            interfaces = self.interfaces.interfaces[1]
            index = interfaces.index(self.paths[pathnumber].interface)
            #received the index of the path.
            #now check if the acceptance is good or not
            intf_accp = stat.intf_accp[index]
            if intf_accp < minacct:
                self.paths[pathnumber].badinterface = True

            #checked if its a bad interface
            if (self.paths[pathnumber].badinterface) and (totpaths%statfreq==0):
                #this is a correction step
                #now find the reason for this bad acceptance
                #possibilities- b interface and AA paths
                logger.info("bad interface")
                logger.info(index)
                aaps = stat.intf_aa[index]
                bps = stat.intf_ba[index] + stat.intf_bb[index]
                #if too many bb paths - beyond transition state
                if (aaps>bps):
                        #find an interface close to transition state
                        maxbapaths = 0
                        goodintf = -1
                        for i in range(len(stat.intf_no)):
                                if (stat.intf_ba[i]>maxbapaths) and (stat.intf_accp>minacct):
                                        maxbapaths = stat.intf_ba[i]
                                        goodintf = stat.intf_no[i]
                        if goodintf==-1:
                                logger.info("No good interface at all!")
                        else:
                                self.paths[pathnumber].shootinterface = goodintf
                                specialshoot=True
                        logger.info("replacement interface")
                        logger.info(goodintf)
        #now shoot - shoot#2
        #adjust shootmove accordingly
        if specialshoot:
                rn = self.shootMove5(dirnumber,pathnumber,sum)
        else:
                rn = self.shootMove2(dirnumber,pathnumber,sum)
        
        return rn

            






    def pickConfigurationsTIS(self,dirnumber):
        for i in self.kernels.kernelPaths:
            self.pickConfigurationLastAccepted(dirnumber,i)  
    
    def pickConfigurationLastAccepted(self,dirnumber,pathnumber):
        """
        A random configuration is taken from a trajectory and the shooting move is performed.
        and copy the files to the workdir of the path.
        """
        
        
        #fdir = os.path.join(self.paths[pathnumber].latemp)
        #flist = self.filesystem.getFileList(fdir, ("path*"+optobj.extension))
        #fsize = len(flist)
        fsize = self.paths[pathnumber].lastAcceptedFullTrajectoryflength
        #bdir = os.path.join(self.paths[pathnumber+1].latemp)
        #blist = self.filesystem.getFileList(bdir, ("path*"+optobj.extension))
        #bsize = len(blist)
        bsize = self.paths[pathnumber].lastAcceptedFullTrajectoryblength
        
        sum = fsize + bsize
        print "sum " ,sum, fsize, bsize
        print "traj ", len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi])
#        if sum != len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi]):
#            self.log.log.info(str(pathnumber) + " " + str(sum) + " " + str( len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi])))
        c = 0 

        #shootmoveselection
        if self.paths[pathnumber].options.runoptions["shootfrominterface"] == 0:
                rn = self.shootMove0(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 1:
                rn = self.shootMove1(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 2:
                rn = self.shootMove2(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 3:
                rn = self.shootMove3(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 4:
                rn = self.shootMove4(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 5:
                rn = self.shootMove5(dirnumber,pathnumber,sum)
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 6:
                rn = self.shootMove6(dirnumber,pathnumber,sum)
        else:
                logger.error("unknown shootfrominterface")
                raise SystemExit()

        
        #SRM: routine to log the shooting point
        trajfor = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==1]
        trajbak = [pos for pos in self.paths[pathnumber].lastAcceptedFullTrajectory[distopi] if pos[2]==0]
        #remove the first slice from forward traj
        trajfor = trajfor[1:]
        #put it back
        fulltraj = trajbak + trajfor
        act = fulltraj[rn][1]
        #the file used for logging shooting points would be called gtpslog.shoot.txt
        shootlog.info(("%d\t%f\t%d")%(pathnumber/2,act,len(fulltraj)))
        #end of routine

        #this needs to be changed
        #here we call the last accepted function as in the new path sampling code
        backwardpart = (rn < bsize)
        if backwardpart:
            sliceno = (bsize-1)-rn
            myfile = os.path.join(os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"backward","traj.dat"))
            #myfile = os.path.join(bdir,("path%d"+optobj.extension) % ((bsize-1)-rn))
        else:
            #myfile = os.path.join(fdir,("path%d"+optobj.extension) % (rn-bsize))
            sliceno = rn-bsize
            myfile = os.path.join(os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"forward","traj.dat"))


        
        self.log.log.debug("pick " +dirnumber + " " + str(pathnumber) +" " + str(rn) +  " " + str( sum) +  " " + str(fsize) + " " +  str(bsize) + " " +  str(len(self.paths[pathnumber].lastAcceptedFullTrajectory)) + " " + myfile)   # JR: was commented out
       
        self.paths[pathnumber].shootingTime = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][0]
        
        destfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,("conforig"+optobj.extension))
        shootfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,("confshoot"+optobj.extension))
        
        if os.path.exists(myfile):    
            #self.filesystem.copyFiles(myfile,destfile)
            self.wrapper.GetSliceFromTraj(myfile,destfile,sliceno)
        else:
#            self.log.log.ERROR("File " + myfile + " not found, this is a fatal error")  #JR original
            self.log.log.error("File " + myfile + " not found, this is a fatal error")
	    #print 'Could not identify shooting point'
	    #print 'rn = ',rn
	    #print 'fsize = ',fsize,' bsize = ',bsize
	    #print myfile
	    #print 'Exiting program...'
	    #sys.exit(1)
        self.helper.shootVelocities(self.paths[pathnumber].options,destfile,shootfile)
        
        for path in [self.paths[pathnumber],self.paths[pathnumber+1]]:
            self.filesystem.copyFiles(shootfile, os.path.join(path.workdir,("conf"+optobj.extension)))
            path.finishedstate = -1
        
        if backwardpart:
            reversepath = pathnumber 
        else:
            reversepath = pathnumber + 1
        self.reverseBackwardFile(reversepath)
        
        self.log.log.debug("shooting from " +str(rn)+ " file " + myfile + " backwardpart  " + str(backwardpart))
        
        #self.filesystem.deleteFileList(fdir, ("path*"+optobj.extension))
        #self.filesystem.deleteFileList(bdir, ("path*"+optobj.extension))


    def reverseBackwardFile(self,pathnumber):
        """
        Reverse all velocities from a lammps dump file.
        """
        conffile = os.path.join(self.paths[pathnumber].workdir,self.paths[pathnumber].options.dumpfiles["conffile"])
        bakfile = os.path.join(self.paths[pathnumber].workdir,self.paths[pathnumber].options.dumpfiles["bakfile"])
        bconffile = os.path.join(self.paths[pathnumber].workdir,self.paths[pathnumber].options.dumpfiles["bconffile"])
        self.helper.reverseVelocities(conffile, bconffile)
        self.filesystem.moveFile(conffile, bakfile)
        self.filesystem.moveFile(bconffile, conffile)


    """
    def gromppPath(self,path):
        workdir = path.workdir
        os.chdir(workdir)
        ok = self.wrapper.checkInputFiles(path.options, path.options.initoptions,                         \
                                     ["-c","-p","-n"],workdir)
        self.log.log.debug(workdir + str(ok) + " grompp executed")
        path.options.writeMdpFile(workdir,"md.mdp")
        
        cmd = self.wrapper.generateCommand(path.options, path.options.initoptions,                        \
                                      ["-c","-o","-p","-f","-n"],                                  \
                                      workdir,"grompp" )
        self.wrapper.executeCommand(cmd)
        os.chdir(self.basedir)
    """       
    
    def shootingQueue(self):
        """
        One of the main tps functions.
        Function that performes the mdrun on each of the paths.
        """
        
        """
        First step is to execute grompp in each of the paths workdirs to get a topol.tpr
        """
        # initialize the normal Paths
        for i in self.kernels.kernelPathsAll:
            self.wrapper.initializeMD(self.paths[i],self.basedir)
        
        # initialize the reverse Paths
        for rpath in self.reversePaths:
            self.wrapper.initializeMD(rpath,self.basedir)            

        
        """
        A process list is generated. It holds all the threads in a dictionary.
        The todo-list is the list of all processes that have to be executed.
        """
        processlist = {}
        reverseprocesslist = {}
        #make the todolist
        todolist = self.kernels.kernelPathsAll[:]
        #make todolist of all reverse paths
#        reversetodolist = range(len(self.reversePaths))
        reversetodolist = range(0)       # DS       
        # fill the cpu with jobs
        toploop = min(len(todolist),self.cores-1)
#	print todolist,toploop   #JR
        for i in range(toploop):
#            print i   #JR
            nextp = todolist.pop()
#	    print self.paths[nextp].workdir #JR
            os.chdir(self.paths[nextp].workdir)
            cmd =self.wrapper.generateMDRunCommand(self.paths[nextp].options, self.paths[nextp].workdir)
#	    tmp = os.getcwd() #JR
#	    print cmd #JR
#	    print tmp #JR
#	    print 'just before subproc',i				#JR
            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)	#JR
#	    print 'and just after subproc',i				#JR
            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
#	    print 'I am here\n'						#JR
            
#        sys.exit(1)    #JR
        finished = False
        """
        loop until all processes in the todo list are finished
        """
        while not(finished):    
            time.sleep(float(self.paths[0].options.runoptions["updatetime"]))
            finished = True
            for i in processlist.keys():
                #poll=None means it is still running
                if processlist[i].poll() == None:
                    finished = False
                    os.chdir(self.paths[i].workdir)
                    pathfinished = True
                    # we check all order parameters that we defined in the orderparameter option file
                    # if all of them give isPathInState = true then the path is finsihed
                    for opi in range(len(self.orderparameters.op)):
#			print 'Checking the order parameter\n'		#JR
#			sys.stdout.flush()   				#JR
                        pathstate,thisfinished,length = self.orderparameters.op[opi].isPathInState(self.paths[i],self.log)
                        #set the length of the path
                        self.paths[i].length = length
                        #SRM : this is the forward part
                        #although forward and backward paths have diffrent random numbers, we choose the random number of the forward part.
                        if ((i%2)==0):
                                rand_no = self.paths[i].length_rand
                        else:
                                rand_no = self.paths[i-1].length_rand
                                self.paths[i].length_rand = rand_no

                        #logging for debugging. Will be removed after checking
                        wdir = self.paths[i].workdir
                        swdir = wdir.split('/')
                        pathtype = swdir[-1]

                        #logger.info(("Pathtype %s pathstate %d")%(pathtype,pathstate))

                        totlength=self.paths[i].length
                        
                        #use the path length criteria only if tis mode is running
                        if self.mode=='tis':
                                #a soft criteria is implemented here. Forward and backward section, each is compared with the maximum path length
                                oldlength = self.paths[i].lastAcceptedFullTrajectoryflength + self.paths[i].lastAcceptedFullTrajectoryblength

                                targetlength = int(oldlength/self.paths[i].length_rand)
                                #logger.info(("Old %d Now %s Target %s")%(oldlength,totlength,targetlength)) 
                                
                                if (totlength > (int(oldlength/self.paths[i].length_rand))) and pathstate==0:
                                        #logger.info("Length exceeded")
                                        pathstate = 0
                                        thisfinished = True
                                #else:
                                        #logger.info("Looks good to go on")
                        
                        self.paths[i].finishedState = pathstate
                        pathfinished = thisfinished and pathfinished
                    if pathfinished:
                        self.log.log.debug("Finished " + str(i) + " " + str(self.paths[i].finishedState))
			print 'terminating process',i			#JR
			sys.stdout.flush()
                        processlist[i].terminate()
                        with open(os.path.join(self.paths[i].workdir,'step.in'),'w') as fout:
                                fout.write("run 1\n")
                        
                        if len(todolist)>0:
                            #if there is still something left to do we start from the todolist
                            nextp = todolist.pop()
                            os.chdir(self.paths[nextp].workdir)
                            cmd =self.wrapper.generateMDRunCommand(self.paths[nextp].options, self.paths[nextp].workdir)
			    
                            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE) #JR
                            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
                        else:
                            # otherwise, we start the equilibrium (called "reverse") paths
                            if len(reversetodolist) >0: 
                                nextr = reversetodolist.pop()
                                self.log.log.debug("Start Reverse Path "+ str(nextr))
                                os.chdir(self.reversePaths[nextr].workdir)
                                cmd =self.wrapper.generateMDRunCommand(self.reversePaths[nextr].options, self.reversePaths[nextr].workdir)

                                reverseprocesslist[nextr] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE) #JR
                                self.log.log.debug("Reversepath " +  str(nextr) + "started")
            
            #this is in case a job ends because the MD programme finished without being killed
            testproclist = True
            for i in processlist.keys():
	        #check if all process in the list are finished
		if processlist[i].poll() == None:
		    testproclist = False

            if len(todolist)>0  and testproclist == True:		#all jobs in processlist are finished but there is still stuff in todolist that needs to be added
                finished = False		
		#if there is still something left to do we start from the todolist
                print 'run next job, previous was finished in MD but not killed'
                nextp = todolist.pop()
                os.chdir(self.paths[nextp].workdir)
                cmd =self.wrapper.generateMDRunCommand(self.paths[nextp].options, self.paths[nextp].workdir)
			    
                processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE) #JR
                self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")


        if finished:
            # finished == True means that all paths from the todolist are finished.
            # now we kill all the remaining reversepaths
            for r in reverseprocesslist.keys():
                if reverseprocesslist[r].poll() == None:
                    reverseprocesslist[r].terminate()
                     
        self.log.log.debug("Shooting is finished and reversepaths killed")
        
        
    def shootingQueueDynamic(self):
        """
        Jeroen: I have added the reverse paths to this class. Here, first all paths from the todo list are executed.
        As soon as the todo list is finished up, but there are still paths running, the reverse paths are executed. 
        If all normal pahts are finished, the reverse paths are killed.  
        """
        
        """ 
        First step is to execute grompp in each of the paths workdirs to get a topol.tpr
        """
        # initialize the normal Paths
        for i in self.kernels.kernelPathsAll:
            self.wrapper.initializeMD(self.paths[i],self.basedir)
        
        # initialize the reverse Paths
        for rpath in self.reversePaths:
            self.wrapper.initializeMD(rpath,self.basedir)
        
        """
        A process list is generated. It holds all the threads in a dictionary.
        The todo-list is the list of all processes that have to be executed.
        """
        processlist = {}
        reverseprocesslist = {}
        #make the todolist
        todolist = self.kernels.kernelPathsAll[:]
        #make todolist of all reverse paths
        reversetodolist = range(len(self.reversePaths))
        
        # fill the cpu with jobs
        toploop = min(len(todolist),self.cores-1)
        for i in range(toploop):
            
            nextp = todolist.pop()
            os.chdir(self.paths[nextp].workdir)
            cmd =self.wrapper.generateMDRunCommand(self.paths[nextp].options, self.paths[nextp].workdir)
            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
            
            
        finished = False
        """
        loop until all processes in the todo list are finished
        """
        while not(finished):    
            time.sleep(float(self.paths[0].options.runoptions["updatetime"]))
            finished = True
            for i in processlist.keys():
                #poll=None means it is still running
                if processlist[i].poll() == None:
                    finished = False
                    os.chdir(self.paths[i].workdir)
                    pathfinished = True
                    # we check all order parameters that we defined in the orderparameter option file
                    # if all of them give isPathInState = true then the path is finsihed
                    for opi in range(len(self.orderparameters.op)):
                        pathstate,thisfinished = self.orderparameters.op[opi].isPathInState(self.paths[i],self.log)
                        self.paths[i].finishedState = pathstate
                        pathfinished = thisfinished and pathfinished
                    if pathfinished:
                        self.log.log.debug("Finished " + str(i) + " " + str(self.paths[i].finishedState))
                        processlist[i].terminate()
                        if len(todolist)>0:
                            #if there is still something left to do we start from the todolist
                            nextp = todolist.pop()
                            os.chdir(self.paths[nextp].workdir)
                            cmd =self.wrapper.generateMDRunCommand(self.paths[nextp].options, self.paths[nextp].workdir)
                            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
                        else:
                            # otherwise, we start the equilibrium (called "reverse") paths
                            if len(reversetodolist) >0: 
                                nextr = reversetodolist.pop()
                                self.log.log.debug("Start Reverse Path "+ str(nextr))
                                os.chdir(self.reversePaths[nextr].workdir)
                                cmd =self.wrapper.generateMDRunCommand(self.reversePaths[nextr].options, self.reversePaths[nextr].workdir)
                                reverseprocesslist[nextr] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                                self.log.log.debug("Reversepath " +  str(nextr) + "started")
        
        
        
        if finished:
            # finished == True means that all paths from the todolist are finished.
            # now we kill all the remaining reversepaths
            for r in reverseprocesslist.keys():
                if reverseprocesslist[r].poll() == None:
                    reverseprocesslist[r].terminate()
                     
        self.log.log.debug("Shooting is finished and reversepaths killed")
        
        


    def writeFinishedFiles(self,dirstring):
        for i in self.kernels.kernelPaths:
            lapath = self.paths[i].nfsladir
            newladir = os.path.join(lapath,dirstring)
            filename = os.path.join(newladir,"finished")
            of = open(filename,"w")
            of.write("finished")
            of.close()
            
    def finalizeShooting(self,dirstring):
        for i in self.kernels.kernelPathsAll:
            workdir = self.paths[i].workdir
            dest = os.path.join(workdir,"paths",dirstring)
            self.filesystem.createDir(dest)
            
#            self.filesystem.moveFile(os.path.join(workdir,"traj.xtc"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"traj.trr"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"topol.tpr"), dest)
             #self.filesystem.moveFile(os.path.join(workdir,"traj.dat"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
            #self.filesystem.copyFiles(os.path.join(workdir,"conf.dump"), dest)
            #self.filesystem.copyFiles(os.path.join(workdir,"md.in"), dest)
            #self.filesystem.deleteFileList(workdir,"#*")
#            self.filesystem.deleteFileList(workdir,"ener.edr")
            for file in optobj.copyfiles2.itervalues():
                self.filesystem.moveFile(os.path.join(workdir,file), dest)
            self.filesystem.deleteFileList(workdir,"#*")
            for file in optobj.deletefiles.itervalues():
                self.filesystem.deleteFileList(workdir,file)

    def deleteOldTrrFiles(self,dirnumber):
        dirlist = []
        for i in self.kernels.kernelPaths:
            #dirlist.append(os.path.join(self.paths[i].nfsalldir,dirnumber,"forward"))
            #dirlist.append(os.path.join(self.paths[i].nfsalldir,dirnumber,"backward"))
            dirlist.append(os.path.join(self.paths[i].nfsladir,dirnumber,"forward"))
            dirlist.append(os.path.join(self.paths[i].nfsladir,dirnumber,"backward"))
        
        for i in self.kernels.kernelPathsAll:
            mpath = self.paths[i]
            dirlist.append(os.path.join(mpath.workdir,"paths",dirnumber))
         
        for directory in dirlist:   
            self.filesystem.deleteFileList(directory,"*.trr")
    

    def checkAllTisPathsAccepted(self):
        """
        This function updates the tisaccepted variable
        """
        
        report = []
        for i in self.kernels.kernelPaths:
            self.getFullTrajectory(i,dirstring = "workdir")
           

            accepted,rejectreason,pathtype = self.paths[i].checkAcceptedTIS(self.orderparameters.op,self.log)
            self.paths[i].tisaccepted = accepted
            self.paths[i+1].tisaccepted = accepted
            
            if accepted or self.mode=="initial":
                self.paths[i].lastAcceptedFullTrajectory = []
                self.paths[i+1].lastAcceptedFullTrajectory = []
                for opi in range(len(self.orderparameters.op)):
                    for ki in range(2):
                        self.paths[i+ki].lastAcceptedFullTrajectory.append(self.paths[i].fullTrajectory[opi][:])
                        self.paths[i+ki].lastAcceptedFullTrajectoryblength = self.paths[i].fullTrajectoryblength
                        self.paths[i+ki].lastAcceptedFullTrajectoryflength = self.paths[i].fullTrajectoryflength
                report.append(str(accepted))
                report.append(' ')
                report.append(str(pathtype))
                report.append(' ')
            elif not accepted:
                report.append(str(accepted))
                report.append(' ')
                report.append(str(pathtype))
                report.append(' ')
                report.append(rejectreason)
        self.log.log.info(" ".join(report))
        
                    
#
#    def updateCrossingHistos(self):
#        for i in self.kernels.kernelPaths:
#            maxv = self.paths[i].getMaxTrajectoryValue()
#            self.paths[i].crossingHisto.addToHisto(maxv)
#    
#    def outputAllCrossingHistos(self):
#        for i in self.kernels.kernelPaths:
#            filename = os.path.join(self.paths[i].nfsladir,"histo.txt")
#            self.paths[i].crossingHisto.outputCrossingHisto(filename)
    
    def finalizeTIS(self,dirstring,newdirstring):
        self.finalizeShooting(dirstring)
        for i in self.kernels.kernelPaths:
            self.finalizeCopyLastAccepted(dirstring, newdirstring, i)
    
    
    
    def finalizeCopyLastAccepted(self,dirstring,newdirstring,pathnumber):
        
        def _copyFiles(fromdir,destdir):
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
             for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)       
            
#	    self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
	    #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#	    self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
            
        
        lapath = self.paths[pathnumber].nfsladir
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        
        #allpath = self.paths[pathnumber].nfsalldir
        #newalldir = os.path.join(allpath,newdirstring)
        #newalldirf = os.path.join(newalldir,"forward")
        #newalldirb = os.path.join(newalldir,"backward")
        dirlist = [lapath,newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)
        
        fdir = os.path.join(self.paths[pathnumber].workdir,"paths",dirstring)
        bdir = os.path.join(self.paths[pathnumber+1].workdir,"paths",dirstring)
        
#        _copyFiles(fdir, newalldirf)
#        _copyFiles(bdir, newalldirb)
#        self.filesystem.deleteFile(os.path.join(newalldirf,"traj.trr"))
#        self.filesystem.deleteFile(os.path.join(newalldirb,"traj.trr"))
        
        
        if self.paths[pathnumber].tisaccepted:
            _copyFiles(fdir, newladirf)
            _copyFiles(bdir, newladirb)
             
#            self.paths[pathnumber].lastAcceptedFullTrajectory = self.paths[pathnumber].fullTrajectory[:]
#            self.paths[pathnumber].lastAcceptedFullTrajectoryblength = self.paths[pathnumber].fullTrajectoryblength
#            self.paths[pathnumber].lastAcceptedFullTrajectoryflength = self.paths[pathnumber].fullTrajectoryflength           
#            
        else:
            oldladir = os.path.join(lapath,dirstring)
            fdir = os.path.join(oldladir,"forward")
            bdir = os.path.join(oldladir,"backward")
            
            _copyFiles(fdir, newladirf)
            _copyFiles(bdir, newladirb)
        

        
    
    
    def getFullTrajectory(self,pathnumber,dirstring):
        self.paths[pathnumber].fullTrajectory = []
        self.paths[pathnumber+1].fullTrajectory = []
        
        if dirstring == "workdir":
            fdir = self.paths[pathnumber].workdir
            bdir = self.paths[pathnumber+1].workdir
        else:
            fdir = os.path.join(self.paths[pathnumber].nfsladir,dirstring,"forward")
            bdir = os.path.join(self.paths[pathnumber+1].nfsladir,dirstring,"backward")

        if (self.paths[0].options.runoptions["parallel"]=='True'):
            parallel = True
        else:
            parallel = False
        
        for opi in range(len(self.orderparameters.op)):
            if parallel==True:
                traj,fl,bl = self.orderparameters.op[opi].getFullTrajectoryParallel(self.paths[pathnumber],self.paths[pathnumber+1],fdir,bdir,self.log,self.paths[0].options.runoptions["queuename"])
            else:
                traj,fl,bl = self.orderparameters.op[opi].getFullTrajectory(self.paths[pathnumber],self.paths[pathnumber+1],fdir,bdir,self.log)
            
            #cut the traj file and adjust path lengths
	    if optobj.runoptions["pathcut"]=="True":
                traj,fl,bl = self.cutAcceptedPaths(traj,fl,bl,opi)
                #now we have to cut the traj file. This is different for each
	        self.wrapper.cutTrajFile(fdir,bdir,fl,bl)
	        #end cutting

            self.paths[pathnumber].fullTrajectory.append(traj[:])
            self.paths[pathnumber+1].fullTrajectory.append(traj[:])
            self.paths[pathnumber].fullTrajectoryblength = bl
            self.paths[pathnumber+1].fullTrajectoryflength = fl
            #reset the length variable to new cut lengths
            self.paths[pathnumber].length = fl+bl
            self.paths[pathnumber+1].length = fl+bl
        
    def cutAcceptedPaths(self,traj,fl,bl,opi):
        """
        SRM : Cut the paths once it reaches the stable state.
        """
        trajfor = [pos for pos in traj if pos[2]==1]
        trajbak = [pos for pos in traj if pos[2]==0]
        #find the stable state location for the forward part
        sstatepos = self.findStableStateCrosssing(trajfor,opi)
        ##trim the traj length
        #logger.info(('forward old length %d new length %d')%(fl,fl-len(trajfor[sstatepos:])))
        fl = fl - len(trajfor[sstatepos:])
        #trim the trajectory
        trajfor = trajfor[:sstatepos]
        #for rev part, reverse first
        trajbakrev = trajbak[::-1]
        sstatepos = self.findStableStateCrosssing(trajbakrev,opi)
        #logger.info(('backward old length %d new length %d')%(bl,bl-len(trajbakrev[sstatepos:])))
        ##trim the traj length
        bl = bl - len(trajbakrev[sstatepos:])
        #trim the trajectory
        trajbakrev= trajbakrev[:sstatepos]
        #now reverse again to make it normal
        trajbak = trajbakrev[::-1]
        traj = trajbak + trajfor
        return traj,fl,bl
        

    def findStableStateCrosssing(self,qtraj,opi):
        """
        SRM: find the location where the traj crosses stable state
        """
        sstate = -1
        #logger.info(self.orderparameters.op[opi].A[1])
        #logger.info(self.orderparameters.op[opi].B[0])

        for i in range(len(qtraj)):
            if ((qtraj[i][1] < self.orderparameters.op[opi].A[1]) or (qtraj[i][1] > self.orderparameters.op[opi].B[0])):
                sstate = i+1
                break
        #logger.info(sstate)
        #added to handle the extra condition
        if sstate==-1:
            sstate=len(qtraj)

        return sstate




    def _outputFullTrajectory(self,directory,pathnumber,dirstring):
        for opi in range(len(self.orderparameters.op)):    
            of = open(os.path.join(directory,"trajectory."+"%02d"% (opi,)+"."+dirstring+".dat"),"w")
            for position in self.paths[pathnumber].lastAcceptedFullTrajectory[opi]:
                of.write("%d %.18f %d\n" % (position[0],position[1],position[2]))
            of.close()
            
            of = open(os.path.join(directory,"trajectory."+"%02d"% (opi,)+"."+dirstring+".trial"),"w")
            for position in self.paths[pathnumber].fullTrajectory[opi]:
                of.write("%d %.18f %d\n" % (position[0],position[1],position[2]))
            of.close()
     
    
    
    def outputAllFullTrajectories(self,dirstring):
        for i in self.kernels.kernelPaths:
            ladir = os.path.join(self.paths[i].nfsladir,dirstring)
            self._outputFullTrajectory(ladir, i,dirstring)
            
    

    def readLastAcceptedTrajectories(self,dirstring):
        for i in self.kernels.kernelPaths:
            ladir = os.path.join(self.paths[i].nfsladir,dirstring)
            self._readFullTrajectory(ladir, i,dirstring)
        
    
    def _readFullTrajectory(self,directory,pathnumber,dirstring):
        traj = []
        self.paths[pathnumber].lastAcceptedFullTrajectoryblength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectoryflength = 0
        self.paths[pathnumber].lastAcceptedFullTrajectory = []
        self.paths[pathnumber+1].lastAcceptedFullTrajectory = []
        
        for opi in range(len(self.orderparameters.op)):    
            for line in open(os.path.join(directory,"trajectory."+"%02d"% (opi,)+"."+dirstring+".dat"),"r"):
                raw = line.split()
                traj.append([int(float(raw[0])),float(raw[1]),int(float(raw[2]))])
                if int(float(raw[2])) == 0:
                    self.paths[pathnumber].lastAcceptedFullTrajectoryblength += 1
                else:
                    self.paths[pathnumber].lastAcceptedFullTrajectoryflength += 1
            
            for ki in range(2):
                self.paths[pathnumber+ki].lastAcceptedFullTrajectory.append(traj[:])
            self.paths[pathnumber+1].lastAcceptedFullTrajectoryflength =self.paths[pathnumber].lastAcceptedFullTrajectoryflength
            self.paths[pathnumber+1].lastAcceptedFullTrajectoryblength =self.paths[pathnumber].lastAcceptedFullTrajectoryblength
            
        
    def _copyLastAcceptedToFull(self,pathnumber):
        self.paths[pathnumber].fullTrajectoryblength = self.paths[pathnumber].lastAcceptedFullTrajectoryblength
        self.paths[pathnumber].fullTrajectoryflength = self.paths[pathnumber].lastAcceptedFullTrajectoryflength
        self.paths[pathnumber].fullTrajectory = []
        for opi in range(len(self.orderparameters.op)): 
            self.paths[pathnumber].fullTrajectory.append(self.paths[pathnumber].lastAcceptedFullTrajectory[opi][:])
        
        
        
 
    """ 
    ************************************************************************
    
        Replica Exchange
    
    ************************************************************************
    """
    
    
    
    def trialReplicaExchange(self,dirstring,newdirstring):
        possiblePairs = [i*2 for i in range(self.npaths)]
        exchange = []
        for i in range(self.npaths*(self.npaths-1)):
            random.shuffle(possiblePairs)
            if self.paths[possiblePairs[0]].forward == self.paths[possiblePairs[1]].forward:
                if self.paths[possiblePairs[0]].forward:
                    interfaces = [self.paths[possiblePairs[0]].interface,self.paths[possiblePairs[1]].interface]
                    maxi = max(interfaces)
                    if self.paths[possiblePairs[0]].getMaxTrajectoryValue() > maxi and self.paths[possiblePairs[1]].getMaxTrajectoryValue() > maxi:
                        exchange = possiblePairs[0:2]
                else:
                    interfaces = [self.paths[possiblePairs[0]].interface,self.paths[possiblePairs[1]].interface]
                    mini = min(interfaces)
                    if self.paths[possiblePairs[0]].getMaxTrajectoryValue() < mini and self.paths[possiblePairs[1]].getMaxTrajectoryValue() < mini:
                        exchange = possiblePairs[0:2]
        
        self.finalizeExchange(dirstring, newdirstring, exchange)
        return exchange
        
        
    def finalizeExchange(self,dirstring,newdirstring,exchange):
        base = [[i*2,i*2] for i in range(self.npaths)]
        if len(exchange) == 2:
            for i,b in enumerate(base):
                if b[1] == exchange[0]:
                    base[i][1] = exchange[1]
                elif b[1] == exchange[1]:
                    base[i][1] = exchange[0]
        for pair in base:
            self.finalizeCopyExchange(dirstring,newdirstring, pair)
        print base
        
    def finalizeCopyExchange(self,dirstring,newdirstring,pair):
        
        def _copyFiles(fromdir,destdir):
            for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)

#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
        
        lapath = self.paths[pair[1]].nfsladir
        #allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        #newalldir = os.path.join(allpath,newdirstring)
        #newalldirf = os.path.join(newalldir,"forward")
        #newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)

        oldlapath = self.paths[pair[0]].nfsladir
        
        fdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
        bdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
        
        _copyFiles(fdir, newladirf)
        _copyFiles(bdir, newladirb)
#        _copyFiles(fdir, newalldirf)
#        _copyFiles(bdir, newalldirb)
        
    """ 
    ************************************************************************
    
        Fast Replica Exchange
    
    ************************************************************************
    """
    
    
    
    def trialReplicaExchangeFast(self,dirstring,newdirstring):
        possiblePairs = [i*2 for i in range(self.npaths)]
        exchange = []
	ex1 = random.randint(0,len(possiblePairs)-2)
	ex2 = ex1 + 1
        if self.paths[possiblePairs[ex1]].forward == self.paths[possiblePairs[ex2]].forward:
 		if self.paths[possiblePairs[ex1]].forward:
                    interfaces = [self.paths[possiblePairs[ex1]].interface,self.paths[possiblePairs[ex2]].interface]
                    maxi = max(interfaces)
                    if self.paths[possiblePairs[ex1]].getMaxTrajectoryValue() > maxi and self.paths[possiblePairs[ex2]].getMaxTrajectoryValue() > maxi:
                        exchange.append(possiblePairs[ex1])
                        exchange.append(possiblePairs[ex2])
            	else:
                    interfaces = [self.paths[possiblePairs[ex1]].interface,self.paths[possiblePairs[ex2]].interface]
                    mini = min(interfaces)
                    if self.paths[possiblePairs[ex1]].getMaxTrajectoryValue() < mini and self.paths[possiblePairs[ex2]].getMaxTrajectoryValue() < mini:
                        exchange.append(possiblePairs[ex1])
                        exchange.append(possiblePairs[ex2])
        
        self.finalizeExchangeFast(dirstring, newdirstring, exchange)
        return exchange
        
        
    def finalizeExchangeFast(self,dirstring,newdirstring,exchange):
        base = [[i*2,i*2] for i in range(self.npaths)]
        if len(exchange) == 2:
            for i,b in enumerate(base):
                if b[1] == exchange[0]:
                    base[i][1] = exchange[1]
                elif b[1] == exchange[1]:
                    base[i][1] = exchange[0]
        for pair in base:
            self.finalizeCopyExchangeFast(dirstring,newdirstring, pair)
        print base
        
    def finalizeCopyExchangeFast(self,dirstring,newdirstring,pair):
        
        def _copyFiles(fromdir,destdir):
            for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)


        def _copyQTrajectory(fromdir,destdir,dirstring,newdirstring):
            fileold = os.path.join(fromdir,"trajectory.00."+dirstring+".dat") 
            fileoldtrial = os.path.join(fromdir,"trajectory.00."+dirstring+".trial") 
            filenew = os.path.join(destdir,"trajectory.00."+newdirstring+".dat")
            filenewtrial = os.path.join(destdir,"trajectory.00."+newdirstring+".trial")
            
	    self.filesystem.copyFiles(fileold, filenew)
	    self.filesystem.copyFiles(fileoldtrial, filenewtrial)


        lapath = self.paths[pair[1]].nfsladir
        #allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        #newalldir = os.path.join(allpath,newdirstring)
        #newalldirf = os.path.join(newalldir,"forward")
        #newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)

        oldlapath = self.paths[pair[0]].nfsladir
        oldladir = os.path.join(oldlapath,dirstring)
        
        fdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
        bdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
        
        _copyFiles(fdir, newladirf)
        _copyFiles(bdir, newladirb)
#        _copyFiles(fdir, newalldirf)
#        _copyFiles(bdir, newalldirb)
	_copyQTrajectory(oldladir,newladir,dirstring,newdirstring)
        
        
    
    """ 
    ************************************************************************
    
        Backward/Forward Exchange
    
    ************************************************************************
    """
    
    def trialBackwardForwardExchange(self,dirstring,newdirstring):
        possiblePairs = [i*2 for i in range(self.npaths)]
        exchange = []
        for i in range(self.npaths*(self.npaths-1)):
            random.shuffle(possiblePairs)
            
            if self.paths[possiblePairs[0]].forward != self.paths[possiblePairs[1]].forward:
                if (self.paths[possiblePairs[0]].checkAcceptedTISAB(self.orderparameters.op)) and \
                   (self.paths[possiblePairs[1]].checkAcceptedTISAB(self.orderparameters.op)):
                    exchange = possiblePairs[0:2]
                    break

        self.finalizeBFExchange(dirstring, newdirstring, exchange)
        return exchange
    
    def finalizeBFExchange(self,dirstring,newdirstring,exchange):
        base = [[i*2,i*2] for i in range(self.npaths)]
        if len(exchange) == 2:
            for i,b in enumerate(base):
                if b[1] == exchange[0]:
                    base[i][1] = exchange[1]
                elif b[1] == exchange[1]:
                    base[i][1] = exchange[0]
        for pair in base:
            self.finalizeBFCopyExchange(dirstring,newdirstring, pair)
        print base
        
    def finalizeBFCopyExchange(self,dirstring,newdirstring,pair):    
        def _copyFiles(fromdir,destdir):
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
             for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)
        
        lapath = self.paths[pair[1]].nfsladir
        #allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        #newalldir = os.path.join(allpath,newdirstring)
        #newalldirf = os.path.join(newalldir,"forward")
        #newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)

        oldlapath = self.paths[pair[0]].nfsladir
        if pair[0] == pair[1]:
            fdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
            bdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
        else:
            bdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
            fdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
            
        _copyFiles(fdir, newladirf)
        _copyFiles(bdir, newladirb)
#        _copyFiles(fdir, newalldirf)
#        _copyFiles(bdir, newalldirb)
        
    
    """ 
    ************************************************************************
    
        Fast Backward/Forward Exchange
    
    ************************************************************************
    """
    
    def trialBackwardForwardExchangeFast(self,dirstring,newdirstring):
        possiblePairs = [i*2 for i in range(self.npaths)]
        
	for i in range(self.npaths*(self.npaths-1)):
		ex1 = random.randint(0,len(possiblePairs)-2)
		ex2 = ex1 + 1
        	exchange = []
            
        	if self.paths[possiblePairs[ex1]].forward != self.paths[possiblePairs[ex2]].forward:
        		if (self.paths[possiblePairs[ex1]].checkAcceptedTISAB(self.orderparameters.op)) and \
                   	(self.paths[possiblePairs[ex2]].checkAcceptedTISAB(self.orderparameters.op)):
                    		exchange.append(possiblePairs[ex1])
                    		exchange.append(possiblePairs[ex2])
                    		break

        self.finalizeBFExchangeFast(dirstring, newdirstring, exchange)
        return exchange
    
    def finalizeBFExchangeFast(self,dirstring,newdirstring,exchange):
        base = [[i*2,i*2] for i in range(self.npaths)]
        if len(exchange) == 2:
            for i,b in enumerate(base):
                if b[1] == exchange[0]:
                    base[i][1] = exchange[1]
                elif b[1] == exchange[1]:
                    base[i][1] = exchange[0]
        for pair in base:
            self.finalizeBFCopyExchangeFast(dirstring,newdirstring, pair)
        print base
        
    def finalizeBFCopyExchangeFast(self,dirstring,newdirstring,pair):    
        def _copyFiles(fromdir,destdir):
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
             for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)            



        def _copyQTrajectory(fromdir,destdir,dirstring,newdirstring):
            fileold = os.path.join(fromdir,"trajectory.00."+dirstring+".dat") 
            fileoldtrial = os.path.join(fromdir,"trajectory.00."+dirstring+".trial") 
            filenew = os.path.join(destdir,"trajectory.00."+newdirstring+".dat")
            filenewtrial = os.path.join(destdir,"trajectory.00."+newdirstring+".trial")
            
	    self.filesystem.copyFiles(fileold, filenew)
	    self.filesystem.copyFiles(fileoldtrial, filenewtrial)

	    
        lapath = self.paths[pair[1]].nfsladir
        #allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        #newalldir = os.path.join(allpath,newdirstring)
        #newalldirf = os.path.join(newalldir,"forward")
        #newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)

        oldlapath = self.paths[pair[0]].nfsladir
        oldladir = os.path.join(oldlapath,dirstring)
        if pair[0] == pair[1]:
            fdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
            bdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
        else:
            bdir = os.path.join(os.path.join(oldlapath,dirstring),"forward")
            fdir = os.path.join(os.path.join(oldlapath,dirstring),"backward")
            
        _copyFiles(fdir, newladirf)
        _copyFiles(bdir, newladirb)
#        _copyFiles(fdir, newalldirf)
#        _copyFiles(bdir, newalldirb)

	if pair[0] == pair[1]:
		_copyQTrajectory(oldladir,newladir,dirstring,newdirstring)
	else:
		self.getFullTrajectory(pair[1],newdirstring)
		self.paths[pair[1]].lastAcceptedFullTrajectory = []
		for opi in range(len(self.orderparameters.op)):
			self.paths[pair[1]].lastAcceptedFullTrajectory.append(self.paths[pair[1]].fullTrajectory[opi][:])
		self._outputFullTrajectory(newladir,pair[1],newdirstring)
	
        
        
        
    """ 
    ************************************************************************
    
        ReversePath
    
    ************************************************************************
    """
        
        
    
    def preperationReverse(self,copyfiles = True):
        """
        Creating structure for the reversePaths
        """            
        for rpath in self.reversePaths:
            initdir = rpath.options.paths["initialpath"]
            workdir = rpath.workdir
            dirlist = [rpath.srcatchbase,rpath.baseworkdir,rpath.workdir,os.path.join(workdir,'paths')]
            #dirlist = [rpath.srcatchbase,rpath.workdir,os.path.join(workdir,'paths')]
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            #dirlist = [os.path.join(initdir,'la'),rpath.nfsbaseworkdir,rpath.nfsladir,rpath.nfsconfigstore]
            dirlist = [os.path.join(initdir,'la'),rpath.nfsladir,rpath.nfsconfigstore]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                #for myfile in rpath.options.standardfiles:
                for myfile in rpath.options.standardfiles.itervalues():
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",myfile), workdir)

	    if rpath.options.runoptions["mddriver"] in ('lammps','lammps_parallel'):
	    	self.wrapper.modifyMDseed(workdir,rpath.options)

                    
            
    def shootingQueueReverse(self):
        """
        Function that performes the mdrun on each of the reverse paths.
        """
        
#      JR: not needed with lammps
#	"""
#        First step is to execute grompp in each of the paths workdirs to get a topol.tpr
#        """
#        for rpaths in self.reversePaths:
#            workdir = rpaths.workdir
#            os.chdir(workdir)
#            ok = self.wrapper.checkInputFiles(rpaths.options, rpaths.options.initoptions,                         \
#                                         ["-c","-p","-n"],workdir)
#            self.log.log.debug(workdir + str(ok) + " grompp executed")
#            rpaths.options.writeMdpFile(workdir,"md.mdp")
#            
#            cmd = self.wrapper.generateCommand(rpaths.options, rpaths.options.initoptions,                        \
#                                          ["-c","-o","-p","-f","-n"],                                  \
#                                          workdir,"grompp" )
#            self.wrapper.executeCommand(cmd)
#            os.chdir(self.basedir)
            
        """
        A process list is generated. It holds all the threads in a dictionary.
        The todo-list is the list of all processes that have to be exsecuted.
        """
        processlist = {}
        todolist = range(len(self.reversePaths))
        
        for i in  range(min(len(self.reversePaths),self.cores-1)):
            nextp = todolist.pop()
            os.chdir(self.reversePaths[nextp].workdir)
            cmd =self.wrapper.generateMDRunCommand(self.reversePaths[nextp].options, self.reversePaths[nextp].workdir)
            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
            
        finished = False
        """
        loop until all processes in the todo list are finished
        """
        while not(finished):   
            time.sleep(float(optobj.runoptions["updatetime"]))
            finished = True
            for i in processlist.keys():
                #poll=None means it is still running
                if processlist[i].poll() == None:
                    finished = False
                    os.chdir(self.reversePaths[i].workdir)

		    q = []						#JR set of order paramters of last slice
		    for opi in range(len(self.orderparameters.op)):
#			print 'Checking the order parameter\n'		#JR
#			sys.stdout.flush()   				#JR
			rlength,qval = self.orderparameters.op[opi].QofLast(self.reversePaths[i]) 
                        q.append(qval)
	

                    #self.log.log.info(str(i) + " " + str(processlist[i].poll()) + " " + str(len(self.distparser.data)) + " " +  str(self.distparser.data[-1]))
#		    pathfinished = self.reversePaths[i].checkFinishedTIS(float(self.distparser.data[-1][self.stablestates.gdistDirection]),self.stablestates.states)
		    #JR only checking for the first order parameter right now!!!
	            if rlength == 0:
		    	pathfinished = False
		    else:
                    	pathfinished = self.reversePaths[i].checkFinishedTIS(q[0])
	            print pathfinished,rlength,q[0]		#JR print for testing
                    if pathfinished or rlength > self.reversePaths[i].options.runoptions["maxlength"]:
                        self.log.log.debug("Finished " + str(i) + " " + str(self.reversePaths[i].finishedState))
                        processlist[i].terminate()
                        if len(todolist)>0:
                            nextp = todolist.pop()
                            os.chdir(self.reversePaths[nextp].workdir)
                            cmd =self.wrapper.generateMDRunCommand(self.reversePaths[nextp].options, self.reversePaths[nextp].workdir)
                            processlist[nextp] = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                            self.log.log.debug(str(todolist) + " " +  str(nextp) + "th thread started")
        self.log.log.debug("Shooting is finished")
     
     
    def finalizeShootingReverse(self,dirstring):
        for rpath in self.reversePaths:
            workdir = rpath.workdir 
            dest = os.path.join(workdir,"paths",dirstring)
            self.filesystem.createDir(dest)
            
#            self.filesystem.moveFile(os.path.join(workdir,"traj.xtc"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"traj.trr"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"topol.tpr"), dest)
             #self.filesystem.moveFile(os.path.join(workdir,"traj.dat"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.dump"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"md.in"), dest)
             #self.filesystem.deleteFileList(workdir,"#*")
#            self.filesystem.deleteFileList(workdir,"ener.edr")
            for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(workdir,file), dest)
            self.filesystem.deleteFileList(workdir,"#*")
            for file in optobj.deletefiles.itervalues():
                self.filesystem.deleteFileList(workdir,file)
    
            
    def finalizeReverse(self,dirstring,newdirstring):
        self.finalizeShootingReverse(dirstring)
        for i in range(len(self.reversePaths)):
            self.finalizeCopyLastAcceptedReverse(dirstring, newdirstring, i)
    
    
    
    def finalizeCopyLastAcceptedReverse(self,dirstring,newdirstring,pathnumber):
        
        def _copyFiles(fromdir,destdir):
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            #self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
            for file in optobj.copyfiles3.itervalues():
                self.filesystem.copyFiles(os.path.join(fromdir,file), destdir)
        
        lapath = self.reversePaths[pathnumber].nfsladir
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        
        dirlist = [lapath,newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)
        
        fdir = os.path.join(self.reversePaths[pathnumber].workdir,"paths",dirstring)
        
        _copyFiles(fdir, newladirf)
        
        

    def getFullTrajectoryReverse(self,pathnumber,dirstring):
        
        if dirstring == "workdir":
            fdir = self.reversePaths[pathnumber].workdir
        else:
            fdir = os.path.join(self.reversePaths[pathnumber].nfsladir,dirstring,"forward")
        
        traj = []
        count = 0
        os.chdir(fdir)

        for opi in range(len(self.orderparameters.op)):
            traj,fl = self.orderparameters.op[opi].getFullTrajectoryReverse(self.reversePaths[pathnumber],fdir,self.log)
            self.reversePaths[pathnumber].fullTrajectory.append(traj[:])
            self.reversePaths[pathnumber].fullTrajectoryflength = fl		#JR: is this okay if there are more than one order parameter?


    def _outputFullTrajectoryReverse(self,directory,pathnumber,dirstring):
	for opi in range(len(self.orderparameters.op)):    
            of = open(os.path.join(directory,"trajectoryZ."+"%02d"% (opi,)+"."+dirstring+".acc"),"w")
            for position in self.reversePaths[pathnumber].lastAcceptedFullTrajectory[opi]:
                of.write("%d %.18f %d\n" % (position[0],position[1],position[2]))
            of.close()
            
    
    def outputAllFullTrajectoriesReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            ladir = os.path.join(self.reversePaths[i].nfsladir,dirstring)
            self._outputFullTrajectoryReverse(ladir, i,dirstring)
            

    def lastReverseToConf(self,dirstring):
        """
        The last reverse trajectory is transformed into pathXX.dump files 
        """
        for i in range(len(self.reversePaths)):
            mpath = self.reversePaths[i]
            workdir = os.path.join(mpath.nfsladir,dirstring,"forward")
            latemp = mpath.latemp								#JR: for 'normal' path this is done in tmp directory...
            for file in optobj.copyfiles4.itervalues():
                self.filesystem.copyFiles(os.path.join(workdir,file) , latemp)
            #self.filesystem.copyFiles(os.path.join(workdir,mpath.options.initoptions["pathfile"][1]) , latemp)
            #self.wrapper.TrajToConf(mpath.options,latemp)
            self.wrapper.TrajToConf(mpath.options,workdir)					#JR: here now in 'regular' directory


    
    def copyStateABFiles(self,dirstring):

        for i in range(len(self.reversePaths)):
            slist = self.filesystem.getFileList(self.reversePaths[i].nfsladir,"path*")
            ns=len(slist)
	    #JR:  in principle there would have to be a loop over the order parameters here!
	    #JR:  only taking the first one right now!
            traj = self.reversePaths[i].fullTrajectory[0]
            if self.reversePaths[i].state == 0:
                for point in traj:
#                    if point[1] > self.stablestates.states[0]:
                    if point[1] > self.reversePaths[i].interface:
                        ns += 1
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward",("path%d"+optobj.extension) % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,("path%d"+optobj.extension) % (ns,))
                        self.filesystem.copyFiles(source, destination)
            else:
                for point in traj:
#                    if point[1] < self.stablestates.states[1]:
                    if point[1] < self.reversePaths[i].interface:
                        ns += 1
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward",("path%d"+optobj.extension) % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,("path%d"+optobj.extension) % (ns,))
                        self.filesystem.copyFiles(source, destination)
            
#            rn = self.reversePaths[i].fullTrajectory[-1][0]
            rn = traj[-1][0]

            source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward",("path%d"+optobj.extension) % rn)
            destination = os.path.join(self.reversePaths[i].nfsladir,dirstring,("conf"+optobj.extension))
            self.filesystem.copyFiles(source, destination)
    
    def pickConfigurationReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            filename = os.path.join(self.reversePaths[i].nfsladir,dirstring,("conf"+optobj.extension))
            self.filesystem.copyFiles(filename, self.reversePaths[i].workdir)

    def writeFinishedFilesReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            lapath = self.reversePaths[i].nfsladir
            newladir = os.path.join(lapath,dirstring)
            filename = os.path.join(newladir,"finished")
            of = open(filename,"w")
            of.write("finished")
            of.close()


    def checkTisPathsAcceptedReverse(self):
        """
        This function updates the tisaccepted variable
        """
        
        report = []
	for i in range(len(self.reversePaths)):
            self.getFullTrajectoryReverse(i,dirstring = "workdir")
           

            accepted,rejectreason = self.reversePaths[i].checkAcceptedTIS(self.orderparameters.op)
            self.reversePaths[i].tisaccepted = accepted
            
            if accepted or self.mode=="initial":
                self.reversePaths[i].lastAcceptedFullTrajectory = []
                for opi in range(len(self.orderparameters.op)):
                    self.reversePaths[i].lastAcceptedFullTrajectory.append(self.reversePaths[i].fullTrajectory[opi][:])
                    self.reversePaths[i].lastAcceptedFullTrajectoryflength = self.reversePaths[i].fullTrajectoryflength
                report.append(str(accepted))
            elif not accepted:
                report.append(str(accepted))
                report.append(' ')
                report.append(rejectreason)

        self.log.log.info(" ".join(report))
    
            
