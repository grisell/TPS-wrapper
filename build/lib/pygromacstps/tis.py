
'''
Created on May 3, 2010

@author: Wolfgang Lechner
'''


import wrappers
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
#python stuff
import random
import multiprocessing
import os
import subprocess
from math import fabs

import time
import sys

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
        self.cores = multiprocessing.cpu_count()
        self.wrapper = wrappers.lammpswrapper()
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

	print self.orderparameters.op[0].A
	print self.orderparameters.op[0].B

	sys.exit(1)
        
        
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
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            dirlist = [os.path.join(initdir,'la'),self.paths[i].nfsbaseworkdir,self.paths[i].nfsladir,os.path.join(self.paths[i].nfsladir,"%07d" % (0))]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                for file in self.paths[i].options.standardfiles:
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",file), workdir)
            
            
    
    def shootingInitialGroFiles(self):
        """
        Perform a shooting move on a gro file.
        """
        for i in self.kernels.kernelPaths:
            conffile = os.path.join(self.paths[i].workdir,"conf.gro")
            bakfile = os.path.join(self.paths[i].workdir,"bakconf.gro")
            bconffile = os.path.join(self.paths[i].workdir,"bconf.gro")
            self.helper.shootVelocities(conffile, bconffile)  
            self.filesystem.moveFile(conffile, bakfile)
            self.filesystem.moveFile(bconffile, conffile)
            dest = self.paths[i+1].workdir
            self.filesystem.copyFiles(conffile, dest)
            self.reverseBackwardGroFile(i+1)
            
            
            
                         
    def finalizeInitial(self):
        
        def _copyFiles(fromdir,destdir):
            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
        
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
        
        for i in range(self.reversePaths):
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
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            dirlist = [os.path.join(initdir,'all'),os.path.join(initdir,'la'),self.paths[i].nfsbaseworkdir,self.paths[i].nfsladir]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                for file in self.paths[i].options.standardfiles:
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",file), workdir)
            
    

    """ 
    ************************************************************************
    
        Shooting procedure
    
    ************************************************************************
    """
    
    def lastAcceptedToGro(self,dirnumber):
        """
        The last accepted trajectory is transformed into gro files using trjconv. 
        """
        for i in self.kernels.kernelPathsAll:
            
            
            mpath = self.paths[i]
            fbdir = mpath.forwardpart and "forward" or "backward"
            nfsworkdir = os.path.join(mpath.nfsladir,dirnumber,fbdir)
            latemp = mpath.latemp
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"traj.trr") , latemp)
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"index.ndx") , latemp)
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"topol.tpr") , latemp)
            cmd = self.wrapper.generateTrjconvCommand(mpath.options,latemp,latemp)
            self.wrapper.executeCommand(cmd, tinput="0")
            
    
    def pickConfigurationsTIS(self,dirnumber):
        for i in self.kernels.kernelPaths:
            self.pickConfigurationLastAccepted(dirnumber,i)  
    
    def pickConfigurationLastAccepted(self,dirnumber,pathnumber):
        """
        A random configuration is taken from a trajectory and the shooting move is performed.
        and copy the files to the workdir of the path.
        """
        
        fdir = os.path.join(self.paths[pathnumber].latemp)
        flist = self.filesystem.getFileList(fdir, "path*.gro")
        fsize = len(flist)
        bdir = os.path.join(self.paths[pathnumber+1].latemp)
        blist = self.filesystem.getFileList(bdir, "path*.gro")
        bsize = len(blist)
        
        sum = fsize + bsize
        print "sum " ,sum, fsize, bsize
        print "traj ", len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi])
#        if sum != len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi]):
#            self.log.log.info(str(pathnumber) + " " + str(sum) + " " + str( len(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi])))
        c = 0
        if self.paths[pathnumber].options.runoptions["shootfrominterface"] == 0:
            rn = sum/2
            while True:
                rn = random.randint(0,sum-2)
                if self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] > self.orderparameters.op[distopi].A[1] and self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] < self.orderparameters.op[distopi].B[0]:
                    break
                c += 1
                if c > 100:
                    print "ERROR picking shootingpoint"
                    break
            
        else:
            mindist = 1000000.0
            rn = sum/2
            for kr in range(0,sum-2):
                act = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][kr][1]
                distance = fabs(act - self.paths[pathnumber].interface)
                accept = (bool(self.paths[pathnumber].forward) == bool(act > self.paths[pathnumber].interface))
                #self.log.log.debug("ptest " + " "+ str(accept)  + " "  + str(act)  + " "  + str(self.paths[pathnumber].interface) ) 
                if distance < mindist and accept:
                    mindist = distance
                    rn = kr
            self.log.log.debug("Picked " +str(rn) + " " + str(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1])+" " + str(self.paths[pathnumber].interface) + "  "+ str(mindist)+  " " + str(pathnumber) )
        
        backwardpart = (rn < bsize)
        if backwardpart:
            file = os.path.join(bdir,"path%d.gro" % (bsize-rn))
        else:
            file = os.path.join(fdir,"path%d.gro" % (rn-bsize))
        
        #self.log.log.debug("pick " +dirnumber + " " + str(pathnumber) +" " + str(rn) +  " " + str( sum) +  " " + str(fsize) + " " +  str(bsize) + " " +  str(len(self.paths[pathnumber].lastAcceptedFullTrajectory)) + " " + file)
       
        self.paths[pathnumber].shootingTime = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][0]
        
        destfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"conforig.gro")
        shootfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"confshoot.gro")
        
        if os.path.exists(file):    
            self.filesystem.copyFiles(file,destfile)
        else:
            self.log.log.ERROR("File " + file + " not found, this is a fatal error")
        
        self.helper.shootVelocities(destfile,shootfile)
        
        for path in [self.paths[pathnumber],self.paths[pathnumber+1]]:
            self.filesystem.copyFiles(shootfile, os.path.join(path.workdir,"conf.gro"))
            path.finishedstate = -1
        
        if backwardpart:
            reversepath = pathnumber 
        else:
            reversepath = pathnumber + 1
        self.reverseBackwardGroFile(reversepath)
        
        self.log.log.debug("shooting from " +str(rn)+ " file " + file + " backwardpart  " + str(backwardpart))
        
        self.filesystem.deleteFileList(fdir, "path*.gro")
        self.filesystem.deleteFileList(bdir, "path*.gro")
        
        
        
    def reverseBackwardGroFile(self,pathnumber):
        """
        Reverse all velocities from a gro file.
        """
        conffile = os.path.join(self.paths[pathnumber].workdir,"conf.gro")
        bakfile = os.path.join(self.paths[pathnumber].workdir,"bakconf.gro")
        bconffile = os.path.join(self.paths[pathnumber].workdir,"bconf.gro")
        self.helper.reverseVelocities(conffile, bconffile)
        self.filesystem.moveFile(conffile, bakfile)
        self.filesystem.moveFile(bconffile, conffile)


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
            self.gromppPath(self.paths[i])
        
        # initialize the reverse Paths
        for rpath in self.reversePaths:
            self.gromppPath(rpath)
        
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
            self.gromppPath(self.paths[i])
        
        # initialize the reverse Paths
        for rpath in self.reversePaths:
            self.gromppPath(rpath)
        
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
            
            self.filesystem.moveFile(os.path.join(workdir,"traj.xtc"), dest)
            self.filesystem.moveFile(os.path.join(workdir,"traj.trr"), dest)
            self.filesystem.moveFile(os.path.join(workdir,"topol.tpr"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
            self.filesystem.deleteFileList(workdir,"#*")
            self.filesystem.deleteFileList(workdir,"ener.edr")

    def deleteOldTrrFiles(self,dirnumber):
        dirlist = []
        for i in self.kernels.kernelPaths:
            dirlist.append(os.path.join(self.paths[i].nfsalldir,dirnumber,"forward"))
            dirlist.append(os.path.join(self.paths[i].nfsalldir,dirnumber,"backward"))
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
            
            
            accepted = self.paths[i].checkAcceptedTIS(self.orderparameters.op,self.log)
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
            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
            
        
        lapath = self.paths[pathnumber].nfsladir
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        
        allpath = self.paths[pathnumber].nfsalldir
        newalldir = os.path.join(allpath,newdirstring)
        newalldirf = os.path.join(newalldir,"forward")
        newalldirb = os.path.join(newalldir,"backward")
        dirlist = [allpath,lapath,newladir,newladirf,newladirb,newalldir,newalldirf,newalldirb]
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
        
        for opi in range(len(self.orderparameters.op)):
            traj,fl,bl = self.orderparameters.op[opi].getFullTrajectory(self.paths[pathnumber],self.paths[pathnumber+1],fdir,bdir,self.log)
            self.paths[pathnumber].fullTrajectory.append(traj[:])
            self.paths[pathnumber+1].fullTrajectory.append(traj[:])
            self.paths[pathnumber].fullTrajectoryblength = bl
            self.paths[pathnumber+1].fullTrajectoryflength = fl
        
    
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
            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
        
        lapath = self.paths[pair[1]].nfsladir
        allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        newalldir = os.path.join(allpath,newdirstring)
        newalldirf = os.path.join(newalldir,"forward")
        newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb,newalldir,newalldirf,newalldirb]
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
            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
        
        lapath = self.paths[pair[1]].nfsladir
        allpath = self.paths[pair[1]].nfsalldir
        
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        newalldir = os.path.join(allpath,newdirstring)
        newalldirf = os.path.join(newalldir,"forward")
        newalldirb = os.path.join(newalldir,"backward")
        
        dirlist = [newladir,newladirf,newladirb,newalldir,newalldirf,newalldirb]
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
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            dirlist = [os.path.join(initdir,'la'),rpath.nfsbaseworkdir,rpath.nfsladir,rpath.nfsconfigstore]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                for file in rpath.options.standardfiles:
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",file), workdir)
                    
            
    def shootingQueueReverse(self):
        """
        Function that performes the mdrun on each of the reverse paths.
        """
        
        """
        First step is to execute grompp in each of the paths workdirs to get a topol.tpr
        """
        for rpaths in self.reversePaths:
            workdir = rpaths.workdir
            os.chdir(workdir)
            ok = self.wrapper.checkInputFiles(rpaths.options, rpaths.options.initoptions,                         \
                                         ["-c","-p","-n"],workdir)
            self.log.log.debug(workdir + str(ok) + " grompp executed")
            rpaths.options.writeMdpFile(workdir,"md.mdp")
            
            cmd = self.wrapper.generateCommand(rpaths.options, rpaths.options.initoptions,                        \
                                          ["-c","-o","-p","-f","-n"],                                  \
                                          workdir,"grompp" )
            self.wrapper.executeCommand(cmd)
            os.chdir(self.basedir)
            
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
            time.sleep(float(self.reversePaths[0].options.runoptions["updatetime"]))
            finished = True
            for i in processlist.keys():
                #poll=None means it is still running
                if processlist[i].poll() == None:
                    finished = False
                    os.chdir(self.reversePaths[i].workdir)
                    cmd = self.wrapper.generateGDistCommand(self.reversePaths[i].workdir,self.reversePaths[i].options)
                    self.wrapper.executeCommand(cmd, self.reversePaths[i].options.runoptions["dist1"]+"\n"+self.reversePaths[i].options.runoptions["dist2"]+"\n" )
                    self.distparser.readDist(os.path.join(self.reversePaths[i].workdir,"dist.xvg"))
                    self.filesystem.deleteFile(os.path.join(self.reversePaths[i].workdir,"dist.xvg"))
                    #self.log.log.info(str(i) + " " + str(processlist[i].poll()) + " " + str(len(self.distparser.data)) + " " +  str(self.distparser.data[-1]))
                    pathfinished = self.reversePaths[i].checkFinishedTIS(float(self.distparser.data[-1][self.stablestates.gdistDirection]),self.stablestates.states)
                    if pathfinished or len(self.distparser.data) > self.reversePaths[i].options.runoptions["maxlength"]:
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
            
            self.filesystem.moveFile(os.path.join(workdir,"traj.xtc"), dest)
            self.filesystem.moveFile(os.path.join(workdir,"traj.trr"), dest)
            self.filesystem.moveFile(os.path.join(workdir,"topol.tpr"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
            self.filesystem.deleteFileList(workdir,"#*")
            self.filesystem.deleteFileList(workdir,"ener.edr")
    
            
    def finalizeReverse(self,dirstring,newdirstring):
        self.finalizeShootingReverse(dirstring)
        for i in range(len(self.reversePaths)):
            self.finalizeCopyLastAcceptedReverse(dirstring, newdirstring, i)
    
    
    
    def finalizeCopyLastAcceptedReverse(self,dirstring,newdirstring,pathnumber):
        
        def _copyFiles(fromdir,destdir):
            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
        
        lapath = self.reversePaths[pathnumber].nfsladir
        newladir = os.path.join(lapath,newdirstring)
        newladirf = os.path.join(newladir,"forward")
        newladirb = os.path.join(newladir,"backward")
        
        dirlist = [lapath,newladir,newladirf,newladirb]
        self.filesystem.createDirList(dirlist, self.log.log)
        
        fdir = os.path.join(self.reversePaths[pathnumber].workdir,"paths",dirstring)
        
        _copyFiles(fdir, newladirf)
        
        

    def getFullTrajectoryReverse(self,pathnumber,dirstring):
        
        fdir = os.path.join(self.reversePaths[pathnumber].nfsladir,dirstring,"forward")
        
        traj = []
        count = 0
        os.chdir(fdir)
        cmd = self.wrapper.generateGDistCommand(fdir,self.reversePaths[1].options)
        self.wrapper.executeCommand(cmd, self.reversePaths[1].options.runoptions["dist1"]+"\n"+self.reversePaths[1].options.runoptions["dist2"]+"\n" )

        self.distparser.readDist(os.path.join(fdir,"dist.xvg"))
        self.filesystem.deleteFile(os.path.join(fdir,"dist.xvg"))
        self.reversePaths[pathnumber].fullTrajectoryflength=0
        for line in self.distparser.data[1:]:
            traj.append([count,float(line[self.stablestates.gdistDirection]),1])
            self.reversePaths[pathnumber].fullTrajectoryflength+=1
            count += 1
        
        self.reversePaths[pathnumber].fullTrajectory = traj[:]
        
    
    def _outputFullTrajectoryReverse(self,directory,pathnumber,dirstring):
        
        of = open(os.path.join(directory,"trajectoryZ."+dirstring+".acc"),"w")
        for position in self.reversePaths[pathnumber].lastAcceptedFullTrajectory:
            of.write("%d %.18f %d\n" % (position[0],position[1],position[2]))
        of.close()
    
    
    def outputAllFullTrajectoriesReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            ladir = os.path.join(self.reversePaths[i].nfsladir,dirstring)
            self._outputFullTrajectoryReverse(ladir, i,dirstring)
            
    def lastReverseToGro(self,dirstring):
        for i in range(len(self.reversePaths)):
            mpath = self.reversePaths[i]
            workdir = os.path.join(mpath.nfsladir,dirstring,"forward")
            cmd = self.wrapper.generateTrjconvCommand(mpath.options,workdir,workdir)
            self.wrapper.executeCommand(cmd, tinput="0")
    
    def copyStateABGroFiles(self,dirstring):

        for i in range(len(self.reversePaths)):
            slist = self.filesystem.getFileList(self.reversePaths[i].nfsladir,"path*")
            ns=len(slist)
            traj = self.reversePaths[i].fullTrajectory
            if self.reversePaths[i].state == 0:
                for point in traj:
                    if point[1] > self.stablestates.states[0]:
                        ns += 1
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.gro" % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,"path%d.gro" % (ns,))
                        self.filesystem.copyFiles(source, destination)
            else:
                for point in traj:
                    if point[1] < self.stablestates.states[1]:
                        ns += 1
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.gro" % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,"path%d.gro" % (ns,))
                        self.filesystem.copyFiles(source, destination)
            
            rn = self.reversePaths[i].fullTrajectory[-1][0]

            source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.gro" % rn)
            destination = os.path.join(self.reversePaths[i].nfsladir,dirstring,"conf.gro")
            self.filesystem.copyFiles(source, destination)
    
    def pickConfigurationReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            filename = os.path.join(self.reversePaths[i].nfsladir,dirstring,"conf.gro")
            self.filesystem.copyFiles(filename, self.reversePaths[i].workdir)

    def writeFinishedFilesReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            lapath = self.reversePaths[i].nfsladir
            newladir = os.path.join(lapath,dirstring)
            filename = os.path.join(newladir,"finished")
            of = open(filename,"w")
            of.write("finished")
            of.close()
    
            
