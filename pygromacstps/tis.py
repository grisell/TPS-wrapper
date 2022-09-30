
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
import helpers_lammps as helpers
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
#        self.cores = multiprocessing.cpu_count()	#JR originial, set to 2 for testing
        self.cores = 2
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

        
        """
        This array holds all the information about the paths. Each trajectory consists of a 
        forward-part and a backward-part. The trajectory can either be a forward trajectory
        or a backward trajectory. Forward trajectroies start in A.
        """
        
        self.paths = []
        self.npaths = 0
        #dvmax=self.paths[0].options.runoptions["dvmax"]
	#print dvmax
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
                for myfile in self.paths[i].options.standardfiles:
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",myfile), workdir)

	    if self.paths[i].options.runoptions["mddriver"] in ('lammps','lammps_parallel'):
	    	self.wrapper.modifyMDseed(workdir,self.paths[i].options)
	
	
            
            
    
    def shootingInitialGroFiles(self):
        """
        Perform a shooting move on a gro file.
        """
        for i in self.kernels.kernelPaths:
            conffile = os.path.join(self.paths[i].workdir,"conf.gro")
            bakfile = os.path.join(self.paths[i].workdir,"bakconf.gro")
            bconffile = os.path.join(self.paths[i].workdir,"bconf.gro")
            self.helper.shootVelocities(self.paths[i].options,conffile, bconffile)  
            self.filesystem.moveFile(conffile, bakfile)
            self.filesystem.moveFile(bconffile, conffile)
            dest = self.paths[i+1].workdir
            self.filesystem.copyFiles(conffile, dest)
            self.reverseBackwardGroFile(i+1)



    def shootingInitialLammpsFiles(self):
        """
        Perform a shooting move on a gro file.
        """
        for i in self.kernels.kernelPaths:
            conffile = os.path.join(self.paths[i].workdir,"conf.dump")
            bakfile = os.path.join(self.paths[i].workdir,"bakconf.dump")
            bconffile = os.path.join(self.paths[i].workdir,"bconf.dump")
            self.helper.shootVelocities(self.paths[i].options,conffile, bconffile)  
            self.filesystem.moveFile(conffile, bakfile)
            self.filesystem.moveFile(bconffile, conffile)
            dest = self.paths[i+1].workdir
            self.filesystem.copyFiles(conffile, dest)
            self.reverseBackwardLammpsFile(i+1)

            
                 
                         
    def finalizeInitial(self):
        
        def _copyFiles(fromdir,destdir):
	#   self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
	#    self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
        
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
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            dirlist = [os.path.join(initdir,'all'),os.path.join(initdir,'la'),self.paths[i].nfsbaseworkdir,self.paths[i].nfsladir]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                for myfile in self.paths[i].options.standardfiles:
                    self.filesystem.copyFiles(os.path.join(initdir,"standardfiles",myfile), workdir)
           
	    if self.paths[i].options.runoptions["mddriver"] in ('lammps','lammps_parallel'):
		self.wrapper.modifyMDseed(workdir,self.paths[i].options)
  
    

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
	    # JR: this means if forwardpart = True -> assign "forward", else assign "backward"
            fbdir = mpath.forwardpart and "forward" or "backward"
            nfsworkdir = os.path.join(mpath.nfsladir,dirnumber,fbdir)
            latemp = mpath.latemp
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"traj.trr") , latemp)
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"index.ndx") , latemp)
            self.filesystem.copyFiles(os.path.join(nfsworkdir,"topol.tpr") , latemp)
            cmd = self.wrapper.generateTrjconvCommand(mpath.options,latemp,latemp)
            self.wrapper.executeCommand(cmd, tinput="0")



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
            self.filesystem.copyFiles(os.path.join(nfsworkdir,mpath.options.initoptions["pathfile"][1]) , latemp)
            self.wrapper.TrajToConf(mpath.options,latemp)

            
    
    def pickConfigurationsTIS(self,dirnumber):
        for i in self.kernels.kernelPaths:
            self.pickConfigurationLastAccepted(dirnumber,i)  
    
    def pickConfigurationLastAccepted(self,dirnumber,pathnumber):
        """
        A random configuration is taken from a trajectory and the shooting move is performed.
        and copy the files to the workdir of the path.
        """
        
        fdir = os.path.join(self.paths[pathnumber].latemp)
        flist = self.filesystem.getFileList(fdir, "path*.dump")
        fsize = len(flist)
        bdir = os.path.join(self.paths[pathnumber+1].latemp)
        blist = self.filesystem.getFileList(bdir, "path*.dump")
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
                rn = random.randint(1,sum-2)
                if self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] > self.orderparameters.op[distopi].A[1] and self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] < self.orderparameters.op[distopi].B[0]:
                    break
                c += 1
                if c > 100:
                    print "ERROR picking shootingpoint"
		    sys.stdout.flush()
                    break
            
        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 1:
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
###################################### DS-pick a slice above the interface #####################################################################################DS
#                if accept:                      
#                    shootpart.append(kr)        
#
#            while True:                         
#                random.shuffle(shootpart)      
#                rn = shootpart[0]
#                if self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] > self.orderparameters.op[distopi].A[1] and self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] < self.orderparameters.op[distopi].B[0]:
#                    break
#                c += 1
#                if c > 100:
#                    print "ERROR picking shootingpoint"
#                    sys.stdout.flush()
#                    break
#            self.log.log.debug("Picked " +str(rn) + " " + str(self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1])+" " + str(self.paths[pathnumber].interface) + "  "+ str(mindist)+  " " + str(pathnumber) )
######################################## DS- pick a slice between fi-1 and fi+1 #################################################################################

        elif self.paths[pathnumber].options.runoptions["shootfrominterface"] == 2:
            shootpart = []     #DS - instead of shooting from the closest interslice to the interface now I choose ramdomly an interslice between fi-1 and fi+1
            rn = sum/2
            for kr in range(1,sum-1):
                act = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][kr][1]

                interfacetot = self.interfaces.interfaces[0] +  self.interfaces.interfaces[1]
                index = interfacetot.index(self.paths[pathnumber].interface)

                if ((index == 0) and (bool(self.paths[pathnumber].forward == 1))):                                                                 #for first interface
                    accept = (bool(self.paths[pathnumber].forward) == bool((act > interfacetot[index]) and (act < interfacetot[index+1])))

                elif ((index == len(interfacetot)-1) and (bool(self.paths[pathnumber].forward == 0))):                                                 #for last interface
                    accept = (bool(self.paths[pathnumber].forward) != bool((act < interfacetot[index]) and (act > interfacetot[index-1])))


                elif ((index == len(self.interfaces.interfaces[0])-1) and (bool(self.paths[pathnumber].forward == 1))):                                #for last f interface
                    accept = (bool(self.paths[pathnumber].forward) == bool((act > interfacetot[index-1]) and (act < interfacetot[index+2])))

                elif ((index == len(self.interfaces.interfaces[0])-1) and (bool(self.paths[pathnumber].forward == 0))):                                 #for first b interface
                    accept = (bool(self.paths[pathnumber].forward) != bool((act > interfacetot[index-1]) and (act < interfacetot[index+2])))

                elif (index == len(self.interfaces.interfaces[0])+1):                                                                                   #for second b interface
                    accept = (bool(self.paths[pathnumber].forward) != bool((act > interfacetot[index-2]) and (act < interfacetot[index+1])))

                elif ((index < len(self.interfaces.interfaces[0])-1) and (index > 0) and (bool(self.paths[pathnumber].forward == 1))):                  #for f interfaces
                    accept = (bool(self.paths[pathnumber].forward) == bool((act > interfacetot[index-1]) and (act < interfacetot[index+1])))

                elif ((index > len(self.interfaces.interfaces[0])+1) and (index < len(interfacetot)-1) and (bool(self.paths[pathnumber].forward == 0))): #for b interfaces
                    accept = (bool(self.paths[pathnumber].forward) != bool((act > interfacetot[index-1]) and (act < interfacetot[index+1])))

                
                if accept:                      
                          shootpart.append(kr) 
                    
            while True:                         
                 random.shuffle(shootpart)      
                 rn = shootpart[0]
                 if self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] > self.orderparameters.op[distopi].A[1] and self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][1] < self.orderparameters.op[distopi].B[0]:
                     break
                 c += 1
                 if c > 100:
                     print "ERROR picking shootingpoint"
                     sys.stdout.flush()
                     break
###############################################################################################################################################################DS
        
        backwardpart = (rn < bsize)
        if backwardpart:
            myfile = os.path.join(bdir,"path%d.dump" % ((bsize-1)-rn))
        else:
            myfile = os.path.join(fdir,"path%d.dump" % (rn-bsize))
        
        self.log.log.debug("pick " +dirnumber + " " + str(pathnumber) +" " + str(rn) +  " " + str( sum) +  " " + str(fsize) + " " +  str(bsize) + " " +  str(len(self.paths[pathnumber].lastAcceptedFullTrajectory)) + " " + myfile)   # JR: was commented out
       
        self.paths[pathnumber].shootingTime = self.paths[pathnumber].lastAcceptedFullTrajectory[distopi][rn][0]
        
        destfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"conforig.dump")
        shootfile = os.path.join(self.paths[pathnumber].nfsladir,dirnumber,"confshoot.dump")
        
        if os.path.exists(myfile):    
            self.filesystem.copyFiles(myfile,destfile)
        else:
#            self.log.log.ERROR("File " + myfile + " not found, this is a fatal error")  #JR original
            self.log.log.error("File " + myfile + " not found, this is a fatal error")
	    print 'Could not identify shooting point'
	    print 'rn = ',rn
	    print 'fsize = ',fsize,' bsize = ',bsize
	    print myfile
	    print 'Exiting program...'
	    sys.exit(1)
        self.helper.shootVelocities(self.paths[pathnumber].options,destfile,shootfile)
        
        for path in [self.paths[pathnumber],self.paths[pathnumber+1]]:
            self.filesystem.copyFiles(shootfile, os.path.join(path.workdir,"conf.dump"))
            path.finishedstate = -1
        
        if backwardpart:
            reversepath = pathnumber 
        else:
            reversepath = pathnumber + 1
        self.reverseBackwardLammpsFile(reversepath)
        
        self.log.log.debug("shooting from " +str(rn)+ " file " + myfile + " backwardpart  " + str(backwardpart))
        
        self.filesystem.deleteFileList(fdir, "path*.dump")
        self.filesystem.deleteFileList(bdir, "path*.dump")
        
        
        
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


    def reverseBackwardLammpsFile(self,pathnumber):
        """
        Reverse all velocities from a lammps dump file.
        """
        conffile = os.path.join(self.paths[pathnumber].workdir,"conf.dump")
        bakfile = os.path.join(self.paths[pathnumber].workdir,"bakconf.dump")
        bconffile = os.path.join(self.paths[pathnumber].workdir,"bconf.dump")
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
#	for i in self.kernels.kernelPathsAll:
#            self.gromppPath(self.paths[i])
        
        # initialize the reverse Paths
#        for rpath in self.reversePaths:
#            self.gromppPath(rpath)
        
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
            #print subprocess.Popen("ls", shell=True, stdout=subprocess.PIPE).stdout.read()
#	    tmp = os.getcwd() #JR
            #print cmd #JR
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
                        pathstate,thisfinished = self.orderparameters.op[opi].isPathInState(self.paths[i],self.log)
                        self.paths[i].finishedState = pathstate
                        pathfinished = thisfinished and pathfinished
                    if pathfinished:
                        self.log.log.debug("Finished " + str(i) + " " + str(self.paths[i].finishedState))
			print 'terminating process',i			#JR
			sys.stdout.flush()
                        processlist[i].terminate()
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
            
#            self.filesystem.moveFile(os.path.join(workdir,"traj.xtc"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"traj.trr"), dest)
#            self.filesystem.moveFile(os.path.join(workdir,"topol.tpr"), dest)
            self.filesystem.moveFile(os.path.join(workdir,"traj.dat"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"conf.dump"), dest)
            self.filesystem.copyFiles(os.path.join(workdir,"md.in"), dest)
            self.filesystem.deleteFileList(workdir,"#*")
#            self.filesystem.deleteFileList(workdir,"ener.edr")

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
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
            
#	    self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
	    self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#	    self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
            
        
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
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
        
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
#            self.filesystem.copyFiles(os.path.join(fromdir,"index.ndx"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"topol.tpr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.gro"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"traj.xtc"), destdir)
#            self.filesystem.moveFile(os.path.join(fromdir,"traj.trr"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"conf.dump"), destdir)
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)


        def _copyQTrajectory(fromdir,destdir,dirstring,newdirstring):
            fileold = os.path.join(fromdir,"trajectory.00."+dirstring+".dat") 
            fileoldtrial = os.path.join(fromdir,"trajectory.00."+dirstring+".trial") 
            filenew = os.path.join(destdir,"trajectory.00."+newdirstring+".dat")
            filenewtrial = os.path.join(destdir,"trajectory.00."+newdirstring+".trial")
            
	    self.filesystem.copyFiles(fileold, filenew)
	    self.filesystem.copyFiles(fileoldtrial, filenewtrial)


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
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
        
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
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)



        def _copyQTrajectory(fromdir,destdir,dirstring,newdirstring):
            fileold = os.path.join(fromdir,"trajectory.00."+dirstring+".dat") 
            fileoldtrial = os.path.join(fromdir,"trajectory.00."+dirstring+".trial") 
            filenew = os.path.join(destdir,"trajectory.00."+newdirstring+".dat")
            filenewtrial = os.path.join(destdir,"trajectory.00."+newdirstring+".trial")
            
	    self.filesystem.copyFiles(fileold, filenew)
	    self.filesystem.copyFiles(fileoldtrial, filenewtrial)

	    
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
            self.filesystem.createDirList(dirlist, self.log.log)            
            
            dirlist = [os.path.join(initdir,'la'),rpath.nfsbaseworkdir,rpath.nfsladir,rpath.nfsconfigstore]
            self.filesystem.createDirList(dirlist, self.log.log)             
            
            if copyfiles:  
                for myfile in rpath.options.standardfiles:
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
            time.sleep(float(self.reversePaths[0].options.runoptions["updatetime"]))
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
            self.filesystem.moveFile(os.path.join(workdir,"traj.dat"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"index.ndx"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.gro"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"conf.dump"), dest)
#            self.filesystem.copyFiles(os.path.join(workdir,"md.in"), dest)
            self.filesystem.deleteFileList(workdir,"#*")
#            self.filesystem.deleteFileList(workdir,"ener.edr")
    
            
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
            self.filesystem.copyFiles(os.path.join(fromdir,"traj.dat"), destdir)
#            self.filesystem.copyFiles(os.path.join(fromdir,"md.in"), destdir)
        
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
            
    def lastReverseToGro(self,dirstring):
        for i in range(len(self.reversePaths)):
            mpath = self.reversePaths[i]
            workdir = os.path.join(mpath.nfsladir,dirstring,"forward")
            cmd = self.wrapper.generateTrjconvCommand(mpath.options,workdir,workdir)
            self.wrapper.executeCommand(cmd, tinput="0")
 

    def lastReverseToConf(self,dirstring):
        """
        The last reverse trajectory is transformed into pathXX.dump files 
        """
        for i in range(len(self.reversePaths)):
            mpath = self.reversePaths[i]
            workdir = os.path.join(mpath.nfsladir,dirstring,"forward")
#            latemp = mpath.latemp								#JR: for 'normal' path this is done in tmp directory...
#            self.filesystem.copyFiles(os.path.join(workdir,mpath.options.initoptions["pathfile"][1]) , latemp)
#            self.wrapper.TrajToConf(mpath.options,latemp)
            self.wrapper.TrajToConf(mpath.options,workdir)					#JR: here now in 'regular' directory


    
    def copyStateABConfFiles(self,dirstring):

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
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.dump" % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,"path%d.dump" % (ns,))
                        self.filesystem.copyFiles(source, destination)
            else:
                for point in traj:
#                    if point[1] < self.stablestates.states[1]:
                    if point[1] < self.reversePaths[i].interface:
                        ns += 1
                        source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.dump" % point[0])
                        destination = os.path.join(self.reversePaths[i].nfsconfigstore,"path%d.dump" % (ns,))
                        self.filesystem.copyFiles(source, destination)
            
#            rn = self.reversePaths[i].fullTrajectory[-1][0]
            rn = traj[-1][0]

            source = os.path.join(self.reversePaths[i].nfsladir,dirstring,"forward","path%d.dump" % rn)
            destination = os.path.join(self.reversePaths[i].nfsladir,dirstring,"conf.dump")
            self.filesystem.copyFiles(source, destination)
    
    def pickConfigurationReverse(self,dirstring):
        for i in range(len(self.reversePaths)):
            filename = os.path.join(self.reversePaths[i].nfsladir,dirstring,"conf.dump")
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
           

            accepted = self.reversePaths[i].checkAcceptedTIS(self.orderparameters.op)
            self.reversePaths[i].tisaccepted = accepted
            
            if accepted or self.mode=="initial":
                self.reversePaths[i].lastAcceptedFullTrajectory = []
                for opi in range(len(self.orderparameters.op)):
                    self.reversePaths[i].lastAcceptedFullTrajectory.append(self.reversePaths[i].fullTrajectory[opi][:])
                    self.reversePaths[i].lastAcceptedFullTrajectoryflength = self.reversePaths[i].fullTrajectoryflength
            report.append(str(accepted))
        self.log.log.info(" ".join(report))
    
            
