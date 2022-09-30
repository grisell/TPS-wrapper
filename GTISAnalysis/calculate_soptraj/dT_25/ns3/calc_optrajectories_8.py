#!/usr/bin/env python

'''
Created Feb 19, 2015

@author: grisell
'''
#this script will create an order parameter trajectory for each path of the TIS ensemble
#you must provide the name and path of the executable of the order parameter in the binary and the name of your output trajectory. 

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
import subprocess as sub

binary = "/home/users/diazlgjj/bin/ns3sp25"
outfilename="ns3sptrajectory"

class getdata(object):
	
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
        
        if len(self.interfaces.reversepaths[0]):
            for rp in self.interfaces.reversepaths[0]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=0))
        if len(self.interfaces.reversepaths[1]):
            for rp in self.interfaces.reversepaths[1]:
                self.reversePaths.append(pathdatareverse.pathdatareverse(0, basedir, mode, forward=True, forwardpart=True, interface=rp, state=1))


        self.kernels.generateKernelLists(self.npaths)

    def getQTrajectory(self,wdir,data):
	    nlines = len(data)
	    natoms = int(data[3])
	    nblock = natoms+9
	    nslices = nlines/nblock

	    if nlines%nblock != 0:
		    print 'Error in getQTrajectory'
		    print 'nlines%nblock != 0'
		    print 'Exiting program...'
		    sys.exit(1)

            qtraj=[]
            q1traj=[]
            q2traj=[]
            q3traj=[]
            q4traj=[]
            q5traj=[]

            #print '#No. of slices = ',nslices
	    for j in range(nslices):
		    start = j*nblock
		    end = (j+1)*nblock
		    tmpfile = os.path.join('/tmp','mytmp_conf')
		    outfile = open(tmpfile,'w')
		    for i in range(start,end):
			    outfile.write(data[i])
		    outfile.flush()
		    outfile.close()
		    cmd = []
		    cmd.append(binary)
		    cmd.append(tmpfile)
		    proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
		    out,err = proc.communicate(input="")
		    proc.wait() 
		    #q = float(out)
                    q=out.split()
                    q0=float(q[0])
                    q1=float(q[1])
                    q2=float(q[2])
                    q3=float(q[3])
                    q4=float(q[4])
                    q5=float(q[5])
                    
                    #q=out
		    #print q

		    qtraj.append(q0)
                    q1traj.append(q1)
                    q2traj.append(q2)
                    q3traj.append(q3)
                    q4traj.append(q4)
                    q5traj.append(q5)

	    os.remove(tmpfile)

	
	    return qtraj, q1traj, q2traj, q3traj, q4traj, q5traj

    def getFullTrajectory(self,basedir,ffile,bfile):
        traj = []
        qtrajback = []
        qtrajfor = []
        
        data = []
        if os.path.exists(bfile):
		infile = open(bfile,"r")
		for line in infile:
			data.append(line)
        if len(data) > 4:
            nlines = len(data)
            natoms = int(data[3])
            nblock = natoms+9
            nslices = nlines/nblock
            if nlines%nblock != 0:
                print 'trajectory not fully written in workdir, removing unfinished slice...'
                newlines = nblock*nslices
                data_tmp = []
                for i in range(newlines):
                    data_tmp.append(data[i])

                    data = []
                    data = data_tmp

                    nlines = len(data)
                    natoms = int(data[3])
                    nblock = natoms+9
                    nslices = nlines/nblock

                    if nlines%nblock != 0:
                        print 'Cannot trajectory file with proper number of lines'
                        print 'Exiting program...'
                        sys.exit(1)
    
        if len(data) > 0:
            qtrajback = self.getQTrajectory(basedir,data)
            qtrajback[0].reverse()
            qtrajback[1].reverse()
            qtrajback[2].reverse()
            qtrajback[3].reverse()
            qtrajback[4].reverse()
            qtrajback[5].reverse()
            blength=0
            count = 0
            for i in range(len(qtrajback[0])):
                traj.append([count,qtrajback[0][i],qtrajback[1][i],qtrajback[2][i],qtrajback[3][i],qtrajback[4][i],qtrajback[5][i],0])
                blength+=1
                count += 1


        data = []
        if os.path.exists(ffile):
            infile = open(ffile,"r")
            for line in infile:
                data.append(line)
                #print data
        if len(data) > 4:
            nlines = len(data)
            natoms = int(data[3])
            nblock = natoms+9
            nslices = nlines/nblock
            if nlines%nblock != 0:
                print 'trajectory not fully written in workdir, removing unfinished slice...'
                newlines = nblock*nslices
                data_tmp = []
                for i in range(newlines):
                    data_tmp.append(data[i])
                             
                             
                    data = data_tmp
                    
                    nlines = len(data)
                    natoms = int(data[3])
                    nblock = natoms+9
                    nslices = nlines/nblock
                        
                    if nlines%nblock != 0:
                        print 'Cannot trajectory file with proper number of lines'
                        print 'Exiting program...'
                        sys.exit(1)
        if len(data) > 0:
            qtrajfor = self.getQTrajectory(basedir,data)
            flength=0
            #print len( qtrajfor[0])
            for i in range(len(qtrajfor[0])):
                traj.append([count,qtrajfor[0][i], qtrajfor[1][i],qtrajfor[2][i], qtrajfor[3][i],qtrajfor[4][i],qtrajfor[5][i],1])
                flength+=1
                count += 1
        return (traj,flength,blength)

if __name__ == '__main__':
     
     basedir = os.path.join(os.getcwd(),"..")
     cdata = getdata(basedir,"tis")
     rstart=0
     rstop=0
     traj = []
     qtrajback = []
     qtrajfor = []
     
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
                
     for dirnumber in range(start,stop):
	     dirstring = "%07d" % dirnumber
	     for i in cdata.kernels.kernelPaths:
                 if i/2 == 8:
                     binfile= os.path.join(cdata.paths[i].nfsladir,dirstring,"backward","traj.dat")
                     finfile= os.path.join(cdata.paths[i].nfsladir,dirstring,"forward","traj.dat")
                     outfile= os.path.join(cdata.paths[i].nfsladir,dirstring,outfilename+"."+dirstring+".dat")
                     traj=cdata.getFullTrajectory(basedir,finfile,binfile)
                     of = open(outfile,"w")
                                        
                     for point in traj[0]:
                         #print point[0],point[1],point[2],point[3],point[4],point[5],point[6], point[7]
                         of.write("%d %.18f  %.18f %.18f %.18f %.18f %.18f %d\n" % (point[0],point[1],point[2],point[3],point[4],point[5],point[6],point[7]))
                     of.close() 

     for dirnumber in range(rstart,rstop):
         dirstring = "%07d" % dirnumber
         rfinfile= os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"forward","traj.dat")
         rbinfile= os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"backward","traj.dat")
         routfile=os.path.join(cdata.reversePaths[0].nfsladir,dirstring,"r"+outfilename+"."+dirstring+".acc")
        
         rtraj=cdata.getFullTrajectory(basedir,rfinfile,rbinfile)
             
         rof = open(routfile,"w")
         for point in rtraj[0]:
             rof.write("%d %.18f %d\n" % (point[0],point[1],point[2]))
         of.close() 
