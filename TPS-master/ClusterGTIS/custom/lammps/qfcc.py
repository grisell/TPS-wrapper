'''
Created on Oct 26, 2010

@author: wolf
'''

#import pygromacstps.wrappers
import pygromacstps.parser
import pygromacstps.filesystem
import math

import os
import sys
import time
import shutil
import wrappers as wrappers

import wrappers as wrappers

# The path where the program water_around is located
#wapath = "/home/cvreede/GTIS/pypwater/watertest"



class TCustomOrderParameter(object):
    def __init__(self,minA,maxA,minB,maxB):
        # The stable states
        self.A = [minA,maxA]
        self.B = [minB,maxB]
        # dircolumn is the column from the output where the interesting coordinated is
	# JR: might not be needed or we might extend this
	#self.dircolumn = 7
        self.wrapper    = wrappers.wrapper()
	
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()

    def _makecalcOP(self,path,workdir):
	cmd = []
	cmd.append(os.path.join(path.options.paths["binpath"],"opcalc.exe"))
	cmd.append(os.path.join(workdir,"tmp_conf"))
	return cmd
        
    
        
    def isPathInState(self,path,log):
        # Does the path end in one of the two stable states?
        # get trajectory and look at the last coordinate
	# get the last slice of current trajectory
	filename = os.path.join(path.workdir,path.options.initoptions["pathfile"][1])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
		
		if len(data) > 4:
			natoms = int(data[3])
		else:
			return 0,False
		#check if file is completely written
		if len(data)%(natoms+9) != 0:
			return 0,False
		if len(data) > natoms:
			#write last slice to tmp file
			tmpfile = os.path.join(path.workdir,"tmp_conf")
			outfile = open(tmpfile,"w")
			for i in range(-(natoms+9),0):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			#calculate order parameter
			#using external program...
			cmd = self._makecalcOP(path,path.workdir)
			output,outerr = self.wrapper.executeCommand(cmd)
			q = float(output)
			os.remove(tmpfile)
		else:
			return 0,False


		if q > self.A[0] and q < self.A[1]:
                	return 1,True
            	elif q > self.B[0] and q < self.B[1]:
                	return 4,True
            	else:
                	return 0,False

	else:
	  	return 0,False

				
				
 
    def getQTrajectory(self,path,wdir,data):
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
	for j in range(nslices):
		start = j*nblock
		end = (j+1)*nblock
		tmpfile = os.path.join(wdir,"tmp_conf")
		outfile = open(tmpfile,"w")
		for i in range(start,end):
			outfile.write(data[i])
		outfile.flush()
		outfile.close()
		cmd = self._makecalcOP(path,wdir)
		output,outerr = self.wrapper.executeCommand(cmd)
		q = float(output)

		qtraj.append(q)

	os.remove(tmpfile)
	
	return qtraj




    def getFullTrajectory(self,fpath,bpath,fdir,bdir,log):
        traj = []
	qtrajback = []
	qtrajfor = []
	flength=0
	blength=0
	#first the backward part
        os.chdir(bdir)
	filename = os.path.join(bdir,bpath.options.initoptions["pathfile"][1])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'backward part of trajectory not fully written in workdir, removing unfinished slice...'
			newlines = nblock*nslices
			tmpfile = os.path.join(bdir,'tmp-traj.dat')
			outfile = open(tmpfile,"w")
			for i in range(newlines):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			infile.close()
			shutil.move(tmpfile,filename)

			#re-read file with proper slices
			data = []
			infile = open(filename,"r")
			for line in infile:
				data.append(line)

			nlines = len(data)
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)

	if len(data) > 0:
		qtrajback = self.getQTrajectory(bpath,bdir,data)
		
        	qtrajback.reverse()
        	blength=0
        	count = 0
#		print 'traback, data',len(qtrajback),len(data)
#		sys.stdout.flush()
        	for i in range(len(qtrajback)):
            		traj.append([count,qtrajback[i],0])
            		blength+=1
            		count += 1
        
	#now the forward part
        os.chdir(fdir)
	filename = os.path.join(fdir,fpath.options.initoptions["pathfile"][1])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'forward part of trajectory not fully written in workdir, removing unfinished slice...'
			newlines = nblock*nslices
			tmpfile = os.path.join(fdir,'tmp-traj.dat')
			outfile = open(tmpfile,"w")
			for i in range(newlines):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			infile.close()
			shutil.move(tmpfile,filename)

			#re-read file with proper slices
			data = []
			infile = open(filename,"r")
			for line in infile:
				data.append(line)

			nlines = len(data)
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)


	

	if len(data) > 0:
		qtrajfor = self.getQTrajectory(fpath,fdir,data)
        
        	flength=0
        	for i in range(len(qtrajfor)):
            		traj.append([count,qtrajfor[i],1])
            		flength+=1
            		count += 1

        if len(traj) > 0:
        	for i in range(len(traj)):
            		traj[i][0] = traj[i][0] - blength + fpath.shootingTime
	else:
		traj.append([0,0.0,-1])
        
        return traj,flength,blength



    def getFullTrajectoryReverse(self,fpath,fdir,log):
        traj = []
	qtrajfor = []
	flength=0
        
	#now the forward part
        os.chdir(fdir)
	filename = os.path.join(fdir,fpath.options.initoptions["pathfile"][1])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
	
	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'reverse trajectory not fully written in workdir, removing unfinished slice...'
			newlines = nblock*nslices
			tmpfile = os.path.join(fdir,'tmp-traj.dat')
			outfile = open(tmpfile,"w")
			for i in range(newlines):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			infile.close()
			shutil.move(tmpfile,filename)

			#re-read file with proper slices
			data = []
			infile = open(filename,"r")
			for line in infile:
				data.append(line)

			nlines = len(data)
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)	

	if len(data) > 0:
		qtrajfor = self.getQTrajectory(fpath,fdir,data)
        
        	flength=0
		count=0
        	for i in range(len(qtrajfor)):
            		traj.append([count,qtrajfor[i],1])
            		flength+=1
            		count += 1

        if len(traj) < 1:
		traj.append([0,0.0,-1])
        
        return traj,flength




    def QofLast(self,path):
        # Q value of the last slice
        # get trajectory and look at the last coordinate
	# get the last slice of current trajectory
	filename = os.path.join(path.workdir,path.options.initoptions["pathfile"][1])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

		if len(data) > 4:
			natoms = int(data[3])
		else:
			return 0,0
		
	        nlines = len(data)
		nblock = natoms+9
		nslices = nlines/nblock
		#check if file is completely written
		if nlines%nblock != 0:
			return 0,0

		if len(data) > natoms:
			#write last slice to tmp file
			tmpfile = os.path.join(path.workdir,"tmp_conf")
			outfile = open(tmpfile,"w")
			for i in range(-(natoms+9),0):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			#calculate order parameter
			#using external program...
			cmd = self._makecalcOP(path,path.workdir)
			output,outerr = self.wrapper.executeCommand(cmd)
			q = float(output)
			os.remove(tmpfile)

			return nslices,q
		else:
			return 0,0

	else:
	  	return 0,0

				

    
    
 
