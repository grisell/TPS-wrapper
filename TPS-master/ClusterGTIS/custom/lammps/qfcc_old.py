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

# The path where the program water_around is located
wapath = "/home/cvreede/GTIS/pypwater/watertest"


class TCustomOrderParameter(object):
    def __init__(self,minA,maxA,minB,maxB):
        # The stable states
        self.A = [minA,maxA]
        self.B = [minB,maxB]
        # dircolumn is the column from the output where the interesting coordinated is
	# JR: might not be needed or we might extend this
        self.dircolumn = 7
        
        self.wrapper    = wrappers.wrapper()
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()
        
    
        
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
		raw=[]
		raw.append(data[8].split())
		ncolumns=len(raw[0])
		for i in range(ncolumns):
			if raw[0][i] == "x":
				xcol = i-2
			if raw[0][i] == "y":
				ycol = i-2
			if raw[0][i] == "z":
				zcol = i-2

		if len(data) > natoms:
			#get box size
			raw=[]
			for i in range(-natoms-4,-natoms-1):
				raw.append(data[i].split())

			boxlength = []
			for i in range(3):
				tmp = float(raw[i][1]) - float(raw[i][0])
				boxlength.append(tmp)
			#get coordinates
			raw=[]
			for i in range(-natoms,0):
				raw.append(data[i].split())
			#copy coordinates of last slice
			coord = [i[:] for i in [[0.0]*3]*natoms]
			#JR debug: check array size
			if len(raw) < natoms:
				print len(raw)
				print natoms
				print raw[0]
				print 'something wrong here, exiting'
				sys.exit(1)
			if len(raw[0]) < zcol:
				print len(raw[0])
				print xcol,ycol,zcol
				print raw[0]
				print 'not enough columns?? Exiting'
				sys.exit(1)
			#JR end debug
			for i in range(natoms):
				coord[i][0] = float(raw[i][xcol])
				coord[i][1] = float(raw[i][ycol])
				coord[i][2] = float(raw[i][zcol])
			q,qatom = self.calculateQ(coord,boxlength)
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

				
				

    def calculateQ(self,coord,boxlength):
	# NN cutoff => must be adjusted to the system!!
	r1_cv = 2.9
	r0_cv = 3.1
	dr_cv = r0_cv - r1_cv

	natoms = len(coord)
	part_c = [i[:] for i in [[0.0]*2]*natoms]
	    
	for i in range(natoms-1):
		for j in range(i+1,natoms):
			dr = self.imageDist(coord[i],coord[j],boxlength) 
#			print dr
			r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]
			r = math.sqrt(r2)
			r4 = r2*r2

			x4 = dr[0]*dr[0]*dr[0]*dr[0]
			y4 = dr[1]*dr[1]*dr[1]*dr[1]
			z4 = dr[2]*dr[2]*dr[2]*dr[2]

			c_alpha = (x4*y4*(1.0 - z4/r4) + y4*z4*(1.0 - x4/r4) + z4*x4*(1.0 - y4/r4))/(r4*r4)
			
			if r <= r1_cv:
				c_r = 1.0
			elif r >= r0_cv:
				c_r = 0.0
			else:
				c_r = ((r-r1_cv)/dr_cv - 1.0) * ((r-r1_cv)/dr_cv - 1.0) * (1.0 + 2.0*(r-r1_cv)/dr_cv)
			
			part_c[i][0] += c_r
			part_c[j][0] += c_r

			part_c[i][1] += c_alpha*c_r
			part_c[j][1] += c_alpha*c_r

	phi = []
	q = 0.0
	for i in range(natoms):
		if part_c[i][0] == 0:
			phiatom = 0.0
		else:
			phiatom = part_c[i][1]/part_c[i][0]
			phiatom = 2288.0/79.0 * phiatom - 64.0/79.0;

		phi.append(phiatom)
		q += phiatom
	
	q /= natoms

	return q,phi



    def imageDist(self,coord_i,coord_j,boxlength):
	dr = []
	for i in range(3):
		tmp = coord_j[i] - coord_i[i]
		if tmp > (boxlength[i]/2.0):
			tmp -= boxlength[i]
		elif tmp < -(boxlength[i]/2.0):
			tmp += boxlength[i]
		
		dr.append(tmp)
	
	return dr



 
    def getQTrajectory(self,data):
	nlines = len(data)
	natoms = int(data[3])
	nblock = natoms+9
	nslices = nlines/nblock

	if nlines%nblock != 0:
		print 'Error in getQTrajectory'
		print 'nlines%nblock != 0'
		print 'Exiting program...'
		sys.exit(1)

	raw=[]
	raw.append(data[8].split())
	ncolumns=len(raw[0])
	for i in range(ncolumns):
		if raw[0][i] == "x":
			xcol = i-2
		if raw[0][i] == "y":
			ycol = i-2
		if raw[0][i] == "z":
			zcol = i-2
	
	qtraj=[]
	for j in range(nslices):
		start = j*nblock
		#get box size
		raw=[]
		for i in range(start+5,start+8):
			raw.append(data[i].split())

		boxlength = []
		for i in range(3):
			tmp = float(raw[i][1]) - float(raw[i][0])
			boxlength.append(tmp)
		#get coordinates
		raw=[]
		for i in range(start+9,start+9+natoms):
			raw.append(data[i].split())
		#copy coordinates of current slice j
		coord = [i[:] for i in [[0.0]*3]*natoms]
		for i in range(natoms):
			coord[i][0] = float(raw[i][xcol])
			coord[i][1] = float(raw[i][ycol])
			coord[i][2] = float(raw[i][zcol])
		q,qatom = self.calculateQ(coord,boxlength)

		qtraj.append(q)
	
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
		qtrajback = self.getQTrajectory(data)
		
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
		qtrajfor = self.getQTrajectory(data)
        
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
		qtrajfor = self.getQTrajectory(data)
        
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
		
		natoms = int(data[3])
	        nlines = len(data)
		nblock = natoms+9
		nslices = nlines/nblock
		#check if file is completely written
		if nlines%nblock != 0:
			return 0,0

		raw=[]
		raw.append(data[8].split())
		ncolumns=len(raw[0])
		for i in range(ncolumns):
			if raw[0][i] == "x":
				xcol = i-2
			if raw[0][i] == "y":
				ycol = i-2
			if raw[0][i] == "z":
				zcol = i-2

		if len(data) > natoms:
			#get box size
			raw=[]
			for i in range(-natoms-4,-natoms-1):
				raw.append(data[i].split())

			boxlength = []
			for i in range(3):
				tmp = float(raw[i][1]) - float(raw[i][0])
				boxlength.append(tmp)
			#get coordinates
			raw=[]
			for i in range(-natoms,0):
				raw.append(data[i].split())
			#copy coordinates of last slice
			coord = [i[:] for i in [[0.0]*3]*natoms]
			#JR debug: check array size
			if len(raw) < natoms:
				print len(raw)
				print natoms
				print raw[0]
				print 'something wrong here, exiting'
				sys.exit(1)
			if len(raw[0]) < zcol:
				print len(raw[0])
				print xcol,ycol,zcol
				print raw[0]
				print 'not enough columns?? Exiting'
				sys.exit(1)
			#JR end debug
			for i in range(natoms):
#				print i,xcol,ycol,zcol   #JR
				coord[i][0] = float(raw[i][xcol])
				coord[i][1] = float(raw[i][ycol])
				coord[i][2] = float(raw[i][zcol])
			q,qatom = self.calculateQ(coord,boxlength)

			return nslices,q
		else:
			return 0,0

	else:
	  	return 0,0

				

    
    
 
