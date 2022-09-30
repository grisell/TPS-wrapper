'''
Created on Mar 2, 2010

@author: wolf
'''
import os
import random
import sys
import math

class helpers(object):
    def __init__(self):
        self.wolf = "helpers"
        #self.options = options.lammpstpsoptions(basedir,mode)        
        #self.options.readOptions(os.path.join(basedir,"options","runoptions.txt"), self.options.runoptions)

    def reverseVelocities(self,filename,outfilename):
	data = []
        if os.path.exists(filename):                
            of = open(outfilename,"w")
            infile = open(filename,"r")
            for line in infile:
		data.append(line)
	    nlines = len(data)
	    for i in range(9):
		of.write(data[i])
	    raw = []
	    raw.append(data[8].split())
	    ncolumns=len(raw[0])
	    for i in range(ncolumns):
	 	if raw[0][i] == "vx":
    			vx=i-2
		if raw[0][i] == "vy":
	 		vy=i-2
       		if raw[0][i] == "vz":
			vz=i-2	


	    raw = []
	    for i in range(9,nlines):
		k=i-9
		raw.append(data[i].split())
		ncolumns=len(raw[k])
		# now modify velocities
		tmp = float(raw[k][vx])
		raw[k][vx] = -tmp 
		tmp = float(raw[k][vy])
		raw[k][vy] = -tmp 
		tmp = float(raw[k][vz])
		raw[k][vz] = -tmp

		for j in range(ncolumns):
			if j==vx or j==vy or j==vz:
				of.write('%.8f' % raw[k][j])
			else:
				of.write(raw[k][j])
			of.write(" ")
		of.write('\n')

            of.close()  




    def shootVelocities(self,options,filename,outfilename):	
        data = []

        dvmax = options.runoptions["dvmax"]	
       
        #print 'dvmax',dvmax
	#sys.stdout.flush() 
	if os.path.exists(filename):                
            of = open(outfilename,"w")
            infile = open(filename,"r")
            for line in infile:
		data.append(line)
	    nlines = len(data)
            natoms = int(data[3])
	    for i in range(9):
		of.write(data[i])
	    raw = []
	    raw.append(data[8].split())
	    ncolumns=len(raw[0])
	    for i in range(ncolumns):
	 	if raw[0][i] == "vx":
    			vx=i-2
		if raw[0][i] == "vy":
	 		vy=i-2
       		if raw[0][i] == "vz":
			vz=i-2	


	    raw = []
	    for i in range(9,nlines):
		raw.append(data[i].split())

	    
	    #old kinetic energy
	    Kold = 0.0
	    for i in range(natoms):
		    Kold += float(raw[i][vx])*float(raw[i][vx])
		    Kold += float(raw[i][vy])*float(raw[i][vy])
		    Kold += float(raw[i][vz])*float(raw[i][vz])
	    # now modify velocities
            for i in range(natoms):
		tmp = float(raw[i][vx])
		raw[i][vx] = tmp + random.gauss(0.0,1.0)*dvmax
		tmp = float(raw[i][vy])
		raw[i][vy] = tmp + random.gauss(0.0,1.0)*dvmax
		tmp = float(raw[i][vz])
		raw[i][vz] = tmp + random.gauss(0.0,1.0)*dvmax

	    #new kinetic energy
	    Knew = 0.0
	    for i in range(natoms):
	  	    Knew += raw[i][vx]*raw[i][vx]	
	  	    Knew += raw[i][vy]*raw[i][vy]	
	  	    Knew += raw[i][vz]*raw[i][vz]	
	    #rescale to have same kinetic energy
	    c = Kold/Knew
	    c = math.sqrt(c)
	    for i in range(natoms):
		    raw[i][vx] *= c
		    raw[i][vy] *= c
		    raw[i][vz] *= c

	    #and now print this out again		
            for i in range(natoms):
		ncolumns=len(raw[i])
		for j in range(ncolumns):
			if j==vx or j==vy or j==vz:
				of.write('%.8f' % raw[i][j])
			else:
				of.write(raw[i][j])
			of.write(" ")
		of.write('\n')

            of.close()  

        
