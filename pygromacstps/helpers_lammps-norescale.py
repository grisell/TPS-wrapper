'''
Created on Mar 2, 2010

@author: wolf
'''
import os
import random


class helpers(object):
    def __init__(self):
        self.wolf = "helpers"
    
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
        dvmax = self.options.runoptions["scdvmax"]
	#dvmax = 0.01
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
		raw[k][vx] = tmp + random.gauss(0.0,1.0)*dvmax
		tmp = float(raw[k][vy])
		raw[k][vy] = tmp + random.gauss(0.0,1.0)*dvmax
		tmp = float(raw[k][vz])
		raw[k][vz] = tmp + random.gauss(0.0,1.0)*dvmax

		for j in range(ncolumns):
			if j==vx or j==vy or j==vz:
				of.write('%.8f' % raw[k][j])
			else:
				of.write(raw[k][j])
			of.write(" ")
		of.write('\n')

            of.close()  

        
