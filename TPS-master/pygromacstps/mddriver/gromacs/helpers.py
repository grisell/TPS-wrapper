'''
Created on Mar 2, 2010

@author: wolf
'''
import os
from random import shuffle
import random
import sys
import math


class helpers(object):
    def __init__(self):
        self.wolf = "helpers"
    
    """
    def reverseVelocities(self,filename,outfilename):
        signpositions = [45,53,61]
        if os.path.exists(filename):                
            of = open(outfilename,"w")
            pos = 0
            infile = open(filename,"r")
            alllines = infile.readlines()
            for line in alllines:    
                raw = line.split()
                ar = []
                if pos > 1 and pos<(len(alllines)-1) and (len(raw) == 9 or len(raw) ==8):
                    for k in range(len(line)):
                        if k in signpositions:
                            if line[k] == '-':
                                ar.append(' ')
                            elif line[k] == ' ':
                                ar.append('-')
                            else:
                                print "ERROR: signs in gro file are expeted to be at positions " + str(signpositions)
                        else:
                            ar.append(line[k])
                    outline = ''.join(ar)
                else:
                    outline = line
                pos += 1
                of.write(outline)
            of.close()
    """

    def reverseVelocities(self,filename,outfilename):
        data = [] 

        if os.path.exists(filename):                
            of = open(outfilename,"w")
            infile = open(filename,"r")
        
        for line in infile:
            data.append(line)
        
        nlines = len(data)
        natoms = int(data[1])
        
        for i in range(2):
            of.write(data[i])

        ncolumns=9
        vx=6
        vy=7
        vz=8
        vel_list = [6,7,8]
        form1_list = [1,3,4,5]

        raw = []
        for i in range(2,natoms+2):
            k=i-2
            raw.append(data[i].split())
            #print raw
            tmp = float(raw[k][vx])
            raw[k][vx] = -tmp 
            tmp = float(raw[k][vy])
            raw[k][vy] = -tmp 
            tmp = float(raw[k][vz])
            raw[k][vz] = -tmp
            for j in range(ncolumns):
                if j==0:
                    of.write(("%7s")%(raw[k][j]))
                elif j in form1_list: 
                    of.write(("%8s")%(raw[k][j]))
                elif j==2:
                    of.write(("%5s")%(raw[k][j]))
                if j in vel_list:
                    of.write('%8.4f' % float(raw[k][j]))
            of.write('\n')

        of.write(data[natoms+2])
        of.close()



    def shootVelocities(self,options,filename,outfilename):

        data = []
        dvmax = options.runoptions["dvmax"] 
        #print filename

        if os.path.exists(filename):                
            of = open(outfilename,"w")
            infile = open(filename,"r")
        #else:
        #    print filename
        
        for line in infile:
            data.append(line)
        
        nlines = len(data)
        natoms = int(data[1])
        
        for i in range(2):
            of.write(data[i])

        ncolumns=9
        vx=6
        vy=7
        vz=8
        vel_list = [6,7,8]
        form1_list = [1,3,4,5]

        raw = []
        for i in range(2,natoms+2):
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
            for j in range(ncolumns):
                if j==0:
                    of.write(("%7s")%(raw[i][j]))
                elif j in form1_list: 
                    of.write(("%8s")%(raw[i][j]))
                elif j==2:
                    of.write(("%5s")%(raw[i][j]))
                if j in vel_list:
                    of.write('%8.4f' % float(raw[i][j]))
            of.write('\n')

        of.write(data[natoms+2])
        of.close()

"""
    def shootVelocities(self,filename,outfilename):
        shootpositions = [49,57,65]
        shoots = [-1,-1,-1,-1,1,1,1,1]
        shootpos = 0
        shuffle(shoots)
        shuffle(shootpositions)
        if os.path.exists(filename):                
            of = open(outfilename,"w")
            pos = 0
            infile = open(filename,"r")
            alllines = infile.readlines()
            for line in alllines:    
                raw = line.split()
                if pos == 1:
                    size = int(float(raw[0]))
                    posrange = range(1,size+1)
                    shuffle(posrange)
                    
                ar = []
                if pos > 1 and pos<(len(alllines)-1) and (len(raw) == 9 or len(raw) ==8):
                        
                    for k in range(len(line)):
                        ar.append(line[k])
                        
                    if pos in posrange[:len(shoots)]:
                        n = int(ar[shootpositions[0]])+ shoots[shootpos]
                        if n >= 0 and n<10:
                            ar[shootpositions[0]] = str(n)
                        shootpos += 1
                    outline = ''.join(ar)
                else:
                    outline = line
                pos += 1
                of.write(outline)
            of.close()  
"""

        
