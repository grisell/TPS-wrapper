'''
Created on Mar 2, 2010

@author: wolf
'''
import os
from random import shuffle


class helpers(object):
    def __init__(self):
        self.wolf = "helpers"
    
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

        
