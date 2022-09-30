'''
Created on May 7, 2010

@author: wolf
'''

import os

class interfaces(object):
    def __init__(self):
        ''' 
        forward and backward interfaces
        read from interfaces.txt in the basedir
        '''
        self.interfaces = [[],[]]
        self.ninterfaces = [0,0]
        self.reversepaths = [[],[]]
        
    def readInterfaces(self,filename):
        fn = 0
        bn = 0
        if os.path.exists(filename):
            for line in open(filename,"r"):
                raw = line.split()
                if len(raw) > 0:
                    if raw[0] == "f":
                        pos = float(raw[2])
                        self.interfaces[0].append(pos)  
                        fn += 1        
                    if raw[0] == "b":
                        pos = float(raw[2])
                        self.interfaces[1].append(pos)
                        bn += 1
                    if raw[0] == "A":
                        pos = float(raw[2])
                        self.reversepaths[0].append(pos)
                    if raw[0] == "B":
                        pos = float(raw[2])
                        self.reversepaths[1].append(pos)
                        
        
        self.ninterfaces[0] = fn
        self.ninterfaces[1] = bn
        self.ninter = fn+bn
        
        
        
                