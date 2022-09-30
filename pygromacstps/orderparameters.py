'''
Created on May 7, 2010

@author: wolf
'''

import os
import orderparameter.dist
from string import upper

class orderparameters(object):
    def __init__(self):
        self.op = []
        self.optionToGDistDirection = {"ABS":2,"X":3,"Y":4,"Z":5}
    
    def readOP(self,filename):
        if os.path.exists(filename):
            for line in open(filename,"r"):
                if not(line.startswith("\\")):
                    raw = line.split()
                    if upper(raw[0]) == "DIST":
                        
                        minA = float(raw[2])
                        maxA = float(raw[3])
                        minB = float(raw[4])
                        maxB = float(raw[5])
                        
                        if not(raw[1].upper() in self.optionToGDistDirection.keys()):
                            print "dir must be x,y,z or abs\n"
                        else:
                            dircolumn = self.optionToGDistDirection[raw[1].upper()]
                        newop = orderparameter.dist.TDistOrderParameter(minA,maxA,minB,maxB,dircolumn)
                        self.op.append(newop)
                        print minA,maxA,minB,maxB
                    if upper(raw[0]) == "CUSTOM":
                        _tempimp = __import__('custom.'+raw[1],globals(), locals(), ['TCustomOrderParameter'], -1)
			
			minA = float(raw[2])
                        maxA = float(raw[3])
                        minB = float(raw[4])
                        maxB = float(raw[5])
    
                        newop = _tempimp.TCustomOrderParameter(minA,maxA,minB,maxB)
                        self.op.append(newop)
                    
                        
                
                        
                
