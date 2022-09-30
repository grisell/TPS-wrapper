'''
Created on May 7, 2010

@author: wolf
'''

import os
import orderparameter.dist
from string import upper
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False   

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
                        try:
                            (float(raw[2]) and float(raw[3]) and float(raw[4]) and float(raw[5]))
                        except ValueError:
                            logging.error("order parameter value should be a number. Program exiting")
                            raise SystemExit()
                        
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
                        #_tempimp = __import__('custom.'+raw[1],globals(), locals(), ['TCustomOrderParameter'], -1)
                        _tempimp = __import__(raw[1],globals(), locals(), ['TCustomOrderParameter'], -1)
                        try:
                            (float(raw[2]) and float(raw[3]) and float(raw[4]) and float(raw[5]))
                        except ValueError:
                            logging.error("order parameter value should be a number. Program exiting")
                            raise SystemExit()
			            
                        minA = float(raw[2])
                        maxA = float(raw[3])
                        minB = float(raw[4])
                        maxB = float(raw[5])
    
                        newop = _tempimp.TCustomOrderParameter(minA,maxA,minB,maxB)
                        self.op.append(newop)
                    
                        
                
                        
                
