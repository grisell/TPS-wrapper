'''
Created on May 7, 2010

@author: wolf
'''

import os
import logging


#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 


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
                        try:
                            float(raw[2])
                        except ValueError:
                            logger.error("interface value should be a number. Program exiting")
                            raise SystemExit()
                        pos = float(raw[2])
                        self.interfaces[0].append(pos)  
                        fn += 1
                    elif raw[0] == "b":
                        try:
                            float(raw[2])
                        except ValueError:
                            logger.error("interface value should be a number. Program exiting")
                            raise SystemExit()
                        pos = float(raw[2])
                        self.interfaces[1].append(pos)
                        bn += 1
                    elif raw[0] == "A":
                        try:
                            float(raw[2])
                        except ValueError:
                            logger.error("interface value should be a number. Program exiting")
                            raise SystemExit()
                        pos = float(raw[2])
                        self.reversepaths[0].append(pos)
                    elif raw[0] == "B":
                        try:
                            float(raw[2])
                        except ValueError:
                            logger.error("interface value should be a number. Program exiting")
                            raise SystemExit()
                        pos = float(raw[2])
                        self.reversepaths[1].append(pos)
                    
                    else:
                        logger.error("Unknown interface keyword found. Correct interfaces.txt")
                        raise SystemExit()
                    

        else:
            logger.error("interfaces.txt not found")
            raise SystemExit()
        
        self.ninterfaces[0] = fn
        self.ninterfaces[1] = bn
        self.ninter = fn+bn
        
        
        
                
