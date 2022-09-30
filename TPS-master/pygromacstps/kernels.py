'''
Created on May 7, 2010

@author: wolf
'''

import os
from math import ceil
from itertools import islice
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 


class kernels(object):
    def __init__(self,kernel):
        self.nkernels = 1
        self.thiskernel = kernel
        self.kernelPaths = []
        self.kernelPathsAll = []
        self.ntps = 1
        self.dynamicshooting = 0
        
    def readKernelOptions(self,filename):
        if os.path.exists(filename):
            for line in open(filename,"r"):
                raw = line.split("=")
                if len(raw) > 0:
                    if raw[0] == "nkernels":
                        try:
                            float(raw[1])
                        except ValueError:
                            logger.error("nkernels value should be a number. Program exiting")
                            raise SystemExit()
                        self.nkernels = int(float(raw[1]))
                    if raw[0] == "dynamicshooting":
                        try:
                            float(raw[1])
                        except ValueError:
                            logger.error("dynamicshooting value should be a number. Program exiting")
                            raise SystemExit()
                        self.dynamicshooting = int(float(raw[1]))
                    if raw[0] == "tpsparallel":
                        try:
                            float(raw[1])
                        except ValueError:
                            logger.error("tpsparallel value should be a number. Program exiting")
                            raise SystemExit()
                        self.ntps = int(float(raw[1]))
        else:
            logger.error("kerneloptions.txt not found")
            raise SystemExit()
                        
        
                
    def splitNodes(self,allPaths, size):
        
        it = iter(allPaths)
        item = list(islice(it, size))
        while item:
            yield item
            item = list(islice(it, size))
    
    def generateKernelLists(self,npaths):
        
        if self.thiskernel == "head":
            self.kernelPathsAll = []
            for i in range(2*npaths):
                self.kernelPathsAll.append(i)
            self.kernelPaths = []
            for i in range(npaths):
                self.kernelPaths.append(2*i)            
            
        elif self.thiskernel == "reverse":
            pass
        else:
            if self.nkernels <= npaths:
                npk = ceil(float(npaths)/float(self.nkernels))
            else:
                npk = 1
            
            tempklist = list(self.splitNodes([i * 2 for i in range(npaths)],npk))
            self.nkernels = len(tempklist)
            self.kernelPaths = tempklist[self.thiskernel]
            self.kernelPathsAll = []
            for p in self.kernelPaths:
                self.kernelPathsAll.append(p)
                self.kernelPathsAll.append(p+1)
                
    
            
        
                
