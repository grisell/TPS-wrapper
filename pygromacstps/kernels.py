'''
Created on May 7, 2010

@author: wolf
'''

import os
from math import ceil
from itertools import islice

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
                        self.nkernels = int(float(raw[1]))
                    if raw[0] == "dynamicshooting":
                        self.dynamicshooting = int(float(raw[1]))
                    if raw[0] == "tpsparallel":
                        self.ntps = int(float(raw[1]))
                        
        
                
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
#	    print "tempklist\n"
#	    print tempklist
            self.nkernels = len(tempklist)
            self.kernelPaths = tempklist[self.thiskernel]
            self.kernelPathsAll = []
            for p in self.kernelPaths:
                self.kernelPathsAll.append(p)
                self.kernelPathsAll.append(p+1)
                
    
            
        
                
