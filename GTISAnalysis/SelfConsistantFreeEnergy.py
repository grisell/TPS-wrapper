'''
Created on Feb 3, 2011

@author: wolf
'''
'''
Created on 06.07.2009

@author: Wolfgang Lechner

This solves the problem of self consistent Histograms 
according to Frenkel&Smit Page 186

'''

import os
from math import log,fabs
from random import random

ospath = "/home/users/diazlgjj/tis_systems/Nickel/tis/dT_20/preparations/setup-01/tis/la/cross_histos"
#ospath = "/home/cpd/wolf/TIS/MD/status"

HISTOSIZE = 2000
HMIN = 0
HMAX = 2000
PRECISION = 1e-10
class THisto(object):
    def __init__(self,index,interface,forward):
        self.forward = forward
        self.interface = interface
        self.index = index
        self.histodata = []
        self.start = 0
        self.Z = 1.0
        self.TZ = 1.0
        for i in range(HISTOSIZE):
            self.histodata.append(0)
        self.readHisto()
        self.M = 1.0
        self.histocut = -1
        integral = 0
        for i in range(HISTOSIZE):
            integral += self.histodata[i]
        
        integral = 0.0
        for i in range(HISTOSIZE):
            integral += float(self.histodata[i])

        for i in range(HISTOSIZE):
            self.histodata[i] /= float(integral)
        
#    
    def readHisto(self):
        filename = self.forward and os.path.join(ospath,"forward/freeenergyhisto.%d.txt" % (self.index*2,)) or os.path.join(ospath,"backward/freeenergyhisto.%d.txt" % (self.index*2,))       
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[1])
            if y > 0.0:
                self.histodata[x] = -log(y)
                if self.start ==0 :
                    self.start = x
            else:
                self.histodata[x] = 0.0
            c+= 1
            
    def printHisto(self):
        ifilename = self.forward and "forward/freeenergyhistolog.%02d.dat" % (self.index*2,) or "backward/freeenergyhistolog.%02d.dat" % (self.index,)
        tfilename = os.path.join(ospath,ifilename)
        of = open(tfilename,"w")
        for i in range(HISTOSIZE):
            of.write("%d %f\n" % (i,self.histodata[i]))
        of.close()
    
    
class TSelfConsistent(object):
    def __init__(self):
        self.fhr = []
        self.bhr = []
        
        for line in open(os.path.join(ospath,"interfaces.txt")):
            raw = line.split()
            if raw[0] == "f":
                self.fhr.append(float(raw[2]))
            else:
                self.bhr.append(float(raw[2]))
         
        self.fhistos = []
        self.bhistos = []
        self.fmab = {}
        self.bmab = {}
        
    def squaredistance(self,i,j):
        sq = 0.0
        sum = 0.0
        for k in range(HISTOSIZE):
            if k >= self.fhistos[i].start and k< self.fhistos[i].start + 100:
                if self.fhistos[i].histodata[k] > 0.0 and self.fhistos[j].histodata[k] > 0.0:
                    sq += (self.fhistos[i].histodata[k] - self.fhistos[j].histodata[k])**2
                    sum += self.fhistos[i].histodata[k] - self.fhistos[j].histodata[k]
        return sq,sum
    
    
    def moveHisto(self,i,diff):
        for k in range(HISTOSIZE):
            if self.fhistos[i].histodata[k]>0.0:
                d = self.fhistos[i].histodata[k] + diff
                self.fhistos[i].histodata[k] = d
     
    def avgCurve(self):
        allavg = []
        for i in range(HISTOSIZE):
            c = 0.0
            avg = 0.0
            for k in range(len(self.fhr)):
                if self.fhistos[k].histodata[i] > 0.0:
                    avg += self.fhistos[k].histodata[i]
                    c += 1.0
            if c >0.0:
                allavg.append(avg/c)
            else:
                allavg.append(0.0)
                
        return allavg
    
if __name__=="__main__":
    sim = TSelfConsistent()

    for i in range(len(sim.fhr)):
        inter = sim.fhr[i]
        sim.fhistos.append(THisto(i,inter,forward=True))
    
        
    for k in range(len(sim.fhr)-1):
        
        osq,sum = sim.squaredistance(k, k+1)
        sumdiff = 0.0
        for i in range(5000):    
            diff = random()/1000.0
            sim.moveHisto(k+1, diff)
            
            sq,sum = sim.squaredistance(k,k+1)
            if osq > sq:
                osq = sq
                sumdiff += diff
            else:
                sim.moveHisto(k+1, -diff)
        print k, sumdiff
    
    avgcurve = sim.avgCurve()
    
    tfilename = os.path.join(ospath,"avgcurve.txt")
    of = open(tfilename,"w")
    for i in range(HISTOSIZE):
        of.write("%d %f\n" % (i,avgcurve[i]+avgcurve[HISTOSIZE-i-1]))
    of.close()

    for i in range(len(sim.fhr)):
        sim.fhistos[i].printHisto()
        
    
       
    

    
