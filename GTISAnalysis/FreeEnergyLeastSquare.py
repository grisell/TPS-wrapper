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
ospath = "/home/users/diazlgjj/tis_systems/Nickel/tis/dT_20/preparations/setup-01/tis/la/cross_histos"
#ospath = "/Users/sopudaww/TPS_wrapper/LammpsTPS-master/example-lammps/test3-tis"
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
        
        ifilename = self.forward and "forward/freeenergyhisto.%02d.dat" % (index*2,) or "backward/freeenergyhisto.%02d.dat" % (index*2,)
        tfilename = os.path.join(ospath,ifilename)
        of = open(tfilename,"w")
        for i in range(HISTOSIZE):
            of.write("%d %f\n" % (i,self.expW(i)*self.histodata[i]/self.M))
        of.close()
        
        
    def expW(self,Q):
        return 1
    
        value = 0.0
        if self.forward:
            if (self.histodata[Q] >= 0):# and (Q>=HMIN) and (Q<=HMAX):
                value = 1.0
        else:
            if (self.histodata[Q] <= 0):# and (Q>=HMIN) and (Q<=HMAX):
                value = 1.0
        return value
#    
    def readHisto(self):
        filename = self.forward and os.path.join(ospath,"forward/freeenergyhisto.%d.txt" % (self.index*2,)) or os.path.join(ospath,"backward/freeenergyhisto.%d.txt" % (self.index*2,))       
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[1])
            self.histodata[x] = y
            c+= 1
            
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
        

    def getInitialFlux(self):
        self.flux = 0.03
        
    def getmab(self):
        fullfilename = os.path.join(ospath,"mabF.dat")
        for line in open(fullfilename,"r"):
            raw = line.split()
            self.fmab[int(raw[0])] = [float(raw[1]),float(raw[2])]

        fullfilename = os.path.join(ospath,"mabB.dat")
        for line in open(fullfilename,"r"):
            raw = line.split()
            
            self.bmab[int(raw[0])] = [float(raw[1]),float(raw[2])]
        
    def Z(self,i,forward):
        mZ = 0

        ihisto = forward and self.fhistos[i] or self.bhistos[i]
        rangehisto = forward and self.fhistos or self.bhistos
        for Q in range(HMIN,HMAX):
            nominator = 0.0
            denominator = 0.0
            #print Q,[x.interface for x in self.fhistos[i-1:i+1]]
            for histo in rangehisto:
                nominator += histo.histodata[Q]
                denominator += histo.expW(Q)*histo.M/histo.Z

            if denominator > 0.0 and ihisto.expW(Q) >0.0:
                mZ += nominator/denominator   
        return mZ
    
    def printP_est(self,forward):
        dirname = forward and "pestF.txt" or "pestB.txt"
        
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        p = []
        for Q in range(HMAX):
            p.append(0.0)
        for Q in range(HMIN,HMAX):
            nominator = 0.0
            denominator = 0.0
            for histo in tr:
                nominator += histo.histodata[Q]
                denominator += histo.expW(Q)*histo.M/histo.Z
            if denominator > 0.0 and nominator > 0.0:
                p_est = nominator/denominator
                p[Q] = p_est
                of.write("%d %.18f %.18f\n" % (Q,p_est,log(p_est)))
            
        for histo in tr:
            dirname = forward and "pesthistoF.%02d.txt" % (histo.interface,) or "pesthistoB.%02d.txt" % (histo.interface,)
            fullfilename = os.path.join(ospath,dirname)
            of = open(fullfilename,"w")
        
            for Q in range(HMIN,HMAX):
                if histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.M>0.0:
                    of.write("%d %.18f %.18f\n" % (Q,log(histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.M),histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.M))
            of.close()
        return p
    
    def printWeights(self,forward):
        dirname = forward and "forwardweights.txt" or "backwardweights.txt"
        
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        for histo in tr:
            of.write("%d %.18f\n" % (histo.interface,histo.Z))
        of.close()
    
    def freeEnergy(self,p0f,p0b):
        ti = 100
        if sim.fmab[ti][0]*p0f[ti] > 0.0 and sim.bmab[ti][0]*p0b[ti] > 0.0:
            ffactor = sim.fmab[ti][1]/(sim.fmab[ti][0]*p0f[ti])
            bfactor = sim.bmab[ti][1]/(sim.bmab[ti][0]*p0b[ti])
            fullfilename = os.path.join(ospath,"freeEnergy.%02d.dat" % (ti,))
            of = open(fullfilename,"w")
            for Q in range(HMAX):
                p = ffactor*p0f[Q] +bfactor*p0b[Q]
                if p > 0.0:
                    of.write("%d %f\n" % (Q,-log(p)))
            of.close()

if __name__=="__main__":
    sim = TSelfConsistent()

    for i in range(len(sim.fhr)):
        inter = sim.fhr[i]
        sim.fhistos.append(THisto(i,i,forward=True))
    
#    for i in range(len(sim.bhr)):
#        inter = sim.bhr[i]
#        if i>1:
#            ninter = sim.bhr[i-1]
#        else:
#            ninter = HMIN
#        sim.bhistos.append(THisto(inter,ninter,forward=False))
#    
#    
    for j in range(1000):
        oldhistoz = sim.fhistos[1].Z
        
        for i,histo in enumerate(sim.fhistos):
            sim.fhistos[i].Z = sim.Z(i,True) 
        sim.fhistos[0].Z =1.0
        ilist = []
        diff = sim.fhistos[1].Z - oldhistoz
#        if fabs(diff) < PRECISION:
#            break
        if j%1 == 0:
            for histo in sim.fhistos:
                ilist.append(histo.Z)
            print fabs(diff),ilist
    p0f = sim.printP_est(forward = True)
    sim.printWeights(forward = True)
    
#    print "Backward"
#    for j in range(100000):
#        oldhistoz = sim.bhistos[-2].Z
#        for i,histo in enumerate(sim.bhistos):
#            if histo.interface < 180:
#                sim.bhistos[i].Z = sim.Z(i,False) 
#             
#        #sim.bhistos[-1].Z =1.0
#        ilist = []
#        diff = sim.bhistos[-2].Z - oldhistoz
#        if fabs(diff) < PRECISION:
#            break
#        if j%100 == 0:
#            for histo in sim.bhistos:
#                ilist.append(histo.Z)
#            print fabs(diff),ilist
#    p0b = sim.printP_est(forward = False)
#    sim.printWeights(forward = False)
#    

    
    

    
