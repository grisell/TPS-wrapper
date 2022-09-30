'''
Created on Jan 5, 2015

@author: wolf version 1.0 (06.07.2009)
@author: grisell version 1.0.1 (05.01.2015)
'''
'''

This code solves the problem of self consistent Histograms 
according to Frenkel&Smit Page 186.


'''

import os
from math import log,fabs, exp
ospath = "/home/users/diazlgjj/tis_systems/Nickel/tis/dT_25/final/setup-final/tis/la/free-energy_ns"
#ospath = "/home/cpd/wolf/TIS/MD/status"

HISTOSIZE = 1000
HMIN = 0.0
HMAX = 5.0
PRECISION = 1e-4
CUTOFF=0.01  
class THisto(object):
    def __init__(self,index,interface,forward):
        self.forward = forward
        self.interface = interface
        self.index = index
        self.histodata = [ 0 ]*HISTOSIZE
        self.lambdahisto = [ 0 ]*HISTOSIZE
        self.Z = 1.0
        self.lnZ = 0.0
        self.newZ = 1.0
        self.lnnewZ = 0.0
        #  for i in range(HISTOSIZE):
        #      self.lambdahisto.append(0)
        #      self.histodata.append(0)
        self.readHisto()
        self.M = 1.0
        #self.histocut = -1
        self.integral = 0.0
        self.max= 0 
       
        for i in range(HISTOSIZE):
            if  self.histodata[i] > self.max :
                self.max = self.histodata[i]
        for i in range(HISTOSIZE):
            if self.histodata[i] < CUTOFF*self.max :
                self.histodata[i] = 0

            self.integral += self.histodata[i]

        #for i in range(HISTOSIZE):
        #    self.histodata[i] /= float(integral)
        
        #ifilename = self.forward and "forward/freeenergyhisto.%02d.dat" % (index*2,) or "backward/freeenergyhisto.%02d.dat" % (index*2,)
        #tfilename = os.path.join(ospath,ifilename)
        #of = open(tfilename,"w")
        #for i in range(HISTOSIZE):
        #    of.write("%d %f\n" % (i,self.expW(i)*self.histodata[i]/self.M))
        #of.close()
        
        
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
            y = float(line.split()[2]) #Get histograms without normalization 
            z = float(line.split()[0]) # Get order parameter 
            self.lambdahisto[x]=z
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
        
    def newZperreplica(self,i,forward):  #GDL: calculate new Z per replica
        newZ = 0

        ihisto = forward and self.fhistos[i] or self.bhistos[i]
        rangehisto = forward and self.fhistos or self.bhistos
        for Q in range(HISTOSIZE):  #GDL: Loop over the bins of each histogram        
          if(ihisto.histodata[Q]>0):
             nominator = 0.0
             denominator = 0.0
             #print Q,[x.interface for x in self.fhistos[i-1:i+1]]
             for histo in rangehisto:   #GDL: Loop over the replicas 
               
		if (histo.histodata[Q]>0):
		 weight= 1.0
	       	else:
		 weight = 0.0

                nominator += histo.histodata[Q]
                denominator += weight*histo.integral/histo.Z

             if denominator > 0.0 :
             #if denominator > 0.0 and ihisto.expW(Q) >0.0:
                newZ += nominator/denominator    
            # print newZ
        return newZ
    def printZvalues(self,forward): #GDL:print values of Z and lnZ when PRECISION is reached 
        dirname = forward and "zvaluesF.txt" or "zvaluesB.txt"
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        rangehistos = forward and self.fhistos or self.bhistos
        for histo in rangehistos:
            #print histo.Z, histo.lnZ 
            of.write("%.18f %.18f\n" % (histo.Z, histo.lnZ))
        of.close()    
        
    def printP_est(self,forward):
        dirname = forward and "pestF.txt" or "pestB.txt"
        
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        #p_list =[]
        p_list = [ 0 ]*HISTOSIZE
        #for i in range(HISTOSIZE):
        #     p_list.append(0)
        norm =0.0
        for Q in range(HISTOSIZE):
            nominator = 0.0
            denominator = 0.0
            for histo in tr:
                if (histo.histodata[Q]>0):
		 weight= 1.0
	       	else:
		 weight = 0.0
                nominator += histo.histodata[Q]
                denominator += weight*histo.integral/histo.Z
            if nominator > 0.0 and denominator > 0.0:
                p_est = nominator/denominator
                p_list[Q] = p_est
                if p_est > norm:
                    norm=p_est
        #print p_list
                #of.write("%d %.18f %.18f\n" % (histo.lambdahisto[Q],p_est,log(p_est)))          
        
        for Q in range(HISTOSIZE):
            #print p_list[Q]
            p_est = p_list[Q]/norm
            if p_est > 0:
                of.write("%f %.18f %.18f\n" % (histo.lambdahisto[Q],p_est,log(p_est)))
        of.close() 
        return p_list
    

    def createallhistos(self,forward): #GDL: create all the histograms
        
        tr = forward and self.fhistos or self.bhistos
        for histo in tr:
            dirname = forward and "pesthistoF.%02d.txt" % (histo.interface,) or "pesthistoB.%02d.txt" % (histo.interface,)
            fullfilename = os.path.join(ospath,dirname)
            of = open(fullfilename,"w")
        
            for Q in range(HISTOSIZE):
                if histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.integral>0.0:
                    of.write("%d %.18f %.18f\n" % (histo.lambdahisto[Q],log(histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.integral),histo.expW(Q)*histo.histodata[Q]*histo.Z/histo.integral))
            of.close()
    
    def printWeights(self,numinter,forward):
        dirname = forward and "forwardweights.txt" or "backwardweights.txt"
        
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        wmax=-1
        w_inter=[0]*numinter
        c=0
        for histo in tr:
            w_inter[c]=histo.Z/histo.integral
            #print w_inter[c]
            if w_inter[c] > wmax:
                wmax= w_inter[c]
            c +=1
        c=0
        for histo in tr:
            w_inter[c] /= wmax
            lnw_inter= log(w_inter[c])
            #print histo.interface, lnw_inter, w_inter[c], histo.integral
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnw_inter,w_inter[c],histo.integral))
            c +=1
        of.close()
    
    def printforwardpathWeights(self,numinter,forward):
        dirname = forward and "fpathweights.txt" or "bpathweights.txt"
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        wmax=-1
        w_inter=[0]*numinter
        w_path=[0]*numinter
        wpathmax=-1
        c=0
        for histo in tr: #calculate weights
            w_inter[c]=histo.Z/histo.integral
            #print w_inter[c]
            c +=1
        isim=0
        #print range(numinter)
        
        for histo in tr:  #calculate path weights
            w_path[isim]=0.0
            #print range(isim,numinter)
            for i in range(isim+1): 
                w_path[isim] +=1.0/w_inter[i]
                #print w_inter[i]
            #print w_path[isim]    
            w_path[isim] = 1.0/w_path[isim]
            #print w_path[isim]
            if w_path[isim] > wpathmax:
                wpathmax = w_path[isim]
            isim +=1
        
        c=0
        for histo in tr:
            w_path[c] /= wpathmax
            lnwpath = log(w_path[c])
            #print w_path[c], lnwpath
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnwpath,w_path[c],histo.integral))
            c +=1
            #print histo.interface, lnw_inter, w_inter[c], histo.integral
            #of.write("%d %.18f %.18f %d\n" % (histo.interface,lnwpath,w_path[c],histo.integral))
        of.close()

    def printBackwardpathWeights(self,numinter,forward):
        dirname = forward and "fpathweights.txt" or "bpathweights.txt"
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        wmax=-1
        w_inter=[0]*numinter
        w_path=[0]*numinter
        wpathmax=-1
        c=0
        for histo in tr: #calculate weights
            w_inter[c]=histo.Z/histo.integral
            c +=1
        isim=0
        
        for histo in tr:  #calculate path weights
            w_path[isim]=0.0
            for i in range(isim,numinter): 
                w_path[isim] +=1.0/w_inter[i]
            w_path[isim] = 1.0/w_path[isim]
            if w_path[isim] > wpathmax:
                wpathmax = w_path[isim]
            isim +=1
        
        c=0
        for histo in tr:
            w_path[c] /= wpathmax
            lnwpath = log(w_path[c])
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnwpath,w_path[c],histo.integral))
            c +=1
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
    print len(sim.fhr)
    for i in range(len(sim.fhr)):
        inter = sim.fhr[i]
        sim.fhistos.append(THisto(i,i,forward=True))
    
    for i in range(len(sim.bhr)-1):        
        inter = sim.bhr[i]
       #if i>1:
        #    ninter = sim.bhr[i-1]
       # else:
        #    ninter = HMIN
        sim.bhistos.append(THisto(i,i,forward=False))
#    

#GDL: Self consistent calculation of Z
    diff=100.0
    iteration = 0
    while (diff > PRECISION)&(iteration < 100000) :
                
            for i,histo in enumerate(sim.fhistos):
                sim.fhistos[i].newZ = sim.newZperreplica(i,True) 
                sim.fhistos[i].lnnewZ = log(sim.fhistos[i].newZ)
                #print  sim.fhistos[i].newZ, sim.fhistos[i].lnnewZ
            
            diff=0.0
	    for i,histo in enumerate(sim.fhistos):
                
                diff += fabs(sim.fhistos[i].lnnewZ - sim.fhistos[i].lnZ)
                sim.fhistos[i].lnZ = (sim.fhistos[i].lnnewZ - sim.fhistos[0].lnnewZ)
                sim.fhistos[i].Z = exp(sim.fhistos[i].lnZ)
                #sim.fhistos[0].Z =1.0
            ilist = []
            for histo in sim.fhistos:
                  ilist.append(histo.Z)
            iteration += 1 
            print diff,ilist, iteration
#GDL:
    p0f = sim.printP_est(forward = True)
    sim.createallhistos(forward= True)
    #sim.printWeights(forward = True)
    sim.printZvalues(forward = True)
    numinter= len(sim.fhr)
    #for i,histo in enumerate(sim.fhistos):
    sim.printWeights(numinter,forward = True)
    #for i,histo in enumerate(sim.fhistos):
    sim.printforwardpathWeights(numinter,forward= True)
    



    diff=100.0
    iteration = 0
    while (diff > PRECISION)&(iteration < 100000) :
                
            for i,histo in enumerate(sim.bhistos):
                sim.bhistos[i].newZ = sim.newZperreplica(i,False) 
                sim.bhistos[i].lnnewZ = log(sim.bhistos[i].newZ)
                #print  sim.fhistos[i].newZ, sim.fhistos[i].lnnewZ
            
            diff=0.0
	    for i,histo in enumerate(sim.bhistos):
                
                diff += fabs(sim.bhistos[i].lnnewZ - sim.bhistos[i].lnZ)
                sim.bhistos[i].lnZ = (sim.bhistos[i].lnnewZ - sim.bhistos[0].lnnewZ)
                sim.bhistos[i].Z = exp(sim.bhistos[i].lnZ)
                #sim.fhistos[0].Z =1.0
            ilistb = []
            for histo in sim.bhistos:
                  ilistb.append(histo.Z)
            iteration += 1 
            print diff,ilistb, iteration
#GDL:
    p0b = sim.printP_est(forward = False)
    sim.createallhistos(forward= False)
    #sim.printWeights(forward = True)
    sim.printZvalues(forward = False)
    

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

    
    

    
