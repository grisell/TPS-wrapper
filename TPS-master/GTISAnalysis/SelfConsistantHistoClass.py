'''
Created on Jan 5, 2015

@author: grisell (05.01.2015) based on wolf version 1.0 (06.07.2009)
'''
'''

This code solves the problem of self consistent Histograms 
according to Frenkel&Smit Page 186.


'''
import sys
import os
from math import log,fabs, exp
import os.path
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.analysis.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

import pygromacstps.kernels as kernels
import pygromacstps.orderparameters as orderparameters

#Run from main directory, results are created in directory 'tis/selfconsistent' 
ospath = os.getcwd()
basedir = os.path.join(os.getcwd())
ospath = os.path.join(ospath,"tis")

class THisto(object):
    def __init__(self,index,interface,histosize,precision,cutoff,temperature,subanalysis,forward): #GDL:initializations
        self.HISTOSIZE = histosize
        self.PRECISION = precision
        self.analysis = subanalysis
        self.CUTOFF = cutoff
        self.temperature = temperature
        self.kernels = kernels.kernels("head")
        self.forward = forward
        self.interface = interface
        self.index = index
        self.histodata = [ 0 ]*self.HISTOSIZE
        self.lambdahisto = [ 0 ]*self.HISTOSIZE
        self.lambdaboundhisto=[ 0 ]*self.HISTOSIZE
        self.boundaryhistodata = [0]*self.HISTOSIZE
        self.Z = 1.0
        self.lnZ = 0.0
        self.newZ = 1.0
        self.lnnewZ = 0.0
        self.readHisto()
        self.readBoundaryHisto()
        self.M = 1.0
        self.integral = 0.0
        self.max= 0
        self.expW=[0.0]*self.HISTOSIZE
       
        for i in range(self.HISTOSIZE): #maximum value of the histogram data
            if  self.histodata[i] > self.max :
                self.max = self.histodata[i]
        for i in range(self.HISTOSIZE):#Compute normalization factor ignoring contribution from histogram values smaller than certain cutoff percentage 
            if self.histodata[i] < self.CUTOFF*self.max :
                self.histodata[i] = 0

            self.integral += self.histodata[i]        
        for i in range(self.HISTOSIZE): #weights for self consistent algorithm
            if self.histodata[i] > 0:
                self.expW[i]=1.0
            else:
                self.expW[i]=0.0
        
    
    def readHisto(self): #GDL:Read histograms to match

        dirstring = self.forward and "%02d" % self.index or "%02d" % self.interface
        
        if self.analysis=="freeenergy":
            filename=self.forward and os.path.join(basedir,"tis","la","f"+dirstring,"freeenergyhisto.%d.txt" % (self.index*2,)) or os.path.join(basedir,"tis","la","b"+dirstring,"freeenergyhisto.%d.txt" % (self.interface*2,))
        elif self.analysis=="crossinghisto":
            filename=self.forward and os.path.join(basedir,"tis","la","f"+dirstring,"crossinghisto.%d.txt" % (self.index*2,)) or os.path.join(basedir,"tis","la","b"+dirstring,"crossinghisto.%d.txt" % (self.interface*2,))
  
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[2]) #GDL:Get histograms without normalization 
            z = float(line.split()[0]) #GDL: Get order parameter 
            self.lambdahisto[x]=z
            self.histodata[x] = y
            c+= 1
    def readBoundaryHisto(self): #GDL:Read boundary histograms to obatin CA and CB factors
        dirstring = self.forward and "%02d" % self.index or "%02d" % self.interface
        filename=self.forward and os.path.join(basedir,"tis","la","f"+dirstring,"boundaryhisto.%d.txt" % (self.index*2,)) or os.path.join(basedir,"tis","la","b"+dirstring,"boundaryhisto.%d.txt" % (self.interface*2,))
  
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[2]) #Get histograms without normalization 
            z = float(line.split()[0]) # Get order parameter 
            self.lambdaboundhisto[x]=z
            self.boundaryhistodata[x] = y
            c+= 1
            
class TSelfConsistent(object):
    def __init__(self,histosize,precision,cutoff,temperature,analysis):
        self.fhr = []
        self.bhr = []
        self.orderparameters = orderparameters.orderparameters()
        self.HISTOSIZE = histosize
        self.PRECISION = precision
        self.CUTOFF = cutoff
        self.temperature = temperature
        self.analysis = analysis

        
        if self.analysis=='freeenergy':
            self.path=os.path.join(ospath,"selfconsistent_fe")
        elif self.analysis=='crossinghisto':
            self.path=os.path.join(ospath,"selfconsistent")
        if not os.path.exists(self.path):
            os.makedirs(self.path)
         #read the stables states from a file
        self.orderparameters.readOP(os.path.join(basedir,"options","orderparameters.txt"))  

        for line in open(os.path.join(basedir,"options","interfaces.txt")):
            raw = line.split()
            if raw[0] == "f":
                self.fhr.append(float(raw[2]))
            else:
                if raw[0] == "b":
                    self.bhr.append(float(raw[2]))
         
        self.fhistos = []
        self.bhistos = []
        self.fmab = {}
        self.bmab = {}
            
    def newZperreplica(self,i,forward):  #GDL: calculate new Z per replica
        newZ = 0

        ihisto = forward and self.fhistos[i] or self.bhistos[i]
        rangehisto = forward and self.fhistos or self.bhistos
        for Q in range(self.HISTOSIZE):  #GDL: Loop over the bins of each histogram        
          if(ihisto.histodata[Q]>0):
             nominator = 0.0
             denominator = 0.0
             
             for histo in rangehisto:   #GDL: Loop over the replicas 
               		
                nominator += histo.histodata[Q]
                denominator += histo.expW[Q]*histo.integral/histo.Z
               
             if denominator > 0.0 :
             #if denominator > 0.0 and ihisto.expW(Q) >0.0:
                newZ += nominator/denominator                
        return newZ

    def printZvalues(self,forward): #GDL:print values of Z and lnZ when PRECISION is reached 
        dirname = forward and "zvaluesF.txt" or "zvaluesB.txt"
        fullfilename = os.path.join(self.path,dirname)
        of = open(fullfilename,"w")
        rangehistos = forward and self.fhistos or self.bhistos
        for histo in rangehistos:
            of.write("%.18f %.18f\n" % (histo.Z, histo.lnZ))
        of.close()    
        
    def printP_est(self,forward): #print final histogram
        dirname = forward and "pestF.txt" or "pestB.txt"
        
        fullfilename = os.path.join(self.path,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        p_list = [ 0 ]*self.HISTOSIZE        
        norm =0.0
        for Q in range(self.HISTOSIZE):
            nominator = 0.0
            denominator = 0.0
            for histo in tr:
                nominator += histo.histodata[Q]
                denominator += histo.expW[Q]*histo.integral/histo.Z
            if nominator > 0.0 and denominator > 0.0:
                p_est = nominator/denominator
                p_list[Q] = p_est
                if p_est > norm:
                    norm=p_est
                
        for Q in range(self.HISTOSIZE):
            p_est = p_list[Q]/norm
            p_list[Q] = p_est
            if p_est > 0:
                of.write("%f %.18f %.18f\n" % (histo.lambdahisto[Q],p_est,log(p_est)))
        of.close() 
        return p_list
    

    def createallhistos(self,forward): #GDL: create all the histograms
        
        tr = forward and self.fhistos or self.bhistos
        for histo in tr:
            dirname = forward and "pesthistoF.%02d.txt" % (histo.interface,) or "pesthistoB.%02d.txt" % (histo.interface,)
            fullfilename = os.path.join(self.path,dirname)
            of = open(fullfilename,"w")
        
            for Q in range(self.HISTOSIZE):                 
                if histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral>0.0:
                    of.write("%d %.18f %.18f\n" % (histo.lambdahisto[Q],log(histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral),histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral))
        of.close()
    
    def printWeights(self,numinter,forward): #print wham weights
        dirname = forward and "forwardweights.txt" or "backwardweights.txt"
        
        fullfilename = os.path.join(self.path,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        wmax=-1
        w_inter=[0]*numinter
        c=0
        for histo in tr:
	    #print histo.Z
            w_inter[c]=histo.Z/histo.integral
            if w_inter[c] > wmax:
                wmax= w_inter[c]
            c +=1
        c=0
        for histo in tr:
            w_inter[c] /= wmax
            lnw_inter= log(w_inter[c])
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnw_inter,w_inter[c],histo.integral))
            c +=1
        of.close()
    
    def printforwardpathWeights(self,numinter,forward): #print path-weights for forward ensemble
        dirname = forward and "fpathweights.txt" or "bpathweights.txt"
        fullfilename = os.path.join(self.path,dirname)
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
        iself=0
        
        for histo in tr:  #calculate path weights
            w_path[iself]=0.0
            for i in range(iself+1): 
                w_path[iself] +=1.0/w_inter[i]
            w_path[iself] = 1.0/w_path[iself]
            if w_path[iself] > wpathmax:
                wpathmax = w_path[iself]
            iself +=1
        
        c=0
        for histo in tr:
            w_path[c] /= wpathmax
            lnwpath = log(w_path[c])
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnwpath,w_path[c],histo.integral))
            c +=1
        of.close()

    def printBackwardpathWeights(self,numinter,forward): #print path weights for backward ensemble
        dirname = forward and "fpathweights.txt" or "bpathweights.txt"
        fullfilename = os.path.join(self.path,dirname)
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
        iself=0
        
        for histo in tr:  #calculate path weights
            w_path[iself]=0.0
            for i in range(iself,numinter): 
                w_path[iself] +=1.0/w_inter[i]
            w_path[iself] = 1.0/w_path[iself]
            if w_path[iself] > wpathmax:
                wpathmax = w_path[iself]
            iself +=1
        
        c=0
        for histo in tr:
            w_path[c] /= wpathmax
            lnwpath = log(w_path[c])
            of.write("%d %.18f %.18f %d\n" % (histo.interface,lnwpath,w_path[c],histo.integral))
            c +=1
        of.close()

    def testcalculatefactors(self,inter,numinter,forward):
        
        dirname = forward and "ffactor.txt" or "bfactor.txt"
        fullfilename = os.path.join(ospath,dirname)
        of = open(fullfilename,"w")
        #ihisto = forward and self.fhistos[i] or self.bhistos[i]
        s_factor=[0.0]*self.HISTOSIZE
        rangehisto = forward and self.fhistos or self.bhistos
        c=0        
        for histo in rangehisto:
            c += 1
            for Q in range(self.HISTOSIZE):
                #if c < 20:
                    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                        s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                #else:
                #    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                #        s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
        for Q in range(self.HISTOSIZE): 
            sfactor=s_factor[Q]
            of.write("%f %f\n" % (histo.lambdahisto[Q], sfactor))
        return s_factor
    
    def calculatefactors(self,plambda,inter,numinter,states,overlap,forward): #calculate scaling factors s_AB and s_BA (JCP,133,174109-11 (2010))               
        dirname = forward and "sfactor_f.dat" or "sfactor_b.dat"
        fullfilename = os.path.join(self.path,dirname)
        of = open(fullfilename,"w")
        s_factor=[0.0]*self.HISTOSIZE
        cfactor=[0.0]*self.HISTOSIZE
	#avc=0.0
        rangehisto = forward and self.fhistos or self.bhistos
        c=0 
        b=0
        wf=0.0
        Nf=0
        novinter=len(overlap)
        #print plambda, inter, numinter, states , overlap, forward
        for histo in rangehisto:
            c += 1
            #qstarf=int(round((float(HISTOSIZE)*(float(inter[c-1]) - float(HISTOMIN)))/(float(HISTOMAX-HISTOMIN))))
            #qstarb=int(round((float(HISTOSIZE)*(float(inter[b]) - float(HISTOMIN)))/(float(HISTOMAX-HISTOMIN))))
            #print qstarf, histo.lambdahisto[qstarf], plambda[qstarf]
            
            for Q in range(self.HISTOSIZE):
                #print c, len(inter)
                if forward:
                    #print plambda[Q]
                    if c < len(inter):
                        if histo.lambdahisto[Q] >= inter[c-1] and  histo.lambdahisto[Q] < inter[c]:
                            if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                                s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                                if plambda[Q]>0.0:
                                #if plambda[qstarf]>0.0:
                                    cfactor[Q]= s_factor[Q]/plambda[Q] 
                                    #cfactor[Q]= s_factor[Q]/plambda[qstarf]
                    else:
                        if histo.lambdahisto[Q] <= states[0].B[0] :
                            if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                                s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                                if plambda[Q]>0.0:
                                #if plambda[qstarf]>0.0:
                                    cfactor[Q]= s_factor[Q]/plambda[Q]
                                    #cfactor[Q]= s_factor[Q]/plambda[qstarf] 
                else:
                    if b>0:
                        if histo.lambdahisto[Q] < inter[b] and histo.lambdahisto[Q] >= inter[b-1]:
                            if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]> 0.0: 
                                s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                                if plambda[Q]>0.0:
                                #if plambda[qstarb]>0.0:
                                    cfactor[Q]= s_factor[Q]/plambda[Q] 
                                    #cfactor[Q]= s_factor[Q]/plambda[qstarb]
                    else:
                        if histo.lambdahisto[Q] >= states[0].A[1]:
                            if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                                s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                                if plambda[Q]>0.0:
                                #if plambda[qstarb]>0.0:
                                    cfactor[Q]= s_factor[Q]/plambda[Q]
                                    #cfactor[Q]= s_factor[Q]/plambda[qstarb]
            b +=1
          
        for Q in range(self.HISTOSIZE): 
            sfactor=s_factor[Q]
            of.write("%f %f\n" % (histo.lambdahisto[Q], sfactor))
        
        return s_factor

    #def calculatefactors(self,inter,numinter,forward): #calculate scaling factors s_AB and s_BA (JCP,133,174109-11 (2010))               
     #   dirname = forward and "ffactor.txt" or "bfactor.txt"
     #   fullfilename = os.path.join(path,dirname)
     #   of = open(fullfilename,"w")
     #   s_factor=[0.0]*HISTOSIZE
     #   rangehisto = forward and self.fhistos or self.bhistos
     #   c=0 
     #   b=0
     #   for histo in rangehisto:
     #       c += 1
     #       for Q in range(HISTOSIZE):
     #           if forward:
     #               if c < len(inter):
     #                   if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0 and histo.lambdahisto[Q] < inter[c]:
     #                       s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                    #else:
                     #   if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                      #      s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
    #            else:
    #                if b>0:
    #                    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0 and histo.lambdahisto[Q] > inter[b-1]: 
    #                        s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                    #else:
                     #    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                      #      s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
    #        b +=1
                       
    #    for Q in range(HISTOSIZE): 
    #        sfactor=s_factor[Q]
    #        of.write("%f %f\n" % (histo.lambdahisto[Q], sfactor))
     #   return s_factor
                
                
    def freeEnergy(self,p0f,p0b,fmab,bmab,overlap,histo): #Compute free energy projected on order parameter used in TIS
        
        ca=[0.0]*self.HISTOSIZE
        cb=[0.0]*self.HISTOSIZE
        p=[0.0]*self.HISTOSIZE
        pfinal=0.0
        maxp=-1.0
        avca=0.0
        avcb=0.0
        wf=0.0
        wb=0.0
        N=0
        novinter=len(overlap)
        eV=float(0.00013806488*self.temperature)/float(1.60217657)
        filename = os.path.join(self.path,"cacb.dat")
        cval = open(filename,"w")
        for Q in range(self.HISTOSIZE):
          if histo.lambdahisto[Q]>=overlap[0] and histo.lambdahisto[Q]<=overlap[novinter-1]:
            if fmab[Q]>0.0 and bmab[Q]>0.0:
                if p0f[Q]>0.0 and p0b[Q]>0.0:
                    ca[Q]=fmab[Q]/p0f[Q]
                    cb[Q]=bmab[Q]/p0b[Q]
                    wf +=ca[Q]
                    wb +=cb[Q]
                    N +=1
                #p = ca[Q]*p0f[Q] +cb[Q]*p0b[Q]
                    avca= wf/N
                    avcb= wb/N
            #if p > 0.0:
                    cval.write("%d %f %f %f %f\n" % (Q,ca[Q],cb[Q],avca,avcb))
        cval.close()
        fullfilename = os.path.join(self.path,"freeEnergy.dat")
        of = open(fullfilename,"w")
        sA = self.orderparameters.op[0].A[1]
        sB = self.orderparameters.op[0].B[0]
        for Q in range(self.HISTOSIZE):
            p[Q]=avca*p0f[Q]+avcb*p0b[Q]
            #p[Q]=avcb*p0b[Q]
        for Q in range(self.HISTOSIZE):
            if p[Q] > maxp:
                maxp=p[Q] 
        for Q in range(self.HISTOSIZE):
            pfinal=p[Q]
            pfinal /= maxp
            if pfinal > 0.0:
                if ((histo.lambdahisto[Q]>=sA) and (histo.lambdahisto[Q]<=sB)):
                    of.write("%f %f\n" % (histo.lambdahisto[Q],-log(pfinal)*eV))
        of.close()

####################################################################################################################

#if __name__=="__main__":
    def run_main(self):
    	
    
    	for i in range(len(self.fhr)): 
    	    self.fhistos.append(THisto(i,i,self.HISTOSIZE,self.PRECISION,self.CUTOFF,self.temperature,self.analysis,forward=True))

    	inter=len(self.fhr)
    	for i in range(len(self.bhr)):        
    	    self.bhistos.append(THisto(i,inter,self.HISTOSIZE,self.PRECISION,self.CUTOFF,self.temperature,self.analysis,forward=False))
    	    inter += 1

    	ovinter=[]
    	for i in range(len(self.fhr)): 
    	    fvalue=self.fhr[i]
    	    for j in range(len(self.bhr)):
    	        bvalue=self.bhr[j]
    	        if fvalue==bvalue:
    	            ovinter.append(fvalue)
    
#GDL: Self consistent calculation of Z (AB)
    	diff=100.0
    	iteration = 0
    	while (diff > self.PRECISION)&(iteration < 100000) :
            for i,histo in enumerate(self.fhistos):
                self.fhistos[i].newZ = self.newZperreplica(i,True) 
                self.fhistos[i].lnnewZ = log(self.fhistos[i].newZ)
            diff=0.0
            for i,histo in enumerate(self.fhistos):
                diff += fabs(self.fhistos[i].lnnewZ - self.fhistos[i].lnZ)
                self.fhistos[i].lnZ = (self.fhistos[i].lnnewZ - self.fhistos[0].lnnewZ)
                self.fhistos[i].Z = exp(self.fhistos[i].lnZ)
            ilist = []
            for histo in self.fhistos:
                ilist.append(histo.Z)
            iteration += 1 
    	        #print diff,ilist, iteration
    	#for i,histo in enumerate(self.fhistos):
    	#        print self.fhistos[i].Z
#GDL	:print final histograms and weights of AB ensemble
    	finter = self.fhr
    	numinter= len(self.fhr)
    	p0f = self.printP_est(forward = True)
    	self.createallhistos(forward= True)
    	self.printZvalues(forward = True)
    	self.printWeights(numinter,forward = True)
    	self.printforwardpathWeights(numinter,forward= True)
    	#fmab=self.calculatefactors(finter,numinter,forward=True)
    	fmab=self.calculatefactors(p0f,finter,numinter,self.orderparameters.op,ovinter,forward=True)

#GDL: Self consistent calculation of Z (BA)
    	diff=100.0
    	iteration = 0
    	while (diff > self.PRECISION)&(iteration < 100000) :
            for i,histo in enumerate(self.bhistos):
                self.bhistos[i].newZ = self.newZperreplica(i,False) 
                self.bhistos[i].lnnewZ = log(self.bhistos[i].newZ)
            diff=0.0
            for i,histo in enumerate(self.bhistos):
                diff += fabs(self.bhistos[i].lnnewZ - self.bhistos[i].lnZ)
                self.bhistos[i].lnZ = (self.bhistos[i].lnnewZ - self.bhistos[0].lnnewZ)
                self.bhistos[i].Z = exp(self.bhistos[i].lnZ)
                #self.fhistos[0].Z =1.0
            ilistb = []
            for histo in self.bhistos:
                ilistb.append(histo.Z)
            iteration += 1 
    	#for histo in self.bhistos:
    	#        print self.bhistos[i].Z	
#GDL	:print final histograms and weights of BA ensemble
    	binter = self.bhr
    	p0b = self.printP_est(forward = False)
    	self.createallhistos(forward= False)
    	self.printZvalues(forward = False)
    	bnuminter= len(self.bhr)
    	self.printWeights(bnuminter,forward = False)
    	self.printBackwardpathWeights(bnuminter,forward= False)
    	#bmab=self.calculatefactors(binter,bnuminter,forward=False)
    	bmab=self.calculatefactors(p0b,binter,bnuminter,self.orderparameters.op,ovinter,forward=False)
	
	
    	self.freeEnergy(p0f,p0b,fmab,bmab,ovinter,histo)
	
	
    	
    	
	
    		
