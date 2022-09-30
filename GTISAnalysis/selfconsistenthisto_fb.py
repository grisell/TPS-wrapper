'''
Created on Jan 5, 2015

@author: grisell (05.01.2015) 
'''
'''

This code solves the problem of self consistent Histograms (according to Frenkel&Smit Page 186),
Additionaly you can compute:
- The scaling factors from matched histograms (J. Chem. Phys. 133, 174109-11 (2010)) 
- The free energy 1d profile from the free energy histograms obtained from TIS code.  

'''

import os
from math import log,fabs, exp
import pygromacstps.kernels as kernels
import os.path
import argparse

#initial and default settings 
DEFAULT_OUTPUT_DIR = "selfconsistent"
HISTOSIZE = 2000
PRECISION = 1e-4
CUTOFF=0.01 

# Parse the input arguments.
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output-dir", help="Output directory containing weights and final histograms", default=DEFAULT_OUTPUT_DIR)
parser.add_argument('-fe', '--print-freeenergy', dest='do_printfe', action='store_true', default=False, help='print free energy and ca, cb matching factors from freeenergy histograms')
args = parser.parse_args()

#Set paths 
ospath = os.getcwd()
basedir = os.path.join(os.getcwd(),"..")
path=os.path.join(ospath, args.output_dir)
if not os.path.exists(path):
    os.makedirs(path)

#def is_true(value):
 #   return value.lower() == 'true' 

do_printfe = args.do_printfe

class THisto(object):
    def __init__(self,index,interface,forward): #GDL:initializations
        self.kernels = kernels.kernels("head")
        self.forward = forward
        self.interface = interface
        self.index = index
        self.histodata = [ 0 ]*HISTOSIZE
        self.lambdahisto = [ 0 ]*HISTOSIZE
        self.lambdaboundhisto=[ 0 ]*HISTOSIZE
        self.boundaryhistodata = [0]*HISTOSIZE
        self.Z = 1.0
        self.lnZ = 0.0
        self.newZ = 1.0
        self.lnnewZ = 0.0
        self.readHisto()
        self.readBoundaryHisto()
        self.M = 1.0
        self.integral = 0.0
        self.max= 0
        self.expW=[0.0]*HISTOSIZE
       
        for i in range(HISTOSIZE): #maximum value of the histogram data
            if  self.histodata[i] > self.max :
                self.max = self.histodata[i]
        for i in range(HISTOSIZE):#Compute normalization factor ignoring contribution from histogram values smaller than certain cutoff percentage 
            if self.histodata[i] < CUTOFF*self.max :
                self.histodata[i] = 0

            self.integral += self.histodata[i]        
        for i in range(HISTOSIZE): #weights for self consistent algorithm
            if self.histodata[i] > 0:
                self.expW[i]=1.0
            else:
                self.expW[i]=0.0
        
    
    def readHisto(self): #GDL:Read histograms to match

        dirstring = self.forward and "%02d" % self.index or "%02d" % self.interface
        if do_printfe:
            filename=self.forward and os.path.join(basedir,"tis","la","f"+dirstring,"freeenergyhisto.%d.txt" % (self.index*2,)) or os.path.join(basedir,"tis","la","b"+dirstring,"freeenergyhisto.%d.txt" % (self.interface*2,))
        else:
            filename=self.forward and os.path.join(basedir,"tis","la","f"+dirstring,"crossinghisto.%d.txt" % (self.index*2,)) or os.path.join(basedir,"tis","la","b"+dirstring,"crossinghisto.%d.txt" % (self.interface*2,))
  
        c = 0
        for line in open( filename , 'r' ):
            x = c
            y = float(line.split()[2]) #Get histograms without normalization 
            z = float(line.split()[0]) #Get order parameter 
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
    def __init__(self):
        self.fhr = []
        self.bhr = []
        
        
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
            
    def newZperreplica(self,i,forward):  #calculate new Z per replica
        newZ = 0

        ihisto = forward and self.fhistos[i] or self.bhistos[i]
        rangehisto = forward and self.fhistos or self.bhistos
        for Q in range(HISTOSIZE):  #Loop over the bins of each histogram        
          if(ihisto.histodata[Q]>0):
             nominator = 0.0
             denominator = 0.0
             
             for histo in rangehisto:   # Loop over the replicas 
               		
                nominator += histo.histodata[Q]
                denominator += histo.expW[Q]*histo.integral/histo.Z
               
             if denominator > 0.0 :
             #if denominator > 0.0 and ihisto.expW(Q) >0.0:
                newZ += nominator/denominator                
        return newZ

    def printZvalues(self,forward): #print values of Z and lnZ when PRECISION is reached 
        dirname = forward and "zvaluesF.txt" or "zvaluesB.txt"
        fullfilename = os.path.join(path,dirname)
        of = open(fullfilename,"w")
        rangehistos = forward and self.fhistos or self.bhistos
        for histo in rangehistos:
            of.write("%.18f %.18f\n" % (histo.Z, histo.lnZ))
        of.close()    
        
    def printP_est(self,forward): #print final histogram
        dirname = forward and "pestF.txt" or "pestB.txt"
        dirname_fe = forward and "pestF_fe.txt" or "pestB_fe.txt"

        fullfilename = os.path.join(path,dirname)
        pfefilename = os.path.join(path,dirname_fe)

       # of = open(fullfilename,"w")
        #ofe = open(pfefilename,"w")

        tr = forward and self.fhistos or self.bhistos
        p_list = [ 0 ]*HISTOSIZE        
        norm =0.0
        for Q in range(HISTOSIZE):
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
        if do_printfe:
            ofe = open(pfefilename,"w")
            for Q in range(HISTOSIZE):
                p_est = p_list[Q]/norm
                p_list[Q] = p_est
                if p_est > 0:
                    ofe.write("%f %.18f %.18f\n" % (histo.lambdahisto[Q],p_est,log(p_est)))
            ofe.close()
        else:
            of = open(fullfilename,"w")
            for Q in range(HISTOSIZE):
                p_est = p_list[Q]/norm
                p_list[Q] = p_est
                if p_est > 0:
                    of.write("%f %.18f %.18f\n" % (histo.lambdahisto[Q],p_est,log(p_est)))
            of.close() 
        return p_list
    

    def createallhistos(self,forward): #create all the histograms
        
        tr = forward and self.fhistos or self.bhistos
        for histo in tr:
            dirname = forward and "pesthistoF.%02d.txt" % (histo.interface,) or "pesthistoB.%02d.txt" % (histo.interface,)
            fullfilename = os.path.join(path,dirname)
            of = open(fullfilename,"w")
        
            for Q in range(HISTOSIZE):                 
                if histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral>0.0:
                    of.write("%d %.18f %.18f\n" % (histo.lambdahisto[Q],log(histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral),histo.expW[Q]*histo.histodata[Q]*histo.Z/histo.integral))
        of.close()
    
    def printWeights(self,numinter,forward): #print wham weights
        dirname = forward and "forwardweights.txt" or "backwardweights.txt"
        
        fullfilename = os.path.join(path,dirname)
        of = open(fullfilename,"w")
        tr = forward and self.fhistos or self.bhistos
        wmax=-1
        w_inter=[0]*numinter
        c=0
        for histo in tr:
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
        fullfilename = os.path.join(path,dirname)
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
            for i in range(isim+1): 
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

    def printBackwardpathWeights(self,numinter,forward): #print path weights for backward ensemble
        dirname = forward and "fpathweights.txt" or "bpathweights.txt"
        fullfilename = os.path.join(path,dirname)
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

    #def calculatefactors(self,inter,numinter,forward):
        
     #   dirname = forward and "ffactor.txt" or "bfactor.txt"
     #   fullfilename = os.path.join(ospath,dirname)
     #   of = open(fullfilename,"w")
     #   #ihisto = forward and self.fhistos[i] or self.bhistos[i]
     #   s_factor=[0.0]*HISTOSIZE
     #   rangehisto = forward and self.fhistos or self.bhistos
     #   c=0        
     #   for histo in rangehisto:
     #       c += 1
     #       for Q in range(HISTOSIZE):
     #           #if c < 20:
     #               if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
     #                   s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
     #           #else:
     #           #    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
     #           #        s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
     #   for Q in range(HISTOSIZE): 
     #       sfactor=s_factor[Q]
     #       of.write("%f %f\n" % (histo.lambdahisto[Q], sfactor))
     #   return s_factor


    def calculatefactors(self,inter,numinter,forward): #calculate scaling factors s_AB and s_BA (JCP,133,174109-11 (2010))               
        dirname = forward and "ffactor.txt" or "bfactor.txt"
        fullfilename = os.path.join(path,dirname)
        of = open(fullfilename,"w")
        s_factor=[0.0]*HISTOSIZE
        rangehisto = forward and self.fhistos or self.bhistos
        c=0 
        b=0
        for histo in rangehisto:
            c += 1
            for Q in range(HISTOSIZE):
                if forward:
                    if c < len(inter):
                        if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0 and histo.lambdahisto[Q] < inter[c]:
                            s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                    #else:
                     #   if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                      #      s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                else:
                    if b>0:
                        if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0 and histo.lambdahisto[Q] > inter[b-1]: 
                            s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
                    #else:
                     #    if histo.boundaryhistodata[Q] > 0.0 and histo.histodata[Q]>0.0:
                      #      s_factor[Q]=histo.histodata[Q]/histo.boundaryhistodata[Q]
            b +=1
                       
        for Q in range(HISTOSIZE): 
            sfactor=s_factor[Q]
            of.write("%f %f\n" % (histo.lambdahisto[Q], sfactor))
        return s_factor
                
                
    def freeEnergy(self,p0f,p0b,fmab,bmab): #Compute free energy projected on order parameter used in TIS
        
        ca=[0.0]*HISTOSIZE
        cb=[0.0]*HISTOSIZE
        p=-1.0
        wf=0.0
        wb=0.0
        N=0
        filename = os.path.join(path,"cacb.dat")
        cval = open(filename,"w")
        for Q in range(HISTOSIZE):
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
        fullfilename = os.path.join(path,"freeEnergy.dat")
        of = open(fullfilename,"w")
        for Q in range(HISTOSIZE):
            p=avca*p0f[Q]+avcb*p0b[Q]
            if p > 0.0:
                of.write("%f %f\n" % (histo.lambdahisto[Q],-log(p)))
        of.close()

##################################### main ###############################################################################

if __name__=="__main__":
    sim = TSelfConsistent()
    
    for i in range(len(sim.fhr)): 
        sim.fhistos.append(THisto(i,i,forward=True))

    inter=len(sim.fhr)

    for i in range(len(sim.bhr)):        
        sim.bhistos.append(THisto(i,inter,forward=False))
        inter += 1
    
#GDL: Self consistent calculation of Z (forward ensemble AB)
    diff=100.0
    iteration = 0
    while (diff > PRECISION)&(iteration < 100000) :
                
            for i,histo in enumerate(sim.fhistos):
                sim.fhistos[i].newZ = sim.newZperreplica(i,True) 
                sim.fhistos[i].lnnewZ = log(sim.fhistos[i].newZ)
                
            
            diff=0.0
            for i,histo in enumerate(sim.fhistos):
                diff += fabs(sim.fhistos[i].lnnewZ - sim.fhistos[i].lnZ)
                sim.fhistos[i].lnZ = (sim.fhistos[i].lnnewZ - sim.fhistos[0].lnnewZ)
                sim.fhistos[i].Z = exp(sim.fhistos[i].lnZ)
                
            ilist = []
            for histo in sim.fhistos:
                ilist.append(histo.Z)
            iteration += 1 
            print diff, ilist, iteration
#GDL:print final results for AB ensemble
    finter = sim.fhr
    numinter= len(sim.fhr)
    p0f = sim.printP_est(forward = True)
    if do_printfe:
        fmab=sim.calculatefactors(finter,numinter,forward=True)
    else:
        sim.createallhistos(forward= True)
        sim.printZvalues(forward = True)
        numinter= len(sim.fhr)
        sim.printWeights(numinter,forward = True)
        sim.printforwardpathWeights(numinter,forward= True)


#Self consistent calculation of Z (backward ensemble BA)
    diff=100.0
    iteration = 0
    while (diff > PRECISION)&(iteration < 100000) :
                
            for i,histo in enumerate(sim.bhistos):
                sim.bhistos[i].newZ = sim.newZperreplica(i,False) 
                sim.bhistos[i].lnnewZ = log(sim.bhistos[i].newZ)
                        
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
            
#print results for BA ensemble
    binter = sim.bhr
    bnuminter= len(sim.bhr)
    p0b = sim.printP_est(forward = False)
#print matching factors 
    if do_printfe:
        bmab=sim.calculatefactors(binter,bnuminter,forward=False)
    else:

        sim.createallhistos(forward= False)
        sim.printZvalues(forward = False)
        sim.printWeights(bnuminter,forward = False)
        sim.printBackwardpathWeights(bnuminter,forward= False)



#print free energy from free energy histograms if flag '-fe' used
    if do_printfe:    
        sim.freeEnergy(p0f,p0b,fmab,bmab)

    

    
    

    
