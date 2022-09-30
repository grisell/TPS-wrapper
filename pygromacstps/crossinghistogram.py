'''
Created on June 16, 2010

@author: wolf
'''

import numpy as np
from math import log

class generichistogram(object):
    """
    A Crossing histogram can be based on a trajectory either forward or backward in time
    histomin is the position of the interface and histomax can be chosen arbitrarily
    """
    def __init__(self,histosize=2000,histomin=-4.0,histomax=4.0,forward=True):
        self.histo = np.zeros(histosize)
        self.fehisto = np.zeros(histosize)
        self.comhisto = np.zeros(histosize)
        self.histomin = histomin
        self.histomax = histomax
        self.histosize = histosize
        self.norm = 0.0
        self.forward = forward
    
    def _getBox(self,value):
        """
        Convert the value to the histogrambox
        """
        return (float(self.histosize)*(float(value) - float(self.histomin)))/(float(self.histomax-self.histomin))
    
    def _getValue(self,hbox):
        """
        Convert the histogram index to the according value
        """
        return (float(hbox)*float(self.histomax-self.histomin) )/float(self.histosize) + float(self.histomin)
    
         
    def addRangeToHisto(self,value):
        """
        Add the value 1 in the range from value to histomax
        """
        hbox = int(round(self._getBox(value)))
        if self.forward:
            if hbox > self.histosize - 1:
                hbox = self.histosize -1
            for i in range(int(hbox)):
                self.histo[i] += 1.0
        else:
            if hbox < 0:
                hbox = 0
            for i in range(self.histosize-int(hbox)):
                self.histo[self.histosize-i-1] += 1.0
        self.norm += 1.0

    def addRangeToHistoFromInterface(self,value,interfacevalue):  #DS+JR for generating histograms from zero
        """
        Add the value 1 in the range from value to histomax
        """
        hbox = int(round(self._getBox(value)))
        minbox = int(round(self._getBox(interfacevalue)))
        if self.forward:
            if hbox > self.histosize - 1:
                hbox = self.histosize -1
            for i in range(int(minbox),int(hbox)+1):
                self.histo[i] += 1.0
        else:
            maxbox = int(round(self._getBox(interfacevalue)))
            if hbox < 0:
                hbox = 0
            #for i in range(int(maxbox),self.histosize-int(hbox)):
            for i in range(int(hbox),int(maxbox)+1):
                #self.histo[self.histosize-i-1] += 1.0
                self.histo[i] += 1.0
        self.norm += 1.0

    
        
    def addPointToHisto(self,value):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self._getBox(value)))
        self.histo[hbox] += 1.0
        self.norm += 1.0

    
    def addweightedPointToHisto(self,value,weight):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self._getBox(value)))
        self.fehisto[hbox] += weight
        #if hbox ==2 :
        #   print value, hbox, weight,  self.fehisto[hbox]
        #self.norm += weight    
        print self.fehisto[hbox]

    def addweightedPointTocommittorHisto(self,value,weight):
        """
        Add a single point to the histogram
        """
        hbox = int(round(self._getBox(value)))
        self.comhisto[hbox] += weight
    
    def outputCrossingHisto(self,filename):
        """
        Print the histogram
        """
        of = open(filename,"w")
        for i in range(self.histosize):
#            of.write("%.18f %.18f\n" % (self._getValue(i),self.histo[i]/self.norm))
            if self.norm != 0:                                                                               #DS to avoid divide by zero
                of.write("%.18f %.18f %.18f\n" % (self._getValue(i),self.histo[i]/self.norm,self.histo[i]))  
            else:
                of.write("%.18f %.18f %.18f\n" % (self._getValue(i),self.histo[i],self.histo[i]))
        of.close()
    
        
    def outputcommitorgHisto(self,filename):
        """
        Print the histogram
        """
        of = open(filename,"w")
        for i in range(self.histosize):
            #print self.fehisto[i]
#            of.write("%.18f %.18f\n" % (self._getValue(i),self.histo[i]/self.norm))
            if self.fehisto[i] != 0.0:                                                                               
                of.write("%.18f %.18f\n" % (self._getValue(i), self.comhisto[i]/self.fehisto[i]))  
                #of.write("%.18f %.18f\n" % (self._getValue(i), self.comhisto[i]))
        of.close()
    
    def outputfreeenergyHisto(self,filename,unit):
        """
        Print the histogram
        """
        max=-1.0
        for i in range(self.histosize):
            if  self.fehisto[i] >  max:
                max = self.fehisto[i]

        of = open(filename,"w")
        for i in range(self.histosize):
#            of.write("%.18f %.18f\n" % (self._getValue(i),self.histo[i]/self.norm))
            if self.fehisto[i]>0.0 :                                                                               
                of.write("%.18f  %.18f  %.18f %.18f\n" % (self._getValue(i),-log(self.fehisto[i]/max)*unit, self.fehisto[i], self.fehisto[i]/max))  
        of.close()
