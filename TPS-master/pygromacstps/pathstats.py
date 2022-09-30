"""
New module for path statistics
by SRM
"""

import os
import logging
#import numpy as np

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.stats.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False 

#we need the information about interfaces.
#this cannot be taken from anywhere. Rather has to be read in.
#use the tis tools approach 

class pathstats(object):
        def __init__(self):
                self.status = True
                self.intf_type = []
                self.intf_no = []
                self.intf_val = []
                self.intf_accp = []
                self.intf_aa = []
                self.intf_ab = []
                self.intf_ba = []
                self.intf_bb = []
                self.intf_un = []
                self.basedir = os.getcwd()
                self.GenerateIntflist()

        def GenerateIntflist(self):
                for line in open(os.path.join(os.getcwd(),"options","interfaces.txt")):
                        raw = line.split()
                        if len(raw)==0:
                                continue
                        self.intf_type.append(raw[0])
                        self.intf_no.append(int(raw[1]))
                        self.intf_val.append(float(raw[2]))

        #function to read the acceptance
        def PathInfo(self):
                
                #set everything to zero
                self.intf_accp = []
                self.intf_aa = []
                self.intf_ab = []
                self.intf_ba = []
                self.intf_bb = []
                self.intf_un = []
                self.tot = 0

                for i in range(len(self.intf_no)):
                        #cycle through each interface
                        instring = (("%03d")%self.intf_no[i])
                        infile = os.path.join(self.basedir,'gtpslog.'+str(instring)+'.txt')
                        
                        ctot = 0
                        ctrue = 0
                        caa = 0
                        cab = 0
                        cba = 0
                        cbb = 0
                        cun = 0

                        for line in open(infile):
                                raw = line.split()
                                #print raw
                                #draw = raw[4].split()
                                if raw[3]=='True':
                                        ctrue+=1
                                if raw[4]=='AA':
                                        caa+=1
                                elif raw[4]=='AB':
                                        cab+=1
                                elif raw[4]=='BA':
                                        cba+=1
                                elif raw[4]=='BB':
                                        cbb+=1
                                elif raw[4]=='UN':
                                        cun+=1
                                ctot+=1

                        #now find statistics
                        if ctot>0:
                                ctrue = float(ctrue)/float(ctot)
                                caa = float(caa)/float(ctot)
                                cab = float(cab)/float(ctot)
                                cba = float(cba)/float(ctot)
                                cbb = float(cbb)/float(ctot)
                                cun = float(cun)/float(ctot)

                        self.intf_accp.append(ctrue)
                        self.intf_aa.append(caa)
                        self.intf_ab.append(cab)
                        self.intf_ba.append(cba)
                        self.intf_bb.append(cbb)
                        self.intf_un.append(cun)
                        self.tot = ctot


        def PrintReport(self):

                self.PathInfo()
                with open('gtpslog.stats.txt','a') as fout:
                        fout.write(("PATH %d\n")%(self.tot))
                        #logger.info(("Intf  accp  AA    AB    BA    BB    Unkn  "))
                        for i in range(len(self.intf_no)):
                                fout.write(("%s %f %f %f %f %f %f\n")%(self.intf_no[i],self.intf_accp[i],self.intf_aa[i],self.intf_ab[i],self.intf_ba[i],self.intf_bb[i],self.intf_un[i]))







