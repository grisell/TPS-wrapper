'''
Created on May 3, 2010

@author: Wolfgang Lechner
@status: beta 
'''
import pygromacstps.tis as gtis
import sys
import os
import random
import time


def shootingKernel(dirnumber):
    for i in range(sim.kernels.nkernels):
        pythonpath= sim.paths[0].options.paths["pythonpath"]
        mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
        mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
        queuename = sim.paths[0].options.runoptions["queuename"]
        
        filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"tis",str(dirnumber),reverse=False,qsubs=mqsubs,walltime=mwalltime,queuename=queuename)
        if mqsubs == "LL":
            os.system("llsubmit " + filename)
        else:
            os.system("qsub " + filename)
        
def shootingKernelReverse(rdirnumber):
    pythonpath= sim.paths[0].options.paths["pythonpath"]
    mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
    mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
    queuename = sim.paths[0].options.runoptions["queuename"]
    
    filename = sim.qsubsystem.writeKernelQsubFile(basedir, 0,pythonpath,"tis",str(rdirnumber),reverse=True,qsubs=mqsubs,walltime=mwalltime,queuename=queuename)
    if mqsubs == "LL":
        os.system("llsubmit " + filename)
    else:
        os.system("qsub " + filename)
    

def allFinished(dirstring):
    allF = True
    for i in sim.kernels.kernelPaths:
        lapath = sim.paths[i].nfsladir
        newladir = os.path.join(lapath,dirstring)
        filename = os.path.join(newladir,"finished")
        if not(os.path.exists(filename)):
            allF = False
    return allF

def allFinishedReverse(dirstring):
    allF = True
    for i in range(len(sim.reversePaths)):
        lapath = sim.reversePaths[i].nfsladir
        newladir = os.path.join(lapath,dirstring)
        filename = os.path.join(newladir,"finished")
        if not(os.path.exists(filename)):
            allF = False
    return allF

def exchangeKernel(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
    exchange = sim.trialReplicaExchange(dirstring,newdirstring)
    for i in sim.kernels.kernelPaths:
        sim.getFullTrajectory(i, newdirstring)
        sim.paths[i].lastAcceptedFullTrajectory =[]
        for opi in range(len(sim.orderparameters.op)):
            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
    sim.outputAllFullTrajectories(newdirstring)
    return exchange

def bfexchangeKernel(dirnumber):
    dirstring = "%07d" % (dirnumber,)
    newdirstring = "%07d" % (dirnumber+1,)
    sim.readLastAcceptedTrajectories(dirstring)
    for i in sim.kernels.kernelPaths:
        sim._copyLastAcceptedToFull(i)
    
    exchange = sim.trialBackwardForwardExchange(dirstring, newdirstring)
    
    for i in sim.kernels.kernelPaths:
        sim.getFullTrajectory(i, newdirstring)
        for opi in range(len(sim.orderparameters.op)):
            sim.paths[i].lastAcceptedFullTrajectory.append(sim.paths[i].fullTrajectory[opi][:])
    sim.outputAllFullTrajectories(newdirstring)
    return exchange

if __name__ == '__main__':

    basedir = os.getcwd()
    dirnumber = 0
    rdirnumber = 0
    runs = 1
    mode = "tis"
    for arg in sys.argv[1:]:
        argument = arg.split("=")
        if argument[0] == "basedir":
            basedir = argument[1]
        elif argument[0] == "mode":
            mode = argument[1]
        elif argument[0] == "startwith":
            dirnumber = int(float(argument[1]))
        elif argument[0] == "rstartwith":
            rdirnumber = int(float(argument[1]))
        elif argument[0] == "runs":
            runs = int(float(argument[1]))
        else:
            print "Arguments should be one of those: basedir, mode, startwith, rstartwith, runs"
            
    if basedir!="":
        if mode == 'initial': 
            sim = gtis.gromacstis(basedir,"initial",kernel="head")
            for i in range(sim.kernels.nkernels):
                pythonpath= sim.paths[0].options.paths["pythonpath"]
                mqsubs = sim.paths[0].options.runoptions["qsubsystem"]
                mwalltime = sim.paths[0].options.runoptions["qsubwalltime"]
    
                filename = sim.qsubsystem.writeKernelQsubFile(basedir, i,pythonpath,"initial",0,qsubs=mqsubs,walltime=mwalltime)
            
                if mqsubs == "LL":
                    os.system("llsubmit " + " " + filename)
                else:
                    os.system("qsub " + " " + filename)
                    
        if mode == 'tis':
            
            dirstring = "%07d" % (dirnumber,)
            newdirstring = "%07d" % (dirnumber+1,)
            
            rdirstring = "%07d" % (rdirnumber,)
            rnewdirstring = "%07d" % (rdirnumber+1,)
            
            os.chdir(basedir)
            sim = gtis.gromacstis(basedir,"tis",kernel="head")
            
            for i in sim.kernels.kernelPaths:
                sim.getFullTrajectory(i, dirstring)
                sim.paths[i].lastAcceptedFullTrajectory = sim.paths[i].fullTrajectory[:]
            sim.outputAllFullTrajectories(dirstring)
            shootingKernel(dirnumber)
            
#            shootingKernelReverse(rdirnumber)
            
            while True:
                time.sleep(60)
                if allFinished(newdirstring):
                    break
#                if allFinishedReverse(rdirstring):
#                    rdirnumber += 1
#                    rdirstring = "%07d" % (rdirnumber,)
#                    rnewdirstring = "%07d" % (rdirnumber+1,)
#                    shootingKernelReverse(rdirnumber)
#                    
            dirnumber += 1
            
            for run in range(runs):
                os.chdir(basedir)
                dirstring = "%07d" % (dirnumber,)
                newdirstring = "%07d" % (dirnumber+1,)
                rn = random.randint(0,2)
                if rn == 0:
                    shootingKernel(dirnumber)
                    ki = 0
                    while True:
                        time.sleep(60)
                        ki += 1
                        if allFinished(newdirstring):
                            break
#                        if allFinishedReverse(rdirstring):
#                            rdirnumber += 1
#                            rdirstring = "%07d" % (rdirnumber,)
#                            rnewdirstring = "%07d" % (rdirnumber+1,)
#                            if rdirnumber < 20:
#                                shootingKernelReverse(rdirnumber)
                    
                    sim.log.log.info(dirstring + " sh")
                elif rn == 1:
                    exchange = exchangeKernel(dirnumber)
                    sim.log.log.info(dirstring + " ex " + str(exchange))  
                elif rn == 2:
                    exchange = bfexchangeKernel(dirnumber)
                    sim.log.log.info(dirstring + " bf " + str(exchange))  
                dirnumber += 1
        if mode=="test":
            sim = gtis.gromacstis(basedir,"initial",kernel="head")
            
            
            
