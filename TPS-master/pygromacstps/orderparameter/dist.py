'''
Created on Oct 18, 2010

@author: wolf
'''

#import pygromacstps.wrappers
import pygromacstps.parser
import pygromacstps.filesystem
import wrappers as wrappers
import os,time

class TDistOrderParameter(object):
    def __init__(self,minA,maxA,minB,maxB,dircolumn):
        self.A = [minA,maxA]
        self.B = [minB,maxB]
        self.dircolumn = dircolumn
        
        self.wrapper    = wrappers.wrapper()
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()
        
    
    def isPathInState(self,path,log):
        cmd = self.wrapper.generateGDistCommand(path.workdir,path.options)
        self.wrapper.executeCommand(cmd, path.options.runoptions["dist1"]+"\n"+path.options.runoptions["dist2"]+"\n" )
        
        if os.path.exists(os.path.join(path.workdir,"dist.xvg")):
            self.distparser.readDist(os.path.join(path.workdir,"dist.xvg"))
            self.filesystem.deleteFileList(os.path.join(path.workdir),"#*.xvg*")
            
            log.log.debug(path.workdir + " " +str(len(self.distparser.data)) + " " +  str(self.distparser.data[0]))
            
            coordinate = float(self.distparser.data[-1][self.dircolumn])
            
            if len(self.distparser.data) > path.options.runoptions["maxlength"]:
                return -1,True
            elif coordinate > self.A[0] and coordinate < self.A[1]:
                return 1,True
            elif coordinate > self.B[0] and coordinate < self.B[1]:
                return 4,True
            else:
                return 0,False
        else:
            return 0,False
            
    
    
    def getFullTrajectory(self,fpath,bpath,fdir,bdir,log):
        traj = []
        os.chdir(bdir)
        cmd = self.wrapper.generateGDistCommand(bdir,bpath.options)
        self.wrapper.executeCommandSTDERR(cmd, bpath.options.runoptions["dist1"]+"\n"+bpath.options.runoptions["dist2"]+"\n" )
        
        filename = os.path.join(bdir,"dist.xvg")
        c = 0
        while not os.path.exists(filename):
            time.sleep(5)
            c += 1
            if c > 20:
                print "ERROR dist.xvg not found"
                break
        self.distparser.readDist(filename)
        #self.filesystem.deleteFile(filename)
        
        self.distparser.data.reverse()
        blength=0
        count = 0
        for line in self.distparser.data:
            traj.append([count,float(line[self.dircolumn]),0])
            blength+=1
            count += 1
        
        os.chdir(fdir)
        cmd = self.wrapper.generateGDistCommand(fdir,fpath.options)
        self.wrapper.executeCommand(cmd, fpath.options.runoptions["dist1"]+"\n"+fpath.options.runoptions["dist2"]+"\n" )
        filename = os.path.join(fdir,"dist.xvg")
        c = 0
        while not os.path.exists(filename):
            time.sleep(5)
            c += 1
            if c > 20:
                print "ERROR dist.xvg not found"
                break
        self.distparser.readDist(filename)
        #self.filesystem.deleteFile(filename)
        
        flength=0
        for line in self.distparser.data[1:]:
            traj.append([count,float(line[self.dircolumn]),1])
            flength+=1
            count += 1

        
        for i in range(len(traj)):
            traj[i][0] = traj[i][0] - blength + fpath.shootingTime
        
        return traj,flength,blength
    
    
            