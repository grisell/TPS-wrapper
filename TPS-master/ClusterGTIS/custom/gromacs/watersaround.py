'''
Created on Oct 26, 2010

@author: wolf
'''

import pygromacstps.wrappers
import pygromacstps.parser
import pygromacstps.filesystem

import os

# The path where the program water_around is located
wapath = "/home/cvreede/GTIS/pypwater/watertest"


class TCustomOrderParameter(object):
    def __init__(self):
        # The stable states of water molecules around the residue
        self.A = [0,4]
        self.B = [10,20]
        # dircolumn is the column from the output where the interesting coordinated is
        self.dircolumn = 7
        
        self.wrapper    = pygromacstps.wrappers.gromacswrapper()
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()
        
    
    def _makeTrajconv(self,path,workdir):
        # Here we generate the following command with absolute paths
        # trjconv -f traj.trr -s conf.gro -o wa.xtc -n wa.ndx -pbc res -ur compact -center >& waxtc.log
        
        cmd = []
        cmd.append(os.path.join(path.options.paths["mdbinpath"],"trjconv"))
        cmd.append("-f")
        cmd.append(os.path.join(workdir,"traj.trr"))
        cmd.append("-o")
        cmd.append(os.path.join(workdir,"wa.xtc"))
        cmd.append("-s")
        cmd.append(os.path.join(workdir,path.options.initoptions["grofile"][1]))
        cmd.append("-n")
        cmd.append(os.path.join(workdir,"wa.ndx"))
        cmd.append("-pbc")
        cmd.append("res")
        cmd.append("-ur")
        cmd.append("compact")
        cmd.append("-center")
        return cmd
        
    def _makeWaterAround(self,path,workdir):
        # This generates the command water_around
        cmd = []
        cmd.append(os.path.join(wapath,"water_around"))
        cmd.append("46GLU")
        cmd.append(os.path.join(workdir,path.options.initoptions["grofile"][1]))
        cmd.append(os.path.join(workdir,"wa.xtc"))
        return cmd
        
    
    def getTraj(self,path,workdir):
        # Get the trajectory of waters around
        # First we need the trajconv to center the pyp
        cmd = self._makeTrajconv(path,workdir)
        self.wrapper.executeCommand(cmd, "1\n0\n")
        # Then execute the water_around
        cmd = self._makeWaterAround(path,workdir)
        # Pipe the output
        output,outerr = self.wrapper.executeCommand(cmd)
        # split at the newlines
        lines = output.split("\n")
        traj = []
        for line in lines:
            raw = line.split()
            if len(raw) >0:
                # and add to the traj
                traj.append([float(raw[0]),float(raw[self.dircolumn])])
        return traj
        
    def isPathInState(self,path,log):
        # Does the path end in one of the two stable states?
        # get trajectory and look at the last coordinate
        traj = self.getTraj(path,path.workdir)
        if len(traj) > 0:
            coordinate = traj[-1][1]
            if coordinate > self.A[0] and coordinate < self.A[1]:
                return 1,True
            elif coordinate > self.B[0] and coordinate < self.B[1]:
                return 4,True
            else:
                return 0,False
        else:
            return 0, False
    
    
    def getFullTrajectory(self,fpath,bpath,fdir,bdir,log):
        fullTraj=[]
        # getpartial Trajectories
        bTraj = self.getTraj(bpath,bdir)
        fTraj = self.getTraj(fpath,fdir)
        flength = len(fTraj)
        blength = len(bTraj)
        
        #reverse the backwardpath
        bTraj.reverse()
        
        #create full trajectory from the backward and forward path
        for i in range(blength):
            fullTraj.append([i-blength,bTraj[i][1]])
        
        for i in range(1,flength):
            fullTraj.append([i-1,fTraj[i][1]])
        
        # subtract the shootingtime
        for i in range(len(fullTraj)):
            fullTraj[i][0] = fullTraj[i][0] + fpath.shootingTime
        return fullTraj,flength,blength
    
    
            
