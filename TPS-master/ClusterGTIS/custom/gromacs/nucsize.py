'''
Created on Apr 22, 2016

@author: gdl
'''

import wrappers as wrappers
import pygromacstps.parser
import pygromacstps.filesystem
import logging

import os

#set up logger
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False

# The path where the program water_around is located
wapath = "/home/cvreede/GTIS/pypwater/watertest"


class TCustomOrderParameter(object):
    def __init__(self,minA,maxA,minB,maxB):
        # The stable states of water molecules around the residue
        self.A = [minA,maxA]
        self.B = [minB,maxB]
        # dircolumn is the column from the output where the interesting coordinated is
        #self.dircolumn = 7
        
        self.wrapper    = wrappers.wrapper()
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()


    def _makeTrajconv(self,path,workdir):
        # Here we generate the following command with absolute paths
        # trjconv -f traj.trr -s conf.gro -o wa.xtc -n wa.ndx -pbc res -ur compact -center >& waxtc.log
        
        cmd = []
        cmd.append(os.path.join(path.options.paths["mdbinpath"],path.options.runoptions["mdbinary"]))
        cmd.append("trjconv")
        cmd.append("-f")
        cmd.append(os.path.join(workdir,"traj.trr"))
        cmd.append("-o")
        cmd.append(os.path.join(workdir,"traj.gro"))
        cmd.append("-n")
        cmd.append(os.path.join(workdir,"index.ndx"))
        return cmd
    
    def getTraj(self,path,workdir):
        # Get the trajectory of waters around
        # First we need the trajconv to center the pyp
        cmd = self._makeTrajconv(path,workdir)
        #self.wrapper.executeCommand(cmd, "1\n0\n")
        output,outerr = self.wrapper.executeCommand(cmd)
        return output,outerr

    def _makecalcOP(self,path,workdir):
        
        cmd = []
        cmd.append(os.path.join(path.options.paths["binpath"],"nucsizegro"))
        cmd.append(os.path.join(workdir,"tmp_conf"))
        return cmd
        
    def isPathInState(self,path,log):
        # Does the path end in one of the two stable states?
        # get trajectory and look at the last coordinate
        # # get the last slice of current trajectory
        # #filename = os.path.join(path.workdir,path.options.initoptions["traj.gro"][1])
        
        filename = os.path.join(path.workdir,"traj.gro")

        if os.path.exists(filename):
            os.system("rm traj.gro")
            #logger.info("traj.gro deleted.")

        output,outerr = self.getTraj(path,path.workdir)
        #logger.info("traj.gro generated using trjconv.")

        #check one more if generation worked
        if not os.path.exists(filename):
            logger.error("traj.gro generation failed. The trjconv output is provided below")
            logger.error(output)
            logger.error(outerr)
            raise SystemExit()


        data = []
        infile = open(filename,"r")
        for line in infile:
            data.append(line)
        
        if len(data) > 2:
            natoms = int(data[1])
        else:
            return 0,False
        
        #check if file is completely written
        if len(data)%(natoms+3) != 0:
            return 0,False
        if len(data) > natoms:
            #write last slice to tmp file
            tmpfile = os.path.join(path.workdir,"tmp_conf")
            outfile = open(tmpfile,"w")
            for i in range(-(natoms+3),0):
                outfile.write(data[i])
            
            outfile.flush()
            outfile.close()
            #calculate order parameter
            #using external program...
            cmd = self._makecalcOP(path,path.workdir)
            output,outerr = self.wrapper.executeCommand(cmd)
            q = float(output)
                        #print q
            os.remove(tmpfile)
        else:
            return 0,False


        if q > self.A[0] and q < self.A[1]:
            return 1,True
        elif q > self.B[0] and q < self.B[1]:
            return 4,True
        else:
            return 0,False

    
    def getQTrajectory(self,path,wdir,data):
        
        nlines = len(data)
        natoms = int(data[1])
        nblock = natoms+3
        nslices = nlines/nblock

        if nlines%nblock != 0:
            logger.error("Error in getQTrajectory")
            logger.error("nlines %% nblock != 0")
            logger.error("Exiting program...")
            raise SystemExit()

        qtraj=[]
        for j in range(nslices):
            start = j*nblock
            end = (j+1)*nblock
            tmpfile = os.path.join(wdir,"tmp_conf")
            outfile = open(tmpfile,"w")
            for i in range(start,end):
                outfile.write(data[i])
            outfile.flush()
            outfile.close()
            cmd = self._makecalcOP(path,wdir)
            output,outerr = self.wrapper.executeCommand(cmd)
            q = float(output)
            qtraj.append(q)

        os.remove(tmpfile)
        return qtraj

    def getFullTrajectory(self,fpath,bpath,fdir,bdir,log):

        dirlist = [bdir,fdir]
        pathlist = [bpath,fpath]
        traj = []
        qtrajback = []
        qtrajfor = []
        flength=0
        blength=0

        #first the backward part
        
        os.chdir(bdir)
        filename = os.path.join(bdir,"traj.gro")

        if os.path.exists(filename):
            os.system("rm traj.gro")
            #logger.info("traj.gro deleted.")

        output,outerr = self.getTraj(bpath,bdir)
        #logger.info("traj.gro generated using trjconv.")

        #check one more if generation worked
        if not os.path.exists(filename):
            logger.error("traj.gro generation failed. The trjconv output is provided below")
            logger.error(output)
            logger.error(outerr)
            raise SystemExit()

        data = []
        infile = open(filename,"r")
        for line in infile:
            data.append(line)


        if len(data) > 2:
            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.info("backward part of trajectory not fully written in workdir, removing unfinished slice.")
                newlines = nblock*nslices
                tmpfile = os.path.join(bdir,'tmp-traj.dat')
                outfile = open(tmpfile,"w")
                for i in range(newlines):
                    outfile.write(data[i])
                outfile.flush()
                outfile.close()
                infile.close()
                shutil.move(tmpfile,filename)

            #re-read file with proper slices
            data = []
            infile = open(filename,"r")
            for line in infile:
                data.append(line)

            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.error("Error in getQTrajectory")
                logger.error("nlines %% nblock != 0")
                logger.error("Exiting program...")
                raise SystemExit()
        
        if len(data) > 0:
            qtrajback = self.getQTrajectory(bpath,bdir,data)
            qtrajback.reverse()
            blength=0
            count = 0
#       print 'traback, data',len(qtrajback),len(data)
#       sys.stdout.flush()
            for i in range(len(qtrajback)):
                traj.append([count,qtrajback[i],0])
                blength+=1
                count += 1
        
    #now the forward part
        os.chdir(fdir)
        filename = os.path.join(fdir,"traj.gro")

        if os.path.exists(filename):
            os.system("rm traj.gro")
            #logger.info("traj.gro deleted.")

        output,outerr = self.getTraj(fpath,fdir)
        #logger.info("traj.gro generated using trjconv.")

        #check one more if generation worked
        if not os.path.exists(filename):
            logger.error("traj.gro generation failed. The trjconv output is provided below")
            logger.error(output)
            logger.error(outerr)
            raise SystemExit()

        data = []
        infile = open(filename,"r")
        for line in infile:
            data.append(line)

        
        if len(data) > 2:
            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.info("backward part of trajectory not fully written in workdir, removing unfinished slice.")
                newlines = nblock*nslices
                tmpfile = os.path.join(fdir,'tmp-traj.dat')
                outfile = open(tmpfile,"w")
                for i in range(newlines):
                    outfile.write(data[i])
                outfile.flush()
                outfile.close()
                infile.close()
                shutil.move(tmpfile,filename)

            #re-read file with proper slices
            data = []
            infile = open(filename,"r")
            for line in infile:
                data.append(line)

            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.error("Error in getQTrajectory")
                logger.error("nlines %% nblock != 0")
                logger.error("Exiting program...")
                raise SystemExit()

        if len(data) > 0:

            qtrajfor = self.getQTrajectory(fpath,fdir,data)
            flength=0
            for i in range(len(qtrajfor)):
                traj.append([count,qtrajfor[i],1])
                flength+=1
                count += 1

            if len(traj) > 0:
                for i in range(len(traj)):
                    traj[i][0] = traj[i][0] - blength + fpath.shootingTime

        else:
            traj.append([0,0.0,-1])
            logger.error("Some generic values are being appended in traj. Now will be the right time to check the code")


        return traj,flength,blength


    def getFullTrajectoryReverse(self,fpath,fdir,log):
        traj = []
        qtrajfor = []
        flength=0        
    
        #now the forward part
        os.chdir(fdir)
        filename = os.path.join(fdir,"traj.gro")
        data = []
        if os.path.exists(filename):
                infile = open(filename,"r")
                for line in infile:
                        data.append(line)  

        if len(data) > 2:
            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.info("reverse trajectory not fully written in workdir, removing unfinished slice...")
                newlines = nblock*nslices
                tmpfile = os.path.join(bdir,'tmp-traj.dat')
                outfile = open(tmpfile,"w")
                for i in range(newlines):
                    outfile.write(data[i])
                outfile.flush()
                outfile.close()
                infile.close()
                shutil.move(tmpfile,filename)

            #re-read file with proper slices
            data = []
            infile = open(filename,"r")
            for line in infile:
                data.append(line)

            nlines = len(data)
            natoms = int(data[1])
            nblock = natoms+3
            nslices = nlines/nblock

            if nlines%nblock != 0:
                logger.error("Error in getQTrajectory")
                logger.error("nlines %% nblock != 0")
                logger.error("Exiting program...")
                raise SystemExit()
        
        if len(data) > 0:
            qtrajfor = self.getQTrajectory(fpath,fdir,data)
            #qtrajback.reverse()
            flength=0
            count = 0
#       print 'traback, data',len(qtrajback),len(data)
#       sys.stdout.flush()
            for i in range(len(qtrajfor)):
                traj.append([count,qtrajfor[i],1])
                flength+=1
                count += 1         

        if len(traj) < 1:
                traj.append([0,0.0,-1])
        
        return traj,flength

    def QofLast(self,path):
        # Q value of the last slice
        # get trajectory and look at the last coordinate
        # get the last slice of current trajectory
        filename = os.path.join(path.workdir,"traj.gro")
        data = []
        if os.path.exists(filename):
                infile = open(filename,"r")
                for line in infile:
                        data.append(line)

                if len(data) > 2:
                        natoms = int(data[1])
                else:
                        return 0,0

                nlines = len(data)
            
                nblock = natoms+3
                nslices = nlines/nblock

                if nlines%nblock != 0:
                        return 0,0

                if len(data) > natoms:
                        #write last slice to tmp file
                        tmpfile = os.path.join(path.workdir,"tmp_conf")
                        outfile = open(tmpfile,"w")
                        for i in range(-(natoms+3),0):
                                outfile.write(data[i])
                        outfile.flush()
                        outfile.close()
                        #calculate order parameter
                        #using external program...
                        cmd = self._makecalcOP(path,path.workdir)
                        output,outerr = self.wrapper.executeCommand(cmd)
                        q = float(output)
                        os.remove(tmpfile)

                        return nslices,q
                else:
                        return 0,0

        else:
                return 0,0