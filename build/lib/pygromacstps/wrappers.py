'''
Created on Mar 2, 2010

@author: wolf
'''
import os
import subprocess as sub

class gromacswrapper(object):
    def __init__(self):
        self.wolf = "test"
    
    def checkInputFiles(self,options,optionlist,selection,commandpath):
        ok = True
        for key,item in optionlist.iteritems():
            if item[0] in selection:
                if not os.path.exists(os.path.join(commandpath,item[1])):
                    ok = False
                    print "Inputfile " + os.path.join(commandpath,item[1]) + " not found" + str(key)
        return ok
    
    def generateCommand(self,options,optionlist,selection,commandpath,command="grompp"):
        cmd = []
        cmd.append(os.path.join(options.paths["gromacspath"],command+options.runoptions["gromacssuffix"]))
        for key,item in optionlist.iteritems():
            if item[0] in selection:
                cmd.append(item[0])
                cmd.append(os.path.join(commandpath,item[1]))
        cmd.append("-maxwarn")
        cmd.append("10")
        return cmd
    
    def executeCommand(self,cmd,tinput=""):
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input=tinput)
        proc.wait()
        return out,err

    def executeCommandSTDERR(self,cmd,tinput=""):
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE)
        out,err = proc.communicate(input=tinput)
        proc.wait()
        return out,err
    
    def generateMDRunCommand(self,options,commandpath):
        cmd = []
        cmd.append(os.path.join(options.paths["gromacspath"],"mdrun"+options.runoptions["gromacssuffix"]))
        cmd.append("-s")
        cmd.append(os.path.join(commandpath,options.initoptions["tprfile"][1]))
        cmd.append("-o")
        cmd.append(os.path.join(commandpath,"traj.trr"))
        cmd.append("-x")
        cmd.append(os.path.join(commandpath,"traj.xtc"))
        cmd.append("-e")
        cmd.append(os.path.join(commandpath,"ener.edr"))
        cmd.append("-dlb")
        cmd.append("no")
        cmd.append("-nt")
        cmd.append("1")
        
        
        return cmd
    
    def generateTrjconvCommand(self,options,commandpath,outputpath):
        cmd = []
        cmd.append(os.path.join(options.paths["gromacspath"],"trjconv"+options.runoptions["gromacssuffix"]))
        cmd.append("-s")
        cmd.append(os.path.join(commandpath,options.initoptions["tprfile"][1]))
        cmd.append("-f")
        cmd.append(os.path.join(commandpath,"traj.trr"))
        cmd.append("-n")
        cmd.append(os.path.join(commandpath,options.initoptions["ndxfile"][1]))
        cmd.append("-o")
        cmd.append(os.path.join(outputpath,"path.gro"))
        cmd.append("-sep")
        return cmd
        
    def generateGDistCommand(self,workdir,options):
        cmd = []
        cmd.append(os.path.join(options.paths["gromacspath"],"g_dist"+options.runoptions["gromacssuffix"]))
        cmd.append("-f")
        cmd.append(os.path.join(workdir,"traj.trr"))
        cmd.append("-s")
        cmd.append(os.path.join(workdir,options.initoptions["tprfile"][1]))
        cmd.append("-n")
        cmd.append(os.path.join(workdir,options.initoptions["ndxfile"][1]))
        cmd.append("-o")
        cmd.append(os.path.join(workdir,"dist.xvg"))
        return cmd
    
    def generateUmbrellaCommand(self,options,i):
        umbrellapath = os.path.join(os.path.join(options.paths["umbrellapath"],"%02d" % i))
        cmd = []
        cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"grompp"+options.runoptions["gromacssuffix"])))
        cmd.append("-c")
        cmd.append(os.path.join(umbrellapath,"conf%d.gro" %(i)))
        cmd.append("-f")
        cmd.append(os.path.join(umbrellapath,"md.mdp"))
        cmd.append("-p")
        cmd.append(os.path.join(umbrellapath,options.initoptions["topfile"][1]))
        cmd.append("-n")
        cmd.append(os.path.join(umbrellapath,options.initoptions["ndxfile"][1]))
        cmd.append("-o")
        cmd.append(os.path.join(umbrellapath,"umbrella%d.tpr" %(i)))
        cmd.append("-maxwarn")
        cmd.append("10")
        return cmd
    
    def generateGDistCommandCustom(self,options,tprfile,ndxfile,xtcfile,outfile):
        cmd = []
        cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"g_dist"+options.runoptions["gromacssuffix"])))
        cmd.append('-s')
        cmd.append(tprfile)
        cmd.append('-n')
        cmd.append(ndxfile)
        cmd.append('-f')
        cmd.append(xtcfile)
        cmd.append('-o')
        cmd.append(outfile)
        return cmd

    def generateGWhamCommandCustom(self,options):
        cmd = []
        cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"g_wham"+options.runoptions["gromacssuffix"])))
        cmd.append('-it')
        cmd.append(os.path.join(options.paths["whampath"],"tprfiles.dat"))
        cmd.append('-if')
        cmd.append(os.path.join(options.paths["whampath"],"pullfiles.dat"))
        cmd.append('-o')
        cmd.append('-hist')
        return cmd
    
    def executeMDRun(self,options,tinput=""):
        workdir = options.paths["initialpath"]
        
        os.chdir(workdir)
        tprfile=os.path.join(workdir,options.initoptions["tprfile"][1])
        cmd = [os.path.join(options.paths["gromacspath"],"mdrun"+options.runoptions["gromacssuffix"])]
        cmd.append('-s')
        cmd.append(tprfile)
        proc = sub.Popen(cmd, stdout=sub.PIPE, stdin=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input=tinput)
        return out,err
    


class lammpswrapper(object):
    def __init__(self):
        self.wolf = "test"

    
