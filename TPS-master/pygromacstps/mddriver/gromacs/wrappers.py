'''
Created on Mar 2, 2010

@author: wolf
'''
import os
import subprocess as sub

class wrapper(object):
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
    
    def generateCommand(self,options,optionlist,commandpath, command="grompp"):
        cmd = []
        #cmd.append(os.path.join(options.paths["gromacspath"],command+options.runoptions["gromacssuffix"]))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append(command)
        cmd.append("-c")
        cmd.append(os.path.join(commandpath,options.initoptions["grofile"]))
        cmd.append("-o")
        cmd.append("-p")
        cmd.append(os.path.join(commandpath,options.initoptions["topfile"]))
        cmd.append("-f")
        cmd.append(os.path.join(commandpath,options.initoptions["mdpfile"]))
        cmd.append("-n")
        cmd.append(os.path.join(commandpath,options.initoptions["ndxfile"]))
        cmd.append("-maxwarn")
        cmd.append("10")
        return cmd
    
    #def generateCommand(self,options,optionlist,selection,commandpath,command="grompp"):
     #   cmd = []
        #cmd.append(os.path.join(options.paths["gromacspath"],command+options.runoptions["gromacssuffix"]))
      #  cmd.append(os.path.join(options.paths["gromacspath"],options.runoptions["mdbinary"]+" "+command))
       # for key,item in optionlist.iteritems():
        #    if item[0] in selection:
       #         cmd.append(item[0])
               # cmd.append(os.path.join(commandpath,item[1]))
       # cmd.append("-maxwarn")
       # cmd.append("10")
       # return cmd
    
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
        #cmd.append(os.path.join(options.paths["gromacspath"],"mdrun"+options.runoptions["gromacssuffix"]))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append("mdrun")
        cmd.append("-s")
        cmd.append(os.path.join(commandpath,options.initoptions["tprfile"]))
        cmd.append("-o")
        cmd.append(os.path.join(commandpath,"traj.trr"))
        cmd.append("-x")
        cmd.append(os.path.join(commandpath,"traj.xtc"))
        cmd.append("-e")
        cmd.append(os.path.join(commandpath,"ener.edr"))
        #cmd.append("-dlb")
        #cmd.append("no")
        #cmd.append("-nt")
        #cmd.append("1")
                
        return cmd
    
    """
    def TrajToConf(self,options,commandpath):
        #generateTrajconvcommand modified
        cmd = []
        #cmd.append(os.path.join(options.paths["gromacspath"],"trjconv"+options.runoptions["gromacssuffix"]))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append("trjconv")
        cmd.append("-s")
        cmd.append(os.path.join(commandpath,options.initoptions["tprfile"]))
        cmd.append("-f")
        cmd.append(os.path.join(commandpath,"traj.trr"))
        cmd.append("-n")
        cmd.append(os.path.join(commandpath,options.initoptions["ndxfile"]))
        cmd.append("-o")
        cmd.append(os.path.join(commandpath,"path.gro"))
        cmd.append("-sep")
        self.executeCommand(cmd, tinput="0")
    """
    
    def TrajToConf(self,options,tmpdirpath):

        os.chdir(tmpdirpath)
        filename = os.path.join(tmpdirpath,'traj.gro')
        data = []
        if os.path.exists(filename):
                infile = open(filename,"r")
                for line in infile:
                        data.append(line)
        if len(data) > 0:
                nlines = len(data)
                natoms = int(data[1])
                nblock = natoms+3
                nslices = nlines/nblock

        for i in range(nslices):
                fout = "path%d.gro" % i
                filename = os.path.join(tmpdirpath,fout)
                wf = open(filename,"w")

                start = i*nblock
                wf.write(data[start])
                #wf.write("0\n")
                for j in range(1,nblock):
                        wf.write(data[start+j])          
        
    def generateGDistCommand(self,workdir,options):
        cmd = []
        #cmd.append(os.path.join(options.paths["gromacspath"],"g_dist"+options.runoptions["gromacssuffix"]))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append("g_dist")
        cmd.append("-f")
        cmd.append(os.path.join(workdir,"traj.trr"))
        cmd.append("-s")
        cmd.append(os.path.join(workdir,options.initoptions["tprfile"]))
        cmd.append("-n")
        cmd.append(os.path.join(workdir,options.initoptions["ndxfile"]))
        cmd.append("-o")
        cmd.append(os.path.join(workdir,"dist.xvg"))
        return cmd
    
    def generateUmbrellaCommand(self,options,i):
        umbrellapath = os.path.join(os.path.join(options.paths["umbrellapath"],"%02d" % i))
        cmd = []
        #cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"grompp"+options.runoptions["gromacssuffix"])))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append("grompp")
        cmd.append("-c")
        cmd.append(os.path.join(umbrellapath,"conf%d.gro" %(i)))
        cmd.append("-f")
        cmd.append(os.path.join(umbrellapath,"md.mdp"))
        cmd.append("-p")
        cmd.append(os.path.join(umbrellapath,options.initoptions["topfile"]))
        cmd.append("-n")
        cmd.append(os.path.join(umbrellapath,options.initoptions["ndxfile"]))
        cmd.append("-o")
        cmd.append(os.path.join(umbrellapath,"umbrella%d.tpr" %(i)))
        cmd.append("-maxwarn")
        cmd.append("10")
        return cmd
    
    def generateGDistCommandCustom(self,options,tprfile,ndxfile,xtcfile,outfile):
        cmd = []
        #cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"g_dist"+options.runoptions["gromacssuffix"])))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append('g_dist')
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
        #cmd.append(os.path.join(os.path.join(options.paths["gromacspath"],"g_wham"+options.runoptions["gromacssuffix"])))
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append('g_wham')
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
        tprfile=os.path.join(workdir,options.initoptions["tprfile"])
        #cmd = [os.path.join(options.paths["gromacspath"],"mdrun"+options.runoptions["gromacssuffix"])]
        cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        cmd.append('mdrun')
        cmd.append('-s')
        cmd.append(tprfile)
        proc = sub.Popen(cmd, stdout=sub.PIPE, stdin=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input=tinput)
        return out,err


    def initializeMD(self,path,basedir):
        workdir = path.workdir
        os.chdir(workdir)
        ok = self.checkInputFiles(path.options, path.options.initoptions,                         \
                                     ["-c","-p","-n"],workdir)
        #self.log.log.debug(workdir + str(ok) + " grompp executed")
       # path.options.writeMdpFile(workdir,"md.mdp")
        
        #cmd = self.wrapper.generateCommand(path.options, path.options.initoptions,                        \
           #                           ["-c","-o","-p","-f","-n"],                                  \
           #                           workdir,"grompp" )
        cmd = self.generateCommand(path.options, path.options.initoptions,                        \
                                      workdir,"grompp" )
        self.executeCommand(cmd)
        os.chdir(basedir)
    
    
    def cutTrajFile(self,fdir,bdir,fl,bl):
        
        filename = os.path.join(fdir,'traj.gro')

        with open(filename) as myfile:
                data = [next(myfile) for x in xrange(4)]
        natoms = int(data[1])
        nblock = natoms + 3

        N = fl*nblock
        with open(filename) as myfile:
                head = [next(myfile) for x in xrange(N)]
        os.remove(filename)
        with open(filename,'w') as fout:
                for line in head:
                        fout.write(line)

        filename = os.path.join(bdir,'traj.gro')
        N = bl*nblock
        with open(filename) as myfile:
                head = [next(myfile) for x in xrange(N)]
        os.remove(filename)
        with open(filename,'w') as fout:
                for line in head:
                        fout.write(line)    
    

    
    
