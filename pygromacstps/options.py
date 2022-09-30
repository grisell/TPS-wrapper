'''
Created on May 3, 2010

@author: wolf
'''
from UserDict import DictMixin
from math import ceil
import os


class gromacstpsoptions(object):
    def __init__(self,basedir=".",mode="initial"):
        if basedir != ".":
            self.basedir = basedir
        else:
            self.basedir = os.getcwd()
        
        self.paths = odict()
        self.paths["gromacspath"] = "/Users/wolf/localApps/gromacs/bin"
        if "TMPDIR" in os.environ:
            self.paths["scratchpath"] = os.environ['TMPDIR']
        else:
            self.paths["scratchpath"] = "/state/partition1/"
        self.paths["pythonpath"] = "."
        if mode == "initial":
            self.paths["initialpath"] = os.path.join(self.basedir,"initial")
        elif mode == "tps":
            self.paths["initialpath"] = os.path.join(self.basedir,"tps")
        elif mode == "tis":
            self.paths["initialpath"] = os.path.join(self.basedir,"tis")
        else:
            print "mode should be initial or tps"
        self._evaluteOptions("paths.txt",self.paths)
        #self._makepaths()
        
        self.standardfiles = []
        for line in open(os.path.join(self.basedir,"options","standardfiles.txt"),"r"):
            self.standardfiles.append(line.rstrip('\n'))
        
        self.runoptions = odict()
        self.runoptions["updatetime"] = 5
        self.runoptions["dvmax"] = 1.0
        self.runoptions["scdvmax"]=0.01
        self.runoptions["gromacssuffix"] = ""
        self.runoptions["shootfrominterface"] = 1
        self.runoptions["queuename"] = "serial"
        
        
        self.runoptions["qsubpath"] = ""
        self.runoptions["maxlength"] = 1000        
        self.runoptions["devmax"]="1.0"
        self.runoptions["dist1"] = "1"
        self.runoptions["dist2"] = "3"      
        self.runoptions["qsubsystem"] = "PBS"   
        self.runoptions["qsubwalltime"] = "100:00:00"   
        self.runoptions["interfacecoordinate"] = "0"
        self.initoptions = odict()
        self.initoptions["grofile"] = ["-c","conf.gro"]
        self.initoptions["topfile"] = ["-p","topol.top"]
        self.initoptions["mdpfile"] = ["-f","md.mdp"]
        self.initoptions["ndxfile"] = ["-n","index.ndx"]
        self.initoptions["tprfile"] = ["-o","topol.tpr"]
        
        
        self.initoptions["pullgroupname"] = ["","DOPC"]
        self.initoptions["pullgroupnumber"] = ["","1"]
        self.initoptions["refgroupname"] = ["","1VAL"]
        self.initoptions["refgroupnumber"] = ["","3"]
        
        self.mdpoptions = odict()
        self.mdpoptions["title"]                 = "TPS" 
        self.mdpoptions["cpp"]                   = "/usr/bin/cpp" 
        self.mdpoptions["integrator"]            = "md" 
        self.mdpoptions["tinit"]                 = "0.0" 
        self.mdpoptions["dt"]                    = "0.0050" 
        self.mdpoptions["nsteps"]                = "90000000" 
        self.mdpoptions["nstcomm"]               = "1" 
        self.mdpoptions["nstxout"]               = "50000" 
        self.mdpoptions["nstvout"]               = "50000" 
        self.mdpoptions["nstfout"]               = "50000" 
        self.mdpoptions["nstlog"]                = "50000" 
        self.mdpoptions["nstenergy"]             = "50000" 
        self.mdpoptions["nstxtcout"]             = "50000" 
        self.mdpoptions["xtc_precision"]         = "1000" 
        self.mdpoptions["nstlist"]               = "10" 
        self.mdpoptions["ns_type"]               = "grid" 
        self.mdpoptions["pbc"]                   = "xyz" 
        self.mdpoptions["rlist"]                 = "1.4" 
        self.mdpoptions["coulombtype"]           = "pme"
        self.mdpoptions["rcoulomb_switch"]       = "0.0" 
        self.mdpoptions["rcoulomb"]              = "1.2" 
        self.mdpoptions["epsilon_r"]             = "15" 
        self.mdpoptions["vdw_type"]              = "Shift"  
        self.mdpoptions["rvdw_switch"]           = "0.9" 
        self.mdpoptions["rvdw"]                  = "1.2" 
        self.mdpoptions["DispCorr"]              = "No" 
        self.mdpoptions["tcoupl"]                = "nose-hoover" 
        self.mdpoptions["tc-grps"]               = "system" 
        self.mdpoptions["tau_t"]                 = "0.3" 
        self.mdpoptions["ref_t"]                 = "323"  
        self.mdpoptions["Pcoupl"]                = "berendsen"  
        self.mdpoptions["Pcoupltype"]            = "semiisotropic" 
        self.mdpoptions["tau_p"]                 = "0.75" 
        self.mdpoptions["compressibility"]       = "10e-5 10e-5" 
        self.mdpoptions["ref_p"]                 = "1.0 1.0" 
        self.mdpoptions["gen_vel"]               = "no" 
        self.mdpoptions["gen_temp"]              = "323" 
        self.mdpoptions["gen_seed"]              = "665" 
        self.mdpoptions["constraints"]           = "none" 
        self.mdpoptions["constraint_algorithm" ] = "Lincs" 
        self.mdpoptions["unconstrained_start"]   = "no" 
        self.mdpoptions["lincs_order"]           = "4" 
        self.mdpoptions["lincs_warnangle"]       = "30" 
#        self.mdpoptions["pull"]                  = "umbrella" 
#        self.mdpoptions["pull_start"]            = "yes" 
#        self.mdpoptions["pull_init1"]            = "0" 
#        self.mdpoptions["pull_geometry"]         = "direction" 
#        self.mdpoptions["pull_group0"]           = "DOPC" 
#        self.mdpoptions["pull_group1"]           = "1VAL" 
#        self.mdpoptions["pull_vec1"]             = "0 0 1"
#        self.mdpoptions["pull_dim"]              = "N N Y"
#        self.mdpoptions["pull_k1"]               = "3000"
#        self.mdpoptions["pull_rate1"]            = "0.001"
#        self.mdpoptions["pull_nstxout"]          = "100"       
#        self.mdpoptions["pull_nstfout"]          = "100"
        

    def writeQsubFileParallel(self,dir,cores):
        nodes = int(ceil(float(self.gdiststeps)/float(cores)))
        print nodes
        
        for i in range(nodes):
            filename = os.path.join(dir,'psubmit.'+str(i)+'.sub')
            f=open(filename, 'w')
            f.write('#!/bin/bash\n')
            f.write('#$ -N pygus'+str(i) +'\n')
            f.write('#$ -S /bin/bash\n')
            f.write('#$ -o pbs.'+str(i) + '.out\n')
            f.write('#$ -e pbs.'+str(i) + '.err\n')
            f.write('#$ -r n\n')
            f.write('#$ -V\n')
            f.write('#$ -q parallel\n')
            f.write('#$ -pe mpi %d\n' % (cores,))
            f.write('#$ -cwd\n')
            for j in range(cores):
                f.write('cd %02d\n' % (i*cores+j,))
                f.write('mdrun -s umbrella%d.tpr &\n' % (i*cores+j,))
                f.write('cd ..\n')
                if i*cores+j >= self.gdiststeps-1:
                    break
            f.close()
        
         
    def readOptions(self,filename,optiondict):
        if os.path.exists(filename):
            for line in open(filename,"r"):
                raw = line.split("=")
                if len(raw) == 2:
                    for key in optiondict.keys():
                        if raw[0] == key:
                            if raw[1][-1:] == '\n':
                                raw[1] = raw[1][:-1]
                            if optiondict[key].__class__ == "".__class__:
                                optiondict[key] = raw[1]
                            elif optiondict[key].__class__ == [].__class__:
                                optiondict[key][1] = raw[1]
                            elif optiondict[key].__class__ == int(1).__class__:
                                optiondict[key] = int(raw[1])
                            elif optiondict[key].__class__ == float(1).__class__:
                                optiondict[key] = float(raw[1])                                       
        else:
            of = open(filename,"w")
            of.write(" ")
            of.close()
            

    def writeOptions(self,filename,optiondict):
        print filename
        of = open(filename,"w")
        for key,item in optiondict.iteritems():
            if item.__class__ == "".__class__:
                of.write("%s=%s\n"% (key,item))
            elif item.__class__ == [].__class__:
                of.write("%s=%s\n"% (key,item[1]))
                
        of.close()
    
    def _evaluteOptions(self,sfilename,optiondict):
        tmpfilename=os.path.join(self.basedir,sfilename)
        if os.path.exists(tmpfilename):
            self.readOptions(tmpfilename,optiondict)
        else:
            self.writeOptions(tmpfilename,optiondict)
            
    def _makepaths(self):
        for key,item in self.paths.iteritems():
            if not (key in ["gromacspath"]):
                if not(os.path.exists(item)):
                    os.mkdir(item)
                     
    def writeMdpFile(self,directory,sfilename):
        filename = os.path.join(directory,sfilename)
        of = open(filename,"w")
        for key,item in self.mdpoptions.iteritems():
            of.write('%-15s = %s\n' % (key, item))
        of.close()
        
    def writeQsubFile(self,directory,sfilename,additional):
        filename = os.path.join(directory,sfilename)
        of = open(filename,"w")
        for line in self.qsuboptions:
            of.write(line + "\n")
        for line in additional:
            of.write(line + "\n")
        of.close()
    
    def mcopyFile(self,destdir,destfile,targetdir,targetfile):
        df = os.path.join(destdir,destfile)
        tf = os.path.join(targetdir,targetfile)
        os.system("cp " + df + " " + tf)
        
       






class lammpstpsoptions(object):
    def __init__(self,basedir=".",mode="initial"):
	if basedir != ".":
            self.basedir = basedir
        else:
            self.basedir = os.getcwd()
        
        self.paths = odict()
        self.paths["binpath"] = "/Users/rogal/bin"
	self.paths["tispath"] = ""
	self.paths["wrapperpath"] = "/Users/rogal/work"
	self.paths["openmpipath"] = " "
        if "TMPDIR" in os.environ:
            self.paths["scratchpath"] = os.environ['TMPDIR']
        else:
            self.paths["scratchpath"] = "/state/partition1/"
        self.paths["pythonpath"] = ""
        if mode == "initial":
            self.paths["initialpath"] = os.path.join(self.basedir,"initial")
        elif mode == "tps":
            self.paths["initialpath"] = os.path.join(self.basedir,"tps")
        elif mode == "tis":
            self.paths["initialpath"] = os.path.join(self.basedir,"tis")
        else:
            print "mode should be initial or tps"
        self._evaluteOptions("paths.txt",self.paths)   #JR resets scratchpath, binpath etc...
        #self._makepaths()
        
        self.standardfiles = []
        for line in open(os.path.join(self.basedir,"options","standardfiles.txt"),"r"):
            self.standardfiles.append(line.rstrip('\n'))
        
        self.runoptions = odict()
	self.runoptions["lammpsbinary"] = "lmp_mac"
        self.runoptions["updatetime"] = 5
	self.runoptions["dvmax"] = 1.0
        self.runoptions["scdvmax"] = 0.01
#        self.runoptions["gromacssuffix"] = ""
        self.runoptions["shootfrominterface"] = 1
        self.runoptions["queuename"] = "serial"        
        self.runoptions["qsubpath"] = ""
        self.runoptions["maxlength"] = 1000        
        self.runoptions["dist1"] = "1"
        self.runoptions["dist2"] = "3"      
        self.runoptions["qsubsystem"] = "PBS"   
        self.runoptions["qsubwalltime"] = "100:00:00"   
        self.runoptions["interfacecoordinate"] = "0"
	#options to select an certain md driver
	self.runoptions["mddriver"] = "lammps"

        self.initoptions = odict()
	self.initoptions["conffile"] = ["","conf.dump"]
	self.initoptions["mdinfile"] = ["","md.in"]
	self.initoptions["pathfile"] = ["","traj.dat"]
#        self.initoptions["grofile"] = ["-c","conf.gro"]
#        self.initoptions["topfile"] = ["-p","topol.top"]
#        self.initoptions["mdpfile"] = ["-f","md.mdp"]
#        self.initoptions["ndxfile"] = ["-n","index.ndx"]
#        self.initoptions["tprfile"] = ["-o","topol.tpr"]
        
        
#        self.initoptions["pullgroupname"] = ["","DOPC"]
#        self.initoptions["pullgroupnumber"] = ["","1"]
#        self.initoptions["refgroupname"] = ["","1VAL"]
#        self.initoptions["refgroupnumber"] = ["","3"]
        
        self.mdpoptions = odict()
        self.mdpoptions["echo"]                  = "both" 
        self.mdpoptions["units"]                 = "metal" 
        self.mdpoptions["atom_style"]            = "atomic" 
        self.mdpoptions["boundary"]              = "p p p" 
        self.mdpoptions["lattice"]               = "fcc 1.0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1"
        self.mdpoptions["region"]                = "box block 0 1 0 1 0 1"
        self.mdpoptions["create_box"]            = "1 box"
        self.mdpoptions["read_dump"]             = "conf.dump 0 x y z vx vy vz scaled no box yes add yes"
        self.mdpoptions["pair_style"]            = "eam/alloy"
        self.mdpoptions["pair_coeff"]            = "* * Cu01.eam.alloy Cu"
        self.mdpoptions["neighbor"]              = "1.0 bin"
        self.mdpoptions["neigh_modify"]          = "every 1 delay 1 check yes"
        self.mdpoptions["fix"]                   = "1 all npt temp 100 100 0.05 iso 0.0 0.0 0.1"
        self.mdpoptions["timestep"]              = "0.002"
        self.mdpoptions["compute"]               = "1 all temp"
        self.mdpoptions["dump"]                  = "1 all custom 5000 traj.dat id type mass x y z vx vy vz"
        self.mdpoptions["thermo_style"]          = "custom step pe c_1 press vol etotal"
	self.mdpoptions["thermo"]                = "10"
	self.mdpoptions["run"]                   = "100000"



 
	
	
    def readOptions(self,filename,optiondict):
        if os.path.exists(filename):
            for line in open(filename,"r"):
                raw = line.split("=")
                if len(raw) == 2:
                    for key in optiondict.keys():
                        if raw[0] == key:
                            if raw[1][-1:] == '\n':
                                raw[1] = raw[1][:-1]
                            if optiondict[key].__class__ == "".__class__:
                                optiondict[key] = raw[1]
                            elif optiondict[key].__class__ == [].__class__:
                                optiondict[key][1] = raw[1]
                            elif optiondict[key].__class__ == int(1).__class__:
                                optiondict[key] = int(raw[1])
                            elif optiondict[key].__class__ == float(1).__class__:
                                optiondict[key] = float(raw[1])                   
        else:
            of = open(filename,"w")
            of.write(" ")
            of.close()
            

    def writeOptions(self,filename,optiondict):
        print filename
        of = open(filename,"w")
        for key,item in optiondict.iteritems():
            if item.__class__ == "".__class__:
                of.write("%s=%s\n"% (key,item))
            elif item.__class__ == [].__class__:
                of.write("%s=%s\n"% (key,item[1]))
                
        of.close()
    
    def _evaluteOptions(self,sfilename,optiondict):
        tmpfilename=os.path.join(self.basedir,sfilename)
        if os.path.exists(tmpfilename):
            self.readOptions(tmpfilename,optiondict)
        else:
            self.writeOptions(tmpfilename,optiondict)

        




	


        
        
class odict(DictMixin):
    
    def __init__(self):
        self._keys = []
        self._data = {}
        
        
    def __setitem__(self, key, value):
        if key not in self._data:
            self._keys.append(key)
        self._data[key] = value
        
        
    def __getitem__(self, key):
        return self._data[key]
    
    
    def __delitem__(self, key):
        del self._data[key]
        self._keys.remove(key)
        
        
    def keys(self):
        return list(self._keys)
    
    
    def copy(self):
        copyDict = odict()
        copyDict._data = self._data.copy()
        copyDict._keys = self._keys[:]
        return copyDict
