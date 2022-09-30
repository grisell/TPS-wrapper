'''
Created on May 3, 2010

@author: wolf
'''
from UserDict import DictMixin
from math import ceil
import os
import logging

#SRM:setting up logger
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False


class tpsoptions(object):
    def __init__(self,basedir=".",mode="initial"):
	if basedir != ".":
	    self.basedir = basedir
	else:
	    self.basedir = os.getcwd()

############check the directory setup first##########################
#
#        if not os.path.exists(os.path.join(os.getcwd(),mode)):
#            if (mode=="initial"):
#                logger.error("Folder initial not found!")
#                raise SystemExit()
#            elif ((mode=="tis") and os.path.exists(os.path.join(os.getcwd(),"initial"))):
#                logger.error("Folder tis not found. But there is an initial folder. Perhaps you forgot to rename it after the initial run?")
#                raise SystemExit()
#            elif (mode=="tis"):
#                logger.error("Folder tis not found!")
#                raise SystemExit()

##########################options paths.txt##########################        
	self.paths = odict()
	self.paths["binpath"] = ""
	self.paths["tispath"] = ""
	self.paths["wrapperpath"] = ""
	self.paths["mdbinpath"] = "/home/users/diazlgjj/programs/gromacs/bin"
	self.paths["openmpipath"] = " "
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
	
	#now reset the paths by reading them from paths.txt
	self.readOptions("paths.txt",self.paths) 

	 #error checking for paths
	if not os.path.exists(self.paths["binpath"]):
	    logger.error("binpath does not exist or is incorrect")
	    raise SystemExit()

	#if scratchpath does not exist,create one
	if not os.path.exists(self.paths["scratchpath"]):
	    if not os.path.exists(os.path.join(os.getcwd(),"scratch")):
		os.mkdir("scratch")
	    self.paths["scratchpath"] = os.path.join(os.getcwd(),"scratch")
	    #logger.info("scratch path incorrect. Scratch folder is set to workingfolder/scratch")

	#wrapperpath
	if not os.path.exists(self.paths["wrapperpath"]):
	    logger.error("wrapperpath does not exist or is incorrect")
	    raise SystemExit()

	#check for tispath. If it doesnt exist, create one
	if not os.path.exists(self.paths["tispath"]):
	    self.paths["tispath"] = self.paths["wrapperpath"] + "/ClusterGTIS"
	    

	if not os.path.exists(self.paths["mdbinpath"]):
	    logger.error("mdbinpath does not exist or is incorrect")
	    raise SystemExit()
	
###########################standard files#######################################
	#self.standardfiles = []
	#for line in open(os.path.join(self.basedir,"options","standardfiles.txt"),"r"):
	#    self.standardfiles.append(line.rstrip('\n'))
	#    SRM edit
	self.standardfiles = odict()
	self.standardfiles["conffile"] = "conf.gro"
	self.standardfiles["mdinfile"] = "md.mdp"
	self.standardfiles["topolfile"] = "topol.top"
	self.standardfiles["indexfile"] = "index.ndx"
	self.standardfiles["customfile1"] = ""
	self.standardfiles["customfile2"] = ""
	self.standardfiles["customfile3"] = ""
	self.standardfiles["customfile4"] = ""
	self.standardfiles["customfile5"] = ""
	self.standardfiles["customfile6"] = ""
	self.standardfiles["customfile7"] = ""
	self.standardfiles["customfile8"] = ""
	self.standardfiles["customfile9"] = ""
	self.standardfiles["customfile10"] = ""
	self.standardfiles["customfile11"] = ""
	self.standardfiles["customfile12"] = ""
	self.standardfiles["customfile13"] = ""
	self.standardfiles["customfile14"] = ""
	self.standardfiles["customfile15"] = ""
	self.standardfiles["customfile16"] = ""
	self.standardfiles["customfile17"] = ""
	self.standardfiles["customfile18"] = ""
	self.standardfiles["customfile19"] = ""
	self.standardfiles["customfile20"] = ""
	self.readOptions("options/standardfiles.txt",self.standardfiles)

	#remove blank entries to prevent key copy errors
	for key,value in self.standardfiles.iteritems():
	    if (value==""):
		try: 
		    del self.standardfiles[key]
		except KeyError:
		    pass

	#error checking for standardfiles
#        for key,value in self.standardfiles.iteritems():
#            if not os.path.exists(os.path.join(os.getcwd(),mode,"standardfiles",value)):
#                logger.error(("%s not found!") % (key))
#                raise SystemExit()

	 #for output files that need to be copied
	self.copyfiles1 = odict()
	self.copyfiles1["trajfile"] = "traj.trr"
	self.copyfiles1["trajfile2"] = "traj.xtc"
	self.copyfiles1["indexfile"] = "index.ndx"
	self.copyfiles1["topolfile"] = "topol.tpr"
	
	#for output files that need to be copied
	self.copyfiles2 = odict()
	self.copyfiles2["trajfile"] = "traj.trr"
	self.copyfiles2["trajfile2"] = "traj.xtc"
	self.copyfiles2["indexfile"] = "index.ndx"
	self.copyfiles2["topolfile"] = "topol.tpr"
	self.copyfiles2["conffile"] = "conf.gro"
	self.copyfiles2["conffile"] = "traj.gro"

	#set III
	self.copyfiles3 = self.copyfiles2.copy()

	#set IV
	self.copyfiles4 = odict()
	self.copyfiles4["trajfile"] = "traj.trr"
	self.copyfiles4["indexfile"] = "index.ndx"
	self.copyfiles4["topolfile"] = "topol.tpr"
	self.copyfiles4["pathfile"] = "traj.gro"


	#set for deleting
	self.deletefiles = odict()
	self.deletefiles["engfile"] = "ener.edr"

	#for file extensions
	self.extension = ".gro"


##########################Run options#########################
	self.runoptions = odict()
	self.runoptions["mdbinary"]= "gmx"
	self.runoptions["updatetime"] = 5
	self.runoptions["dvmax"] = 1.0
	self.runoptions["scdvmax"] = 0.01
	#self.runoptions["gromacssuffix"] = ""
	self.runoptions["shootfrominterface"] = 1
	self.runoptions["queuename"] = "serial"      
	self.runoptions["qsubpath"] = ""
	self.runoptions["maxlength"] = 1000        
	self.runoptions["dist1"] = "1"
	self.runoptions["dist2"] = "3"      
	self.runoptions["qsubsystem"] = "local"    
	self.runoptions["qsubwalltime"] = "100:00:00"   
	self.runoptions["interfacecoordinate"] = "0"
	self.runoptions["mddriver"] = ""

	self.runoptions["shoot"] = 45
	self.runoptions["swap"] = 45
	self.runoptions["bfswap"] = 10
	self.runoptions["initialopeval"] = 'True'
	self.runoptions["jobname"] = ''
	
	self.readOptions("options/runoptions.txt",self.runoptions)

	#fix swapping moves and shooting moves
	self.runoptions["shoot"] = int(self.runoptions["shoot"])
	self.runoptions["swap"] = int(self.runoptions["swap"])
	self.runoptions["bfswap"] = int(self.runoptions["bfswap"])

	#error checking for runoptions
	#mdbinary check
	if not os.path.exists(os.path.join(self.paths["mdbinpath"],self.runoptions["mdbinary"])):
	    logger.error("md binary not found!")
	    raise SystemExit()

	qsubsystemlist = ['Standard','PBS','lisa','LL','local','vulcan','vulcan-parallel']

	if self.runoptions["qsubsystem"] not in qsubsystemlist:
	    logger.error("unknown qsubsystem entry found!")
	    raise SystemExit() 


	

##################SRM addition###############################################
	self.dumpfiles = odict()
	self.dumpfiles["conffile"] = "conf.gro"
	self.dumpfiles["bakfile"] = "backconf.gro"
	self.dumpfiles["bconffile"] = "bconf.gro"


####################MD init options ########################
	#GDLcomment: names for inputs of MD should go to MD driver specific
	self.initoptions = odict()
	self.initoptions["grofile"] = "conf.gro"
	self.initoptions["topfile"] = "topol.top"
	self.initoptions["mdpfile"] = "md.mdp"
	self.initoptions["ndxfile"] = "index.ndx"
	self.initoptions["tprfile"] = "topol.tpr"
	self.initoptions["pullgroupname"] = "DOPC"
	self.initoptions["pullgroupnumber"] = "1"
	self.initoptions["refgroupname"] = "1VAL"
	self.initoptions["refgroupnumber"] = "3"

#GDL: taking away mdrun options to use a file        
	self.mdpoptions = odict()
	self.mdpoptions["title"]                 = "TPS" 
	self.mdpoptions["cpp"]                   = "/usr/bin/cpp" 
	self.mdpoptions["integrator"]            = "md" 
	self.mdpoptions["tinit"]                 = "0.0" 
	self.mdpoptions["dt"]                    = "0.0050" 
	self.mdpoptions["nsteps"]                = "90000000" 
	self.mdpoptions["nstcomm"]               = "1" 
	self.mdpoptions["comm-mode"]               = "Linear"
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
	self.mdpoptions["domain-decomposition"]  = "no"
	self.mdpoptions["coulombtype"]           = "pme"
	self.mdpoptions["rcoulomb_switch"]       = "0.0" 
	self.mdpoptions["rcoulomb"]              = "1.2" 
	self.mdpoptions["epsilon_r"]             = "15" 
	self.mdpoptions["vdw_type"]              = "Shift"  
	self.mdpoptions["rvdw_switch"]           = "0.9" 
	self.mdpoptions["rvdw"]                  = "1.2" 
	self.mdpoptions["DispCorr"]              = "No" 
	self.mdpoptions["tcoupl"]                = "nose-hoover" 
	self.mdpoptions["nsttcouple"]            = "-1"
	self.mdpoptions["nh-chain-length"]       = "10"
	self.mdpoptions["tc-grps"]               = "System" 
	self.mdpoptions["tau_t"]                 = "0.3" 
	self.mdpoptions["ref_t"]                 = "323"  
	self.mdpoptions["Pcoupl"]                = "berendsen"  
	self.mdpoptions["Pcoupltype"]            = "semiisotropic"
	self.mdpoptions["nstpcouple"]            = "-1"
	self.mdpoptions["tau_p"]                 = "0.75" 
	self.mdpoptions["compressibility"]       = "10e-5 10e-5" 
	self.mdpoptions["ref_p"]                 = "1.0 1.0" 
	self.mdpoptions["gen_vel"]               = "no" 
	self.mdpoptions["gen_temp"]              = "323" 
	self.mdpoptions["gen_seed"]              = "665" 
	self.mdpoptions["constraints"]           = "none" 
	self.mdpoptions["constraint_algorithm" ] = "Lincs" 
	self.mdpoptions["continuation"]   = "no" 
	self.mdpoptions["lincs_order"]           = "4" 
	self.mdpoptions["lincs_warnangle"]       = "30" 
	self.mdpoptions["morse"]                 = "no"
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
 
	

    #def writeQsubFileParallel(self,dir,cores):
    #    nodes = int(ceil(float(self.gdiststeps)/float(cores)))
    #    print nodes
    #    
    #    for i in range(nodes):
    #        filename = os.path.join(dir,'psubmit.'+str(i)+'.sub')
    #        f=open(filename, 'w')
    #        f.write('#!/bin/bash\n')
    #        f.write('#$ -N pygus'+str(i) +'\n')
    #        f.write('#$ -S /bin/bash\n')
    #        f.write('#$ -o pbs.'+str(i) + '.out\n')
    #        f.write('#$ -e pbs.'+str(i) + '.err\n')
    #        f.write('#$ -r n\n')
    #        f.write('#$ -V\n')
    #        f.write('#$ -q parallel\n')
    #        f.write('#$ -pe mpi %d\n' % (cores,))
    #        f.write('#$ -cwd\n')
    #        for j in range(cores):
    #            f.write('cd %02d\n' % (i*cores+j,))
    #            f.write('mdrun -s umbrella%d.tpr &\n' % (i*cores+j,))
    #            f.write('cd ..\n')
    #            if i*cores+j >= self.gdiststeps-1:
    #                break
    #        f.close()
	
 #########GDL: functions read/write/evaluate options#########################        

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
	    
    #def _makepaths(self):
    #    for key,item in self.paths.iteritems():
    #        if not (key in ["gromacspath"]):
    #            if not(os.path.exists(item)):
    #                os.mkdir(item)
		     
    #def writeMdpFile(self,directory,sfilename):
    #    filename = os.path.join(directory,sfilename)
    #    of = open(filename,"w")
    #    for key,item in self.mdpoptions.iteritems():
    #        of.write('%-15s = %s\n' % (key, item))
    #    of.close()
	
    #def writeQsubFile(self,directory,sfilename,additional):
    #    filename = os.path.join(directory,sfilename)
    #    of = open(filename,"w")
    #    for line in self.qsuboptions:
    #        of.write(line + "\n")
    #    for line in additional:
    #        of.write(line + "\n")
    #    of.close()
    
    #def mcopyFile(self,destdir,destfile,targetdir,targetfile):
    #    df = os.path.join(destdir,destfile)
    #    tf = os.path.join(targetdir,targetfile)
    #    os.system("cp " + df + " " + tf)
	
	
	
	
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
