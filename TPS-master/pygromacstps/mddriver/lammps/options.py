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
	self.paths = odict()
        self.paths["binpath"] = "/Users/rogal/bin"
        self.paths["tispath"] = ""
        self.paths["wrapperpath"] = "/Users/rogal/work"
        self.paths["mdbinpath"] = "/home/users/diazlgjj/bin"
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
        #old code
        #self.standardfiles = []
        #for line in open(os.path.join(self.basedir,"options","standardfiles.txt"),"r"):
        #    self.standardfiles.append(line.rstrip('\n'))
        
        #    SRM edit
        self.standardfiles = odict()
        self.standardfiles["conffile"] = "conf.dump"
        self.standardfiles["mdinfile"] = "md.in"
        self.standardfiles["potentialfile"] = "Ni_u3.eam"
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


        #for output files that need to be copied:set I
        self.copyfiles1 = odict()
        self.copyfiles1["trajfile"] = "traj.dat"
        self.copyfiles1["conffile"] = "conf.dump"
        self.copyfiles1["mdinfile"] = "md.in"

        #for output files that need to be copied:set II
        self.copyfiles2 = self.copyfiles1.copy()

        #set III
        self.copyfiles3 = odict()
        self.copyfiles3["trajfile"] = "traj.dat"

        #set IV
        self.copyfiles4 = self.copyfiles3.copy()

        #set for deleting
        self.deletefiles = odict()

        #for file extensions
        self.extension = ".dump"



##########################Run options##############################
        self.runoptions = odict()
        self.runoptions["mdbinary"] = ""
        self.runoptions["updatetime"] = 5
        self.runoptions["dvmax"] = 1.0
        self.runoptions["scdvmax"] = 0.01
#       self.runoptions["gromacssuffix"] = ""
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
        self.runoptions["parallel"] = "False"

        self.runoptions["shoot"] = 45
        self.runoptions["swap"] = 45
        self.runoptions["bfswap"] = 10
        self.runoptions["initialopeval"] = 'True'
        self.runoptions["jobname"] = 'ktis'
        self.runoptions["pathcut"] = 'False'
        self.runoptions["statfrequency"] = 1
        self.runoptions["sigma"] = 0.2
        self.runoptions["minprob"] = 0.001
        self.runoptions["minacceptance"] = 0.40

        self.readOptions("options/runoptions.txt",self.runoptions)


        #error checking for runoptions
        #mdbinary check
        if not os.path.exists(os.path.join(self.paths["mdbinpath"],self.runoptions["mdbinary"])):
            logger.error("md binary not found!")
            raise SystemExit()

        qsubsystemlist = ['Standard','PBS','lisa','LL','local','vulcan','vulcan-parallel','slurm','slurm-parallel']

        if self.runoptions["qsubsystem"] not in qsubsystemlist:
            logger.error("unknown qsubsystem entry found!")
            raise SystemExit()

        #fix swapping moves and shooting moves
        



####################MD init options ########################   
        #GDLcomment: names for inputs of MD should go to MD driver specific
        self.initoptions = odict()
        self.initoptions["conffile"] = "conf.dump"
        self.initoptions["mdinfile"] = "md.in"
        self.initoptions["pathfile"] = "traj.dat"
        
        
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


##################SRM addition###############################################
        self.dumpfiles = odict()
        self.dumpfiles["conffile"] = "conf.dump"
        self.dumpfiles["bakfile"] = "backconf.dump"
        self.dumpfiles["bconffile"] = "bconf.dump"

 
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
