"""
Test set of functions.
by SRM
"""


import os
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.runrecord.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False   


#SRM:function to read mddriver becuase we need it initially
def  read_options(filename,option):
    if os.path.exists(filename):
        for line in open(filename,"r"):
            raw = line.split("=")
            if len(raw) == 2:
                if raw[0] == option:
                    return raw[1]
 
#SRM: function to process paths
def find_paths():
    #find out paths here
    #opens paths.txt
    basepath = os.getcwd()

    masterdir_path = os.path.join(basepath,"paths.txt")
    
    #check if masterdir_path is okay
    if not os.path.exists(masterdir_path):
        logger.error("paths.txt not found!")
        raise SystemExit()

    #check options
    if not os.path.exists(os.path.join(os.getcwd(),"options")):
        logger.error("folder options not found!")
        raise SystemExit()

    #opens runoptions.txt
    ropt_path = os.path.join(basepath,"options","runoptions.txt")

    if not os.path.exists(ropt_path):
        logger.error("runoptions.txt not found!")
        raise SystemExit()

    md_driver = read_options(ropt_path,"mddriver")
    
    #SRM: errorcheck for mddriver
    checkd = md_driver.strip()
    print checkd
    
    if ((checkd != "lammps") and (checkd!="gromacs")):
        logger.error("mddriver is not lammps or gromacs. Please check the entry in runoptions.txt")
        raise SystemExit()
    
    #finds path of the master code
    md_path = read_options(masterdir_path,"wrapperpath")
    md_path = md_path.strip()
    #check for md path
    if not os.path.exists(md_path):
        logger.error("wrapper path is wrong. Please check the entry in paths.txt")
        raise SystemExit()
    
    #modify the path to suit us
    module_path = os.path.join(md_path,"pygromacstps","mddriver")
    custom_path = os.path.join(md_path,"ClusterGTIS","custom")
    module_path = module_path.strip()
    custom_path = custom_path.strip()
    
    #print module_path
    module_path = os.path.join(module_path,md_driver)
    custom_path = os.path.join(custom_path,md_driver)
    module_path = module_path.strip()
    custom_path = custom_path.strip()
    
    #print module_path
    #now we can pass module_path to whichever functions and 
    #open the necessary stuff as we wish
    return module_path,custom_path