
import os
import sys
import logging

#SRM:set up logger for the general error messages
logger = logging.getLogger(__name__)
handler = logging.FileHandler("gtpslog.analysis.txt")
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.propagate = False

#SRM:import the module for selecting lammps/gromacs
import pygromacstps.findmd as md

#SRM:find the module path
module_path,custom_path = md.find_paths()

#SRM:check if module path is okay and add it to list of paths
if os.path.exists(module_path):
    sys.path.append(module_path)
else:
    logger.error("lammps/gromacs specific modules not found in the code directory")
    raise SystemExit()

#SRM:check if custom path is okay and add it to list of paths
if os.path.exists(custom_path):
    sys.path.append(custom_path)
else:
    logger.error("lammps/gromacs specific order parameter files not found in the code directory")
    raise SystemExit()

#imports here
import GTISAnalysis.CreateCrossingHistogramsClass as cc
import GTISAnalysis.SelfConsistantHistoClass as sh


# a function to read values
def readOptions(filename,optiondict):
	if os.path.exists(filename):
		for line in open(filename,"r"):
			raw = line.split("=")
			if len(raw) == 2:
				for key,value in optiondict.iteritems():
					if raw[0] == key:
						optiondict[key]=raw[1].strip()



if __name__=="__main__":
	
	options = empty_dict = dict.fromkeys(['analysis','subanalysis','startpath','stoppath','histosize','histomax','histomin','precision','cutoff','temperature'])
	readOptions("options/analysisoptions.txt",options)

	start = int(options['startpath'])
	stop = int(options['stoppath'])
	histosize = int(options['histosize'])
	analysis = options['analysis']
	subanalysis = options['subanalysis']

	logger.info(("Analysis selected : %s %s")%(analysis,subanalysis))

	if (options["analysis"]=="wham"):
		if (options["subanalysis"]=="freeenergy"):
			histomin = float(options['histomin'])
			histomax = float(options['histomax'])
			precision = float(options['precision'])
			cutoff = float(options['cutoff'])
			temperature = float(options['temperature'])
			basedir = os.path.join(os.getcwd())
			crossinghisto = cc.crossingHistoData(histosize,histomax,histomin,start,stop,basedir,"tis")
			crossinghisto.run_main()
			sim = sh.TSelfConsistent(histosize,precision,cutoff,temperature,subanalysis)
			sim.run_main()
		elif (options["subanalysis"]=="crossinghisto"):
			histomin = float(options['histomin'])
			histomax = float(options['histomax'])
			precision = float(options['precision'])
			cutoff = float(options['cutoff'])
			temperature = float(options['temperature'])
			basedir = os.path.join(os.getcwd())
			crossinghisto = cc.crossingHistoData(histosize,histomax,histomin,start,stop,basedir,"tis")
			crossinghisto.run_main()
			sim = sh.TSelfConsistent(histosize,precision,cutoff,temperature,subanalysis)
			sim.run_main()


	elif (options["analysis"]=="histogram"):
		histomin = float(options['histomin'])
		histomax = float(options['histomax'])
		basedir = os.path.join(os.getcwd())
		crossinghisto = cc.crossingHistoData(histosize,histomax,histomin,start,stop,basedir,"tis")
		crossinghisto.run_main()

	

    #stuff if plot is true
    	
    	






