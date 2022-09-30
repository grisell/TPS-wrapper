
import os
import subprocess as sub
import random
import shutil
import sys




class wrapper(object):
    def __init__(self):
        self.wolf = "test"



    def generateMDRunCommand(self,options,commandpath):
	if options.runoptions["mddriver"] == "lammps":
		if options.runoptions["parallel"]=="True":
			if options.runoptions["qsubsystem"] == "vulcan-parallel":
				cmd = []	
				cmd.append(os.path.join(options.paths["openmpipath"],'mpirun'))
				cmd.append('-np')
				#cmd.append(os.environ['NSLOTS'])
				#srm edit to fix double number of jobs running
				slots = os.environ['NSLOTS']
				slots = int(slots)/2
				cmd.append(str(slots))
				cmd.append('-machinefile')
				cmd.append(os.path.join(os.environ['TMPDIR'],'machines'))
				cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        			cmd.append('-in')
        			cmd.append('md.in')
#				cmd.append(options.initoptions["mdinfile"][1]) 
				cmd.append('-screen')
				cmd.append('lammscreen.log')

                        
                        elif options.runoptions["qsubsystem"] == "slurm-parallel":
                                cmd = []        
                                cmd.append('mpirun')
                                cmd.append('-np')
                                dum = options.runoptions["queuename"].split('-')
                                nslots = int(dum[1])
                                slots = nslots/2
                                cmd.append(str(slots))
                                #cmd.append('-machinefile')
                                #cmd.append(os.path.join(os.environ['TMPDIR'],'machines'))
                                cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
                                cmd.append('-in')
                                cmd.append('md.in')
#                               cmd.append(options.initoptions["mdinfile"][1]) 
                                cmd.append('-screen')
                                cmd.append('lammscreen.log')

                        elif options.runoptions["qsubsystem"] == "local":
                                cmd = []        
                                cmd.append('mpirun')
                                cmd.append('-np')
                                dum = options.runoptions["queuename"].split('-')
                                nslots = int(dum[1])
                                slots = nslots/2
                                cmd.append(str(slots))
                                #cmd.append('-machinefile')
                                #cmd.append(os.path.join(os.environ['TMPDIR'],'machines'))
                                cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
                                cmd.append('-in')
                                cmd.append('md.in')
#                               cmd.append(options.initoptions["mdinfile"][1]) 
                                cmd.append('-screen')
                                cmd.append('lammscreen.log')


			else:
				print 'Parallel lammps run works only for sge on vulcan'
				print 'qsubsystem must be vulcan-parallel'
				print 'Exiting program...'
				sys.exit(1)
		else:
        		cmd = []
        		cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
        		cmd.append('-in')
        		cmd.append('md.in')
#			cmd.append(options.initoptions["mdinfile"][1]) 
			cmd.append('-screen')
			cmd.append('lammscreen.log')

	elif options.runoptions["mddriver"] == "ljmd":
		cmd = []
		cmd.append(os.path.join(options.paths["mdbinpath"],options.runoptions["mdbinary"]))
		cmd.append('md.in')


	return cmd



    def TrajToConf(self,options,tmpdirpath):
	os.chdir(tmpdirpath)
	filename = os.path.join(tmpdirpath,options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
	if len(data) > 0:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock

	for i in range(nslices):
		fout = "path%d.dump" % i
		filename = os.path.join(tmpdirpath,fout)
		wf = open(filename,"w")

		start = i*nblock
		wf.write(data[start])
		wf.write("0\n")
		for j in range(2,nblock):
			wf.write(data[start+j])
    

    def GetSliceFromTraj(self,filename,outfile,rn):
        """
        Extract a prticular slice from the trajectory

        Parameters
        ----------
        filename : string 
            name of the input file

        rn : int
            slice number

        gzip : Trur or False
            True if the file is compressed, False otherwise
        """
        natoms = 0
        gz = open(filename,'rb')
        f = gz

        
        count=0
        for line in f:
            if count==3:
                natoms = int(line.strip())
                break
            count+=1
        gz.close()

        #now re read with necessary info
        startread= (natoms+9)*rn
        endread = (natoms+9)*(rn+1)


        fout = open(outfile,'w')
        gz = open(filename,'rb')
        f = gz

        count=0
        
        for line in f:
            if count>=startread:
                if count==startread+1:
                    fout.write("0\n")
                else:
                    fout.write(line)
                    #print line
            
            count+=1
            if count==endread:
                break

        gz.close()
        fout.close()

    def modifyMDseed(self,wdir,options):
	    filename = os.path.join(wdir,options.initoptions["mdinfile"])
	    data = []
	    if os.path.exists(filename):
		    infile = open(filename,'r')
		    for line in infile:
			    data.append(line)
		    infile.close()
		    
		    newfile = os.path.join(wdir,'newmd.in')
		    outfile = open(newfile,'w')
		    nlines = len(data)
		    raw = []
		    for i in range(nlines):
			    raw.append(data[i].split())

		    for i in range(nlines):
			if(len(raw[i]) > 3):
				if raw[i][3] == 'langevin':
					ncolumn = len(raw[i])
					seed = random.randint(0,1000000)
					for k in range(7):
						outfile.write(raw[i][k])
						outfile.write(" ")
					outfile.write(str(seed))
					outfile.write(" ")
					if ncolumn > 8:
			 			for k in range(8,ncolumn):
							outfile.write(raw[i][k])
							outfile.write(" ")
					outfile.write("\n")
                                elif raw[i][3] == 'temp/bussinvt':                      #DS    
                                        ncolumn = len(raw[i])
                                        seed = random.randint(-1000000,-1)
                                        for k in range(6):
						outfile.write(raw[i][k])
						outfile.write(" ")
                                        outfile.write(str(seed))
                                        outfile.write(" ")
                                        if ncolumn > 7:
			 			for k in range(7,ncolumn):
							outfile.write(raw[i][k])
							outfile.write(" ")
					outfile.write("\n")   
                                    
				else:
			   		outfile.write(data[i])
			else:
				outfile.write(data[i])
	
		    outfile.flush()
		    outfile.close()
		  
		    shutil.move(newfile,filename)




    def executeCommand(self,cmd,tinput=""):
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        out,err = proc.communicate(input=tinput)
        proc.wait()
        return out,err


    def initializeMD(self,path,basedir):
    	#nothing to do here
    	return 0

    
    def cutTrajFile(self,fdir,bdir,fl,bl):
    	
    	filename = os.path.join(fdir,'traj.dat')

    	with open(filename) as myfile:
    		data = [next(myfile) for x in xrange(4)]
    	natoms = int(data[3])
    	nblock = natoms + 9

    	N = fl*nblock
    	with open(filename) as myfile:
    		head = [next(myfile) for x in xrange(N)]
    	os.remove(filename)
    	with open(filename,'w') as fout:
    		for line in head:
    			fout.write(line)

    	filename = os.path.join(bdir,'traj.dat')
    	N = bl*nblock
    	with open(filename) as myfile:
    		head = [next(myfile) for x in xrange(N)]
    	os.remove(filename)
    	with open(filename,'w') as fout:
    		for line in head:
    			fout.write(line)

