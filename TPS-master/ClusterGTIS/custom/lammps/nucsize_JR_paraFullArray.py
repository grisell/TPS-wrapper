'''
Created on Oct 26, 2010

@author: wolf
'''

#import pygromacstps.wrappers
import pygromacstps.parser
import pygromacstps.filesystem
import math

import os
import sys
import time
import shutil
import subprocess as sub
import multiprocessing as mp
import copy_reg
import types
import wrappers as wrappers

# The path where the program water_around is located
#wapath = "/home/cvreede/GTIS/pypwater/watertest"

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

def return_result(result):
    return result

copy_reg.pickle(types.MethodType, _pickle_method)


class TCustomOrderParameter(object):
    def __init__(self,minA,maxA,minB,maxB):
        # The stable states
        self.A = [minA,maxA]
        self.B = [minB,maxB]
	        # dircolumn is the column from the output where the interesting coordinated is
	# JR: might not be needed or we might extend this
	#self.dircolumn = 7
	#SRM code changes        
	
	    
        self.wrapper    = wrappers.wrapper()
	
        self.distparser = pygromacstps.parser.gdistparser()
        self.filesystem = pygromacstps.filesystem.filesystem()



    def _makecalcOP(self,path,workdir):
	cmd = []
	cmd.append(os.path.join(path.options.paths["binpath"],"main"))
	cmd.append(os.path.join(workdir,"tmp_conf"))
	return cmd

    def _makecalcOP_file(self,path,tmpfile):
	cmd = []
	cmd.append(os.path.join(path.options.paths["binpath"],"main"))
	cmd.append(tmpfile)
	return cmd
        
    
        
    def isPathInState(self,path,log):
        # Does the path end in one of the two stable states?
        # get trajectory and look at the last coordinate
	# get the last slice of current trajectory
	filename = os.path.join(path.workdir,path.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
		
		if len(data) > 4:
			natoms = int(data[3])
		else:
			return 0,False
		#check if file is completely written
		if len(data)%(natoms+9) != 0:
			return 0,False
		if len(data) > natoms:
			#write last slice to tmp file
			tmpfile = os.path.join(path.workdir,"tmp_conf")
			outfile = open(tmpfile,"w")
			for i in range(-(natoms+9),0):
				outfile.write(data[i])
			outfile.flush()
			outfile.close()
			#calculate order parameter
			#using external program...
			cmd = self._makecalcOP(path,path.workdir)
			output,outerr = self.wrapper.executeCommand(cmd)
			q = float(output)
			os.remove(tmpfile)
		else:
			return 0,False


		if q > self.A[0] and q < self.A[1]:
                	return 1,True
            	elif q > self.B[0] and q < self.B[1]:
                	return 4,True
            	else:
                	return 0,False

	else:
	  	return 0,False

				
				
 
    def getQTrajectory(self,path,wdir,data):
	nlines = len(data)
	natoms = int(data[3])
	nblock = natoms+9
	nslices = nlines/nblock

	if nlines%nblock != 0:
		print 'Error in getQTrajectory'
		print 'nlines%nblock != 0'
		print 'Exiting program...'
		sys.exit(1)

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
		#print cmd
		output,outerr = self.wrapper.executeCommand(cmd)
		q = float(output)

		qtraj.append(q)

	os.remove(tmpfile)
	
	return qtraj




    def getFullTrajectory(self,fpath,bpath,fdir,bdir,log):
        #JR add logger information
	log.log.debug("JR: Getting full trajectory serial")    
	#JR end
        traj = []
	qtrajback = []
	qtrajfor = []
	flength=0
	blength=0
	#first the backward part
        os.chdir(bdir)
	filename = os.path.join(bdir,bpath.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'backward part of trajectory not fully written in workdir, removing unfinished slice...'
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
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)

	if len(data) > 0:
		qtrajback = self.getQTrajectory(bpath,bdir,data)
		
        	qtrajback.reverse()
        	blength=0
        	count = 0
#		print 'traback, data',len(qtrajback),len(data)
#		sys.stdout.flush()
        	for i in range(len(qtrajback)):
            		traj.append([count,qtrajback[i],0])
            		blength+=1
            		count += 1
        
	#now the forward part
        os.chdir(fdir)
	filename = os.path.join(fdir,fpath.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'forward part of trajectory not fully written in workdir, removing unfinished slice...'
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
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)


	

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
        
        return traj,flength,blength



    
    #JR : parallel evaluation of op
    def getQTrajectoryParallel(self,data,instance,kslices,path,wdir):
	#JR add logger
	#print "JR: now in getQtrajectoryParallel"
	#JR end add logger
	nlines = len(data)
	natoms = int(data[3])
	nblock = natoms + 9
	nslices = nlines/nblock

	if nlines%nblock != 0:
		print 'Error in getQTrajectoryParallel'
		print 'nlines%nblock != 0'
		print 'Exiting program...'
		sys.exit(1)
	qtraj = []
	for j in range(kslices[0],kslices[1]):
		start = j*nblock
		end = (j+1)*nblock
		#JR print
		#print start, end
		#sys.stdout.flush()
		#JR end print
		tmpfile = os.path.join(wdir,"my_tmp_conf" + str(instance))
		#JR print
		#print tmpfile
		#sys.stdout.flush()
		#JR end print
		outfile = open(tmpfile, "w")
		for i in range(start,end):
			outfile.write(data[i])
		outfile.flush()
		outfile.close()
		cmd = self._makecalcOP_file(path,tmpfile)
		#JR print
		#print cmd
		#sys.stdout.flush()
		#JR end print
		output,outerr = self.wrapper.executeCommand(cmd)
		q = float(output)
		#JR print
		#print q
		#sys.stdout.flush()
		#JR end print
		qtraj.append(q)

	os.remove(tmpfile)

	return instance,qtraj


    def getFullTrajectoryParallel(self,fpath,bpath,fdir,bdir,log):
	#JR add logger information
	log.log.debug("JR: Getting full trajectory parallel")    
	#JR end
        traj = []
	qtrajback = []
	qtrajfor = []
	flength=0
	blength=0
	#JR: this needs to be modified to use proper number of cores!
	#cores=mp.cpu_count()-1
	#JR get number of cores from the queuing system on vulcan, determine cores if 'local', exit otherwise
	if bpath.options.runoptions["qsubsystem"] == "vulcan-parallel":
		cores = int(os.environ['NSLOTS'])
	elif bpath.options.runoptions["qsubsystem"] == "local":
		cores = mp.cpu_count()-1
	else:
		print 'Error in getFullTrajectoryParallel:'
		print 'Parallel evaluation of the order parameter only works for qsubsystem vulcan-parallel and local'
		print 'Exiting progam...'
		sys.exit(1)
	#JR add logger information
	log.log.debug("JR: number of cores:" + str(cores))    
	#JR end

	#first the backward part
        os.chdir(bdir)
	results = []
        kernel_list = [ 0 for x in range(cores)]
	
	filename = os.path.join(bdir,bpath.options.initoptions["pathfile"])
	data = []
	#JR add logger
	#log.log.debug("JR: bw filename:" + filename)
	#JR end
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

	#JR add logger
	#log.log.debug("JR: bw length of data array:" + str(len(data)))	
	#JR end
	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'backward part of trajectory not fully written in workdir, removing unfinished slice...'
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
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)

	if len(data) > 0:
    		per_core = nslices/cores
		#JR: in case nslices < ncore
		if (per_core == 0):
			kernel_list = [ 0 for x in range(nslices)]

		#JR add logger
		log.log.debug("JR: bw: nslices = " + str(nslices) + " cores = " + str(cores) + " per_core = " + str(per_core))
		#JR end logger	
		for j in range(len(kernel_list)):
			kernel_list[j] = per_core
		for j in range(nslices%cores):
			kernel_list[j] += 1
		#kernel_list contains no. of slices for each core
		#assigen start and end for the loop on each core in kernel_slices
		start = 0
		end = 0
		kernel_slices = []
		for j in range(len(kernel_list)):
			end += kernel_list[j]
			kernel_slices.append([start,end])
			start += kernel_list[j]
		#JR add logger
		for j in range(len(kernel_slices)):
			log.log.debug("JR: kernel_slices[" + str(j) + "] = " + str(kernel_slices[j]))
		#JR end logger	
        	
		#now launch process
		pool = mp.Pool(processes=len(kernel_list))
		for x in range(len(kernel_list)):
			results.append(pool.apply_async(self.getQTrajectoryParallel, args=(data,x,kernel_slices[x],bpath,bdir,)))
		pool.close()
		pool.join()
		output = [p.get() for p in results]
		#sort the results
		output.sort()
		#trim the unnecessary variable
		output = [out[1] for out in output]
		#join output
		output = sum(output,[])

  		#give it to qtrajback
  		qtrajback = output
  	      	qtrajback.reverse()
        	blength=0
        	count = 0
#		print 'traback, data',len(qtrajback),len(data)
#		sys.stdout.flush()
        	for i in range(len(qtrajback)):
            		traj.append([count,qtrajback[i],0])
            		blength+=1
            		count += 1

	#now the forward part
        os.chdir(fdir)
	results=[]
	kernel_list = [ 0 for x in range(cores)]

	filename = os.path.join(fdir,fpath.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'forward part of trajectory not fully written in workdir, removing unfinished slice...'
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
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)

	if len(data) > 0:
    		per_core = nslices/cores
		#JR: in case nslices < ncore
		if (per_core == 0):
			kernel_list = [ 0 for x in range(nslices)]

		#JR add logger
		log.log.debug("JR: fw: nslices = " + str(nslices) + " cores = " + str(cores) + " per_core = " + str(per_core))
		#JR end logger	
		for j in range(len(kernel_list)):
			kernel_list[j] = per_core
		for j in range(nslices%cores):
			kernel_list[j] += 1
		#kernel_list contains no. of slices for each core
		#assigen start and end for the loop on each core in kernel_slices
		start = 0
		end = 0
		kernel_slices = []
		for j in range(len(kernel_list)):
			end += kernel_list[j]
			kernel_slices.append([start,end])
			start += kernel_list[j]
		#JR add logger
		for j in range(len(kernel_slices)):
			log.log.debug("JR: kernel_slices[" + str(j) + "] = " + str(kernel_slices[j]))
		#JR end logger	
        	
		#now launch process
		pool = mp.Pool(processes=len(kernel_list))
		for x in range(len(kernel_list)):
			results.append(pool.apply_async(self.getQTrajectoryParallel, args=(data,x,kernel_slices[x],fpath,fdir,)))
		pool.close()
		pool.join()
		output = [p.get() for p in results]
		#sort the results
		output.sort()
		#trim the unnecessary variable
		output = [out[1] for out in output]
		#join output
		output = sum(output,[])

  		#give it to qtrajfor
  		qtrajfor = output
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
        
        return traj,flength,blength



    def getFullTrajectoryReverse(self,fpath,fdir,log):
        traj = []
	qtrajfor = []
	flength=0
        
	#now the forward part
        os.chdir(fdir)
	filename = os.path.join(fdir,fpath.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)
	
	if len(data) > 4:
		nlines = len(data)
		natoms = int(data[3])
		nblock = natoms+9
		nslices = nlines/nblock
		if nlines%nblock != 0:
			print 'reverse trajectory not fully written in workdir, removing unfinished slice...'
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
			natoms = int(data[3])
			nblock = natoms+9
			nslices = nlines/nblock

			if nlines%nblock != 0:
				print 'Cannot write trajectory file with proper number of lines'
				print filename
				print 'Exiting program...'
				sys.exit(1)	

	if len(data) > 0:
		qtrajfor = self.getQTrajectory(fpath,fdir,data)
        
        	flength=0
		count=0
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
	filename = os.path.join(path.workdir,path.options.initoptions["pathfile"])
	data = []
	if os.path.exists(filename):
		infile = open(filename,"r")
		for line in infile:
			data.append(line)

		if len(data) > 4:
			natoms = int(data[3])
		else:
			return 0,0
		
	        nlines = len(data)
		nblock = natoms+9
		nslices = nlines/nblock
		#check if file is completely written
		if nlines%nblock != 0:
			return 0,0

		if len(data) > natoms:
			#write last slice to tmp file
			tmpfile = os.path.join(path.workdir,"tmp_conf")
			outfile = open(tmpfile,"w")
			for i in range(-(natoms+9),0):
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

				

    
    
 
