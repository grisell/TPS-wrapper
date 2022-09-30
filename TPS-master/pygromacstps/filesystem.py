'''
Created on Mar 2, 2010

@author: wolf
'''
import os
import glob
import time

class filesystem(object):
    def __init__(self):
        self.data = []
        
    
    def createDirList(self,directories,log):
        for directory in directories:
            if not os.path.exists(directory):
                os.system("mkdir " + directory)
                #log.debug(directory+ " created")
        
                
                
    def createDir(self,directory):
        if not os.path.exists(directory):
            os.system("mkdir " + directory)
            
    
    def copyFiles(self,source,destination):
        if os.path.exists(source):
            os.system("cp " + source + " " + destination)
        else:
            print source + " doesn't exist."
    
    
    def copyFileList(self,source,destination,wildcards):
        if os.path.exists(source):
            if os.path.exists(destination):
                flist = glob.glob(os.path.join(source,wildcards))
                #os.system("cp " + sourcef + " " + destination)
                for sourcef in flist:
                    os.system("cp " + sourcef + " " + destination)
            else:
                print destination + " doesn't exist."
        else:
            print source + " doesn't exist."
    
    def deleteFileList(self,source,wildcards):
        if os.path.exists(source):
            flist = glob.glob(os.path.join(source,wildcards))
            for sourcef in flist:
                os.system("rm " + sourcef)
        else:
            print source + " doesn't exist."
    
    def deleteFile(self,filename):
        if os.path.exists(filename):
            os.remove(filename)
    
            
    def getFileList(self,source,wildcards):
        flist = []
        if os.path.exists(source):
            flist = glob.glob(os.path.join(source,wildcards))
        else:
            print source + " doesn't exist"
        return flist

    def moveFile(self,source,destination):
        if os.path.exists(source):
            os.system("mv " + source + " " + destination)
        else:
            print source + " doesn't exist."
        
            
    def clearTempFiles(self,directory):
        os.system('rm '+ directory+'/\#*')
        
    
    def waitForFile(self,filename):
        while True:
            if os.path.exists(filename):
                break
            else:
                time.sleep(20)
    
         
        
        
        