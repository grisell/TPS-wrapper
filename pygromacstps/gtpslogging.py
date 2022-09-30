'''
Created on May 6, 2010

@author: wolf
'''
import logging
import os

class log(object):
    def __init__(self,levelname,basedir,kernel):

            
        LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}
        self.level = LEVELS.get(levelname, logging.NOTSET)
        self.log = logging
        if kernel =="head":
            self.filename = os.path.join(basedir,"gtpslog.head.txt")
        elif kernel == "reverse":
            self.filename = os.path.join(basedir,"gtpslog.reverse.txt")
        elif kernel == "analysis":
            self.filename = os.path.join(basedir,"gtpslog.analysis.txt")
        else:            
            self.filename = os.path.join(basedir,"gtpslog.%03d.txt" % kernel)
        self.log.basicConfig(level=self.level,filename=self.filename,format="'%(asctime)s %(levelname)s %(message)s'")
        
        
        