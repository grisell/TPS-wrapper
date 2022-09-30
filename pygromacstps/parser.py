'''
Created on Mar 2, 2010

@author: wolf
'''


class gdistparser(object):
    def __init__(self):
        self.data = []
        
        
    def readDist(self,filename):
        index = 0
        self.data = []
        for line in open(filename,'r'):
            if not line[0] in ['@','#']:
                wlist = [index]
                wlist.extend(line.split())
                self.data.append(wlist)
                index += 1
        
        
        
        