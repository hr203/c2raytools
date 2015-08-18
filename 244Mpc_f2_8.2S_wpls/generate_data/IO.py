import numpy as np
import pylab as pl
#import redshifts as rs
import sys
sys.path.append('../../src/')
import c2raytools as c2t
import matplotlib.pyplot as plt

mesh=250

def writedata(data,filename):
    file = open(filename,'w')
    for it in range(len(data)):
        file.write(str(data[it])+'\n')
    print "Written"

def readmap(name):
    map = np.zeros(mesh**3).reshape(mesh,mesh,mesh)
    file = open('data/map_'+name+'.dat', 'r')
    c=0
    for line in file:
        x = int((c)/(mesh**2))
        y = int(c%(mesh)/2.0) 
        z = c%3
        map[x,y,z] = float(line)
        c=c+1
    return map

def writemap(data,filename):
    file = open(filename, 'w')
    l = len(data[:,1,1])
    for i in range(l):
        for j in range(l):
            for k in range(l):
                file.write(str(data[i,j,k])+'\n')
    file.close()

