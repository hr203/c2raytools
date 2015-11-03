import numpy as np
import pylab as pl
import setup_dirs 
#import redshifts as rs
import sys

sys.path.append('../../src/')
import c2raytools as c2t

#sys.path.append('/src/')
#import c2raytools as c2t

import matplotlib.pyplot as plt

mesh=250

def writedata(data,filename):
    file = open(filename,'w')
    for it in range(len(data)):
        if str(data[it])!='nan':
            file.write(str(data[it])+'\n')
        else:
            print "deleted nan bin"
    print "Written to " + str(filename)


def write2data(data,data2,filename,filename2):
    file = open(filename,'w')
    file2 = open(filename2,'w')
    for it in range(len(data)):
        if str(data[it])!='nan':
            file.write(str(data[it])+'\n')
            file2.write(str(data2[it])+'\n')
        else:
            print "deleted nan bin"
    print "Written to " + str(filename)


def readmap(name,resultdir=setup_dirs.resultsdir()):
    map = np.zeros(mesh**3).reshape(mesh,mesh,mesh)
    #if 
    file = open('../generate_data/'+resultdir+'map_'+name+'.dat', 'r')
    c=0
    for line in file:
        x = int((c)/(mesh**2))
        y = int(c/mesh)%(mesh)
       # y = int(c%(mesh/2.0)) 
        z = c%mesh
        map[x,y,z] = float(line)
        c=c+1
        #if (True):
        #    print x+1,y+1,z+1
    #print map
    print 'Read map from ../generate_data/'+setup_dirs.resultsdir()+'map_'+name+'.dat'
    return map

def readoned(name):
    file = open('../generate_data/'+setup_dirs.resultsdir()+name+'.dat','r')
    cnt=0
    for line in file:
        cnt = cnt+1
    data=np.zeros(cnt)
    file.close()
    file = open('../generate_data/'+setup_dirs.resultsdir()+name+'.dat','r')
    c=-1
    for line in file:
        if (c!=-1):
            data[c] = float(line)
        c=c+1
    file.close()
    print 'read data from ../generate_data/'+setup_dirs.resultsdir()+name+'.dat'
    return data


def writemap(data,filename):
#    print data
    file = open(filename, 'w')
    l = len(data[:,1,1])
    for i in range(l):
        for j in range(l):
            for k in range(l):
                file.write(str(data[i,j,k])+'\n')
    print "Plotted data from " + filename
    file.close()

