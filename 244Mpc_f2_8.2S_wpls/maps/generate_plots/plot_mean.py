import numpy as np
import pylab as pl
#import redshifts as rs
import sys
sys.path.append('../../')
import redshifts as rs
sys.path.append('../../../src/')
import c2raytools as c2t
import matplotlib.pyplot as plt

redshifts = rs.read_redshifts("../../red_ori2.dat")
mesh=250
maxi = 0
mini = 0

def map(name):
#   redshifts = rs.read_redshifts("../red_ori2.dat")
#    redshifts = rs.read_redshifts("../../red_ori.dat")
    map = np.zeros(mesh**3).reshape(mesh,mesh,mesh)
    file = open('../generate_data/data/map_'+name+'.dat', 'r')
    c=0
    for line in file:
        x = int((c)/(mesh**2))
        y = int(c%(mesh)/2.0) 
        z = c%3
        map[x,y,z] = float(line)
        c=c+1
    return map

def findmax(data):
    return np.max(data)

def plot(dataslice,it,title,cmap='hot'):
    plt.figure()
    plt.title(str(title) + ", Redshift: " + str(rs))
    plt.imshow(dataslice,cmap=cmap,vmin=mini,vmax=maxi,origin='lower')
    plt.colorbar(orientation='vertical')
    plt.title(title+", Redshift: "+str(redshifts[it]))
    plt.savefig("plots/"+title+"_"+str(it)+"_"+str(redshifts[it])+".png")
    plt.close()

for i in range(len(redshifts)-1,-1,-1):
    temperature = map("temper_"+str('%.3f' % redshifts[i]))
    if (i == (len(redshifts)-1)):
        maxi = findmax(temperature)
    plot(temperature[mesh/2,:,:],i,"Temperature")
    

#plot_mean(mean("temp"),"Tempererature","Temperature (K)")

#means = np.zeros(len(redshifts)*3).reshape(len(redshifts),3)
#means[:,0] = mean('xfrac')
#means[:,1] = mean('xfracHe1')
#means[:,2] = mean('xfracHe2')
##plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")
