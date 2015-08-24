import numpy as np
import pylab as pl
import sys
sys.path.append('../')
import redshifts as rs
#sys.path.append('../')
import IO
sys.path.append('../../src/')
import c2raytools as c2t
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

redshifts = rs.read_redshifts("../red_ori2.dat")
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

def plot(dataslice,it,title,name,mini,maxi,cmap='hot',type='lin'):
    plt.figure()
    print mini,maxi
    plt.title(str(title) + ", Redshift: " + str(rs))
    if (type=='lin'):
        plt.imshow(dataslice,cmap=cmap,vmin=mini,vmax=maxi,origin='lower')
    else:
        plt.imshow(dataslice,cmap=cmap,norm=LogNorm(),vmin=mini,vmax=maxi,origin='lower')
    plt.colorbar(orientation='vertical')
    plt.title(title+", Redshift: "+str(redshifts[it]))
    plt.savefig("plots/"+name+"_"+str(it+10)+"_"+str(redshifts[it])+".png")
    plt.close()

def plottemp():
    for i in range(len(redshifts)-1,-1,-1):
        temperature = IO.readmap("temper_"+str('%.3f' % redshifts[i]))
        mini=0.0
        if (i == (len(redshifts)-1)):
            maxi = findmax(temperature)/2.0
        #print mini,maxi
        #print temperature[mesh/2,:,:]
        plot(temperature[mesh/2,:,:],i,"Temperature","temp",mini,maxi)

def plotdbt():
    dbt = IO.readmap("dbt_"+str('%.3f' % redshifts[len(redshifts)-8]))
    maxi = findmax(dbt)
    dbt = IO.readmap("dbt_"+str('%.3f' % redshifts[0]))
    mini = np.min(dbt)
    for i in range(len(redshifts)):
        dbt = IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
        plot(dbt[mesh/2,:,:],i,"Differential Brightness Temperature","dbt",'gnuplot2')
    print "Complete"    

def plotxfrac(id):
    for i in range(len(redshifts)):
        xfrac=IO.readmap('xfrac'+id+'_'+str('%.3f' % redshifts[i]))
        maxi = 1.0
#        if (i==0):
#            mini=np.min(xfrac)
        mini = 0.000199999994948
        plot(xfrac[mesh/2,:,:],i,id + "Ionised Fraction","xfrac"+id,mini,maxi,cmap='Blues_r',type='log')

#plottemp()
plotxfrac('He1')
plotxfrac('He2')
#plot_mean(mean("temp"),"Tempererature","Temperature (K)")
#print plotdbt()
#means = np.zeros(len(redshifts)*3).reshape(len(redshifts),3)
#means[:,0] = mean('xfrac')
#means[:,1] = mean('xfracHe1')
#means[:,2] = mean('xfracHe2')
##plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")
