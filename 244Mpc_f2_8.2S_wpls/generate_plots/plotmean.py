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

def mean(name):
#   redshifts = rs.read_redshifts("../red_ori2.dat")
#    redshifts = rs.read_redshifts("../../red_ori.dat")
    mean = np.zeros(len(redshifts))
    file = open('../generate_data/data/mean_'+name+'.dat', 'r')
    i=0
    for line in file:
        mean[i] = float(line)
        i=i+1
    return mean

def plot_mean(data,title,ylabel):

#    if type(data.shape==(len(redshifts),)):
    plt.figure()
    if(data.shape==(len(redshifts),)):  
        plt.plot(redshifts,data)
    elif(data.shape==(len(redshifts),2)):
        plt.plot(redshifts,data[:,0])
        plt.plot(redshifts,data[:,1])
    elif (data.shape==(len(redshifts),3)):
        plt.plot(redshifts,data[:,0])
        plt.plot(redshifts,data[:,1])
        plt.plot(redshifts,data[:,2])
    else:
        print "Data is in wrong format"

     #   plt.errorbar(rds,mn,yerr=er,fmt='o')
     #     plt.errorbar(rds,mn2,yerr=er2,fmt='o')
     #   plt.errorbar(rds,mn3,yerr=er3,fmt='o')
    plt.xlim(redshifts[1],redshifts[len(redshifts)-1])
    plt.title(title)
    plt.xlabel("Redshift")
    plt.ylabel(ylabel)
    #plt.legend(loc=2)
    plt.savefig("plots/"+title+".png")
    plt.close()

plot_mean(mean("temp"),"Tempererature","Temperature (K)")

means = np.zeros(len(redshifts)*3).reshape(len(redshifts),3)
means[:,0] = mean('xfrac')
means[:,1] = mean('xfracHe1')
means[:,2] = mean('xfracHe2')
plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")
