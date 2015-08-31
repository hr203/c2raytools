import numpy as np
import pylab as pl
import sys
#sys.path.append('../../../src/')
#import c2raytools as c2t
import matplotlib.pyplot as plt
sys.path.append('../')
import setup

redshifts = setup.read_redshifts()

def mean(name):
#   redshifts = rs.read_redshifts("../red_ori2.dat")
#    redshifts = rs.read_redshifts("../../red_ori.dat")
    mean = np.zeros(len(redshifts))
    file = open('../generate_data/'+setup.resultsdir()+'/mean_'+name+'.dat', 'r')
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
    plt.savefig(setup.plotsdir()+title+".png")
    plt.close()


def plot_powerspectra(fr,data,title):
    plt.figure()
    #if(data.shape==(len(redshifts),)):
    print "plotting"#fr[1],data[1]
    plt.plot(fr,data)
    #elif(data.shape==(len(redshifts),2)):
    #    plt.plot(fr,data[:,0])
    #    plt.plot(fr,data[:,1])
    #elif (data.shape==(len(redshifts),3)):
    ##    plt.plot(fr,data[:,0])
    #    plt.plot(fr,data[:,1])
    #    plt.plot(fr,data[:,2])
    #else:
    #    print "Data is in wrong format"
    plt.xlim(0,1)
    plt.ylim(0,1e6)
    plt.title(title)
    plt.xlabel("k")
    plt.ylabel("P(k)")
    plt.savefig(setup.plotsdir()+title+".png")
    plt.close()

def plot_histogram(data,title,name,xlabel,nobins=5000,type='log'):
    fig = plt.figure()

#method 1
    ax = fig.add_subplot(111)
    ax.set_yscale('log',basey=10)
    weights = np.ones_like(data)/float(len(data))
    ax.hist(np.log10(data),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='red',edgecolor='red')

#method 2 
#    results, edges = np.histogram(data,range=(0,10000),bins=100000,normed=1)#,range=[0,1000])
#    binWidth = edges[1] - edges[0]
#    plt.bar(edges[:-1],results*binWidth, binWidth,color='red',edgecolor='red')
#    if (type=='log'):
#        plt.xscale('log')
#    plt.title(title)
    if type == 'log':
        plt.xlabel('log10('+xlabel+')')
    else:
        plt.xlabel(xlabel)
    plt.ylabel("Probability")
    plt.savefig(setup.plotsdir()+name+".png")
    plt.close()
    print "histogram complete"

def plot_histogram3(data,data2,data3,title,name,xlabel,nobins=5000,type='log'):
    fig = plt.figure()
    #print np.amin(data3)
    ax = fig.add_subplot(111)
    ax.set_yscale('log',basey=10)
    if (type=='log'):
        weights = np.ones_like(data)/len(data2)
        ax.hist(np.log10(data),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',alpha=0.4,color='red',edgecolor='none',label='HII')
        ax.hist(np.log10(data2),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',alpha=0.4,color='yellow',edgecolor='none',label='HeII')
        ax.hist(np.log10(data3),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',alpha=0.4,color='blue',edgecolor='none',label='HeIII')
        

    else:
        ax.hist(data,nobins,color='red',alpha=0.01)
        ax.hist(data,nobins,color='green',alpha=0.01)
    plt.title(title)
    plt.legend()
    if type=='log':
        plt.xlabel('log10('+xlabel+')')
        plt.yscale=('log')
    else:
        plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.savefig(setup.plotsdir()+name+".png")
    plt.close()
    print "histogram complete"


#means = np.ones(len(redshifts)*3).reshape(len(redshifts),3)
#means[:,0] = mean('xfrac')
#means[:,1] = mean('xfracHe1')
#means[:,2] = mean('xfracHe2')
#print means
#plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")
#
#mean = mean('temp')
#plot_mean(mean,'temp', 'Temperature (K)')

mean = mean('dbt')
plot_mean(mean,'dbt','Differential Brightness Temperature')
