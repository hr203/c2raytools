import numpy as np
import pylab as pl
import sys
sys.path.append('../')
import IO
import setup_dirs 

sys.path.append('../../src/')
import c2raytools as c2t

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
mpl.rc('xtick', labelsize=22) 
mpl.rc('ytick', labelsize=22) 
mpl.rc('font',family='serif')
fontsize=24
numberfontsize=22
tickwidth=1.5

redshifts = setup_dirs.read_redshifts()
nosteps=setup_dirs.nosteps()


start=55#49
mesh=250
maxi = 0
mini = 0

def map(name):
    map = np.zeros(mesh**3).reshape(mesh,mesh,mesh)
    file = open('../generate_data/'+setup_dirs.resultsdir()+'/map_'+name+'.dat', 'r')
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

def plot_scatter_contour(x,y,name,ylabel,xlabel,it,nbins=200):
    #plt.figure()
    fig=plt.figure(figsize=(9, 7.25), dpi= 300)
    plt.subplot('111', axisbg='black')
    ax = plt.gca()
    #general tick parameters
    ax.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')
    ax.tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')

    #y-axis

    #ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2),useOffset=False)
    m0 = MultipleLocator(5000)
    ax.yaxis.set_major_locator(m0)
    m0 = MultipleLocator(1000)
    ax.yaxis.set_minor_locator(m0)
    #x-axis
    m0 = MultipleLocator(1)
    ax.xaxis.set_major_locator(m0)
    m0 = MultipleLocator(0.2)
    ax.xaxis.set_minor_locator(m0)

    norm = 1.0
    if name =="equation_of_state":
        h =0.7
        omega_m = 0.27
        omega_delta = 0.73
        z =  redshifts[it]
        G = 6.67408e-11
        Mpc_in_metres=3.086e+22

        H0 = h/Mpc_in_metres*1e5
        pcrit = 3.0*H0**2/(8*np.pi*G)#*(omega_m*(0.0+z)**3+omega_delta)
        pcrit = pcrit/1.0e4
        norm = pcrit


    #plot data
    cax =plt.hist2d(x/norm, y, bins=nbins,norm=LogNorm(),normed=True,cmap='Blues_r')
    cbar = plt.colorbar(orientation='vertical')
    #cbar.set_ticks([maxi,maxi/10.0,maxi/100.0,maxi/1000.0,maxi/10000.0,maxi/1000000.0])
    #cbar.set_ticklabels(['$10**0$','$10**{-1}$','$10**{-2}$','$10**{-3}$','$10**{-4}$','$10**{-5}$'])
    plt.ylabel(ylabel,size=fontsize)
    plt.ylim(0,25100)
    plt.xlabel(xlabel,size=fontsize)
    plt.tight_layout()   
    plt.savefig(setup_dirs.plotsdir()+name+"_"+str(it+10)+"_"+str(redshifts[it])+".png")
    plt.close()
    
#    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)
#
#    axes[0, 0].set_title('Scatterplot')
#    axes[0, 0].plot(x, y, 'ko',alpha=0.01,color='orange')
#
#    axes[0, 1].set_title('Hexbin plot')
#    axes[0, 1].hexbin(x, y, gridsize=nbins)
#
#    axes[1, 0].set_title('2D Histogram')
#    axes[1, 0].hist2d(x, y, bins=nbins,norm=LogNorm())
#
# #   axes[1, 1].set_title('Combination')
 #   axes[1, 1].scatter(x,y,alpha=0.1,s=0.1,color='orange')
 #   fig.tight_layout()
#    plt.savefig(setup.plotsdir()+name+"_"+str(it+10)+"_"+str(redshifts[it])+".png")



def plot(dataslice,it,label,name,mini,maxi,cmap='hot',type='lin'):
    plt.figure()
    print mini,maxi
#    plt.title(str(title) + ", Redshift: " + str(redshifts[it]))
    if (type=='lin'):
        plt.imshow(dataslice,cmap=cmap,vmin=mini,vmax=maxi,origin='lower')
    else:
        plt.imshow(dataslice,cmap=cmap,norm=LogNorm(),vmin=mini,vmax=maxi,origin='lower')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label(label,size=19)

    cnv=250.0/244.0
    tick_locs=[0,50*cnv,100*cnv,150*cnv,200*cnv,250*cnv]
    tick_lbls=[0,50,100,150,200,250]
    plt.xlabel("Distance (Mpc/h)",size=18)
    plt.ylabel("Distance (Mpc/h)",size=18)
    #plt.title(title+", Redshift: "+str(redshifts[it]))
    plt.savefig(setup_dirs.plotsdir()+name+"_"+str(it+10)+"_"+str(redshifts[it])+".png")
    plt.close()

def plottemp():
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        temperature = c2t.TemperFile(temp_filename).temper
        mini=0.0
        maxi = 14300.0
#        if (i == (len(redshifts)-25)):
 #           maxi = findmax(temperature)/2.0
        plot(temperature[mesh/2,:,:],i,"Temperature (K)","temp",mini,maxi)

def plotdbt():
    mini = -300
    maxi=300
    print mini,maxi
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        dbt = IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
        plot(dbt[mesh/2,:,:],i,"$\delta T_b$","dbt",mini,maxi,'gnuplot2')
    print "Complete"    

def plotxfrac(id):
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        xfrac_filename = setup_dirs.path()+'xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac = c2t.XfracFile(xfrac_filename).xi
        maxi = 1.0
        mini = 0.000199999994948
        plot(xfrac[mesh/2,:,:],i,id + "Ionised Fraction","xfrac"+id,mini,maxi,cmap='Blues_r',type='log')

def scatter_contour():
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        xfrac_filename = setup_dirs.path()+'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac = c2t.XfracFile(xfrac_filename).xi.flatten()
        print "Read ionized fraction"
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        temperature = c2t.TemperFile(temp_filename).temper.flatten()
        print "Read Temperature"
        print "Generating plots..."
        plot_scatter_contour(xfrac,temperature,"contour_scatter_xfrac_temperature","Temperature (K)","Ionised Fraction",i)

def equation_of_state():
    print redshifts
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        if nosteps == 2 or i<nosteps:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-i%nosteps]) + 'n_all.dat'
            print "density filename: "+ str('%.3f' % redshifts[i-(i)%nosteps]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-i%nosteps-1]) + 'n_all.dat'
            print "density filename: "+ str('%.3f' % redshifts[i-(i)%nosteps-1]) + 'n_all.dat'
        density = c2t.DensityFile(density_filename).cgs_density.flatten()
        print "Read density file"
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        temperature = c2t.TemperFile(temp_filename).temper.flatten()
        print "Read temperature"
        print "Generating plots..."
        plot_scatter_contour(density,temperature,"equation_of_state","Temperature (K)","Density ($\mathrm{\Omega}$)",i)


#scatter()
#plottemp()
#plotxfrac('')
#plotxfrac('He1')
#plotxfrac('He2')
#plotdbt()
#plot_mean(mean("temp"),"Tempererature","Temperature (K)")
#print plotdbt()
#means = np.zeros(len(redshifts)*3).reshape(len(redshifts),3)
#means[:,0] = mean('xfrac')
#means[:,1] = mean('xfracHe1')
#means[:,2] = mean('xfracHe2')
#plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")
