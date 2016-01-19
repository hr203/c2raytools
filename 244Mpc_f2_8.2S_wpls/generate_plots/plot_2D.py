from matplotlib.font_manager import FontProperties
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
mpl.rc('xtick', labelsize=25)#26) 
mpl.rc('ytick', labelsize=25)#26) 
mpl.rc('font',family='serif')
fontsize=26
numberfontsize=25
tickwidth=1.5

redshifts = setup_dirs.read_redshifts()
nosteps=setup_dirs.nosteps()

ss=121
start=0
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
    fig=plt.figure(figsize=(9, 10), dpi= 300)
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

def plot_lightcone(dataslice,it,label,name,mini,maxi,redshifts,cmap='hot',type='lin'):
    if name=='dbt' or name=='dbt_lightcone':
        print "using dbt colourmap"
        cmapdbt = {'red':   ((0.0, 0.0, 0.0),
                            (0.846, 0.0, 1.0),
                            (1.0, 1.0, 1.0)),

                   'green': ((0.0, 0.0, 0.0),
                            (0.846,0.0, 0.0),
                            (1.0, 1.0, 1.0)),

                   'blue':  ((0.0, 0.0, 0.0),
                            (0.846, 1.0, 0.0),
                            (1.0, 1.0, 1.0))
                  }
        plt.register_cmap(name='dbtmap',data=cmapdbt)
        cmap = plt.get_cmap('dbtmap')

    fig,axes =plt.subplots(1,1,figsize=(26,6))
    if (type=='lin'):
        im = plt.imshow(dataslice,cmap=cmap,vmin=mini,vmax=maxi,origin='lower')
    else:
        im = plt.imshow(dataslice,cmap=cmap,norm=LogNorm(),vmin=mini,vmax=maxi,origin='lower')
    font = FontProperties()
    font.set_weight('bold')

    cnv=250.0/244.0
    ytick_locs=[0,50*cnv,100*cnv,150*cnv,200*cnv,250*cnv]
    ytick_lbls=[0,50,100,150,200,250]
    
    axes.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')
    m0=MultipleLocator(200)
    axes.xaxis.set_major_locator(m0)

    m0=MultipleLocator(100)
    axes.xaxis.set_minor_locator(m0)

    labels = np.zeros(len(redshifts)/200+2)
    print len(labels-1)
    for i in range(len(dataslice[1,:])):
        if (i+1)%200==0: 
            if i<len(redshifts) and i/100<len(labels):
                labels[i/100] = str('%.2f'%redshifts[i])
    labels[len(labels)-2]=str('%.2f'%redshifts[len(redshifts)-2])
    print labels
    axes.set_xticklabels(labels[::-1])
    plt.xlabel("z",size=fontsize)
    plt.ylabel("Distance (Mpc)",size=fontsize) 
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.25, 0.02, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[-400,-300,-200,-100,0,100,200,300])
    cbar.ax.set_yticklabels
    cbar.set_label("$\delta T_b \ (mK)$",size=fontsize)

#    plt.tight_layout()
    print "saving as..." + setup_dirs.plotsdir()+name+".png"
    plt.savefig(setup_dirs.plotsdir()+name+".png")
    plt.close()



def plot(dataslice,it,label,name,mini,maxi,cmap='hot',type='lin'):
    
    if name=='dbt' or name=='dbt_lightcone_':
        print "using dbt colourmap"
        cmapdbt = {'red':   ((0.0, 0.0, 0.0),
                            (0.846, 0.0, 1.0),
                            (1.0, 1.0, 1.0)),
                   
                   #'green':   ((0.1, 0.1, 0.1),
                   #         (0.846, 0.1, 0.0),
                   #         (0.0, 0.0, 0.0)),
    
                   'green': ((0.0, 0.0, 0.0),
                            (0.846,0.0, 0.0),
                            (1.0, 1.0, 1.0)),

                   'blue':  ((0.0, 0.0, 0.0),
                            (0.846, 1.0, 0.0),
                            (1.0, 1.0, 1.0))

                   #'white': ((0.0, 0.0, 0.0),
                   #         (0.846, 1.0, 0.0),
                   #         (1.0, 0.1, 1.0)),

                  }


        plt.register_cmap(name='dbtmap',data=cmapdbt)
        cmap = plt.get_cmap('dbtmap')

    #plt.figure()
    plt.figure(figsize=(9,7.25),dpi=300)
    #plt.subplots_adjust(right=0.3)
    #plt.subplots_adjust(bottom=1.3)

    print mini,maxi
#    plt.title(str(title) + ", Redshift: " + str(redshifts[it]))
    if (type=='lin'):
        im = plt.imshow(dataslice,cmap=cmap,vmin=mini,vmax=maxi,origin='lower')
    else:
        im = plt.imshow(dataslice,cmap=cmap,norm=LogNorm(),vmin=mini,vmax=maxi,origin='lower')
    font = FontProperties()
    font.set_weight('bold')
    if name == "xfrac" : 
        plt.text(20,220,"HII",fontsize=fontsize,color='white')
    if name == "xfracHe1" : 
        plt.text(20,220,"HeII",fontsize=fontsize,color='white')
    if name == "xfracHe2" : 
        plt.text(20,220,"HeIII",fontsize=fontsize,color='white')
        plt.text(150,20,"z = "+str(redshifts[it]),fontsize=fontsize,color='red')
    if name == "Temperature":
        plt.text(20,220, "Temperature",fontsize=fontsize,color='white')
    
    if name!='dbt':
        print "not dbt..."
        cbar = plt.colorbar(im,orientation='vertical',fraction=0.046,pad=0.04)
        cbar.set_label(label,size=fontsize)
    #cbar.ax.set_yticklabels
        if name!="Temperature" :
            cbar.set_ticks([1.,1./10.0,1./100.0,1./1000.0,1./10000.0,1./1000000.0])
            cbar.set_ticklabels(['$ \ 10 ^0 $  ','$ \ 10  ^{-1}  $  ','$ \ 10  ^{-2}  $  ','$ \ 10  ^{-3}  $  ','$ \ 10  ^{-4}  $  ','$ \ 10 ^{-5}  $  '])
        else:
            cbar.set_ticks([2000,4000,6000,8000,10000,12000,14000])
            cbar.set_ticklabels(['  2  ','  4  ','  6  ','  8  ','  10  $ \ $','  12  ' , '  14   '])

    cnv=250.0/244.0
    tick_locs=[0,50*cnv,100*cnv,150*cnv,200*cnv,250*cnv]
    tick_lbls=[0,50,100,150,200,250]
    #plt.xlabel("Distance (Mpc/h)",size=18)
    #plt.ylabel("Distance (Mpc/h)",size=18)
    #plt.title(title+", Redshift: "+str(redshifts[it]))i

    #plt.tight_layout()
#    if name == "Temperature":
    #plt.gcf().subplots_adjust(left=0.5)
    print "saving as..." + setup_dirs.plotsdir()+name+"_"+str(it+10)+"_"+str(redshifts[it])+".png"
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
        plot(temperature[len(temperature[:,0,0])/2,:,:],i,"Temperature (k K)","Temperature",mini,maxi)

def plotdbt():
   
    mini = -457
    maxi=83
    print mini,maxi
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        dbt = np.load("../generate_data/"+setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin")
#IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
        print np.min(dbt)
        plot(dbt[len(dbt[:,0,0])/2,:,:],i,"$\delta T_b$","dbt",mini,maxi,'seismic')
    print "Complete"    

def plotdbt_lightcone():

    mini = -457
    maxi=83
    print mini,maxi
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        dbt = IO.readmap("dbt_lightcone_"+str('%.3f' % redshifts[i]))
        print np.min(dbt)
        plot(dbt[len(dbt[:,1,1])/2,:,:],i,"$\delta T_b$","dbt_lightcone",mini,maxi,'seismic')
    print "Complete"


def plotxfrac(id):
    for i in range(start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        xfrac_filename = setup_dirs.path()+'xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac = c2t.XfracFile(xfrac_filename).xi
        maxi = 1.0
        mini = 0.000199999994948
        plot(xfrac[len(xfrac[:,1,1])/2,:,:],i,"Ionised Fraction","xfrac"+id,mini,maxi,cmap='Blues_r',type='log')

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
        #density = np.ones(ss**3).reshape(ss,ss,ss)*1.981e-10*(1+redshifts[i])**3
        print "Read density file"
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        temperature = c2t.TemperFile(temp_filename).temper.flatten()
        print "Read temperature"
        print "Generating plots..."
        plot_scatter_contour(density,temperature,"equation_of_state","Temperature (K)","Density ($\mathrm{\Omega}$)",i)

def lightcone():
    filenames = ["" for x in range(len(redshifts))]
#    filenames[0] = len(redshifts)
    mini = -457
    maxi=83 
    for i in range(len(filenames)):
        #filenames[i] = setup_dirs.path() + 'lightcone_temp/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        filenames[i-20] = '../generate_data/'+setup_dirs.resultsdir()+'dbt/map_dbt_'+str('%.3f' % redshifts[i]) + '.bin'
    im,z=c2t.make_lightcone(filenames,redshifts[len(redshifts)-1],redshifts[0]) 
    im=np.asarray(im)
    print(im[125,:,::-1])
    print(im[1,1,1])
    plot_lightcone(im[125,:,::-1],i,"$\delta$ T","dbt_lightcone",mini,maxi,z,cmap='Blues_r')
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
