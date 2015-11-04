import numpy as np
import plot_1D
import plot_2D
import sys

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
mpl.rc('xtick', labelsize=20)
mpl.rc('ytick', labelsize=20)
mpl.rc('font',family='serif')
fontsize=20
numberfontsize=20
tickwidth=1.5


sys.path.append('../')
import setup_dirs
import IO



redshifts_wstars = setup_dirs.read_redshifts('wstars')
redshifts_wpls = setup_dirs.read_redshifts('wpls')
start=0
def compare_mean_temp():
    wstars_temp = plot_1D.mean('temp','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('temp','data_wpls/',redshifts_wpls)

    length = min(len(wstars_temp),len(wpls_temp))

    means = np.zeros(2*length).reshape(length,2)
    for i in range(length):
        means[i,0] = wstars_temp[i]
    for i in range(length):
        means[i,1] = wpls_temp[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_wstars_wpls/',redshifts_wpls,"Stellar only","Stellar and X-ray")

def compare_mean_dbt():
    wstars_temp = plot_1D.mean('dbt','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('dbt','data_wpls/',redshifts_wpls)

    length = min(len(wstars_temp),len(wpls_temp))/2
    means = np.zeros(2*length).reshape(length,2)
    redshifts2 = np.zeros(length)
    for i in range(length*2):
        if i%2==0:
            means[i/2,0] = wstars_temp[i]
            means[i/2,1] = wpls_temp[i]
            redshifts2[i/2]=redshifts_wpls[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"dbt","$\delta T_b$",'compare_wstars_wpls/',redshifts2,"Stellar only","Stellar and X-ray")



def compare_mean_xfrac():
    wstars_xfrac = plot_1D.mean('xfrac','data_wstars/',redshifts_wstars)
    wpls_xfrac = plot_1D.mean('xfrac','data_wpls/',redshifts_wpls)
    wstars_xfracHe1 = plot_1D.mean('xfracHe1','data_wstars/',redshifts_wstars)
    wpls_xfracHe1 = plot_1D.mean('xfracHe1','data_wpls/',redshifts_wpls)
    wstars_xfracHe2 = plot_1D.mean('xfracHe2','data_wstars/',redshifts_wstars)
    wpls_xfracHe2 = plot_1D.mean('xfracHe2','data_wpls/',redshifts_wpls)

    length = min(len(wstars_xfrac),len(wpls_xfrac))

    means = np.zeros(6*length).reshape(length,6)
    for i in range(length):
        means[i,0] = wpls_xfrac[i]
        means[i,1] = wpls_xfracHe1[i]
        means[i,2] = wpls_xfracHe2[i]
        means[i,3] = wstars_xfrac[i]
        means[i,4] = wstars_xfracHe1[i]
        means[i,5] = wstars_xfracHe2[i]

    means[0,2] = means[1,2]
    means[0,5] = means[1,5]  
   #print means[:,2], means[:,5]
    #print redshifts_wstars[i], redshifts_wpls[i]

    plot_1D.plot_mean(means,"xfrac","log$_{10}$ (Ionised Fraction)",'compare_wstars_wpls/',redshifts_wpls,"Stellar and X-ray binaries","Stellar only")


def mjfFormatter(x,pos):
    conv=244.0/250.0
    print x
    x = int(conv*x)
    if x!=100 and x!=200:
        x = ''
    print x
    return x #"{0:.0f}".format(x)

def compare_hisograms_temperature():
    for i in range(12,len(redshifts)):
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(temp_filename).temper
        data = data.flatten() #need to import other dataset and put in right format
#        plot.1D.plot_log_histogram(data,"Temperature, Redshift:" +str('%.3f' % redshifts[i]),"loghist_temper_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)")


def compare_dbtmaps():
    mini = -300
    maxi=300
    redshifts = [17.525,15.360,13.733]
#    mainfig = plt.figure()
#    mainfig.text(0.45,0.04,"Distance (Mpc)",ha='center',va='center')
#    mainfig.text(0.04,0.5,"Distance (Mpc)",ha='center',va='center',rotation='vertical')

    fig, axes = plt.subplots(2,3,figsize=(12,8))
    fig.text(0.5,0.04,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    fig.text(0.03,0.5,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)

#    plt.setp(axes[0,:],xticks=[])
#    plt.setp(axes[:,1:3],yticks=[])
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.02,hspace=0.02)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        dbt_wpls = IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')

        
       # ax=plt.subplot(gs1[i])
       # plt.axis('on')
        im = axes[0,i].imshow(dbt_wpls[250/2,:,:],cmap='gnuplot2',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
       # ax.set_xticklabels([])
       # ax.set_yticklabels([])
       # ax.set_aspect('equal')
        #plt.subp

        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
       
       # ax2=plt.subplot(gs1[i+3])
       # plt.axis('on')

        axes[1,i].imshow(dbt_wstars[250/2,:,:],cmap='gnuplot2',vmin=mini,vmax=maxi,origin='lower')
#        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
       # ax2.set_xticklabels([])
       # ax2.set_yticklabels([])
       # ax2.set_aspect('equal')
        #plt.subp

    #fig = plt.figure() 

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.9, 0.12, 0.03, 0.35])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[-300,-150,0,150,300])
#    cbar.ax.set_yticklabels
    cbar.set_label("$\delta T_b (mK)$",size=fontsize)

    cbar_ax = fig.add_axes([0.9,0.52,0.03,0.35])
    cbar = fig.colorbar(im,cax=cbar_ax,ticks=[-300,-150,0,150,300])
    cbar.set_label("$\delta T_b (mK)$",size=fontsize)

    cnv=250.0/244.0
    for i in range(len(axes[:,0])):
        for j in range(len(axes[0,:])):
                axes[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')

                m0 = MultipleLocator(50*cnv)
                axes[i,j].yaxis.set_major_locator(m0)
                axes[i,j].xaxis.set_major_locator(m0)

                axes[i,j].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                axes[i,j].xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))

                m0 = MultipleLocator(10*cnv)
                axes[i,j].yaxis.set_minor_locator(m0)
                axes[i,j].xaxis.set_minor_locator(m0)
        
#                axes[i,j].set_xticklabels([])

#                axes[i,j].tick_locs=[0,50*cnv,100*cnv,150*cnv,200*cnv]
#                axes[i,j].tick_lbls=['0','','100','','200']


    axes[0,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off')
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')


    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:3],yticks=[])

#    plt.setp(axes[:,2],yticks=[])
#    axes[0,1].set_xticks
    

#    f.subplots_adjust(right=0.8)
#    cbar_ax = fig.add_axes([0.85,0.15,0.05,0.7])
#    f.colorbar(im,cax=cbar_ax)
#    cbar = fig=.colorbar(orientation='vertical')
#    cbar.set_label("$\delta T_b$",size=fontsize)

    #plt.tight_layout()
    plt.savefig("compare_wstars_wpls/map_dbt.png")
    plt.close()


compare_dbtmaps()

#compare_mean_xfrac()
#compare_mean_dbt()
#compare_mean_temp()
