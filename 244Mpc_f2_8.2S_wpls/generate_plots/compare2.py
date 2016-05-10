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
mpl.rc('xtick', labelsize=16)
mpl.rc('ytick', labelsize=16)
mpl.rc('font',family='serif')
fontsize=16
numberfontsize=16
tickwidth=1.5


sys.path.append('../')
import setup_dirs
import IO

sys.path.append('../../src/')
import c2raytools as c2t

results_dir='/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_wstars'


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

#def compare_mean_temp():
#    wstars_temp = plot_1D.mean('temp','data_wstars/',redshifts_wstars)
#    wpls_temp = plot_1D.mean('temp','data_wpls/',redshifts_wpls)
#
#    length = max(len(wstars_temp),len(wpls_temp))
#    
#    tcmb_wstars = 2.725*(np.ones(len(wstars_temp)+redshifts_wstars))
#    tcmb_wpls = 2.725 *(np.ones(len(wpls_temp)+redshifts_wpls))
#
#    wstars_test = (wstars_temp-tcmb_wstars)/wstars_temp
#    wspls_test = (wpls_temp-tcmb_wpls)/wpls_temp
#
#    means = np.zeros(2*length).reshape(length,2)
#    for i in range(length):
#        means[i,0] = wstars_test[i]
#    for i in range(length):
#        means[i,1] = wpls_test[i]
#        #print means[i,1]
#    #print redshifts_wstars[i], redshifts_wpls[i]
#    plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_wstars_wpls/',redshifts_wpls,"Stellar only","Stellar and X-ray")



def compare_mean_test():
    wstars_temp = plot_1D.mean('temp','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('temp','data_wpls/',redshifts_wpls)

    length = max(len(wstars_temp),len(wpls_temp))
    
    wstars_tcmb = 2.725*(np.ones(len(redshifts_wstars))+redshifts_wstars)
    wpls_tcmb = 2.725*(np.ones(len(redshifts_wpls))+redshifts_wpls)

    wstars_test = (wstars_temp-wstars_tcmb)/wstars_temp
    wpls_test = (wpls_temp-wpls_tcmb)/wpls_temp


    means = np.zeros(2*length).reshape(length,2)
    for i in range(length):
        means[i,0] = wstars_test[i]
    for i in range(length):
        if (i<len(wpls_temp)):
            means[i,1] = wpls_test[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"Test","$(T_S - T_{CMB}) / T_S$",'compare_wstars_wpls/',redshifts_wstars,"Stellar only","Stellar and X-ray")



def compare_mean_dbt():
    wstars_temp = plot_1D.mean('dbt','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('dbt','data_wpls/',redshifts_wpls)

    length = max(len(wstars_temp),len(wpls_temp))/2
    means = np.zeros(2*length).reshape(length,2)
    redshifts2 = np.zeros(length)
    for i in range(length*2):
        if i%2==0:
            means[i/2,0] = wstars_temp[i]
            if (i<len(wpls_temp)):
                means[i/2,1] = wpls_temp[i]
            redshifts2[i/2]=redshifts_wstars[i]
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
    x = int(conv*x)
    if x!=100 and x!=199:
        x = ''
    if x==199:
        x=x+1
    #if x==0:
    #    print "x=0"
    #    x = 0
    #if x == 1:
    #    print "x=1"
    #    x=0
    return x #"{0:.0f}".format(x)

def compare_hisograms_temperature():
    for i in range(12,len(redshifts)):
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(temp_filename).temper
        data = data.flatten() #need to import other dataset and put in right format
#        plot.1D.plot_log_histogram(data,"Temperature, Redshift:" +str('%.3f' % redshifts[i]),"loghist_temper_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)")


def compare_dbtmaps():

    cmapdbt = {'red':   ((0.0, 0.0, 0.0),
                    (0.85, 0.0, 1.0),
                    (1.0, 1.0, 1.0)),

               'green': ((0.0, 0.0, 0.0),
                        (0.85,0.0, 0.0),
                        (1.0, 1.0, 1.0)),

               'blue':  ((0.0, 0.0, 0.0),
                        (0.85, 1.0, 0.0),
                        (1.0, 1.0, 1.0))
              }


    plt.register_cmap(name='dbtmap',data=cmapdbt)

    mini = -457
    maxi=83
    redshifts = [14.912,14.294,13.733]

    fig, axes = plt.subplots(2,3,figsize=(12,7))
    fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.02,hspace=0.02)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        dbt_wpls = IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')
        im = axes[0,i].imshow(dbt_wpls[250/2,:,:],cmap='dbtmap',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
        axes[1,i].imshow(dbt_wstars[250/2,:,:],cmap='dbtmap',vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.17, 0.03, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[-450,-400,-350,-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300])
    cbar.ax.set_yticklabels
    cbar.set_label("$\delta T_b \ (mK)$",size=fontsize)

    #cbar_ax = fig.add_axes([0.87,0.54,0.03,0.35])
    #cbar = fig.colorbar(im,cax=cbar_ax,ticks=[-300,-150,0, 150, 300])
    #cbar.ax.set_yticklabels
    #cbar.set_label("$\delta T_b \ (mK)$",size=fontsize)

    cnv=250.0/244.0
    for i in range(len(axes[:,0])):
        for j in range(len(axes[0,:])):
                axes[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')
                axes[i,j].set_xlim([0,250])
                axes[i,j].set_ylim([0,250])

                m0 = MultipleLocator(50.0*cnv)
                axes[i,j].yaxis.set_major_locator(m0)
                axes[i,j].xaxis.set_major_locator(m0)
                axes[i,j].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                axes[i,j].xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                m0 = MultipleLocator(10.0*cnv)
                axes[i,j].yaxis.set_minor_locator(m0)
                axes[i,j].xaxis.set_minor_locator(m0)

    axes[0,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off')
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:3],yticks=[])

    plt.savefig("compare_wstars_wpls/map_dbt.png")
    plt.close()

#########################################################################
def compare_xfracmaps(id):
    maxi = 1.0 
    mini= -5.00
    redshifts = [17.525,15.360,13.733]

    fig, axes = plt.subplots(2,3,figsize=(12,7))
    fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.02,hspace=0.02)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        #dbt_wpls = #IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_w'+'pls/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wpls = c2t.XfracFile(xfrac_filename).xi

        im = axes[0,i].imshow(np.log10(dbt_wpls[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
#        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_w'+'stars/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wstars = c2t.XfracFile(xfrac_filename).xi
        axes[1,i].imshow(np.log10(dbt_wstars[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.33])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[-300,-150,0,150,300])
    cbar.ax.set_yticklabels
    cbar.set_label("$log_{10}$ (xfrac)",size=fontsize)

    cbar_ax = fig.add_axes([0.87,0.54,0.03,0.33])
    cbar = fig.colorbar(im,cax=cbar_ax,ticks=[-300,-150,0, 150, 300])
    cbar.ax.set_yticklabels
    cbar.set_label("$\delta T_b \ (mK)$",size=fontsize)

    cnv=250.0/244.0
    for i in range(len(axes[:,0])):
        for j in range(len(axes[0,:])):
                axes[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')
                axes[i,j].set_xlim([0,250])
                axes[i,j].set_ylim([0,250])

                m0 = MultipleLocator(50.0*cnv)
                axes[i,j].yaxis.set_major_locator(m0)
                axes[i,j].xaxis.set_major_locator(m0)
                axes[i,j].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                axes[i,j].xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                m0 = MultipleLocator(10.0*cnv)
                axes[i,j].yaxis.set_minor_locator(m0)
                axes[i,j].xaxis.set_minor_locator(m0)

    axes[0,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off')
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:3],yticks=[])

    plt.savefig("compare_wstars_wpls/map_xfrac"+id+".png")
    plt.close()
###################################################################################
def compare_tempmaps():
    maxi =14000.0
    mini= 0.00
    redshifts = [17.525,15.360,13.733]

    fig, axes = plt.subplots(2,3,figsize=(12,7))
    fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.02,hspace=0.02)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        #dbt_wpls = #IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_w'+'pls/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wpls = c2t.TemperFile(xfrac_filename).temper

        im = axes[0,i].imshow(np.log10(dbt_wpls[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
#        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_w'+'stars/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wstars = c2t.TemperFile(xfrac_filename).temper
        axes[1,i].imshow(np.log10(dbt_wstars[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.33])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[0,1000,2000,10000])
    cbar.ax.set_yticklabels
    cbar.set_label("$Temperature \ (K)",size=fontsize)

    cbar_ax = fig.add_axes([0.87,0.54,0.03,0.33])
    cbar = fig.colorbar(im,cax=cbar_ax,ticks=[-300,-150,0, 150, 300])
    cbar.ax.set_yticklabels
    cbar.set_label("$Temperature \ (K)$",size=fontsize)

    cnv=250.0/244.0
    for i in range(len(axes[:,0])):
        for j in range(len(axes[0,:])):
                axes[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='out',pad=14.0,top='off',right='off')
                axes[i,j].set_xlim([0,250])
                axes[i,j].set_ylim([0,250])

                m0 = MultipleLocator(50.0*cnv)
                axes[i,j].yaxis.set_major_locator(m0)
                axes[i,j].xaxis.set_major_locator(m0)
                axes[i,j].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                axes[i,j].xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjfFormatter))
                m0 = MultipleLocator(10.0*cnv)
                axes[i,j].yaxis.set_minor_locator(m0)
                axes[i,j].xaxis.set_minor_locator(m0)

    axes[0,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off')
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:3],yticks=[])

    plt.savefig("compare_wstars_wpls/map_temp.png")
    plt.close()




#compare_dbtmaps()
#compare_xfracmaps('')
#compare_xfracmaps('He1')
#compare_xfracmaps('He2')
#compare_tempmaps()

#compare_mean_xfrac()
#compare_mean_dbt()
#compare_mean_temp()
compare_mean_test()
