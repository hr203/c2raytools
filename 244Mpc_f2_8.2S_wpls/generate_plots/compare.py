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
mpl.rc('xtick', labelsize=30)
mpl.rc('ytick', labelsize=30)
mpl.rc('font',family='serif')
fontsize=30
numberfontsize=24
tickwidth=1.5


sys.path.append('../')
import setup_dirs
import IO

sys.path.append('../../src/')
import c2raytools as c2t

results_dir='/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_wstars'


#redshifts_wstars = setup_dirs.read_redshifts('wstars')
#redshifts_wpls = setup_dirs.read_redshifts('wpls')
#redshifts_wqs = setup_dirs.read_redshifts('wquasars')
start=0

flag1='wstars'
flag2='wquasars'

redshifts1=setup_dirs.read_redshifts(flag1)
redshifts2=setup_dirs.read_redshifts(flag2)

if flag1=='wstars'or flag2=='wstars':
    label1="Stellar Only"        
if flag1=='wpls'or flag2=='wpls':
    if flag1 =='wstars' or flag2=='wstars':
        label2="X-Ray and Stellar"
    else:
        label1="X-Ray and Stellar"
if flag1=='wquasars'or flag2=='wquasars':
    label2="Quasar and Stellar"

def compare_mean_temp():
    
    minmax = "min"

    data1 = plot_1D.mean('temp','data_'+flag1+'/',redshifts1)
    data2 = plot_1D.mean('temp','data_'+flag2+'/',redshifts2)

    if minmax=="max":
        length = max(len(data1),len(data2))

        means = np.zeros(2*length).reshape(length,2)
        for i in range(len(data1)):
            means[i,0] = data1[i]
        for i in range(len(data2)):
            means[i,1] = data2[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
        plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_'+flag1+'_'+flag2+'/',redshifts1,redshifts1,label1,label2)

    elif minmax=="min":
        length = min(len(data1),len(data2))
        means = np.zeros(2*length).reshape(length,2)
        for i in range(length):
            means[i,0] = data1[i]
        for i in range(length):
            means[i,1] = data2[i]
        print means[i,1]
        print means[i,0]
    #prin redshifts_wstars[i], redshifts_wpls[i]
        plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_'+flag1+'_'+flag2+'/',redshifts2,redshifts2,label1,label2)


def compare_rms():

    length=min(len(redshifts1),len(redshifts2))
    if length!=0:
        length=length-1
    rmsdata=np.zeros(length).reshape((length)/2,2)
    redshifts_short=np.zeros((length)/2)
    for i in range(0,length,2):
        redshifts_short[i/2]=redshifts2[i] 
        print "opening ../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin and ../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin" 
        data1=np.load("../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin")
        data2=np.load("../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin")
        mean1=np.mean(data1)
        mean2=np.mean(data2)
      #  length=len(data)
        sum=0.0
        print  "summing..."
        for x in range(length):
            for y in range(length):
                for z in range(length):
                    rmsdata[i/2,0]=rmsdata[i/2,0]+(data1[x,y,z]-mean1)**2
                    rmsdata[i/2,1]=rmsdata[i/2,1]+(data2[x,y,z]-mean2)**2
        rmsdata[i/2,:]=np.sqrt(rmsdata[i/2,:]/(length**3)/2.0)
#        rmsdata[:,i]=np.sqrt(rmsdata[:,i]/(length**3)/2.0)
    plot_1D.plot_mean(rmsdata,"rms_100b","rms (mK)","compare_"+flag1+"_"+flag2+"/",redshifts_short,redshifts_short,label1,label2)


#    for i in range(length): 
#        fr_wstars=np.loadtxt("../generate_data/data_wstars/powerSpectraFrequencies_dbt_"+str('%.3f' % redshifts[i]))
#        data_wstars = np.loadtxt("../generate_data/data_wstars/powerSpectra_"+str('%.3f' % redshifts[i]))*1.0e-3
#        rmsdata[0,i]=sci.simps(data_wstars,fr_wstars)/1.0e-3
#
#        fr_wpls=np.loadtxt("../generate_data/data_wstars/powerSpectraFrequencies_dbt_"+str('%.3f' % redshifts[i]))
#        data_wpls = np.loadtxt("../generate_data/data_wstars/powerSpectra_"+str('%.3f' % redshifts[i]))*1.0e-3
#        rmsdata[1,i]=sci.simps(data_pls,fr_wpls)/1.0e-3
#
#    plot_1D.plot_mean(rmsdata,"rms","rms (mK)",'compare_wstars_wpls/',redshifts_wstars,redshifts_wstars,"Stellar only","Stellar and X-ray")
#
#    wstars_temp = plot_1D.mean('temp','data_wstars/',redshifts_wstars)
#    wpls_temp = plot_1D.mean('temp','data_wpls/',redshifts_wpls)
#
#    length = max(len(wstars_temp),len(wpls_temp))
#    print length
#    means = np.zeros(2*length).reshape(length,2)
#    for i in range(length):
#        means[i,0] = wstars_temp[i]
#    for i in range(length):
#        means[i,1] = wpls_temp[i]
#        #print means[i,1]
#    #print redshifts_wstars[i], redshifts_wpls[i]
#    plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_wstars_wpls/',redshifts_wstars,redshifts_wstars,"Stellar only","Stellar and X-ray")

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
    plot_1D.plot_mean(means,"Test","$(T_S - T_{CMB}) / T_S$",'compare_wstars_wpls/',redshifts_wstars,redshifts_wpls,"Stellar only","Stellar and X-ray")



def compare_mean_dbt():
    wstars_temp = plot_1D.mean('dbt','data_'+flag1+'/',redshifts1)
    wpls_temp = plot_1D.mean('dbt','data_'+flag2+'/',redshifts2)
    length = max(len(wstars_temp),len(wpls_temp))/2
    if (length%2!=0):
        length=length-1
    means = np.zeros(2*length).reshape(length,2)
    redshifts3 = np.zeros(length)
    redshifts4 = np.zeros(len(wpls_temp)/2)
    for i in range(length*2):
        if i%2==0:
            means[i/2,0] = wstars_temp[i]
            if (i<len(wpls_temp)):
                means[i/2,1] = wpls_temp[i]
            redshifts3[i/2]=redshifts1[i]
    for i in range(length): 
        redshifts4[i/2]=redshifts2[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"dbt","$\delta T_b$",'compare_'+flag1+'_'+flag2+'/',redshifts3,redshifts4,label1,label2)



def compare_mean_xfrac():
    wstars_xfrac = plot_1D.mean('xfrac','data_'+flag1+'/',redshifts1)
    wpls_xfrac = plot_1D.mean('xfrac','data_'+flag2+'/',redshifts2)
    wstars_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag1+'/',redshifts1)
    wpls_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag2+'/',redshifts2)
    wstars_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag1+'/',redshifts1)
    wpls_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag2+'/',redshifts2)

    length = max(len(wstars_xfrac),len(wpls_xfrac))

    means = np.zeros(6*length).reshape(length,6)
    for i in range(length):
        if i<len(wpls_xfrac):
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

    plot_1D.plot_mean(means,"xfrac","log$_{10}$ (Ionised Fraction)",'compare_'+flag1+'_'+flag2+'/',redshifts1,redshifts2,"Stellar & X-ray binaries","Stellar only")



def mjfFormatter(x,pos):
    conv=244.0/250.0
    x = int(conv*x)
    if x!=100 and x!=199 and x!=0:
        x = ''
    if x==199:
        x=x+1
#    if x==0:
    #    print "x=0"
 #       x=0
    #if x == 1:
    #    print "x=1"
    #    x=0
    return x #"{0:.0f}".format(x)

def compare_hisograms_temperature():
    for i in range(0,len(redshifts2)):
         data = np.zeros(2*250**3).reshape(250**3,2)
         xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag2+'/Temper3D_'+str('%.3f' % redshifts2[i]) + '.bin'
         dbt_wpls = c2t.TemperFile(xfrac_filename).temper
         data[:,0]=dbt_wpls.flatten()
#        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
         xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag1+'/Temper3D_'+str('%.3f' % redshifts1[i]) + '.bin'
         dbt_wstars = c2t.TemperFile(xfrac_filename).temper
         data[:,1] = dbt_wstars.flatten() #need to import other dataset and put in right format
         plot_1D.plot_log_histogram(data,"Temperature, Redshift:" +str('%.3f' % redshifts1[i]),"loghist_temper_"+str(i+10)+'_'+str('%.3f' % redshifts2[i]),"Temperature(K)",5000,'log',"compare_"+flag1+"_"+flag2+"/")


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
    if flag1=='wpls' or flag2=='wpls':
        redshifts = [14.912,14.294,13.733,13.221]
    else:
        redshifts = [14.294,13.733,12.603,11.546]

    fig, axes = plt.subplots(2,len(redshifts),figsize=(13,6))
    #fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    #fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
#    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.00001)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        f1 = open("../generate_data/data_"+flag2+"/dbt/map_dbt_"+str('%.3f' % redshifts[i])+'.bin')
        dbt_wpls = np.load(f1)
        f1.close()
        im = axes[0,i].imshow(dbt_wpls[250/2,:,:],cmap='dbtmap',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize,color='white')
#        dbt_wstars= IO.readbin("dbt_"+str('%.3f' % redshifts[i]),'data_'+flag1+'/')
        f2=open("../generate_data/data_"+flag1+"/dbt/map_dbt_"+str('%.3f' % redshifts[i])+'.bin')
        dbt_wstars= np.load(f2)
        f2.close()
        axes[1,i].imshow(dbt_wstars[250/2,:,:],cmap='dbtmap',vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.01, 0.7])
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
    axes[0,0].set_ylabel("X-Rays & Stellar",fontsize=numberfontsize)
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[1,0].set_ylabel("Stellar",fontsize=numberfontsize)
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,3].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:len(redshifts)],yticks=[])
    print "saveing fig as " + "compare_"+flag1+"_"+flag2+"/map_dbt.png"
    plt.savefig("compare_"+flag1+"_"+flag2+"/map_dbt.png")
    plt.close()

#########################################################################
def compare_xfracmaps(id):
    maxi = 0.0
    mini= -14.00

    if flag1=='wpls' or flag2=='wpls':
        redshifts = [14.912,14.294,13.733,13.221]
    else:
        redshifts = [14.294,13.733,13.221,11.918]

    fig, axes = plt.subplots(2,len(redshifts),figsize=(13,6))

    #fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    #fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
#    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.00001)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        #dbt_wpls = #IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag2+'/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wpls = c2t.XfracFile(xfrac_filename).xi

        im = axes[0,i].imshow(np.log10(dbt_wpls[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,220,"z = "+str(redshifts[i]),fontsize=fontsize-2,color='white')
#        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag1+'/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wstars = c2t.XfracFile(xfrac_filename).xi
        axes[1,i].imshow(np.log10(dbt_wstars[250/2,:,:]),cmap='Blues_r',vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.01, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[-14,-12,-10,-8,-6,-4,-2,0])
    cbar.ax.set_yticklabels
    cbar.set_label("$Ionised Fraction$",size=fontsize)
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
    axes[0,0].set_ylabel(label2,fontsize=numberfontsize)
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[1,0].set_ylabel(label1,fontsize=numberfontsize)
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,3].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:len(redshifts)],yticks=[])

    plt.savefig("compare_"+flag1+"_"+flag2+"/map_xfrac"+id+".png")
    plt.close()
  


###################################################################################
def compare_tempmaps():
    maxi =13000
    mini= 0.00

    tcmap = {'red':   ((0.0, 0.0, 0.0),
                    (0.85, 0.0, 1.0),
                    (1.0, 1.0, 1.0)),

               'green': ((0.0, 0.0, 0.0),
                        (0.85,0.0, 0.0),
                        (1.0, 1.0, 1.0)),

               'blue':  ((0.0, 0.0, 0.0),
                        (0.85, 1.0, 0.0),
                        (1.0, 1.0, 1.0))
            }


#    plt.register_cmap(name='cmapt',data=tcmap)
#    cmapt = plt.get_cmap('cmapt')

    cmapt='PuRd'

    if flag1=='wpls' or flag2=='wpls':
        redshifts = [14.912,14.294,13.733,13.221]
    else:
        redshifts = [14.294,13.733,13.221,11.918]

    fig, axes = plt.subplots(2,len(redshifts),figsize=(13,6))
    #fig.text(0.5,0.02,"Distance (Mpc)",ha='center',va='center',fontsize=fontsize)
    #fig.text(0.06,0.55,"Distance (Mpc)",ha='center',va='center',rotation='vertical',fontsize=fontsize)
    plt.gcf().subplots_adjust(right=0.15)
#    plt.gcf().subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.00001)
    for i in range(len(redshifts)):#start,len(redshifts)):
        print "Doing redshift " + str(redshifts[i])+"..."
        #dbt_wpls = #IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wpls/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag2+'/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wpls = c2t.TemperFile(xfrac_filename).temper

        im = axes[0,i].imshow(dbt_wpls[250/2,:,:],cmap=cmapt,vmin=mini,vmax=maxi,origin='lower')
        axes[0,i].text(20,210,"z = "+str('%.2f'%redshifts[i]),size=fontsize,color='blue')#,bbox=dict(facecolor='white',edgecolor='none'))
#        dbt_wstars= IO.readmap("dbt_"+str('%.3f' % redshifts[i]),'data_wstars/')
        xfrac_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag1+'/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        dbt_wstars = c2t.TemperFile(xfrac_filename).temper
#        print dbt_wstars[125,125,125], dbt_wpls[125,125,125]
        axes[1,i].imshow(dbt_wstars[250/2,:,:],cmap=cmapt,vmin=mini,vmax=maxi,origin='lower')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.01, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax,ticks=[0,100,200,300,400,500,600,700,800,900,1000])
    cbar.ax.set_yticklabels
    cbar.set_label("Temperature (K)",size=fontsize)

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
    axes[0,0].set_ylabel(label2,size=numberfontsize)
    axes[0,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')
    axes[1,0].set_ylabel(label1,size=numberfontsize)
    axes[0,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',bottom='off',left='off')

    axes[1,0].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off')
    axes[1,1].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,2].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')
    axes[1,3].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='out',top='off',right='off',left='off')

    plt.setp(axes[0,:],xticks=[])
    plt.setp(axes[:,1:len(redshifts)],yticks=[])

    plt.savefig("compare_"+flag1+"_"+flag2+"/map_temp.png")
    plt.close()






#compare_dbtmaps()
#compare_xfracmaps('')
#compare_xfracmaps('He1')
#compare_xfracmaps('He2')
compare_tempmaps()

#compare_mean_xfrac()
#compare_mean_dbt()
#compare_mean_temp()
#compare_rms()

#compare_hisograms_temperature()
