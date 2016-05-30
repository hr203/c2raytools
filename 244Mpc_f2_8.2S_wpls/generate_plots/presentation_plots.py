import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
#import compare
import sys
import plot_1D
sys.path.append('../')
import setup_dirs
sys.path.append('../../src')
import c2raytools as c2t
c2t.set_sim_constants(244)
plt.rc('xtick', labelsize=26)#30)
plt.rc('ytick', labelsize=26)#30)
plt.rc('font',family='serif')
fontsize=26#30
numberfontsize=24#24#30
tickwidth=1.5

size=26#30
lw = 2.5
lc_zs=np.loadtxt('../generate_data/data_wstars/smooth_zs_.dat')

flag1="wstars"#"wpls"
flag2="wpls"#"wquasars"#"wpls"#quasars"
redshifts1=setup_dirs.read_redshifts(flag1)
redshifts2=setup_dirs.read_redshifts(flag2)
if flag1=='wstars'or flag2=='wstars':
    label1="Stellar"
if flag1=='wpls'or flag2=='wpls':
    if flag1 =='wstars' or flag2=='wstars':
        label2="X-ray & Stellar"
    else:
        label1="X-ray & Stellar"

    if flag2=="wquasars_wpls":
        label2="X-rays & Quasar"
if flag1=='wquasars'or flag2=='wquasars':
    label2="Quasars &  Stellar"

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
#    subax.xaxis.set_tick_params(labelsize=x_labelsize)
#    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def plot1():
    print "Comparing: " + flag1 + " and "+ flag2    
    minmax = "min"
    #load temperature means
    #data1 = plot_1D.mean('temp','data_'+flag1+'/',redshifts1)
    #data2 = plot_1D.mean('temp','data_'+flag2+'/',redshifts2)
    data1 = np.loadtxt('../generate_data/data_'+flag1+'/mean_temp.dat')
    data2 = np.loadtxt('../generate_data/data_'+flag2+'/mean_temp.dat')
    data3 = np.loadtxt('../generate_data/data_'+flag1+'/median_temp.dat')
    data4 = np.loadtxt('../generate_data/data_'+flag2+'/median_temp.dat')
    
#    print data1

    #data3 = plot_1D.median('temp','data_'+flag1+'/',redshifts1)
    #data3 = plot_1D.median('temp','data_'+flag2+'/',redshifts2)
   
    #load xfrac means
    wstars_xfrac = plot_1D.mean('xfrac','data_'+flag1+'/',redshifts1)
    wpls_xfrac = plot_1D.mean('xfrac','data_'+flag2+'/',redshifts2)
    wstars_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag1+'/',redshifts1)
    wpls_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag2+'/',redshifts2)
    wstars_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag1+'/',redshifts1)
    wpls_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag2+'/',redshifts2)

    length = min(len(data1),len(data2))
    if minmax=="max":
        length = min(len(data1),len(data2))
    #print redshifts2 
    meantemp = np.zeros(4*length).reshape(length,4)
    meanxfrac = np.zeros(6*length).reshape(length,6)
    if minmax=="max":
        for i in range(len(data1)):
            meantemp[i,0] = data1[i]
            meantemp[i,2] = data3[i]
            meanxfrac[i,3] = wstars_xfrac[i]
            meanxfrac[i,4] = wstars_xfracHe1[i]
            meanxfrac[i,5] = wstars_xfracHe2[i]

        for i in range(len(data2)):
            meanxfrac[i,0] = wpls_xfrac[i]
            meanxfrac[i,1] = wpls_xfracHe1[i]
            meanxfrac[i,2] = wpls_xfracHe2[i]
            meantemp[i,1] = data2[i]
            meantemp[i,3] = data4[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
#        plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_'+flag1+'_'+flag2+'/',redshifts1,redshifts1,label1,label2)

    elif minmax=="min":
        for i in range(length):
            meantemp[i,0] = data1[i]
            meantemp[i,1] = data2[i]
            meantemp[i,2] = data3[i]
            meantemp[i,3] = data4[i]

            meanxfrac[i,0] = wpls_xfrac[i]
            meanxfrac[i,1] = wpls_xfracHe1[i]
            meanxfrac[i,2] = wpls_xfracHe2[i]
            meanxfrac[i,3] = wstars_xfrac[i]
            meanxfrac[i,4] = wstars_xfracHe1[i]
            meanxfrac[i,5] = wstars_xfracHe2[i]
    meanxfrac[0,5]=meanxfrac[0,4]
    meanxfrac[0,2]=meanxfrac[0,1]
   
    #print len(redshifts2), len(meantemp[:,1])
#    print meantemp[:,1]
#    print meantemp[:,0]


    T0=2.725
    Tcmb = np.zeros(len(redshifts1))
    for z in range(len(redshifts1)):
        Tcmb[z] = T0*(1.0+redshifts1[z])

    lw = 2.5

    #print np.log10(meanxfrac[:,5])
    f, (ax1, ax2) = plt.subplots(2, figsize=(15,15),sharex=True)
    ax1.set_xlim(redshifts2[1],redshifts2[len(redshifts2)-1])
    ax1.plot(redshifts2[1:length],np.log10(meanxfrac[1:,0]),color="Red",label="Stellar & X-ray: HII",linewidth=lw)
    ax1.plot(redshifts2[1:length],np.log10(meanxfrac[1:,1]),color="Red",label="Stellar & X-ray: HeII",linewidth=lw,linestyle='--')
    ax1.plot(redshifts2[1:length],np.log10(meanxfrac[1:,2]),color="Red",label="Stellar & X-ray: HeIII",linewidth=lw,linestyle=':')
    ax1.plot(redshifts1[1:length],np.log10(meanxfrac[1:,3]),color="Blue",label="Stellar: HII",linewidth=lw)
    ax1.plot(redshifts1[1:length],np.log10(meanxfrac[1:,4]),color="Blue",label="Stellar: HeII",linewidth=lw,linestyle='--')
    ax1.plot(redshifts1[1:length],np.log10(meanxfrac[1:,5]),color="Blue",label="Stellar: HeIII",linewidth=lw,linestyle=':')
    ax1.set_ylabel(r'log(x$_\mathrm{v}$)',size=size)
    ax1.set_ylim(-17.5,0)
    
    ax2.set_ylim(0,800)
    ax2.plot(redshifts1,Tcmb,label='$T_{cmb}$',color='black',linestyle=':',linewidth=lw)
    ax2.plot(redshifts1[1:length],meantemp[1:,0],color='Blue',label=label1+' mean',linewidth=lw,linestyle='--')
    ax2.plot(redshifts2[1:length],meantemp[1:,1],color='Red',label=label2+' mean',linewidth=lw,linestyle='--')
    ax2.plot(redshifts1[1:length],meantemp[1:,2],color='Blue',label=label1+' median',linewidth=lw)
    ax2.plot(redshifts2[1:length],meantemp[1:,3],color='Red',label=label2+' median',linewidth=lw)
    ax2.set_ylabel("Temperature [K] \n",size=size)
    ax2.set_xlabel("Redshift",size=size+1)
    #ax2.set_ylim(-250,80)

    ax1.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
    ax1.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
    ax2.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11,direction='in',top='off',right='off')
    ax2.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5,direction='in',top='off',right='off')


    lg1 = ax1.legend(bbox_to_anchor=(0.04, 0.4), loc=2, ncol=2, borderaxespad=0.,prop={'size':size-1})
    lg1.draw_frame(False)
    lg = ax2.legend(bbox_to_anchor=(0.04, 0.8), loc=2, borderaxespad=0.,prop={'size':size-1})
    lg.draw_frame(False)

#    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    #plt.tight_layout()
    f.subplots_adjust(hspace=0)


#    fig = plt.figure(figsize=(12,10))
#    ax = fig.add_subplot(111)
#    ax.plot(redshifts2,meantemp[:,0],color='Blue',label=label1,linewidth=lw)
#    ax.plot(redshifts2,meantemp[:,1],color='Red',label=label2,linewidth=lw)
#    ax.plot(redshifts1,Tcmb,label='$T_{cmb}$',color='black',linestyle='--',linewidth=lw)
#    ax.set_ylabel("Temperature [K]",fontsize=size)
#    ax.set_xlabel("Redshift",fontsize=size)
#    ax.xaxis.set_tick_params(labelsize=size)
#    ax.yaxis.set_tick_params(labelsize=size)
    #lg = plt.legend(bbox_to_anchor=(0.55, 0.97), loc=2,prop={'size':size-1})
    #lg.draw_frame(False)
#    plt.xlim(redshifts2[0],redshifts1[len(redshifts1)-1])
#    ax.text(22.7,1100,"(b)",fontsize=size)
#
#    rect = [0.15,0.2,0.55,0.55]
#    ax1 = add_subplot_axes(ax,rect)
#    ax1.set_ylabel("log$_{10}$(xfrac)",fontsize=size-2)
##    ax1.set_xlabel("z",fontsize=size-2)
#    ax1.xaxis.set_tick_params(labelsize=size-2)
#    ax1.yaxis.set_tick_params(labelsize=size-2)  
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,0]),color="Red",label="HII",linewidth=lw)
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,1]),color="Red",label="HeII",linewidth=lw,linestyle='--')
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,2]),color="Red",label="HeIII",linewidth=lw,linestyle=':')
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,3]),color="Blue",label="HII",linewidth=lw)
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,4]),color="Blue",label="HeII",linewidth=lw,linestyle='--')
#    ax1.plot(redshifts2,np.log10(meanxfrac[:,5]),color="Blue",label="HeIII",linewidth=lw,linestyle=':')
#    ax1.text(22,-2,"(a)",fontsize=size)
#    ax2.text(22,700,"(b)",fontsize=size)
#    lg =plt.legend(loc=4,ncol=2,prop={'size':size-3})
#    lg.draw_frame(False)
#    plt.ylim(-16,0)
    #plt.xlim(redshifts2[0],13.221)
 


    print "saving as presenationplots/plot1.png"
    #plt.tight_layout()
    plt.savefig("presentationplots/plot1.png")

def example2():
    fig = plt.figure(figsize=(15,16))
    axes = []

def plot3(type):
    length=0
    if type=='coeval':
        length=min(len(redshifts1),len(redshifts2))
    elif type=='lightcone':
        length=len(lc_zs)
    redshifts_short=np.zeros(length)
    wstars_smooth=np.zeros((length,2))
    wpls_smooth=np.zeros((length,2))
    wstars=np.zeros((length,2))
    wpls=np.zeros((length,2))

    means=np.zeros((len(redshifts1),3))

    means[:,0]=np.loadtxt('../generate_data/data_'+flag1+'/mean_dbt.dat')
    means[:,1]=np.loadtxt('../generate_data/data_'+flag2+'/mean_dbt.dat')
    means[:,2]=np.loadtxt('../generate_data/data_'+flag2+'/mean_dbthightemp.dat')

    wstars_smooth[:,1]=np.loadtxt('../generate_data/data_'+flag1+'/smooth_rms_'+type+'.dat')
    wpls_smooth[:,1]=np.loadtxt('../generate_data/data_'+flag2+'/smooth_rms_'+type+'.dat')
    wpls_hightemp=np.loadtxt('../generate_data/data_'+flag2+'/smooth_rms_hightemp'+type+'.dat')

    z=np.zeros(len(redshifts1))
    L=len(redshifts2)-1
    f, (ax1, ax2) = plt.subplots(2, figsize=(15,15),sharex=True) 
    ax1.set_xlim(redshifts2[1],redshifts2[L]) 
    ax1.plot(redshifts1,z,linewidth=lw+2,color="Black",linestyle=':')
    ax1.plot(redshifts1[0:L:2],means[0:L:2,0],label=label1,linestyle='-',linewidth=lw,color='Blue')
    ax1.plot(redshifts1[0:L:2],means[0:L:2,1],label=label2,linestyle='-',linewidth=lw,color='Red')
    ax1.plot(redshifts1[0:L:2],means[0:L:2,2],linestyle='--',linewidth=lw,color='Orange',label='High temperature limit')
    ax1.set_ylabel(r'$\bar{\delta T}$ [mK]',size=size+1)

    if type=='lightcone':
        ax2.plot(lc_zs,wstars_smooth[:,1],linestyle='-',linewidth=lw,color='Blue',label=label1)
        ax2.plot(lc_zs,wpls_smooth[:,1],linestyle='-',linewidth=lw,color='Red',label=label2)
        ax2.plot(lc_zs,wstars[:,1],linestyle='--',linewidth=lw,color='Orange',label='High temperature limit')
    else:
        ax2.plot(redshifts1[0:L:2],wstars_smooth[0:L:2,1],linestyle='-',linewidth=lw,color='Blue',label=label1)
        ax2.plot(redshifts1[0:L:2],wpls_smooth[0:L:2,1],linestyle='-',linewidth=lw,color='Red',label=label2)
        ax2.plot(redshifts1[0:L:2],wpls_hightemp[0:L:2],linestyle='--',linewidth=lw,color='Orange',label='High temperature limit')
        
    ax2.set_ylabel("\n\n$\delta$T RMS [mK] \n",size=size+1)
    ax2.set_xlabel("Redshift",size=size+1)

    ax1.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
    ax1.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
    ax2.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11,direction='in',top='off',right='off')
    ax2.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5,direction='in',top='off',right='off')

    ax1.set_ylim(-298,200)
    lg1 = ax1.legend(loc=2,prop={'size':size})#(bbox_to_anchor=(0.04, 0.64), loc=2, borderaxespad=0.,prop={'size':size})
    lg1.draw_frame(False)
    #lg = ax2.legend(loc=2,prop={'size':size})#(bbox_to_anchor=(0.04, 0.8), loc=2, borderaxespad=0.,prop={'size':size})
    #lg.draw_frame(False)

    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    plt.tight_layout()
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saving fig as presentationplots/plot3_"+type+".png"
    plt.savefig("presentationplots/plot3_"+type+".png")
    

def plot6(type):

    if type=='coeval':
        length=min(len(redshifts1),len(redshifts2))
        L=length-1
    elif type=='lightcone':
        length=len(lc_zs)
    else:
        print "unknown type"
    redshifts_short=np.zeros(length)
    wstars_smooth=np.zeros((length,2))
    wpls_smooth=np.zeros((length,2))
    wstars=np.zeros((length,2))
    wpls=np.zeros((length,2))
    rmsname=''
    rmsnameht=''
#    wstarsnp.loadtxt('../generate_data/data_'+flag1+'/smooth_skewness.dat')
    wstars_smooth[:,0]=np.loadtxt('../generate_data/data_'+flag1+'/smooth_skewness_'+type+'.dat')
    wpls_smooth[:,0]=np.loadtxt('../generate_data/data_'+flag2+'/smooth_skewness_'+type+'.dat')
    wstars_smooth[:,1]=np.loadtxt('../generate_data/data_'+flag1+'/smooth_kurtosis_'+type+'.dat')
    wpls_smooth[:,1]=np.loadtxt('../generate_data/data_'+flag2+'/smooth_kurtosis_'+type+'.dat')


    #wstars[:,0]=np.loadtxt('../generate_data/data_'+flag1+'/smoothskewnesshightemp.dat')

    wpls[:,0]=np.loadtxt('../generate_data/data_'+flag2+'/smooth_skewness_hightemp'+type+'.dat')

#    wstars[:,1]=np.loadtxt('../generate_data/data_'+flag1+'/smooth_kurtosis_hightemp'+type+'.dat')

    wpls[:,1]=np.loadtxt('../generate_data/data_'+flag2+'/smooth_kurtosis_hightemp'+type+'.dat')



    f, (ax1, ax2) = plt.subplots(2,figsize=(15,15), sharex=True)
#    ax1.plot(redshifts2[0:L:2],wstars[0:L:2,0],linewidth=lw,linestyle="--",color='Blue')
    if type=='lightcone':
        ax1.plot(lc_zs,wstars_smooth[:,0],label=label1,linewidth=lw,linestyle="-",color='Blue')
        ax1.plot(lc_zs,wpls_smooth[:,0],label=label2,linewidth=lw,linestyle="-",color='Red')
        ax1.plot(lc_zs,wpls[:,0],linewidth=lw,linestyle="--",color='Orange',label='High temperature limit')
    else:
        ax1.plot(redshifts2[0:L:2],wstars_smooth[0:L:2,0],linewidth=lw,linestyle="-",color='Blue')
        ax1.plot(redshifts2[0:L:2],wpls_smooth[0:L:2,0],linewidth=lw,linestyle="-",color='Red')
        ax1.plot(redshifts2[0:L:2],wpls[0:L:2,1],linewidth=lw,linestyle="--",color='Orange',label='High temperature limit')
    ax1.set_ylabel("\n\nSkewness",size=size+1)
    ax1.set_xlim(redshifts2[0],redshifts2[len(redshifts2)-1])
    if type=='lightcone':
        ax2.plot(lc_zs,wstars_smooth[:,1],label=label1,linewidth=lw,linestyle="-",color='Blue')
        ax2.plot(lc_zs,wpls_smooth[:,1],label="X-ray & Stellar",linewidth=lw,linestyle="-",color='Red')
        ax2.plot(lc_zs,wpls[:,1],linewidth=lw,linestyle="--",color='Orange',label='High temperature limit')
    else:
	ax2.plot(redshifts1[0:L:2],wpls_smooth[0:L:2,1],label="X-ray & Stellar",linewidth=lw,linestyle="-",color='Red')
        ax2.plot(redshifts1[0:L:2],wstars_smooth[0:L:2,1],label="Stellar Only",linewidth=lw,linestyle="-",color='Blue')
        ax2.plot(redshifts1[0:L:2],wpls[0:L:2,1],linewidth=lw,linestyle="--",color='Orange',label='High temperature limit')

    ax2.set_ylabel("\n\nKurtosis",size=size+1)
    ax2.set_xlabel("Redshift",size=size+1) 

    ax1.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
    ax1.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
    ax2.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11,direction='in',top='off',right='off')
    ax2.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5,direction='in',top='off',right='off')

#    ax1.text(14,1.75,"(c)",fontsize=size)
#    ax2.text(14,4.5,"(d)",fontsize=size)

    #lg = ax1.legend(loc=2,prop={'size':size})#bbox_to_anchor=(0.13, 0.35),loc=2,prop={'size':size})
    #lg.draw_frame(False)
    #lg2 = ax2.legend(loc=2,prop={'size':size})#bbox_to_anchor=(0.13, 1.0),loc=2,prop={'size':size})
    #lg2.draw_frame(False)

    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    plt.tight_layout()
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saving figure as presentationplots/plot6_"+type+'.png'
    plt.savefig("presentationplots/plot6_"+type+".png")

def plot1b():
#redshifts4 = [23.268,18.910,15.596,13.733]
    nobins=100
    redshifts_short = np.ones(4).reshape(2,2)
    redshifts_short[0,0]=16.359
    redshifts_short[0,1]=14.493
    redshifts_short[1,0]= 13.733
    redshifts_short[1,1]=13.221
    print redshifts_short.shape
    length=min(len(redshifts1),len(redshifts2))

    f, ax = plt.subplots(2,2, figsize=(15,15),sharex=True,sharey=True)
    for i in range(2):
        for j in range(2):
            data = np.zeros(2*250**3).reshape(250**3,2)
            temp_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag2+'/Temper3D_'+str('%.3f' % redshifts_short[i,j]) + '.bin'
            dbt_wpls = c2t.TemperFile(temp_filename).temper
            data[:,0]=dbt_wpls.flatten()
            temp_filename = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+flag1+'/Temper3D_'+str('%.3f' % redshifts_short[i,j]) + '.bin'
            dbt_wstars = c2t.TemperFile(temp_filename).temper
            data[:,1] = dbt_wstars.flatten() #need to import other dataset and put in right format
            weights = np.ones_like(data[:,1])/len(data)
            ax[i,j].hist(np.log10(data[:,1]),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='blue',edgecolor='none',alpha=0.5,label=label1)

            ax[i,j].hist(np.log10(data[:,0]),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='Red',edgecolor='none',alpha=0.5,label="X-ray & \n Stellar")
            ax[i,j].text(1.75,0.8,"z = "+str(redshifts_short[i,j]),fontsize=fontsize,color='Purple')

            ax[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
            ax[i,j].tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
            ax[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11,direction='in',top='off',right='off')
            ax[i,j].tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5,direction='in',top='off',right='off')
#
           # ax[i,j].set_ylim(0,1)
           # ax[i,j].set_xlim(0,4)

#            m0 = MultipleLocator(0.1)
#            ax[i,j].yaxis.set_major_locator(m0)
#            m0 = MultipleLocator(0.02)
#            ax[i,j].yaxis.set_minor_locator(m0)
#    
#            m0 = MultipleLocator(50)
#            ax[i,j].yaxis.set_major_locator(m0)
#            m0 = MultipleLocator(10)
#            ax[i,j].yaxis.set_minor_locator(m0)
#
            m0 = MultipleLocator(1.0)
            ax[i,j].xaxis.set_major_locator(m0)
            m0 = MultipleLocator(0.5)
            ax[i,j].xaxis.set_minor_locator(m0)
#    
#            m0 = MultipleLocator(2)
#            ax[i,j].xaxis.set_major_locator(m0)
#            m0 = MultipleLocator(1)
#            ax[i,j].xaxis.set_minor_locator(m0)

            if i==0 and j==0:
                lg1 = ax[i,j].legend(bbox_to_anchor=(0.35, 0.75), loc=2, borderaxespad=0.,prop={'size':size})
                lg1.draw_frame(False)
                lg = ax[i,j].legend(bbox_to_anchor=(0.35, 0.75), loc=2, borderaxespad=0.,prop={'size':size})
                lg.draw_frame(False)

            ax[i,j].set_xlim(0,3.9)
            ax[i,j].set_ylim(0,0.95)
    ax[1,1].set_xlabel('log$_{10}$(T) [K]',size=size)
    ax[1,0].set_xlabel("log$_{10}$(T) [K]",size=size)
    ax[0,0].set_ylabel("\n\nProbability\n",size=size)
    ax[1,0].set_ylabel("\n\nProbability\n",size=size)

#    ax[0,0].text(22,0.8,"(c)",fontsize=size)
  #labels1=["","0.0","","0.2","","0.4","","0.6","","0.8"]
#    labels12=["","0.0","","0.2","","0.4","","0.6","","0.8","","1.0"]
#    labels2=["","0.0","","1.0","","2.0","","3.0",""]
#    labels22=["","0.0","","1.0","","2.0","","3.0","","4.0"]
    nolabels=[" "," "," "," "," "," "," "," "," "]
 #   ax[0,0].set_yticklabels(labels12,size=size)
 #   ax[1,0].set_yticklabels(labels1,size=size)
 ##   ax[1,1].set_yticklabels(nolabels,size=size)
 ##   ax[0,1].set_yticklabels(nolabels,size=size)
##    ax[1,0].set_yticklabels(labels1,size=size)
 #   ax[1,1].set_xticklabels(labels22,size=size)
 #   ax[1,0].set_xticklabels(labels2,size=size)
 ##   ax[0,0].set_xticklabels(nolabels,size=size)
 ##   ax[0,1].set_xticklabels(nolabels,size=size)

#    plt.tight_layout()
#    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    f.subplots_adjust(hspace=0,wspace=0)
#    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saving fig as presentationplots/plot1b.png"
    plt.savefig("presentationplots/plot1b.png")

def kplot():
    #k=0.1 is in positiion 6
    print "hello"
    ks_wstars=np.zeros((len(redshifts1),3))
    ks_wpls=np.zeros((len(redshifts1),3))
    ks_ht=np.zeros((len(redshifts1),3))

    for z in range(len(redshifts1)):
        ps_wstars=np.loadtxt("../generate_data/data_"+flag1+'/powerSpectra_100b_'+str('%.3f'%redshifts1[z]+'.dat'))
        ps_wpls=np.loadtxt("../generate_data/data_"+flag2+'/powerSpectra_100b_'+str('%.3f'%redshifts1[z]+'.dat'))
        ps_ht=np.loadtxt("../generate_data/data_"+flag2+'/powerSpectra_100b_hightemp'+str('%.3f'%redshifts1[z]+'.dat'))
        fr=np.loadtxt("../generate_data/data_"+flag2+"/powerSpectraFrequencies_dbt_100b_"+str('%.3f'%redshifts1[z]+'.dat'))

        ks_wstars[z,0]=ps_wstars[6]*fr[6]**3/(4.*np.pi**2)
        ks_wpls[z,0]=ps_wpls[6]*fr[6]**3/(4.*np.pi**2)
        ks_ht[z,0]=ps_ht[6]*fr[6]**3/(4.*np.pi**2)

        ks_wstars[z,1]=ps_wstars[48]*fr[48]**3/(4.*np.pi**2)
        ks_wpls[z,1]=ps_wpls[48]*fr[48]**3/(4.*np.pi**2)
        ks_ht[z,1]=ps_ht[48]*fr[48]**3/(4.*np.pi**2)

        ks_wstars[z,2]=ps_wstars[35]*fr[35]**3/(4.*np.pi**2)
        ks_wpls[z,2]=ps_wpls[35]*fr[35]**3/(4.*np.pi**2)
        ks_ht[z,2]=ps_ht[35]*fr[35]**3/(4.*np.pi**2)

    ks_wstars=np.sqrt(ks_wstars)
    ks_wpls=np.sqrt(ks_wpls)
    ks_ht=np.sqrt(ks_ht)

    L=len(redshifts1)-1

    fig,ax=plt.subplots(3,1,figsize=(15,15),sharex=True,sharey=True)
    #lg.draw_frame(False)#

    #k=0.1
    ax[0].plot(redshifts1[0:L:2],np.log10(ks_wpls[0:L:2,0]),label="X-ray & Stellar",color='red',linewidth=lw+2,linestyle='--')
    ax[0].plot(redshifts1[0:L:2],np.log10(ks_wstars[0:L:2,0]),label="Stellar",color='blue',linewidth=lw+2,linestyle='--')
    ax[0].plot(redshifts1[0:L:2],np.log10(ks_ht[0:L:2,0]),label="High T limit",color='orange',linewidth=lw+2,linestyle='--')
    ax[0].text(16,1.75,"k = 0.1 Mpc$^{-1}$",fontsize=size+5)

    #k=0.5
    ax[1].plot(redshifts1[0:L:2],np.log10(ks_wpls[0:L:2,2]),color='red',linewidth=lw+2,linestyle='--')
    ax[1].plot(redshifts1[0:L:2],np.log10(ks_wstars[0:L:2,2]),color='blue',linewidth=lw+2,linestyle='--')
    ax[1].plot(redshifts1[0:L:2],np.log10(ks_ht[0:L:2,2]),color='orange',linewidth=lw+2,linestyle='--')
    ax[1].text(16,1.75,"k = 0.5 Mpc$^{-1}$",fontsize=size+5)

    #k=1.0
    ax[2].plot(redshifts1[0:L:2],np.log10(ks_wpls[0:L:2,1]),color='red',linewidth=lw+2,linestyle='--')
    ax[2].plot(redshifts1[0:L:2],np.log10(ks_wstars[0:L:2,1]),color='blue',linewidth=lw+2,linestyle='--')
    ax[2].plot(redshifts1[0:L:2],np.log10(ks_ht[0:L:2,1]),color='orange',linewidth=lw+2,linestyle='--')
    ax[2].text(16,1.75,"k = 1.0 Mpc$^{-1}$",fontsize=size+5)

    lg=ax[0].legend(loc=2,prop={'size':size+5})
    lg.draw_frame(False)
    plt.xlabel("Redshifts",size=size+5)
    ax[0].set_ylabel("log$_{10}$($\Delta_{21cm}$)",size=size+5)
    ax[1].set_ylabel("log$_{10}$($\Delta_{21cm}$)",size=size+5)
    ax[2].set_ylabel("log$_{10}$($\Delta_{21cm}$)",size=size+5)
    plt.ylim(0.0,2.2)
    plt.xlim(redshifts1[0],redshifts1[L])
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig("presentationplots/kplot5.png")

kplot()
#plot1()
#plot1b()
#plot3("coeval")
#plot6("coeval")
#plot3("lightcone")
#plot6("lightcone")
