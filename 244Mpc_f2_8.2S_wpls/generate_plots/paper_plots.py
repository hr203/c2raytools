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

plt.rc('xtick', labelsize=24)
plt.rc('ytick', labelsize=24)
plt.rc('font',family='serif')
fontsize=24#3dd0
numberfontsize=24
tickwidth=1.5

size=24#26
lw = 2.5


flag1="wstars"
flag2="wpls"
redshifts1=setup_dirs.read_redshifts(flag1)
redshifts2=setup_dirs.read_redshifts(flag2)
if flag1=='wstars'or flag2=='wstars':
    label1="Stellar"
if flag1=='wpls'or flag2=='wpls':
    if flag1 =='wstars' or flag2=='wstars':
        label2="X-Ray & Stellar"
    else:
        label1="X-Ray & Stellar"
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
    data1 = plot_1D.mean('temp','data_'+flag1+'/',redshifts1)
    data2 = plot_1D.mean('temp','data_'+flag2+'/',redshifts2)
   
    #load xfrac means
    wstars_xfrac = plot_1D.mean('xfrac','data_'+flag1+'/',redshifts1)
    wpls_xfrac = plot_1D.mean('xfrac','data_'+flag2+'/',redshifts2)
    wstars_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag1+'/',redshifts1)
    wpls_xfracHe1 = plot_1D.mean('xfracHe1','data_'+flag2+'/',redshifts2)
    wstars_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag1+'/',redshifts1)
    wpls_xfracHe2 = plot_1D.mean('xfracHe2','data_'+flag2+'/',redshifts2)

    if minmax=="max":
        length = max(len(data1),len(data2))
    elif minmax=="min":
        length = min(len(data1),len(data2))
    print redshifts2 
    meantemp = np.zeros(2*length).reshape(length,2)
    meanxfrac = np.zeros(6*length).reshape(length,6)
    if minmax=="max":
        for i in range(len(data1)):
            meantemp[i,0] = data1[i]
            meanxfrac[i,3] = wstars_xfrac[i]
            meanxfrac[i,4] = wstars_xfracHe1[i]
            meanxfrac[i,5] = wstars_xfracHe2[i]

        for i in range(len(data2)):
            meanxfrac[i,0] = wpls_xfrac[i]
            meanxfrac[i,1] = wpls_xfracHe1[i]
            meanxfrac[i,2] = wpls_xfracHe2[i]
            meantemp[i,1] = data2[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
#        plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_'+flag1+'_'+flag2+'/',redshifts1,redshifts1,label1,label2)

    elif minmax=="min":
        for i in range(length):
            meantemp[i,0] = data1[i]
            meantemp[i,1] = data2[i]

            meanxfrac[i,0] = wpls_xfrac[i]
            meanxfrac[i,1] = wpls_xfracHe1[i]
            meanxfrac[i,2] = wpls_xfracHe2[i]
            meanxfrac[i,3] = wstars_xfrac[i]
            meanxfrac[i,4] = wstars_xfracHe1[i]
            meanxfrac[i,5] = wstars_xfracHe2[i]
    meanxfrac[0,5]=meanxfrac[0,4]
    meanxfrac[0,2]=meanxfrac[0,1]
   
    print len(redshifts2), len(meantemp[:,1])
    print meantemp[:,1]
    print meantemp[:,0]


    T0=2.725
    Tcmb = np.zeros(len(redshifts1))
    for z in range(len(redshifts1)):
        Tcmb[z] = T0*(1.0+redshifts1[z])

    size=26
    lw = 2.5

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    ax.plot(redshifts2,meantemp[:,0],color='Blue',label=label1,linewidth=lw)
    ax.plot(redshifts2,meantemp[:,1],color='Red',label=label2,linewidth=lw)
    ax.plot(redshifts1,Tcmb,label='$T_{cmb}$',color='black',linestyle='--',linewidth=lw)
    ax.set_ylabel("Temperature [K]",fontsize=size)
    ax.set_xlabel("Redshift",fontsize=size)
    ax.xaxis.set_tick_params(labelsize=size)
    ax.yaxis.set_tick_params(labelsize=size)
    lg = plt.legend(bbox_to_anchor=(0.55, 0.97), loc=2,prop={'size':size-1})
    lg.draw_frame(False)
    plt.xlim(redshifts2[0],redshifts2[len(redshifts2)-1])
    ax.text(22.7,1100,"(b)",fontsize=size)

    rect = [0.15,0.2,0.55,0.55]
    ax1 = add_subplot_axes(ax,rect)
    ax1.set_ylabel("log$_{10}$(xfrac)",fontsize=size-2)
#    ax1.set_xlabel("z",fontsize=size-2)
    ax1.xaxis.set_tick_params(labelsize=size-2)
    ax1.yaxis.set_tick_params(labelsize=size-2)  
    ax1.plot(redshifts2,np.log10(meanxfrac[:,0]),color="Red",label="HII",linewidth=lw)
    ax1.plot(redshifts2,np.log10(meanxfrac[:,1]),color="Red",label="HeII",linewidth=lw,linestyle='--')
    ax1.plot(redshifts2,np.log10(meanxfrac[:,2]),color="Red",label="HeIII",linewidth=lw,linestyle=':')
    ax1.plot(redshifts2,np.log10(meanxfrac[:,3]),color="Blue",label="HII",linewidth=lw)
    ax1.plot(redshifts2,np.log10(meanxfrac[:,4]),color="Blue",label="HeII",linewidth=lw,linestyle='--')
    ax1.plot(redshifts2,np.log10(meanxfrac[:,5]),color="Blue",label="HeIII",linewidth=lw,linestyle=':')
    ax1.text(22,-2,"(a)",fontsize=size-2)
    lg =plt.legend(loc=4,ncol=2,prop={'size':size-3})
    lg.draw_frame(False)
    plt.ylim(-16,0)
    plt.xlim(redshifts2[0],redshifts2[len(redshifts2)-1])
 



    plt.tight_layout()
    plt.savefig("exampe.png")

def example2():
    fig = plt.figure(figsize=(10,10))
    axes = []

def plot3():
    length=min(len(redshifts1),len(redshifts2))
#    if length!=0:
#        length=length-1
    wstars=np.zeros(length).reshape(length/2,2)
    wpls=np.zeros(length).reshape(length/2,2)
    redshifts_short=np.zeros((length)/2)

    file=open('../generate_data/data_'+flag1+'/mean_dbt.dat')
    print 'Read mean from ../generate_data/data'+flag1+'/mean_dbt.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wstars[i/2,0]=float(line)
        i+=1
    file.close

    file=open('../generate_data/data_'+flag2+'/mean_dbt.dat')
    print 'Read mean from ../generate_data/data_'+flag2+'/mean_dbt.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wpls[i/2,0]=float(line)
            redshifts_short[i/2]=redshifts2[i]
        i+=1    
    file.close

    file=open('../generate_data/data_'+flag1+'/rms.dat')
    print 'Read mean from ../generate_data/data_'+flag1+'/rms.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wstars[i/2,1]=float(line)
        i+=1
    file.close

    file=open('../generate_data/data_'+flag2+'/rms.dat')
    print 'Reading mean from ../generate_data/data_'+flag2+'/rms.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wpls[i/2,1]=float(line)
        i+=1
    file.close

    z=np.zeros(len(redshifts_short))

    f, (ax1, ax2) = plt.subplots(2, figsize=(10,10),sharex=True) 
    ax1.set_xlim(redshifts_short[0],redshifts2[len(redshifts2)-1]) 
    ax1.plot(redshifts_short,z,linewidth=lw,color="Black",linestyle='--')
    ax1.plot(redshifts_short,wstars[:,0],label=label1,linestyle='-',linewidth=lw,color='Blue')
    ax1.plot(redshifts_short,wpls[:,0],label=label2,linestyle='-',linewidth=lw,color='Red')
    ax1.set_ylabel(r'$\bar{\delta T}$ [mK]',size=size)
    ax2.plot(redshifts_short,wstars[:,0],linestyle='-',linewidth=lw,color='Blue',label=label1)
    ax2.plot(redshifts_short,wpls[:,1],linestyle='-',linewidth=lw,color='Red',label=label2)
    ax2.set_ylabel("$\delta$T RMS [mK]",size=size)
    ax2.set_xlabel("Redshifts",size=size)
    ax2.set_ylim(-250,80)

    ax1.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
    ax1.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
    ax2.tick_params(axis='both', which='major', width = tickwidth, length = 11,direction='in',top='off',right='off')
    ax2.tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='in',top='off',right='off')

    m0 = MultipleLocator(50)
    ax1.yaxis.set_major_locator(m0)
    m0 = MultipleLocator(10)
    ax1.yaxis.set_minor_locator(m0)

    m0 = MultipleLocator(50)
    ax2.yaxis.set_major_locator(m0)
    m0 = MultipleLocator(10)
    ax2.yaxis.set_minor_locator(m0)

    m0 = MultipleLocator(2)
    ax1.xaxis.set_major_locator(m0)
    m0 = MultipleLocator(1)
    ax1.xaxis.set_minor_locator(m0)

    m0 = MultipleLocator(2)
    ax2.xaxis.set_major_locator(m0)
    m0 = MultipleLocator(1)
    ax2.xaxis.set_minor_locator(m0)

    labels1=["","","-200","","-100","","0",""]
    labels2=["","","-200","","-100","","0",""]
    ax1.set_yticklabels(labels1,size=size)
    ax2.set_yticklabels(labels2,size=size)

    lg1 = ax1.legend(bbox_to_anchor=(0.04, 0.65), loc=2, borderaxespad=0.,prop={'size':size})
    lg1.draw_frame(False)
    lg = ax2.legend(bbox_to_anchor=(0.04, 0.65), loc=2, borderaxespad=0.,prop={'size':size})
    lg.draw_frame(False)

    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    #plt.tight_layout()
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saving fig as paperplots/plot2.png"
    plt.savefig("paperplots/plot3.png")
    

def plot6():

    length=min(len(redshifts1),len(redshifts2))
    #if length!=0:
    #    length=length-1
    wstars=np.zeros(length).reshape(length/2,2)
    wpls=np.zeros(length).reshape(length/2,2)
    redshifts_short=np.zeros(length/2)


    file=open('../generate_data/data_'+flag1+'/skewness.dat')
    print 'Read mean from ../generate_data/data'+flag1+'/skewness.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wstars[i/2,0]=float(line)
        i+=1
#        print len(wstars), i, line

    file.close

    file=open('../generate_data/data_'+flag2+'/skewness.dat')
    print 'Read mean from ../generate_data/data_'+flag2+'/skewness.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wpls[i/2,0]=float(line)
            redshifts_short[i/2]=redshifts2[i]
        i+=1
    file.close

    file=open('../generate_data/data_'+flag1+'/kurtosis.dat')
    print 'Read mean from ../generate_data/data_'+flag1+'/kurtosis.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wstars[i/2,1]=float(line)
        i+=1
    file.close

    file=open('../generate_data/data_'+flag2+'/kurtosis.dat')
    print 'Reading mean from ../generate_data/data_'+flag2+'/kurtosis.dat'
    i=-1
    for line in file:
        if i!=-1 and i<length and i%2==0:
            wpls[i/2,1]=float(line)
        i+=1



#    for i in range(0,length,2):
#        print "opening ../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin and ../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin"
#        data1=np.load("../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin")
#        data2=np.load("../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin")
#        skewnessdata[i,0] = c2t.statistics.skewness(data1)
#        skewnessdata[i,1] = c2t.statistics.skewness(data2)
#        kurtosisdata[i,0] = c2t.statistics.kurtosis(data1)
#        kurtosisdata[i,1] = c2t.statistics.kurtosis(data2)
    f, (ax1, ax2) = plt.subplots(2,figsize=(10,10), sharex=True)
    ax1.plot(redshifts_short,wstars[:,0],label=label1,linewidth=lw,linestyle="-",color='Blue')
    ax1.plot(redshifts_short,wpls[:,0],label=label2,linewidth=lw,linestyle="-",color='Red')
    ax1.set_ylabel("Skewness",size=size)
    ax1.set_xlim(redshifts_short[0],redshifts2[len(redshifts2)-1])
    ax2.plot(redshifts_short,wstars[:,1],label=label1,linewidth=lw,linestyle="-",color='Blue')
    ax2.plot(redshifts_short,wpls[:,1],label="X-Ray & \n Stellar",linewidth=lw,linestyle="-",color='Red')
    ax2.set_ylabel("Kurtosis",size=size)
    ax2.set_xlabel("Redshift",size=size) 
    ax2.set_ylim(0,19)
    ax1.set_ylim(-0.9,2.5)

    ax1.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
    ax1.tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
    ax2.tick_params(axis='both', which='major', width = tickwidth, length = 11,direction='in',top='off',right='off')
    ax2.tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='in',top='off',right='off')


    labels2=["","0$ \ $","5$ \ $","10$ \ $","15$ \ $"]
    labels1=["","","0.0","","1.0","","2.0",""]
    ax1.set_yticklabels(labels1,size=size)
    ax2.set_yticklabels(labels2,size=size)


    lg = ax1.legend(bbox_to_anchor=(0.4, 0.4),loc=2,prop={'size':size})
    lg.draw_frame(False)

    lg2 = ax2.legend(bbox_to_anchor=(0.5, 1.0),loc=2,prop={'size':size})
    lg2.draw_frame(False)

    m0 = MultipleLocator(0.5)
    ax1.yaxis.set_major_locator(m0)
    m0 = MultipleLocator(0.1)
    ax1.yaxis.set_minor_locator(m0)

    m0 = MultipleLocator(5)
    ax2.yaxis.set_major_locator(m0)
    m0 = MultipleLocator(1)
    ax2.yaxis.set_minor_locator(m0)

    m0 = MultipleLocator(2)
    ax2.xaxis.set_major_locator(m0)
    m0 = MultipleLocator(1)
    ax2.xaxis.set_minor_locator(m0)



    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
#    ax3.scatter(x, 2 * y ** 2 - 1, color='r')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saveing figure as paperplots/plot6"
    plt.savefig("paperplots/plot6.png")

def plot1b():
    nobins=100
    redshifts_short = np.ones(4).reshape(2,2)
    redshifts_short[0,0]=22.666
    redshifts_short[0,1]=18.910
    redshifts_short[1,0]= 15.360
    redshifts_short[1,1]=13.060
    print redshifts_short.shape
    length=min(len(redshifts1),len(redshifts2))

    f, ax = plt.subplots(2,2, figsize=(12,10))#,sharex=True,sharey=True)
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

            ax[i,j].hist(np.log10(data[:,0]),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='Red',edgecolor='none',alpha=0.5,label="X-Ray & \n Stellar")
            ax[i,j].text(1.75,0.8,"z = "+str(redshifts_short[i,j]),fontsize=fontsize,color='Purple')

            ax[i,j].tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, direction='in',pad=14.0,top='off',right='off')
            ax[i,j].tick_params(axis='both', which='minor', labelsize=numberfontsize, width = tickwidth, length = 5.5, direction='in',pad=14.0,top='off',right='off')
            ax[i,j].tick_params(axis='both', which='major', width = tickwidth, length = 11,direction='in',top='off',right='off')
            ax[i,j].tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,direction='in',top='off',right='off')
#
            ax[i,j].set_ylim(0,1)
            ax[i,j].set_xlim(0,4)

            m0 = MultipleLocator(0.1)
            ax[i,j].yaxis.set_major_locator(m0)
            m0 = MultipleLocator(0.02)
            ax[i,j].yaxis.set_minor_locator(m0)
#    
#            m0 = MultipleLocator(50)
#            ax[i,j].yaxis.set_major_locator(m0)
#            m0 = MultipleLocator(10)
#            ax[i,j].yaxis.set_minor_locator(m0)
#
            m0 = MultipleLocator(0.5)
            ax[i,j].xaxis.set_major_locator(m0)
            m0 = MultipleLocator(0.1)
            ax[i,j].xaxis.set_minor_locator(m0)
#    
#            m0 = MultipleLocator(2)
#            ax[i,j].xaxis.set_major_locator(m0)
#            m0 = MultipleLocator(1)
#            ax[i,j].xaxis.set_minor_locator(m0)


            lg1 = ax[i,j].legend(bbox_to_anchor=(0.3, 0.75), loc=2, borderaxespad=0.,prop={'size':size})
            lg1.draw_frame(False)
            lg = ax[i,j].legend(bbox_to_anchor=(0.3, 0.75), loc=2, borderaxespad=0.,prop={'size':size})
            lg.draw_frame(False)

#    ax1.set_xlim(redshifts_short[0],redshifts2[len(redshifts2)-1])
    ax[1,1].set_xlabel('log$_{10}$(Temperature) [K]',size=size)
    ax[1,0].set_xlabel("log$_{10}$(Temperature) [K]",size=size)
    ax[0,0].set_ylabel("Probability",size=size)
    ax[1,0].set_ylabel("Probability",size=size)


    labels1=["","0.0","","0.2","","0.4","","0.6","","0.8"]
    labels12=["","0.0","","0.2","","0.4","","0.6","","0.8","","1.0"]
    labels2=["","0.0","","1.0","","2.0","","3.0",""]
    labels22=["","0.0","","1.0","","2.0","","3.0","","4.0"]
    nolabels=[" "," "," "," "," "," "," "," "," "]
    ax[0,0].set_yticklabels(labels12,size=size)
    ax[1,0].set_yticklabels(labels1,size=size)
    ax[1,1].set_yticklabels(nolabels,size=size)
    ax[0,1].set_yticklabels(nolabels,size=size)
#    ax[1,0].set_yticklabels(labels1,size=size)
    ax[1,1].set_xticklabels(labels22,size=size)
    ax[1,0].set_xticklabels(labels2,size=size)
    ax[0,0].set_xticklabels(nolabels,size=size)
    ax[0,1].set_xticklabels(nolabels,size=size)

    plt.tight_layout()
#    plt.gcf().subplots_adjust(left=0.2,bottom=0.15)
    f.subplots_adjust(hspace=0,wspace=0)
#    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    print "saving fig as paperplots/plot2.png"
    plt.savefig("paperplots/plot1b.png")



plot1b()
#plot3()
#plot6()
