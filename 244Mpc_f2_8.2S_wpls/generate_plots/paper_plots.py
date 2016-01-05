import matplotlib.pyplot as plt
import numpy as np
#import compare
import sys
import plot_1D
sys.path.append('../')
import setup_dirs
sys.path.append('../../src')
import c2raytools as c2t

flag1="wstars"
flag2="wpls"
redshifts1=setup_dirs.read_redshifts(flag1)
redshifts2=setup_dirs.read_redshifts(flag2)
if flag1=='wstars'or flag2=='wstars':
    label1="Stellar Only"
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

def plot2():
    print "calculating mean"
    length=min(len(redshifts1),len(redshifts2))
    if length!=0:
        length=length-1
    rmsdata=np.zeros(length+1).reshape((length+1)/2,2)
    meandata=np.zeros(length*2).reshape(length,2)
    redshifts_short=np.zeros((length)/2)
    for i in range(0,length):
        print "opening ../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin and ../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin"
        data1=np.load("../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin")
        data2=np.load("../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin")
        meandata[i,0]=np.mean(data1)
        meandata[i,1]=np.mean(data2)
        print "calculating rms"
        if i%2==0 and i>len(redshifts_short):
            redshifts_short[i/2]=redshifts2[i]
            mean2=np.mean(data2)
            mean1=np.mean(data1)
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



   # x = np.linspace(0, 2 * np.pi, 400)
   # y = np.sin(x ** 2)
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    ax1.plot(redshifts2, meandata[:,0],label=label1,color='Blue')
    ax1.plot(redshifts2, meandata[:,1],label=label2,color='Red')
    ax1.set_title('Sharing both axes')
    ax2.scatter(redshifts2, rmsdata[:,0],label=label1,color='Blue')
    ax2.scatter(redshifts2, rmsdata[:,1],label=label2,color='Red')
#    ax3.scatter(x, 2 * y ** 2 - 1, color='r')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig("plot2.png")
    
    for i in range(len(meandata)):
        print meandata[i]

    for i in range(len(rmsdata)):
        print meanrms[i]

def plot6():

    length=min(len(redshifts1),len(redshifts2))
    #if length!=0:
    #    length=length-1
    skewnessdata=np.zeros(length*2).reshape((length),2)
    kurtosisdata=np.zeros(length*2).reshape(length,2)
    redshifts_short=np.zeros((length)/2)
    for i in range(0,length,2):
        print "opening ../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin and ../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin"
        data1=np.load("../generate_data/data_"+flag1+"/map_dbt_"+str('%.3f' % redshifts1[i])+".bin")
        data2=np.load("../generate_data/data_"+flag2+"/map_dbt_"+str('%.3f' % redshifts2[i])+".bin")
        skewnessdata[i,0] = c2t.statistics.skewness(data1)
        skewnessdata[i,1] = c2t.statistics.skewness(data2)
        kurtosisdata[i,0] = c2t.statistics.kurtosis(data1)
        kurtosisdata[i,1] = c2t.statistics.kurtosis(data2)
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    ax1.plot(redshifts2, skewnessdata[:,0],label=label1,color='Blue')
    ax1.plot(redshifts2, skewnessdata[:,1],label=label2,color='Red')
    ax1.set_title('Sharing both axes')
    ax2.scatter(redshifts2, kurtosisdata[:,0],label=label1,color='Blue')
    ax2.scatter(redshifts2, kurtosisdata[:,1],label=label2,color='Red')
#    ax3.scatter(x, 2 * y ** 2 - 1, color='r')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig("plot6.png")

plot6()
