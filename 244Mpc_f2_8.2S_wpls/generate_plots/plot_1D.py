import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('../../src/')
import c2raytools as c2t
sys.path.append('../')
import setup_dirs
import IO
import scipy.integrate as sci
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
mpl.rc('xtick', labelsize=25)
mpl.rc('ytick', labelsize=25)
mpl.rc('font',family='serif')
fontsize=26
numberfontsize=25
tickwidth=1.5
start=0

#print plotsdir()
redshifts = setup_dirs.read_redshifts()
path_data = setup_dirs.resultsdir()
path_plots = setup_dirs.plotsdir()

def mean(name,pathd = path_data,redshifts1=redshifts):
    mean = np.zeros(len(redshifts1))
    print 'Reading mean from ../generate_data/'+pathd+'mean_'+name+'.dat'
    file = open('../generate_data/'+pathd+'mean_'+name+'.dat', 'r')
    i=0
    for line in file:
        mean[i] = float(line)
        i=i+1
    return mean

def plot_mean(data,title,ylabel,pathp = path_plots, redshifts1=redshifts, redshifts2=redshifts,legend1='',legend2='',legend3='',legend4=''):
    #plt.figure()
    fig=plt.figure(figsize=(10.25, 10.25), dpi= 300)
#    plt.subplot('111')
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=numberfontsize, width = tickwidth, length = 11, pad=14.0,top='off',right='off')
    ax.tick_params(axis='both', which='minor', width = tickwidth, length = 5.5,top='off',right='off')
    print data.shape
    print len(redshifts1)
    if(data.shape==(len(redshifts1),)):  
        plt.plot(redshifts1,data)
    elif(data.shape==(len(redshifts1),2) or data.shape==((len(redshifts1))/2,2) or  data.shape==((len(redshifts1)-1)/2,2) or data.shape==(len(redshifts2),2) or data.shape==(len(redshifts2)/2,2) or data.shape==((len(redshifts2)-1)/2,2)):
        if legend1 == '' or legend2 =='':
            print "WARNING: Specify legend for plotting two data sets"
        plt.gcf().subplots_adjust(left=0.15)
        l = 2.0 #linewidth value
        linestyle2='--'
        #plot CMB temperature
        if title == "dbt":
            #yaxis
            m0 = MultipleLocator(50)
            ax.yaxis.set_major_locator(m0)
            m0 = MultipleLocator(10)
            ax.yaxis.set_minor_locator(m0)
            zeros = np.zeros(len(redshifts1))
            plt.plot(redshifts1,zeros,color='black',linestyle='--',linewidth=l)
            #plt.ylim(-250,100)
        if title=="Temperature (K)":
            linestyle2 ='-'
            T0=2.725
            Tcmb = np.zeros(len(redshifts1))
            for z in range(len(redshifts1)):
                Tcmb[z] = T0*(1.0+redshifts1[z])
            plt.plot(redshifts1,Tcmb,label='$T_{cmb}$',color='black',linestyle='--',linewidth=l)
                    #plt.legend(loc=4,prop={'size':18})
            #yaxis
            m0 = MultipleLocator(100)
            ax.yaxis.set_major_locator(m0)
            m0 = MultipleLocator(20)
            ax.yaxis.set_minor_locator(m0)
        #x-axis
        m0 = MultipleLocator(1)
        ax.xaxis.set_major_locator(m0)
        m0 = MultipleLocator(0.5)
        ax.xaxis.set_minor_locator(m0)

#        print data[0 
#        print legend1, legend2 
        plt.plot(redshifts1,data[0:len(redshifts1),0],label=legend1,color='blue',linestyle=linestyle2,linewidth=l)
        plt.plot(redshifts1,data[0:len(redshifts1),1],label=legend2,color='red',linestyle=linestyle2,linewidth=l)
        lg = plt.legend(loc=2,prop={'size':fontsize})
        lg.draw_frame(False)

    elif (data.shape==(len(redshifts1),3) or data.shape==(len(redshifts2),3)):
        plt.plot(redshifts1,data[:,0],label=legend1)
        plt.plot(redshifts1,data[:,1],label=legend2)
        plt.plot(redshifts1,data[:,2],label=legend3)

        if legend1 == '' or legend2=='' or legend3=='':
            print "WARNING: Specify legend for plotting three data sets"
        else:
            lg = plt.legend(loc=4,prop={'size':fontsize})
            lg.draw_frame(False)

    elif(data.shape==(len(redshifts1),4) or data.shape==((len(redshifts1))/2,4) or  data.shape==((len(redshifts1)-1)/2,4) or data.shape==(len(redshifts2),4) or data.shape==(len(redshifts2)/2,4) or data.shape==((len(redshifts2)-1)/2,4)):
        if legend1 == '' or legend2 =='' or legend3 == '' or legend4=='':
            print "WARNING: Specify legend for plotting two data sets"
        plt.gcf().subplots_adjust(left=0.15)
        l = 2.0 #linewidth value
        linestyle2='--'
        #plot CMB temperature
        if title == "dbt":
            #yaxis
            m0 = MultipleLocator(50)
            ax.yaxis.set_major_locator(m0)
            m0 = MultipleLocator(10)
            ax.yaxis.set_minor_locator(m0)
            zeros = np.zeros(len(redshifts1))
            plt.plot(redshifts1,zeros,color='black',linestyle='--',linewidth=l)
            #plt.ylim(-250,100)
        if title=="Temperature (K)":
            linestyle2 ='-'
            T0=2.725
            Tcmb = np.zeros(len(redshifts1))
            for z in range(len(redshifts1)):
                Tcmb[z] = T0*(1.0+redshifts1[z])
            plt.plot(redshifts1,Tcmb,label='$T_{cmb}$',color='black',linestyle=':',linewidth=l)
                    #plt.legend(loc=4,prop={'size':18})
            #yaxis
            m0 = MultipleLocator(100)
            ax.yaxis.set_major_locator(m0)
            m0 = MultipleLocator(20)
            ax.yaxis.set_minor_locator(m0)
        #x-axis
        m0 = MultipleLocator(1)
        ax.xaxis.set_major_locator(m0)
        m0 = MultipleLocator(0.5)
        ax.xaxis.set_minor_locator(m0)
        print "hello"
#        print data[0 
#        print legend1, legend2 
        plt.plot(redshifts1,data[0:len(redshifts1),0],label=legend1,color='blue',linestyle=linestyle2,linewidth=l)
        plt.plot(redshifts1,data[0:len(redshifts1),1],label=legend2,color='red',linestyle=linestyle2,linewidth=l)
        plt.plot(redshifts1,data[0:len(redshifts1),2],label=legend2,color='blue',linestyle='--',linewidth=l)
        plt.plot(redshifts1,data[0:len(redshifts1),3],label=legend2,color='red',linestyle='--',linewidth=l)
        lg = plt.legend(loc=2,prop={'size':fontsize})
        lg.draw_frame(False)


    elif (data.shape==(len(redshifts1),6) or data.shape==(len(redshifts2),6)): #plotting for comparing xfracs
        if legend1 == '' or legend2=='':
            print "WARNING: Specify legend for plotting three data sets"
        a = 0.5 #alpha value
        l = 2.0 #linewidth value
        plt.plot(redshifts1,np.log10(data[:,0]),color='red',label=legend1,alpha=a,linewidth=l)
        plt.plot(redshifts1,np.log10(data[:,1]),linestyle='--',color='red',alpha=a,linewidth=l+1)
        plt.plot(redshifts1,np.log10(data[:,2]),linestyle=':',color='red',alpha=a,linewidth=l+2)
        plt.plot(redshifts1,np.log10(data[:,3]),color='blue',label=legend2,alpha=a,linewidth=l)
        plt.plot(redshifts1,np.log10(data[:,4]),linestyle='--',color='blue',alpha=a,linewidth=l+1)
        plt.plot(redshifts1,np.log10(data[:,5]),linestyle=':',color='blue',alpha=a,linewidth=l+2)
        #plt.legend(loc=4,prop={'size':18})
        lg = plt.legend(loc=4,prop={'size':fontsize})
        lg.draw_frame(False)
        #yaxis
        m0 = MultipleLocator(2)
        ax.yaxis.set_major_locator(m0)
        m0 = MultipleLocator(1)
        ax.yaxis.set_minor_locator(m0)
        #x-axis
        m0 = MultipleLocator(1)
        ax.xaxis.set_major_locator(m0)
        m0 = MultipleLocator(0.5)
        ax.xaxis.set_minor_locator(m0)

    elif (data.shape==(len(redshifts1),4)):
        plt.plot(redshifts1,data[:,0],label=legend1)
        plt.plot(redshifts1,data[:,1],label=legend2)
        plt.plot(redshifts1,data[:,2],label=legend3)
        plt.plot(redshifts1,data[:,3],label=legend4)

        if legend1 == '' or legend2=='' or legend3=='' or legend4=='':
            print "WARNING: Specify legend for plotting three data sets"
        else:
            lg = plt.legend(loc=4,prop={'size':fontsize})
            lg.draw_frame(False)
 
    else:
        print "Data is in wrong format"
     #   plt.errorbar(rds,mn,yerr=er,fmt='o')
     #     plt.errorbar(rds,mn2,yerr=er2,fmt='o')
     #   plt.errorbar(rds,mn3,yerr=er3,fmt='o')
    plt.xlim(redshifts1[1],redshifts1[len(redshifts1)-1])
#    plt.title(title)
    plt.xlabel("Redshift",size=fontsize)
    plt.ylabel(ylabel,size=fontsize)
    #plt.legend(loc=2)i
    #plt.tight_layout()
    print "Saving file as " + pathp+title+".png ..."
    plt.savefig(pathp+title+".png")
    plt.close()

def plot_powerspectra(fr,data,title,i):
    fig,axes=plt.subplots(1,1,figsize=(10,6))
    print "plotting"
    if i!='null':
        plt.plot(np.log10(fr),np.log10(data))
        fig.text(0.2,0.85,"z = "+str(redshifts[i]),fontsize=fontsize,color='red')
    else:
        print "plotting redshifts:"
        for z in range(0,len(redshifts),10):
            print redshifts[z]
            plt.plot(np.log10(fr),np.log10(data[z,:]),label="z = " +str(redshifts[z]))
            lg=plt.legend(loc=2,prop={'size':fontsize-1},ncol=2)
            lg.draw_frame(False)
    plt.ylim(0,2.0)
    m0=MultipleLocator(0.5)
    axes.xaxis.set_major_locator(m0)

    m0=MultipleLocator(0.1)
    axes.xaxis.set_minor_locator(m0)

    #plt.ylim(0,1e6)
#    plt.title(title)
    plt.ylabel("log$_{10}$($\Delta_{21cm}$)  [mK] ",size=fontsize)
    plt.xlabel("log$_{10}$(k) [h Mpc$^{-1}$]",size=fontsize)
    plt.tight_layout()
    print "saing as: "+setup_dirs.plotsdir()+title+".png" 
    plt.savefig(setup_dirs.plotsdir()+title+".png")
    plt.close()

def plot_dbt_histogram(data,title,name,xlabel,nobins=500,type='log'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yscale('log',basey=10)
    weights = np.ones_like(data)/float(len(data))
    #logbins = np.zeros(nobins)
    #logbins[nobins/2:nobins] = np.logspace(-10,np.log10(100),nobins/2)
    #temp = -1*np.logspace(-10,np.log10(100),nobins/2)
    #logbins[0:nobins/2] = temp[::-1]
    #print logbins
    #print logbins

    #if(data.shape==(len(redshifts1),))
    ax.hist(data,weights=weights,bins=500,histtype='stepfilled',color='blue',edgecolor='blue')

#    elif(data.shape==(len(redshifts1),2)):
#        if legend1 == '' or legend2 =='':
#            print "WARNING: Specify legend for plotting two data sets"
#        ax.hist(data[:,0],weights=weights,bins=logbins,histtype='stepfilled',color='blue',edgecolor='red')
#        ax.hist(data[:,1],weights=weights,bins=logbins,histtype='stepfilled',color='red',edgecolor='red')
#

    plt.xlabel(xlabel,size=30)
    plt.xlim(-457,83)
    plt.ylabel("Probability",size=30)
    plt.savefig(setup_dirs.plotsdir()+name+".png")
    plt.close()
    print "histogram complete"

def plot_log_histogram(data,name,xlabel,label1='Stellar Sources',label2='Stellar sources and X-Rays',type='log',filedir='',nobins=90):
    fig = plt.figure()
    ax=plt.gca()
    ax = fig.add_subplot(111)
    ax.set_yscale('log',basey=10)
    if(data.shape==(250**3,)):
        weights = np.ones_like(data)/float(len(data))
        ax.hist(np.log10(data),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='red',edgecolor='red')
    elif(data.shape==(250**3,2)):
        print "Plotting 2 data sets..."+label1+" and "+label2
        weights = np.ones_like(data[:,0])/float(len(data))
        ax.hist(np.log10(data[:,0]),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='red',edgecolor='none',alpha=0.5,label=label1)
        ax.hist(np.log10(data[:,1]),weights=weights,bins=10**np.linspace(-1, 0.6, nobins),histtype='stepfilled',color='blue',edgecolor='none',alpha=0.5,label=label2)
        lg = plt.legend(loc=2,prop={'size':fontsize})
        lg.draw_frame(False)

#    lg=plt.legend(loc=2,prop={'size':14})
#    lg.draw_frame(False)
    plt.xlabel('log$_{10}$('+xlabel+')',size=fontsize)
    plt.ylabel("Probability",size=fontsize)
    plt.tight_layout()
    if filedir=='':filedir=setup_dirs.plotsdir()
    plt.savefig(filedir+name+".png")
    plt.close()
    print "histogram complete, saving as "+filedir+name+".png"

def plot_3_log_histograms(dataHI,dataHeI,dataHeII,title,name,xlabel,nobins=5000,type='log'):
    fig = plt.figure()
    #print np.amin(data3)
    dataHI = dataHI.flatten()
    dataHeI = dataHeI.flatten()
    dataHeII = dataHeII.flatten()
    ax = fig.add_subplot(111)
    ax.set_yscale('log',basey=10)
    if (type=='log'):
        weights = np.ones_like(dataHI)/len(dataHI)
        ax.hist(np.log10(dataHI),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',edgecolor='none',alpha=0.4,color='red',label='HII')
        ax.hist(np.log10(dataHeI),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',alpha=0.4,color='yellow',edgecolor='none',label='HeII')
        ax.hist(np.log10(dataHeII),weights=weights,bins=np.linspace(-16, 0, nobins),histtype='stepfilled',alpha=0.4,color='blue',edgecolor='none',label='HeIII')
        
    else:
        ax.hist(data,nobins,color='red',alpha=0.01)
        ax.hist(data,nobins,color='green',alpha=0.01)
    plt.title(title)
    plt.legend(loc=2)
    if type=='log':
        plt.xlabel('log10('+xlabel+')')
        plt.yscale=('log')
    else:
        plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.savefig(setup_dirs.plotsdir()+name+".png")
    plt.close()
    print "histogram complete"

#================================================================================

def xfrac_histogram():
    '''plot hisograms'''
    for i in range(1,len(redshifts)):#max(1,start),len(redshifts)):
        print redshifts[i]
        xfrac_filename = setup_dirs.path()+'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_HI = c2t.XfracFile(xfrac_filename).xi
        xfracHe1_filename = setup_dirs.path()+'xfrac3dHe1_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_He1 = c2t.XfracFile(xfracHe1_filename).xi
        xfracHe2_filename = setup_dirs.path()+'xfrac3dHe2_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_He2 = c2t.XfracFile(xfracHe2_filename).xi

        plot_3_log_histograms(xfrac_HI[len(xfrac_HI[:,0,0])/2,:,:],xfrac_He1[len(xfrac_HI[:,0,0])/2,:,:],xfrac_He2[len(xfrac_HI[:,0,0])/2,:,:],"Ionised Fraction, Redshift:" +str('%.3f' % redshifts[i]),"loghist_xfrac_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Ionised Fraction")


def temp_histrogram():
    for i in range(start,len(redshifts)):
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(temp_filename).temper
        data = data.flatten()
        print len(data), 250**3
        plot_log_histogram(data,"loghist_temper_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)","Redshift:" +str('%.3f' % redshifts[i]),"loghist_temper_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)")

def dbt_histogram():
    for i in range(len(redshifts)):
        print "Doing redshift " + str(redshifts[i]) + "..."
        dbt = IO.readbin("map_dbt_"+str('%.3f' % redshifts[i])+'.bin')
        data = dbt.flatten()
        print "Generating histogram..."
        plot_dbt_histogram(data,"Differential Brightness Temperature, Redshift:" +str('%.3f' % redshifts[i]),"loghist_dbt_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Differential Brightness Temperature(K)")

def powerSpec():
    ''' plotmap.plotdbt() '''
    for i in range(len(redshifts)):
        fr=IO.readoned("powerSpectraFrequencies_dbt_100b_"+str('%.3f' % redshifts[i]))
        data = IO.readoned("powerSpectra_100b_"+str('%.3f' % redshifts[i]))
        ps=data*fr**3./(4.*np.pi**2.)
        ps=np.sqrt(ps)
        plot_powerspectra(fr,ps,"ps_dbt_100b_notsquare"+str(i+10)+'_'+str(redshifts[i]),i)

def allPowerSpec():
    fr=IO.readoned("powerSpectraFrequencies_dbt_100b_"+str('%.3f' % redshifts[0]))
    ps=np.zeros(len(redshifts)*len(fr)).reshape(len(redshifts),len(fr))
    for i in range(0,len(redshifts),10):
        ps[i,:] = IO.readoned("powerSpectra_100b_"+str('%.3f' % redshifts[i]))
        #print len(data), len(fr)
    ps=ps*fr**3./(4.*np.pi**2.)
    ps=np.sqrt(ps)
    plot_powerspectra(fr,ps,"ps_dbt_100b_notsquare_all",'null')


def rms():
    rmsdata=np.zeros(len(redshifts))
    for i in range(len(redshifts)):
        print "opening ../generate_data/"+setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin"
        data=np.load("../generate_data/"+setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin")
        mean=np.mean(data)
        length=len(data)
        sum=0.0
        print  "summing..."
        for x in range(length):
            for y in range(length):
                for z in range(length):
                    rmsdata[i]=rmsdata[i]+(data[x,y,z]-mean)**2
        rmsdata[i]=np.sqrt(rmsdata[i]/length**3)
    plot_mean(rmsdata,"rms_500b","rms (mK)")

 


