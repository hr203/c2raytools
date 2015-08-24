#import plotmap
#import plotmean
#import power_spectrum as ps
import sys
sys.path.append('../')
import redshifts as rs
import IO
import plot_1D

redshifts = rs.read_redshifts("../red_ori2.dat")
mesh=250


def plotmap():
    ''' plotmap.plotdbt() '''
    for i in range(len(redshifts)-3):
        fr = IO.readoned("powerSpectraFrequencies_dbt_"+str('%.3f' % redshifts[i]))
        data = IO.readoned("powerSpectra_"+str('%.3f' % redshifts[i]))
        print data
        #print ' ---------------------------------------------'
        #print fr
        print '=============================================='
        #print data
        #print fr
        plot_1D.plot_powerspectra(fr,data,"ps_dbt_more_bins_"+str(redshifts[i]))


def histogram():
    '''plot hisograms'''
    for i in range(len(redshifts)-3):
        data = IO.readoned("map_temper_"+str('%.3f' % redshifts[i]))
        plot_1D.plot_histogram(data,"Temperature, Redshift:" +str('%.3f' % redshifts[i]),"loghist_temper_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)",nobins=200) 

#create powerspectra datafile from dbt
#for i in range(len(redshifts)):
#    data=IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
#    powerspec=ps.power_spectrum_1d(data)
#    IO.write2data(powerspec[0],powerspec[1],'data/powerSpectra_'+str('%.3f' % redshifts[i])+'.dat','data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')
#    #IO.writedata(powerspec[1],'data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')

histogram()
