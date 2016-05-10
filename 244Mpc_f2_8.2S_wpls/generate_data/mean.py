'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
import scipy
import power_spectrum as ps
import sys
sys.path.append('../')
import setup_dirs
#import redshifts as rs
sys.path.append('../../src/')
import c2raytools as c2t
c2t.set_sim_constants(244)
sys.path.append('../')
import IO 

ss=setup_dirs.get_res()
start = 0#48#27  
redshifts = setup_dirs.read_redshifts()
c2t.set_verbose(True)
c2t.set_sim_constants(boxsize_cMpc = 244)

def mean_temp():
    file = open(setup_dirs.resultsdir()+'mean_temp.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        file.write(str(np.mean(data.temper)) + '\n')
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'mean_temp.dat'

def mean_dbt():
    file=open(setup_dirs.resultsdir()+'mean_dbt.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_dbt_'+str('%.3f' % redshifts[i]) + '.bin'
        data = np.load(filename)
        file.write(str(np.mean(data)) + '\n')
        print np.mean(data) 
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'mean_dbt.dat'

def mean_dbt_smooth():
    file=open(setup_dirs.resultsdir()+'smooth_mean_dbt.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.resultsdir()+'smoothed_map_dbt_'+str('%.3f' % redshifts[i]) + '.bin'
        data = np.load(filename)
        file.write(str(np.mean(data)) + '\n')
        print np.mean(data)
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'smooth_mean_dbt.dat'


def median_temp():
    file = open(setup_dirs.resultsdir()+'median_temp.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        file.write(str(np.median(data.temper)) + '\n')
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'median_temp.dat'

def mean_temp_igm():
    file = open(setup_dirs.resultsdir()+'mean_temp_igm.dat', 'w')
    file2= open(setup_dirs.resultsdir()+'fraction_ionised.dat','w')
    for i in range(len(redshifts)):
        filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        xfrac_filename = setup_dirs.path()+'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        igm_data=[]
        ionisedpoints=0
        for i in range(len(data.temper[1,:,1])):
            for j in range(len(data.temper[1,:,1])):
                for k in range(len(data.temper[1,:,1])):
                    if xfile.xi[i,j,k] < 1e3:
                        igm_data.append(data.temper[i,j,k])
                    else:
                        ionisedpoints=ionisedpoints+1
        file.write(str(np.mean(igm_data)) + '\n')
        file.write(str(ionisedpoints/250**3))
        
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'mean_temp_igm.dat'

def min_temp():
    file = open(setup_dirs.resultsdir()+'min_temp.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        file.write(str(np.amin(data.temper)) + '\n')
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'min_temp.dat'


def mean_xfrac(id):
    file = open(setup_dirs.resultsdir()+'mean_xfrac'+id+'.dat', 'w')
    for i in range(len(redshifts)):
        xfrac_filename = setup_dirs.path()+'xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        file.write(str(np.mean(xfile.xi)) + '\n') 
    file.close()
    print "Written mean xfrac to "+setup_dirs.resultsdir()+'mean_xfrac'+id+'.dat'

def xfrac_power_spectrum(id):
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=c2t.XfracFile(setup_dirs.path()+"xfrac3d"+id+'_'+str('%.3f' % redshifts[i])+".bin").xi
        powerspec=ps.power_spectrum_1d(data,100)
        IO.write2data(powerspec[0],powerspec[1],setup_dirs.resultsdir()+'/powerSpectra_xfrac_'+id+'_'+str('%.3f' % redshifts[i])+'.dat',setup_dirs.resultsdir()+'/powerSpectraFrequencies_dbt_100b_'+str('%.3f' % redshifts[i])+'.dat')

def temp_power_spectrum():
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=c2t.TemperFile(setup_dirs.path()+"Temper3D_"+str('%.3f' % redshifts[i])+".bin").temper
        powerspec=ps.power_spectrum_1d(data,100)
        IO.write2data(powerspec[0],powerspec[1],setup_dirs.resultsdir()+'/powerSpectra_temp_'+str('%.3f' % redshifts[i])+'.dat',setup_dirs.resultsdir()+'/powerSpectraFrequencies_dbt_100b_'+str('%.3f' % redshifts[i])+'.dat')

#ean_dbt():
#    file = open(setup_dirs.resultsdir()+'mean_dbt.dat','w')
#    for i in range(len(redshifts)):
#        filename = setup_dirs.resultsdir()+'map_dbt_'+str('%.3f' % redshifts[i])+'.dat'
#        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
##        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
#        if i%2==0:
 #           density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
 #       else:
 #           density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'
 #       tfile = c2t.TemperFile(temp_filename)
 #       xfile =  c2t.XfracFile(xfrac_filename).xi
 #       if ss!=' ':
 #           dfile = np.ones(ss**3).reshape(ss,ss,ss)*1.981e-10*(1+redshifts[i])**3
 #       else:
 #           dfile = c2t.DensityFile(density_filename).cgs_density

 #       dT_box = c2t.calc_dt_full(xfile, tfile, dfile, redshifts[i]) #returned in micro kelvin! update plots
 #       file.write(str(np.mean(dT_box))+'\n')
 #   print "Written mean dbt to " + setup_dirs.resultsdir()+'mean_dbt.dat'

def mean_dbt_hightemp():
    file = open(setup_dirs.resultsdir()+'mean_dbt_hightemp.dat','w')
    for i in range(len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_dbt_hightemp'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if i%2==0:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'
        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename).xi
        if ss!=' ':
            dfile = np.ones(ss**3).reshape(ss,ss,ss)*1.981e-10*(1+redshifts[i])**3
        else:
            dfile = c2t.DensityFile(density_filename).cgs_density

        dT_box = c2t.calc_dt(xfile, dfile, redshifts[i]) #returned in micro kelvin! update plots
        file.write(str(np.mean(dT_box))+'\n')
    print "Written mean dbt to " + setup_dirs.resultsdir()+'mean_dbt_hightemp.dat'

def median_temp_hightemp():
    file = open(setup_dirs.resultsdir()+'median_temp_hightemp.dat', 'w')
    for i in range(len(redshifts)):
        filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        file.write(str(np.median(data.temper)) + '\n')
    file.close()
    print "Written mean Temperature to "+setup_dirs.resultsdir()+'median_temp_hightemp.dat'

def power_spectrum(id=''):
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=np.load(setup_dirs.resultsdir()+"map_dbt_"+id+str('%.3f' % redshifts[i])+".bin")
        powerspec=ps.power_spectrum_1d(data,100)
        IO.write2data(powerspec[0],powerspec[1],setup_dirs.resultsdir()+'/powerSpectra_100b_'+id+str('%.3f' % redshifts[i])+'.dat',setup_dirs.resultsdir()+'/powerSpectraFrequencies_dbt_100b_'+str('%.3f' % redshifts[i])+'.dat')
    #    IO.writedata(powerspec[1],'data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')

def rms():
    file = open(setup_dirs.resultsdir()+'rms.dat','w')
#    file2 = open(setup_dirs.resultsdir()+'rms_density.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        print "opening "+setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin"
        data=np.load(setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin")
        mean=np.mean(data)
        print "mean dbt: " + str(mean)
        length=len(data)
        rmsdata=0.0
        print  "summing..."
#        for x in range(length):
#            for y in range(length):
#                for z in range(length):
#                    rmsdata=rmsdata+(data[x,y,z]-mean)**2
#                    #rmsdata=rmsdata+(data[x,y,z])**2
        rmsdata=np.sum((data-mean)**2)
        rmsdata=np.sqrt(rmsdata/length**3)
        print rmsdata, np.sqrt(np.var(data)) 
        #rmsdata=(rmsdata/length**3)
        file.write(str(rmsdata)+'\n')
        print(str(rmsdata)+'\n')
    print "Written rms to " + setup_dirs.resultsdir()+" rms.dat"

def smooth_rms(id=''):
    file = open(setup_dirs.resultsdir()+'smooth_rms'+id+'.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        print "opening "+setup_dirs.resultsdir()+"smoothed_map_dbt_"+id+str('%.3f' % redshifts[i])+".bin"
        data=np.load(setup_dirs.resultsdir()+"smoothed_map_dbt_"+id+str('%.3f' % redshifts[i])+".bin")
        mean=np.mean(data)
        print "mean dbt: " + str(mean)
        length=len(data)
        rmsdata=0.0
        print  "summing..."
        for x in range(length):
            for y in range(length):
#                for z in range(length):
		#rmsdata=rmsdata+(data[x,y]-mean)**2
                rmsdata=rmsdata+(data[x,y])**2
        rmsdata=np.sqrt(rmsdata/length**2)
        #rmsdata=(rmsdata/length**2)
        file.write(str(rmsdata)+'\n')
        print(str(rmsdata)+'\n')
    print "Written rms to " + setup_dirs.resultsdir()+"smooth_rms"+id+".dat"
 
def skewness():
    file = open(setup_dirs.resultsdir()+'skewness.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        print setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin"
        data=np.load(setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin")
        data=data.flatten()
        print scipy.stats.skew(data), c2t.statistics.skewness(data)
        #skewness=scipy.stats.skew(data)        
        skewness=c2t.statistics.skewness(data)

        file.write(str(skewness)+'\n')
    print "Written skewness to " + setup_dirs.resultsdir()+"skewness.dat"

def kurtosis():
    file = open(setup_dirs.resultsdir()+'kurtosis.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=np.load(setup_dirs.resultsdir()+"map_dbt_"+str('%.3f' % redshifts[i])+".bin")
        data=data.flatten()
#        kurtosis=scipy.stats.kurtosis(data)
        kurtosis=c2t.statistics.kurtosis(data)
        file.write(str(kurtosis)+'\n')
    print "Written kurtosis to " + setup_dirs.resultsdir()+"kurtosis.dat"



def smooth_skewness(id=''):
    file = open(setup_dirs.resultsdir()+'smooth_skewness'+id+'.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        print setup_dirs.resultsdir()+"smoothed_map_dbt_"+id+str('%.3f' % redshifts[i])+".bin"
        data=np.load(setup_dirs.resultsdir()+"smoothed_map_dbt_"+id+str('%.3f' % redshifts[i])+".bin")
        data=data.flatten()
        skewness=scipy.stats.skew(data)
#        skewness=c2t.statistics.skewness(data)
        file.write(str(skewness)+'\n')
    print "Written skewness to " + setup_dirs.resultsdir()+"smooth_skewness"+id+".dat" 

def smooth_kurtosis(id=''):
    file = open(setup_dirs.resultsdir()+'smooth_kurtosis'+id+'.dat','w')
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=np.load(setup_dirs.resultsdir()+"smoothed_map_dbt_"+id+str('%.3f' % redshifts[i])+".bin")
        data=data.flatten()
        kurtosis=scipy.stats.kurtosis(data)
#        kurtosis=c2t.statistics.kurtosis(data)
        file.write(str(kurtosis)+'\n')
    print "Written kurtosis to " + setup_dirs.resultsdir()+"smooth_kurtosis"+id+".dat"
