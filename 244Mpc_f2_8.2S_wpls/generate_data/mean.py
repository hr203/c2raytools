'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
import power_spectrum as ps
import sys
sys.path.append('../')
import setup_dirs
#import redshifts as rs
sys.path.append('../../src/')
import c2raytools as c2t
sys.path.append('../')
import IO 

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

def mean_xfrac(id):
    file = open(setup_dirs.resultsdir()+'mean_xfrac'+id+'.dat', 'w')
    for i in range(len(redshifts)):
        xfrac_filename = setup_dirs.path()+'xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        file.write(str(np.mean(xfile.xi)) + '\n') 
    file.close()
    print "Written mean xfrac to "+setup_dirs.resultsdir()+'mean_xfrac'+id+'.dat'

def mean_dbt():
    file = open(setup_dirs.resultsdir()+'mean_dbt.dat','w')
    print len(redshifts)
    for i in range(len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_dbt_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if i%2==0:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'
        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename)
        dfile = c2t.DensityFile(density_filename)

        dT_box = c2t.calc_dt_full(xfile, tfile, dfile, redshifts[i]) #returned in micro kelvin! update plots
        file.write(str(np.mean(dT_box))+'\n')
    print "Written mean dbt to " + setup_dirs.resultsdir()+'mean_dbt.dat'

def power_spectrum():
    for i in range(len(redshifts)):
        print "Doing redshift: " + str(redshifts[i])
        data=IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
        powerspec=ps.power_spectrum_1d(data,500)
        IO.write2data(powerspec[0],powerspec[1],setup_dirs.resultsdir()+'/powerSpectra_'+str('%.3f' % redshifts[i])+'.dat',setup_dirs.resutlsdir()+'/powerSpectraFrequencies_dbt_'+str('%.3f' % redshifts[i])+'.dat')
    #    IO.writedata(powerspec[1],'data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')
 
