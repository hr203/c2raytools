'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
#import IO
import sys
sys.path.append('../')
import redshifts as rs
sys.path.append('../')
import IO
sys.path.append('../../src/')
import c2raytools as c2t

redshifts = rs.read_redshifts("../red_ori2.dat")
print redshifts

c2t.set_verbose(True)
c2t.set_sim_constants(boxsize_cMpc = 244)

base_path = '/lustre/scratch/astro/hr203/RESULTS/'

def map_temp():
    for i in range(len(redshifts)):
        filename = 'data/map_temper_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        tfile = c2t.TemperFile(temp_filename)
        IO.writemap(tfile.temper,filename)
    print "Done"

def map_xfrac(id):
    for i in range(len(redshifts)):
        filename = 'data/map_xfrac'+id+'_'+str('%.3f' % redshifts[i])+'.dat'
        xfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        IO.writemap(xfile.xi,filename)
    print "Done"

def map_dbt():
    for i in range(len(redshifts)):
        filename = 'data/map_dbt_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if (i%2==0):
#            density_filename = '/lustre/scratch/astro/hr203/244Mpc_f2_8.2S_wstars/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'

        elif(i%2!=0):
#            density_filename = '/lustre/scratch/astro/hr203/244Mpc_f2_8.2S_wstars/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'


        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename)
        dfile = c2t.DensityFile(density_filename)

        dT_box = c2t.calc_dt_full(xfile, tfile, dfile, redshifts[i])

        IO.writemap(dT_box,filename)
    print "Done"

print map_temp()