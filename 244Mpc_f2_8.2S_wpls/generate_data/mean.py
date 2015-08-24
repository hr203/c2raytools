'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
import sys
sys.path.append('../')
import redshifts as rs
sys.path.append('../../src/')
import c2raytools as c2t

redshifts = rs.read_redshifts("../red_ori2.dat")
c2t.set_verbose(True)
c2t.set_sim_constants(boxsize_cMpc = 244)

def mean_temp():
    file = open('data/mean_temp.dat', 'w')
    base_path = '/lustre/scratch/astro/hr203/RESULTS/'
    for i in range(len(redshifts)):
        filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(filename)
        file.write(str(np.mean(data.temper)) + '\n')
    file.close()
    print "Done"

def mean_xfrac(id):
    file = open('data/mean_xfrac'+id+'.dat', 'w')
    base_path = '/lustre/scratch/astro/hr203/RESULTS/'

    for i in range(len(redshifts)):
        xfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        file.write(str(np.mean(xfile.xi)) + '\n') 
    file.close()
    print "Done"
