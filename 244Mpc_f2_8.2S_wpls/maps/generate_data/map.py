'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
import sys
sys.path.append('../../')
import redshifts as rs
sys.path.append('../../../src/')
import c2raytools as c2t

redshifts = rs.read_redshifts("../../red_ori2.dat")
c2t.set_verbose(True)
c2t.set_sim_constants(boxsize_cMpc = 244)

base_path = '/lustre/scratch/astro/hr203/RESULTS/'

def printf(data,filename):
    file = open(filename, 'w')
    l = len(data[:,1,1])      
    for i in range(l):
        for j in range(l):
            for k in range(l):
                file.write(str(data[i,j,k])+'\n')
    file.close()

def mean_temp():
   # file = open('data/map_temp.dat', 'w')
    for i in range(len(redshifts)):
        filename = 'data/map_temper_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        tfile = c2t.TemperFile(temp_filename)
        printf(tfile.temper,filename)
    return "Done"

def mean_xfrac(id):
    for i in range(len(redshifts)):
        filename = 'data/map_xfrac'+id+'_'+str('%.3f' % redshifts[i])+'.dat'
        xfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        printf(xfile.xi,filename)
#        l = len(xfile.xi[:,1,1])
#        for i in range(l):
#            for j in range(l):
#                for k in range(l):                    
#                    file.write(str(xfile.xi[i,j,k])+'\n')
                     
#        file.close()
    return "Done"

print mean_temp()
#print mean_xfrac('')
print mean_xfrac('He1')
#print mean_xfrac('He2')
