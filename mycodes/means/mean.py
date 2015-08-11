'''
This file illustrates some of the most basic usage for the c2raytools package.
It reads some data files and prints and plots some statistics. 
To run this script, you must have c2raytools installed in a location where
Python can find it. To do this, see the Readme file, or http://ttt.astro.su.se/~hjens/c2raytools/
Second, you must modify the path names to point to files that you actually have access to.

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

#import c2raytools as c2t
import numpy as np
import pylab as pl
#import redshifts as rs
import sys
sys.path.append('../')
import redshifts as rs
sys.path.append('../src/c2raytools/')
import c2raytools as c2t



def mean_xfrac():

    redshifts = rs.read_redshifts("../red_ori2.dat")
    redshifts2 = rs.read_redshifts("../red_ori.dat")
    meanxfrac = np.zeros(len(redshifts))
    
    file = open('data/xfrac_stars.dat', 'w')
    file.write("Stars only: Mean Xfrac")

    file_pls = open('data/xfrac_pls.dat','w')
    file_pls.write("Power laws: Mean xfrac")

    file_quasars = open('data/xfrac_quasars.dat', 'w')
    file_quasars.write("Quasars: Mean Xfrac")

    file_quasars_pls = open('data/xfrac_quasars_pls.dat', 'w')
    file_quasars_pls.write("Power laws and Quasars: Mean Xfrac")

    for i in range(len(redshifts)):

        #Some path and file names. Modify these as needed.
        base_path = '/lustre/scratch/astro/hr203/RESULTS/'

        #for i in range(noRedshifts):  
        density_filename = '/lustre/scratch/astro/hr203/244Mpc_f2_8.2S_H250_wstars/coarser_densities/'

        xfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wstars/xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        qxfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wquasars/xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if (i%2 ==0):
            pxfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wpls/xfrac3d_'+str('%.3f' % redshifts[i/2]) + '.bin'
        qpxfrac_filename = base_path + '/244Mpc_f2_8.2S_H250_wquasars_wpls/xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'

        #Enable the printing of various messages
        c2t.set_verbose(True)

        #We are using the 114/h Mpc simulation box, so set all the proper conversion factors.
        #Be sure to always set this before loading anything. Otherwise, c2raytools will
        #not know how to convert densities and velocities to physical units!
        #c2t.set_sim_constants(boxsize_cMpc = 50)

        #Read an ionized fractions file and store it as an XfracFile object
        xfile = c2t.XfracFile(xfrac_filename)
        if (i%2==0):
            xfile_pls = c2t.XfracFile(pxfrac_filename)
        xfile_qs = c2t.XfracFile(pxfrac_filename)
        xfile_qs_pls = c2t.XfracFile(qpxfrac_filename)

       #The most important property of an XfracFile object is xi, which
       #is a numpy array containing the ionized fraction
       #print 'The ionized fraction in point (10,223,45) is: ', xfile.xi[10,223,45]
       #print 'The volume-averaged mean ionized fraction is: ', xfile.xi.mean()
        file.write(xfile.xi.mean())
        if (i%2):
            file_pls.write(xfile_pls.xi.mean())
        file.write(xfilei_qs.xi.mean())
        file.write(xfile_qs_pls.xi.mean())
        return "Complete"

print mean_xfrac()
