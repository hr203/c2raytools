'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''

import numpy as np
import pylab as pl
#import IO
import sys
#sys.path.append('../')
#import redshifts as rs
sys.path.append('../')
import IO
sys.path.append('../../src/')
import c2raytools as c2t
sys.path.append('../')
c2t.set_sim_constants(244)
import setup_dirs

ss=setup_dirs.get_res()
start=0#48#0
redshifts = setup_dirs.read_redshifts()
print redshifts
print len(redshifts)
c2t.set_verbose(True)
c2t.set_sim_constants(boxsize_cMpc = 244)

def map_temp():
    for i in range(start,len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_temper_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        tfile = c2t.TemperFile(temp_filename)
        IO.writemap(tfile.temper,filename)
    print "Writen map to "+filename

def map_xfrac(id):
#    for i in range(len(redshifts)-6,len(redshifts)-4):
    print len(redshifts)
    for i in range(start,len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_xfrac'+id+'_'+str('%.3f' % redshifts[i])+'.dat'
        xfrac_filename = setup_dirs.path()+'xfrac3d'+id+'_'+str('%.3f' % redshifts[i]) + '.bin'
        xfile = c2t.XfracFile(xfrac_filename)
        IO.writemap(xfile.xi,filename)
    print "Written map to " + filename

def map_dbt():
    for i in range(start,len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_dbt_'+str('%.3f' % redshifts[i])+'.bin'
        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if i%2==0:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'

        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename).xi
        
        if ss!=' ':
            abu_he=0.074
            G_grav=6.674e-11
            pc=3.086e16
            Mpc=1e6*pc 
            H0=0.7*100.0*1.0e5/Mpc
            rho_crit_0=3.0*H0*H0/(8.0*np.pi*G_grav)
            Omega_B=0.044
            mu=(1.0-abu_he)+4.0*abu_he
            m_p=1.672661e-24

            dfile = np.ones(ss**3).reshape(ss,ss,ss)
            dfile = dfile*rho_crit_0*Omega_B/(mu*m_p)*(1.0+redshifts[i])**3
            dfile = dfile*1.23581719037e-35

        else:
            dfile = c2t.DensityFile(density_filename).cgs_density
        dT_box = c2t.calc_dt_full(xfile, tfile, dfile, redshifts[i])
#        IO.writemap(dT_box,filename)
        IO.writebin(dT_box,filename)
        print "Written map to "+filename


def map_dbt_hightemp():
    for i in range(start,len(redshifts)):
        filename = setup_dirs.resultsdir()+'map_dbt_hightemp'+str('%.3f' % redshifts[i])+'.bin'
        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if i%2==0:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'

        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename).xi

        if ss!=' ':
            abu_he=0.074
            G_grav=6.674e-11
            pc=3.086e16
            Mpc=1e6*pc
            H0=0.7*100.0*1.0e5/Mpc
            rho_crit_0=3.0*H0*H0/(8.0*np.pi*G_grav)
            Omega_B=0.044
            mu=(1.0-abu_he)+4.0*abu_he
            m_p=1.672661e-24

            dfile = np.ones(ss**3).reshape(ss,ss,ss)
            dfile = dfile*rho_crit_0*Omega_B/(mu*m_p)*(1.0+redshifts[i])**3
            dfile = dfile*1.23581719037e-35

        else:
            dfile = c2t.DensityFile(density_filename).cgs_density
        dT_box = c2t.calc_dt(xfile, dfile, redshifts[i])
#        IO.writemap(dT_box,filename)
        IO.writebin(dT_box,filename)
        print "Written map to "+filename



def smoothed_map_dbt(id=''):
    for i in range(start,len(redshifts)):
        filename = setup_dirs.resultsdir()+'smoothed_map_dbt_'+id+str('%.3f' % redshifts[i])+'.bin'
        dT_box=np.load(setup_dirs.resultsdir()+"map_dbt_"+id+str('%.3f' % redshifts[i])+".bin")
        wl=0.21*(1+redshifts[i])
        c=299792458.
        bw_r=wl/(2.e3) #radians
        bw = bw_r*3437.74677#arcminutesi
        print "Wavelength of 21-cm from redshift " + str(redshifts[i]) + " is " +str(wl)
        print "At redshift: "+ str(redshifts[i]) + " smoothing with a "+ str(bw) + " arc minute beam"
        dT_box = c2t.beam_convolve(dT_box[len(dT_box[:,0,0])/2,:,:],redshifts[i],244.,beam_w=bw)
        #IO.writebin(dT_box,filename)
        H0_SI = 2.2685e-18 
        Nc=bw_r*redshifts[i]*c*250/(H0_SI*244*3.085677581e+22)
        print "This corresponds to: " + str(Nc) + " cells"
#        print "Written map to "+filename




def map_dbt_lightcone():
    for i in range(0,len(redshifts)):
       # print 'filename: setup_dirs.resultsdir()+'map_dbt_lightcone_'+str('%.3f' % redshifts[i])+'.dat'

        filename = setup_dirs.resultsdir()+'map_dbt_lightcone_'+str('%.3f' % redshifts[i])+'.dat'
        temp_filename = setup_dirs.path()+'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        xfrac_filename = setup_dirs.path() +'xfrac3d_'+str('%.3f' % redshifts[i]) + '.bin'
        if i%2==0:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i]) + 'n_all.dat'
        else:
            density_filename='/research/prace/244Mpc_RT/244Mpc_f2_8.2pS_250/coarser_densities/' + str('%.3f' % redshifts[i-1]) + 'n_all.dat'

        tfile = c2t.TemperFile(temp_filename)
        xfile =  c2t.XfracFile(xfrac_filename).xi
        if ss==' ':
            dfile = c2t.DensityFile(density_filename).cgs_density
        else:
            dfile = np.ones(ss**3).reshape(ss,ss,ss)*1.981e-10*(1+redshifts[i])**3
        
        dT_box = c2t.calc_dt_full_lightcone(xfile, tfile, dfile, redshifts[i])
        IO.writemap(dT_box,filename)
        print "Written map to "+filename

