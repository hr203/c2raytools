'''
This file saves the mean ionization fraction for each redshift to data/xfrac.dat

For more information, see the full documentation at 
http://ttt.astro.su.se/~hjens/c2raytools/
'''
from scipy.integrate import quad
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
#import lightone_smooth

ss=setup_dirs.get_res()
start=27#48#0
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

def dbt_lightcone(ht=''):
    filenames=["" for x in range(len(redshifts))]
    for i in range(0,len(redshifts)):
        filenames[i]=setup_dirs.resultsdir()+"dbt"+ht+"/map_dbt_"+ht+str('%.3f' % redshifts[i])+".bin"
    output = setup_dirs.resultsdir()+'dbt_lightcone'+ht+'.bin'
    output2 = setup_dirs.resultsdir()+'dbt_lightcone_redshifts.bin'
    lightcone,z=c2t.make_lightcone(filenames,redshifts[len(redshifts)-1],redshifts[0],interpolation='linear')
    lightcone2=np.asarray(lightcone)
    z2=np.asarray(z)
    IO.writebin(lightcone2,output)
    IO.writebin(z2,output2)
    print "written lightcone to " + output

def smoothed_lightcone_dbt(id='',ht=''):
    infile = open(setup_dirs.resultsdir()+'dbt_lightcone'+ht+'.bin','rb')
    outfile = setup_dirs.resultsdir()+'dbt_lightcone_smooth'+ht+'.bin'
    zfile = open(setup_dirs.resultsdir()+'dbt_lightcone_redshifts.bin','rb')
    dT_box=np.load(infile)
    zs=np.load(zfile)

    dT_box3 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    dT_box2 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    for z in range(len(dT_box[1,1,:])-1,-1,-1):
        wl=0.21*(1+zs[z])
        c=299792.458
        omm=0.27
        oml=0.73
        omr=0.0
        omk=1.0-omm-oml-omr
        H0=70.

        def integrand(x,omm,oml,omr,omk): return 1.0/np.sqrt(omr*((1+x)**4)+omm*((1+x)**3)+oml+omk*((1+x)**2)) 
        def dc(z,omm,oml,omr,omk): return quad(integrand,0,z,args=(omm,oml,omr,omk))[0]#/(1.0+z)
        vec_dc = np.vectorize(dc)

        bw_r=wl/(2.e3) #radians
        bw = bw_r*3437.74677#arcminutesi
        print "Wavelength of 21-cm from redshift " + str(zs[z]) + " is " +str(wl) + "m"
        print "At redshift: "+ str(zs[z]) + " smoothing with a "+ str(bw) + " arc minute beam"

        rc = bw_r*c/H0*vec_dc(zs[z],omm,oml,omr,omk) #comoving Mpc
        Hz= H0*np.sqrt(omr*(1+z)**4+omm*(1+z)**3+oml+omk*(1+z)**2)

#        dnu=nu0*Hz*rc/(c*(1+zs[z])**2)
        dz=rc*Hz/c
        print "$\Delta$ z = " +str(dz)
        ncs=dz*zs[len(zs)-1]/len(zs)
        print str(ncs) + " cells in the z direction"
#        rc_h=rc/0.7 # comoving Mpc/h
#        print "This corresponds to "+str(rc_h)+"Mpc/h on the sky"

        dT_box2[:,:,z] = c2t.beam_convolve(dT_box[:,:,z],zs[z],244.,beam_w=bw)

        if z>ncs and z+ncs<zs[len(zs)-1]:
            for x in range(len(dT_box2)):
                for y in range(len(dT_box2)):
                    dT_box3[x,y,z]=np.mean(dT_box2[x,y,z-ncs/2.:z+ncs/2.])
        else:
            print "..."
            dT_box3[:,:,z] = dT_box2[:,:,z]

        IO.writebin(dT_box3[:,:,z],setup_dirs.resultsdir()+'smoothed_map_dbt_'+ht+str('%.3f'%zs[z])+'.bin')

    IO.writebin(dT_box3,outfile)
    print "Written map to "+outfile

def smoothed_lightcone_dbt(id='',ht=''):
    infile = open(setup_dirs.resultsdir()+'dbt_lightcone'+ht+'.bin','rb')
    outfile = setup_dirs.resultsdir()+'dbt_lightcone_smooth'+ht+'.bin'
    zfile = open(setup_dirs.resultsdir()+'dbt_lightcone_redshifts.bin','rb')
    dT_box=np.load(infile)
    log=open("stats.dat","w")
    zs=np.load(zfile)

    dT_box3 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    dT_box2 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    for z in range(len(dT_box[1,1,:])-1,-1,-1):
        wl=0.21*(1+zs[z])
        c=299792.458
        omm=0.27
        oml=0.73
        omr=0.0
        omk=1.0-omm-oml-omr
        H0=70.

        def integrand(x,omm,oml,omr,omk): return 1.0/np.sqrt(omr*((1+x)**4)+omm*((1+x)**3)+oml+omk*((1+x)**2))
        def dc(z,omm,oml,omr,omk): return quad(integrand,0,z,args=(omm,oml,omr,omk))[0]#/(1.0+z)
        vec_dc = np.vectorize(dc)

        bw_r=wl/(2.e3) #radians
        bw = bw_r*3437.74677#arcminutesi
        log.write( "Wavelength of 21-cm from redshift " + str(zs[z]) + " is " +str(wl) + "m\n")
        log.write("At redshift: "+ str(zs[z]) + " smoothing with a "+ str(bw) + " arc minute beam.\n")

        rc = bw_r*c/H0*vec_dc(zs[z],omm,oml,omr,omk) #comoving Mpc
        Hz= H0*np.sqrt(omr*(1+zs[z])**4+omm*(1+zs[z])**3+oml+omk*(1+zs[z])**2)

        log.write("rc = " + str(rc)+'\n')
#        dnu=nu0*Hz*rc/(c*(1+zs[z])**2)
        dz=rc*Hz/c
        log.write("$\Delta$ z = " +str(dz)+'\n')
        ncs=rc*250.*0.7/244.
        log.write(str(ncs) + " cells in the z direction\n")
        log.write('\n')
#        rc_h=rc/0.7 # comoving Mpc/h
#        print "This corresponds to "+str(rc_h)+"Mpc/h on the sky"

        dT_box2[:,:,z] = c2t.beam_convolve(dT_box[:,:,z],zs[z],244.,beam_w=bw)

        if z>ncs and z+ncs<zs[len(zs)-1]:
            for x in range(len(dT_box2)):
                for y in range(len(dT_box2)):
                    dT_box3[x,y,z]=np.mean(dT_box2[x,y,z-ncs/2.:z+ncs/2.])
        else:
            print "..."
            dT_box3[:,:,z] = dT_box2[:,:,z]

        IO.writebin(dT_box3[:,:,z],setup_dirs.resultsdir()+'smoothed_map_dbt_'+ht+str('%.3f'%zs[z])+'.bin')

    IO.writebin(dT_box3,outfile)
    print "Written map to "+outfile

def smoothed_coeval_box(ht=''):
    log2=open("coeval_stats.dat","w")
    for z in range(len(redshifts)):
        infile = 'map_dbt_'+ht+str('%0.3f'%redshifts[z])+'.bin'
        dT_box=IO.readbin(infile)
        dT_box3 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
        dT_box2 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
        wl=0.21*(1+redshifts[z])
        c=299792.458
        omm=0.27
        oml=0.73
        omr=0.0
        omk=1.0-omm-oml-omr 
        H0=70.
        def integrand(x,omm,oml,omr,omk): return 1.0/np.sqrt(omr*((1+x)**4)+omm*((1+x)**3)+oml+omk*((1+x)**2))
        def dc(z,omm,oml,omr,omk): return quad(integrand,0,z,args=(omm,oml,omr,omk))[0]#/(1.0+z)
        vec_dc = np.vectorize(dc)
        bw_r=wl/(2.e3) #radians
        bw = bw_r*3437.74677#arcminutesi
        log2.write( "Wavelength of 21-cm from redshift " + str(redshifts[z]) + " is " +str(wl) + "m\n")
        log2.write("At redshift: "+ str(redshifts[z]) + " smoothing with a "+ str(bw) + " arc minute beam.\n")
    
        rc = bw_r*c/H0*vec_dc(redshifts[z],omm,oml,omr,omk) #comoving Mpc
        Hz= H0*np.sqrt(omr*(1+redshifts[z])**4+omm*(1+redshifts[z])**3+oml+omk*(1+redshifts[z])**2)
        log2.write("rc = " + str(rc)+'\n')
        dz=rc*Hz/c
        log2.write("$\Delta$ z = " +str(dz)+'\n')
        ncs=rc*250.*0.7/244.
        log2.write(str(ncs) + " cells in the z direction\n")
        log2.write('\n')
        for k in range(len(dT_box2[1,1,:])):
            dT_box2[:,:,k] = c2t.beam_convolve(dT_box[:,:,k],redshifts[z],244.,beam_w=bw)
        for i in range(0):#len(dT_box2[:,1,1])):
            for j in range(0):#len(dT_box2[1,:,1])):
                for k  in range(0):#len(dT_box2[1,1,:])):
                    L=len(dT_box2[1,1,:])
                    if k-ncs/2>0 and k+ncs/2<L:
                       dT_box3[i,j,k]=np.mean(dT_box2[i,j,k-ncs/2.:k+ncs/2.])
                       if k+ncs/2.-k+ncs/2.!=ncs:
                            print "Wrong smoothing width for tophat function!!"
                    elif (k-ncs)/2<0:
                       if k+ncs/2.+k-ncs/2!=ncs:
                            print "Wrong smoothing width for tophat function!!"
                       dT_box3[i,j,k]=(np.sum(dT_box2[i,j,0:(k+ncs/2.)])+np.sum(dT_box2[i,j,L-ncs/2+k:L]))/ncs
                    elif (k+ncs)/2>10:
                       if L-k+ncs/2+L-k+ncs/2!=ncs:
                            print "Wrong smoothing width for tophat function!!"
                       dT_box3[i,j,k]=(np.sum(dT_box2[i,j,k-ncs/2:L])+np.sum(dT_box2[i,j,L-k+ncs/2:L]))/ncs

        IO.writebin(dT_box2,setup_dirs.resultsdir()+'smooth_coevalmap_dbt_'+ht+str('%.3f'%redshifts[z])+'.bin')
        print np.mean(dT_box2)
        print "Written map"


def smoothed_lightcone_dbt2(id='',ht=''):
    print "hello" 
    infile = open(setup_dirs.resultsdir()+'dbt_lightcone'+ht+'.bin','rb')
    outfile = setup_dirs.resultsdir()+'dbt_lightcone_smooth'+ht+'.bin'
    zfile = open(setup_dirs.resultsdir()+'dbt_lightcone_redshifts.bin','rb')
    dT_box=np.load(infile)
    log=open("stats.dat",'w')
    zs=np.load(zfile)

    dT_box3 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    dT_box2 = np.zeros((len(dT_box[:,1,1])/2,len(dT_box[1,:,1])/2,len(dT_box[1,1,:])))
    for z in range(len(dT_box[1,1,:])):
        wl=0.21*(1+zs[z])
        c=299792.458
        omm=0.27
        oml=0.73
        omr=0.0
        omk=1.0-omm-oml-omr
        H0=70.
        print 1
        def integrand(x,omm,oml,omr,omk): return 1.0/np.sqrt(omr*((1+x)**4)+omm*((1+x)**3)+oml+omk*((1+x)**2))
        def dc(z,omm,oml,omr,omk): return quad(integrand,0,z,args=(omm,oml,omr,omk))[0]#/(1.0+z)
        vec_dc = np.vectorize(dc)

        print 2
        bw_r=wl/(2.e3) #radians
        bw = bw_r*3437.74677#arcminutesi
        stats.write("Wavelength of 21-cm from redshift " + str(zs[z]) + " is " +str(wl) + "m")
        stats.write("At redshift: "+ str(zs[z]) + " smoothing with a "+ str(bw_r) + " radian beam")
        stats.write("At redshift: "+ str(zs[z]) + " smoothing with a "+ str(bw) + " arc minute beam")
        stats.write("\n")

        print 3
        dT_box2[:,:,z] = c2t.beam_convolve(dT_box[:,:,z],zs[z],244.,beam_w=bw)

        print 4
        dT_box3=lightcone_smooth.frequency_direction_same(dT_box2,zs[0],244.,2)
        print 5
        IO.writebin(dT_box3[:,:,z],setup_dirs.resultsdir()+'smoothed_map_dbt_'+ht+str('%.3f'%zs[z])+'.bin')

    IO.writebin(dT_box3,outfile)
    print "Written map to "+outfile




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

