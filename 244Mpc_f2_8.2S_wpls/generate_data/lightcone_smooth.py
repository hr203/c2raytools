import numpy as np
import sys
sys.path.append('../../src/')
import c2raytools as c2t
import c2raytools.cosmology as cm
import c2raytools.helper_functions as hf

def angular_direction(lightcone, sigma=1.0):
	nx, ny, nz = lightcone.shape
	smooth = np.zeros((nx,ny,nz))
	for k in xrange(nz):
		smooth[:,:,k] = c2t.smooth_gauss(lightcone[:,:,k], sigma=sigma)
	return smooth

def frequency_direction_(lightcone, z_low, box_size_mpc, dnu):

	#Figure out the comoving distances of the slices in input lightcone
	cell_size = box_size_mpc/lightcone.shape[0]
	distances = cm.z_to_cdist(z_low) + np.arange(lightcone.shape[2])*cell_size
	input_redshifts = cm.cdist_to_z(distances)
	input_frequencies = cm.z_to_nu(input_redshifts)
	nu1 = input_frequencies[0]
	nu2 = input_frequencies[-1]
	output_frequencies = np.arange(nu1, nu2, -dnu) - dnu/2.
	
	output_lightcone = np.zeros((lightcone.shape[0], lightcone.shape[1], \
                                 len(output_frequencies)))
    
	#Bin in frequencies by smoothing and indexing    
	for i in xrange(output_lightcone.shape[2]-1):
		max_cell_size = cm.nu_to_cdist(output_frequencies[i])-cm.nu_to_cdist(output_frequencies[i+1])
		smooth_scale = np.round(max_cell_size/cell_size)
		if smooth_scale < 1:
			smooth_scale = 1
		hf.print_msg('Smooth along LoS with scale %f' % smooth_scale)
		tophat3d = np.ones((1,1,smooth_scale))
		tophat3d /= np.sum(tophat3d)
		lightcone_smoothed = hf.fftconvolve(lightcone, tophat3d)
		nu = output_frequencies[i]
		idx = hf.find_idx(input_frequencies, nu)
		output_lightcone[:,:,i] = lightcone_smoothed[:,:,idx]

	return output_lightcone, output_frequencies

def frequency_direction(lightcone, z_low, box_size_mpc, dnu):
	cell_size = box_size_mpc/lightcone.shape[0]
	distances = cm.z_to_cdist(z_low) + np.arange(lightcone.shape[2])*cell_size
	input_redshifts = cm.cdist_to_z(distances)
	input_frequencies = cm.z_to_nu(input_redshifts)
	nu1 = input_frequencies[0]
	nu2 = input_frequencies[-1]
	output_frequencies = np.arange(nu1, nu2, -dnu) - dnu/2.
	
	output_lightcone = np.zeros((lightcone.shape[0], lightcone.shape[1], \
                                 len(output_frequencies)))

	for i in xrange(output_lightcone.shape[2]):
		nu = output_frequencies[i]
		z_out_low  = cm.nu_to_z(nu+dnu/2.)
		z_out_high = cm.nu_to_z(nu-dnu/2.)
		if i==0:idx_low = np.ceil(hf.find_idx(input_redshifts, z_out_low))
		else:	idx_low = idx_high
		idx_high = np.ceil(hf.find_idx(input_redshifts, z_out_high))
		output_lightcone[:,:,i] = np.mean(lightcone[:,:,idx_low:idx_high], axis=2)

	return output_lightcone, output_frequencies

#def sigma_to_dnu(sigma):
	
def frequency_direction_same(lightcone, z_low, box_size_mpc, baseline):
    cell_size = 1.0*box_size_mpc/lightcone.shape[0]   
    distances = cm.z_to_cdist(z_low) + np.arange(lightcone.shape[2])*cell_size
    input_redshifts = cm.cdist_to_z(distances)
    #output_frequencies = cm.z_to_nu(input_redshifts)
    output_dtheta  = (1+input_redshifts)*21e-5/baseline
    output_ang_res = output_dtheta*cm.z_to_cdist(input_redshifts)
    output_dz      = output_ang_res/c2t.const.c
    for i in xrange(len(output_dz)):
        output_dz[i] = output_dz[i] * hubble_parameter(input_redshifts[i])

    output_lightcone = np.zeros(lightcone.shape)
    for i in xrange(output_lightcone.shape[2]):
        z_out_low  = input_redshifts[i]-output_dz[i]/2
        z_out_high = input_redshifts[i]+output_dz[i]/2
        if i==0:idx_low = np.ceil(hf.find_idx(input_redshifts, z_out_low))
        else:    idx_low = idx_high
        idx_high = np.ceil(hf.find_idx(input_redshifts, z_out_high))
        output_lightcone[:,:,i] = np.mean(lightcone[:,:,idx_low:idx_high+1], axis=2)

    return output_lightcone, input_redshifts


def hubble_parameter(z):
    part = np.sqrt(c2t.const.Omega0*(1.+z)**3+c2t.const.lam)
    return c2t.const.H0 * part
