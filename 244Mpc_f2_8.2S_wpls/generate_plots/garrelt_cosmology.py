# Import the integration module
from scipy.integrate import quad
#Import numpy for some mathematical operation
import numpy as np
#Import pylab for plotting
import pylab as pl

# Define some astronomical units
Mpc=3.086e24 # in cm
kms=1e5      # in cm/s
year=3.156e7 # year in seconds

# Define speed of light [km/s]
c=3e5

# Define the integrand
def integrand(x,omm,oml,omr,omk): return 1.0/np.sqrt(omr*((1+x)**4)+omm*((1+x)**3)+oml+omk*((1+x)**2))

# Define the comoving distance as the integral.
# The answer is in units of a Hubble radius (c/H0)
def dc(z,omm,oml,omr,omk): return quad(integrand,0,z,args=(omm,oml,omr,omk))[0]#/(1.0+z)

# Define the redshift array (here from 0 to 20 with steps of 0.1)
z=np.arange(6.0,12.0,0.1)

# Set the Omegas you want to use
omm=0.27
oml=0.73
omr=0.0
omk=1.0-omm-oml-omr

# Set the Hubble constant and calculate the Hubble time
H0=70.0 # in km/sec/Mpc
H0cgs=H0*kms/Mpc # inverse Hubble time in sec

#R0=c/H0/sqrt(-omk)

#-------------------------------------------------------------------
# Part a: making a plot of Da versus z
import matplotlib.pyplot as pl
# Vectorize the age function (so you can give it a vector z as input)
vec_dc = np.vectorize(dc)

# Plot redshift against age
pl.plot(z,vec_dc(z,omm,oml,omr,omk))

# Add labels to the axes
pl.xlabel(r'$z$') # redshift
pl.ylabel(r'$\left(\frac{c}{H_0}\right)d_c$')

#pl.show()
pl.savefig("da_vs_z.png")

#-------------------------------------------------------------------
# Part b: Considering a object of 10 cMpc

# Define comoving size
s_cMpc=10.0

# Calculate the angular size for a range of z [rad]
theta_s=s_cMpc/(c/H0*vec_dc(z,omm,oml,omr,omk))

# Convert radians to arcmin
rad2arcmin=60.*180./np.pi

# Plot angular size [arcsec] as a function of z
pl.plot(z,theta_s*rad2arcmin)
pl.savefig("angularSize_vs_z.png")
#-------------------------------------------------------------------
# Part d: frequency

# 21cm
nu0=1.42e3

# H(z)
Hz=H0*np.sqrt(omr*((1+z)**4)+omm*((1+z)**3)+oml+omk*((1+z)**2))

# Frequency range corresponding to distance s_cMpc
dnu=nu0*Hz*s_cMpc/(c*(1.0+z)**2)

#import matplotlib.pyplot as plt

#plt.plot(z,dnu)
pl.savefig("frequency_vs_z.png")
