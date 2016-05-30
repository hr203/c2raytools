import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../')
import setup_dirs

redshifts = setup_dirs.read_redshifts()

two_wstars = np.loadtxt('data_wstars/2e4/mean_dbt.dat')
onehalf_wstars = np.loadtxt('data_wstars/1.5e4/mean_dbt.dat')
one_wstars =np.loadtxt('data_wstars/1e4/mean_dbt.dat')
plt.figure()
plt.plot(redshifts,two_wstars,label='THII = 2e4')
plt.plot(redshifts,one_wstars,label='THII = 1e4')
plt.plot(redshifts,onehalf_wstars,label='THII = 1e4')
plt.legend()
plt.ylabel("Mean")
plt.xlabel("z")
plt.savefig('compare_mean.png')

two_wstars = np.loadtxt('data_wstars/2e4/_rms_.dat')
one_wstars =np.loadtxt('data_wstars/1e4/_rms_.dat')
onehalf_wstars =np.loadtxt('data_wstars/1.5e4/_rms_.dat')
plt.plot(redshifts,one_wstars,label='THII = 1e4')
plt.figure()
plt.plot(redshifts,two_wstars,label='THII = 2e4')
plt.plot(redshifts,onehalf_wstars,label='THII = 1.5e4')
plt.plot(redshifts,one_wstars,label='THII = 1e4')
plt.legend()
plt.ylabel("rms")
plt.xlabel("z")
plt.savefig('compare_rms.png')

two_wstars = np.loadtxt('data_wstars/2e4/skewness.dat')
one_wstars =np.loadtxt('data_wstars/1e4/skewness.dat')
onehalf_wstars =np.loadtxt('data_wstars/1.5e4/skewness.dat')
plt.figure()
plt.plot(redshifts,two_wstars,label='THII = 2e4')
plt.plot(redshifts,onehalf_wstars,label='THII = 1.5e4')
plt.plot(redshifts,one_wstars,label='THII = 1e4')
plt.legend()
plt.ylabel("Skewness")
plt.xlabel("z")
plt.savefig('compare_skewness.png')


two_wstars = np.loadtxt('data_wstars/2e4/kurtosis.dat')
one_wstars =np.loadtxt('data_wstars/1e4/kurtosis.dat')
onehalf_wstars =np.loadtxt('data_wstars/1.5e4/kurtosis.dat')
plt.figure()
plt.plot(redshifts,two_wstars,label='THII = 2e4')
plt.plot(redshifts,onehalf_wstars,label='THII = 1.5e4')
plt.plot(redshifts,one_wstars,label='THII = 1e4')
plt.legend()
plt.ylabel("Kurtosis")
plt.xlabel("z")
plt.savefig('compare_kurtosis.png')


