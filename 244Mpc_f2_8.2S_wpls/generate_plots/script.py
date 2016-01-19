import numpy as np
import plot_1D
import plot_2D
import sys
sys.path.append('../')
import setup_dirs
#sys.path.append('../../src/')
#import c2raytools as c2t

redshifts = setup_dirs.read_redshifts()
#mesh=250

plot_1D.temp_histrogram()
#plot_1D.xfrac_histogram()
#
plot_2D.scatter_contour()
plot_2D.plottemp()
plot_2D.plotxfrac('')
plot_2D.plotxfrac('He1')
plot_2D.plotxfrac('He2')
#plot_2D.equation_of_state()


#these need information from generate_data so leave until last

#plot_1D.dbt_histogram()
#plot_2D.plotdbt()
plot_2D.plotdbt_lightcone()
plot_2D.lightcone()

##plot_1D.plot_mean(plot_1D.mean("temp"),"Tempererature","Temperature (K)")

means = np.zeros(len(redshifts)*3).reshape(len(redshifts),3)
means[:,0] = plot_1D.mean('xfrac')
means[:,1] = plot_1D.mean('xfracHe1')
means[:,2] = plot_1D.mean('xfracHe2')
plot_1D.plot_mean(np.log10(means),"x-frac","Log (Ionised Fraction)")

meantemp=plot_1D.mean('temp')
plot_1D.plot_mean(meantemp,"Temper","Temperature (K)")

#CMBtemp=2.725*(np.ones(len(redshifts))+redshifts)
#meantest=(meantemp-CMBtemp)/meantemp
#plot_1D.plot_mean(meantest,'test','Test')

meandbt = plot_1D.mean('dbt')
plot_1D.plot_mean(meandbt,'dbt','Differential Brightness Temperature')

plot_1D.powerSpec()
plot_1D.allPowerSpec()
plot_1D.rms()
