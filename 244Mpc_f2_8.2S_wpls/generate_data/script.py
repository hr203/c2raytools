
import map
import mean

#map.map_temp()
#map.map_xfrac('')
#map.map_xfrac('He1')
#map.map_xfrac('He2')
map.map_dbt()

#mean.mean_temp()
#mean.mean_xfrac('')
#mean.mean_xfrac('He1')
#mean.mean_xfrac('He2')
#mean.mean_dbt()

#create powerspectra datafile from dbt
#for i in range(len(redshifts)):
#    data=IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
#    powerspec=ps.power_spectrum_1d(data,500)
#    IO.write2data(powerspec[0],powerspec[1],'data/powerSpectra_'+str('%.3f' % redshifts[i])+'.dat',paths.resutlsdir()+'/powerSpectraFrequencies_dbt_'+str('%.3f' % redshifts[i])+'.dat')
#    #IO.writedata(powerspec[1],'data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')
