import IO
import map
import mean
import power_spectrum as ps
import sys
sys.path.append('../')
import redshifts as rs

redshifts = rs.read_redshifts("../red_ori2.dat")
mesh=250



#create powerspectra datafile from dbt
for i in range(len(redshifts)):
    data=IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
    powerspec=ps.power_spectrum_1d(data,500)
    IO.write2data(powerspec[0],powerspec[1],'data/powerSpectra_'+str('%.3f' % redshifts[i])+'.dat','data/powerSpectraFrequencies_dbt_'+str('%.3f' % redshifts[i])+'.dat')
    #IO.writedata(powerspec[1],'data/powerSpectraFrequencies_'+str('%.3f' % redshifts[i])+'.dat')
