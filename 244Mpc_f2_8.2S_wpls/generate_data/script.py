import IO
import map
import mean
import power_spectrum as ps
import sys
sys.path.append('../')
import redshifts as rs

redshifts = rs.read_redshifts("../red_ori2.dat")
mesh=250



#create powerspectra datafile
for i in range(len(redshifts)):
    data=IO.readmap("dbt_"+str('%.3f' % redshifts[i]))
    powerspec=ps.power_spectrum_1d(data)
    IO.writedata(powerspec,'data/powerSpectra_'+str('%.3f' % redshifts[i])+'.dat')

