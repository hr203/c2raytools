import numpy as np
import plot_1D
import plot_2D
import sys
sys.path.append('../')
import setup_dirs

redshifts_wstars = setup_dirs.read_redshifts('wstars')
redshifts_wpls = setup_dirs.read_redshifts('wpls')

def compare_mean_temp():
    wstars_temp = plot_1D.mean('temp','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('temp','data_wpls/',redshifts_wpls)

    length = min(len(wstars_temp),len(wpls_temp))

    means = np.zeros(2*length).reshape(length,2)
    for i in range(length):
        means[i,0] = wstars_temp[i]
    for i in range(length):
        means[i,1] = wpls_temp[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"Temperature","Temperature (K)",'compare_wstars_wpls/',redshifts_wpls,"Stellar only","Stellar and X-ray")

def compare_mean_dbt():
    wstars_temp = plot_1D.mean('dbt','data_wstars/',redshifts_wstars)
    wpls_temp = plot_1D.mean('dbt','data_wpls/',redshifts_wpls)

    length = min(len(wstars_temp),len(wpls_temp))/2
    means = np.zeros(2*length).reshape(length,2)
    redshifts2 = np.zeros(length)
    for i in range(length*2):
        if i%2==0:
            means[i/2,0] = wstars_temp[i]
            means[i/2,1] = wpls_temp[i]
            redshifts2[i/2]=redshifts_wpls[i]
        #print means[i,1]
    #print redshifts_wstars[i], redshifts_wpls[i]
    plot_1D.plot_mean(means,"dbt","$\delta T_b$",'compare_wstars_wpls/',redshifts2,"Stellar only","Stellar and X-ray")



def compare_mean_xfrac():
    wstars_xfrac = plot_1D.mean('xfrac','data_wstars/',redshifts_wstars)
    wpls_xfrac = plot_1D.mean('xfrac','data_wpls/',redshifts_wpls)
    wstars_xfracHe1 = plot_1D.mean('xfracHe1','data_wstars/',redshifts_wstars)
    wpls_xfracHe1 = plot_1D.mean('xfracHe1','data_wpls/',redshifts_wpls)
    wstars_xfracHe2 = plot_1D.mean('xfracHe2','data_wstars/',redshifts_wstars)
    wpls_xfracHe2 = plot_1D.mean('xfracHe2','data_wpls/',redshifts_wpls)

    length = min(len(wstars_xfrac),len(wpls_xfrac))

    means = np.zeros(6*length).reshape(length,6)
    for i in range(length):
        means[i,0] = wpls_xfrac[i]
        means[i,1] = wpls_xfracHe1[i]
        means[i,2] = wpls_xfracHe2[i]
        means[i,3] = wstars_xfrac[i]
        means[i,4] = wstars_xfracHe1[i]
        means[i,5] = wstars_xfracHe2[i]

    means[0,2] = means[1,2]
    means[0,5] = means[1,5]  
   #print means[:,2], means[:,5]
    #print redshifts_wstars[i], redshifts_wpls[i]

    plot_1D.plot_mean(means,"xfrac","Ionised Fraction (%)",'compare_wstars_wpls/',redshifts_wpls,"Stellar and X-ray binaries","Stellar only")

def compare_hisograms_temperature():
    for i in range(12,len(redshifts)):
        temp_filename = setup_dirs.path() + 'Temper3D_'+str('%.3f' % redshifts[i]) + '.bin'
        data = c2t.TemperFile(temp_filename).temper
        data = data.flatten() #need to import other dataset and put in right format
#        plot.1D.plot_log_histogram(data,"Temperature, Redshift:" +str('%.3f' % redshifts[i]),"loghist_temper_shortrange_"+str(i+10)+'_'+str('%.3f' % redshifts[i]),"Temperature(K)")


compare_mean_xfrac()
compare_mean_dbt()
compare_mean_temp()
