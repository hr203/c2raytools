'''Module for setting up paths to results directories, data directories, plot directories and redshifts files'''

import numpy


FLAG = 'wpls' #for power laws only
#FLAG = 'wquasars'  #for quasars only
#FLAG = 'wquasars_wpls' #for both

#sim broken
#FLAG = 'wstars' #for neither

basepath = '/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+FLAG+'/' 
results_dir= 'data_'+FLAG+'/'
plots_dir = 'plots_'+FLAG+'/'

def path():
    return basepath

def resultsdir():
#    print results_dir
    return results_dir

def plotsdir():
    return plots_dir

def read_redshifts():

    if FLAG == 'wpls':
        file = "../red_wpls.dat" 
    elif FLAG == 'wquasars':
        file = "../red_wquasars.dat"
    elif FLAG == "wquasars_wpls":
        file = "../red_wquasars_wpls.dat"
    elif FLAG == "wstars":
        file = "../red_wstars.dat"
    else:
        print "Flag incorrect"

    noRedshifts = -1
    with open(file, "r") as f:
        for line in f:
            noRedshifts=noRedshifts+1
    f.close()

    redshifts = numpy.zeros(noRedshifts)

    i=-1
    with open(file,"r") as f2:
        for line in f2:
            if (i!=-1):
                redshifts[i] = line
            i=i+1
    f2.close()
    return redshifts

#def get_dens_redshift(it):
#    if FLAG == 'wpls':
#        it = it 
#    elif FLAG == 'wquasars' or FLAG == "wquasars_wpls" or FLAG == "wstars":
#        if it%2 == 0:
#            it = it
##        else:
#            it = it - 1
#    else:
 #       print "Flag incorrect"
 #   return it



