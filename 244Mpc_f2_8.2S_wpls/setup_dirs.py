'''Module for setting up paths to results directories, data directories, plot directories and redshifts files'''

import numpy


#FLAG = 'wpls' #for power laws only
#FLAG = 'wquasars_wpls' #for both
#FLAG = 'wquasars_0.8'
FLAG = 'wquasars_1.6'

#sims redone on apollo#FLAG = 'wstars' #for neither
#FLAG = 'wstars'
#FLAG = 'wstars_equalphotons' #for neither

#tests
#FLAG = 'kylsim'
#FLAG = 'test'
#res=30
#steps = 2

if FLAG!='test':
    results_dir= 'data_'+FLAG+'/'
    plots_dir = 'plots_'+FLAG+'/'
    basepath ='/lustre/scratch/astro/hr203/RESULTS/244Mpc_f2_8.2S_H250_'+FLAG+'_2/' 
else:
    results_dir='tests/data_'+str(res)+'_2source/'
    plots_dir='tests/data_'+str(res)+'_2source/'
    basepath = '/lustre/scratch/astro/hr203/RESULTS/tests/'+str(res)+'_xbox_1.6/'


def get_flag():
    return FLAG


def get_res():
    if FLAG=='test':
        return res
    else:
        return ' ' 

def path():
    return basepath

def resultsdir():
#    print results_dir
    return results_dir

def plotsdir():
    return plots_dir

def read_redshifts(flag = FLAG):

    file = '../red_' + flag +'.dat'

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


def nosteps():
    steps =2
    if FLAG == 'wpls_higheff':
        steps=40
    return steps

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



