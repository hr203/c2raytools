import numpy
import paths



def read_redshifts(file):

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
