


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kde

np.random.seed(1977)

# Generate 200 correlated x,y points
data = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], 200)
x, y = data.T

nbins = 20

fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)

axes[0, 0].set_title('Scatterplot')
axes[0, 0].plot(x, y, 'ko')

axes[0, 1].set_title('Hexbin plot')
axes[0, 1].hexbin(x, y, gridsize=nbins)

axes[1, 0].set_title('2D Histogram')
axes[1, 0].hist2d(x, y, bins=nbins)

# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
k = kde.gaussian_kde(data.T)
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

axes[1, 1].set_title('Gaussian KDE')
axes[1, 1].pcolormesh(xi, yi, zi.reshape(xi.shape))

fig.tight_layout()
plt.savefig('fig2.png')











#X,Y,Z = grid(x,y,z)
#plt.contourf(X,Y,Z)
#
#l=100
#x = np.arange(0.0, 22000.0, 22000.0/l)
#y = np.arange(0.0, 1.0, 1.0/l)
#X, Y = np.meshgrid(x, y)
#
#Z = np.zeros(l**2).reshape(l,l)
#for i in range(l):
#    for j in range(l):
#        for t in range(len(temp)-1):
#               if (xfrac[t] > float(i-1)/l and xfrac[t]<float(i)/l):
#                   if(temp[t] > 22000.0*float(j-1)/l and temp[t]< 22000.0*float(j)/l):
#                      Z[i,j] = Z[i,j]+1
#
#
#norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
#cmap = cm.PRGn
#
#levels = np.arange(0.0, l, 0.4)
#
#fig, axes = plt.subplots(1, sharey=True)
##for ax, zord in zip(axes, [1, -1]):
#axes.contourf(X, Y, Z, levels,
#            cmap=cm.get_cmap(cmap, len(levels)-1),
#            norm=norm)
##           zorder=1)#zord)
##axes.scatter([0.0,11000,20000,22000],[0.0,0.4,0.1,0.9])
#print x,y
#axes.scatter(x,y)
#axes.set_title('Scatter with Contour')#.format(zord))
#plt.savefig('fig.png')
