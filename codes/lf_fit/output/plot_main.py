import corner
import matplotlib.pyplot as plt
import numpy as np

import sys
fname=sys.argv[1]

#samples = np.load("chain_main.npy")
data = np.genfromtxt(fname+".dat")
id = data[:,1]>1000
samples = data[id,2:]
print samples.shape

best_fit = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
          zip(*np.percentile(samples, [16, 50, 84],axis=0)))
best_fit = np.array(best_fit)
print best_fit.shape
for i in range(best_fit.shape[0]):
	print "p{0:d}:".format(i), best_fit[i,0], "+{0:f}".format(best_fit[i,1]),"-{0:f}".format(best_fit[i,2])
fig = corner.corner(samples,top_ticks=True)#range=np.ones(samples.shape[1])*0.96, top_ticks=True)
fig.savefig("triangle.png",format='png')
