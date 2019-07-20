import corner
import matplotlib.pyplot as plt
import numpy as np

import sys
pid = int(sys.argv[1])
fname="chain"

data = np.genfromtxt(fname+".dat")
nburn = 1000
nwalker = 100

for i in range(100):
	select = data[:,0]==i
	samples = data[select,2:]

	plt.plot(samples[:,pid])
#plt.axvline(8700, linestyle='--')
plt.show()

