import corner
import matplotlib.pyplot as plt
import numpy as np

samples = np.load("chain_main.npy")
print samples.shape
fig = corner.corner(samples)
plt.show()
