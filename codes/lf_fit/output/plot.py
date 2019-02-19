import corner
import matplotlib.pyplot as plt
import numpy as np

samples = np.load("chain.npy")
print samples.shape
fig = corner.corner(samples, labels=["g1", "g2", "Phi" , "L"])
plt.show()
