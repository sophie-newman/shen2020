import numpy as np
import os, errno
import __main__

#COSMOLOGY
hubble      = 0.6774
Omega0      = 0.3089
OmegaBaryon = 0.0486
OmegaLambda = 0.6911

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
