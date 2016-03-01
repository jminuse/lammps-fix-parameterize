from merlin import *

import matplotlib.pyplot as plt
import numpy as np

eps = 1.0
sig = 3.0
lj = lambda r: 4*eps*( (sig/r)**12 - (sig/r)**6 )
r_cutoff = 4.
cut = lambda r: (r**6 + r_cutoff**6)**(1./6)

lj2 = lambda r: 4*eps*( (sig/cut(r))**12 - (sig/cut(r))**6 )

xx = np.arange(2.,10.,0.1)
Es = [lj(r) for r in xx]
E2s = [lj2(r) for r in xx]

plt.plot(xx, Es)
plt.plot(xx, E2s)
plt.ylabel('E')
plt.ylim(-10.0, 10.0)
plt.show()


