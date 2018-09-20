from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

Za =
delta_x = 17 #[fm]
m = 3728.401 #[MeV/c^2]
c = 2.998e23 #[fm/s]
hbar = 6.582e-22 #[MeVs]

def T(E):
    return 1/(1 + V0**2 / (4*E*(V0 - E))* np.sinh(delta_x * np.sqrt(2*m/(c**2) * (V0 - E))/hbar)**2)

E = np.linspace(0.1, 10, 100000)
plt.plot(E, T(E))
plt.xlabel("E [MeV]")
plt.ylabel("T(E)")
plt.savefig("oppg_d.png")
plt.show()