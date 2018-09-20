from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

x0 = 5 #[fm]
a = 1 #[fm]
k = 1.38 #[1/fm]

def psi(x):
    return np.sqrt(1/(2*np.pi*a**2))*np.exp(-2*(x-x0)**2/(4*a**2))

x = np.linspace(0, 10, 1000)
plt.plot(x, psi(x))
plt.xlabel("x [fm]")
plt.ylabel("$|\Psi(x, 0)|^2$ [1/fm]")
plt.savefig("oppg_a.png")
plt.show()