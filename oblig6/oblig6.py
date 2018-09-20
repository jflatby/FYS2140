from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")


h = 6.626e-34
hbar = h/(2*np.pi)
k = 1
m = (hbar)**2/k
omega = np.sqrt(k/m)
gamma = m*omega/hbar

x = np.linspace(-5, 5, 10000)

def psi_0(x):
    return (gamma/np.pi)**(1/4)*np.exp(-(gamma/2)*x**2)

def psi_1(x):
    return (gamma/np.pi)**(1/4) * np.sqrt(2*gamma)* x *np.exp(-(gamma/2)*x**2)

def psi_2(x):
    return 1/np.sqrt(2) * (gamma/np.pi)**(1/4)*(2*gamma*x**2-1)*np.exp(-(gamma/2)*x**2)

plt.plot(x, psi_0(x), x, psi_1(x), x, psi_2(x))
plt.xlabel("x [nm]")
plt.ylabel("Psi(x) [nm^(-1/2)]")
plt.legend(["Psi_0", "Psi_1", "Psi_2"])
plt.savefig("oblig6.png")
plt.show()