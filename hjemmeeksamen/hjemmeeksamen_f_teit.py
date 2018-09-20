from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")
# Tools for sparse matrices
import scipy.sparse as sparse
import scipy.sparse.linalg

"""Physical constants"""
E0 = 4*938.27       # alpha-particle rest energy [MeV]
hbarc = 197.32      # [MeV fm]
c = 2.998e23      # Speed of light [fm / s]

deltax = 17 #[fm]
V0 = 34 #[MeV]

x_tot = 50   #[fm]
dx = 0.1         #??????????
nx = int(x_tot/dx)
x = np.linspace(0, x_tot, nx)

T = 2e-20 #[s]
dt = 1e-23 #[s]
t = 0
num_time_steps = int(T/dt)


def psi0(x):
    #Initial state for a travelling niggerfaggot

    x0 = 5    #[fm]
    a = 1     #[fm]

    A = (1/(2 * np.pi * a**2))**0.25
    K1 = np.exp(-(x-x0)**2 / (4*a**2))

    return A*K1

def potential(x):
    pot_height = V0
    pot_width = deltax

    potential = np.zeros(len(x))
    potential[int(7.3/dx):int((7.3 + deltax)/dx)] = pot_height

    return potential

k1 = -(1j * hbarc * c)/(2*E0)
k2 = (1j * c)/hbarc

#init psi
psi = psi0(x)

#krerifisering av matrisen som konteiner sentrale differenser
data = np.ones((3, nx))
data[1] = -2 * data[1]
diags = [-1, 0, 1]
D2 = k1/dx**2 * sparse.spdiags(data, diags, nx, nx)

#Identitetsneger
I = sparse.identity(nx)

#Krer diagonaaaal matrifixen som konteiner potensialet
V_data = potential(x)
V_diags = [0]
V = k2*sparse.spdiags(V_data, V_diags, nx, nx)

fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot( x, abs(psi)**2, label='$|\Psi(x,t)|^2$' )

ax.plot(x, V_data, label="$V(x)$")
#ax.set_ylim(0, 0.01)
ax.set_xlabel("$x$ [fm]")
ax.set_ylabel("$|\Psi(x, t)|^2$ [1/fm] (og $V(x)$ [MeV])")
plt.draw()
while t<T:

    A = (I - dt/2*(D2 + V))
    b = (I + dt/2 * (D2 + V)) * psi

    psi = sparse.linalg.spsolve(A, b)

    t += dt

    line.set_ydata(abs(psi)**2)
    plt.draw()
    plt.pause(0.001)

plt.show()