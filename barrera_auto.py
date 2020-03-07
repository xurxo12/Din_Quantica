import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import numfor as nf
from datetime import datetime

def psi_0(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k*x)

def R():
    return np.sum(prob[:int(n/2)])*dx

def T():
    return np.sum(prob[int(n/2):])*dx

# BC = int(input('Condicions de contorn periòdiques (0), potencial infinit (1) o obertes (2) ? '))
BC=1
a = 5
n = 1000
k = 20
w = 2

m  = 1
L=5
x0 =-L/2
V = np.zeros(n)
ample=int(0.5*(n/2)/20)
x=np.linspace(-L,L,n)
dx=(2*L)/(n-1)
dt=0.9/(2/(m*dx**2)+np.max(V))

R_arr=np.empty(50)
T_arr=np.empty(50)
alcada = np.linspace(0,1.5,50)
for i in range(50):
    psi=psi_0(x) ; psi_old = psi_0(x)
    V[int(n/2)-ample:int(n/2)+ample]=alcada[i]*k**2/(2*m)
    nf.euler_b(psi=psi, psi_old=psi_old, steps=int(0.35/dt), dt=dt, m=m, v=V, dx=dx, bc=BC)
    prob = np.abs(psi)**2
    R_arr[i]=R()
    T_arr[i]=T()

plt.figure()
plt.plot(alcada,R_arr, 'C0.-', label='R')
plt.plot(alcada,T_arr, 'C1.-', label='T')
plt.legend()
plt.show()


# fig, axs = plt.subplots(1, figsize=(6,6))
# axs.plot(x, V/max(V), 'r-')
# axs.set(ylim=(-2,2))
# real_line, = axs.plot(x, np.real(psi), lw=1, label='Real')
# imag_line, = axs.plot(x, np.imag(psi), lw=1, label='Imaginary')
# abs_line,  = axs.plot(x, prob, 'k-', label='Probability')
# plt.tight_layout()
# plt.show()
