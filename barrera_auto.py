import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numfor as nf
from datetime import datetime

def psi_0(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k*x)

def R():
    return np.sum(prob[:int(n/2)])*dx

def T():
    return np.sum(prob[int(n/2):])*dx

a = 5
n = 1000
k = 20
w = 2

m  = 1
L  = 5
x0 =-L/2
V  = np.zeros(n)
x  = np.linspace(-L, L, n)
dx = (2*L)/(n-1)
dt = 2/(2/(m*dx**2)+np.max(V))

R_arr, T_arr = np.empty(50), np.empty(50)
ample=int(0.5*(n/2)/20)
alcada = np.linspace(0,1.5,50)
for i in range(50):
    psi=psi_0(x)
    V  = np.zeros(n)
    V[int(n/2)-ample:int(n/2)+ample]=alcada[i]*k**2/(2*m)
    nf.euler_rk4(psi=psi, steps=int(0.35/dt), dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()


plt.figure(dpi=300)
plt.plot(alcada,R_arr, 'C0.-', label='R')
plt.plot(alcada,T_arr, 'C1.-', label='T')
plt.xlabel(r'Energia de la barrera partit per $\frac{k^2}{2m}$')
plt.legend()
plt.tight_layout()
plt.savefig('barrera_energia.png')


ample=np.linspace(0,0.5*(n/2)/20,50).astype(int)
for i in range(50):
    psi=psi_0(x)
    V  = np.zeros(n)
    V[int(n/2)-ample[i]:int(n/2)+ample[i]]=0.9*k**2/(2*m)
    nf.euler_rk4(psi=psi, steps=int(0.35/dt), dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()

plt.figure(dpi=300)
plt.plot(2*ample*dx/L,R_arr, 'C0.-', label='R')
plt.plot(2*ample*dx/L,T_arr, 'C1.-', label='T')
plt.xlabel(r'Amplada de la barrera partit per L')
plt.xticks([0, 0.025, 0.05], ['0', 'L/40', 'L/20'])
plt.legend()
plt.tight_layout()
plt.savefig('barrera_amplada.png')
