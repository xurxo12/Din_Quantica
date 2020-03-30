import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft as FFT, ifft as IFFT
import numfor as nf
from datetime import datetime

def psi_0(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k0*x)

def split_operator(steps):
    global psi
    x_operator = np.exp(-0.5j*dt*V)
    k_operator = np.exp(-0.5j*dt*(np.fft.fftfreq(n)*(2*np.pi/dx))**2/m)
    for i in range(steps):
        psi =IFFT(FFT(psi*x_operator)*k_operator)*x_operator



a  = 5
n  = 1000
k0 = 20
m  = 1
L  = 8
x0 = -L/4

E = k0**2/(2*m)


dx=(2*L)/(n-1)

x=np.linspace(-L,L,n)

dt=1.e-6
N=1000
steps=100

SE = np.empty(N)
ME = np.empty(N)
RK = np.empty(N)
SS = np.empty(N)
time = np.empty(N)

V  = np.zeros(n)
V[int(n/2)-10:int(n/2)+10]=0.9*E

plt.figure(dpi=300, figsize=(6,3))
plt.grid()
psi1=psi_0(x)
psi2=psi_0(x)
psi3=psi_0(x)
psi4=psi_0(x)
psi=psi_0(x)
for i in range(N):
    nf.euler(psi=psi1, n_steps=steps, v=V, dt=1.e-6, m=m, dx=dx)
    nf.euler_b(psi=psi2, psi_old=psi3, n_steps=steps, dt=1.e-4, m=m, v=V, dx=dx)
    nf.euler_rk4(psi=psi4, n_steps=steps, dt=3.55e-4, m=m, v=V, dx=dx)
    split_operator(steps)

    SE[i]=np.sum(np.abs(psi1)**2)*dx
    ME[i]=np.sum(np.abs(psi2)**2)*dx
    RK[i]=np.sum(np.abs(psi4)**2)*dx
    SS[i]=np.sum(np.abs(psi)**2)*dx
    time[i]=(i+1)*steps

plt.plot(time, ME, 'C1', label='Modified Euler', lw=1)
plt.plot(time, SE, 'C0', label='Simple Euler')
plt.plot(time, RK, 'C2', label='RK')
plt.plot(time, SS, 'C3', label='Split-Step')
plt.ylim(0.993,1.007)
plt.legend(ncol=2)
plt.xlabel('Steps')
plt.xlim(0,100000)
plt.ylabel('Norm')
plt.tight_layout()
plt.savefig('stability.png')
