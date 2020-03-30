'''
Programa de càlcul dels coeficients de transmissió i reflexió per diferents alçades
i amples d'una barrera de potencial.
El programa calcula "steps" passos i fa una mesura de R i T per cada tipus de barrera
i per cada un dels tres mètodes de càlcul (Modified Euler, RK4 i Split-Step method)
i ho ploteja.
'''
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

def R():
    return np.sum(prob[:round(n/2)])*dx
def T():
    return np.sum(prob[round(n/2):])*dx

# DEFINIM LES VARIABLES
a  = 5
n  = 1000
k0 = 20
m  = 1
L  = 8
x0 = -L/4

# DEFINIM EL NOSTRE SISTEMS
E = k0**2/(2*m)
dx=(2*L)/(n-1)
x=np.linspace(-L,L,n)
dt=1.e-4
steps=int(0.2/dt) # quantita de passos a calcular
R_arr, T_arr = np.empty(50), np.empty(50)
ample=int(0.5*(n/2)/20)
ample=10
alcada = np.linspace(0,1.5,50)

# INICIALITZEM LA FIGURA I COMENCEM A CALCULAR
plt.figure(dpi=300, figsize=(5,5))
plt.subplot(2,1,1)
plt.yticks(np.linspace(0,1,5))
plt.xticks(np.linspace(0,1.5,7))
plt.grid()

# CÀLCUL DE DIFERENTS ALÇADES AMB MODIFIED EULER
for i in range(50):
    psi=psi_0(x) ; psi_old=psi_0(x)
    V  = np.zeros(n)
    V[int(n/2)-ample:int(n/2)+ample]=alcada[i]*E
    nf.euler_b(psi=psi, psi_old=psi_old, n_steps=steps, dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(alcada,R_arr, 'C1.-', label='Modified Euler', lw=1)
plt.plot(alcada,T_arr, 'C1.-', lw=1)

# CÀLCUL DE DIFERENTS ALÇADES AMB RUNGE-KUTTA 4
for i in range(50):
    psi=psi_0(x)
    V  = np.zeros(n)
    V[int(n/2)-ample:int(n/2)+ample]=alcada[i]*E
    nf.euler_rk4(psi=psi, n_steps=steps, dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(alcada,R_arr, 'C2.-', label='RK-4', lw=1)
plt.plot(alcada,T_arr, 'C2.-', lw=1)

# CÀLCUL DE DIFERENTS ALÇADES AMB SPLIT-STEP
for i in range(50):
    psi=psi_0(x)
    V  = np.zeros(n)
    V[int(n/2)-ample:int(n/2)+ample]=alcada[i]*E
    split_operator(steps)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(alcada,R_arr, 'C3.-', label='Split Operator', lw=1)
plt.plot(alcada,T_arr, 'C3.-', lw=1)

plt.xlabel(r'$E_{wall}/E_{wave}$')
plt.ylim(-0.05,1.05)
plt.xlim(0,1.5)
plt.legend()



R_arr, T_arr = np.empty(20), np.empty(20)
ample=np.array(range(20))
plt.subplot(2,1,2)
plt.ylim(-0.05,1.05)
plt.xlim(-0.01,0.26)

# CÀLCUL DE DIFERENTS AMPLADES AMB MODIFIED EULER
for i in range(20):
    psi=psi_0(x); psi_old=psi_0(x)
    V  = np.zeros(n)
    if ample[i]%2==0:
        V[int(n/2)-ample[i]//2:int(n/2)+ample[i]//2]=0.9*k0**2/(2*m)
    else:
        V[int(n/2)-(ample[i]+1)//2:int(n/2)+(ample[i]-1)//2]=0.9*k0**2/(2*m)
    nf.euler_b(psi=psi, psi_old=psi_old, n_steps=steps, dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(ample*dx,R_arr, 'C1.-', lw=1)
plt.plot(ample*dx,T_arr, 'C1.-', lw=1)

# CÀLCUL DE DIFERENTS AMPLADES AMB RUNGE-KUTTA 4
for i in range(20):
    psi=psi_0(x)
    V  = np.zeros(n)
    if ample[i]%2==0:
        V[int(n/2)-ample[i]//2:int(n/2)+ample[i]//2]=0.9*k0**2/(2*m)
    else:
        V[int(n/2)-(ample[i]+1)//2:int(n/2)+(ample[i]-1)//2]=0.9*k0**2/(2*m)
    nf.euler_rk4(psi=psi, n_steps=steps, dt=dt, m=m, v=V, dx=dx)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(ample*dx,R_arr, 'C2.-', label='RK4', lw=1)
plt.plot(ample*dx,T_arr, 'C2.-', lw=1)

# CÀLCUL DE DIFERENTS AMPLADES AMB SPLIT-STEP OPERATOR
for i in range(20):
    psi=psi_0(x)
    V  = np.zeros(n)
    if ample[i]%2==0:
        V[int(n/2)-ample[i]//2:int(n/2)+ample[i]//2]=0.9*k0**2/(2*m)
    else:
        V[int(n/2)-(ample[i]+1)//2:int(n/2)+(ample[i]-1)//2]=0.9*k0**2/(2*m)
    split_operator(steps)
    prob = np.abs(psi)**2
    R_arr[i], T_arr[i] = R(), T()
plt.plot(ample*dx,R_arr, 'C3.-', label='RK4', lw=1)
plt.plot(ample*dx,T_arr, 'C3.-', lw=1)


plt.xlabel(r'Wall width (r.u.)')
plt.xticks(np.linspace(0,0.25,6))
plt.yticks(np.linspace(0,1,5))
plt.grid()
plt.tight_layout()
plt.savefig('barrera.png')
