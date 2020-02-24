import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numfor as nf
from copy import copy

def complex_fun(x,t):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-a**2*(x-x0)**2/4)*np.exp(1j*(k0*x-w*t))

a  = 5
n  = 500
k0 = 20

w  = 1.0
x0 = -3.0
m=10
V = np.zeros(n)
V[int(1.5*n/2):]=20.0
L=5
x=np.linspace(-L,L,n)
dx=(2*L)/n
dt=1/(2/(m*dx**2)+np.max(V))
dt=2.e-4
print(dt)
t=0.0
psi=np.array(complex_fun(x, 0), order='F')
psi[0] =0
psi[-1]=0
psi_old = copy(psi)
freq = np.fft.fftfreq(n)*(2*np.pi/dx)

# fig = plt.figure(figsize=(6,4), dpi=120)
fig, axs = plt.subplots(2)
fig.set(dpi=130)
axs[0].set(xlim=(-L, L), ylim=(-1.5*(a**2/(2*np.pi))**(1/2), 1.5*(a**2/(2*np.pi))**(1/2)))
axs[1].set(xlim=(-30,30),  ylim=(0,7))
plt.tight_layout()
axs[0].plot(x, 0.4*V, 'r-')
# title_text = axs[0].text(0.02, 0.95, '', transform=axs[0].transAxes)
abs_line,  = axs[0].plot(x, x, 'k-')
real_line, = axs[0].plot(x, x, lw=1)
imag_line, = axs[0].plot(x, x, lw=1)
fft_line, = axs[1].plot(freq, freq)

def animate(i):
    global t
    nf.euler_b(psi=psi, psi_old=psi_old, steps=20, dt=dt, m=m, v=V, dx=dx)
    # prop = np.abs(psi**2)
    # if i%500 == 0:
    #     print(np.sum(prop)*dx)
    fft = np.abs(np.fft.fft(psi, norm="ortho"))**2
    real_line.set_ydata(np.real(psi))
    imag_line.set_ydata(np.imag(psi))
    abs_line.set_ydata(np.abs(psi)**2)
    fft_line.set_ydata(fft)
    # title_text.set_text(('Polaritzation: %.4f' % t))
    return real_line, imag_line, abs_line, fft_line,



ani = animation.FuncAnimation(fig, animate, frames=1000, interval=0, blit=True, repeat=True)
plt.show()
