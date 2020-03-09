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
    return np.sum(prob[int(n/2)+ample:])*dx

# high = float(input('Alçada de la barrera de potencial ( de l.ordre de 15.0)? '))
# BC = int(input('Condicions de contorn periòdiques (0), potencial infinit (1) o obertes (2) ? '))
BC=1
a  = 10
n  = 1000
k = 50

m  = 1
L=10
x0 =-L/2
E = k**2/(2*m)
V = np.zeros(n)
dx=(2*L)/(n-1)
ample=int(0.5*(n/2)/20)
V[int(n/2)-ample:int(n/2)+ample]=E*1.5
x=np.linspace(-L,L,n)
dt=0.5/(2/(m*dx**2)+np.max(V))
t = 0.0
sig = np.sqrt(2)/a
psi=psi_0(x) ; psi_old = psi_0(x)
freq = np.fft.fftshift(np.fft.fftfreq(n)*(2*np.pi/dx))



fig, axs = plt.subplots(2, figsize=(6,6))
# fig.set(dpi=150)
axs[0].set(xlim=(-L, L), ylim=(-1.2*(a**2/(2*np.pi))**(1/4), 1.2*(a**2/(2*np.pi))**(1/4)), xlabel='x', ylabel='$\Psi$')
axs[1].set(xlim=(-40,40),  ylim=(0,15), xlabel='k', ylabel='FFT')

axs[0].plot(x, V/max(V), 'r-')
real_line, = axs[0].plot(x, x, lw=1, label='Real')
imag_line, = axs[0].plot(x, x, lw=1, label='Imaginary')
abs_line,  = axs[0].plot(x, x, 'k-', label='Probability')
fft_line,  = axs[1].plot(freq, freq)
title_text = axs[0].text(0.02, 0.92, '', transform=axs[0].transAxes)
norm_text  = axs[0].text(0.02, 0.84, '', transform=axs[0].transAxes)
axs[0].legend(ncol=3, loc=8)
plt.tight_layout()

def animate(i):
    global t, prob
    # nf.euler_b(psi=psi, psi_old=psi_old, steps=30, dt=dt, m=m, v=V, dx=dx, bc=BC)
    nf.euler_rk4(psi=psi, n=n, steps=30, dt=dt, m=m, v=V, dx=dx)
    t += 30*dt
    fft = np.abs(np.fft.fftshift(np.fft.fft(psi, norm="ortho")))**2
    prob = np.abs(psi)**2
    norm = np.sum(prob)*dx
    real_line.set_ydata(np.real(psi))
    imag_line.set_ydata(np.imag(psi))
    abs_line.set_ydata(prob)
    fft_line.set_ydata(fft)
    title_text.set_text(('Time: %.3f' % t))
    norm_text.set_text((" R={0:.0%}".format(R()))+
                       (" T={0:.0%}".format(T()))+
                       (" Norma=%.3f" % norm))
    return real_line, imag_line, abs_line, fft_line, title_text, norm_text,

ani = animation.FuncAnimation(fig, animate, frames=350, interval=0, blit=True, repeat=True)
plt.show()

# print(datetime.now())
# plt.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'
# ani = animation.FuncAnimation(fig, animate, frames=1000, blit=True, repeat=False)
# ani.save('entrega.mp4', writer='ffmpeg', fps=30)
# print(datetime.now())
