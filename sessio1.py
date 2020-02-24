import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import numfor as nf
from datetime import datetime

def complex_fun(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k0*x)

high = float(input('Alçada de la barrera de potencial ( de l.ordre de 15.0)? '))

a  = 3
n  = 1000
k0 = 30

x0 =-3.0
m  = 50
L=5
V = np.zeros(n)
V[int(n/2):]=high
x=np.linspace(-L,L,n)
dx=(2*L)/n
dt=0.5/(2/(m*dx**2)+np.max(V))
t = 0.0
psi=complex_fun(x) ; psi_old = complex_fun(x)
freq = np.fft.fftshift(np.fft.fftfreq(n)*(2*np.pi/dx))

fig = plt.figure(dpi=200)
ax = fig.gca(projection='3d')
ax.set_xlim(-4.5,-1.5)
ax.plot(x, np.real(psi), np.imag(psi), label='$\Psi$')
plt.tight_layout()
plt.savefig('plot3d.png')
plt.close()


fig, axs = plt.subplots(2, figsize=(5,5))
fig.set(dpi=170)
axs[0].set(xlim=(-L, L), ylim=(-1.2*(a**2/(2*np.pi))**(1/2), 1.2*(a**2/(2*np.pi))**(1/2)), xlabel='x', ylabel='$\Psi$', title=('Barrera potencial de %.1f' % high))
axs[1].set(xlim=(-40,40),  ylim=(0,15), xlabel='k', ylabel='FFT')

axs[0].plot(x, V/max(V), 'r-')
abs_line,  = axs[0].plot(x, x, 'k-', label='Probability')
real_line, = axs[0].plot(x, x, lw=1, label='Real')
imag_line, = axs[0].plot(x, x, lw=1, label='Imaginary')
fft_line,  = axs[1].plot(freq, freq)
title_text = axs[0].text(0.02, 0.92, '', transform=axs[0].transAxes)
norm_text = axs[0].text(0.02, 0.84, '', transform=axs[0].transAxes)
axs[0].legend(ncol=3, loc=8)
plt.tight_layout()

def animate(i):
    global t
    nf.euler_b(psi=psi, psi_old=psi_old, steps=30, dt=dt, m=m, v=V, dx=dx)
    # nf.euler(psi=psi, steps=20, dt=dt, m=m, v=V, dx=dx)
    t += 30*dt
    fft = np.abs(np.fft.fftshift(np.fft.fft(psi, norm="ortho")))**2
    prob = np.abs(psi)**2
    norm = np.sum(prob)*dx
    real_line.set_ydata(np.real(psi))
    imag_line.set_ydata(np.imag(psi))
    abs_line.set_ydata(prob)
    fft_line.set_ydata(fft)
    title_text.set_text(('Time: %.3f' % t))
    norm_text.set_text(('Norma: %.3f' % norm))
    return real_line, imag_line, abs_line, fft_line, title_text, norm_text
for i in range(500): animate(0)
plt.savefig('screenshot.png')

ani = animation.FuncAnimation(fig, animate, frames=100, interval=0, blit=True, repeat=True)
plt.show()

# print(datetime.now())
# plt.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'
# ani = animation.FuncAnimation(fig, animate, frames=1000, blit=True, repeat=False)
# ani.save('entrega.mp4', writer='ffmpeg', fps=30)
# print(datetime.now())
