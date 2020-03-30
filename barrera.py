import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy.fft import fft as FFT, ifft as IFFT
import numfor as nf

def psi_0(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k0*x)

def split_operator(steps):
    global psi
    for i in range(steps):
        psi =IFFT(FFT(psi*x_operator)*k_operator)*x_operator

def R():
    return np.sum(prob[:int(n/2)])*dx

def T():
    return np.sum(prob[int(n/2)+ample:])*dx


a  = 5
n  = 1000
k0 = 20
m  = 1
L  = 8
x0 = -L/4

E = k0**2/(2*m)

V = np.zeros(n)
dx=(2*L)/(n-1)
ample=int(0.5*(n/2)/40)
V[int(n/2)-ample:int(n/2)+ample]=E*0.9
x=np.linspace(-L,L,n)
# dt=0.5/(2/(m*dx**2)+np.max(V))
dt=1.e-4
t = 0.0
psi=psi_0(x) ; psi_old = psi_0(x)
freq = np.fft.fftshift(np.fft.fftfreq(n)*(2*np.pi/dx))

x_operator = np.exp(-0.5j*dt*V)
k_operator = np.exp(-0.5j*dt*(np.fft.fftfreq(n)*(2*np.pi/dx))**2/m)

fig, axs = plt.subplots(2, figsize=(6,6), dpi=200)
axs[0].set(xlim=(-L, L), ylim=(-1.2*(a**2/(2*np.pi))**(1/4), 1.2*(a**2/(2*np.pi))**(1/4)), xlabel='x', ylabel='$\Psi$')
axs[1].set(xlim=(freq[400],freq[-400]),  ylim=(0,5), xlabel='k', ylabel='FFT')
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
    # nf.euler_b(psi=psi, psi_old=psi_old, n_steps=30, dt=dt, m=m, v=V, dx=dx)
    # nf.euler(psi=psi, n_steps=30, v=V, dt=0.05*dt, m=m, dx=dx)
    # nf.euler_rk4(psi=psi, n_steps=20, dt=dt, m=m, v=V, dx=dx)
    split_operator(100)
    t += 20*dt
    fft  = np.abs(np.fft.fftshift(np.fft.fft(psi, norm="ortho")))**2
    prob = np.abs(psi)**2
    norm = np.sum(prob)*dx
    norm = np.sum(fft)*2*np.pi/dx
    real_line.set_ydata(np.real(psi))
    imag_line.set_ydata(np.imag(psi))
    abs_line.set_ydata(prob)
    fft_line.set_ydata(fft)
    title_text.set_text(('Time: %.3f' % t))
    norm_text.set_text(("R={0:.0%}".format(R()))+
                      (" T={0:.0%}".format(T()))+
                      (" Norma=%.8f" % norm))
    return real_line, imag_line, abs_line, fft_line, title_text, norm_text,

ani = animation.FuncAnimation(fig, animate, frames=150, interval=50, blit=True, repeat=False)

plt.show()
# ani.save('animation1.mp4')

# print(datetime.now())
# plt.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'
# ani = animation.FuncAnimation(fig, animate, frames=1000, blit=True, repeat=False)
# ani.save('entrega.mp4', writer='ffmpeg', fps=30)
# print(datetime.now())
