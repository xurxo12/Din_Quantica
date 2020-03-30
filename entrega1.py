import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numfor as nf

def psi_0(x):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*(a*(x-x0))**2)*np.exp(1j*k0*x)

# BLOC DE DEFINICIÓ DE VARIABLES
a  = 5
n  = 1000
k0 = 20
m  = 1
L  = 8
x0 = -L/2
E = k0**2/(2*m)

# BLOC DE CONSTRUCCIO DEL SISTEMA
V = np.zeros(n)
V[int(n/2):]=E
x=np.linspace(-L,L,n)
dx=(2*L)/(n-1)
dt=1/(2/(m*dx**2)+np.max(V))
t = 0.0
psi=psi_0(x) ; psi_old = psi_0(x)
freq = np.fft.fftshift(np.fft.fftfreq(n))*2*np.pi/dx
some_real_array=x

# BLOC DE INICIALITZACIO DE L'ANIMACIÓ
fig, axs = plt.subplots(2, figsize=(6,6))
axs[0].set(xlim=(-L, L), ylim=(-1.3*np.max(np.abs(psi)), 1.3*np.max(np.abs(psi))), xlabel='x', ylabel='$\psi$', title='Barrera potencial de E')
axs[1].set(xlim=(-40,40),  ylim=(0,5), xlabel='k', ylabel='FFT')
axs[0].plot(x, V/max(V), 'r-')
real_line, = axs[0].plot(x, some_real_array, lw=1, label='Real')
imag_line, = axs[0].plot(x, some_real_array, lw=1, label='Imaginary')
abs__line, = axs[0].plot(x, some_real_array, 'k-', label='Probability')
fft__line, = axs[1].plot(freq, some_real_array)
titl_text  = axs[0].text(0.02, 0.92, '', transform=axs[0].transAxes)
norm_text  = axs[0].text(0.02, 0.84, '', transform=axs[0].transAxes)
axs[0].legend(ncol=3, loc=8)
plt.tight_layout()

def animate(i): # FUNCIÓ D'ANIMACIÓ
    global t
    nf.euler_b(psi=psi, psi_old=psi_old, n_steps=30, dt=dt, m=m, v=V, dx=dx)
    t += 30*dt
    fft  = np.abs(np.fft.fftshift(np.fft.fft(psi, norm="ortho")))**2
    prob = np.abs(psi)**2
    norm = np.sum(prob)*dx
    real_line.set_ydata(np.real(psi))
    imag_line.set_ydata(np.imag(psi))
    abs__line.set_ydata(prob)
    fft__line.set_ydata(fft)
    titl_text.set_text(('Time: %.3f' % t))
    norm_text.set_text(('Norma: %.5f' % norm))
    return real_line, imag_line, abs__line, fft__line, titl_text, norm_text,

ani = animation.FuncAnimation(fig, animate, frames=100, interval=10, blit=True, repeat=True)
plt.show()
# ani.save('video_entrega1.mp4')
