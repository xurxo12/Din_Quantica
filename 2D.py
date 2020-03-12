import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numfor as nf

def psi_0(x,y):
    return (a**2/(2*np.pi))**(1/4)*np.exp(-0.25*a**2*((x-x0)**2+(y-y0)**2))*np.exp(1j*(kx*x+ky*y))

BC = 1
a  = 1
n  = 250
kx = 0
ky = 40
steps=20

m  = 1
L  = 10
x0 = 0.0
y0 = 0.0
E = (kx**2+ky**2)/(2*m)
V = np.zeros((n,n))
dx=(2*L)/(n-1)
ample=int(0.5*(n/2)/20)
x=np.linspace(-L,L,n)
y=np.linspace(-L,L,n)
dt=1.e-3
t = 0.0
X = np.linspace(-L/2, L/2, n)
Y = np.linspace(-L/2, L/2, n)
X, Y = np.meshgrid(X, Y)
psi=np.array(psi_0(X, Y), order='F')
psi /= np.sqrt(np.sum(np.abs(psi)**2)*dx**2)


fig, ax = plt.subplots(1)
surf = ax.imshow( np.abs(psi)**2, cmap='inferno')
title_text = ax.text(0.02, 0.92, '', transform=ax.transAxes, color='w')
norm_text  = ax.text(0.02, 0.85, '', transform=ax.transAxes, color='w')
ax.set_xticks([]) ; ax.set_yticks([]) ; plt.tight_layout()

def animate(i):
    global t, prob
    nf.euler_rk4_3d(psi=psi, steps=steps, dt=dt, m=m, dx=dx)
    t += steps*dt
    prob = np.abs(psi)**2
    surf.set_data(prob)
    title_text.set_text(('Time: %.3f' % t))
    norm_text.set_text((" Norma=%.5f" % (np.sum(prob)*dx**2)))
    return surf, title_text, norm_text

ani = animation.FuncAnimation(fig, animate, frames=350, interval=1, blit=True, repeat=True)

plt.show()
