import matplotlib.pylab as plt
import numpy as np

x0,x1,x2,dx,theta = np.genfromtxt('v.log', delimiter='\t', unpack=True, skip_footer=1)
n = len(x0)
dt = 0.1
t = dt * np.linspace(0,n-1,n)
window = 250

smooth = np.convolve(theta, np.ones(window)/window)[(window-1):]

mask = (t < dt*(n-window))

plt.plot(t, x0, linewidth=3, color='k')
plt.xlabel(r'$t$',fontsize=28)
plt.ylabel(r'$x_{tj}$',fontsize=28)
plt.savefig("trijunction.png", dpi=600, bbox_inches='tight')
plt.close()

plt.plot(t, theta, linewidth=5, color='k')
plt.plot(t[mask], smooth[mask], linewidth=3, color='r')
plt.xlabel(r'$t$',fontsize=28)
plt.ylabel(r'$\theta$',fontsize=28)
plt.savefig("theta.png", dpi=600, bbox_inches='tight')
plt.close()
