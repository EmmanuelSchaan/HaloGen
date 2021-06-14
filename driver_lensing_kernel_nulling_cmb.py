import universe
reload(universe)
from universe import *

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import p2d
reload(p2d)
from p2d import *


##################################################################################

u = UnivPlanck15()

##################################################################################
# Nulling low-z signal for single-source planes

# single source planes
z = np.array([5., 6., 1100.])
w0 = WeightLensSingle(u, z_source=z[0], name="limz5lens", nameLatex=r'$\kappa_{\text{LIM } z=5}$')
w1 = WeightLensSingle(u, z_source=z[1], name="limz6lens", nameLatex=r'$\kappa_{\text{LIM } z=6}$')
w2 = WeightLensSingle(u, z_source=z[2], name="cmblens", nameLatex=r'$\kappa_\text{CMB}$')

# linear combination to null the low-z signal
w_combined = WeightLensSingle(u, z_source=z[2], name="nulllens", nameLatex=r'$\kappa_\text{Null}$')
chi = u.bg.comoving_distance(z)
alpha = (1./chi[2] - 1./chi[0]) / (1./chi[0] - 1./chi[1])
w_combined.f = lambda a:  w2.f(a) + alpha * w1.f(a) - (1.+alpha) * w0.f(a)

# Quick plot
#w_combined.plotW(ws=[w0, w1, w2, w_combined], ls=['--', '--', '--', '-'], zMax=1100.)


##################################################################################
# Plots

Na = 501
A = np.linspace(1.-1.e-5, 1./(1.+1100.), Na)
Z = 1./A - 1.
ComovDistToObs = np.array( map(u.bg.comoving_distance, Z ) )
H_A = u.hubble(1./A-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)


# W(chi)
fig=plt.figure(-1)
ax=plt.subplot(111)
#
# CMB 
W = np.array( map( lambda a: w2.f(a), A ) )
ax.plot(ComovDistToObs, W/H_A, lw=1, ls='--', label=w2.nameLatex)
#
# LIM z=6
W = np.array( map( lambda a: w1.f(a), A ) )
ax.plot(ComovDistToObs, W/H_A, lw=1, ls='--', label=w1.nameLatex)
#
# LIM z=5
W = np.array( map( lambda a: w0.f(a), A ) )
ax.plot(ComovDistToObs, W/H_A, lw=1, ls='--', label=w0.nameLatex)
#
# Null combination
W = np.array( map( lambda a: w_combined.f(a), A ) )
#ax.plot(ComovDistToObs, W/H_A, lw=2, ls='-', label=w_combined.nameLatex)
ax.fill_between(ComovDistToObs, 0.*W/H_A, W/H_A, facecolor='r', edgecolor='', label=w_combined.nameLatex)
#
#ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$\chi$ [Mpc/$h$]', fontsize=22)
ax.set_ylabel(r'$W(z)$', fontsize=22)
ax.set_xlim((0., u.bg.comoving_distance(1100.)))
ax.set_ylim((0., 0.5))
#
#fig.savefig("./figures/weight/W_cmblens.pdf")
plt.show()




##################################################################################
# What about source bins with wider redshift distribution?


##################################################################################

save = False

p2d_lim1lens = P2dAuto(u, u, w0, nProc=3, name='halofit', save=save)
p2d_lim2lens = P2dAuto(u, u, w1, nProc=3, name='halofit', save=save)
p2d_cmblens = P2dAuto(u, u, w2, nProc=3, name='halofit', save=save)
p2d_nulllens = P2dAuto(u, u, w_combined, nProc=3, name='halofit', save=save)

p2d_lim1cmblens = P2dCross(u, u, w2, w0, nProc=3, name='halofit', save=save)
p2d_lim2cmblens = P2dCross(u, u, w2, w1, nProc=3, name='halofit', save=save)
p2d_nulllenscmblens = P2dCross(u, u, w2, w_combined, nProc=3, name='halofit', save=save)


##################################################################################

'''
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(p2d_cmblens.L, p2d_cmblens.Ptot, label=r'$\kappa_\text{CMB} \kappa_\text{CMB}$')
ax.plot(p2d_nulllenscmblens.L, p2d_nulllenscmblens.Ptot, label=r'$\kappa_\text{CMB} \kappa_\text{Null}$')
ax.plot(p2d_nulllens.L, p2d_nulllens.Ptot, label=r'$\kappa_\text{Null} \kappa_\text{Null}$')
#
ax.legend(loc=1)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'L')
ax.set_ylabel(r'$C_L$')

plt.show()


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(p2d_cmblens.L, p2d_nulllenscmblens.Ptot / p2d_cmblens.Ptot)
#ax.plot(p2d_cmblens.L, p2d_nulllens.Ptot / p2d_cmblens.Ptot)
#
ax.set_ylim((1.e-2, 1.))
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'L')
ax.set_ylabel(r'$C_L^{\kappa_\text{CMB} \kappa_\text{Null}} / C_L^{\kappa_\text{CMB}}$')

plt.show()
'''

