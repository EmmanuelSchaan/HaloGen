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
z = np.array([1., 1.5, 5.])
w0 = WeightLensSingle(u, z_source=z[0], name="z1lens", nameLatex=r'$\kappa_{\text{gal } z=1}$')
w1 = WeightLensSingle(u, z_source=z[1], name="z1p5lens", nameLatex=r'$\kappa_{\text{gal } z=1.5}$')
w2 = WeightLensSingle(u, z_source=z[2], name="limz5lens", nameLatex=r'$\kappa_{\text{LIM } z=5}$')

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



