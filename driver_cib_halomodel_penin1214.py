import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *

import profile
reload(profile)
from profile import *

import i_halo_model
reload(i_halo_model)
from i_halo_model import *

import pn_3d
reload(pn_3d)
from pn_3d import *

import weight
reload(weight)
from weight import *

import pn_2d
reload(pn_2d)
from pn_2d import *

import hsvpn_2d
reload(hsvpn_2d)
from hsvpn_2d import *

import covp_2d
reload(covp_2d)
from covp_2d import *

import cmb
reload(cmb)
from cmb import *

import cmb_lensing_rec
reload(cmb_lensing_rec)
from cmb_lensing_rec import *


##################################################################################

u = UnivPlanck15()
#massFunc = MassFuncPS(u, save=False)
#massFunc = MassFuncST(u, save=False)
massFunc = MassFuncTinker(u, save=False)
iHaloModel = IHaloModel(u, massFunc)

##################################################################################
##################################################################################
# CIB halo model from Penin+12, 14, using Bethermin+12 flux counts
# HOD from Penin+12,14, from Tinker Wetzel 10
# Penin+12,14 halo model for CIB
# No: flux cuts for Planck, from table 2 in Penin+12
# Yes: flux cuts in table 1 of Planck13 XXX. These two don't match!
'''
# 217GHz
profHODPenin12_217 = ProfHODPenin12(u, massFunc, nu=217)
p3d_galPenin12_217 = P3dAuto(u, iHaloModel, profHODPenin12_217, fPnoise=profHODPenin12_217.fPshotNoise, fTnoise=profHODPenin12_217.fTshotNoise, name="galpenin217", doT=True, save=False)
w_cib217 = WeightCIBPenin12(u, nu=217.e9, fluxCut=225.e-3, name='cibpenin12')
p2d_cib217 = P2dAuto(u, p3d_galPenin12_217, w_cib217, fPnoise=w_cib217.fPshotNoise, fTnoise=w_cib217.fTshotNoise, doT=True, save=False, nProc=1)
'''

# 353GHz
profHODPenin12_353 = ProfHODPenin12(u, massFunc, nu=353)
p3d_galPenin12_353 = P3dAuto(u, iHaloModel, profHODPenin12_353, fPnoise=profHODPenin12_353.fPshotNoise, fTnoise=profHODPenin12_353.fTshotNoise, name="galpenin353", doT=True, save=False)
w_cib353 = WeightCIBPenin12(u, nu=353.e9, fluxCut=315.e-3, name='cibpenin12')
p2d_cib353 = P2dAuto(u, p3d_galPenin12_353, w_cib353, fPnoise=w_cib353.fPshotNoise, fTnoise=w_cib353.fTshotNoise, doT=True, save=False, nProc=3)


# 545GHz
profHODPenin12_545 = ProfHODPenin12(u, massFunc, nu=545)
p3d_galPenin12_545 = P3dAuto(u, iHaloModel, profHODPenin12_545, fPnoise=profHODPenin12_545.fPshotNoise, fTnoise=profHODPenin12_545.fTshotNoise, name="galpenin545", doT=True, save=False)
w_cib545 = WeightCIBPenin12(u, nu=545.e9, fluxCut=350.e-3, name='cibpenin12')
#w_cib545.aMax = 1./(1.+1.)
p2d_cib545 = P2dAuto(u, p3d_galPenin12_545, w_cib545, fPnoise=w_cib545.fPshotNoise, fTnoise=w_cib545.fTshotNoise, doT=True, save=False, nProc=3)

'''
# 857GHz
profHODPenin12_857 = ProfHODPenin12(u, massFunc, nu=857)
p3d_galPenin12_857 = P3dAuto(u, iHaloModel, profHODPenin12_857, fPnoise=profHODPenin12_857.fPshotNoise, fTnoise=profHODPenin12_857.fTshotNoise, name="galpenin857", doT=True, save=False)
w_cib857 = WeightCIBPenin12(u, nu=857.e9, fluxCut=710.e-3, name='cibpenin12')
p2d_cib857 = P2dAuto(u, p3d_galPenin12_857, w_cib857, fPnoise=w_cib857.fPshotNoise, fTnoise=w_cib857.fTshotNoise, doT=True, save=False, nProc=1)
'''


##################################################################################
# CIB monopole z-dependence at different frequencies,
# for Penin+12, 14, using Bethermin+12 flux counts
'''
# inverse hubble length: H/c in (h Mpc^-1)
# used to convert the kernel from chi to z
H = u.Hubble(w_cib217.A) / 3.e5
# projection kernel
W_cib217 = np.array(map(w_cib217.f, w_cib217.A))
W_cib353 = np.array(map(w_cib353.f, w_cib353.A))
W_cib545 = np.array(map(w_cib545.f, w_cib545.A))
W_cib857 = np.array(map(w_cib857.f, w_cib857.A))

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(w_cib217.Z, (W_cib217/H) / np.max(W_cib217/H), c=plt.cm.winter(0.), lw=2, label=r'217 GHz')
ax.plot(w_cib353.Z, (W_cib353/H) / np.max(W_cib353/H), c=plt.cm.winter(0.33), lw=2, label=r'353 GHz')
ax.plot(w_cib545.Z, (W_cib545/H) / np.max(W_cib545/H), c=plt.cm.winter(0.66), lw=2, label=r'545 GHz')
ax.plot(w_cib857.Z, (W_cib857/H) / np.max(W_cib857/H), c=plt.cm.winter(1.), lw=2, label=r'857 GHz')
#
ax.legend(loc=1)
ax.set_xlim((0., 6.))
ax.set_ylim((0., 1.1))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'CIB intensity kernel, arbitrary units')
#ax.set_xscale('log', nonposx='clip')
#ax.set_yscale('log', nonposy='clip')
#
#path="./figures/cib_penin12/emissivity_kernels_penin12.pdf"
#fig.savefig(path, bbox_inches='tight')

plt.show()
'''

##################################################################################
# CIB power spectrum z-dependence
# for Penin+12, 14, using Bethermin+12 flux counts
'''
# values of a/z to plot
Na = 101
A = np.linspace(1./(1.+5.1), 1./(1.+1.e-5), Na)
Z = 1./A - 1.

# CIB source distribution: power spectrum
l = 100.
f = lambda a: p2d_cib545.integrand(a, p2d_cib545.Pn.fP, l) * a**2
dcldz_cibpenin14_545_l100 = np.array(map(f, A))
#
l = 500.
f = lambda a: p2d_cib545.integrand(a, p2d_cib545.Pn.fP, l) * a**2
dcldz_cibpenin14_545_l500 = np.array(map(f, A))
#
l = 2000.
f = lambda a: p2d_cib545.integrand(a, p2d_cib545.Pn.fP, l) * a**2
dcldz_cibpenin14_545_l2000 = np.array(map(f, A))
#
l = 10000.
f = lambda a: p2d_cib545.integrand(a, p2d_cib545.Pn.fP, l) * a**2
dcldz_cibpenin14_545_l10000 = np.array(map(f, A))


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CIB source distributions: power spectrum
ax.plot(Z, dcldz_cibpenin14_545_l500/np.max(dcldz_cibpenin14_545_l500), c=plt.cm.autumn(0.), lw=1, label=r'Penin+14: $C^0_{\ell=500}$')
ax.plot(Z, dcldz_cibpenin14_545_l2000/np.max(dcldz_cibpenin14_545_l2000), c=plt.cm.autumn(0.25), lw=1, label=r'Penin+14: $C^0_{\ell=2000}$')
ax.plot(Z, dcldz_cibpenin14_545_l10000/np.max(dcldz_cibpenin14_545_l10000), c=plt.cm.autumn(0.5), lw=1, label=r'Penin+14: $C^0_{\ell=10000}$')
#
# arbitrary single source
ax.axvline(2., color='m', ymax=0.95, lw=2, ls='--', label=r'Single source')
#
ax.legend(loc=1, fontsize=14)
ax.set_xlim((0., 6.))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$d(\text{CIB}) / dz$')

plt.show()
'''

##################################################################################
# Effect of trispectrum on reconstruction noise
'''
# trispectrum for Penin+12, 14, using Bethermin+12 flux counts
t3d_galPenin12 = T3dAuto(u, iHaloModel, profHODPenin12, noiseBias=profHODPenin12.fTshotNoise)

t2d_cib217 = T2dAuto(u, t3d_galPenin12, w_cib217, noiseBias=w_cib217.fTshotNoise, save=False)
t2d_cib353 = T2dAuto(u, t3d_galPenin12, w_cib353, noiseBias=w_cib353.fTshotNoise, save=False)
t2d_cib545 = T2dAuto(u, t3d_galPenin12, w_cib545, noiseBias=w_cib545.fTshotNoise, save=False)
t2d_cib857 = T2dAuto(u, t3d_galPenin12, w_cib857, noiseBias=w_cib857.fTshotNoise, save=False)

cib = CIB(beam=5., noise=47., nu1=217.e9, nu2=217.e9)
cib.lMaxT = 1.e4
cib.funlensedTT = p2d_cib217.fPinterp
cib.ftotalTT = lambda l: p2d_cib217.fPinterp(l) + cib.fdetectorNoise(l)

cibLensRec = CMBLensRec(cib, save=False)
Result = cibLensRec.trispectrumCorrection(t2d_cib217.fTnondiag, L=20.)
#cibLensRec.plotSnrDensity_TT()
'''

##################################################################################
##################################################################################
# CIB monopole z-dependence: comparison Penin vs world
'''
# assuming CIB is a single source
w_ciblens = WeightLensSingle(u, z_source=2., name="ciblens")

# z-dist for CIB monopole, Schmidt+15
w_ciblensschmidt15_353 = WeightLensCIBSchmidt15(u, z0=1.25, alpha=1.24, name="ciblensSchmidt15_353")
w_ciblensschmidt15_545 = WeightLensCIBSchmidt15(u, z0=1.34, alpha=1.50, name="ciblensSchmidt15_545")
w_ciblensschmidt15_857 = WeightLensCIBSchmidt15(u, z0=1.08, alpha=1.52, name="ciblensSchmidt15_857")

# z-dist for CIB monopole, Pullen+17
w_ciblenspullen17_353 = WeightLensCIBPullen17(u, nu=353, name="ciblensPullent17_353")
w_ciblenspullen17_545 = WeightLensCIBPullen17(u, nu=545, name="ciblensPullent17_545")
w_ciblenspullen17_857 = WeightLensCIBPullen17(u, nu=857, name="ciblensPullent17_857")

# kernels from Penin+12, 14, using Bethermin+12 flux counts
w_cib353 = WeightCIBPenin12(u, nu=353.e9, fluxCut=315.e-3, name="cibpenin12_353")
w_cib545 = WeightCIBPenin12(u, nu=545.e9, fluxCut=350.e-3, name="cibpenin12_545")
w_cib857 = WeightCIBPenin12(u, nu=857.e9, fluxCut=710.e-3, name="cibpenin12_857")

# kappa_CIB kernel
W_single = np.array( map( lambda a: w_ciblens.f(a), A ) )
W_single /= H_A
W_single *= (Z<=w_ciblens.z_source)
W_ciblensschmidt15_545 = np.array( map( lambda a: w_ciblensschmidt15_545.f(a), A ) )
W_ciblensschmidt15_545 /= H_A
W_ciblenspullen17_545 = np.array( map( lambda a: w_ciblenspullen17_545.f(a), A ) )
W_ciblenspullen17_545 /= H_A

# values of a/z to plot
Na = 101
A = np.linspace(1./(1.+5.1), 1./(1.+1.e-5), Na)
Z = 1./A - 1.

# to convert d/dchi to d/dz
H_A = u.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)

# CIB source distribution: monopole
dmonodz_cibschmidt15_545 = np.array(map(w_ciblensschmidt15_353.fdpdz, Z))
#
dmonodz_cibpullen17_545 = np.array(map(w_ciblenspullen17_353.fdpdz, Z))
#
dmonodz_cibpenin14_545 = np.array(map(w_cib545.f, A))
dmonodz_cibpenin14_545 /= H_A


# plot CIB source distributions: monopole
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(Z, dmonodz_cibschmidt15_545/np.max(dmonodz_cibschmidt15_545), c=plt.cm.winter(0.), lw=1, label=r'Schmidt+15: monopole')
ax.plot(Z, dmonodz_cibpenin14_545/np.max(dmonodz_cibpenin14_545), c=plt.cm.winter(0.25), lw=1, label=r'Bethermin+12: monopole')
ax.plot(Z, dmonodz_cibpullen17_545/np.max(dmonodz_cibpullen17_545), c=plt.cm.winter(0.5), lw=1, label=r'Pullen+17: monopole')
#
# arbitrary single source
ax.axvline(2., color='m', ymax=0.95, lw=2, ls='--', label=r'Single source')
#
ax.legend(loc=1, fontsize=14)
ax.set_xlim((0., 6.))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$d(\text{CIB}) / dz$')

plt.show()


# CIB lensing kernel: Penin vs world
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# corresponding kappa_CIB kernels
ax.plot(Z, W_single, 'g-', lw=2, label=r'$z_\text{source} = 2$')
ax.plot(Z, W_ciblensschmidt15_545, 'r-', lw=2, label=r'Schmidt+15 monopole')
ax.plot(Z, W_ciblenspullen17_545, 'b-', lw=2, label=r'Pullen+17 monopole')
#
ax.legend(loc=1, fontsize=14)
ax.set_xlim((0., 6.))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$W_{\kappa_\text{CIB}}$')

plt.show()
'''

##################################################################################
##################################################################################
# CIB trispectrum from the halo model and data
'''
# read the measured CIB trispectrum from Planck GNILC maps at 545GHz
data = np.genfromtxt("./input/cib_gnilc_545ghz_fullsky/trispectrum_gnilc_545ghz_collapsed_fullsky.txt")
L = data[:,0]
Tmeas = data[:,1]
sTmeas = data[:,2]

factor = 1.#p2d_cib545.L*(p2d_cib545.L+1.)/(2.*np.pi)

# T
fig=plt.figure(0)
ax=plt.subplot(111)
#
ax.plot(p2d_cib545.L, factor*p2d_cib545.Ptot**2, 'k--', lw=3, label=r'$C^2$')
#plot(p2d_cib545.L, factor*p2d_cib545.Ptot**2 * 4. * 7.5e-8, 'k--', lw=2, label=r'$C^2 / N_\text{modes}$')
#
#ax.fill_between(p2d_cib545.L, 1.e-10, 21, facecolor='b', alpha=0.05, label=r'T upper bound')
ax.plot(p2d_cib545.L, factor*p2d_cib545.Ttot, 'k', lw=3, label=r'$\mathcal{T}^\text{total}$')
ax.plot(p2d_cib545.L, factor*p2d_cib545.T1h, 'r', lw=1.5, label=r'$\mathcal{T}^{1h}$')
ax.plot(p2d_cib545.L, factor*p2d_cib545.T2h, 'orange', lw=1.5, label=r'$\mathcal{T}^{2h}$')
ax.plot(p2d_cib545.L, factor*p2d_cib545.T4h, 'gold', lw=1.5, label=r'$\mathcal{T}^{4h}$')
ax.plot(p2d_cib545.L, factor*p2d_cib545.Tssv, 'limegreen', lw=1.5, label=r'$\mathcal{T}^\text{SSV}$')
ax.plot(p2d_cib545.L, factor*p2d_cib545.Tnoise, 'royalblue', lw=1.5, label=r'$\mathcal{T}^\text{shot}$')
#
# measured trispectrum
ax.scatter(L[:-2], Tmeas[:-2]+2.*sTmeas[:-2], s=30, marker='.', facecolors='k', linewidths=1.5, edgecolors='k', label=r'$\mathcal{T}$ upper bound')
ax.quiver(L[:-2], Tmeas[:-2]+2.*sTmeas[:-2], 0.*np.ones_like(L[:-2]), -0.5 * np.ones_like(L[:-2]))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_ylim((1.e-1, 1.e11))
##ax.set_ylim((1.e-2, 2.e8))
ax.set_xlabel(r'\ell')
ax.set_xlim((50., 5.e4))
#ax.set_ylabel(r'$T(\ell)$')
ax.set_ylabel(r'$\langle \mathcal{T}^0_{\boldsymbol{\ell}, \boldsymbol{L}-\boldsymbol{\ell}, \boldsymbol{\ell}, \boldsymbol{-L}-\boldsymbol{\ell}} \rangle$     [Jy$^4$/sr]')
path = "./figures/pn2d/t_"+str(p2d_cib545.Weight)+"_test.pdf"
#fig.savefig(path, bbox_inches='tight')

plt.show()
'''

##################################################################################
# Effect of zMin on the CIB trispectrum
'''
# values of zMin to try
ZMin = np.array([0., 0.5, 1., 1.5])
Color = ['k', 'dimgray', 'darkgray', 'gainsboro']

factor = 1.#p2d_cib545.L*(p2d_cib545.L+1.)/(2.*np.pi)

# T
fig=plt.figure(0)
ax=plt.subplot(111)
#
ax.plot(p2d_cib545.L, factor*p2d_cib545.Ptot**2, 'k--', lw=3, label=r'$C^2$')
#
for iZMin in range(len(ZMin)):
   zMin = ZMin[iZMin]
   w_cib545.aMax = 1./(1.+zMin)
   p2d_cib545 = P2dAuto(u, p3d_galPenin12_545, w_cib545, fPnoise=w_cib545.fPshotNoise, fTnoise=w_cib545.fTshotNoise, doT=True, save=True, nProc=3)
   #
   ax.plot(p2d_cib545.L, factor*p2d_cib545.Ttot, c=Color[iZMin], lw=2, label=r'$\mathcal{T}^\text{total}$, $z_\text{min}=$'+str(zMin))
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_ylim((1.e-1, 1.e11))
##ax.set_ylim((1.e-2, 2.e8))
ax.set_xlabel(r'\ell')
ax.set_xlim((50., 5.e4))
#ax.set_ylabel(r'$T(\ell)$')
ax.set_ylabel(r'$\langle \mathcal{T}^0_{\boldsymbol{\ell}, \boldsymbol{L}-\boldsymbol{\ell}, \boldsymbol{\ell}, \boldsymbol{-L}-\boldsymbol{\ell}} \rangle$     [Jy$^4$/sr]')
path = "./figures/pn2d/t_"+str(p2d_cib545.Weight)+"_zMinDpdce.pdf"
#fig.savefig(path, bbox_inches='tight')

plt.show()

# clean up
w_cib545.aMax = 1./(1.+0.)
p2d_cib545 = P2dAuto(u, p3d_galPenin12_545, w_cib545, fPnoise=w_cib545.fPshotNoise, fTnoise=w_cib545.fTshotNoise, doT=True, save=True, nProc=3)
'''





















