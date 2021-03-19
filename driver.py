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

import p3d
reload(p3d)
from p3d import *

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import p2d
reload(p2d)
from p2d import *


##################################################################################

# Basic functions and parameters, for background and fluctuations
u = UnivPlanck15()


##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=False)
#massFunc = MassFuncST(u, save=False)
massFunc = MassFuncTinker(u, save=False)

'''
massFunc.testInterp(z=0.)
massFunc.plotMassFunc()
massFunc.plotIntegralConstraintsZ()
massFunc.plotMassCounterTerms()
massFunc.plotMassConstraintMMin()
massFunc.plotBiasConstraintMMin()
'''

##################################################################################

# 3d profiles for matter density, Compton y, and HOD for CIB galaxies
profNFW = ProfNFW(u)
profY = ProfY(u)
profHODPenin12_353 = ProfHODPenin12(u, massFunc, nu=353) # for CIB

# plot the Fourier transform of any profile:
# here, the NFW profile
'''
profNFW.plotU()
'''

##################################################################################
# 3d power spectra, auto and cross

import p3d
reload(p3d)
from p3d import *

# Matter density, Compton y, and CIB galaxies
p3d_d = P3dAuto(u, massFunc, profNFW, save=False)
p3d_y = P3dAuto(u, massFunc, profY, save=False)
p3d_dy = P3dCross(u, massFunc, profNFW, profY, save=False)
p3d_galPenin12_353 = P3dAuto(u, massFunc, profHODPenin12_353, pNoise=profHODPenin12_353.fPshotNoise, name="galpenin353", save=False)

'''
p3d_d.plotPInt(z=0.)
p3d_d.plotP1hMMinDependence(k=0.01)
p3d_d.plotP2hMMinDependence(k=0.01)

p3d_y.plotPInt(z=0.)
p3d_y.plotP1hMMinDependence(k=0.01)
p3d_y.plotP2hMMinDependence(k=0.01)

p3d_dy.plotPInt(z=0.)
p3d_dy.plotP1hMMinDependence(k=0.01)
p3d_dy.plotP2hMMinDependence(k=0.01)

p3d_galPenin12_353.plotPInt(z=0.)
p3d_galPenin12_353.plotP1hMMinDependence(k=0.01)
p3d_galPenin12_353.plotP2hMMinDependence(k=0.01)
'''


##################################################################################

# projection kernels to get 2d power spectra
# from 3d power spectra,
# for various tracers, for lensing convergence, for tSZ
w_cmblens = WeightLensSingle(u, z_source=1100., name="cmblens")
w_gallens = WeightLensSingle(u, z_source=1., name="gallens")
w_y = WeightY(u)
w_cmass = WeightTracerCMASS(u)
w_wise = WeightTracerWISE(u)
w_cib353 = WeightCIBPenin12(u, nu=353.e9, fluxCut=315.e-3, name='cibpenin12')
w_qso = WeightTracerDESIQSO(u)
w_lsstgold = WeightTracerLSSTGold(u)

'''
# plot any projection kernel:
# here, the galaxy lensing kernel, and the CMB lensing kernel
w_gallens.plotW()
w_cmblens.plotW()
'''

##################################################################################
# 2d power spectra, in the Limber and flat sky approximations

import p2d
reload(p2d)
from p2d import *


# auto
p2d_cmblens = P2dAuto(u, p3d_d, w_cmblens, nProc=3, save=False)
p2d_y = P2dAuto(u, p3d_y, w_y, nProc=3, save=False)
p2d_cmblens = P2dAuto(u, p3d_d, w_cmblens, nProc=3, save=False)
p2d_gallens = P2dAuto(u, p3d_d, w_gallens, nProc=3, save=False)
p2d_cmass = P2dAuto(u, p3d_d, w_cmass, pNoise=lambda l:1./w_cmass.ngal, nProc=3, save=False)
p2d_wise = P2dAuto(u, p3d_d, w_wise, pNoise=lambda l:1./w_wise.ngal, nProc=3, save=False)
p2d_lsstgold = P2dAuto(u, p3d_d, w_lsstgold, pNoise=lambda l:1./w_lsstgold.ngal, nProc=3, save=False)
p2d_qso = P2dAuto(u, p3d_d, w_qso, pNoise=lambda l:1./w_qso.ngal, nProc=3, save=False)
p2d_cib353 = P2dAuto(u, p3d_galPenin12_353, w_cib353, pNoise=w_cib353.fPshotNoise, save=False, nProc=3)
# cross
p2d_cmblensgallens = P2dCross(u, p3d_d, w_cmblens, w_gallens, nProc=3, save=False)
p2d_cmasscmblens = P2dCross(u, p3d_d, w_cmass, w_cmblens, nProc=3, save=False)
p2d_cmassgallens = P2dCross(u, p3d_d, w_cmass, w_gallens, nProc=3, save=False)
p2d_wisegallens = P2dCross(u, p3d_d, w_wise, w_gallens, nProc=3, save=False)
p2d_wisecmblens = P2dCross(u, p3d_d, w_wise, w_cmblens, nProc=3, save=False)
p2d_lsstgoldcmblens = P2dCross(u, p3d_d, w_lsstgold, w_cmblens, nProc=3, save=False)


# plot any 2d power spectrum:
# here, the CMB lensing power spectrum
p2d_cmblens.plotP()
# here, the CIB power spectrum at 353GHz
p2d_cib353.plotP()






