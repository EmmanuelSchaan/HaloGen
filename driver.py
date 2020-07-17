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

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import pn_2d
reload(pn_2d)
from pn_2d import *


##################################################################################

# Basic functions and parameters, for background and fluctuations
u = UnivPlanck15()


##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=False)
#massFunc = MassFuncST(u, save=False)
massFunc = MassFuncTinker(u, save=False)


##################################################################################

# 3d profiles for matter density, Compton y, and HOD for CIB galaxies
profNFW = ProfNFW(u)
profY = ProfY(u)
profHODPenin12_353 = ProfHODPenin12(u, massFunc, nu=353) # for CIB

# plot the Fourier transform of any profile:
# here, the NFW profile
profNFW.plotU()

##################################################################################

# halo model integrals
iHaloModel = IHaloModel(u, massFunc)


##################################################################################

# 3d power spectra, auto and cross
# of matter density, Compton y, and CIB galaxies
p3d_d = P3dAuto(u, iHaloModel, profNFW, doT=False, save=False)
p3d_y = P3dAuto(u, iHaloModel, profY, doT=False, save=False)
p3d_dy = P3dCross(u, iHaloModel, profNFW, profY, doT=False, save=False)
p3d_galPenin12_353 = P3dAuto(u, iHaloModel, profHODPenin12_353, fPnoise=profHODPenin12_353.fPshotNoise, fTnoise=profHODPenin12_353.fTshotNoise, name="galpenin353", doT=False, save=False)

# plot any 3d power spectrum:
# here, compare halo model matter power spectrum
# to linear power spectrum
p3d_d.plotP()

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

# plot any projection kernel:
# here, the galaxy lensing kernel, and the CMB lensing kernel
w_gallens.plotW()
w_cmblens.plotW()

##################################################################################
# 2d power spectra, in the Limber and flat sky approximations

# auto
p2d_y = P2dAuto(u, p3d_y, w_y, doT=False, nProc=3, save=False)
p2d_cmblens = P2dAuto(u, p3d_d, w_cmblens, doT=False, nProc=3, save=False)
p2d_gallens = P2dAuto(u, p3d_d, w_gallens, doT=False, nProc=3, save=False)
p2d_cmass = P2dAuto(u, p3d_d, w_cmass, fPnoise=lambda l:1./w_cmass.ngal, doT=False, nProc=3, save=False)
p2d_wise = P2dAuto(u, p3d_d, w_wise, fPnoise=lambda l:1./w_wise.ngal, doT=False, nProc=3, save=False)
p2d_lsstgold = P2dAuto(u, p3d_d, w_lsstgold, fPnoise=lambda l:1./w_lsstgold.ngal, doT=False, nProc=3, save=False)
p2d_qso = P2dAuto(u, p3d_d, w_qso, fPnoise=lambda l:1./w_qso.ngal, doT=False, nProc=3, save=False)
p2d_cib353 = P2dAuto(u, p3d_galPenin12_353, w_cib353, fPnoise=w_cib353.fPshotNoise, fTnoise=w_cib353.fTshotNoise, doT=False, save=False, nProc=3)
# cross
p2d_cmblensgallens = P2dCross(u, p3d_d, w_cmblens, w_gallens, doT=False, nProc=3, save=False)
p2d_cmasscmblens = P2dCross(u, p3d_d, w_cmass, w_cmblens, doT=False, nProc=3, save=False)
p2d_cmassgallens = P2dCross(u, p3d_d, w_cmass, w_gallens, doT=False, nProc=3, save=False)
p2d_wisegallens = P2dCross(u, p3d_d, w_wise, w_gallens, doT=False, nProc=3, save=False)
p2d_wisecmblens = P2dCross(u, p3d_d, w_wise, w_cmblens, doT=False, nProc=3, save=False)
p2d_lsstgoldcmblens = P2dCross(u, p3d_d, w_lsstgold, w_cmblens, doT=False, nProc=3, save=False)

# plot any 2d power spectrum:
# here, the CMB lensing power spectrum
p2d_cmblens.plotP()
# here, the CIB power spectrum at 353GHz
p2d_cib353.plotP()






