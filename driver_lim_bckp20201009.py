import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *

import sfr
reload(sfr)
from sfr import *

import luminosity_function
reload(luminosity_function)
from luminosity_function import *

import profile
reload(profile)
from profile import *

import i_halo_model
reload(i_halo_model)
from i_halo_model import *

import pn_3d
reload(pn_3d)
from pn_3d import *


##################################################################################

# Basic functions and parameters, for background and fluctuations
u = UnivPlanck15()


##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=False)
massFunc = MassFuncST(u, save=False)
#massFunc = MassFuncTinker(u, save=False)


##################################################################################

sfr = Sfr(u, massFunc)
#sfr.testSfrd()


##################################################################################
##################################################################################
# Luminosity functions

# SPHEREx lines
#lineNames = np.array(['halpha', 'hbeta', 'o3_5007', 'lyalpha'])

lfHaSobral12 = LFHaSobral12(u)
lfHaColbert13 = LFHaColbert13(u)
lfOiiiColbert13 = LFOiiiColbert13(u)
lfHaCochrane17 = LFHaCochrane17(u)
lfOiiiMehta15 = LFOiiiMehta15(u)


# EGG: many lines
#lineNames = np.array(['c2_157', 'n2_205', 'c1_609', 'co10', 'co21', 'co32', 'co43',
#           'co54', 'co65', 'co76', 'halpha', 'hbeta', 'hgamma', 'hdelta',
#           'n2_6583', 'n2_6548', 'o3_5007', 'o3_4959', 'o2_3727',
#           'lyalpha'])
lfHaEgg = LFEGG(u, lineName='halpha')
lfOiii5007Egg = LFEGG(u, lineName='o3_5007')
lfOiii4959Egg = LFEGG(u, lineName='o3_4959')


# lf = lfHaSobral12
#lf.plotnGal()
#lf.plotMeanIntensity()
#lf.plotLf()
#lf.plotnGalEff()
#lf.plotShotNoise()
#lf.plotS2()

#print lfOiiiMehta15.computeCorrCoeffHaOIII()
#lfOiiiMehta15.plotBivariateLf()



##################################################################################
##################################################################################
# LIM from LF

profLIMLFHaSobral12 = ProfLIMLF(u, sfr, lfHaSobral12, trunc=4.)
profLIMLFOiiiMehta15 = ProfLIMLF(u, sfr, lfOiiiMehta15, trunc=4.)
profLIMLFHaColbert13 = ProfLIMLF(u, sfr, lfHaColbert13, trunc=4.)
profLIMLFOiiiColbert13 = ProfLIMLF(u, sfr, lfOiiiColbert13, trunc=4.)
profLIMLFHaCochrane17 = ProfLIMLF(u, sfr, lfHaCochrane17, trunc=4.)


# EGG
profLimLfHaEgg = ProfLIMLF(u, sfr, lfHaEgg, trunc=4.)
profLimLfOiii5007Egg = ProfLIMLF(u, sfr, lfOiii5007Egg, trunc=4.)
profLimLfOiii4959Egg = ProfLIMLF(u, sfr, lfOiii4959Egg, trunc=4.)


##################################################################################
##################################################################################

# halo model integrals
iHaloModel = IHaloModel(u, massFunc)


##################################################################################
##################################################################################
# Power spectrum


# Sobral+12 Ha
p3d_limsobral12 = P3dAuto(u, iHaloModel, profLIMLFHaSobral12, fPnoise = lambda k,z: profLIMLFSobral12.Pshot(z), doT=False, save=True)
# Mehta+15 OIII
p3d_limmehta15 = P3dAuto(u, iHaloModel, profLIMLFOiiiMehta15, fPnoise = lambda k,z: profLIMLFMehta15.Pshot(z), doT=False, save=True)
p3d_limhacolbert13 = P3dAuto(u, iHaloModel, profLIMLFHaColbert13, fPnoise = lambda k,z: profLIMLFHaColbert13.Pshot(z), doT=False, save=True)
p3d_limoiiicolbert13 = P3dAuto(u, iHaloModel, profLIMLFOiiiColbert13, fPnoise = lambda k,z: profLIMLFOiiiColbert13.Pshot(z), doT=False, save=True)
p3d_limhacochrane17 = P3dAuto(u, iHaloModel, profLIMLFHaCochrane17, fPnoise = lambda k,z: profLIMLFHaCochrane17.Pshot(z), doT=False, save=True)


# EGG
p3d_limhaegg = P3dAuto(u, iHaloModel, profLimLfHaEgg, fPnoise = lambda k,z: profLimLfHaEgg.Pshot(z), doT=False, save=True)
p3d_limoiii5007egg = P3dAuto(u, iHaloModel, profLimLfOiii5007Egg, fPnoise = lambda k,z: profLimLfOiii5007Egg.Pshot(z), doT=False, save=True)
p3d_limoiii4959egg = P3dAuto(u, iHaloModel, profLimLfOiii4959Egg, fPnoise = lambda k,z: profLimLfOiii4959Egg.Pshot(z), doT=False, save=True)


##p3d_limhaegg.plotP(z=1.)
##p3d_limhaegg.plotBEff()

#p3d_limoiii5007egg.plotP(z=1.)
#p3d_limoiii4959egg.plotP(z=1.)

##profLimLfHaEgg.plotP3dGong17(p3d_limegg.fPtotinterp, lineName='halpha')




