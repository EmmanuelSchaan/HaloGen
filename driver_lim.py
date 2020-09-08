import universe
reload(universe)
from universe import *

import mass_function
reload(mass_function)
from mass_function import *

import sfr
reload(sfr)
from sfr import *

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
# EGG: many lines

profLIMEGG = ProfLIMEGG(u, sfr, lineName='halpha', trunc=4., unit='dI/I')
#profLIMEGG = ProfLIM(u, sfr, lineName='halpha', trunc=4., unit='Jy/sr')

#profLIMEGG.plotU(z=1.)
#profLIMEGG.plotNgal()
#profLIMEGG.plotnGal()
#profLIMEGG.plotMeanIntensity()


##################################################################################
# Sobral+12: H alpha

profLIMSobral12 = ProfLIMSobral12(u, sfr, trunc=4., unit='dI/I')

#profLIMSobral12.plotLF()
#profLIMSobral12.plotnGal(profs=[profLIMSobral12, profLIMEGG])
#profLIMSobral12.plotMeanIntensity(profs=[profLIMSobral12, profLIMEGG])
profLIMSobral12.plotShotNoiseUncertainty(compareProfs=[profLIMEGG])


##################################################################################

# halo model integrals
iHaloModel = IHaloModel(u, massFunc)

##################################################################################

'''
# All lines available
#lineNames = np.array(['c2_157', 'n2_205', 'c1_609', 'co10', 'co21', 'co32', 'co43',
#           'co54', 'co65', 'co76', 'halpha', 'hbeta', 'hgamma', 'hdelta',
#           'n2_6583', 'n2_6548', 'o3_5007', 'o3_4959', 'o2_3727',
#           'lyalpha'])

# SPHEREx lines
lineNames = np.array(['halpha', 'hbeta', 'o3_5007', 'lyalpha'])


p3d_limegg = P3dAuto(u, iHaloModel, profLIMEGG, fPnoise = lambda k,z: profLIMEGG.Pshot(z), doT=False, save=False)
p3d_limsobral12 = P3dAuto(u, iHaloModel, profLIMSobral12, fPnoise = lambda k,z: profLIMSobral12.Pshot(z), doT=False, save=True)

#p3d_lim = P3dAuto(u, iHaloModel, profLIM, fPnoise = lambda k,z: 1./profLIM.nGal(z), doT=False, save=False)
#p3d_limlim = P3dCross(u, iHaloModel, profLIM, profLIM, doT=False, save=False)





p3d_limegg.plotP(z=1.)
#p3d_limegg.plotBEff()

#profLIMEGG.plotP3dGong17(p3d_limegg.fPtot, lineName='halpha')
#profLIMEGG.plotP3dGong17(p3d_limegg.fPtotinterp, lineName='h_alpha')
'''
