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

profLIM = ProfLIM(u, sfr, lineName='halpha', trunc=4.)

#profLIM.plotU(z=1.)
#profLIM.plotNgal()
#profLIM.plotnGal()
profLIM.plotMeanIntensity()

##################################################################################

# halo model integrals
iHaloModel = IHaloModel(u, massFunc)

##################################################################################


#p3d_lim = P3dAuto(u, iHaloModel, profLIM, fPnoise = lambda k,z: profLIM.Pshot(z), doT=False, save=False)
p3d_lim = P3dAuto(u, iHaloModel, profLIM, fPnoise = lambda k,z: 1./profLIM.nGal(z), doT=False, save=False)
#p3d_limlim = P3dCross(u, iHaloModel, profLIM, profLIM, doT=False, save=False)

#p3d_lim.plotP(z=1.)
#p3d_lim.plotBEff()

