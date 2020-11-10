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

#import pn_3d
#reload(pn_3d)
#from pn_3d import *

import p3d_rsd
reload(p3d_rsd)
from p3d_rsd import *


##################################################################################

# Basic functions and parameters, for background and fluctuations
u = UnivPlanck15()

#u.plotSigma2V1d()
#u.plotSigma2DispFog()
#u.plotKMaxParaSpectroRes()
u.plotKMaxPerpPsf()
#u.plotKFPerp()
#u.plotTradeOffNModes()



##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=False)
massFunc = MassFuncST(u, save=False)
#massFunc = MassFuncTinker(u, save=False)


##################################################################################

#sfr = SfrFonseca16(u, massFunc)
#sfr = SfrMoster13(u, massFunc)

sfr = SfrMoster13Speagle14(u, massFunc, scatter=True, nProc=3, save=False)

#sfr.plotSfr()
#sfr.plotSfrd()
#sfr.plotnHEff()
#sfr.plotNHEffSpherex()
#sfr.plotBEff()
sfr.plotdbEff2dlnm()
sfr.plotdP1hdlnm()
sfr.plotdlnMeanIntensitydlnm()

##################################################################################
##################################################################################
# Luminosity functions

# SPHEREx lines
#lineNames = np.array(['halpha', 'hbeta', 'o3_5007', 'lyalpha'])

# EGG: available lines
#lineNames = np.array(['c2_157', 'n2_205', 'c1_609', 'co10', 'co21', 'co32', 'co43',
#           'co54', 'co65', 'co76', 'halpha', 'hbeta', 'hgamma', 'hdelta',
#           'n2_6583', 'n2_6548', 'o3_5007', 'o3_4959', 'o2_3727',
#           'lyalpha'])



#import pandas as pd
#df = pd.DataFrame(columns=['lineName', 'ref', 'lf', 'profLimLf', 'p3d'])
#df = df.append({'lineName': 'halpha', 'ref': 'Sobral12', 'lf': LFHaSobral12(u)}, ignore_index=True)


lfHa = {}
lfHa['Sobral12'] = LFHaSobral12(u)
lfHa['Colbert13'] = LFHaColbert13(u)
lfHa['Cochrane17'] = LFHaCochrane17(u)
#lfHa['Egg'] = LFEGG(u, lineName='halpha')

lfOiii = {}
lfOiii['Colbert13'] = LFOiiiColbert13(u)
lfOiii['Mehta15'] = LFOiiiMehta15(u)
#lfOiii['5007Egg'] = LFEGG(u, lineName='o3_5007')
##lfOiii['4959Egg'] = LFEGG(u, lineName='o3_4959')


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
# Plot: luminosity functions

#for key in lfHa.keys():
#   try: lfHa[key].plotLf()
#   except: pass
#
#for key in lfOiii.keys():
#   try: lfOiii[key].plotLf()
#   except: pass

'''
lfHa['Cochrane17'].plotLf(lfs=[lfHa[key] for key in lfHa.keys()])

lfOiii['Colbert13'].plotLf(lfs=[lfOiii[key] for key in lfOiii.keys()])
'''

##################################################################################
# Plot: mean intensity

'''
lfHa['Sobral12'].plotMeanIntensity(lfs=[lfHa[key] for key in lfHa.keys()])

lfOiii['Colbert13'].plotMeanIntensity(lfs=[lfOiii[key] for key in lfOiii.keys()])
'''

##################################################################################
# Plot: nGalEff

'''
lfHa['Sobral12'].plotnGalEff(lfs=[lfHa[key] for key in lfHa.keys()])

lfOiii['Colbert13'].plotnGalEff(lfs=[lfOiii[key] for key in lfOiii.keys()])
'''

##################################################################################
# Plot: shot noise
'''
lfHa['Sobral12'].plotShotNoise(lfs=[lfHa[key] for key in lfHa.keys()])
lfOiii['Colbert13'].plotShotNoise(lfs=[lfOiii[key] for key in lfOiii.keys()])
'''

##################################################################################
# Contributions from each luminosity
'''
lfHa['Sobral12'].plotdlnMeanIntensitydlnL()
lfHa['Sobral12'].plotdlnPshotdlnL()
'''

##################################################################################
##################################################################################
# LIM from LF

profLimLfHa = {}
for key in lfHa.keys():
   profLimLfHa[key] = ProfLIMLF(u, sfr, lfHa[key], trunc=4.)

profLimLfOiii = {}
for key in lfOiii.keys():
   profLimLfOiii[key] = ProfLIMLF(u, sfr, lfOiii[key], trunc=4.)



#profLIMLFHaSobral12 = ProfLIMLF(u, sfr, lfHaSobral12, trunc=4.)
#profLIMLFOiiiMehta15 = ProfLIMLF(u, sfr, lfOiiiMehta15, trunc=4.)
#profLIMLFHaColbert13 = ProfLIMLF(u, sfr, lfHaColbert13, trunc=4.)
#profLIMLFOiiiColbert13 = ProfLIMLF(u, sfr, lfOiiiColbert13, trunc=4.)
#profLIMLFHaCochrane17 = ProfLIMLF(u, sfr, lfHaCochrane17, trunc=4.)
#
#
## EGG
#profLimLfHaEgg = ProfLIMLF(u, sfr, lfHaEgg, trunc=4.)
#profLimLfOiii5007Egg = ProfLIMLF(u, sfr, lfOiii5007Egg, trunc=4.)
#profLimLfOiii4959Egg = ProfLIMLF(u, sfr, lfOiii4959Egg, trunc=4.)


##################################################################################
##################################################################################

# halo model integrals
iHaloModel = IHaloModel(u, massFunc)


##################################################################################
##################################################################################
# Power spectrum

# obsolete, see previous versions

##################################################################################
##################################################################################
# RSD power spectrum

import p3d_rsd
reload(p3d_rsd)
from p3d_rsd import *


pRsdHa = {}
for key in lfHa.keys():
   pRsdHa[key] = P3dRsdAuto(u, profLimLfHa[key], massFunc, nProc=3)

pRsdOiii = {}
for key in lfOiii.keys():
   pRsdOiii[key] = P3dRsdAuto(u, profLimLfOiii[key], massFunc, nProc=3)

key = 'Cochrane17'
p = pRsdHa[key]


#p.plotPMuDpdce(lfHa[key].Z[0])
#p.plotP(lfHa[key].Z[0])
#p.plotBEff()
#p.save(z=p.Prof.Lf.Z[0])
#p.load(z=p.Prof.Lf.Z[0])


##################################################################################
# Compute the mass build up of the power spectrum
'''
p.plotCumulMassContributionP(mu=0)
p.plotMassContributionP(mu=0)

p.plotCumulMassContributionP(mu=0.5)
p.plotMassContributionP(mu=0.5)
'''

##################################################################################
# Show the 2h/1h/shot noise terms
'''
p.plotPTermsZ(mu=0.)
p.plotPTermsZ(mu=0.5)
'''
'''
for key in lfHa.keys():
   pRsdHa[key].plotPTermsZ(mu=0.)
   pRsdHa[key].plotPTermsZ(mu=0.5)
'''


##################################################################################
# Fisher forecast for f

'''
# SPHEREx specs: fractional uncertainty on f
z = p.Z[0]
dz = 0.5 # 1.
R = 40.  # 150.
fSky = 100. * (np.pi/180.)**2 / (4.*np.pi)   # convert [sq deg] to [sr]
fwhmPsf = 6. * np.pi/(180.*3600.)   # convert [arcsec] to [rad]
print p.sFOverFFisher(z, R, fwhmPsf, fSky, dz)
'''

p.plotRequiredAreaToDetectF()
