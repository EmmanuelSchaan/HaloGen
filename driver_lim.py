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

import lim_specs
reload(lim_specs)
from lim_specs import *


##################################################################################

# Basic functions and parameters, for background and fluctuations
u = UnivPlanck15()

#u.plotSigma2V1d()
#u.plotSigma2DispFog()
#u.plotKMaxParaSpectroRes()
#u.plotKMaxPerpPsf()
#u.plotKFPerp()
#u.plotTradeOffNModes()



##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=True)
massFunc = MassFuncST(u, save=False)
#massFunc = MassFuncTinker(u, save=True)


##################################################################################

#sfr = SfrFonseca16(u, massFunc)
#sfr = SfrMoster13(u, massFunc)

sfr = SfrMoster13Speagle14(u, massFunc, scatter=False, nProc=3, save=False)

'''
sfr.plotSfr()
sfr.plotSfrd()
sfr.plotnHEff()
sfr.plotNHEffSparsitySummary()
#sfr.plotNHEffSparsity(exp='SPHEREx')
#sfr.plotNHEffSparsity(exp='COMAP')
#sfr.plotNHEffSparsity(exp='CONCERTO')
sfr.plotBEff()
sfr.plotdlnAlldlnm()
#sfr.plotdlnMeanIntensitydlnm()
#sfr.plotdlnbEff2dlnm()
#sfr.plotdlnP1hdlnm()
'''

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

lfCii = {}
lfCii['Popping16'] = LFCiiPopping16(u)

lfCO = {}
lfCO['Popping16'] = LFCOPopping16(u, 1)   # CO 1-0 transition

lfLya = {}
lfLya['Cassata11'] = LFLyaCassata11(u)

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

# check that I reproduce fig 9 in Popping+16
#lfCii['Popping16'].plotLf(xLim=(1.e6, 1.e9), yLim=(1.e-6, 0.1), unit='Lsun')
lfCii['Popping16'].plotLf(xLim=(1.e6 * 3.846e33, 1.e9 * 3.846e33), yLim=(1.e-6, 0.1))

#unit = lfCO['Popping16'].convertLumUnit('cgs') / lfCO['Popping16'].convertLumUnit('Jy*km/s*Mpc^2')
#lfCO['Popping16'].plotLf(xLim=(1.e5*unit, 1.e10*unit), yLim=(10.**(-4.5), 10.**(-0.5)))
lfCO['Popping16'].plotLf(xLim=(1.e36, 1.e40), yLim=(10.**(-4.5), 10.**(-0.5)))

lfLya['Cassata11'].plotLf(xLim=(1.e41, 1.e44), yLim=(1.e-5, 0.1))
'''



##################################################################################
# Plot: mean intensity
'''
unit = 'Jy/sr'
#lfHa['Sobral12'].plotMeanIntensity(lfs=[lfHa[key] for key in lfHa.keys()], unit=unit)
#lfOiii['Colbert13'].plotMeanIntensity(lfs=[lfOiii[key] for key in lfOiii.keys()], unit=unit)
lfCii['Popping16'].plotMeanIntensity(unit=unit)
#lfCO['Popping16'].plotMeanIntensity(unit=unit)
#lfLya['Cassata11'].plotMeanIntensity(unit=unit)
'''

##################################################################################
# Plot: nGalEff

'''
lfHa['Sobral12'].plotnGalEff(lfs=[lfHa[key] for key in lfHa.keys()])
lfOiii['Colbert13'].plotnGalEff(lfs=[lfOiii[key] for key in lfOiii.keys()])
lfCii['Popping16'].plotnGalEff()
lfCO['Popping16'].plotnGalEff()
lfLya['Cassata11'].plotnGalEff()
'''

##################################################################################
# Plot: NGalEff

'''
lfHa['Sobral12'].plotNGalEffSparsitySummary()
#lfHa['Sobral12'].plotNGalEffSparsity(lfs=[lfHa[key] for key in lfHa.keys()], exp='SPHEREx')
#lfOiii['Colbert13'].plotNGalEffSparsity(lfs=[lfOiii[key] for key in lfOiii.keys()], exp='SPHEREx')
#lfCii['Popping16'].plotNGalEffSparsity(exp='CONCERTO')
#lfCO['Popping16'].plotNGalEffSparsity(exp='COMAP')
#lfLya['Cassata11'].plotNGalEffSparsity(exp='SPHEREx')
'''

##################################################################################
# Plot: shot noise
'''
lfHa['Sobral12'].plotShotNoise(lfs=[lfHa[key] for key in lfHa.keys()])
lfOiii['Colbert13'].plotShotNoise(lfs=[lfOiii[key] for key in lfOiii.keys()])
lfCii['Popping16'].plotShotNoise()
lfCO['Popping16'].plotShotNoise()
lfLya['Cassata11'].plotShotNoise()
'''

##################################################################################
# Contributions from each luminosity
'''
lfHa['Sobral12'].plotdlnAlldlnL()
#lfHa['Sobral12'].plotdlnMeanIntensitydlnL()
#lfHa['Sobral12'].plotdlnPshotdlnL()
'''

##################################################################################
##################################################################################
# LIM from LF

profLimLfHa = {}
for key in lfHa.keys():
   profLimLfHa[key] = ProfLIMLF(u, sfr, lfHa[key], trunc=4., a=1.)

profLimLfOiii = {}
for key in lfOiii.keys():
   profLimLfOiii[key] = ProfLIMLF(u, sfr, lfOiii[key], trunc=4., a=1.)

profLimLfCii = {}
for key in lfCii.keys():
   profLimLfCii[key] = ProfLIMLF(u, sfr, lfCii[key], trunc=4., a=1.)

profLimLfCO = {}
for key in lfCO.keys():
   profLimLfCO[key] = ProfLIMLF(u, sfr, lfCO[key], trunc=4., a=0.6)

profLimLfLya = {}
for key in lfLya.keys():
   profLimLfLya[key] = ProfLIMLF(u, sfr, lfLya[key], trunc=4., a=1.)



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
# RSD auto power spectrum

import p3d_rsd
reload(p3d_rsd)
from p3d_rsd import *

pRsdHa = {}
for key in lfHa.keys():
   pRsdHa[key] = P3dRsdAuto(u, profLimLfHa[key], massFunc, nProc=3)

pRsdOiii = {}
for key in lfOiii.keys():
   pRsdOiii[key] = P3dRsdAuto(u, profLimLfOiii[key], massFunc, nProc=3)

pRsdCii = {}
for key in lfCii.keys():
   pRsdCii[key] = P3dRsdAuto(u, profLimLfCii[key], massFunc, nProc=3)

pRsdCO = {}
for key in lfCO.keys():
   pRsdCO[key] = P3dRsdAuto(u, profLimLfCO[key], massFunc, nProc=3)

pRsdLya = {}
for key in lfLya.keys():
   pRsdLya[key] = P3dRsdAuto(u, profLimLfLya[key], massFunc, nProc=3)



key = 'Cochrane17'
p = pRsdHa[key]


#p.plotP(lfHa[key].Z[0])
#p.plotBEff()
#p.save(z=p.Prof.Lf.Z[0])
#p.load(z=p.Prof.Lf.Z[0])


##################################################################################
# Scales probed by the various experiments
'''
p.plotPMuDpdce(lfHa[key].Z[0], exp='SPHEREx')
pRsdCO['Popping16'].plotPMuDpdce(2., exp='COMAP')
pRsdCii['Popping16'].plotPMuDpdce(6., exp='CONCERTO')

p.plotFourierModes()
'''


##################################################################################
# LIM vs galaxy surveys
'''
# in 3d RSD
pRsdHa['Cochrane17'].plotSigmaLumMatchedFilter(exp='SPHEREx')
pRsdOiii['Colbert13'].plotSigmaLumMatchedFilter(exp='SPHEREx')
pRsdCO['Popping16'].plotSigmaLumMatchedFilter(exp='COMAP')
pRsdCii['Popping16'].plotSigmaLumMatchedFilter(exp='CONCERTO')
'''
# in 2d
pRsdCii['Popping16'].plotSigmaLumMatchedFilter(exp='CCAT-P')


##################################################################################
# Compare references

'''
pRsdHa['Cochrane17'].compareP(ps=[pRsdHa[key] for key in pRsdHa.keys()])
pRsdOiii['Colbert13'].compareP(ps=[pRsdOiii[key] for key in pRsdOiii.keys()])
pRsdCii['Popping16'].compareP()
pRsdCO['Popping16'].compareP()
pRsdLya['Cassata11'].compareP()
'''


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


## SPHEREx specs: fractional uncertainty on f
#z = p.Z[0]
#dz = 0.5 # 1.
#R = 40.  # 150.
#fSky = 100. * (np.pi/180.)**2 / (4.*np.pi)   # convert [sq deg] to [sr]
#fwhmPsf = 6. * np.pi/(180.*3600.)   # convert [arcsec] to [rad]
#print p.sFOverFFisher(z, R, fwhmPsf, fSky, dz)


'''
# Halpha with SPHEREx
##p.plotRequiredAreaToDetectBeta()
#p.plotRequiredAreaToDetectA(2)
#p.plotRequiredAreaToDetectA(0)
#p.plotRequiredAreaToDetectAUnmarginalized(0, kMax=0.1, exp='SPHEREx')
#p.plotRequiredAreaToDetectAUnmarginalized(2, kMax=0.1, exp='SPHEREx')
p.plotRequiredAreaToDetectA(kMax=0.1, exp='SPHEREx', marg=False)
p.plotRequiredAreaToDetectA(kMax=0.1, exp='SPHEREx', marg=True)
'''
'''
# CO with COMAP
pRsdCO['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='COMAP', marg=False)
pRsdCO['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='COMAP', marg=True)
'''
'''
# [CII] with CONCERTO
pRsdCii['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='CONCERTO', marg=False)
pRsdCii['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='CONCERTO', marg=True)
'''


##################################################################################
##################################################################################
# RSD cross power spectrum

'''
key = 'Cochrane17'
profHa = ProfLIMLF(u, sfr, lfHa[key], trunc=4., a=1.1)
key = 'Mehta15'
profOiii = ProfLIMLF(u, sfr, lfOiii[key], trunc=4., a=0.8)


pHaOiii = P3dRsdCross(u, profHa, profOiii, massFunc, r=0.65, nProc=3)


pHaOiii.plotCorrCoeff(Z=[1., 2.])
'''

'''
# Example with high correlation coefficients: Ha - Oiii
# line correlation coefficient from Mehta+15
pHaOiii = P3dRsdCross(u, profLimLfHa['Cochrane17'], profLimLfOiii['Mehta15'], massFunc, r=0.65, nProc=3)

pHaOiii.plotCorrCoeff(Z=[1., 2.])
'''
'''
# Example with low correlation coefficient: Lya - CO
# line correlation coefficient from EGG
pLyaCO = P3dRsdCross(u, profLimLfLya['Cassata11'], profLimLfCO['Popping16'], massFunc, r=0.087, nProc=3)

pLyaCO.plotCorrCoeff(Z=[3., 4.])
'''


##################################################################################
# Experimental specs

'''
import lim_specs
reload(lim_specs)
from lim_specs import *

spherexSpecs = LimSpecs(u, exp='SPHEREx')
print spherexSpecs.whiteNoisePower(1.)

comapSpecs = LimSpecs(u, exp='COMAP')
print comapSpecs.whiteNoisePower(1.)

concertoSpecs = LimSpecs(u, exp='CONCERTO')
print concertoSpecs.whiteNoisePower(5.)
'''



