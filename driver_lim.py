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

#import i_halo_model
#reload(i_halo_model)
#from i_halo_model import *

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
'''
u.plotKMaxPerpPsf()
u.plotKFPerp()
'''
#u.plotTradeOffNModes()



##################################################################################

# Several mass functions implemented: Press-Schechter, Sheth-Tormen, Tinker
#massFunc = MassFuncPS(u, save=True)
massFunc = MassFuncST(u, save=False)
#massFunc = MassFuncTinker(u, save=True)


##################################################################################
# Experimental specs


#import lim_specs
#reload(lim_specs)
#from lim_specs import *

spherexSpecs = LimSpecs(u, exp='SPHEREx')
'''
print floatExpForm(spherexSpecs.whiteNoisePower(0.8))
print floatExpForm(spherexSpecs.whiteNoisePower(1.))
print floatExpForm(spherexSpecs.whiteNoisePower(2.))
'''
comapSpecs = LimSpecs(u, exp='COMAP')
'''
print floatExpForm(comapSpecs.whiteNoisePower(1.))
print floatExpForm(comapSpecs.whiteNoisePower(2.))
print floatExpForm(comapSpecs.whiteNoisePower(4.))
print floatExpForm(comapSpecs.whiteNoisePower(6.))
'''
concertoSpecs = LimSpecs(u, exp='CONCERTO')
'''
print floatExpForm(concertoSpecs.whiteNoisePower(1.))
print floatExpForm(concertoSpecs.whiteNoisePower(2.))
print floatExpForm(concertoSpecs.whiteNoisePower(3.))
print floatExpForm(concertoSpecs.whiteNoisePower(4.))
print floatExpForm(concertoSpecs.whiteNoisePower(6.))
'''
hetdexSpecs = LimSpecs(u, exp='HETDEX')
'''
print floatExpForm(hetdexSpecs.whiteNoisePower(2.475))
print floatExpForm(hetdexSpecs.whiteNoisePower(3.775))
print floatExpForm(hetdexSpecs.whiteNoisePower(5.575))
'''
cdimSpecs = LimSpecs(u, exp='CDIM')



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
sfr.plotdlnAlldlnm(alpha=1.)
sfr.plotdlnAlldlnm(alpha=0.6)
sfr.plotdlnAlldlnm(alpha=1.1)
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

#lf = lfHa['Cochrane17']
#lf.plotnGal()
#lf.plotMeanIntensity()
#lf.plotLf()
#lf.plotnGalEff()
#lf.plotShotNoise()
#lf.plotS2()
#lf.plotSchechterProperties()

#print lfOiiiMehta15.computeCorrCoeffHaOIII()
#lfOiiiMehta15.plotBivariateLf()

##################################################################################
# General properties of Schechter functions
'''
lf = LFHaCochrane17(u)
lf.plotSchechterProperties()
'''

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
lfHa['Cochrane17'].plotLf(lfs=[lfHa[key] for key in lfHa.keys()], yLim=(1.e-10, 1.))
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

unit = 'Jy/sr'
'''
lfHa['Sobral12'].plotMeanIntensity(lfs=[lfHa[key] for key in lfHa.keys()], unit=unit)
lfOiii['Colbert13'].plotMeanIntensity(lfs=[lfOiii[key] for key in lfOiii.keys()], unit=unit)
lfCii['Popping16'].plotMeanIntensity(unit=unit)
lfCO['Popping16'].plotMeanIntensity(unit=unit)
lfLya['Cassata11'].plotMeanIntensity(unit=unit)
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
lfHa['Sobral12'].plotNGalEffSparsity(lfs=[lfHa[key] for key in lfHa.keys()], exp='SPHEREx', sfr=sfr)
lfOiii['Colbert13'].plotNGalEffSparsity(lfs=[lfOiii[key] for key in lfOiii.keys()], exp='SPHEREx', sfr=sfr)
lfCii['Popping16'].plotNGalEffSparsity(exp='CONCERTO', sfr=sfr)
lfCO['Popping16'].plotNGalEffSparsity(exp='COMAP', sfr=sfr)
lfLya['Cassata11'].plotNGalEffSparsity(exp='SPHEREx', sfr=sfr)
'''
'''
lfHa['Sobral12'].plotNGalEffSparsitySummary()
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


## EGG
#profLimLfHaEgg = ProfLIMLF(u, sfr, lfHaEgg, trunc=4.)
#profLimLfOiii5007Egg = ProfLIMLF(u, sfr, lfOiii5007Egg, trunc=4.)
#profLimLfOiii4959Egg = ProfLIMLF(u, sfr, lfOiii4959Egg, trunc=4.)


##################################################################################
# Plot the halo mass-luminosity relation

'''
profLimLfHa['Cochrane17'].plotLuminosityMassRelation()
profLimLfOiii['Colbert13'].plotLuminosityMassRelation()
profLimLfCii['Popping16'].plotLuminosityMassRelation()
profLimLfCO['Popping16'].plotLuminosityMassRelation()
profLimLfLya['Cassata11'].plotLuminosityMassRelation()
'''

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


#p = pRsdLya['Cassata11']
#z = p.Prof.Lf.Z[-1] 
#p.load(z=z)
##p.plotP(z=z, unit='cgs')
#p.plotP(z=z, unit='Lsun/(Mpc/h)^2/sr/Hz')
#p.plotP(z=z, unit='Jy/sr')




##################################################################################
# Scales probed by the various experiments
'''
p.plotPMuDpdce(lfHa[key].Z[0], exp='SPHEREx')
pRsdCO['Popping16'].plotPMuDpdce(2., exp='COMAP')
pRsdCii['Popping16'].plotPMuDpdce(6., exp='CONCERTO')
'''
'''
p.plotFourierModes()
'''

##################################################################################
# LIM vs galaxy surveys
'''
# in 3d RSD
pRsdHa['Cochrane17'].plotSigmaLumMatchedFilter(specs=cdimSpecs)
pRsdLya['Cassata11'].plotSigmaLumMatchedFilter(specs=hetdexSpecs)

pRsdHa['Cochrane17'].plotSigmaLumMatchedFilter(specs=spherexSpecs)
pRsdOiii['Colbert13'].plotSigmaLumMatchedFilter(specs=spherexSpecs)

pRsdCO['Popping16'].plotSigmaLumMatchedFilter(specs=comapSpecs)

pRsdCii['Popping16'].plotSigmaLumMatchedFilter(specs=concertoSpecs)
'''
'''
# summary plots
pList = [pRsdHa['Cochrane17'], pRsdHa['Cochrane17'], pRsdOiii['Colbert13'], pRsdLya['Cassata11'],  pRsdCO['Popping16'], pRsdCii['Popping16']]
specsList = [spherexSpecs, cdimSpecs, spherexSpecs, hetdexSpecs, comapSpecs, concertoSpecs]
pRsdHa['Cochrane17'].plotLimVsGalDet(pList, specsList)
'''

#pList = [pRsdLya['Cassata11']]
#specsList = [hetdexSpecs]
#pRsdLya['Cassata11'].plotLimVsGalDet(pList, specsList)


'''
# in 2d
pRsdCii['Popping16'].plotSigmaLumMatchedFilter2d(exp='CCAT-P')
'''

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

'''
# Halpha with SPHEREx
##p.plotRequiredAreaToDetectBeta()
#p.plotRequiredAreaToDetectA(2)
#p.plotRequiredAreaToDetectA(0)
#p.plotRequiredAreaToDetectAUnmarginalized(0, kMax=0.1, exp='SPHEREx')
#p.plotRequiredAreaToDetectAUnmarginalized(2, kMax=0.1, exp='SPHEREx')
p.plotRequiredAreaToDetectA(kMax=0.1, exp='SPHEREx', marg=False)
p.plotRequiredAreaToDetectA(kMax=0.1, exp='SPHEREx', marg=True)


# CO with COMAP
pRsdCO['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='COMAP', marg=False)
pRsdCO['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='COMAP', marg=True)


# [CII] with CONCERTO
pRsdCii['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='CONCERTO', marg=False)
pRsdCii['Popping16'].plotRequiredAreaToDetectA(kMax=0.1, exp='CONCERTO', marg=True)
'''


##################################################################################
##################################################################################
# RSD cross power spectrum

'''
# Example with high correlation coefficients: Ha - Oiii
# line correlation coefficient from Mehta+15
pHaOiii = P3dRsdCross(u, profLimLfHa['Cochrane17'], profLimLfOiii['Mehta15'], massFunc, r=0.65, nProc=3)
pHaOiii.plotCorrCoeff(Z=[1., 2.])

# Example with low correlation coefficient: Lya - CO
# line correlation coefficient from EGG
pLyaCO = P3dRsdCross(u, profLimLfLya['Cassata11'], profLimLfCO['Popping16'], massFunc, r=0.088, nProc=3)
pLyaCO.plotCorrCoeff(Z=[3., 4.])
'''


##################################################################################
##################################################################################
# Bispectrum

import b3d_rsd
reload(b3d_rsd)
from b3d_rsd import *
'''
# Matter density profile
profMatter = ProfNFW(u)

# the bispectrum between matter - Ha - Lya
#b3d = B3dRsdCross(u, profLimLfHa['Cochrane17'], profLimLfLya['Cassata11'], profMatter, massFunc, nProc=3)
b3d = B3dRsdCross(u, profLimLfHa['Cochrane17'], profLimLfHa['Cochrane17'], profMatter, massFunc, nProc=3)

print "b1h"
print b3d.b1h(1.e-5, 1.e-5, 1.e-5, 2.)
'''

'''
# Compute the amplitude of the matter - LIM - LIM bispectrum
# 1-halo term only
print "amplitude"
print profLimLfHa['Cochrane17'].b1hMatLimLimAmplitude(2.)
'''

