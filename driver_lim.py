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
lfHa['Egg'] = LFEGG(u, lineName='halpha')

lfOiii = {}
lfOiii['Colbert13'] = LFOiiiColbert13(u)
lfOiii['Mehta15'] = LFOiiiMehta15(u)
lfOiii['5007Egg'] = LFEGG(u, lineName='o3_5007')
lfOiii['4959Egg'] = LFEGG(u, lineName='o3_4959')


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
'''
for key in lfHa.keys():
   try: lfHa[key].plotLf()
   except: pass

for key in lfOiii.keys():
   try: lfOiii[key].plotLf()
   except: pass
'''


##################################################################################
# Plot: mean intensity

'''
lfHa['Sobral12'].plotMeanIntensity(lfs=[lfHa[key] for key in lfHa.keys()])

lfOiii['Colbert13'].plotMeanIntensity(lfs=[lfOiii[key] for key in lfOiii.keys()])
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

p3dHa = {}
for key in lfHa.keys():
   p3dHa[key] = P3dAuto(u, iHaloModel, profLimLfHa[key], doT=False, save=True)

p3dOiii = {}
for key in lfOiii.keys():
   p3dOiii[key] = P3dAuto(u, iHaloModel, profLimLfOiii[key], doT=False, save=True)


##p3d_limhaegg.plotP(z=1.)
##p3d_limhaegg.plotBEff()

#p3d_limoiii5007egg.plotP(z=1.)
#p3d_limoiii4959egg.plotP(z=1.)

##profLimLfHaEgg.plotP3dGong17(p3d_limegg.fPtotinterp, lineName='halpha')




##################################################################################
# Plot: power spectrum


for key in lfHa.keys():

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for z in lfHa[key].Z:    
      f = lambda k: p3dHa[key].fPtotinterp(k, z)
      pTot = np.array(map(f, p3dHa[key].K)) * lfHa[key].convertPowerSpectrumUnit('Jy/sr')
      f = lambda k: p3dHa[key].fP1hinterp(k, z)
      p1h = np.array(map(f, p3dHa[key].K)) * lfHa[key].convertPowerSpectrumUnit('Jy/sr')
      #f = lambda k: p3dHa[key].fP2hinterp(k, z)
      #p2h = np.array(map(f, p3dHa[key].K)) * lfHa[key].convertPowerSpectrumUnit('Jy/sr')
      #print p2h
      f = lambda k: p3dHa[key].fPnoise(k, z)
      pShot = np.array(map(f, p3dHa[key].K)) * lfHa[key].convertPowerSpectrumUnit('Jy/sr')
      #
      plot=ax.loglog(p3dHa[key].K, pTot, label=r'$z=$'+str(round(z, 2)))
      ax.loglog(p3dHa[key].K, p1h, ls='--', c=plot[0].get_color(), lw=1)
      #ax.loglog(p3dHa[key].K, p2h, ls='-', c=plot[0].get_color())
      ax.loglog(p3dHa[key].K, pShot, ls=':', c=plot[0].get_color(), lw=1)
   #
   ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$k$ [$h$/Mpc]')
   ax.set_ylabel(r'$P(k,z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
   ax.set_title(lfHa[key].name.replace("_", ""))
   #
   path = './figures/pn_3d/p3d_'+lfHa[key].name+'.pdf'
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()



for key in lfOiii.keys():

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   for z in lfOiii[key].Z:    
      f = lambda k: p3dOiii[key].fPtotinterp(k, z)
      pTot = np.array(map(f, p3dOiii[key].K)) * lfOiii[key].convertPowerSpectrumUnit('Jy/sr')
      f = lambda k: p3dOiii[key].fP1hinterp(k, z)
      p1h = np.array(map(f, p3dOiii[key].K)) * lfOiii[key].convertPowerSpectrumUnit('Jy/sr')
      #f = lambda k: p3dOiii[key].fP2hinterp(k, z)
      #p2h = np.array(map(f, p3dOiii[key].K)) * lfOiii[key].convertPowerSpectrumUnit('Jy/sr')
      #print p2h
      f = lambda k: p3dOiii[key].fPnoise(k, z)
      pShot = np.array(map(f, p3dOiii[key].K)) * lfOiii[key].convertPowerSpectrumUnit('Jy/sr')
      #
      plot=ax.loglog(p3dOiii[key].K, pTot, label=r'$z=$'+str(round(z, 2)))
      ax.loglog(p3dOiii[key].K, p1h, ls='--', c=plot[0].get_color(), lw=1)
      #ax.loglog(p3dOiii[key].K, p2h, ls='-', c=plot[0].get_color())
      ax.loglog(p3dOiii[key].K, pShot, ls=':', c=plot[0].get_color(), lw=1)
   #
   ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$k$ [$h$/Mpc]')
   ax.set_ylabel(r'$P(k,z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
   ax.set_title(lfOiii[key].name.replace("_", ""))
   #
   path = './figures/pn_3d/p3d_'+lfOiii[key].name+'.pdf'
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
   #plt.show()

