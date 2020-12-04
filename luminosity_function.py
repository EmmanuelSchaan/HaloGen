from headers import *


class LF(object):
   '''Generic luminosity function class,
   for a given line / band.
   Default luminosity unit: [Lsun]
   Default intensity unit: [Lsun/(Mpc/h)^2/sr/Hz]
   '''

   def __init__(self, U):
      self.U = U

      # required attributes
      #self.name
      #self.nameLatex   # which line and which paper reference
      #self.refLatex # paper reference
      #self.lineName # which line this is
      #self.lambdaMicrons = 656.28e-3   # [mu]
      #self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]
      #self.lMin
      #self.lMax
      #self.zMin
      #self.zMax
      #self.Z
      #self.phi(l, z)#  [number / (Mpc/h)^3 / Lsun]


      # Interpolate for speed
      x = np.array(map(self.nGal, self.Z))
      self.nGalInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   # [(Mpc/h)^-3]
      #
      x = np.array(map(self.lumDensity, self.Z))
      self.lumDensityInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   #[Lsun (Mpc/h)^-3]
      #
      x = np.array(map(self.sqLumDensity, self.Z))
      self.sqLumDensityInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   #[Lsun^2 (Mpc/h)^-3]
      #
      x = np.array(map(self.meanIntensity, self.Z))
      self.meanIntensityInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   #[Lsun / (Mpc/h)^2 / sr / Hz]
      #
      x = np.array(map(self.meanGalLum, self.Z))
      self.meanGalLumInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   #[Lsun]
      #
      x = np.array(map(self.pShot, self.Z))
      self.pShotInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   # [(intensity unit)^2 (Mpc/h)^3] = [Lsun^2 / (Mpc/h) / sr^2 / Hz^2]
      #
      x = np.array(map(self.nGalEff, self.Z))
      self.nGalEffInterp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   # [(Mpc/h)^-3]
      #
      x = np.array(map(self.s2, self.Z))
      self.s2Interp = interp1d(self.Z, x, kind='linear', bounds_error=False, fill_value=0.)   # [(Mpc/h)^-3]



   def __str__(self):
      return self.name


   ##################################################################################
   # Moments of the luminosity function

   def lumMoment(self, z, n, lMin=0., lMax=np.inf):
      ''' [Lsun^n (Mpc/h)^-3]
      This integral of the LF formally diverges at low luminosities,
      so it is highly cutoff dependent. Use with caution.
      '''
      def integrand(lnl):
         l = np.exp(lnl)
         result = self.phi(z, l)
         result *= l**(n+1)
         return result
      # integration bounds
      lMin = max(lMin, self.lMin)
      lMax = min(lMax, self.lMax)
      result = integrate.quad(integrand, np.log(lMin), np.log(lMax), epsabs=0., epsrel=1.e-3)[0]
      return result


   def nGal(self, z, lMin=0., lMax=np.inf):
      '''Mean number density of galaxies [(Mpc/h)^{-3}]
      This quantity may be luminosity-cutoff dependent, and may not converge,
      depending on the LF
      '''
      return self.lumMoment(z, 0, lMin=lMin, lMax=lMax)


   def lumDensity(self, z, lMin=0., lMax=np.inf):
      '''Mean luminosity density [Lsun (Mpc/h)^-3]
      '''
      return self.lumMoment(z, 1, lMin=lMin, lMax=lMax)


   def sqLumDensity(self, z, lMin=0., lMax=np.inf):
      '''Mean squared luminosity density [Lsun^2 (Mpc/h)^-3]
      '''
      return self.lumMoment(z, 2, lMin=lMin, lMax=lMax)


   ##################################################################################
   # Quantities derived from the moments of the luminosity function

   def meanIntensity(self, z, lMin=0., lMax=np.inf):
      '''[Lsun / (Mpc/h)^2 / sr / Hz]
      '''
      result = self.lumDensity(z, lMin=lMin, lMax=lMax)  # [Lsun / (Mpc/h)^3]
      result *= 3.e5 / self.U.hubble(z)   # *[Mpc/h]
      result /= 4. * np.pi * self.nuHz # *[/sr/Hz]
      return result


   def meanGalLum(self, z, lMin=0., lMax=np.inf):
      '''Mean galaxy luminosity [Lsun]
      This quantity may be luminosity-cutoff dependent,
      depending on the LF
      '''
      result = self.lumDensity(z, lMin=lMin, lMax=lMax)  # [Lsun / (Mpc/h)^3]
      result /= self.nGal(z, lMin=lMin, lMax=lMax) # / [(Mpc/h)^{-3}] = [Lsun]
      return result


   def pShot(self, z, lMin=0., lMax=np.inf):
      '''[(intensity unit)^2 (Mpc/h)^3] = [Lsun^2 / (Mpc/h) / sr^2 / Hz^2]
      '''
      result = self.sqLumDensity(z, lMin=lMin, lMax=lMax)  # [Lsun^2 / (Mpc/h)^3]
      result *= (3.e5 / self.U.hubble(z))**2   # *[(Mpc/h)^2]
      result /= (4. * np.pi * self.nuHz)**2 # *[1/sr^2/Hz^2]
      return result


   def nGalEff(self, z, lMin=0., lMax=np.inf):
      '''Effective galaxy number density for the shot noise [(Mpc/h)^{-3}]
      '''
      result = self.meanIntensity(z, lMin=lMin, lMax=lMax)**2
      result /= self.pShot(z, lMin=lMin, lMax=lMax)
      return result

   
   def s2(self, z, lMin=0., lMax=np.inf):
      '''Fractional luminosity variance s2 = var(L_Ha) / mean(L_Ha)^2 [dimless]
      This quantity may be luminosity-cutoff dependent,
      depending on the LF
      '''
      result = self.nGal(z, lMin=lMin, lMax=lMax) 
      result /= self.nGalEff(z, lMin=lMin, lMax=lMax) 
      result -= 1.
      return result


   def dlnLumMomentdlnL(self, z, l, n, lMin=0., lMax=np.inf):
      '''l is galaxy luminosity in Lsun
      Result is [dimless]
      '''
      result = self.phi(z, l) * l**(n+1)
      result /= self.lumMoment(z, n, lMin=lMin, lMax=lMax)
      result *= (l>=lMin) * (l<=lMax)
      return result


   ##################################################################################

   def plotdlnMeanIntensitydlnL(self):
      L_cgs = np.logspace(np.log10(1.e35), np.log10(1.e44), 101, 10.) # [erg/s]
      L_Lsun = L_cgs / self.convertLumUnit('cgs')

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(self.Z)):
         z = self.Z[iZ]
         f = lambda l: self.dlnLumMomentdlnL(z, l, 1)
         y = np.array(map(f, L_Lsun))
         ax.plot(L_cgs, y, c=plt.cm.cool(iZ/(len(self.Z)-1.)), label=r'$z=$'+str(z))
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$L$ [erg/s]')
      ax.set_ylabel(r'$d \text{ln} I / d\text{ln} L$')
      #
      path = './figures/lf/dlnmeanintensitydlnl_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()

   
   def plotdlnPshotdlnL(self):
      L_cgs = np.logspace(np.log10(1.e39), np.log10(1.e44), 101, 10.) # [erg/s]
      L_Lsun = L_cgs / self.convertLumUnit('cgs')

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(self.Z)):
         z = self.Z[iZ]
         f = lambda l: self.dlnLumMomentdlnL(z, l, 2)
         y = np.array(map(f, L_Lsun))
         ax.plot(L_cgs, y, c=plt.cm.cool(iZ/(len(self.Z)-1.)), label=r'$z=$'+str(z))
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$L$ [erg/s]')
      ax.set_ylabel(r'$d \text{ln} P_\text{shot} / d\text{ln} L$')
      #
      path = './figures/lf/dlnpshotdlnl_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   ##################################################################################

   def convertLumUnit(self, unit):
      '''Choice of luminosity unit:
      'Lsun'
      'cgs' for [erg/s]
      'Jy*m^2*Hz'
      Returns the factor by which to multiply to convert to the requested unit
      '''
      if unit=='Lsun':
         return 1.
      elif unit=='cgs':
         return 3.839e33  # [Lsun] to [erg/s]
      elif unit=='Jy*m^2*Hz':
         result = 3.827e26   # [Lsun] to [W]
         result /= 1.e-26  # [W] to [Jy*m^2*Hz]
         return result

   def convertLfUnit(self, unit):
      return 1. / self.convertLumUnit(unit)

   def convertIntensityUnit(self, unit):
      '''Choice of intensity units:
      'Lsun/(Mpc/h)^2/sr/Hz'
      'cgs' for [erg/s/cm^2/sr/Hz]
      'Jy/sr'
      '''
      if unit=='Lsun/(Mpc/h)^2/sr/Hz':
         return 1.
      if unit=='Jy/sr':
         result = 3.827e26   # [Lsun] to [W]
         result /= (3.086e22 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [m^{-2}]
         result /= 1.e-26  # [W/m^2/Hz/sr] to [Jy/sr]
         return result
      elif unit=='cgs':
         result = 3.839e33  # [Lsun] to [erg/s]
         result /= (3.086e24 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [cm^{-2}]
         return result

   def convertPowerSpectrumUnit(self, unit):
      return self.convertIntensityUnit(unit)**2


   ##################################################################################

   def plotLf(self, lfs=None):
      
      if lfs is None:
         lfs = [self]

      L_cgs = np.logspace(np.log10(1.e40), np.log10(1.e44), 501, 10.) # [erg/s]
      L_Lsun = L_cgs / self.convertLumUnit('cgs')
      
      # find min and max redshifts in all the LFs to plot
      zMin = np.min([lf.zMin for lf in lfs])
      zMax = np.max([lf.zMax for lf in lfs])

      # "Observed" LF,
      # ie affected by dust attenuation
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      lineStyles = np.array(['-', '--', ':'])
      legendItems = []
      for iLf in range(len(lfs)):
         lf = lfs[iLf]
         ls = lineStyles[iLf]
         for iZ in range(lf.nZ):
            z = lf.Z[iZ]
            c = plt.cm.cool((z-zMin)/(zMax-zMin))
            #c = plt.cm.rainbow((z-zMin)/(zMax-zMin))
            #
            # Schechter fit
            if hasattr(lf, 'phi'):
               y = lf.phi(z, L_Lsun) * L_Lsun * np.log(10.) * self.U.bg.h**3
               #ax.plot(L_cgs, y, label=r'$z=$'+str(round(z,2))+' '+lf.nameLatex)
               #line, =ax.plot(L_cgs, y, c=c, ls=ls,label=r'$z=$'+str(round(z,2))+' '+lf.refLatex)
               line, =ax.plot(L_cgs, y, c=c, ls=ls)
               legendItems.append((z, line, r'$z=$'+str(round(z,2))+' '+lf.refLatex))
      #
      legendItems.sort()
      #ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
      ax.legend([x[1] for x in legendItems], [x[2] for x in legendItems], loc=3, fontsize='x-small', labelspacing=0.1)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((10.**(40.), 10.**(44)))
      ax.set_ylim((1.e-8, 1.))
      ax.set_xlabel(r'$L$ [erg/s]')
      ax.set_ylabel(r'$\text{log}_{10} \left( \Phi \times \text{Mpc}^3 \right)$')
      ax.set_title(self.lineNameLatex+' Luminosity function')
      #
      path = './figures/lf/lf_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()
      #plt.close()
      
      


   def plotnGal(self, lfs=None):
      if lfs is None:
         lfs = [self]
      
      print("nGal is highly cutoff dependent (formally divergent at low luminosity)")
      print("so its value is meaningless")
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for lf in lfs:
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         nGal = np.array(map(lf.nGal, Z))
         #ax.semilogy(Z, nGal, label=str(lf).replace('_', ' '))
         ax.semilogy(Z, nGal, label=lf.nameLatex)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{n}_\text{gal}$ [(Mpc/h)$^{-3}$]')

      plt.show()


   def plotMeanIntensity(self, lfs=None):
      if lfs is None:
         lfs = [self]

      # reproduce Fig 3 in Gong Cooray Silva + 17
      # in Jy/sr
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # compare to Gong+17
      # in Jy/sr
      if self.lineName=='halpha':
         # center values
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_center.csv"
         data = np.genfromtxt(path, delimiter=', ')
         plt.plot(data[:,0], data[:,1], 'b', label=r'HB06')
         # error band
         #path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_low.csv"
         #low = np.genfromtxt(path, delimiter=', ')
         #path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_high.csv"
         #high = np.genfromtxt(path, delimiter=', ')
         #plt.fill(np.append(low[:,0], high[::-1,0]), np.append(low[:,1], high[::-1,1]), facecolor='b', alpha=0.5)
      #
      for lf in lfs:
         f = lambda z: lf.meanIntensity(z) * self.convertIntensityUnit('Jy/sr') 
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         meanIntensity = np.array(map(f, Z))
         #ax.plot(Z, meanIntensity, label=str(lf).replace('_', ' '))
         ax.plot(Z, meanIntensity, label=lf.refLatex)
      #
      ax.set_title(self.lineNameLatex+' Mean intensity')
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{I}(z)$ [Jy/sr]')
      ax.set_xlim((0., 5.))
      ax.set_ylim((0., 40.))
      #
      path = './figures/lf/mean_intensity_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()
      #plt.close()

      '''
      # reproduce Fig 3 in Fonseca+16
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      f = lambda z: self.meanIntensity(z, unit='cgs')
      meanIntensity = np.array(map(f, Z))
      ax.semilogy(Z, self.nuHz * meanIntensity)
      #
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\nu \bar{I}(z)$ [erg/s/cm$^2$/sr]')

      plt.show()
      '''


   def plotnGalEff(self, lfs=None):
      if lfs is None:
         lfs = [self]
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for lf in lfs:
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         nGalEff = np.array(map(lf.nGalEff, Z))
         #ax.semilogy(Z, nGalEff, label=str(lf).replace('_', ' '))
         ax.semilogy(Z, nGalEff, label=lf.refLatex)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{n}_\text{gal eff}$ [(Mpc/h)$^{-3}$]')
      ax.set_title(self.lineNameLatex)
      #
      path = './figures/lf/ngaleff_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for lf in lfs:
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         # SPHEREx voxel size
         # the spectral resolution power is R=40 for the lower z, and 150 for the high z
         R = 40.
         # hence the redshift size of the voxel
         dz = (1. + Z) / R
         # and the comoving depth of the voxel
         dChi = dz * 3.e5 / self.U.hubble(Z)   # [Mpc/h]
         # angular pixel size: 6.2 arcsec
         thetaPix = 6.2 * np.pi/(180.*3600.)
         # hence the voxel comoving volume
         vVoxSpherex = (self.U.bg.comoving_distance(Z) * thetaPix)**2 * dChi  # [(Mpc/h)^3]
         #print "vVoxSpherex=", vVoxSpherex, "(Mpc/h)^3"

         # effective number of galaxies
         n = np.array(map(lf.nGalEff, Z))
         n *= vVoxSpherex
         ax.plot(Z, n, label=lf.refLatex)
      #
      ax.legend(loc='center right', fontsize='x-small', labelspacing=0.1)
      #ax.set_xlim((np.min(self.Z), np.max(self.Z)))
      #ax.set_xlim((0., 7.))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{N}^\text{gal eff}$ per SPHEREx voxel')
      ax.set_title(self.lineNameLatex)
      #
      path = './figures/lf/galaxy_sparsity_spherex_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()





   def plotShotNoise(self, lfs=None):
      if lfs is None:
         lfs = [self]

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for lf in lfs:
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         #f = lambda z: 1. / lf.nGalEff(z) # shot noise power of dI/I
         f = lambda z: lf.pShot(z)
         pShot = np.array(map(f, Z)) * self.convertPowerSpectrumUnit('Jy/sr')
         plt.plot(Z, pShot, label=lf.refLatex)
         
      #
      ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$P_\text{shot}$ [(Jy/sr)$^2$(Mpc/h)$^3$]')
      ax.set_title(self.lineNameLatex+' Shot noise')
      #
      path = './figures/lf/shot_noise_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      plt.show()


   def plotS2(self, lfs=None):
      if lfs is None:
         lfs = [self]
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for lf in lfs:
         if hasattr(lf, 'Z'):
            Z = lf.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         s2 = np.array(map(lf.s2, Z))
         #ax.semilogy(Z, s2, label=str(lf))
         ax.semilogy(Z, s2, label=lf.nameLatex)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$s^2 \equiv \text{var}(L) / \langle L \rangle^2$ [dimless]')

      plt.show()




##################################################################################
##################################################################################


class LFHaSobral12(LF):
   '''Halpha luminosity function from Sobral+12 arXiv:1202.3436v2
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'sobral12halpha'
      self.nameLatex = r'S13 H$\alpha$'
      self.refLatex = 'S13'
      self.lineName = 'halpha'
      self.lineNameLatex = r'H$\alpha$'
      self.lambdaMicrons = 656.28e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 5
      path = './input/LIM_literature/Sobral+12/table5/Sobral+12_table5.txt'
      data = np.genfromtxt(path)
      z = data[:,0]
      #
      self.Z = z
      self.nZ = len(z)
      self.zMin = np.min(z)
      self.zMax = np.max(z)
      #
      LStar = 10.**data[:,1] / 3.839e33  # convert from [erg/s] to [Lsun]
      #LStarHigh = 10.**(data[:,1]+data[:,2]) / 3.839e33  # convert from [erg/s] to [Lsun]
      #LStarLow = 10.**(data[:,1]+data[:,3]) / 3.839e33  # convert from [erg/s] to [Lsun]
      #
      PhiStar = 10.**data[:,4] / U.bg.h**3  # [(Mpc/h)^-3]
      #PhiStarHigh = 10.**(data[:,4]+data[:,5]) / U.bg.h**3  # [(Mpc/h)^-3]
      #PhiStarLow = 10.**(data[:,4]+data[:,6]) / U.bg.h**3  # [(Mpc/h)^-3]
      #
      # for alpha, they recommend using -1.6, rather than the best fit...
      #Alpha = -1.6 * np.ones_like(z) # [dimless]
      Alpha = data[:,7] # [dimless]
      #AlphaHigh = data[:,7]+data[:,8] # [dimless]
      #AlphaLow = data[:,7]+data[:,9] # [dimless]
      

      # interpolate the Schechter fits (Eq 2)
      self.lStar = interp1d(z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #self.lStarHigh = interp1d(z, LStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.lStarLow = interp1d(z, LStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.phiStar = interp1d(z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
      #self.phiStarHigh = interp1d(z, PhiStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.phiStarLow = interp1d(z, PhiStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.alpha = interp1d(z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #self.alphaHigh = interp1d(z, AlphaHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.alphaLow = interp1d(z, AlphaLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      # Intrinsic LF, with dust attenuation removed
      self.phiInt = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]
      # Observed LF, where the luminosities experience dust extinction
      self.phiObs = lambda z,l: 100.**(1./5.) * self.phiInt(z, l * 100.**(1./5.))   # [(Mpc/h)^-3 / Lsun]

      self.phi = self.phiObs 


      super(LFHaSobral12, self).__init__(U)

      self.loadMeasuredLF()



   def loadMeasuredLF(self):
      # Measured luminosity functions
      # table 4
      self.Z = np.array([0.4, 0.84, 1.47, 2.23])
      zString = ['0.40', '0.84', '1.47', '2.23']
      self.nZ = len(self.Z)
      pathIn = './input/LIM_literature/Sobral+12/table4/'
      # intrinsic luminosities [erg/s], ie corrected to undo dust extinction
      self.LInt = {}
      self.dLInt = {}
      self.dLog10LInt = {}

      # "observed" luminosity [erg/s], ie affected by dust extinction
      self.LObs = {} 
      self.dLObs = {}
      self.dLog10LObs = {}

      # luminosity functions and uncertainties [(Mpc/h)^-3 (erg/s)^-1]
      self.PhiInt = {}
      self.PhiIntHigh = {}
      self.PhiIntLow = {}
      #
      self.PhiObs = {}
      self.PhiObsHigh = {}
      self.PhiObsLow = {}

      # interpolated LFs [(Mpc/h)^-3 (erg/s)^-1]
      self.fPhiObsMeas = {}
      self.fPhiObsHighMeas = {}
      self.fPhiObsLowMeas = {}

      for iZ in range(self.nZ):
         z = self.Z[iZ]
         path = pathIn + 'Sobral+12_table5_z'+zString[iZ]+'.txt'
         data = np.genfromtxt(path)
         # intrinsic luminosities
         self.LInt[iZ] = 10.**data[:,0]   # [erg/s]
         self.dLInt[iZ] = 10.**(data[:,0] + data[:,1]) - 10.**(data[:,0] - data[:,1])   # [erg/s]
         self.dLog10LInt[iZ] = 2. * data[:,1]  # log10(L/(erg/s))
         # observed luminosities
         self.LObs[iZ] = self.LInt[iZ]  / (100.**(1./5.))  # adding 1 magnitude of dust extinction [erg/s]
         self.dLObs[iZ] = self.dLInt[iZ] / (100.**(1./5.))
         self.dLog10LObs[iZ] = self.dLog10LInt[iZ]  # unchanged by multiplicative dust extinction
         # intrinsic LF 
         self.PhiInt[iZ] = 10.**data[:,5] / self.U.bg.h**3 / self.LInt[iZ] / np.log(10.)
         self.PhiIntHigh[iZ] = self.PhiInt[iZ] * 10.**data[:,6]
         self.PhiIntLow[iZ] = self.PhiInt[iZ] / 10.**data[:,6]

         # observed LF
         self.PhiObs[iZ] = 10.**data[:,5] / self.U.bg.h**3 / self.LObs[iZ] / np.log(10.)
         self.PhiObsHigh[iZ] = self.PhiObs[iZ] * 10.**data[:,6]
         self.PhiObsLow[iZ] = self.PhiObs[iZ] / 10.**data[:,6]

         # interpolate
         self.fPhiObsMeas[iZ] = interp1d(self.LObs[iZ], self.PhiObs[iZ], kind='linear', bounds_error=False, fill_value=0.)
         self.fPhiObsHighMeas[iZ] = interp1d(self.LObs[iZ], self.PhiObsHigh[iZ], kind='linear', bounds_error=False, fill_value=0.)
         self.fPhiObsLowMeas[iZ] = interp1d(self.LObs[iZ], self.PhiObsLow[iZ], kind='linear', bounds_error=False, fill_value=0.)
      


   def plotLf(self, lfs=None):
      '''Reproduces fig 8 in Sobral+12.
      '''
      if lfs is None:
         lfs = [self]



      L_cgs = np.logspace(np.log10(1.e40), np.log10(1.e44), 501, 10.) # [erg/s]
      L_Lsun = L_cgs / self.convertLumUnit('cgs')

      if self.name=='sobral12halpha':
         # Intrinsic LF,
         # ie corrected for the dust attenuation
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         colors = ['gray', 'g', 'r', 'b']
         for iZ in range(self.nZ):
            z = self.Z[iZ]
            factor = self.LInt[iZ] * np.log(10.) * self.U.bg.h**3
            #
            # measured LF
            yerr = np.vstack(((self.PhiInt[iZ]-self.PhiIntLow[iZ])*factor, (self.PhiIntHigh[iZ]-self.PhiInt[iZ])*factor))
            ax.errorbar(self.LInt[iZ], self.PhiInt[iZ]*factor, yerr=yerr, fmt='o', c=colors[iZ], label=r'$z=$'+str(z))
            #
            # error band
            ax.fill_between(self.LInt[iZ], self.PhiIntLow[iZ]*factor, self.PhiIntHigh[iZ]*factor, edgecolor='',facecolor=colors[iZ], alpha=0.3)
            #
            # Schechter fit
            y = self.phiInt(z, L_Lsun) * self.convertLfUnit('cgs')  * L_cgs * np.log(10.) * self.U.bg.h**3
            ax.plot(L_cgs, y, c=colors[iZ])
         #
         ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((10.**(40.5), 10.**(44)))
         ax.set_ylim((10.**(-5.5), 10.**(-0.5)))
         ax.set_xlabel(r'$L_{H_\alpha}$ [erg/s]')
         ax.set_ylabel(r'$\text{log}_{10} \left( \Phi \times \text{Mpc}^3 \right)$')
         ax.set_title(r'Intrinsic luminosity function')
         #
         #path = './figures/profile/Sobral12/'+'lf_intrinsic_sobral12_fig8.pdf'
         path = './figures/lf/lf_halpha_intrinsic_sobral12_fig8.pdf'
         fig.savefig(path, bbox_inches='tight')
         fig.clf()
         #plt.show()
         #plt.close()



      # "Observed" LF,
      # ie affected by dust attenuation
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      colors = ['gray', 'g', 'r', 'b']
      for iZ in range(self.nZ):
         z = self.Z[iZ]
         factor = self.LObs[iZ] * np.log(10.) * self.U.bg.h**3
         #
         # measured LF
         yerr = np.vstack(((self.PhiObs[iZ]-self.PhiObsLow[iZ])*factor, (self.PhiObsHigh[iZ]-self.PhiObs[iZ])*factor))
         ax.errorbar(self.LObs[iZ], self.PhiObs[iZ]*factor, yerr=yerr, fmt='o', c=colors[iZ], label=r'$z=$'+str(z))
         #
         # error band
         ax.fill_between(self.LObs[iZ], self.PhiObsLow[iZ]*factor, self.PhiObsHigh[iZ]*factor, edgecolor='',facecolor=colors[iZ], alpha=0.3)
         #
         # Schechter fit
         y = self.phi(z, L_Lsun) * self.convertLfUnit('cgs')  * L_cgs * np.log(10.) * self.U.bg.h**3
         ax.plot(L_cgs, y, c=colors[iZ])
         #
#         # my high and low curves
#         #y = self.phiObs(z, L) + self.fPhiObsHighMeas[iZ](L) - self.fPhiObsMeas[iZ](L)
#         #ax.plot(L, y * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         #y = self.phiObs(z, L) + self.fPhiObsLowMeas[iZ](L) - self.fPhiObsMeas[iZ](L)
#         #ax.plot(L, y * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         #
#         # varying the Schechter fits within 1 sigma
#         ax.plot(L, self.phi(z, L, obs=True, lStar='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         ax.plot(L, self.phi(z, L, obs=True, lStar='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         ax.plot(L, self.phi(z, L, obs=True, phiStar='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         ax.plot(L, self.phi(z, L, obs=True, phiStar='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         ax.plot(L, self.phi(z, L, obs=True, alpha='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
#         ax.plot(L, self.phi(z, L, obs=True, alpha='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((10.**(40.5), 10.**(44)))
      ax.set_ylim((1.e-7, 1.))
      ax.set_xlabel(r'$L_{H_\alpha}$ [erg/s]')
      ax.set_ylabel(r'$\text{log}_{10} \left( \Phi \times \text{Mpc}^3 \right)$')
      ax.set_title(r'``Observed" luminosity function')
      #
      #path = './figures/profile/Sobral12/'+'lf_observed.pdf'
      path = './figures/lf/lf_halpha_sobral12.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()
      #plt.close()



##################################################################################
##################################################################################


class LFHaCochrane17(LF):
   '''Halpha luminosity function from Cochrane+17
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'cochrane17halpha'
      self.nameLatex = r'C17 H$\alpha$'
      self.refLatex = 'C17'
      self.lineName = 'halpha'
      self.lineNameLatex = r'H$\alpha$'
      self.lambdaMicrons = 656.28e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 2
      #
      self.Z = np.array([0.810, 1.466, 2.231])
      self.nZ = len(self.Z)
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      #
      LStar = 10.**np.array([42.12, 42.56, 42.87]) / 3.839e33  # convert from [erg/s] to [Lsun]
      LStarLow = 10.**np.array([42.12 - 0.02, 42.56 - 0.05, 42.87 - 0.06]) / 3.839e33
      LStarHigh = 10.**np.array([42.12 + 0.03, 42.56 + 0.06, 42.87 + 0.08]) / 3.839e33
      #
      PhiStar = 10.**np.array([-2.31, -2.61, -2.78]) / U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarLow = 10.**np.array([-2.31 - 0.05, -2.61 - 0.09, -2.78 - 0.09]) / U.bg.h**3
      PhiStarHigh = 10.**np.array([-2.31 + 0.04, -2.61 + 0.08, -2.78 + 0.08]) / U.bg.h**3
      #
      Alpha = np.array([-1.6, -1.62, -1.59]) # [dimless]
      AlphaLow = np.array([-1.6 - 0.2, -1.62 - 0.29, -1.59 - 0.13])
      AlphaHigh = np.array([-1.6 + 0.2, -1.62 + 0.25, -1.59 + 0.12])

#      # Schechter fits to the intrinsic luminosity functions
#      # table 5
#      path = './input/LIM_literature/Sobral+12/table5/Sobral+12_table5.txt'
#      data = np.genfromtxt(path)
#      z = data[:,0]
#      #
#      self.Z = z
#      self.zMin = np.min(z)
#      self.zMax = np.max(z)
#      #
#      LStar = 10.**data[:,1] / 3.839e33  # convert from [erg/s] to [Lsun]
#      #LStarHigh = 10.**(data[:,1]+data[:,2]) / 3.839e33  # convert from [erg/s] to [Lsun]
#      #LStarLow = 10.**(data[:,1]+data[:,3]) / 3.839e33  # convert from [erg/s] to [Lsun]
#      #
#      PhiStar = 10.**data[:,4] / U.bg.h**3  # [(Mpc/h)^-3]
#      #PhiStarHigh = 10.**(data[:,4]+data[:,5]) / U.bg.h**3  # [(Mpc/h)^-3]
#      #PhiStarLow = 10.**(data[:,4]+data[:,6]) / U.bg.h**3  # [(Mpc/h)^-3]
#      #
#      # for alpha, they recommend using -1.6, rather than the best fit...
#      #Alpha = -1.6 * np.ones_like(z) # [dimless]
#      Alpha = data[:,7] # [dimless]
#      #AlphaHigh = data[:,7]+data[:,8] # [dimless]
#      #AlphaLow = data[:,7]+data[:,9] # [dimless]

      # interpolate the Schechter fits (Eq 2)
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #self.lStarHigh = interp1d(self.Z, LStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.lStarLow = interp1d(self.Z, LStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
      #self.phiStarHigh = interp1d(self.Z, PhiStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.phiStarLow = interp1d(self.Z, PhiStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #self.alphaHigh = interp1d(self.Z, AlphaHigh, kind='linear', bounds_error=False, fill_value=0.)
      #self.alphaLow = interp1d(self.Z, AlphaLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      # Intrinsic LF, with dust attenuation removed
      self.phiInt = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]
      # Observed LF, where the luminosities experience dust extinction
      self.phiObs = lambda z,l: 100.**(1./5.) * self.phiInt(z, l * 100.**(1./5.))   # [(Mpc/h)^-3 / Lsun]

      self.phi = self.phiObs


      super(LFHaCochrane17, self).__init__(U)






##################################################################################
##################################################################################


class LFOiiiMehta15(LF):
   '''[OIII] luminosity function from Mehta+15 arXiv:1505.07843v2
   This is really the total luminosity in the [OIII] doublet
   4959A and 5007A.
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'mehta15oiii'
      self.nameLatex = r'M15 O{\sc iii}'
      self.refLatex = 'M15'
      self.lineName = 'oiii'
      self.lineNameLatex = r'[O{\sc iii}]'
      self.lambdaMicrons = 495.9e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the observed, dust attenuated LF
 
      # OIII LF: tables 2 and 3
      zMin = np.array([0.8, 1.85])  # lower bin edge
      zMax = np.array([1.2, 2.2])   # higher bin edge
      self.Z = 0.5 * (zMin + zMax)
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      #
      Alpha = np.array([-1.42, -1.57]) 
      AlphaLow = np.array([-1.42 - 0.43, -1.57 - 0.77]) 
      AlphaHigh = np.array([-1.42 + 0.23, -1.57 + 0.28]) 
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #
      LStar = 10.**np.array([42.21, 42.55]) / 3.839e33  # convert from [erg/s] to [Lsun]
      LStarLow = 10.**np.array([42.21 - 0.18, 42.55 - 0.19]) / 3.839e33  # convert from [erg/s] to [Lsun]
      LStarHigh = 10.**np.array([42.21 + 0.22, 42.55 + 0.28]) / 3.839e33  # convert from [erg/s] to [Lsun]
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #
      PhiStar = 10.**np.array([-3.17, -2.69]) / U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarLow = 10.**np.array([-3.17 - 0.39, -2.69 - 0.51]) / U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarHigh = 10.**np.array([-3.17 + 0.27, -2.69 + 0.31]) / U.bg.h**3  # [(Mpc/h)^-3]
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      # Observed LF, where the luminosities experience dust extinction
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFOiiiMehta15, self).__init__(U)



   def phiJ(self, lHa, lOIII): 
      '''Joint LF for Ha and OIII
      lHa, lOIII [Lsun]
      [number / (Mpc/h)^3 / Lsun^2]
      '''
      # convert from Lun to erg/sec
      lHa *= 3.839e33
      lOIII *= 3.839e33
      # take the log10
      llHa = np.log10(lHa)
      llOIII = np.log10(lOIII)

      # Joint bivariate OIII - Ha LF: table 1
      alphaJ = -1.5 
      alphaJLow = -1.5 - 0.2 
      alphaJHigh = -1.5 + 0.5 
      #
      lStarJ = 10.**42.1   # [erg/s]
      lStarJLow = 10.**(42.1 - 0.2) # [erg/s]
      lStarJHigh = 10.**(42.1 + 0.1) # [erg/s]
      #
      phiStarJ = 10.**-2.95   # [Mpc^-3]
      phiStarJLow = 10.**(-2.95 - 0.18)   # [Mpc^-3]
      phiStarJHigh = 10.**(-2.95 + 0.33)  # [Mpc^-3]
      #
      beta = 1.13
      betaLow = 1.13 - 0.26
      betaHigh = 1.13 + 0.06
      #
      r = 0.28
      rLow = 0.28 - 0.20
      rHigh = 0.28 + 0.35
      # 
      # sLog10LHa = sLnLHa / np.log(10.)
      sLnLHa = 0.92
      sLnLHaLow = 0.92 - 0.11
      sLnLHaHigh = 0.92 + 0.08
      #
      l0 = 10.**40. # [erg/s]

      # mean Ha luminosity, given the OIII luminosity
      # Eq 3 in Mehta+15. Not used, and inconsistent with Eq 5
      # because of a ln VS log10 issue.
      #lHaMean = lambda lOIII: l0 * r * (lOIII / l0)**beta   # [erg/s]
      # Eq 5 in Mehta+15
      llHaMean = lambda llOII: np.log10(l0) + np.log10(r) + beta * (llOIII - np.log10(l0))

      # bivariate LF
      # [/Mpc^3 / log10(erg/s)^2]
      result = np.log(10.) 
      # the phiStarJ is missing from Eq 4, but is quoted in table 1
      # and seems needed by dimensional analysis,
      # and to reproduce fig 2.
      result *= phiStarJ
      result *= (lOIII / lStarJ)**(alphaJ+1.)
      result *= np.exp(-lOIII / lStarJ)
      result *= np.log(10.) / sLnLHa / np.sqrt(2. * np.pi)
      # here I assume that Mehta+15 did not realize that log<lHa> <> <log lHa>
      #result *= np.exp(- 0.5 * (llHa - np.log10(lHaMean(lOIII)))**2 * np.log(10)**2 / sLnLHa**2)
      result *= np.exp(- 0.5 * (llHa - llHaMean(lOIII))**2 * np.log(10)**2 / sLnLHa**2)
      # convert to [/(Mpc/h)^3 / Lsun^2]
      result /= self.U.bg.h**3
      result *= 3.839e33**2   # convert from [/(erg/s)^2] to [/Lsun^2]
      result /= np.log(10.) * lHa 
      result /= np.log(10.) * lOIII 
      
      return result


   def computeCorrCoeffHaOIII(self):
      '''Integrate the bivariate LF
      to compute the correlation coefficient between the Ha and OIII
      line luminosities.
      This is a relevant quantity for LIM cross-correlations.
      [dimless]
      '''

      # moments of the joint bivariate LF
      def integrand(lnlHa, lnlOIII, nHa=0., nOIII=0.):
         '''input luminosities [Lsun]
         '''
         lHa = np.exp(lnlHa)
         lOIII = np.exp(lnlOIII)
         result = self.phiJ(lHa, lOIII)   # [/(Mpc/h)^3 / Lsun^2]
         result *= lHa**(1. + nHa) * lOIII**(1. + nOIII)   # [Lsun^(nHa + nOIII) / (Mpc/h)^3]
         return result

      # normalization
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=0., nOIII=0.)
      nGal = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      print("nGal = "+str(nGal)+" Mpc/h^-3")

      # <lHa * lOIII>
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=1., nOIII=1.)
      meanLHaLOIII = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      meanLHaLOIII /= nGal
      print("<lHa * lOIII> = "+str(meanLHaLOIII)+" Lsun^2")

      # <lHa>
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=1., nOIII=0.)
      meanLHa = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      meanLHa /= nGal
      print("<lHa> = "+str(meanLHa)+" Lsun")

      # <lOIII>
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=0., nOIII=1.)
      meanLOIII = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      meanLOIII /= nGal
      print("<lOIII> = "+str(meanLOIII)+" Lsun")

      # cov[lHa, lOIII] = <lHa * lOIII> - <lHa> * <lOIII>
      covLHaLOIII = meanLHaLOIII - meanLHa * meanLOIII
      print("cov[lHa, lOIII] = "+str(covLHaLOIII)+" Lsun^2")

      # <lHa**2>
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=2., nOIII=0.)
      meanSqLHa = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      meanSqLHa /= nGal
      print("<lHa^2> = "+str(meanSqLHa)+" Lsun^2")

      # var[lHa]
      varLHa = meanSqLHa - meanLHa**2
      print("var[lHa] = "+str(varLHa)+" Lsun^2")

      # <lOIII**2>
      f = lambda lnlHa, lnlOIII: integrand(lnlHa, lnlOIII, nHa=0., nOIII=2.)
      meanSqLOIII = integrate.dblquad(f, np.log(self.lMin), np.log(self.lMax), lambda x: np.log(self.lMin), lambda x: np.log(self.lMax), epsabs=0., epsrel=1.e-3)[0]
      meanSqLOIII /= nGal
      print("<lOIII^2> = "+str(meanSqLOIII)+" Lsun^2")

      # var[lOIII]
      varLOIII = meanSqLOIII - meanLOIII**2
      print("var[lOIII] = "+str(varLOIII)+" Lsun^2")

      # r = cov[lHa, lOIII] / sqrt(var[lHa] * var[lOIII])
      rLHaLOIII = covLHaLOIII / np.sqrt(varLHa * varLOIII)
      print("r[lHa, lOIII] = "+str(rLHaLOIII))

      return rLHaLOIII


   ##################################################################################
   
   def plotBivariateLf(self):
      '''Reproduce fig 2 in Mehta+15
      '''

      # log10(lOIII / (erg/s))
      xMin = 40.5
      xMax = 43.5
      nX = 100  # nb of cells. nb of edges is nX+1
      dX = (xMax-xMin)/nX
      xCenters = np.linspace(xMin, xMax, nX)
      xEdges = np.linspace(xMin-0.5*dX, xMax+0.5*dX, nX+1)
      # luminosities in Lsun
      l_Lsun = 10.**xCenters / 3.839e33

      xx,yy = np.meshgrid(xEdges, xEdges, indexing='ij') # log10(lOIII / (erg/s)), log10(lHa / (erg/s))
      lOIII, lHa = np.meshgrid(l_Lsun, l_Lsun, indexing='ij')  # [Lsun]
      
      # compute bivariate LF
      data = self.phiJ(lHa, lOIII)  # [number / (Mpc/h)^3 / Lsun^2]
      # convert to [number / Mpc^3 / (log10(L/(erg/s)))^2]
      data *= self.U.bg.h**3
      data /= 3.839e33**2   # convert from [/Lsun^2] to [/(erg/s)^2]
      data *= np.log(10.) * lHa 
      data *= np.log(10.) * lOIII 
      # take the log10
      data = np.log10(data)

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      cp=ax.pcolormesh(xx, yy, data, linewidth=0, rasterized=True, cmap=plt.cm.jet)
      cp.set_clim(-5., -2.)
      fig.colorbar(cp, label=r'log$_{10}\left(\text{nb/Mpc}^3/\text{log}_{10}(L/(\text{erg/s}))^2\right)$')
      #
      plt.axis('scaled')
      ax.set_xlabel(r'$\text{log}_{10}\left( L_\text{OIII} / \text{(erg/s)} \right)$')
      ax.set_ylabel(r'$\text{log}_{10}\left( L_{\text{H}\alpha} / \text{(erg/s)} \right)$')


      plt.show()



##################################################################################
##################################################################################


class LFHaColbert13(LF):
   '''Halpha luminosity functions from Colbert+13 
   These LFs have been corrected for completeness, but are the observed,
   dust-extincted ones, as needed for LIM.
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'colbert13halpha'
      self.nameLatex = r'C13 H$\alpha$'
      self.refLatex = 'C13'
      self.lineName = 'halpha'
      self.lineNameLatex = r'H$\alpha$'
      self.lambdaMicrons = 656.28e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 4
      zMin = np.array([0.3, 0.9])  # lower bin edge
      zMax = np.array([0.9, 1.5])   # higher bin edge
      self.Z = 0.5 * (zMin + zMax)
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      #
      Alpha = np.array([-1.27, -1.43]) 
      AlphaLow = np.array([-1.27 - 0.12, -1.43 - 0.17]) 
      AlphaHigh = np.array([-1.27 + 0.12, -1.43 + 0.17]) 
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #
      LStar = 10.**np.array([41.72, 42.18]) / 3.839e33  # convert from [erg/s] to [Lsun]
      LStarLow = 10.**np.array([41.72 - 0.09, 42.18 - 0.10]) / 3.839e33
      LStarHigh = 10.**np.array([41.72 + 0.09, 42.18 + 0.10]) / 3.839e33
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #
      PhiStar = 10.**np.array([-2.51, -2.70]) / U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarLow = 10.**np.array([-2.51 - 0.11, -2.70 - 0.12]) / U.bg.h**3
      PhiStarHigh = 10.**np.array([-2.51 + 0.11, -2.70 + 0.12]) / U.bg.h**3
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      # Observed LF, where the luminosities experience dust extinction
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFHaColbert13, self).__init__(U)

      

##################################################################################
##################################################################################


class LFOiiiColbert13(LF):
   '''[OIII] luminosity functions from Colbert+13 
   Unclear from the paper, but I believe this is for the sum
   of the luminosities in the [OIII] doublet at 4959A and 5007A.
   These LFs have been corrected for completeness, but are the observed,
   dust-extincted ones, as needed for LIM.
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'colbert13oiii'
      self.nameLatex = r'C13 O{\sc iii}'
      self.refLatex = 'C13'
      self.lineName = 'oiii'
      self.lineNameLatex = r'[O{\sc iii}]'
      self.lambdaMicrons = 495.9e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 4
      zMin = np.array([0.7, 1.5])  # lower bin edge
      zMax = np.array([1.5, 2.3])   # higher bin edge
      self.Z = 0.5 * (zMin + zMax)
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      #
      Alpha = np.array([-1.40, -1.67]) 
      AlphaLow = np.array([-1.40 - 0.15, -1.67 - 0.78]) 
      AlphaHigh = np.array([-1.40 + 0.15, -1.67 + 0.78]) 
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #
      LStar = 10.**np.array([42.34, 42.91]) / 3.839e33  # convert from [erg/s] to [Lsun]
      LStarLow = 10.**np.array([42.34 - 0.06, 42.91 - 0.37]) / 3.839e33
      LStarHigh = 10.**np.array([42.34 + 0.06, 42.91 + 0.37]) / 3.839e33
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #
      PhiStar = 10.**np.array([-3.19, -3.74]) / U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarLow = 10.**np.array([-3.19 - 0.09, -3.74 - 0.43]) / U.bg.h**3
      PhiStarHigh = 10.**np.array([-3.19 + 0.09, -3.74 + 0.43]) / U.bg.h**3
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      # Observed LF, where the luminosities experience dust extinction
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFOiiiColbert13, self).__init__(U)




##################################################################################
##################################################################################


class LFCiiPopping16(LF):
   '''[Cii] luminosity functions from Popping+16, table 2.
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'popping16cii'
      self.nameLatex = r'P16 [C{\sc ii}]'
      self.refLatex = 'P16'
      self.lineName = 'cii'
      self.lineNameLatex = r'[C{\sc ii}]'
      self.lambdaMicrons = 158.  # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 2
      self.Z = np.array([0., 1., 2., 3. 4., 6.])
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      #
      Alpha = np.array([-1.25, -1.43, -1.52, -1.41, -1.53, -1.77]) 
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #
      LStar = 10.**np.array([7.47, 7.66, 7.81, 7.80, 7.85, 7.80])  # [Lsun]
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #
      PhiStar = 10.**np.array([-2.33, -2.15, -2.20, -2.12, -2.37, -2.95]) / U.bg.h**3  # [(Mpc/h)^-3]
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      # Observed LF
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFCiiPopping16, self).__init__(U)



##################################################################################
##################################################################################


class LFCOPopping16(LF):
   '''[CO j->j-1] luminosity functions from Popping+16, table 1.
   '''

   def __init__(self, U, j):
      # choose the line CO j->j-1
      self.j = j
      # required attributes
      self.name = 'popping16co'+str(j)+'-'+str(j-1)
      self.nameLatex = r'P16 [C0'+str(j)+'-'+str(j-1)+']'
      self.refLatex = 'P16'
      self.lineName = 'co'+str(j)+'-'+str(j-1)
      self.lineNameLatex = r'[C0'+str(j)+'-'+str(j-1)+']'
      # frequency of j->j-1 transition:
      self.nuHz = 115.271208e9 * j   # [Hz]
      self.lambdaMicrons = 299792458. / self.nuHz * 1.e6 # [mu]

      # Luminosity bounds 
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 1
      self.Z = np.array([0., 1., 2., 3. 4., 6.])
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      
      if self.j==1:
         Alpha = np.array([-1.36, -1.49, -1.52, -1.71, -1.94]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([6.97, 7.25, 7.30, 7.26, 6.99])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.85, -2.73, -2.63, -2.94, -3.46]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)

      if self.j==2:
         Alpha = np.array([-1.35, -1.47, -1.52, -1.75, -2.00]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([7.54, 7.84, 7.92, 7.89, 7.62])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.85, -2.72, -2.66, -3.00, -3.56]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      if self.j==3:
         Alpha = np.array([-1.29, -1.47, -1.53, -1.76, -2.00]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([7.83, 8.23, 8.36, 8.26, 7.95])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.81, -2.79, -2.78, -3.11, -3.60]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      if self.j==4:
         Alpha = np.array([-1.29, -1.45, -1.51, -1.80, -2.03]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([8.16, 8.50, 8.64, 8.70, 8.23])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.93, -2.84, -2.85, -3.45, -3.78]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)

      if self.j==5:
         Alpha = np.array([-1.20, -1.47, -1.45, -1.76, -1.95]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([8.37, 8.80, 8.74, 8.73, 8.30])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.94, -3.03, -2.80, -3.34, -3.67]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)

      if self.j==6:
         Alpha = np.array([-1.15, -1.41, -1.43, -1.73, -1.93]) 
         self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
         #
         LStar = 10.**np.array([8.38, 8.74, 8.77, 8.84, 8.38])  # [Jy * km/s * Mpc^2]
         LStar *= self.nuHz / 299792458.e-3  # convert to [Jy * Hz * Mpc^2] = [1.e-23 * erg/s /cm^2 * Mpc^2]
         LStar *= 3.086e24**2 # convert to [1.e-23 * erg/s]
         LStar *= 1.e-23 / 3.846e33   # convert to [Lsun]
         self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
         #
         PhiStar = 10.**np.array([-2.92, -2.92, -2.80, -3.40, -3.72]) / U.bg.h**3  # [(Mpc/h)^-3]
         self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)

      # Observed LF
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFCOPopping16, self).__init__(U)



##################################################################################
##################################################################################


class LFLyaCassata11(LF):
   '''Lyman-alpha luminosity functions from Cassata+11, table 2.
   '''

   def __init__(self, U):
      # required attributes
      self.name = 'popping16cii'
      self.nameLatex = r'P16 [C{\sc ii}]'
      self.refLatex = 'P16'
      self.lineName = 'cii'
      self.lineNameLatex = r'[C{\sc ii}]'
      self.lambdaMicrons = 121.567e3  # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]

      # Luminosity bounds 
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]


      # Schechter fits to the intrinsic luminosity functions
      # table 2
      zMin = np.array([1.95, 3., 4.55])  # lower bin edge
      zMax = np.array([3., 4.55, 6.6])   # higher bin edge
      self.Z = 0.5 * (zMin + zMax)
      self.zMin = np.min(self.Z)
      self.zMax = np.max(self.Z)
      self.nZ = len(self.Z)
      #
      Alpha = np.array([-1.6, -1.78, -1.69]) 
      self.alpha = interp1d(self.Z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      #
      # luminosities after IGM absorption
      LStar = 10.**np.array([42.70, 42.70, 42.72]) / 3.839e33  # convert from [erg/s] to [Lsun]
      self.lStar = interp1d(self.Z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      #
      PhiStar = 10.**np.array([7.1, 4.8, 9.2]) / U.bg.h**3  # [(Mpc/h)^-3]
      self.phiStar = interp1d(self.Z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
     
      # Observed LF
      self.phi = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 / Lsun]

      super(LFLyaCassata11, self).__init__(U)









##################################################################################
##################################################################################



class LFEGG(LF):
   '''Various line luminosity functions from Sobral+12 arXiv:1202.3436v2
   '''

   def __init__(self, U, lineName='halpha'):
      self.U = U

      # required attributes
      self.lineName = lineName
      self.name = 'egg' + self.lineName
      self.refLatex = 'EGG'
      self.nameLatex = ('EGG ' + self.lineName).replace('_', ' ')

      # Luminosity bounds 
      # Arbitrary here, copied from Sobral+12
      # the measurements only span 1.e40 to 1.e44 [erg/sec]
      self.lMin = 1.e35 / 3.839e33  # convert from [erg/sec] to [Lsun]
      self.lMax = 1.e44 / 3.839e33  # convert from [erg/sec] to [Lsun]
      #self.phi(l, z) = #  [number / (Mpc/h)^3 / (luminosity unit)]



      # line setup: read from egg
      # read line parameters from egg
      pathIn = './input/egg_results/'
      # read redshift and mean number density of galaxies [(mpc/h)^{-3}]
      d1 = np.loadtxt(pathIn+'ngal.txt')
      self.Z = d1[:,0]
      self.nZ = len(self.Z)
      self.nGal = interp1d(self.Z, d1[:,1], kind='linear', bounds_error=False, fill_value=0.)
      # find the corresponding line number
      lineNames = np.loadtxt(pathIn+'line_names.txt', dtype='str')
      nLines = len(lineNames)
      iLine = np.where(lineNames==self.lineName)[0][0]
      # find the line wavelength and frequency
      self.lambdaMicrons = np.loadtxt(pathIn+'line_lambda_microns.txt')[iLine] # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [hz]
      # read mean galaxy luminosity [lsun]
      d2 = np.loadtxt(pathIn+'mean_gal_lum.txt')
      self.meanGalLum = interp1d(self.Z, d2[:,iLine], kind='linear', bounds_error=False, fill_value=0.)
      # read total galaxy luminosity density # [lsun / (mpc/h)^3]
      d3 = np.loadtxt(pathIn+'total_gal_lum_density.txt')
      self.lumDensity = interp1d(self.Z, d3[:,iLine], kind='linear', bounds_error=False, fill_value=0.)

      # read fractional covariance of galaxy line luminosities [dimless]
      d4 = np.loadtxt(pathIn+'s2ij.txt')
      d4 = d4.reshape((self.nZ, nLines, nLines))
      self.s2ij = interp1d(self.Z, d4, axis=0, kind='linear', bounds_error=False, fill_value=0.)      
   
      #
      self.s2 = interp1d(self.Z, d4[:, iLine, iLine], kind='linear', bounds_error=False, fill_value=0.)      
      self.nGalEff = lambda z: self.nGal(z) / (1. + self.s2(z))
      #self.sqLumDensity = 
      self.pShot = lambda z: self.meanIntensity(z)**2 / self.nGalEff(z)

      #super(LFEGG, self).__init__(U)


   ##################################################################################
   # Quantities derived from the moments of the luminosity function
   # !!! had to be repeated here to remove the lMin and lMax arguments inside,
   # which break everything

   def meanIntensity(self, z, lMin=0., lMax=np.inf):
      '''[Lsun / (Mpc/h)^2 / sr / Hz]
      '''
      result = self.lumDensity(z)  # [Lsun / (Mpc/h)^3]
      result *= 3.e5 / self.U.hubble(z)   # *[Mpc/h]
      result /= 4. * np.pi * self.nuHz # *[/sr/Hz]
      return result


   def meanGalLum(self, z, lMin=0., lMax=np.inf):
      '''Mean galaxy luminosity [Lsun]
      This quantity may be luminosity-cutoff dependent,
      depending on the LF
      '''
      result = self.lumDensity(z)  # [Lsun / (Mpc/h)^3]
      result /= self.nGal(z, lMin=lMin, lMax=lMax) # / [(Mpc/h)^{-3}] = [Lsun]
      return result


   def pShot(self, z, lMin=0., lMax=np.inf):
      '''[(intensity unit)^2 (Mpc/h)^3] = [Lsun^2 / (Mpc/h) / sr^2 / Hz^2]
      '''
      result = self.sqLumDensity(z)  # [Lsun^2 / (Mpc/h)^3]
      result *= (3.e5 / self.U.hubble(z))**2   # *[(Mpc/h)^2]
      result /= (4. * np.pi * self.nuHz)**2 # *[1/sr^2/Hz^2]
      return result


   def nGalEff(self, z, lMin=0., lMax=np.inf):
      '''Effective galaxy number density for the shot noise [(Mpc/h)^{-3}]
      '''
      result = self.meanIntensity(z)**2
      result /= self.pShot(z)
      return result


