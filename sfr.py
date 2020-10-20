from headers import *

class Sfr(object):

   def __init__(self, U, MassFunc):
      self.U = U
      self.MassFunc = MassFunc

      # required attributes:
      # self.name
      # self.Z
      # self.sfr(m [Msun/h], z) [Msun/yr]

      # SFRD interpolation
      #z = np.logspace(np.log10(1.e-3), np.log10(10.), 501, 10.)
      sfrd = np.array(map(self.sfrdForInterp, self.Z))
      self.sfrd = interp1d(self.Z, sfrd, kind='linear', bounds_error=False, fill_value=0.)


   def sfrdForInterp(self, z, alpha=1, bias=False):
      '''SFR density:
      \int dm dn/dm SFR(m)^alpha [(Msun/yr)^alpha / (Mpc/h)^3]
      '''
      def integrand(lnm):
         '''the mass in integral is in [Msun/h]
         '''
         m = np.exp(lnm)
         result =  m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m, z)**alpha
         if bias:
            result *= self.MassFunc.fb1(m, 1./(1.+z))
         return result

      # [Msun /yr / (Mpc/h)^3]
      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
      return result

   
   def sfrdBehroozi13(self, z):
      '''[Msun/yr/(Mpc/h)^3]
      '''
      # from Fonseca+16
      # [Msun/yr/Mpc^3]
      result = 0.180 / (10.**(-0.997*(z-1.243)) + 10.**(0.241*(z-1.243))) 
      # convert from [Msun/yr/Mpc^3] to [Msun/yr/(Mpc/h)^3]
      result /= self.U.bg.h**3
      return result

   def plotSfrd(self):
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      Sfrd = np.array(map(self.sfrd, self.Z))
      Sfrd *= self.U.bg.h**3  # convert to [Msun/yr/Mpc]
      ax.plot(self.Z, Sfrd, 'b--', label=self.name)
      #
      # Behroozi+13 parametrization
      z = np.linspace(0., 5., 101)
      ax.plot(z, self.sfrdBehroozi13(z) * self.U.bg.h**3, 'r-', label=r'Behroozi+13')
      #
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc='lower center')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'SFRD$(z)$/[M$_\odot$/yr / Mpc$^3$]')
      #
      #fig.savefig(pathFig + "sfrd.pdf", bbox_inches='tight')
      plt.show()
      fig.clf()


   def nHEff(self, z, alpha1=1, alpha2=1):
      '''Effective halo number density [(Mpc/h)^-3]
      which determines the amplitude of the halo shot noise,
      ie the 1-halo term
      '''
      result = self.sfrdForInterp(z, alpha=alpha1) * self.sfrdForInterp(z, alpha=alpha2)
      result /= self.sfrdForInterp(z, alpha=alpha1+alpha2)
      return result



   def plotnHEff(self):

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # range of reasonable scalings for L = SFR^alpha
      A = [0.8, 1.0, 1.1]
      for a in A:
         f = lambda z: self.nHEff(z, alpha1=a, alpha2=a)
         n = np.array(map(f, self.Z))
         ax.plot(self.Z, n, label=r'$\alpha_i=\alpha_j=$'+str(round(a, 1)))
      #
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(self.Z), np.max(self.Z)))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{n}^h_\text{eff}$ [(Mpc/$h$)$^{-3}$]')
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      #
      path = './figures/sfr/nheff_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()

      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # range of reasonable scalings for L = SFR^alpha
      A = [0.8, 1.0, 1.1]
      for a in A:
         f = lambda z: self.nHEff(z, alpha1=a, alpha2=a)
         n = np.array(map(f, self.Z))
         ax.plot(self.Z, 1./n, label=r'$\alpha_i=\alpha_j=$'+str(round(a, 1)))
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(self.Z), np.max(self.Z)))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$P^\text{1h}$ [(Mpc/$h$)$^3$]')
      #
      path = './figures/sfr/p1h_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()




   def plotNHEffSpherex(self):
      # SPHEREx voxel size
      # the spectral resolution power is R=40 for the lower z, and 150 for the high z
      R = 40.
      # hence the redshift size of the voxel
      dz = (1. + self.Z) / R
      # and the comoving depth of the voxel
      dChi = dz * 3.e5 / self.U.hubble(self.Z)   # [Mpc/h]
      # angular pixel size: 6.2 arcsec
      thetaPix = 6.2 * np.pi/(180.*3600.)
      # hence the voxel comoving volume
      vVoxSpherex = (self.U.bg.comoving_distance(self.Z) * thetaPix)**2 * dChi  # [(Mpc/h)^3]
      print "vVoxSpherex=", vVoxSpherex, "(Mpc/h)^3"


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # range of reasonable scalings for L = SFR^alpha
      A = [0.8, 1.0, 1.1]
      for a in A:
         f = lambda z: self.nHEff(z, alpha1=a, alpha2=a) 
         n = np.array(map(f, self.Z)) * vVoxSpherex
         n *= vVoxSpherex
         ax.plot(self.Z, n, label=r'$\alpha_i=\alpha_j=$'+str(round(a, 1)))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlim((np.min(self.Z), np.max(self.Z)))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{N}^h_\text{eff}$ per SPHEREx voxel')
      #
      path = './figures/sfr/halo_sparsity_spherex_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()

   

   def bEff(self, z, alpha=1.):
      '''Effective halo bias [dimless]
      bEff = int dm n(m) SFR(m)^alpha b(m) / SFRD,
      SFRD = int dm n(m) SFR(m)
      '''
      result = self.sfrdForInterp(z, alpha=alpha, bias=True)
      result /= self.sfrdForInterp(z, alpha=1, bias=False)
      return result


   def plotBEff(self):

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # range of reasonable scalings for L = SFR^alpha
      A = [0.8, 1.0, 1.1]
      for a in A:
         f = lambda z: self.bEff(z, alpha=a) 
         b = np.array(map(f, self.Z))
         ax.plot(self.Z, b, label=r'$\alpha$='+str(round(a, 1)))
      #
      ax.set_xlim((np.min(self.Z), np.max(self.Z)))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$b^h_\text{eff}$')
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      #
      path = './figures/sfr/bheff_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def dbEff2dlnm(self, m, z, alpha=1.):
      '''d(b_eff^2)/dlnM [dimless]
      This gives the contribution of each halo mass to 
      the 2halo power spectrum
      '''
      result = 2. * m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m, z)**alpha
      result *= self.MassFunc.fb1(m, 1./(1.+z))
      result /= self.sfrdForInterp(z, alpha=alpha)
      result *= self.bEff(z, alpha=alpha)
      return result

   def plotdbEff2dlnm(self):

      M = np.logspace(np.log10(1.e10), np.log10(2.e14), 101, 10.) # [Msun/h]
      Z = np.array([0.001, 1., 2., 5.])


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         z = Z[iZ]
         f = lambda m: self.dbEff2dlnm(m, z)
         y = np.array(map(f, M))
         ax.plot(M, y, c=plt.cm.autumn_r(iZ/(len(Z)-1.)), label=r'$z=$'+str(int(z)))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'Halo mass $m$ [$M_\odot/h$]')
      ax.set_ylabel(r'$\frac{d b^h_\text{eff}(z)}{d\text{ln} m}$', fontsize=14)
      #
      path = './figures/sfr/dbheff2dlnm_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()




   def dP1hdlnm(self, m, z, alpha1=1, alpha2=1):
      '''dP1h/dlnm [(Mpc/h)^3]
      '''
      result = m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m, z)**(alpha1+alpha2)
      result /= self.sfrdForInterp(z, alpha=alpha1) * self.sfrdForInterp(z, alpha=alpha2)
      return result


   def plotdP1hdlnm(self):
      
      M = np.logspace(np.log10(1.e10), np.log10(2.e14), 101, 10.) # [Msun/h]
      Z = np.array([0.001, 1., 2., 5.])


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         z = Z[iZ]
         f = lambda m: self.dP1hdlnm(m, z)
         y = np.array(map(f, M))
         ax.plot(M, y, c=plt.cm.autumn_r(iZ/(len(Z)-1.)), label=r'$z=$'+str(int(z)))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'Halo mass $m$ [$M_\odot/h$]')
      ax.set_ylabel(r'$\frac{dP^\text{1h}(z)}{d\text{ln} m}$  [(Mpc/h)$^3$]', fontsize=14)
      #
      path = './figures/sfr/dp1hdlnm_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


##################################################################################


class SfrFonseca16(Sfr):

   def __init__(self, U, MassFunc):
      self.name = 'Fonseca16'

      # from Fonseca+16 (1607.05288v2).
      # Table I
      Z = np.array([0., 1., 2., 3., 4., 5.])
      self.Z = Z.copy()
      M0 = np.array([3.0e-10, 1.7e-9, 4.0e-9, 1.1e-8, 6.6e-8, 7.0e-7])  # [Msun/yr]
      Mb = np.array([6.0e10, 9.0e10, 7.0e10, 5.0e10, 5.0e10, 6.0e10])   # [Msun]
      Mc = np.array([1.0e12, 2.0e12, 2.0e12, 3.0e12, 2.0e12, 2.0e12])   # [Msun]
      a = np.array([3.15, 2.9, 3.1, 3.1, 2.9, 2.5])
      b = np.array([-1.7, -1.4, -2.0, -2.1, -2.0, -1.6])
      c = np.array([-1.7, -2.1, -1.5, -1.5, -1.0, -1.0])
      # interpolate values
      self.fM0 = interp1d(Z, M0, kind='linear', bounds_error=False, fill_value=0.) # [Msun/yr]
      self.fMb = interp1d(Z, Mb, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
      self.fMc = interp1d(Z, Mc, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
      self.fa = interp1d(Z, a, kind='linear', bounds_error=False, fill_value=0.)
      self.fb = interp1d(Z, b, kind='linear', bounds_error=False, fill_value=0.)
      self.fc = interp1d(Z, c, kind='linear', bounds_error=False, fill_value=0.)
      # below Eq 11
      self.Ma= 1.e8   # [Msun]

      super(SfrFonseca16, self).__init__(U, MassFunc)


   def sfr(self, m, z):
      ''' SFR [Msun/yr] as a function of halo mass [Msun/h] and redshift
      from Fonseca+16 (1607.05288v2), Eq 11.
      I inferred the units from Eq 9 and 11.
      input halo mass m is Mvir [Msun/h]
      '''
      # bounds for the fitting function
      if z<0. or z>5.:
         return 0.

      # convert from [Msun/h] to [Msun]
      m /= self.U.bg.h

      # Eq 11
      result = self.fM0(z)
      result *= (m/self.Ma)**self.fa(z)
      result *= (1.+m/self.fMb(z))**self.fb(z)
      result *= (1.+m/self.fMc(z))**self.fc(z)

      return result


##################################################################################
#!!!!! WARNING:
#This SFR(m,z) relation from Moster+13 is somewhat sketchy,
#in that m represents the "halo mass at z=0",
#not the halo mass at that redshift.
#So it is difficult for me to use, without knowing how to convert the halo mass
#at a given redshift into a halo mass at z=0

#class SfrMoster13(Sfr):
#
#   def __init__(self, U, MassFunc):
#      self.name = 'Moster13'
#      self.Z = np.linspace(0., 9., 101)
#
#      super(SfrMoster13, self).__init__(U, MassFunc)
#
#
#
#   def sfr(self, m, z):
#      ''' SFR [Msun/yr] as a function of halo mass [Msun/h] and redshift
#      from Moster+13 Eq 17-20
#      input halo mass m is Mvir [Msun/h]
#      '''
#      # bounds for the fitting function
#      if z<0. or z>100.:
#         return 0.
#      
#      # convert from Mvir to M200c [Msun/h]
#      # in Eq 17-20, they call their halo mass Mvir,
#      # but at the start of the paper they explain 
#      # that they actually mean M200crit
#      #m, r = self.U.massRadiusConversion(m, z, value=200, ref='c')
#      m, r = self.U.massRadiusConversion(m, 0., value=200, ref='c')
#
#      # convert from [Msun/h] to [Msun]
#      m /= self.U.bg.h
#
#      f10 = 2.658 # [Msun/yr]
#      f11 = 11.937
#      f12 = 0.424
#      #
#      f20 = 5.507
#      f21 = 2.437
#      #
#      f30 = -0.915
#      f31 = 0.696
#      f32 = -0.159
#
#      f1 = f10 * np.exp(-0.5 * (np.log10(m)-f11)**2 / f12**2)  # [Msun/yr]
#      f2 = f20 + f21 * np.log10(m/1.e12)
#      f3 = 10.**(f30 + f31 * (m/1.e12)**f32)
#
#      result = f1 / (1.+z)**f2 * np.exp(f2 / f3 * z / (1.+z))
#      return result
#
#
#   def sfrdForInterp(self, z, alpha=1):
#      '''SFR density:
#      \int dm dn/dm SFR(m)^alpha [(Msun/yr)^alpha / (Mpc/h)^3]
#      '''
#      def integrand(lnm):
#         '''the mass in integral is in [Msun/h]
#         '''
#         m = np.exp(lnm)
#         #return m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m / self.U.bg.h, z)
#         return m * self.MassFunc.fmassfunc(m, 1./(1.+0.)) * self.sfr(m, z)**alpha
#      # [Msun /yr / (Mpc/h)^3]
#      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
#      # [Msun /yr / Mpc^3]
#      #result *= self.U.bg.h**3 
#      return result





##################################################################################


class SfrMoster13Speagle14(Sfr):
   '''Uses the halo mass to stellar mass relation in Moster+13,
   and the star forming main sequence from Speagle+14
   '''

   def __init__(self, U, MassFunc, scatter=True, nProc=3, save=False):
      self.name = 'Moster13Speagle14'
      self.nProc = nProc
      self.U = U
      # values of z to precompute SFR
      self.zMin = 0.
      self.zMax = 10.
      self.Nz = 101
      self.Z = np.linspace(self.zMin, self.zMax, self.Nz)
      
      # decide whether to take into account
      # the scatters in mVir-mStar and mStar-SFR relations
      if not scatter:
         self.sfr = self.sfrNoScatter
      else:

         # path to precompute SFR
         self.pathSfr = './output/sfr/sfr_moster13speagle14.txt'
         
         # integration bounds, for the scatter in mStar and sfr
         self.mStarMin = 1.e8 # [Msun]
         self.mStarMax = 1.e13   # [Msun
         self.sfrMin = 1.e-2  # [Msun/yr]
         self.sfrMax = 1.e4   # [Msun/yr]

         # values of mVir to precompute SFR 
         self.mMin = 1.e10  # in (h^-1 solarM)
         self.mMax = 1.e16  # in (h^-1 solarM)
         self.Nm = 201 # nb of m pts
         self.M = np.logspace(np.log10(self.mMin), np.log10(self.mMax), self.Nm, 10.) # masses in h^-1 solarM

         if save:
            self.save()
         self.load()

      super(SfrMoster13Speagle14, self).__init__(U, MassFunc)



# I believe that Minf is M_{200crit} in the Moster+13 notation

# According to Moster+13, sec 3.3.2,
# The stellar mass in the central galaxies is obtained 
# when inputing M_{200crit} for the halo mass.
# The stellar mass of the satellite galaxies is obtained
# when inputing the infall mass and infall redshift for the subhalo

#   def mStar(Minf,z):
#       """
#       moster(Minf,z):
#       Returns the stellar mass (M*) given Minf and z from Table 1 and
#       Eq. (2,11-14) of Moster++13 [1205.5807].
#       This version works in terms of Msun units, not Msun/h.
#       To get "true" stellar mass, add 0.15 dex of lognormal scatter.
#       To get "observed" stellar mass, add between 0.1-0.45 dex extra scatter.
#       """
#       zzp1  = z/(1+z)
#       M1    = 10.0**(11.590+1.195*zzp1)
#       mM    = 0.0351 - 0.0247*zzp1
#       beta  = 1.376  - 0.826*zzp1
#       gamma = 0.608  + 0.329*zzp1
#       Mstar = 2*mM / ( (Minf/M1)**(-beta) + (Minf/M1)**gamma )
#       Mstar*= Minf
#       return Mstar

   def meanMStar(self, m, z):
      """
      Moster+13
      m: halo virial mass [Msun/h]
      z: redshift
      returns mStar, the stellar mass of the central galaxy [Msun]
      From Table 1 and Eq. (2,11-14) of Moster++13 [1205.5807].
      To get "true" stellar mass, add 0.15 dex of lognormal scatter.
      To get "observed" stellar mass, add between 0.1-0.45 dex extra scatter.
      """
      # convert virial mass to M_{200c}, as in Moster+13
      # they call their halo mass Mvir,
      # but at the start of the paper they explain 
      # that they actually mean M200crit
      m, r = self.U.massRadiusConversion(m, z, value=200, ref='c')
      # convert from [Msun/h] to [Msun]
      m /= self.U.bg.h

      zzp1  = z/(1+z)
      M1    = 10.0**(11.590+1.195*zzp1)
      mM    = 0.0351 - 0.0247*zzp1
      beta  = 1.376  - 0.826*zzp1
      gamma = 0.608  + 0.329*zzp1
      Mstar = 2*mM / ( (m/M1)**(-beta) + (m/M1)**gamma )
      Mstar*= m  # [Msun]
      return Mstar


   def pMStarGivenMVir(self, mStar, mVir, z):
      '''Moster+13
      proba(mStar | mVir, z)
      As suggested in Moster+13,
      assume a lognormal distribution for stellar mass
      given halo mass, with a scatter of 0.15 dex
      '''
      sLog10MStar = 0.15
      result = np.exp(-0.5 * (np.log10(mStar) - np.log10(self.meanMStar(mVir, z)))**2 / sLog10MStar**2)
      result /= np.sqrt(2. * np.pi) * sLog10MStar
      result /= mStar * np.log(10.)
      return result


   def meanSfr(self, mStar, z):
      '''Speagle+14
      returns the mean SFR [Msun/yr]
      given the stellar mass of the central galaxy [Msun]
      The scatter around this mean is 
      0.2 dex lognormal for the true SFR
      0.3 dex log normal for the observed SFR
      '''
      # age of universe at z
      t = self.U.bg.time(z)   # [Gyr]
      # fitting function
      result = (0.84 - 0.026 * t) * np.log10(mStar)
      result -= 6.51 - 0.11 * t
      result = 10.**result
      return result
      


   def pSfrGivenMStar(self, sfr, mStar, z):
      '''Speagle+14
      proba(sfr + mStar, z)
      The scatter around the mean SFR, at a fixed stellar mass, is 
      0.2 dex lognormal for the true SFR
      (0.3 dex log normal for the observed SFR)
      '''
      sLog10Sfr = 0.2
      result = np.exp(-0.5 * (np.log10(sfr) - np.log10(self.meanSfr(mStar, z)))**2 / sLog10Sfr**2)
      result /= np.sqrt(2. * np.pi) * sLog10Sfr
      result /= sfr * np.log(10.)
      return result


   def sfrNoScatter(self, mVir, z):
      '''SFR [Msun/yr] in a halo
      with virial mass mVir [Msun/h]
      at redshift z,
      ignoring the scatter in the mVir-mStar and mStar-SFR relations
      '''
      return meanSfr(meanMStar(mVir, z), z)


   def sfrForInterp(self, mVir, z, alpha=1.):
      '''Mean SFR [Msun/yr] for a halo of mass mVir [Msun/h] at z,
      including the scatter in the mVir-mStar relation
      and the scatter in the mStar-sfr relation
      '''
      def integrand(lnSfr, lnMStar):
         sfr = np.exp(lnSfr)
         mStar = np.exp(lnMStar)
         result = sfr * mStar
         result *= self.pMStarGivenMVir(mStar, mVir, z)
         result *= self.pSfrGivenMStar(sfr, mStar, z)
         result *= sfr**alpha
         return result

      result = integrate.dblquad(integrand, np.log(self.mStarMin), np.log(self.mStarMax), lambda x: np.log(self.sfrMin), lambda x: np.log(self.sfrMax), epsabs=0., epsrel=1.e-3)[0]
      return result



   def save(self):
      print "Computing SFR(mVir, z)"
      Sfr = np.zeros((self.Nm, self.Nz))

      # Evaluate in parallel
      with sharedmem.MapReduce(np=self.nProc) as pool:
         for iz in range(self.Nz):
            z = self.Z[iz]
            # compute mass function
            f = lambda m: self.sfrForInterp(m, z, alpha=1.)
            Sfr[:, iz] = np.array(pool.map(f, self.M))
            print "- done "+str(iz+1)+" of "+str(self.Nz)

      np.savetxt(self.pathSfr, Sfr)
      return


   def load(self):
      print "Loading SFR(mVir, z)"
      self.Sfr = np.genfromtxt(self.pathSfr)

      # interpolate
#      interp_sfr = RectBivariateSpline(np.log(self.M), self.Z, np.log(self.Sfr), s=0)
#      self.sfr = lambda m, z: (z>=self.zMin and z<=self.zMax) * np.exp( interp_sfr(np.log(m), z)[0,0] )
      interp_sfr = RectBivariateSpline(np.log(self.M), self.Z, self.Sfr, s=0, kx=1, ky=1)
      self.sfr = lambda m, z: (z>=self.zMin and z<=self.zMax) * interp_sfr(np.log(m), z)[0,0]


   def sfrdForInterpNoScatter(self, z, alpha=1):
      '''SFRD calculation that neglects the scatters
      in mVir-mStar and mStar-SFR.
      Used as a check of the importance of the scatter
      \int dm dn/dm SFR(m)^alpha [(Msun/yr)^alpha / (Mpc/h)^3]
      '''
      def integrand(lnm):
         '''the mass in integral is in [Msun/h]
         '''
         m = np.exp(lnm)
         result = m * self.MassFunc.fmassfunc(m, 1./(1.+z)) 
         result *= self.meanSfr(self.meanMStar(m, z), z)**alpha
         return result
      # [Msun /yr / (Mpc/h)^3]
      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
      return result


   ##################################################################################

   
   def plotSfrd(self):
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # Including scatter
      Sfrd = np.array(map(self.sfrd, self.Z))
      Sfrd *= self.U.bg.h**3  # convert to [Msun/yr/Mpc]
      ax.plot(self.Z, Sfrd, 'b-', label=self.name)
      #
      # neglecting scatter in mVir-mStar and mStar-SFR
      SfrdNoScatter = np.array(map(self.sfrdForInterpNoScatter, self.Z))
      SfrdNoScatter *= self.U.bg.h**3  # convert to [Msun/yr/Mpc]
      ax.plot(self.Z, SfrdNoScatter, 'b--', label=self.name+' no scatter')
      #
      # Behroozi+13 parametrization
      z = np.linspace(0., 5., 101)
      ax.plot(z, self.sfrdBehroozi13(z) * self.U.bg.h**3, 'r-', label=r'Behroozi+13')
      #
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc='lower center')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'SFRD$(z)$/[M$_\odot$/yr / Mpc$^3$]')
      #
      #fig.savefig(pathFig + "sfrd.pdf", bbox_inches='tight')
      plt.show()
      fig.clf()


   def plotSfr(self):
      '''Reproduce fig 8 in Speagle+14
      '''

      Z = np.array([0., 0.25, 0.5, 1., 2., 4.])
      MStar = np.logspace(9.6, 11.2, 101, 10.)

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         f = lambda mStar: self.meanSfr(mStar, z)
         sfr = np.array(map(f, MStar))
         plt.plot(MStar, sfr)
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$M_\star$ [$M_\odot$]')
      ax.set_ylabel(r'SFR [$M_\odot$/yr]')

      plt.show()

