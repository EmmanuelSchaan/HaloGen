from headers import *

##################################################################################
# All radii and volumes are comoving
# Densities are per unit comoving volume


class Profile(object):
   """takes a Universe class U as input
   required variable: use_correction_factor, rMax for the truncation radius, mMin, mMax
   required function: u
   optional functions I would like: rho3d, rho2d, rho3dFourier, rho3dFourierNorm
   """
   
   def __init__(self, U):
      # copy U
      self.U = U

   def uN(self, k, m, z, n=1):
      return self.u(k,m,z)**n

   def plotU(self, z=0., func=None):
      if func is None:
         func = self.u
      M = np.logspace(np.log10(1.e11), np.log10(1.e16), 6, 10.)
      K = np.logspace(np.log10(3.e-2), np.log10(3.e3), 501, 10.)
      
      plt.figure(0)
      ax = plt.subplot(111)
      #
      for im in range(len(M))[::-1]:
         m = M[im]
         f = lambda k: func(k, m, z) #self.u(k, m, z)
         U = np.array(map(f, K))
         ax.loglog(K, U, lw=2, label=r'$M=$'+floatExpForm(m, round=2)+'$M_\odot/h$')
      #
      ax.legend(loc=3, framealpha=0.05)
      ax.set_xlim((min(K), max(K)))
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$u(k)$ [dimensionless]')
      
      plt.show()
   
   ##################################################################################

   def rho2d(self, R, m, z):
      """project the 3d into the 2d profile
      """
      f = lambda r: self.rho3d(r, m, z) * 2.*r / np.sqrt(r**2 - R**2)
      result = integrate.quad(f, R, np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result

   def convolveGauss(self, R, m, z, fwhm):
      """convolve the 2d profile with a Gaussian
      the fwhm should be in Mpc/h comoving, same as R
      """
      sigmaBeam = fwhm / np.sqrt(8.*np.log(2.))
      # do the smoothing
      def f(x):
         if np.isfinite(i0(x*X/sigmaBeam**2)):
            result = x * self.rho2d(x) * np.exp(-0.5*(x**2+X**2)/sigmaBeam**2) * i0(x*X/sigmaBeam**2) / sigmaBeam**2
         else:
            result = 0.
         return result
      result = integrate.quad(f, 0., np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result

   def rhoFourier(self, k, m, z):
      """Fourier transform of 3d profile (same as 2d profile)
      """
      # do Fourier transform
      integrand = lambda r: 4.*np.pi*r**2 * self.rho3d(r, m, z) * j0(k * r)
      result = integrate.quad(integrand, 0., np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result


   def fourier(self, func, k):
      # do Fourier transform
      integrand = lambda r: 4.*np.pi*r**2 * func(r) * j0(k * r)
      result = integrate.quad(integrand, 0., np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result

   def inverseFourier(self, func, r):
      integrand = lambda k: k**2 * func(k) * j0(k * r) / np.pi
      result = integrate.quad(integrand, 0., np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result



##################################################################################
##################################################################################

class ProfNFW(Profile):

   def __init__(self, U, trunc=4.):
      super(ProfNFW, self).__init__(U)
      self.LoadNonLinMass()
      self.trunc = trunc   # truncation radius, in units of rVir
      self.use_correction_factor = True
      self.mMin = 0.
      self.mMax = np.inf
   
   def __str__(self):
      return "nfwdensity"
   
   def LoadNonLinMass(self):
      """precompute the nonlinear mass at z=0.
      """
      print "Loading non-lin mass at z=0"
      z = 0.
      self.m_nonlin = self.U.nonLinMass(z)


   def rS_rhoS_c(self, m, z):
      """comoving scale radius for NFW profile
      in Mpc/h
      """
      Rvir = self.U.frvir(m, z)
      # concentration parameter
      #c = 10./(1.+z) * (m / self.m_nonlin)**(-0.2)   # from Takada & Jain 2002
      c = 9./(1.+z) * (m / self.m_nonlin)**(-0.13) # Takada & Jain 2003
      # scale radius
      RS = Rvir / c  # in Mpc/h
      # normalize the mass within rVir to be mVir
      rhoS = m / (4.*np.pi*RS**3)
      rhoS /= np.log(1.+c) - c/(1.+c)  # (Msun/h) / (Mpc/h)^3
      return RS, rhoS, c
   
   
   def totalMass(self, trunc=None):
      """total mass within truncation radius
      trunc in units of rVir
      mass in Msun/h
      if trunc=infinity, the mass is infinite
      """
      if trunc is None:
         trunc = self.trunc
      rVir = self.U.frvir(m, z)
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      # truncation radius over scale radius
      xMax = trunc * rVir/rS
      result = 4./3. * np.pi * rS**3 * rhoS
      result = xMax - np.log(1 + xMax)
      return result
   

   def rho3d(self, r, m, z, trunc=None):
      """comoving 3d density, truncated at trunc*rVir
      in (Msun/h) / (Mpc/h)^3
      r comoving in Mpc/h
      m is Mvir in Mpc/h
      """
      if trunc is None:
         trunc = self.trunc
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      x = r / rS  # dimensionless
      # 3d profile
      result = rhoS
      if x <= trunc * c:
         result /= x*(1.+x)**2
      else:
         result = 0.
      return result


#!!!!!!!!!!! if the truncation radius is infinite, there is an analytical expression. I should use it.
   def rho2d(self, R, m, z, trunc=None):
      """comoving 2d density
      where the 3d NFW is truncated at trunc*rVir
      in (Msun/h) / (Mpc/h)^2
      R in Mpc/h
      x = R/RS = theta/thetaS, dimensionless
      """
      if trunc is None:
         trunc = self.trunc
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      x = r / rS  # dimensionless
      # projected density in (Msun/h) / (Mpc/h)^2
      result = 2.*rhoS*RS
      # truncate the 3d profile at rvir
      if X>trunc*c:
         result *= 0.
      else:
         f = lambda x: x * 1./(x*(1.+x)**2) / np.sqrt(x**2 - X**2)
         result *= integrate.quad(f, X, trunc*c, epsabs=0., epsrel=1.e-3)[0]
      return result


   def nfw(self, k, m, z):
      """Fourier transform of NFW density profile,
      normalized such that u(k=0, m, z) = 1
      ie rhoNFW(k,m,z) = m * nfw(k,m,z)
      truncation radius is taken to be infinite (unphysical)
      k in h Mpc^-1
      m in h^-1 solarM
      """
      RS, rhoS, c = self.rS_rhoS_c(m, z)
      #
      result = np.sin(k * RS) * (  Si((1+c) * k * RS) - Si(k * RS)  )
      result += - np.sin(c * k * RS) / ((1+c) * k * RS)
      result += np.cos(k * RS) * (  Ci((1+c) * k * RS) - Ci(k * RS)  )
      result /= (np.log(1+c) - c/(1+c))
      return result
   
   
   def rho3dFourier(self, k, m, z, trunc=None):
      """approximate Fourier transform:
      take the exact fourier transform for trunc=infty
      and just rescale by the total mass for finite trunc
      """
      if trunc is None:
         trunc = self.trunc
      result = self.nfw(k, m, z)
      result *= self.totalMass(trunc=trunc)
      return result


   def u(self, k, m, z):
      """returns m/\bar{rho} * NFW u, in (h^-1 Mpc)^3
      k in h Mpc^-1, m in h^-1 solarM
      this is the quantity typically used in the halo model
      """
      return self.nfw(k, m, z) * m / self.U.rho_m(z)
   

   def phiOfR(self, r, m, z):
      """gravitational potential
      function of comoving radius in Mpc/h
      """
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      x = r/rS
      result = -4.*np.pi*self.U.G
      result *= rhoS * rS**2
      result *= np.log(1.+x) / x
      return result


   def test_nfw(self):
      
      M = np.logspace(np.log10(1.e11), np.log10(1.e16), 6, 10.)
      K = np.logspace(np.log10(3.e-2), np.log10(3.e3), 501, 10.)
      z=0.
      
      plt.figure(0)
      ax = plt.subplot(111)
      #
      for im in range(len(M)):
         m = M[im]
         #f = lambda k: self.u(k, m, z) / (m/self.U.rho_m(z))
         f = lambda k: self.nfw(k, m, z)
         Y = np.array(map(f, K))
         ax.loglog(K, Y, label=r'$M=$'+floatExpForm(m)+r'$M_\odot/h$')
      #
      ax.legend(loc=3)
      ax.set_xlim((min(K), max(K)))
      ax.set_ylim((1.e-3, 2.5))
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'u(k), normalized to 1 for k=0')
      ax.grid()
      
      plt.show()


##################################################################################
##################################################################################

class ProfNFWPlusGaussian(Profile):
   """DM in NFW profile,
   gas in Gaussian profile with std dev = Rvir,
   consistent with the kSZ detection in Schaan Ferraro +16
   """
   
   def __init__(self, U, fGauss=0.17):
      super(ProfNFWPlusGaussian, self).__init__(U)
      # instantiate profiles to combine
      self.ProfNFW = ProfNFW(U)
      # fraction of the total mass in the Gaussian profile
      self.fGauss = fGauss
      self.mMin = 0.
      self.mMax = np.inf
   
   def __str__(self):
      return "nfwplusgaussdensity"
   
   
   def gauss(self, k, m, z):
      """Gaussian density profile
      k in h Mpc^-1, m in h^-1 solarM
      """
      Rvir = self.U.frvir(m, z)
      result = np.exp(-0.5 * Rvir**2 * k**2)
      return result
   
   def u(self, k, m, z):
      """returns m/\bar{rho} * NFW u, in (h^-1 Mpc)^3
      k in h Mpc^-1, m in h^-1 solarM;
      """
      result = self.fGauss * self.gauss(k, m, z)
      result += (1. - self.fGauss) * self.nfw(k, m, z)
      result *= m / self.U.rho_m(z)
      return result


##################################################################################
##################################################################################

class ProfHOD(Profile):
   """generic HOD
   """
   
   def __init__(self, U, MassFunc):
      super(ProfHOD, self).__init__(U)
      self.MassFunc = MassFunc
      self.use_correction_factor = False
      self.ProfNFW = ProfNFW(U)
      
      # interpolate nBarGal for speed
      A = self.MassFunc.A.copy()
      nBarGal = np.array(map(self.nBarGalForInterp, A))
      f = UnivariateSpline(A, nBarGal, k=1, s=0)
      self.nBarGal = lambda a: f(a) * (a>=np.min(A))*(a<=np.max(A))

   
   ##################################################################################
   
   
   def Ncen(self, m):
      """number of central galaxies per halo
      """
      pass
   
   def Nsat(self, m):
      """number of satellite galaxies per halo
      """
      pass
   
   def Ngal(self, m):
      """total number of galaxies per halo
      """
      return self.Ncen(m) + self.Nsat(m)
   
   def plotHOD(self):
      M = np.logspace(np.log10(1.e11), np.log10(1.e15), 201, 10.)
      Ncen = np.array(map(self.Ncen, M))
      Nsat = np.array(map(self.Nsat, M))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(M, Ncen+Nsat, 'k', lw=4, label=r'$N_\text{total}$')
      ax.plot(M, Ncen, 'b', lw=2, label=r'$N_\text{cen}$')
      ax.plot(M, Nsat, 'r', lw=2, label=r'$N_\text{sat}$')
      #
      ax.legend(loc=4)
      ax.set_xlim((1.e11, 1.e16))
      ax.set_ylim((1.e-1, 1.e2))
      ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'halo mass [$M_\odot$/h]')
      ax.set_ylabel(r'mean number of galaxies per halo $\langle N_\text{gal} \rangle$')
      #
      #fig.savefig("./figures/cib_penin12/hod.pdf", bbox_inches='tight')
      
      plt.show()
   
   ##################################################################################
   
   
   def nBarGalForInterp(self, a):
      """number of galaxies per unit volume
      """
      integrand = lambda lnm: np.exp(lnm)*self.MassFunc.fmassfunc(np.exp(lnm), a)*self.Ngal(np.exp(lnm))
      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0, epsrel=1.e-2)[0]
      return result
   
   def plotNBarGalIntegrand(self):
      """checking for convergence of mass integral: very low masses required
      """
      Z = np.array([0., 0.5, 1., 2., 5.])
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         a = 1./(1.+z)
         integrand = lambda lnm: np.exp(lnm)*self.MassFunc.fmassfunc(np.exp(lnm), a)*self.Ngal(np.exp(lnm))
         M = self.MassFunc.M.copy()
         Integrand = np.array(map(integrand, np.log(M)))
         #
         ax.plot(M, Integrand, lw=2, label=r'$z=$'+str(z))
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_xlabel(r'halo mass $M$ [$M_\odot/h$]')
      ax.set_ylabel(r'$\frac{d \bar{n}_\text{gal}}{d \ln M}$ [(Mpc/h)$^{-3}$]')
      #
      #fig.savefig("./figures/cib_penin12/dngal_dlnm.pdf", bbox_inches='tight')
      
      plt.show()
   
   
   def plotNBarGal(self):
      Z = np.linspace(0., 10., 201)
      #Z = np.logspace(np.log10(0.01), np.log10(10.), 101, 10.)
      nBarGal = np.array(map(self.nBarGal, 1./(1.+Z)))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, nBarGal, lw=3)
      #
      #ax.set_xscale('log')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'redshift $z$')
      ax.set_ylabel(r'$\bar{n}_\text{gal}$ [(Mpc/h)$^{-3}$]')
      #
      #fig.savefig("./figures/cib_penin12/ngal.pdf", bbox_inches='tight')
      
      plt.show()
   
   ##################################################################################
   
   def fPshotNoise(self, k, z):
      """shot noise for the galaxy number density power spectrum
      i.e. 1./ngal
      """
      result = 1./self.nBarGal(1./(1.+z))
      return result
   
   def fTshotNoise(self, k, z):
      """shot noise for the galaxy number density trispectrum
      i.e. 1./ngal^3
      """
      result = 1./self.nBarGal(1./(1.+z))**3
      return result



##################################################################################
##################################################################################

class ProfHODPenin12(ProfHOD):
   """HOD for IR galaxies, used for CIB halo model
   from Penin Dore Lagache Bethermin 2012
   DOI: 10.1051/0004-6361/201117489
   nu should be 217, 353, 545, 857 (in GHz)
   """
   
   def __init__(self, U, MassFunc, nu=217):
      self.nu = nu
      self.U = U
      self.mMin = 0.
      self.mMax = np.inf
      
#      # HOD params, Penin+12
#      self.mMinHod = 10.**11.5 * self.U.bg.h # in Msun/h
#      self.sLogM = 0.65#np.sqrt(0.65)
#      self.mSat = 10. * self.mMinHod # in Msun/h
#      self.aSat = 1.4
      
      # best fit to Planck CIB. Table 1 in Penin+14
      if self.nu==217:
         self.mMinHod = 10.**12. * self.U.bg.h # in Msun/h
         self.aSat = 1.5
      if self.nu==353:
         self.mMinHod = 10.**12.2 * self.U.bg.h # in Msun/h
         self.aSat = 1.7
      if self.nu==545:
         self.mMinHod = 10.**12.6 * self.U.bg.h # in Msun/h
         self.aSat = 1.9
      if self.nu==857:
         self.mMinHod = 10.**12.8 * self.U.bg.h # in Msun/h
         self.aSat = 1.7
      
      # HOD params, Penin+14
      self.sLogM = 0.65
      self.mSat = 10. * self.mMinHod # in Msun/h
   
      super(ProfHODPenin12, self).__init__(U, MassFunc)

   def __str__(self):
      return "hodpenin12"
   
   ##################################################################################
   
   def Ncen(self, m):
      """number of central galaxies per halo, between 0 and 1
      from Penin+12, from Tinker Wetzel 10
      """
      result = np.log10(m) - np.log10(self.mMinHod)
      result /= self.sLogM
      result = 0.5 * (1. + special.erf(result))
      return result
   
   def Nsat(self, m):
      """number of satellite galaxies per halo
      from Penin+12, from Tinker Wetzel 10
      """
      result = np.log10(m) - np.log10(2.*self.mMinHod)
      result /= self.sLogM
      result = 0.5 * (1. + special.erf(result))
      result *= (m/self.mSat)**self.aSat
      return result

   ##################################################################################

   def u(self, k, m, z):
      """Penin+12 effectively assumes that both centrl and satellites
      are NFW distributed, which is weird.
      """
      result = self.ProfNFW.nfw(k, m, z) * self.Ngal(m) / self.nBarGal(1./(1.+z))
      return result
   
   def uN(self, k, m, z, n=1):
      """Penin+12 wants to use
      <Ngal> for the 2-halo and
      <Ngal*(Ngal-1)> for the 1-halo term
      """
      if n==1:
         result = self.u(k, m, z)
      # for the 1-halo profile with an HOD,
      # you don't just square the 2-halo profile,
      # although this is only different at very low mass
      # see Lacasa+13
      elif n==2:
         Ncen = self.Ncen(m)
         Nsat = self.Nsat(m)
         result = Nsat**2 + 2.*Nsat*Ncen
         result *= self.ProfNFW.nfw(k, m, z)**2 / self.nBarGal(1./(1.+z))**2
      else:
         result = self.u(k, m, z)**n
      return result


##################################################################################
##################################################################################

class ProfHODAlam16(ProfHOD):
   """HOD for CMASS from Alam+16.
   Uses halo model of More+13
   """

   def __init__(self, U, MassFunc):
      # HOD params, Alam+16
      self.mCut = 1.77e13  # Msun/h
      self.m1 = 1.51e14 # Msun/h
      self.sigma = 0.8897  # std dev of lnM, dimless
      self.kappa = 0.137
      self.alpha = 1.151
      self.mMin = 0.
      self.mMax = np.inf
   
      super(ProfHODAlam16, self).__init__(U, MassFunc)
   
   def __str__(self):
      return "hodalam16"
   
   ##################################################################################
   
   def Ncen(self, m):
      """number of central galaxies per halo, between 0 and 1
      """
      result = np.log(m/self.mCut)
      result /= np.sqrt(2.) * self.sigma
      result = 0.5 * (1. + special.erf(result))
      return result
   
   def Nsat(self, m):
      """number of satellite galaxies per halo
      """
      result = (m - self.kappa * self.mCut)
      if result>0.:
         result /= self.m1
         result **= self.alpha
         result *= self.Ncen(m)
      else:
         result = 0.
      return result

   ##################################################################################

   def uCen(self, k, m, z):
      """Distribution of central galaxies, no miscentering.
      Here I am guessing, Alam+16 doesn't say
      """
      result = self.Ncen(m) / self.nBarGal(1./(1.+z))
      return result

   def uSat(self, k, m, z):
      """Distribution of satellite galaxies
      following NFW.
      Here I am guessing, Alam+16 doesn't say
      """
      result = self.ProfNFW.nfw(k, m, z)
      result *= self.Nsat(m) / self.nBarGal(1./(1.+z))
      return result

   def u(self, k, m, z):
      """Total galaxy profile,
      including central and satellites
      """
      result = self.uCen(k, m, z)
      result += self.uSat(k, m, z)
      return result

   def uN(self, k, m, z, n=1):
      """More+13,15 convention
      """
      if n==1:
         result = self.u(k, m, z)
      # for the 1-halo profile with an HOD,
      # you don't just square the 2-halo profile,
      # although this is only different at very low mass
      elif n==2:
         uCen = self.uCen(k, m, z)
         uSat = self.uSat(k, m, z)
         result = uSat**2 + 2.*uCen*uSat
      else:
         result = self.u(k, m, z)**n
      return result



##################################################################################
##################################################################################

class ProfHODMore15(ProfHOD):
   """HOD for CMASS
   from More Miyatake +15
   """
   
   def __init__(self, U, MassFunc):
      # HOD params, More+15, average of the 3 columns in Table 1
      self.alphaInc = 0.51
      self.log10mInc = 13.84
      self.log10mMin = 13.42
      self.sLog10m = np.sqrt(0.49)
      self.kappa = 1.10
      self.mMinHod = 10.**(13.42)
      self.m1 = 10.**(14.43)
      self.alpha = 1.09
      self.pOff = 0.357
      self.rOff= 2.3
      
      self.mMin = 0.
      self.mMax = np.inf
      
      super(ProfHODMore15, self).__init__(U, MassFunc)
   
   def __str__(self):
      return "hodmore15"
   
   ##################################################################################
   
   def fInc(self, m):
      """More+15 assume that for a given halo mass,
      a fixed fraction of CMASS galaxies are seen
      how physical is this?
      """
      result = np.min([1., 1.+self.alphaInc*(np.log10(m) - self.log10mInc)])
      result = np.max([0., result])
      return result
   
   def Ncen(self, m):
      """number of central galaxies per halo, between 0 and 1
      """
      result = np.log10(m) - self.log10mMin
      result /= self.sLog10m
      result = 0.5 * (1. + special.erf(result))
      result *= self.fInc(m)
      return result
   
   def Nsat(self, m):
      """number of satellite galaxies per halo
      """
      result = (m - self.kappa * self.mMinHod)
      if result>0.:
         result /= self.m1
         result **= self.alpha
         result *= self.Ncen(m)
      else:
         result = 0.
      return result

   
   ##################################################################################
   
   def uOff(self, k, m, z):
      """Gaussian describing the off-centering
      between matter and the central galaxy
      """
      rS, rhoS, c = self.ProfNFW.rS_rhoS_c(m, z)
      result = k * rS * self.rOff
      result = np.exp(-0.5 * result**2)
      return result
   
   def uCen(self, k, m, z):
      """Distribution of central galaxies,
      to account for miscentering
      """
      result = (1.-self.pOff) + self.pOff * self.uOff(k, m, z)
      result *= self.Ncen(m) / self.nBarGal(1./(1.+z))
      return result
   
   def uSat(self, k, m, z):
      """Distribution of satellite galaxies
      following NFW
      """
      result = self.ProfNFW.nfw(k, m, z)
      result *= self.Nsat(m) / self.nBarGal(1./(1.+z))
      return result
   
   def u(self, k, m, z):
      """Total galaxy profile,
      including central and satellites
      """
      result = self.uCen(k, m, z)
      result += self.uSat(k, m, z)
      return result

   def uN(self, k, m, z, n=1):
      """More+13,15 convention
      """
      if n==1:
         result = self.u(k, m, z)
      # for the 1-halo profile with an HOD,
      # you don't just square the 2-halo profile,
      # although this is only different at very low mass
      elif n==2:
         uCen = self.uCen(k, m, z)
         uSat = self.uSat(k, m, z)
         result = uSat**2 + 2.*uCen*uSat
      else:
         result = self.u(k, m, z)**n
      return result




##################################################################################
##################################################################################

class ProfY(Profile):

   def __init__(self, U):
      super(ProfY, self).__init__(U)
      #
      self.profile_id = "compton_y"
      self.use_correction_factor = False
      self.mMin = 5.e11
      self.mMax = 5.e15
      # parameters for electron pressure profile, from Battaglia et al 2012
      self.alpha = 1.
      self.gamma = -0.3
      #
      self.P00 = 18.1
      self.P0am = 0.154
      self.P0az = -0.758
      #
      self.xc0 = 0.497
      self.xcam = -0.00865
      self.xcaz = 0.731
      #
      self.beta0 = 4.35
      self.betaam = 0.0393
      self.betaaz = 0.415
#      # (Duffy et al 2008)
#      self.cNFW0 = 5.71
#      self.cNFWam = -0.084
#      self.cNFWaz = -0.47
      # link between Pe and Pth
      self.xH = 0.76 # primordial hydrogen mass fraction
      self.Pth_to_Pe = 2.*(self.xH+1.)/(5.*self.xH+3)

   def __str__(self):
      return "pressure"
   
#   # used for the y3d profile
#   def m200c_r200c(self, m, z):
#      # from Duffy et al 2008: different pivot mass
#      cNFW = self.cNFW0 * (m / (2.e12))**self.cNFWam * (1.+z)**self.cNFWaz
#
#      # physical virial radius and scale radius in h^-1 Mpc
#      Rvir = ( 3.*m / (4*np.pi*self.U.rho_crit(z) * self.U.Deltacrit_z(z)) )**(1./3.)
#      Rs = Rvir / cNFW
#
#      # NFW scale density (physical)
#      rhoS = m / (4.*np.pi*Rs**3) / (np.log(1.+cNFW) - cNFW/(1.+cNFW))
#
#      # get R200 and M200
#      f = lambda x: -1. + 1./(1.+x) + np.log(1.+x) - 200./3.*(self.U.rho_crit(z)/rhoS)*x**3
#      x = optimize.brentq(f , 0.1, 100.)
#      R200 = x * Rs  # physical
#      M200 = 4./3.*np.pi*R200**3 * self.U.rho_crit(z) * 200.
#      R200 *= (1.+z) # comoving
#
#      return M200, R200

   
   # dimless normalized thermal3 pressure
   # x=r/r_200c is dimless, m=m_200c in h^-1 solarM
   def P_PDelta_x(self, x, m, z):
      # parameters of the fit, from Battaglia et al 2012
      P0 = self.P00 * (m / (1.e14 * self.U.bg.h))**self.P0am * (1.+z)**self.P0az
      xc = self.xc0 * (m / (1.e14 * self.U.bg.h))**self.xcam * (1.+z)**self.xcaz
      beta = self.beta0 * (m / (1.e14 * self.U.bg.h))**self.betaam * (1.+z)**self.betaaz
      # normalized pressure profile from Battaglia et al. 2012, x = R/Rvir, with R and Rvir comoving
      P_PDelta_x = P0 * (x/xc)**self.gamma / ( 1. + (x/xc)**self.alpha )**beta
      return P_PDelta_x
   
  # y_3d profile
   # k in h Mpc^-1, m=m_vir in Msun/h, y3D(k) in (Mpc/h)^2 (y3D(x) in (Mpc/h)^-1)
   def u(self, k, m, z, cut=False):
      
      # cut from Hill & Pajer 2013
      if (k<1.e-4)or(k>3.):
         return 0.

      # get M200c, R200c
      M200, R200 = self.U.massRadiusConversion(m, z, value=200., ref="c")#self.m200c_r200c(m, z)
      
      # P200, in s^-2 kg / (Mpc/h)
      msun = 1.988435e30 # 1 solar mass in kg
      pc = 3.08567758e16   # 1 parsec in meters
      P200 = self.U.G * (M200*msun)/self.U.bg.h
      P200 *= 200.
      P200 *= self.U.rho_crit(z) * self.U.bg.h**2 * msun/(1.e6*pc)**3 # critical density in kg m^-3 (comoving)
      P200 *= (1.+z)**3 # physical
      P200 *= self.U.bg.Omega0_b/self.U.bg.Omega0_m # baryon fraction
      P200 /= 2. * R200/(1.+z)* 1.e6*pc/self.U.bg.h  # twice R200 in m (physical)
      P200 *= 1.e6*pc/self.U.bg.h
      
      # y3d(x), in (Mpc/h)^-1
      factor = 8.125532e-16   # sigma_T/(me * c**2) in s^2 kg^-1
      y3d_x = lambda x: factor * P200 * self.Pth_to_Pe * self.P_PDelta_x(x, M200, z)
      
      # bound for integral over x=r/r200c
      Rvir = self.U.frvir(m, z)
      xmax = 4. * Rvir / R200
      
      # compute Fourier transform
      # Hill Pajer 13, Eq A9
      f = lambda x: y3d_x(x) * x * np.sin(k * x * R200)

      # if cut=True (default), stop the integral at 4 cvir like in Battaglia et al 2012
      if cut:
         result = integrate.quad(f, 0., xmax, epsrel=1.e-2, epsabs=0)[0]
      else:
         result = integrate.quad(f, 0., np.inf, epsrel=1.e-2, epsabs=0)[0]
      result *= 4.*np.pi* R200**2 / k
      
      return result


   ##################################################################################

   # reproduces fig2 in Battaglia et al 2012
   def testP_PDelta_x(self):
      
      # values of m_200
      M = np.array([6.2e13, 9.5e13, 1.4e14, 2.2e14, 3.45e14, 5.35e14, 8.3e14, 1.29e15])
      # values of redshift
      Z = np.array([0., 0.3, 0.5, 0.7, 1., 1.5])
      # values of r/r_200
      X = np.logspace(np.log10(0.03), np.log10(3.), 50, 10.)
      
      # normalized pressure profile for various masses
      z = 0.
      plt.figure(0)
      ax = plt.subplot(111)
      for im in range(len(M)):
         m = M[im]
         name = '{:.1e}'.format(float(m))
         Y = np.array(map(lambda x: self.P_PDelta_x(x, m, z), X ))
         Y *= X**3
         ax.loglog(X, Y, label=r'$m_{200} = $'+name)
      
      ax.set_xlim((0.03, 3.))
      ax.set_ylim((8.e-4, 0.25))
      ax.grid()
      ax.legend(loc=4)
      ax.set_xlabel(r'$x= r/r_{200}$', fontsize=20)
      ax.set_ylabel(r'$x^3 P_{th}(x) / P_{200}$', fontsize=20)
      
      # normalized pressure profile for various redshifts
      m = 1.4e14
      plt.figure(1)
      ax = plt.subplot(111)
      for iz in range(len(Z)):
         z = Z[iz]
         Y = np.array(map(lambda x: self.P_PDelta_x(x, m, z), X))
         Y *= X**3
         ax.loglog(X, Y, label=r'z = '+str(z))
      ax.set_xlim((0.03, 3.))
      ax.set_ylim((8.e-4, 0.25))
      ax.legend(loc=4)
      ax.grid()
      ax.set_xlabel(r'$x= r/r_{200}$', fontsize=20)
      ax.set_ylabel(r'$x^3 P_{th}(x) / P_{200}$', fontsize=20)
      
      plt.show()
      return


   def testy3d_profile(self):
      M = np.logspace(np.log10(1.e11), np.log10(1.e16), 6, 10.)
      K = np.logspace(np.log10(1.e-3), np.log10(1.e1), 101, 10.)
      z=0.
      
      plt.figure(0)
      ax = plt.subplot(111)
      
      for im in range(len(M)):
         m = M[im]
         f = lambda k: self.u(k, m, z, cut=True)
         tStart = time()
         Y = np.array(map(f, K))
         tStop = time()
         print "took ", tStop-tStart, "sec"
         ax.loglog(K, Y, 'b')
      
      for im in range(len(M)):
         m = M[im]
         f = lambda k: self.u(k, m, z, cut=False)
         tStart = time()
         Y = np.array(map(f, K))
         tStop = time()
         print "took ", tStop-tStart, "sec"
         ax.loglog(K, Y, 'r')
      
      plt.show()


   def compareColin(self):
#!!!!! the disagreement is because I am setting to zero the k>3 h/Mpc, as in Hill Pajer 2013
      # data from Colin
      data = np.genfromtxt("./input/tests/Colin_tSZ/y_l.txt")
      m_vir = data[0, 0]   # [Msun/h]
      z = data[0, 1]
      L = data[:, 2] # dimless
      Y = data[:, 3] # dimless
      
      # my calculation
      a = 1./(1.+z)
      da = self.U.bg.comoving_distance(z)
      K = L / da
      f = lambda k: self.u(k, m_vir, z, cut=True) * a/ da**2.
      Y_me = np.array(map(f, K))
      
      fig=plt.figure(0)
      ax=plt.subplot(211)
      ax.loglog(L, L*Y, 'r', label='Colin')
      ax.loglog(L, L*Y_me, 'b.', label='Manu')
      ax.legend(loc=3)
      ax.set_ylabel(r'$\ell \, y(\ell)$')
      ax.grid()
      #
      ax=plt.subplot(212)
      ax.semilogx(L, Y_me/Y-1., 'b.')
      ax.semilogx(L, np.zeros_like(L), 'k')
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'Manu/Colin-1')
      ax.grid()
      #fig.savefig("./figures/tests/compare_yl_Colin.pdf")
      
      plt.show()



##################################################################################
##################################################################################

class ProfCIBPlanck15(ProfNFW):
   """CIB profile: from Planck XXIII 2015, Planck XXX 2013
   this is an NFW profile, with an extra mass dependence,
   coming from the mass-emissivity relation
   !!!!!!!!!! Doesn't match the plots in the paper... argh
   """
   
   def __init__(self, U, nu=150.e9):
      self.nu = nu
      super(ProfCIBPlanck15, self).__init__(U)

      # constants
      self.c = 3.e8  # m/s
      self.h = 6.63e-34 # SI
      self.kB = 1.38e-23   # SI
      self.Tcmb = 2.726 # K
      
      # typical dust temperature
      # table I of Planck XXIII 2015
      self.Td0 = 24.4
      self.aCIB = 0.36
      self.Td = lambda a: self.Td0 * a**(-self.aCIB)
      
      # SED power index
      # for typical galaxy SED
      # table I of Planck XXIII 2015
      self.bCIB = 1.75  # at low freq
      self.cCIB = 1.70  # at high freq
      
      # L500-M500 relation
      # table I of Planck XXIII 2015
      # mass exponent
      self.eCIB = 1.0
      # z-dependent factor
      self.dCIB = 3.2
      self.psi = lambda a: a**(-self.dCIB)
      # normalization
      # not quoted in the Planck papers, for being only "a normalization parameter" (Planck15 XXII),
      # "which being not physically meaningful will not be further discussed" (Planck13 XXX)
#!!!!!!!!!!!!!!!! factor set to 1 for now!!!
      self.L0 = 2.9e18
#      self.L0 = 7.4e-4

      # interpolating the sed normalization for speed
      Td = np.linspace(1., 100., 201)
      Norm = np.array(map(self.sedNormalization, Td))
      self.sedNormalizationInterp = UnivariateSpline(Td, Norm, k=1, s=0)

      ##################################################################################
      # not needed for this halo model, just for curiosity
      
      # effective linear bias of CIB-emitting galaxies
      # table 8 of Planck XXX 2013
      self.b0 = 0.82
      self.b1 = 0.34
      self.b2 = 0.31
      self.bias = lambda z: self.b0 + self.b1*z + self.b2*z**2
      
      # values of SFR density
      # table 8 of Planck XXX
      # in Msun / yr / Mpc^3
      self.rhoSFR0 = 1.88e-2
      self.rhoSFR1 = 16.07e-2
      self.rhoSFR2 = 16.61e-2
      self.rhoSFR4 = 4.0e-2
      
      # Kennicutt constant, Msun/year
      # from Planck XXIII 2015, from Kennicutt 98
      # in (Msun/h) / yr
      self.K = 1.7e-10 * self.U.bg.h

   
   def __str__(self):
      return "cibplanck15"+str(int(self.nu/1.e9))
   
   def plotTd(self):
      Z = np.linspace(0., 10., 501)
      A = 1./(1.+Z)
      Td = np.array(map(self.Td, A))
   
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, Td, 'k', lw=2)
      #
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'dust temperature $T_d$ [K]')
      #
      #fig.savefig("./figures/cmb/cib_td_planck15xxiii.pdf", bbox_inches='tight')
      
      plt.show()
   
   
   def blackbody(self, nu, T):
      """blackbody function
      input and output in SI
      """
      x = self.h*nu/(self.kB*T)
      result = 2.*self.h*nu**3 /self.c**2
      result /= np.exp(x) - 1.
      return result
   

   def dLnBlackdLnNu(self, nu, T):
      """some derivative of blackbody
      """
      x = self.h*nu/(self.kB*T)
      result = 3.
      if x>=10.:
         result -= x
      else:
         result -= x*np.exp(x) / (np.exp(x)-1.)
      return result
   

   def nu0(self, Td):
      """transition frequency for the typical galaxy SED
      found by matching the log slopes of the SED
      output in Hz
      """
      f = lambda nu: self.bCIB + self.cCIB + self.dLnBlackdLnNu(nu, Td)
      nuMin = 1.e1 * 1.e9  # in Hz
      nuMax = 1.e5 * 1.e9  # in Hz
      '''
      Nu = np.linspace(nuMin, nuMax, 201)
      F = np.array(map(f, Nu))
      plt.semilogx(Nu, F)
      plt.semilogx(Nu, 0.*Nu)
      plt.show()
      '''
      result = optimize.brentq(f , nuMin, nuMax)
      return result
   
   
   def sedRaw(self, nu, Td):
      """typical SED of galaxy contributing to CIB
      called Theta in Planck XXIII 2015
      the normalization of this function is ambiguous from the Planck paper,
      but is very important for the amplitude of the profile, and may introduce a z-dpdce
      I am going to assume it is normalized to unit integral
      nu is the rest-frame frequency
      """
      nu0 = self.nu0(Td)
      result = (nu<nu0) * nu**self.bCIB * self.blackbody(nu, Td)
      result += (nu>=nu0) * (nu/nu0)**(-self.cCIB) * nu0**self.bCIB*self.blackbody(nu0, Td)
      return result
   
   def sedNormalization(self, Td):
      nuMin = 1.  # Hz
      nuMax = 1e5 * 1.e9 # Hz
      result = integrate.quad(lambda nu: self.sedRaw(nu, Td), nuMin, nuMax, epsabs=0, epsrel=1.e-2)[0]
      return result
   
   def sed(self, nu, Td):
      """sed, normalized so that int dnu sed(nu) = 1
      """
      return self.sedRaw(nu, Td) / self.sedNormalizationInterp(Td)
   
#   # here I try a Td-dpdt normalization,
#   # which changes the nu-dpdce of the CIB power spectrum... but still not good
#   def sed(self, nu, Td):
#      nu0 = self.nu0(Td)
#      result = (nu<nu0) * (nu/nu0)**self.bCIB * self.blackbody(nu, Td)
#      result += (nu>=nu0) * (nu/nu0)**(-self.cCIB) * self.blackbody(nu0, Td)
#      return result

   def plotSED(self):
      Nu = np.logspace(np.log10(1.e2), np.log10(1.e4), 501, 10.)*1.e9   # in Hz
      Z = np.linspace(0., 5., 6)
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         a = 1./(1.+z)
         # rest frame sed, normalized to integrate to unity
         f = lambda nu: self.sed(nu, self.Td(a))
         sed = np.array(map(f, Nu))
         ax.plot(Nu/1.e9, sed, c=plt.cm.rainbow((z-min(Z))/(max(Z)-min(Z))), lw=2, label=r'z='+str(z))
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'frequency $\nu$ [GHz]')
      ax.set_ylabel(r'rest frame specific intensity [arbitrary units]')
      #
      #fig.savefig("./figures/cib_planck15/cib_rf_gal_sed_planck15xxiii.pdf", bbox_inches='tight')
      
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         a = 1./(1.+z)
         # observed sed, normalized to unity then redshifted... the relative amplitudes are meaningless now
         f = lambda nu: self.sed(nu/a, self.Td(a))
         sed = np.array(map(f, Nu))
         ax.plot(Nu/1.e9, sed, c=plt.cm.rainbow((z-min(Z))/(max(Z)-min(Z))), lw=2, label=r'z='+str(z))
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'frequency $\nu$ [GHz]')
      ax.set_ylabel(r'redshifted specific intensity [arbitrary units]')
      #
      #fig.savefig("./figures/cib_planck15/cib_obs_gal_sed_planck15xxiii.pdf", bbox_inches='tight')


      plt.show()
   
   
   # M500 in Msun
#!!!!!!!!!!!!!!!!!!!!!!!! value of self.L0? Units?
   # unit of output?
   def L500(self, nu, m500, a):
      result = self.L0 * (m500/1.e14*self.U.bg.h)**self.eCIB
      result *= self.psi(a)
      result *= self.sed(nu/a, self.Td(a))
      return result
   
   # S500
   def S500(self, nu, m500, a):
      chi = self.U.bg.comoving_distance(1./a-1.)
      result = self.L500(nu, m500, a)
      result *= a/(4.*np.pi*chi**2)
      return result
   
   def u(self, k, m, z):
      """k in h/Mpc, m is Mvir in Msun/h
      """
      a = 1./(1.+z)
      m500 = self.U.massRadiusConversion(m, z, value=500., ref="m")[0]
      result = self.S500(self.nu, m500, a)
      result *= self.nfw(k, m, z)
      return result

#   # k in h/Mpc, m is Mvir in Msun/h
## this is my interpretation of what the halo model profile should be...
#   def u(self, k, m, z):
#      a = 1./(1.+z)
#      m500 = self.U.massRadiusConversion(m, z, value=500., ref="m")[0]
#      result = self.L500(self.nu, m500, a)
#      result *= self.nfw(k, m, z)
##      result /= 4.*np.pi
#      return result

   
   
   ##################################################################################
   # Not needed for halo model, just for curiosity
   
   def rhoSFR(self, a):
      """SFR density, in Msun/yr/Mpc^3
      """
      z = 1./a-1.
      if (z<0):
         print "error: asking for rhoSFR at z<0!"
         return 0.
      elif (z<=1):
         z0 = 0.
         z1 = 1.
         rhoSFR0 = self.rhoSFR0
         rhoSFR1 = self.rhoSFR1
      elif (z<=2):
         z0 = 1.
         z1 = 2.
         rhoSFR0 = self.rhoSFR1
         rhoSFR1 = self.rhoSFR2
      else:
         z0 = 2.
         z1 = 4.
         rhoSFR0 = self.rhoSFR2
         rhoSFR1 = self.rhoSFR4
      exponent = np.log(rhoSFR1/self.rhoSFR0)
      exponent /= np.log((1.+z1)/(1.+z0))
      result = rhoSFR0 * (1.+z)**exponent
      return result
   
   
   def j(self, nu, a):
      """emissivity of galaxies in CIB
      """
      result = self.rhoSFR(a)
      result /= a*self.K
      result *= self.sed(nu, a)
      chi = self.comoving_distance(1./a-1.)
      result *= chi**2
      return result
   
   def I(nu):
      """CIB intensity
      """
      f = lambda a: a*self.j(nu, a) / (self.Hubble(a)/3.e5 * a**2)
      aMin = 1./(1.+10.)
      aMax = 1.
      result = integrate.quad(f, aMin, aMax, epsabs=0, epsrel=1.e-5)[0]
      return result


##################################################################################
##################################################################################

class ProfGasBattaglia16(Profile):
   """Tau profile, from Battaglia 2016
   watch typos in paper. This code is correct.
   """
   
   def __init__(self, U, trunc=2.):
      super(ProfGasBattaglia16, self).__init__(U)
      self.use_correction_factor = False
      self.mMin = 0.
      self.mMax = np.inf
      self.trunc = trunc
   
      # parameters from Battaglia 17
      self.xc = 0.5
      self.gamma = -0.2
      
      # m is M200c in Msun/h, z is redshift
      self.rho0 = lambda m,z: 4.e3 * ((m/self.U.bg.h)/1.e14)**0.29 * (1.+z)**(-0.66)
      self.alpha = lambda m,z: 0.88 * ((m/self.U.bg.h)/1.e14)**(-0.03) * (1.+z)**0.19
      self.beta = lambda m,z: 3.83 * ((m/self.U.bg.h)/1.e14)**0.04 * (1.+z)**(-0.025)


   def __str__(self):
      return "taubattaglia17"
   
   
   def rhoFit3d(self, x, m, z):
      """3d scaled gas density profile
      dimensionless
      """
      result = 1. + (x/self.xc)**self.alpha(m, z)
      result **= -(self.beta(m,z)+self.gamma) / self.alpha(m,z)
      result *= (x/self.xc)**(self.gamma)
      result *= self.rho0(m, z)
      return result
   
   def rhoFit2d(self, X, m, z):
      """2d scaled gas density profile
      dimensionless
      """
      f = lambda x:  self.rhoFit3d(x, m, z) * 2.*x / np.sqrt(x**2 - X**2)
      result = integrate.quad(f, X, np.inf, epsabs=0., epsrel=1.e-2)[0]
      return result
   
   
   
   def rho3d(self, r, m, z):
      """3d physical gas density profile
      in (Msun/h) / (Mpc/h)^3
      m is mVir in Msun/h
      r is comoving radius in Mpc/h
      """
      # convert mass from Mvir to M200c (both in [Msun/h])
      # and get comoving R200c [Mpc/h]
      m200c, r200c = self.U.massRadiusConversion(m, z, value=200, ref="c")
      # fitting function from Battaglia 17
      result = self.rhoFit3d(x, m200c, z)
      # rescale with comoving critical density, in (Msun/h) / (Mpc/h)^3
      result *= self.U.rho_crit(z)
      # rescale with the baryon fraction
      result *= self.U.bg.Omega0_b/self.U.bg.Omega0_m
      return result
   
   
   
   def rho3dFourier(self, k, m, z):
      """Fourier transform of 3d density profile,
      k in h Mpc^-1, m in h^-1 solarM
      """
      # truncation radius
      rVir = self.U.frvir(m, z)
      m200c, r200c = self.U.massRadiusConversion(m, z, value=200, ref="c")
      xTrunc = self.trunc * rVir / r200c
      # do Fourier transform
      integrand = lambda x: 4.*np.pi*x**2 * self.rhoFit3d(x, m, z) * j0(k * r200c * x)
      result = integrate.quad(integrand, 0., xTrunc, epsabs=0., epsrel=1.e-3)[0]
      result *= r200c**2
      result *= self.U.rho_crit(z) * self.U.bg.Omega0_b/self.U.bg.Omega0_m
      return result
   
   def u(self, k, m, z):
      return self.rho3dFourier(k, m, z)


   def rho3dFourierNorm(self, k, m, z):
      """Fourier transform of 3d density profile,
      normalized such that rho3dFourierNorm(k=0, m, z) = 1,
      k in h Mpc^-1, m in h^-1 solarM;
      output dimensionless
      """
      # truncation radius
      rVir = self.U.frvir(m, z)
      m200c, r200c = self.U.massRadiusConversion(m, z, value=200, ref="c")
      xTrunc = self.trunc * rVir / r200c
      # do Fourier transform
      integrand = lambda x: 4.*np.pi*x**2 * self.rhoFit3d(x, m, z) * j0(k * r200c * x)
      result = integrate.quad(integrand, 0., xTrunc, epsabs=0., epsrel=1.e-2)[0]
      # compute normalization
      integrand = lambda x: 4.*np.pi*x**2 * self.rhoFit3d(x, m, z)
      result /= integrate.quad(integrand, 0., xTrunc, epsabs=0., epsrel=1.e-2)[0]
      return result
   

   def rho2d(self, R, m, z):
      """2d physical gas density profile
      in (Msun/h) / (Mpc/h)^2
      m is mVir in Msun/h
      R is comoving projected radius in Mpc/h
      """
      # convert mass from Mvir to M200c (both in [Msun/h])
      # and get comoving R200c [Mpc/h]
      m200c, r200c = self.U.massRadiusConversion(m, z, value=200, ref="c")
      X = R / r200c
      # projected fitting function from Battaglia 17
      result = self.rhoFit2d(X, m200c, z)
      # rescale by comoving critical density, in (Msun/h) / (Mpc/h)^3
      result *= self.U.rho_crit(z)
      # rescale by comoving radius, to get (Msun/h) / (Mpc/h)^2
      result *= r200c
      return result


   def ne3d(self, r, m, z):
      """3d physical electron number density profile
      in 1/(Mpc/h)^3
      assuming fully ionized gas and primordial He abundance
      m is mVir in Msun/h
      r is comoving radius in Mpc/h
      """
      result = self.rho3d(r, m, z)  # (Msun/h) / (Mpc/h)^3

      # convert from baryon mass to electron number
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      factor = (me + nH_ne*mH + nHe_ne*mHe)  # in kg
      msun = 1.989e30   # solar mass (kg)
      factor /= msun # in Msun
      factor *= self.U.bg.h   # in Msun/h
      
      # get the physical 3d electron number density
      result /= factor  # in (Mpc/h)^(-3)

      return result
   
   
   
   def ne2d(self, R, m, z):
      """2d physical electron number density profile
      in 1/(Mpc/h)^2
      assuming fully ionized gas and primordial He abundance
      m is mVir in Msun/h
      R is comoving projected radius in Mpc/h
      """
      result = self.rho2d(R, m, z)  # (Msun/h) / (Mpc/h)^2
   
      # convert from baryon mass to electron number
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      factor = (me + nH_ne*mH + nHe_ne*mHe)  # in kg
      msun = 1.989e30   # solar mass (kg)
      factor /= msun # in Msun
      factor *= self.U.bg.h   # in Msun/h
   
      # get the physical 2d electron number density
      result /= factor  # in (Mpc/h)^(-2)
      
      return result


   def tau3d(self, r, m, z):
      """Thompson scattering optical depth "3d profile"
      in 1/(Mpc/h) comoving
      ie you get the 2d tau profile by projecting this profile
      along the physical (not comoving) radial coordinate
      assuming fully ionized gas and primordial He abundance
      m is mVir in Msun/h
      r is comoving radius in Mpc/h
      """
      result = self.ne3d(r, m, z)
      
      # multiply by Thompson cross section (physical)
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T *= (self.U.bg.h/mpc)**2  # in (Mpc/h)^2
      result *= sigma_T # dimensionless
      
      return result


   def tau(self, R, m, z):
      """Thompson scattering optical depth profile
      dimensionless
      assuming fully ionized gas and primordial He abundance
      m is mVir in Msun/h
      R is comoving projected radius in Mpc/h
      """
      result = self.ne2d(R, m, z)
   
      # multiply by Thompson cross section (physical)
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T *= (self.U.bg.h/mpc)**2  # in (Mpc/h)^2
      
      result *= sigma_T # dimensionless
      return result
   
   
   def tauBeam(self, R, m, z, fwhm):
      """Thompson scattering optical depth profile,
      convolved with a Gaussian beam
      dimensionless
      assuming fully ionized gas and primordial He abundance
      m is mVir in Msun/h
      R is comoving projected radius in Mpc/h
      fwhm in arcmin
      """
      # comoving beam size at redshift z
      fwhm *= np.pi / (180.*60.) # radians
      sigmaBeam = fwhm / np.sqrt(8.*np.log(2.)) # radians
      sigmaBeam *= self.U.bg.comoving_distance(z) # in comoving Mpc/h
      # beam, function of comoving projected radius
      fbeam = lambda R: np.exp(-0.5*R**2/sigmaBeam**2) / (2.*np.pi*sigmaBeam**2)
      
      # do the smoothing
      f = lambda r: r * self.tau(r, m, z) * np.exp(-0.5*(r**2+R**2)/sigmaBeam**2) / sigmaBeam**2 * i0(r*R/sigmaBeam**2)
      result = integrate.quad(f, 0., np.inf, epsabs=0., epsrel=1.e-2)[0]
      print "convolved with beam"
      return result
   
   
   def tauBeamDiskRing(self, m, z, fwhm, R0):
      """convolve tau profile with beam,
      then apply equal area disk-ring filter
      output is dimensionless
      m is mVir in Msun/h
      fwhm in arcmin
      R0 comoving projected radius of disk in Mpc/h
      """
      f = lambda r: 2.*np.pi*r * self.tauBeam(r, m, z, fwhm)
      result = integrate.quad(f, 0., R0, epsabs=0, epsrel=1.e-2)[0]
      result -= integrate.quad(f, R0, np.sqrt(2.)*R0, epsabs=0, epsrel=1.e-2)[0]
      result /= np.pi*R0**2
      print "wahoo!"
      return result
   

   ##################################################################################

   def testFig5Battaglia16(self):
      # radii to evaluate
      X= np.logspace(np.log10(7.e-2), np.log10(4.), 101, 10.)
      
      # masses to evaluate
      M200c = np.array([7.5e13, 8.5e13, 9.5e13, 1.5e14, 2.5e14, 1.15e15])  # Msun
      M200c *= self.U.bg.h   # Msun/h
      
      # redshifts to evaluate
      Z = np.array([0., 0.5, 1.])


      # vary mass
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iM in range(len(M200c)):
         m = M200c[iM]
         z = 0.
         f = lambda x: self.rhoFit3d(x, m, z)
         Rho = np.array(map(f, X))
         ax.loglog(X, X**2 * Rho, label=r'$M_\text{200c}=$'+str(m)+r'$M_\odot$')
      #
      ax.legend(loc=3)
      #ax.set_ylim((1.e1, 1.e2))
      ax.set_xlabel(r'$x\equiv R/R_\text{200c}$')
      ax.set_ylabel(r'$x^2 \rho(x) / f_b \rho_\text{crit}(z)$')
      #
      #fig.savefig("/Users/Emmanuel/Desktop/fig5a_battaglia16.pdf", bbox_inches='tight')


      # vary redshift
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         m = 2.5e14
         z = Z[iZ]
         f = lambda x: self.rhoFit3d(x, m, z)
         Rho = np.array(map(f, X))
         ax.loglog(X, X**2 * Rho, label=r'$z=$'+str(z))
      #
      ax.legend(loc=3)
      #ax.set_ylim((1.e1, 1.e2))
      ax.set_xlabel(r'$x\equiv R/R_\text{200c}$')
      ax.set_ylabel(r'$x^2 \rho(x) / f_b \rho_\text{crit}(z)$')
      #
      #fig.savefig("/Users/Emmanuel/Desktop/fig5b_battaglia16.pdf", bbox_inches='tight')

      plt.show()


##################################################################################
##################################################################################

class ProfNFWPlusBattaglia16(Profile):
   """Matter density profile, with:
   DM following NFW
   baryons following Battaglia 16
   """
   
   def __init__(self, U):
      super(ProfNFWPlusBattaglia16, self).__init__(U)
      self.use_correction_factor = True
      self.mMin = 0.
      self.mMax = np.inf
   
      # instantiate the two profiles to combine
      self.ProfNFW = ProfNFW(U)
      self.ProfGasBattaglia16 = ProfGasBattaglia16(U)
   
      # baryon fraction
      self.fb = self.U.bg.Omega0_b / self.U.bg.Omega0_m
   
   
   def __str__(self):
      return "nfwplusbattagliadensity"

   def u(self, k, m, z):
      """returns m/\bar{rho} * NFW u, in (h^-1 Mpc)^3
      k in h Mpc^-1, m in h^-1 solarM
      """
      result = (1. - self.fb) * self.ProfNFW.nfw(k, m, z)
      result += self.fb * self.ProfGasBattaglia16.rho3dFourierNorm(k, m, z)
      result *= m / self.U.rho_m(z)
      return result


##################################################################################
##################################################################################

class ProfHODAlam16GasBattaglia16(Profile):
   """Mean gas profile for CMASS galaxies:
   gas from Battaglia+16
   CMASS HOD from Alam+16
   """
   
   def __init__(self, U, MassFunc, save=False):
      super(ProfHODAlam16GasBattaglia16, self).__init__(U)
      self.U = U
      self.MassFunc = MassFunc
      self.ProfGasBattaglia16 = ProfGasBattaglia16(U)
      self.ProfHODAlam16 = ProfHODAlam16(U, MassFunc)
      self.mMin = 0.
      self.mMax = np.inf
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      
      if save:
         self.saveAll()
      self.loadAll()
   
   def __str__(self):
      return "alam16battaglia16gas"

   def saveAll(self):
      # compute the convolution in Fourier space
      data = np.zeros((len(self.K), 2))
      data[:,0] = self.K
      data[:,1] = np.array(map(self.rho3dFourier, self.K))
      np.savetxt("./output/profile/profhodalam16gasbattaglia16_k.txt", data)
      self.rho3dFourierInterp = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
      
      # compute the real space convolution
      R = np.logspace(np.log10(0.001), np.log10(10.), 101, 10.)
      data = np.zeros((len(R), 2))
      data[:,0] = R
      f = lambda r: self.inverseFourier(self.rho3dFourierInterp, r)
      data[:,1] = np.array(map(f, R))
      np.savetxt("./output/profile/profhodalam16gasbattaglia16_r.txt", data)
   
   
   def loadAll(self):
      # interpolate the Fourier profile (after convolution)
      data = np.genfromtxt("./output/profile/profhodalam16gasbattaglia16_k.txt")
      forRho3dFourierInterp = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
      self.rho3dFourierInterp = lambda k, m, z: forRho3dFourierInterp(k)
   
      # interpolate the real space 2d profile
      data = np.genfromtxt("./output/profile/profhodalam16gasbattaglia16_r.txt")
      forRho2dInterp = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
      self.rho2dInterp = lambda r, m, z: forRho2dInterp(r)
   
   
   def u(self, k, m, z):
      result = self.ProfHODAlam16.u(k, m, z)
      result *= self.ProfGasBattaglia16.rho3dFourier(k, m, z)
      return result

   def rho3dFourier(self, k):
      z = 0.57
      a = 1./(1.+z)
      integrand = lambda lnm: np.exp(lnm)*self.MassFunc.fmassfunc(np.exp(lnm), a) * self.u(k, np.exp(lnm), z)
      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0, epsrel=1.e-2)[0]
      print k, result
      return result


   def plotHODSmoothing(self):
      """Compare the gas profile of CMASS galaxies in these scenarios:
      - CMASS galaxies are all central with Mhalo=2.e13, z=0.57
      - CMASS galaxies obey an HOD, at z=0.57
      """
      
      R = np.logspace(np.log10(0.001), np.log10(10.), 101, 10.)

      # simple CMASS gas profile

      # more sophisticated CMASS gas profile



##################################################################################
##################################################################################
##################################################################################






##################################################################################
##################################################################################

class ProfLIMEGG(Profile):

   def __init__(self, U, Sfr, lineName='halpha', trunc=4., unit='dI/I'):
      '''Choice of unit for intensity:
      unit=='dI/I' for delta I / I [dimless]
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      super(ProfLIMEGG, self).__init__(U)
      self.Sfr = Sfr

      # choose to compute dI/I or I
      self.unit = unit

      # NFW profile setup
      self.LoadNonLinMass()
      self.trunc = trunc   # truncation radius, in units of rVir
      self.use_correction_factor = False#False
      self.mMin = 0.
      self.mMax = np.inf

      # Line setup: read from EGG
      self.lineName = lineName
      # read line parameters from EGG
      pathIn = './input/EGG_results/'
      # read redshift and mean number density of galaxies [(Mpc/h)^{-3}]
      d1 = np.loadtxt(pathIn+'ngal.txt')
      z = d1[:,0]
      nZ = len(z)
      self.nGal = interp1d(z, d1[:,1], kind='linear', bounds_error=False, fill_value=0.)
      # find the corresponding line number
      lineNames = np.loadtxt(pathIn+'line_names.txt', dtype='str')
      nLines = len(lineNames)
      iLine = np.where(lineNames==self.lineName)[0][0]
      # find the line wavelength and frequency
      self.lambdaMicrons = np.loadtxt(pathIn+'line_lambda_microns.txt')[iLine] # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]
      # read mean galaxy luminosity [Lsun]
      d2 = np.loadtxt(pathIn+'mean_gal_lum.txt')
      self.meanGalLum = interp1d(z, d2[:,iLine], kind='linear', bounds_error=False, fill_value=0.)
      # read total galaxy luminosity density # [Lsun / (Mpc/h)^3]
      d3 = np.loadtxt(pathIn+'total_gal_lum_density.txt')
      self.meanLumDensity = interp1d(z, d3[:,iLine], kind='linear', bounds_error=False, fill_value=0.)

      # read fractional covariance of galaxy line luminosities [dimless]
      d4 = np.loadtxt(pathIn+'s2ij.txt')
      d4 = d4.reshape((nZ, nLines, nLines))
      self.s2ij = interp1d(z, d4, axis=0, kind='linear', bounds_error=False, fill_value=0.)      
   
   def __str__(self):
      return "limegg"+self.lineName
   
   def LoadNonLinMass(self):
      """precompute the nonlinear mass at z=0.
      """
      print "Loading non-lin mass at z=0"
      z = 0.
      self.m_nonlin = self.U.nonLinMass(z)


   def rS_rhoS_c(self, m, z):
      """comoving scale radius for NFW profile
      in Mpc/h
      """
      Rvir = self.U.frvir(m, z)
      # concentration parameter
      #c = 10./(1.+z) * (m / self.m_nonlin)**(-0.2)   # from Takada & Jain 2002
      c = 9./(1.+z) * (m / self.m_nonlin)**(-0.13) # Takada & Jain 2003
      # scale radius
      RS = Rvir / c  # in Mpc/h
      # normalize the mass within rVir to be mVir
      rhoS = m / (4.*np.pi*RS**3)
      rhoS /= np.log(1.+c) - c/(1.+c)  # (Msun/h) / (Mpc/h)^3
      return RS, rhoS, c
   
   
   def totalMass(self, trunc=None):
      """total mass within truncation radius
      trunc in units of rVir
      mass in Msun/h
      if trunc=infinity, the mass is infinite
      """
      if trunc is None:
         trunc = self.trunc
      rVir = self.U.frvir(m, z)
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      # truncation radius over scale radius
      xMax = trunc * rVir/rS
      result = 4./3. * np.pi * rS**3 * rhoS
      result = xMax - np.log(1 + xMax)
      return result
   

   def nfw(self, k, m, z):
      """Fourier transform of NFW density profile,
      normalized such that u(k=0, m, z) = 1
      ie rhoNFW(k,m,z) = m * nfw(k,m,z)
      truncation radius is taken to be infinite (unphysical)
      k in h Mpc^-1
      m in h^-1 solarM
      """
      RS, rhoS, c = self.rS_rhoS_c(m, z)
      #
      result = np.sin(k * RS) * (  Si((1+c) * k * RS) - Si(k * RS)  )
      result += - np.sin(c * k * RS) / ((1+c) * k * RS)
      result += np.cos(k * RS) * (  Ci((1+c) * k * RS) - Ci(k * RS)  )
      result /= (np.log(1+c) - c/(1+c))
      return result



   def Ngal(self, m, z):
      '''Mean number of galaxies in one halo [dimless]
      of mass m [Msun/h] at redshift z.
      Assumed propto SFR(m)
      '''
      result = self.nGal(z)  # [(Mpc/h)^{-3}]
      result *= self.Sfr.sfr(m, z)   # [Msun/yr]
      result /= self.Sfr.sfrd(z) # [(Msun/yr) (Mpc/h)^{-3}]
#!!!!! test
#      result *= m / self.U.rho_m(z)
      return result  # [dimless]


   def meanHaloLum(self, m, z):
      '''Mean luminosity of one halo [Lsun],
      of mass m [Msun/h] at redshift z.
      '''
      return self.meanGalLum(z) * self.Ngal(m, z)

   
   def u(self, k, m, z, unit=None):
      '''Effective profile for the halo model integrals
      Unit is [(Mpc/h)^3], multiplied by the intensity unit:
      unit=='dI/I' for delta I / I
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      if unit is None:
         unit = self.unit
      result = self.nfw(k, m, z)
      result *= self.meanHaloLum(m, z)
      result /= self.meanLumDensity(z)
      result *= self.meanIntensity(z, unit=unit)
      #result *= np.exp(- 0.5 * sigma**2 * k**2 * mu**2)
      if not np.isfinite(result):
         result = 0.
      return result
   

   def meanIntensity(self, z, unit=None):
      '''Mean intensity in 
      unit=='dI/I' for delta I / I [dimless]
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      if unit is None:
         unit = self.unit

      if unit=='dI/I':
         return 1.

      result = self.meanLumDensity(z)  # [Lsun / (Mpc/h)^3]
      result *= 3.e5 / self.U.hubble(z)   # *[Mpc/h]
      result /= 4. * np.pi * self.nuHz # *[/sr/Hz]
      if unit=='Lsun/(Mpc/h)^2/sr/Hz':
         return result
      if unit=='Jy/sr':
         result *= 3.827e26   # [Lsun] to [W]
         result /= (3.086e22 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [m^{-2}]
         result /= 1.e-26  # [W/m^2/Hz/sr] to [Jy/sr]
      elif unit=='cgs':
         result *= 3.839e33  # [Lsun] to [erg/s]
         result /= (3.086e24 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [cm^{-2}]
      return result


   def nGalEff(self, z, lineName1=None, lineName2=None):
      '''Effective galaxy number density for shot noise
      [(Mpc/h)^{-3}]
      '''
      if lineName1 is None:
         lineName1 = self.lineName
      if lineName2 is None:
         lineName2 = self.lineName
      # find the corresponding line numbers
      pathIn = './input/EGG_results/'
      lineNames = np.loadtxt(pathIn+'line_names.txt', dtype='str')
      iLine1 = np.where(lineNames==lineName1)[0][0]
      iLine2 = np.where(lineNames==lineName2)[0][0]
      # get the interpolated fractional covariance
      # of galaxy line luminosities [dimless]
      s2 = self.s2ij(z)[iLine1, iLine2]

      return self.nGal(z) / (1. + s2)

   def Pshot(self, z, lineName1=None, lineName2=None, unit=None):
      '''shot noise power spectrum in [(Mpc/h)^3] multiplied
      by the square of the intensity unit
      '''
      if unit is None:
         unit = self.unit
      if lineName1 is None:
         lineName1 = self.lineName
      if lineName2 is None:
         lineName2 = self.lineName

      result = 1. / self.nGalEff(z, lineName1, lineName2)

      if unit=='dI/I':
         return result
      else:
#!!!! manuwaring: this is super super slow. Fix it
         if lineName1==self.lineName:
            p1 = self
         else:
            p1 = ProfLIM(self.U, self.Sfr, lineName=lineName1, trunc=self.trunc, unit=unit)
         if lineName2==self.lineName:
            p2 = self
         else:
            p2 = ProfLIM(self.U, self.Sfr, lineName=lineName2, trunc=self.trunc, unit=unit)
         result *= p1.meanIntensity(z, unit=unit)
         result *= p2.meanIntensity(z, unit=unit)
         return result
      


   ####################################################

   def plotNgal(self):
      Z = np.linspace(0.,5.,6)
      M = np.logspace(np.log10(1.e10), np.log10(1.e16), 101, 10.) # masses in h^-1 solarM


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         f = lambda m: self.Ngal(m, z)
         Ngal = np.array(map(f, M))
         ax.semilogx(M, Ngal, label=r'$z=$'+str(round(z,2)))
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'halo mass $M$  [$M_\odot$/h]')
      ax.set_ylabel(r'$N_\text{gal}(M, z)$')

      plt.show()


   def plotnGal(self):
      Z = np.linspace(0.71, 6.,101)

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      nGal = np.array(map(self.nGal, Z))
      ax.semilogy(Z, nGal)
      #
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{n}_\text{gal}$ [(Mpc/h)$^{-3}$]')

      plt.show()


   def plotMeanIntensity(self):
      Z = np.linspace(0.71, 6.,101)

      # reproduce Fig 3 in Gong Cooray Silva + 17
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # compare to Gong+17
      if self.lineName=='halpha':
         # center values
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_center.csv"
         data = np.genfromtxt(path, delimiter=', ')
         plt.plot(data[:,0], data[:,1], 'b', label=r'Hopkins Beacom 06')
         # error band
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_low.csv"
         low = np.genfromtxt(path, delimiter=', ')
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_high.csv"
         high = np.genfromtxt(path, delimiter=', ')
         plt.fill(np.append(low[:,0], high[::-1,0]), np.append(low[:,1], high[::-1,1]), facecolor='b', alpha=0.5)
      #
      # results from my calculation (EGG)
      f = lambda z: self.meanIntensity(z, unit='Jy/sr')
      meanIntensity = np.array(map(f, Z))
      ax.semilogy(Z, meanIntensity, 'k', label=r'from EGG')
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{I}(z)$ [Jy/sr]')

      plt.show()

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



   def plotP3dGong17(self, fp3d, lineName=None):
      '''Compares power spectrum to Gong+17 fig5.
      Power spectrum function should be in [(Jy/sr)^3(Mpc/h)^3]
      '''
      if lineName is None:
         lineName = self.lineName
      
      zNames = np.array(['1.0', '1.4', '1.8', '2.2', '2.7', '3.3', '4.0', '4.8'])
      Z = zNames.astype(float)
      nZ = len(Z)
      colors = plt.cm.cool(np.arange(nZ)/(nZ-1.))


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(nZ):
         zName = zNames[iZ]
         z = Z[iZ]

         # read Gong+17 fig 5 power spectrum
         path = "./input/LIM_literature/Gong+17/fig5/Gong+17_fig5_"+lineName+"_z_"+zName+".csv"
         data = np.genfromtxt(path)
         k = data[:,0]  # h/Mpc
         p = data[:,1]  # k^3 P(k) / (2pi^2) [(Jy/sr)^2]
         #ax.loglog(k, p, c=colors[iZ], ls='-', label=zName)
         p /= k**3 / (2. * np.pi**2) # convert to [(Jy/sr)^2 (Mpc/h)^3]
#         ax.loglog(k, p, c=colors[iZ], ls='-', label=zName)

         # convert to [(Mpc/h)^3]
         # by dividing out the factor of intensity squared 
         # from Gong+17 fig3, model by Hopkins Beacom 06
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_center.csv"
         data = np.genfromtxt(path, delimiter=', ')
         fMeanIntensityGong17 = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
         p /= fMeanIntensityGong17(z)**2
         


         # my calculation
         f = lambda k: fp3d(k, z)
         myp = np.array(map(f, k))
         #ax.loglog(k, myp, c=colors[iZ], ls='--', label=zName)
         ax.loglog(k, myp/p, c=colors[iZ], ls='--', label=zName)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$P(k, z)$ [(Jy/sr)$^2$(Mpc/h)$^3$]')

      plt.show()




##################################################################################
##################################################################################


class ProfLIMSobral12(Profile):
   '''Halpha luminosity functions from Sobral+12 arXiv:1202.3436v2
   '''

   def __init__(self, U, Sfr, trunc=4., unit='dI/I'):
      '''Choice of unit for intensity:
      unit=='dI/I' for delta I / I [dimless]
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      super(ProfLIMSobral12, self).__init__(U)
      self.Sfr = Sfr

      # choose to compute dI/I or I
      self.unit = unit

      # NFW profile setup
      self.LoadNonLinMass()
      self.trunc = trunc   # truncation radius, in units of rVir
      self.use_correction_factor = False#False
      self.mMin = 0.
      self.mMax = np.inf

      # Line properties
      self.lineName = 'halpha'
      self.lambdaMicrons = 656.28e-3   # [mu]
      self.nuHz = 299792458. / self.lambdaMicrons * 1.e6 # [Hz]
      

      # Measured luminosity functions
      # table 4
      self.Z = [0.4, 0.84, 1.47, 2.23]
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
      

      # Schechter fits to the intrinsic luminosity functions
      # table 5
      path = './input/LIM_literature/Sobral+12/table5/Sobral+12_table5.txt'
      data = np.genfromtxt(path)
      z = data[:,0]
      #
      LStar = 10.**data[:,1]  # L^star_Halpha [erg/s]
      LStarHigh = 10.**(data[:,1]+data[:,2])  # L^star_Halpha [erg/s]
      LStarLow = 10.**(data[:,1]+data[:,3])  # L^star_Halpha [erg/s]
      #
      PhiStar = 10.**data[:,4] / self.U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarHigh = 10.**(data[:,4]+data[:,5]) / self.U.bg.h**3  # [(Mpc/h)^-3]
      PhiStarLow = 10.**(data[:,4]+data[:,6]) / self.U.bg.h**3  # [(Mpc/h)^-3]
      #
      # for alpha, they recommend using -1.6, rather than the best fit
      #Alpha = -1.6 * np.ones_like(z) # [dimless]
      Alpha = data[:,7] # [dimless]
      AlphaHigh = data[:,7]+data[:,8] # [dimless]
      AlphaLow = data[:,7]+data[:,9] # [dimless]
      

      # interpolate the Schechter fits (Eq 2)
      self.lStar = interp1d(z, LStar, kind='linear', bounds_error=False, fill_value=0.)
      self.lStarHigh = interp1d(z, LStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      self.lStarLow = interp1d(z, LStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.phiStar = interp1d(z, PhiStar, kind='linear', bounds_error=False, fill_value=0.)
      self.phiStarHigh = interp1d(z, PhiStarHigh, kind='linear', bounds_error=False, fill_value=0.)
      self.phiStarLow = interp1d(z, PhiStarLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.alpha = interp1d(z, Alpha, kind='linear', bounds_error=False, fill_value=0.)
      self.alphaHigh = interp1d(z, AlphaHigh, kind='linear', bounds_error=False, fill_value=0.)
      self.alphaLow = interp1d(z, AlphaLow, kind='linear', bounds_error=False, fill_value=0.)
      #
      self.phiInt = lambda z,l: self.phiStar(z) * (l/self.lStar(z))**self.alpha(z) * np.exp(-l/self.lStar(z)) / self.lStar(z) # [(Mpc/h)^-3 (erg/s)^-1]
      self.phiObs = lambda z,l: 100.**(1./5.) * self.phiInt(z, l * 100.**(1./5.))

      '''
      self.nGal(z)
      self.meanGalLum(z)
      self.meanLumDensity(z)
      self.s2ij(z)
      '''

   def phi(self, z, l, obs=True, lStar=None, phiStar=None, alpha=None):
      '''output in [(Mpc/h)^-3 (erg/s)^-1]
      luminosity in [erg/s]
      obs = True for dust extincted luminosities
      obs = False for intrinsic luminosities
      '''
      if obs:
         lStarCorr = 100.**(1./5.)
      else:
         lStarCorr = 1.

      if lStar is None:
         lStar = self.lStar(z) / lStarCorr
      elif lStar is 'high':
         lStar = self.lStarHigh(z) / lStarCorr
      elif lStar is 'low':
         lStar = self.lStarLow(z) / lStarCorr


      if phiStar is None:
         phiStar = self.phiStar(z)
      elif phiStar is 'high':
         phiStar = self.phiStarHigh(z)
      elif phiStar is 'low':
         phiStar = self.phiStarLow(z)

      if alpha is None:
         alpha = self.alpha(z)
      elif alpha is 'high':
         alpha = self.alphaHigh(z)
      elif alpha is 'low':
         alpha = self.alphaLow(z)

      result = phiStar * (l/lStar)**alpha * np.exp(-l/lStar) / lStar

      return result



   def nGal(self, z, lStar=None, phiStar=None, alpha=None):
      ''' [(Mpc/h)^-3]
      This integral of the LF formally diverges at low luminosities,
      so it is highly cutoff dependent.
      '''
      def integrand(lnl):
         l = np.exp(lnl)
         result = self.phi(z, l, obs=True, lStar=lStar, phiStar=phiStar, alpha=alpha)
         result *= l
         return result
      #result = integrate.quad(integrand, np.log(1.e30), np.log(1.e44), epsabs=0., epsrel=1.e-3)[0]
      result = integrate.quad(integrand, np.log(1.e30), np.log(1.e44), epsabs=0., epsrel=1.e-3)[0]
      return result


   def meanLumDensity(self, z, lStar=None, phiStar=None, alpha=None):
      '''[Lsun / (Mpc/h)^3]
      '''
      def integrand(lnl):
         l = np.exp(lnl)
         result = self.phi(z, l, obs=True, lStar=lStar, phiStar=phiStar, alpha=alpha)
         result *= l**2
         return result
      result = integrate.quad(integrand, np.log(1.e30), np.log(1.e44), epsabs=0., epsrel=1.e-3)[0]
      # convert from [erg/s] to [Lsun]
      result /= 3.839e33
      return result
      

   def meanGalLum(self, z, lStar=None, phiStar=None, alpha=None):
      '''[Lsun]
      '''
      result = self.meanLumDensity(z, lStar=lStar, phiStar=phiStar, alpha=alpha)
      result /= self.nGal(z, lStar=lStar, phiStar=phiStar, alpha=alpha)
      return result


   def s2ij(self, z, lStar=None, phiStar=None, alpha=None):
      '''var(L_Ha) / mean(L_Ha)^2 [dimless]
      '''
      def integrand(lnl):
         l = np.exp(lnl)
         result = self.phi(z, l, obs=True, lStar=lStar, phiStar=phiStar, alpha=alpha)
         result *= l**3
         return result
      result = integrate.quad(integrand, np.log(1.e30), np.log(1.e44), epsabs=0., epsrel=1.e-3)[0]
      result /= self.nGal(z, lStar=lStar, phiStar=phiStar, alpha=alpha)
      # convert from [(erg/s)^2] to [Lsun^2]
      result /= 3.839e33**2
      result /= self.meanGalLum(z, lStar=lStar, phiStar=phiStar, alpha=alpha)**2
      result -= 1.
      return result
   



   def nGalEff(self, z, lStar=None, phiStar=None, alpha=None):
      '''Effective galaxy number density for shot noise
      [(Mpc/h)^{-3}]
      '''
      nGal = self.nGal(z, lStar=lStar, phiStar=phiStar, alpha=alpha)
      s2 = self.s2ij(z, lStar=lStar, phiStar=phiStar, alpha=alpha)
      return nGal / (1. + s2)


   def meanIntensity(self, z, unit=None):
      '''Mean intensity in
      unit=='dI/I' for delta I / I [dimless]
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      if unit is None:
         unit = self.unit

      if unit=='dI/I':
         return 1.

      result = self.meanLumDensity(z)  # [Lsun / (Mpc/h)^3]
      result *= 3.e5 / self.U.hubble(z)   # *[Mpc/h]
      result /= 4. * np.pi * self.nuHz # *[/sr/Hz]
      if unit=='Lsun/(Mpc/h)^2/sr/Hz':
         return result
      if unit=='Jy/sr':
         result *= 3.827e26   # [Lsun] to [W]
         result /= (3.086e22 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [m^{-2}]
         result /= 1.e-26  # [W/m^2/Hz/sr] to [Jy/sr]
      elif unit=='cgs':
         result *= 3.839e33  # [Lsun] to [erg/s]
         result /= (3.086e24 / self.U.bg.h)**2  # [(Mpc/h)^{-2}] to [cm^{-2}]
      return result


   def Pshot(self, z, unit=None, lStar=None, phiStar=None, alpha=None):
      '''shot noise power spectrum in [(Mpc/h)^3] multiplied
      by the square of the intensity unit
      '''
      if unit is None:
         unit = self.unit

      result = 1. / self.nGalEff(z, lStar=lStar, phiStar=phiStar, alpha=alpha)

      if unit=='dI/I':
         return result
      else:
         result *= self.meanIntensity(z, unit=unit)
         result *= self.meanIntensity(z, unit=unit)
         return result


   def plotShotNoiseUncertainty(self, compareProfs=None):
      '''Vary the Schechter fit parameters to get an idea of the uncertainty 
      on the shot noise power spectrum
      '''

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      colors = ['gray', 'g', 'r', 'b']
      for iZ in range(self.nZ):
         z = self.Z[iZ]
         y = [self.Pshot(z, lStar=None, phiStar=None, alpha=None, unit='dI/I'),
              self.Pshot(z, lStar='high', phiStar=None, alpha=None, unit='dI/I'),
              self.Pshot(z, lStar='low', phiStar=None, alpha=None, unit='dI/I'),
              self.Pshot(z, lStar=None, phiStar='high', alpha=None, unit='dI/I'),
              self.Pshot(z, lStar=None, phiStar='low', alpha=None, unit='dI/I'),
              self.Pshot(z, lStar=None, phiStar=None, alpha='high', unit='dI/I'),
              self.Pshot(z, lStar=None, phiStar=None, alpha='low', unit='dI/I')]
         #ax.axhline(y[0], xmin=1.*iZ/self.nZ, xmax=(iZ+1.)/self.nZ, color=colors[iZ], label=r'$z=$'+str(z) )
         #ax.axhspan(np.min(y), np.max(y), xmin=1.*iZ/self.nZ, xmax=(iZ+1.)/self.nZ, color=colors[iZ], alpha=0.3)
         
         ax.errorbar([z], [y[0]], yerr=[0.5*(np.max(y) - np.min(y))], c='b', fmt='o-')
      ax.errorbar([], [], c='b', label=r'Sobral+12')



      #
      for prof in compareProfs:
         if hasattr(prof, 'Z'):
            Z = prof.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         f = lambda z: prof.Pshot(z, unit='dI/I')
         pShot = np.array(map(f, Z))
         plt.plot(Z, pShot, label=str(prof))
         
      #
      ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
      #ax.axes.xaxis.set_visible(False)
      #ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_ylabel(r'$P_\text{shot}$ [(Mpc/h)$^3$]')
      #ax.set_xlabel(r'$k$ [h/Mpc]')
      #
      path = './figures/profile/Sobral12/'+'shot_noise_uncertainty.pdf'
      fig.savefig(path, bbox_inches='tight')
      plt.show()
































   def plotLF(self):
      '''Reproduces fig 8 in Sobral+12.
      '''

      L = np.logspace(np.log10(1.e40), np.log10(1.e44), 501, 10.)

      
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
         ax.plot(L, self.phiInt(z, L) * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ])
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
      path = './figures/profile/Sobral12/'+'lf_intrinsic_sobral12_fig8.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      plt.show()
      

      
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
         ax.plot(L, self.phiObs(z, L) * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ])
         #
         # my high and low curves
         #y = self.phiObs(z, L) + self.fPhiObsHighMeas[iZ](L) - self.fPhiObsMeas[iZ](L)
         #ax.plot(L, y * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         #y = self.phiObs(z, L) + self.fPhiObsLowMeas[iZ](L) - self.fPhiObsMeas[iZ](L)
         #ax.plot(L, y * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         #
         # varying the Schechter fits within 1 sigma
         ax.plot(L, self.phi(z, L, obs=True, lStar='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         ax.plot(L, self.phi(z, L, obs=True, lStar='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         ax.plot(L, self.phi(z, L, obs=True, phiStar='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         ax.plot(L, self.phi(z, L, obs=True, phiStar='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         ax.plot(L, self.phi(z, L, obs=True, alpha='high') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
         ax.plot(L, self.phi(z, L, obs=True, alpha='low') * L * np.log(10.) * self.U.bg.h**3, c=colors[iZ], ls='--')
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
      path = './figures/profile/Sobral12/'+'lf_observed.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      plt.show()
      
      
















   def __str__(self):
      return "limsobral12"+self.lineName

   def LoadNonLinMass(self):
      """precompute the nonlinear mass at z=0.
      """
      print "Loading non-lin mass at z=0"
      z = 0.
      self.m_nonlin = self.U.nonLinMass(z)


   def rS_rhoS_c(self, m, z):
      """comoving scale radius for NFW profile
      in Mpc/h
      """
      Rvir = self.U.frvir(m, z)
      # concentration parameter
      #c = 10./(1.+z) * (m / self.m_nonlin)**(-0.2)   # from Takada & Jain 2002
      c = 9./(1.+z) * (m / self.m_nonlin)**(-0.13) # Takada & Jain 2003
      # scale radius
      RS = Rvir / c  # in Mpc/h
      # normalize the mass within rVir to be mVir
      rhoS = m / (4.*np.pi*RS**3)
      rhoS /= np.log(1.+c) - c/(1.+c)  # (Msun/h) / (Mpc/h)^3
      return RS, rhoS, c


   def totalMass(self, trunc=None):
      """total mass within truncation radius
      trunc in units of rVir
      mass in Msun/h
      if trunc=infinity, the mass is infinite
      """
      if trunc is None:
         trunc = self.trunc
      rVir = self.U.frvir(m, z)
      rS, rhoS, c = self.rS_rhoS_c(m, z)
      # truncation radius over scale radius
      xMax = trunc * rVir/rS
      result = 4./3. * np.pi * rS**3 * rhoS
      result = xMax - np.log(1 + xMax)
      return result


   def nfw(self, k, m, z):
      """Fourier transform of NFW density profile,
      normalized such that u(k=0, m, z) = 1
      ie rhoNFW(k,m,z) = m * nfw(k,m,z)
      truncation radius is taken to be infinite (unphysical)
      k in h Mpc^-1
      m in h^-1 solarM
      """
      RS, rhoS, c = self.rS_rhoS_c(m, z)
      #
      result = np.sin(k * RS) * (  Si((1+c) * k * RS) - Si(k * RS)  )
      result += - np.sin(c * k * RS) / ((1+c) * k * RS)
      result += np.cos(k * RS) * (  Ci((1+c) * k * RS) - Ci(k * RS)  )
      result /= (np.log(1+c) - c/(1+c))
      return result



   def Ngal(self, m, z):
      '''Mean number of galaxies in one halo [dimless]
      of mass m [Msun/h] at redshift z.
      Assumed propto SFR(m)
      '''
      result = self.nGal(z)  # [(Mpc/h)^{-3}]
      result *= self.Sfr.sfr(m, z)   # [Msun/yr]
      result /= self.Sfr.sfrd(z) # [(Msun/yr) (Mpc/h)^{-3}]
#!!!!! test
#      result *= m / self.U.rho_m(z)
      return result  # [dimless]


   def meanHaloLum(self, m, z):
      '''Mean luminosity of one halo [Lsun],
      of mass m [Msun/h] at redshift z.
      '''
      return self.meanGalLum(z) * self.Ngal(m, z)


   def u(self, k, m, z, unit=None):
      '''Effective profile for the halo model integrals
      Unit is [(Mpc/h)^3], multiplied by the intensity unit:
      unit=='dI/I' for delta I / I
      unit=='Lsun/(Mpc/h)^2/sr/Hz' for I
      unit=='Jy/sr' for I
      unit=='cgs' for I [erg/s/cm^2/sr/Hz]
      '''
      if unit is None:
         unit = self.unit
      result = self.nfw(k, m, z)
      result *= self.meanHaloLum(m, z)
      result /= self.meanLumDensity(z)
      result *= self.meanIntensity(z, unit=unit)
      #result *= np.exp(- 0.5 * sigma**2 * k**2 * mu**2)
      if not np.isfinite(result):
         result = 0.
      return result




   ####################################################

   def plotNgal(self):
      Z = np.linspace(0.,5.,6)
      M = np.logspace(np.log10(1.e10), np.log10(1.e16), 101, 10.) # masses in h^-1 solarM


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         f = lambda m: self.Ngal(m, z)
         Ngal = np.array(map(f, M))
         ax.semilogx(M, Ngal, label=r'$z=$'+str(round(z,2)))
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'halo mass $M$  [$M_\odot$/h]')
      ax.set_ylabel(r'$N_\text{gal}(M, z)$')

      plt.show()


   def plotnGal(self, profs=None):
      if profs is None:
         profs = [self]
      
      print("nGal is highly cutoff dependent (formally divergent at low luminosity)")
      print("so its value is meaningless")
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for prof in profs:
         if hasattr(prof, 'Z'):
            Z = prof.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         nGal = np.array(map(prof.nGal, Z))
         ax.semilogy(Z, nGal, label=str(prof))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{n}_\text{gal}$ [(Mpc/h)$^{-3}$]')

      plt.show()


   def plotMeanIntensity(self, profs=None):
      if profs is None:
         profs = [self]

      # reproduce Fig 3 in Gong Cooray Silva + 17
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # compare to Gong+17
      if self.lineName=='halpha':
         # center values
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_center.csv"
         data = np.genfromtxt(path, delimiter=', ')
         plt.plot(data[:,0], data[:,1], 'b', label=r'Hopkins Beacom 06')
         # error band
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_low.csv"
         low = np.genfromtxt(path, delimiter=', ')
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_high.csv"
         high = np.genfromtxt(path, delimiter=', ')
         plt.fill(np.append(low[:,0], high[::-1,0]), np.append(low[:,1], high[::-1,1]), facecolor='b', alpha=0.5)
      #
      for prof in profs:
         f = lambda z: prof.meanIntensity(z, unit='Jy/sr')
         if hasattr(prof, 'Z'):
            Z = prof.Z
         else:
            Z = np.linspace(0.71, 6.,101)
         meanIntensity = np.array(map(f, Z))
         ax.semilogy(Z, meanIntensity, label=str(prof))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\bar{I}(z)$ [Jy/sr]')

      plt.show()

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




   def plotP3dGong17(self, fp3d, lineName=None):
      '''Compares power spectrum to Gong+17 fig5.
      Power spectrum function should be in [(Jy/sr)^3(Mpc/h)^3]
      '''
      if lineName is None:
         lineName = self.lineName

      zNames = np.array(['1.0', '1.4', '1.8', '2.2', '2.7', '3.3', '4.0', '4.8'])
      Z = zNames.astype(float)
      nZ = len(Z)
      colors = plt.cm.cool(np.arange(nZ)/(nZ-1.))


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(nZ):
         zName = zNames[iZ]
         z = Z[iZ]

         # read Gong+17 fig 5 power spectrum
         path = "./input/LIM_literature/Gong+17/fig5/Gong+17_fig5_"+lineName+"_z_"+zName+".csv"
         data = np.genfromtxt(path)
         k = data[:,0]  # h/Mpc
         p = data[:,1]  # k^3 P(k) / (2pi^2) [(Jy/sr)^2]
         #ax.loglog(k, p, c=colors[iZ], ls='-', label=zName)
         p /= k**3 / (2. * np.pi**2) # convert to [(Jy/sr)^2 (Mpc/h)^3]
#         ax.loglog(k, p, c=colors[iZ], ls='-', label=zName)

         # convert to [(Mpc/h)^3]
         # by dividing out the factor of intensity squared
         # from Gong+17 fig3, model by Hopkins Beacom 06
         path = "./input/LIM_literature/Gong+17/fig3/mean_intensity_"+self.lineName+"_hopkinsbeacom06_center.csv"
         data = np.genfromtxt(path, delimiter=', ')
         fMeanIntensityGong17 = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
         p /= fMeanIntensityGong17(z)**2



         # my calculation
         f = lambda k: fp3d(k, z)
         myp = np.array(map(f, k))
         #ax.loglog(k, myp, c=colors[iZ], ls='--', label=zName)
         ax.loglog(k, myp/p, c=colors[iZ], ls='--', label=zName)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$P(k, z)$ [(Jy/sr)$^2$(Mpc/h)$^3$]')

      plt.show()


