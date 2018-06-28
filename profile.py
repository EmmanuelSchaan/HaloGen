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
         ax.loglog(K, U, lw=2, label=r'$M=$'+str(m)+'$M_\odot/h$')
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

   def __init__(self, U, trunc=1.):
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
      return self.nfw(k, m, z) * m / self.U.rho_z(z)
   

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
         #f = lambda k: self.u(k, m, z) / (m/self.U.rho_z(z))
         f = lambda k: self.nfw(k, m, z)
         Y = np.array(map(f, K))
         ax.loglog(K, Y, 'b')
      #
      ax.set_xlim((min(K), max(K)))
      ax.set_ylim((1.e-3, 2.5))
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
      result *= m / self.U.rho_z(z)
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
#      self.mMinHod = 10.**11.5 * self.U.h # in Msun/h
#      self.sLogM = 0.65#np.sqrt(0.65)
#      self.mSat = 10. * self.mMinHod # in Msun/h
#      self.aSat = 1.4
      
      # best fit to Planck CIB. Table 1 in Penin+14
      if self.nu==217:
         self.mMinHod = 10.**12. * self.U.h # in Msun/h
         self.aSat = 1.5
      if self.nu==353:
         self.mMinHod = 10.**12.2 * self.U.h # in Msun/h
         self.aSat = 1.7
      if self.nu==545:
         self.mMinHod = 10.**12.6 * self.U.h # in Msun/h
         self.aSat = 1.9
      if self.nu==857:
         self.mMinHod = 10.**12.8 * self.U.h # in Msun/h
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
#      Rvir = ( 3.*m / (4*np.pi*self.U.rhocrit_z(z) * self.U.Deltacrit_z(z)) )**(1./3.)
#      Rs = Rvir / cNFW
#
#      # NFW scale density (physical)
#      rhoS = m / (4.*np.pi*Rs**3) / (np.log(1.+cNFW) - cNFW/(1.+cNFW))
#
#      # get R200 and M200
#      f = lambda x: -1. + 1./(1.+x) + np.log(1.+x) - 200./3.*(self.U.rhocrit_z(z)/rhoS)*x**3
#      x = optimize.brentq(f , 0.1, 100.)
#      R200 = x * Rs  # physical
#      M200 = 4./3.*np.pi*R200**3 * self.U.rhocrit_z(z) * 200.
#      R200 *= (1.+z) # comoving
#
#      return M200, R200

   
   # dimless normalized thermal3 pressure
   # x=r/r_200c is dimless, m=m_200c in h^-1 solarM
   def P_PDelta_x(self, x, m, z):
      # parameters of the fit, from Battaglia et al 2012
      P0 = self.P00 * (m / (1.e14 * self.U.h))**self.P0am * (1.+z)**self.P0az
      xc = self.xc0 * (m / (1.e14 * self.U.h))**self.xcam * (1.+z)**self.xcaz
      beta = self.beta0 * (m / (1.e14 * self.U.h))**self.betaam * (1.+z)**self.betaaz
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
      P200 = self.U.G * (M200*msun)/self.U.h
      P200 *= 200.
      P200 *= self.U.rhocrit_z(z) * self.U.h**2 * msun/(1.e6*pc)**3 # critical density in kg m^-3 (comoving)
      P200 *= (1.+z)**3 # physical
      P200 *= self.U.OmB/self.U.OmM # baryon fraction
      P200 /= 2. * R200/(1.+z)* 1.e6*pc/self.U.h  # twice R200 in m (physical)
      P200 *= 1.e6*pc/self.U.h
      
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
      da = self.U.ComovDist(1./(1.+z), 1.)
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















