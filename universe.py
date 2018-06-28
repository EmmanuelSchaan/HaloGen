from headers import *
##################################################################################


class Universe(object):

   def __init__(self):
      """Required variables: h, OmC, OmB, OmM, path_lin_matter_power, a_obs
      """
      
      # physical constants
      self.G = 6.67e-11   # in SI
      self.rhoCRIT = 2.7617e11 # critical density today, in (h^-1 solarM) (h^-1 Mpc)^-3
      
      # background universe
      self.rho0 = self.OmM * self.rhoCRIT # matter density today, in (h^-1 solarM) (h^-1 Mpc)^-3

      # linear power spectrum from CAMB
      readPlin0 = np.genfromtxt(self.path_lin_matterpower)
      K0 = readPlin0[:, 0]  # comoving k in units h * Mpc^-1
      Plin0 = readPlin0[:,1]  # linear power spectrum in (h^-1 Mpc)^3, at redshift 0
      # number of k points to compute sigma2 (1500 gives the same sigma8 as CAMB)
      self.Nk = 1500
      # fit to extend linear power spectrum to higher k than CAMB
      self.k_eq = 0.0155 # turning point of the linear Plin in (h Mpc^-1)
      self.k_fit_min = 5.e-4 # in (h Mpc^-1), fit in k
      self.k_fit_max = 10. # in (h Mpc^-1), fit in k^-3 ln(k)^2
      # interpolate k, and create k array to compute sigma2
      fK = UnivariateSpline(np.linspace(0., 1., len(K0)),K0,k=1,s=0)
      self.K = fK(np.linspace(0., 1., self.Nk)) # new k array
      # interpolate linear power spectrum
      self.fPlin0 = UnivariateSpline(K0,Plin0,k=1,s=0) # function for power spectrum
      self.Plin = self.fPlin0(self.K) # new array for power spectrum

      # read non-linear power spectrum (from CAMB with halofit)
      #self.path_nonlin_matterpower = "./input/nonlin_matterpower_HillPajer13.dat"
      #readPnonlin = np.genfromtxt(self.path_nonlin_matterpower)
      #self.Knonlin = readPnonlin[:, 0]
      #self.Pnonlin = readPnonlin[:, 1]
      
      # fitting function for the correlation coefficient
      # between the reconstructed and true velocity
      # from Mariana's simulations
      self.fr = lambda k: min(0.92, 0.54/k**0.165)
      #self.fr = lambda k: min(1., 0.54/k**0.165)


   ##################################################################################
   # background universe
   ##################################################################################

   def rho_z(self, z):
      """Comoving matter density at redshift z, in (h^-1 solarM) (h^-1 Mpc)^-3
      """
      return self.rho0


   def rhocrit_z(self, z):
      """Comoving critical density at redshift z, in (h^-1 solarM) (h^-1 Mpc)^-3
      """
      # critical density per unit physical volume
      result = self.rhoCRIT * ( self.OmM * (1.+z)**3 + (1. - self.OmM) )
      # critical density per unit comoving volume
      result /= (1.+z)**3
      return result
   
   
   def Hubble(self, a):
      """a dimless, Hubble in (km s^-1 (h Mpc^-1))
      """
      return 100 * ( self.OmM/a**3 + (1. - self.OmM) )**0.5
   
   
   def ComovDist(self, a_min, a_max):
      """Comoving distance along light cone from a=a_min to a=a_max, a dimless, ComovDist(a) in (h^-1 Mpc)
      """
      f = lambda a: 1./ (self.Hubble(a)/3.e5 * a**2)
      result = integrate.quad(f, a_min, a_max, epsabs=0, epsrel=1.e-5)[0]
      return result
   
   
   def Lookback(self, a_min, a_max):
      """Lookback time in Gyr/h
      """
      Mpc_to_km = 3.08567758e19
      f = lambda a: 1./ (self.Hubble(a) * a)
      result = integrate.quad(f, a_min, a_max, epsabs=0, epsrel=1.e-5)[0]
      result *= Mpc_to_km / (3600.*24.*365.*1.e9)
      return result
   
   
   def LinGrowth(self, a):
      """a dimless, LinGrowth(a) dimless, normalize to 1 for a=1
      """
      f = lambda a: 1./(self.Hubble(a) * a)**3
      result = integrate.quad(f, 0., a, epsabs=0, epsrel=1.e-5)[0]
      result /= integrate.quad(f, 0., 1., epsabs=0, epsrel=1.e-5)[0] # normalize to 1 for a=1
      result *= self.Hubble(a)/100.
      return result
   

   ##################################################################################
   # spherical collapse
   ##################################################################################


   def deltaC_z(self, z):
      """critical density for spherical collapse at redshift z
      from Henry 2000, from Nakamura & Suto 1997
      usual 3.*(12.*pi)**(2./3.) / 20. = 1.686 if OmM=1.
      """
      x = ( 1./self.OmM - 1. )**(1./3.)
      x /= 1.+z
      dc = 3.*(12.*np.pi)**(2./3.) / 20.
      dc *= 1. - 0.0123* np.log( 1. + x**3 )
      return dc
   

   '''
   def DeltaVir(self, z):
      """Overdensity wrt mean for virialized halo at z
      from Bullock et al 2001, from Bryan & Norman 1998
      usual 18*pi**2 if OmM=1.
      gives 337 at z=0 for OmM= 0.3
      Omega = rhocrit(z)/rho_matter(z)
      """
      Omega = self.OmM*(1.+z)**3
      Omega /= self.OmM*(1.+z)**3 + (1. - self.OmM)
      x = Omega - 1.
      Dvir = 18*np.pi**2 + 82.*x - 39.*x**2
      # convert between rho_m(z) and rho_crit(z)
      Dvir /= Omega
      return Dvir
   '''
   '''
   def DeltaVir(self, z):
      """Overdensity wrt mean for virialized halo at z
      from Henry 2000, from Nakamura & Suto 1997
      usual 18*pi**2 if OmM=1.
      """
      x = ( 1./self.OmM - 1. )**(1./3.)
      x /= 1.+z
      Dvir = 18*np.pi**2 * ( 1. + 0.4093* x**2.71572 )
      return Dvir
   '''


   def Deltacrit_z(self, z):
      """ratio of virialized density to critical density at collapse (dimless).
      from Bullock et al 2001, from Bryan & Norman 1998
      usual 18*pi**2 if OmM=1.
      Omega = rhocrit(z)/rho_matter(z).
      """
      f = self.OmM * (1.+z)**3 / ( self.OmM * (1.+z)**3 + (1. - self.OmM) )
      return 18.*np.pi**2 + 82.*(f-1.) - 39.*(f-1.)**2
   
   
   def frvir(self, m, z):
      """Comoving virial and scale radii (Mpc/h)
      input mass is mvir (Msun/h)
      """
      Rvir = ( 3.*m / (4*np.pi*self.rhocrit_z(z)*self.Deltacrit_z(z)) )**(1./3.)  # in h^-1 Mpc
      return Rvir

   ##################################################################################
   # density fluctuations
   ##################################################################################


   def Plin_z(self, z):
      """array for linear power spectrum at redshift z
      """
      return self.Plin * self.LinGrowth(1./(1.+z))**2


   def fPlin_z(self, k, z):
      """Linear power spectrum, with a fit to extrapolate to high k
      output is (h^-1 Mpc)^3
      """
      if (k < self.k_fit_min):
         result = self.fPlin0(self.k_fit_min) * (k/self.k_fit_min)
      elif (k > self.k_fit_max):
         f = k**(-3) * np.log(k/(8.*self.k_eq))**2
         result = f * self.fPlin0(self.k_fit_max) / ( self.k_fit_max**(-3) * np.log(self.k_fit_max/(8.*self.k_eq))**2 )
      else:
         result = self.fPlin0(k)
      growth = self.LinGrowth(1./(1.+z))**2
      result *= growth
      return result
   
   
   def fdPlindK_z(self, k, z):
      """derivative of linear power spectrum, with a fit to extrapolate to high k
      output is (Mpc/h)^4
      """
      if (k < self.k_fit_min):
         result = self.fPlin0(self.k_fit_min) / self.k_fit_min
      elif (k > self.k_fit_max):
         f = np.log(k/self.k_eq)/k**4 * (2. - 3.*np.log(k/self.k_eq))
         result = f * self.fPlin0(self.k_fit_max) / ( self.k_fit_max**(-3) * np.log(self.k_fit_max/self.k_eq)**2 )
      else:
         result = self.fPlin0.derivatives(k)[1]
      growth = self.LinGrowth(1./(1.+z))**2
      result *= growth
      return result


   def Sigma2(self, R, z, W3d):
      """variance of delta on an isotropic 3d domain,
      defined by W3d
      R in h^-1 Mpc, comoving scale, output is dimless
      """
      F = (self.K**3) * self.Plin_z(z) * ( np.array(map(W3d, self.K * R))**2 ) / (2* np.pi**2)  # dimensionless
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      return Itrap
   
   
   def nonLinMass(self, z):
      """nonlin mass at z, in Msun/h
      from Takada and Jain 2002/2003
      """
      # bounds for looking for m_nonlinin, in (h^-1 solarM) (h^-1 Mpc)^-3
      ma = 1.e11
      mb = 1.e13
      # solve for sigma(m)^2 = deltaC**2
      R = lambda m: ( 3.* m / (4*np.pi*self.rho_z(z)) )**(1./3.)   # in h^-1 Mpc
      f = lambda m: self.Sigma2(R(m), z, W3d_sth) - self.deltaC_z(z)**2
      # find mass such that nu(m, z) = 1
      result = optimize.brentq(f , ma, mb)
      return result



   def Sigma2_2d(self, r, z, W2d):
      """this is Dchi * sigma^2,
      where sigma^2 = <( delta averaged on bin Dchi )^2>
      r is comoving scale in h^-1 Mpc
      """
      F = self.K * self.Plin_z(z) * ( np.array(map(W2d, self.K * r))**2 ) / (2* np.pi)
      dK = self.K[1:]-self.K[0:-1]
      Itrap = np.sum( dK * ( F[:-1] + F[1:] ) ) * 0.5
      return Itrap



   def Sigma2_cone(self, aMin, aMax, theta, W2d):
      """Variance of delta on a section of a cone,
      with fixed angular radius theta
      from aMin to aMax
      output dimless
      """
      z = lambda a: 1./a-1.
      r = lambda a: self.ComovDist(a, 1.) * theta
      integrand = lambda a: 3.e5/(self.Hubble(a) * a**2) * self.Sigma2_2d(r(a), z(a), W2d)
   
      result = integrate.quad(integrand, aMin, aMax, epsabs=0., epsrel=1.e-3)[0]
      chiMin = self.ComovDist(aMax, 1.)
      chiMax = self.ComovDist(aMin, 1.)
      result /= (chiMax-chiMin)**2
      return result



   def dlnSigma2_dlnR(self, R, z):
      """R in h^-1 Mpc, comoving scale, output is dimless
      dln(sigma2) / dln(R)
      """
      F = (self.K**3) * self.Plin_z(z) * 2.*np.array(map(W3d_sth, self.K*R))*np.array(map(dW3d_sth, self.K * R))*(self.K * R) / (2* np.pi**2)  # dimensionless
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      result = Itrap / self.Sigma2(R, z, W3d_sth)
      return result


   def fdlnSigma_dlnM(self, m, z):
      """dln(sigma)/dln(m)
      """
      R = (3.*m / (4.*np.pi*self.rho_z(z)))**(1./3.)
      result = self.dlnSigma2_dlnR(R, z) /6.
      return result



   def fnu(self, m, z):
      """nu = dc**2/sigma2(m, z)
      """
      r = (3.*m / (4.*np.pi*self.rho_z(z)))**(1./3.)
      s2 = self.Sigma2(r, z, W3d_sth)
      nu = self.deltaC_z(z)**2 / s2
      return nu


   def fdlnnu_dlnm(self, m, z):
      """dln(nu)/dln(m)
      """
      return -2.* self.fdlnSigma_dlnM(m, z)


   ##################################################################################
   # halo mass conversion


   def massRadiusConversion(self, m, z, value=200, ref="m"):
      """converts virial mass (Msun/h)
      into M_value,ref (Msun/h),
      and corresponding comoving R_value,ref (Mpc/h).
      Assumes concentration from Duffy et al 2008.
      Used by Tinker mass function
      """
      
      # concentration params from Duffy et al 2008, used by Tinker
      cNFW0 = 5.71
      cNFWam = -0.084
      cNFWaz = -0.47
      
      # from Duffy et al 2008: different pivot mass
      cNFW = cNFW0 * (m/2.e12)**cNFWam * (1.+z)**cNFWaz
      # comoving virial radius and scale radius in h^-1 Mpc
      Rvir = ( 3.*m / (4*np.pi*self.rhocrit_z(z) * self.Deltacrit_z(z)) )**(1./3.)
      Rs = Rvir / cNFW
      # NFW scale density (comoving)
      rhoS = m / (4.*np.pi*Rs**3) / (np.log(1.+cNFW) - cNFW/(1.+cNFW))
      
      # comoving reference density
      if ref=="m":   # ie wrt mean density
         rhoRef = self.rho_z(z)
      elif ref=="c": # ie wrt critical density
         rhoRef = self.rhocrit_z(z)

      # get R200 and M200
      f = lambda x: -1. + 1./(1.+x) + np.log(1.+x) - value/3.*(rhoRef/rhoS)*x**3
      x = optimize.brentq(f , 0.1, 100.)
      Rnew = x * Rs
      Mnew = 4./3.*np.pi*Rnew**3 * rhoRef * value
      
      return Mnew, Rnew


   ##################################################################################
   # response of power spectrum to local overdensity
   
   def fdlnPdDelta(self, k, z):
      result = 68./21.
      result -= 1.
      result -= k/self.fPlin_z(k, z) * self.fdPlindK_z(k, z) / 3.
      return result
   
   
   
   def testdLnPlindDelta(self, z=0.):
      """for the linear power spectrum, dlnP/ddelta is independent of z
      but dP/ddelta is not
      """
      K = np.logspace(np.log10(1.e-5), np.log10(1.e2), 1001, 10.)
      #
      f = lambda k: self.fPlin_z(k, z)
      P = np.array(map(f, K))
      #
      f = lambda k: self.fdlnPdDelta(k, z)
      dP = np.array(map(f, K))
      
      # dlnP/ddelta
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.semilogx(K, dP, 'b', label=r'$\frac{68}{21} - \frac{1}{3}\frac{d\ln k^3P_\text{lin}}{d\ln k}$')
      ax.axhline(68./21., color='r', label=r'$\frac{68}{21}$')
      #
      ax.legend(loc=4)
      ax.grid()
      #ax.set_xlim((0., 0.4))
      #ax.set_ylim((2.2, 3.2))
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$d\ln P_\text{lin}/d\delta$ [(Mpc/h)$^3$]')
      #fig.savefig('./figures/response_powerspectrum/dLnPlindDelta.pdf', bbox_inches='tight')
      
      # dP/ddelta
      fig = plt.figure(1)
      ax = plt.subplot(111)
      #
      ax.loglog(K, P*dP, 'b', label=r'$P(k)\left[\frac{68}{21} - \frac{1}{3}\frac{d\ln k^3P_\text{lin}}{d\ln k} \right]$')
      ax.loglog(K, P*68./21., 'r', label=r'$P(k)\left[\frac{68}{21}\right]$')
      ax.loglog(K, P, 'k', label=r'$P(k)$')
      #
      ax.legend(loc=3)
      ax.grid()
      #ax.set_xlim((0., 0.4))
      #ax.set_ylim((2.2, 3.2))
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$dP_\text{lin}/d\delta$ [(Mpc/h)$^3$]')
      #fig.savefig('./figures/response_powerspectrum/dPlindDelta_loglog.pdf', bbox_inches='tight')
      
      plt.show()

   
   ##################################################################################
   # Velocity fluctuations
   
   def RMSVelocity(self, R, z, W3d):
      """R in h^-1 Mpc, comoving scale, output is the rms velocity in (km/s)
      """
      F = (self.K**3) * ( np.array(map(W3d, self.K * R))**2 ) / (2* np.pi**2)  # dimensionless
      F *= self.Plin_z(z) / self.K**2
      F *= ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
      F *= ( self.Hubble(1./(1.+z))/(1.+z) )**2
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      return np.sqrt(Itrap)


   def RMSErrorRecVel(self, R, z, W3d):
      """R in h^-1 Mpc, comoving scale, output is the rms error on reconstructed velocity in (km/s)
      """
      CorrCoeff = np.array(map(self.fr, self.K))
      F = (self.K**3) * ( np.array(map(W3d, self.K * R))**2 ) / (2* np.pi**2)  # dimensionless
      F *= self.Plin_z(z) / self.K**2
      F *= 2.* (1. - CorrCoeff)
      F *= ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
      F *= ( self.Hubble(1./(1.+z))/(1.+z) )**2
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      return np.sqrt(Itrap)

   def growthLogDerivativeF(self, z):
      return ( self.OmM*(1.+z)**3 / (self.OmM*(1.+z)**3+(1. - self.OmM)) )**(5./9.)

   def convertVelToDisp(self, z):
      """multiply a velocity at redshift z by this factor to get the displacement at that redshift
      v in km/s, displacement in Mpc/h
      """
      factor = self.Hubble(z)/(1.+z)*self.growthLogDerivativeF(z)
      return 1./factor
   
   def RMSDisplacement(self, R, z, W3d):
      """R in h^-1 Mpc, comoving scale, output is the rms displacement in (Mpc/h)
      """
      result = self.RMSVelocity(R, z, W3d)
      result *= self.convertVelToDisp(z)
      return result
   
   # variance of vRadial_RMS^2, the average vRadial^2 on a cylindrical volume
   # with area "area" in deg2, and between redshifts zMin and zMax
#   def varRMSVRadial(self, area, zMin, zMax):


#   # variance of vRadial_RMS^2, the average vRadial^2 on a spherical volume
#   # with radius R in comoving Mpc/h
#   def varRMSVRadialSpherical(self, R, z, W3d):
#      
#      def integrand(pars):
#         k = pars[0]
#         K = pars[1]
#         mu = pars[2]
#         x = K/k
#         # Put the density instead of velocity!!!
#         result = self.fPlin_z(k, z) / k**2
#         result *= self.fPlin_z(k*np.sqrt(1 - 2.*x*mu + x**2), z) / (k*np.sqrt(1 - 2.*x*mu + x**2))**2
#         
#         factor = ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
#         factor *= ( self.Hubble(1./(1.+z))/(1.+z) )**2 # (a*H)**2
#         
#         result *= factor**2
#         result *= W3d(k*x*R)**2
#         result *= x**2 * k**5
#         result *= k
#         result *= 4./ 3. / (2.*np.pi)**4
#         return result
#   
#      integ = vegas.Integrator([[1.e-4, 1.e1], [1.e-4, 1.e1], [-1., 1.]])
#      integ(integrand, nitn=4, neval=1000)
#      result = integ(integrand, nitn=8, neval=1000)
#      print result.sdev/result.mean
#      return result.mean
##      return result.summary()


   def varRMSVSpherical(self, R, z, W3d):
      """variance of v_RMS^2, the average v^2 on a spherical volume
      with radius R in comoving Mpc/h
      """

      def integrand(pars):
         k = pars[0]
         K = pars[1]
         mu = pars[2]
#         x = K/k
         # Put the density instead of velocity!!!
         result = self.fPlin_z(k, z) / k**2
         result *= self.fPlin_z(np.sqrt(k**2 - 2.*k*K*mu + K**2), z) / (k**2 - 2.*k*K*mu + K**2)
         
         factor = ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
         factor *= ( self.Hubble(1./(1.+z))/(1.+z) )**2 # (a*H)**2
         
         result *= factor**2
         result *= W3d(K*R)**2
         result *= K**2 * k**2
         result *= 4. / (2.*np.pi)**4
         return result
      
      integ = vegas.Integrator([[1.e-4, 1.e1], [1.e-4, 1.e1], [-1., 1.]])
      integ(integrand, nitn=8, neval=1000)
      result = integ(integrand, nitn=8, neval=1000)
      print result.sdev/result.mean
#      print result.summary()
      return result.mean
   

   def plotVarRMSVSpherical(self, z=0.):
      vMin = 1.e-2 * 1.e9  # (Mpc/h)^3
      vMax = 1.e1 * 1.e9
      V = np.logspace(np.log10(vMin), np.log10(vMax), 11, 10.)
      R = ( 3.*V/(4.*np.pi) )**(1./3.)
      
      # We want the v_RMS^2, not the square of the average v,
      # this is why the scale is 0 Mpc/h
      f = lambda r: self.RMSVelocity(r*0., z, W3d_sth)**2
      v2 = np.array(map(f, R))
      # We want the cosmic variance on v_RMS^2;
      # now this quantity will depend on the scale r
      f = lambda r: self.varRMSVSpherical(r, z, W3d_sth)
      varV2 = np.array(map(f, R))
      
      # volume of D56
      volumeD56 = 700.*(np.pi/180.)**2 # area in sr
      volumeD56 *= self.ComovDist(1./(1.+0.57), 1.)**2 # area in (Mpc/h)^2
      volumeD56 *= self.ComovDist(1./(1.+0.7), 1./(1.+0.4)) # volume in (Mpc/h)^3
      # volume of D56+BOSS N
      volumeD56BN = volumeD56 * 2700./700.
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(V / 1.e9, np.sqrt(varV2)/v2, 'b-')
      #
      ax.axvline(volumeD56 / 1.e9, c='k', linestyle='--')
      ax.axvline(volumeD56BN / 1.e9, c='k', linestyle='--')
      #
      ax.set_xscale('log')
      ax.set_xlabel(r'Volume [Gpc/h]$^3$')
      ax.set_ylabel(r'$\sigma_{v_\text{RMS}^2} / v_\text{RMS}^2$')
      #
#      fig.savefig("./figures/velocities/variance_of_vrms2_spherical.pdf", bbox_inches='tight')

      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.plot(V / 1.e9, np.sqrt(v2), 'b-')
      ax.plot(V / 1.e9, np.sqrt(np.sqrt(varV2)), 'b-')
      #
      ax.set_xscale('log')
      ax.set_xlabel(r'Volume [Gpc/h]$^3$')
      ax.set_ylabel(r'$v_\text{RMS}$')


      plt.show()
      return


   ##################################################################################
   # show properties
   ##################################################################################
   
   
   def testLinearPowerSpectrum(self):
      
      # show linear power spectrum; test interpolation
      K = np.logspace(np.log10(1.e-5), np.log10(1.e2), 10001, 10.)
      z=0.
      f = lambda k: self.fPlin_z(k, z)
      Plin = np.array(map(f, K))
      fig0 = plt.figure(0)
      ax = plt.subplot(111)
      ax.loglog(self.K, self.Plin, 'k.')
      ax.loglog(K, Plin, 'r')
      ax.grid()
#      ax.set_xlim((1.e-3, 1.e0))
#      ax.set_ylim((1.e1, 1.e5))
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$P(k)$ [(Mpc/h)$^3$]')
      #fig0.savefig('./figures/Plin.pdf')
      
      # compute variance of matter overdensity as a function of scale
      z = 0.
      fR  = lambda m: (3.*m / (4.*np.pi*self.rho_z(z)))**(1./3.)
      fS2 = lambda m: self.Sigma2(fR(m), z, W3d_sth)
      M = np.logspace(np.log10(1.e10), np.log10(1.e16), 101, 10.) # masses in h^-1 solarM
      R = np.array(map(fR, M))
      S2 = np.array(map( fS2, M ))
      # Sigma2 = f(m)
      fig1 = plt.figure(1)
      ax = plt.subplot(111)
      ax.loglog(M, S2, 'r')
      ax.loglog(M, M/M, 'k')
      ax.grid()
      ax.set_xlabel(r'mass $m$ [$M_{sun}/h$]')
      ax.set_ylabel(r'variance of $\delta_m$, smoothed on a scale $m$')
      #fig1.savefig('./figures/S2_m.pdf')
      # Sigma2 = f(r)
      fig2 = plt.figure(2)
      ax = plt.subplot(111)
      ax.loglog(R, S2, 'b')
      ax.loglog(R, R/R, 'k')
      ax.grid()
      ax.set_xlabel(r'scale $R$ [Mpc/h]')
      ax.set_ylabel(r'variance of $\delta_m$, smoothed on a scale $R$')
      #fig2.savefig('./figures/S2_r.pdf')
      
      plt.show()
   
   
   def testRMSErrorRecVel(self):
      # input
      R = 0.   # smoothing scale
      z = 0.   # redshift
      W3d = W3d_sth  # choice of window function

      # correlation coefficient
      CorrCoeff = np.array(map(self.fr, self.K))
      
      # velocity power spectrum, normalized to its maximum
      F = self.Plin_z(z) / self.K**2
      #F *= 2.* (1. - CorrCoeff)
      F *= ( self.OmM*(1.+z)**3 )**(2.*5./9.)   # f**2, where f = Omega_matter**5/9
      F *= self.LinGrowth(1./(1.+z))**2
      F *= ( self.Hubble(1./(1.+z))/(1.+z) )**2
      F /= np.max(F)
      
      # integrand for RMS velocity, normalized to its maximum
      G = F * (self.K**3) * ( np.array(map(W3d, self.K * R))**2 ) / (2* np.pi**2)  # dimensionless
      G /= np.max(G)
      
      # integrand for RMS velocity error, normalized to its maximum
      H = G * 2.* (1. - CorrCoeff)
   
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.loglog(self.K, F, 'b', label=r'$P_{vv}$')
      ax.loglog(self.K, G, 'g', label=r'$k^3 P_{vv}$')
      ax.loglog(self.K, (1-CorrCoeff), 'r', label=r'$(1-r)$')
      ax.loglog(self.K, H, 'm', label=r'$(1-r(k))k^3 P_{vv}$')
      #
      ax.grid()
      ax.legend(loc=4)
      ax.set_ylim((1.e-5, 1.))
      ax.set_xlabel(r'$k$ h/Mpc')
      #fig.savefig("./figures/velocities/rms_error.pdf")
   
      plt.show()
  

   ##################################################################################
  
   
   def plotDistances(self):
      
      Z = np.linspace(0., 10., 501)
      # comoving distance, i.e. comoving angular diameter distance (flat universe)
      fcomov = lambda z: self.ComovDist(1./(1.+z), 1.)
      Comov = np.array(map( fcomov, Z ))
      # lookback time
      ftime = lambda z: self.Lookback(1./(1.+z), 1.)
      Time = np.array(map( ftime, Z ))
      # proper distance, i.e. physical angular diameter distance (flat universe)
      Angular = Comov / (1.+Z)
      # luminosity distance, and distance modulus
      Lumi = Comov * (1.+Z)
      DistMod = 5.*( np.log10(Lumi*1.e6*self.h) - 1. )
      # linear growth factor
      fLinGrowth = lambda z: self.LinGrowth(1./(1.+z))
      LinGrowth = np.array(map( fLinGrowth, Z ))
      
      plt.figure(0)
      ax=plt.subplot(111)
      ax.plot(Z, Comov)
      ax.grid()
      ax.set_title('Comoving distance, i.e.comov. ang. diam. dist. (flat universe)')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('distance [Mpc/h]')
      
      plt.figure(1)
      ax=plt.subplot(111)
      ax.plot(Z, Angular)
      ax.grid()
      ax.set_title('Physical ang. diam. dist. (flat universe)')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('distance [Mpc/h]')
      
      plt.figure(2)
      ax=plt.subplot(111)
      ax.plot(Z, Lumi)
      ax.grid()
      ax.set_title('Luminosity distance')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('distance [Mpc/h]')
      
      plt.figure(3)
      ax=plt.subplot(111)
      ax.plot(Z, DistMod)
      ax.grid()
      ax.set_title(r'Distance modulus $\mu=5 (log_{10}(d/pc) - 1)$')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('distance modulus')
      
      plt.figure(4)
      ax=plt.subplot(111)
      ax.plot(Z, Time)
      ax.grid()
      ax.set_title(r'Lookback time')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('time [Gyr/h]')
      
      plt.figure(5)
      ax=plt.subplot(111)
      ax.plot(Z, LinGrowth)
      ax.grid()
      ax.set_title(r'Linear growth factor')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('linear growth factor')
      
      plt.show()


   ##################################################################################

   def plotHSV_2d_scaling(self):
      # range of redshift
      Z = np.linspace(0.005, 4., 8.)
      Z = np.array([0.005, 0.01, 0.02, 0.1, 0.2, 1.,  2., 3., 4.])
      A = 1./(1.+Z)
      # values of fsky to compute
      Fsky = np.logspace(np.log10(1.e-2), np.log10(1.), 51, 10.)

      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      for a in A:
         fscaling = lambda fsky: self.Sigma2_2d(self.ComovDist(a, self.a_obs) * 2.*np.sqrt(fsky), 1./a-1., W2d_cth) * fsky
         HSV = np.array(map(fscaling, Fsky))
         ax.loglog(Fsky, HSV, label=r'$z=$'+str(round(1./a-1., 2)))
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$f_{sky}$', fontsize=18)
      ax.set_ylabel(r'$\sigma^2(z, f_{sky})  f_{sky}$', fontsize=18)
      #fig.savefig('./figures/hsv_2d_scaling.pdf')

      plt.show()
      
      
      
   def plotHSV_3d_scaling(self):
      # survey volumes
      Vmin = 1.e-3 * (1.e3)**3
      Vmax = 10. * (1.e3)**3
      NV = 51
      V = np.logspace(np.log10(Vmin), np.log10(Vmax), NV, 10.)
      
      # variance of the matter overdensity smoothed over the survey (only used for HSV for 3d cov)
      R = (3.*V/(4.*np.pi))**(1./3.)   # for a spherical survey
      HSV = np.array( map(lambda r: self.Sigma2(r, 1./self.a_obs-1., W3d_sth), R) ) * V
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.loglog(V/(1.e3)**3, HSV)
      #ax.loglog(V/(1.e3)**3, V**(-1./3.)* 1.e7)
      #ax.loglog(V/(1.e3)**3, V)
      #
      ax.set_xlabel(r'survey volume [(Gpc/h)$^3$]', fontsize=18)
      ax.set_ylabel(r'$\sigma^2 V$', fontsize=18)
      fig.savefig('./figures/hsv_3d_scaling.pdf')
      
      plt.show()


   ##################################################################################
   
   
   def plotRMSVel(self, z=0.):
      Z = np.linspace(0., 1., 101)
      
      f = lambda z: self.RMSVelocity(0., z, W3d_sth)
      vRMS = np.array(map(f, Z))
   
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, vRMS, 'b', lw=2)
      #
      ax.axvspan(0.1, 0.3, alpha=0.3, color='g', label=r'maxBCG')
      ax.axvline(0.2, color='g')
      ax.axvspan(0.08, 0.55, alpha=0.3, color='r', label=r'redMaPPer')
      ax.axvline(0.35, color='r')
      ax.axvspan(0.4, 0.7, alpha=0.3, color='c', label=r'CMASS')
      ax.axvline(0.57, color='c')
      #
      ax.grid()
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\sqrt{ \langle v^2 \rangle }$ [km/s]')
      #fig.savefig("./figures/velocities/vrms_z.pdf")
   
      plt.show()


   ##################################################################################

   def plotSizeACTDeep(self):
      # typical size of the ACT Deep footprint
      theta = 10.*(np.pi/180.)
      Z = np.linspace(0., 1., 21)
      
      f = lambda z: self.ComovDist(1./(1.+z), 1.) * theta
      Size = np.array(map(f, Z))
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.axvspan(0.1, 0.3, alpha=0.3, color='g', label=r'maxBCG')
      ax.axvline(0.2, color='g')
      ax.axvspan(0.08, 0.55, alpha=0.3, color='r', label=r'redMaPPer')
      ax.axvline(0.35, color='r')
      ax.axvspan(0.4, 0.7, alpha=0.3, color='c', label=r'CMASS')
      ax.axvline(0.57, color='c')
      #
      ax.plot(Z, Size, 'k', lw=2)
      #
      ax.legend(loc=4)
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'Comoving size [Mpc/h]')
      #fig.savefig("./figures/velocities/transsize_act_deep.pdf")
      
      fig=plt.figure(1)
      ax=plt.subplot(111)
      #
      ax.axhspan(self.ComovDist(1./(1.+0.1), 1.), self.ComovDist(1./(1.+0.3), 1.), alpha=0.3, color='g', label=r'maxBCG')
      ax.axhline(self.ComovDist(1./(1.+0.2), 1.), color='g')
      ax.axhspan(self.ComovDist(1./(1.+0.08), 1.), self.ComovDist(1./(1.+0.55), 1.), alpha=0.3, color='r', label=r'redMaPPer')
      ax.axhline(self.ComovDist(1./(1.+0.35), 1.), color='r')
      ax.axhspan(self.ComovDist(1./(1.+0.4), 1.), self.ComovDist(1./(1.+0.7), 1.), alpha=0.3, color='c', label=r'CMASS')
      ax.axhline(self.ComovDist(1./(1.+0.57), 1.), color='c')
      #
      ax.set_ylim((0., 2000.))
      ax.legend(loc=4)
      ax.set_ylabel(r'Comoving size [Mpc/h]')
      #fig.savefig("./figures/velocities/longsize_act_deep.pdf")
      
      plt.show()


   ##################################################################################


   def testVelCorrCoeff(self):
      """check my fitting function for the correlation coefficient
      between the true and reconstructed velocities
      """
   
      # read correlation coefficient between the reconstructed and true velocity
      # from Mariana's simulations
      data = np.genfromtxt("./input/plot_mariana/r.txt")

      fr = lambda k: min(0.92, 0.54/k**0.165)
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.semilogx(data[:,0], data[:,1], 'ko', label=r'from Mariana')
      ax.semilogx(self.K, 0.92*np.ones_like(self.K), 'r', label=r'$0.92$')
      ax.semilogx(self.K, 0.54/self.K**0.165, 'r', label=r'$0.54 / k^{0.165}$')
      ax.semilogx(self.K, np.array(map(fr, self.K)), 'b', lw=2, label=r'fitting function')
      #
      ax.grid()
      ax.legend(loc=1)
      ax.set_xlim((7.e-3, 2.))
      ax.set_ylim((0.5, 1.))
      ax.set_xlabel(r'$k$ Mpc/h')
      ax.set_ylabel(r'correlation coefficient')
      ax.set_title(r'smoothing on $5$Mpc/h')
      #fig.savefig("./figures/velocities/corr_coeff.pdf")

      plt.show()



   def fvelocityCorrelation(self, R, z):
      """R in h^-1 Mpc, comoving scale, output is the rms velocity in (km/s)
      """
      f = lambda k: special.jv(0,k*R)
      F = np.array(map(f, self.K))
      
      F *= (self.K**3) / (2* np.pi**2)
      F *= self.Plin_z(z) / self.K**2
      F *= ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
      F *= ( self.Hubble(1./(1.+z))/(1.+z) )**2
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      return Itrap


   def fdensityCorrelation(self, R, z):
      """R in h^-1 Mpc, comoving scale, output is in (km/s)^2
      """
      f = lambda k: special.jv(0,k*R)
      F = np.array(map(f, self.K))
      
      F *= (self.K**3) / (2* np.pi**2)
      F *= self.Plin_z(z)
      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
      return Itrap


   def plotCorrelationFunction(self):
      z = 0.
      R = np.logspace(np.log10(1.), np.log10(1.e3), 101, 10.)
      
      f = lambda r: self.fvelocityCorrelation(r, z)
      VelCorr = np.array(map(f, R))
      
      f = lambda r: self.fdensityCorrelation(r, z)
      DenCorr = np.array(map(f, R))
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.semilogx(R, VelCorr, 'b')
      
      fig=plt.figure(1)
      ax=plt.subplot(111)
      #
      ax.semilogx(R, DenCorr, 'b')
      
      plt.show()




##################################################################################

class UnivHillPajer13(Universe):
   
   def __init__(self):
      
      # redshift of observer
      self.a_obs = 1.
      
      # background universe
      self.h = 0.697  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmC = 0.236
      self.OmB = 0.0461
      self.OmM = self.OmC + self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/lin_matterpower_HillPajer13.dat"

      super(UnivHillPajer13, self).__init__()


##################################################################################

class UnivTinkerEtAl08(Universe):
   
   def __init__(self):
   
      # redshift of observer
      self.a_obs = 1.
 
      # background universe
      self.h = 0.7  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmC = 0.26
      self.OmB = 0.04
      self.OmM = self.OmC + self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/lin_matterpower_wmap1_TinkerEtAl08.dat"
      
      super(UnivTinkerEtAl08, self).__init__()


##################################################################################

class UnivSchaanEtAl14(Universe):
   
   def __init__(self):
   
      # redshift of observer
      self.a_obs = 1.
      
      # background universe
      self.h = 0.732  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmM = 0.238
      self.OmB = 0.042
      self.OmC = self.OmM - self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/lin_matterpower_SchaanEtAl14.dat"
      
      super(UnivSchaanEtAl14, self).__init__()


##################################################################################

# uses WMAP9 + eCMB parameters
class UnivHandEtAl13(Universe):
   
   def __init__(self):
      
      # redshift of observer
      self.a_obs = 1.
      
      # background universe
      self.h = 0.705  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmC = 0.227
      self.OmB = 0.0449
      self.OmM = self.OmC + self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/universe_HandEtAl13/lin_wmap9eCMB_HandEtAl13_matterpower.dat"
      
      super(UnivHandEtAl13, self).__init__()


##################################################################################

# uses the OmM and h that Mariana used for her reconstructed velocities
class UnivMariana(Universe):
   
   def __init__(self):
      
      # redshift of observer
      self.a_obs = 1.
      
      # background universe
      self.h = 0.70  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmM = 0.29
      self.OmB = 0.0458571
      self.OmC = self.OmM - self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/universe_Mariana/lin_Mariana_matterpower.dat"
      
      super(UnivMariana, self).__init__()


##################################################################################

class UnivPlanck15(Universe):
   
   def __init__(self):
      
      # redshift of observer
      self.a_obs = 1.
      
      # background universe
      self.h = 0.6712  # H0 = 100.*h km s^-1 Mpc^-1
      self.OmB = 0.0493
      self.OmC = 0.267
      self.OmM = self.OmC + self.OmB
      
      # read linear power spectrum from CAMB
      self.path_lin_matterpower = "./input/universe_Planck15/camb/matterpower_z0_lin.dat"
      
      super(UnivPlanck15, self).__init__()
