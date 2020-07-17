from headers import *
##################################################################################

class Projection(object):

   def __init__(self, U, name="", nProc=1):
      # copy U
      self.U = U
      self.name=name
      self.nProc=nProc
      #self.aMin
      #self.aMax

      # interpolate the projection kernel, for speed
      nA = 501
      A = np.linspace(self.aMin, self.aMax, nA)
      F = np.array(map(self.fForInterp, A))

#      print "nProc=", self.nProc
#      tStart = time()
#      with sharedmem.MapReduce(np=self.nProc) as pool:
#         F = np.array(pool.map(self.fForInterp, A))
#      tStop = time()
#      print "para took", tStop-tStart

      self.f = interp1d(A, F, kind='linear', bounds_error=False, fill_value=0.)

   def __str__(self):
      return self.name
   
   
   # projection kernel, such that, e.g.:
   # kappa = int dchi f(a) delta
   # to be interpolated for speed
   def fForInterp(self, a):
      pass



   ##################################################################################

   def dWddelta(self, aMin, aMax):
      """response of the projected quantity (kappa, y2d, ...)
      to the mean overdensity along the line of sight,
      measured between aMin and aMax
      """
      integrand = lambda a: 3.e5/(self.U.hubble(1./a-1.) * a**2) * self.f(a)
      result = integrate.quad(integrand, aMin, aMax, epsabs=0., epsrel=1.e-2)[0]
      return result

   ##################################################################################
   
   def plotW(self):
      
      # choice of weight function
      n = 2 # n of n-pt function
      nW = n
      nd = 2*n-2.
      
      # range to plot
      Na = 501
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
      ComovDistToObs = np.array( map( self.U.bg.comoving_distance, Z ) )
      W = np.array( map( lambda a: self.f(a), A ) )
      
      H_A = self.U.hubble(1./A-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
      # compute weight for n-pt function
      F_a = 1. / (H_A * A**2)
      F_a *= W**nW / ComovDistToObs**nd
      #
      F_z = 1. / H_A
      F_z *= W**nW / ComovDistToObs**nd
      #
      F_chi = W**nW / ComovDistToObs**nd
      
      #return Z, W/H_A
      
      # projection kernel W
      fig=plt.figure(-1)
      ax=plt.subplot(111)
      ax.plot(Z, W/H_A, 'b', lw=2)
      ax.set_xlabel(r'$z$', fontsize=22)
      ax.set_ylabel(r'$W (z)$', fontsize=22)
      #fig.savefig("./figures/weight/W_cmblens.pdf")
      '''
      # projection kernel W
      fig=plt.figure(0)
      ax=plt.subplot(111)
      ax.plot(Z, W, 'b')
      ax.set_xlabel(r'$z$', fontsize=18)
      ax.set_ylabel(r'$W_\chi (z)$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      
      # weight for n-pt function...
      # per z interval
      fig=plt.figure(1)
      ax=plt.subplot(111)
      ax.loglog(Z, F_z, 'b')
      ax.set_xlabel(r'$z$', fontsize=18)
      ax.set_ylabel(r'$ \frac{c}{H} \frac{ W^{'+str(nW)+'} }{ \chi^{'+str(int(nd))+'} }$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      fig.savefig('./figures/plotW.pdf')
      
      # per a interval
      fig=plt.figure(2)
      ax=plt.subplot(111)
      ax.plot(A, F_a, 'b')
      ax.set_xlabel(r'$a$', fontsize=18)
      ax.set_ylabel(r'$ \frac{c}{H a^2} \frac{ W^{'+str(nW)+'} }{ \chi^{'+str(int(nd))+'} }$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      
      # per radial comoving distance interval
      fig=plt.figure(3)
      ax=plt.subplot(111)
      ax.plot(ComovDistToObs, F_chi, 'b')
      ax.set_xlabel(r'$\chi$', fontsize=18)
      ax.set_ylabel(r'$\frac{ W^{'+str(nW)+'} }{ \chi^{'+str(int(nd))+'} }$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      
      # per logk interval
      L = np.logspace(log10(1.), log10(1.e4), 11, 10.)
      #
      fig=plt.figure(4)
      ax=plt.subplot(111)
      ax.set_xlabel(r'$k$ [h/Mpc]', fontsize=18)
      ax.set_ylabel(r'$ \frac{l}{k^2}  \frac{ W^{'+str(nW)+'} }{ \chi^{'+str(int(nd))+'} }$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      ax.grid()
      #
      for il in range(len(L)):
         l = L[il]
         #print 'l=', l, 'k=', l/self.U.bg.comoving_distance(0.5)
         K = l / ComovDistToObs
         #
         F_k = l /K**2
         F_k *= W**nW / ComovDistToObs**nd
         #
         ax.loglog(K, F_k, label=r'l='+str(round(l, 1)))
      ax.legend(loc=1)
      
      
      # per logk interval, normalized to have a max of 1
      L = np.logspace(log10(1.), log10(1.e4), 11, 10.)
      #
      fig=plt.figure(5)
      ax=plt.subplot(111)
      ax.set_xlabel(r'$k$ [h/Mpc]', fontsize=18)
      ax.set_ylabel(r'$ \frac{l}{k^2}  \frac{ W^{'+str(nW)+'} }{ \chi^{'+str(int(nd))+'} }$', fontsize=18)
      ax.set_title(r'weight for '+str(n)+'-point function', fontsize=18)
      ax.grid()
      #
      for il in range(len(L)):
         l = L[il]
         #print 'l=', l, 'k=', l/self.U.bg.comoving_distance(0.5)
         K = l / ComovDistToObs
         #
         F_k = l /K**2
         F_k *= W**nW / ComovDistToObs**nd
         #
         ax.semilogx(K, F_k/max(F_k), label=r'l='+str(round(l, 1)))
      ax.legend(loc=4)
      '''
      
      plt.show()



##################################################################################
##################################################################################

class WeightY(Projection):
   """Compton-y projection
   """
   
   def __init__(self, U, name='y'):
      self.aMin = 0.2   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
      super(WeightY, self).__init__(U, name=name)
   
   def fForInterp(self, a):
      """Compton y projection kernel
      """
      return a


##################################################################################
##################################################################################

class WeightKsz(Projection):
   """Projection for kSZ, to convert the radial momentum q_r into dT/T
   """
   
   def __init__(self, U, name='ksz'):
      self.aMin = 1./(1.+10.)   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
      super(WeightKsz, self).__init__(U, name=name)

   def fForInterp(self, a):
      """kSZ projection kernel:
      bar{n}_e^0 * sigma_T / a^2
      """
      # comoving baryon mass density (ie physical baryon mass density today)
      result = self.U.rho_b(1./a-1.)   # [(Msun/h) / (Mpc/h)^3]

      # convert to comoving electron number density (ie physical number density today)
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      msun = 1.989e30   # solar mass (kg)
      factor = (me + nH_ne*mH + nHe_ne*mHe) / (msun/self.U.bg.h)   # [Msun/h]
      result /= factor  # [(Mpc/h)^3]

      # multiply by Thomson cross section (physical)
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      sigma_T *= (self.U.bg.h/mpc)**2  # convert to (Mpc/h)^2
      result *= sigma_T

      result /= a**2
      return result


##################################################################################
##################################################################################

class WeightLensSingle(Projection):
   """Lensing projection: single source. The default is z_source = 1.
   """
   
   def __init__(self, U, z_source=1., name='lens'):
      # a for mass func, biases, and projection
      self.z_source = z_source
      self.a_source = 1./(1.+z_source)
      #
      self.aMin = max(self.a_source, 1./11.)   # don't go further than z=10
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
      super(WeightLensSingle, self).__init__(U, name=name)
   
   
   def fForInterp(self, a):
      """lensing projection kernel
      for a single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
#      # Expression for flat cosmology
#      d_a = self.U.bg.comoving_distance(1./a-1.)
#      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m * d_a / a
#      wlensing *= 1. - d_a/self.dist_source

      # Expression for flat/curved cosmology
      z = 1./a - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(self.z_source)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(self.z_source)
      d_LS = self.U.fK(chi_S - chi_L)
      #
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing
   
   def testHandEtAl13_fig1(self):
      # fig 1 from Hand et al 2013
      Data = np.genfromtxt('./input/tests/HandEtAl13/HandEtAl13_fig1_cmb.txt')
      Z = Data[:, 0]
      A = 1./(1.+Z)
      Wgal_ref = Data[:, 1]
      # from my code
      Zme = np.linspace(0., 10., 501)
      Ame = 1./(1.+Zme)
      fW = lambda a: self.f(a) * (3.e5/self.U.hubble(1./a-1.))
      Wgal_me = np.array(map(fW, Ame))
      Wgal_me /= np.max(Wgal_me)

      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.plot(Z, Wgal_ref, 'b-', label=r'Hand et al 2013')
      ax.plot(Zme, Wgal_me, 'r-', label=r'my code')
      #
      ax.legend(loc=4)
      ax.set_xlabel(r'$z$', fontsize=18)
      ax.set_ylabel(r'$W_{CMB}$ such that $\kappa = \int dz W(z)$', fontsize=18)
      ax.set_title(r'Fig1 from Hand et al 2013')
      #ax.set_xlim((0., 4.))
      ax.set_ylim((0., 1.))
      #fig.savefig('./figures/tests/HandEtAl13_fig1_cmb.pdf')

      plt.show()


##################################################################################
##################################################################################

class WeightLensOguriTakada11(Projection):
   """lensing projection: source distribution from Oguri & Takada 2011
   """
   
   def __init__(self, U, z_source=1., name='lens'):
      # parameters for source distribution (Oguri & Takada 2011)
      self.z0 = z_source/3. # should be <z_source>/3
      self.nz0 = 20.   # width of integral over sources
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + self.nz0 * self.z0 ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
      super(WeightLensOguriTakada11, self).__init__(U, name=name)


   def fdpdz(self, z):
      """source distribution from Oguri & Takada 2011
      int_0^inf dz dPdz = 1
      """
      result = 0.5 * z**2 / self.z0**3
      result *= exp( -z / self.z0)
      return result
   
   
   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing


   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result


##################################################################################
##################################################################################

class WeightLensHandEtAl13(Projection):
   """lensing projection: source distribution from Hand et al 2013
   """
   
   def __init__(self, U, name='lens'):
      # source distribution (eq6 from Hand et al 2013)
      A = 0.688
      a = 0.531
      b = 7.810
      c = 0.517
      fdpdz_nonorm = lambda z: A*(z**a + z**(a*b))/(z**b + c)
      norm = integrate.quad(fdpdz_nonorm, 0., np.inf, epsabs=0, epsrel=1.e-2)[0]
      self.fdpdz = lambda z: fdpdz_nonorm(z) / norm
      
      # a for mass func, biases, and projection
      self.aMin = 1./(1.+10.)  # arbitrary for now
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
   
      super(WeightLensHandEtAl13, self).__init__(U, name=name)

   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing

   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def testHandEtAl13_fig1(self):
      # fig 1 from Hand et al 2013
      Data = np.genfromtxt('./input/tests/HandEtAl13/HandEtAl13_fig1_gal.txt')
      Z = Data[:, 0]
      A = 1./(1.+Z)
      Wgal_ref = Data[:, 1]
      # from my code
      fW = lambda a: self.f(a) * (3.e5/self.U.hubble(1./a-1.))
      Wgal_me = np.array(map(fW, A))
      Wgal_me /= np.max(Wgal_me)
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.plot(Z, Wgal_ref, 'b-', label=r'Hand et al 2013')
      ax.plot(Z, Wgal_me, 'r-', label=r'my code')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=18)
      ax.set_ylabel(r'$W_{gal}$ such that $\kappa = \int dz W(z)$', fontsize=18)
      ax.set_title(r'Fig1 from Hand et al 2013')
      #fig.savefig('./figures/tests/HandEtAl13_fig1_gal.pdf')
      
      plt.show()


   def testHandEtAl13_fig2(self):
      # fig 2 from Hand et al 2013
      Data = np.genfromtxt('./input/tests/HandEtAl13/HandEtAl13_fig2.txt')
      # my interpolation
      Z = np.linspace(0., 3., 501)
      Me = np.array(map(self.fdpdz, Z))
      
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.plot(Data[:,0], Data[:, 1]/np.max(Data[:, 1]), 'b-', label=r'data from fig2')
      ax.plot(Z, Me/np.max(Me), 'r-', label=r'fit from eq6')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=18)
      ax.set_ylabel(r'$dn/dz$ normalized by its maximum value', fontsize=18)
      ax.set_title(r'Fig2 from Hand et al 2013')
      #fig.savefig('./figures/tests/HandEtAl13_fig2.pdf')
      
      plt.show()


##################################################################################
##################################################################################

class WeightLensDasEtAl13(Projection):
   """lensing projection: source distribution from Das Errard Spergel 2013
   """
   
   def __init__(self, U, name='lens'):
      # source distribution (eq11 from Das Errard Spergel 2013)
      z0 = 0.69   # ie median z is 1
      self.fdpdz = lambda z: 1.5 * z**2/z0**3 * np.exp(-(z/z0)**1.5)
      # a for mass func, biases, and projection
      self.aMin = 1./(1.+10.)  # arbitrary for now
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
      super(WeightLensDasEtAl13, self).__init__(U, name=name)
   
   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing

   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def plot(self):
      Na = 501
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
   
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.hubble(1./a-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
      fig=plt.figure(-1)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dPdz/np.max(dPdz), 'r', lw=2, label=r'source')
      ax.plot(Z, (W/H_A)/np.max(W/H_A), 'b', lw=2, label=r'$\kappa$')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=22)
      ax.set_ylabel(r'$W (z)$', fontsize=22)
      #fig.savefig("./figures/weight/W_cmblens.pdf")
      
      plt.show()

##################################################################################
##################################################################################

class WeightLensCustom(Projection):
   
   def __init__(self, U, fdndz, m=lambda z: 0., zMinG=1.e-4, zMaxG=2., name='lens', nProc=1):
      '''Here zMinG and zMaxG are those of the galaxy sample,
      not those of the resulting lensing kernel.
      m is the shear multiplicative bias.
      '''
      # bounds for the tracer sample (dn/dz)
      self.aMinG = 1./(1.+zMaxG)
      self.aMaxG = 1./(1.+zMinG)
      
      # bounds for the lensing kernel
      self.aMin = self.aMinG
      self.aMax = 1.-0.005
      
      # shear multiplicative bias m(z)
      self.m = m
      
      # fdndz doesn't need to be normalized to anything, for lensing purposes
      self.fdndz = lambda z: fdndz(z) * (z>=zMinG) * (z<=zMaxG)
      self.ngal = integrate.quad(self.fdndz, zMinG, zMaxG, epsabs=0., epsrel=1.e-2)[0]
      # dpdz normalized such that int dz dpdz = 1
      self.fdpdz = lambda z: self.fdndz(z) / self.ngal
   
      super(WeightLensCustom, self).__init__(U, name=name, nProc=nProc)

   
#   def fForInterp(self, a):
#      """Only valid in the absence of curvature.
#      """
#      if a<self.aMinG:
#         return 0.
#      else:
#         d_a = self.U.bg.comoving_distance(1./a-1.)
#         result = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m * d_a / a
#         integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.bg.comoving_distance(1./a_s-1.))
#         result *= integrate.quad(integrand, self.aMinG, a, epsabs=0, epsrel=1.e-2)[0]
#         result *= (1.+self.m(1./a-1.))  # shear multiplicative bias
#         return result

   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source,
      valid even with curvature.
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing


   def fForInterp(self, a):
      """Valid even with curvature.
      """
      if a<self.aMinG:
         return 0.
      else:
         integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
         result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
         result *= (1.+self.m(1./a-1.))  # shear multiplicative bias
         return result

   def plot(self):
      Na = 501
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
   
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.hubble(1./A-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
      fig=plt.figure(-1)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dPdz/np.max(dPdz), 'r', lw=2, label=r'source')
      ax.plot(Z, (W/H_A)/np.max(W/H_A), 'b', lw=2, label=r'$\kappa$')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=22)
      ax.set_ylabel(r'$W (z)$', fontsize=22)
      #fig.savefig("./figures/weight/W_cmblens.pdf")
      
      plt.show()

##################################################################################
##################################################################################

class WeightLensCIBSchmidt15(Projection):
   """lensing projection: CIB source distribution from Schmidt Menard Scranton+15
   the values of z0 and alpha are in table 2
   approximate calculation: assumes that the CIB monopole redshift distribution
   is the relevant source distribution for CIB lensing
   """
   
   def __init__(self, U, z0=1., alpha=1., name='lens'):
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + 10. ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
      #
      # source distribution from Schmidt Menard Scranton+15
      fdpdzNonNormalized = lambda z: z**alpha * np.exp( -(z / z0)**alpha)
      # normalize to have int_zMin^zMax dz dPdz = 1
      norm = integrate.quad(fdpdzNonNormalized, 1./self.aMax, 1./self.aMin, epsabs=0, epsrel=1.e-4)[0]
      self.fdpdz = lambda z: fdpdzNonNormalized(z) / norm
      #
      super(WeightLensCIBSchmidt15, self).__init__(U, name=name)
   
   
   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing


   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result
   

   def plot(self):
      Na = 501
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
      
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.hubble(1./a-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
      fig=plt.figure(-1)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dPdz/np.max(dPdz), 'r', lw=2, label=r'source')
      ax.plot(Z, (W/H_A)/np.max(W/H_A), 'b', lw=2, label=r'$\kappa$')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=22)
      ax.set_ylabel(r'$W (z)$', fontsize=22)
      #fig.savefig("./figures/weight/W_cmblens.pdf")
      
      plt.show()


##################################################################################
##################################################################################

class WeightLensCIBPullen17(Projection):
   """lensing projection: CIB source distribution from Pullen+17
   these have been digitized from figure 6
   nu should be 353, 545 or 857 (in GHz)
   approximate calculation: assumes that the CIB monopole redshift distribution
   is the relevant source distribution for CIB lensing
   """
   
   def __init__(self, U, nu=353, name='lens'):
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + 5. ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
      
      # read digitized values
      path = "./input/cib_zdist_pullen17/Pullen+17_"+str(nu)+".txt"
      data = np.genfromtxt(path)
      fdpdzNonNormalized = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
      # normalize to have int_zMin^zMax dz dPdz = 1
      norm = integrate.quad(fdpdzNonNormalized, 1./self.aMax, 1./self.aMin, epsabs=0, epsrel=1.e-4)[0]
      self.fdpdz = lambda z: fdpdzNonNormalized(z) / norm
      
      super(WeightLensCIBPullen17, self).__init__(U, name=name)
      
   
   def fSingleSource(self, a, aSource):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      # Expression valid for curved cosmology
      z = 1./a - 1.
      zSource = 1./aSource - 1.
      d_L = self.U.bg.comoving_transverse_distance(z)
      d_S = self.U.bg.comoving_transverse_distance(zSource)
      chi_L = self.U.bg.comoving_distance(z)
      chi_S = self.U.bg.comoving_distance(zSource)
      d_LS = self.U.fK(chi_S - chi_L)
      
      wlensing = 1.5 * (100./3.e5)**2 * self.U.bg.Omega0_m
      wlensing *= d_L/a * d_LS / d_S
      return wlensing


   def fForInterp(self, a):
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * self.fSingleSource(a, a_s)
      result = integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def plot(self):
      Na = 501
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
      
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.hubble(1./a-1.) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
      fig=plt.figure(-1)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dPdz/np.max(dPdz), 'r', lw=2, label=r'source')
      ax.plot(Z, (W/H_A)/np.max(W/H_A), 'b', lw=2, label=r'$\kappa$')
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$z$', fontsize=22)
      ax.set_ylabel(r'$W (z)$', fontsize=22)
      #fig.savefig("./figures/weight/W_cmblens.pdf")
      
      plt.show()




##################################################################################
##################################################################################

class WeightCIBPlanck15(Projection):
   """CIB projection: from Planck XXIII 2015, Planck XXX 2013
   the split between projection kernel and profile is
   somewhat arbitrary here, given that there isn't a well-defined 3d quantity.
   I could have put all the z-dependence in the profile,
   and make the projection kernel trivial
   """
   
   def __init__(self, U, name='cibplanck'):
      # a for mass func, biases, and projection
      self.aMin = 1./(1.+10.)  # arbitrary for now
      epsilon = 1.e-5
      self.aMax = 1.-epsilon
   
      super(WeightCIBPlanck15, self).__init__(U, name=name)
   
   def fForInterp(self, a):
      return self.U.bg.comoving_distance(1./a-1.)**2


##################################################################################
##################################################################################

class WeightTracer(Projection):
   """Projected density field of tracers
   requires defining b(z) and dn/dz(z)
   """
   
   def __init__(self, U, name='d'):
      self.name = name
      # normalization of dn/dz, ie number of galaxies per unit steradian
      self.ngal = integrate.quad(self.dndz, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-2)[0]
      # convert to number of galaxies per square arcmin
      self.ngal_per_arcmin2 = self.ngal * (np.pi/180./60.)**2
      # print number of gal per square arcmin
      print(self.name+": "+str(np.round(self.ngal_per_arcmin2, 2))+"gal/arcmin^2")
      
      super(WeightTracer, self).__init__(U, name=name)

   def b(self, z):
      """tracer bias
      """
      pass
   
   def dndz(z):
      """normalized such that int dz dn/dz = ngal,
      ie the number of gals per unit steradian
      """
      pass

   def fForInterp(self, a):
      """projection kernel
      """
      z = 1./a - 1.
      result = self.U.hubble(1./a-1.) / 3.e5
      result *= self.dndz(z)
      result /= self.ngal
      result *= self.b(z)
      return result


   def plotDndz(self, wArr=None, logx=False, logy=False):
      if wArr is None:
         wArr = [self]
      
      # plot dN / dz / dOmega [arcmin^-2]
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for w in wArr:
         Z = np.linspace(1./w.aMax-1., 1./w.aMin-1., 501)
         Dndz = np.array(map(w.dndz, Z))
         # normalize such that int dz dn/dz = ngal in arcmin^-2
         Dndz /= (180.*60. / np.pi)**2
         # plot it
         ax.plot(Z, Dndz, label=w.name)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      if logx:
         ax.set_xscale('log', nonposx='clip')
      if logy:
         ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$dN/dz d\Omega$ [arcmin$^{-2}$]')


      # plot dN / dV [(Mpc/h)^-3]
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      for w in wArr:
         Z = np.linspace(1./w.aMax-1., 1./w.aMin-1., 501)
         Dndz = np.array(map(w.dndz, Z))  # [sr^-2]
         Dndz /= self.U.bg.comoving_distance(Z)**2 * (3.e5 / self.U.hubble(Z))
         # plot it
         ax.plot(Z, Dndz, label=w.name)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      if logx:
         ax.set_xscale('log', nonposx='clip')
      if logy:
         ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$dN/dV$ [(Mpc/h)$^{-3}$]')


      # plot dN / dz [dimless], assuming full sky
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      for w in wArr:
         Z = np.linspace(1./w.aMax-1., 1./w.aMin-1., 501)
         Dndz = np.array(map(w.dndz, Z))
         # normalize such that int dz dn/dz = ngal in arcmin^-2
         Dndz /= (180.*60. / np.pi)**2
         # convert to dN/dz, assumng full sky is covered
         Dndz *= 4.*np.pi * (180.*60./np.pi)**2
         # plot it
         ax.plot(Z, Dndz, label=w.name)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      if logx:
         ax.set_xscale('log', nonposx='clip')
      if logy:
         ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$dN/dz$ for $f_\text{sky}=1$')

      plt.show()


   def splitBins(self, nBins):
      '''Splits the dn/dz into nBins bins,
      with equal number of objects.
      '''
      nBins = np.int(nBins)
      
      def cumulZDist(zMax):
         f = lambda z: self.dndz(z) / self.ngal
         result = integrate.quad(f, 1./self.aMax-1., zMax, epsabs=0., epsrel=1.e-2)[0]
         return result
      
      # fill an array with the z-bounds of the bins
      zBounds = np.zeros(nBins+1)
      zBounds[0] = 1./self.aMax-1.
      zBounds[-1] = 1./self.aMin-1.
      for iBin in range(nBins-1):
         # force equal number of objects per bins
         f = lambda zMax: cumulZDist(zMax) - (iBin+1.)/nBins
         zBounds[iBin+1] = optimize.brentq(f , zBounds[0], zBounds[-1])
      return zBounds

   def bMean(self):
      '''Computes the dn/dz-weighted galaxy bias.
      '''
      f = lambda z: self.dndz(z) / self.ngal * self.b(z)
      result = integrate.quad(f, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-2)[0]
      # the normalization should be 1, but just in case:
      f = lambda z: self.dndz(z) / self.ngal
      result /= integrate.quad(f, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-2)[0]
      return result

   def zMean(self):
      '''Computes the dn/dz-weighted mean redshift.
      '''
      f = lambda z: self.dndz(z) / self.ngal * z
      result = integrate.quad(f, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-2)[0]
      # the normalization should be 1, but just in case:
      f = lambda z: self.dndz(z) / self.ngal
      result /= integrate.quad(f, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-2)[0]
      return result


##################################################################################
##################################################################################

class WeightTracerCMASS(WeightTracer):
   """Projected number density of CMASS DR12 galaxies
   """

   def __init__(self, U, name='cmass'):
      self.aMin = 1./(1.+0.7)   # min bound for integral over a
      self.aMax = 1./(1.+0.4)   # max bound for integral over a

      # tracer bias
      self.b = lambda z: 2.

      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      data = np.genfromtxt("./input/dndz/cmass_dndz.txt")
      Z = data[:,0]
      Dndz = data[:,1]
      f = UnivariateSpline(Z, Dndz, k=1, s=0)
      self.dndz = lambda z: f(z) * (z>=np.min(Z)) * (z<=np.max(Z))

      super(WeightTracerCMASS, self).__init__(U, name=name)


##################################################################################
##################################################################################

class WeightTracerWISE(WeightTracer):
   """Projected number density of WISE galaxies
   """
   
   def __init__(self, U, name='wise'):
      self.aMin = 1./(1.+1.)   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
   
      # tracer bias
      self.b = lambda z: 1.2   #1.
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      data = np.genfromtxt("./input/dndz/wise_dndz.txt")
      Z = data[:,0]
      Dndz = data[:,1]
      f = UnivariateSpline(Z, Dndz, k=1, s=0)
      self.dndz = lambda z: f(z) * (z>=np.min(Z)) * (z<=np.max(Z))

      super(WeightTracerWISE, self).__init__(U, name=name)

##################################################################################
##################################################################################

class WeightTracerLSSTGold(WeightTracer):
   """Projected number density of LSST gold galaxies,
   From the LSST Science book, chapter 3 and 13.
   """
   
   def __init__(self, U, name='lsstgold', iLim=25.3):
      self.aMin = 1./(1.+3.)   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
      
      self.iLim = iLim  # limiting i-band magnitude
      
      # tracer bias
      self.b = lambda z: 1 + 0.84*z
      
      self.ngal_per_arcmin2 = 46.*10**(0.31*(iLim-25.)) # galaxies per squared arcmin
      self.ngal = self.ngal_per_arcmin2 / (np.pi/180./60.)**2
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.z0 = 0.0417*iLim - 0.744
      # the normalization to ngal below is approximate, but correct to better than 1%
      self.dndz = lambda z: self.ngal * (z/self.z0)**2 * np.exp(-z/self.z0) / (2.*self.z0)

      super(WeightTracerLSSTGold, self).__init__(U, name=name)


   def magnificationBias(self, z):
      '''Implements the magnification bias from Joachimi Bridle 2010.
      Computes alpha = dn/dS, such that:
      W = W_g + 2(alpha-1)*W_kappa.
      '''
      b11 = 0.44827
      b12 = -1.
      b13 = 0.05617
      b14 = 0.07704
      b15 = -11.3768
      result = b11 + b12*(b13*self.iLim - b14)**b15
      #
      b21 = 0.
      b22 = 1.
      b23 = 0.19658
      b24 = 3.31359
      b25 = -2.5028
      result += (b21 + b22*(b23*self.iLim - b24)**b25 ) * z
      #
      b31 = 0.
      b32 = 1.
      b33 = 0.18107
      b34 = 3.05213
      b35 = -2.5027
      result += (b31 + b32*(b33*self.iLim - b34)**b35 ) * z**2
      return result



##################################################################################
##################################################################################

class WeightTracerLSSTSourcesSchaanKrauseEifler16(WeightTracer):
   """Projected number density of LSST source galaxies,
   used for shear measurements,
   as in Schaan Krause Eifler +16.
   """
   
   def __init__(self, U, name='lsstsources', zMin=0.005, zMax=4., ngal=26.):
      self.zMin = zMin
      self.zMax = zMax
      self.aMin = 1./(1.+self.zMax)   # min bound for integral over a
      self.aMax = 1./(1.+self.zMin)   # max bound for integral over a
      
      # tracer bias
      # copied from the LSST gold sample!
      self.b = lambda z: 1 + 0.84*z
      
      self.ngal_per_arcmin2 = ngal # galaxies per squared arcmin
      self.ngal = self.ngal_per_arcmin2 / (np.pi/180./60.)**2  # per steradian
      

      # dn/dz, non-normalized
      self.z0 = 0.5
      self.alpha = 1.27
      self.beta = 1.02
      f = lambda z: z**self.alpha * np.exp(-(z/self.z0)**self.beta)
      # normalization
      norm = integrate.quad(f, self.zMin, self.zMax, epsabs=0., epsrel=1.e-2)[0]
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = lambda z: self.ngal * f(z) / norm

      super(WeightTracerLSSTSources, self).__init__(U, name=name)



   def magnificationBias(self, z):
      '''Implements the magnification bias from Joachimi Bridle 2010.
      Computes alpha = dn/dS, such that:
      W = W_g + 2(alpha-1)*W_kappa.
      '''
      b11 = 0.44827
      b12 = -1.
      b13 = 0.05617
      b14 = 0.07704
      b15 = -11.3768
      result = b11 + b12*(b13*self.iLim - b14)**b15
      #
      b21 = 0.
      b22 = 1.
      b23 = 0.19658
      b24 = 3.31359
      b25 = -2.5028
      result += (b21 + b22*(b23*self.iLim - b24)**b25 ) * z
      #
      b31 = 0.
      b32 = 1.
      b33 = 0.18107
      b34 = 3.05213
      b35 = -2.5027
      result += (b31 + b32*(b33*self.iLim - b34)**b35 ) * z**2
      return result


##################################################################################
##################################################################################

class WeightTracerLSSTSourcesDESCSRDV1(WeightTracer):
   """Projected number density of LSST source galaxies,
   used for shear measurements,
   as in the DESC SRD v1.
   """
   
   def __init__(self, U, name='lsstsources', zMin=0.005, zMax=4., ngal=27.):
      self.zMin = zMin
      self.zMax = zMax
      self.aMin = 1./(1.+self.zMax)   # min bound for integral over a
      self.aMax = 1./(1.+self.zMin)   # max bound for integral over a
      
      # tracer bias
      # copied from the LSST gold sample, from LSST Science book
      self.b = lambda z: 0.95 / self.U.bg.scale_independent_growth_factor(z) #1 + 0.84*z
      
      self.ngal_per_arcmin2 = ngal # galaxies per squared arcmin
      self.ngal = self.ngal_per_arcmin2 / (np.pi/180./60.)**2  # per steradian
      

      # dn/dz, non-normalized
      self.z0 = 0.11
      self.alpha = 0.68
      f = lambda z: z**2. * np.exp(-(z/self.z0)**self.alpha)
      # normalization
      norm = integrate.quad(f, self.zMin, self.zMax, epsabs=0., epsrel=1.e-2)[0]
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = lambda z: self.ngal * f(z) / norm

      super(WeightTracerLSSTSourcesDESCSRDV1, self).__init__(U, name=name)



   def magnificationBias(self, z):
      '''Implements the magnification bias from Joachimi Bridle 2010.
      Computes alpha = dn/dS, such that:
      W = W_g + 2(alpha-1)*W_kappa.
      '''
      b11 = 0.44827
      b12 = -1.
      b13 = 0.05617
      b14 = 0.07704
      b15 = -11.3768
      result = b11 + b12*(b13*self.iLim - b14)**b15
      #
      b21 = 0.
      b22 = 1.
      b23 = 0.19658
      b24 = 3.31359
      b25 = -2.5028
      result += (b21 + b22*(b23*self.iLim - b24)**b25 ) * z
      #
      b31 = 0.
      b32 = 1.
      b33 = 0.18107
      b34 = 3.05213
      b35 = -2.5027
      result += (b31 + b32*(b33*self.iLim - b34)**b35 ) * z**2
      return result




##################################################################################
##################################################################################

class WeightTracerDESIBGS(WeightTracer):
   """Projected number density of DESI BGS,
   dn/dz from Rongpu Zhou Feb 19 2020.
   Bias b(z) = 1.34 / D(z) from DESI FDR arxiv:1611.00036v2, section 2.4.2
   BGS: expect 2.5e13 Msun halos, from DESI Final Design Report, fig 3.4
   """

   def __init__(self, U, name='desibgs'):
      
      # read dn/dz
      # nGalBin [gal / sq deg]
      path = "./input/dndz/dndz_DESI_BGS.txt"
      data = np.genfromtxt(path)
      zLow = data[:,0]
      zHigh = data[:,1]
      nGalBin = data[:,2]
      nGalBin /= (np.pi/180.)**2 # convert [gal/sq deg] to [gal/sr]

      # bounds for a integrals
      self.aMin = 1./(1.+np.max(zHigh))
      self.aMax = 1./(1.+1.e-3)
      
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      zMean = np.zeros(len(zLow) + 2)
      zMean[0] = zLow[0]
      zMean[1:-1] = 0.5 * (zLow + zHigh)
      zMean[-1] = zHigh[-1]
      #
      nGal = np.zeros(len(zLow)+2)
      nGal[1:-1] = nGalBin
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = interp1d(zMean, nGal, kind='linear', bounds_error=False, fill_value=0.)

      # tracer bias
      self.b = lambda z: 1.34 / U.bg.scale_independent_growth_factor(z)

      super(WeightTracerDESIBGS, self).__init__(U, name=name)


##################################################################################

class WeightTracerDESILRG(WeightTracer):
   """Projected number density of DESI LRG,
   dn/dz from Rongpu Zhou Feb 19 2020.
   Bias b(z) = 1.7 / D(z) from DESI FDR arxiv:1611.00036v2, section 2.4.2
   LRG: expect 2.e13 Msun halos from CMASS galaxies
   """

   def __init__(self, U, name='desilrg'):
      
      # read dn/dz
      # nGalBin [gal / sq deg]
      path = "./input/dndz/dndz_DESI_LRG.txt"
      data = np.genfromtxt(path)
      zLow = data[:,0]
      zHigh = data[:,1]
      nGalBin = data[:,2]
      nGalBin /= (np.pi/180.)**2 # convert [gal/sq deg] to [gal/sr]

      # bounds for a integrals
      self.aMin = 1./(1.+np.max(zHigh))
      self.aMax = 1./(1.+1.e-3)
      
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      zMean = np.zeros(len(zLow) + 2)
      zMean[0] = zLow[0]
      zMean[1:-1] = 0.5 * (zLow + zHigh)
      zMean[-1] = zHigh[-1]
      #
      nGal = np.zeros(len(zLow)+2)
      nGal[1:-1] = nGalBin
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = interp1d(zMean, nGal, kind='linear', bounds_error=False, fill_value=0.)

      # tracer bias
      self.b = lambda z: 1.7 / U.bg.scale_independent_growth_factor(z)

      super(WeightTracerDESILRG, self).__init__(U, name=name)


##################################################################################

class WeightTracerDESIELG(WeightTracer):
   """Projected number density of DESI ELG,
   dn/dz from Rongpu Zhou Feb 19 2020.
   Bias b(z) = 0.84 / D(z) from DESI FDR arxiv:1611.00036v2, section 2.4.2
   ELG: expect 2.e11 Msun halos; DESI Cosmo Sim Req Doc mentions 1.e11Msun as minimum mass
   """

   def __init__(self, U, name='desielg'):
      
      # read dn/dz
      # nGalBin [gal / sq deg]
      path = "./input/dndz/dndz_DESI_ELG.txt"
      data = np.genfromtxt(path)
      zLow = data[:,0]
      zHigh = data[:,1]
      nGalBin = data[:,2]
      nGalBin /= (np.pi/180.)**2 # convert [gal/sq deg] to [gal/sr]

      # bounds for a integrals
      self.aMin = 1./(1.+np.max(zHigh))
      self.aMax = 1./(1.+1.e-3)
      
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      zMean = np.zeros(len(zLow) + 2)
      zMean[0] = zLow[0]
      zMean[1:-1] = 0.5 * (zLow + zHigh)
      zMean[-1] = zHigh[-1]
      #
      nGal = np.zeros(len(zLow)+2)
      nGal[1:-1] = nGalBin
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = interp1d(zMean, nGal, kind='linear', bounds_error=False, fill_value=0.)

      # tracer bias
      self.b = lambda z: 0.84 / U.bg.scale_independent_growth_factor(z)

      super(WeightTracerDESIELG, self).__init__(U, name=name)

##################################################################################

class WeightTracerDESIQSO(WeightTracer):
   """Projected number density of DESI QSO,
   dn/dz from Rongpu Zhou Feb 19 2020.
   Bias b(z) = 1.2 / D(z) from DESI FDR arxiv:1611.00036v2, section 2.4.2
   QSO: guess 5.e12 Msun halos.
   """

   def __init__(self, U, name='desiqso'):
      
      # read dn/dz
      # nGalBin [gal / sq deg]
      path = "./input/dndz/dndz_DESI_QSO.txt"
      data = np.genfromtxt(path)
      zLow = data[:,0]
      zHigh = data[:,1]
      nGalBin = data[:,2]
      nGalBin /= (np.pi/180.)**2 # convert [gal/sq deg] to [gal/sr]

      # bounds for a integrals
      self.aMin = 1./(1.+np.max(zHigh))
      self.aMax = 1./(1.+1.e-3)
      
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      zMean = np.zeros(len(zLow) + 2)
      zMean[0] = zLow[0]
      zMean[1:-1] = 0.5 * (zLow + zHigh)
      zMean[-1] = zHigh[-1]
      #
      nGal = np.zeros(len(zLow)+2)
      nGal[1:-1] = nGalBin
      
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = interp1d(zMean, nGal, kind='linear', bounds_error=False, fill_value=0.)

      # tracer bias
      self.b = lambda z: 1.2 / U.bg.scale_independent_growth_factor(z)

      super(WeightTracerDESIQSO, self).__init__(U, name=name)


##################################################################################
##################################################################################

class WeightTracerSPHEREx(WeightTracer):
   """SPHEREx, from text file from Olivier.
   iSample selects the sample, depending on the redshift quality sigma_z/(1+z).
   sigma_z/(1+z) in [0., 0.003], [0.003, 0.01], [0.01, 0.03], [0.03, 0.1], [0.1, 0.2]
   for the sample 0, 1, 2, 3, 4.
   Ex: iSample = [0], iSample = [4], iSample = [0,1,2]
   """

   def __init__(self, U, name='spherex', ISample=0):
      
      # sample chosen
      self.ISample = ISample
      # update the name
      name += ''.join([str(iSample) for iSample in ISample])

      # galaxy samples
      # sz = sigma(z)/(1+z)
      nSample = 5
      sZMin = np.array([0., 0.003, 0.01, 0.03, 0.1])
      sZMax = np.array([0.003, 0.01, 0.03, 0.1, 0.2])
      # keep min and max redshift uncertainties
      self.sZMin = np.min(sZMin[ISample])
      self.sZMax = np.max(sZMax[ISample])
      

      #redshift bins
      zLow = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.6, 2.2, 2.8, 3.4, 4.0])
      zHigh = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.6, 2.2, 2.8, 3.4, 4.0, 4.6])
      nZBin = len(zLow)
      
      # bounds for a integrals
      self.aMin = 1./(1.+np.max(zHigh))
      self.aMax = 1./(1.+1.e-3)
      
      # mean nb densities
      # nbar [(Mpc/h)^-3]
      nBarSamples = np.zeros((nSample, nZBin))
      nBarSamples[0,:] = np.array([6.210e-03, 1.540e-03, 4.930e-05, 1.430e-06, 1.000e-09, 1.000e-09, 1.000e-09, 1.000e-09, 1.000e-09, 2.210e-07, 1.000e-09])
      nBarSamples[1,:] = np.array([8.200e-03, 4.580e-03, 7.180e-04, 1.190e-04, 2.310e-05, 5.250e-07, 2.170e-07, 4.300e-07, 1.000e-09, 1.000e-09, 1.000e-09])
      nBarSamples[2,:] = np.array([8.110e-03, 6.150e-03, 2.360e-03, 1.440e-03, 6.550e-04, 4.200e-05, 3.910e-06, 3.700e-06, 1.000e-09, 4.420e-07, 1.000e-09])
      nBarSamples[3,:] = np.array([1.340e-02, 1.030e-02, 4.320e-03, 2.840e-03, 2.120e-03, 3.060e-04, 4.240e-05, 2.400e-05, 1.040e-05, 6.640e-07, 4.850e-07])
      nBarSamples[4,:] = np.array([1.150e-02, 6.390e-03, 2.580e-03, 2.930e-03, 2.330e-03, 5.050e-04, 1.480e-04, 9.150e-05, 4.630e-05, 2.430e-06, 1.160e-06])
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      zMean = np.zeros(len(zLow) + 2)
      zMean[0] = zLow[0]
      zMean[1:-1] = 0.5 * (zLow + zHigh)
      zMean[-1] = zHigh[-1]
      #
      nBar = np.zeros(len(zLow)+2)
      nBar[1:-1] = np.sum([nBarSamples[iSample,:] for iSample in ISample], axis=0)   # [(Mpc/h)^-3]
      # convert to [gal/arcmin^2]
      nBar *= U.bg.comoving_distance(zMean)**2 * (3.e5 / U.hubble(zMean))
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = interp1d(zMean, nBar, kind='linear', bounds_error=False, fill_value=0.)

      # tracer bias
      biasSamples = np.zeros((nSample, nZBin))
      biasSamples[0,:] = np.array([0.860, 1.100, 2.300, 4.400, 1.000, 1.000, 1.000, 1.000, 1.000, 7.300, 1.000])
      biasSamples[1,:] = np.array([0.790, 0.910, 1.300, 2.000, 2.900, 5.300, 6.200, 6.200, 1.000, 1.000, 1.000])
      biasSamples[2,:] = np.array([0.760, 0.840, 1.100, 1.300, 1.600, 2.900, 4.500, 4.900, 1.000, 6.700, 1.000])
      biasSamples[3,:] = np.array([0.730, 0.790, 0.940, 1.100, 1.200, 2.000, 3.200, 3.900, 4.800, 6.300, 7.200])
      biasSamples[4,:] = np.array([0.720, 0.780, 0.910, 1.000, 1.100, 1.700, 2.600, 3.200, 3.900, 5.700, 6.500])
      # replace the piecewise constant dn/dz by an approximate piecewise linear function
      # get the mean bias, weighted by number density
      bias = np.zeros(len(zLow)+2)
      bias[1:-1] = np.sum([biasSamples[iSample,:] * nBarSamples[iSample,:]  for iSample in ISample], axis=0)
      bias[1:-1] /= np.sum([nBarSamples[iSample,:]  for iSample in ISample], axis=0)
      # interpolate tracer bias
      self.b = interp1d(zMean, bias, kind='linear', bounds_error=False, fill_value=0.)

      super(WeightTracerSPHEREx, self).__init__(U, name=name)



##################################################################################
##################################################################################

class WeightTracerCustom(WeightTracer):
   """Projected density field of tracers
   requires defining functions b(z) and dn/dz(z)
   """
   
   def __init__(self, U, b, dndz, zMin=1.e-4, zMax=2., name='dcustom'):
      self.aMin = 1./(1.+zMax)
      self.aMax = 1./(1.+zMin)
      
      # bias as a function of z
      self.b = b
      # dndz as a function of z
      # normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = dndz
      
      super(WeightTracerCustom, self).__init__(U, name=name)


##################################################################################
##################################################################################

class WeightCIBPenin12(Projection):
   """Projection kernel for CIB, for P2h and P1h (Pshot treated separately)
   where P3d is the power spectrum of number density of IR galaxies
   from Penin Dore Lagache Bethermin 2012
   DOI: 10.1051/0004-6361/201117489
   Uses galaxy luminosity functions from Bethermin+12,
   available at http://irfu.cea.fr/Sap/Phocea/Page/index.php?id=537
   """
   
   def __init__(self, U, nu=217.e9, fluxCut=160.e-3, name='cibpenin12'):
      self.U = U
      self.nu = nu   # in Hz
      self.fluxCut = fluxCut  # in Jy
      
      # read the Bethermin+12 flux number counts
      # dNdSnudzdOmega in gal/Jy/sr,
      # improperly called dN / dSnu dz in Bethermin+12 and Penin+12
      self.Z = np.genfromtxt("./input/cib_bethermin12_2sfm/converted/z.txt")
      self.A = 1./(1.+self.Z)
      self.Snu = np.genfromtxt("./input/cib_bethermin12_2sfm/converted/Snu.txt") # in Jy
      self.dNdSnudzdOmega = np.genfromtxt("./input/cib_bethermin12_2sfm/converted/dNdSnudz_Planck"+str(int(nu/1.e9))+"GHz.txt")
      # put A in growing order
      self.A = self.A[::-1]
      self.Z = self.Z[::-1]
      self.dNdSnudzdOmega = self.dNdSnudzdOmega[::-1,:]
      # implement the flux cut
      iFluxCut = np.argmin(np.abs(self.Snu - self.fluxCut))
#      print self.fluxCut, iFluxCut, len(self.Snu), self.Snu[iFluxCut], self.Snu[-1]
#      print len(self.Snu)
      self.Snu = self.Snu[:iFluxCut]
      self.dNdSnudzdOmega = self.dNdSnudzdOmega[:, :iFluxCut]
#      print len(self.Snu)

      # convert flux number counts to gal/Jy/(Mpc/h)^3
      # dNdSnudV in gal/Jy/(Mpc/h)^3
      Chi = np.array(map(lambda a: self.U.bg.comoving_distance(1./a-1.), self.A))
      Hubble = np.array(map(lambda a: self.U.hubble(1./a-1.), self.A))
      dV_dzdOmega = Chi**2 * (3.e5/Hubble)
      self.dNdSnudV = self.dNdSnudzdOmega / dV_dzdOmega[:,np.newaxis]
      
   
      '''
      # Beware of 2d interpolations!!!!!! They are the devil
      # convert flux number counts to gal/Jy/(Mpc/h)^3
      # interpolate flux number count
      print("interpolating Bethermin+12 flux number count")
      fordNdSnudV = RectBivariateSpline(np.log(A), np.log(Snu), dNdSnudz, s=0)
      self.dNdSnudV = lambda snu, a: fordNdSnudV(np.log(a), np.log(snu)) *\
                                    (a>=np.min(A))*(a<=np.max(A))*\
                                    (snu>=np.min(Snu))*(snu<=np.max(Snu))*\
                                    4.*np.pi*self.U.bg.comoving_distance(1./a-1.)**2* 3.e5/self.U.hubble(1./a-1.)
      '''

      
      # compute emissivity tables, in Jy^p / (Mpc/h)^3
      # mean flux^p per unit volume, for projection kernel and shot noises
      self.JNu1 = np.trapz(self.dNdSnudV * self.Snu[np.newaxis, :], self.Snu, axis=1)
      self.JNu2 = np.trapz(self.dNdSnudV * self.Snu[np.newaxis, :]**2, self.Snu, axis=1)
      self.JNu3 = np.trapz(self.dNdSnudV * self.Snu[np.newaxis, :]**3, self.Snu, axis=1)
      self.JNu4 = np.trapz(self.dNdSnudV * self.Snu[np.newaxis, :]**4, self.Snu, axis=1)
      
      # z bounds
      self.aMin = np.min(self.A)
      self.aMax = np.max(self.A)
      # flux bounds
      self.snuMin = np.min(self.Snu)
      self.snuMax = np.max(self.Snu)
      
      # interpolate emissivities
      forjNu1 = UnivariateSpline(self.A, self.JNu1, k=1, s=0)
      self.jNu1 = lambda a: forjNu1(a) * (a>=self.aMin)*(a<=self.aMax)
      #
      forjNu2 = UnivariateSpline(self.A, self.JNu2, k=1, s=0)
      self.jNu2 = lambda a: forjNu2(a) * (a>=self.aMin)*(a<=self.aMax)
      #
      forjNu3 = UnivariateSpline(self.A, self.JNu3, k=1, s=0)
      self.jNu3 = lambda a: forjNu3(a) * (a>=self.aMin)*(a<=self.aMax)
      #
      forjNu4 = UnivariateSpline(self.A, self.JNu4, k=1, s=0)
      self.jNu4 = lambda a: forjNu4(a) * (a>=self.aMin)*(a<=self.aMax)
   
      super(WeightCIBPenin12, self).__init__(U, name=name+'_'+str(int(nu/1.e9))+'GHZ')

   def fForInterp(self, a):
      """projection kernel
      """
      result = self.U.bg.comoving_distance(1./a-1.)**2
      result *= self.jNu1(a)
      return result
   
   def fdPshotNoise_da(self, a, l):
      """contribution of each scale factor to the shot noise
      """
      result = (3.e5/self.U.hubble(1./a-1.)) / a**2
      result *= self.U.bg.comoving_distance(1./a-1.)**2 * self.jNu2(a)
      return result
   
   def fPshotNoise(self, l):
      """Shot noise for CIB 2d power spectrum
      """
      integrand = lambda a: (3.e5/self.U.hubble(1./a-1.)) / a**2 *\
                           self.U.bg.comoving_distance(1./a-1.)**2 * self.jNu2(a)
      result = integrate.quad(integrand, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      return result

   def fTshotNoise(self, l):
      """shot noise for CIB 2d trispectrum
      """
      integrand = lambda a: (3.e5/self.U.hubble(1./a-1.)) / a**2 *\
                           self.U.bg.comoving_distance(1./a-1.)**2 * self.jNu4(a)
      result = integrate.quad(integrand, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      return result
   
   def plotFig1Penin14(self):
      """Successfully reproduces fig 1 in Penin+14
      their j_nu quantity is weirdly defined,
      with a useless factor of a:
      their j_nu = my j_nu * chi^2 / a
      without the factor of a, this would be the relevant projection kernel for CIB
      """
      Chi = np.array(map(lambda a: self.U.bg.comoving_distance(1./a-1.), self.A))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(self.Z, self.JNu1/self.A*Chi**2, 'b', lw=2, label=str(int(self.nu/1.e9))+' GHz')
      #
      ax.legend(loc=1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'their $j_\nu=$ my $j_\nu \chi^2 / a = $ proj. kernel$/a$')
      #
      plt.show()
   
   
   def plotFig3Bethermin12(self):
      """Sort of successfully reproduces fig 3 in Bethermin+12
      they plot a meaningless quantity related to counts,
      presumably to reduce the span of the y-axis...
      this test is not too important, since I can reproduce fig 1 in Penin+14
      """
      Quantity = np.trapz(self.dNdSnudzdOmega * self.Snu[np.newaxis, :]**(5./2.), self.Z, axis=0)
      Quantity = np.abs(Quantity)
      #print Quantity
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.loglog(self.Snu*1.e3, Quantity)
      #
      ax.set_xlabel(r'$S_\nu$ [mJy]')
      ax.set_ylabel(r'$S_\nu^{5/2} dN / d\Omega dS_\nu$')
      #
      plt.show()
   
   
   def plotdNdSnudV(self):
      Z = np.linspace(0., 5., 6)
      A = 1./(1.+Z)
      IA = np.array([ np.argmin((self.A-a)**2) for a in A ])
      
      # luminosity functions
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         iA = IA[iZ]
         z = Z[iZ]
         ax.plot(self.Snu, self.Snu * self.dNdSnudV[iA, :], lw=2, label=r'$z=$'+str(int(z)))
      #
      ax.legend(loc=3)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$S_\nu$ at $\nu=$'+str(np.int(self.nu/1.e9))+'GHz, [Jy]')
      ax.set_ylabel(r'$dN/d\text{ln}S_{\nu}/dV$ [gal/(Mpc/h)$^3$]')
      #
      #fig.savefig("./figures/cib_penin12/bethermin12_dNdlnSnudV.pdf", bbox_inches='tight')
      
      # shot noise power spectrum density
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         iA = IA[iZ]
         z = Z[iZ]
         ax.plot(self.Snu, self.Snu**2 * self.Snu * self.dNdSnudV[iA, :], lw=2, label=r'$z=$'+str(int(z)))
      #
      ax.legend(loc=3)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$S_\nu$ at $\nu=$'+str(np.int(self.nu/1.e9))+'GHz, [Jy]')
      ax.set_ylabel(r'$\frac{dP^\text{shot}}{d\text{ln}S_{\nu}} = S_{\nu}^2 \; dN/d\text{ln}S_{\nu}/dV$ [gal/(Mpc/h)$^3$]')
      #
      #fig.savefig("./figures/cib_penin12/bethermin12_Snu2dNdlnSnudV.pdf", bbox_inches='tight')
      
      # shot noise trispectrum density
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         iA = IA[iZ]
         z = Z[iZ]
         ax.plot(self.Snu, self.Snu**4 * self.Snu * self.dNdSnudV[iA, :], lw=2, label=r'$z=$'+str(int(z)))
      #
      ax.legend(loc=3)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$S_\nu$ at $\nu=$'+str(np.int(self.nu/1.e9))+'GHz, [Jy]')
      ax.set_ylabel(r'$\frac{d\mathcal{T}^\text{shot}}{d\text{ln}S_{\nu}} = S_{\nu}^4 \; dN/d\text{ln}S_{\nu}/dV$ [gal/(Mpc/h)$^3$]')
      #
      #fig.savefig("./figures/cib_penin12/bethermin12_Snu4dNdlnSnudV.pdf", bbox_inches='tight')
      
      plt.show()


   def plotJnu1(self):
      Chi = np.array(map(lambda a: self.U.bg.comoving_distance(1./a-1.), self.A))
      # inverse hubble length: H/c in (h Mpc^-1)
      # used to convert the kernel from chi to z
      H = self.U.hubble(1./self.A-1.) / 3.e5
      # projection kernel
      W = np.array(map(self.f, self.A))
      # factor to keep jnu in the plot
      factorJNu1 = 2.*self.JNu1[np.argmin(np.abs(self.Z-2.))]

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(self.Z, (Chi**2/H) / np.max(Chi**2/H), 'g', lw=1.5, label=r'$dV / dzd\Omega = \chi^2 c /H$')
      ax.plot(self.Z, (self.JNu1) / factorJNu1, 'r', lw=1.5, label=r'$\bar{j}_\nu = \int dS_\nu \frac{dN_\text{gal}}{dVdS_\nu} S_\nu$')
      ax.plot(self.Z, (W/H) / np.max(W/H), 'b', lw=3, label=r'$W_z =  \bar{j}_\nu \times dV/dzd\Omega$')
      #
      ax.legend(loc=1)
      ax.set_ylim((0., 1.1))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'arbitrary units')
      #ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      #
      #path="./figures/cib_penin12/emissivity_kernel"+str(int(self.nu/1.e9))+"GHz_penin12.pdf"
      #fig.savefig(path, bbox_inches='tight')
      
      plt.show()







