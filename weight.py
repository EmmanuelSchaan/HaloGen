from headers import *
##################################################################################

class Projection(object):
   
   # Required variables: a_min, a_max, name
   
   def __str__(self):
      return self.name
   
   # projection kernel, such that, e.g.:
   # kappa = int dchi f(a) delta
   def f(self, a):
      pass

   def __init__(self, U, name=''):
      # copy U
      self.U = U
      self.name=name

   ##################################################################################

   def dWddelta(self, aMin, aMax):
      """response of the projected quantity (kappa, y2d, ...)
      to the mean overdensity along the line of sight,
      measured between aMin and aMax
      """
      integrand = lambda a: 3.e5/(self.U.Hubble(a) * a**2) * self.f(a)
      result = integrate.quad(integrand, aMin, aMax, epsabs=0., epsrel=1.e-3)[0]
      return result

   ##################################################################################
   
   def plotW(self):
      
      # choice of weight function
      n = 2 # n of n-pt function
      nW = n
      nd = 2*n-2.
      
      # range to plot
      Na = 101
      A = np.linspace(self.aMin, self.aMax, Na)
      ComovDistToObs = np.array( map( lambda a: self.U.ComovDist(a, self.U.a_obs), A ) )
      W = np.array( map( lambda a: self.f(a), A ) )
      Z = 1./A - 1.
      H_A = self.U.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
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
         #print 'l=', l, 'k=', l/self.U.ComovDist(1./(1.+0.5), self.U.a_obs)
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
         #print 'l=', l, 'k=', l/self.U.ComovDist(1./(1.+0.5), self.U.a_obs)
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
      super(WeightY, self).__init__(U, name=name)
      #
      self.aMin = 0.2   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
   
   def f(self, a):
      """Compton y projection kernel
      """
      return a


##################################################################################
##################################################################################

class WeightLensSingle(Projection):
   """Lensing projection: single source. The default is z_source = 1.
   """
   
   def __init__(self, U, z_source=1., name='lens'):
      super(WeightLensSingle, self).__init__(U, name=name)
      # a for mass func, biases, and projection
      a_source = 1./(1.+z_source)
      #
      self.z_source = z_source
      self.dist_source = self.U.ComovDist(a_source, self.U.a_obs)
      self.aMin = max(a_source, 1./11.)   # don't go further than z=10
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)
   
   def f(self, a):
      """lensing projection kernel
      for single source
      a is dimless, Wlensing(a) in (h Mpc^-1)
      """
      d_a = self.U.ComovDist(a, self.U.a_obs)
      wlensing = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      wlensing *= 1. - d_a/self.dist_source
      return wlensing
   
   def testHandEtAl13_fig1(self):
      # fig 1 from Hand et al 2013
      Data = np.genfromtxt('./input/tests/HandEtAl13/HandEtAl13_fig1_cmb.txt')
      Z = Data[:, 0]
      A = 1./(1.+Z)
      Wgal_ref = Data[:, 1]
      # from my code
      Zme = np.linspace(0., 10., 101)
      Ame = 1./(1.+Zme)
      fW = lambda a: self.f(a) * (3.e5/self.U.Hubble(a))
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
      super(WeightLensOguriTakada11, self).__init__(U, name=name)
      #
      # parameters for source distribution (Oguri & Takada 2011)
      self.z0 = z_source/3. # should be <z_source>/3
      self.nz0 = 20.   # width of integral over sources
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + self.nz0 * self.z0 ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)


   def fdpdz(self, z):
      """source distribution from Oguri & Takada 2011
      int_0^inf dz dPdz = 1
      """
      result = 0.5 * z**2 / self.z0**3
      result *= exp( -z / self.z0)
      return result
   
   
   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result


##################################################################################
##################################################################################

class WeightLensHandEtAl13(Projection):
   """lensing projection: source distribution from Hand et al 2013
   """
   
   def __init__(self, U, name='lens'):
      super(WeightLensHandEtAl13, self).__init__(U, name=name)
      
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
      self.aMax = self.U.a_obs*(1.-epsilon)

   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def testHandEtAl13_fig1(self):
      # fig 1 from Hand et al 2013
      Data = np.genfromtxt('./input/tests/HandEtAl13/HandEtAl13_fig1_gal.txt')
      Z = Data[:, 0]
      A = 1./(1.+Z)
      Wgal_ref = Data[:, 1]
      # from my code
      fW = lambda a: self.f(a) * (3.e5/self.U.Hubble(a))
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
      Z = np.linspace(0., 3., 101)
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
      super(WeightLensDasEtAl13, self).__init__(U, name=name)
      
      # source distribution (eq11 from Das Errard Spergel 2013)
      z0 = 0.69   # ie median z is 1
      self.fdpdz = lambda z: 1.5 * z**2/z0**3 * np.exp(-(z/z0)**1.5)
      
      # a for mass func, biases, and projection
      self.aMin = 1./(1.+10.)  # arbitrary for now
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)
   
   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def plot(self):
      Na = 101
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
   
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
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
   
   def __init__(self, U, fdndz, zMin=1.e-4, zMax=2., name='lens'):
      super(WeightLensCustom, self).__init__(U, name=name)
      
      self.aMin = 1./(1.+zMax)
      self.aMax = 1./(1.+zMin)
      
      # fdndz doesn't need to be normalized to anything
      self.fdndz = fdndz
      self.ngal = integrate.quad(self.fdndz, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-3)[0]
      # dpdz normalized such that int dz dpdz = 1
      self.fdpdz = lambda z: fdndz(z) / self.ngal

   
   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def plot(self):
      Na = 101
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
   
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
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
      super(WeightLensCIBSchmidt15, self).__init__(U, name=name)
      #
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + 10. ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)
      #
      # source distribution from Schmidt Menard Scranton+15
      fdpdzNonNormalized = lambda z: z**alpha * np.exp( -(z / z0)**alpha)
      # normalize to have int_zMin^zMax dz dPdz = 1
      norm = integrate.quad(fdpdzNonNormalized, 1./self.aMax, 1./self.aMin, epsabs=0, epsrel=1.e-4)[0]
      self.fdpdz = lambda z: fdpdzNonNormalized(z) / norm
   
   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result
   

   def plot(self):
      Na = 101
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
      
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
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
      super(WeightLensCIBPullen17, self).__init__(U, name=name)
      #
      # a for mass func, biases, and projection
      self.aMin = 1./( 1. + 5. ) # integrate far enough to get all the sources
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)
      #
      # read digitized values
      path = "./input/cib_zdist_pullen17/Pullen+17_"+str(nu)+".txt"
      data = np.genfromtxt(path)
      fdpdzNonNormalized = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)
      # normalize to have int_zMin^zMax dz dPdz = 1
      norm = integrate.quad(fdpdzNonNormalized, 1./self.aMax, 1./self.aMin, epsabs=0, epsrel=1.e-4)[0]
      self.fdpdz = lambda z: fdpdzNonNormalized(z) / norm
      
   
   def f(self, a):
      d_a = self.U.ComovDist(a, self.U.a_obs)
      result = 1.5 * (100./3.e5)**2 * self.U.OmM * d_a / a
      integrand = lambda a_s: self.fdpdz(1./a_s-1.) /a_s**2 * (1. - d_a/self.U.ComovDist(a_s, self.U.a_obs))
      result *= integrate.quad(integrand, self.aMin, a, epsabs=0, epsrel=1.e-2)[0]
      return result

   def plot(self):
      Na = 101
      A = np.linspace(self.aMin, self.aMax, Na)
      Z = 1./A - 1.
      
      # source distribution
      dPdz = np.array(map(self.fdpdz, Z))
      
      # lensing kernel
      W = np.array( map( lambda a: self.f(a), A ) )
      H_A = self.U.Hubble(A) / 3.e5   # inverse hubble length: H/c in (h Mpc^-1)
      
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
      super(WeightCIBPlanck15, self).__init__(U, name=name)
   
      # a for mass func, biases, and projection
      self.aMin = 1./(1.+10.)  # arbitrary for now
      epsilon = 1.e-5
      self.aMax = self.U.a_obs*(1.-epsilon)
   
   def f(self, a):
      return self.U.ComovDist(a, self.U.a_obs)**2


##################################################################################
##################################################################################

class WeightTracer(Projection):
   """Projected density field of tracers
   requires defining b(z) and dn/dz(z)
   """
   
   def __init__(self, U, name='d'):
      super(WeightTracer, self).__init__(U, name=name)

      # normalization of dn/dz, ie number of galaxies per unit steradian
      self.ngal = integrate.quad(self.dndz, 1./self.aMax-1., 1./self.aMin-1., epsabs=0., epsrel=1.e-3)[0]
      # convert to number of galaxies per square arcmin
      self.ngal_per_arcmin2 = self.ngal * (np.pi/180./60.)**2


   def b(self, z):
      """tracer bias
      """
      pass
   
   def dndz(z):
      """normalized such that int dz dn/dz = ngal,
      ie the number of gals per unit steradian
      """
      pass

   def f(self, a):
      """projection kernel
      """
      z = 1./a - 1.
      result = self.U.Hubble(a) / 3.e5
      result *= self.dndz(z)
      result /= self.ngal
      result *= self.b(z)
      return result

   def plotDndz(self):

      Z = np.linspace(1./self.aMax-1., 1./self.aMin-1., 101)
      Dndz = np.array(map(self.dndz, Z))
      # normalize such that int dz dn/dz = ngal in arcmin^-2
      Dndz /= (180.*60. / np.pi)**2

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, Dndz)
      #
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$dn/dz$ [arcmin$^{-2}$]')
      
      plt.show()


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


##################################################################################
##################################################################################

class WeightTracerLSSTSources(WeightTracer):
   """Projected number density of LSST source galaxies,
   used for shear measurements,
   as in Schaan Krause Eifler +16.
   """
   
   def __init__(self, U, name='lsstsources'):
      self.aMin = 1./(1.+4.)   # min bound for integral over a
      self.aMax = 1.-0.005   # max bound for integral over a
      self.zMin = 1./self.aMax-1.
      self.zMax = 1./self.aMin-1.
      
      # tracer bias
      # copied from the LSST gold sample!
      self.b = lambda z: 1 + 0.84*z
      
      self.ngal_per_arcmin2 = 26. # galaxies per squared arcmin
      self.ngal = self.ngal_per_arcmin2 / (np.pi/180./60.)**2  # per steradian
      

      # dn/dz, non-normalized
      self.z0 = 0.5
      self.alpha = 1.27
      self.beta = 1.02
      f = lambda z: z**self.alpha * np.exp(-(z/self.z0)**self.beta)
      # normalization
      norm = integrate.quad(f, self.zMin, self.zMax, epsabs=0., epsrel=1.e-3)[0]
      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      self.dndz = lambda z: self.ngal * f(z) / norm

      super(WeightTracerLSSTSources, self).__init__(U, name=name)




##################################################################################
##################################################################################

class WeightTracerDESIQSO(WeightTracer):
   """Projected number density of CMASS DR12 galaxies
   """

   def __init__(self, U, name='desiqso'):
      self.aMin = 1./(1.+2.1)   # min bound for integral over a
      self.aMax = 1./(1.+0.9)   # max bound for integral over a

      # tracer bias
      self.b = lambda z: 2.

      # dn/dz, normalized such that int dz dn/dz = ngal
      # where ngal = number of gals per unit steradian
      ngal = 180. # per deg^2
      ngal /= (np.pi/180.)**2
      self.dndz = lambda z: ngal * (z>=0.9) * (z<=2.1) / (2.1-0.9)

      super(WeightTracerDESIQSO, self).__init__(U, name=name)


##################################################################################
##################################################################################

class WeightTracerCustom(Projection):
   """Projected density field of tracers
   requires defining b(z) and dn/dz(z)
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
      super(WeightCIBPenin12, self).__init__(U, name=name+'_'+str(int(nu/1.e9))+'GHZ')
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
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), self.A))
      Hubble = np.array(map(lambda a: self.U.Hubble(a), self.A))
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
                                    4.*np.pi*self.U.ComovDist(a, 1.)**2* 3.e5/self.U.Hubble(a)
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

   def f(self, a):
      """projection kernel
      """
      result = self.U.ComovDist(a, 1.)**2
      result *= self.jNu1(a)
      return result
   
   def fdPshotNoise_da(self, a, l):
      """contribution of each scale factor to the shot noise
      """
      result = (3.e5/self.U.Hubble(a)) / a**2
      result *= self.U.ComovDist(a, 1.)**2 * self.jNu2(a)
      return result
   
   def fPshotNoise(self, l):
      """Shot noise for CIB 2d power spectrum
      """
      integrand = lambda a: (3.e5/self.U.Hubble(a)) / a**2 *\
                           self.U.ComovDist(a, 1.)**2 * self.jNu2(a)
      result = integrate.quad(integrand, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      return result

   def fTshotNoise(self, l):
      """shot noise for CIB 2d trispectrum
      """
      integrand = lambda a: (3.e5/self.U.Hubble(a)) / a**2 *\
                           self.U.ComovDist(a, 1.)**2 * self.jNu4(a)
      result = integrate.quad(integrand, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      return result
   
   def plotFig1Penin14(self):
      """Successfully reproduces fig 1 in Penin+14
      their j_nu quantity is weirdly defined,
      with a useless factor of a:
      their j_nu = my j_nu * chi^2 / a
      without the factor of a, this would be the relevant projection kernel for CIB
      """
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), self.A))
      
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
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), self.A))
      # inverse hubble length: H/c in (h Mpc^-1)
      # used to convert the kernel from chi to z
      H = self.U.Hubble(self.A) / 3.e5
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







