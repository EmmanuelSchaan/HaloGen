from headers import *
##################################################################################

class MassFunction(object):

   def __init__(self, U, save=False):
   
      # required variables: path_MassFunc, path_B1, path_B2
      
      # masses to compute mass function and biases
      self.mMin = 1.e10  # in (h^-1 solarM)
      self.mMax = 1.e16  # in (h^-1 solarM)
      self.Nm = 201 # nb of m pts
      self.M = np.logspace(np.log10(self.mMin), np.log10(self.mMax), self.Nm, 10.) # masses in h^-1 solarM
      
      # values of a to be computed
      self.aMin = 0.08  #1./(1.+10.)
      self.aMax = 1.
      self.Na = 101
      self.A = np.linspace(self.aMin, self.aMax, self.Na)
   
      # copy U and Proj
      self.U = U
      
      # compute things if necessary
      if save==True:
         self.Save()
      
      self.Load()
   
   ##################################################################################
   
   def massfunc(self, m, z):
      pass
   
   def b(m, z):
      pass
 
   def Save(self):
      # mass function and biases
      print "Computing mass function and biases"
      Massfunc = np.zeros((self.Nm, self.Na))
      B1 = np.zeros((self.Nm, self.Na))
      B2 = np.zeros((self.Nm, self.Na))
      
      for ia in range(self.Na):
         z = 1./self.A[ia] - 1.
         for im in range(self.Nm):
            m = self.M[im]
            # make sure not to get 0
            Massfunc[im, ia] = max( self.massfunc(m, z), 1.e-300 )
            B1[im, ia], B2[im, ia] = self.b(m, z)
         print "- done "+str(ia+1)+" of "+str(self.Na)
      
      np.savetxt(self.path_Massfunc, Massfunc)
      np.savetxt(self.path_B1, B1)
      np.savetxt(self.path_B2, B2)
      return
   
   def Load(self):
      # load sigma, nu, mass function and biases
      print "Loading mass function and biases"
      self.Massfunc = np.genfromtxt(self.path_Massfunc)
      self.B1 = np.genfromtxt(self.path_B1)
      self.B2 = np.genfromtxt(self.path_B2)
      
      # interpolate
      interp_massfunc = RectBivariateSpline(np.log(self.M), self.A, np.log(self.Massfunc), s=0)
      self.fmassfunc = lambda m, a: (a>=self.aMin and a<=self.aMax) * np.exp( interp_massfunc(np.log(m), a)[0,0] )
      #
      interp_b1 = RectBivariateSpline(np.log(self.M), self.A, np.log(self.B1), s=0)
      self.fb1 = lambda m, a: (a>=self.aMin and a<=self.aMax) * np.exp( interp_b1(np.log(m), a)[0,0] )
      #
      interp_b2 = RectBivariateSpline(np.log(self.M), self.A,  self.B2, s=0)
      self.fb2 = lambda m, a: (a>=self.aMin and a<=self.aMax) * interp_b2(np.log(m), a)[0,0]


   def testInterp(self, z=0.):
      
      # closest computed redshift
      a = 1./(1.+z)
      iaBest = np.where(abs(self.A-a)==np.min(abs(self.A-a)))[0][0]
      
      # interpolated values
      M = np.logspace(np.log10(self.mMin), np.log10(self.mMax), 5*self.Nm, 10.)
      #
      f = lambda m: self.fmassfunc(m, a)
      RecMassfunc = np.array(map(f, self.M))
      #
      f = lambda m: self.fb1(m, a)
      RecB1 = np.array(map(f, self.M))
      #
      f = lambda m: self.fb2(m, a)
      RecB2 = np.array(map(f, self.M))
      
      # mass function
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.loglog(self.M, self.M * self.Massfunc[:, iaBest], 'r-', lw=3, label=r'closest computed: $z=$'+str(round(1./self.A[iaBest]-1.,1)))
      for ia in range(self.Na):
         ax.loglog(self.M, self.M * self.Massfunc[:, ia], 'r-', alpha=0.2)
      ax.loglog(self.M, self.M * RecMassfunc, 'b', lw=3, label=r'interpolated: $z=$'+str(round(z,1)))
      #
      ax.legend(loc=3, numpoints=1)
      ax.set_xlabel(r'M [M$_\odot$/h]', fontsize=18)
      ax.set_ylabel(r'$dn/d\ln(m)$', fontsize=18)
      
      # b1
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.loglog(self.M, self.B1[:, iaBest], 'r-', lw=3, label=r'closest computed: $z=$'+str(round(1./self.A[iaBest]-1.,1)))
      for ia in range(self.Na):
         ax.loglog(self.M, self.B1[:, ia], 'r-', alpha=0.2)
      ax.loglog(self.M, RecB1, 'b', lw=3, label=r'interpolated: $z=$'+str(round(z,1)))
      #
      ax.legend(loc=3, numpoints=1)
      ax.set_xlabel(r'M [M$_\odot$/h]', fontsize=18)
      ax.set_ylabel(r'$b_1(m)$', fontsize=18)
      
      # b2
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      ax.semilogx(self.M, self.B2[:, iaBest], 'r-', lw=3, label=r'closest computed: $z=$'+str(round(1./self.A[iaBest]-1.,1)))
      for ia in range(self.Na):
         ax.semilogx(self.M, self.B2[:, ia], 'r-', alpha=0.2)
      ax.semilogx(self.M, RecB2, 'b', lw=3, label=r'interpolated: $z=$'+str(round(z,1)))
      #
      ax.legend(loc=3, numpoints=1)
      ax.set_xlabel(r'M [M$_\odot$/h]', fontsize=18)
      ax.set_ylabel(r'$b_2(m)$', fontsize=18)

      plt.show()


##################################################################################


class MassFuncPS(MassFunction):

   def f(self, nu, z):
      """returns f(nu) for the Sheth-Tormen mass function
      """
      f = np.exp(-nu/2.)
      f /= np.sqrt(2.*np.pi*nu)
      return f
   
   def b(self, m, z):
      """only b1 implemented
      """
      nu = self.U.fnu(m, z)
      return 1 + (nu-1) / self.U.deltaC_z(z), 0.
   
   def massfunc(self, m, z):
      """dn/dm
      """
      nu = self.U.fnu(m, z)
      f = self.U.rho_z(z)/m**2
      f *= nu * self.f(nu, z)
      f *= self.U.fdlnnu_dlnm(m, z)
      return f
   
   def __init__(self, U, save=False):
      # paths for saving/loading data
      self.path_Massfunc = "./output/massfunc/PS_Massfunc.txt"
      self.path_B1 = "./output/massfunc/PS_B1.txt"
      self.path_B2 = "./output/massfunc/PS_B2.txt"
      #
      super(MassFuncPS, self).__init__(U, save=save)



##################################################################################

class MassFuncST(MassFunction):
   
   def f(self, nu, z):
      """returns f(nu) for the Sheth-Tormen mass function
      """
      f = self.A_p
      f *= 1 + (self.q*nu)**(-self.p)
      f *= np.sqrt( self.q*nu / (2.*np.pi) )
      f *= exp( -0.5 * self.q*nu )
      f /= nu  # to get f(nu) and not nu*f(nu)
      return f

   def b(self, m, z):
      """returns bias1 and bias2 for the Sheth-Tormen mass function
      """
      nu = self.U.fnu(m, z)
      dc = self.U.deltaC_z(z)
      #
      b1 = 1. + (self.q*nu - 1.)/dc
      b1 += 2.*self.p / ( dc * (1. + (self.q*nu)**self.p) )
      #
      b2 = 8./21. * (b1 - 1.)
      b2 += nu * (nu - 3.) / dc**2
      
      return b1, b2

   def massfunc(self, m, z):
      """dn/dm
      """
      nu = self.U.fnu(m, z)
      f = self.U.rho_z(z)/m**2
      f *= nu * self.f(nu, z)
      f *= self.U.fdlnnu_dlnm(m, z)
      return f

   def __init__(self, U, save=False):
      # params for Sheth-Tormen
      # different from Sheth & Tormen 1999, q = 0.707
      self.p = .3
      #self.q = .707 # from Takada & Jain 2002
      self.q = .75  # from Takada & Jain 2003
      # amplitude, function of p only, to have "int f(nu) dnu = 1"
      self.A_p = 1./( 1. + special.gamma(0.5-self.p)/(2**self.p * np.sqrt(np.pi)) )

      # paths for saving/loading data
      self.path_Massfunc = "./output/massfunc/ST_Massfunc.txt"
      self.path_B1 = "./output/massfunc/ST_B1.txt"
      self.path_B2 = "./output/massfunc/ST_B2.txt"
      #
      super(MassFuncST, self).__init__(U, save=save)


##################################################################################


class MassFuncTinker(MassFunction):
   
   def __init__(self, U, save=False):
      # copy U and Proj
      self.U = U
      
      # params for Tinker
      self.A0 = 0.186
      self.Az = -0.14
      self.a0 = 1.47
      self.az = -0.06
      self.b0 = 2.57
      #self.bz = - np.exp(  - (0.75/np.log(200./0.75))**1.2  )
      self.bz = - 10.**(  - (0.75/np.log10(200./0.75))**1.2  )
      self.c0 = 1.19
      # concentration params from Duffy et al 2008 for Tinker
      self.cNFW0 = 5.71
      self.cNFWam = -0.084
      self.cNFWaz = -0.47
      # paths for saving/loading data
      self.path_Massfunc = "./output/massfunc/Tinker_Massfunc.txt"
      self.path_B1 = "./output/massfunc/Tinker_B1.txt"
      self.path_B2 = "./output/massfunc/Tinker_B2.txt"
      #
      super(MassFuncTinker, self).__init__(U, save=save)


   def f(self, sigma, z):
      """returns f(sigma) for the Tinker mass function (Tinker et al 2008)
      notice that here m stands for M200d (the jacobian is included in Save())
      """
      z = min(z, 2.5)   # recommended by Tinker et al 2008
      A = self.A0 * (1.+z)**self.Az
      a = self.a0 * (1.+z)**self.az
      b = self.b0 * (1.+z)**self.bz
      c = self.c0
      # get sigma from our convention for nu
      f = A * ( (sigma/b)**(-a) + 1. ) * np.exp(-c/sigma**2)
      return f
   
   def Tinker_massfunc_m200d(self, m, z):
      """dn/dm200d; input m=m200d here
      """
      r = (3.*m / (4.*np.pi*self.U.rho_z(z)))**(1./3.)
      sigma = np.sqrt( self.U.Sigma2(r, z, W3d_sth) )
      #
      f = self.f(sigma, z)
      f *= self.U.rho_z(z)/m**2
      f *= abs(self.U.fdlnSigma_dlnM(m, z))
      return f
   
   def massfunc(self, m, z):
      """dn/dm, including the jacobian dm200_d/dm; m=m_vir here
      """
      m200d = self.U.massRadiusConversion(m, z, 200., 'm')[0]
      
      f = self.Tinker_massfunc_m200d(m200d, z)
      #f = self.Tinker_massfunc_m200d(m, z)
      
      # jacobian for mass conversion
      epsilon = 1.e-3
      dm200_dm = self.U.massRadiusConversion(m*(1.+epsilon), z, 200., 'm')[0] - self.U.massRadiusConversion(m*(1.-epsilon), z, 200., 'm')[0]
      dm200_dm /= 2.*m*epsilon
      
      f *= dm200_dm
      return f
   
   def Tinker_b1_nu(self, nu, z):
      """Tinker's convention: nu = dc/sigma, with dc=1.686
      """
      # params for Tinker bias
      dc = 1.686  # from Tinker et al 2010
      y = np.log10(200.)
      A = 1. + 0.24*y*np.exp( -(4./y)**4 )
      a = 0.44*y - 0.88
      B = 0.183
      b = 1.5
      C = 0.019 + 0.107*y + 0.19*np.exp( -(4./y)**4 )
      c = 2.4
      # first order bias
      b1 = 1. - A * 1./(1. + (dc/nu)**a) + B*nu**b + C*nu**c
      return b1
   
   
   def b(self, m, z):
      # Tinker's convention: nu = dc/sigma
      R = (3.*m / (4.*np.pi*self.U.rho_z(z)))**(1./3.)
      dc = 1.686  # from Tinker et al 2010
      nu = dc / np.sqrt( self.U.Sigma2(R, z, W3d_sth) )
      # first order bias
      b1 = self.Tinker_b1_nu(nu, z)
      # second order bias
      b2 = 0.
      return b1, b2


   def plotTinkerMassFunc(self):
      """reproduces fig6 and fig5 in Tinker et al 2008
      requires massfunc_id = "ST"
      """
      z = 0.
      
      # fig6, from my code
      Sigma = np.logspace(-0.4, 0.6, 51, 10.)
      f = lambda sigma: self.f(sigma, z)
      Tinker_F = np.array(map(f, Sigma))
      # fig6, from Tinker et al 2008
      Test = np.genfromtxt('./input/tests/Tinkeretal08_fig6.txt')
      X = Test[:,0]
      Y = Test[:, 1]
      #
      fig6 = plt.figure(6)
      ax = plt.subplot(111)
      ax.loglog(Sigma, Tinker_F, 'g', label='my code')
      ax.loglog(X, Y, 'r', label='Tinker et al 2008')
      ax.grid()
      ax.set_xlabel(r'$\sigma$')
      ax.set_ylabel(r'$f(\sigma)$')
      ax.set_title('fig6 from Tinker et al 2008')
      ax.legend(loc=4)
      #fig6.savefig("./figures/tests/Tinkeretal08_fig6.pdf")
      
      # fig5 from my code
      M = np.logspace(np.log10(1.e10), np.log10(1.e16), 51, 10.) # masses in h^-1 solarM
      f = lambda m: self.Tinker_massfunc_m200d(m, z)
      Tinker_massfunc = np.array(map(f, M))
      # fig5 from Tinker et al 2008
      Test = np.genfromtxt('./input/tests/Tinkeretal08_fig5.txt')
      X = Test[:,0]
      Y = Test[:, 1]
      #
      fig5 = plt.figure(5)
      ax = plt.subplot(111)
      ax.loglog(M, M**2/self.U.rho_z(z) * Tinker_massfunc, 'g', label='my code')
      ax.loglog(X, Y, 'r', label='Tinker et al 2008')
      ax.grid()
      ax.set_xlabel(r'$M_{200,d}$ [$M_{sun}$/h]')
      ax.set_ylabel(r'$(m^2/\rho)dn/dm$')
      ax.set_title('fig5 from Tinker et al 2008')
      ax.legend(loc=4)
      #fig5.savefig("./figures/tests/Tinkeretal08_fig5.pdf")
      
      plt.show()


   def plotTinkerBias(self):
      """reproduces fig1 in Tinker et al 2010
      """
      # from my code
      TinkerNu = np.logspace(-0.5, 0.6, 101, 10.) # log in base 10
      z = 0.
      f = lambda nu: self.Tinker_b1_nu(nu, z)
      B = np.array(map( f, TinkerNu ))
      # from Tinker et al 2010, fig1
      Test = np.genfromtxt('./input/tests/Tinkeretal10_fig1.txt')
      X = Test[:,0]
      Y = Test[:, 1]
      #
      fig=plt.figure(0)
      ax=plt.subplot(111)
      ax.semilogx(TinkerNu, B, 'g', label='my code')
      ax.semilogx(X, Y, 'r', label='Tinker et al 2010')
      ax.grid()
      ax.set_xlim((10**(-0.5), 10**(0.6)))
      ax.set_ylim((0., 9.))
      ax.set_title('fig1 from Tinker et al 2010')
      ax.set_xlabel(r'Tinker $\nu$')
      ax.set_ylabel(r'$b(\nu)$')
      #fig.savefig("./figures/tests/Tinkeretal10_fig1.pdf")
      
      plt.show()


   def plotMasses(self):
      # Mvir and redshifts
      Mvir = self.M  # in Msun/h
      Z = np.linspace(0., 3., 5)
      
      # M200d
      plt.figure(0)
      ax=plt.subplot(111)
      for iz in range(len(Z)):
         z = Z[iz]
         fm200d = lambda m: self.U.massRadiusConversion(m, z, 200., 'm')[0]
         M200d = np.array(map(fm200d, Mvir))
         ax.semilogx(Mvir, M200d/Mvir, label=r'$z=$'+str(z))
      ax.set_xlabel('$M_{vir}$ $[M_{sun}/h]$')
      ax.set_ylabel('$M_{200d}/M_{vir}$')
      ax.legend(loc=4)
      
      # M200c
      plt.figure(1)
      ax=plt.subplot(111)
      for iz in range(len(Z)):
         z = Z[iz]
         fm200c = lambda m: self.U.massRadiusConversion(m, z, 200, 'c')[0]
         M200c = np.array(map(fm200c, Mvir))
         ax.semilogx(Mvir, M200c/Mvir, label=r'$z=$'+str(z))
      ax.set_xlabel('$M_{vir}$ $[M_{sun}/h]$')
      ax.set_ylabel('$M_{200c}/M_{vir}$')
      ax.legend(loc=4)
      
      plt.show()



   def compareMassesColin(self):
      # fixed redshift
      z = 0.2
      # data from Colin
      data = np.genfromtxt("./input/tests/Colin_tSZ/mvir_m200_duffy_z0.2.txt")
      Mvir = data[:, 0] # [Msun/h]
      M200c_Colin = data[:, 1]   # [Msun/h]
      # my calculation
      f = lambda m: self.U.massRadiusConversion(m, z, 200., 'c')[0]
      M200c = np.array(map(f, Mvir))
      #
      fig=plt.figure(0)
      ax=plt.subplot(211)
      ax.semilogx(Mvir, M200c_Colin/Mvir, 'r', label='Colin')
      ax.semilogx(Mvir, M200c/Mvir, 'b', label='me')
      ax.legend(loc=3)
      ax.set_ylabel(r'$M_{200c}/M_{vir}$')
      ax.grid()
      #
      ax=plt.subplot(212, sharex=ax)
      ax.semilogx(Mvir, M200c/M200c_Colin - 1., 'b')
      ax.grid()
      ax.set_xlabel(r'$M_{vir}$ [$M_{sun}$/h]')
      ax.set_ylabel(r'me/Colin -1')
      
      #fig.savefig("./figures/tests/m200c.pdf")
      
      plt.show()
