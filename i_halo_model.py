from headers import *
##################################################################################

class IHaloModel():

   def __init__(self, U, MassFunc, save=0):
      # copy classes
      self.U = U
      self.MassFunc = MassFunc

   ##################################################################################

   def f(self, i, Profiles, z, nbias=0, mMin=None, test=False):
      """integrals Mij, and their modified versions when nbias is included
      from theory: nbias = 0 or 1
      K is the size j array of moduli [k1, ..., kj]
      """
      # find closest redshift to z in Data.A
      a = 1./(z+1.)
      
      # mass function
      massfunc = lambda m: self.MassFunc.fmassfunc(m, a)
      # bias if necessary
      if i==1:
         bias = lambda m: self.MassFunc.fb1(m, a)
      elif i==2:
         bias = lambda m: self.MassFunc.fb2(m, a)
      else :
         bias = lambda m: 1.
      # extra bias if necessary
      if nbias<>0:
         extrabias = lambda m: self.MassFunc.fb1(m, a)**nbias
      else :
         extrabias = lambda m: 1.
      # profiles
      profiles = lambda m: reduce(lambda x,y: x*y, [ Profiles[i][0].uN(Profiles[i][2], m, z, n=Profiles[i][1]) for i in range(len(Profiles)) ])
      
      # integrand in lnm, for speed
      def integrand(lnm):
         m = np.exp(lnm)
         result = massfunc(m) * profiles(m) * extrabias(m) * bias(m)
         result *= m # because integrating in lnm and not m
         return result
      
      if test:
         # plot integrand to check for convergence
         M = self.MassFunc.M.copy()
         Integrand = np.array(map(integrand, np.log(M)))
         #
         #
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(M, Integrand, lw=2, label=r'$z=$'+str(z))
         #
         ax.legend(loc=1)
         ax.set_xscale('log')
         ax.set_xlabel(r'halo mass $M$ [$M_\odot/h$]')
         ax.set_ylabel(r'$\frac{d I}{d \ln M}$')
         #
         plt.show()

      # integration bounds
      mMinIntegral = np.max([ Profiles[j][0].mMin for j in range(len(Profiles)) ])
      mMinIntegral = np.max([mMinIntegral, mMin, self.MassFunc.mMin])
      mMaxIntegral = np.min([ Profiles[j][0].mMax for j in range(len(Profiles)) ])
      mMaxIntegral = np.min([mMaxIntegral, self.MassFunc.mMax])
      # compute integral
      integral = integrate.quad(integrand, np.log(mMinIntegral), np.log(mMaxIntegral), epsabs=0, epsrel=1.e-2)[0]

      # correction factor?
      if i==1 and len(Profiles)==1 and Profiles[0][0].use_correction_factor==1:
         integral /= self.correction_factor_m11(a)
      return integral


   ##################################################################################

   def correction_factor_m11(self, a):
      # integrand in lnm, for speed
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.fmassfunc(m, a) * self.MassFunc.fb1(m, a) * m / self.U.rho_z(1./a-1.)
         result *= m # because integrating in lnm and not m
         return result

      return integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0, epsrel=1.e-2)[0]


   def plotCorrectionFactor(self):
      A = np.copy(self.MassFunc.A)
      CorrFactor = np.array(map(self.correction_factor_m11, A))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(1./A-1., CorrFactor, 'b')
      
      plt.show()



