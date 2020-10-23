from headers import *

##################################################################################
##################################################################################

class P3dRsdAuto(object):
   
   def __init__(self, U, Prof, MassFunc, name=""):
      # copy classes
      self.U = U
      #self.IHaloModel = IHaloModel
      self.MassFunc = MassFunc
      self.Prof = Prof
      
      # mass bounds for integrals
      self.mMin = max(self.Prof.mMin, self.MassFunc.mMin)
      self.mMax = min(self.Prof.mMax, self.MassFunc.mMax)


      #self.fPnoise = lambda k,z: self.Prof.Pshot(z)
      #self.fPnoise = fPnoise
      #self.fTnoise = fTnoise
      #self.Vs = Vs   # in (Mpc/h)^3
      self.name = str(self.Prof) + name
      #self.nProc = nProc
      
      # values of k to evaluate
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      self.Ke = np.genfromtxt("./input/K.txt") # edges of the bins for k
      self.dK = np.genfromtxt("./input/dK.txt") # widths of the bins for k
      # redshifts to evaluate
      #self.Z = np.linspace(0., 10., 11)
      self.Z = self.Prof.Lf.Z
   
      # create folder if needed
      directory = "./output/p_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)

   def __str__(self):
      return self.name


   ##################################################################################
   # Power spectrum ingredients


   def p1h(self, k, z, mu=0.):
      '''
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, 1./(1.+z))
         result *= self.Prof.u(k, m, z, mu)**2
         result *= m # because integrating in lnm and not m
         return result

      result = integrate.quad(integrand, np.log(self.mMin), np.log(self.mMax), epsabs=0, epsrel=1.e-3)[0]
      return result



   def correctionFactorI11(self, z):
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z) 
         result *= self.MassFunc.b1(m, z) 
         result *= m / self.U.rho_m(z)
         result *= m # because integrating in lnm and not m
         return result

      return integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0, epsrel=1.e-2)[0]



   def bEff(self, k, z, mu=0.):
      '''
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof.u(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result

      result = integrate.quad(integrand, np.log(self.mMin), np.log(self.mMax), epsabs=0, epsrel=1.e-3)[0]
      if self.Prof.use_correction_factor==1:
         result /= self.correctionFactorI11(z)
      return result

   def fEff(self, k, z, mu=0.):
      '''
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.uMat(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result

      result = integrate.quad(integrand, np.log(self.mMin), np.log(self.mMax), epsabs=0, epsrel=1.e-3)[0]
      result *= self.U.bg.scale_independent_growth_rate(z)
      # fix the halo model total mass problem
      result /= self.correctionFactorI11(z)
      return result

   def p2h(self, k, z, mu=0.):
      '''
      '''
      result = self.bEff(k, z, mu=mu) + self.fEff(k, z, mu=mu) * mu**2
      result = result**2 * self.U.pLin(k, z)
      return result

   def pShot(self, z):
      '''
      '''
      return self.Prof.Pshot(z)

   def p(self, k, z, mu=0.):
      '''
      '''
      return self.p1h(k, z, mu=mu) + self.p2h(k, z, mu=mu)

   def pTot(self, k, z, mu=0.):
      '''
      '''
      return self.p1h(k, z, mu=mu) + self.p2h(k, z, mu=mu) + self.pShot(z)

   
   ##################################################################################


   def plotPMuDpdce(self, z=0.):
      
      Mu = np.array([0., 0.5, 1.])

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      for mu in Mu:
         f = lambda k: self.pTot(k, z, mu=mu)
         p = np.array(map(f, self.K))
         #
         ax.loglog(self.K, p, lw=2, label=r'$\mu=$'+str(round(mu, 1)))
      #
      ax.grid()
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$P(k)$')
      #path = "./figures/pn3d/p_"+self.name+"z_"+str(z)+".pdf"
      #fig.savefig(path, bbox_inches='tight')

      plt.show()



   def plotBEff(self):
      Z = np.linspace(0., 10., 101)
      bEff = np.array(map(lambda z: self.bEff(1.e-4, z), Z))

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, bEff, lw=2)
      #
      ax.set_xlabel(r'redshift $z$')
      ax.set_ylabel(r'$\sqrt{P^\text{2h}(k=0) / P^\text{lin}(k=0)}$')
      ax.set_title(r'Effective bias')
      #
      #fig.savefig("./figures/pn3d/beff_"+self.name+".pdf", bbox_inches='tight')

      plt.show()




##################################################################################
##################################################################################


