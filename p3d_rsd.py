from headers import *

##################################################################################
##################################################################################

class P3dRsdAuto(object):
   
   def __init__(self, U, Prof, MassFunc, name="", nProc=1):
      # copy classes
      self.U = U
      #self.IHaloModel = IHaloModel
      self.MassFunc = MassFunc
      self.Prof = Prof
      self.nProc = nProc
      self.name = str(self.Prof) + name
      
      # mass bounds for integrals
      self.mMin = max(self.Prof.mMin, self.MassFunc.mMin)
      self.mMax = min(self.Prof.mMax, self.MassFunc.mMax)


      # k moduli, to precompute
      self.K = np.logspace(np.log10(1.e-3), np.log10(1.), 31, 10.)
      self.nK = len(self.K)
      self.kMin = np.min(self.K)
      self.kMax = np.max(self.K)
         
      # mu values, to precompute
      self.Mu = np.linspace(0., 1., 11)
      self.nMu = len(self.Mu)
      self.muMin = np.min(self.Mu)
      self.muMax = np.max(self.Mu)

      ## values of k to evaluate
      #self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      #self.Ke = np.genfromtxt("./input/K.txt") # edges of the bins for k
      #self.dK = np.genfromtxt("./input/dK.txt") # widths of the bins for k

      # redshifts to evaluate, from the luminosity function
      self.Z = self.Prof.Lf.Z
   
      # create folder if needed
      directory = "./output/p3d_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)

   def __str__(self):
      return self.name


   ##################################################################################
   # Power spectrum ingredients


   def p1h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.u(k, m, z, mu)**2
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result



   def correctionFactorI11(self, z, mMin=0., mMax=np.inf):
      '''Used for the matter 2-halo term, or in fEff,
      to correct the fact that 
      int dm n(m) m/rho = 1
      is not satisfied for usual halo mass functions.
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z) 
         result *= self.MassFunc.b1(m, z) 
         result *= m / self.U.rho_m(z)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      return integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]



   def bEff(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof.u(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      if self.Prof.use_correction_factor==1:
         result /= self.correctionFactorI11(z, mMin=mMin, mMax=mMax)
      return result

   def fEff(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Effective growth rate of structure
      Converges to f when k-->0
      Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.uMat(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      result *= self.U.bg.scale_independent_growth_rate(z)
      # fix the halo model total mass problem
      result /= self.correctionFactorI11(z, mMin=mMin, mMax=mMax)
      return result

   def p2h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes RSD: Kaiser effect and FOGs
      '''
      result = self.bEff(k, z, mu=mu, mMin=mMin, mMax=mMax) 
      result += self.fEff(k, z, mu=mu, mMin=mMin, mMax=mMax) * mu**2
      result = result**2 * self.U.pLin(k, z)
      return result

   def pShot(self, z, mMin=0., mMax=np.inf):
      '''(RSD has no effect on the shot noise power spectrum)
      '''
      result = self.Prof.pShot(z)
      # build up with mass,
      # assuming phi(L|m) = N(m) * phi(L)
      result *= self.Prof.Sfr.sfrdForInterp(z, mMin=mMin, mMax=mMax) / self.Prof.Sfr.sfrd(z)
      return result

   def p(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Sum of 1h and 2h terms
      '''
      result = self.p1h(k, z, mu=mu, mMin=mMin, mMax=mMax) 
      result += self.p2h(k, z, mu=mu, mMin=mMin, mMax=mMax)
      return result

   def pTot(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Sum of 1h, 2h and shot noise terms
      '''
      result = self.p1h(k, z, mu=mu, mMin=mMin, mMax=mMax) 
      result += self.p2h(k, z, mu=mu, mMin=mMin, mMax=mMax) 
      result += self.pShot(z, mMin=mMin, mMax=mMax)
      return result

   
   ##################################################################################


   def save(self, z=1):
      '''At the requested redshift, precompute on a (k,mu) grid:
      bEff, fEff, P2h, P1h, Pshot, Ptot
      and save to file
      '''
      print("Precompute the RSD power spectrum at z="+str(z))
      # arrays to compute
      bEff = np.zeros((self.nK, self.nMu))
      fEff = np.zeros((self.nK, self.nMu))
      p2h = np.zeros((self.nK, self.nMu))
      p1h = np.zeros((self.nK, self.nMu))
      pShot = np.zeros((self.nK, self.nMu))
      pTot = np.zeros((self.nK, self.nMu))

      # Evaluate in parallel
      with sharedmem.MapReduce(np=self.nProc) as pool:
         for iMu in range(self.nMu):
            mu = self.Mu[iMu]
            # compute for all k
            f = lambda k: self.bEff(k, z, mu)
            bEff[:,iMu] = np.array(pool.map(f, self.K))
            f = lambda k: self.fEff(k, z, mu)
            fEff[:,iMu] = np.array(pool.map(f, self.K))
            f = lambda k: self.p2h(k, z, mu)
            p2h[:,iMu] = np.array(pool.map(f, self.K))
            f = lambda k: self.p1h(k, z, mu)
            p1h[:,iMu] = np.array(pool.map(f, self.K))
            pShot[:,iMu] = self.pShot(z) * np.ones(self.nK)
            pTot[:,iMu] = p1h[:,iMu] + p2h[:,iMu] + pShot[:,iMu]
            print("- done "+str(iMu+1)+" of "+str(self.nMu))

      # save to files
      np.savetxt('./output/p3d_rsd/beff_'+self.name+'_z'+str(z)+'.txt', bEff)
      np.savetxt('./output/p3d_rsd/feff_'+self.name+'_z'+str(z)+'.txt', fEff)
      np.savetxt('./output/p3d_rsd/p1h_'+self.name+'_z'+str(z)+'.txt', p1h)
      np.savetxt('./output/p3d_rsd/p2h_'+self.name+'_z'+str(z)+'.txt', p2h)
      np.savetxt('./output/p3d_rsd/pshot_'+self.name+'_z'+str(z)+'.txt', pShot)
      np.savetxt('./output/p3d_rsd/ptot_'+self.name+'_z'+str(z)+'.txt', pTot)


   def load(self, z=1):
      print("Load the precomputed RSD power spectrum at z="+str(z))
      # read files
      bEff = np.genfromtxt('./output/p3d_rsd/beff_'+self.name+'_z'+str(z)+'.txt')
      fEff = np.genfromtxt('./output/p3d_rsd/feff_'+self.name+'_z'+str(z)+'.txt')
      p1h = np.genfromtxt('./output/p3d_rsd/p1h_'+self.name+'_z'+str(z)+'.txt')
      p2h = np.genfromtxt('./output/p3d_rsd/p2h_'+self.name+'_z'+str(z)+'.txt')
      pShot = np.genfromtxt('./output/p3d_rsd/pshot_'+self.name+'_z'+str(z)+'.txt')
      pTot = np.genfromtxt('./output/p3d_rsd/ptot_'+self.name+'_z'+str(z)+'.txt')
      
      # interpolate
      self.bEffInt = {}
      self.fEffInt = {}
      self.p1hInt = {}
      self.p2hInt = {}
      self.pShotInt = {}
      self.pTotInt = {}
      #
      self.bEffInt[z] = interp2d(self.K, self.Mu, bEff.T, kind='linear', bounds_error=False, fill_value=0.)
      self.fEffInt[z] = interp2d(self.K, self.Mu, fEff.T, kind='linear', bounds_error=False, fill_value=0.)
      self.p1hInt[z] = interp2d(self.K, self.Mu, p1h.T, kind='linear', bounds_error=False, fill_value=0.)
      self.p2hInt[z] = interp2d(self.K, self.Mu, p2h.T, kind='linear', bounds_error=False, fill_value=0.)
      self.pShotInt[z] = interp2d(self.K, self.Mu, pShot.T, kind='linear', bounds_error=False, fill_value=0.)
      self.pTotInt[z] = interp2d(self.K, self.Mu, pTot.T, kind='linear', bounds_error=False, fill_value=0.)




   ##################################################################################
   # Old crappy RSD forecast

   def sBetaOverBetaFisher(self, z, R, fwhmPsf, fSky, dz):
      '''Relative uncertainty on RSD parameter beta = f/b
      at redshift z
      From Fisher forecast, without marginalizing over
      any other parameter.
      R spectral resolving power [dimless]
      fwhmPsf [rad]
      fSky [dimless]
      dz width of redshift slice
      '''
      # precompute the RSD power spectrum at the requested redshift
      try:
         self.load(z=z)
      except:
         self.save(z=z)
         self.load(z=z)


      # survey volume
      dChi = self.U.c_kms/self.U.hubble(z) * dz
      volume = 4.*np.pi*fSky  # sky area in [srd]
      volume *= self.U.bg.comoving_distance(z)**2 * dChi # volume [(Mpc/h)^3]
      # k perp max
      kPerpMax = self.U.kMaxPerpPsf(fwhmPsf, z)
      # k para max
      kParaMax = self.U.kMaxParaSpectroRes(R, z)

      def integrand(lnKPara, lnKPerp):
         kPerp = np.exp(lnKPerp)
         kPara = np.exp(lnKPara)
         k = np.sqrt(kPerp**2 + kPara**2)
         mu = kPara / k
         b = self.bEffInt[z](k, mu)
         f = self.fEffInt[z](k, mu)
         pTot = self.pTotInt[z](k, mu)
         #print k, mu, b, f, pTot
         #
         # derivative wrt f, at fixed b
         #result = (2.*f* mu**2 * (b+f*mu**2) * self.U.pLin(k, z) / pTot)**2
         # derivative wrt beta, at fixed b
         result = (2.*b**2 * mu**2 * (1.+(f/b)*mu**2) * self.U.pLin(k, z) / pTot)**2
         result *= kPerp / (2.*np.pi)**2
         result *= kPerp * kPara   # because int wrt ln
         result *= volume / 2.
         result *= 4.   # symmetry factor, since we reduce the integration domain
         return result
      
      # Fisher matrix element
      result = integrate.dblquad(integrand, np.log(self.kMin), np.log(kPerpMax), lambda x: np.log(self.kMin), lambda x: np.log(kParaMax), epsabs=0., epsrel=1.e-2)[0]
      # Unmarginalized uncertainty
      result = 1./np.sqrt(result)
      return result


   def plotRequiredAreaToDetectBeta(self):
      '''For a given spectral resolution R,
      survey depth Delta z, what sky area is required
      to give a $10\%$ measurement of the growth rate f?
      '''
      RR = np.array([40., 150., 300.])
      #Z = np.linspace(0., 7., 11)
      Z = self.Z.copy()

      def fSkyReq(z, R, dz=1., fwhmPsf=6.*np.pi/(180.*3600.), target=0.1):
         sBetaOverBeta = self.sBetaOverBetaFisher(z, R, fwhmPsf, 1., dz)
         result = (sBetaOverBeta / target)**2
         return result

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for R in RR:
         f = lambda z: fSkyReq(z, R)
         fSky = np.array(map(f, Z))
         ax.plot(Z, fSky, label=r'$\mathcal{R}=$'+str(np.int(R)))
      #
      # SPHEREx: 100 deg^2
      ax.axhline(2.*100. * (np.pi/180.)**2 / (4.*np.pi), ls='--', label=r'SPHEREx deep fields')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(Z), np.max(Z)))
      #ax.set_ylim((1.e-2, 1.e-1))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'Required sky fraction $f_\text{sky}$')
      #
      ax2=ax.twinx()
      ylim = ax.get_ylim()
      ax2.set_ylim((ylim[0] * 4.*np.pi*(180./np.pi)**2, ylim[1] * 4.*np.pi*(180./np.pi)**2))
      ax2.set_yscale('log', nonposy='clip')
      ax2.set_ylabel(r'Required sky area [deg$^2$]')
      #
      fig.savefig('./figures/p3d_rsd/fsky_tradeoff_beta.pdf', bbox_inches='tight')
      plt.show()


   ##################################################################################
   # old forecast for power spectrum monopole and quadrupole:
   # only gives the unmarginalized constraints


   def sAOverAFisher(self, i, z, R, fwhmPsf, fSky, dz, kMax=np.inf):
      '''Relative uncertainty on A_i in the multipole expansion
      P(k, mu) = A_0 P_0(k) + A_2 P_2(k) 1/2*(3mu^2-1) + ...
      Fiducial:
      P_0(k) = (1 + 13*beta/15) * I^2 * b^2 * Plin(k)
      P_2(k) = 4*beta/3 * I^2 * b^2 * Plin(k)
      at redshift z
      From Fisher forecast.
      Since the angular dependences are orthogonal,
      no need to marginalize over other A_j.
      R spectral resolving power [dimless]
      fwhmPsf [rad]
      fSky [dimless]
      dz width of redshift slice
      '''
      # precompute the RSD power spectrum at the requested redshift
      try:
         self.load(z=z)
      except:
         self.save(z=z)
         self.load(z=z)


      # survey volume
      dChi = self.U.c_kms/self.U.hubble(z) * dz
      volume = 4.*np.pi*fSky  # sky area in [srd]
      volume *= self.U.bg.comoving_distance(z)**2 * dChi # volume [(Mpc/h)^3]
      # k perp max, impose additional cut if requested
      kPerpMax = min(kMax, self.U.kMaxPerpPsf(fwhmPsf, z))
      # k para max
      kParaMax = min(kMax, self.U.kMaxParaSpectroRes(R, z))
      print "kPerpMax", kPerpMax
      print "kParaMax", kParaMax

      def integrand(lnKPara, lnKPerp):
         kPerp = np.exp(lnKPerp)
         kPara = np.exp(lnKPara)
         k = np.sqrt(kPerp**2 + kPara**2)
         if k>=kMax:
            return 0.
         mu = kPara / k
         b = self.bEffInt[z](k, mu)
         f = self.fEffInt[z](k, mu)
         beta = f / b
         pTot = self.pTotInt[z](k, mu)
         #print k, mu, b, f, pTot
         #
         if i==0:
            result = (1. + 2./3.*beta + beta**2/5.)
         elif i==2:
            result = 4./3.*beta + 4./7.*beta**2
            result *= 0.5 * (3 * mu**2 - 1.)
         result *= b**2 * self.U.pLin(k, z) 
         result = result**2 / pTot**2
         result *= kPerp / (2.*np.pi)**2
         result *= kPerp * kPara   # because int wrt ln
         result *= volume / 2.
         result *= 2.   # symmetry factor, since we reduce the integration domain
         return result
      
      # Fisher matrix element
      result = integrate.dblquad(integrand, np.log(self.kMin), np.log(kPerpMax), lambda x: np.log(self.kMin), lambda x: np.log(kParaMax), epsabs=0., epsrel=1.e-2)[0]
      # Unmarginalized uncertainty
      result = 1./np.sqrt(result)
      return result


   def plotRequiredAreaToDetectAUnmarginalized(self, i, kMax=np.inf, exp='SPHEREx'):
      '''For a given spectral resolution R,
      survey depth Delta z, what sky area is required
      to give a $10\%$ measurement of the power spectrum multipole amplitude A_i?
      '''

      if exp=='SPHEREx':
         RR = np.array([40., 150., 300.])
         fwhmPsf = 6.*np.pi/(180.*3600.)  # 6'' in [rad]
         Z = self.Z.copy()
         dz = 1.
         fSkyExp = 2. * 100. * (np.pi/180.)**2 / (4.*np.pi) # 2 * 100 deg2 deep fields
      elif exp=='COMAP':
         RR = np.array([800.])
         fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
         Z = self.Z.copy()
         dz = 1.
         fSkyExp = 2.5 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
      elif exp=='CONCERTO':
         RR = np.array([300.])
         fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]
         Z = self.Z.copy()
         dz = 1.
         fSkyExp = 2. * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2


      def fSkyReq(z, R, dz=dz, fwhmPsf=fwhmPsf, target=0.1):
         sAOverA = self.sAOverAFisher(i, z, R, fwhmPsf, 1., dz, kMax=kMax)
         result = (sAOverA / target)**2
         return result

      # Naive mode counting, to check the Fisher forecast
      # for the monopole power spectrum
      #def fSkyReqNaiveModeCounting(z, R, dz=0.5, nModes=200., kPerpMax=0.1):
      def fSkyReqNaiveModeCounting(z, R, dz=dz, fwhmPsf=fwhmPsf, target=0.1):
         nModes = 2. / target**2
         kPerpMax = kMax   # here loosely equating these two
         result = nModes / self.U.bg.comoving_distance(z)**2 / kPerpMax**2
         result *= np.pi * (1.+z) / (dz * R)
         return result

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # Sky area of experiment considered
      ax.axhline(fSkyExp, ls='--', label=exp)
      #
      for R in RR:
         f = lambda z: fSkyReq(z, R)
         fSky = np.array(map(f, Z))
         plot=ax.plot(Z, fSky, label=r'$\mathcal{R}=$'+str(np.int(R)))

         # if power spectrum monopole, compare with
         # naive Nmode forecast
         if i==0:
            f = lambda z: fSkyReqNaiveModeCounting(z, R)
            fSky = np.array(map(f, Z))
            ax.plot(Z, fSky, c=plot[0].get_color(), ls=':', alpha=0.5)
      # add legend entry for the naive mode counting
      if i==0:
         ax.plot([], [], c='gray', ls=':', alpha=0.5, label=r'Na\"ive mode counting')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(Z), np.max(Z)))
      #ax.set_ylim((1.e-2, 1.e-1))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'Required sky fraction $f_\text{sky}$')
      if i==0:
         ax.set_title(r'Power spectrum monopole')
      elif i==2:
         ax.set_title(r'Power spectrum quadrupole')
      #
      ax2=ax.twinx()
      ylim = ax.get_ylim()
      ax2.set_ylim((ylim[0] * 4.*np.pi*(180./np.pi)**2, ylim[1] * 4.*np.pi*(180./np.pi)**2))
      ax2.set_yscale('log', nonposy='clip')
      ax2.set_ylabel(r'Required sky area [deg$^2$]')
      #
      fig.savefig('./figures/p3d_rsd/fsky_tradeoff_a'+str(i)+'_'+self.name+'_'+exp+'_kmax'+str(kMax)+'.pdf', bbox_inches='tight')
      #plt.show()
      fig.clf()




   ##################################################################################
   # RSD forecast: unmarginalized and marginalized uncertainties
   # on the power spectrum monopole and quadrupole

   def fisherAij(self, i, j, z, R, fwhmPsf, fSky, dz, kMax=np.inf):
      '''Fisher matrix F_{Ai, Aj} for the amplitudes A_i:
      P(k, mu) = A_0 P_0(k) + A_2 P_2(k) 1/2*(3mu^2-1) + ...
      Fiducial:
      P_0(k) = (1 + 13*beta/15) * I^2 * b^2 * Plin(k)
      P_2(k) = 4*beta/3 * I^2 * b^2 * Plin(k)
      at redshift z
      R spectral resolving power [dimless]
      fwhmPsf [rad]
      fSky [dimless]
      dz width of redshift slice
      '''
      # precompute the RSD power spectrum at the requested redshift
      try:
         self.load(z=z)
      except:
         self.save(z=z)
         self.load(z=z)

      # survey volume
      dChi = self.U.c_kms/self.U.hubble(z) * dz
      volume = 4.*np.pi*fSky  # sky area in [srd]
      volume *= self.U.bg.comoving_distance(z)**2 * dChi # volume [(Mpc/h)^3]
      # k perp max, impose additional cut if requested
      kPerpMax = min(kMax, self.U.kMaxPerpPsf(fwhmPsf, z))
      # k para max
      kParaMax = min(kMax, self.U.kMaxParaSpectroRes(R, z))
      print "kPerpMax", kPerpMax
      print "kParaMax", kParaMax

      def integrand(lnKPara, lnKPerp):
         kPerp = np.exp(lnKPerp)
         kPara = np.exp(lnKPara)
         k = np.sqrt(kPerp**2 + kPara**2)
         if k>=kMax:
            return 0.
         mu = kPara / k
         b = self.bEffInt[z](k, mu)
         f = self.fEffInt[z](k, mu)
         beta = f / b
         pTot = self.pTotInt[z](k, mu)
         #print k, mu, b, f, pTot

         # integrand
         result = b**2 * self.U.pLin(k, z)
         result = result**2 / pTot**2
         result *= kPerp / (2.*np.pi)**2
         result *= kPerp * kPara   # because int wrt ln
         result *= volume / 2.
         result *= 2.   # symmetry factor, since we reduce the integration domain
         #
         if i==0:
            result *= (1. + 2./3.*beta + beta**2/5.)
         elif i==2:
            result *= 4./3.*beta + 4./7.*beta**2
            result *= 0.5 * (3 * mu**2 - 1.)
         #
         if j==0:
            result *= (1. + 2./3.*beta + beta**2/5.)
         elif j==2:
            result *= 4./3.*beta + 4./7.*beta**2
            result *= 0.5 * (3 * mu**2 - 1.)
         return result

      # Fisher matrix element
      result = integrate.dblquad(integrand, np.log(self.kMin), np.log(kPerpMax), lambda x: np.log(self.kMin), lambda x: np.log(kParaMax), epsabs=0., epsrel=1.e-2)[0]
      return result


   def sAOverA(self, z, R, fwhmPsf, fSky, dz, kMax=np.inf, marg=True):
      '''Compute marginalized/unmarginalized uncertainties 
      on A0, A2, such that:
      P(k, mu) = A_0 P_0(k) + A_2 P_2(k) 1/2*(3mu^2-1) + ...
      Fiducial:
      P_0(k) = (1 + 13*beta/15) * I^2 * b^2 * Plin(k)
      P_2(k) = 4*beta/3 * I^2 * b^2 * Plin(k)
      at redshift z
      R spectral resolving power [dimless]
      fwhmPsf [rad]
      fSky [dimless]
      dz width of redshift slice
      '''
      
      # build Fisher matrix
      F = np.zeros((2,2))
      F[0,0] = self.fisherAij(0, 0, z, R, fwhmPsf, fSky, dz, kMax=kMax)
      F[1,1] = self.fisherAij(2, 2, z, R, fwhmPsf, fSky, dz, kMax=kMax)
      
      # unmarginalized uncertainties
      if not marg:
         sA0 = 1. / np.sqrt(F[0,0])
         sA2 = 1. / np.sqrt(F[1,1])
         return sA0, sA2

      # marginalized uncertainties
      else:
         F[0,1] = self.fisherAij(0, 2, z, R, fwhmPsf, fSky, dz, kMax=kMax)
         F[1,0] = F[0,1]
         # invert to get marginalized uncertainties
         invF = np.linalg.inv(F)
         sA0 = np.sqrt(invF[0,0])
         sA2 = np.sqrt(invF[1,1])
         return sA0, sA2


   def plotRequiredAreaToDetectA(self, kMax=np.inf, exp='SPHEREx', marg=True):
      '''For a given spectral resolution R,
      survey depth Delta z, what sky area is required
      to give a $10\%$ measurement of the power spectrum multipole amplitude A_i?
      '''

      if exp=='SPHEREx':
         #RR = np.array([40.])
         RR = np.array([40., 150., 300.])
         fwhmPsf = 6.*np.pi/(180.*3600.)  # 6'' in [rad]
         #Z = np.linspace(self.Z.min(), self.Z.max(), 2)
         dz = 1.
         fSkyExp = 2. * 100. * (np.pi/180.)**2 / (4.*np.pi) # 2 * 100 deg2 deep fields
      elif exp=='COMAP':
         RR = np.array([800.])
         fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
         dz = 1.
         fSkyExp = 2.5 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
      elif exp=='CONCERTO':
         RR = np.array([300.])
         fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]
         dz = 1.
         fSkyExp = 2. * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2

      Z = self.Z.copy()
      # avoid z=0
      if Z[0]<0.01:
         Z = Z[1:]

      def fSkyReq(z, R, dz=dz, fwhmPsf=fwhmPsf, target=0.1):
         # relative uncertainties on power spectrum monopole and quadrupole
         s0, s2 = self.sAOverA(z, R, fwhmPsf, 1., dz, kMax=kMax, marg=marg)
         result0 = (s0 / target)**2
         result2 = (s2 / target)**2
         return result0, result2

      # Naive mode counting, to check the Fisher forecast
      # for the monopole power spectrum
      #def fSkyReqNaiveModeCounting(z, R, dz=0.5, nModes=200., kPerpMax=0.1):
      def fSkyReqNaiveModeCounting(z, R, dz=dz, fwhmPsf=fwhmPsf, target=0.1):
         nModes = 2. / target**2
         kPerpMax = kMax   # here loosely equating these two
         result = nModes / self.U.bg.comoving_distance(z)**2 / kPerpMax**2
         result *= np.pi * (1.+z) / (dz * R)
         return result
      
      # Monopole
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # Sky area of experiment considered
      ax.axhline(fSkyExp, ls='--', label=exp)


      # Quadrupole
      fig2=plt.figure(2)
      ax2=fig2.add_subplot(111)
      #
      # Sky area of experiment considered
      ax2.axhline(fSkyExp, ls='--', label=exp)


      for R in RR:
         f = lambda z: fSkyReq(z, R)
         fSky = np.array(map(f, Z))
         print fSky.shape

         plot=ax.plot(Z, fSky[:,0], label=r'$\mathcal{R}=$'+str(np.int(R)))
         plot2=ax2.plot(Z, fSky[:,1], label=r'$\mathcal{R}=$'+str(np.int(R)))

         # if power spectrum monopole, compare with
         # naive Nmode forecast
         f = lambda z: fSkyReqNaiveModeCounting(z, R)
         fSky = np.array(map(f, Z))
         ax.plot(Z, fSky, c=plot[0].get_color(), ls=':', alpha=0.2)

      # add legend entry for the naive mode counting
      ax.plot([], [], c='gray', ls=':', alpha=0.5, label=r'Na\"ive mode counting')
      

      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(Z), np.max(Z)))
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'Required sky fraction $f_\text{sky}$')
      ax.set_title(r'Power spectrum monopole')
      #
      # have the ticks in scientific format 
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      axr=ax.twinx()
      ylim = ax.get_ylim()
      axr.set_ylim((ylim[0] * 4.*np.pi*(180./np.pi)**2, ylim[1] * 4.*np.pi*(180./np.pi)**2))
      axr.set_yscale('log', nonposy='clip')
      axr.set_ylabel(r'Required sky area [deg$^2$]')
      #
      # have the ticks in scientific format 
      axr.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      axr.yaxis.set_major_locator(LogLocator(numticks=15))
      axr.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      if marg:
         fig.savefig('./figures/p3d_rsd/fsky_tradeoff_a0_marg_'+self.name+'_'+exp+'_kmax'+str(kMax)+'.pdf', bbox_inches='tight')
      else:
         fig.savefig('./figures/p3d_rsd/fsky_tradeoff_a0_unmarg_'+self.name+'_'+exp+'_kmax'+str(kMax)+'.pdf', bbox_inches='tight')
      #plt.show()
      fig.clf()


      ax2.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax2.set_yscale('log', nonposy='clip')
      ax2.set_xlim((np.min(Z), np.max(Z)))
      ax2.set_xlabel(r'$z$')
      ax2.set_ylabel(r'Required sky fraction $f_\text{sky}$')
      ax2.set_title(r'Power spectrum quadrupole')
      #
      # have the ticks in scientific format 
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      ax2r=ax2.twinx()
      ylim = ax2.get_ylim()
      ax2r.set_ylim((ylim[0] * 4.*np.pi*(180./np.pi)**2, ylim[1] * 4.*np.pi*(180./np.pi)**2))
      ax2r.set_yscale('log', nonposy='clip')
      ax2r.set_ylabel(r'Required sky area [deg$^2$]')
      #
      # have the ticks in scientific format 
      axr.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      axr.yaxis.set_major_locator(LogLocator(numticks=15))
      axr.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      if marg:
         fig2.savefig('./figures/p3d_rsd/fsky_tradeoff_a2_marg_'+self.name+'_'+exp+'_kmax'+str(kMax)+'.pdf', bbox_inches='tight')
      else:
         fig2.savefig('./figures/p3d_rsd/fsky_tradeoff_a2_unmarg_'+self.name+'_'+exp+'_kmax'+str(kMax)+'.pdf', bbox_inches='tight')
      #plt.show()
      fig2.clf()





   ##################################################################################



   def plotP(self, z=None, mu=0., unit='Lsun/(Mpc/h)^2/sr/Hz'):
      '''Choice of intensity units:
      'Lsun/(Mpc/h)^2/sr/Hz'
      'cgs' for [erg/s/cm^2/sr/Hz]
      'Jy/sr'
      '''
      unitConversion = self.Prof.Lf.convertPowerSpectrumUnit(unit) 

      if z is None:
         z = self.Z[0]

      # Plin
      f = lambda k: self.U.pLin(k, z)
      Plin = np.array(map(f, self.K)) * unitConversion
      # P1h
      f = lambda k: self.p1h(k, z, mu)
      #f = lambda k: self.p1hInt[z](k, mu)
      P1h = np.array(map(f, self.K)) * unitConversion
      # P2h
      f = lambda k: self.p2h(k, z, mu)
      #f = lambda k: self.p2hInt[z](k, mu)
      P2h = np.array(map(f, self.K)) * unitConversion
      # noise bias
      f = lambda k: self.pShot(z)
      #f = lambda k: self.pShotInt[z](k, mu)
      Pnoise = np.array(map(f, self.K)) * unitConversion
      # Ptot
      P = P1h + P2h + Pnoise

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, P, 'k', lw=4, label=r'$P_\text{tot}$')
      ax.loglog(self.K, P2h, 'b-', lw=2, label=r'$P_\text{2h}$')
      #ax.loglog(self.K, self.Prof.bias(1./(1.+z))**2 * Plin, 'b--', lw=2, label=r'$b_\text{eff}^2 P_\text{lin}$')
      #ax.loglog(self.K, Plin, 'k--', lw=2, label=r'$P_\text{lin}$')
      ax.loglog(self.K, P1h, 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.K, Pnoise, 'g-', lw=2, label=r'$P_\text{noise}$')

      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k)$ [$'+unit+r'$(Mpc/$h$)$^3$')
      #path = "./figures/pn3d/p_"+self.name+"z_"+str(z)+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      plt.show()




   def plotFourierModes(self):
      # width of redshift slice
      dz = 1.

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      currentAxis = plt.gca()
      #
      # k_perp = k_para line
      kk = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
      ax.plot(kk, kk, c='gray', alpha=0.3, label=r'$k_\perp = k_\parallel$')
      #
      # iso |k| contours
      for k in [1.e-3, 1.e-2, 1.e-1, 1., 10., 100.]:
         kk = np.logspace(np.log10(1.e-3), np.log10(k), 101, 10.)
         ax.plot(kk, np.sqrt(k**2-kk**2), c='gray', alpha=0.2, lw=1.)
      #
      # iso-mu lines
      kk = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
      #for theta in np.linspace(0., np.pi/2., 11):
      #   ax.plot(kk * np.tan(theta), kk, c='gray', alpha=0.2, lw=1.)
      for mu in np.linspace(0., 1., 11):
         theta = np.arccos(mu)
         ax.plot(kk * np.tan(theta), kk, c='gray', alpha=0.2, lw=1.)


      # Specs
      exp = 'SPHEREx'
      z = 0.81
      RR = np.array([40., 150., 300.])
      R = RR[0]
      fwhmPsf = 6.*np.pi/(180.*3600.)  # 6'' in [rad]
      fSkyExp = 2. * 100. * (np.pi/180.)**2 / (4.*np.pi) # 2 * 100 deg2 deep fields
      # scales probed
      kFPerp = self.U.kFPerp(z, fSkyExp)
      kMaxPerp = self.U.kMaxPerpPsf(fwhmPsf, z)
      kFPara = self.U.kFPara(z, dz=dz)
      kMaxPara = self.U.kMaxParaSpectroRes(R, z)
      #
      currentAxis.add_patch(Rectangle((kFPerp, kFPara), kMaxPerp-kFPerp, kMaxPara-kFPara,
                      alpha=1., facecolor='none', edgecolor='r', linewidth=2,
                      label=exp+r' $z=$'+str(round(z, 1))))


      # Specs
      exp = 'COMAP'
      z = 2.
      RR = np.array([800.])
      R = RR[0]
      fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
      dz = 1.
      fSkyExp = 2.5 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
      # scales probed
      kFPerp = self.U.kFPerp(z, fSkyExp)
      kMaxPerp = self.U.kMaxPerpPsf(fwhmPsf, z)
      kFPara = self.U.kFPara(z, dz=dz)
      kMaxPara = self.U.kMaxParaSpectroRes(R, z)
      #
      currentAxis.add_patch(Rectangle((kFPerp, kFPara), kMaxPerp-kFPerp, kMaxPara-kFPara,
                      alpha=1., facecolor='none', edgecolor='b', linewidth=2,
                      label=exp+r' $z=$'+str(round(z, 1))))

      
      # Specs
      exp = 'CONCERTO'
      z = 6.
      RR = np.array([300.])
      R = RR[0]
      fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]
      dz = 1.
      fSkyExp = 2. * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
      # scales probed
      kFPerp = self.U.kFPerp(z, fSkyExp)
      kMaxPerp = self.U.kMaxPerpPsf(fwhmPsf, z)
      kFPara = self.U.kFPara(z, dz=dz)
      kMaxPara = self.U.kMaxParaSpectroRes(R, z)
      #
      currentAxis.add_patch(Rectangle((kFPerp, kFPara), kMaxPerp-kFPerp, kMaxPara-kFPara,
                      alpha=1., facecolor='none', edgecolor='g', linewidth=2,
                      label=exp+r' $z=$'+str(round(z, 1))))


      
      ax.set_title('Fourier plane coverage')
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((1.e-3, 100.))
      ax.set_ylim((1.e-3, 100.))
      ax.set_aspect('equal', adjustable='box')
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k_\perp$ [$h$/Mpc]')
      ax.set_ylabel(r'$k_\parallel$ [$h$/Mpc]')
      path = "./figures/p3d_rsd/fourier_plane_coverage.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


















   def plotPMuDpdce(self, z=None, exp='SPHEREx'):

      if exp=='SPHEREx':
         #RR = np.array([40.])
         RR = np.array([40., 150., 300.])
         fwhmPsf = 6.*np.pi/(180.*3600.)  # 6'' in [rad]
         Z = self.Z.copy()
         #Z = np.linspace(self.Z.min(), self.Z.max(), 2)
         dz = 1.
         fSkyExp = 2. * 100. * (np.pi/180.)**2 / (4.*np.pi) # 2 * 100 deg2 deep fields
      elif exp=='COMAP':
         RR = np.array([800.])
         fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
         Z = self.Z.copy()
         dz = 1.
         fSkyExp = 2.5 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
      elif exp=='CONCERTO':
         RR = np.array([300.])
         fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]
         Z = self.Z.copy()
         dz = 1.
         fSkyExp = 2. * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2

      if z is None:
         z = self.Z[0]
         print("z="+str(z))
      Mu = np.array([0., 0.5, 1.])
      K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      for mu in Mu:
         f = lambda k: self.pTot(k, z, mu=mu) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         p = np.array(map(f, K))
         #
         plot=ax.loglog(K, p, lw=2, label=r'$\mu=$'+str(round(mu, 1)), c=plt.cm.cool(mu))
         #
         # Show the scales probed by SPHEREx
         if mu==0.:
            # k_perp
            kFPerp = self.U.kFPerp(z, fSkyExp)
            kMaxPerp = self.U.kMaxPerpPsf(fwhmPsf, z)
            print "perp k_min, k_max =", kFPerp, kMaxPerp
            #ax.hlines(5.e3, xmin=kFPerp, xmax=kMaxPerp, colors=plot[0].get_color())
            ax.axvspan(kFPerp, kMaxPerp, fc=plot[0].get_color(), ec=None, alpha=0.1)
         if mu==1.:
            # k_para
            R = RR[0]
            kFPara = self.U.kFPara(z, dz=dz)
            kMaxPara = self.U.kMaxParaSpectroRes(R, z)
            print "para k_min, k_max =", kFPara, kMaxPara
            #ax.hlines(5.e3, xmin=kFPara, xmax=kMaxPara, colors=plot[0].get_color())
            ax.axvspan(kFPara, kMaxPara, fc=plot[0].get_color(), ec=None, alpha=0.1)
      #
      ax.set_title(self.Prof.Lf.lineNameLatex+' with '+exp+r' $(z=$'+str(round(z, 1))+r'$)$')
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k, \mu, z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      #
      path = "./figures/p3d_rsd/p_mudpdce_"+self.name+"z_"+str(z)+"_"+exp+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()



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


   def compareP(self, ps=None):
      if ps is None:
         ps = [self]

      # find min and max redshifts in all the LFs to plot
      zMin = np.min([p.Prof.Lf.zMin for p in ps])
      zMax = np.max([p.Prof.Lf.zMax for p in ps])

      fig=plt.figure(0, figsize=(8,8))
      ax=fig.add_subplot(111)
      #
      lineStyles = np.array(['-', '--', ':'])
      legendItems = []
      for iP in range(len(ps)):
         print "working on", p.name
         p = ps[iP]
         ls = lineStyles[iP]
         for z in p.Prof.Lf.Z:
            c = plt.cm.cool((z-zMin)/(zMax-zMin))
            if z<=6:
               f = lambda k: p.pTot(k, z) #, mu=0.5)
               pTot = np.array(map(f, p.K)) * p.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
               line, =ax.loglog(p.K, pTot, ls=ls, c=c)
               label=r'$z=$'+str(round(z, 2))+' '+p.Prof.Lf.refLatex
               legendItems.append((z, line, label))
      #
      #ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      legendItems.sort()
      ax.legend([x[1] for x in legendItems], [x[2] for x in legendItems], loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k,z, \mu=0)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      #ax.set_ylim((1.e3, 1.e8))
      ax.set_title(self.Prof.Lf.lineNameLatex+' Power spectrum')
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      path = './figures/p3d_rsd/ptot_'+self.name+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      #fig.clf()
      plt.show()


   def plotCumulMassContributionP(self, z=None, mu=0.):
      if z is None:
         z = self.Z[0]
      M = np.logspace(np.log10(5.e10), np.log10(4.e14), 6, 10.) # [Msun/h]
      K = np.logspace(np.log10(1.e-3), np.log10(40.), 101, 10.) # [h/Mpc]

      # Compute the build up of the power spectrum,
      # as the maximum mass increases
      P = np.zeros((len(M)+1, len(K)))
      for iM in range(len(M)):
         m = M[iM]
         f = lambda k: self.pTot(k, z, mu, mMin=None, mMax=m) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         P[iM+1, :] = np.array(map(f, K))
         print P[iM+1,:]
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iM in range(len(M))[::-1]:
         m = M[iM]
         #
         ax.fill_between(K, P[iM,:], P[iM+1,:], facecolor=plt.cm.YlOrRd(1.*iM/(len(M)-1.)), edgecolor='', label=r'$m\leqslant$'+floatSciForm(m, round=1)+r' $M_\odot/h$')

      #
      ax.legend(loc=3, fontsize=12, labelspacing=0.1, handlelength=0.5)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(K), np.max(K)))
      #ax.set_ylim((1.e-14, 1.e-7))
      ax.set_ylim((10., 5.e7))
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      #ax.set_ylabel(r'$P(k, \mu='+str(mu)+', z='+str(z)+')$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_ylabel(r'$P(k, \mu, z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_title(r'$z='+str(z)+', \mu='+str(mu)+'$')
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      path = './figures/p3d_rsd/ptot_'+self.name+'_m_cumul_mu'+str(mu)+'_z'+str(z)+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def plotMassContributionP(self, z=None, mu=0.):
      if z is None:
         z = self.Z[0]
      M = np.logspace(np.log10(5.e10), np.log10(7.e13), 7, 10.) # [Msun/h]
      K = np.logspace(np.log10(1.e-3), np.log10(40.), 101, 10.) # [h/Mpc]

      # Compute the build up of the power spectrum,
      # as the maximum mass increases
      P = np.zeros((len(M)-1, len(K)))
      for iM in range(len(M)-1):
         f = lambda k: self.pTot(k, z, mu, mMin=M[iM], mMax=M[iM+1]) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         P[iM, :] = np.array(map(f, K))
         print P[iM,:]


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # full power spectrum
      f = lambda k: self.p(k, z, mu)
      p = np.array(map(f, K)) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
      ax.plot(K, p, 'k-', label=r'Total')
      #
      for iM in range(len(M)-1)[::-1]:
         m = M[iM]
         #
         ax.plot(K, P[iM,:], ls='--', c=plt.cm.YlOrRd((iM+1.)/(len(M)-1.)), label=floatSciForm(M[iM], round=1)+r'$\leqslant m<$'+floatSciForm(M[iM+1], round=1)+r' $M_\odot/h$')
      #
      ax.legend(loc=3, fontsize=12, labelspacing=0.1, handlelength=0.5)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((np.min(K), np.max(K)))
      #ax.set_ylim((1.e-14, 1.e-7))
      ax.set_ylim((0.1, 5.e7))
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      #ax.set_ylabel(r'$P(k, \mu='+str(mu)+', z='+str(z)+')$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_ylabel(r'$P(k, \mu, z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_title(r'$z='+str(z)+', \mu='+str(mu)+'$')
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      # to get more tick marks on the y axis
      ax.yaxis.set_major_locator(LogLocator(numticks=15))
      ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      path = './figures/p3d_rsd/ptot_'+self.name+'_m_mu'+str(mu)+'_z'+str(z)+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def plotPTermsZ(self, mu=0.):

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in self.Prof.Lf.Z:    
         if z<=5:
            f = lambda k: self.pTot(k, z, mu) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            pTot = np.array(map(f, self.K)) 
            f = lambda k: self.p1h(k, z, mu) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            p1h = np.array(map(f, self.K))
            #f = lambda k: self.fP2hinterp(k, z) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            #p2h = np.array(map(f, self.K))
            f = lambda k: self.pShot(z) * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            pShot = np.array(map(f, self.K))
            #
            plot=ax.loglog(self.K, pTot, label=r'$z=$'+str(round(z, 2)))
            ax.loglog(self.K, p1h, ls=(0, (10, 3)), c=plot[0].get_color(), lw=1)
            #ax.loglog(self.K, p2h, ls='-', c=plot[0].get_color())
            ax.loglog(self.K, pShot, ls=':', c=plot[0].get_color(), lw=1)
      #
      ax.plot([], [], ls='-', c='gray', alpha=0.5, label=r'total')
      ax.plot([], [], ls=(0, (10, 3)), c='gray', alpha=0.5, label=r'1h')
      ax.plot([], [], ls=':', c='gray', alpha=0.5, label=r'shot')
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k,z, \mu)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_ylim((4.e3, 5.e7))
      ax.set_title(r'$\mu='+str(mu)+'$')
      #ax.set_title(self.Prof.Lf.nameLatex)
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      path = './figures/p3d_rsd/p1h2hshot_'+self.name+'_mu'+str(mu)+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()









   ##################################################################################
   # 3d matched filter for point sources in redshift-space

   def sigmaFluxMatchedFilter(self, detNoisePower, R, fwhmPsf, z):
      '''Computes the uncertainty on the flux [Lsun / (Mpc/h)^2]
      of a point source in the map,
      given the detector noise power spectrum detNoisePower.
      detNoisePower: default unit [(Lsun/(Mpc/h)^2/sr/Hz)^2 (Mpc/h)^3],
      ie [Lsun^2/(Mpc/h)/sr^2/Hz^2].
      This assumes a matched filter is used, with the point source profile
      determined by the PSF and SPSF of the instrument.
      R: spectral resolving power [dimless]
      fwhmPsf: [rad]
      '''
      # precompute the RSD power spectrum at the requested redshift
      try:
         self.load(z=z)
      except:
         self.save(z=z)
         self.load(z=z)

      def integrand(lnKPara, lnKPerp):
         kPerp = np.exp(lnKPerp)
         kPara = np.exp(lnKPara)
         k = np.sqrt(kPerp**2 + kPara**2)
         mu = kPara / k
         
         # keep default unit [(Lsun/(Mpc/h)^2/sr/Hz)^2 (Mpc/h)^3]
         # ie [Lsun^2/(Mpc/h)/sr^2/Hz^2]
         pTot = self.pTotInt[z](k, mu)
         if pTot==0.:
            print "watch out ", k, mu, pTot
         
         psf = self.U.psfF(kPerp, fwhmPsf, z)   # [dimless]
         spsf = self.U.spectralPsfF(kPara, R, z)   # [dimless]
         w = psf * spsf

         result = w**2   # [dimless]
         result /= 2. * (w**2 * pTot + detNoisePower)   # [Lsun^-2 * (Mpc/h) * sr^2 * Hz^2]
         result *= kPerp / (2.*np.pi)**2  # [Lsun^-2 * sr^2 * Hz^2]
         result *= kPerp * kPara   # [Lsun^-2 * sr^2 * Hz^2 / (Mpc/h)^2] because int wrt ln
         result *= 2.   # symmetry factor, since we reduce the integration domain
         return result  # [Lsun^-2 * sr^2 * Hz^2 / (Mpc/h)^2]

      # compute 2d integral
      # integration bounds: if modulus of k goes above self.kMax, 
      # pTot will be zero, and sigmaMatchedFilter will be zero
      # if det noise is zero. --> avoid this!
      result = integrate.dblquad(integrand, np.log(self.kMin), np.log(self.kMax/np.sqrt(2.)), lambda x: np.log(self.kMin), lambda x: np.log(self.kMax/np.sqrt(2.)), epsabs=0., epsrel=1.e-2)[0]
      result = 1. / np.sqrt(result) # [Lsun / sr / Hz * (Mpc/h)]
      result *= self.U.hubble(z) * self.Prof.Lf.nuHz / (1.+z)**2 / self.U.c_kms / self.U.bg.comoving_distance(z)**2  # * [Hz * sr / (Mpc/h)^3] = [Lsun / (Mpc/h)^2]
      return result


   def sigmaLumMatchedFilter(self, detNoisePower, R, fwhmPsf, z):
      '''Return the matched filter uncertainty in terms of luminosity [Lsun]
      '''
      # get the flux
      result = self.sigmaFluxMatchedFilter(detNoisePower, R, fwhmPsf, z)   # [Lsun / (Mpc/h)^2]
      # convert to luminosity
      result *= 4.*np.pi * (1.+z)**2 * self.U.bg.comoving_distance(z)**2   # [Lsun]
      return result



   def plotSigmaLumMatchedFilter(self, specs):

      for z in self.Z:

         # experimental specs
         exp = specs.exp
         R = specs.R
         fwhmPsf = specs.fwhmPsf
         detNoiseFid = specs.whiteNoisePower(z)
         fSkyExp = specs.fSkyExp

         # shot noise level
         pShot = self.pShot(z)   # default units involving [Lsun]
         pShot *= self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')  # [(Jy/sr)^2 (Mpc/h)^3]
         

         # detector noise values to explore
         # in default power units, converted later when plotting
         if exp=='SPHEREx':
            DetNoisePower = np.logspace(np.log10(1.e-16), np.log10(1.e-7), 21, 10.)
         elif exp=='COMAPPathfinder':
            DetNoisePower = np.logspace(np.log10(1.e-17), np.log10(1.e-9), 21, 10.)
         elif exp=='COMAP':
            DetNoisePower = np.logspace(np.log10(1.e-17), np.log10(1.e-8), 21, 10.)
         elif exp=='CONCERTO':
            DetNoisePower = np.logspace(np.log10(1.e-16), np.log10(1.e-6), 21, 10.)
         if exp=='CDIM':
            DetNoisePower = np.logspace(np.log10(1.e-16), np.log10(1.e-7), 21, 10.)
         if exp=='HETDEX':
            DetNoisePower = np.logspace(np.log10(1.e-16), np.log10(1.e-7), 21, 10.)



         # min luminosity detectable [Lsun]: 5 sigma
         f = lambda detNoisePower: 5. * self.sigmaLumMatchedFilter(detNoisePower, R, fwhmPsf, z)
         LMin = np.array(map(f, DetNoisePower))
         print "LMin", LMin

         # Convert the minimum detected luminosity lMin
         # into a minimum detected halo mass
         # get the Kennicutt-Schmidt constant
         K = self.Prof.Sfr.kennicuttSchmidtConstant(z, self.Prof.Lf, alpha=self.Prof.a)
         print "KS constant", K
         # use it to convert lMin to mMin
         f = lambda l: self.Prof.Sfr.massFromLum(l, z, K, alpha=self.Prof.a)
         MMin = np.array(map(f, LMin))

         

         # Plot minimum luminosity detectable
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         # fiducial detector noise
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
         #
         # compare with shot noise
         ax.axvline(pShot, c='gray', ls=':', alpha=0.5, label=r'Shot noise')
         #
         # lStar of the LF, to compare
         lStar = self.Prof.Lf.lStar(z) # [Lsun]
         lStar *= self.Prof.Lf.convertLumUnit('cgs')  # [cgs] = [erg/s]
         ax.axhline(lStar, ls='--', c='gray', alpha=0.5, label=r'$L^\star$')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         ax.plot(x, LMin*self.Prof.Lf.convertLumUnit('cgs'))
         #
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'$L_\text{min}$ [erg/s]')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z='+str(round(z,1))+r'$')
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         #ax2=ax.twinx()
         #ax2.set_yscale('log', nonposy='clip')
         #lTicks = ax.get_yticks()
         #mTicks = np.array(map(f, lTicks))
         #ax2.set_yticks(mTicks)
         #ax2.set_ylabel(r'$M_\text{min}$ [$M_\odot/h$]')
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_lmin_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()




         # Plot minimum halo mass detected
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         ax.plot(x, MMin)
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'$m_\text{min}$ [$M_\odot/h$]')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         ax.set_ylim((1.e11, 1.e19))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_mmin_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         ###################################

         
         # fraction of mean intensity from undetected sources
         def f(lMin):
            result = self.Prof.Lf.meanIntensity(z, lMin=0., lMax=lMin)
            result /= self.Prof.Lf.meanIntensityInterp(z)
            return result
         fracMeanIntensity = np.array(map(f, LMin))
         #
         # fraction of shot noise from undetected sources
         def f(lMin):
            result = self.Prof.Lf.pShot(z, lMin=0., lMax=lMin)
            result /= self.Prof.Lf.pShotInterp(z)
            return result
         fracShotNoise = np.array(map(f, LMin))


#         # Check: fraction of mean intensity from undetected sources
#         # from mass cut instead of luminosity cut
#         def f(mMin):
#            result = self.Prof.Sfr.sfrdForInterp(z, alpha=1, bias=False, mMin=0., mMax=mMin)
#            result /= self.Prof.Sfr.sfrdForInterp(z, alpha=1, bias=False, mMin=0., mMax=np.inf)
#            return result
#         fracMeanIntensityM = np.array(map(f, MMin))
         
         # fraction of 2h from undetected sources
         # at a fiducial k
         def f(mMin):
            k = 0.01
            result = self.p2h(k, z, mu=0., mMin=0., mMax=mMin)
            result /= self.p2h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac2h = np.array(map(f, MMin))
         #
         # fraction of 1h from undetected sources
         # at a fiducial k
         def f(mMin):
            k = 0.1
            result = self.p1h(k, z, mu=0., mMin=0., mMax=mMin)
            result /= self.p1h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac1h = np.array(map(f, MMin))

#         # Check: fraction of 2h from undetected sources,
#         # quick independent estimate
#         def f(mMin):
#            result = self.Prof.Sfr.sfrdForInterp(z, alpha=1, bias=True, mMin=0., mMax=mMin)**2
#            result /= self.Prof.Sfr.sfrdForInterp(z, alpha=1, bias=True, mMin=0., mMax=np.inf)**2
#            return result
#         frac2hApprox = np.array(map(f, MMin))



         # fraction of LIM from undetected sources
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid.')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         ax.plot(x, fracMeanIntensity, label=r'Mean')
#         ax.plot(x, fracMeanIntensityM, label=r'Mean intensity M')
         ax.plot(x, fracShotNoise, label=r'Shot')
         ax.plot(x, frac2h, label=r'2h')
#         ax.plot(x, frac2hApprox, label=r'2h ($k=0.01h/$Mpc) Approx')
         ax.plot(x, frac1h, label=r'1h')
         #
         #ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.))
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'Fraction from undetected sources')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_fracundetected_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()
         

         ###################################
         # SNR on fraction of LIM from undetected sources
         
         
         # Width of the redshift slice
         Dz = 0.5
         fSky = 1. * fSkyExp
         DetNoisePowerJy = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')

         print "SNR on observables"
         print "fsky =", fSky

         # Mean intensity
         # uncertainty
         f = lambda whiteNoisePower: specs.uncertaintyMeanIntensity(z, Dz, fSky, whiteNoisePower)
         uncertaintyMeanIntensity = np.array(map(f, DetNoisePowerJy))   # [Jy/sr]
         uncertaintyMeanIntensity /= self.Prof.Lf.convertIntensityUnit('Jy/sr')  # default intensity units
         # signal
         meanIntensity = self.Prof.Lf.meanIntensityInterp(z)
         # SNR on undetected fraction
         snrFracMeanIntensity = fracMeanIntensity * meanIntensity / uncertaintyMeanIntensity


         # 2h term
         k = 0.01
         # uncertainty
         kPerpMin = k / 2.
         kPerpMax = k * 2.
         kParaMin = 0.
         kParaMax = self.U.kMaxParaSpectroRes(R, z)
         f = lambda whiteNoisePower: specs.uncertaintyPowerAmplitude(z, Dz, kPerpMax, kParaMax, fSky, kPerpMin, kParaMin, whiteNoisePower)
         uncertaintyP2h = np.array(map(f, DetNoisePowerJy))   # [(Jy/sr)^2 * (Mpc/h)^3]
         uncertaintyP2h /= self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')  # default power units
         print "Test uncertainty p2h"
         print uncertaintyP2h
         # signal
         p2h = self.p2h(k, z, mu=0., mMin=0., mMax=np.inf)
         # SNR on undetected fraction
         snrFrac2h = frac2h * p2h / uncertaintyP2h


         # 1h term
         k = 0.1
         # uncertainty
         kPerpMin = k / 2.
         kPerpMax = k * 2.
         kParaMin = 0.
         kParaMax = self.U.kMaxParaSpectroRes(R, z)
         f = lambda whiteNoisePower: specs.uncertaintyPowerAmplitude(z, Dz, kPerpMax, kParaMax, fSky, kPerpMin, kParaMin, whiteNoisePower=whiteNoisePower)
         uncertaintyP1h = np.array(map(f, DetNoisePowerJy))
         uncertaintyP1h /= self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')  # default power units
         # signal
         p1h = self.p1h(k, z, mu=0., mMin=0., mMax=np.inf)
         # SNR on undetected fraction
         snrFrac1h = frac1h * p1h / uncertaintyP1h


         # Shot noise
         # uncertainty
         kPerpMin = self.U.kMaxPerpPsf(fwhmPsf, z) / 2.
         kPerpMax = self.U.kMaxPerpPsf(fwhmPsf, z) * 2.
         kParaMin = 0.
         kParaMax = self.U.kMaxParaSpectroRes(R, z)
         f = lambda whiteNoisePower: specs.uncertaintyPowerAmplitude(z, Dz, kPerpMax, kParaMax, fSky, kPerpMin, kParaMin, whiteNoisePower=whiteNoisePower)
         uncertaintyShotNoise = np.array(map(f, DetNoisePowerJy))
         uncertaintyShotNoise /= self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')  # default power units
         # signal
         shotNoise = self.Prof.Lf.pShotInterp(z)
         # SNR on undetected fraction
         snrFracShotNoise = fracShotNoise * shotNoise / uncertaintyShotNoise


         # SNR on fraction of LIM from undetected sources
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid.')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         ax.plot(x, snrFracMeanIntensity, label=r'Mean')
         ax.plot(x, snrFracShotNoise, label=r'Shot')
         ax.plot(x, snrFrac2h, label=r'2h')
         ax.plot(x, snrFrac1h, label=r'1h')
         #
         #ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         #ax.set_ylim((0., 1.))
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'SNR on fraction from undetected sources')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_snr_fracundetected_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()
         


         # fraction of LIM from undetected sources
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid.')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         #
         # Mean intensity
         I = np.where(snrFracMeanIntensity>=1.)[0]
         J = np.where(snrFracMeanIntensity<1.)[0]
         plot=ax.plot(x[I], fracMeanIntensity[I], label=r'Mean')
         ax.plot(x[J], fracMeanIntensity[J], ls='--', c=plot[0].get_color())
         #
         # Shot noise
         I = np.where(snrFracShotNoise>=1.)[0]
         J = np.where(snrFracShotNoise<1.)[0]
         plot=ax.plot(x[I], fracShotNoise[I] - 0.005, label=r'Shot')
         ax.plot(x[J], fracShotNoise[J] - 0.005, ls='--', c=plot[0].get_color())
         #
         # 2-halo
         I = np.where(snrFrac2h>=1.)[0]
         J = np.where(snrFrac2h<1.)[0]
         plot=ax.plot(x[I], frac2h[I] - 0.01, label=r'2h')
         ax.plot(x[J], frac2h[J] - 0.01, ls='--', c=plot[0].get_color())
         #
         # 1-halo
         I = np.where(snrFrac1h>=1.)[0]
         J = np.where(snrFrac1h<1.)[0]
         plot=ax.plot(x[I], frac1h[I] - 0.015, label=r'1h')
         ax.plot(x[J], frac1h[J] - 0.015, ls='--', c=plot[0].get_color())
         #
         #ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.1))
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel('Fraction of observable \n from undetected sources')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_fracundetected_snr_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()
         










         ###################################
         
         # Bias of LIM
         # full LIM
         bLimFull = self.Prof.Sfr.bEff(z, alpha=1.)
         # masked LIM
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=0., mMax=mMin)
         bLimLowM = np.array(map(f, MMin))
         # detected galaxies, weighted by luminosity
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=mMin, mMax=np.inf)
         bLimHighM = np.array(map(f, MMin))
         
         # Bias of detected galaxies
         # assume Ngal propto SFR, as for LIM
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=mMin, mMax=np.inf)
         bGalM = np.array(map(f, MMin))



         # Compare bias: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         #
         ax.axhline(bLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, bGalM, label=r'Galaxies $(M>M_\text{min})$')
         #
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((1., 50.))
         #ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'Effective bias $b$')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_bias_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


#         # Compare bias: unmasked vs masked LIM
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
#         #
#         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
#         #
#         ax.axhline(bLimFull, c='k', ls='--', label=r'LIM')
#         ax.plot(x, bLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
#         #
#         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
#         ax.set_xscale('log', nonposx='clip')
#         #ax.set_yscale('log', nonposy='clip')
#         ax.set_xlim((np.min(x), np.max(x)))
#         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
#         ax.set_ylabel(r'Effective bias $b$')
#         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
#         # to get more tick marks on the x axis
#         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
#         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
#         #
#         fig.savefig('./figures/p3d_rsd/'+exp+'_bias_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
#         plt.clf()
#         #plt.show()


#         # Compare bias: bright gal, number or luminosity weighted
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
#         #
#         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
#         #
#         ax.plot(x, bGalM, label=r'Number weighted galaxies $(M>M_\text{min})$')
#         ax.plot(x, bLimHighM, '--', label=r'Luminosity-weighted galaxies $(M>M_\text{min})$')
#         #
#         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
#         ax.set_xscale('log', nonposx='clip')
#         #ax.set_yscale('log', nonposy='clip')
#         ax.set_xlim((np.min(x), np.max(x)))
#         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
#         ax.set_ylabel(r'Effective bias $b$')
#         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
#         # to get more tick marks on the x axis
#         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
#         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
#         #
#         fig.savefig('./figures/p3d_rsd/'+exp+'_bias_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
#         plt.clf()
#         #plt.show()
         


         ###################################      
         
         # nGalEff of LIM
         # full LIM
         nGalEffLimFull = self.Prof.Lf.nGalEff(z)
         # masked LIM
         def f(lMin): return self.Prof.Lf.nGalEff(z, lMin=0., lMax=lMin)
         nGalEffLimLowL = np.array(map(f, LMin))
         # detected galaxies, weighted by luminosity
         def f(lMin): return self.Prof.Lf.nGalEff(z, lMin=lMin, lMax=np.inf)
         nGalEffLimHighL = np.array(map(f, LMin))


         # nGalEff of detected galaxies = nGal
         # from luminosity cut
         def f(lMin): 
            return self.Prof.Lf.nGal(z, lMin=lMin, lMax=np.inf)
         nGalEffGalL = np.array(map(f, LMin))


         # Compare nGalEff: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         #
         ax.axhline(nGalEffLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, nGalEffGalL, label=r'Galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         #ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
         #
         # have the ticks in scientific format 
         ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the x axis
         ax.xaxis.set_major_locator(LogLocator(numticks=15))
         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_neff_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


#         # Compare nGalEff: unmasked vs masked LIM
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
#         #
#         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
#         #
#         ax.axhline(nGalEffLimFull, c='k', ls='--', label=r'LIM')
#         ax.plot(x, nGalEffLimLowL, label=r'Masked LIM $(L<L_\text{min})$')
#         #
#         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
#         ax.set_xscale('log', nonposx='clip')
#         ax.set_yscale('log', nonposy='clip')
#         ax.set_xlim((np.min(x), np.max(x)))
#         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
#         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
#         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
#         # to get more tick marks on the x axis
#         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
#         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
#         #
#         fig.savefig('./figures/p3d_rsd/'+exp+'_neff_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
#         plt.clf()
#         #plt.show()


#         # Compare nGalEff: bright gal, number or luminosity weighted
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
#         #
#         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
#         #
#         ax.plot(x, nGalEffGalL, label=r'Number weighted galaxies $(L>L_\text{min})$')
#         ax.plot(x, nGalEffLimHighL, label=r'Luminosity weighted galaxies $(L>L_\text{min})$')
#         #
#         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
#         ax.set_xscale('log', nonposx='clip')
#         ax.set_yscale('log', nonposy='clip')
#         ax.set_xlim((np.min(x), np.max(x)))
#         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
#         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
#         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
#         # to get more tick marks on the x axis
#         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
#         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
#         #
#         fig.savefig('./figures/p3d_rsd/'+exp+'_neff_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
#         plt.clf()
#         #plt.show()
         


         ###################################

         # values of k to test
         K = np.array([0.01, 0.1, 1.]) # [h/Mpc]
         for iK in range(len(K)):
            k = K[iK]

            # Compare nP with 1 for detected galaxies
            # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
            #k = 0.1  # [h/Mpc]
            pLin = self.U.pLin(k, z)   # [(Mpc/h)^3]
            snrGal = bGalM**2 * pLin / (bGalM**2 * pLin + 1./nGalEffGalL)

            # compare p2h/(pShot+Ndet) with 1 for LIM
            # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
            # Get the beam value at that k
            psf = self.U.psfF(k, fwhmPsf, z)   # [dimless]
            # assume the Fourier mode is radial
            spsf = 1.   #self.U.spectralPsfF(kPara, R, z)   # [dimless]
            w = psf * spsf

            # LIM full
            snrLimFull = self.p2h(k,z) / (self.p2h(k,z) + self.pShot(z) + DetNoisePower / w**2)
            # masked LIM
            def f(mMin):
               return self.p2h(k, z, mMin=0., mMax=mMin)
            p2hLimLowM = np.array(map(f, MMin))
            def f(mMin):
               return self.pShot(z, mMin=0., mMax=mMin)
            pShotLimLowM = np.array(map(f, MMin))
            snrLimLowM = p2hLimLowM / (p2hLimLowM + pShotLimLowM + DetNoisePower / w**2)
            # detected galaxies, weighted by luminosity,
            # hence no detector noise
            def f(mMin):
               return self.p2h(k, z, mMin=mMin, mMax=np.inf)
            p2hLimHighM = np.array(map(f, MMin))
            def f(mMin):
               return self.pShot(z, mMin=mMin, mMax=np.inf)
            pShotLimHighM = np.array(map(f, MMin))
            # For the bright galaxies, one would generate the LIM
            # from the individually measured luminosities,
            # so there is no detector noise
            #snrLimHighM = p2hLimHighM / (p2hLimHighM + pShotLimHighM)
            snrLimHighM = bLimHighM**2 * pLin / (bLimHighM**2 * pLin + 1. / nGalEffLimHighL)

            # Compare nP: unmasked LIM vs bright galaxies
            fig=plt.figure(0)
            ax=fig.add_subplot(111)
            #
            ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
            #
            x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            #
            ax.plot(x, snrLimFull, c='k', ls='--', label=r'LIM')
            ax.plot(x, snrGal, label=r'Galaxies $(L>L_\text{min})$')
            #
            ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
            ax.set_xscale('log', nonposx='clip')
            #ax.set_yscale('log', nonposy='clip')
            ax.set_ylim((0., 1.1))
            #ax.set_xlim((np.min(x), np.max(x)))
            ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
            #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
            ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k_\perp='+str(k)+r'h/\text{Mpc})$')
            ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
            #
            # have the ticks in scientific format 
            ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
            # to get more tick marks on the x axis
            ax.xaxis.set_major_locator(LogLocator(numticks=15))
            ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
            #
            fig.savefig('./figures/p3d_rsd/'+exp+'_snr_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'_k'+str(k)+'.pdf')
            plt.clf()
            #plt.show()


   #         # Compare nP: unmasked vs masked LIM
   #         fig=plt.figure(0)
   #         ax=fig.add_subplot(111)
   #         #
   #         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
   #         #
   #         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
   #         #
   #         ax.plot(x, snrLimFull, c='k', ls='--', label=r'LIM')
   #         ax.plot(x, snrLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
   #         #
   #         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
   #         ax.set_xscale('log', nonposx='clip')
   #         #ax.set_yscale('log', nonposy='clip')
   #         ax.set_ylim((0., 1.1))
   #         ax.set_xlim((np.min(x), np.max(x)))
   #         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
   #         #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
   #         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k=0.1 h/\text{Mpc})$')
   #         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
   #         # to get more tick marks on the x axis
   #         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
   #         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
   #         #
   #         fig.savefig('./figures/p3d_rsd/'+exp+'_snr_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
   #         plt.clf()
   #         #plt.show()


   #         # Compare nP: bright gal, number or luminosity weighted
   #         fig=plt.figure(0)
   #         ax=fig.add_subplot(111)
   #         #
   #         ax.axvline(detNoiseFid, c='gray', alpha=0.5, label=r'Fid. detector noise')
   #         #
   #         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
   #         #
   #         ax.plot(x, snrGal, label=r'Number weighted galaxies $(L>L_\text{min})$')
   #         ax.plot(x, snrLimHighM, label=r'Luminosity weighted galaxies $(L>L_\text{min})$')
   #         #
   #         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   #         ax.set_xscale('log', nonposx='clip')
   #         #ax.set_yscale('log', nonposy='clip')
   #         ax.set_ylim((0., 1.1))
   #         ax.set_xlim((np.min(x), np.max(x)))
   #         ax.set_xlabel(r'Detector noise power [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
   #         #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
   #         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k=0.1 h/\text{Mpc})$')
   #         ax.set_title(exp+' '+self.Prof.Lf.lineNameLatex+r' $z=$'+str(round(z,1)))
   #         # to get more tick marks on the x axis
   #         ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
   #         ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
   #         #
   #         fig.savefig('./figures/p3d_rsd/'+exp+'_snr_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
   #         plt.clf()
   #         #plt.show()




   ##################################################################################



   def plotLimVsGalDet(self, P3d, Specs):
      '''Given a list of p3d_rsd objects
      and a list of experimental specs for each of them,
      generate the summary plots:
      Lmin = f(z)
      Mmin = f(z)
      frac undetected = f(z)
      '''

      # Lmin plot
      fig0=plt.figure(-1)
      ax0=fig0.add_subplot(111)

      # Mmin plot
      fig1=plt.figure(1)
      ax1=fig1.add_subplot(111)


      
      # Compute for each experiment
      for iExp in range(len(P3d)):
         p = P3d[iExp]
         specs = Specs[iExp]

         # remove z=0
         if p.Z[0]<0.1:
            Z = p.Z[1:]
         else:
            Z = p.Z.copy()

         print "#################################"
         print specs.exp, p.Prof.Lf.lineName

         # experimental specs
         exp = specs.exp
         R = specs.R
         fwhmPsf = specs.fwhmPsf
         def fDetNoise(z): 
            result = specs.whiteNoisePower(z)  # [(Jy/sr)^2 (Mpc/h)^3]
            # convert to fiducial unit [(Lsun/(Mpc/h)^2/sr/Hz)^2 (Mpc/h)^3]
            result /= p.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
            return result
         DetNoiseFid = np.array(map(fDetNoise, Z))
         


         # min luminosity detectable [Lsun]: 5 sigma
         f = lambda z: 5. * p.sigmaLumMatchedFilter(fDetNoise(z), R, fwhmPsf, z)
         LMin = np.array(map(f, Z))
         # interpolate it for next steps
         fLMin = interp1d(Z, LMin, kind='linear', bounds_error=False, fill_value=0.)

         # plot lMin
         plot=ax0.plot(Z, LMin*p.Prof.Lf.convertLumUnit('cgs'), label=exp+' '+p.Prof.Lf.lineNameLatex)
         # compare with lStar for the given line
         lStar = np.array(map(p.Prof.Lf.lStar, Z))  # [Lsun]
         lStar *= p.Prof.Lf.convertLumUnit('cgs')  # [cgs] = [erg/s]
         ax0.plot(Z+iExp*0.05, lStar, c=plot[0].get_color(), ls='--')


         # test
         print "Lmin [Lsun]"
         print LMin
         print "Lmin [erg/s]"
         print LMin*self.Prof.Lf.convertLumUnit('cgs')
         print "Lstar [erg/s]"
         print lStar



         # Convert the minimum detected luminosity lMin
         # into a minimum detected halo mass
         # get the Kennicutt-Schmidt constant
         fK = lambda z: p.Prof.Sfr.kennicuttSchmidtConstant(z, p.Prof.Lf, alpha=self.Prof.a)
         # use it to convert lMin to mMin
         f = lambda z: p.Prof.Sfr.massFromLum(fLMin(z), z, fK(z), alpha=self.Prof.a)
         MMin = np.array(map(f, Z))
         # interpolate it for next steps
         fMMin = interp1d(Z, MMin, kind='linear', bounds_error=False, fill_value=0.)

         # plot mMin,
         # except for COMAP and CONCERTO,
         # where the minimum mass is too large for any halo in the Universe
         if not (specs.exp=='COMAP' or specs.exp=='CONCERTO'):
            ax1.plot(Z, MMin, label=exp+' '+p.Prof.Lf.lineNameLatex)

         # test
         print "Mmin [Msun/h]"
         print MMin


         ###################################
         # fraction undetected

         # fraction of mean intensity from undetected sources
         def f(z):
            result = p.Prof.Lf.meanIntensity(z, lMin=0., lMax=fLMin(z))
            result /= p.Prof.Lf.meanIntensityInterp(z)
            return result
         fracMeanIntensity = np.array(map(f, Z))
         #
         # fraction of shot noise from undetected sources
         def f(z):
            result = p.Prof.Lf.pShot(z, lMin=0., lMax=fLMin(z))
            result /= p.Prof.Lf.pShotInterp(z)
            return result
         fracShotNoise = np.array(map(f, Z))


         # fraction of 2h from undetected sources
         # at a fiducial k
         def f(z):
            k = 0.01
            result = p.p2h(k, z, mu=0., mMin=0., mMax=fMMin(z))
            result /= p.p2h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac2h = np.array(map(f, Z))
         #
         # fraction of 1h from undetected sources
         # at a fiducial k
         def f(z):
            k = 0.1
            result = p.p1h(k, z, mu=0., mMin=0., mMax=fMMin(z))
            result /= p.p1h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac1h = np.array(map(f, Z))

         # test
         print "frac undetected"
         print fracMeanIntensity
         print fracShotNoise
         print frac2h
         print frac1h

         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(Z, fracMeanIntensity, label=r'Mean intensity')
         ax.plot(Z, fracShotNoise, label=r'Shot noise')
         ax.plot(Z, frac2h, label=r'2h ($k=0.01h/$Mpc)')
         ax.plot(Z, frac1h, label=r'1h ($k=0.1h/$Mpc)')
         #
         ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
         #ax.set_xscale('log', nonposx='clip')
         ax.set_ylim((0., 1.1))
         ax.set_xlabel(r'$z$')
         ax.set_ylabel(r'Fraction from undetected sources')
         ax.set_title(exp+' '+p.Prof.Lf.lineNameLatex)
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_fracundetected_z_'+p.name+'.pdf')
         fig.clf()
         #plt.show()
         


         ###################################
         # bias comparison

         # Bias of LIM
         # full LIM
         f = lambda z: p.Prof.Sfr.bEff(z, alpha=1.)
         bLimFull = np.array(map(f, Z))
         # masked LIM
         f = lambda z: p.Prof.Sfr.bEff(z, alpha=1, mMin=0., mMax=fMMin(z))
         bLimLowM = np.array(map(f, Z))
         # detected galaxies, weighted by luminosity
         f = lambda z: p.Prof.Sfr.bEff(z, alpha=1, mMin=fMMin(z), mMax=np.inf)
         bLimHighM = np.array(map(f, Z))
         
         # Bias of detected galaxies
         # assume Ngal propto SFR, as for LIM
         f = lambda z: p.Prof.Sfr.bEff(z, alpha=1, mMin=fMMin(z), mMax=np.inf)
         bGalM = np.array(map(f, Z))



         # Compare bias: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(Z, bGalM, label=r'Galaxies $(M>M_\text{min})$')
         ax.plot(Z, bLimFull, c='k', ls='--', label=r'LIM')
         #
         ax.legend(loc='center right', fontsize='x-small', labelspacing=0.1)
         ax.set_yscale('log', nonposy='clip')
         #ax.set_ylim((1., 50.))
         ax.set_xlabel(r'$z$')
         ax.set_ylabel(r'Effective bias $b$')
         ax.set_title(exp+' '+p.Prof.Lf.lineNameLatex)
         #
         # have the ticks in scientific format 
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_bias_z_lim_vs_brightgal_'+p.name+'.pdf')
         fig.clf()
         #plt.show()


         ###################################
         # nGalEff plots


         # nGalEff of LIM
         # full LIM
         f = lambda z: p.Prof.Lf.nGalEff(z)
         nGalEffLimFull = np.array(map(f, Z))
         # masked LIM
         def f(z): return p.Prof.Lf.nGalEff(z, lMin=0., lMax=fLMin(z))
         nGalEffLimLowL = np.array(map(f, Z))
         # detected galaxies, weighted by luminosity
         def f(z): return p.Prof.Lf.nGalEff(z, lMin=fLMin(z), lMax=np.inf)
         nGalEffLimHighL = np.array(map(f, Z))


         # nGalEff of detected galaxies = nGal
         # from luminosity cut
         def f(z): 
            return self.Prof.Lf.nGal(z, lMin=fLMin(z), lMax=np.inf)
         nGalEffGalL = np.array(map(f, Z))

         # test
         print "nGalEff"
         print nGalEffLimFull
         print nGalEffGalL


         # Compare nGalEff: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(Z, nGalEffLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(Z, nGalEffGalL, ls='-', label=r'Galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(r'$z$')
         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
         ax.set_title(exp+' '+p.Prof.Lf.lineNameLatex)
         #
         # have the ticks in scientific format 
         ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
         # to get more tick marks on the y axis
         ax.yaxis.set_major_locator(LogLocator(numticks=15))
         ax.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_neff_z_lim_vs_brightgal_'+p.name+'.pdf')
         fig.clf()
         #plt.show()


         ###################################
         # SNR plot


         # Compare nP with 1 for detected galaxies
         # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
         k = 0.1  # [h/Mpc]
         pLin = p.U.pLin(k, Z)   # [(Mpc/h)^3]
         snrGal = bGalM**2 * pLin / (bGalM**2 * pLin + 1./nGalEffGalL)

         # compare p2h/(pShot+Ndet) with 1 for LIM
         # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
         # Get the beam value at that k
         f = lambda z: p.U.psfF(k, fwhmPsf, z)   # [dimless]
         psf = np.array(map(f, Z))   # [dimless]
         # assume the Fourier mode is across the LOS
         spsf = 1.   #self.U.spectralPsfF(kPara, R, z)   # [dimless]
         w = psf * spsf

         # LIM full
         f = lambda z: p.p2h(k, z)
         p2h = np.array(map(f, Z))
         f = lambda z: p.pShot(k, z)
         pShot = np.array(map(f, Z))
         snrLimFull = p2h / (p2h + pShot + DetNoiseFid / w**2)
         # masked LIM
         f = lambda z: p.p2h(k, z, mMin=0., mMax=fMMin(z))
         p2hLimLowM = np.array(map(f, Z))
         f = lambda z: p.pShot(z, mMin=0., mMax=fMMin(z))
         pShotLimLowM = np.array(map(f, Z))
         snrLimLowM = p2hLimLowM / (p2hLimLowM + pShotLimLowM + DetNoiseFid / w**2)
         # detected galaxies, weighted by luminosity,
         # hence no detector noise
         f = lambda z: p.p2h(k, z, mMin=fMMin(z), mMax=np.inf)
         p2hLimHighM = np.array(map(f, Z))
         f = lambda z: p.pShot(z, mMin=fMMin(z), mMax=np.inf)
         pShotLimHighM = np.array(map(f, Z))
         # For the bright galaxies, one would generate the LIM
         # from the individually measured luminosities,
         # so there is no detector noise
         #snrLimHighM = p2hLimHighM / (p2hLimHighM + pShotLimHighM)
         snrLimHighM = bLimHighM**2 * pLin / (bLimHighM**2 * pLin + 1. / nGalEffLimHighL)

         # test
         print "SNR"
         print snrLimFull
         print snrGal


         # Compare nP: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(Z, snrLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(Z, snrGal, label=r'Galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
         ax.set_ylim((0., 1.1))
         ax.set_xlabel(r'$z$')
         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k_\perp=0.1 h/\text{Mpc})$')
         ax.set_title(exp+' '+p.Prof.Lf.lineNameLatex)
         #
         fig.savefig('./figures/p3d_rsd/'+exp+'_snr_z_lim_vs_brightgal_'+p.name+'.pdf')
         fig.clf()
         #plt.show()


         ###################################


      # Lmin plot
      ax.plot([], [], c='gray', ls='--', alpha=0.5, label=r'$L^\star$')
      #
      ax0.legend(loc=4, fontsize='x-small', labelspacing=0.1)
      ax0.set_yscale('log', nonposy='clip')
      ax0.set_xlabel(r'$z$')
      ax0.set_ylabel(r'$L_\text{min}$ [erg/s]')
      #
      # have the ticks in scientific format 
      ax0.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      ax0.yaxis.set_major_locator(LogLocator(numticks=15))
      ax0.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      fig0.savefig('./figures/p3d_rsd/summary_lmin_z_'+self.name+'.pdf')
      fig0.clf()
      #plt.show()

      # Mmin plot
      ax1.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax1.set_yscale('log', nonposy='clip')
      ax1.set_xlabel(r'$z$')
      ax1.set_ylabel(r'$m_\text{min}$ [$M_\odot/h$]')
      ax.set_ylim((1.e11, 1.e19))
      #
      # have the ticks in scientific format
      ax1.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the y axis
      ax1.yaxis.set_major_locator(LogLocator(numticks=15))
      ax1.yaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      fig1.savefig('./figures/p3d_rsd/summary_mmin_z_'+self.name+'.pdf')
      fig1.clf()
      #plt.show()






















   ##################################################################################
   # 2d matched filter for point sources in projected 2d maps

   def sigmaFluxMatchedFilter2d(self, detNoisePower, R, fwhmPsf, z, dz):
      '''Computes the uncertainty on the flux [Lsun / (Mpc/h)^2]
      of a point source in the projected 2d map,
      given the 3d detector noise power spectrum detNoisePower.
      detNoisePower: default unit [(Lsun/(Mpc/h)^2/sr/Hz)^2 (Mpc/h)^3],
      ie [Lsun^2/(Mpc/h)/sr^2/Hz^2].
      This assumes a matched filter is used, with the point source profile
      determined by the PSF of the instrument.
      R: spectral resolving power [dimless]
      fwhmPsf: [rad]
      dz: width of the redshift slice
      '''
      # precompute the RSD power spectrum at the requested redshift
      try:
         self.load(z=z)
      except:
         self.save(z=z)
         self.load(z=z)

      # volume conversion from p3d to p2d
      # in the thin redshift slice
      chi = self.U.bg.comoving_distance(z)   # [Mpc/h]
      dchi = dz * self.U.c_kms/self.U.hubble(z) # [Mpc/h]
      dVdOmega = chi**2 * dchi   # [(Mpc/h)^3/sr]

      # frequency range corresponding to the depth
      # of the thin redshift slice
      dnu = self.Prof.Lf.nuHz * dz / (1.+z)**2  # [Hz]


      def integrand(lnL):
         l = np.exp(lnL)
         kPerp = (l + 0.5) / chi
         k = kPerp
         mu = 0.

         # keep default unit [(Lsun/(Mpc/h)^2/sr/Hz)^2 (Mpc/h)^3]
         # ie [Lsun^2/(Mpc/h)/sr^2/Hz^2]
         pTot = self.pTotInt[z](k, mu)
         if pTot==0.:
            print "watch out ", k, mu, pTot

         psf = self.U.psfF(kPerp, fwhmPsf, z)   # [dimless]

         result = (psf/dnu)**2   # [1/Hz^2]
         result /= 2. * (psf**2 * pTot + detNoisePower) / dVdOmega   # [(Lsun/(Mpc/h)^2/sr)^-2 / sr]
         result *= l / (2.*np.pi)
         result *= l # because we integrate wrt lnL [(Lsun/(Mpc/h)^2)^-2]
         result *= 2.   # symmetry factor, since we reduce the integration domain
         return result

      # compute 2d integral
      # integration bounds: if modulus of k goes above self.kMax,
      # pTot will be zero, and sigmaMatchedFilter will be zero
      # if det noise is zero. --> avoid this!
      result = integrate.quad(integrand, np.log(self.kMin*chi), np.log(self.kMax*chi/np.sqrt(2.)), epsabs=0., epsrel=1.e-2)[0]
      result = 1. / np.sqrt(result) # [Lsun/(Mpc/h)^2]
      return result


   def sigmaLumMatchedFilter2d(self, detNoisePower, R, fwhmPsf, z, dz):
      '''Return the matched filter uncertainty in terms of luminosity [Lsun]
      '''
      # get the flux
      result = self.sigmaFluxMatchedFilter2d(detNoisePower, R, fwhmPsf, z, dz)   # [Lsun / (Mpc/h)^2]
      # convert to luminosity
      result *= 4.*np.pi * (1.+z)**2 * self.U.bg.comoving_distance(z)**2   # [Lsun]
      return result



   def plotSigmaLumMatchedFilter2d(self, exp='SPHEREx'):
      #if z is None:
      #   z = self.Z[0]

      # width of the redshift slice
      dz = 1.


      for z in self.Z:

         # volume conversion from p3d to p2d
         # in the thin redshift slice
         chi = self.U.bg.comoving_distance(z)   # [Mpc/h]
         dchi = dz * self.U.c_kms/self.U.hubble(z) # [Mpc/h]
         dVdOmega = chi**2 * dchi   # [(Mpc/h)^3/sr]


         # default power units, converted later when plotting
         DetNoisePower = np.logspace(np.log10(1.e-12), np.log10(1.e-4), 11, 10.)

         if exp=='SPHEREx':
            # SPHEREx specs
            R = 40.
            fwhmPsf = 6. * np.pi / (180.*3600.)
         elif exp=='CCAT-P':
            R = 1.
            fwhmPsf = 1. * np.pi / (180. * 60.)
         elif exp=='COMAP':
            R = 800.
            fwhmPsf = 3. * np.pi / (180. * 60.)
         elif exp=='CONCERTO':
            R = 300.
            fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]

         # min luminosity detectable [Lsun]: 5 sigma
         f = lambda detNoisePower: 5. * self.sigmaLumMatchedFilter2d(detNoisePower, R, fwhmPsf, z, dz)
         LMin = np.array(map(f, DetNoisePower))
         print "LMin", LMin

         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         ax.plot(x, LMin*self.Prof.Lf.convertLumUnit('cgs'))
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'$L_\text{min}$ [erg/s]')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_lmin_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Convert the minimum detected luminosity lMin
         # into a minimum detected halo mass
         # get the Kennicutt-Schmidt constant
         K = self.Prof.Sfr.kennicuttSchmidtConstant(z, self.Prof.Lf, alpha=1.)
         print "KS constant", K
         # use it to convert lMin to mMin
         f = lambda l: self.Prof.Sfr.massFromLum(l, z, K, alpha=1.)
         MMin = np.array(map(f, LMin))

         # Plot minimum halo mass detected
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         ax.plot(x, MMin)
         #
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'$M_\text{min}$ [$M_\odot/h$]')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_mmin_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         ###################################


         # fraction of mean intensity from undetected sources
         def f(lMin):
            result = self.Prof.Lf.meanIntensity(z, lMin=0., lMax=lMin)
            result /= self.Prof.Lf.meanIntensityInterp(z)
            return result
         fracMeanIntensity = np.array(map(f, LMin))
         #
         # fraction of shot noise from undetected sources
         def f(lMin):
            result = self.Prof.Lf.pShot(z, lMin=0., lMax=lMin)
            result /= self.Prof.Lf.pShotInterp(z)
            return result
         fracShotNoise = np.array(map(f, LMin))


         # fraction of 2h from undetected sources
         # at a fiducial k
         def f(mMin):
            k = 0.01
            result = self.p2h(k, z, mu=0., mMin=0., mMax=mMin)
            result /= self.p2h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac2h = np.array(map(f, MMin))
         #
         # fraction of 1h from undetected sources
         # at a fiducial k
         def f(mMin):
            k = 0.1
            result = self.p1h(k, z, mu=0., mMin=0., mMax=mMin)
            result /= self.p1h(k, z, mu=0., mMin=0., mMax=np.inf)
            return result
         frac1h = np.array(map(f, MMin))


         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         ax.plot(x, fracMeanIntensity, label=r'Mean intensity')
         ax.plot(x, fracShotNoise, label=r'Shot noise')
         ax.plot(x, frac2h, label=r'2h ($k=0.01h/$Mpc)')
         ax.plot(x, frac1h, label=r'1h ($k=0.1h/$Mpc)')
         #
         #ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'Fraction from undetected sources')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_fracundetected_detnoise_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         ###################################

         # Bias of LIM
         # full LIM
         bLimFull = self.Prof.Sfr.bEff(z, alpha=1.)
         # masked LIM
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=0., mMax=mMin)
         bLimLowM = np.array(map(f, MMin))
         # detected galaxies, weighted by luminosity
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=mMin, mMax=np.inf)
         bLimHighM = np.array(map(f, MMin))

         # Bias of detected galaxies
         # assume Ngal propto SFR, as for LIM
         f = lambda mMin: self.Prof.Sfr.bEff(z, alpha=1, mMin=mMin, mMax=np.inf)
         bGalM = np.array(map(f, MMin))



         # Compare bias: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.axhline(bLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, bGalM, label=r'Galaxies $(M>M_\text{min})$')
         #
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'Effective bias $b$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_bias_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare bias: unmasked vs masked LIM
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.axhline(bLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, bLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
         #
         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'Effective bias $b$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_bias_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare bias: bright gal, number or luminosity weighted
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.plot(x, bGalM, label=r'Number weighted galaxies $(M>M_\text{min})$')
         ax.plot(x, bLimHighM, '--', label=r'Luminosity-weighted galaxies $(M>M_\text{min})$')
         #
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'Effective bias $b$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_bias_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


   #      # Compare bias
   #      fig=plt.figure(0)
   #      ax=fig.add_subplot(111)
   #      #
   #      x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
   #      x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
   #      #
   #      ax.axhline(bLimFull, c='k', ls='--', label=r'LIM')
   #      ax.plot(x, bLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
   #      ax.plot(x, bGalM, label=r'Galaxies $(M>M_\text{min})$')
   #      ax.plot(x, bLimHighM, label=r'LIM-weighted galaxies $(M>M_\text{min})$')
   #      #
   #      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   #      ax.set_xscale('log', nonposx='clip')
   #      #ax.set_yscale('log', nonposy='clip')
   #      ax.set_xlim((np.min(x), np.max(x)))
   #      ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
   #      ax.set_ylabel(r'Effective bias $b$')
   #      #
   #      plt.show()


         ###################################

         # nGalEff of LIM
         # full LIM
         nGalEffLimFull = self.Prof.Lf.nGalEff(z)
         # masked LIM
         def f(lMin): return self.Prof.Lf.nGalEff(z, lMin=0., lMax=lMin)
         nGalEffLimLowL = np.array(map(f, LMin))
         # detected galaxies, weighted by luminosity
         def f(lMin): return self.Prof.Lf.nGalEff(z, lMin=lMin, lMax=np.inf)
         nGalEffLimHighL = np.array(map(f, LMin))


         # nGalEff of detected galaxies = nGal
         # from luminosity cut
         def f(lMin):
            return self.Prof.Lf.nGal(z, lMin=lMin, lMax=np.inf)
         nGalEffGalL = np.array(map(f, LMin))


         # Compare nGalEff: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.axhline(nGalEffLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, nGalEffGalL, label=r'Galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_neff_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare nGalEff: unmasked vs masked LIM
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.axhline(nGalEffLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, nGalEffLimLowL, label=r'Masked LIM $(L<L_\text{min})$')
         #
         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_neff_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare nGalEff: bright gal, number or luminosity weighted
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.plot(x, nGalEffGalL, label=r'Number weighted galaxies $(L>L_\text{min})$')
         ax.plot(x, nGalEffLimHighL, label=r'Luminosity weighted galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         ax.set_yscale('log', nonposy='clip')
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_neff_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()



   #      # Compare nGalEff
   #      fig=plt.figure(0)
   #      ax=fig.add_subplot(111)
   #      #
   #      x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
   #      x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
   #      #
   #      ax.axhline(nGalEffLimFull, c='k', ls='--', label=r'LIM')
   #      ax.plot(x, nGalEffLimLowL, label=r'Masked LIM $(L<L_\text{min})$')
   #      ax.plot(x, nGalEffGalL, label=r'Galaxies $(L>L_\text{min})$')
   #      ax.plot(x, nGalEffLimHighL, label=r'LIM-weighted galaxies $(L>L_\text{min})$')
   #      #
   #      ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   #      ax.set_xscale('log', nonposx='clip')
   #      ax.set_yscale('log', nonposy='clip')
   #      ax.set_xlim((np.min(x), np.max(x)))
   #      ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
   #      ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}$ [(Mpc/$h$)$^{-3}$]')
   #      #
   #      plt.show()

         ###################################

         # Compare nP with 1 for detected galaxies
         # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
         k = 0.1  # [h/Mpc]
         pLin = self.U.pLin(k, z)   # [(Mpc/h)^3]
         snrGal = bGalM**2 * pLin / (bGalM**2 * pLin + 1./nGalEffGalL)

         # compare p2h/(pShot+Ndet) with 1 for LIM
         # at k = 0.1 h/Mpc, where we are clearly 2-halo dominated
         # Get the beam value at that k
         psf = self.U.psfF(k, fwhmPsf, z)   # [dimless]
         # assume the Fourier mode is radial
         spsf = 1.   #self.U.spectralPsfF(kPara, R, z)   # [dimless]
         w = psf * spsf

         # LIM full
         snrLimFull = self.p2h(k,z) / (self.p2h(k,z) + self.pShot(z) + DetNoisePower / w**2)
         # masked LIM
         def f(mMin):
            return self.p2h(k, z, mMin=0., mMax=mMin)
         p2hLimLowM = np.array(map(f, MMin))
         def f(mMin):
            return self.pShot(z, mMin=0., mMax=mMin)
         pShotLimLowM = np.array(map(f, MMin))
         snrLimLowM = p2hLimLowM / (p2hLimLowM + pShotLimLowM + DetNoisePower / w**2)
         # detected galaxies, weighted by luminosity,
         # hence no detector noise
         def f(mMin):
            return self.p2h(k, z, mMin=mMin, mMax=np.inf)
         p2hLimHighM = np.array(map(f, MMin))
         def f(mMin):
            return self.pShot(z, mMin=mMin, mMax=np.inf)
         pShotLimHighM = np.array(map(f, MMin))
         # For the bright galaxies, one would generate the LIM
         # from the individually measured luminosities,
         # so there is no detector noise
         #snrLimHighM = p2hLimHighM / (p2hLimHighM + pShotLimHighM)
         snrLimHighM = bLimHighM**2 * pLin / (bLimHighM**2 * pLin + 1. / nGalEffLimHighL)

         # Compare nP: unmasked LIM vs bright galaxies
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.plot(x, snrLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, snrGal, label=r'Galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.1))
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k=0.1 h/\text{Mpc})$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_snr_detnoise_lim_vs_brightgal_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare nP: unmasked vs masked LIM
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.plot(x, snrLimFull, c='k', ls='--', label=r'LIM')
         ax.plot(x, snrLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
         #
         ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.1))
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k=0.1 h/\text{Mpc})$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_snr_detnoise_masked_vs_unmasked_lim_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


         # Compare nP: bright gal, number or luminosity weighted
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
         x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
         #
         ax.plot(x, snrGal, label=r'Number weighted galaxies $(L>L_\text{min})$')
         ax.plot(x, snrLimHighM, label=r'Luminosity weighted galaxies $(L>L_\text{min})$')
         #
         ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
         ax.set_xscale('log', nonposx='clip')
         #ax.set_yscale('log', nonposy='clip')
         ax.set_ylim((0., 1.1))
         ax.set_xlim((np.min(x), np.max(x)))
         ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
         #ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
         ax.set_ylabel(r'$\text{SNR}_{P_\text{lin}}(k=0.1 h/\text{Mpc})$')
         #
         fig.savefig('./figures/p3d_rsd/2d_'+exp+'_snr_detnoise_brightgal_weighting_'+self.name+'_z'+str(z)+'.pdf')
         plt.clf()
         #plt.show()


   #      fig=plt.figure(0)
   #      ax=fig.add_subplot(111)
   #      #
   #      x = DetNoisePower * self.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
   #      x /= dVdOmega  # convert 3d power to 2d power [(Jy/sr)^2 * sr]
   #      #
   #      ax.plot(x, nPLimFull, c='k', ls='--', label=r'LIM')
   #      ax.plot(x, nPLimLowM, label=r'Masked LIM $(M<M_\text{min})$')
   #      #ax.plot(x, p2hLimM / pShotLimM, 'k--')
   #      #ax.plot(x, p2hLimM / DetNoisePower, 'ko')
   #      ax.plot(x, nPLimHighM, label=r'LIM-weighted galaxies $(M>M_\text{min})$')
   #      #
   #      ax.plot(x, nPGal, label=r'Galaxies $(L>L_\text{min})$')
   #      #
   #      ax.legend(loc='center left', fontsize='x-small', labelspacing=0.1)
   #      ax.set_xscale('log', nonposx='clip')
   #      ax.set_yscale('log', nonposy='clip')
   #      ax.set_xlim((np.min(x), np.max(x)))
   #      ax.set_xlabel(r'Detector noise power [Jy$^2$/sr]')
   #      ax.set_ylabel(r'$\bar{n}_\text{gal}^\text{eff}\ b^2 P_\text{lin}$')
   #      #
   #      plt.show()

      



##################################################################################
##################################################################################
##################################################################################
##################################################################################

class P3dRsdCross(P3dRsdAuto):
   
   def __init__(self, U, Prof, Prof2, MassFunc, r=1., name="", nProc=1):
      # copy classes
      self.U = U
      #self.IHaloModel = IHaloModel
      self.MassFunc = MassFunc
      self.Prof = Prof
      self.Prof2 = Prof2
      self.r = r
      self.nProc = nProc
      self.name = str(self.Prof) + '_' + str(self.Prof2) + name
      
      # mass bounds for integrals
      self.mMin = np.max([self.Prof.mMin, self.Prof2.mMin, self.MassFunc.mMin])
      self.mMax = min([self.Prof.mMax, self.Prof2.mMax, self.MassFunc.mMax])


      # k moduli, to precompute
      self.K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
      self.nK = len(self.K)
      self.kMin = np.min(self.K)
      self.kMax = np.max(self.K)
         
      # mu values, to precompute
      self.Mu = np.linspace(0., 1., 11)
      self.nMu = len(self.Mu)
      self.muMin = np.min(self.Mu)
      self.muMax = np.max(self.Mu)

      # redshifts to evaluate, from the luminosity function
      self.Z = self.Prof.Lf.Z
   
      # create folder if needed
      directory = "./output/p3d_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)


   ##################################################################################
   # Power spectrum ingredients


   def p1h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.u(k, m, z, mu) * self.Prof2.u(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result



   def fEff2(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Effective growth rate of structure
      Converges to f when k-->0
      Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof2.uMat(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      result *= self.U.bg.scale_independent_growth_rate(z)
      # fix the halo model total mass problem
      result /= self.correctionFactorI11(z, mMin=mMin, mMax=mMax)
      return result



   def bEff2(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof2.u(k, m, z, mu)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      if self.Prof.use_correction_factor==1:
         result /= self.correctionFactorI11(z, mMin=mMin, mMax=mMax)
      return result


   def p2h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes RSD: Kaiser effect and FOGs
      '''
      b1 = self.bEff(k, z, mu=mu, mMin=mMin, mMax=mMax)
      b2 = self.bEff2(k, z, mu=mu, mMin=mMin, mMax=mMax)
      f1 = self.fEff(k, z, mu=mu, mMin=mMin, mMax=mMax)
      f2 = self.fEff2(k, z, mu=mu, mMin=mMin, mMax=mMax)
      return (b1 + f1*mu**2) * (b2 + f2*mu**2) * self.U.pLin(k, z)


   def pShot(self, z, mMin=0., mMax=np.inf):
      '''(RSD has no effect on the shot noise power spectrum)
      '''
      result = self.r * np.sqrt(self.Prof.pShot(z) * self.Prof2.pShot(z))
      # build up with mass,
      # assuming phi(L|m) = N(m) * phi(L)
      result *= self.Prof.Sfr.sfrdForInterp(z, mMin=mMin, mMax=mMax) / self.Prof.Sfr.sfrd(z)
      return result

   
   ##################################################################################

#   def plotCorrCoeff(self, Z=None):
#      if Z is None:
#         Z = self.Z
#      P12 = self
#      P11 = P3dRsdAuto(self.U, self.Prof, self.MassFunc)
#      P22 = P3dRsdAuto(self.U, self.Prof2, self.MassFunc)
#
#      #K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 7, 10.)
#      K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 51, 10.)
#
#      Mu = np.array([0., 0.5, 1.])
#      Ls = ['-', '--', ':']
#
#      fig=plt.figure(0)
#      ax=fig.add_subplot(111)
#      #
#      for z in Z:
#         plot=ax.plot([],[], ls='-', label=r'$z=$'+str(round(z, 2)))
#         c = plot[0].get_color()
#         for iMu in range(len(Mu)):
#            mu = Mu[iMu]
#            f = lambda k: P12.pTot(k, z, mu)
#            p12 = np.array(map(f, K))
#            print "done 12"
#            print p12
#            f = lambda k: P11.pTot(k, z, mu)
#            p11 = np.array(map(f, K))
#            print "done 11"
#            print p11
#            f = lambda k: P22.pTot(k, z, mu)
#            p22 = np.array(map(f, K))
#            print "done 22"
#            print p22
#            #
#            plot=ax.plot(K, p12 / np.sqrt(p11 * p22), c=c, ls=Ls[iMu])
#      #
#      for iMu in range(len(Mu)):
#         ax.plot([], [], c='k', ls=Ls[iMu], label=r'$\mu=$'+str(Mu[iMu]))
#      #
#      ax.axhline(self.r, c='gray', alpha=0.3)
#      #
#      ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
#      ax.set_xlim((0.01, 1.e2))
#      ax.set_ylim((0., 1.))
#      ax.set_xscale('log', nonposx='clip')
#      ax.set_xlabel(r'$k$ [$h/$Mpc]')
#      ax.set_ylabel(r'$r_{12}(k, \mu, z)$')
#      ax.set_title(self.Prof.Lf.lineNameLatex+r' -- '+self.Prof2.Lf.lineNameLatex)
#      #
#      fig.savefig('./figures/p3d_rsd/r_'+str(self.Prof)+'_'+str(self.Prof2)+'.pdf', bbox_inches='tight')
#      fig.clf()
#      #plt.show()



   def plotCorrCoeff(self, Z=None):
      if Z is None:
         Z = self.Z
      P12 = self
      P11 = P3dRsdAuto(self.U, self.Prof, self.MassFunc)
      P22 = P3dRsdAuto(self.U, self.Prof2, self.MassFunc)

      #K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 5, 10.)
      K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 51, 10.)

      Mu = np.array([0., 0.5, 1.])
      Ls = ['-', '--', ':']

      # plot r
      fig=plt.figure(0)
      ax=fig.add_subplot(111)

      # plot 1-r**2
      fig2=plt.figure(1)
      ax2=fig2.add_subplot(111)



      # compute the correlation coefficients
      for z in Z:
         plot=ax.plot([],[], ls='-', label=r'$z=$'+str(round(z, 2)))
         plot2=ax2.plot([],[], ls='-', label=r'$z=$'+str(round(z, 2)))
         c = plot[0].get_color()
         for iMu in range(len(Mu)):
            mu = Mu[iMu]
            f = lambda k: P12.pTot(k, z, mu)
            p12 = np.array(map(f, K))
            print "done 12"
            print p12
            f = lambda k: P11.pTot(k, z, mu)
            p11 = np.array(map(f, K))
            print "done 11"
            print p11
            f = lambda k: P22.pTot(k, z, mu)
            p22 = np.array(map(f, K))
            print "done 22"
            print p22
            #
            r = p12 / np.sqrt(p11 * p22)
            plot=ax.plot(K, r, c=c, ls=Ls[iMu])
            plot2=ax2.plot(K, 1. - r**2, c=c, ls=Ls[iMu])
      

      for iMu in range(len(Mu)):
         ax.plot([], [], c='gray', ls=Ls[iMu], label=r'$\mu=$'+str(Mu[iMu]))
         ax2.plot([], [], c='gray', ls=Ls[iMu], label=r'$\mu=$'+str(Mu[iMu]))

      # expected r corr coeff in the shot noise regime
      ax.axhline(self.r, c='r', alpha=0.1)
      ax2.axhline(1. - self.r**2, c='r', alpha=0.1)
      

      # Clean up r plot
      ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
      ax.set_xlim((0.01, 1.e2))
      ax.set_ylim((0., 1.))
      ax.set_xscale('log', nonposx='clip')
      ax.set_xlabel(r'$k$ [$h/$Mpc]')
      ax.set_ylabel(r'$r_{12}(k, \mu, z)$')
      ax.set_title(self.Prof.Lf.lineNameLatex+r' -- '+self.Prof2.Lf.lineNameLatex)
      #
      # have the ticks in scientific format 
      ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15))
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      fig.savefig('./figures/p3d_rsd/r_'+str(self.Prof)+'_'+str(self.Prof2)+'.pdf', bbox_inches='tight')
      fig.clf()
      #plt.show()

      # Clean up 1-r**2 plot
      ax2.legend(loc=4, fontsize='x-small', labelspacing=0.1)
      ax2.set_xlim((0.01, 1.e2))
      ax2.set_ylim((0., 1.))
      ax2.set_xscale('log', nonposx='clip')
      ax2.set_xlabel(r'$k$ [$h/$Mpc]')
      ax2.set_ylabel(r'$1 - r^2_{12}(k, \mu, z)$')
      ax2.set_title(self.Prof.Lf.lineNameLatex+r' -- '+self.Prof2.Lf.lineNameLatex)
      #
      # have the ticks in scientific format 
      ax2.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
      # to get more tick marks on the x axis
      ax2.xaxis.set_major_locator(LogLocator(numticks=15))
      ax2.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10)))
      #
      fig2.savefig('./figures/p3d_rsd/r2_'+str(self.Prof)+'_'+str(self.Prof2)+'.pdf', bbox_inches='tight')
      fig2.clf()
      #plt.show()















