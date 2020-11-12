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
      self.K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
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
      directory = "./output/p_rsd/"
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
         result = self.MassFunc.massfunc(m, 1./(1.+z))
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
      result = self.Prof.Pshot(z)
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

   def sFOverFFisher(self, z, R, fwhmPsf, fSky, dz):
      '''Relative uncertainty on growth of structure f
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
         result = (2.*f* mu**2 * (b+f*mu**2) * self.U.pLin(k, z) / pTot)**2
         result *= kPerp / (2.*np.pi)**2
         result *= kPerp * kPara   # because int wrt ln
         result *= volume / 2.
         result *= 4.   # symmetry factor, since we reduce the integration domain
         return result
      
      result = integrate.dblquad(integrand, np.log(self.kMin), np.log(kPerpMax), lambda x: np.log(self.kMin), lambda x: np.log(kParaMax), epsabs=0., epsrel=1.e-2)[0]
      result = 1./np.sqrt(result)
      return result


   ##################################################################################



   def plotP(self, z=None, mu=0.):
      if z is None:
         z = self.Z[0]

      # Plin
      f = lambda k: self.U.pLin(k, z)
      Plin = np.array(map(f, self.K))
      # P1h
      f = lambda k: self.p1h(k, z)
      P1h = np.array(map(f, self.K))
      # P2h
      f = lambda k: self.p2h(k, z)
      P2h = np.array(map(f, self.K))
      # noise bias
      f = lambda k: self.pShot(z)
      Pnoise = np.array(map(f, self.K))
      # Ptot
      P = P1h + P2h + Pnoise

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, P, 'k', lw=4, label=r'$P_\text{tot}$')
      ax.loglog(self.K, P2h, 'b-', lw=2, label=r'$P_\text{2h}$')
      #ax.loglog(self.K, self.Prof.bias(1./(1.+z))**2 * Plin, 'b--', lw=2, label=r'$b_\text{eff}^2 P_\text{lin}$')
      ax.loglog(self.K, Plin, 'k--', lw=2, label=r'$P_\text{lin}$')
      ax.loglog(self.K, P1h, 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.K, Pnoise, 'g-', lw=2, label=r'$P_\text{noise}$')

      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$P(k)$')
      #path = "./figures/pn3d/p_"+self.name+"z_"+str(z)+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      plt.show()



   def plotPMuDpdce(self, z=None):
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
            fSky = 100.*(np.pi/180.)**2 / (4.*np.pi)
            fwhmPsf = 6. * np.pi/(180.*3600.)   # [arcsec] to [rad]
            kFPerp = self.U.kFPerp(z, fSky)
            kMaxPerp = self.U.kMaxPerpPsf(fwhmPsf, z)
            print "perp", kFPerp, kMaxPerp
            #ax.hlines(5.e3, xmin=kFPerp, xmax=kMaxPerp, colors=plot[0].get_color())
            ax.axvspan(kFPerp, kMaxPerp, fc=plot[0].get_color(), ec=None, alpha=0.1)
         if mu==1.:
            # k_para
            R = 40.
            dz = 0.5
            kFPara = self.U.kFPara(z, dz=dz)
            kMaxPara = self.U.kMaxParaSpectroRes(R, z)
            print "para", kFPara, kMaxPara
            #ax.hlines(5.e3, xmin=kFPara, xmax=kMaxPara, colors=plot[0].get_color())
            ax.axvspan(kFPara, kMaxPara, fc=plot[0].get_color(), ec=None, alpha=0.1)
      #
      #
      # to get more tick marks on the x axis
      ax.xaxis.set_major_locator(LogLocator(numticks=15)) #(1)
      ax.xaxis.set_minor_locator(LogLocator(numticks=15,subs=np.arange(2,10))) #(2)
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k, \mu, z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      path = "./figures/p3d_rsd/p_mudpdce_"+self.name+"z_"+str(z)+".pdf"
      fig.savefig(path, bbox_inches='tight')

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


   def compareP(self, ps=None):
      if ps is None:
         ps = [self]

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for p in ps:
         for z in p.Prof.Lf.Z:
            if z<=5:
               f = lambda k: p.pTot(k, z)
               pTot = np.array(map(f, p.K)) * p.Prof.Lf.convertPowerSpectrumUnit('Jy/sr')
               plot=ax.loglog(p.K, pTot, label=r'$z=$'+str(round(z, 2))+' '+p.Prof.Lf.nameLatex)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k,z)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      #
      #path = './figures/pn_3d/p3d_'+lfHa[key].name+'.pdf'
      #fig.savefig(path, bbox_inches='tight')
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
         ax.fill_between(K, P[iM,:], P[iM+1,:], facecolor=plt.cm.YlOrRd(1.*iM/(len(M)-1.)), edgecolor='', label=r'$M\leqslant$'+floatSciForm(m, round=1)+r' $M_\odot/h$')

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
         ax.plot(K, P[iM,:], ls='--', c=plt.cm.YlOrRd((iM+1.)/(len(M)-1.)), label=floatSciForm(M[iM], round=1)+r'$\leqslant M<$'+floatSciForm(M[iM+1], round=1)+r' $M_\odot/h$')
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
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'$P(k,z, \mu)$ [(Jy/sr)$^2$ (Mpc/$h$)$^3$]')
      ax.set_ylim((4.e3, 5.e7))
      ax.set_title(r'$\mu='+str(mu)+'$')
      #ax.set_title(self.Prof.Lf.nameLatex)
      #
      path = './figures/p3d_rsd/p1h2hshot_'+self.name+'_mu'+str(mu)+'.pdf'
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def plotRequiredAreaToDetectF(self):
      '''For a given spectral resolution R,
      survey depth Delta z, what sky area is required
      to give a $10\%$ measurement of the growth rate f?
      '''
      RR = np.array([40., 150., 300.])
      #Z = np.linspace(0., 7., 11)
      Z = self.Z.copy()

      def fSkyReq(z, R, dz=1., fwhmPsf=6.*np.pi/(180.*3600.), target=0.1):
         sFOverF = self.sFOverFFisher(z, R, fwhmPsf, 1., dz)
         result = (sFOverF / target)**2
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
      ax.axhline(100. * (np.pi/180.)**2 / (4.*np.pi), ls='--', label=r'SPHEREx deep fields')
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
      fig.savefig('./figures/p3d_rsd/fsky_tradeoff_f.pdf', bbox_inches='tight')
      plt.show()



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
      directory = "./output/p_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)


   ##################################################################################
   # Power spectrum ingredients


   def p1h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, 1./(1.+z))
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
      result = self.r * np.sqrt(self.Prof.Pshot(z) * self.Prof2.Pshot(z))
      # build up with mass,
      # assuming phi(L|m) = N(m) * phi(L)
      result *= self.Prof.Sfr.sfrdForInterp(z, mMin=mMin, mMax=mMax) / self.Prof.Sfr.sfrd(z)
      return result

   
   ##################################################################################

   def plotCorrCoeff(self, Z=None):
      if Z is None:
         Z = self.Z
      P12 = self
      P11 = P3dRsdAuto(self.U, self.Prof, self.MassFunc)
      P22 = P3dRsdAuto(self.U, self.Prof2, self.MassFunc)

      K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 51, 10.)

      Mu = np.array([0., 0.5, 1.])
      Ls = ['-', '--', ':']

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for z in Z:
         plot=ax.plot([],[], ls='-', label=r'$z=$'+str(round(z, 2)))
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
            plot=ax.plot(K, p12 / np.sqrt(p11 * p22), c=c, ls=Ls[iMu])
      #
      for iMu in range(len(Mu)):
         ax.plot([], [], c='k', ls=Ls[iMu], label=r'$\mu=$'+str(Mu[iMu]))
      #
      ax.axhline(self.r, c='gray', alpha=0.3)
      #
      ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
      ax.set_xlim((0.01, 1.e2))
      ax.set_ylim((0., 1.))
      ax.set_xscale('log', nonposx='clip')
      ax.set_xlabel(r'$k$ [$h/$Mpc]')
      ax.set_ylabel(r'$r_{12}(k, \mu, z)$')
      plt.show()

