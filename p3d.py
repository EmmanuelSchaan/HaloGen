from headers import *

##################################################################################
##################################################################################

class P3dAuto(object):
   
   def __init__(self, U, MassFunc, Prof, pNoise=lambda k,z: 0., name="", Vs=1., save=False):
      # copy classes
      self.U = U
      #self.IHaloModel = IHaloModel
      self.MassFunc = MassFunc
      self.Prof = Prof
      self.pNoise = pNoise
      self.Vs = Vs   # in (Mpc/h)^3
      self.name = str(self.Prof) + name

      # mass bounds for integrals
      self.mMin = max(self.Prof.mMin, self.MassFunc.mMin)
      self.mMax = min(self.Prof.mMax, self.MassFunc.mMax)
      self.nM = 101  # for mass integrals

#      # k moduli, to precompute
#      self.K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
#      self.nK = len(self.K)
#      self.kMin = np.min(self.K)
#      self.kMax = np.max(self.K)
      
      # values of k to evaluate
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      self.Ke = np.genfromtxt("./input/K.txt") # edges of the bins for k
      self.dK = np.genfromtxt("./input/dK.txt") # widths of the bins for k

      # redshifts to evaluate
      self.Z = np.linspace(0., 10., 51)
      #self.Z = np.linspace(0., 10., 5)
      #self.Z = np.linspace(0., 10., 3)
   
      # create folders if needed
      self.pathOut = "./output/p3d/p3dauto_"+self.name+"/"
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      self.pathFig = "./figures/p3d/p3dauto_"+self.name+"/"
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)

      # power spectrum
      if save:
         self.saveP()
      self.loadP()
      

   ##################################################################################
   ##################################################################################
   # Precompute P3d in batch


   def massIntegrals(self, z, fOfkmz):
      '''Computes
      int dm n(m,z) fOfkm(k,m)
      at the z requested
      for all k in self.K
      '''
      m = np.logspace(np.log10(self.mMin), np.log10(self.mMax), self.nM, 10.)

      fVect = np.vectorize(self.MassFunc.massfunc)
      integrand = fVect(m, z)
      
      fVect = np.vectorize(fOfkmz)
      f = lambda k: fVect(k, m, z)
      integrand = integrand[None, :] * np.array(map(f, self.K))

      result = np.trapz(integrand, m, axis=-1)
      return result


   def computeP1hBare(self, z):
      f = lambda k,m,z: self.Prof.u(k, m, z)**2
      return self.massIntegrals(z, f)

   def computeP1hCounterTerm(self, z):
      result = self.mMin 
      result *= self.MassFunc.nMinForInterp(z, self.mMin, self.mMax)
      f = lambda k: self.Prof.u(k, self.mMin, z)**2
      result *= np.array(map(f, self.K))
      return result
   

   def computeBEffBare(self, z):
      f = lambda k,m,z: self.MassFunc.b1(m, z) * self.Prof.u(k, m, z)
      return self.massIntegrals(z, f)

   def computeBEffCounterTerm(self, z):
      result = self.mMin
      result *= self.MassFunc.nMinForInterp(z, self.mMin, self.mMax) 
      result *= self.MassFunc.b1MinForInterp(z, self.mMin, self.mMax)
      f = lambda k: self.Prof.u(k, self.mMin, z)
      result *= np.array(map(f, self.K))
      return result
      
   def computeBEffCorrected(self, z):
      result = self.computeBEffBare(z)
      result *= self.MassFunc.correctionFactorBias(z)
      return result


   ##################################################################################

   def saveP(self):
      print "Precomputing p3d "+self.name
      #print "FAST VERSION"
      tStart = time()
      nZ = len(self.Z)
      nK = len(self.K)
      
      # power spectra
      P1h = np.zeros((nK, nZ))
      P2h = np.zeros((nK, nZ))
      PNoise = np.zeros((nK, nZ))
      #dP = np.zeros((nK, nZ))
      P1hBare = np.zeros((nK, nZ))
      P1hCounterTerm = np.zeros((nK, nZ))
      BEffBare = np.zeros((nK, nZ))
      BEffCounterTerm = np.zeros((nK, nZ))
      BEffCorrected = np.zeros((nK, nZ))


      for iZ in range(nZ):
         z = self.Z[iZ]
         P1hBare[:,iZ] = self.computeP1hBare(z)
         P1hCounterTerm[:,iZ] = self.computeP1hCounterTerm(z)
         #
         BEffBare[:,iZ] = self.computeBEffBare(z)
         BEffCounterTerm[:,iZ] = self.computeBEffCounterTerm(z)
         BEffCorrected[:,iZ] = self.computeBEffCorrected(z)
         print "done z "+str(iZ)+" of "+str(nZ)

      P1h = P1hBare + P1hCounterTerm
      #
      BEff = BEffBare + BEffCounterTerm
      pLinVect = np.vectorize(self.U.pLin)
      f = lambda k: pLinVect(k, self.Z)
      PLin = np.array(map(f, self.K))
      #
      P2hBare = BEffBare**2 * PLin
      P2hCounterTerm = BEffCounterTerm**2 * PLin
      P2h = BEff**2 * PLin
      P2hCorrected = BEffCorrected**2 * PLin

      pNoiseVect = np.vectorize(self.pNoise)
      f = lambda k: pNoiseVect(k, self.Z)
      PNoise = np.array(map(f, self.K))

      # save power spectra
      path = self.pathOut+"p3d_"+self.name
      np.savetxt(path+"_z.txt", self.Z)
      np.savetxt(path+"_k.txt", self.K)
      np.savetxt(path+"_1h.txt", P1h)
      np.savetxt(path+"_2h.txt", P2h)
      np.savetxt(path+"_noise.txt", PNoise)
      #np.savetxt(path+"_dPdd.txt", dP)
      np.savetxt(path+"_1hbare.txt", P1hBare)
      np.savetxt(path+"_1hcounterterm.txt", P1hCounterTerm)
      np.savetxt(path+"_2hbare.txt", P2hBare)
      np.savetxt(path+"_2hcounterterm.txt", P2hCounterTerm)
      np.savetxt(path+"_2hcorrected.txt", P2hCorrected)

      tStop = time()
      print "took ", (tStop-tStart)/60., "min"
      return


   def loadP(self):
      # read the precomputed power spectra
      path = self.pathOut+"p3d_"+self.name
      Z = np.genfromtxt(path+"_z.txt")
      K = np.genfromtxt(path+"_k.txt")
      self.P1h = np.genfromtxt(path+"_1h.txt")
      self.P2h = np.genfromtxt(path+"_2h.txt")
      self.PNoise = np.genfromtxt(path+"_noise.txt")
      self.P = self.P1h + self.P2h
      self.PTot = self.P1h + self.P2h + self.PNoise
      #self.dP = np.genfromtxt(path+"_dPdd.txt")
      self.P1hBare = np.genfromtxt(path+"_1hbare.txt")
      self.P1hCounterTerm = np.genfromtxt(path+"_1hcounterterm.txt")
      self.P2hBare = np.genfromtxt(path+"_2hbare.txt")
      self.P2hCounterTerm = np.genfromtxt(path+"_2hcounterterm.txt")
      self.P2hCorrected = np.genfromtxt(path+"_2hcorrected.txt")
      
      # interpolate them
      self.p1hInt = interp2d(K, Z, self.P1h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p2hInt = interp2d(K, Z, self.P2h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.pNoiseInt = interp2d(K, Z, self.PNoise.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.pInt = interp2d(K, Z, self.P.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.pTotInt = interp2d(K, Z, self.PTot.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      #self.fdPinterp = interp2d(K, Z, self.dP.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p1hBareInt = interp2d(K, Z, self.P1hBare.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p1hCounterTermInt = interp2d(K, Z, self.P1hCounterTerm.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p2hBareInt = interp2d(K, Z, self.P2hBare.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p2hCounterTermInt = interp2d(K, Z, self.P2hCounterTerm.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.p2hCorrectedInt = interp2d(K, Z, self.P2hCorrected.transpose(), kind='linear', bounds_error=False, fill_value=0.)
   

   ##################################################################################
   ##################################################################################
   # 1-halo term

   def p1hBare(self, k, z, mMin=0., mMax=np.inf):
      '''P1h, without the mass counterterm
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.u(k, m, z)**2
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      if not np.isfinite(result):
         result = 0.
      return result


   def p1hCounterTerm(self, k, z, mMin=0., mMax=np.inf):
      '''Mass counter term for P1h,
      only evaluated at the fiducial mass cut
      '''
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      # counter term
      result = self.MassFunc.nMinForInterp(z, mMin, mMax) 
      result *= mMin
      result *= self.Prof.u(k, mMin, z)**2
      return result
      

   def p1h(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of the bare and counter term
      '''
      result = self.p1hBare(k, z, mMin, mMax)
      result += self.p1hCounterTerm(k, z, mMin, mMax)
      return result


   def plotP1hMMinDependence(self, k=0.01):
      '''Select a scale k in [h/Mpc]
      '''
      MMin = np.logspace(np.log10(1.e6), np.log10(1.e16), 101, 10.) # [Msun/h]
      Z = np.array([0.001, 1., 2., 3., 4.])

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot([], [], '-', c='gray', label=r'Bare + Counter term')
      ax.plot([], [], '--', c='gray', label=r'Bare')
      ax.plot([], [], '-.', c='gray', label=r'Counter term')
      #
      for iZ in range(len(Z)):
         z = Z[iZ]
         # bare + counter term
         f = lambda mMin: self.p1h(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         plot=ax.plot(MMin, y, c=plt.cm.cool(iZ/(len(Z)-1.)), ls='-', label=r'$z=$'+str(int(z)))
         # bare  
         f = lambda mMin: self.p1hBare(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         ax.plot(MMin, y, c=plot[0].get_color(), ls='--')
         # counter term  
         f = lambda mMin: self.p1hCounterTerm(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         ax.plot(MMin, y, c=plot[0].get_color(), ls='-.')
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((MMin.min(), MMin.max()))
      #ax.set_ylim((1.e-5, 1.e5))
      ax.set_xlabel(r'Minimum halo mass $m_\text{min}$ [$M_\odot/h$]')
      ax.set_ylabel(r'$P^{1-\text{halo}}(k='+str(round(k, 3))+r' h/\text{Mpc})$')
      #
      fig.savefig(self.pathFig+"p1h_mmin_dependence.pdf", bbox_inches='tight')
      fig.clf()
      #plt.show()


   ##################################################################################
   # 2-halo term

   def bEffBare(self, k, z, mMin=0., mMax=np.inf):
      '''Without mass counter term
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof.u(k, m, z)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result


   def bEffCounterTerm(self, k, z, mMin=0., mMax=np.inf):
      '''Mass counter term,
      such that the matter bias is indeed unity.
      '''
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      # counter term
      result = self.MassFunc.nMinForInterp(z, mMin, mMax) 
      result *= self.MassFunc.b1MinForInterp(z, mMin, mMax)
      result *= mMin
      result *= self.Prof.u(k, mMin, z)
      return result


   def bEff(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of the bare and counter term
      '''
      result = self.bEffBare(k, z, mMin, mMax)
      result += self.bEffCounterTerm(k, z, mMin, mMax)
      return result

   def bEffCorrected(self, k, z, mMin=0., mMax=np.inf):
      '''With the correction factor
      '''
      result = self.bEffBare(k, z, mMin, mMax)
      result *= self.MassFunc.correctionFactorBias(z)
      return result


   def p2hBare(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEffBare(k, z, mMin=mMin, mMax=mMax)**2
      result *= self.U.pLin(k, z)
      return result

   def p2hCounterTerm(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEffCounterTerm(k, z, mMin=mMin, mMax=mMax)**2
      result *= self.U.pLin(k, z)
      return result


   def p2h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEff(k, z, mMin=mMin, mMax=mMax)**2 
      result *= self.U.pLin(k, z)
      return result


   def plotP2hMMinDependence(self, k=0.01):
      '''Select a scale k in [h/Mpc]
      '''
      MMin = np.logspace(np.log10(1.e6), np.log10(1.e16), 101, 10.) # [Msun/h]
      Z = np.array([0.001, 1., 2., 3., 4.])

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot([], [], '-', c='gray', label=r'Bare + Counter term')
      ax.plot([], [], '--', c='gray', label=r'Bare')
      ax.plot([], [], '-.', c='gray', label=r'Counter term')
      #
      for iZ in range(len(Z)):
         z = Z[iZ]
         # bare + counter term
         f = lambda mMin: self.p2h(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         plot=ax.plot(MMin, y, c=plt.cm.cool(iZ/(len(Z)-1.)), ls='-', label=r'$z=$'+str(int(z)))
         # bare  
         f = lambda mMin: self.p2hBare(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         ax.plot(MMin, y, c=plot[0].get_color(), ls='--')
         # counter term  
         f = lambda mMin: self.p2hCounterTerm(k, z, mMin=mMin)
         y = np.array(map(f, MMin))
         ax.plot(MMin, y, c=plot[0].get_color(), ls='-.')
      #
      ax.legend(loc=4, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((MMin.min(), MMin.max()))
      #ax.set_ylim((1.e1, 1.e5))
      ax.set_xlabel(r'Minimum halo mass $m_\text{min}$ [$M_\odot/h$]')
      ax.set_ylabel(r'$P^{2-\text{halo}}(k='+str(round(k, 3))+r' h/\text{Mpc})$')
      #
      fig.savefig(self.pathFig+"p2h_mmin_dependence.pdf", bbox_inches='tight')
      fig.clf()
      #plt.show()


   def plotP2hCorrections(self):
      '''2-halo term: compare the bare, counter term and corrected values
      as a function of k.
      '''
      Z = np.array([0., 1., 3., 5.])

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot([], [], '-', c='gray', label=r'Bare + Counter term')
      ax.plot([], [], '--', c='gray', label=r'Bare')
      ax.plot([], [], '-.', c='gray', label=r'Counter term')
      ax.plot([], [], ':', c='gray', label=r'Corrected')
      #
      for iZ in range(len(Z)):
         z = Z[iZ]

         yBare = self.p2hBareInt(self.K, z)
         yCounterTerm = self.p2hCounterTermInt(self.K, z)
         yFull = self.p2hInt(self.K, z)
         yCorrected = self.p2hCorrectedInt(self.K, z)
         
         plot=ax.plot([], [], c=plt.cm.cool(iZ/(len(Z)-1.)), ls='-', label=r'$z=$'+str(int(z)))
         ax.plot(self.K, yBare / yFull, c=plot[0].get_color(), ls='--')
         ax.plot(self.K, yCounterTerm / yFull, c=plot[0].get_color(), ls='-.')
         ax.plot(self.K, yCorrected / yFull, c=plot[0].get_color(), ls=':')
      #
      ax.axhline(1., c='gray', ls='-')
      #
      ax.legend(loc=3, fontsize='x-small', labelspacing=0.2)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      #ax.set_xlim((MMin.min(), MMin.max()))
      #ax.set_ylim((1.e-3, 1.2))
      ax.set_xlabel(r'$k$ [$h$/Mpc]')
      ax.set_ylabel(r'Fraction of $P^\text{2-halo}(k, z)$')
      #
      fig.savefig(self.pathFig + "counterterms_2halo_"+self.name+".pdf", bbox_inches='tight')
      fig.clf()
      #plt.show()

   ##################################################################################


   def plotPInt(self, z=0.):
      # Plin
      f = lambda k: self.U.pLin(k,z)
      PLin = np.array(map(f, self.K))
#      # P1h
#      f = lambda k: self.p1hInt(k, z)
#      P1h = np.array(map(f, self.K))
#      print P1h
#      # P2h
#      f = lambda k: self.p2hInt(k, z)
#      P2h = np.array(map(f, self.K))
#      # noise bias
#      f = lambda k: self.pNoise(k, z)
#      PNoise = np.array(map(f, self.K))
#      # P1h+P2h
#      P = P1h + P2h + PNoise

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, self.pTotInt(self.K, z), 'k', lw=4, label=r'$P_\text{tot}$')
      ax.loglog(self.K, self.p2hInt(self.K, z), 'b-', lw=2, label=r'$P_\text{2h}$')
      ax.loglog(self.K, PLin, 'k--', lw=2, label=r'$P_\text{lin}$')
      ax.loglog(self.K, self.p1hInt(self.K, z), 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.K, self.pNoiseInt(self.K, z), 'g-', lw=2, label=r'$P_\text{noise}$')
      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$P(k)$')
      #
      path = self.pathFig + "pint_"+self.name+"z_"+str(z)+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def plotP(self, z=0.):
      # Plin
      f = lambda k: self.U.pLin(k,z)
      Plin = np.array(map(f, self.K))
      # P1h
      f = lambda k: self.p1h(k, z)
      P1h = np.array(map(f, self.K))
      # P2h
      f = lambda k: self.p2h(k, z)
      P2h = np.array(map(f, self.K))
      # noise bias
      f = lambda k: self.pNoise(k, z)
      PNoise = np.array(map(f, self.K))
      # P1h+P2h
      P = P1h + P2h + PNoise

      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, P, 'k', lw=4, label=r'$P_\text{tot}$')
      ax.loglog(self.K, P2h, 'b-', lw=2, label=r'$P_\text{2h}$')
      ax.loglog(self.K, Plin, 'k--', lw=2, label=r'$P_\text{lin}$')
      ax.loglog(self.K, P1h, 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.K, PNoise, 'g-', lw=2, label=r'$P_\text{noise}$')
      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$P(k)$')
      #
      path = self.pathFig + "p_"+self.name+"z_"+str(z)+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()


   def testPInt(self):
      # read the precomputed power spectra
      path = "./output/pn_3d/p3d_"+self.name
      Z = np.genfromtxt(path+"_z.txt")
      nZ = len(Z)
      K = np.genfromtxt(path+"_k.txt")
      P1h = np.genfromtxt(path+"_1h.txt")
      P2h = np.genfromtxt(path+"_2h.txt")
      Pshot = np.genfromtxt(path+"_shot.txt")
      P = P1h + P2h + Pshot
   
      # show precomputed values
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      for iZ in range(len(Z)):
         ax.plot(K, P[:, iZ], lw=0.5, c=plt.cm.jet(float(iZ)/nZ))
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'precomputed $P(k)$')
      
      plt.show()
      
      # show interpolated values
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      Knew = np.logspace(np.log10(1.e-3), np.log10(1.e1), 201, 10.)
      for iZ in range(len(Z)):
         z = Z[iZ]
         f = lambda k: self.fPinterp(k, z)
         Pinterp = np.array(map(f, Knew))
         ax.plot(Knew, Pinterp, lw=0.5, c=plt.cm.jet(float(iZ)/nZ))
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'interpolated $P(k)$')
      
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
      ax.set_ylabel(r'effective bias $P^\text{2h}(k=0) / P^\text{lin}(k=0)$')
      #
      path = self.pathFig + "beff_"+self.name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      #plt.show()



##################################################################################
##################################################################################


class P3dCross(P3dAuto):
   
   def __init__(self, U, MassFunc, Prof1, Prof2, pNoise=lambda k,z: 0., name="", Vs=1., save=False):
      # copy classes
      self.U = U
      self.MassFunc = MassFunc
      self.Prof1 = Prof1
      self.Prof2 = Prof2
      self.pNoise = pNoise
      self.Vs = Vs   # in (Mpc/h)^3
      self.name = str(self.Prof1) + str(self.Prof2) + name
   
      # mass bounds for integrals
      self.mMin = np.max([self.Prof1.mMin, self.Prof2.mMin, self.MassFunc.mMin])
      self.mMax = np.min([self.Prof1.mMax, self.Prof2.mMax, self.MassFunc.mMax])
      self.nM = 101  # for mass integrals

#      # k moduli, to precompute
#      self.K = np.logspace(np.log10(1.e-3), np.log10(1.e2), 101, 10.)
#      self.nK = len(self.K)
#      self.kMin = np.min(self.K)
#      self.kMax = np.max(self.K)
      
      # values of k to evaluate
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      self.Ke = np.genfromtxt("./input/K.txt") # edges of the bins for k
      self.dK = np.genfromtxt("./input/dK.txt") # widths of the bins for k

      # redshifts to evaluate
      self.Z = np.linspace(0., 10., 101)
   
      # create folders if needed
      self.pathOut = "./output/p3d/p3dcross_"+self.name+"/"
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      self.pathFig = "./figures/p3d/p3dcross_"+self.name+"/"
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)

      # power spectrum
      if save:
         self.saveP()
      self.loadP()


   ##################################################################################
   ##################################################################################
   # Precompute P3d in batch

   def computeP1hBare(self, z):
      f = lambda k,m,z: self.Prof1.u(k, m, z) * self.Prof2.u(k, m, z)
      return self.massIntegrals(z, f)

   def computeP1hCounterTerm(self, z):
      result = self.mMin 
      result *= self.MassFunc.nMinForInterp(z, self.mMin, self.mMax)
      f = lambda k: self.Prof1.u(k, self.mMin, z) * self.Prof2.u(k, self.mMin, z)
      result *= np.array(map(f, self.K))
      return result
   

   def computeBEffBare1(self, z):
      f = lambda k,m,z: self.MassFunc.b1(m, z) * self.Prof1.u(k, m, z)
      return self.massIntegrals(z, f)

   def computeBEffBare2(self, z):
      f = lambda k,m,z: self.MassFunc.b1(m, z) * self.Prof2.u(k, m, z)
      return self.massIntegrals(z, f)

   def computeBEffCounterTerm1(self, z):
      result = self.mMin
      result *= self.MassFunc.nMinForInterp(z, self.mMin, self.mMax) 
      result *= self.MassFunc.b1MinForInterp(z, self.mMin, self.mMax)
      f = lambda k: self.Prof1.u(k, self.mMin, z)
      result *= np.array(map(f, self.K))
      return result
      
   def computeBEffCounterTerm2(self, z):
      result = self.mMin
      result *= self.MassFunc.nMinForInterp(z, self.mMin, self.mMax) 
      result *= self.MassFunc.b1MinForInterp(z, self.mMin, self.mMax)
      f = lambda k: self.Prof2.u(k, self.mMin, z)
      result *= np.array(map(f, self.K))
      return result
      
   def computeBEffCorrected1(self, z):
      result = self.computeBEffBare1(z)
      result *= self.MassFunc.correctionFactorBias(z)
      return result

   def computeBEffCorrected2(self, z):
      result = self.computeBEffBare2(z)
      result *= self.MassFunc.correctionFactorBias(z)
      return result


   def saveP(self):
      print "Precomputing p3d "+self.name
      #print "FAST VERSION"
      tStart = time()
      nZ = len(self.Z)
      nK = len(self.K)
      
      # power spectra
      P1h = np.zeros((nK, nZ))
      P2h = np.zeros((nK, nZ))
      PNoise = np.zeros((nK, nZ))
      #dP = np.zeros((nK, nZ))
      P1hBare = np.zeros((nK, nZ))
      P1hCounterTerm = np.zeros((nK, nZ))
      BEffBare1 = np.zeros((nK, nZ))
      BEffCounterTerm1 = np.zeros((nK, nZ))
      BEffCorrected1 = np.zeros((nK, nZ))
      BEffBare2 = np.zeros((nK, nZ))
      BEffCounterTerm2 = np.zeros((nK, nZ))
      BEffCorrected2 = np.zeros((nK, nZ))


      for iZ in range(nZ):
         z = self.Z[iZ]
         P1hBare[:,iZ] = self.computeP1hBare(z)
         P1hCounterTerm[:,iZ] = self.computeP1hCounterTerm(z)
         #
         BEffBare1[:,iZ] = self.computeBEffBare1(z)
         BEffCounterTerm1[:,iZ] = self.computeBEffCounterTerm1(z)
         BEffCorrected1[:,iZ] = self.computeBEffCorrected1(z)
         BEffBare2[:,iZ] = self.computeBEffBare2(z)
         BEffCounterTerm2[:,iZ] = self.computeBEffCounterTerm2(z)
         BEffCorrected2[:,iZ] = self.computeBEffCorrected2(z)
         print "done z "+str(iZ)+" of "+str(nZ)

      P1h = P1hBare + P1hCounterTerm
      #
      BEff1 = BEffBare1 + BEffCounterTerm1
      BEff2 = BEffBare2 + BEffCounterTerm2
      pLinVect = np.vectorize(self.U.pLin)
      f = lambda k: pLinVect(k, self.Z)
      PLin = np.array(map(f, self.K))
      #
      P2hBare = BEffBare1 * BEffBare2 * PLin
      P2hCounterTerm = BEffCounterTerm1 * BEffCounterTerm2 * PLin
      P2h = BEff1 * BEff2 * PLin
      P2hCorrected = BEffCorrected1 * BEffCorrected2 * PLin

      pNoiseVect = np.vectorize(self.pNoise)
      f = lambda k: pNoiseVect(k, self.Z)
      PNoise = np.array(map(f, self.K))

      # save power spectra
      path = self.pathOut+"p3d_"+self.name
      np.savetxt(path+"_z.txt", self.Z)
      np.savetxt(path+"_k.txt", self.K)
      np.savetxt(path+"_1h.txt", P1h)
      np.savetxt(path+"_2h.txt", P2h)
      np.savetxt(path+"_noise.txt", PNoise)
      #np.savetxt(path+"_dPdd.txt", dP)
      np.savetxt(path+"_1hbare.txt", P1hBare)
      np.savetxt(path+"_1hcounterterm.txt", P1hCounterTerm)
      np.savetxt(path+"_2hbare.txt", P2hBare)
      np.savetxt(path+"_2hcounterterm.txt", P2hCounterTerm)
      np.savetxt(path+"_2hcorrected.txt", P2hCorrected)

      tStop = time()
      print "took ", (tStop-tStart)/60., "min"
      return




   ##################################################################################
   ##################################################################################
   # 1-halo term


   def p1hBare(self, k, z, mMin=0., mMax=np.inf):
      '''P1h, without the mass counterterm
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof1.u(k, m, z)
         result *= self.Prof2.u(k, m, z)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      if not np.isfinite(result):
         result = 0.
      return result


   def p1hCounterTerm(self, k, z, mMin=0., mMax=np.inf):
      '''Mass counter term for P1h,
      only evaluated at the fiducial mass cut
      '''
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      # counter term
      result = self.MassFunc.nMinForInterp(z, mMin, mMax) 
      result *= mMin
      result *= self.Prof1.u(k, mMin, z)
      result *= self.Prof2.u(k, mMin, z)
      return result
      

   ##################################################################################
   # 2-halo term


   def bEffBare1(self, k, z, mMin=0., mMax=np.inf):
      '''Without mass counter term
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof1.u(k, m, z)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result

   def bEffBare2(self, k, z, mMin=0., mMax=np.inf):
      '''Without mass counter term
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.MassFunc.b1(m, z)
         result *= self.Prof2.u(k, m, z)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result


   def bEffCounterTerm1(self, k, z, mMin=0., mMax=np.inf):
      '''Mass counter term,
      such that the matter bias is indeed unity.
      '''
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      # counter term
      result = self.MassFunc.nMinForInterp(z, mMin, mMax) 
      result *= self.MassFunc.b1MinForInterp(z, mMin, mMax)
      result *= mMin
      result *= self.Prof1.u(k, mMin, z)
      return result

   def bEffCounterTerm2(self, k, z, mMin=0., mMax=np.inf):
      '''Mass counter term,
      such that the matter bias is indeed unity.
      '''
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      # counter term
      result = self.MassFunc.nMinForInterp(z, mMin, mMax) 
      result *= self.MassFunc.b1MinForInterp(z, mMin, mMax)
      result *= mMin
      result *= self.Prof2.u(k, mMin, z)
      return result


   def bEff1(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of the bare and counter term
      '''
      result = self.bEffBare1(k, z, mMin, mMax)
      result += self.bEffCounterTerm1(k, z, mMin, mMax)
      return result

   def bEff2(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of the bare and counter term
      '''
      result = self.bEffBare2(k, z, mMin, mMax)
      result += self.bEffCounterTerm2(k, z, mMin, mMax)
      return result


   def p2hBare(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEffBare1(k, z, mMin=mMin, mMax=mMax)
      result *= self.bEffBare2(k, z, mMin=mMin, mMax=mMax)
      result *= self.U.pLin(k, z)
      return result


   def p2hCounterTerm(self, k, z, mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEffCounterTerm1(k, z, mMin=mMin, mMax=mMax)
      result *= self.bEffCounterTerm2(k, z, mMin=mMin, mMax=mMax)
      result *= self.U.pLin(k, z)
      return result


   def p2h(self, k, z, mu=0., mMin=0., mMax=np.inf):
      '''Sum of bare and counter term
      '''
      result = self.bEff1(k, z, mMin=mMin, mMax=mMax) 
      result *= self.bEff2(k, z, mMin=mMin, mMax=mMax) 
      result *= self.U.pLin(k, z)
      return result


   ##################################################################################



