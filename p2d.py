from headers import *

##################################################################################

class P2dAuto(object):

   def __init__(self, U, P3dAuto, Weight, name="", pNoise=lambda l: 0., save=False, nProc=1):
      # copy classes
      self.U = U
      self.Pn = P3dAuto
      self.Weight = Weight
      self.name = str(Weight) + name
      self.pNoise = pNoise
      self.nProc = nProc
      
      # bounds for z integrals
      self.aMin = self.Weight.aMin
      self.aMax = self.Weight.aMax
      
      # values of ell to evaluate
      self.L = np.genfromtxt("./input/Lc.txt") # center of the bins for l

      # create directory if needed
      self.pathOut = "./output/p2d/p2dauto_"+self.name+"/"
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      self.pathFig = "./figures/p2d/p2dauto_"+self.name+"/"
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      # power spectrum
      if save:
         self.saveP()
      self.loadP()
      
      
   ##################################################################################
      
#   def saveP(self):
#      print "precomputing p2d "+self.name
#      data = np.zeros((len(self.L), 4))
#      data[:,0] = self.L.copy()
#      pool = Pool(ncpus=self.nProc)
#      data[:,1] = np.array(map(self.fP_1h, self.L))
#      data[:,2] = np.array(map(self.fP_2h, self.L))
#      data[:,3] = np.array(map(self.fPnoise, self.L))
#      np.savetxt(self.pathOut+"p2d_"+self.name+".txt", data)
#
#   def loadP(self):
#      data = np.genfromtxt(self.pathOut+"p2d_"+self.name+".txt")
#      self.L = data[:,0]
#      self.P1h = data[:,1]
#      self.P2h = data[:,2]
#      self.Pnoise = data[:,3]
#      self.P = self.P1h + self.P2h
#      self.Ptot = self.P1h + self.P2h + self.Pnoise
#      #
#      # interpolate power spectra
#      forP1h = UnivariateSpline(self.L,self.P1h,k=1,s=0)
#      self.fP1hinterp = lambda l: forP1h(l)*(l>=min(self.L))*(l<=max(self.L))
#      forP2h = UnivariateSpline(self.L,self.P2h,k=1,s=0)
#      self.fP2hinterp = lambda l: forP2h(l)*(l>=min(self.L))*(l<=max(self.L))
#      forP = UnivariateSpline(self.L,self.P,k=1,s=0)
#      self.fPinterp = lambda l: forP(l)*(l>=min(self.L))*(l<=max(self.L))
#      forPtot = UnivariateSpline(self.L,self.Ptot,k=1,s=0)
#      self.fPtotinterp = lambda l: forPtot(l)*(l>=min(self.L))*(l<=max(self.L))
   

   ##################################################################################

   def computeP(self, fp3d):
      '''Compute P2d for all self.L at once,
      given the 3d power spectrum fp3d.
      '''
      z = self.Pn.Z.copy()
      if z[0]==0:
         z = z[1:]
      a = 1. / (1. + z)
      chi = self.U.bg.comoving_distance(z)

      integrand = 3.e5/( self.U.hubble(z) * a**2 )
      integrand *= self.Weight.f(a)**2 / chi**2

      fp3dVect = np.vectorize(fp3d)
      f = lambda l: fp3dVect((l + 0.5)/chi, z)
      integrand = integrand[None,:] * np.array(map(f, self.L))
      integrand *= -1.  # because a is in decreasing order
      
      result = np.trapz(integrand, a, axis=-1)
      return result

   
   def saveP(self):
      print "precomputing p2d "+self.name
      data = np.zeros((len(self.L), 8))
      data[:,0] = self.L.copy()
      #pool = Pool(ncpus=self.nProc)
      data[:,1] = self.computeP(self.Pn.p1hInt)
      data[:,2] = self.computeP(self.Pn.p2hInt)
      data[:,3] = np.array(map(self.pNoise, self.L))

      print "precomputing bare vs counter terms"
      data[:,4] = self.computeP(self.Pn.p1hBareInt)
      data[:,5] = self.computeP(self.Pn.p1hCounterTermInt)
      data[:,6] = self.computeP(self.Pn.p2hBareInt)
      data[:,7] = self.computeP(self.Pn.p2hCounterTermInt)

      np.savetxt(self.pathOut+"p2d_"+self.name+".txt", data)


   def loadP(self):
      data = np.genfromtxt(self.pathOut+"p2d_"+self.name+".txt")
      self.L = data[:,0]
      self.P1h = data[:,1]
      self.P2h = data[:,2]
      self.Pnoise = data[:,3]
      self.P = self.P1h + self.P2h
      self.Ptot = self.P1h + self.P2h + self.Pnoise
      #
      # interpolate power spectra
      forP1h = UnivariateSpline(self.L,self.P1h,k=1,s=0)
      self.fP1hinterp = lambda l: forP1h(l)*(l>=min(self.L))*(l<=max(self.L))
      forP2h = UnivariateSpline(self.L,self.P2h,k=1,s=0)
      self.fP2hinterp = lambda l: forP2h(l)*(l>=min(self.L))*(l<=max(self.L))
      forP = UnivariateSpline(self.L,self.P,k=1,s=0)
      self.fPinterp = lambda l: forP(l)*(l>=min(self.L))*(l<=max(self.L))
      forPtot = UnivariateSpline(self.L,self.Ptot,k=1,s=0)
      self.fPtotinterp = lambda l: forPtot(l)*(l>=min(self.L))*(l<=max(self.L))



   ##################################################################################
   # power spectrum

   def integrandP(self, a, fP, l):
      '''dP2d/da, to be integrated wrt a
      '''
      z = 1./a-1.
      chi = self.U.bg.comoving_distance(z)
      #
      result = 3.e5/( self.U.hubble(z) * a**2 )
      result *= self.Weight.f(a)**2
      result /= chi**2
      result *= fP(l/chi, z)
      return result
   
   def fP_1h(self, l):
      f = lambda a: self.integrandP(a, self.Pn.p1hInt, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result

   def fP_2h(self, l):
      f = lambda a: self.integrandP(a, self.Pn.p2hInt, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result

   def fP(self, l):
      result = self.fP_1h(l) + self.fP_2h(l)
      return result

   ##################################################################################
   

   def fdPdz_1h(self, l, z):
      a = 1./(1.+z)
      result = self.integrandP(a, self.Pn.fP1hinterp, l)
      result *= a**2
      return result

   def fdPdz_2h(self, l, z):
      a = 1./(1.+z)
      result = self.integrandP(a, self.Pn.fP2hinterp, l)
      result *= a**2
      return result

   def fdPdz(self, l, z):
      result = self.fdPdz_1h(l, z) + self.fdPdz_2h(l, z)
      return result

   def fdPnoisedz(self, l, z):
      a = 1./(1.+z)
      result = self.Weight.fdPshotNoise_da(a, l)
      result *= a**2
      return result



   ##################################################################################
   # trispectrum

   def integrandT(self, a, fP, l):
      z = 1./a-1.
      chi = self.U.ComovDist(a, self.U.a_obs)
      #
      result = 3.e5/( self.U.Hubble(a) * a**2 )
      result *= self.Weight.f(a)**4
      result /= chi**6
      result *= fP(l/chi, z)
      return result

   def fT_1h(self, l):
      f = lambda a: self.integrandT(a, self.Pn.fT1hinterp, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result

   def fT_2h(self, l):
      f = lambda a: self.integrandT(a, self.Pn.fT2hinterp, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result

   def fT_4h(self, l):
      f = lambda a: self.integrandT(a, self.Pn.fT4hinterp, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result

   def fT(self, l):
      result = self.fT_1h(l)
      return result


   def fT_ssv(self, l, L=1.e2):
      """This is the difference between the almost-squeezed
      and the exactly squeezed trispectra:
      T(l, -l+L, l, -l-L) = T(l,-l,l,-l) + T_ssv,
      where L << l.
      """
      g = lambda k,z: self.Pn.fT_ssv(k, k, k*L/l, z)
      f = lambda a: self.integrandT(a, g, l)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l
      return result



   ##################################################################################
   # trispectrum non-diagonal

   def integrandTNonDiag(self, a, fP, l1, l2):
      z = 1./a-1.
      chi = self.U.ComovDist(a, self.U.a_obs)
      #
      result = 3.e5/( self.U.Hubble(a) * a**2 )
      result *= self.Weight.f(a)**4
      result /= chi**6
      result *= fP(l1/chi, l2/chi, z)
      return result

   def fTnondiag(self, l1, l2):
      f = lambda a: self.integrandTNonDiag(a, self.Pn.fTnondiag, l1, l2)
      result = integrate.quad(f, self.aMin, self.aMax, epsabs=0, epsrel=1.e-2)[0]
      print "done ell=",l1, l2
      return result


   ##################################################################################

   def plotP(self):
      
      # P
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.L, self.P1h+self.P2h+self.Pnoise, 'k', lw=4, label=r'$P_\text{total}$')
      ax.loglog(self.L, self.P2h, 'b-', lw=2, label=r'$P_\text{2h}$')
      ax.loglog(self.L, self.P1h, 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.L, self.Pnoise, 'g-', lw=2, label=r'$P_\text{noise}$')
      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$P(\ell)$')
      path = "./figures/pn2d/p_"+self.name+".pdf"
      #fig.savefig(path)
      
      # l*(l+1)*P
      fig = plt.figure(1)
      ax = plt.subplot(111)
      #
      factor = self.L*(self.L+1.)/(2.*np.pi)
      ax.loglog(self.L, factor*(self.P1h+self.P2h+self.Pnoise), 'k', lw=4, label=r'$P_\text{total}$')
      ax.loglog(self.L, factor*self.P2h, 'b-', lw=2, label=r'$P_\text{2h}$')
      ax.loglog(self.L, factor*self.P1h, 'r-', lw=2, label=r'$P_\text{1h}$')
      ax.loglog(self.L, factor*self.Pnoise, 'g-', lw=2, label=r'$P_\text{noise}$')
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$\ell(\ell+1)P(\ell) / (2\pi)$')
      path = "./figures/pn2d/l2p_"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      
      plt.show()
   
   
   def plotdP(self):
      
      # P
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      factor = self.L*(self.L+1.)/(2.*np.pi)
      ax.loglog(L, factor*self.dP1h, '--', label=r'1h')
      ax.loglog(L, factor*self.dP2h, '--', label=r'2h')
      ax.loglog(L, factor*(self.dP1h+self.dP2h), 'k', label=r'tot')
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'l')
      ax.set_ylabel(r'$\ell(\ell+1)\frac{dP}{d\delta}(\ell)$')
      ax.set_title(r'Power spectrum response')
      path = "./figures/pn2d/dp_"+self.name+".pdf"
      #      fig.savefig(path, bbox_inches='tight')
      
      # P
      fig = plt.figure(1)
      ax = plt.subplot(111)
      #
      ax.semilogx(self.L, self.dP / self.P, 'b-')
      ax.grid()
      ax.set_xlabel(r'l')
      ax.set_ylabel(r'$\frac{dlnP}{d\delta}(\ell)$')
      ax.set_title(r'Power spectrum response')
      
      plt.show()



   def plotdPdz(self, l=1.e3):
      A = np.linspace(self.aMin, self.aMax, 201)
      Z = 1./A-1.
      print Z
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), A))
      H = np.array(map(lambda a: self.U.Hubble(a), A))
      W = np.array(map(self.Weight.f, A))
      dChidA = 3.e5 / (H*A**2)
      dChidZ = 3.e5 / H
      
      # redshift contributions for P1h and P2h
      f = lambda a: self.integrandP(a, self.Pn.fP_1h, l)
      dP1h_da = np.array(map(f, A))
      f = lambda a: self.integrandP(a, self.Pn.fP_2h, l)
      dP2h_da = np.array(map(f, A))
      #
      dP1h_dz = dP1h_da * A**2
      dP2h_dz = dP2h_da * A**2
      
      # redshift contributions for Pshot
      if hasattr(self.Weight, 'fdPshotNoise_da'):
         f = lambda a: self.Weight.fdPshotNoise_da(a, l)
         dPshot_da = np.array(map(f, A))
         dPshot_dz = dPshot_da * A**2
      
      '''
      def f(a):
         z = 1./a-1.
         chi = self.U.ComovDist(a, 1.)
         return self.Pn.fP_1h(l/chi, z)
      P3d_1h = np.array(map(f, A))
      
      def f(a):
         z = 1./a-1.
         chi = self.U.ComovDist(a, 1.)
         return self.Pn.fP_2h(l/chi, z)
      P3d_2h = np.array(map(f, A))
      
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(A, A* dChidA * W**2/Chi**2 / np.max(A* dChidA * W**2/Chi**2), 'k', lw=2, label=r'kernel')
      #
      ax.plot(A, P3d_1h / np.max(P3d_1h), 'b--', lw=2, label=r'$P_\text{3d}^\text{1h}$')
      ax.plot(A, A*dP1h_da/np.max(A*dP1h_da), 'b', lw=2, label=r'integrand 1h')
      #
      ax.plot(A, P3d_2h / np.max(P3d_2h), 'g--', lw=2, label=r'$P_\text{3d}^\text{2h}$')
      ax.plot(A, A*dP2h_da/np.max(A*dP2h_da), 'g', lw=2, label=r'integrand 2h')
      #
      ax.legend(loc=3)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'scale factor $a$')
      ax.set_ylabel(r'$d C_{\ell='+str(int(l))+'} / d\ln a$')
      #
      path = "./figures/pn2d/dp2d_dlna_"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      '''

      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # factors in the integrands
      #ax.plot(Z, dChidZ * W**2/Chi**2 / np.max(dChidZ * W**2/Chi**2), 'k', lw=2, label=r'kernel')
      #ax.plot(Z, P3d_1h / np.max(P3d_1h), 'b--', lw=2, label=r'$P_\text{3d}^\text{1h}$')
      #ax.plot(Z, P3d_2h / np.max(P3d_2h), 'g--', lw=2, label=r'$P_\text{3d}^\text{2h}$')
      #
      # normalized integrands
#      ax.plot(Z, dP1h_dz / np.max(dP1h_dz), 'b', lw=2, label=r'1h')
#      ax.plot(Z, dP2h_dz / np.max(dP2h_dz), 'g', lw=2, label=r'2h')
#      if hasattr(self.Weight, 'fdPshotNoise_da'):
#         ax.plot(Z, dPshot_dz / np.max(dPshot_dz), 'r', lw=2, label=r'shot')
#         ax.plot(Z, (dP1h_dz+dP2h_dz+dPshot_dz) / np.max(dP1h_dz+dP2h_dz+dPshot_dz), 'k', lw=2, label=r'1h+2h+shot')
#      else:
#         ax.plot(Z, (dP1h_dz+dP2h_dz) / np.max(dP1h_dz+dP2h_dz), 'k', lw=2, label=r'1h+2h')
      #
      # non-normalized ingredients
      ax.plot(Z, dP1h_dz, 'b', lw=2, label=r'1h')
      ax.plot(Z, dP2h_dz, 'g', lw=2, label=r'2h')
      if hasattr(self.Weight, 'fdPshotNoise_da'):
         ax.plot(Z, dPshot_dz, 'r', lw=2, label=r'shot')
         ax.plot(Z, dP1h_dz+dP2h_dz+dPshot_dz, 'k', lw=2, label=r'1h+2h+shot')
      else:
         ax.plot(Z, dP1h_dz+dP2h_dz, 'k', lw=2, label=r'1h+2h')
      #
      ax.legend(loc=4)
      #ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'redshift $z$')
      ax.set_ylabel(r'$d C_{\ell='+str(int(l))+'} / dz$')
      #
      path = "./figures/pn2d/dp2d_dz"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
   
      '''
      fig=plt.figure(2)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dP1h_dz / np.max(dP1h_dz), 'b', lw=2, label=r'1h')
      #
      ax.plot(Z, dP2h_dz / np.max(dP2h_dz), 'g', lw=2, label=r'2h')
      #
      ax.legend(loc=4)
      #ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'redshift $z$')
      ax.set_ylabel(r'$d C_{\ell='+str(int(l))+'} / dz$ [arbitrary unit]')
      #
      path = "./figures/pn2d/dp2d_dz_summary"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      '''

      plt.show()


   def plotdPdz_color(self):
   
      # redshifts to evaluate
      nZ = 51
      zMin = 1./self.aMax-1.
      zMax = 5.   #1./self.Weight.aMin-1.
      dZ = (zMax-zMin)/nZ
      Z = np.linspace(zMin, zMax, nZ)
      zEdges = np.linspace(zMin-0.5*dZ, zMax+0.5*dZ, nZ+1)

      A = 1./(1.+Z)
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), A))
      H = np.array(map(lambda a: self.U.Hubble(a), A))
      W = np.array(map(self.Weight.f, A))
      dChidA = 3.e5 / (H*A**2)
      dChidZ = 3.e5 / H
      
      # multipoles to evaluate
      nL = 51  #51
      lnlMin = np.log10(10.)
      lnlMax = np.log10(1.e4)
      dlnl = (lnlMax-lnlMin)/nL
      lnL = np.linspace(lnlMin, lnlMax, nL)
      lnlEdges = np.linspace(lnlMin-0.5*dlnl, lnlMax+0.5*dlnl, nL+1)
      L = 10.**lnL
      lEdges = 10.**lnlEdges

      '''
      # 1h
      dP1hdz = np.zeros((nZ, nL))
      for iL in range(nL):
         l = L[iL]
         f = lambda a: self.integrandP(a, self.Pn.fP_1h, l)
         dP1hdz[:,iL] = np.array(map(f, A))
         dP1hdz[:,iL] *= A**2
         #dP1hdz[:,iL] /= np.trapz(Z, dP1hdz[:,iL])
      print "done 1h"
      
      # 2h
      dP2hdz = np.zeros((nZ, nL))
      for iL in range(nL):
         l = L[iL]
         f = lambda a: self.integrandP(a, self.Pn.fP_2h, l)
         dP2hdz[:,iL] = np.array(map(f, A))
         dP2hdz[:,iL] *= A**2
         #dP2hdz[:,iL] /= np.trapz(Z, dP2hdz[:,iL])
      print "done 2h"

      # shot noise
      if hasattr(self.Weight, 'fdPshotNoise_da'):
         dPshotdz = np.zeros((nZ, nL))
         for iL in range(nL):
            l = L[iL]
            f = lambda a: self.integrandP(a, self.Weight.fdPshotNoise_da, l)
            dPshotdz[:,iL] = np.array(map(f, A))
            dPshotdz[:,iL] *= A**2
            #dPshotdz[:,iL] /= np.trapz(Z, dPshotdz[:,iL])
         print "done shot noise"
      '''
   
      # total
      dPdz = np.zeros((nZ, nL))
      for iL in range(nL):
         l = L[iL]
         f = lambda a: self.integrandP(a, self.Pn.fP, l)
         dPdz[:,iL] = np.array(map(f, A))
         dPdz[:,iL] *= A**2
#         # normalize so int dz dP/dz = 1 for all ell
#         dPdz[:,iL] /= np.trapz(Z, dPdz[:,iL])
      dPdz = np.abs(dPdz)
      print "done total"
      
      # show the 2d color plot
      zz,ll = np.meshgrid(zEdges, lEdges, indexing='ij')

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      cp=ax.pcolormesh(zz, ll, np.log(dPdz), linewidth=0, rasterized=True, cmap=plt.cm.YlOrRd_r)
      #
      cp.set_clim(0., 13.)
      cb=fig.colorbar(cp)
      cb.ax.set_title(r'$\text{ln}\left(\frac{dC^0_\ell}{dz}\right)$')
      ax.set_yscale('log')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'$\ell$')
      #ax.set_title(r'$\frac{dC^0_\ell}{dz}$, normalized to $\int dz\; \frac{dC^0_\ell}{dz} = 1$', fontsize=18)
      
      plt.show()


   def plotPlanckCIB(self):
      
      # load Planck data points (Planck 13 XXX)
      nu = self.Weight.nu
      name = str(int(nu/1.e9))
      p13 = Planck13CIBData()
      L = p13.PlanckPCIB['ell']
      P = p13.PlanckPCIB[name]
      sP = p13.PlanckPCIB[name+'_error']
      Pshot = p13.PlanckPCIBShot[name]
      
      
      # sensitivity in Jy/rad
      # beam in arcmin
      def fdetectorNoise(l, sensitivity, beam):
         beam *= np.pi/180./60.  # convert arcmin to rad
         sigma_beam = beam / np.sqrt(8.*np.log(2.))   # convert fwhm to sigma
         return sensitivity**2 * np.exp(l**2 * sigma_beam**2)
      
      # compare to some noise levels
      # !!!!! these are only relevant for \sim 545GHz, for Planck and CCAT
      f = lambda l: fdetectorNoise(l, sensitivity=13.5, beam=4.8)
      noisePlanck = np.array(map(f, self.L))
      f = lambda l: fdetectorNoise(l, sensitivity=1.2, beam=0.5)
      noiseCCAT = np.array(map(f, self.L))
      
      
      # P
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      # halo model
#      ax.plot(self.L, self.P2h, c=plt.cm.autumn(0.4), lw=1.5, label=r'2-halo')
#      ax.plot(self.L, self.P1h, c=plt.cm.autumn(0.7), lw=1.5, label=r'1-halo')
#      ax.plot(self.L, self.Pnoise, c=plt.cm.autumn(0.9), lw=1.5, label=r'1-galaxy')
#      ax.plot(self.L, (self.P1h+self.P2h+self.Pnoise), 'r', lw=3, label=r'total')
      #
      ax.plot(self.L, self.P2h, c='r', lw=1.5, label=r'2-halo')
      ax.plot(self.L, self.P1h, c='g', lw=1.5, label=r'1-halo')
      ax.plot(self.L, self.Pnoise, c='y', lw=1.5, label=r'1-galaxy')
      ax.plot(self.L, (self.P1h+self.P2h+self.Pnoise), 'b', lw=3, label=r'total')
      
      #
      # Planck data points
      ax.errorbar(L, P, yerr=sP, fmt='.', c='k', label='Planck13 XXX')
      #
      # noise levels
      ax.plot(self.L, noisePlanck, c='gray', ls='--', lw=1, label=r'Planck noise')
      ax.plot(self.L, noiseCCAT, c='grey', ls='-.', lw=1, label=r'CCAT noise')
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((10., 5.e4))
      ax.set_ylim((1.e-1, 1.e6))
      ax.legend(loc=3, numpoints=1, fontsize=14, framealpha=1)
      ax.set_xlabel(r'$\ell$')
      ax.set_ylabel(r'$C_\ell^{545\text{GHz}}$ [Jy$^2$/sr]')
      #
      #path="./figures/pn2d/p_"+str(self.Weight)+".pdf"
      #path="./figures/cib_penin12/planck13model1_vs_planck"+name+".pdf"
      path="./figures/cib_penin12/penin1214_"+name+"GHz.pdf"
      #path="./figures/cib_planck13/planck13_"+name+"GHz.pdf"
      fig.savefig(path, bbox_inches='tight')

      plt.show()


   ##################################################################################

   def plotT(self, func=None):
      factor = 1.#self.L*(self.L+1.)/(2.*np.pi)
      
      # T
      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.plot(self.L, factor*self.Ptot**2, 'k--', lw=3, label=r'$C^2$')
      #
      ax.plot(self.L, factor*self.Ttot, 'k', lw=3, label=r'$T^\text{total}$')
      ax.plot(self.L, factor*self.T1h, 'r', lw=1, label=r'$T^{1h}$')
      ax.plot(self.L, factor*self.T2h, 'orange', lw=1, label=r'$T^{2h}$')
      ax.plot(self.L, factor*self.T4h, 'gold', lw=1, label=r'$T^{4h}$')
      ax.plot(self.L, factor*self.Tssv, 'm', lw=1, label=r'$T^\text{SSV}$')
      ax.plot(self.L, factor*self.Tnoise, 'fuchsia', lw=1, label=r'$T^\text{shot}$')
      #
      if func is not None:
         F = np.array(map(func, self.L))
         ax.plot(self.L, factor*F**2, 'b', lw=2)
      #
      ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'\ell')
      #ax.set_xlim((50., 5.e4))
      #ax.set_ylabel(r'$T(\ell)$')
      ax.set_ylabel(r'$\mathcal{T}^0_{\ell, L-\ell, \ell, -L-\ell}$')
      path = "./figures/pn2d/t_"+self.name+"_test.pdf"
      fig.savefig(path, bbox_inches='tight')
      
      plt.show()



   def plotIntegrandT(self, l=5.e2):
      A = np.linspace(self.aMin, self.aMax, 101)
      Z = 1./A-1.
      Chi = np.array(map(lambda a: self.U.ComovDist(a, 1.), A))
      H = np.array(map(lambda a: self.U.Hubble(a), A))
      W = np.array(map(self.Weight.f, A))
      dChidA = 3.e5 / (H*A**2)
      dChidZ = 3.e5 / H
      
      def f(a):
         z = 1./a-1.
         chi = self.U.ComovDist(a, 1.)
         return self.Pn.fT1hinterp(l/chi, z)
      T3d_1h = np.array(map(f, A))
      
      #      print T3d_1h
      
      # integrand
      f = lambda a: self.integrand(a, self.Pn.fT1hinterp, l)
      dT1h_da = np.array(map(f, A))
      
      dT1h_dz = dT1h_da * A**2
      
      #      print dT1h_dz
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(A, A* dChidA * W**4/Chi**6 / np.max(A* dChidA * W**4/Chi**6), 'k', lw=2, label=r'kernel')
      ax.plot(A, T3d_1h / np.max(T3d_1h), 'b--', lw=2, label=r'$T_\text{1h}^\text{3d}$')
      ax.plot(A, A*dT1h_da/np.max(A*dT1h_da), 'b', lw=2, label=r'integrand for T1h')
      #
      ax.legend(loc=2)
      ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'scale factor $a$')
      ax.set_ylabel(r'$d T_{\ell='+str(int(l))+'} / d\ln a$')
      #
      path = "./figures/pn2d/dt2d_dlna_"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.plot(Z, dChidZ * W**4/Chi**6 / np.max(dChidZ * W**4/Chi**6), 'k', lw=2, label=r'kernel')
      ax.plot(Z, T3d_1h / np.max(T3d_1h), 'b--', lw=2, label=r'$T_\text{1h}^\text{3d}$')
      ax.plot(Z, dT1h_dz/np.max(dT1h_dz), 'b', lw=2, label=r'integrand for T1h')
      #
      ax.legend(loc=2)
      #ax.set_xscale('log', nonposx='clip')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(r'redshift $z$')
      ax.set_ylabel(r'$d T_{\ell='+str(int(l))+'} / dz$')
      #
      path = "./figures/pn2d/dt2d_dz_"+self.name+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      
      
      plt.show()






##################################################################################
##################################################################################

class P2dCross(P2dAuto):
   
   def __init__(self, U, P3dCross, Weight1, Weight2, name="", pNoise=lambda l: 0., save=False, nProc=1):
      # copy classes
      self.U = U
      self.Pn = P3dCross
      self.Weight1 = Weight1
      self.Weight2 = Weight2
      self.name = str(self.Weight1) + str(self.Weight2) + name
      self.pNoise = pNoise
      self.nProc = nProc
      
      # bounds for z integrals
      self.aMin = max(self.Weight1.aMin, self.Weight2.aMin)
      self.aMax = min(self.Weight1.aMax, self.Weight2.aMax)
      
      # values of ell to evaluate
      self.L = np.genfromtxt("./input/Lc.txt") # center of the bins for l
      
      # create directory if needed
      self.pathOut = "./output/p2d/p2dauto_"+self.name+"/"
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      self.pathFig = "./figures/p2d/p2dauto_"+self.name+"/"
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      # power spectrum
      if save:
         self.saveP()
      self.loadP()


   ##################################################################################


   def computeP(self, fp3d):
      '''Compute P2d for all self.L at once,
      given the 3d power spectrum fp3d.
      '''
      z = self.Pn.Z.copy()
      if z[0]==0:
         z = z[1:]
      a = 1. / (1. + z)
      chi = self.U.bg.comoving_distance(z)

      integrand = 3.e5/( self.U.hubble(z) * a**2 )
      integrand *= self.Weight1.f(a)
      integrand *= self.Weight2.f(a)
      integrand /= chi**2

      fp3dVect = np.vectorize(fp3d)
      f = lambda l: fp3dVect((l + 0.5)/chi, z)
      integrand = integrand[None,:] * np.array(map(f, self.L))
      integrand *= -1.  # because a is in decreasing order
      
      result = np.trapz(integrand, a, axis=-1)
      return result





   ##################################################################################
   
   def integrandP(self, a, fP, l):
      z = 1./a-1.
      chi = self.U.ComovDist(a, self.U.a_obs)
      #
      result = 3.e5/( self.U.Hubble(a) * a**2 )
      result *= self.Weight1.f(a) * self.Weight2.f(a)
      result /= chi**2
      result *= fP(l/chi, z)
      return result

   def integrandT(self, a, fP, l):
      z = 1./a-1.
      chi = self.U.ComovDist(a, self.U.a_obs)
      #
      result = 3.e5/( self.U.Hubble(a) * a**2 )
      result *= self.Weight1.f(a)**2 * self.Weight2.f(a)**2
      result /= chi**6
      result *= fP(l/chi, z)
      return result








##################################################################################
##################################################################################

class Planck13CIBData(object):
   """Measured CIB power spectra from Planck13 XXX
   """

   # here name has to be the frequency of the CIB map
   def __init__(self):

      # Shot noise on CIB power spectra
      # table 9 from Planck13 XXX on CIB
      # unit is Jy^2 / sr
      self.PlanckPCIBShot = {}
      self.PlanckPCIBShot['857'] = 5364.
      self.PlanckPCIBShot['857_545'] = 2701.
      self.PlanckPCIBShot['857_353'] = 953.
      self.PlanckPCIBShot['857_217'] = 181.
      self.PlanckPCIBShot['545'] = 1690.
      self.PlanckPCIBShot['545_353'] = 626.
      self.PlanckPCIBShot['545_217'] = 121.
      self.PlanckPCIBShot['353'] = 262.
      self.PlanckPCIBShot['353_217'] = 54.
      self.PlanckPCIBShot['217'] = 21.

      # Measured CIB power spectra
      # table D2 from Planck13 XXX on CIB
      # unit is Jy^2/sr
      self.PlanckPCIB = {}

      # ell values
      self.PlanckPCIB['ell'] = np.array([53., 114., 187., 320., 502., 684., 890., 1158., 1505., 1956., 2649.])

      # Auto and cross spectra of the CIB: these maps are CMB-free, Galactic dust-free, corrected for SZ contamination and CIB contamination from CMB template
      # The values of the spectra at the first two ell values are only upper limits
      self.PlanckPCIB['857'] = np.array([0., 0., 2.87e5, 1.34e5, 7.20e4, 4.38e4, 3.23e4, 2.40e4, 1.83e4, 1.46e4, 1.16e4])
      self.PlanckPCIB['857_545'] = np.array([0., 0., 1.30e5, 6.36e4, 3.53e4, 2.21e4, 1.63e4, 1.22e4, 9.31e3, 7.38e3, 5.91e3])
      self.PlanckPCIB['857_353'] = np.array([0., 0., 4.30e4, 2.20e4, 1.25e4, 7.99e3, 5.88e3, 4.25e3, 3.24e3, 2.54e3, 0.])
      self.PlanckPCIB['857_217'] = np.array([0., 0., 9.70e3, 5.26e3, 3.03e3, 1.88e3, 1.31e3, 9.18e2, 7.00e2, 5.38e2, 0.])
      self.PlanckPCIB['857_143'] = np.array([0., 0., 1.84e3, 1.06e3, 6.52e2, 3.86e2, 2.55e2, 1.76e2, 1.23e2, 1.03e2, 0.])
      #
      self.PlanckPCIB['545'] = np.array([0., 0., 6.63e4, 3.34e4, 1.91e4, 1.25e4, 9.17e3, 6.83e3, 5.34e3, 4.24e3, 3.42e3])
      self.PlanckPCIB['545_353'] = np.array([0., 0., 2.22e4, 1.19e4, 6.93e3, 4.61e3, 3.39e3, 2.50e3, 1.93e3, 1.52e3, 0.])
      self.PlanckPCIB['545_217'] = np.array([0., 0., 4.97e3, 2.79e3, 1.65e3, 1.06e3, 7.41e2, 5.38e2, 4.30e2, 3.30e2, 0.])
      self.PlanckPCIB['545_143'] = np.array([0., 0., 1.01e3, 5.98e2, 3.77e2, 2.29e2, 1.54e2, 1.03e2, 7.09e1, 5.89e1, 0.])
      #
      self.PlanckPCIB['353'] = np.array([0., 0., 7.88e3, 4.35e3, 2.60e3, 1.74e3, 1.29e3, 9.35e2, 7.45e2, 6.08e2, 0.])
      self.PlanckPCIB['353_217'] = np.array([0., 0., 1.75e3, 1.02e3, 6.21e2, 3.97e2, 2.87e2, 1.99e2, 1.59e2, 1.35e2, 0.])
      self.PlanckPCIB['353_143'] = np.array([0., 0., 3.61e2, 2.32e2, 1.48e2, 9.42e1, 6.33e1, 4.56e1, 2.77e1, 3.53e1, 0.])
      #
      self.PlanckPCIB['217'] = np.array([0., 0., 4.17e2, 2.62e2, 1.75e2, 1.17e2, 8.82e1, 6.42e1, 3.34e1, 4.74e1, 0.])
      self.PlanckPCIB['217_143'] = np.array([0., 0., 1.04e2, 7.49e1, 5.87e1, 3.93e1, 2.64e1, 2.21e1, 1.07e1, 1.45e1, 0.])
      #
      self.PlanckPCIB['143'] = np.array([0., 0., 3.64e1, 3.23e1, 2.81e1, 2.27e1, 1.84e1, 1.58e1, 1.25e1, 0., 0.])

      # Error bars on the power spectra:
      self.PlanckPCIB['857_error'] = np.array([2.76e6, 7.99e5, 0.37e5, 0.08e5, 0.26e4, 0.18e4, 0.09e4, 0.05e4, 0.03e4, 0.02e4, 0.01e4])
      self.PlanckPCIB['857_545_error'] = np.array([9.73e5, 3.23e5, 0.13e5, 0.30e4, 0.10e4, 0.07e4, 0.04e4, 0.02e4, 0.11e3, 0.07e3, 0.06e3])
      self.PlanckPCIB['857_353_error'] = np.array([2.91e5, 1.05e5, 0.41e4, 0.11e4, 0.06e4, 0.39e3, 0.27e3, 0.17e3, 0.10e3, 0.07e3, 0.])
      self.PlanckPCIB['857_217_error'] = np.array([6.43e4, 2.49e4, 1.22e3, 0.53e3, 0.32e3, 0.22e3, 0.16e3, 0.87e2, 0.23e2, 0.12e2, 0.])
      self.PlanckPCIB['857_143_error'] = np.array([1.81e4, 5.12e3, 0.45e3, 0.12e3, 0.59e2, 0.41e2, 0.30e2, 0.23e2, 0.16e2, 0.15e2, 0.])
      #
      self.PlanckPCIB['545_error'] = np.array([3.74e5, 1.45e5, 0.51e4, 0.12e4, 0.04e4, 0.03e4, 0.17e3, 0.10e3, 0.06e3, 0.04e3, 0.04e3])
      self.PlanckPCIB['545_353_error'] = np.array([1.15e5, 4.84e4, 0.16e4, 0.05e4, 0.23e3, 0.16e3, 0.11e3, 0.07e3, 0.04e3, 0.03e3, 0.])
      self.PlanckPCIB['545_217_error'] = np.array([2.58e4, 1.16e4, 0.48e3, 0.21e3, 0.12e3, 0.09e3, 0.63e2, 0.35e2, 0.12e2, 0.07e2, 0.])
      self.PlanckPCIB['545_143_error'] = np.array([7.04e3, 2.53e3, 0.19e3, 0.67e2, 0.39e2, 0.27e2, 0.19e2, 0.14e2, 1.80e1, 1.73e1, 0.])
      #
      self.PlanckPCIB['353_error'] = np.array([3.68e4, 1.66e4, 0.53e3, 0.18e3, 0.10e3, 0.07e3, 0.05e3, 0.33e2, 0.22e2, 0.16e2, 0.])
      self.PlanckPCIB['353_217_error'] = np.array([8.01e3, 3.82e3, 0.15e3, 0.06e3, 0.38e2, 0.27e2, 0.20e2, 0.14e2, 0.10e2, 0.05e2, 0.])
      self.PlanckPCIB['353_143_error'] = np.array([2.05e3, 8.26e2, 0.62e2, 0.24e2, 0.14e2, 1.06e1, 0.83e1, 0.91e1, 1.11e1, 0.69e1, 0.])
      #
      self.PlanckPCIB['217_error'] = np.array([1.78e3, 8.47e2, 0.47e2, 0.20e2, 0.13e2, 0.10e2, 0.89e1, 1.61e1, 2.15e1, 0.65e1, 0.])
      self.PlanckPCIB['217_143_error'] = np.array([4.74e2, 1.89e2, 0.19e2, 0.81e1, 0.58e1, 0.50e1, 0.52e1, 1.19e1, 1.65e1, 0.54e1, 0.])
      #
      self.PlanckPCIB['143_error'] = np.array([1.55e2, 6.41e1, 0.73e1, 0.35e1, 0.30e1, 0.29e1, 0.35e1, 0.91e1, 1.28e1, 0., 0.])



   ##################################################################################
   
   def plotPCIB(self, name='353'):
      L = self.PlanckPCIB['ell']
      P = self.PlanckPCIB[name]
      sP = self.PlanckPCIB[name+'_error']
      Pshot = self.PlanckPCIBShot[name]
      
      # P
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.errorbar(L, P, yerr=sP, fmt='.', c='k', label=name)
      ax.plot(L, Pshot*np.ones_like(L), 'b--', label=r'quoted shot noise')
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      ax.legend(loc=1, numpoints=1)
      ax.set_xlabel(r'\ell')
      ax.set_ylabel(r'$C_\ell^\text{CIB}$')
      path = "./figures/pn2d/planck15_cib_"+name+".pdf"
      #fig.savefig(path)
      
      plt.show()
