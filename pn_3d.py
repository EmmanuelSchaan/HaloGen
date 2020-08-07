from headers import *

##################################################################################
##################################################################################

class P3dAuto(object):
   
   def __init__(self, U, IHaloModel, Prof, fPnoise=lambda k,z: 0., fTnoise=lambda k,z: 0., name="", doT=False, Vs=1., nProc=1, save=False):
      # copy classes
      self.U = U
      self.IHaloModel = IHaloModel
      self.Prof = Prof
      self.fPnoise = fPnoise
      self.fTnoise = fTnoise
      self.Vs = Vs   # in (Mpc/h)^3
      self.name = str(self.Prof) + name
      self.nProc = nProc
      
      # values of k to evaluate
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      self.Ke = np.genfromtxt("./input/K.txt") # edges of the bins for k
      self.dK = np.genfromtxt("./input/dK.txt") # widths of the bins for k
      # redshifts to evaluate
      self.Z = np.linspace(0., 10., 101)
   
      # create folder if needed
      directory = "./output/pn_3d/"
      if not os.path.exists(directory):
         os.makedirs(directory)

      # power spectrum
      if (save==True):
         self.SaveP()
      self.LoadP()
      
      # trispectrum, if needed
      if doT:
         if (save==True):
            self.SaveT()
         self.LoadT()

   def __str__(self):
      return self.name

   def SaveP(self):
      nZ = len(self.Z)
      nK = len(self.K)
      
      # power spectra
      P1h = np.zeros((nK, nZ))
      P2h = np.zeros((nK, nZ))
      Pshot = np.zeros((nK, nZ))
      dP = np.zeros((nK, nZ))
   
      # precompute the polyspectra
      print "precomputing p3d "+self.name
      with sharedmem.MapReduce(np=self.nProc) as pool:
         for iZ in range(nZ):
            z = self.Z[iZ]
            f = lambda k: self.fP_1h(k, z)
            P1h[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fP_2h(k, z)
            P2h[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fPnoise(k, z)
            Pshot[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fdP(k, z)
            dP[:, iZ] = np.array(pool.map(f, self.K))
            print "done z "+str(iZ)+" of "+str(nZ)

      # save power spectra
      path = "./output/pn_3d/p3d_"+self.name
      np.savetxt(path+"_z.txt", self.Z)
      np.savetxt(path+"_k.txt", self.K)
      np.savetxt(path+"_1h.txt", P1h)
      np.savetxt(path+"_2h.txt", P2h)
      np.savetxt(path+"_shot.txt", Pshot)
      np.savetxt(path+"_dPdd.txt", dP)
      return


   def LoadP(self):
      # read the precomputed power spectra
      path = "./output/pn_3d/p3d_"+self.name
      Z = np.genfromtxt(path+"_z.txt")
      K = np.genfromtxt(path+"_k.txt")
      self.P1h = np.genfromtxt(path+"_1h.txt")
      self.P2h = np.genfromtxt(path+"_2h.txt")
      self.Pshot = np.genfromtxt(path+"_shot.txt")
      self.P = self.P1h + self.P2h
      self.Ptot = self.P1h + self.P2h + self.Pshot
      self.dP = np.genfromtxt(path+"_dPdd.txt")
      
      # interpolate them
      self.fP1hinterp = interp2d(K, Z, self.P1h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fP2hinterp = interp2d(K, Z, self.P2h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fPshotinterp = interp2d(K, Z, self.Pshot.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fPinterp = interp2d(K, Z, self.P.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fPtotinterp = interp2d(K, Z, self.Ptot.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fdPinterp = interp2d(K, Z, self.dP.transpose(), kind='linear', bounds_error=False, fill_value=0.)


   def SaveT(self):
      nZ = len(self.Z)
      nK = len(self.K)

      # trispectra
      T1h = np.zeros((nK, nZ))
      T2h = np.zeros((nK, nZ))
      T4h = np.zeros((nK, nZ))
      Tshot = np.zeros((nK, nZ))
      
      # precompute the polyspectra
      print "precomputing t3d equilateral "+self.name
      with sharedmem.MapReduce(np=self.nProc) as pool:
         for iZ in range(nZ):
            z = self.Z[iZ]
            f = lambda k: self.fT_1h(k, z)
            T1h[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fT_2h(k, z)
            T2h[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fT_4h(k, z)
            T4h[:, iZ] = np.array(pool.map(f, self.K))
            f = lambda k: self.fTnoise(k, z)
            Tshot[:, iZ] = np.array(pool.map(f, self.K))
            print "done z "+str(iZ)+" of "+str(nZ)
   
      # save trispectra
      path = "./output/pn_3d/t3d_"+self.name
      np.savetxt(path+"_z.txt", self.Z)
      np.savetxt(path+"_k.txt", self.K)
      np.savetxt(path+"_1h.txt", T1h)
      np.savetxt(path+"_2h.txt", T2h)
      np.savetxt(path+"_4h.txt", T4h)
      np.savetxt(path+"_shot.txt", Tshot)
      return


   def LoadT(self):
      # read the precomputed trispectra
      path = "./output/pn_3d/t3d_"+self.name
      Z = np.genfromtxt(path+"_z.txt")
      K = np.genfromtxt(path+"_k.txt")
      self.T1h = np.genfromtxt(path+"_1h.txt")
      self.T2h = np.genfromtxt(path+"_2h.txt")
      self.T4h = np.genfromtxt(path+"_4h.txt")
      self.Tshot = np.genfromtxt(path+"_shot.txt")
      self.T = self.T1h + self.T2h + self.T4h
      self.Ttot = self.T + self.Tshot
      
      # interpolate them
      self.fT1hinterp = interp2d(K, Z, self.T1h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fT2hinterp = interp2d(K, Z, self.T2h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fT4hinterp = interp2d(K, Z, self.T4h.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fTshotinterp = interp2d(K, Z, self.Tshot.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fTinterp = interp2d(K, Z, self.T.transpose(), kind='linear', bounds_error=False, fill_value=0.)
      self.fTtotinterp = interp2d(K, Z, self.Ttot.transpose(), kind='linear', bounds_error=False, fill_value=0.)

   ##################################################################################
   # Power spectrum covariance
   
   def LoadCovP(self, z=0.):
      # find closest precomputed redshift
      iZ = np.argmin((self.Z-z)**2)
      
      # number of pairs
      self.Npairs = self.Vs * self.K**2 * self.dK / (2.*np.pi**2)
      
      # variance of mean overdensity, needed for HSV
      R = (3.*self.Vs/(4.*np.pi))**(1./3.)
      self.s2 = self.U.Sigma2(R, z, W3d_sth)
      
      # Diagonal covariances
      self.s2P_PP = 2. * self.P[:,iZ]**2 / self.Npairs
      self.s2P_T = self.T[:,iZ] / self.Vs
      self.s2P_HSV = self.dP[:,iZ]**2 * self.s2
      self.s2P = self.s2P_PP + self.s2P_T + self.s2P_HSV

      # Non-diagonal covariances
      Nk = len(self.K)
      self.CovP = np.zeros((Nk, Nk))
      for ik1 in range(Nk):
         k1 = self.K[ik1]
         for ik2 in range(Nk):
            k2 = self.K[ik2]
            self.CovP[ik1, ik2] = (ik1==ik2) * self.s2P_PP[ik1]
            self.CovP[ik1, ik2] += self.fTnondiag(k1, k2, z) / self.Vs
            self.CovP[ik1, ik2] += self.dP[ik1,iZ] * self.dP[ik2,iZ] * self.s2
         print "done", ik1, "of", Nk
      # invert the matrix
      #self.InvCovP = np.linalg.inv(self.CovP)
      
      
   def plotSNRCumul(self, z=0., plot=True):
      """Cumulative SNR^2 on P
      """
      # load the covariance matrix
      self.LoadCovP(z=z)
      
      # compute power spectrum
      f = lambda k: self.fPinterp(k, z)[0]  # reduce the dim 0 array to a float
      P = np.array(map(f, self.K))
      
      # compute SNR2 as a function of kMax
      Snr2 = np.zeros(len(self.K))
      for ikmax in range(len(self.K)):
         invCovP = np.linalg.inv(self.CovP[:ikmax+1, :ikmax+1])
         Snr2[ikmax] = np.dot( np.dot(P[:ikmax+1], invCovP), P[:ikmax+1])
      
      if plot:
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.loglog(self.K, np.cumsum(self.Npairs)/2., 'k', lw=2)
         ax.loglog(self.K, Snr2, 'r', lw=4)
         #ax.loglog(K, K**3 * self.Vs / (6.*np.pi**2), 'r')
         #ax.semilogx(K, Snr2 / (np.cumsum(self.Npairs)/2.), 'k', lw=3)
         #ax.loglog(K, 4.e3*(K/0.1)**0.7, 'r--', lw=1)
         #
         ax.set_xlabel(r'$k_\text{max}$ [h/Mpc]', fontsize=22)
         ax.set_ylabel(r'Cumulative SNR$^2$', fontsize=22)
         
         plt.show()
         
      return Snr2


   def plotCovP(self, z=0.):
      # load the covariance matrix
      self.LoadCovP(z=z)
      
      # compute power spectrum
      f = lambda k: self.fPinterp(k, z)[0]  # reduce the dim 0 array to a float
      P = np.array(map(f, self.K))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.loglog(self.K, np.sqrt(self.s2P_PP)/P, 'g', lw=2, label=r'PP')
      ax.loglog(self.K, np.sqrt(self.s2P_T)/P, 'b', lw=2, label=r'T')
      ax.loglog(self.K, np.sqrt(self.s2P_HSV)/P, 'r', lw=2, label=r'HSV')
      ax.loglog(self.K, np.sqrt(self.s2P)/P, 'k', lw=2, label=r'total')
      #
      ax.grid()
      ax.legend()
      ax.set_xlabel(r'$k$ [h/Mpc]')
      ax.set_ylabel(r'$\sigma_P / P$')
      
      plt.show()


   ##################################################################################
   # Power spectrum
   
   def fP_1h(self, k, z):
      Profiles = [[self.Prof, 2, k]]
      result = self.IHaloModel.f(0, Profiles, z)
      if not np.isfinite(result):
         result = 0.
      return result

   def fP_2h(self, k, z):
      Profiles = [[self.Prof, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)**2
      result *= self.U.fPlin(k, z)
      if not np.isfinite(result):
         result = 0.
      return result

   def fP(self, K, z):
      result = self.fP_1h(K, z) + self.fP_2h(K, z)
      if not np.isfinite(result):
         result = 0.
      return result

   def fPtot(self, K, z):
      result = self.fP_1h(K, z) + self.fP_2h(K, z) + self.fPnoise(K, z)
      if not np.isfinite(result):
         result = 0.
      return result

   
   # effective bias, such that P2h = beff^2 * Plin
   def bEff(self, k, z):
      Profiles = [[self.Prof, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      return result
   
   
   ##################################################################################
   # Power spectrum response: overdensity
   
   def fdP_1h(self, k, z):
      Profiles = [[self.Prof, 2, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      return result
   
   def fdP_2h(self, k, z):
      Profiles = [[self.Prof, 1, k]]
      I1 = self.IHaloModel.f(1, Profiles, z)
      I2 = self.IHaloModel.f(2, Profiles, z)
      
      result = 2*I1*I2
      result += I1**2 * self.U.dlnPlindDelta(k,z)
      result *= self.U.fPlin(k, z)
      return result
   
   def fdP(self, k, z):
      result = self.fdP_1h(k, z) + self.fdP_2h(k, z)
      return result


   ##################################################################################
   # Power spectrum response: number of halos

   def fndnP_1h(self, k, z, mMin=1.e12):
      Profiles = [[self.Prof, 2, k]]
      result = self.IHaloModel.f(0, Profiles, z, mMin=mMin)
      return result

   def fndnP_2h(self, k, z, mMin=1.e12):
      Profiles = [[self.Prof, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      result *= self.IHaloModel.f(1, Profiles, z, mMin=mMin)
      result *= self.U.fPlin(k, z)
      result *= 2.
      return result

   def fndnP(self, k, z, mMin=1.e12):
      result = self.fndnP_1h(k, z) + self.fndnP_2h(k, z)
      return result


   ##################################################################################
   # Number of halos, and response to overdensity

   def fn(self, z, mMin=1.e12):
      result = self.IHaloModel.f(0, [[self.Prof, 0, 0.01]], 0., mMin=mMin)
      return result

   def fdn(self, z, mMin=1.e12):
      result = self.IHaloModel.f(1, [[self.Prof, 0, 0.01]], 0., mMin=mMin)
      return result

   def fb(self, z, mMin=1.e12):
      result = self.IHaloModel.f(1, [[self.Prof, 0, 0.01]], 0., mMin=mMin)
      result /= self.IHaloModel.f(0, [[self.Prof, 0, 0.01]], 0., mMin=mMin)
      return result


   ##################################################################################
   # Trispectrum: equilateral

   def fT_1h(self, k, z):
      Profiles = [[self.Prof, 4, k]]
      result = self.IHaloModel.f(0, Profiles, z)
      if not np.isfinite(result):
         result = 0.
      return result

   def fT_2h(self, k, z):
      ''' This is the squeezed, azimuthally averaged, diagonal, 2h trispectrum,
      ie T2h(k, -k, k, -k) azimuthally averaged.
      Following Takada Hu 14
      '''
      # term 2h 13
      Profiles = [[self.Prof, 1, k]]
      result13 = self.IHaloModel.f(1, Profiles, z)
      Profiles = [[self.Prof, 3, k]]
      result13 *= self.IHaloModel.f(1, Profiles, z)
      result13 *= self.U.fPlin(k, z)
      result13 *= 4. # permutations
      # term 2h 22
      Profiles = [[self.Prof, 2, k]]
      result22 = self.IHaloModel.f(1, Profiles, z)
      result22 **= 2.
      result22 *= self.U.fPlin(k*np.sqrt(2.), z)  # the sqrt(2.) is from the azimuthal average
      result22 *= 2. # non-zero permutations
      # sum
      result = result13 + result22
      if not np.isfinite(result):
         result = 0.
      return result

   # the 3h term is zero in the collapsed limit,
   # because the squeezed matter bispectrum is zero,
   # and T3h collapsed \propto Bmatter squeezed
#!!! I think this is not exactly true though... argh!

   def fT_4h(self, k, z):
      ''' This is the squeezed, azimuthally averaged, diagonal, 4h trispectrum,
      ie T4h(k, -k, k, -k) azimuthally averaged.
      Following Takada Hu 14
      '''
      Profiles = [[self.Prof, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      result **= 4.
      # from Scoccimarro Zaldarriaga Hui 99:
      # T(k,-k,k,-k) azimuthally averaged in PT \simeq 232/441 * P(k)**3
      result *= 232./441. * self.U.fPlin(k, z)**3
      if not np.isfinite(result):
         result = 0.
      return result

   def fT(self, K, z):
      result = self.fT_1h(K, z)
      if not np.isfinite(result):
         result = 0.
      return result

   def fTnondiag(self, k1, k2, z):
      Profiles = [[self.Prof, 2, k1], [self.Prof, 2, k2]]
      result = self.IHaloModel.f(0, Profiles, z)
      if not np.isfinite(result):
         result = 0.
      return result


   def fT_ssv(self, k1, k2, K, z):
      """This is the difference between the almost-squeezed
      and the exactly squeezed trispectra:
      T(k1, -k1+K, k2, -k2-K) = T(k1,-k1,k2,-k2) + T_ssv,
      where K << k1, k2.
      T_ssv = (dPk1/ddelta) * (dPk2/ddelta) * PK.
      """
      result = self.fdPinterp(k1, z)
      result *= self.fdPinterp(k2, z)
      result *= self.U.fPlin(k, z)
      if not np.isfinite(result):
         result = 0.
      return result


   ##################################################################################


   def plotPInterp(self):
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
   


   def plotP(self, z=0., p3d=None):
   
      # Plin
      f = lambda k: self.U.fPlin(k, z)
      Plin = np.array(map(f, self.K))
      # P1h
      f = lambda k: self.fP_1h(k, z)
      P1h = np.array(map(f, self.K))
      # P2h
      f = lambda k: self.fP_2h(k, z)
      P2h = np.array(map(f, self.K))
      # noise bias
      f = lambda k: self.fPnoise(k, z)
      Pnoise = np.array(map(f, self.K))
      # P1h+P2h
      P = P1h + P2h + Pnoise
      # other power spectrum to be compared
      if p3d is not None:
         f = lambda k: p3d.fPinterp(k, z)
         P3d = np.array(map(f, self.K))
      
      
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, P, 'k', lw=4, label=r'$P_\text{tot}$')
      #
      if p3d is not None:
         ax.loglog(self.K, P3d, 'k--', lw=2, label=str(p3d))
      #
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


      if p3d is not None:
         fig=plt.figure(1)
         ax=fig.add_subplot(111)
         #
         ax.semilogx(self.K, P/P3d, 'r')
         #
         ax.axhline(0.9)
         ax.axhline(1.)
         ax.axhline(1.1)
         #
         ax.set_xlabel(r'$k$ [h/Mpc]')
         ax.set_ylabel(r'ratio with '+str(p3d))

      plt.show()


   def plotdP(self, z=0.):
      
      # dP1h
      f = lambda k: self.fdP_1h(k, z)
      dP1h = np.array(map(f, self.K))
      # dP2h
      f = lambda k: self.fdP_2h(k, z)
      dP2h = np.array(map(f, self.K))
      # dPtot
      dP = dP1h + dP2h
      # P
      f = lambda k: self.fP(k, z)
      P = np.array(map(f, self.K))
      
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.semilogx(self.K, dP/P, 'k', lw=2, label=r'tot')
      ax.semilogx(self.K, dP1h/P, 'r-', label=r'1h')
      ax.semilogx(self.K, dP2h/P, 'b-', label=r'2h')
      #
      ax.grid()
      ax.legend(loc=3)
      ax.set_xlabel(r'k [Mpc/h]')
      ax.set_ylabel(r'$d\ln P(k) / d\delta$')
      ax.set_title(r'Power spectrum')
      #path = "./figures/pn3d/dp_"+self.name+".pdf"
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


   def plotT(self, z=0.):
      
      # 1h
      f = lambda k: self.fT1hinterp(k, z)
      T1h = np.array(map(f, self.K))
      # 2h
      f = lambda k: self.fT2hinterp(k, z)
      T2h = np.array(map(f, self.K))
      # 4h
      f = lambda k: self.fT4hinterp(k, z)
      T4h = np.array(map(f, self.K))
      # noise
      f = lambda k: self.fTshotinterp(k, z)
      Tnoise = np.array(map(f, self.K))
      
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, T1h, 'r-', lw=2, label=r'$T_\text{1h}$')
      ax.loglog(self.K, T2h, 'b-', lw=2, label=r'$T_\text{2h}$')
      ax.loglog(self.K, T4h, 'm-', lw=2, label=r'$T_\text{4h}$')
      ax.loglog(self.K, Tnoise, 'g-', lw=2, label=r'$T_\text{noise}$')
      #
      ax.legend(loc=3)
      ax.set_xlabel(r'k [Mpc/h]')
      ax.set_ylabel(r'$T(k)$ [(Mpc/h)$^3$]')
      #path = "./figures/pn3d/t_"+self.name+"z_"+str(z)+".pdf"
      #fig.savefig(path, bbox_inches='tight')
      
      plt.show()



##################################################################################
##################################################################################


class P3dCross(P3dAuto):
   
   def __init__(self, U, IHaloModel, Prof1, Prof2, fPnoise=lambda k,z: 0., fTnoise=lambda k,z: 0., name="", doT=False, Vs=1., nProc=1, save=False):
      # copy classes
      self.U = U
      self.IHaloModel = IHaloModel
      self.Prof1 = Prof1
      self.Prof2 = Prof2
      self.fPnoise = fPnoise
      self.fTnoise = fTnoise
      self.Vs = Vs   # in (Mpc/h)^3
      self.name = str(self.Prof1) + str(self.Prof2) + name
      self.nProc = nProc
   
      # values of k to evaluate
      self.K = np.genfromtxt("./input/Kc.txt") # center of the bins for k
      # redshifts to evaluate
      self.Z = np.linspace(0., 10., 101)
      
      # create folder if needed
      directory = "./output/pn_3d/"
      if not os.path.exists(directory):
         os.makedirs(directory)
      
      # power spectrum
      if (save==True):
         self.SaveP()
      self.LoadP()
      
      # trispectrum, if needed
      if doT:
         if (save==True):
            self.SaveT()
         self.LoadT()

   def __str__(self):
      return self.name

   def SaveT(self):
      nZ = len(self.Z)
      nK = len(self.K)

      # trispectra
      T1h = np.zeros((nK, nZ))
      T2h = np.zeros((nK, nZ))
      T4h = np.zeros((nK, nZ))
      Tshot = np.zeros((nK, nZ))
      
      # precompute the polyspectra
      print "precomputing t3d "+self.name
      for iZ in range(nZ):
         z = self.Z[iZ]
         for iK in range(nK):
            k = self.K[iK]
            # equilateral trispectra
            T1h[iK, iZ] = self.fT_1h(k, z)
            T2h[iK, iZ] = 0.  #self.fT_2h(k, z)
            T4h[iK, iZ] = 0.  #self.fT_4h(k, z)
            Tshot[iK, iZ] = self.fTnoise(k, z)
         #print "done k "+str(iK)+" of "+str(nK)
         print "done z "+str(iZ)+" of "+str(nZ)
      
      # save trispectra
      path = "./output/pn_3d/t3d_"+self.name
      np.savetxt(path+"_z.txt", self.Z)
      np.savetxt(path+"_k.txt", self.K)
      np.savetxt(path+"_1h.txt", T1h)
      np.savetxt(path+"_2h.txt", T2h)
      np.savetxt(path+"_4h.txt", T4h)
      np.savetxt(path+"_shot.txt", Tshot)
      return

   ##################################################################################
   
   def fP_1h(self, k, z):
      Profiles = [[self.Prof1, 1, k], [self.Prof2, 1, k]]
      result = self.IHaloModel.f(0, Profiles, z)
      return result
   
   def fP_2h(self, k, z):
      Profiles = [[self.Prof1, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      Profiles = [[self.Prof2, 1, k]]
      result *= self.IHaloModel.f(1, Profiles, z)
      result *= self.U.fPlin(k, z)
      return result
   
   def fP(self, K, z):
      result = self.fP_1h(K, z) + self.fP_2h(K, z)
      return result

   ##################################################################################
   # Response of the power spectrum to overdensity

   def fdP_1h(self, k, z):
      Profiles = [[self.Prof1, 1, k], [self.Prof2, 1, k]]
      result = self.IHaloModel.f(1, Profiles, z)
      return result

   def fdP_2h(self, k, z):
      Profiles = [[self.Prof1, 1, k]]
      I1_1 = self.IHaloModel.f(1, Profiles, z)
      I2_1 = self.IHaloModel.f(2, Profiles, z)
      
      Profiles = [[self.Prof2, 1, k]]
      I1_2 = self.IHaloModel.f(1, Profiles, z)
      I2_2 = self.IHaloModel.f(2, Profiles, z)
      
      result = I2_1*I1_2 + I1_1*I2_2
      result += I1_1*I1_2 * self.U.dlnPlindDelta(k,z)
      result *= self.U.fPlin(k, z)
      return result

   def fdP(self, K, z):
      result = self.fdP_1h(K, z) + self.fdP_2h(K, z)
      return result

   ##################################################################################
   
   def fT_1h(self, k, z):
      Profiles = [[self.Prof1, 2, k], [self.Prof2, 2, k]]
      result = self.IHaloModel.f(0, Profiles, z)
      return result
   
   def fT(self, K, z):
      result = self.fT_1h(K, z)
      return result

   ##################################################################################














































#
#
#   def fdP_1h(self, K, z):
#      k = K[0]
#      P1h = self.Mij.f(0, 2, [k, k], z, nbias=1)
#      return P1h
#
#   def fdP_2h(self, K, z):
#      k = K[0]
#      P2h = self.Mij.f(1, 1, [k], z)
#      P2h *= self.Mij.f(1, 1, [k], z, nbias=1)
#      P2h *= self.U.fPlin(k, z)
#      P2h *= 2.
#      return P2h
#
#   def fdP(self, K, z):
#      return self.fdP_1h(K, z) + self.fdP_2h(K, z)
#
#
#   def fndnP_1h(self, K, z):
#      k = K[0]
#      P1h = self.Mij.f(0, 2, [k, k], z, mMin=self.Mij.mMin)
#      return P1h
#
#   def fndnP_2h(self, K, z):
#      k = K[0]
#      P2h = self.Mij.f(1, 1, [k], z, mMin=self.Mij.mMin)
#      P2h *= self.Mij.f(1, 1, [k], z)
#      P2h *= self.U.fPlin(k, z)
#      P2h *= 2.
#      return P2h
#
#   def fndnP(self, K, z, Data):
#      return self.fndnP_1h(K, z) + self.fndnP_2h(K, z)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#   ##################################################################################
#   # Principal function
#   ##################################################################################
#
#
#   # K has to be a valid input for the corresponding n-pt func
#   def f(self, n, K, z, ihalo=0, id=''):
#      
#      # int_{m>mMin} dn/dm * dm
#      if n==0:
#         if id=="d":
#            return self.Mij.f(1, 0, [], z, mMin=self.Mij.mMin)
#         else:
#            return self.Mij.f(0, 0, [], z, mMin=self.Mij.mMin)
#      
#      # power spectrum
#      if n==2:
#         # linear
#         if id=="linear":
#            return self.U.fPlin(K[0], z)
#         # non-linear
#         if id=='':
#            if ihalo==0:
#               return self.fP(K, z)
#            elif ihalo==1:
#               return self.fP_1h(K, z)
#            elif ihalo==2:
#               return self.fP_2h(K, z)
#         # dP/d(delta_b)
#         elif id=="d":
#            if ihalo==0:
#               return self.fdP(K, z)
#            elif ihalo==1:
#               return self.fdP_1h(K, z)
#            elif ihalo==2:
#               return self.fdP_2h(K, z)
#         # n*dP/dn
#         elif id=="ndn":
#            if ihalo==0:
#               return self.fndnP(K, z)
#            elif ihalo==1:
#               return self.fndnP_1h(K, z)
#            elif ihalo==2:
#               return self.fndnP_2h(K, z)
#
#      # bispectrum
#      if n==3:
#         # non-linear
#         if id=='':
#            if ihalo==0:
#               return self.fB(K, z)
#            elif ihalo==1:
#               return self.fB_1h(K, z)
#            elif ihalo==2:
#               return self.fB_2h(K, z)
#            elif ihalo==3:
#               return self.fB_3h(K, z)
#         # dB/d(delta_b)
#         if id=="d":
#            if ihalo==0:
#               return self.fdB(K, z)
#            elif ihalo==1:
#               return self.fdB_1h(K, z)
#            elif ihalo==2:
#               return self.fdB_2h(K, z)
#            elif ihalo==3:
#               return self.fdB_3h(K, z)
#         # n*dB/dn
#         elif id=="ndn":
#            if ihalo==0:
#               return self.fndnB(K, z)
#            elif ihalo==1:
#               return self.fndnB_1h(K, z)
#            elif ihalo==2:
#               return self.fndnB_2h(K, z)
#            elif ihalo==3:
#               return self.fndnB_3h(K, z)
#
#      
#      # trispectrum
#      elif n==4:
#         return self.fT(K, z)
#      
#      # P5
#      elif n==5:
#         return self.fP5(K, z)
#      
#      # P6
#      elif n==6:
#         return self.fP6(K, z)
#
#      return 
#
#
#   ##################################################################################
#   # Power Spectrum functions
#   ##################################################################################
#
#   def fP_1h(self, K, z):
#      k = K[0]
#      P1h = self.Mij.f(0, 2, [k, k], z)
#         #if k > 3. or k < 1.e-4: # WATCH OUT !!!
#         #P1h = 0.
#      return P1h
#   
#   def fP_2h(self, K, z):
#      k = K[0]
#      P2h = self.Mij.f(1, 1, [k], z)**2
#      P2h *= self.U.fPlin(k, z)
#         #if k > 3. or k < 1.e-4: # WATCH OUT !!!
#         #P2h = 0.
#      return P2h
#
#   def fP(self, K, z):
#      P = self.fP_1h(K, z) + self.fP_2h(K, z)
#         #if k > 3. or k < 1.e-4: # WATCH OUT !!!
#         #P = 0.
#      return P
#
#
#   def fdP_1h(self, K, z):
#      k = K[0]
#      P1h = self.Mij.f(0, 2, [k, k], z, nbias=1)
#      return P1h
#
#   def fdP_2h(self, K, z):
#      k = K[0]
#      P2h = self.Mij.f(1, 1, [k], z)
#      P2h *= self.Mij.f(1, 1, [k], z, nbias=1)
#      P2h *= self.U.fPlin(k, z)
#      P2h *= 2.
#      return P2h
#
#   def fdP(self, K, z):
#      return self.fdP_1h(K, z) + self.fdP_2h(K, z)
#
#
#   def fndnP_1h(self, K, z):
#      k = K[0]
#      P1h = self.Mij.f(0, 2, [k, k], z, mMin=self.Mij.mMin)
#      return P1h
#   
#   def fndnP_2h(self, K, z):
#      k = K[0]
#      P2h = self.Mij.f(1, 1, [k], z, mMin=self.Mij.mMin)
#      P2h *= self.Mij.f(1, 1, [k], z)
#      P2h *= self.U.fPlin(k, z)
#      P2h *= 2.
#      return P2h
#   
#   def fndnP(self, K, z, Data):
#      return self.fndnP_1h(K, z) + self.fndnP_2h(K, z)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#   ##################################################################################
#   # Bispectrum: equilateral config
#   ##################################################################################
#
#   # second order PT kernel for delta
#   # here assumes OmM=1
#   def F2(self, K):
#      q1 = K[0]
#      q2 = K[1]
#      theta1_2 = K[2]
#      
#      f2 = 5./7.
#      f2 += 1./2. * (1./q1**2 + 1./q2**2) * q1*q2*cos(theta1_2)
#      f2 += 2./7. * cos(theta1_2)**2
#      return f2
#
#
#   def fB_1h(self, K, z):
#      k1=k2=k3= K[0]
#      B1h = self.Mij.f(0, 3, [k1, k2, k3], z)
#      return B1h
#   
#   def fB_2h(self, K, z):
#      k1=k2=k3= K[0]
#      B2h = self.Mij.f(1, 1, [k1], z)
#      B2h *= self.Mij.f(1, 2, [k2, k3], z)
#      B2h *= self.U.fPlin(k1, z)
#      B2h *= 3.
#      return B2h
#   
#   def fB_3h(self, K, z):
#      k1=k2=k3= K[0]
#      B3h_PT = 2. * self.F2([k1, k2, 2.*np.pi/3.])
#      B3h_PT *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_PT *= self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_PT *= 3.
#      B3h_b2 = self.Mij.f(2, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_b2 *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_b2 *= 3.
#      B3h = B3h_PT + B3h_b2
#      return B3h
#   
#   def fB(self, K, z):
#      return self.fB_1h(K, z) + self.fB_2h(K, z) + self.fB_3h(K, z)
#   
#
#
#   def fdB_1h(self, K, z):
#      k1=k2=k3= K[0]
#      # B1h
#      B1h = self.Mij.f(0, 3, [k1, k2, k3], z, nbias=1)
#      return B1h
#   
#   def fdB_2h(self, K, z):
#      k1=k2=k3= K[0]
#      B2h = self.Mij.f(1, 1, [k1], z, nbias=1) * self.Mij.f(1, 2, [k2, k3], z)
#      B2h += self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 2, [k2, k3], z, nbias=1)
#      B2h *= self.U.fPlin(k1, z)
#      B2h *= 3.
#      return B2h
#   
#   def fdB_3h(self, K, z):
#      k1=k2=k3= K[0]
#      B3h_PT =  self.Mij.f(1, 1, [k1], z, nbias=1) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_PT +=  self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 1, [k2], z, nbias=1) * self.Mij.f(1, 1, [k3], z)
#      B3h_PT +=  self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z, nbias=1)
#      B3h_PT *= 2. * self.F2([k1, k2, 2.*np.pi/3.])
#      B3h_PT *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_PT *= 3.
#      B3h_b2 = self.Mij.f(2, 1, [k1], z, nbias=1) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_b2 += self.Mij.f(2, 1, [k1], z) * self.Mij.f(1, 1, [k2], z, nbias=1) * self.Mij.f(1, 1, [k3], z)
#      B3h_b2 += self.Mij.f(2, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z, nbias=1)
#      B3h_b2 *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_b2 *= 3.
#      B3h = B3h_PT + B3h_b2
#      return B3h
#   
#   def fdB(self, K, z):
#      return self.fdB_1h(K, z) + self.fdB_2h(K, z) + self.fdB_3h(K, z)
#   
#
#
#   def fndnB_1h(self, K, z):
#      k1=k2=k3= K[0]
#      B1h = self.Mij.f(0, 3, [k1, k2, k3], z, mMin=self.Mij.mMin)
#      return B1h
#   
#   def fndnB_2h(self, K, z):
#      k1=k2=k3= K[0]
#      B2h = self.Mij.f(1, 1, [k1], z, mMin=self.Mij.mMin) * self.Mij.f(1, 2, [k2, k3], z)
#      B2h += self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 2, [k2, k3], z, mMin=self.Mij.mMin)
#      B2h *= self.U.fPlin(k1, z)
#      B2h *= 3.
#      return B2h
#   
#   def fndnB_3h(self, K, z):
#      k1=k2=k3= K[0]
#      B3h_PT =  self.Mij.f(1, 1, [k1], z, mMin=self.Mij.mMin) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_PT +=  self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 1, [k2], z, mMin=self.Mij.mMin) * self.Mij.f(1, 1, [k3], z)
#      B3h_PT +=  self.Mij.f(1, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z, mMin=self.Mij.mMin)
#      B3h_PT *= 2. * self.F2([k1, k2, 2.*np.pi/3.])
#      B3h_PT *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_PT *= 3.
#      B3h_b2 = self.Mij.f(2, 1, [k1], z, mMin=self.Mij.mMin) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z)
#      B3h_b2 += self.Mij.f(2, 1, [k1], z) * self.Mij.f(1, 1, [k2], z, mMin=self.Mij.mMin) * self.Mij.f(1, 1, [k3], z)
#      B3h_b2 += self.Mij.f(2, 1, [k1], z) * self.Mij.f(1, 1, [k2], z) * self.Mij.f(1, 1, [k3], z, mMin=self.Mij.mMin)
#      B3h_b2 *= self.U.fPlin(k1, z) * self.U.fPlin(k2, z)
#      B3h_b2 *= 3.
#      B3h = B3h_PT + B3h_b2
#      return B3h
#   
#   def fndnB(self, K, z):
#      return self.fndnB_1h(K, z) + self.fndnB_2h(K, z) + self.fndnB_3h(K, z)
#
#
#   ##################################################################################
#   # T, P5, P6: 1 halo term only
#   # independent on config
#   ##################################################################################
#
#   def fT(self, K, z):
#      T1h = self.Mij.f(0, 4, K, z)
#      return T1h
#      
#   def fP5(self, K, z):
#      P5_1h = self.Mij.f(0, 5, K, z)
#      return P5_1h
#
#   def fP6(self, K, z):
#      P6_1h = self.Mij.f(0, 6, K, z)
#      return P6_1h
#   
#   ##################################################################################
#   # Name, Path, Nhalos, to save things later
#   ##################################################################################
#   
#   def Name(self, n, ihalo=0, id='', nondiag=0):
#      # tag for n
#      stringn=""
#      if n==0:
#         stringn = "N"
#      elif n==2:
#         stringn = "P"
#      elif n==3:
#         stringn = "B"
#      elif n==4:
#         stringn = "T"
#      else:
#         stringn = "P"+str(n)
#      # tag for 1h, ..., nh contribution
#      stringhalo=""
#      if ihalo<>0:
#         stringhalo = "_"+str(ihalo)+"h"
#      # tag for id
#      prefix=''
#      suffix=''
#      if id=='d' or id=='ndn':
#         prefix = id
#      else:
#         suffix = id
#      # tag for diag/non_diag
#      stringdiag=""
#      if nondiag==1:
#         stringdiag = "_nondiag"
#      # name
#      name = prefix+stringn+suffix+stringhalo+stringdiag
#      return name
#
#
#   def Path(self, n, ihalo=0, id='', nondiag=0):
#      name = self.Name(n, ihalo, id=id, nondiag=nondiag)
#      path = "./output/"+name+"_3d.txt"
#      return path
#   
#   # gives halo contributions to be computed (1h, 2h, ...)
#   def Nhalos(self, n):
#      if n in [2, 3]:
#         return range(1, n+1)
#      elif n in [0, 4, 5, 6]:
#         return [0]
#
#
#
#
#   # computes and save all the halo contributions to the 3d n-pt function
#   # at z=0
#   def SavePn3d(self, n, id='', nondiag=0):
#      # redshift of observer
#      z = 0.  #1./self.U.a_obs - 1.
#      K0 = np.ones(n)
#      
#      print "- computing "+self.Name(n, ihalo=0, id=id, nondiag=nondiag)
#      
#      if nondiag==0:
#         pn_tot = np.zeros(self.Nk)
#         for ihalo in self.Nhalos(n):
#            f = lambda k: self.f(n, k*K0, z, ihalo=ihalo, id=id)
#            pn = np.array(map(f, self.K))
#            # N needs an extra factor of Vs
#            if n==0:
#               pn *= self.Vs
#            np.savetxt(self.Path(n, ihalo=ihalo, id=id), pn)
#            pn_tot += pn
#         np.savetxt(self.Path(n, id=id), pn_tot)
#   
#      elif nondiag==1:
#         pn_tot = np.zeros((self.Nk, self.Nk))
#         pn = np.zeros((self.Nk, self.Nk))
#         for ihalo in self.Nhalos(n):
#   
#            for ik1 in range(self.Nk):
#               k1 = self.K[ik1]
#               for ik2 in range(self.Nk):
#                  k2 = self.K[ik2]
#                  
#                  # arguments for T, P5, P6
#                  K1 = ones(int(n/2.))
#                  K2 = ones(n-int(n/2.))
#                  K = np.append( k1*K1, k2*K2 )
#                  
#                  pn[ik1, ik2] = self.f(n, K, z, ihalo=ihalo, id=id)
#                  # N needs an extra factor of Vs
#                  if n==0:
#                     pn *= self.Vs
#
#            np.savetxt(self.Path(n, ihalo=ihalo, id=id, nondiag=1), pn)
#            pn_tot += pn
#         np.savetxt(self.Path(n, id=id, nondiag=1), pn_tot)
#      
#      return
#   
#   
#   
#   
#   ##################################################################################
#   # class constructor
#   ##################################################################################
#
#
#      
#
#
#   ##################################################################################
#   # tests
#   ##################################################################################
#
#   def plotPn(self):
#      
#      # P
#      fig0 = plt.figure(0)
#      ax = plt.subplot(111)
#      ax.loglog(self.U.K, self.U.Plin, 'g--', label=r'lin')
#      ax.loglog(self.K, self.P_1h, 'r--', label=r'1h')
#      ax.loglog(self.K, self.P_2h, 'b--', label=r'2h')
#      ax.loglog(self.K, self.P, 'k', label=r'tot')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'$P(k)$')
#      ax.set_title(r'Power spectrum')
#      '''
#      # B
#      fig1 = plt.figure(1)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.B_1h, '--', label=r'1h')
#      ax.loglog(self.K, self.B_2h, '--', label=r'2h')
#      ax.loglog(self.K, self.B_3h, '--', label=r'3h')
#      ax.loglog(self.K, self.B, 'k', label=r'tot')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'equilateral $B(k)$')
#      ax.set_title(r'Equilateral bispectrum')
#      
#      # T
#      fig2 = plt.figure(2)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.T, 'r--', label=r'1h')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'square $T(k)$')
#      ax.set_title(r'Square trispectrum')
#
#      # P5
#      fig3 = plt.figure(3)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.P5, 'r--', label=r'1h')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'diagonal $P_5(k)$')
#      ax.set_title(r'Diagonal $P_5(k)$')
#      
#      # P6
#      fig3 = plt.figure(4)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.P6, 'r--', label=r'1h')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'diagonal $P_6(k)$')
#      ax.set_title(r'Diagonal $P_6(k)$')
#      
#      
#      # dP
#      fig0 = plt.figure(5)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.dP_1h, '--', label=r'1h')
#      ax.loglog(self.K, self.dP_2h, '--', label=r'2h')
#      ax.loglog(self.K, self.dP, 'k', label=r'tot')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'$dP(k) / d\delta_b$')
#      ax.set_title(r'deriv. power spectrum wrt $\delta_b$')
#      
#      # dB
#      fig0 = plt.figure(6)
#      ax = plt.subplot(111)
#      ax.loglog(self.K, self.dB_1h, '--', label=r'1h')
#      ax.loglog(self.K, self.dB_2h, '--', label=r'2h')
#      ax.loglog(self.K, self.dB_3h, '--', label=r'3h')
#      ax.loglog(self.K, self.dB, 'k', label=r'tot')
#      ax.grid()
#      ax.legend(loc=3)
#      ax.set_xlabel(r'k [Mpc/h]')
#      ax.set_ylabel(r'$dB(k) / d\delta_b$')
#      ax.set_title(r'deriv. equi. bispectrum wrt $\delta_b$')
#      '''
#      plt.show()
#   
#   
#   
#   
#   
#   
#   
#   def plotN_volume(self):
#      
#      z = 0.
#      # survey volumes [(Mpc/h)^3]
#      V = np.logspace(log10(0.01), log10(10.), 51, 10.) * (1.e3)**3
#      # nb of halos of mass > mMin
#      f = lambda v: v*self.f(0, [0.], z)
#      N = np.array(map(f, V))
#      # HSV on N
#      r = lambda v: (3.*v/(4.*np.pi))**(1./3.)
#      f = lambda v: self.U.Sigma2(r(v), z, W3d_sth) * (v*self.f(0, [0.], z, id='d'))**2
#      HSV = np.array(map(f, V))
#      
#      fig=plt.figure(0)
#      ax=plt.subplot(111)
#      ax.loglog(V/(1.e3)**3, np.sqrt(N)/N, 'k--', label=r'Poisson $\propto 1/\sqrt{V}$')
#      ax.loglog(V/(1.e3)**3, np.sqrt(HSV)/N, 'b--', label=r'HSV')
#      ax.loglog(V/(1.e3)**3, np.sqrt(N+HSV)/N, 'r', label=r'total')
#      #ax.loglog(V/(1.e3)**3, np.sqrt(1./V) * np.sqrt(N[0])/N[0]/np.sqrt(1./V)[0], 'k.', label=r'$1/\sqrt{V}$')
#      ax.legend()
#      ax.set_xlabel(r'survey volume [(Gpc/h)$^3$]')
#      ax.set_ylabel(r'$\sigma_N/N$')
#      plt.savefig('./figures/cov_3d/errorN_volume.pdf')
#
#      plt.show()


