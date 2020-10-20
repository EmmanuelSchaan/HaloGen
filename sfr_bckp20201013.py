from headers import *

class Sfr(object):

   def __init__(self, U, MassFunc):
      self.U = U
      self.MassFunc = MassFunc

      # from Fonseca+16 (1607.05288v2).
      # Table I
      Z = np.array([0., 1., 2., 3., 4., 5.])
      self.Z = Z.copy()
      M0 = np.array([3.0e-10, 1.7e-9, 4.0e-9, 1.1e-8, 6.6e-8, 7.0e-7])  # [Msun/yr]
      Mb = np.array([6.0e10, 9.0e10, 7.0e10, 5.0e10, 5.0e10, 6.0e10])   # [Msun]
      Mc = np.array([1.0e12, 2.0e12, 2.0e12, 3.0e12, 2.0e12, 2.0e12])   # [Msun]
      a = np.array([3.15, 2.9, 3.1, 3.1, 2.9, 2.5])
      b = np.array([-1.7, -1.4, -2.0, -2.1, -2.0, -1.6])
      c = np.array([-1.7, -2.1, -1.5, -1.5, -1.0, -1.0])
      # interpolate values
      self.fM0 = interp1d(Z, M0, kind='linear', bounds_error=False, fill_value=0.) # [Msun/yr]
      self.fMb = interp1d(Z, Mb, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
      self.fMc = interp1d(Z, Mc, kind='linear', bounds_error=False, fill_value=0.) # [Msun]
      self.fa = interp1d(Z, a, kind='linear', bounds_error=False, fill_value=0.)
      self.fb = interp1d(Z, b, kind='linear', bounds_error=False, fill_value=0.)
      self.fc = interp1d(Z, c, kind='linear', bounds_error=False, fill_value=0.)
      # below Eq 11
      self.Ma= 1.e8   # [Msun]

      # SFRD interpolation
      #z = np.logspace(np.log10(1.e-3), np.log10(10.), 501, 10.)
      z = np.linspace(0., 10., 501)
      sfrd = np.array(map(self.sfrdForInterp, z))
      self.sfrd = interp1d(z, sfrd, kind='linear', bounds_error=False, fill_value=0.)


   def sfr(self, m, z):
      ''' SFR [Msun/yr] as a function of halo mass [Msun] and redshift
      from Fonseca+16 (1607.05288v2), Eq 11.
      I inferred the units from Eq 9 and 11.
      input halo mass m [Msun/h]
      '''
      # bounds for the fitting function
      if z<0. or z>5.:
         return 0.

      # convert from [Msun/h] to [Msun]
      m /= self.U.bg.h

      # Eq 11
      result = self.fM0(z)
      result *= (m/self.Ma)**self.fa(z)
      result *= (1.+m/self.fMb(z))**self.fb(z)
      result *= (1.+m/self.fMc(z))**self.fc(z)

      return result


   def sfrdForInterp(self, z, alpha=1):
      '''SFR density:
      \int dm dn/dm SFR(m)^alpha [(Msun/yr)^alpha / (Mpc/h)^3]
      '''
      def integrand(lnm):
         '''the mass in integral is in [Msun/h]
         '''
         m = np.exp(lnm)
         #return m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m / self.U.bg.h, z)
         return m * self.MassFunc.fmassfunc(m, 1./(1.+z)) * self.sfr(m, z)**alpha
      # [Msun /yr / (Mpc/h)^3]
      result = integrate.quad(integrand, np.log(self.MassFunc.mMin), np.log(self.MassFunc.mMax), epsabs=0., epsrel=1.e-3)[0]
      # [Msun /yr / Mpc^3]
      #result *= self.U.bg.h**3 
      return result

   
   def sfrdBehroozi13(self, z):
      # from Fonseca+16
      # [Msun/yr/Mpc^3]
      result = 0.180 / (10.**(-0.997*(z-1.243)) + 10.**(0.241*(z-1.243))) 
      # convert from [Msun/yr/Mpc^3] to [Msun/yr/(Mpc/h)^3]
      result /= self.U.bg.h**3
      return result

   def testSfrd(self):
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # SMill fitting function
      Sfrd = np.array(map(self.sfrd, self.Z))
      ax.plot(self.Z, np.log10(Sfrd), 'b--', label=r'Fonseca+16')
      #
      # Behroozi+13 parametrization
      z = np.linspace(0., 5., 101)
      ax.plot(z, np.log10(self.sfrdBehroozi13(z)), 'r-', label=r'Behroozi+13')
      #
      ax.legend(loc='lower center')
      ax.set_xlabel(r'$z$')
      ax.set_ylabel(r'log$_{10}$(SFRD$(z)$/[M$_\odot$/yr / (Mpc/h)$^3$])')
      #
      #fig.savefig(pathFig + "sfrd.pdf", bbox_inches='tight')
      plt.show()
      fig.clf()


   def nHEff(self, z):
      '''Effective halo number density [(Mpc/h)^-3]
      which determines the amplitude of the halo shot noise,
      ie the 1-halo term
      '''
      result = self.sfrdForInterp(z, alpha=1)**2
      result /= self.sfrdForInterp(z, alpha=2)
      return result

   def plotnHEff(self):

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      n = np.array(map(self.nHEff, self.Z))
      ax.plot(self.Z, n)
      #
      plt.show()

