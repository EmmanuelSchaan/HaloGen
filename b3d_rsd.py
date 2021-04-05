from headers import *

##################################################################################
##################################################################################

class B3dRsdAuto(object):
   
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
      directory = "./output/b3d_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)

   def __str__(self):
      return self.name


   ##################################################################################
   # Power spectrum ingredients


   def b1h(self, k1, k2, k3, z, mu1=0., mu2=0., mu3=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.u(k1, m, z, mu1)
         result *= self.Prof.u(k2, m, z, mu1)
         result *= self.Prof.u(k3, m, z, mu1)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result



##################################################################################
##################################################################################
##################################################################################
##################################################################################

class B3dRsdCross(P3dRsdAuto):

   def __init__(self, U, Prof, Prof2, Prof3, MassFunc, name="", nProc=1):
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
      directory = "./output/b3d_rsd/"
      if not os.path.exists(directory):
         os.makedirs(directory)


   ##################################################################################
   # Power spectrum ingredients


   def b1h(self, k1, k2, k3, z, mu1=0., mu2=0., mu3=0., mMin=0., mMax=np.inf):
      '''Includes FOGs
      '''
      def integrand(lnm):
         m = np.exp(lnm)
         result = self.MassFunc.massfunc(m, z)
         result *= self.Prof.u(k1, m, z, mu1) 
         result *= self.Prof2.u(k2, m, z, mu2)
         result *= self.Prof3.u(k3, m, z, mu3)
         result *= m # because integrating in lnm and not m
         return result
      # integration bounds
      mMin = max(self.mMin, mMin)
      mMax = min(self.mMax, mMax)
      result = integrate.quad(integrand, np.log(mMin), np.log(mMax), epsabs=0, epsrel=1.e-2)[0]
      return result


