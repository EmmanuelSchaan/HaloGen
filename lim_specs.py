from headers import *



class LimSpecs(object):
   
   def __init__(self, U, exp='SPHEREx'):

      self.U = U
      self.exp = exp

      if exp=='SPHEREx':
         self.RR = np.array([40., 150., 300.])
         self.fwhmPsf = 6.*np.pi/(180.*3600.)  # 6'' in [rad]
         self.fSkyExp = 2. * 100. * (np.pi/180.)**2 / (4.*np.pi) # 2 * 100 deg2 deep fields
         self.pixelAngularArea = (6.*np.pi/(180.*3600.))**2  # 6''*6'' in [sr]


      elif exp=='COMAPPathfinder':
         self.RR = np.array([800.])
         self.fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
         self.fSkyExp = 2.5 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
         self.pixelAngularArea = (self.fwhmPsf / 2.)**2 # assume two pixels per PSF FWHM

      elif exp=='COMAP':
         self.RR = np.array([3000.])   #np.array([800.])
         self.fwhmPsf = 3.*np.pi/(180.*60.) # 3' in [rad]
         self.fSkyExp = 6.25 * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
         self.pixelAngularArea = (self.fwhmPsf / 2.)**2 # assume two pixels per PSF FWHM



      elif exp=='CONCERTO':
         self.RR = np.array([300.])
         self.fwhmPsf = 0.24 * np.pi/(180.*60.) # 3' in [rad]
         self.fSkyExp = 2. * (np.pi/180.)**2 / (4.*np.pi) # 2.5 deg^2
         self.pixelAngularArea = (self.fwhmPsf / 2.)**2 # assume two pixels per PSF FWHM


      self.R = self.RR[0]
      self.surveyAngularArea = self.fSkyExp * 4.*np.pi


   def voxelComovingDepth(self, z, R=None):
      '''Compute dchi [Mpc/h]
      the comoving depth of the voxel
      at the requested redshift
      and spectral resolving factor
      '''
      if R is None:
         R = self.R
      result = self.U.c_kms / self.U.hubble(z)
      result *= (1.+z) / R
      return result
      

   def pixelComovingArea(self, z):
      '''Pixel comoving area [(Mpc/h)^2]
      at the requested redshift
      '''
      return self.pixelAngularArea * self.U.bg.comoving_distance(z)**2


   def voxelComovingVolume(self, z, R=None):
      '''Voxel comoving volume [(Mpc/h)^3]
      '''
      result = self.voxelComovingDepth(z, R=R)
      result *= self.pixelComovingArea(z)
      return result


#   def dBdT(self, nu, T):
#      '''d(blackbody)/dT, such that
#      dI = d(blackbody)/dT * dT
#      input: nu [Hz], T thermo temperature of the black body [K]
#      output in SI: [W / Hz / m^2 / sr / K]
#      '''
#      x = self.h*nu/(self.kB*T)
#      result = 2.*self.h**2*nu**4
#      result /= self.kB*T**2*self.c**2
#      result *= np.exp(x) / (np.exp(x) - 1.)**2
#      return result


   def whiteNoisePower(self, z, R=None):
      '''White noise power spectrum [(Jy/sr)^2 * (Mpc/h)^3]
      '''
      if self.exp=='SPHEREx':
         # Inferred from SPHEREx science book
         mAB5Sigma = 22 # 5sigma lim mag for point source (Cheng+18)
         f5Sigma = 10.**((8.9-mAB5Sigma)/2.5)   # 5sigma lim point source flux [Jy]
         sigmaFSource = f5Sigma / 5. # 1sigma lim point source flux [Jy]
         # This point source flux is the output of a spatial matched filter
         # for one frequency element.
         # It needs to be converted to pixel flux.
         # The SPHEREx doc, fig 9, gives the conversion using 
         # the effective number of pixels covered by the PSF
         nPixEff = 3.   # 2-5 in fig 9 of SPHEREx doc
         sigmaFPixel = sigmaFSource / np.sqrt(nPixEff)   # [Jy]
         sigmaIPixel = sigmaFPixel / self.pixelAngularArea  # [Jy/sr]
         # convert from pixel variance [(Jy/sr)^2]
         # to white noise power spectrum [(Jy/sr)^2 * (Mpc/h)^3]
         result = sigmaIPixel**2 * self.voxelComovingVolume(z, R=R)

      elif self.exp=='COMAP':
         # Full configuration, from Li Wechsler+16
         tSys = 35   # system temperature [K]
         nFeed = 500.  # number of feeds
         dnu = 10.e6 # spectral element width [Hz]
         nuCenter = 32.e9  # 30-34GHz is the COMAP band [Hz]
         tObs = 2250 * 3600.  # total observing time [s]
         # observing time per pixel [s]
         tPixel = tObs * self.pixelAngularArea / self.surveyAngularArea
         # radiometer equation (App C1 in Li Wechsler+16)
         # giving the pixel noise standard deviation [K]
         sigmaIPixel = tSys / np.sqrt(nFeed * dnu * tPixel) # [K]
         # convert from [K] to [Jy/sr]
         # using the Rayleigh-Jeans temperature definition
         kB = 1.38e-23  # [SI] = [m^2*kg/s^2/K]
         sigmaIPixel *= 2. * nuCenter**2 * kB / 299792458.**2  # [W/m^2/sr/Hz]
         # convert to [Jy/sr]
         sigmaIPixel /= 1.e-26   # [Jy/sr]
         # convert from pixel variance [(Jy/sr)^2]
         # to white noise power spectrum [(Jy/sr)^2 * (Mpc/h)^3]
         result = sigmaIPixel**2 * self.voxelComovingVolume(z, R=R)

      elif self.exp=='CONCERTO':
         # From Dumitru+19, from Serra+16
         lambdaCii = 158.e-6  # rest wavelength[m]
         lambdaObsCii = lambdaCii * (1.+z)   # obs wavelength [m]
         D = 12.  # telescope aperture diameter [m]
         airyDiskRadius = 1.22 * lambdaObsCii / D  # [rad]; fwhm would be 1.00 lambda/D
         beamSolidAngle = 2.*np.pi * (airyDiskRadius/2.355)**2 # [sr]
         # pixel noise [Jy/sr*sqrt(s)]
         sigmaPixel = 155.e-3 / beamSolidAngle # [Jy/sr*sqrt(s)]
         # compute observing time per pixel
         tSurvey = 1500.*3600.   # [s]
         nPixel = 1500.
         surveyAngularArea = self.fSkyExp * 4.*np.pi  # [sr]
         # I don't understand the formula below (scaling with nPix...)
         # but it is in Serra+16 and Dumitru+19
         tPixel = tSurvey * nPixel * beamSolidAngle / surveyAngularArea   # [s]
         # compute white noise power spectrum
         result = sigmaPixel**2 / tPixel * self.voxelComovingVolume(z, R=R)

      return result

