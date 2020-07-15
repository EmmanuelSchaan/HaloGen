from headers import *

##################################################################################
#  Mathematical functions



def W3d_sth(x):
   """Fourier transform of a 3d spherical top hat window function.
   Use x = k*R as input,
   where R is the tophat radius and k the wave vector.
   Input and output are dimensionless.
   """
   if x < 1.e-3:  # for small x, replace by expansion, for numerical stability
      f = 1. - 0.1* x**2 + 0.00357143* x**4
   else:
      f = (3./(x**3)) * ( np.sin(x) - x * np.cos(x) )
   return f


def dW3d_sth(x):
   """Derivative of the FT of the top hat.
   Input and output are dimensionless.
   """
   f = 3. * (3. * x * np.cos(x) - 3. * np.sin(x) + (x**2) * np.sin(x)) / (x**4)
   return f


def W2d_cth(x):
   """FT of a 2d circular top hat window function.
   Input and output are dimensionless.
   """
   return 2.*special.jn(1, x) / x

def W1d_th(x):
   """FT of a 1d tophat
   normalized to unity at k=0 (ie real space integral is 1)
   Input and output are dimensionless.
   """
   return sinc(x/2.)
   
def Si(x):
   return special.sici(x)[0]

def Ci(x):
   return special.sici(x)[1]

def sinc(x):
   return special.sph_jn(0, x)[0][0]

def j0(x):
   """relevant for isotropic Fourier transform in 2d
   """
   return special.jn(0, x)


def i0(x):
   """Modified Bessel function of the first kind
   """
   return special.iv(0, x)


##################################################################################
# formatting numbers

def intExpForm(input):
   """
   clean scientific notation for file names
   removes trailing decimal point if not needed
   """
   a = '%e' % np.float(input)
   # mantissa: remove trailing zeros
   # then remove dot if no decimal digits
   mantissa = a.split('e')[0].rstrip('0').rstrip('.')
   # exponent: remove + sign if there, and leading zeros
   exponent = np.int(a.split('e')[1])
   exponent = np.str(exponent)
   if exponent=='0':
      return mantissa
   else:
      return mantissa + 'e' + exponent



def floatExpForm(input):
   """same as intExpForm, except always leaves the decimal point
   """
   a = '%e' % np.float(input)
   # mantissa: remove trailing zeros
   # then remove dot if no decimal digits
   mantissa = a.split('e')[0].rstrip('0')
   # exponent: remove + sign if there, and leading zeros
   exponent = np.int(a.split('e')[1])
   exponent = np.str(exponent)
   if exponent=='0':
      return mantissa
   else:
      return mantissa + 'e' + exponent


##################################################################################
# Matrix inversion for ill-conditioned matrices, with SVD

def invertMatrixSvdTruncated(matrix, epsilon=1.e-5, keepLow=True):
   '''Invert a matrix by inverting its SVD,
   and setting to zero the singular values that are too small/large.
   epsilon sets the tolerance for discarding singular values.
   keepLow=True: for inverting a cov matrix, want to keep the modes with lowest variance.
   keepLow=False: for inverting a Fisher maytrix, want to keep the modes with highest Fisher information.
   '''
   # Perform SVD on matrix
   U, s, Vh = scipy.linalg.svd(matrix)
   # invert the singular values
   sInv = 1./s
   # remove the super poorly constrained modes, that lead to numerical instabilities
   if keepLow:
      sInvMax = np.max(sInv)
      sInv[sInv<=sInvMax*epsilon] = 0.
   else:
      sInvMin = np.min(sInv)
      sInv[sInv>=sInvMin/epsilon] = 0.
   # create diagonal matrix
   sInv = np.diag(sInv)
   # invert the hermitian matrices
   V = np.conj(Vh.transpose())
   Uh = np.conj(U.transpose())
   # generate the inverse
   result = np.dot(V, np.dot(sInv, Uh))
   return result


##################################################################################
# Generating ell bins with constant number of modes


def generateEllBins(lMin, lMax, nL, fsky=1.):
   '''Generates nL bins between lMin and lMax,
   such that the number of 2d modes in each bin is identical.
   Returns the bin centers, the bin edges, the bin widths, and the number of modes per bin.
   '''
   # area in ell space between lMin and l,
   # normalized to 1 when l=lMax
   farea = lambda l: (l**2 - lMin**2) / (lMax**2 - lMin**2)

   # find the bin edges,
   # such that each bin has equal number of modes
   Le = np.zeros(nL+1)
   for iL in range(nL+1):
      f = lambda l: farea(l) - float(iL) / nL
      Le[iL] = optimize.brentq(f , lMin, lMax)
   
   # use the average ell in the bin, weighted by number of modes, as bin center
   Lc =  2./3. * (Le[1:]**3 - Le[:-1]**3) / (Le[1:]**2 - Le[:-1]**2)
   # bin spacing
   dL = Le[1:] - Le[:-1]
   # Number of modes
   lF = 2.*np.pi / np.sqrt(4. * np.pi * fsky)
   nModesTotal = np.pi * (lMax**2 - lMin**2) / lF**2
   Nmodes = nModesTotal / nL * np.ones(nL)   
   
   return Lc, dL, Nmodes, Le


##################################################################################
#  Extract non-mask data vector and cov matrix

def extractMaskedMat(cov, mask=None, I=None):
   '''cov: large matrix
   mask: 1d array, 0 if unmasked, anything else if masked
   I: indices of the large matrix to keep, pre-masking
   '''
   # convert mask to 0 and 1
   mask = mask.astype(bool)
   # extract indices of interest, if needed
   if I is not None:
      mask = mask[I]
      J = np.ix_(I, I)
      cov = cov[J]
   if mask is not None:
      # nb of unmasked rows
      nNew = np.int(np.sum(1 - mask))
      # mask cov matrix
      mask = np.diag(mask)
      cov = ma.masked_array(cov, mask=mask)
      cov = ma.mask_rowcols(cov)
      # extract the non-masked elements
      cov = cov.compressed().reshape((nNew, nNew))
   return cov

def extractMaskedVec(vec, mask=None, I=None):
   '''vec: large vector
   mask: 1d array, 0 if unmasked, anything else if masked
   I: indices of the large vector to keep, pre-masking
   '''
   # convert mask to 0 and 1
   mask = mask.astype(bool)
   # extract indices of interest, if needed
   if I is not None:
      mask = mask[I]
      vec = vec[I]
   if mask is not None:
      # mask vec matrix
      vec = ma.masked_array(vec, mask=mask)
      # extract the non-masked elements
      vec = vec.compressed()
   return vec

##################################################################################
# Measuring RAM usage of the current process

def currentRssMB(pid=None):
   """Returns the RSS (resident set size, ie portion of RAM occupied by a process).
   Output in MB.
   """
   memory = dict(psutil.Process(pid=None).memory_info()._asdict())
   result = memory['rss']  # in Bytes
   return result / 1.e6

##################################################################################

# !!! works only on scalars, not on arrays
#def divide(x, y, exceptOut=0.):
#   '''Returns 0. or any requested value
#   when dividing by zero.
#   '''
#   try: return x/y
#   except ZeroDivisionError: return exceptOut


def divide(x, y, exceptOut=0.):
   '''Returns 0. or any requested value
   when dividing by zero.
   Works both on scalars and arrays.
   '''
   # if both inputs are scalar, make one an array
   inputScalar = False
   if np.shape(x)==() and np.shape(y)==():
      inputScalar = True
      x = np.array([x])
   # do the division
   result = x / y
   # where it went wrong, replace value with exceptOut
   result[np.where(np.isfinite(result)==False)] = exceptOut
   # if inputs were both scalar, make the result a scalar
   if inputScalar:
      result = result[0]
   return result


##################################################################################

def myHistogram(X, nBins=71, lim=None, S2Theory=[], path=None, plot=False, nameLatex=r'$x$', semilogx=False, semilogy=False, doGauss=False):
   """Generic histogram plotter.
   Flattens the input array X first thing.
   """
   # Flatten the array in case
   X = X.flatten()
   # value limits for the histogram
   if lim is None:
      lim = (np.min(X), np.max(X))
   # Bin edges
   if semilogx:
      Bins = np.logspace(np.log10(lim[0]), np.log10(lim[1]), nBins, 10.)
   else:
      Bins = np.linspace(lim[0], lim[1], nBins)
   binwidth = Bins[1:] - Bins[:-1]

   # Data histogram
   histX = np.histogram(X, Bins)[0]
   histX = histX.astype(float)

   # Plot
   fig = plt.figure(0)
   ax = fig.add_subplot(111)
   #
   # Histogram from data
   ax.bar(Bins[:-1], histX, binwidth, color='b', alpha=0.5, label=r'Data')
   #
   # histogram for a Gaussian with the variance from the data
   if doGauss:
      mean = np.mean(X)
      std = np.std(X)
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
      histGaussFit = np.array(map(g, range(nBins-1)))
      histGaussFit *= len(X)
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
   #
   # Theory histogram
   for s2Theory in S2Theory:
      av = 0.
      fPDF = lambda x: (2*np.pi*s2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*s2Theory))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
      histTheory = np.array(map(g, range(nBins-1)))
      histTheory *= len(X)
      ax.step(Bins[:-1], histTheory, color='r', lw=3, where='post')#, label=r'Theory')
   #
   ax.legend(loc=1)
   ax.set_xlim((lim[0], lim[1]))
   if semilogx:
      ax.set_xscale('log', nonposx='clip')
   if semilogy:
      ax.set_yscale('log', nonposy='clip')
   ax.set_ylim((0.5*np.min(histX[histX>0]), 2.*np.max(histX)))
   ax.set_xlabel(nameLatex)
   ax.set_ylabel(r'number of objects')
   if path is not None:
      fig.savefig(path, bbox_inches='tight')
   if plot:
      plt.show()
   else:
      fig.clf()



def my2dHistogram(X, Y, nBins=(71, 71), limx=None, limy=None, limc=None, fTheory=[], path=None, plot=False, nameLatexX=r'$x$', nameLatexY=r'$y$', logx=False, logy=False, logColor=False, cmap=plt.cm.jet):
   """Generic 2d histogram plotter.
   """
   # limits for bins and colors
   if limx is None:
      limx = (np.min(X), np.max(X))
   if limy is None:
      limy = (np.min(Y), np.max(Y))
   if limc is None:
      limc = (None, None)
   # x-bin edges
   if logx:
      BinsX = np.logspace(np.log10(limx[0]), np.log10(limx[1]), nBins[0], 10.)
   else:
      BinsX = np.linspace(limx[0], limx[1], nBins[0])
   # y-bin edges
   if logy:
      BinsY = np.logspace(np.log10(limy[0]), np.log10(limy[1]), nBins[1], 10.)
   else:
      BinsY = np.linspace(limy[0], limy[1], nBins[1])

   # Plot
   fig = plt.figure(0)
   ax = fig.add_subplot(111)
   #
   # 2d histogram
   H, xEdges, yEdges = np.histogram2d(X, Y, bins=[BinsX, BinsY], range=[[limx[0], limx[1]],[limy[0], limy[1]]])
   H = H.T
   x, y = np.meshgrid(xEdges, yEdges)
   if logColor:
      im = ax.pcolormesh(x, y, H, cmap=cmap, linewidth=0, rasterized=True, norm=LogNorm(limc[0], limc[1]))
   else:
      im = ax.pcolormesh(x, y, H, cmap=cmap, linewidth=0, rasterized=True, vmin=limc[0], vmax=limc[1])
   #
   # theory curves, if any
   for f in fTheory:
      Ytheory = np.array(map(f, X))
      ax.plot(X, Ytheory, 'y')
   #
   plt.colorbar(im)
   if logx:
      ax.set_xscale('log', nonposx='clip')
   if logy:
      ax.set_yscale('log', nonposy='clip')
   ax.set_xlim((limx[0], limx[1]))
   ax.set_ylim((limy[0], limy[1]))
   ax.set_xlabel(nameLatexX)
   ax.set_ylabel(nameLatexY)
   if path is not None:
      fig.savefig(path, bbox_inches='tight')
   if plot:
      plt.show()
   else:
      fig.clf()
