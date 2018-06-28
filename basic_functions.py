from headers import *

##################################################################################
#  Mathematical functions


# FT of a 3d spherical top hat window function
def W3d_sth(x):
   if x < 1.e-3:  # for small x, replace by expansion, for numerical stability
      f = 1. - 0.1* x**2 + 0.00357143* x**4
   else:
      f = (3./(x**3)) * ( np.sin(x) - x * np.cos(x) )
   return f

# derivative of the FT of the top hat
def dW3d_sth(x):
   f = 3. * (3. * x * np.cos(x) - 3. * np.sin(x) + (x**2) * np.sin(x)) / (x**4)
   return f

# FT of a 2d circular top hat window function
def W2d_cth(x):
   return 2.*special.jn(1, x) / x

# FT of a 1d tophat
# normalized to unity at k=0 (ie real space integral is 1)
def W1d_th(x):
   return sinc(x/2.)
   
def Si(x):
   return special.sici(x)[0]

def Ci(x):
   return special.sici(x)[1]

def sinc(x):
   return special.sph_jn(0, x)[0][0]

# relevant for isotropic Fourier transform in 2d
def j0(x):
   return special.jn(0, x)

# Modified Bessel function of the first kind
def i0(x):
   return special.iv(0, x)


##################################################################################
# formatting numbers

# clean scientific notation for file names
# removes trailing decimal point if not needed
def intExpForm(input):
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


# same as intExpForm, except always leaves the decimal point
def floatExpForm(input):
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


