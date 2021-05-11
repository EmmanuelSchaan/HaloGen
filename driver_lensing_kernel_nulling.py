import universe
reload(universe)
from universe import *

import projection_kernel
reload(projection_kernel)
from projection_kernel import *

import p2d
reload(p2d)
from p2d import *


##################################################################################

u = UnivPlanck15()

##################################################################################
# Nulling low-z signal for single-source planes

# single source planes
z = np.array([5., 6., 1100])
w0 = WeightLensSingle(u, z_source=z[0], name="limz5lens")
w1 = WeightLensSingle(u, z_source=z[1], name="limz6lens")
w2 = WeightLensSingle(u, z_source=z[2], name="cmblens")

# linear combination to null the low-z signal
w_combined = WeightLensSingle(u, z_source=z[2], name="null combi")
chi = u.bg.comoving_distance(z)
alpha = (1./chi[2] - 1./chi[0]) / (1./chi[0] - 1./chi[1])
w_combined.f = lambda a:  w2.f(a) + alpha * w1.f(a) - (1.+alpha) * w0.f(a)
w_combined.plotW(ws=[w0, w1, w2, w_combined])


##################################################################################
# What about source bins with wider redshift distribution?


##################################################################################


p2d_lim1lens = P2dAuto(u, u, w0, nProc=3, name='halofit', save=True)
p2d_lim2lens = P2dAuto(u, u, w1, nProc=3, name='halofit', save=True)
p2d_cmblens = P2dAuto(u, u, w2, nProc=3, name='halofit', save=True)

p2d_lim1cmblens = P2dCross(u, u, w2, w0, nProc=3, name='halofit', save=True)
p2d_lim2cmblens = P2dCross(u, u, w2, w1, nProc=3, name='halofit', save=True)



