#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from magnoh_solve import find_v0

steps=300000
eta_st=1.000001e5
eta_end=1e-5

chi=1.5
gam=1.6666667
rho0=0.6

ratio=1.0
Bphi0=1.0
Bz0=1.0
MA=np.sqrt(2)/np.sqrt(4*np.pi*rho0)

U0min=-100
U0max=-0.1

out=np.zeros((1,1))

for gam in np.arange(1.1,1.55,0.05):
    v0=find_v0(U0min, U0max, Bphi0, ratio, rho0, chi, gam, 0, 0, steps, eta_st, eta_end, out)
    print("%.3f %.3f %.3f %.5f %.7e"%(ratio,MA,chi,gam,-v0))
