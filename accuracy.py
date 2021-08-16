#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import argparse
from magnoh_solve import solve,find_v0,find_bphi

parser=argparse.ArgumentParser()
parser.add_argument("chi",help="main exponent")
parser.add_argument("gamma",help="gas gamma")
parser.add_argument("bphi",help="toroidal b")
parser.add_argument("bz",help="b along axis")
parser.add_argument("rho",help="density")
parser.add_argument("--rot",help="rotation fraction",default=0)
parser.add_argument("--beta",help="initial plasma pressure in units of Bphi pressure",default=0)
parser.add_argument("--v0solve",action='store_true',help="assume zero inflowing velocity")
parser.add_argument("--u0min",help="min bracket multiplyer for vr",default=0.1)
parser.add_argument("--u0max",help="max bracket multiplyer for vr",default=25.0)
parser.add_argument("--bmax",help="max bracket multiplyer to solve for Bphi",default=20)
parser.add_argument("--no_plot",action='store_true',help="do not generate plot")


args = parser.parse_args()

chi=float(args.chi)
gam=float(args.gamma)
Bphi0=float(args.bphi)
Bz0=float(args.bz)
rho0=float(args.rho)
rot=float(args.rot)
beta=float(args.beta)
ratio=Bz0/Bphi0

args1=(300000, 1.000001e5,1e-5)
args2=(600000, 1.000001e5,1e-5)
out1=np.zeros((args1[0],8))# 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F -?? must be different
out2=np.zeros((args2[0],8))

if not args.v0solve: 
    vA=np.sqrt(Bphi0**2+Bz0**2)/np.sqrt(4*np.pi*rho0)
    U0min=-vA*float(args.u0max)
    U0max=-vA*float(args.u0min)
    print("iteratively finding v0 between %.2e and %.2e to obtain BC v(r=0)=0..."%(U0max,U0min))
    U0=find_v0(U0min, U0max, Bphi0, ratio, rho0, chi, gam, rot, beta, *args1, out1)
    print('found solution for v0=%.5e'%U0)
    solve(U0, Bphi0, ratio, rho0, chi, gam, rot, beta, *args1, out1, True)
    print("iteratively finding v0 between %.2e and %.2e to obtain BC v(r=0)=0..."%(U0max,U0min))
    U0=find_v0(U0min, U0max, Bphi0, ratio, rho0, chi, gam, rot, beta, *args2, out2)
    print('found solution for v0=%.5e'%U0)
    solve(U0, Bphi0, ratio, rho0, chi, gam, rot, beta, *args2, out2, True)
else: # finding Bphi for v(0)=0 case
    rot=0
    bmin=Bphi0
    bmax=Bphi0*float(args.bmax)
    print("iteratively finding b_phi between %.2e and %.2e to obtain BC v(r=0)=0..."%(bmin,bmax))
    bphi=find_bphi(bmin, bmax, ratio, rho0, chi, gam, 0, beta, *args1, out1)
    print('found solution for bphi=%.5e'%bphi)
    solve(0, bphi, ratio, rho0, chi, gam, 0, beta, *args1, out1, True)
    print("iteratively finding b_phi between %.2e and %.2e to obtain BC v(r=0)=0..."%(bmin,bmax))
    bphi=find_bphi(bmin, bmax, ratio, rho0, chi, gam, 0, beta, *args2, out2)
    print('found solution for bphi=%.5e'%bphi)
    solve(0, bphi, ratio, rho0, chi, gam, 0, beta, *args2, out2, True)


# 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F -?? must be different

eta1=out1[:,0]
eta2=out2[:,0]
out1[:,5]=np.abs(eta1*out1[:,5])
out2[:,5]=np.abs(eta2*out2[:,5])
out1[:,6]=np.abs(eta1*out1[:,6])
out2[:,6]=np.abs(eta2*out2[:,6])

#out1=out1[eta1>1e-4,:]
#out2=out2[eta2>1e-4,:]


vars=["$\eta$", "$S$", "$H_\phi$", "$H_z$", "$N$", "$-\eta U$", "$|\eta W|$"]
vars2=["eta ", "S   ", "Hphi ", "Hz  ", "N   ", "etaU", "etaW"]

if rot==0: max_iter=6
else: max_iter=7

if not args.no_plot:
    for i in range(1,max_iter):
        plt.figure(i-1)
        plt.title(vars[i])
        plt.yscale('log')
        plt.xscale('log')
        plt.plot(out2[:,0],out2[:,i],'g-')
        plt.plot(out1[:,0],out1[:,i],'r-')
    plt.show()

eta1=np.log(out1[::-1,0])
s1=out1[::-1,1]
out1=np.log(out1[::-1,2:max_iter+1])
mask1=eta1<0
s1=np.log(s1[mask1])


eta2=np.log(out2[::-1,0])
s2=out2[::-1,1]
out2=np.log(out2[::-1,2:max_iter+1])
mask2=eta2<0
s2=np.log(s2[mask2])

out3=np.zeros_like(out1)
s3=np.zeros_like(s1)

s3=np.interp(eta1[mask1],eta2[mask2],s2)

for i in range(4):
    out3[:,i]=np.interp(eta1,eta2,out2[:,i])


print( "mean error "+vars2[1]+":",np.nanmean(np.abs(s3-s1)))

for i in range(2,max_iter):
    print( "mean error "+vars2[i]+":",np.nanmean(np.abs(out3[:,i-2]-out1[:,i-2])))

