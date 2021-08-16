#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import os
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
mpl.use('Agg')
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
parser.add_argument("--u0min",help="min bracket multiplyer for vr",default=0.1)
parser.add_argument("--u0max",help="max bracket multiplyer for vr",default=25.0)
parser.add_argument("--v0solve",action='store_true',help="assume zero inflowing velocity")
parser.add_argument("--bmax",help="max bracket multiplyer to solve for Bphi",default=20)


args = parser.parse_args()

chi=float(args.chi)
gam=float(args.gamma)
Bphi0=float(args.bphi)
Bz0=float(args.bz)
rho0=float(args.rho)
rot=float(args.rot)
vshock=1.0
beta=0.0

ratio=Bz0/Bphi0
MA=np.sqrt(Bphi0**2+Bz0**2)/np.sqrt(4*np.pi*rho0)

print("trying to find the solution with initial conditions:")
print("rho    = %.5e*r^%.2f [g/cm^3] "%(rho0,chi*2))
print("B_phi  = %.5e*r^%.2f   [gauss]"%(Bphi0,chi))
print("B_z    = %.5e*r^%.2f   [gauss]"%(Bz0,chi))


steps=100000
out=np.zeros((steps,8))# 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F -?? must be different
eta_max=1.000001e5
eta_min=1e-5

if not args.v0solve: 
    vA=np.sqrt(Bphi0**2+Bz0**2)/np.sqrt(4*np.pi*rho0)
    U0min=-vA*float(args.u0max)
    U0max=-vA*float(args.u0min)
    print("iteratively finding v0 between %.2e and %.2e to obtain BC v(r=0)=0..."%(U0max,U0min))
    U0=find_v0(U0min, U0max, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_max, eta_min, out)
    v0=U0
    print('found solution for v0=%.5e [cm/s]'%v0)
    solve(U0, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_max, eta_min, out, True)
else: # finding Bphi for v(0)=0 case
    rot=0
    #rho0=1
    bmin=Bphi0
    bmax=Bphi0*float(args.bmax)
    print("iteratively finding b_phi between %.2e and %.2e to obtain BC v(r=0)=0..."%(bmin,bmax))
    bphi=find_bphi(bmin, bmax, ratio, rho0, chi, gam, 0, beta, steps, eta_max, eta_min, out)
    MA=np.sqrt(bphi**2+(bphi*ratio)**2)/np.sqrt(4*np.pi)
    print('found solution for bphi=%.5e'%bphi)
    solve(0, bphi, ratio, rho0, chi, gam, 0, beta, steps, eta_max, eta_min, out, True)
    Bphi0=bphi
    Bz0=bphi*ratio
    #rho0=1
    v0=0
    


#   out=np.zeros((steps,8)) # 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F

out[:,1]*=(1-out[:,5])*out[:,4]*(out[:,0]**2)/gam
out[:,5]*=-1

#   out=np.zeros((steps,8)) # 0eta, 1P, 2H, 3H2, 4N, 5-U, 6W, 7F

#out[:,6]=-(np.roll(out[:,3],1)-np.roll(out[:,3],-1))/(np.log(np.roll(out[:,0],1))-np.log(np.roll(out[:,0],-1)))/out[:,0]
#out[0,6]=np.nan
#out[-1,6]=np.nan

titles=["$\eta$","$P$","$H_\phi$","$H_z$","$N$","$-\eta U$","$j_\phi$"]


plt.figure(figsize=(17,8.))

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['axes.titlesize']= 20
plt.rcParams['axes.labelsize']= 20
plt.rcParams['xtick.labelsize']= 15
plt.rcParams['ytick.labelsize']= 15

plot_layout=231
nplots=5


eta=out[:,0]
eta_min=eta[-1]
eta_max=eta[0]
#print(eta[0],eta[-1])
mask1=eta<0.3 #max to plot
eta1=eta[mask1]
#print(eta1[0],eta1[-1])

sigma=2*chi*(2-gam)/(chi+2)
omega=2*chi/(chi+2)
epsilon=gam*out[-1,1]/out[-1,3]**2
print("sigma, epsilon:",sigma, epsilon)

out2=np.zeros((len(eta1),8))
out2[:,0]=eta1[:]

U0=-chi/2
U1=-chi*(chi+2)*(2-gam)/2/gam/(3*chi-gam*chi+2)
N1=-((4-gam)*chi+2)/gam/((3-gam)*chi+2)
P1=-((3-gam)*chi+1)/gam/((3-gam)*chi+2)
Hp1=-((3-gam)*chi+1)/gam/((3-gam)*chi+2)
Hz1=-1/gam

#   out2=np.zeros((steps,8)) # 0eta, 1P, 2H, 3H2, 4N, 5-U, 6W, 7F
#out2[:,5]=-eta1*(U0+epsilon*U1*eta1**sigma)
out2[:,5]=-(U0+epsilon*U1*eta1**sigma)

out2[:,4]=(1+epsilon*N1*eta1**sigma)*eta1**(2*chi/(chi+2))
out2[:,1]=(1+epsilon*P1*eta1**sigma)*eta1**sigma
out2[:,2]=(1+epsilon*Hp1*eta1**sigma)*eta1**(chi/(chi+2))
out2[:,3]=(1+epsilon*Hz1*eta1**sigma)
#normalization
out2[:,4]*=out[-1,4]/out2[-1,4]
out2[:,1]*=out[-1,1]/out2[-1,1]
out2[:,2]*=out[-1,2]/out2[-1,2]
out2[:,3]*=out[-1,3]/out2[-1,3]

mask3=eta>3 # max to plot
eta3=eta[mask3]
out3=np.zeros((len(eta3),8))
out3[:,0]=eta3
#   out=np.zeros((steps,8)) # 0eta, 1P, 2H, 3H2, 4N, 5-U, 6W, 7F
Uinf=-out[0,5]*out[0,0] # etaU at inf
Azinf=(out[0,3]**2)*(1+out[0,5])/out[0,4]
Apinf=(out[0,2]**2)*(1+out[0,5])/out[0,4]
print("Uinf,Azinf,Apinf:",Uinf,Azinf,Apinf)

out3[:,5]=Uinf-((chi+1)*Apinf+chi*Azinf)/eta3-((3*chi+2)*Apinf+(chi+1)*Azinf)*Uinf/2/eta3**2
out3[:,5]/=-eta3 #  make -U

out3[:,4]=(1-(2*chi+1)*Uinf/eta3+chi*((chi+1)*Apinf+chi*Azinf+(2*chi+1)*Uinf**2)/eta3**2)*eta3**(2*chi)
out3[:,2]=(1-chi*Uinf/eta3+(chi-1)*((chi+1)*Apinf+chi*Azinf+chi*Uinf**2)/2/eta3**2)*eta3**chi
out3[:,3]=(1-(chi+1)*Uinf/eta3+chi*((chi+1)*Apinf+chi*Azinf+(chi+1)*Uinf**2)/2/eta3**2)*eta3**chi
# normalization
out3[:,4]*=out[0,4]/out3[0,4]
out3[:,2]*=out[0,2]/out3[0,2]
out3[:,3]*=out[0,3]/out3[0,3]


# plot etaU
out[:,5]*=eta
out2[:,5]*=eta1
out3[:,5]*=eta3


for i in range(1,nplots+1):
    ax=plt.subplot(plot_layout+i-1)
    plt.title(titles[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(out[:,0],out[:,i],'g-')
    plt.plot(out2[:,0],out2[:,i],'r-')
    plt.plot(out3[:,0],out3[:,i],'b-')
    plt.xlim([1e-4,1e4])


plt.figtext(0.8, 0.35, '$r_b=$%1.1f,  $M_a=$%1.2f, $v_\phi=$%.2f$v_r$'%(ratio, MA, rot), ha='center', va='center',fontsize=20)
plt.figtext(0.8, 0.25, '$\chi=$%1.1f,  $\gamma=$%1.2f'%(chi, gam), ha='center', va='center',fontsize=20)

#plt.show()

fpicture2='asympt1.eps'
plt.savefig(fpicture2,dpi=50,transparent=True,bbox_inches="tight",pad_inches = 0)

#plt.savefig(fpicture,dpi=100,bbox_inches="tight", pad_inches = 0)
