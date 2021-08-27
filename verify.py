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
parser.add_argument("--beta",help="initial plasma pressure in units of Bphi pressure",default=0)
parser.add_argument("--v0solve",action='store_true',help="assume zero inflowing velocity")
parser.add_argument("--nr",help="Athena domain size in r",default=4000)
parser.add_argument("--u0min",help="min bracket multiplyer for vr",default=0.1)
parser.add_argument("--u0max",help="max bracket multiplyer for vr",default=25.0)
parser.add_argument("--bmax",help="max bracket multiplyer to solve for Bphi",default=20)
parser.add_argument("--vshock",help="shock speed",default=1e7)
parser.add_argument("--no_athena",action='store_true',help="Do not run Athena - output self-similar solution and stop")

args = parser.parse_args()

nr=4*int(int(args.nr)/4.0)
chi=float(args.chi)
gam=float(args.gamma)
Bphi0=float(args.bphi)
Bz0=float(args.bz)
rho0=float(args.rho)
rot=float(args.rot)
beta=float(args.beta)
vshock=float(args.vshock)

beta_athena=beta/(8*np.pi*(1+(Bz0/Bphi0)**2))
ratio=Bz0/Bphi0
MA=np.sqrt(Bphi0**2+Bz0**2)/np.sqrt(4*np.pi*rho0)/vshock

print("trying to find the solution with initial conditions:")
print("rho    = %.5e*r^%.2f [g/cm^3] "%(rho0,chi*2))
print("B_phi  = %.5e*r^%.2f   [gauss]"%(Bphi0,chi))
print("B_z    = %.5e*r^%.2f   [gauss]"%(Bz0,chi))

rho_coeff=1e4
B_coeff=np.sqrt(rho_coeff)/vshock

rho0*=rho_coeff
Bphi0*=B_coeff #1e-5?
Bz0*=B_coeff

steps=300000
out=np.zeros((steps,8))# 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F -?? must be different
eta_max=1.000001e5
eta_min=1e-5

if not args.v0solve: 
    vA=np.sqrt(Bphi0**2+Bz0**2)/np.sqrt(4*np.pi*rho0)
    U0min=-vA*float(args.u0max)
    U0max=-vA*float(args.u0min)
    print("iteratively finding v0 between %.2e and %.2e to obtain BC v(r=0)=0..."%(U0max*vshock,U0min*vshock))
    U0=find_v0(U0min, U0max, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_max, eta_min, out)
    v0=U0*vshock
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
    
#file name base
sig="r%02dma%03dchi%03dgam%03drot%03d"%(int(ratio*10+0.5),int(MA*100+0.5),int(chi*100+0.5),int(gam*100+0.5),int(rot*100+0.5))
fout="solutions/"+sig+".txt"

print("output to file ",fout)
f1=open(fout,'w')
preamble='''\
###########################################################################
# Magnetized Noh self-similar solution for
# gamma= %.5f
#  initial conditions:
# density = %.5e*r^%.2f [g/cm^3] 
# vshock  = %.5e        [cm/s]
# vr      = %.5e        [cm/s]
# vphi    = %.5e        [cm/s]
# B_phi   = %.5e*r^%.2f [gauss]
# B_z     = %.5e*r^%.2f [gauss]
# P       = %.5e*r^%.2f [erg/cm^3] 
###########################################################################
# step radius     density      velocity     B_phi          Bz        Pressure Vphi\
'''%(gam,rho0/rho_coeff,chi*2,vshock,v0,-v0*rot,Bphi0/B_coeff,chi,Bz0/B_coeff,chi,beta*(Bphi0/B_coeff)**2/8/np.pi,chi*2)
print(preamble,file=f1)

#for vshock=1e7
#Bx10^5
#rho/10^4
#vr*10^7

t=0.3/vshock
tm=1.0/vshock
t_fac=(t/tm)**chi
t_fac2=t_fac**2
pm=1.0/B_coeff/B_coeff
rhom=1.0/rho_coeff
bm=np.sqrt(4*np.pi)/B_coeff
vs=vshock

# 0eta, 1S, 2H, 3H2, 4N, 5U, 6F
r=out[:,0]*vs*t
v=vs*out[:,0]*out[:,5]
vphi=vs*out[:,0]*out[:,6]
rho=rhom*t_fac2*out[:,4]
B_phi=bm*t_fac*out[:,2]
B_z=bm*t_fac*out[:,3]
Pr=pm*t_fac2*out[:,1]*out[:,0]**2*out[:,4]*(1-out[:,5])/gam

for i in range(0,steps,100): print("%6d %.5e  %.5e  %.5e  %.5e  %.5e  %.5e   %.5e"%(i,r[i],rho[i],v[i],B_phi[i],B_z[i],Pr[i],vphi[i]),file=f1)
    
f1.close()

if args.no_athena: quit()

print("generating Athena++ input file")
f1=open('athinput.template')
template=f1.read()
f1.close()
if -v0>20*vshock: print("ATTENTION: v0 is too large, Athena result will be bogus"); print("ATTENTION: edit athinput.template and increase x1max")

template=template.replace("problem_id=magnoh","problem_id="+sig)
template=template.replace("nx1        = 2000","nx1=%d"%nr)
template=template.replace("x1min=3e-5","x1min=%.3e"%(0.25/nr))
template=template.replace("nx1        = 500","nx1=%d"%(int(nr/4)))
template=template.replace("gamma=g1111","gamma=%.5f"%gam)
template=template.replace("alpha=a1111","alpha=%.5f"%(chi*2))
template=template.replace("beta=b1111",  "beta=%.5f"%chi)
template=template.replace("pcoeff=p1111","pcoeff=%.5e"%(beta_athena))
template=template.replace("d=d1111","d=%.5e"%(rho0/rho_coeff))
template=template.replace("vr=v1111","vr=%.5e"%(v0))
template=template.replace("vphi=0","vphi=%.5e"%(-v0*rot))
template=template.replace("bphi=b1111","bphi=%.5e"%(Bphi0/B_coeff))
template=template.replace("bz=b2222","bz=%.5e"%(Bz0/B_coeff))
finput='solutions/'+sig+'input.txt'
f1=open(finput,"w")
f1.write(template)
f1.close()
print("running Athena with %d domain points"%nr)
os.system("mpirun --oversubscribe -np 4 ./athena -i "+finput)

print("plotting self-similar vs Athena solution")

a=np.loadtxt(sig+'.block0.out2.00001.tab')
a[:,8]*=np.sqrt(4*np.pi)
a[:,9]*=np.sqrt(4*np.pi)

if rot!=0: 
    variables=['B_phi','rho','Pr','B_z','v','vphi']
    titles=["$B_\phi$","$\\rho$","$P$","$B_z$","$v_r$","$v_\phi$"]
    v2index=[8,2,3,9,4,5]
else:
    variables=['rho','v','B_phi','B_z','Pr']
    titles=["$\\rho$","$v$","$B_\phi$","$B_z$","$P$"]
    v2index=[2,4,8,9,3]


if rot==0: plt.figure(figsize=(24.,3.5))
else: plt.figure(figsize=(17,6.))

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['axes.titlesize']= 20
plt.rcParams['axes.labelsize']= 20
plt.rcParams['xtick.labelsize']= 15
plt.rcParams['ytick.labelsize']= 15

if rot!=0: 
    plot_layout=231
    nplots=6
else: 
    plot_layout=151
    nplots=5

mask1=r<0.6
r2=a[:,1]
mask2=r2<0.6
    
for i in range(nplots):
    ax=plt.subplot(plot_layout+i)
    var=globals()[variables[i]]
    plt.title(titles[i])#,fontsize=18)
    if rot==0: plt.xlabel('$r$')
    #if rot==0: plt.ylabel(titles[i],rotation='0')
    plt.xlim([0,0.6])
    #ax.ticklabel_format(axis='y',style='sci',useOffset=False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    var2=a[:,v2index[i]]
    if i==5: plt.ylim([0,0.5*np.max(var2[mask2])])
    if i==4 and rot!=0: plt.ylim([np.min(var2[mask2]),-0.1*np.min(var2[mask2])])
    plt.plot(r2[mask2],var2[mask2],'go')
    plt.plot(r[mask1],var[mask1],'r-')
    if rot!=0: plt.subplots_adjust(wspace=0.3,hspace=0.3)
    else: plt.subplots_adjust(wspace=0.3)

if rot!=0:
    plt.figtext(0.8, 0.4, '$r_b=$%1.1f,  $M_a=$%1.2f, $v_\phi=$%.2f$v_r$'%(ratio, MA, rot), ha='center', va='center',fontsize=20)
    plt.figtext(0.8, 0.3, '$\chi=$%1.1f,  $\gamma=$%1.2f'%(chi, gam), ha='center', va='center',fontsize=20)
else:
    plt.figtext(0.85, 0.4, '$r_b=$%1.1f,  $M_a=$%1.2f'%(ratio, MA), ha='center', va='center',fontsize=20)
    plt.figtext(0.85, 0.3, '$\chi=$%1.1f,  $\gamma=$%1.2f'%(chi, gam), ha='center', va='center',fontsize=20)

fpicture='solutions/'+sig+'.png'
fpicture2='solutions/'+sig+'.eps'
print("saving into "+fpicture)
plt.savefig(fpicture,dpi=100,bbox_inches="tight", pad_inches = 0)
plt.savefig(fpicture2,dpi=50,transparent=True,bbox_inches="tight",pad_inches = 0)

print("removing Athena output files")
os.system("rm "+sig+".hst")
os.system("rm "+sig+".*.tab")
print("Done")
