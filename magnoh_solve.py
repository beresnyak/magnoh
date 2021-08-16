#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from numba import jit,b1,i8,f8

#remove "from numba import" and @jit if you don't have numba installed, expect much slower code

@jit(f8(f8,f8,f8,f8,f8,f8,f8,f8,i8,f8,f8,f8[:,::1],b1),nopython=True)
def solve(v0, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out, record):
    gam2=(gam-1)/gam
    if gam>2: Umin=-chi/gam
    else: Umin=-chi/2
    dleta=(np.log(eta_end)-np.log(eta_st))/steps
    U=v0/eta_st
    if rot==0: W=1e-5*np.abs(U)
    else: W=rot*np.abs(U)
    #log initials
    #(eta,S,A,A2,N,W)=np.exp((leta,lS,lA,lA2,lN,lW)) 
    a=np.zeros(6,dtype=np.float64)
    la=np.zeros(6,dtype=np.float64)
    dla=np.zeros(6,dtype=np.float64)
    la[0]=np.log(eta_st)
    la[2]=np.log(Bphi0**2/(4*np.pi*rho0)/(1-U))-2*la[0]
    if beta==0: la[1]=-1000.0
    else: la[1]=la[2]+np.log(gam*beta*0.5) #S=0.5*gamma beta Aphi
    Bz0=ratio*Bphi0
    la[3]=  np.log(Bz0**2/(4*np.pi*rho0)/(1-U))-2*la[0]
    la[4]=np.log(rho0)+2*chi*la[0]
    la[5]=np.log(W)
    dla[0]=dleta

#   out=np.zeros((steps,8)) # 0eta, 1S, 2H, 3H2, 4N, 5U, 6W, 7F
#   integrate ODE
    for i in np.arange(steps): 
        if(la[0]>0 and la[0]+dla[0]<0): 
# begin shock jump condition
            if beta==0:
                MAsq=(1-U)/(a[2]+a[3])
                mu=(gam2+1.0/MAsq+np.sqrt((gam2+1.0/MAsq)**2+4*(1-2*gam2)*(2-gam2)/MAsq))/2.0/(2-gam2)
                MSsq_inv=gam*(1-mu-(1/mu**2-1)/2/MAsq)
                la[1]=np.log(MSsq_inv*(1-U))
            else:
                MAsq_inv=(a[2]+a[3])/(1-U)
                MSsq_inv=(a[1])/(1-U)
                minus_b=gam-1+2*MSsq_inv+gam*MAsq_inv
                aaa=gam+1
                ccc=(gam-2)*MAsq_inv
                mu=(minus_b+np.sqrt(minus_b**2-4*aaa*ccc))/(2*aaa)
                la[1]+=np.log(1+(1-0.5*(1+mu)*MAsq_inv/mu/mu)*gam*(1-mu)/MSsq_inv)                    
            la[2]-=2*np.log(mu)
            la[3]-=2*np.log(mu)
            la[4]-=np.log(mu)
            U=1-mu*(1-U)
# end shock jump condition
        a=np.exp(la) 
        F=U+a[1]+a[2]+a[3]-1.0

        if(la[0]<0.0 and F<0.0): return -10.0
        
        if record:    
            out[i,0]=a[0]
            out[i,2]=(a[2]*a[0]**2*a[4]*(1-U))**0.5
            out[i,3]=(a[3]*a[0]**2*a[4]*(1-U))**0.5
            out[i,4]=a[4]
            out[i,1]=a[1]
            out[i,5]=U
            out[i,6]=a[5]
            out[i,7]=F# M_f=(F/(1-U)+1)**(-0.5)

        if rot==0: WW=0
        else: WW=a[5]*a[5]
        dU=dla[0]*(-(chi+1.0)*a[2]-(chi+2.0*U)*a[3]-2.0*a[1]*(U+chi/gam)-U*(U-1.0)+WW)/F
        U1inv=1.0/(1.0-U)
        dU2=2.0*dU*U1inv
        dla[1]=(gam*dU+dla[0]*2.0*(gam*U-1.0))*U1inv
        dla[2]=dU2-2.0*dla[0]
        dla[3]=dU2+(4.0*U-2.0)*dla[0]*U1inv
        dla[4]=(dU+dla[0]*(2.0*U+2.0*chi))*U1inv
        dla[5]=(2*U-1)*dla[0]*U1inv
        U+=dU
        la+=dla
        
        if(U>0.0): return 1.0
        
    if record:
        if rot==0: out[:,6]=0

    return U-Umin

def solvev0(Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out):
    return solve(0, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out, False)

import scipy.optimize as opt

def find_v0(v0min, v0max, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out):
    return opt.bisect(solve,v0min,v0max,args=(Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out, False))

def find_bphi(bphi_min, bphi_max, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out):
    return opt.bisect(solvev0,bphi_min,bphi_max,args=(ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end, out))    



#    
#    amin=solve(v0min, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end)
#    amax=solve(v0max, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end)
#    if amin*amax>0: raise ValueError('Error in finding v0: brackets result in same sign, increase the interval [vmin,vmax]')
#    while True:
#        middle=(v0min+v0max)/2
#        amiddle=solve(middle, Bphi0, ratio, rho0, chi, gam, rot, beta, steps, eta_st, eta_end)
#        if(amiddle==0 or np.abs(v0min-v0max)<1e-6+1e-10*np.abs(v0max)): return middle
#        if(amiddle*amin>0):
#            v0min=middle
#            amin=amiddle
#        else:
#            v0max=middle
#            amax=amiddle
                
             

