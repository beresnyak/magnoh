//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
/*============================================================================*/
/*! \file magnoh_cyl.cpp -- with many domains in phi
    \brief Magnetized Noh with perturbation in B_phi
    A. Beresnyak */
/*============================================================================*/

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"

static Real alpha,beta,rho0,P0,pcoeff,vr,perturb;
static Real r_perturb;
static double grids_per_cm;
static double x1min;
static Real bphi0,bz;
const int nzones=5;
static int mphi[nzones];
const int ngap=5;
int nx2_new, nx2_dom;
int perturb_on;
const double perturb_time=4e-9;



FILE* stab_output;
FILE* rad_output;
AthenaArray<double> cosk,sink;

Real DivergenceB(MeshBlock *pmb, int iout);

void ix2_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void ox2_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(1);
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // initialize global variables
   //nu_iso = pin->GetOrAddReal("problem","nu_iso",0.0);
   //eta_Ohm = pin->GetOrAddReal("problem","eta_Ohm",0.0);

   alpha  = pin->GetReal("problem", "alpha");
   beta  = pin->GetReal("problem", "beta");
   pcoeff= pin->GetReal("problem", "pcoeff");
   rho0  = pin->GetReal("problem", "d");
   vr =  pin->GetReal("problem", "vr");
   bphi0=pin->GetReal("problem", "bphi")/sqrt(4*M_PI); // convert from CGS to Athena Heaviside units
   bz=pin->GetReal("problem", "bz")/sqrt(4*M_PI);
   P0=4*M_PI*pcoeff*(bphi0*bphi0+bz*bz);
   perturb=pin->GetReal("problem","perturb");
   perturb_on=pin->GetInteger("problem","perturb_on"); // 1 if at t=0, 0 if at t=perturb_time
   r_perturb=pin->GetOrAddReal("problem","r_perturb",0.5*(1e7-vr)*perturb_time); // for t=0
   grids_per_cm=mesh_size.nx1/(mesh_size.x1max-mesh_size.x1min);  
   x1min=mesh_size.x1min;
   mphi[0]=0;char var_str[6]="mphi1";
 for(int dom=1; dom<nzones; ++dom)
   {
   mphi[dom]=pin->GetInteger("problem",var_str);
   var_str[4]+=1; //"mphi2" etc
   };
   
   if (Globals::my_rank==0 && NGHOST!=2) 
      { printf("NGHOST=%d !=2\n",NGHOST); throw std::runtime_error("something went wrong");}
   
   EnrollUserBoundaryFunction(INNER_X2, ix2_bc);
   EnrollUserBoundaryFunction(OUTER_X2, ox2_bc);

   AllocateUserHistoryOutput(1);
   EnrollUserHistoryOutput(0, DivergenceB, "divB");
   
if(Globals::my_rank==0){
 stab_output=fopen("etot.txt","w");
 fprintf(stab_output,"# total energy diff\n");
 fprintf(stab_output,"# time E_tot  ");
 for(int dom=1; dom<nzones; ++dom)
 fprintf(stab_output,"  mode=%d   ",mphi[dom]);
 fprintf(stab_output,"\n");        
 fflush(stab_output); 
 rad_output=fopen("erad.txt","w");
 fprintf(rad_output,"# r<r_sh energy diff\n");
 fprintf(rad_output,"# time E_r   ");
 for(int dom=1; dom<nzones; ++dom)
 fprintf(rad_output,"  mode=%d   ",mphi[dom]);
 fprintf(rad_output,"\n");        
 fflush(rad_output); 
}

  return;
}



//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for zpinch problem
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

if ((je-js+1-ngap*(nzones-1))%nzones!=0) 
    throw std::runtime_error("nx2-ngap*(nzones-1)) must be divisible by nzones");

nx2_new=(je-js+1-ngap*(nzones-1))/nzones;  //size of individual zone
nx2_dom=nx2_new+ngap;

//printf("%d %e %e %e\n",Globals::my_rank,pcoord->dx1f(is),pcoord->dx2f(js),pcoord->dx3f(ks));


if (MAGNETIC_FIELDS_ENABLED)
{
// initialize vector potential for inflowing B
  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx2,nx1);


cosk.NewAthenaArray(nzones,nx2_new);
sink.NewAthenaArray(nzones,nx2_new);
for(int dom=0; dom<nzones; ++dom)
for(int j=0; j<nx2_new; ++j)
{    
double phase=mphi[dom]*2.0*M_PI*j/nx2_new;
cosk(dom,j)=cos(phase);
sink(dom,j)=sin(phase);
}

if(Globals::my_rank==0)
 {
  fprintf(stderr,"zone size: %d\n",nx2_new);
  //fprintf(stderr,"%d %d\n",ks,ke);
 for(int dom=0; dom<nzones; ++dom)
  {
  int js_new=js+dom*(nx2_new+ngap);
  int je_new=js_new+nx2_new-1;
  if(fabs(-pcoord->x2f(js_new)/M_PI+pcoord->x2f(je_new+1)/M_PI-2.0)>1e-5){
    fprintf(stderr,"\n error: mismatched zone boundaries (in PI): ");
    fprintf(stderr,"%.6f %.6f\n",pcoord->x2f(js_new)/M_PI,pcoord->x2f(je_new+1)/M_PI);
#ifdef MPI_PARALLEL
    MPI_Abort(MPI_COMM_WORLD,1);
#else        
    exit(1);
#endif      
  }   
  }
  fprintf(stderr,"\n");
 } 
  
  for (int j=js; j<=je+1; ++j) {
  int ndom=(j-js)/nx2_dom; Real pert; 
  if (ndom<0 || ndom>nzones-1) ndom=0; 
  if (ndom==0) pert=0; else pert=perturb;
  int j_dom=(j-js)%nx2_dom;
  //if(j_dom>nx2_new) pert=0; 
  double phase=mphi[ndom]*2*M_PI*double(j_dom)/nx2_new;
  for (int i=is; i<=ie+1; ++i){
    if(perturb_on && pcoord->x1f(i)<3*r_perturb)
      az(j,i)=(bphi0/(beta+1))*pow(pcoord->x1f(i),beta+1)
                *(1+pert*exp(-2.*SQR(pcoord->x1f(i)/r_perturb))*cos(phase));       
    else   
      az(j,i)=(bphi0/(beta+1))*pow(pcoord->x1f(i),beta+1);
  }}

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j)/pcoord->x1f(i);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x3f(k,j,i) = bz*pow(pcoord->x1v(i),beta);
  }}}
  az.DeleteAthenaArray();
} //MAGNETIC_FIELDS_ENABLED

  // initialize conserved variables
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
             // Volume centered coordinates
             Real  rho=rho0*pow(pcoord->x1v(i),alpha);
             Real P  = P0*pow(pcoord->x1v(i),2*beta);
             phydro->u(IDN,k,j,i) = rho;
             phydro->u(IM1,k,j,i) = rho*vr;
             phydro->u(IM2,k,j,i) = 0.0;
             phydro->u(IM3,k,j,i) = 0.0;
             phydro->u(IEN,k,j,i) = P/(peos->GetGamma()-1.0) + 0.5*rho*SQR(vr);
          if (NON_BAROTROPIC_EOS && MAGNETIC_FIELDS_ENABLED) {
             phydro->u(IEN,k,j,i)+=
                0.5*0.25*(SQR(pfield->b.x1f(k,j,i)+pfield->b.x1f(k,j,i+1))
                        + SQR(pfield->b.x2f(k,j,i)+pfield->b.x2f(k,j+1,i))
                        + SQR(pfield->b.x3f(k,j,i)+pfield->b.x3f(k+1,j,i)));
          } // NON_BAROTROPIC_EOS && MAGNETIC_FIELDS_ENABLED
  }}}; //main loop

  return;
}

void project_array(int var, AthenaArray<Real> &arr, int is, int ie, int js, int je)
{
  //project to m mode    
for(int dom=0; dom<nzones; ++dom) //domain
 {
  int js_new=js+dom*(nx2_new+ngap);
  int je_new=js_new+nx2_new-1;     
  for (int i=is; i<=ie; ++i)
  {
    Real a0=0,a1=0,a2=0;
    for (int j=js_new; j<=je_new; ++j)
    {
      a0+=arr(var,0,j,i);
      a1+=arr(var,0,j,i)*cosk(dom,j-js_new);
      a2+=-arr(var,0,j,i)*sink(dom,j-js_new);    
    }    
    a0/=nx2_new; a1/=nx2_new; a2/=nx2_new;    
    for (int j=js_new; j<=je_new; ++j)
    {
      arr(var,0,j,i)=a0+(2*a1*cosk(dom,j-js_new)-2*a2*sink(dom,j-js_new))*(mphi[dom]!=0); 
    }
  }
 }
    
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
if(Globals::my_rank==0){fclose(stab_output);} return;
}

void MeshBlock::UserWorkInLoop(void)
{

if(Globals::my_rank!=0) return;
    

//introduce perturbation at t=perturb_time
if(!perturb_on && pmy_mesh->time>perturb_time)
{
  perturb_on=1;
  double shock_pos=pmy_mesh->time*1e7;
  double ipos=(shock_pos-x1min)*grids_per_cm;
  int i_shock=int(ipos);
  i_shock+=is;

  AthenaArray<Real> az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  az.NewAthenaArray(nx2,nx1);

  for (int j=js; j<=je+1; ++j) {
  int ndom=(j-js)/nx2_dom; Real pert; 
  if (ndom<0 || ndom>nzones-1) ndom=0; 
  if (ndom==0) pert=0; else pert=perturb;
  int j_dom=(j-js)%nx2_dom;
  //if(j_dom>nx2_new) pert=0; 
  double phase=mphi[ndom]*2*M_PI*double(j_dom)/nx2_new;
  double az0=0;
  for (int i=is; i<=ie+1; ++i){
    az(j,i)=az0*pert*((pcoord->x1f(i)>shock_pos)?0:(shock_pos-pcoord->x1f(i)))*cos(phase)/shock_pos;
    az0+=-pfield->b.x2f(ks,js,i)*pcoord->dx1f(i);       
  }}

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfield->b.x1f(k,j,i)+= (az(j+1,i) - az(j,i))/pcoord->dx2f(j)/pcoord->x1f(i);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfield->b.x2f(k,j,i)+= (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
  }}}

  az.DeleteAthenaArray();

//printf("%d %e\n",Globals::my_rank,DivergenceB(this,0));  
}  //introduce perturbation   


for(int var=0; var<=4; ++var)
  project_array(var, phydro->u, is, ie, js, je); 

project_array(0, pfield->b.x1f, is, ie, js, je);
project_array(0, pfield->b.x2f, is, ie, js, je);
project_array(0, pfield->b.x3f, is, ie, js, je);


if(pmy_mesh->ncycle%20==0 && pmy_mesh->time*1e7>x1min+1.0/grids_per_cm)
 {
  double shock_pos=pmy_mesh->time*1e7;
  double ipos=(shock_pos-x1min)*grids_per_cm;
  int i_shock=int(ipos);
  double i_d=ipos-i_shock;
  i_shock+=is;
  fprintf(stab_output,"%.4e ",pmy_mesh->time);       
  fprintf(rad_output, "%.4e ",pmy_mesh->time);       
  for(int dom=0; dom<nzones; ++dom) //domain
  {
  int jump=dom*(nx2_new+ngap);    
  double diff=0; 
  double diff2=0;
  double en=0; 
  double en2=0; 

  for (int k=ks; k<=ke; k++) 
    for (int j=js; j<js+nx2_new; j++)
      for (int i=is; i<=ie; i++)
      {
  double E=0.5*pcoord->x1f(i)*SQR(pfield->b.x1f(k,j,i))
          +0.5*pcoord->x1v(i)*(
                 SQR(pfield->b.x2f(k,j,i))
                +SQR(pfield->b.x3f(k,j,i))
                              )
          +0.5*pcoord->x1v(i)*(
                 SQR(phydro->u(IM1,k,j,i))
                +SQR(phydro->u(IM2,k,j,i))       
                +SQR(phydro->u(IM3,k,j,i))       
                              )/phydro->u(IDN,k,j,i);
               
  double Ed=0.5*pcoord->x1f(i)*SQR(pfield->b.x1f(k,j,i)-pfield->b.x1f(k,j+jump,i))
           +0.5*pcoord->x1v(i)*SQR(pfield->b.x2f(k,j,i)-pfield->b.x2f(k,j+jump,i))
           +0.5*pcoord->x1v(i)*SQR(pfield->b.x3f(k,j,i)-pfield->b.x3f(k,j+jump,i))    
           +0.5*pcoord->x1v(i)*(
               SQR(phydro->u(IM1,k,j,i)-phydro->u(IM1,k,j+jump,i))
              +SQR(phydro->u(IM2,k,j,i)-phydro->u(IM2,k,j+jump,i))       
              +SQR(phydro->u(IM3,k,j,i)-phydro->u(IM3,k,j+jump,i))       
               )/phydro->u(IDN,k,j,i);
      en+=E;
      diff+=Ed;
      if(i<i_shock)
        {
          en2+=E;
          diff2+=Ed;
        }      
      if(i==i_shock)
        {
          en2+=E*i_d;
          diff2+=Ed*i_d;              
        }  
      } 
   if(dom==0) fprintf(stab_output,"%.4e ",en);
         else fprintf(stab_output,"%.4e ",diff);       
   if(dom==0) fprintf(rad_output,"%.4e ",en2);
         else fprintf(rad_output,"%.4e ",diff2);       
  }//dom
  fprintf(stab_output,"\n");
  fprintf(rad_output, "\n");
// }

 }//if
 
 
} //end UserWorkInLoop

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
    {
 //User output, div B
  Real dx1=pcoord->dx1f(is);
  Real dx2=pcoord->dx2f(js);
 for(int k=ks; k<=ke; k++)
    for(int j=js; j<=je; j++)
        for(int i=is; i<=ie; i++)
       {
         user_out_var(0,k,j,i) =
          ((pcoord->x1f(i+1)*pfield->b.x1f(k,j,i+1)-pcoord->x1f(i)*pfield->b.x1f(k,j,i))/dx1
         + (                 pfield->b.x2f(k,j+1,i)-               pfield->b.x2f(k,j,i))/dx2
          )/pcoord->x1v(i);
       }
        
    }


void PeriodicInnerX2(AthenaArray<Real> &prim, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
//if (Globals::my_rank==0 ) printf("call InnerX2 J=%d,%d\n",js,je);
  for (int n=0; n<(NHYDRO); ++n)
    for (int k=ks; k<=ke; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie; ++i)
       prim(n,k,js-j,i) = prim(n,k,je-(j-1),i);

  if (MAGNETIC_FIELDS_ENABLED){
    for (int k=ks; k<=ke; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie+1; ++i)
        b.x1f(k,js-j,i) = b.x1f(k,je-(j-1),i);

    for (int k=ks; k<=ke; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie; ++i)
        b.x2f(k,js-j,i) = b.x2f(k,je-(j-1),i); 

    for (int k=ks; k<=ke+1; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie; ++i)
        b.x3f(k,js-j,i) = b.x3f(k,je-(j-1),i);
  } //MAGNETIC_FIELDS_ENABLED

  return;
}

void PeriodicOuterX2(AthenaArray<Real> &prim, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke, int ngh) {
//if (Globals::my_rank==0 ) printf("call OuterX2 J=%d,%d\n",js,je);
  for (int n=0; n<(NHYDRO); ++n)
    for (int k=ks; k<=ke; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie; ++i)
       prim(n,k,je+j,i) = prim(n,k,js+(j-1),i);

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie+1; ++i)
        b.x1f(k,je+j,i) =  b.x1f(k,js+(j-1),i);

    for (int k=ks; k<=ke; ++k)
    for (int j=2; j<=ngh+1; ++j)
    for (int i=is; i<=ie; ++i)
        b.x2f(k,je+j,i) = b.x2f(k,js+(j-1),i);

    for (int k=ks; k<=ke+1; ++k)
    for (int j=1; j<=ngh; ++j)
    for (int i=is; i<=ie; ++i)
        b.x3f(k,je+j,i) = b.x3f(k,js+(j-1),i);
  } //MAGNETIC_FIELDS_ENABLED

  return;
}


void ix2_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{       
 for(int j=0; j<nzones; ++j)
  {
  int js_new=js+j*(nx2_new+ngap);
  int je_new=js_new+nx2_new-1;     
  PeriodicInnerX2(prim, b, is, ie, js_new, je_new, ks, ke, ngh);
  }
}

void ox2_bc(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{                    
 for(int j=0; j<nzones; ++j)
  {
  int js_new=js+j*(nx2_new+ngap);
  int je_new=js_new+nx2_new-1;     
  PeriodicOuterX2(prim, b, is, ie, js_new, je_new, ks, ke, ngh);
  }
}


Real DivergenceB(MeshBlock *pmb, int iout)
{
  if(Globals::my_rank!=0) return 0.0;

  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  FaceField &b = pmb->pfield->b;
  Coordinates* pco=pmb->pcoord;

  Real dx1=pco->dx1f(is);
  Real dx2=pco->dx2f(js);
  Real dx3=pco->dx3f(ks);

  Real measure_sq=0; int ncells=0;
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
     if(((j-js)%nx2_dom)<nx2_new) 
      for(int i=is+1; i<=ie; i++) {
        Real divb=((pco->x1f(i+1)*b.x1f(k,j,i+1)-pco->x1f(i)*b.x1f(k,j,i))/dx1
                  +(              b.x2f(k,j+1,i)-            b.x2f(k,j,i))/dx2)/pco->x1v(i);
        measure_sq+=4*SQR(divb*dx1)/(SQR(b.x2f(k,j,i))+SQR(b.x2f(k,j+1,i))
                                    +SQR(b.x1f(k,j,i))+SQR(b.x1f(k,j,i+1)));    
        ncells++;
      }
    }
  }

  return sqrt(measure_sq/ncells);
}
