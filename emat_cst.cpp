/*
 *  emat_cst.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 8/2/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */

#include <math.h>
#include <iostream>

#include "delfem/vector3d.h"
#include "emat_cst.h"

void MakePositiveDefinite_Sim22(const double s2[3],double s3[3])
{
  const double b = (s2[0]+s2[1])*0.5;
//  const double c = s2[0]*s2[1]-s2[2]*s2[2];
  const double d = (s2[0]-s2[1])*(s2[0]-s2[1])*0.25 + s2[2]*s2[2];
  const double e = sqrt(d);
  if( b-e > 1.0e-20 ){ 
    s3[0] = s2[0];
    s3[1] = s2[1];
    s3[2] = s2[2];
    return;
  }
  if( b+e < 0 ){
    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;    
    return;
  }
  const double l = b+e;
  double t0[2] = { s2[0]-l, s2[2]   };
  double t1[2] = { s2[2],   s2[1]-l };
  //  std::cout << t0[0]*t1[1]-t0[1]*t1[0] << std::endl;
  const double sqlen_t0 = t0[0]*t0[0]+t0[1]*t0[1];
  const double sqlen_t1 = t1[0]*t1[0]+t1[1]*t1[1]; 
  if( sqlen_t0 > sqlen_t1 ){
    if( sqlen_t0 < 1.0e-20 ){ 
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;    
      return;      
    }      
    const double invlen_t0 = 1.0/sqrt(sqlen_t0);
    t0[0] *= invlen_t0;
    t0[1] *= invlen_t0; 
    s3[0] = l*t0[0]*t0[0];
    s3[1] = l*t0[1]*t0[1];
    s3[2] = l*t0[0]*t0[1];    
  }
  else{    
    if( sqlen_t1 < 1.0e-20 ){ 
      s3[0] = 0;
      s3[1] = 0;
      s3[2] = 0;    
      return;      
    }      
    const double invlen_t1 = 1.0/sqrt(sqlen_t1);
    t1[0] *= invlen_t1;
    t1[1] *= invlen_t1; 
    s3[0] = l*t1[0]*t1[0];
    s3[1] = l*t1[1]*t1[1];
    s3[2] = l*t1[0]*t1[1];        
  }
  return;
}


void GetMatRes_CST_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3],
 double lambda, double myu,
 double dt, double damp_coeff,
 double& strain_energy)
{	
	double Kmat0[3][3][3][3];
	const double c[3][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
		{ C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
	GetKmatRes_CST(Kmat0,Res, C,c, lambda,myu, strain_energy);
	for(unsigned int i=0;i<3*3*3*3;i++){ (&Kmat[0][0][0][0])[i] = dt*(dt+damp_coeff)*(&Kmat0[0][0][0][0])[i]; }
	for(unsigned int ino=0;ino<3;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    Res[ino][idim] *= dt;
    for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      Res[ino][idim] += damp_coeff*Kmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt;
      Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt*dt;          
    }
    }
  }
	}  
}


void GetMatRes_CST_NewmarkBeta
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3], const double a[3][3],
 double lambda, double myu,
 double dt, double gamma_n, double beta_n, bool isinit)
{
	double Kmat0[3][3][3][3];
	const double c[3][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
		{ C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
  double se = 0;
	GetKmatRes_CST(Kmat0,Res, C,c, lambda,myu,se);
	for(unsigned int i=0;i<3*3*3*3;i++){ (&Kmat[0][0][0][0])[i] = dt*dt*beta_n*(&Kmat0[0][0][0][0])[i]; }
	if( isinit ){
		for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jno=0;jno<3;jno++){
      for(unsigned int jdim=0;jdim<3;jdim++){
        Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt;
        Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*a[jno][jdim]*dt*dt*0.5;
      }
      }
    }
		}
	}
}

void GetKmatResdRdC_CST
(double Kmat[3][3][3][3], double Res[3][3], double dRdC[3][3][3][3],
 const double C[3][3], const double c[3][3],
 double lambda, double myu)
{
  double Gd[3][3] = {
		{ C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
		{ C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  Com::UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
	
	double Gu[2][3];
	{
    double invArea = 0.5/Area;
    Com::Cross3D(Gu[0], Gd[1], Gd[2]);
		Gu[0][0] *= invArea;	Gu[0][1] *= invArea;	Gu[0][2] *= invArea;
		////
    Com::Cross3D(Gu[1], Gd[2], Gd[0]);
		Gu[1][0] *= invArea ;	Gu[1][1] *= invArea;	Gu[1][2] *= invArea;
	}
	
	const double gd[2][3] = { 
		{ c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
		{ c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
  const double E2[3] = {  // engineering strain
		0.5*( Com::Dot3D(gd[0],gd[0]) - Com::Dot3D(Gd[0],Gd[0]) ),
		0.5*( Com::Dot3D(gd[1],gd[1]) - Com::Dot3D(Gd[1],Gd[1]) ),
		1.0*( Com::Dot3D(gd[0],gd[1]) - Com::Dot3D(Gd[0],Gd[1]) ) };      
  
  const double E2dC[3][3][3] = { // E2dC[a][i][k] : derivative of E2[k] w.r.t. C[a][i]
    { {+Gd[0][0],+Gd[1][0],+Gd[0][0]+Gd[1][0]},
      {+Gd[0][1],+Gd[1][1],+Gd[0][1]+Gd[1][1]},
      {+Gd[0][2],+Gd[1][2],+Gd[0][2]+Gd[1][2]} },
    { {-Gd[0][0],0,        -Gd[1][0]},
      {-Gd[0][1],0,        -Gd[1][1]},
      {-Gd[0][2],0,        -Gd[1][2]} },
    { {0,        -Gd[1][0],-Gd[0][0]},
      {0,        -Gd[1][1],-Gd[0][1]},
      {0,        -Gd[1][2],-Gd[0][2]} } };
  
  double DAdC[3][3];
  {
    const double D[3][3] = {
      { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] },
      { C[0][0]-C[2][0], C[0][1]-C[2][1], C[0][2]-C[2][2] },
      { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] } };
    Com::Cross3D(DAdC[0], Gd[2], D[0]);
    Com::Cross3D(DAdC[1], Gd[2], D[1]);
    Com::Cross3D(DAdC[2], Gd[2], D[2]);
  }
  
  double GuDC[3][3][2][3];  // GuDC[a][i][b][j] : derivative of Gu[b][j] w.r.t. C[a][i]
  {
    const double invSqDArea = 0.25/(Area*Area);
    double G101[3]; Com::Cross3D(G101,  Gd[1], Gd[2]);
    double G010[3]; Com::Cross3D(G010,  Gd[2], Gd[0]);    
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){      
        GuDC[0][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[0][idim] - 2.0*Gd[1][idim]*Gd[0][jdim] + (Gd[0][idim]+Gd[1][idim])*Gd[1][jdim] );
        GuDC[1][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[1][idim] - Gd[1][idim]*Gd[1][jdim] );
        GuDC[2][idim][0][jdim] = invSqDArea*( -2.0*G101[jdim]*DAdC[2][idim] + 2.0*Gd[1][idim]*Gd[0][jdim] - Gd[0][idim]*Gd[1][jdim] );
        GuDC[0][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[0][idim] - 2.0*Gd[0][idim]*Gd[1][jdim] + (Gd[0][idim]+Gd[1][idim])*Gd[0][jdim] );
        GuDC[1][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[1][idim] + 2.0*Gd[0][idim]*Gd[1][jdim] - Gd[1][idim]*Gd[0][jdim] );
        GuDC[2][idim][1][jdim] = invSqDArea*( -2.0*G010[jdim]*DAdC[2][idim] - Gd[0][idim]*Gd[0][jdim] );
      }
      GuDC[0][idim][0][idim] += invSqDArea*( -Com::Dot3D(Gd[1],Gd[1]) + Com::Dot3D(Gd[0],Gd[1]) );
      GuDC[1][idim][0][idim] += invSqDArea*( +Com::Dot3D(Gd[1],Gd[1]) );
      GuDC[2][idim][0][idim] += invSqDArea*( -Com::Dot3D(Gd[1],Gd[0]) );
      GuDC[0][idim][1][idim] += invSqDArea*( -Com::Dot3D(Gd[0],Gd[0]) + Com::Dot3D(Gd[0],Gd[1]) );
      GuDC[1][idim][1][idim] += invSqDArea*( -Com::Dot3D(Gd[1],Gd[0]) );
      GuDC[2][idim][1][idim] += invSqDArea*( +Com::Dot3D(Gd[0],Gd[0]) );
    }
  }
  
  double GuGu2dC[3][3][3];  // GuDC[a][i][k] : derivative of GuGu[k] w.r.t. C[a][i]
  for(unsigned int jno=0;jno<3;jno++){
  for(unsigned int jdim=0;jdim<3;jdim++){
    GuGu2dC[jno][jdim][0] = 2.0*Com::Dot3D(GuDC[jno][jdim][0],Gu[0]);
    GuGu2dC[jno][jdim][1] = 2.0*Com::Dot3D(GuDC[jno][jdim][1],Gu[1]);
    GuGu2dC[jno][jdim][2] = Com::Dot3D(GuDC[jno][jdim][0],Gu[1]) + Com::Dot3D(GuDC[jno][jdim][1],Gu[0]);
  }
  }
  const double GuGu2[3] = { Com::Dot3D(Gu[0],Gu[0]), Com::Dot3D(Gu[1],Gu[1]), Com::Dot3D(Gu[1],Gu[0]) }; 
  
  const double Cons2[6] = {
    (lambda+2*myu)*GuGu2[0]*GuGu2[0],
    (lambda+2*myu)*GuGu2[1]*GuGu2[1],
    (lambda+  myu)*GuGu2[2]*GuGu2[2] + myu*GuGu2[0]*GuGu2[1],
    (lambda+2*myu)*GuGu2[1]*GuGu2[2],
    (lambda+2*myu)*GuGu2[0]*GuGu2[2],
    lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]) };   
  
  double Cons2dC[3][3][6];
  {
    for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      Cons2dC[jno][jdim][0] = 2*(lambda+2*myu)*GuGu2[0]*GuGu2dC[jno][jdim][0];
      Cons2dC[jno][jdim][1] = 2*(lambda+2*myu)*GuGu2[1]*GuGu2dC[jno][jdim][1];
      Cons2dC[jno][jdim][2] = 2*(lambda+  myu)*GuGu2[2]*GuGu2dC[jno][jdim][2] + myu*( GuGu2[1]*GuGu2dC[jno][jdim][0] + GuGu2[0]*GuGu2dC[jno][jdim][1] );
      Cons2dC[jno][jdim][3] =   (lambda+2*myu)*( GuGu2[1]*GuGu2dC[jno][jdim][2] + GuGu2[2]*GuGu2dC[jno][jdim][1] );
      Cons2dC[jno][jdim][4] =   (lambda+2*myu)*( GuGu2[0]*GuGu2dC[jno][jdim][2] + GuGu2[2]*GuGu2dC[jno][jdim][0] );
      Cons2dC[jno][jdim][5] = lambda*( GuGu2[1]*GuGu2dC[jno][jdim][0] + GuGu2[0]*GuGu2dC[jno][jdim][1] ) + 4*myu*GuGu2[2]*GuGu2dC[jno][jdim][2];
    }
    }    
  }

  const double S2[3] = {  // stress
    Cons2[0]*E2[0] + Cons2[5]*E2[1] + Cons2[4]*E2[2], 
    Cons2[5]*E2[0] + Cons2[1]*E2[1] + Cons2[3]*E2[2], 
    Cons2[4]*E2[0] + Cons2[3]*E2[1] + Cons2[2]*E2[2] };  
  
  const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };    
  for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      Res[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }  
  
  double S2dC[3][3][3];
  for(unsigned int jno=0;jno<3;jno++){
  for(unsigned int jdim=0;jdim<3;jdim++){
    S2dC[jno][jdim][0]
    = Cons2[0]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][0]*E2[0]
    + Cons2[5]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][5]*E2[1]
    + Cons2[4]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][4]*E2[2];
    S2dC[jno][jdim][1]
    = Cons2[5]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][5]*E2[0]
    + Cons2[1]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][1]*E2[1]
    + Cons2[3]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][3]*E2[2];    
    S2dC[jno][jdim][2]
    = Cons2[4]*E2dC[jno][jdim][0] + Cons2dC[jno][jdim][4]*E2[0]
    + Cons2[3]*E2dC[jno][jdim][1] + Cons2dC[jno][jdim][3]*E2[1]
    + Cons2[2]*E2dC[jno][jdim][2] + Cons2dC[jno][jdim][2]*E2[2];        
  }    
  }

  for(unsigned int ino=0;ino<3;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      dRdC[ino][jno][idim][jdim] = 
      +(Area*S2dC[jno][jdim][0]+0.5*DAdC[jno][jdim]*S2[0])*gd[0][idim]*dNdr[ino][0]
      +(Area*S2dC[jno][jdim][2]+0.5*DAdC[jno][jdim]*S2[2])*gd[0][idim]*dNdr[ino][1]
      +(Area*S2dC[jno][jdim][2]+0.5*DAdC[jno][jdim]*S2[2])*gd[1][idim]*dNdr[ino][0]
      +(Area*S2dC[jno][jdim][1]+0.5*DAdC[jno][jdim]*S2[1])*gd[1][idim]*dNdr[ino][1];
    }
    }
  }
  }    
  
  //////////////////////////
  // Calc stiffness matrix 
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
	for(unsigned int ino=0;ino<3;ino++){		
    for(unsigned int jno=0;jno<ino+1;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){	
          double dtmp0 = 0;				  
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[0]+dNdr[ino][1]*Cons2[4])*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[5]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[1][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[5])*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[3]+dNdr[ino][1]*Cons2[1])*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][0];
          Kmat[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
       +S3[2]*dNdr[ino][0]*dNdr[jno][1]
       +S3[2]*dNdr[ino][1]*dNdr[jno][0]
       +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
      Kmat[ino][jno][0][0] += dtmp1;
      Kmat[ino][jno][1][1] += dtmp1;
      Kmat[ino][jno][2][2] += dtmp1;
    }
	}
  
  for(unsigned int ino=0;ino<3;ino++){		
    for(unsigned int jno=ino+1;jno<3;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Kmat[ino][jno][idim][jdim] = Kmat[jno][ino][jdim][idim];
        }
      }
    }
  }                        
}

static double TriArea2D(const double v1[2], const double v2[2], const double v3[2]){
  double z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return z*0.5;
}


void GetKmatRes_CST
(double Kmat[3][3][3][3], double res[3][3],
 const double C[3][3], const double c[3][3], 
 const double lambda, const double myu,
 double& strain_energy)
{
	double Gd[3][3] = {
		{ C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
		{ C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  Com::UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
	
	double Gu[2][3];
	{
    const double invDArea = 0.5/Area;
    Com::Cross3D(Gu[0], Gd[1], Gd[2]);
		Gu[0][0] *= invDArea;	Gu[0][1] *= invDArea;	Gu[0][2] *= invDArea;
		////
    Com::Cross3D(Gu[1], Gd[2], Gd[0]);
		Gu[1][0] *= invDArea;	Gu[1][1] *= invDArea;	Gu[1][2] *= invDArea;
	}
	
	const double gd[2][3] = { 
		{ c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
		{ c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
	    	
	const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
        
  const double E2[3] = {  // engineering strain
		0.5*( Com::Dot3D(gd[0],gd[0]) - Com::Dot3D(Gd[0],Gd[0]) ),
		0.5*( Com::Dot3D(gd[1],gd[1]) - Com::Dot3D(Gd[1],Gd[1]) ),
		1.0*( Com::Dot3D(gd[0],gd[1]) - Com::Dot3D(Gd[0],Gd[1]) ) };    
  const double GuGu2[3] = { Com::Dot3D(Gu[0],Gu[0]), Com::Dot3D(Gu[1],Gu[1]), Com::Dot3D(Gu[1],Gu[0]) };
  const double Cons2[6] = {
    (lambda+2*myu)*GuGu2[0]*GuGu2[0],
    (lambda+2*myu)*GuGu2[1]*GuGu2[1],
    (lambda+  myu)*GuGu2[2]*GuGu2[2] + myu*GuGu2[0]*GuGu2[1],
    (lambda+2*myu)*GuGu2[1]*GuGu2[2],
    (lambda+2*myu)*GuGu2[0]*GuGu2[2],
    lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]) };
  const double S2[3] = {  // stress
    Cons2[0]*E2[0] + Cons2[5]*E2[1] + Cons2[4]*E2[2], 
    Cons2[5]*E2[0] + Cons2[1]*E2[1] + Cons2[3]*E2[2], 
    Cons2[4]*E2[0] + Cons2[3]*E2[1] + Cons2[2]*E2[2] };
  
  strain_energy += 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);

  ////
  for(unsigned int ino=0;ino<3;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    res[ino][idim] = Area*
    (+S2[0]*gd[0][idim]*dNdr[ino][0]
     +S2[2]*gd[0][idim]*dNdr[ino][1]
     +S2[2]*gd[1][idim]*dNdr[ino][0]
     +S2[1]*gd[1][idim]*dNdr[ino][1]);
	}
  }
  
  double S3[3] = { S2[0], S2[1], S2[2] };
  MakePositiveDefinite_Sim22(S2,S3);
	for(unsigned int ino=0;ino<3;ino++){		
  for(unsigned int jno=0;jno<ino+1;jno++){
    for(unsigned int idim=0;idim<3;idim++){
    for(unsigned int jdim=0;jdim<3;jdim++){	
      double dtmp0 = 0;				  
      dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[0]+dNdr[ino][1]*Cons2[4])*gd[0][jdim]*dNdr[jno][0];
      dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[5]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][1];
      dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[0][jdim]*dNdr[jno][1];
      dtmp0 += gd[0][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[2])*gd[1][jdim]*dNdr[jno][0];
      dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[4]+dNdr[ino][1]*Cons2[5])*gd[0][jdim]*dNdr[jno][0];
      dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[3]+dNdr[ino][1]*Cons2[1])*gd[1][jdim]*dNdr[jno][1];
      dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[0][jdim]*dNdr[jno][1];
      dtmp0 += gd[1][idim]*(dNdr[ino][0]*Cons2[2]+dNdr[ino][1]*Cons2[3])*gd[1][jdim]*dNdr[jno][0];
      Kmat[ino][jno][idim][jdim] = dtmp0*Area;
    }
    }
    const double dtmp1 = Area*
    (+S3[0]*dNdr[ino][0]*dNdr[jno][0]
     +S3[2]*dNdr[ino][0]*dNdr[jno][1]
     +S3[2]*dNdr[ino][1]*dNdr[jno][0]
     +S3[1]*dNdr[ino][1]*dNdr[jno][1]);
    Kmat[ino][jno][0][0] += dtmp1;
    Kmat[ino][jno][1][1] += dtmp1;
    Kmat[ino][jno][2][2] += dtmp1;
  }
	}
  
  for(unsigned int ino=0;ino<3;ino++){		
  for(unsigned int jno=ino+1;jno<3;jno++){
    for(unsigned int idim=0;idim<3;idim++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      Kmat[ino][jno][idim][jdim] = Kmat[jno][ino][jdim][idim];
    }
    }
  }
  }                      
}

/*
void GetKmatRes_CST
(double Kmat[3][3][3][3], double res[3][3],
 const double C[3][3], const double c[3][3], 
 const double lambda, const double myu,
 double& strain_energy)
{
	double Gd[3][3] = {
		{ C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
		{ C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
	UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
	
	double Gu[3][2];
	{
		Cross3D(Gu[0], Gd[1], Gd[2]);
		const double invtmp1 = 1.0/Dot3D(Gu[0],Gd[0]);
		Gu[0][0] *= invtmp1;	Gu[0][1] *= invtmp1;	Gu[0][2] *= invtmp1;
		////
		Cross3D(Gu[1], Gd[2], Gd[0]);
		const double invtmp2 = 1.0/Dot3D(Gu[1],Gd[1]);
		Gu[1][0] *= invtmp2;	Gu[1][1] *= invtmp2;	Gu[1][2] *= invtmp2;
	}
	
	const double gd[2][3] = { 
		{ c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] },
		{ c[2][0]-c[0][0], c[2][1]-c[0][1], c[2][2]-c[0][2] } };
  
	const double dNdr[3][2] = { {-1.0, -1.0}, {+1.0, +0.0}, {+0.0, +1.0} };
  
  const double E2[3] = {  // engineering strain
		0.5*( Dot3D(gd[0],gd[0]) - Dot3D(Gd[0],Gd[0]) ),
		0.5*( Dot3D(gd[1],gd[1]) - Dot3D(Gd[1],Gd[1]) ),
		1.0*( Dot3D(gd[0],gd[1]) - Dot3D(Gd[0],Gd[1]) ) };    
  const double GuGu2[3] = { Dot3D(Gu[0],Gu[0]), Dot3D(Gu[1],Gu[1]), Dot3D(Gu[1],Gu[0]) };
  const double Cons2[3][3] = {
    { lambda*GuGu2[0]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[0]),
      lambda*GuGu2[0]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[0]*GuGu2[2] + 2*myu*(GuGu2[0]*GuGu2[2]) },
    { lambda*GuGu2[1]*GuGu2[0] + 2*myu*(GuGu2[2]*GuGu2[2]),
      lambda*GuGu2[1]*GuGu2[1] + 2*myu*(GuGu2[1]*GuGu2[1]),
      lambda*GuGu2[1]*GuGu2[2] + 2*myu*(GuGu2[2]*GuGu2[1]) },
    { lambda*GuGu2[2]*GuGu2[0] + 2*myu*(GuGu2[0]*GuGu2[2]),
      lambda*GuGu2[2]*GuGu2[1] + 2*myu*(GuGu2[2]*GuGu2[1]),
      lambda*GuGu2[2]*GuGu2[2] + 1*myu*(GuGu2[0]*GuGu2[1] + GuGu2[2]*GuGu2[2]) } };
  const double S2[3] = {  // stress
    Cons2[0][0]*E2[0] + Cons2[0][1]*E2[1] + Cons2[0][2]*E2[2], 
    Cons2[1][0]*E2[0] + Cons2[1][1]*E2[1] + Cons2[1][2]*E2[2], 
    Cons2[2][0]*E2[0] + Cons2[2][1]*E2[1] + Cons2[2][2]*E2[2] };
  
  strain_energy += 0.5*Area*(E2[0]*S2[0] + E2[1]*S2[1] + E2[2]*S2[2]);
  
  for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      res[ino][idim] = Area*
      (+S2[0]*gd[0][idim]*dNdr[ino][0]
       +S2[2]*gd[0][idim]*dNdr[ino][1]
       +S2[2]*gd[1][idim]*dNdr[ino][0]
       +S2[1]*gd[1][idim]*dNdr[ino][1]);
    }
  }
	
	for(unsigned int ino=0;ino<3;ino++){		
    for(unsigned int jno=0;jno<ino+1;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){	
          double dtmp0 = 0;				  
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][0]*Cons2[0][2]*gd[1][jdim]*dNdr[jno][0];                                        
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][1]*Cons2[1][2]*gd[1][jdim]*dNdr[jno][0];                                        
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[0][idim]*dNdr[ino][1]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];                            
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][0]*gd[0][jdim]*dNdr[jno][0];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][1]*gd[1][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[0][jdim]*dNdr[jno][1];
          dtmp0 += gd[1][idim]*dNdr[ino][0]*Cons2[2][2]*gd[1][jdim]*dNdr[jno][0];                            
          Kmat[ino][jno][idim][jdim] = dtmp0*Area;
        }
      }
      const double dtmp1 = Area*
      (+S2[0]*dNdr[ino][0]*dNdr[jno][0]
       +S2[2]*dNdr[ino][0]*dNdr[jno][1]
       +S2[2]*dNdr[ino][1]*dNdr[jno][0]
       +S2[1]*dNdr[ino][1]*dNdr[jno][1]);
      Kmat[ino][jno][0][0] += dtmp1;
      Kmat[ino][jno][1][1] += dtmp1;
      Kmat[ino][jno][2][2] += dtmp1;
    }
	}
  for(unsigned int ino=0;ino<3;ino++){		
    for(unsigned int jno=ino+1;jno<3;jno++){
      for(unsigned int idim=0;idim<3;idim++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Kmat[ino][jno][idim][jdim] = Kmat[jno][ino][jdim][idim];
        }
      }
    }
  }                    
}

*/


void GetMatRes_MassCST_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3],
 const double C[3][3], const double u[3][3], const double v[3][3],
 const double lambda, const double myu,
 const double rho, double gx, double gy, double gz,
 const double dt, const double damp_coeff,
 double& kinetic_energy,
 double& strain_energy,
 double& potential_energy )
{
  double c[3][3] = {
    { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
    { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
    { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
  GetKmatRes_CST(Kmat,Res, C,c, lambda, myu,strain_energy);
	for(unsigned int i=0;i<3*3;    i++){ (&Res[0][0])[i] *= dt; }
  for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jno=0;jno<3;jno++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Res[ino][idim] += Kmat[ino][jno][idim][jdim]*v[jno][jdim]*dt*dt;
        }
      }
    }
  }
	for(unsigned int i=0;i<3*3*3*3;i++){ (&Kmat[0][0][0][0])[i] *= dt*(dt+damp_coeff); }
  ////
  const double Area = Com::TriArea3D(C[0], C[1], C[2]);
  const double tmp_r = Area/3.0*rho*dt;
  Res[0][0] -= gx*tmp_r;
  Res[0][1] -= gy*tmp_r;
  Res[0][2] -= gz*tmp_r;
  Res[1][0] -= gx*tmp_r;
  Res[1][1] -= gy*tmp_r;
  Res[1][2] -= gz*tmp_r;
  Res[2][0] -= gx*tmp_r;
  Res[2][1] -= gy*tmp_r;
  Res[2][2] -= gz*tmp_r;
  const double tmp_m = Area/3.0*rho;
  for(unsigned int ino=0;ino<3;ino++){
    Kmat[ino][ino][0][0] += tmp_m;
    Kmat[ino][ino][1][1] += tmp_m;
    Kmat[ino][ino][2][2] += tmp_m;
  }
  const double sqv0 = v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2];
  const double sqv1 = v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2];
  const double sqv2 = v[2][0]*v[2][0] + v[2][1]*v[2][1] + v[2][2]*v[2][2];  
  kinetic_energy += 0.5*tmp_m*(sqv0+sqv1+sqv2);
  potential_energy -= tmp_m*( gx*(u[0][0]+u[1][0]+u[2][0]) +  gy*(u[0][1]+u[1][1]+u[2][1]) +  gz*(u[0][2]+u[1][2]+u[2][2]) );
}





