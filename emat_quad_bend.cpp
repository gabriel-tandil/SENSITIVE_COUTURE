/*
 *  emat_quad_bend.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 11/30/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <iostream>
#include <math.h>

#include "delfem/vector3d.h"
#include "emat_quad_bend.h"

void GetKmatRes_QuadBend
(double Kmat[4][4][3][3], double res[4][3],
 const double C[4][3], const double c[4][3],
 double stiff,
 double& strain_energy)
{
  const double A0 = Com::TriArea3D(C[0],C[2],C[3]);
  const double A1 = Com::TriArea3D(C[1],C[3],C[2]);
  const double L0 = Com::Distance3D(C[2],C[3]);
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };  
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };  
  double cot023, cot032;
  {
    const double r2 = -Com::Dot3D(e02,e23);
    const double r3 = +Com::Dot3D(e03,e23);    
    cot023 = r2/H0;
    cot032 = r3/H0;    
  }
  double cot123, cot132;
  {
    const double r2 = -Com::Dot3D(e12,e23);
    const double r3 = +Com::Dot3D(e13,e23);
    cot123 = r2/H1;
    cot132 = r3/H1;    
  }  
  const double tmp0 = stiff/((A0+A1)*L0*L0);  
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  for(unsigned int i=0;i<4*4*3*3;i++){ (&Kmat[0][0][0][0])[i] = 0; }
  for(unsigned int ino=0;ino<4;ino++){
  for(unsigned int jno=0;jno<4;jno++){
    const double tmp = K[ino]*K[jno]*tmp0;
    Kmat[ino][jno][0][0] = tmp;
    Kmat[ino][jno][1][1] = tmp;
    Kmat[ino][jno][2][2] = tmp;
  }    
  }
  for(unsigned int ino=0;ino<4;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    res[ino][idim] = 0;
    for(unsigned int jno=0;jno<4;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      res[ino][idim] += Kmat[ino][jno][idim][jdim]*c[jno][jdim];  
    }
    }
    strain_energy += res[ino][idim]*c[ino][idim];
  }
  }  
}


void GetMatRes_QuadBend_BackwardEular
(double Kmat[4][4][3][3], double Res[4][3], 
 const double C[4][3], 
 const double u[4][3], const double v[4][3],
 double stiff,
 double dt,
 double& strain_energy)
{
	const double c[4][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
		{ C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] },
		{ C[3][0]+u[3][0], C[3][1]+u[3][1], C[3][2]+u[3][2] }, };
	GetKmatRes_QuadBend(Kmat,Res, C,c, stiff, strain_energy);
	for(unsigned int i=0;i<4*3;    i++){ (&Res[0][0])[i] *= dt; }
	////    
  for(unsigned int ino=0;ino<4;ino++){
  for(unsigned int idim=0;idim<3;idim++){
    for(unsigned int jno=0;jno<4;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      Res[ino][idim] += Kmat[ino][jno][idim][jdim]*v[jno][jdim]*dt*dt;
    }
    }
  }
	}
	for(unsigned int i=0;i<4*4*3*3;i++){ (&Kmat[0][0][0][0])[i] *= dt*dt; }   
}

static void DerDoubleAreaTri3D
(double dAdC[3][3],
 const double c0[3], const double c1[3], const double c2[3],
 const double un[3])	// unit normal
{
  const double c[3][3] = {
    { c2[0]-c1[0], c2[1]-c1[1], c2[2]-c1[2] },
    { c0[0]-c2[0], c0[1]-c2[1], c0[2]-c2[2] },
    { c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2] } };
  Com::Cross3D(dAdC[0], un, c[0]);
  Com::Cross3D(dAdC[1], un, c[1]);
  Com::Cross3D(dAdC[2], un, c[2]);
}

inline void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
	n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
	const double invlen = 0.5/a;
	n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

void GetKmatResdRdC_QuadBend
(double Kmat[4][4][3][3], double Res[4][3], double dRdC[4][4][3][3],
 const double C[4][3], const double c[4][3],
 double stiff)
{
  double A0,UN0[3];  UnitNormalAreaTri3D(UN0, A0, C[0], C[2], C[3]);
  double A1,UN1[3];  UnitNormalAreaTri3D(UN1, A1, C[1], C[3], C[2]);
  const double L0 = Com::Distance3D(C[2],C[3]);
  double coeff  = stiff/((A0+A1)*L0*L0);
  double dcoeffdC[4][3];
  {
    double dA0dC023[3][3]; DerDoubleAreaTri3D(dA0dC023,C[0],C[2],C[3], UN0);  // 023    
    double dA1dC132[3][3]; DerDoubleAreaTri3D(dA1dC132,C[1],C[3],C[2], UN1);  // 132      
    const double t0 = stiff/((A0+A1)*(A0+A1)*L0*L0*L0*L0);
    const double u23[3] =  { (C[3][0]-C[2][0])/L0, (C[3][1]-C[2][1])/L0, (C[3][2]-C[2][2])/L0 };
    for(unsigned int i=0;i<3;i++){    
      dcoeffdC[0][i] = -dA0dC023[0][i]*L0*L0*0.5*t0;
      dcoeffdC[1][i] = -dA1dC132[0][i]*L0*L0*0.5*t0;
      dcoeffdC[2][i] = -((dA0dC023[1][i]+dA1dC132[2][i])*L0*L0*0.5 - (A0+A1)*2.0*L0*u23[i])*t0;
      dcoeffdC[3][i] = -((dA0dC023[2][i]+dA1dC132[1][i])*L0*L0*0.5 + (A0+A1)*2.0*L0*u23[i])*t0;
    }    
  }  /////////
  const double H0 = A0*2.0/L0;
  const double H1 = A1*2.0/L0;  
  ////
  const double e23[3] = { C[3][0]-C[2][0], C[3][1]-C[2][1], C[3][2]-C[2][2] };
  const double e02[3] = { C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] };
  const double e03[3] = { C[3][0]-C[0][0], C[3][1]-C[0][1], C[3][2]-C[0][2] };
  const double e12[3] = { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] };  
  const double e13[3] = { C[3][0]-C[1][0], C[3][1]-C[1][1], C[3][2]-C[1][2] };    
  double dKdC[4][3][4];
  double cot023, cot032, cot123, cot132;
  {
    const double d02 = -Com::Dot3D(e02,e23);
    const double d03 = +Com::Dot3D(e03,e23);    
    cot023 = d02/H0;        
    cot032 = d03/H0;    
    const double r02 = d03/(d02+d03);
    const double invH0 = 1.0/H0;
    double uvH0[3] =  { 
      ( r02*C[2][0]+(1-r02)*C[3][0] - C[0][0] )*invH0,
      ( r02*C[2][1]+(1-r02)*C[3][1] - C[0][1] )*invH0,
      ( r02*C[2][2]+(1-r02)*C[3][2] - C[0][2] )*invH0 };
//    std::cout << uvH0[0]*uvH0[0] + uvH0[1]*uvH0[1] + uvH0[2]*uvH0[2] << std::endl;
    double dcotdC023[3][3] = {
      { +e23[0]        *invH0 +         d02*uvH0[0]*invH0*invH0,
        +e23[1]        *invH0 +         d02*uvH0[1]*invH0*invH0,
        +e23[2]        *invH0 +         d02*uvH0[2]*invH0*invH0 },
      { (e02[0]-e23[0])*invH0 - r02    *d02*uvH0[0]*invH0*invH0,
        (e02[1]-e23[1])*invH0 - r02    *d02*uvH0[1]*invH0*invH0,
        (e02[2]-e23[2])*invH0 - r02    *d02*uvH0[2]*invH0*invH0 }, 
      { -e02[0]        *invH0 - (1-r02)*d02*uvH0[0]*invH0*invH0,
        -e02[1]        *invH0 - (1-r02)*d02*uvH0[1]*invH0*invH0,
        -e02[2]        *invH0 - (1-r02)*d02*uvH0[2]*invH0*invH0 } 
    };
    double dcotdC032[3][3] = {
      { -e23[0]        *invH0 +         d03*uvH0[0]*invH0*invH0,
        -e23[1]        *invH0 +         d03*uvH0[1]*invH0*invH0,
        -e23[2]        *invH0 +         d03*uvH0[2]*invH0*invH0 },
      { (e23[0]+e03[0])*invH0 - (1-r02)*d03*uvH0[0]*invH0*invH0,
        (e23[1]+e03[1])*invH0 - (1-r02)*d03*uvH0[1]*invH0*invH0,
        (e23[2]+e03[2])*invH0 - (1-r02)*d03*uvH0[2]*invH0*invH0 }, 
      { -e03[0]        *invH0 - r02    *d03*uvH0[0]*invH0*invH0,
        -e03[1]        *invH0 - r02    *d03*uvH0[1]*invH0*invH0,
        -e03[2]        *invH0 - r02    *d03*uvH0[2]*invH0*invH0 } 
    };    
    const double d12 = -Com::Dot3D(e12,e23);
    const double d13 = +Com::Dot3D(e13,e23);
    cot123 = d12/H1;
    cot132 = d13/H1; 
    const double r12 = d13/(d12+d13);
    const double invH1 = 1.0/H1;
    double uvH1[3] =  { 
      ( r12*C[2][0]+(1-r12)*C[3][0] - C[1][0] )*invH1,
      ( r12*C[2][1]+(1-r12)*C[3][1] - C[1][1] )*invH1,
      ( r12*C[2][2]+(1-r12)*C[3][2] - C[1][2] )*invH1 };
//    std::cout << uvH0[0]*uvH0[0] + uvH0[1]*uvH0[1] + uvH0[2]*uvH0[2] << std::endl;
    double dcotdC123[3][3] = {
      { +e23[0]        *invH1 +         d12*uvH1[0]*invH1*invH1,       
        +e23[1]        *invH1 +         d12*uvH1[1]*invH1*invH1,
        +e23[2]        *invH1 +         d12*uvH1[2]*invH1*invH1 },
      { (e12[0]-e23[0])*invH1 - r12    *d12*uvH1[0]*invH1*invH1,
        (e12[1]-e23[1])*invH1 - r12    *d12*uvH1[1]*invH1*invH1,
        (e12[2]-e23[2])*invH1 - r12    *d12*uvH1[2]*invH1*invH1 }, 
      { -e12[0]        *invH1 - (1-r12)*d12*uvH1[0]*invH1*invH1,
        -e12[1]        *invH1 - (1-r12)*d12*uvH1[1]*invH1*invH1,
        -e12[2]        *invH1 - (1-r12)*d12*uvH1[2]*invH1*invH1 } 
    };
    double dcotdC132[3][3] = {
      { -e23[0]        *invH1 +         d13*uvH1[0]*invH1*invH1,
        -e23[1]        *invH1 +         d13*uvH1[1]*invH1*invH1,
        -e23[2]        *invH1 +         d13*uvH1[2]*invH1*invH1 },
      { (e23[0]+e13[0])*invH1 - (1-r12)*d13*uvH1[0]*invH1*invH1,
        (e23[1]+e13[1])*invH1 - (1-r12)*d13*uvH1[1]*invH1*invH1,
        (e23[2]+e13[2])*invH1 - (1-r12)*d13*uvH1[2]*invH1*invH1 }, 
      { -e13[0]        *invH1 - r12    *d13*uvH1[0]*invH1*invH1,
        -e13[1]        *invH1 - r12    *d13*uvH1[1]*invH1*invH1,
        -e13[2]        *invH1 - r12    *d13*uvH1[2]*invH1*invH1 } 
    };            
    for(unsigned int i=0;i<3;i++){
      dKdC[0][i][0] = -dcotdC023[0][i]-dcotdC032[0][i];
      dKdC[1][i][0] = 0;
      dKdC[2][i][0] = -dcotdC023[1][i]-dcotdC032[2][i];
      dKdC[3][i][0] = -dcotdC023[2][i]-dcotdC032[1][i];
      ///
      dKdC[0][i][1] = 0;
      dKdC[1][i][1] = -dcotdC123[0][i]-dcotdC132[0][i];
      dKdC[2][i][1] = -dcotdC123[1][i]-dcotdC132[2][i];
      dKdC[3][i][1] = -dcotdC123[2][i]-dcotdC132[1][i];
      ///
      dKdC[0][i][2] = +dcotdC032[0][i];
      dKdC[1][i][2] = +dcotdC132[0][i];
      dKdC[2][i][2] = +dcotdC032[2][i]+dcotdC132[2][i];
      dKdC[3][i][2] = +dcotdC032[1][i]+dcotdC132[1][i];      
      ///
      dKdC[0][i][3] = +dcotdC023[0][i];
      dKdC[1][i][3] = +dcotdC123[0][i];
      dKdC[2][i][3] = +dcotdC023[1][i]+dcotdC123[1][i];
      dKdC[3][i][3] = +dcotdC023[2][i]+dcotdC123[2][i]; 
    }
  }
  const double K[4] = { -cot023-cot032, -cot123-cot132, cot032+cot132, cot023+cot123 };
  for(unsigned int i=0;i<4*4*3*3;i++){ (&Kmat[0][0][0][0])[i] = 0; }
  for(unsigned int ino=0;ino<4;ino++){
    for(unsigned int jno=0;jno<4;jno++){
      const double tmp = K[ino]*K[jno]*coeff;
      Kmat[ino][jno][0][0] = tmp;
      Kmat[ino][jno][1][1] = tmp;
      Kmat[ino][jno][2][2] = tmp;
    }    
  }
  for(unsigned int ino=0;ino<4;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      Res[ino][idim] = 0;
      for(unsigned int jno=0;jno<4;jno++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Res[ino][idim] += Kmat[ino][jno][idim][jdim]*c[jno][jdim];  
        }
      }
    }
  }  
  {    
    const double Kc[3] = {
      K[0]*c[0][0] + K[1]*c[1][0] + K[2]*c[2][0] + K[3]*c[3][0],
      K[0]*c[0][1] + K[1]*c[1][1] + K[2]*c[2][1] + K[3]*c[3][1],
      K[0]*c[0][2] + K[1]*c[1][2] + K[2]*c[2][2] + K[3]*c[3][2] };
    for(unsigned int ino=0;ino<4;ino++){
    for(unsigned int jno=0;jno<4;jno++){      
      for(unsigned int idim=0;idim<3;idim++){
      for(unsigned int jdim=0;jdim<3;jdim++){        
        dRdC[ino][jno][idim][jdim] = 
          dcoeffdC[jno][jdim]*K[ino]*Kc[idim] 
        + coeff*dKdC[jno][jdim][ino]*Kc[idim]
        + coeff*K[ino]*(dKdC[jno][jdim][0]*c[0][idim]+
                        dKdC[jno][jdim][1]*c[1][idim]+
                        dKdC[jno][jdim][2]*c[2][idim]+
                        dKdC[jno][jdim][3]*c[3][idim]);
        
      }    
      }
    }  
    }
  }
  
}