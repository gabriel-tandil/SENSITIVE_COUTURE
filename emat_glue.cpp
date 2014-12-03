/*
 *  emat_glue.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <math.h>

#include "delfem/vector3d.h"
#include "emat_glue.h"

void GetMatRes_Glue_Lagrange_NewmarkBeta
(double mat[2][2][3][3], double mat_ll[3][3],  
 double mat_ld[2][3][3], double mat_dl[2][3][3],
 double Res[2][3],  double Resl[3],
 double C[2][3], double u[2][3], double v[2][3], double a[2][3], 
 double ul[3], double vl[3], double al[3], 
 double dt, double gamma_n, double beta_n,
 bool isinit)
{
	const double c[2][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] } };
	double eps = -0.0e-0;
	double Kmat_ll0[3][3] = { {eps,0,0}, {0,eps,0}, {0,0,eps} };
	const double k = 1.0e+0;
	double Kmat_dl0[2][3][3] = { { {+k,0,0}, {0,+k,0}, {0,0,+k} }, { {-k,0,0}, {0,-k,0}, {0,0,-k} } };
	double Kmat_ld0[2][3][3] = { { {+k,0,0}, {0,+k,0}, {0,0,+k} }, { {-k,0,0}, {0,-k,0}, {0,0,-k} } };
	const double d = 0.0;
	double Cmat_dl0[2][3][3] = { { {+d,0,0}, {0,+d,0}, {0,0,+d} }, { {-d,0,0}, {0,-d,0}, {0,0,-d} } };
	double Cmat_ld0[2][3][3] = { { {+d,0,0}, {0,+d,0}, {0,0,+d} }, { {-d,0,0}, {0,-d,0}, {0,0,-d} } };
	const double m = 0.0e-0;
	double Mmat_dl0[2][3][3] = { { {+m,0,0}, {0,+m,0}, {0,0,+m} }, { {-m,0,0}, {0,-m,0}, {0,0,-m} } };
	double Mmat_ld0[2][3][3] = { { {+m,0,0}, {0,+m,0}, {0,0,+m} }, { {-m,0,0}, {0,-m,0}, {0,0,-m} } };
	/*	Res[0][0] = ( + ul[0]*k + vl[0]*d + al[0]*m);
	 Res[0][1] = ( + ul[1]*k + vl[1]*d + al[1]*m);
	 Res[0][2] = ( + ul[2]*k + vl[2]*d + al[2]*m);
	 Res[1][0] = -Res[0][0];
	 Res[1][1] = -Res[0][1];
	 Res[1][2] = -Res[0][2];
	 Resl[0] = + ( (c[0][0]-c[1][0])*k + (v[0][0]-v[1][0])*d + (a[0][0]-a[1][0])*m );
	 Resl[1] = + ( (c[0][1]-c[1][1])*k + (v[0][1]-v[1][1])*d + (a[0][1]-a[1][1])*m );
	 Resl[2] = + ( (c[0][2]-c[1][2])*k + (v[0][2]-v[1][2])*d + (a[0][2]-a[1][2])*m );*/
	Res[0][0] = + ul[0]*k;
	Res[0][1] = + ul[1]*k;
	Res[0][2] = + ul[2]*k;
	Res[1][0] = - ul[0]*k;
	Res[1][1] = - ul[1]*k;
	Res[1][2] = - ul[2]*k;
	Resl[0] = (c[0][0]-c[1][0])*k;
	Resl[1] = (c[0][1]-c[1][1])*k;
	Resl[2] = (c[0][2]-c[1][2])*k;
	/*	{
	 double len = (c[0][0]-c[1][0])*(c[0][0]-c[1][0])+(c[0][1]-c[1][1])*(c[0][1]-c[1][1])+(c[0][2]-c[1][2])*(c[0][2]-c[1][2]);
	 std::cout << "Len : " << len << "   ---   " << ul[0] << " " << ul[1] << " " << ul[2] << std::endl;
	 std::cout << "      " << Resl[0] << " " << Resl[1] << " " << Resl[2] << std::endl;
	 }*/
	for(unsigned int i=0;i<2*2*3*3;i++){ (&mat[0][0][0][0])[i] = 0; }
	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_ld[0][0][0])[i] = dt*dt*beta_n*(&Kmat_ld0[0][0][0])[i] + dt*gamma_n*(&Cmat_ld0[0][0][0])[i] + (&Mmat_ld0[0][0][0])[i]; }
	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_dl[0][0][0])[i] = dt*dt*beta_n*(&Kmat_dl0[0][0][0])[i] + dt*gamma_n*(&Cmat_dl0[0][0][0])[i] + (&Mmat_dl0[0][0][0])[i]; }
	for(unsigned int i=0;i<3*3;    i++){ (&mat_ll[0][0]   )[i] = dt*dt*beta_n*(&Kmat_ll0[0][0]   )[i]; }
	//	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_ld[0][0][0])[i] = dt*dt*beta_n*(&Kmat_ld0[0][0][0])[i]; }
	//	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_dl[0][0][0])[i] = dt*dt*beta_n*(&Kmat_dl0[0][0][0])[i]; }
	//	for(unsigned int i=0;i<3*3;    i++){ (&mat_ll[0][0]   )[i] = dt*dt*beta_n*(&Kmat_ll0[0][0]   )[i]; }	
	if( isinit ){
		for(unsigned int ino=0;ino<2;ino++){
			for(unsigned int idim=0;idim<3;idim++){
				for(unsigned int jdim=0;jdim<3;jdim++){
					Res[ino][idim] += Kmat_dl0[ino][idim][jdim]*( vl[jdim]*dt+al[jdim]*dt*dt*0.5 )
					+ Cmat_dl0[ino][idim][jdim]*( al[jdim]*dt );
				}
			}
		}
		for(unsigned int idim=0;idim<3;idim++){
			for(unsigned int jno=0;jno<2;jno++){
				for(unsigned int jdim=0;jdim<3;jdim++){
					Resl[idim] += Kmat_ld0[jno][idim][jdim]*( v[jno][jdim]*dt+a[jno][jdim]*dt*dt*0.5 )
					+ Cmat_ld0[jno][idim][jdim]*( a[jno][jdim]*dt );
				}
			}
			for(unsigned int jdim=0;jdim<3;jdim++){
				Resl[idim] += Kmat_ll0[idim][jdim]*( vl[jdim]*dt+al[jdim]*dt*dt*0.5 );
			}
		}
	}
}

void GetKmatRes_GluePointPoint_Penalty
(double Kmat[2][2][3][3], double res[2][3],
 const double C[2][3], const double u[2][3],
 double stiff)
{
	const double c[2][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] } };
	const double vtmp0[3] = { c[1][0]-c[0][0], c[1][1]-c[0][1], c[1][2]-c[0][2] };
	res[0][0] = -stiff*vtmp0[0];
	res[0][1] = -stiff*vtmp0[1];
	res[0][2] = -stiff*vtmp0[2];
	res[1][0] = -res[0][0];
 	res[1][1] = -res[0][1];
	res[1][2] = -res[0][2];
	Kmat[0][0][0][0] = Kmat[1][1][0][0] = stiff;
	Kmat[0][0][0][1] = Kmat[1][1][0][1] = 0;
	Kmat[0][0][0][2] = Kmat[1][1][0][2] = 0;
	Kmat[0][0][1][0] = Kmat[1][1][1][0] = 0;
	Kmat[0][0][1][1] = Kmat[1][1][1][1] = stiff;
	Kmat[0][0][1][2] = Kmat[1][1][1][2] = 0;
	Kmat[0][0][2][0] = Kmat[1][1][2][0] = 0;
	Kmat[0][0][2][1] = Kmat[1][1][2][1] = 0;
	Kmat[0][0][2][2] = Kmat[1][1][2][2] = stiff;
	Kmat[0][1][0][0] = Kmat[1][0][0][0] = -Kmat[0][0][0][0];
	Kmat[0][1][0][1] = Kmat[1][0][0][1] = -Kmat[0][0][0][1];
	Kmat[0][1][0][2] = Kmat[1][0][0][2] = -Kmat[0][0][0][2];
	Kmat[0][1][1][0] = Kmat[1][0][1][0] = -Kmat[0][0][1][0];
	Kmat[0][1][1][1] = Kmat[1][0][1][1] = -Kmat[0][0][1][1];
	Kmat[0][1][1][2] = Kmat[1][0][1][2] = -Kmat[0][0][1][2];
	Kmat[0][1][2][0] = Kmat[1][0][2][0] = -Kmat[0][0][2][0];
	Kmat[0][1][2][1] = Kmat[1][0][2][1] = -Kmat[0][0][2][1];
	Kmat[0][1][2][2] = Kmat[1][0][2][2] = -Kmat[0][0][2][2];
}


void GetMatRes_Glue_Penalty_NewmarkBeta
(double Kmat[2][2][3][3], double Res[2][3], 
 const double C[2][3], 
 const double u[2][3], const double v[2][3], const double a[2][3],
 double stiff, 
 double dt, double gamma_n, double beta_n, bool isinit)
{
	double Kmat0[2][2][3][3];
	GetKmatRes_GluePointPoint_Penalty(Kmat0,Res, C,u, stiff);
	////////////////
  double damp1 = 0.0, damp2 = 0;
	double Cmat0[2][2][3][3];
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Cmat0[0][0][0][0])[i] = damp1*(&Kmat0[0][0][0][0])[i]; }
	for(unsigned int ino=0;ino<2;ino++){
		Cmat0[ino][ino][0][0] += damp2;
		Cmat0[ino][ino][1][1] += damp2;
		Cmat0[ino][ino][2][2] += damp2;
	}
	for(unsigned int ino=0;ino<2;ino++){
		for(unsigned int idim=0;idim<3;idim++){
			for(unsigned int jno=0;jno<2;jno++){
				for(unsigned int jdim=0;jdim<3;jdim++){
					Res[ino][idim] += Cmat0[ino][jno][idim][jdim]*v[jno][jdim];
				}
			}
		}
	}
	////////////////
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Kmat[0][0][0][0])[i] = dt*dt*beta_n*(&Kmat0[0][0][0][0])[i]; }
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Kmat[0][0][0][0])[i] +=  dt*gamma_n*(&Cmat0[0][0][0][0])[i]; }
	////////////////
	if( isinit ){
		for(unsigned int ino=0;ino<2;ino++){
			for(unsigned int idim=0;idim<3;idim++){
				for(unsigned int jno=0;jno<2;jno++){
					for(unsigned int jdim=0;jdim<3;jdim++){
						Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt;
						Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*a[jno][jdim]*dt*dt*0.5;
						Res[ino][idim] += Cmat0[ino][jno][idim][jdim]*a[jno][jdim]*dt;
					}
				}
			}
		}
	}
}



void GetMatRes_GluePointPoint_Penalty_BackwardEular
(double Kmat[2][2][3][3], double Res[2][3], 
 const double C[2][3], 
 const double u[2][3], const double v[2][3],
 double stiff, 
 double dt, bool isinit)
{
	double Kmat0[2][2][3][3];
	GetKmatRes_GluePointPoint_Penalty(Kmat0,Res, C,u, stiff);
	////////////////
  double damp1 = 0.0001*stiff;
	double Cmat0[2][2][3][3];
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Cmat0[0][0][0][0])[i] = damp1*(&Kmat0[0][0][0][0])[i]; }
	for(unsigned int ino=0;ino<2;ino++){
		for(unsigned int idim=0;idim<3;idim++){
			Res[ino][idim] *= dt;
			for(unsigned int jno=0;jno<2;jno++){
				for(unsigned int jdim=0;jdim<3;jdim++){
          Res[ino][idim] += Cmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt;
				}
			}
		}
	}
	////////////////
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Kmat[0][0][0][0])[i] = dt*dt*(&Kmat0[0][0][0][0])[i]; }
	for(unsigned int i=0;i<2*2*3*3;i++){ (&Kmat[0][0][0][0])[i] +=  dt*(&Cmat0[0][0][0][0])[i]; }
	////////////////
	if( isinit ){
		for(unsigned int ino=0;ino<2;ino++){
			for(unsigned int idim=0;idim<3;idim++){
				for(unsigned int jno=0;jno<2;jno++){
					for(unsigned int jdim=0;jdim<3;jdim++){
						Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*v[jno][jdim]*dt*dt;
						//						Res[ino][idim] += Kmat0[ino][jno][idim][jdim]*a[jno][jdim]*dt*dt*0.5;
						//						Res[ino][idim] += Cmat0[ino][jno][idim][jdim]*a[jno][jdim]*dt;
					}
				}
			}
		}
	}
} 

void GetKmatRes_GluePointLine_Penalty
(double Kmat[3][3][3][3], double res[3][3],
 const double C[3][3], const double c[3][3], double r,
 double stiff,
 double& se)
{
  double m[3] = {
    c[1][0]*(1-r)+c[2][0]*r,
    c[1][1]*(1-r)+c[2][1]*r,
    c[1][2]*(1-r)+c[2][2]*r };
	const double vtmp0[3] = { m[0]-c[0][0], m[1]-c[0][1], m[2]-c[0][2] };
  se += 0.5*stiff*(vtmp0[0]*vtmp0[0] + vtmp0[1]*vtmp0[1] + vtmp0[2]*vtmp0[2]);
	res[0][0] = -stiff*vtmp0[0];
	res[0][1] = -stiff*vtmp0[1];
	res[0][2] = -stiff*vtmp0[2];
	res[1][0] = -res[0][0]*(1-r);
 	res[1][1] = -res[0][1]*(1-r);
	res[1][2] = -res[0][2]*(1-r);
	res[2][0] = -res[0][0]*r;
 	res[2][1] = -res[0][1]*r;
	res[2][2] = -res[0][2]*r;
  for(unsigned int i=0;i<81;i++){ (&Kmat[0][0][0][0])[i] = 0; }
  Kmat[0][0][0][0] = stiff;
  Kmat[0][0][1][1] = stiff;
  Kmat[0][0][2][2] = stiff;      
  Kmat[1][1][0][0] = stiff*(1-r)*(1-r);
  Kmat[1][1][1][1] = stiff*(1-r)*(1-r);
  Kmat[1][1][2][2] = stiff*(1-r)*(1-r);
  Kmat[2][2][0][0] = stiff*r*r;
  Kmat[2][2][1][1] = stiff*r*r;
  Kmat[2][2][2][2] = stiff*r*r; 
  Kmat[0][1][0][0] = Kmat[1][0][0][0] = -stiff*(1-r);
	Kmat[0][1][1][1] = Kmat[1][0][1][1] = -stiff*(1-r);
	Kmat[0][1][2][2] = Kmat[1][0][2][2] = -stiff*(1-r);
  Kmat[0][2][0][0] = Kmat[2][0][0][0] = -stiff*r;
	Kmat[0][2][1][1] = Kmat[2][0][1][1] = -stiff*r;
	Kmat[0][2][2][2] = Kmat[2][0][2][2] = -stiff*r;
  Kmat[1][2][0][0] = Kmat[2][1][0][0] = +stiff*r*(1-r);
	Kmat[1][2][1][1] = Kmat[2][1][1][1] = +stiff*r*(1-r);
	Kmat[1][2][2][2] = Kmat[2][1][2][2] = +stiff*r*(1-r); 
}

void GetMatRes_GluePointLine_Penalty_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3], double r,
 double stiff, 
 double dt,
 double damp_coeff,
 double& strain_energy)
{
  const double c[3][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },    
		{ C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };  
	double Kmat0[3][3][3][3];
	GetKmatRes_GluePointLine_Penalty(Kmat0,Res, C,c,r, stiff, strain_energy);
	for(unsigned int i=0;i<3*3*3*3;i++){ (&Kmat[0][0][0][0])[i] = dt*(dt+damp_coeff)*(&Kmat0[0][0][0][0])[i]; }  
	for(unsigned int ino=0;ino<3;ino++){
    for(unsigned int idim=0;idim<3;idim++){
      Res[ino][idim] *= dt;
      for(unsigned int jno=0;jno<3;jno++){
        for(unsigned int jdim=0;jdim<3;jdim++){
          Res[ino][idim] += dt*(damp_coeff+dt)*Kmat0[ino][jno][idim][jdim]*v[jno][jdim];
        }
      }
    }
	}
}

