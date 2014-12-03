
/*
 *  emat_cst.h
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 8/2/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */

#if !defined(EMAT_CST_H)
#define EMAT_CST_H

void GetMatRes_CST_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3],
 double lambda, double myu,
 double dt, double damp_coeff,
 double& strain_energy);

void GetMatRes_MassCST_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3],
 const double C[3][3], const double u[3][3], const double v[3][3],
 const double lambda, const double myu,
 const double rho, double gx, double gy, double gz,
 const double dt, const double damp_coeff,
 double& kinetic_energy,
 double& strain_energy,
 double& potential_energy );

void GetMatRes_CST_NewmarkBeta
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3], const double a[3][3],
 double lambda, double myu,
 double dt, double gamma_n, double beta_n, bool isinit);

void GetKmatRes_CST
(double Kmat[3][3][3][3], double res[3][3],
 const double C[3][3], const double c[3][3], 
 const double lambda, const double myu,
 double& strain_energy);

void GetKmatResdRdC_CST
(double Kmat[3][3][3][3], double Res[3][3], double dRdC[3][3][3][3],
 const double C[3][3], const double c[3][3],
 double lambda, double myu);

void GetMatRes_CST_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3],
 double lambda, double myu,
 double dt );

#endif
