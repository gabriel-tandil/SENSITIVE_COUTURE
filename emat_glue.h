/*
 *  emat_glue.h
 *  ds
 *
 *  Created by Nobuyui Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */

#ifndef EMAT_GLUE_H
#define EMAT_GLUE_H

void GetKmatRes_GluePointLine_Penalty
(double Kmat[3][3][3][3], double res[3][3],
 const double C[3][3], const double c[3][3], double r,
 double stiff,
 double& se);

void GetMatRes_Glue_Penalty_NewmarkBeta
(double Kmat[2][2][3][3], double Res[2][3], 
 const double C[2][3], 
 const double u[2][3], const double v[2][3], const double a[2][3],
 double stiff, 
 double dt, double gamma_n, double beta_n, bool isinit);


void GetMatRes_GluePointLine_Penalty_BackwardEular
(double Kmat[3][3][3][3], double Res[3][3], 
 const double C[3][3], 
 const double u[3][3], const double v[3][3], double r,
 double stiff, 
 double dt,
 double damp_coeff,
 double& strain_energy);

void GetMatRes_Glue_Lagrange_NewmarkBeta
(double mat[2][2][3][3], double mat_ll[3][3],  
 double mat_ld[2][3][3], double mat_dl[2][3][3],
 double Res[2][3],  double Resl[3],
 double C[2][3], double u[2][3], double v[2][3], double a[2][3], 
 double ul[3], double vl[3], double al[3], 
 double dt, double gamma_n, double beta_n,
 bool isinit);



#endif