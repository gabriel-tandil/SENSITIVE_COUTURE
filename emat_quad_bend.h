/*
 *  emat_quad_bend.h
 *  ds
 *
 *  Created by Nobuyuki Umetani on 11/30/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */

#if !defined(EMAT_QUAD_BEND_H)
#define EMAT_QUAD_BEND_H


void GetKmatRes_QuadBend
(double Kmat[4][4][3][3], double res[4][3],
 const double C[4][3], const double c[4][3],
 double stiff,
 double& strain_energy);


void GetMatRes_QuadBend_BackwardEular
(double Kmat[4][4][3][3], double Res[4][3], 
 const double C[4][3], 
 const double u[4][3], const double v[4][3],
 double stiff,
 double dt,
 double& strain_energy);


void GetKmatResdRdC_QuadBend
(double Kmat[4][4][3][3], double Res[4][3], double dRdC[4][4][3][3],
 const double C[4][3], const double c[4][3],
 double bend_stiff);


#endif