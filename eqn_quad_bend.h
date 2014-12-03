#if !defined(EQN_QUAD_BEND_H)
#define EQN_QUAD_BEND_H


/*
 *  eqn_quadric_bend.h
 *  ds
 *
 *  Created by Nobuyuki Umetani on 12/1/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */


#include "delfem/linearsystem_interface_eqnsys.h"


bool AddLinSys_QuadBend_CST_BackWardEular
(double dt, double damp_coeff,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_bend, double lambda, double myu,
 double rho, double g_x, double g_y, double g_z, 
 const unsigned int id_field_disp, const unsigned int id_field_hinge,
 const Fem::Field::CFieldWorld& world, 
 double& kinetic_energy,
 double& strain_energy,
 double& potential_energy);

bool AddLinSys_QuadBend_CST_Sensitivity_FictitousBending
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_bend, double lambda, double myu,
 double rho, double g_x, double g_y, double g_z,
 const unsigned int id_field_disp, const unsigned int id_field_hinge,
 const unsigned int id_field_senseX, const unsigned int id_field_senseY,
 const unsigned int id_field_lamX, const unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world);

#endif
