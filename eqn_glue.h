#ifndef EQN_GLUE_H
#define EQN_GLUE_H

/*
 *  eqn_glue.h
 *  ds
 *
 *  Created by Nobuyuki Umetani on 7/27/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */
#include "delfem/femls/linearsystem_field.h"

bool AddLinSys_Glue_Lagrange_NewmarkBeta
(double dt, double gamma_newmark, double beta_newmark,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 unsigned int id_field_disp, unsigned int id_field_lambda,
 const Fem::Field::CFieldWorld& world,
 bool is_initial);

bool AddLinSys_Glue_Lagrange_BackwardEular
(double dt, 
 Fem::Eqn::ILinearSystem_Eqn& ls,
 unsigned int id_field_disp, unsigned int id_field_lambda,
 const Fem::Field::CFieldWorld& world,
 bool is_initial);

bool AddLinSys_Glut_Penalty_NewmarkBeta
(double dt, double gamma_newmark, double beta_newmark,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world,
 bool is_initial);

bool AddLinSys_Glut_Penalty_BackwardEular
(double dt, double damp_coeff,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world,
 double& strain_energy);

bool AddLinSys_Glut_Penalty_Sensitivity
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world);



unsigned int MakeField_GlueEdge_Lambda
(Fem::Field::CFieldWorld& world,
 unsigned int id_field, unsigned int id_ea1, unsigned int id_ea2, const int derivative_type );

unsigned int MakePartialField_GlueEdge_Penalty
(Fem::Field::CFieldWorld& world,
 unsigned int id_field, unsigned int id_ea1, unsigned int id_ea2, bool is_same_dir,
 bool& is_reordering_needed);



#endif
