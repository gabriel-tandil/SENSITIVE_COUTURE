/*
 *  eqn_contact3d.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#if !defined(EQN_CONTACT_3D_H)
#define EQN_CONTACT_3D_H

class CContactTarget3D;

class CFrictionPoint
{
public:
  CFrictionPoint(){
    is_pin = false;
  }
public:
	double pd;
	double aloc[3];
	int itype_contact;	// 0:not_contact 1:static contact 2:dynamic contact
  bool is_pin;
};

void Update_FrictionalContact
(const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp,
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint );


bool AddLinSys_FrictionalContact_Penalty_NonStatic_BackwardEular
(double dt,
 Fem::Ls::CLinearSystem_Field& ls, 
 const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp, 
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint );

bool AddLinSys_FrictionalContact_Penalty_NonStatic_Sensitivity
(Fem::Ls::CLinearSystem_Field& ls, 
 const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp, 
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint );


#endif
