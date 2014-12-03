/*
DelFEM (Finite Element Analysis)
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// AUTHOR
// Nobuyuki Umetani

// DESCRIPTION
// This file declar the functions for linear solid equation.
// The function builds element matrices and marges them to linear system.
// The definition of each functions can be found in eqn_linear_solid.cpp

#if !defined(EQN_LINEAR_SOLID_2D_H)
#define EQN_LINEAR_SOLID_2D_H


#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace MatVec
{
	class CVector_Blk;
}

namespace Fem
{
namespace Ls
{
	class CLinearSystem_Field;
	class CLinearSystem_Save;
	class CLinearSystem_SaveDiaM_NewmarkBeta;
	class CLinearSystem_Eigen;
	class CPreconditioner;
}

namespace Field
{
	class CField;
	class CFieldWorld;
}

namespace Eqn{

	class ILinearSystem_Eqn;

	////////////////////////////////////////////////////////////////
	// 2-dimensional euqations

	// 2D static linear elastic solid 
	bool AddLinSys_LinearSolid2D_Static
	(Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu, double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp,
	 unsigned int id_ea = 0 );

	// 2D static linear elastic solid with thermal-stress
	bool AddLinSys_LinearSolidThermalStress2D_Static
	(Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu, double rho, double g_x, double g_y, 	double thermoelastic, 
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp,
	 unsigned int id_ea = 0 );

	// 2D dynamic linear elastic solid
	// this function use backward eular time integration method
	bool AddLinSys_LinearSolid2D_NonStatic_BackwardEular
	(double dt, Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu, double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
	 const MatVec::CVector_Blk& velo_pre,
	 bool is_initial = true,	 
	 unsigned int id_ea = 0 );
	
	// 2D dynamic linear elastic solid
	// this function use newmark-beta time integration method
	bool AddLinSys_LinearSolid2D_NonStatic_NewmarkBeta
	(double dt, double gamma, double beta, Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu, double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
	 bool is_initial = true,
	 unsigned int id_ea = 0 );

	// 2D dynamic lienar solid with thermal-stress
	bool AddLinSys_LinearSolidThermalStress2D_NonStatic_NewmarkBeta
	(double dt, double gamma, double beta, Fem::Ls::CLinearSystem_Field& ls,
	 double lambda, double myu, double rho, double g_x, double g_y, double thermoelastic, 
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, unsigned int id_field_temp, 
	 bool is_inital = true,
	 unsigned int id_ea = 0);

	// static linear elastic solid (saving stiffness matrix)
	bool AddLinSys_LinearSolid2D_Static_SaveStiffMat
	(Fem::Ls::CLinearSystem_Save& ls,
	 double lambda, double myu, double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
	 unsigned int id_ea = 0 );

	// dynamic linear elastic solid (saving stiffness matrix)
	bool AddLinSys_LinearSolid2D_NonStatic_Save_NewmarkBeta
	(Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls,
	 double lambda, double myu, double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
	 unsigned int id_ea = 0 );

	// buidling matirx for eigennalysis
	bool AddLinSys_LinearSolid2D_Eigen
	(Fem::Ls::CLinearSystem_Eigen& ls,
	 double lambda, double myu, double rho,
	 const Fem::Field::CFieldWorld& world, unsigned int id_field_disp, 
	 unsigned int id_ea = 0 );
	
}
}

#endif
