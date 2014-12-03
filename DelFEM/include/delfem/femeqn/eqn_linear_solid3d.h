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
// The definition of each functions can be found in eqn_linear_solid_3d.cpp

#if !defined(EQN_LINEAR_SOLID_3D_H)
#define EQN_LINEAR_SOLID_3D_H

#include <vector>

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace Fem{
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

	// lienar elastic solid static
	bool AddLinSys_LinearSolid3D_Static
	(Fem::Ls::CLinearSystem_Field& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp );

	// linear elastic solid dynamic with newmark-beta time integration
	bool AddLinSys_LinearSolid3D_NonStatic_NewmarkBeta
	(double dt, double gamma, double beta,
	 Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp );

	// linear elastic solid static (saving stiffness matrix)
	bool AddLinSys_LinearSolid3D_Static_SaveStiffMat
	(Fem::Ls::CLinearSystem_Save& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp );

	// buld linear system for eigenanalysis of dynamic linear elastic solid
	bool AddLinSys_LinearSolid3D_Eigen
	(Fem::Ls::CLinearSystem_Eigen& ls,
	 double lambda, double myu, double rho,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp );
}
}

#endif
