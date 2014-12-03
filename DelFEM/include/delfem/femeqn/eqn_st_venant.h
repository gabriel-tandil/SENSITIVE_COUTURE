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
// This file declar the functions for St.Venant-Kirchhoff solid equation.
// The function builds element matrices and marges them to linear system.
// The definition of each functions can be found in eqn_linear_solid.cpp

#if !defined(EQN_ST_VENSNT_H)
#define EQN_ST_VENANT_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include "delfem/linearsystem_interface_eqnsys.h"

namespace Fem{
namespace Ls{
	class CLinearSystem_Field;
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{
	// marge 2D static St.Venant-Kirchhoff material euqtion
	// ls [in,out] : linear system
	// lambda [in] : first lame constant
	// myu [in] : second lame constant
	// rho [in] : mass density
	// g_x [in] : gravity in x axes
	// g_y [in] : gravity in y axes	
	bool AddLinSys_StVenant2D_Static
	(Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double rho, double f_x, double f_y,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp,
	 unsigned int id_ea = 0);
	
	// marge 2D dynamic St.Venant-Kirchhoff material euqtion
	// newmark-beta time integration
	// ls [in,out] : linear system
	// lambda [in] : first lame constant
	// myu [in] : second lame constant
	// rho [in] : mass density
	// g_x [in] : gravity in x axes
	// g_y [in] : gravity in y axes
	bool AddLinSys_StVenant2D_NonStatic_NewmarkBeta
	(double dt, double gamma_newmark, double beta, 
	 Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp, 
	 bool is_initial,
	 unsigned int id_ea = 0 );
	
	bool AddLinSys_StVenant2D_NonStatic_BackwardEular
	(double dt, 
	 Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double rho, double g_x, double g_y,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp, 
	 const MatVec::CVector_Blk& velo_pre,
	 bool is_initial,
	 unsigned int id_ea = 0 );
	
	// marge 3D static St.Venant-Kirchhoff material euqtion
	// ls [in,out] : linear system
	// lambda [in] : first lame constant
	// myu [in] : second lame constant
	// rho [in] : mass density
	// g_x [in] : gravity in x axes
	// g_y [in] : gravity in y axes
	// g_z [in] : gravity in z axes
	bool AddLinSys_StVenant3D_Static
	(Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp,
	 unsigned int id_ea = 0 );
	

	// marge 3D dynamic St.Venant-Kirchhoff material euqtion
	// newmark-beta time integration
	// ls [in,out] : linear system
	// lambda [in] : first lame constant
	// myu [in] : second lame constant
	// rho [in] : mass density
	// g_x [in] : gravity in x axes
	// g_y [in] : gravity in y axes
	// g_z [in] : gravity in z axes
	bool AddLinSys_StVenant3D_NonStatic_NewmarkBeta
	(double dt, double gamma, double beta,
	 Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_field_disp,
	 bool is_initial,
	 unsigned int id_ea = 0);
}
}

#endif
