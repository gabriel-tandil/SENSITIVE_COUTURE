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


/*! @file
@brief DKT薄板方程式オブジェクトのインターフェース
@author Nobuyuki Umetani
*/


#if !defined(EQN_DKT_H)
#define EQN_DKT_H

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace Fem{
namespace Ls{
	class CLinearSystem_Field;
	class CLinearSystem_Save;	
	class CLinearSystem_SaveDiaM_NewmarkBeta;
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{
	bool AddLinearSystem_DKT2D_Static(
			Fem::Ls::CLinearSystem_Field& ls,
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );

	bool AddLinearSystem_DKT3D_Linear_Static(
			Fem::Ls::CLinearSystem_Field& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
		
	bool AddLinearSystem_DKT3D_Linear_Static_Save(
			Fem::Ls::CLinearSystem_Save& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
		
	bool AddLinearSystem_DKT3D_Linear_NonStatic(
			double dt, double gamma_newmark, double beta_newmark,
			Fem::Ls::CLinearSystem_Field& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
		
	bool AddLinearSystem_DKT3D_Linear_NonStatic_Save(
			Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
		
	bool AddLinearSystem_DKT3D_NonLinear_Static(
			Fem::Ls::CLinearSystem_Field& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
	
	bool AddLinearSystem_DKT3D_NonLinear_NonStatic(
			double dt, double gamma_newmark, double beta_newmark, bool is_first_itr,
			Fem::Ls::CLinearSystem_Field& ls,
			double young, double poisson, double thickness, double arearho,
			double g_x, double g_y, double g_z, double press, 
			const Fem::Field::CFieldWorld& world,
			const unsigned int id_field_deflect, unsigned int id_field_rot, 
			unsigned int id_ea = 0 );
		
}
}

#endif

