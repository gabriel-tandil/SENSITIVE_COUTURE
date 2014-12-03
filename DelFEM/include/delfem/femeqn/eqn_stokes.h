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
@brief ストークス方程式の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
*/


#if !defined(EQN_STOKES_H)
#define EQN_STOKES_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace Fem{
namespace Ls{
	class CLinearSystem_Field;
	class CLinearSystem_LinEqn;
	class CPreconditioner;
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn{

/*! @defgroup eqn_stokes ストークス方程式を連立一次方程式にマージする関数群
@ingroup FemEqnMargeFunction


静的ストークス方程式
@f$ -\nabla v^2 = f, \nabla\cdot v = 0 @f$

動的ストークス方程式
@f$ \rho\frac{\partial v}{\partial t} = \nabla^2 v + f, \nabla\cdot v = 0 @f$
*/
//! @{

/*!
@brief２次元の静的stokes方程式のマージ
@param [in,out] ls 連立一次方程式
@param [in] alpha 粘性係数 @f$ \alpha @f$
@param [in] f_x x方向の体積力
@param [in] f_y y方向の体積力
*/
bool AddLinSys_Stokes2D_Static(
	double alpha, 
	double f_x, double f_y,
	Fem::Ls::CLinearSystem_Field& ls,
	const unsigned int id_field_velo, unsigned int id_field_press, const Fem::Field::CFieldWorld& world,
	unsigned int id_ea = 0 );

/*!
@brief２次元の動的stokes方程式のマージ
@param [in,out] ls 連立一次方程式
@param [in] rho 密度 @f$ \rho @f$
@param [in] alpha 粘性係数 @f$ \alpha @f$
@param [in] f_x x方向の体積力
@param [in] f_y y方向の体積力
*/
bool AddLinSys_Stokes2D_NonStatic_Newmark(
	double rho, double alpha, 
	double g_x, double g_y,
	double gamma, double dt,
	Fem::Ls::CLinearSystem_Field& ls,
	const unsigned int id_field_velo, unsigned int id_field_press, const Fem::Field::CFieldWorld& world,
	unsigned int id_ea = 0);

bool AddLinSys_Stokes3D_Static(
		double alpha, 
		double rho, double g_x, double g_y, double g_z,
		Fem::Ls::CLinearSystem_Field& ls,
		const unsigned int id_field_velo, unsigned int id_field_press, const Fem::Field::CFieldWorld& world );

//! @}
}
}

#endif

