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
@brief ヘルムホルツ方程式の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
@sa http://ums.futene.net/wiki/FEM/46454D20666F7220536F756E64204669656C6420416E616C79736973.html
*/

#if !defined(EQN_HELMHOLTZ_H)
#define EQN_HELMHOLTZ_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace Fem
{

namespace Ls{
	class CZLinearSystem;
	class CZLinearSystem_GeneralEigen;
	class CPreconditioner;
}
namespace Field{
	class CField;
	class CFieldWorld;
}

namespace Eqn{

/*!
@brief Helmholtz方程式のマージ
@param[in] wave_length 波長
*/
bool AddLinSys_Helmholtz(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

/*!
@brief Helmholtz方程式のマージ
@param[in] wave_length 波長
*/
bool AddLinSys_Helmholtz_AxalSym(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

/*!
@brief Helmholtz方程式の放射境界条件のマージ
@param[in] wave_length 波長
*/
bool AddLinSys_SommerfeltRadiationBC(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

/*!
@brief Helmholtz方程式の放射境界条件のマージ
@param[in] wave_length 波長
*/
bool AddLinSys_SommerfeltRadiationBC_AxalSym(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

/*!
@brief Helmholtz方程式のマージ
@param[in] wave_length 波長
*/
bool AddLinSys_MassMatrixEigen_AxalSym(
		Fem::Ls::CZLinearSystem_GeneralEigen& ls,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

}
}

#endif
