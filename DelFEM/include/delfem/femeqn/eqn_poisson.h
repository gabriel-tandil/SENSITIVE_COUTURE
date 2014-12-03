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
@brief ポアソン方程式の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
@sa http://ums.futene.net/wiki/FEM/46454D20666F7220506F6973736F6E2773204571756174696F6E.html
*/

#if !defined(EQN_POISSON_H)
#define EQN_POISSON_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace Fem
{

namespace Ls
{
class CLinearSystem_Field;
class CLinearSystem_Save;
class CLinearSystem_Eigen;
class CPreconditioner;
}

namespace Field
{
class CField;
class CFieldWorld;
}

namespace Eqn
{

/*! @defgroup eqn_poisson ポアソン方程式をマージする関数群です
@ingroup FemEqnMargeFunction
　
次のような方程式を解きます．
\f$ \alpha \nabla^2 \phi = f \f$
*/
//! @{

/*!
@brief ポアソン方程式のマージ
@param [in,out] ls 連立一次方程式
@param [in] alpha 拡散係数 @f$ \alpha @f$
@param [in] source ソース項
@param [out] world 場管理クラス
@param [in] id_field_val 値場のID
*/
bool AddLinSys_Poisson(
		Fem::Ls::CLinearSystem_Field& ls,
		double alpha, double source,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて

/*!
@brief ポアソン方程式(剛性行列保存)のマージ
@param [in,out] ls 剛性行列を保存する連立一次方程式
@param [in] alpha 拡散係数 @f$ \alpha @f$
@param [in] source ソース項
@param [out] world 場管理クラス
@param [in] id_field_val 値場のID
*/
bool AddLinSys_Poisson(
		Fem::Ls::CLinearSystem_Save& ls,
		double alpha, double source,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea = 0 );	// 0だとid_field_valすべてについて
		
////////////////
bool AddLinearSystem_Wave(
		Fem::Ls::CLinearSystem_Eigen& ls,
		double rho, double alpha, double source,
		const Fem::Field::CFieldWorld& world,
		const unsigned int id_field_val );
// @}
}
}

#endif
