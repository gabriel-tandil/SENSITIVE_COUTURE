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
@brief ナビア・ストークス方程式の要素剛性作成部のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(EQN_NAVIER_STOKES_H)
#define EQN_NAVIER_STOKES_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace Fem{
namespace Ls{
class CLinearSystem_Field;
}

namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn{

/*! @defgroup eqn_ns ナビア・ストークス方程式の連立一次方程式へのマージ関数群
@ingroup FemEqnMargeFunction
*/

/*! 
@brief ナビア・ストークス方程式の連立一次方程式を作成するクラス
@param [in] rho 質量密度
@param [in] alpha 粘性
@param [in] g_x ｘ方向の体積力
@param [in] g_y ｙ方向の体積力
@param [in] gamma Newmark法のパラメータ
@param [in] dt 時間刻み
@param [in] id_field_velo 流速場のID
@param [in] id_field_press 圧力場のID
*/
bool AddLinSys_NavierStokes2D_NonStatic_Newmark(
	double dt, double gamma, Fem::Ls::CLinearSystem_Field& ls, 
	double rho, double alpha,      double g_x, double g_y,
	const unsigned int id_field_velo, unsigned int id_field_press, 
    const Fem::Field::CFieldWorld& world,
	unsigned int id_ea = 0 );

/*! 
@brief ALE座標におけるナビア・ストークス方程式の連立一次方程式を作成するクラス
@param [in] rho 質量密度
@param [in] alpha 粘性
@param [in] g_x ｘ方向の体積力
@param [in] g_y ｙ方向の体積力
@param [in] gamma Newmark法のパラメータ
@param [in] dt 時間刻み
@param [in] id_field_velo 流速場のID
@param [in] id_field_press 圧力場のID
*/
bool AddLinSys_NavierStokesALE2D_NonStatic_Newmark(
	double dt, double gamma, Fem::Ls::CLinearSystem_Field& ls, 
	double rho, double alpha,      double g_x, double g_y,
	const unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_msh_velo,
	const Fem::Field::CFieldWorld& world );

/*! 
@brief 熱の浮力をうけるナビア・ストークス方程式の連立一次方程式を作成するクラス
@param [in] rho 質量密度
@param [in] alpha 粘性
@param [in] g_x ｘ方向の体積力
@param [in] g_y ｙ方向の体積力
@param [in] gamma Newmark法のパラメータ
@param [in] dt 時間刻み
@param [in] id_field_velo 流速場のID
@param [in] id_field_press 圧力場のID
*/	
bool AddLinSys_NavierStokes2DThermalBuoy_NonStatic_Newmark(
	double rho, double alpha, double g_x, double g_y,
	double gamma, double dt,
	Fem::Ls::CLinearSystem_Field& ls, 
	const unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_temp,
	const Fem::Field::CFieldWorld& world );

}
}

#endif
