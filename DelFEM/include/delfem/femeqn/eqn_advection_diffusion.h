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
@brief 移流拡散方程式の要素剛性作成関数のインターフェース
@author Nobuyuki Umetani
@sa http://ums.futene.net/wiki/FEM/46454D20666F7220416476656374696F6E2D446966667573696F6E204571756174696F6E.html
*/


#if !defined(EQN_ADVECTION_H)
#define EQN_ADVECTION_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace Fem{
namespace Ls
{
	class CLinearSystem_Field;
	class CLinearSystem_Save;
	class CLinearSystem_SaveDiaM_Newmark;
}
namespace Field
{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{

/*! @defgroup eqn_advection_diffusion 移流拡散方程式の連立一次方程式へのマージする関数群
@ingroup FemEqnMargeFunction

　　
静的な移流拡散方程式
@f$ \rho\frac{\partial\phi}{\partial t} = -\mu \nabla^2 \phi + f @f$

動的な移流拡散方程式
@f$\rho\frac{\partial\phi}{\partial t}+\rho v \cdot(\nabla\phi) = -\mu \nabla^2 \phi + f @f$

移流拡散方程式は，ある値が流れにしたがって移流されるとともに，拡散するような現象を表す方程式です．
例えば「川の流れにインクが落ち広がっていく」，「煙突から出る煙が風によって流される」などの様子を表します．

@pre この方程式の性質として，流れの上流に固定境界条件が設定されていないと連立一次方程式を解くことができません．

*/
/*!@{*/


/*!
@brief 静的な２次元移流拡散方程式の連立一次方程式へのマージ
@param [in,out] ls マージされる連立一次方程式
@param [in] myu 拡散係数
@param [in] source 面積に比例するソース項
@param [in] world 場管理クラス
@param [in] id_field_val 値場のID
@param [in] id_field_velo 速度場のID
@pre 値場はスカラー場でなければなりません．速度場は２次ベクトル場でなければなりません．
*/
bool AddLinSys_AdvectionDiffusion_Static(
    Fem::Ls::CLinearSystem_Field& ls,
	double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea = 0 );

/*!
@brief 静的な移流拡散方程式の連立一次方程式へのマージ(剛性行列保存)
@param [in,out] マージされる連立一次方程式
@param [in] myu 拡散係数
@param [in] source 面積に比例するソース項
@param [in] world 場管理クラス
@param [in] id_field_val 値場のID
@param [in] id_field_velo 速度場のID
@pre 値場はスカラー場でなければなりません．速度場は２次ベクトル場でなければなりません．
*/
bool AddLinSys_AdvectionDiffusion_Static(
	Fem::Ls::CLinearSystem_Save& ls,
	double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea = 0 );

/*!
@brief 動的な移流拡散方程式の連立一次方程式へのマージ
@param [in,out] マージされる連立一次方程式
@param [in] rho 慣性項 @f$\rho@f$
@param [in] myu 拡散係数 @f$\mu@f$
@param [in] source 面積に比例するソース項
@param [in] id_field_val 値場のID
@param [in] id_field_velo 速度場のID
@pre 値場はスカラー場でなければなりません．速度場は２次ベクトル場でなければなりません．
*/	
bool AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
	double dt, double gamma, 
	Fem::Ls::CLinearSystem_Field& ls,
	double rho, double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea = 0 );

/*!
@brief 動的な移流拡散方程式の連立一次方程式へのマージ(剛性行列保存)
@param [in,out] マージされる連立一次方程式(剛性行列保存型)
@param [in] rho 慣性項
@param [in] myu 拡散係数
@param [in] source 面積に比例するソース項
@param [in] id_field_val 値場のID
@param [in] id_field_velo 速度場のID
@pre 値場はスカラー場でなければなりません．速度場は２次ベクトル場でなければなりません．
*/	
bool AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
	Fem::Ls::CLinearSystem_SaveDiaM_Newmark& ls,
	double rho, double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea = 0 );

/*!@}*/
}
}

#endif
