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
@brief 固体の方程式クラスのインターフェース
@author Nobuyuki Umetani
*/

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#if !defined(EQN_SYS_SOLID_H)
#define EQN_SYS_SOLID_H

#include <vector>
#include <iostream>

#include "delfem/eqnsys.h"

namespace Fem{

namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{

/*! 
@brief ３次元固体方程式クラス
(Linearと言っているのは物性地が線形の意味．幾何学的非線形を意味しているのではない．紛らわしいので将来的に変更予定)
@ingroup FemEqnObj
*/
class CEqn_Solid3D_Linear : public CEqnSystem
{
public:
	CEqn_Solid3D_Linear();
	CEqn_Solid3D_Linear(Fem::Field::CFieldWorld& world);
	CEqn_Solid3D_Linear(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
			
	// 仮想化可能関数
	virtual bool SetDomain_Field(unsigned int id_field, Fem::Field::CFieldWorld& world);
	virtual bool Solve(Fem::Field::CFieldWorld& world);

	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	unsigned int GetIdField_Disp() const { return m_IdFieldDisp; }

	// 非virtual
	void SetYoungPoisson( double young, double poisson ){
		m_lambda = young*poisson / ( (1.0+poisson)*(1-2.0*poisson) );
		m_myu = young / ( 2.0*(1.0+poisson) );
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	//! 重力加速度設定
	void SetGravitation( double g_x, double g_y, double g_z ){
		m_g_x = g_x;
		m_g_y = g_y;
		m_g_z = g_z;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	void SetGeometricalNonLinear(){
		if( m_IsGeomNonlin ) return;
		m_IsGeomNonlin = true;
		m_IsSaveStiffMat = false;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void UnSetGeometricalNonLinear(){
		if( !m_IsGeomNonlin ) return;
		m_IsGeomNonlin = false;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	void SetSaveStiffMat(){
		if( m_IsSaveStiffMat || m_IsGeomNonlin || !m_IsStationary ) return;
		m_IsSaveStiffMat = true;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void UnSetSaveStiffMat(){
	if( !m_IsSaveStiffMat ) return;
		m_IsSaveStiffMat = false;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	void SetStationary(){
		if( m_IsStationary ) return;
		m_IsStationary = true;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void UnSetStationary(){
		if( !m_IsStationary ) return;
		m_IsSaveStiffMat = false;
		m_IsStationary = false;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
		// TODO : 加速度場を０に設定しないとダメだね．
	}
private:
	// 仮想化可能関数
	double MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_initial);
	bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
private:
	bool m_IsGeomNonlin;
	bool m_IsSaveStiffMat;
	bool m_IsStationary;
	double m_lambda, m_myu, m_rho;
	double m_g_x, m_g_y, m_g_z;
	
	unsigned int m_IdFieldDisp;
	////////////////
	// 境界条件について
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*! 
@brief ２次元固体方程式クラス
@ingroup FemEqnObj
*/
class CEqn_Solid2D
{
public:
	CEqn_Solid2D(unsigned int id_ea, unsigned int id_field_disp)
		: m_id_ea(id_ea), m_IdFieldDisp(id_field_disp)
	{
		m_lambda = 0.0; m_myu = 1.0; m_rho = 1.0; m_g_x = 0.0; m_g_y = 0.0;
		m_IsGeomNonlin = false;
		m_IdFieldTemperature = 0;
	}
	CEqn_Solid2D()
		: m_id_ea(0), m_IdFieldDisp(0)
	{
		m_lambda = 0.0; m_myu = 1.0; m_rho = 1.0; m_g_x = 0.0; m_g_y = 0.0;
		m_IsGeomNonlin = false;
      m_IdFieldTemperature = 0;
  }
   CEqn_Solid2D( const CEqn_Solid2D& eqn ){
       this->m_id_ea = eqn.m_id_ea;
       this->m_IdFieldDisp = eqn.m_IdFieldDisp;
       this->m_IdFieldTemperature = eqn.m_IdFieldTemperature;
       this->CopyParameters(eqn);
   }

	unsigned int GetIdField_Disp() const { return m_IdFieldDisp; }
	void CopyParameters( const CEqn_Solid2D& eqn ){
		m_lambda = eqn.m_lambda; 
		m_myu    = eqn.m_myu;
		m_rho    = eqn.m_rho;
		m_g_x    = eqn.m_g_x;
		m_g_y    = eqn.m_g_y;
		m_IdFieldTemperature = eqn.m_IdFieldTemperature;
		m_IsGeomNonlin = eqn.m_IsGeomNonlin;
	}
	// Setメソッド
	void SetYoungPoisson( double young, double poisson, bool is_plane_stress ){
		m_lambda = young*poisson / ( (1.0+poisson)*(1-2.0*poisson) );
		m_myu = young / ( 2.0*(1.0+poisson) );
		if( !is_plane_stress ){
			m_lambda = 2*m_lambda*m_myu / (m_lambda+2*m_myu);
		}
	}
	//! ラメ定数の取得
	void GetLambdaMyu( double& lambda, double& myu) const{
		lambda = m_lambda;
		myu    = m_myu;
	}
	//! 質量密度の取得
	double GetRho() const{ return m_rho; }
	void GetYoungPoisson( double& young, double& poisson) const{
//		if( is_plane_stress ){
			poisson = m_lambda*0.5/(m_lambda+m_myu);
			young   = 2*m_myu*(1+poisson);
/*		}
		else{

		}*/
	}
	//! 質量密度の設定
	void SetRho(  double rho){ m_rho = rho; }
	//! 体積力の設定
	void SetGravitation( double g_x, double g_y ){
		m_g_x = g_x;
		m_g_y = g_y;
	}
	void GetGravitation( double& g_x, double& g_y ){
		g_x = m_g_x;
		g_y = m_g_y;
	}
	void SetIdFieldDisp(unsigned int id_field_disp){ m_IdFieldDisp = id_field_disp; }
	void SetGeometricalNonlinear(bool is_nonlin){   m_IsGeomNonlin = is_nonlin;  }
	void SetThermalStress(unsigned int id_field_temp){ m_IdFieldTemperature = id_field_temp; }
	// Getメソッド
	unsigned int GetIdEA() const { return m_id_ea; }
	bool IsTemperature()  const { return m_IdFieldTemperature!=0;  }
	bool IsGeometricalNonlinear() const { return m_IsGeomNonlin; }

	// 連立一次方程式マージメソッド
	bool AddLinSys_NewmarkBetaAPrime( double dt, double gamma, double beta, bool is_inital, Ls::CLinearSystem_Field& ls, const Field::CFieldWorld& world );
	bool AddLinSys_NewmarkBetaAPrime_Save( Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls, const Field::CFieldWorld& world );
	bool AddLinSys( Ls::CLinearSystem_Field& ls, const Field::CFieldWorld& world );
	bool AddLinSys_Save( Ls::CLinearSystem_Save& ls, const Field::CFieldWorld& world );
private:
	unsigned int m_id_ea;
	unsigned int m_IdFieldDisp;
	unsigned int m_IdFieldTemperature;

	bool m_IsGeomNonlin;
	double m_lambda, m_myu, m_rho;
	double m_g_x, m_g_y;
};


//2Dと3Dは分ける．外力の次元が違うし、平面応力平面歪の設定があるから．


/*!
@brief ２次元固体，連成方程式クラス
@ingroup FemEqnSystem
*/
class CEqnSystem_Solid2D : public CEqnSystem
{
public:
	//! デフォルトコンストラクタ
	CEqnSystem_Solid2D();
	CEqnSystem_Solid2D(Fem::Field::CFieldWorld& world);
	CEqnSystem_Solid2D(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
			
	// 仮想化可能関数
	virtual bool UpdateDomain_Field(unsigned int id_field, Fem::Field::CFieldWorld& world);
	virtual bool SetDomain_FieldEA(unsigned int id_field, unsigned int id_ea, Fem::Field::CFieldWorld& world);
	//! 方程式を解く
	virtual bool Solve(Fem::Field::CFieldWorld& world);
	virtual void Clear(){
		CEqnSystem::Clear();
		m_IdFieldDisp = 0;
		m_aIdFixField.clear();
		m_IsSaveStiffMat = false;
		m_IsStationary = true;
		m_aEqn.clear();
		m_young_back = 1.0;
		m_poisson_back = 0;
		m_is_plane_stress_back = true;
		m_rho_back = 1;
		m_is_geom_nonlin_back = false;
		m_aLoad.clear();
		m_num_iter = 100;
		m_conv_ratio = 1.0e-6;
	}

	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	unsigned int GetIdField_Disp() const { return m_IdFieldDisp; }

	bool ToplogicalChangeCad_InsertLoop(Fem::Field::CFieldWorld& world, 
		unsigned int id_l_back, unsigned id_l_ins);

	////////////////////////////////
	// 非virtual
	
	//! 方程式の取得，設定
	//! @{
	bool SetEquation( const CEqn_Solid2D& eqn );
   CEqn_Solid2D GetEquation(unsigned int id_ea) const;
	//! @}

	////////////////
	// 各領域の方程式のパラメータを一度に変えるオプション

	//! @{
	void SetYoungPoisson( double young, double poisson, bool is_plane_stress ); //!< 剛性パラメータ設定
	void SetRho( double rho ); //!< 密度設定
	void SetGravitation( double g_x, double g_y ); //!< 重力加速度設定
	void SetThermalStress(unsigned int id_field_temp);	//!< 熱応力を考慮する０を代入するとで解除される
	void SetGeometricalNonlinear(bool is_nonlin);	//!< 幾何学的非線形性を考慮する
	void SetStationary(bool is_stationary);	//!< 静的問題に設定
	void SetSaveStiffMat(bool is_save);		//!< 剛性行列を保存する
	
	bool IsStationary() const { return m_IsStationary; }
	bool IsSaveStiffMat() const { return m_IsSaveStiffMat; }
	void GetYoungPoisson(double& young, double& poisson, bool& is_plane_str) const{
		young = m_young_back;
		poisson = m_poisson_back;
		is_plane_str = m_is_plane_stress_back;
	}
	bool IsGeometricalNonlinear() const{ return m_is_geom_nonlin_back; }
	double GetRho(){ return m_rho_back; }

	//! @}

	//! 荷重の設定
	void SetLoad(double load, unsigned int id_ea, Fem::Field::CFieldWorld& world);
	void ClearLoad(unsigned int id_ea);
	void ClearLoad(){ m_aLoad.clear(); }

	//! 変位から応力の値を計算して場(ID:id_field)にセットする
	bool SetEquivStressValue(unsigned int id_field, Fem::Field::CFieldWorld& world);
	bool SetStressValue(unsigned int id_field, Fem::Field::CFieldWorld& world);

private:
	// 仮想化可能関数
	bool EqnationProperty(bool& is_nonlin);

	double MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_inital);
	bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
	bool InitializePreconditioner();
	bool MakePreconditioner();
private:
	////////////////
	// 境界条件について
	unsigned int m_IdFieldDisp;
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

	bool m_IsSaveStiffMat;
	bool m_IsStationary;

	std::vector<CEqn_Solid2D> m_aEqn;
    double m_young_back;  // 新しく方程式要素が出来たときはこの値で初期化
    double m_poisson_back;  // 新しく方程式要素が出来たときはこの値で初期化
    bool m_is_plane_stress_back;  // 新しく方程式要素が出来たときはこの値で初期化
    double m_rho_back;  // 新しく方程式要素が出来たときはこの値で初期化
	bool m_is_geom_nonlin_back;

	std::vector< std::pair<unsigned int,double> > m_aLoad;

	unsigned int m_num_iter;
	double m_conv_ratio;
};

}
}

#endif

