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
@brief 流体の方程式クラス(Fem::Eqn::CEqnSystem_Fluid2D, Fem::Eqn::CEqn_Fluid2D, Fem::Eqn::CEqn_Fluid3D,)のインターフェース
@author Nobuyuki Umetani
*/



#if !defined(EQN_OBJ_FLUID_H)
#define EQN_OBJ_FLUID_H

#include <vector>
#include <cassert>
#include "delfem/eqnsys.h"

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace Fem{

namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn
{

/*! 
@brief ２次元流体方程式クラス
@ingroup FemEqnObj
*/
class CEqn_Fluid2D
{
public:
	CEqn_Fluid2D(unsigned int id_ea, unsigned int id_velo, unsigned int id_press){
		m_rho = 1.0;
		m_myu = 1.0;
		m_g_x = 0.0;
		m_g_y = 0.0;
		this->m_id_ea = id_ea;
		this->m_IdFieldVelo  = id_velo;
		this->m_IdFieldPress = id_press;
		m_IsStokes = true;
		m_IsNavierStokes = false;
		m_IsNavierStokesALE = false;
	}
    bool IsNavierStokes() const
    {
		if( m_IsStokes ) return false;
		if( m_IsNavierStokes ) return true;
		if( m_IsNavierStokesALE ) return true;
		assert(0);
		return false;
	}
	//! 密度を設定する
	void SetRho(double rho){ m_rho = rho; }
    double GetRho() const { return m_rho; }
	//! 粘性係数を設定する
	void SetMyu(double myu){ m_myu = myu; }
    double GetMyu() const { return m_myu; }
	//! 外力項
	void SetBodyForce(double g_x, double g_y ){
		m_g_x = g_x;
		m_g_y = g_y;
	}
	void SetStokes(){
		m_IsStokes = true;
		m_IsNavierStokesALE = false;
		m_IsNavierStokes = false;
	}
	void SetNavierStokes(){
		m_IsNavierStokes = true;
		m_IsStokes = false;
		m_IsNavierStokesALE = false;
	}
	void SetNavierStokesALE(unsigned int id_field_msh_velo){
		m_IdFieldMshVelo = id_field_msh_velo;
		m_IsNavierStokesALE = true;
		m_IsNavierStokes = false;
		m_IsStokes = false;
	}
	bool AddLinSys(Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world );
	bool AddLinSys_NewmarkBetaAPrime( double dt, double gamma, double beta, bool is_initial, 
		Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world );

	unsigned int GetIdEA() const { return m_id_ea; }
	void SetIdEA(unsigned int id_ea){ m_id_ea = id_ea; }
	void SetIdFieldVelocity(unsigned int id_field_velo){ m_IdFieldVelo = id_field_velo; }
	void SetIdFieldPressure(unsigned int id_field_press){ m_IdFieldPress = id_field_press; }
private:
	unsigned int m_id_ea;

	unsigned int m_IdFieldVelo;
	unsigned int m_IdFieldPress;
	unsigned int m_IdFieldMshVelo;

	double m_rho, m_myu;
	double m_g_x, m_g_y;

	bool m_IsNavierStokes;
	bool m_IsNavierStokesALE;
	bool m_IsStokes;
};

/*! 
@brief ２次元流体,連成方程式クラス
@ingroup FemEqnSystem
*/
class CEqnSystem_Fluid2D : public CEqnSystem
{
public:
	CEqnSystem_Fluid2D();
	CEqnSystem_Fluid2D(Fem::Field::CFieldWorld& world);
	CEqnSystem_Fluid2D(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);

	virtual ~CEqnSystem_Fluid2D(){}
			
	////////////////
	// virtual関数
	virtual bool UpdateDomain_Field(unsigned int id_field, Fem::Field::CFieldWorld& world);
	virtual bool UpdateDomain_FieldVeloPress(unsigned int id_base_field_velo, unsigned int id_field_press, Fem::Field::CFieldWorld& world);
	virtual bool UpdateDomain_FieldElemAry(unsigned int id_field, unsigned int id_ea, Fem::Field::CFieldWorld& world);

	////////////////////////////////
	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	virtual bool Solve(Fem::Field::CFieldWorld& world);

	////////////////
	// 非virtual
	unsigned int GetIdField_Velo(){ return m_id_velo; }
	unsigned int GetIdField_Press(){ return m_id_press; }
	// 静的か動的かを定める
	void SetIsStationary(bool is_stat);
	// 静的か動的かを得る
	bool GetIsStationary(){ return this->m_IsStationary; }
	// 方程式を全てクリアする
	virtual void Clear()
	{ 
		CEqnSystem::Clear();
		this->m_aIdFixField.clear();
		m_id_press = 0;
		m_id_velo = 0;
		m_id_msh_velo = 0;
        m_id_force = 0;
		m_aEqn.clear(); 
		m_rho_back = 1.0;
		m_myu_back = 1.0;
		m_IsStationary = true;
		m_IsntInterpolationBubble = true;
	//	m_IsntCombine = true;
		m_IsntCombine = false;
	}
	// 方程式の要素IDの配列を得る
	std::vector<unsigned int> GetAry_EqnIdEA(){
		std::vector<unsigned int> aId;
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){
			aId.push_back( m_aEqn[ieqn].GetIdEA() );
		}
		return aId;
	}
	// 要素IDがid_eaの方程式を得る
   CEqn_Fluid2D GetEquation(unsigned int id_ea) const;
	// 要素単位で方程式を定める．要素が存在しなければ方程式を追加する
	bool SetEquation( const CEqn_Fluid2D& eqn );
	// 要素IDを更新する
	void UpdateIdEA( const std::vector<unsigned int>& map_ea2ea ){
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){
			const unsigned int id_ea0 = m_aEqn[ieqn].GetIdEA();
			if( id_ea0 >= map_ea2ea.size() ) continue;
			const unsigned int id_ea1 = map_ea2ea[id_ea0];
			m_aEqn[ieqn].SetIdEA(id_ea1);
		}
	}
	////////////////
	void SetNavierStokes(){
        m_is_stokes_back = false;
        this->SetIsStationary(false);
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){
			m_aEqn[ieqn].SetNavierStokes();
		}
	}	
    void SetNavierStokesALE(unsigned int id_field_msh_velo)
    {
        m_is_stokes_back = false;
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){
			m_aEqn[ieqn].SetNavierStokesALE(id_field_msh_velo);
		}
		m_id_msh_velo = id_field_msh_velo;
    }
    void SetForceField(unsigned int id_field_force){
        this->m_id_force = id_field_force;
    }
	void SetInterpolationBubble(){
		if( !m_IsntInterpolationBubble ) return;
		m_IsntInterpolationBubble = false;
		this->m_IsntCombine = true;
		this->ClearLinearSystemPreconditioner();
	}
	void UnSetInterpolationBubble(){
		if( m_IsntInterpolationBubble ) return;
		m_IsntInterpolationBubble = true;
		this->ClearLinearSystemPreconditioner();
	}
    ////////////////////////////////
	void SetRho(double rho){ 	
        m_rho_back = rho;
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){ m_aEqn[ieqn].SetRho(rho); }
	}
	void SetMyu(double myu){ 
        m_myu_back = myu;
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){ m_aEqn[ieqn].SetMyu(myu); }
	}	
	void SetStokes(){			
        m_is_stokes_back = true;
		for(unsigned int ieqn=0;ieqn<this->m_aEqn.size();ieqn++){ m_aEqn[ieqn].SetStokes(); }
	}
private:
	// 仮想化可能関数
	virtual double MakeLinearSystem(const Fem::Field::CFieldWorld& world);
	virtual bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
	bool EqnationProperty(bool& is_asym){
		is_asym = false;
		for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
			if( m_aEqn[ieqn].IsNavierStokes() ){
				is_asym = true;
				break;
			}
		}
		return true;
	}
private:
	unsigned int m_id_press;
	unsigned int m_id_velo;
	unsigned int m_id_msh_velo;
    unsigned int m_id_force;

	////////////////
	// 境界条件について
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

    ////////////////////////////////
    // 物性値について
	std::vector<CEqn_Fluid2D> m_aEqn;
    double m_rho_back;  // 新しく方程式要素が出来たときはこの値で初期化
    double m_myu_back;  // 新しく方程式要素が出来たときはこの値で初期化
    bool m_is_stokes_back;  // 新しく方程式要素が出来たときはこの値で初期化
    ////////////////
	bool m_IsStationary;
	bool m_IsntInterpolationBubble;
	bool m_IsntCombine;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/*! 
@brief ３次元流体方程式クラス
@ingroup FemEqnObj
*/
class CEqn_Fluid3D : public CEqnSystem
{
public:
	CEqn_Fluid3D();
	CEqn_Fluid3D(Fem::Field::CFieldWorld& world);
	CEqn_Fluid3D(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
			
	// 仮想化可能関数
	virtual bool SetDomain(unsigned int id_base, Fem::Field::CFieldWorld& world);
	virtual bool Solve(Fem::Field::CFieldWorld& world);

	////////////////////////////////
	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	void SetGravitation(double g_x, double g_y, double g_z){
		m_g_x = g_x;
		m_g_y = g_y;
		m_g_z = g_z;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	// 非virtual
	unsigned int GetIDField_Velo(){ return m_id_velo; }
	unsigned int GetIDField_Press(){ return m_id_press; }

			
private:
	// 仮想化可能関数
	virtual double MakeLinearSystem(const Fem::Field::CFieldWorld& world);
	virtual bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);

private:
	////////////////
	// 境界条件について
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

	unsigned int m_id_press;
	unsigned int m_id_velo;

    double m_myu, m_rho;
	double m_g_x, m_g_y, m_g_z;
};

}
}

#endif

