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
@brief シェルの方程式オブジェクトのインターフェース
@author Nobuyuki Umetani
*/


#if !defined(EQN_OBJ_SHELL_H)
#define EQN_OBJ_SHELL_H

#include <vector>
#include "delfem/eqnsys.h"

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

namespace Fem{
namespace Field
{
class CField;
class CFieldWorld;
}
namespace Eqn{
//! DKTシェル要素方程式クラス
class CEqnSystem_DKT : public CEqnSystem
{
public:
	CEqnSystem_DKT();
	CEqnSystem_DKT(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
			
	// 仮想化可能関数
	virtual bool SetDomain_FieldElemAry(unsigned int id_base, unsigned int id_ea, Fem::Field::CFieldWorld& world);
	virtual bool SetDomain_Field(unsigned int id_base, Fem::Field::CFieldWorld& world);
	virtual bool Solve(Fem::Field::CFieldWorld& world);

	virtual bool         AddFixField(                    unsigned int id_field, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry(                  unsigned int id_ea,    Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	// 非virtual
	unsigned int GetIDField_Disp(){ return m_id_disp; }
	unsigned int GetIDField_Rot(){ return m_id_rot; }

	// 非virtual
	void SetYoungPoisson( double young, double poisson ){
		m_young = young;
		m_poisson = poisson;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void SetThickness( double thick ){
		m_thick = thick;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void SetArealDensity( double arho ){
		m_arealrho = arho;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void SetPressure( double press ){
		m_press = press;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	// 重力加速度設定
	void SetGravitation( double g_x, double g_y, double g_z ){
		m_g_x = g_x;
		m_g_y = g_y;
		m_g_z = g_z;
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	void SetGeometricalNonlinear(bool is_nl){
		if( m_IsGeomNonlin == is_nl ) return;
		m_IsGeomNonlin = is_nl;
		if( m_IsSaveStiffMat ){
			m_IsSaveStiffMat = false;
			this->ClearLinearSystemPreconditioner();
		}
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}
	void SetSaveStiffMat(bool is_save){
		if( m_IsSaveStiffMat == is_save ) return;
		if( m_IsGeomNonlin && is_save ) return;
		m_IsSaveStiffMat = is_save;
		this->ClearLinearSystemPreconditioner();
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
	}

	void SetStationary(bool is_stat){
		if( m_IsStationary == is_stat) return;
		m_IsStationary = is_stat;
		this->ClearLinearSystemPreconditioner();
		this->m_is_cleared_value_ls = true;
		this->m_is_cleared_value_prec = true;
		// TODO : 加速度場を０に設定しないとダメだね．
	}
private:
	// 仮想化可能関数
	double MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_initial);
	bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
private:
	unsigned int m_id_disp;
	unsigned int m_id_rot;
    std::vector< std::pair<unsigned int,int> > m_aIdFixField;
    ////////////////
    // 物性値
	bool m_IsGeomNonlin;
	bool m_IsStationary;
	bool m_IsSaveStiffMat;
	double m_young, m_poisson, m_arealrho, m_thick;
	double m_g_x, m_g_y, m_g_z;
	double m_press;
	bool isnt_direct;
};

}	// Eqn
}	// Fem

#endif
