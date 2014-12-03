/*
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
@brief 連立一次方程式クラス(Fem::Ls::CLinearSystem_Field)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(LINEAR_SYSTEM_FIELD_H)
#define LINEAR_SYSTEM_FIELD_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <assert.h>

#include "delfem/field.h"
#include "delfem/linearsystem_interface_eqnsys.h"
#include "delfem/ls/linearsystem_interface_solver.h"
#include "delfem/ls/linearsystem.h"

namespace MatVec{
class CVector_Blk;
class CZVector_Blk;
class CMat_BlkCrs;
class CMatDia_BlkCrs;
class CDiaMat_Blk;
class CBCFlag;
}

namespace Fem{
namespace Field{
	class CFieldWorld;
}
namespace Ls{
  
void BoundaryCondition
(unsigned int id_field, const Fem::Field::ELSEG_TYPE& elseg_type,  
 MatVec::CBCFlag& bc_flag, const Fem::Field::CFieldWorld& world,  
 unsigned int ioffset=0);  
  
void BoundaryCondition
(unsigned int id_field, const Fem::Field::ELSEG_TYPE& elseg_type, unsigned int idofns, 
 MatVec::CBCFlag& bc_flag, const Fem::Field::CFieldWorld& world);

  
/*! 
@brief 連立一次方程式クラス
@ingroup FemLs
*/
class CLinearSystem_Field : public LsSol::ILinearSystem_Sol, public Eqn::ILinearSystem_Eqn
{
public:
  CLinearSystem_Field(){}	//! default constructor
  virtual ~CLinearSystem_Field(){ this->Clear(); }	//! destructor
	
	////////////////
  virtual void Clear(){ m_ls.Clear(); m_aSegField.clear(); }

	////////////////////////////////
	// function for marge element matrix

	//! Get residual vector for node location (elseg_type) in field (id_field)
  virtual MatVec::CVector_Blk& GetResidual(unsigned int id_field, Field::ELSEG_TYPE elseg_type, const Field::CFieldWorld& world);
  virtual MatVec::CVector_Blk& GetUpdate(  unsigned int id_field, Field::ELSEG_TYPE elseg_type, const Field::CFieldWorld& world);
	//! Get square sub-matrix from diagonal part of full linear system
  virtual MatVec::CMatDia_BlkCrs& GetMatrix
  (unsigned int id_field, Fem::Field::ELSEG_TYPE elseg_type, const Fem::Field::CFieldWorld& world);
	//! Get non-square sub-matrix from off-diagonal part of full linear system
  virtual MatVec::CMat_BlkCrs& GetMatrix
  (unsigned int id_field_col, Field::ELSEG_TYPE elseg_type_col,
   unsigned int id_field_row, Field::ELSEG_TYPE elseg_type_row,
   const Field::CFieldWorld& world);
	
	////////////////////////////////
	// function for marge
	
	//! Initialize before marge (set zero value to residual & matrix)	
    virtual void InitializeMarge(){ m_ls.InitializeMarge(); }
	//! Finalization before marge (set boundary condition, return residual square norm)
    virtual double FinalizeMarge(){ return m_ls.FinalizeMarge(); }

    ////////////////////////////////
    // パターン初期化関数

	//! fieldで初期化する、fieldの中の非ゼロパターンを作る
	virtual bool AddPattern_Field(const unsigned int id_field, const Field::CFieldWorld& world);
	//! fieldで初期化する、fieldとfield-field2の中の非ゼロパターンを作る
	virtual bool AddPattern_Field(unsigned int id_field, unsigned int id_field2, const Field::CFieldWorld& world);
	//! fieldとfield2がパターンが同じだとして，ブロックが結合された一つの行列を作る
	virtual bool AddPattern_CombinedField(unsigned id_field, unsigned int id_field2, const Field::CFieldWorld& world);

	////////////////////////////////
	// function for fixed boundary condition

	void ClearFixedBoundaryCondition();
	//! set fix boundary condition to dof (idof) in field (id_field)
	bool SetFixedBoundaryCondition_Field( unsigned int id_field, unsigned int idofns, const Field::CFieldWorld& world );
	//! set fix boundary condition to all dof in field (id_field)
	bool SetFixedBoundaryCondition_Field( unsigned int id_field, const Field::CFieldWorld& world );

	////////////////////////////////
	// function for update solution

	//! 残差を作る(LinearSystemSaveと宣言を一致させるためのダミーの関数)
	virtual double MakeResidual(const Fem::Field::CFieldWorld& world){ return 0.0; }

	//! fieldの値を更新する
	virtual bool UpdateValueOfField( unsigned int id_field, Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt );
	virtual bool UpdateValueOfField_RotCRV( unsigned int id_field, Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt );
	bool UpdateValueOfField_Newmark(double gamma, double dt, unsigned int id_field_val, Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt, bool IsInitial );
	bool UpdateValueOfField_NewmarkBeta(double gamma, double beta, double dt, unsigned int id_field_val, Field::CFieldWorld& world, bool IsInitial );
	bool UpdateValueOfField_BackwardEular(double dt, unsigned int id_field_val, Field::CFieldWorld& world, bool IsInitial );

	// 存在しないなら-1を返す
	int FindIndexArray_Seg(
		unsigned int id_field, const Field::ELSEG_TYPE& type, const Field::CFieldWorld& world );
    unsigned int GetNLynSysSeg() const { return m_aSegField.size(); }

	////////////////////////////////
	// function for linear solver
	// v=-1:residual    v=-2:update

	//! ソルバに必要な作業ベクトルの数を得る
	virtual unsigned int GetTmpVectorArySize() const { return m_ls.GetTmpVectorArySize(); }
	//! ソルバに必要な作業ベクトルの数を設定
    virtual bool ReSizeTmpVecSolver(unsigned int size_new){ return m_ls.ReSizeTmpVecSolver(size_new); }
    virtual double DOT(int iv1, int iv2){ return m_ls.DOT(iv1,iv2); }//!< ベクトルの内積 (return {v1} * {v2})
    virtual bool COPY(int iv1, int iv2){ return m_ls.COPY(iv1,iv2); } //!< ベクトルのコピー ({v2} := {v1})
    virtual bool SCAL(double alpha, int iv1){ return m_ls.SCAL(alpha,iv1); } //!< ベクトルのスカラー倍 ({v1} := alpha * {v1})
    virtual bool AXPY(double alpha, int iv1, int iv2){ return m_ls.AXPY(alpha,iv1,iv2); }//!< ベクトルの足し算({v2} := alpha*{v1} +　{v2})	
    virtual bool MATVEC(double alpha, int iv1, double beta, int iv2){ return m_ls.MATVEC(alpha,iv1,beta,iv2); } //!< 行列ベクトル積 ({v2} := alpha*[MATRIX]*{v1} + beta*{v2})
protected:
	class CLinSysSeg_Field{
	public:
		unsigned int id_field;	// parent_fieldでなければならない
		unsigned int id_field2;	// combined fieldの時は０以外の値を持つ
		Field::ELSEG_TYPE node_config;
		////////////////
		unsigned int nnode;
		unsigned int len;
	};
protected:
	int AddLinSysSeg( const CLinSysSeg_Field& seg );
	int AddLinSysSeg_Field(unsigned int id_field, Fem::Field::ELSEG_TYPE es_type, const Fem::Field::CFieldWorld& world);
  bool AddMat_NonDia(unsigned int ils_col, unsigned int ils_row, const Com::CIndexedArray& crs ){
    return m_ls.AddMat_NonDia(ils_col,ils_row,crs);
  }
  bool AddMat_Dia(unsigned int ils, Com::CIndexedArray& crs ){
    return m_ls.AddMat_Dia(ils,crs);
  }  
protected:
	std::vector< CLinSysSeg_Field > m_aSegField;
public:
  LsSol::CLinearSystem m_ls;
};


}	// end namespace Ls
}	// end namespace Fem

#endif
