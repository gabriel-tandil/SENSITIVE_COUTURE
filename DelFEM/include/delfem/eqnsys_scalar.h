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
@brief スカラー型の方程式クラス(Fem::Eqn::CEqnSystem_Scalar2D, Fem::Eqn::CEqn_Scalar2D, Fem::Eqn::CEqn_Scalar3D)のインターフェース
@author Nobyuki Umetani

@li ポアソン方程式
@li 拡散方程式
@li 静的移流拡散方程式
@li 動的移流拡散方程式
*/



#if !defined(EQN_OBJ_SCALAR_H)
#define EQN_OBJ_SCALAR_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include "delfem/eqnsys.h"

namespace Fem{

namespace Ls{
	class CLinearSystem_Field;
	class CLinearSystem_SaveDiaM_Newmark;
	class CPreconditioner;
}

namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn{

/*! 
@brief 2Dのスカラー型方程式クラス
@ingroup FemEqnObj
*/
class CEqn_Scalar2D
{
public:
	CEqn_Scalar2D(unsigned int id_ea, unsigned int id_field_val)
		: m_id_ea(id_ea), m_IdFieldVal(id_field_val)
	{
		m_alpha = 1.0; 
		m_capa = 1.0;
		m_source = 0.0;
		m_IdFieldAdvec = 0;
	}
	// Setメソッド
	//! @{
	void SetAlpha(   double alpha ){ m_alpha  = alpha;  }
	void SetCapacity(double capa  ){ m_capa   = capa;   }
	void SetSource(  double source){ m_source = source; }
	void SetAdvection(unsigned int id_field_advec){ m_IdFieldAdvec = id_field_advec; }
	//! @}
        double GetAlpha() const { return m_alpha; }
        double GetCapacity() const { return m_capa; }
        double GetSource() const { return m_source; }
	// Getメソッド
	unsigned int GetIdEA() const { return m_id_ea; }
	bool IsAdvection()  const { return m_IdFieldAdvec != 0; }

	// 連立一次方程式マージメソッド
	bool AddLinSys( Ls::CLinearSystem_Field& ls, const Field::CFieldWorld& world );
	bool AddLinSys_Newmark( double dt, double gamma, Ls::CLinearSystem_Field& ls, 
		bool is_ax_sym,
		const Field::CFieldWorld& world );
	bool AddLinSys_Save( Ls::CLinearSystem_Save& ls, const Field::CFieldWorld& world );
	bool AddLinSys_SaveKDiaC( Ls::CLinearSystem_SaveDiaM_Newmark& ls, const Field::CFieldWorld& world );
private:
	unsigned int m_id_ea;
	unsigned int m_IdFieldVal;
	unsigned int m_IdFieldAdvec;

	double m_alpha;
	double m_capa;
	double m_source;
};
	
////////////////////////////////////////////////////////////////

/*! 
@ingroup FemEqnSystem
@brief 2Dのスカラー型,連成方程式クラス

解くことができる方程式
@li ポアソン方程式
@li 拡散方程式
@li 静的移流拡散方程式
@li 動的移流拡散方程式
*/
class CEqnSystem_Scalar2D : public CEqnSystem
{
public:
	//! デフォルト・コンストラクタの設定
	CEqnSystem_Scalar2D();
	CEqnSystem_Scalar2D(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
	// vritual 関数
	virtual bool SetDomain_Field(unsigned int id_field, Fem::Field::CFieldWorld& world);
	virtual bool SetDomain_FieldElemAry(unsigned int id_field, unsigned int id_ea, Fem::Field::CFieldWorld& world);

	//! 方程式を解く
	virtual bool Solve(Fem::Field::CFieldWorld& world);
	
	////////////////////////////////
	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	virtual unsigned int GetIdField_Value() const { return m_IdFieldVal; }
	
	//! 方程式を設定する
	bool SetEquation( const CEqn_Scalar2D& eqn );
	//! 方程式を取得する
   CEqn_Scalar2D GetEquation(unsigned int id_ea) const;

	//! 拡散係数の設定
	void SetAlpha(double alpha);
	//! 容量の設定
	void SetCapacity(double capa);
	//! 面積に比例するソース項の設定
	void SetSource(double source);
	/*! 
	@brief 静的問題
	@param [in] is_stat 静的問題かどうか
	*/
	void SetStationary(bool is_stat);	
	/*! 
	@brief 移流項の設定
	@param [in] id_field_advec 値を移流させる速度場のID
	@remark もしもid_field_advecが０または存在しないIDの場合は移流項は無視される
	*/
	void SetAdvection(unsigned int id_field_advec);
	void SetAxialSymmetry(bool is_ax_sym);
	
	void SetSaveStiffMat(bool is_save){
		m_IsSaveStiffMat = is_save;
		this->ClearLinearSystemPreconditioner();
	}

private:
	virtual double MakeLinearSystem(const Fem::Field::CFieldWorld& world);
	virtual bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
	virtual bool MatrixProperty(bool& is_c, bool& is_m, bool& is_asym );
private:
	unsigned int m_IdFieldVal;
	////////////////
	// 境界条件について
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

	std::vector<CEqn_Scalar2D> m_aEqn;
	bool m_IsAxialSymmetry;
	bool m_IsStationary;
	bool m_IsSaveStiffMat;
};

////////////////

/*! 
@brief 3Dのスカラー型方程式クラス
@ingroup FemEqnObj
*/
class CEqn_Scalar3D : public CEqnSystem
{
public:
	//! デフォルト・コンストラクタ
	CEqn_Scalar3D();
	CEqn_Scalar3D(const unsigned int id_field_val, Fem::Field::CFieldWorld& world);
	
	virtual bool SetDomain(unsigned int id_base, Fem::Field::CFieldWorld& world);
	//! 解く
	virtual bool Solve(Fem::Field::CFieldWorld& world);

	////////////////////////////////
	// 固定境界条件を追加&削除する
	virtual bool         AddFixField(   unsigned int id_field,                  Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( unsigned int id_ea,                     Fem::Field::CFieldWorld& world, int idof = -1);
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1);
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world );
	virtual void ClearFixElemAry();

	virtual unsigned int GetIdField_Value() const { return m_IdFieldVal; }

	//! 拡散係数を設定する
	void SetAlpha(double alpha){ 
		m_alpha = alpha; 
		this->m_is_cleared_value_ls   = true;
		this->m_is_cleared_value_prec = true;
	}
	//! 体積に比例するソース項を設定する
	void SetSource(double source){ 
		m_source = source; 
		this->m_is_cleared_value_ls   = true;
		this->m_is_cleared_value_prec = true;
	}
private:
	virtual double MakeLinearSystem(const Fem::Field::CFieldWorld& world);
	virtual bool InitializeLinearSystem(const Fem::Field::CFieldWorld& world);
private:
	unsigned int m_IdFieldVal;
	////////////////
	// 境界条件について
	std::vector< std::pair<unsigned int,int> > m_aIdFixField;

	double m_alpha;
	double m_source;
};

}
}

#endif
