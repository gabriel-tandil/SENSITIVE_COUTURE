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
@brief 複素数有限要素法連立一次方程式クラスのインターフェース
@author Nobuyuki Umetani
*/


#if !defined(ZLINEAR_SYSTEM_H)
#define ZLINEAR_SYSTEM_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <assert.h>

#include "delfem/field.h"

namespace MatVec{
    class CZVector_Blk;
    class CZMatDia_BlkCrs;
    class CDiaMat_Blk;
    class CZMat_BlkCrs;
    class CBCFlag;
}

namespace Fem{
namespace Field{
	class CFieldWorld;
}
namespace Ls{

/*!
@brief 複素数有限要素法連立一次方程式クラス
@ingroup FemLs
*/
class CZLinearSystem
{
public:

	CZLinearSystem();
	virtual ~CZLinearSystem();

	////////////////
	virtual void Clear();

	// fieldで初期化する、fieldの中の非ゼロパターンを作る
	virtual bool AddPattern_Field(const unsigned int id_field, const Field::CFieldWorld& world);

	bool SetVector_fromField(unsigned int iv, unsigned int id_field_val, const Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt);
	bool NormalizeVector(int iv);

	////////////////////////////////
	// function for marge element matrix

	// 残差ベクトルをゲットする
    MatVec::CZVector_Blk* GetResidualPtr(unsigned int id_field, const Field::ELSEG_TYPE& elseg_type, const Field::CFieldWorld& world);
	// 更新ベクトルをゲットする
    MatVec::CZVector_Blk* GetUpdatePtr(unsigned int id_field, const Field::ELSEG_TYPE& elseg_type, const Field::CFieldWorld& world);
	// 対角行列をゲットする
    MatVec::CZMatDia_BlkCrs* GetMatrixPtr(unsigned int id_field, const Field::ELSEG_TYPE& elseg_type, const Field::CFieldWorld& world);
	// 非対角行列をゲットする
    MatVec::CZMat_BlkCrs* GetMatrixPtr(unsigned int id_field_col,const Field::ELSEG_TYPE& elseg_type_col,
		unsigned int id_field_row,const Field::ELSEG_TYPE& elseg_type_row,
		const Field::CFieldWorld& world);
	

	////////////////////////////////
	// function for marge

	// マージに必要なバッファサイズを得る
	unsigned int GetTmpBufferSize();
	// マージ前の初期化
	virtual void InitializeMarge();
	// マージ後の処理（残差ノルムを返す)
	virtual double FinalizeMarge(); 

	////////////////////////////////
	// function for fixed boundary condition

	// 固定境界条件を全解除
	void ClearFixedBoundaryCondition();
	// 固定境界条件の設定
	// idof:固定する自由度
	bool SetFixedBoundaryCondition_Field( unsigned int id_field, unsigned int idofns, const Field::CFieldWorld& world );
	// 固定境界条件の設定---field中の全自由度固定
	bool SetFixedBoundaryCondition_Field( unsigned int id_field, const Field::CFieldWorld& world );

	////////////////////////////////
	// function for update solution

	// 残差を作る
	virtual double MakeResidual(const unsigned int id_field_val, const Fem::Field::CFieldWorld& world){ return 0.0; }
	virtual bool UpdateValueOfField( unsigned int id_field, Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt );

	////////////////////////////////
	// function for linear solver
	// v=-1:residual    v=-2:update

	// ソルバに必要な作業ベクトルの数を得る
	unsigned int GetTmpVectorArySize() const { return m_TmpVectorArray.size(); }
	// ソルバに必要な作業ベクトルの数を設定
	bool ReSizeTmpVecSolver(unsigned int size_new);
	// ベクトルの共役を取る
	bool Conjugate(int iv1);
	// ベクトルの掛け算
	// return {v1} * {v2}
	Com::Complex DOT(int iv1, int iv2);
	// ベクトルの内積
	// return {v1} * {v2}^H
	Com::Complex INPROCT(int iv1, int iv2);
	// ベクトルのコピー
	// {v2} := {v1}
	bool COPY(int iv_from, int iv_to);
	// ベクトルのスカラー倍
	// {v1} := alpha * {v1}
	bool SCAL(Com::Complex alpha, int iv1);
	// ベクトルの足し算
	// {v2} := alpha*{v1} +　{v2}
	bool AXPY(Com::Complex alpha, int iv1, int iv2);
	// 行列ベクトル積
	// {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
	bool MatVec(double alpha, int iv1, double beta, int iv2);
	// 行列のエルミートとベクトル積
	// {v2} := alpha*[MATRIX]^H*{v1} + beta*{v2}
	bool MatVec_Hermitian(double alpha, int iv1, double beta, int iv2);

	////////////////////////////////
	// function for preconditioner
	
	unsigned int GetNLynSysSeg() const { return this->m_aSeg.size(); }
	// 存在しないなら-1を返す
	int FindIndexArray_Seg( unsigned int id_field, const Field::ELSEG_TYPE& type, const Field::CFieldWorld& world );

	bool IsMatrix(unsigned int ilss, unsigned int jlss) const {
		const unsigned int nlss = m_aSeg.size();
		assert( ilss < nlss );
		assert( jlss < nlss );
		if( ilss == jlss ) return true;
		if( this->m_Matrix_NonDia[ilss][jlss] == 0 ) return false;
		return true;
	}
    const MatVec::CZMatDia_BlkCrs& GetMatrix(unsigned int ilss) const {
		const unsigned int nlss = m_aSeg.size();
		assert( ilss < nlss );
		assert( m_Matrix_Dia[ilss] != 0 );
		return *m_Matrix_Dia[ilss];
	}
    const MatVec::CZMat_BlkCrs& GetMatrix(unsigned int ilss, unsigned int jlss) const {
		const unsigned int nlss = m_aSeg.size();
		assert( ilss < nlss );
		assert( jlss < nlss );
		assert( m_Matrix_NonDia[ilss][jlss] != 0 );
		return *m_Matrix_NonDia[ilss][jlss];
	}
    MatVec::CZVector_Blk& GetVector(int iv, unsigned int ilss){
		assert( iv < (int)this->m_TmpVectorArray.size() );
		if( iv == -1 ){
			return *m_Residual[ilss];
		}
		else if( iv == -2 ){
			return *m_Update[ilss];
		}
		const unsigned int nlss = m_aSeg.size();
		assert( ilss < nlss );
		return *m_TmpVectorArray[iv][ilss];
	}
private:
	class CLinSysSeg{
	public:
		unsigned int id_field;	// parent_fieldでなければならない
		Field::ELSEG_TYPE node_config;
		unsigned int nnode;
		unsigned int len;
	};
protected:
	int AddLinSysSeg( const CLinSysSeg& seg );
    bool AddMat_NonDia(unsigned int ils_col, unsigned int ils_row, const Com::CIndexedArray& crs );
	bool AddMat_Dia(unsigned int ils, const Field::CElemAry& ea, unsigned int id_es);
protected:
	std::vector< CLinSysSeg > m_aSeg;
    std::vector< std::vector< MatVec::CZMat_BlkCrs* > > m_Matrix_NonDia;
    std::vector< MatVec::CZMatDia_BlkCrs* > m_Matrix_Dia;
    std::vector< MatVec::CZVector_Blk* > m_Residual, m_Update;
    std::vector< MatVec::CBCFlag* > m_BCFlag;	// Boundary Condition Flag
    std::vector< std::vector< MatVec::CZVector_Blk* > > m_TmpVectorArray;	// Working Buffer for Linear Solver
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CZLinearSystem_GeneralEigen : public CZLinearSystem
{
public:
	virtual void Clear();
	// fieldで初期化する、fieldの中の非ゼロパターンを作る
	virtual bool AddPattern_Field(const unsigned int id_field, const Field::CFieldWorld& world);
	virtual void InitializeMarge();
    MatVec::CDiaMat_Blk* GetDiaMassMatrixPtr(unsigned int id_field, const Field::ELSEG_TYPE& elseg_type, const Field::CFieldWorld& world);
    /*!
    @brief 場の値をベクトルに設定
    @param[in] iv 0以上なら作業ベクトル, -1なら残差ベクトル, -2なら更新ベクトル
    */
    bool SetVector_fromField(int iv, unsigned int id_field_val, const Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt );
	bool DecompMultMassMatrix();
	void OffsetDiagonal(double lambda);
	bool MultUpdateInvMassDecomp();
	bool MultVecMassDecomp(int ivec);
	void RemoveConstant(int iv);
private:
    std::vector< MatVec::CDiaMat_Blk* > m_DiaMassMatrix;
};


}	// Ls
}	// Fem

#endif
