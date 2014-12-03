/*
Copyright (C) 2009  Nobuyuki Umetani    n.umetani@gmail.com

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distriｑbuted in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*! @file
@brief 連立一次方程式クラス(LsSol::CLinearSystem)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(LINEAR_SYSTEM_H)
#define LINEAR_SYSTEM_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <assert.h>

#include "delfem/ls/linearsystem_interface_solver.h"
#include "delfem/indexed_array.h"

namespace MatVec{
class CVector_Blk;
class CZVector_Blk;
class CMat_BlkCrs;
class CMatDia_BlkCrs;
class CDiaMat_Blk;
class CBCFlag;
}

namespace LsSol{

/*! 
@brief 連立一次方程式クラス
@ingroup LsSol
*/
class CLinearSystem : public LsSol::ILinearSystem_Sol
{
public:
    //! default constructer
    CLinearSystem(){}
    //! destructor
    virtual ~CLinearSystem(){ this->Clear(); }

	////////////////
	virtual void Clear();	//!< 全てのデータのクリア(行列のサイズは保持，非ゼロパターンは消去)
    void ClearFixedBoundaryCondition();

	////////////////////////////////
	// function for marge
	//! マージ前の初期化(行列，ベクトルに０を設定する)
	virtual void InitializeMarge();
	//! マージ後の処理（境界条件を設定し，残差ノルムを返す)
	virtual double FinalizeMarge(); 

	////////////////////////////////
	// function for linear solver
	// v=-1:residual    v=-2:update

	//! ソルバに必要な作業ベクトルの数を得る
	virtual unsigned int GetTmpVectorArySize() const { return m_TmpVectorArray.size(); }
	//! ソルバに必要な作業ベクトルの数を設定
	virtual bool ReSizeTmpVecSolver(unsigned int size_new);

	virtual double DOT(int iv1, int iv2); //!< ベクトルの内積 (return {v1} * {v2})
	virtual bool COPY(int iv1, int iv2); //!< ベクトルのコピー ({v2} := {v1})
	virtual bool SCAL(double alpha, int iv1); //!< ベクトルのスカラー倍 ({v1} := alpha * {v1})
	virtual bool AXPY(double alpha, int iv1, int iv2); //!< ベクトルの足し算({v2} := alpha*{v1} +　{v2})	
	virtual bool MATVEC(double alpha, int iv1, double beta, int iv2); //!< 行列ベクトル積 ({v2} := alpha*[MATRIX]*{v1} + beta*{v2})

	////////////////////////////////
	// function for preconditioner

	int AddLinSysSeg( unsigned int nnode, unsigned int len );
    int AddLinSysSeg( unsigned int nnode, const std::vector<unsigned int>& aLen );
	unsigned int GetNLinSysSeg() const { return this->m_Matrix_Dia.size(); }
    const unsigned int GetBlkSizeLinSeg(unsigned int ilss) const{
        assert( ilss < this->GetNLinSysSeg() );
        return m_aSeg[ilss].nnode;
    }

	bool IsMatrix(unsigned int ilss, unsigned int jlss) const {
		const unsigned int nlss = m_Matrix_Dia.size();
		assert( ilss < nlss );
		assert( jlss < nlss );
        if( ilss == jlss ){
            if( this->m_Matrix_Dia[ilss] == 0 ) return false;
            return true;
        }
		if( this->m_Matrix_NonDia[ilss][jlss] == 0 ) return false;
		return true;
	}
    const MatVec::CMatDia_BlkCrs& GetMatrix(unsigned int ilss) const {
		const unsigned int nlss = m_Matrix_Dia.size();
		assert( ilss < nlss );
		assert( m_Matrix_Dia[ilss] != 0 );
		return *m_Matrix_Dia[ilss];
	}
    MatVec::CMatDia_BlkCrs& GetMatrix(unsigned int ilss){
		const unsigned int nlss = m_Matrix_Dia.size();
		assert( ilss < nlss );
		assert( m_Matrix_Dia[ilss] != 0 );
		return *m_Matrix_Dia[ilss];
	}
    const MatVec::CMat_BlkCrs& GetMatrix(unsigned int ilss, unsigned int jlss) const {
		const unsigned int nlss = this->GetNLinSysSeg();
		assert( ilss < nlss );
		assert( jlss < nlss );
		assert( m_Matrix_NonDia[ilss][jlss] != 0 );
		return *m_Matrix_NonDia[ilss][jlss];
	}
    MatVec::CMat_BlkCrs& GetMatrix(unsigned int ilss, unsigned int jlss){
		const unsigned int nlss = this->GetNLinSysSeg();
		assert( ilss < nlss );
		assert( jlss < nlss );
		assert( m_Matrix_NonDia[ilss][jlss] != 0 );
		return *m_Matrix_NonDia[ilss][jlss];
	}
    MatVec::CVector_Blk& GetVector(int iv, unsigned int ilss){
		assert( iv < (int)this->m_TmpVectorArray.size() );
        assert( ilss < this->GetNLinSysSeg() );
        if(      iv == -1 ){ return *m_Residual[ilss]; }
        else if( iv == -2 ){ return *m_Update[  ilss]; }
		return *m_TmpVectorArray[iv][ilss];
	}
    const MatVec::CVector_Blk& GetVector(int iv, unsigned int ilss) const {
		assert( iv < (int)this->m_TmpVectorArray.size() );
        assert( ilss < this->GetNLinSysSeg() );
        if(      iv == -1 ){ return *m_Residual[ilss]; }
        else if( iv == -2 ){ return *m_Update[ilss];   }
		return *m_TmpVectorArray[iv][ilss];
	}
    MatVec::CBCFlag& GetBCFlag(unsigned int ilss){
		assert( ilss < this->GetNLinSysSeg() );
        return *m_BCFlag[ilss];
    }
public:
    bool AddMat_NonDia(unsigned int ils_col, unsigned int ils_row, const Com::CIndexedArray& crs );
	bool AddMat_Dia(unsigned int ils, const Com::CIndexedArray& crs );
private:
	class CLinSysSeg{
    public:
        CLinSysSeg(unsigned int n, unsigned int l): nnode(n), len(l){}
        CLinSysSeg(unsigned int n, std::vector<unsigned int> al): nnode(n), len(-1), aLen(al){
            assert( al.size() == n );
        }
	public:
		unsigned int nnode;
		int len;
        std::vector<unsigned int> aLen;
	};
	std::vector< CLinSysSeg > m_aSeg;
public:
    std::vector< std::vector< MatVec::CMat_BlkCrs* > > m_Matrix_NonDia;
    std::vector< MatVec::CMatDia_BlkCrs* > m_Matrix_Dia;
    std::vector< MatVec::CVector_Blk* > m_Residual, m_Update;
private:
    std::vector< std::vector< MatVec::CVector_Blk* > > m_TmpVectorArray;	// Working Buffer for Linear Solver
    std::vector< MatVec::CBCFlag* > m_BCFlag;	// Boundary Condition Flag
};

}

#endif

