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
@brief 複素数前処理クラスのインターフェース
@author Nobuyuki Umetani
*/

#if !defined(ZPRECONDITIONER_H)
#define ZPRECONDITIONER_H

#include <assert.h>

#include "delfem/femls/zlinearsystem.h"

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/zmatdiafrac_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/bcflag_blk.h"

namespace Fem{
namespace Ls{

/*!
@brief 複素数抽象前処理クラス
@ingroup FemLs
*/
class CZPreconditioner
{
public:
	virtual ~CZPreconditioner(){};
	virtual void SetLinearSystem(const CZLinearSystem& ls) = 0;
	virtual bool SetValue(const CZLinearSystem& ls) = 0;
	virtual bool SolvePrecond(CZLinearSystem& ls, unsigned int iv) = 0;
};

/*!
@brief 複素数ILU前処理クラス
@ingroup FemLs
*/
class CZPreconditioner_ILU : public CZPreconditioner
{
public:
	CZPreconditioner_ILU(){
		m_nlev = 0;
	}
	CZPreconditioner_ILU(const CZLinearSystem& ls, unsigned int nlev = 0){ 
		m_nlev = nlev;
		this->SetLinearSystem(ls); 
	}
	virtual ~CZPreconditioner_ILU(){
		this->Clear();
	}
	void Clear(){
		for(unsigned int i=0;i<m_Matrix_Dia.size();i++){
			if( m_Matrix_Dia[i] != 0 ) delete m_Matrix_Dia[i];
		}
		m_Matrix_Dia.clear();
	}
	void SetFillInLevel(unsigned int nlev){ m_nlev = nlev; }
	// ILU(0)のパターン初期化
	virtual void SetLinearSystem(const CZLinearSystem& ls)
	{
		const unsigned int nlss = ls.GetNLynSysSeg();
		assert( nlss == 1 );
		m_Matrix_Dia.clear();
        m_Matrix_Dia.push_back( new MatVec::CZMatDiaFrac_BlkCrs(m_nlev,ls.GetMatrix(0)) );
	}

	// 値を設定してILU分解を行う関数
	// ILU分解が成功しているかどうかはもっと詳細なデータを返したい
	virtual bool SetValue(const CZLinearSystem& ls){
		unsigned int nlss = ls.GetNLynSysSeg();
		assert( nlss==1 );
		assert( m_Matrix_Dia.size() == 1 );
		////////////////
		// 値をコピー
		m_Matrix_Dia[0]->SetValue_Initialize( ls.GetMatrix(0) );
		////////////////
		// ILU分解
		m_Matrix_Dia[0]->DoILUDecomp();
		return true;
	}

	// Solve Preconditioning System
	virtual bool SolvePrecond(CZLinearSystem& ls, unsigned int iv){
		const int nlss = ls.GetNLynSysSeg();
		assert( nlss == 1 );
		m_Matrix_Dia[0]->ForwardSubstitution( ls.GetVector(iv,0));
		m_Matrix_Dia[0]->BackwardSubstitution(ls.GetVector(iv,0));
		return true;
	}
private:
	unsigned int m_nlev;
//	std::vector< std::vector< MatVec::CMatFrac_BlkCrs* > > m_Matrix_NonDia;
    std::vector< MatVec::CZMatDiaFrac_BlkCrs* > m_Matrix_Dia;
};

}	// end namespace Ls
}	// end namespace Fem

#endif
