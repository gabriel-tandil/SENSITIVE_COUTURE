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
@brief 複素数ILU分解対角ブロックCRS行列クラス(MatVec::CZMatDiaFrac_BlkCrs)クラスのインターフェース
@author Nobuyuki Umetani
*/

#if !defined(ZMAT_DIA_CRS_FRAC)
#define ZMAT_DIA_CRS_FRAC

#include <vector>

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/zmatprecond_blk.h"

namespace MatVec{

/*! 
@brief 複素数ILU分解対角ブロックCRS行列クラス
@ingroup MatVec
*/
class CZMatDiaFrac_BlkCrs : public MatVec::CZMatDia_BlkCrs, public MatVec::CZMatPrecond_Blk
{
private:
	struct SRowLevNext{
		unsigned int row;
		unsigned int lev;
		int next;
	};
	struct SRowLev{
		unsigned int row;
		unsigned int lev;
	};
public:
	CZMatDiaFrac_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow);
	CZMatDiaFrac_BlkCrs(const CZMatDia_BlkCrs& rhs);
	CZMatDiaFrac_BlkCrs(const unsigned int lev_fill, const CZMatDia_BlkCrs& rhs);
	virtual ~CZMatDiaFrac_BlkCrs();

	bool SetValue_Initialize(const CZMatDia_BlkCrs& rhs){
		CZMatDia_BlkCrs::SetValue(rhs,true);
		return true;
	}

	//　この関数はILU分解もやってしまう
	virtual bool SetValue(const CZMatDia_BlkCrs& rhs);

	bool AddFracPtn(const int lev_fill);
	
	bool FracInitialize(const CZMatDia_BlkCrs& rhs);
	bool FracFinalize();
	
	bool Solve(CZVector_Blk& vec) const{
		this->ForwardSubstitution(vec);
		this->BackwardSubstitution(vec);
		return true;
	}

	bool ForwardSubstitution(CZVector_Blk& vec) const;
	bool BackwardSubstitution(CZVector_Blk& vec) const;

	virtual bool SolvePrecond(const CZMatDia_BlkCrs& mat, CZVector_Blk& vec) const{
		return this->Solve(vec);
	}

	bool DoILUDecomp();
private:

private:
	// 0 for not pattern 
	// 1 for during pattern
	// 2 for after pattern
	unsigned int m_ConditionFlag;

	unsigned int* m_DiaInd;
	std::vector<SRowLev>* m_pRowLev;
};

}	// namespace Ls

#endif 
