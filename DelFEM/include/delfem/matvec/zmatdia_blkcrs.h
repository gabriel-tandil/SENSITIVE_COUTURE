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
@brief 複素数対角ブロックCRS行列クラス(MatVec::CZMatDia_BlkCrs)クラスのインターフェース
@author Nobuyuki Umetani
*/

#if !defined(ZMAT_DIA_BLK_CRS_H)
#define ZMAT_DIA_BLK_CRS_H

#include <string>

#include "delfem/complex.h"
#include "delfem/matvec/bcflag_blk.h"
#include "delfem/matvec/zmat_blkcrs.h"

namespace MatVec{

class CZVector_Blk;
/*! 
@brief 複素数対角ブロックCRS行列クラス
@ingroup MatVec
*/
class CZMatDia_BlkCrs  : public CZMat_BlkCrs
{
public:
	CZMatDia_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow);
	CZMatDia_BlkCrs(const std::string& file_path);
	CZMatDia_BlkCrs(const CZMatDia_BlkCrs& rhs, const bool is_value, const bool isnt_trans, const bool isnt_conj);
	virtual ~CZMatDia_BlkCrs();

	virtual bool DeletePattern();

    virtual bool AddPattern(const Com::CIndexedArray& crs);
	bool AddPattern(const CZMatDia_BlkCrs& rhs, const bool isnt_trans);
	bool AddPattern(const CZMat_BlkCrs& m1, const CZMatDia_BlkCrs& m2, const CZMat_BlkCrs& m3);

	bool SetValue(const CZMatDia_BlkCrs& rhs, const bool isnt_trans);
	bool SetValue(const CZMat_BlkCrs& m1, const CZMatDia_BlkCrs& m2, const CZMat_BlkCrs& m3);

	virtual bool SetZero();
	virtual bool Mearge(
		unsigned int nblkel_col, const unsigned int* blkel_col,
		unsigned int nblkel_row, const unsigned int* blkel_row,
		unsigned int blksize, const Com::Complex* emat, int* tmp_buffer);
	void AddUnitMatrix(const Com::Complex& epsilon);

	bool SetBoundaryCondition(const CBCFlag& bc_flag);

	bool MatVec(double alpha, const CZVector_Blk& rhs, double beta, CZVector_Blk& lhs) const;
	bool MatVec_Hermitian(double alpha, const CZVector_Blk& rhs, double beta, CZVector_Blk& lhs) const;

	////////////////////////////////

	const Com::Complex* GetPtrValDia(const unsigned int ipoin) const	{ return &m_valDia_Blk[ipoin]; }
	Com::Complex* GetPtrValDia(const unsigned int ipoin){ return &m_valDia_Blk[ipoin]; }

protected:
	Com::Complex* m_valDia_Blk;
};

}	// namespace Ls

#endif // MATDIA_CRS_H
