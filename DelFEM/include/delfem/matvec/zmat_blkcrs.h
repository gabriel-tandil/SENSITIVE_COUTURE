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
@brief 複素数ブロックCRS行列クラス(MatVec::CZMat_BlkCrs)クラスのインターフェース
@author Nobuyuki Umetani
*/

#if !defined(ZMAT_BLK_CRS_H)
#define ZMAT_BLK_CRS_H

#include <assert.h>

#include "delfem/complex.h"

namespace Com{	
    class CIndexedArray;
}

namespace MatVec{

class CBCFlag;
class CZVector_Blk;
/*! 
@brief 複素数ブロックCRS行列クラス
@ingroup MatVec
*/
class CZMat_BlkCrs  
{
public:
	CZMat_BlkCrs(const unsigned int nblk_col, const unsigned int len_col, 
		const unsigned int nblk_row, const unsigned int len_row );
	CZMat_BlkCrs(const CZMat_BlkCrs& rhs, const bool is_value, const bool isnt_trans, const bool isnt_conj);
	CZMat_BlkCrs(); // Default Constructor
	virtual ~CZMat_BlkCrs();

	void Initialize(const unsigned int nblk_col, const unsigned int len_col, 
		const unsigned int nblk_row, const unsigned int len_row );

	unsigned int NBlkMatCol() const { return m_nblk_MatCol; }
	unsigned int NBlkMatRow() const { return m_nblk_MatRow; }
	
	unsigned int LenBlkCol() const { return m_len_BlkCol; }
	unsigned int LenBlkRow() const { return m_len_BlkRow; }

	virtual bool SetZero();
	virtual bool DeletePattern();
	
    virtual bool AddPattern(const Com::CIndexedArray& crs);
    virtual bool AddPattern(const MatVec::CZMat_BlkCrs& rhs, const bool isnt_trans);

    virtual bool SetValue(const MatVec::CZMat_BlkCrs& rhs, const bool isnt_trans, bool isnt_conj);
	
    virtual bool MatVec(double alpha, const MatVec::CZVector_Blk& x, double beta, MatVec::CZVector_Blk& b, const bool isnt_trans) const;

	bool SetBoundaryCondition_Row(const MatVec::CBCFlag& bc_flag);
	bool SetBoundaryCondition_Colum(const MatVec::CBCFlag& bc_flag);
	
	////////////////////////////////////////////////
	// Crs Original Function ( Non-virtual )
	////////////////////////////////////////////////
	const unsigned int NCrs() const { return m_ncrs_Blk; }
	const unsigned int* GetPtrIndPSuP(const unsigned int ipoin, unsigned int& npsup) const
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_rowPtr_Blk[icrs];
	}
	const Com::Complex* GetPtrValPSuP(const unsigned int ipoin, unsigned int& npsup) const
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_valCrs_Blk[icrs];
	}
	Com::Complex* GetPtrValPSuP(const unsigned int ipoin, unsigned int& npsup)
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_valCrs_Blk[icrs];
	}

protected:
	unsigned int m_nblk_MatCol;
	unsigned int m_nblk_MatRow;

	unsigned int  m_ncrs_Blk;
	unsigned int* m_colInd_Blk;
	unsigned int* m_rowPtr_Blk;

	unsigned int m_len_BlkCol;
	unsigned int m_len_BlkRow;
	Com::Complex* m_valCrs_Blk;
};

}	// namespace LS

#endif 
