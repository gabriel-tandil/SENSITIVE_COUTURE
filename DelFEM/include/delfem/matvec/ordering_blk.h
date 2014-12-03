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
@brief interface of node ordering class (MatVec::COrdering_Blk)
@author Nobuyuki Umetani
*/

#if !defined(ORDERING_BLK_H)
#define ORDERING_BLK_H

#include <assert.h>
#include <vector>

namespace MatVec{

class CMatDia_BlkCrs;
class CVector_Blk;

/*! 
@brief node ordering class
@ingroup MatVec
*/
class COrdering_Blk{
public:
	COrdering_Blk(){
		m_pOrder = 0;
		m_pInvOrder = 0;
	}
	~COrdering_Blk(){
		if( m_pOrder != 0 ) delete[] m_pOrder;
		if( m_pInvOrder != 0 ) delete[] m_pInvOrder;
	}
  void SetOrdering(const std::vector<int>& ord);
	void MakeOrdering_RCM(const CMatDia_BlkCrs& mat);
	void MakeOrdering_RCM2(const CMatDia_BlkCrs& mat);
	void MakeOrdering_AMD(const CMatDia_BlkCrs& mat);
	unsigned int NBlk() const { return m_nblk; }
	int NewToOld(unsigned int iblk_new) const { return m_pOrder[iblk_new]; }
	int OldToNew(unsigned int iblk_old) const { return m_pInvOrder[iblk_old]; }
	void OrderingVector_NewToOld(CVector_Blk& vec_to, const CVector_Blk& vec_from);
	void OrderingVector_OldToNew(CVector_Blk& vec_to, const CVector_Blk& vec_from);
private:
	unsigned int m_nblk;
	int* m_pOrder;
	int* m_pInvOrder;
};

}

#endif
