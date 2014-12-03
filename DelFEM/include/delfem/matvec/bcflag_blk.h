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
@brief 境界条件フラグクラス(MatVec::CBCFlag)のインターフェイス
*/


#if !defined(BC_FLAG_H)
#define BC_FLAG_H

#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/zvector_blk.h"

namespace MatVec{
/*!
@brief 境界条件フラグ
@ingroup MatVec
*/
class CBCFlag
{
public:
	CBCFlag(unsigned int nBlk, unsigned int lenBlk)
  : m_lenBlk(lenBlk), m_nBlk(nBlk){
    m_DofPtr = 0;
    m_Flag = new int [m_lenBlk*m_nBlk];
      this->SetAllFlagZero();
	}
  CBCFlag(unsigned int nBlk, const std::vector<unsigned int>& aLen )
  : m_lenBlk(-1), m_nBlk(nBlk)
  {
        m_DofPtr = new unsigned int [m_nBlk+1];
        m_DofPtr[0] = 0;
        for(unsigned int iblk=0;iblk<m_nBlk;iblk++){
            m_DofPtr[iblk+1] = m_DofPtr[iblk] + aLen[iblk];
        }
        const unsigned int ni = m_DofPtr[m_nBlk];
		m_Flag = new int [ni];
		this->SetAllFlagZero();
	}
	virtual ~CBCFlag(){
        if( m_DofPtr != 0 ){ delete[] m_DofPtr; }
		if( m_Flag   != 0 ){ delete[] m_Flag;   }
	}

	////////////////
	int LenBlk() const { return m_lenBlk; }
	unsigned int LenBlk(unsigned int iblk) const { 
        if( m_lenBlk != -1 ){
            return m_lenBlk; 
        }
        assert( m_DofPtr != 0 );
        return m_DofPtr[iblk+1] - m_DofPtr[iblk];
    }
	unsigned int NBlk() const { return m_nBlk; }
	int GetBCFlag(unsigned int iblk, unsigned int idofblk) const{
		assert( iblk < m_nBlk );
        if( m_lenBlk != -1 ){
            assert( m_DofPtr == 0 );
            assert( (int)idofblk < m_lenBlk );
		    return m_Flag[iblk*m_lenBlk+idofblk];
        }
        assert( m_DofPtr != 0 );
        assert( idofblk < m_DofPtr[iblk+1] - m_DofPtr[iblk] );
        const unsigned int ipos = m_DofPtr[iblk];
        return m_Flag[ipos+idofblk];
	}

	////////////////

	// フラグを全てリセット
	void SetAllFlagZero(){
        unsigned int ndof;
        if( m_lenBlk == -1 ){
            assert( m_DofPtr != 0 );
            ndof = m_DofPtr[m_nBlk];
        }
        else{ 
            assert( m_DofPtr == 0 );
            ndof = m_lenBlk*m_nBlk; 
        }
		for(unsigned int idof=0;idof<ndof;idof++){ m_Flag[idof] = 0; }
	}
	
	void SetZeroToBCDof(CVector_Blk& vec){
		assert( vec.NBlk() == this->m_nBlk );
		assert( vec.Len() == this->m_lenBlk );
        const unsigned int nblk = m_nBlk;
        if( this->m_lenBlk == -1 ){
            assert( m_DofPtr != 0 );
		    for(unsigned int iblk=0;iblk<nblk;iblk++){
                const unsigned int ipos0 = m_DofPtr[iblk];
                const unsigned int len_i = m_DofPtr[iblk+1] - m_DofPtr[iblk];
                assert( len_i == vec.Len(iblk) );
		        for(unsigned int ilen=0;ilen<len_i;ilen++){
		            if( m_Flag[ipos0+ilen] == 0 ) continue;
			        vec.SetValue(iblk,ilen,0);
			    }
		    }
            return;
        }
        ////////////////
        assert( m_DofPtr == 0 );
        const unsigned int nlen = m_lenBlk;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
		for(unsigned int ilen=0;ilen<nlen;ilen++){
		    if( m_Flag[iblk*nlen+ilen] == 1 ){
			    vec.SetValue(iblk,ilen,0);
			}
		}
        }
	}
    void SetZeroToBCDof(MatVec::CZVector_Blk& vec)
    {
        assert( this->m_lenBlk >= 0 );
        assert( vec.BlkVecLen() == m_nBlk );
        assert( vec.BlkLen() == (unsigned int)this->m_lenBlk );
        const unsigned int nblk = m_nBlk;
        if( this->m_lenBlk == -1 ){
            assert(0);
            return;
        }
        ////////////////
        assert( m_DofPtr == 0 );
        const unsigned int nlen = m_lenBlk;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			if( m_Flag[iblk*nlen+ilen] == 1 ){
				vec.SetValue(iblk,ilen,0);
			}
		}
		}
	}

	// iblk番目のブロックのidofblk番目の自由度を固定する
	inline bool SetBC(unsigned int iBlk, unsigned int iDofBlk){
		assert( iBlk < m_nBlk );
        if( this->m_lenBlk == -1 ){
            assert(0);
            return false;
        }
        assert( (int)iDofBlk < m_lenBlk );
		m_Flag[iBlk*m_lenBlk+iDofBlk] = 1;
		return true;
	}
private:
    const int m_lenBlk; // Flexサイズなら-1
    const unsigned int m_nBlk;
    unsigned int* m_DofPtr; // Flexサイズじゃないなら０
	int* m_Flag;
};

}	// end namespace 'Ls'

#endif
