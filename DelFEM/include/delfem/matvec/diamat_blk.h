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
@brief CMatDia_BlkCrs クラスのインターフェイス
@author Nobuyuki Umetani
*/


#if !defined(DIAMAT_BLK_H)
#define DIAMAT_BLK_H

#ifndef for 
#define for if(0); else for
#endif

#include <math.h>

#include "delfem/matvec/vector_blk.h"

namespace MatVec{


/*! 
@brief diagonal matrix class
@ingroup MatVec
*/
class CDiaMat_Blk
{
public:
	CDiaMat_Blk(const unsigned int nblk, const unsigned int nlen){
		m_nblk = nblk;
		m_nlen = nlen;
		this->m_valDia_Blk = new double [nblk*nlen*nlen];
	}
	CDiaMat_Blk();
	virtual ~CDiaMat_Blk(){
		delete[] m_valDia_Blk;
	}

	unsigned int NBlk(){
		return m_nblk;
	}
	unsigned int LenBlk(){
		return m_nlen;
	}

	//! 行列に０をセットするための関数
	virtual void SetZero(){
		for(unsigned int i=0;i<m_nblk*m_nlen*m_nlen;i++){ m_valDia_Blk[i] = 0.0; }
	}
	//! マージするための関数
	virtual bool Mearge(unsigned int iblk, unsigned int blksize, const double* emat ){
		if( iblk >= m_nblk ){ assert(0); return false; }
		if( blksize != m_nlen*m_nlen ){ assert(0); return false; }
		for(unsigned int i=0;i<m_nlen*m_nlen;i++){
			m_valDia_Blk[iblk*blksize+i] += emat[i];
		}
		return true;
	}

	bool CholeskyDecomp(){
		if( m_nlen == 1 ){
			for(unsigned int iblk=0;iblk<m_nblk;iblk++){
				if( m_valDia_Blk[iblk] >= 0.0 ){
//					std::cout << iblk << " " << m_valDia_Blk[iblk] << std::endl;
					m_valDia_Blk[iblk] = 1.0/sqrt( m_valDia_Blk[iblk] );
				}
				else{ assert(0); }
			}
		}
		else if( m_nlen == 2 ){
			for(unsigned int iblk=0;iblk<m_nblk;iblk++){
				assert( m_valDia_Blk[iblk*4+3] >= 0.0 );
				assert( m_valDia_Blk[iblk*4  ] >= 0.0 );
				m_valDia_Blk[iblk*4  ] = 1.0/sqrt( m_valDia_Blk[iblk*4  ] );
				m_valDia_Blk[iblk*4+3] = 1.0/sqrt( m_valDia_Blk[iblk*4+3] );
			}
		}
		else if( m_nlen == 3 ){
			for(unsigned int iblk=0;iblk<m_nblk;iblk++){
				assert( m_valDia_Blk[iblk*9  ] >= 0.0 );
				assert( m_valDia_Blk[iblk*9+4] >= 0.0 );
				assert( m_valDia_Blk[iblk*9+8] >= 0.0 );
				m_valDia_Blk[iblk*9  ] = 1.0/sqrt( m_valDia_Blk[iblk*9  ] );
				m_valDia_Blk[iblk*9+4] = 1.0/sqrt( m_valDia_Blk[iblk*9+4] );
				m_valDia_Blk[iblk*9+8] = 1.0/sqrt( m_valDia_Blk[iblk*9+8] );
			}
		}
		else{ assert(0); }
		return true;
	}

	/*!
	@brief 行列ベクトル積
	{lhs} = beta*{lhs} + alpha*[A]{rhs}
	*/
    virtual bool MatVec(double alpha, const MatVec::CVector_Blk& rhs, double beta, MatVec::CVector_Blk& lhs) const{
        if( rhs.Len() < 0 || lhs.Len() < 0 ){
			assert(0);
            return false;
        }
        const unsigned int nblk = m_nblk;
        const unsigned int nlen  = m_nlen;
        assert( rhs.Len() == (int)nlen );
        assert( rhs.NBlk() == nblk );
        assert( lhs.Len() == (int)nlen );
        assert( lhs.NBlk() == nblk );
		assert( &rhs != &lhs );
		if( beta != 1 ){ lhs *= beta; }

		if( nlen == 1 ){
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				double val0 = rhs.GetValue(iblk,0)*alpha;
				val0 *= m_valDia_Blk[iblk];
				lhs.AddValue(iblk,0,val0);
			}
		}
		else if( nlen == 2 ){
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				const double val0_a = rhs.GetValue(iblk,0)*alpha;
				const double val1_a = rhs.GetValue(iblk,1)*alpha;
				const double* pMat = &m_valDia_Blk[iblk*4];
				const double val0_b = pMat[0]*val0_a + pMat[1]*val1_a;
				const double val1_b = pMat[2]*val0_a + pMat[3]*val1_a;
				lhs.AddValue(iblk,0,val0_b);
				lhs.AddValue(iblk,1,val1_b);
			}
		}
		else{
			double* pval_vec_r = new double[nlen];
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				for(unsigned int ilen=0;ilen<nlen;ilen++){
					pval_vec_r[ilen] = rhs.GetValue(iblk,ilen)*alpha;
				}
				const double* pMat = &m_valDia_Blk[iblk*nlen*nlen];
				for(unsigned int ilen=0;ilen<nlen;ilen++){
					double dtmp1 = 0;
					for(unsigned int jlen=0;jlen<nlen;jlen++){
						dtmp1 += pMat[ilen*nlen+jlen]*pval_vec_r[jlen];
					}
					lhs.AddValue(iblk,ilen,dtmp1);
				}
			}
			delete[] pval_vec_r;
		}

		return true;
	}

	const double* GetPtrValDia(unsigned int iblk) const {
		assert( iblk < m_nblk );
		return &m_valDia_Blk[iblk*m_nlen*m_nlen];
	}

protected:
	unsigned int m_nblk;
	unsigned int m_nlen;
	double* m_valDia_Blk;
};

}	// end namespace 'Ls'

#endif // MATDIA_CRS_H
