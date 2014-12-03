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
@brief 複素数ブロックベクトルクラス(MatVec::CZVector_Blk)のインターフェイス
@author Nobuyuki Umetani
*/

#if !defined(ZVECTOR_BLK_H)
#define ZVECTOR_BLK_H

#include "delfem/complex.h"

namespace MatVec{

class CZMat_BlkCrs;
class CZMatDia_BlkCrs;
class CZMatDiaInv_BlkDia;

/*!
@brief 複素数ベクトルクラス
@ingroup MatVec
*/
class CZVector_Blk  
{
	friend class CZMat_BlkCrs;
	friend class CZMatDia_BlkCrs;
	friend class CZMatDiaInv_BlkDia;
    friend Com::Complex operator*(const CZVector_Blk& lhs, const CZVector_Blk& rhs);
	friend Com::Complex InnerProduct(const CZVector_Blk& lhs, const CZVector_Blk& rhs);
public:
	CZVector_Blk(const unsigned int& iblkveclen, const unsigned int& iblklen)
		: m_BlkVecLen(iblkveclen), m_BlkLen(iblklen){
		m_Value = new Com::Complex [m_BlkLen*m_BlkVecLen];
	}
	virtual ~CZVector_Blk(){ if( m_Value!=0) delete[] m_Value; }

	CZVector_Blk& operator*=(const double& d0);	// Scaler Product
	CZVector_Blk& operator*=(const Com::Complex& c0);	// Scaler Product
	CZVector_Blk& operator=(const CZVector_Blk& rhs); // Substitue Vector
	CZVector_Blk& operator+=(const CZVector_Blk& rhs); // Add 


	////////////////////////////////
	// Menber function

	CZVector_Blk& AXPY(const double& alpha, const CZVector_Blk& rhs); // Add scaler scaled Vector
	CZVector_Blk& AXPY(const Com::Complex& alpha, const CZVector_Blk& rhs); // Add scaler scaled Vector

	void SetVectorZero();	// Set 0 to Value
	double GetSquaredVectorNorm() const;	// return norm of this vector
	void SetVectorConjugate();

    inline unsigned int BlkLen() const { return m_BlkLen; }	// Size of One Block
    inline unsigned int BlkVecLen() const { return m_BlkVecLen; }	// Size of Bocks Vector have

	const Com::Complex& GetValue(const unsigned int& iblk, const unsigned int& idofblk) const{
		assert( iblk<m_BlkVecLen ); assert( idofblk<m_BlkLen );
		return m_Value[iblk*m_BlkLen+idofblk]; 
	}
	void SetValue(const unsigned int& iblk, const unsigned int& idofblk, const Com::Complex& val){
		assert( iblk<m_BlkVecLen ); assert( idofblk<m_BlkLen );
		m_Value[iblk*m_BlkLen+idofblk] = val; 	
	}
	void AddValue(const unsigned int& iblk, const unsigned int& idofblk, const Com::Complex& val){
		assert( iblk<m_BlkVecLen ); assert( idofblk<m_BlkLen );
		m_Value[iblk*m_BlkLen+idofblk] += val; 	
	}
private:
    const unsigned int m_BlkVecLen;
    const unsigned int m_BlkLen;
	Com::Complex* m_Value;
};

}	// end namespace Ls

#endif // VEC_H
