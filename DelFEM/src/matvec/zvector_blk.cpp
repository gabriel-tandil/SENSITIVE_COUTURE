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

//////////////////////////////////////////////////////////////////////
// Vector_Blk.cpp: implementation of vector class CVector_Blk
//////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <cstdlib> //(abort)

#include "delfem/complex.h"
#include "delfem/matvec/zvector_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 非メンバのフレンドのオペレータ
//////////////////////////////////////////////////////////////////////

namespace MatVec{

Com::Complex operator*(const CZVector_Blk& lhs, const CZVector_Blk& rhs){	
	assert( lhs.m_BlkLen == rhs.m_BlkLen );
	assert( lhs.m_BlkVecLen == rhs.m_BlkVecLen );
  if( lhs.BlkLen() < 0 ){
    std::cout << "Error!-->Not Implemented" << std::endl;
    assert(0);
    abort();
  }
  const unsigned int ndof = (unsigned int)lhs.BlkLen() * lhs.BlkVecLen();
	double real0 = 0;
	double imag0 = 0;
	for(unsigned int idof=0;idof<ndof;idof++){ 
		double d0 = lhs.m_Value[idof].Real();
		double d1 = lhs.m_Value[idof].Imag();
		double d2 = rhs.m_Value[idof].Real();
		double d3 = rhs.m_Value[idof].Imag();
		real0 += d0*d2-d1*d3;
		imag0 += d0*d3+d1*d2;
	}
	return Com::Complex(real0,imag0);
}


Com::Complex InnerProduct(const CZVector_Blk& lhs, const CZVector_Blk& rhs){
	assert( lhs.m_BlkLen == rhs.m_BlkLen );
	assert( lhs.m_BlkVecLen == rhs.m_BlkVecLen );
    if( lhs.BlkLen() < 0 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        abort();
    }
    const unsigned int ndof = (unsigned int)lhs.BlkLen() * lhs.BlkVecLen();
	double real0 = 0;
	double imag0 = 0;
	for(unsigned int idof=0;idof<ndof;idof++){ 
		double d0 = lhs.m_Value[idof].Real();
		double d1 = lhs.m_Value[idof].Imag();
		double d2 = rhs.m_Value[idof].Real();
		double d3 = rhs.m_Value[idof].Imag();
		real0 += d0*d2+d1*d3;
		imag0 += d0*d3-d1*d2;
	}
	return Com::Complex(real0,imag0);
}

}

////////////////////////////////////////////////
// メンバーオペレータ
////////////////////////////////////////////////

CZVector_Blk& CZVector_Blk::operator*=(const double& d0){	// Scaler Product
	const int ndof = m_BlkLen * m_BlkVecLen;
	for(int idof=0;idof<ndof;idof++){ m_Value[idof] *= d0; }
	return *this; 
}

CZVector_Blk& CZVector_Blk::operator*=(const Com::Complex& c0){	// Scaler Product
	const int ndof = m_BlkLen * m_BlkVecLen;
	for(int idof=0;idof<ndof;idof++){ m_Value[idof] *= c0; }
	return *this; 
}

CZVector_Blk& CZVector_Blk::operator=(const CZVector_Blk& rhs){ // Substitue Vector
	assert( m_BlkLen == rhs.m_BlkLen ); 
	assert( m_BlkVecLen == rhs.m_BlkVecLen );
	const int ndof = m_BlkLen*m_BlkVecLen;
	for(int idof=0;idof<ndof;idof++){ m_Value[idof] = rhs.m_Value[idof]; }
	return *this; 
}

CZVector_Blk& CZVector_Blk::operator+=(const CZVector_Blk& rhs){ // Add 
	assert( m_BlkLen == rhs.m_BlkLen ); 
	assert( m_BlkVecLen == rhs.m_BlkVecLen );
	const int ndof = m_BlkLen*m_BlkVecLen;
	for(int idof=0;idof<ndof;idof++){ m_Value[idof] += rhs.m_Value[idof]; }
	return *this; 
}

////////////////////////////////////////////////
// メンバ関数
////////////////////////////////////////////////

CZVector_Blk& CZVector_Blk::AXPY(const double& alpha, const CZVector_Blk& rhs){ // Add scaler scaled Vector
	assert( m_BlkLen == rhs.m_BlkLen ); 
	assert( m_BlkVecLen == rhs.m_BlkVecLen );
	const int ndof = m_BlkLen * m_BlkVecLen;
	for(int idof=0;idof<ndof;idof++){ m_Value[idof] += alpha * rhs.m_Value[idof]; }
	return *this;
}


CZVector_Blk& CZVector_Blk::AXPY(const Com::Complex& alpha, const CZVector_Blk& rhs){ // Add scaler scaled Vector
	assert( m_BlkLen == rhs.m_BlkLen ); 
	assert( m_BlkVecLen == rhs.m_BlkVecLen );
    if( m_BlkLen < 0 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        abort();
    }
    const unsigned int ndof = (unsigned int)m_BlkLen * m_BlkVecLen;
	for(unsigned int idof=0;idof<ndof;idof++){ m_Value[idof] += alpha * rhs.m_Value[idof]; }
	return *this;
}

void CZVector_Blk::SetVectorZero(){	// Set 0 to Value
	const unsigned int ndof = m_BlkLen * m_BlkVecLen;
	for(unsigned int idof=0;idof<ndof;idof++){ m_Value[idof]=0.0; } 
}

void CZVector_Blk::SetVectorConjugate(){	// Set 0 to Value
	const unsigned int ndof = m_BlkLen * m_BlkVecLen;
	for(unsigned int idof=0;idof<ndof;idof++){ m_Value[idof] = Com::Conjugate(m_Value[idof]); } 
}

double CZVector_Blk::GetSquaredVectorNorm() const {
	const unsigned int ndof = m_BlkLen*m_BlkVecLen;
	double tmp_value = 0.0;
	for(unsigned int idof=0;idof<ndof;idof++){
		tmp_value += Com::SquaredNorm(m_Value[idof]);
	}
	return tmp_value;
}
