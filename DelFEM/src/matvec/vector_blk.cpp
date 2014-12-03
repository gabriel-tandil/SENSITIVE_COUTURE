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

////////////////////////////////////////////////////////////////
// Vector_Blk.h : implentation of vector class (CVector_Blk)
////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <cstring> //(memcpy)

#include "delfem/matvec/vector_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 非メンバのフレンドのオペレータ
//////////////////////////////////////////////////////////////////////

namespace MatVec{

double operator*(const CVector_Blk& lhs, const CVector_Blk& rhs){	
	assert( lhs.Len() == rhs.Len() );
	assert( lhs.NBlk() == rhs.NBlk() );
    assert( lhs.GetTotalDofSize() == rhs.GetTotalDofSize() );
	double dot = 0.0;
	const double* prhs = rhs.m_Value;
	const double* plhs = lhs.m_Value;
	const unsigned int ndof = lhs.GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ dot += plhs[idof]*prhs[idof]; }
	return dot;
}

}

////////////////////////////////////////////////
// メンバーオペレータ
////////////////////////////////////////////////

CVector_Blk& CVector_Blk::operator*=(double d0){	// Scaler Product
	double* plhs = this->m_Value;
	const unsigned int ndof = this->GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ plhs[idof] *= d0; }
	return *this; 
}

CVector_Blk& CVector_Blk::operator=(const CVector_Blk& rhs){ // Substitue Vector
	assert( this->Len() == rhs.Len() ); 
	assert( this->NBlk() == rhs.NBlk() );
    assert( this->GetTotalDofSize() == rhs.GetTotalDofSize() );
	const void* prhs = (void*)rhs.m_Value;
	void* plhs = (void*)(this->m_Value);
	const unsigned int nbit = this->GetTotalDofSize()*8;
	memcpy(plhs,prhs,nbit);
	return *this; 
}

CVector_Blk& CVector_Blk::operator+=(const CVector_Blk& rhs){ // Add 
	assert( this->Len() == rhs.Len() ); 
	assert( this->NBlk() == rhs.NBlk() );
	assert( this->GetTotalDofSize() == rhs.GetTotalDofSize() );
	const double* prhs = rhs.m_Value;
	double* plhs = this->m_Value;
	const unsigned int ndof = this->GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ plhs[idof] += prhs[idof]; }
	return *this; 
}

////////////////////////////////////////////////
// メンバ関数
////////////////////////////////////////////////

CVector_Blk& CVector_Blk::AXPY(const double& alpha, const CVector_Blk& rhs){ // Add scaler scaled Vector
	assert( this->Len() == rhs.Len() ); 
	assert( this->NBlk() == rhs.NBlk() );
	assert( this->GetTotalDofSize() == rhs.GetTotalDofSize() );
	const double* prhs = rhs.m_Value;
	double* plhs = this->m_Value;
	const unsigned int ndof = this->GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ plhs[idof] += alpha*prhs[idof]; }
	return *this;
}


void CVector_Blk::SetVectorZero(){	// Set 0 to Value
	double* plhs = this->m_Value;
	const unsigned int ndof = this->GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ plhs[idof]=0.0; } 
}

double CVector_Blk::GetSquaredVectorNorm() const{
	double dtmp1 = 0.0;
	const double* plhs = this->m_Value;
	const unsigned int ndof = this->GetTotalDofSize();
	for(unsigned int idof=0;idof<ndof;idof++){ dtmp1 += plhs[idof]*plhs[idof]; }
	return dtmp1;
}
