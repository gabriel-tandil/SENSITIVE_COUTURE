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

////////////////////////////////////////////////////////////////
// MatDiaCrsFrac.cpp: implementation of ILU factorization class CMatDiaCrsFrac
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )   // C4786なんて表示すんな( ﾟДﾟ)ｺﾞﾙｧ
#endif
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <vector>
#include <set>
#include <math.h>

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"
#include "delfem/matvec/matfrac_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/ordering_blk.h"

//#include "ker_mat.h"

using namespace MatVec;

static void CalcMatPr(double* out, const double* d, double* tmp,
			          const unsigned int ni, const unsigned int nj )
{
	unsigned int i,j,k;
	for(i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			tmp[i*nj+j] = 0.0;
			for(k=0;k<ni;k++){
				tmp[i*nj+j] += d[i*ni+k]*out[k*nj+j];
			}
		}
	}
	for(i=0;i<ni*nj;i++){
		out[i] = tmp[i];
	}
}

static void CalcSubMatPr(double* out, const double* a, const double* b,
				  const unsigned int ni, const unsigned int nk, const unsigned int nj )
{
	unsigned int i,j,k;
	for(i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			for(k=0;k<nk;k++){
				out[i*nj+j] -= a[i*nk+k]*b[k*nj+j];
			}
		}
	}
}
/*
static void CalcInvMat(double* a, const unsigned int& n, int& info )
{
	double tmp1;

	info = 0;
	unsigned int i,j,k;
	for(i=0;i<n;i++){
		if( fabs(a[i+i*n]) < 1.0e-30 ){
			info = 1;
			return;
		}
		if( a[i+i*n] < 0.0 ){
			info--;
		}
		tmp1 = 1.0 / a[i+i*n];
		a[i+i*n] = 1.0;
		for(k=0;k<n;k++){
			a[i+k*n] *= tmp1;
		}
		for(j=0;j<n;j++){
			if( j!=i ){
				tmp1 = a[j+i*n];
				a[j+i*n] = 0.0;
				for(k=0;k<n;k++){
					a[j+k*n] -= tmp1*a[i+k*n];
				}
			}
		}
	}
}
*/

static void CalcInvMat(double* a, const unsigned int n, int& info )
{
	double tmp1;

	info = 0;
	unsigned int i,j,k;
	for(i=0;i<n;i++){
		if( fabs(a[i*n+i]) < 1.0e-30 ){
			info = 1;
			return;
		}
		if( a[i*n+i] < 0.0 ){
			info--;
		}
		tmp1 = 1.0 / a[i*n+i];
		a[i*n+i] = 1.0;
		for(k=0;k<n;k++){
			a[i*n+k] *= tmp1;
		}
		for(j=0;j<n;j++){
			if( j!=i ){
				tmp1 = a[j*n+i];
				a[j*n+i] = 0.0;
				for(k=0;k<n;k++){
					a[j*n+k] -= tmp1*a[i*n+k];
				}
			}
		}
	}
}

// tはtmp_bufferのこと
static inline void CalcInvMat3(double a[], double t[] )
{
	const double det = a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
		             - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
	const double inv_det = 1.0/det;

    for(int i=0;i<9;i++){ t[i] = a[i]; }

	a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
	a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
	a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);

	a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
	a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
	a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);

	a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
	a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
	a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
}
//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CMatDiaFrac_BlkCrs::CMatDiaFrac_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow)
	:CMatDia_BlkCrs(nblk_colrow,len_colrow)
{
	m_ConditionFlag = 0;
	m_DiaInd = 0;
	m_pRowLev = 0;
}

CMatDiaFrac_BlkCrs::CMatDiaFrac_BlkCrs(const CMatDia_BlkCrs& rhs)
{
	m_ConditionFlag = -1;
	m_DiaInd = 0;
	m_pRowLev = 0;

    // サイズを設定する
    if( rhs.LenBlkCol() >= 0 ){
        CMatDia_BlkCrs::Initialize(rhs.NBlkMatCol(),rhs.LenBlkCol());
    }
    else{
        std::vector<unsigned int> aLen;
        aLen.resize(rhs.NBlkMatCol());
        for(unsigned int iblk=0;iblk<rhs.NBlkMatCol();iblk++){
            aLen[iblk] = rhs.LenBlkCol(iblk);
        }
        CMatDia_BlkCrs::Initialize(rhs.NBlkMatCol(),aLen);
    }
    
    m_ConditionFlag = 0;

	CMatDia_BlkCrs::AddPattern(rhs,true);
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( LenBlkRow() == LenBlkCol() );
	{	// Assertion Routine
		const unsigned int nblk = NBlkMatCol();
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				if( icrs != m_colInd_Blk[iblk] ){	// check if row is in ascending order
					assert( jblk0 > m_rowPtr_Blk[icrs-1] );
				}
				assert( jblk0 < NBlkMatRow() );
			}
		}
	}
	{	// Make m_DiaInd
		const unsigned int nblk = NBlkMatCol();
		m_DiaInd = new unsigned int [nblk];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < NBlkMatRow() );
				if( jblk0 > iblk ){
					m_DiaInd[iblk] = icrs;
					break;
				}
			}
		}
	}
	m_ConditionFlag = 2;
}

CMatDiaFrac_BlkCrs::CMatDiaFrac_BlkCrs(const int lev_fill, const CMatDia_BlkCrs& rhs)
{	
	m_ConditionFlag = -1;
	m_DiaInd = 0;
	m_pRowLev = 0;

    // サイズを設定する
    if( rhs.LenBlkCol() >= 0 ){
        CMatDia_BlkCrs::Initialize(rhs.NBlkMatCol(),rhs.LenBlkCol());
    }
    else{
        std::vector<unsigned int> aLen;
        aLen.resize(rhs.NBlkMatCol());
        for(unsigned int iblk=0;iblk<rhs.NBlkMatCol();iblk++){
            aLen[iblk] = rhs.LenBlkCol(iblk);
        }
        CMatDia_BlkCrs::Initialize(rhs.NBlkMatCol(),aLen);
    }

    m_ConditionFlag = 0;

	if( lev_fill == 0 ){
		CMatDia_BlkCrs::AddPattern(rhs,true);
		{	// Assertion Routine
			const unsigned int nblk = NBlkMatCol();
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					if( icrs != m_colInd_Blk[iblk] ){	// check if row is in ascending order
						assert( jblk0 > m_rowPtr_Blk[icrs-1] );
					}
					assert( jblk0 < NBlkMatRow() );
				}
			}
		}
		{	// Make m_DiaInd
			const unsigned int nblk = NBlkMatCol();
			m_DiaInd = new unsigned int [nblk];
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < NBlkMatRow() );
					if( jblk0 > iblk ){
						m_DiaInd[iblk] = icrs;
						break;
					}
				}
			}
		}
		m_ConditionFlag = 2;
	}
	else{
		MakePattern_Initialize(rhs);
		AddFracPtn(lev_fill);
		MakePatternFinalize();
		assert( m_ConditionFlag == 2 );
	}
}

CMatDiaFrac_BlkCrs::CMatDiaFrac_BlkCrs(
        const int lev_fill, const CMatDia_BlkCrs& rhs, const COrdering_Blk& order)
:CMatDia_BlkCrs(rhs.NBlkMatCol(),rhs.LenBlkCol())
{	
	m_ConditionFlag = 0;
	m_DiaInd = 0;
	m_pRowLev = 0;
	if( lev_fill == 0 ){
		CMatDia_BlkCrs::AddPattern(rhs,order);
		{	// Make m_DiaInd
			const unsigned int nblk = NBlkMatCol();
			m_DiaInd = new unsigned int [nblk];
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < NBlkMatRow() );
					if( jblk0 > iblk ){
						m_DiaInd[iblk] = icrs;
						break;
					}
				}
			}
		}
		const unsigned int nlen = rhs.LenBlkCol();
		m_valCrs_Blk = new double [m_ncrs_Blk*nlen*nlen];
		m_ConditionFlag = 2;
	}
	else{
		MakePattern_Initialize(rhs,order);
		AddFracPtn(lev_fill);
		MakePatternFinalize();
		assert( m_ConditionFlag == 2 );
	}
}

CMatDiaFrac_BlkCrs::~CMatDiaFrac_BlkCrs()
{
	if( m_DiaInd != 0 ){ delete[] m_DiaInd; }
	if( m_pRowLev != 0 ){ delete m_pRowLev; }
}

bool CMatDiaFrac_BlkCrs::ForwardSubstitution( CVector_Blk& vec ) const
{

	assert( NBlkMatRow() == NBlkMatCol() );
	assert( vec.NBlk() == NBlkMatCol() );
	assert( vec.Len() == LenBlkCol() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        if( LenBlkCol() >= 0 || LenBlkRow() >= 0 ){
            std::cout << "Error!-->Not Implemented" << std::endl;
            assert(0);
            return false;
        }
        const unsigned int nblk = this->NBlkMatCol();
		double* pTmpVec;
        {
            unsigned int maxlen = 0;
            for(unsigned int iblk=0;iblk<nblk;iblk++){
                maxlen = ( maxlen > this->LenBlkCol(iblk) ) ? maxlen : this->LenBlkCol(iblk);
            }
            pTmpVec = new double [maxlen];
        }
		for(unsigned int iblk=0;iblk<nblk;iblk++){
            const unsigned int lencol = this->LenBlkCol(iblk); assert( vec.Len(iblk) == lencol );
			for(unsigned int idof=0;idof<lencol;idof++){
				pTmpVec[idof] = vec.GetValue(iblk,idof);
			}
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_DiaInd[iblk];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<iblk);
                const unsigned int lenrow = this->LenBlkRow(jblk0); assert( vec.Len(jblk0) == lenrow );
				const double* pVal_ij = &m_valCrs_Blk[ m_ValPtr[ijcrs] ];
				for(unsigned int idof=0;idof<lencol;idof++){
					for(unsigned int jdof=0;jdof<lenrow;jdof++){
						pTmpVec[idof] -= pVal_ij[idof*lenrow+jdof]*vec.GetValue(jblk0,jdof);
					}
				}
			}
			const double* pVal_ii = &m_valDia_Blk[ m_DiaValPtr[iblk] ];
			for(unsigned int idof=0;idof<lencol;idof++){
				double dtmp1 = 0.0;
				for(unsigned int jdof=0;jdof<lencol;jdof++){
					dtmp1 += pVal_ii[idof*lencol+jdof]*pTmpVec[jdof];
				}
				vec.SetValue(iblk,idof,dtmp1);
			}
		}
		delete[] pTmpVec;
        return true;
    }
    
    assert( LenBlkCol() >= 0 || LenBlkRow() >= 0 );

	if( LenBlkCol() == 1 ){	// ブロック行列のサイズが1×１な場合
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* matval_nd = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		for(unsigned int iblk=0;iblk<nblk;iblk++){	// 前進消去
			double lvec_i = vec.GetValue(iblk,0);
			for(unsigned int ijcrs=colind[iblk];ijcrs<diaind[iblk];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				lvec_i -= matval_nd[ijcrs]*vec.GetValue(jblk0,0);
			}
			vec.SetValue(iblk,0,matval_dia[iblk]*lvec_i);
		}
	}
	else if( LenBlkCol() == 2 ){
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* matval_nd = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		double pTmpVec[2];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			pTmpVec[0] = vec.GetValue(iblk,0);
			pTmpVec[1] = vec.GetValue(iblk,1);
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = diaind[iblk];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				const double* pVal_ij = &matval_nd[ijcrs*4];
                const double valj0 = vec.GetValue(jblk0,0);
				const double valj1 = vec.GetValue(jblk0,1);
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1;
				pTmpVec[1] -= pVal_ij[2]*valj0+pVal_ij[3]*valj1;
			}
			const double* pVal_ii = &matval_dia[iblk*4];
			const double dtmp1 = pVal_ii[0]*pTmpVec[0]+pVal_ii[1]*pTmpVec[1];
			vec.SetValue(iblk,0,dtmp1);
			const double dtmp2 = pVal_ii[2]*pTmpVec[0]+pVal_ii[3]*pTmpVec[1];
			vec.SetValue(iblk,1,dtmp2);
		}
	}
	else if( LenBlkCol() == 3 ){
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* matval_nd = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		double pTmpVec[3];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			pTmpVec[0] = vec.GetValue(iblk,0);
			pTmpVec[1] = vec.GetValue(iblk,1);
			pTmpVec[2] = vec.GetValue(iblk,2);
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = diaind[iblk];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0<iblk );
				const double* pVal_ij = &matval_nd[ijcrs*9];
				const double valj0 = vec.GetValue(jblk0,0);
				const double valj1 = vec.GetValue(jblk0,1);
				const double valj2 = vec.GetValue(jblk0,2);
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
				pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
				pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
			}
			const double* pVal_ii = &matval_dia[iblk*9];
			const double dtmp1 = pVal_ii[0]*pTmpVec[0]+pVal_ii[1]*pTmpVec[1]+pVal_ii[2]*pTmpVec[2];
			vec.SetValue(iblk,0,dtmp1);
			const double dtmp2 = pVal_ii[3]*pTmpVec[0]+pVal_ii[4]*pTmpVec[1]+pVal_ii[5]*pTmpVec[2];
			vec.SetValue(iblk,1,dtmp2);
			const double dtmp3 = pVal_ii[6]*pTmpVec[0]+pVal_ii[7]*pTmpVec[1]+pVal_ii[8]*pTmpVec[2];
			vec.SetValue(iblk,2,dtmp3);
		}
	}
	else{	// ブロック行列のサイズが3×3より大きい場合
		const unsigned int BlkLen = LenBlkCol();
		const unsigned int BlkSize = BlkLen*BlkLen;
		double* pTmpVec = new double [BlkLen];	// 作業用の小さな配列
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			for(unsigned int idof=0;idof<BlkLen;idof++){
				pTmpVec[idof] = vec.GetValue(iblk,idof);
			}
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_DiaInd[iblk];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0<iblk );
				const double* pVal_ij = &m_valCrs_Blk[ijcrs*BlkSize];
				for(unsigned int idof=0;idof<BlkLen;idof++){
					for(unsigned int jdof=0;jdof<BlkLen;jdof++){
						pTmpVec[idof] -= pVal_ij[idof*BlkLen+jdof]*vec.GetValue(jblk0,jdof);
					}
				}
			}
			const double* pVal_ii = &m_valDia_Blk[iblk*BlkSize];
			for(unsigned int idof=0;idof<BlkLen;idof++){
				double dtmp1 = 0.0;
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){
					dtmp1 += pVal_ii[idof*BlkLen+jdof]*pTmpVec[jdof];
				}
				vec.SetValue(iblk,idof,dtmp1);
			}
		}
		delete[] pTmpVec;
	}
	return true;
}

bool CMatDiaFrac_BlkCrs::BackwardSubstitution( CVector_Blk& vec ) const
{
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( vec.NBlk() == NBlkMatCol() );
	assert( vec.Len() == LenBlkCol() );
	assert( vec.Len() == LenBlkCol() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        if( LenBlkCol() >= 0 || LenBlkRow() >= 0 ){
            std::cout << "Error!-->Not Implemented" << std::endl;
            assert(0);
            return false;
        }
        const unsigned int nblk = this->NBlkMatCol();
		double* pTmpVec;
        {
            unsigned int maxlen = 0;
            for(unsigned int iblk=0;iblk<nblk;iblk++){
                maxlen = ( maxlen > this->LenBlkCol(iblk) ) ? maxlen : this->LenBlkCol(iblk);
            }
            pTmpVec = new double [maxlen];
        }
        assert( nblk > 0 );
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (unsigned int)iblk < NBlkMatCol() );
            const unsigned int lencol = this->LenBlkCol(iblk);
            assert( lencol == vec.Len(iblk) );
			for(unsigned int idof=0;idof<lencol;idof++){
				pTmpVec[idof] = vec.GetValue(iblk,idof);
			}
            for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; 
				assert( jblk0>(unsigned int)iblk && jblk0<NBlkMatRow() );
                const unsigned int lenrow = this->LenBlkRow(jblk0); assert( lenrow == vec.Len(jblk0) );
				const double* pVal_ij = &m_valCrs_Blk[ this->m_ValPtr[ijcrs] ];
				for(unsigned int idof=0;idof<lencol;idof++){
					for(unsigned int jdof=0;jdof<lenrow;jdof++){
						pTmpVec[idof] -= pVal_ij[idof*lenrow+jdof]*vec.GetValue(jblk0,jdof);
					}
				}
			}
			for(unsigned int idof=0;idof<lencol;idof++){
				vec.SetValue(iblk,idof,pTmpVec[idof]);
			}
		}
		delete[] pTmpVec;
        return true;
    }

	if( LenBlkCol() == 1 ){	// ブロック行列のサイズが1×１な場合
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* matval_nd = m_valCrs_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(int iblk=nblk-1;iblk>=0;iblk--){	// 後退代入
			assert( (unsigned int)iblk < nblk );
			double lvec_i = vec.GetValue(iblk,0);
			for(unsigned int ijcrs=diaind[iblk];ijcrs<colind[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0>(unsigned int)iblk && jblk0<nblk );
				lvec_i -=  matval_nd[ijcrs]*vec.GetValue(jblk0,0);
			}
			vec.SetValue(iblk,0,lvec_i);
		}
	}
	else if( LenBlkCol() == 2 ){
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* valmat_nd = m_valCrs_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		double pTmpVec[2];
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (unsigned int)iblk < nblk );
			pTmpVec[0] = vec.GetValue(iblk,0);
			pTmpVec[1] = vec.GetValue(iblk,1);
			const unsigned int icrs0 = diaind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0>(unsigned int)iblk && jblk0<nblk );
				const double* pVal_ij = &valmat_nd[ijcrs*4];
				const double valj0 = vec.GetValue(jblk0,0);
				const double valj1 = vec.GetValue(jblk0,1);
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1;
				pTmpVec[1] -= pVal_ij[2]*valj0+pVal_ij[3]*valj1;
			}
			vec.SetValue(iblk,0,pTmpVec[0]);
			vec.SetValue(iblk,1,pTmpVec[1]);
		}
	}
	else if( LenBlkCol() == 3 ){
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		const double* valmat_nd = m_valCrs_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		double pTmpVec[3];
		for(int iblk=nblk-1;iblk>=0;iblk--){
			assert( (unsigned int)iblk < nblk );
			pTmpVec[0] = vec.GetValue(iblk,0);
			pTmpVec[1] = vec.GetValue(iblk,1);
			pTmpVec[2] = vec.GetValue(iblk,2);
			const unsigned int icrs0 = diaind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int ijcrs=icrs0;ijcrs<icrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs];
				assert( jblk0>(unsigned int)iblk && jblk0<nblk );
				const double* pVal_ij = &valmat_nd[ijcrs*9];
				const double valj0 = vec.GetValue(jblk0,0);
				const double valj1 = vec.GetValue(jblk0,1);
				const double valj2 = vec.GetValue(jblk0,2);
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
				pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
				pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
			}
			vec.SetValue(iblk,0,pTmpVec[0]);
			vec.SetValue(iblk,1,pTmpVec[1]);
			vec.SetValue(iblk,2,pTmpVec[2]);
		}
	}
	else{
		const unsigned int BlkLen = LenBlkCol();
		const unsigned int BlkSize = BlkLen*BlkLen;
		double* pTmpVec = new double [BlkLen];	// 作業用の小さな配列
		for(int iblk=NBlkMatCol()-1;iblk>=0;iblk--){
			assert( (unsigned int)iblk < NBlkMatCol() );
			for(unsigned int idof=0;idof<BlkLen;idof++){
				pTmpVec[idof] = vec.GetValue(iblk,idof);
			}
			for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0>(unsigned int)iblk && jblk0<NBlkMatRow() );
				const double* pVal_ij = &m_valCrs_Blk[ijcrs*BlkSize];
				for(unsigned int idof=0;idof<BlkLen;idof++){
					for(unsigned int jdof=0;jdof<BlkLen;jdof++){
						pTmpVec[idof] -= pVal_ij[idof*BlkLen+jdof]*vec.GetValue(jblk0,jdof);
					}
				}
			}
			for(unsigned int idof=0;idof<BlkLen;idof++){
				vec.SetValue(iblk,idof,pTmpVec[idof]);
			}
		}
		delete[] pTmpVec;
	}
	return true;
}

bool CMatDiaFrac_BlkCrs::DoILUDecompLowUp( 
	const CMatFrac_BlkCrs& mat_low, const CMatFrac_BlkCrs& mat_up )
{
	assert( NBlkMatRow() == NBlkMatCol() );
    assert( mat_up.NBlkMatCol() == mat_low.NBlkMatRow() );
	assert( mat_low.NBlkMatCol() == NBlkMatCol() );
	assert( mat_up.NBlkMatRow() == NBlkMatRow() );
	assert( LenBlkRow() == LenBlkCol() );

	const unsigned int nblk = NBlkMatCol();
	const unsigned int nblk_dia = mat_up.NBlkMatCol();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
	    int* row2crs_f = new int [nblk];
	    for(unsigned  jblk=0;jblk<nblk;jblk++){ row2crs_f[jblk] = -1; }
	    for(unsigned int iblk=0;iblk<nblk;iblk++){
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<this->m_ncrs_Blk );
			    const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
			    row2crs_f[jblk0] = ijcrs;
		    }
            const unsigned int len_i = this->LenBlkCol(iblk);
            assert( mat_low.LenBlkCol(iblk) == len_i );
            for(unsigned int ikcrs=mat_low.m_colInd_Blk[iblk];ikcrs<mat_low.m_colInd_Blk[iblk+1];ikcrs++){ assert( ikcrs<mat_low.m_ncrs_Blk );
			    const double* ikvalue = mat_low.GetValCrsPtr(ikcrs);
			    const unsigned int kblk0 = mat_low.m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk_dia );
                const unsigned int len_k = mat_up.LenBlkCol(kblk0);
                assert( mat_low.LenBlkRow(kblk0) == len_k );
                for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<mat_up.m_ncrs_Blk ); 
				    const double* kjvalue = mat_up.GetValCrsPtr(kjcrs);
				    const unsigned int jblk0 = mat_up.m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk );
                    const unsigned int len_j = this->LenBlkRow(jblk0);
                    assert( mat_up.LenBlkRow(jblk0) == len_j );
				    if( iblk != jblk0 ){
					    const int ijcrs0 = row2crs_f[jblk0];
                        if( ijcrs0 == -1 ){ continue; }
                        assert( ijcrs0>=0 && ijcrs0<(int)m_ncrs_Blk );
					    double* ijvalue = &m_valCrs_Blk[ m_ValPtr[ijcrs0] ];
                        CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_i,len_k,len_j);	// Aij = Aij-Aik*Akk*Akj
				    }
				    else{
                        assert( len_i == len_j );
					    double* iivalue = &this->m_valDia_Blk[ m_DiaValPtr[iblk] ];
					    CalcSubMatPr(iivalue,ikvalue,kjvalue,len_i,len_k,len_j);	// Aii = Aii-Aik*Akk*Aki
				    }
			    }		
		    }
            for(unsigned ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
			    const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
			    row2crs_f[jblk0] = -1;
		    }
		    ////////////////
	    }	// inode
	    delete[] row2crs_f;
        return true;
    }

//	assert( mat_low.LenBlkCol() == LenBlkCol() );
//  assert( mat_up.LenBlkRow() == LenBlkRow() );

	const unsigned int len = LenBlkCol();

	int* row2crs_f = new int [nblk];
	for(unsigned  jblk=0;jblk<nblk;jblk++){ row2crs_f[jblk] = -1; }
	for(unsigned int iblk=0;iblk<nblk;iblk++){
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<this->m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
			row2crs_f[jblk0] = ijcrs;
		}
        assert( mat_low.LenBlkCol(iblk) == len );
        for(unsigned int ikcrs=mat_low.m_colInd_Blk[iblk];ikcrs<mat_low.m_colInd_Blk[iblk+1];ikcrs++){ assert( ikcrs<mat_low.m_ncrs_Blk );
			const double* ikvalue = mat_low.GetValCrsPtr(ikcrs);
			const unsigned int kblk0 = mat_low.m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk_dia );
            const unsigned int len_dia = mat_low.LenBlkRow(kblk0);
            assert( mat_up.LenBlkCol( kblk0) == len_dia );
            for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<mat_up.m_ncrs_Blk ); 
				const double* kjvalue = mat_up.GetValCrsPtr(kjcrs);
				const unsigned int jblk0 = mat_up.m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk );
                assert( mat_up.LenBlkRow(jblk0) == len );
				if( iblk != jblk0 ){
					const int ijcrs0 = row2crs_f[jblk0];
                    if( ijcrs0 == -1 ){ continue; }
                    assert( ijcrs0>=0 && ijcrs0<(int)m_ncrs_Blk );
					double* ijvalue = &m_valCrs_Blk[ijcrs0*len*len];
					CalcSubMatPr(ijvalue,ikvalue,kjvalue,len,len_dia,len);	// Aij = Aij-Aik*Akk*Akj
				}
				else{
					double* iivalue = &this->m_valDia_Blk[iblk*len*len];
					CalcSubMatPr(iivalue,ikvalue,kjvalue,len,len_dia,len);	// Aii = Aii-Aik*Akk*Aki
				}
			}		
		}
        for(unsigned ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
			row2crs_f[jblk0] = -1;
		}
		////////////////
	}	// inode
	delete[] row2crs_f;

	return true;
}

// numerical factorization
bool CMatDiaFrac_BlkCrs::DoILUDecomp()
{
  
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( this->LenBlkRow() == LenBlkCol() );
	const unsigned int nmax_sing = 10;	// ブロックILU分解に失敗していい数の上限
	unsigned int icnt_sing = 0;
  
  if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
    if( LenBlkCol() >= 0 || LenBlkRow() >= 0 ){
      std::cout << "Error!-->Not Implemented" << std::endl;
      assert(0);
      return false;
    }
    assert( this->m_ValPtr != 0 );
    const unsigned int nblk = this->NBlkMatCol();
    
    double* pTmpBlk;
    {
      unsigned int maxlen = 0;
      for(unsigned int iblk=0;iblk<nblk;iblk++){
        maxlen = ( maxlen > this->LenBlkCol(iblk) ) ? maxlen : this->LenBlkCol(iblk);
      }
      pTmpBlk = new double [maxlen*maxlen];
    }
    int* row2crs_f = new int [nblk];
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      const unsigned int lenblk_i = this->LenBlkCol(iblk);
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
        const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
        row2crs_f[jblk0] = ijcrs;
      }
      for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_DiaInd[iblk];ikcrs++){
        const unsigned int kblk0 = m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk );
        const unsigned int lenblk_k = this->LenBlkCol(kblk0); 
        const double* pVal_ik = &this->m_valCrs_Blk[ m_ValPtr[ikcrs] ];
        for(unsigned int kjcrs=m_DiaInd[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){
          const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk && jblk0>kblk0 );
          const unsigned int lenblk_j = this->LenBlkRow(jblk0);
          const double* pVal_kj = &m_valCrs_Blk[ m_ValPtr[kjcrs] ]; assert( pVal_kj != 0 );
          double* pVal_ij = 0;
          if( jblk0 != iblk ){
            const int ijcrs0 = row2crs_f[jblk0];
            if( ijcrs0 == -1 ){ continue; }
            pVal_ij = &m_valCrs_Blk[ m_ValPtr[ijcrs0] ];
          }
          else{ pVal_ij = &m_valDia_Blk[ m_DiaValPtr[iblk] ]; }
          assert( pVal_ij != 0 );
          CalcSubMatPr(pVal_ij,pVal_ik,pVal_kj, lenblk_i,lenblk_k,lenblk_j);
        }
      }
      {
        double* pVal_ii = &m_valDia_Blk[ m_DiaValPtr[iblk] ];
        int info = 0;
        CalcInvMat(pVal_ii,lenblk_i,info);
        if( info==1 ){
          //				    std::cout << "frac false" << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){ 
						delete[] row2crs_f;
						delete[] pTmpBlk;
						return false; 
					}
        }
      }
      for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
        double* pVal_ij = &m_valCrs_Blk[ m_ValPtr[ijcrs] ];
        const double* pVal_ii = &m_valDia_Blk[ m_DiaValPtr[iblk] ];
        const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
        const unsigned int lenblk_j = this->LenBlkRow(jblk0);
        CalcMatPr(pVal_ij,pVal_ii,pTmpBlk,  lenblk_i,lenblk_j);
      }
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
        const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk );
        row2crs_f[jblk0] = -1;
      }
    }	// end iblk
    delete[] row2crs_f;
    delete[] pTmpBlk;
    return true;
  }
  
	const unsigned int BlkLen = LenBlkCol();
	const unsigned int BlkSize = BlkLen*BlkLen;
  
  //	int info= 0;
	int* row2crs_f = new int [NBlkMatRow()];
	for(unsigned int jblk=0;jblk<NBlkMatRow();jblk++){ row2crs_f[jblk] = -1; }
	if( BlkLen == 1 ){
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			// 非ゼロフラグの初期化
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = ijcrs;
			}
			// [L] * [D^-1*U] の計算
			for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_DiaInd[iblk];ikcrs++){
				const unsigned int kblk = m_rowPtr_Blk[ikcrs]; assert( kblk<NBlkMatCol() );
				const double ikvalue = m_valCrs_Blk[ikcrs];
				for(unsigned int kjcrs=m_DiaInd[kblk];kjcrs<m_colInd_Blk[kblk+1];kjcrs++){
					const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<NBlkMatRow() );
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs_f[jblk0];
						if( ijcrs0 == -1 ) continue;	// もしパターンがなければ戻る
						m_valCrs_Blk[ijcrs0] -= ikvalue*m_valCrs_Blk[kjcrs];
					}
					else{ m_valDia_Blk[iblk] -= ikvalue*m_valCrs_Blk[kjcrs]; }
				}
			}
			// 対角の逆数を計算
			double iivalue = m_valDia_Blk[iblk];
			if( fabs(iivalue) > 1.0e-30 ){
				m_valDia_Blk[iblk] = 1.0 / iivalue; 
			}
			else{ 
				std::cout << "frac false" << iblk << std::endl;
				icnt_sing++;
				if( icnt_sing > nmax_sing ){
					delete[] row2crs_f;
					return false;
				}
			}
			// 対角の逆数を計算して上三角行列に掛ける
      for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				m_valCrs_Blk[ijcrs] = m_valCrs_Blk[ijcrs] * m_valDia_Blk[iblk];
			}
			// 非ゼロフラグをもとに戻す
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = -1;
			}
		}	// end iblk
	}
	////////////////////////////////////////////////////////////////
	else if( BlkLen == 2 ){
		// 高速化のためにメンバ変数をローカルにコピー
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const unsigned int* diaind = m_DiaInd;
		double* matval_nd = m_valCrs_Blk;
		double* matval_dia = m_valDia_Blk;
		const unsigned int nblk = this->NBlkMatCol();
		double TmpBlk[4];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			// 非ゼロフラグの初期化
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = ijcrs;
			}
			// i行目について分解を行う。[L] * [D^-1*U] を計算する。
			for(unsigned int ikcrs=colind[iblk];ikcrs<diaind[iblk];ikcrs++){
				const unsigned int kblk = rowptr[ikcrs]; assert( kblk<nblk );
				const double* pVal_ik = &matval_nd[ikcrs*BlkSize];
				for(unsigned int kjcrs=diaind[kblk];kjcrs<colind[kblk+1];kjcrs++){
					const unsigned int jblk0 = rowptr[kjcrs]; assert( jblk0<nblk );
					double* pVal_kj = &matval_nd[kjcrs*4]; assert( pVal_kj != 0 );
					double* pVal_ij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs_f[jblk0];
						if( ijcrs0 == -1 ) continue;	// もしパターンがなければ戻る
						pVal_ij = &matval_nd[ijcrs0*4];
					}
          else{ pVal_ij = &matval_dia[iblk *4]; }
          assert( pVal_ij != 0 );
					pVal_ij[0] -= pVal_ik[0]*pVal_kj[0]+pVal_ik[1]*pVal_kj[2];
					pVal_ij[1] -= pVal_ik[0]*pVal_kj[1]+pVal_ik[1]*pVal_kj[3];
					pVal_ij[2] -= pVal_ik[2]*pVal_kj[0]+pVal_ik[3]*pVal_kj[2];
					pVal_ij[3] -= pVal_ik[2]*pVal_kj[1]+pVal_ik[3]*pVal_kj[3];
				}
			}
			{	// 対角の逆数を計算
				double* pVal_ii = &matval_dia[iblk*4];
				const double det = pVal_ii[0]*pVal_ii[3]-pVal_ii[1]*pVal_ii[2];
				if( fabs(det) > 1.0e-30 ){
					const double inv_det = 1.0/det;
					double dtmp1 = pVal_ii[0];
					pVal_ii[0] =  inv_det*pVal_ii[3];
					pVal_ii[1] = -inv_det*pVal_ii[1];
					pVal_ii[2] = -inv_det*pVal_ii[2];
					pVal_ii[3] =  inv_det*dtmp1;
				}
				else{
					std::cout << "frac false" << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
						delete[] row2crs_f;
						return false;
					}
				}
			}
			// 対角の逆数を計算して上三角行列に掛ける。[U] = [1/D][U]
      for(unsigned int ijcrs=diaind[iblk];ijcrs<colind[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				double* pVal_ij = &matval_nd[ijcrs*4];
				const double* pVal_ii = &matval_dia[iblk*4];
				for(unsigned int i=0;i<4;i++){ TmpBlk[i] = pVal_ij[i]; }
				pVal_ij[0] = pVal_ii[0]*TmpBlk[0] + pVal_ii[1]*TmpBlk[2];
				pVal_ij[1] = pVal_ii[0]*TmpBlk[1] + pVal_ii[1]*TmpBlk[3];
				pVal_ij[2] = pVal_ii[2]*TmpBlk[0] + pVal_ii[3]*TmpBlk[2];
				pVal_ij[3] = pVal_ii[2]*TmpBlk[1] + pVal_ii[3]*TmpBlk[3];
			}
			// 非ゼロフラグをもとに戻す
      for(unsigned int ijcrs=colind[iblk];ijcrs<colind[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[ijcrs]; assert( jblk0<nblk );
				row2crs_f[jblk0] = -1;
			}
		}	// end iblk
	}
	////////////////////////////////////////////////////////////////
	else if( BlkLen == 3 ){	// lenBlk >= 3
		double tmpBlk[9];
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			// 非ゼロフラグの初期化
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = ijcrs;
			}
			// i行目について分解を行う。[L] * [D^-1*U] を計算する。
			for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_DiaInd[iblk];ikcrs++){
				const unsigned int kblk = m_rowPtr_Blk[ikcrs]; assert( kblk<NBlkMatCol() );
				const double* pVal_ik = &m_valCrs_Blk[ikcrs*BlkSize];
				for(unsigned int kjcrs=m_DiaInd[kblk];kjcrs<m_colInd_Blk[kblk+1];kjcrs++){
					const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<NBlkMatRow() );
					double* pVal_kj = &m_valCrs_Blk[kjcrs *BlkSize]; assert( pVal_kj != 0 );
					double* pVal_ij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs_f[jblk0];
            if( ijcrs0 == -1 ){ continue; }// もしパターンがなければ戻る
						pVal_ij = &m_valCrs_Blk[ijcrs0*BlkSize];
					}
          else{ pVal_ij = &m_valDia_Blk[iblk *BlkSize]; }
					assert( pVal_ij != 0 );
          for(unsigned int i=0;i<3;i++){
            pVal_ij[i*3+0] -= pVal_ik[i*3+0]*pVal_kj[0] + pVal_ik[i*3+1]*pVal_kj[3] + pVal_ik[i*3+2]*pVal_kj[6];
            pVal_ij[i*3+1] -= pVal_ik[i*3+0]*pVal_kj[1] + pVal_ik[i*3+1]*pVal_kj[4] + pVal_ik[i*3+2]*pVal_kj[7];
            pVal_ij[i*3+2] -= pVal_ik[i*3+0]*pVal_kj[2] + pVal_ik[i*3+1]*pVal_kj[5] + pVal_ik[i*3+2]*pVal_kj[8];
          }
				}
			}
			{   // 対角の逆数を計算
				double* a_ii = &m_valDia_Blk[iblk*BlkSize];
				const double det = a_ii[0]*a_ii[4]*a_ii[8] + a_ii[3]*a_ii[7]*a_ii[2] + a_ii[6]*a_ii[1]*a_ii[5]
        - a_ii[0]*a_ii[7]*a_ii[5] - a_ii[6]*a_ii[4]*a_ii[2] - a_ii[3]*a_ii[1]*a_ii[8];				
				if( fabs(det) > 1.0e-30 ){
					CalcInvMat3(a_ii,tmpBlk);
				}
				else{
					std::cout << "frac false 3 " << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
						delete[] row2crs_f;
            std::cout << "ilu frac false exceeds tolerance" << std::endl;
						return false;
					}
				}
			}
			// 対角の逆数を計算して上三角行列に掛ける。[U] = [1/D][U]
      for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				double* pVal_ij = &m_valCrs_Blk[ijcrs*BlkSize];
				const double* pVal_ii = &m_valDia_Blk[iblk *BlkSize];
        for(unsigned int i=0;i<9;i++){ tmpBlk[i] = pVal_ij[i]; }
        for(unsigned int i=0;i<3;i++){
          pVal_ij[i*3+0] = pVal_ii[i*3+0]*tmpBlk[0] + pVal_ii[i*3+1]*tmpBlk[3] + pVal_ii[i*3+2]*tmpBlk[6];
          pVal_ij[i*3+1] = pVal_ii[i*3+0]*tmpBlk[1] + pVal_ii[i*3+1]*tmpBlk[4] + pVal_ii[i*3+2]*tmpBlk[7];
          pVal_ij[i*3+2] = pVal_ii[i*3+0]*tmpBlk[2] + pVal_ii[i*3+1]*tmpBlk[5] + pVal_ii[i*3+2]*tmpBlk[8];
        }
			}
			// 非ゼロフラグをもとに戻す
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = -1;
			}
		}	// end iblk
	}
  ////////////////////////////////////////////////////////////////
	else{	// lenBlk >= 4
		double* pTmpBlk = new double [BlkSize];
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			// 非ゼロフラグの初期化
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = ijcrs;
			}
			// i行目について分解を行う。[L] * [D^-1*U] を計算する。
			for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_DiaInd[iblk];ikcrs++){
				const unsigned int kblk = m_rowPtr_Blk[ikcrs]; assert( kblk<NBlkMatCol() );
				const double* pVal_ik = &m_valCrs_Blk[ikcrs*BlkSize];
				for(unsigned int kjcrs=m_DiaInd[kblk];kjcrs<m_colInd_Blk[kblk+1];kjcrs++){
					const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<NBlkMatRow() );
					double* pVal_kj = &m_valCrs_Blk[kjcrs *BlkSize]; assert( pVal_kj != 0 );
					double* pVal_ij = 0;
					if( jblk0 != iblk ){
						const int ijcrs0 = row2crs_f[jblk0];
            if( ijcrs0 == -1 ){ continue; }// もしパターンがなければ戻る
						pVal_ij = &m_valCrs_Blk[ijcrs0*BlkSize];
					}
          else{ pVal_ij = &m_valDia_Blk[iblk *BlkSize]; }
					assert( pVal_ij != 0 );
          CalcSubMatPr(pVal_ij,pVal_ik,pVal_kj, BlkLen,BlkLen,BlkLen);
				}
			}
			// 対角の逆数を計算
			{
				double* pVal_ii = &m_valDia_Blk[iblk*BlkSize];
				int info = 0;
				CalcInvMat(pVal_ii,BlkLen,info);
				if( info==1 ){
					std::cout << "frac false" << iblk << std::endl;
					icnt_sing++;
					if( icnt_sing > nmax_sing ){
						delete[] pTmpBlk;
						delete[] row2crs_f;
						return false;
					}
				}
			}
			// 対角の逆数を計算して上三角行列に掛ける。[U] = [1/D][U]
      for(unsigned int ijcrs=m_DiaInd[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				double* pVal_ij = &m_valCrs_Blk[ijcrs*BlkSize];
				const double* pVal_ii = &m_valDia_Blk[iblk *BlkSize];
        CalcMatPr(pVal_ij,pVal_ii,pTmpBlk,  BlkLen,BlkLen);
			}
			// 非ゼロフラグをもとに戻す
      for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<NBlkMatRow() );
				row2crs_f[jblk0] = -1;
			}
		}	// end iblk
		delete[] pTmpBlk;
	}
  
	delete[] row2crs_f;
  
	return true;  
}


bool CMatDiaFrac_BlkCrs::AddFracPtnLowUp( const int lev_fill, const CMatFrac_BlkCrs& mat_low, const CMatFrac_BlkCrs& mat_up )
{
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( mat_low.NBlkMatCol() == NBlkMatCol() );
	assert( mat_up.NBlkMatRow() == NBlkMatRow() );
	assert( mat_up.NBlkMatCol() == mat_low.NBlkMatRow() );

	assert( m_ConditionFlag != -1 );
	if(  m_ConditionFlag == -1 ) return false;
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

	if( m_ncrs_Blk == NBlkMatCol() * ( NBlkMatRow() - 1 ) ){
		std::cout << "   Full Matrix " << std::endl;
		return true;
	}

	if( lev_fill == 0 ){
		std::cout << "    Fill Lev 0 " << std::endl;
		return true;
	}

    ////////////////

	if( m_ConditionFlag == 2 ){
		if( m_valCrs_Blk    != 0 ){ 	// 値配列が確保されている場合は確保し直す
            std::cout << "deallocate value list " << std::endl;
            delete[] m_valCrs_Blk; 
            m_valCrs_Blk = 0; 
        }
        assert( m_DiaInd != 0 );
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < NBlkMatRow() );
				if( jblk0 > iblk ){
					m_DiaInd[iblk] = icrs;
					break;
				}
			}
		}

		m_pRowLev = new std::vector<CRowLev>;
		m_pRowLev->resize( m_ncrs_Blk );
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
			(*m_pRowLev)[icrs].lev = 0;
		}
		delete[] m_rowPtr_Blk;
		m_rowPtr_Blk = 0;

		m_ConditionFlag = 1;
	}

    ////////////////////////////////////////////////////////////////
	const unsigned int* ColInd_pre = this->m_colInd_Blk;
	const unsigned int* DiaInd_pre = this->m_DiaInd;
	const std::vector<CRowLev>* pRowLev_pre = this->m_pRowLev;

	const unsigned int nnode = NBlkMatCol();
	const unsigned int ncrs_pre = NCrs();

	m_colInd_Blk = new unsigned int [nnode+1];
	m_colInd_Blk[0] = 0;
	m_DiaInd = new unsigned int [nnode];
	m_DiaInd[0] = 0;
	assert( m_rowPtr_Blk == 0 );
	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*2);

    ////////////////////////////////
    CListRLN nonzero(NBlkMatRow());
	for(unsigned int iblk=0;iblk<nnode;iblk++){
        nonzero.Initialize((*pRowLev_pre), ColInd_pre[iblk], ColInd_pre[iblk+1] );
		if( mat_low.m_ConditionFlag == 1 && mat_up.m_ConditionFlag == 1 ){
            for(unsigned int ikcrs_low=mat_low.m_colInd_Blk[iblk];ikcrs_low<mat_low.m_colInd_Blk[iblk+1];ikcrs_low++){
				assert( ikcrs_low>=0 && ikcrs_low<mat_low.NCrs() );
				const int kblk0 = (*mat_low.m_pRowLev)[ikcrs_low].row;
                assert( kblk0>=0 && kblk0<(int)mat_low.NBlkMatRow() );
				const int ik_lev0 = (*mat_low.m_pRowLev)[ikcrs_low].lev;
				if( ik_lev0+1>lev_fill && lev_fill!=-1 ) continue;
                int jnz_cur = nonzero.ipos_entry;
			    for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){
				    const int kj_lev0 = (*mat_up.m_pRowLev)[kjcrs].lev;
				    if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				    const unsigned int jblk0 = (*mat_up.m_pRowLev)[kjcrs].row;
				    assert( jblk0<NBlkMatRow() );
				    if( jblk0 == iblk ) continue;   // 対角にＣＲＳ成分は無いとは思うけど一応
				    const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
			    }
		    }
		}
		else if( mat_low.m_ConditionFlag == 2 && mat_up.m_ConditionFlag == 2 ){
            for(unsigned int ikcrs_low=mat_low.m_colInd_Blk[iblk];ikcrs_low<mat_low.m_colInd_Blk[iblk+1];ikcrs_low++){
				assert( ikcrs_low>=0 && ikcrs_low<mat_low.NCrs() );
                const unsigned int kblk0 = mat_low.m_rowPtr_Blk[ikcrs_low];
				assert( kblk0>=0 && kblk0<mat_low.NBlkMatRow() );
				const int ik_lev0 = 0;
				if( ik_lev0+1>lev_fill && lev_fill!=-1 ) continue;
                int jnz_cur = nonzero.ipos_entry;
			    for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){
				    const int kj_lev0 = 0;
				    if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				    const unsigned int jblk0 = mat_up.m_rowPtr_Blk[kjcrs];
				    assert( jblk0<NBlkMatRow() );
				    if( jblk0 == iblk ) continue;   // 対角にＣＲＳ成分は無いとは思うけど一応
				    const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
			    }
		    }
		}
		else{
            std::cout << "Error!-->Not Implemente" << std::endl;
			assert(0);
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk,m_DiaInd);
	}
/*
	for(int inode=0;inode<nnode;inode++){
		std::cout << inode << "-->";
		for(int ijcrs=m_BlkColInd[inode];ijcrs<m_BlkDiaInd[inode];ijcrs++){
			const int jnode0 = (*m_BlkRowLevPtr)[ijcrs].row;
			const int lev0 = (*m_BlkRowLevPtr)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
		std::cout << "  \\" << inode << "\\  ";
		for(int ijcrs=m_BlkDiaInd[inode];ijcrs<m_BlkColInd[inode+1];ijcrs++){
			const int jnode0 = (*m_BlkRowLevPtr)[ijcrs].row;
			const int jlev0 = (*m_BlkRowLevPtr)[ijcrs].lev;
			std::cout << jnode0 << "-" << jlev0 << "   ";
		}
		std::cout << std::endl;
	}
*/
	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete[] const_cast<unsigned int*>(DiaInd_pre);
	delete pRowLev_pre;
	return true;
}


bool CMatDiaFrac_BlkCrs::AddFracPtn(const int lev_fill)
{
	assert( NBlkMatCol() == NBlkMatRow() );

	assert( m_ConditionFlag != -1 );
	if(  m_ConditionFlag == -1 ) return false;
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

	if( m_ncrs_Blk == 0 ){
		std::cout << "	Dia Center " << std::endl;
		return true;
	}

	if( m_ncrs_Blk == NBlkMatCol() * ( NBlkMatRow() - 1 ) ){
		std::cout << "   Full Matrix " << m_ncrs_Blk << " " << NBlkMatCol() << " " << NBlkMatRow() << std::endl;
        std::cout << m_ConditionFlag << std::endl;
		return true;
	}

	if( lev_fill == 0 ){
		std::cout << "    Fill Lev 0 " << std::endl;
		return true;
	}

    ////////////////

	if( m_ConditionFlag == 2 ){
		if( m_valCrs_Blk    != 0 ){ 	// 値配列が確保されている場合は確保し直す
            std::cout << "deallocate value list " << std::endl;
            delete[] m_valCrs_Blk; 
            m_valCrs_Blk = 0; 
        }
        assert( m_DiaInd != 0 );
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < NBlkMatRow() );
				if( jblk0 > iblk ){
					m_DiaInd[iblk] = icrs;
					break;
				}
			}
		}

		m_pRowLev = new std::vector<CRowLev>;
		m_pRowLev->resize( m_ncrs_Blk );
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
			(*m_pRowLev)[icrs].lev = 0;
		}
		delete[] m_rowPtr_Blk;
		m_rowPtr_Blk = 0;

		m_ConditionFlag = 1;
	}

	////////////////
	const unsigned int* ColInd_pre = m_colInd_Blk;  // もとのColInd
	const unsigned int* DiaInd_pre = m_DiaInd;      // もとのDiaInd
	const std::vector<CRowLev>* pRowLev_pre = m_pRowLev;

	const unsigned int ncrs_pre = m_ncrs_Blk;

	m_colInd_Blk = new unsigned int [NBlkMatCol()+1];   // 新しいColInd
	m_colInd_Blk[0] = 0;
	m_DiaInd = new unsigned int [NBlkMatCol()]; // 新しいDiaInd

	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*4);

    CListRLN nonzero(NBlkMatRow());
	nonzero.list.reserve(NBlkMatCol());
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
        nonzero.Initialize((*pRowLev_pre), ColInd_pre[iblk], ColInd_pre[iblk+1] );
		int knz_cur = nonzero.ipos_entry;
        if( knz_cur == -1 ){
            nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk,m_DiaInd);
            continue;
        }
		for(;;){
			const unsigned int kblk0 = nonzero.list[knz_cur].row;
			assert( kblk0<NBlkMatCol() );
			const int ik_lev0 = nonzero.list[knz_cur].lev;
			if( ik_lev0+1>lev_fill && lev_fill!=-1 ){
				knz_cur = nonzero.list[knz_cur].next;
				if( knz_cur == -1 ) break;
				continue;
			}
			if( kblk0 >= iblk ) break;
			int jnz_cur = knz_cur;
			for(unsigned int kjcrs=m_DiaInd[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){
				const int kj_lev0 = (*m_pRowLev)[kjcrs].lev;
				if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				const int unsigned jblk0 = (*m_pRowLev)[kjcrs].row;
				assert( jblk0>kblk0 && jblk0<NBlkMatRow() );
				assert( nonzero.list[jnz_cur].row < jblk0 );
				if( jblk0 == iblk ) continue;   // 対角にＣＲＳ成分は無いとは思うけど一応
				const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
				// check if this is fill in
                nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
			}
			knz_cur = nonzero.list[knz_cur].next;
			assert( (knz_cur>=0&&knz_cur<(int)NBlkMatRow()) || knz_cur==-1 );
			if( knz_cur == -1 ) break;
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk,m_DiaInd);
		////////////////
	}	// end iblk

	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete[] const_cast<unsigned int*>(DiaInd_pre);
	delete const_cast< std::vector<CRowLev>* >(pRowLev_pre);
    std::cout << "hoge" << m_ncrs_Blk << std::endl;
	/*
	for(unsigned int inode=0;inode<this->NBlkMatCol();inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = (*m_pRowLev)[ijcrs].row;
			const int lev0 = (*m_pRowLev)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
        std::cout << std::endl;
	}
    */
	return true;
}

bool CMatDiaFrac_BlkCrs::AddFracPtn(const int lev_fill, const std::vector<unsigned int>& aBlkFill){
	assert( NBlkMatCol() == NBlkMatRow() );

	assert( m_ConditionFlag != -1 );
	if(  m_ConditionFlag == -1 ) return false;
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

    ////////////////

	if( m_ncrs_Blk == 0 ){
		std::cout << "	Dia Center " << std::endl;
		return true;
	}

	if( m_ncrs_Blk == NBlkMatCol() * ( NBlkMatRow() - 1 ) ){
		std::cout << "   Full Matrix " << m_ncrs_Blk << " " << NBlkMatCol() << " " << NBlkMatRow() << std::endl;
        std::cout << m_ConditionFlag << std::endl;
		return true;
	}

	if( lev_fill == 0 && aBlkFill.empty() ){
		std::cout << "    Fill Lev 0 " << std::endl;
		return true;
	}

    ////////////////

	if( m_ConditionFlag == 2 ){
		if( m_valCrs_Blk    != 0 ){ 	// 値配列が確保されている場合は確保し直す
            std::cout << "deallocate value list " << std::endl;
            delete[] m_valCrs_Blk; 
            m_valCrs_Blk = 0; 
        }
        assert( m_DiaInd != 0 );
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < NBlkMatRow() );
				if( jblk0 > iblk ){
					m_DiaInd[iblk] = icrs;
					break;
				}
			}
		}

		m_pRowLev = new std::vector<CRowLev>;
		m_pRowLev->resize( m_ncrs_Blk );
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
			(*m_pRowLev)[icrs].lev = 0;
		}
		delete[] m_rowPtr_Blk;
		m_rowPtr_Blk = 0;

		m_ConditionFlag = 1;
	}

    const unsigned int nblk = this->NBlkMatCol();

    ////////////////
    std::vector<int> aFlgFill;
    aFlgFill.resize( nblk, lev_fill );
    for(unsigned int i=0;i<aBlkFill.size();i++){
        const unsigned int iblkfill = aBlkFill[i];
        assert( iblkfill < nblk );
        aFlgFill[iblkfill] = -1;
    }

	////////////////
	const unsigned int* ColInd_pre = m_colInd_Blk;  // もとのColInd
	const unsigned int* DiaInd_pre = m_DiaInd;      // もとのDiaInd
	const std::vector<CRowLev>* pRowLev_pre = m_pRowLev;

	const unsigned int ncrs_pre = m_ncrs_Blk;

	m_colInd_Blk = new unsigned int [nblk+1];   // 新しいColInd
	m_colInd_Blk[0] = 0;
	m_DiaInd = new unsigned int [nblk]; // 新しいDiaInd

	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*4);

    CListRLN nonzero(NBlkMatRow());
	nonzero.list.reserve(NBlkMatCol());
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
        const int levfill_i = aFlgFill[iblk];
        nonzero.Initialize((*pRowLev_pre), ColInd_pre[iblk], ColInd_pre[iblk+1] );
		int knz_cur = nonzero.ipos_entry;
        if( knz_cur == -1 ){
            nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk,m_DiaInd);
            continue;
        }
		for(;;){
			const unsigned int kblk0 = nonzero.list[knz_cur].row;
			assert( kblk0<NBlkMatCol() );
			const int ik_lev0 = nonzero.list[knz_cur].lev;
			if( kblk0 >= iblk ) break;
			int jnz_cur = knz_cur;
			for(unsigned int kjcrs=m_DiaInd[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){
				const int kj_lev0 = (*m_pRowLev)[kjcrs].lev;
				const int unsigned jblk0 = (*m_pRowLev)[kjcrs].row;
				assert( jblk0>kblk0 && jblk0<NBlkMatRow() );
				assert( nonzero.list[jnz_cur].row < jblk0 );
                const int levfill_j = aFlgFill[jblk0];
				if( jblk0 == iblk ) continue;   // 対角にＣＲＳ成分は無いとは思うけど一応
				const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                if( ((int)max_lev0+1>levfill_i || (int)max_lev0+1>levfill_j)
                    && (levfill_i!=-1 && levfill_j!=-1) ) continue;
//                std::cout << max_lev0 << " " << levfill_i << " " << levfill_j << std::endl;
				// check if this is fill in
                nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
			}
			knz_cur = nonzero.list[knz_cur].next;
			assert( (knz_cur>=0&&knz_cur<(int)NBlkMatRow()) || knz_cur==-1 );
			if( knz_cur == -1 ) break;
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk,m_DiaInd);
		////////////////
	}	// end iblk

	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete[] const_cast<unsigned int*>(DiaInd_pre);
	delete const_cast< std::vector<CRowLev>* >(pRowLev_pre);
    std::cout << "hoge : " << m_ncrs_Blk << std::endl;
/*	for(unsigned int inode=0;inode<this->NBlkMatCol();inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = (*m_pRowLev)[ijcrs].row;
			const int lev0 = (*m_pRowLev)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
        std::cout << std::endl;
	}*/
    return true;
}

bool CMatDiaFrac_BlkCrs::SetValue(const CMatDia_BlkCrs& rhs)
{
	assert( NBlkMatCol() == rhs.NBlkMatCol() );
	assert( NBlkMatRow() == rhs.NBlkMatRow() );
	assert( NBlkMatRow() == NBlkMatCol() );

	if( m_ConditionFlag == 0 ){	// no pattern & no value
		// make lev_fill 0 pattern
		CMatDia_BlkCrs::AddPattern(rhs,true);
        CMatDia_BlkCrs::SetZero();  // 値も確保しておく
		m_ConditionFlag = 2;
	}
	else if( m_ConditionFlag == 1 ){	// fill_lev & pattern & no value
		MakePatternFinalize();
		assert( m_ConditionFlag == 2 );
	}
    // 自分の行列の値に０をセットして，相手の行列の値をセットする
	this->SetValue_Initialize(rhs);
	this->DoILUDecomp();
	return true;
}

// 値をコピーしてＩＬＵ分解する
bool CMatDiaFrac_BlkCrs::SetValue(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order)
{
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( rhs.NBlkMatCol() == rhs.NBlkMatRow()   );
	assert( order.NBlk() == rhs.NBlkMatCol()   );
	assert( order.NBlk()   == this->NBlkMatCol() );

	assert( this->LenBlkCol() == this->LenBlkRow() );
	assert( rhs.LenBlkCol()   == rhs.LenBlkRow()   );
	assert( this->LenBlkCol() == rhs.LenBlkCol()   );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        return false;
    }

	if( m_ConditionFlag == 0 ){
		// make lev_fill 0 pattern
		CMatDia_BlkCrs::AddPattern(rhs,order);
		m_valCrs_Blk = new double [m_ncrs_Blk * LenBlkCol()*LenBlkRow()];
		m_ConditionFlag = 2;
	}
	else if( m_ConditionFlag == 1 ){
		MakePatternFinalize();
	}
    // 自分の行列の値に０をセットして，相手の行列の値をセットする
    this->SetValue_Initialize(rhs,order);
/*
	for(unsigned int inode=0;inode<this->NBlkMatCol();inode++){
        std::cout << inode << "-->" << std::endl;
        const unsigned int len = this->LenBlkCol();
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jblk0 = m_rowPtr_Blk[ijcrs];
			std::cout << "   " << jblk0 << "  ";
            for(unsigned int i=0;i<len*len;i++){
                std::cout << m_valCrs_Blk[ijcrs*len*len+i] << std::endl;
            }
            std::cout << std::endl;
		}
	}
*/
	DoILUDecomp();
	return true;
}

bool CMatDiaFrac_BlkCrs::MakePattern_Initialize(const CMatDia_BlkCrs& rhs)
{
    if( m_ConditionFlag == -1 ){
        if( rhs.LenBlkCol() == -1 ){
            const unsigned int nblk = rhs.NBlkMatCol();
            std::vector<unsigned int> aBlkSize;
            aBlkSize.resize(nblk);
            for(unsigned int iblk=0;iblk<nblk;iblk++){
                aBlkSize[iblk] = rhs.LenBlkCol(iblk);
            }
            CMatDia_BlkCrs::Initialize(nblk,aBlkSize);
        }
        else{
            CMatDia_BlkCrs::Initialize(rhs.NBlkMatCol(),rhs.LenBlkCol());
        }
        m_ConditionFlag = 0;
    }
    if( !CMatDia_BlkCrs::AddPattern(rhs,true) ){
        assert(0);
        return false;
    }
	assert( NBlkMatCol() == rhs.NBlkMatCol() );
	assert( NBlkMatRow() == rhs.NBlkMatRow() );
	assert( NBlkMatCol() == NBlkMatRow() );

	m_DiaInd = new unsigned int [NBlkMatCol()];
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
		m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			assert( icrs < m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[icrs];
			assert( jblk0 < NBlkMatRow() );
			if( jblk0 > iblk ){
				m_DiaInd[iblk] = icrs;
				break;
			}
		}
	}

	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->resize( m_ncrs_Blk );
	for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
		(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
		(*m_pRowLev)[icrs].lev = 0;
	}
	delete[] m_rowPtr_Blk;
	m_rowPtr_Blk = 0;

	m_ConditionFlag = 1;

	return true;
}

bool CMatDiaFrac_BlkCrs::MakePattern_Initialize(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order)
{
	assert( NBlkMatRow() == NBlkMatCol() );
	assert( rhs.NBlkMatCol() == rhs.NBlkMatRow()   );
	assert( order.NBlk() == rhs.NBlkMatCol()   );

	assert( this->LenBlkCol() == this->LenBlkRow() );
	assert( rhs.LenBlkCol()   == rhs.LenBlkRow()   );
	assert( this->LenBlkCol() == rhs.LenBlkCol()   );

	CMatDia_BlkCrs::AddPattern(rhs,order);

	m_DiaInd = new unsigned int [NBlkMatCol()];
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
		m_DiaInd[iblk] = m_colInd_Blk[iblk+1];
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			assert( icrs < m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[icrs];
			assert( jblk0 < NBlkMatRow() );
			if( jblk0 > iblk ){
				m_DiaInd[iblk] = icrs;
				break;
			}
		}
	}

	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->resize( m_ncrs_Blk );
	for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
		(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
		(*m_pRowLev)[icrs].lev = 0;
	}
	delete[] m_rowPtr_Blk;
	m_rowPtr_Blk = 0;

	m_ConditionFlag = 1;

	return true;
}

bool CMatDiaFrac_BlkCrs::MakePatternFinalize(){
	assert( NBlkMatRow() == NBlkMatCol() );

	assert( m_ConditionFlag != 0 );
    if(  m_ConditionFlag == 0 ){ return false; }

    assert( m_ConditionFlag != -1 );
    if(  m_ConditionFlag == -1 ){ return false; }

	if( m_ConditionFlag == 2 ){
		assert( m_pRowLev == 0 );
//		assert( m_rowPtr_Blk != 0 );
		return true;
	}

	m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
	assert( m_pRowLev->size() == m_ncrs_Blk );
	for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
		m_rowPtr_Blk[icrs] = (*m_pRowLev)[icrs].row;
	}

	delete m_pRowLev;
	m_pRowLev = 0;

    ////////////////
    // 値の確保
    
    unsigned int ntotdof = 0;
    if( this->LenBlkCol() >= 0 && this->LenBlkRow() >= 0 ){
	    ntotdof = NCrs()*LenBlkCol()*LenBlkRow();
    }
    else{   // 可変ブロックサイズの場合
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
        assert( m_DofPtrCol != 0 );
        assert( m_DofPtrRow != 0 );
        if( m_ValPtr == 0 ){    // m_ValPtrが作成されて無い場合は作る．
            m_ValPtr = new unsigned int [m_ncrs_Blk+1];
            for(unsigned int iblk=0;iblk<this->NBlkMatCol();iblk++){
                const unsigned int lencol = m_DofPtrCol[iblk+1] - m_DofPtrCol[iblk];
                for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
                    const unsigned int jblk0 = m_rowPtr_Blk[icrs];
                    const unsigned int lenrow = m_DofPtrRow[jblk0+1] - m_DofPtrRow[jblk0];
                    m_ValPtr[icrs] = ntotdof;
                    ntotdof += lencol*lenrow;
                }
            }
            m_ValPtr[m_ncrs_Blk] = ntotdof;
        }
	    ntotdof = m_ValPtr[m_ncrs_Blk];
    }
    
    if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; }
	// まだ値の領域をを確保していない場合
	m_valCrs_Blk = new double [ntotdof];

	m_ConditionFlag = 2;

	return true;
}
