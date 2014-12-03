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
// MatDiaCrsFrac.cpp: CMatDiaCrsFrac クラスのインプリメンテーション
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

#include "delfem/matvec/matfrac_blkcrs.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
//#include "ker_mat.h"

using namespace MatVec;


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




//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CMatFrac_BlkCrs::CMatFrac_BlkCrs(unsigned int nblk_col, unsigned int len_col,
								 unsigned int nblk_row, unsigned int len_row)
	:CMat_BlkCrs(nblk_col,len_col,  nblk_row,len_row)
{
    assert( m_colInd_Blk != 0 );
    assert( m_ncrs_Blk == 0 );
    m_ConditionFlag = 2;
	m_pRowLev = 0;
}

CMatFrac_BlkCrs::CMatFrac_BlkCrs(const CMat_BlkCrs& rhs)
{
    if( rhs.LenBlkCol() >= 0 && rhs.LenBlkRow() >= 0 ){
        this->Initialize(rhs.NBlkMatCol(),rhs.LenBlkCol(),  rhs.NBlkMatRow(),rhs.LenBlkRow());
    }
    else{
        const unsigned int nblk_col = rhs.NBlkMatCol();
        std::vector<unsigned int> aLen_col;
        aLen_col.resize(nblk_col);
        for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            aLen_col[iblk] = rhs.LenBlkCol(iblk);
        }
        const unsigned int nblk_row = rhs.NBlkMatRow();
        std::vector<unsigned int> aLen_row;
        aLen_row.resize(nblk_row);
        for(unsigned int iblk=0;iblk<nblk_row;iblk++){
            aLen_row[iblk] = rhs.LenBlkRow(iblk);
        }
        this->Initialize(nblk_col,aLen_col,  nblk_row,aLen_row);
    }

    ////////////////
	m_pRowLev = 0;

	CMat_BlkCrs::AddPattern(rhs,true);
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
	m_valCrs_Blk = 0;
	m_ConditionFlag = 2;
}





bool CMatFrac_BlkCrs::MakePatternInitialize(const CMat_BlkCrs& rhs)
{
	CMat_BlkCrs::AddPattern(rhs,true);

	assert( NBlkMatCol() == rhs.NBlkMatCol() );
	assert( NBlkMatRow() == rhs.NBlkMatRow() );

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

bool CMatFrac_BlkCrs::AddFracPtnLowUp(const CMatFrac_BlkCrs& mat_low, const CMatFrac_BlkCrs& mat_up,  const int lev_fill )
{
	assert( mat_low.NBlkMatCol() == NBlkMatCol() );
	assert( mat_up.NBlkMatRow() == NBlkMatRow() );
	assert( mat_up.NBlkMatCol() == mat_low.NBlkMatRow() );

	assert( m_ConditionFlag != -1 );
    if(  m_ConditionFlag == -1 ) return false;

	if( m_ncrs_Blk == NBlkMatCol() * NBlkMatRow() ){
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
		m_pRowLev = new std::vector<CRowLev>;
		m_pRowLev->resize( m_ncrs_Blk );
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
			(*m_pRowLev)[icrs].lev = 0;
		}
		if( m_rowPtr_Blk != 0 ) delete[] m_rowPtr_Blk;
		m_rowPtr_Blk = 0;

		m_ConditionFlag = 1;
	}

    ////////////////////////////////////////////////////////////////
	const unsigned int* ColInd_pre = this->m_colInd_Blk;
	const std::vector<CRowLev>* pRowLev_pre = this->m_pRowLev;

	const unsigned int nnode = NBlkMatCol();
	const unsigned int ncrs_pre = NCrs();

	m_colInd_Blk = new unsigned int [nnode+1];
	m_colInd_Blk[0] = 0;
	assert( m_rowPtr_Blk == 0 );
	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*2);

    ////////////////////////////////
    CListRLN nonzero(NBlkMatRow());
	for(unsigned int iblk=0;iblk<nnode;iblk++){
        nonzero.Initialize((*pRowLev_pre), ColInd_pre[iblk], ColInd_pre[iblk+1] );
		for(unsigned int ikcrs_low=mat_low.m_colInd_Blk[iblk];ikcrs_low<mat_low.m_colInd_Blk[iblk+1];ikcrs_low++){
			assert( ikcrs_low>=0 && ikcrs_low<mat_low.NCrs() );
            const unsigned int kblk0 = mat_low.GetRow(ikcrs_low);
			assert( kblk0>=0 && kblk0<mat_low.NBlkMatRow() );
			const int ik_lev0 = mat_low.GetLevel(ikcrs_low);
			if( ik_lev0+1>lev_fill && lev_fill!=-1 ) continue;
            int jnz_cur = nonzero.ipos_entry;
		    for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){
			    const int kj_lev0 = mat_up.GetLevel(kjcrs);
			    if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
			    const unsigned int jblk0 = mat_up.GetRow(kjcrs);
			    assert( jblk0<NBlkMatRow() );
			    const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
		    }
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk);
	}

/*
	for(unsigned int inode=0;inode<nnode;inode++){
        std::cout << inode << "  " << m_colInd_Blk[inode] << std::endl;
    }
	for(unsigned int inode=0;inode<nnode;inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = (*m_pRowLev)[ijcrs].row;
			const int lev0 = (*m_pRowLev)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
		std::cout << std::endl;
	}
*/
	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete pRowLev_pre;
	return true;
}

// 右上の非対角ブロックの行列にパターンを追加
bool CMatFrac_BlkCrs::AddFracPtnUp(const CMatDiaFrac_BlkCrs& mat_dia, const int lev_fill)
{
    ////////////////
	assert( mat_dia.NBlkMatCol() == NBlkMatCol() );
	assert( mat_dia.LenBlkCol() == LenBlkCol() );

	if( this->NCrs() == NBlkMatCol() * NBlkMatRow() ){
		std::cout << "   Full Matrix " << std::endl;
		return true;
	}

	if( lev_fill == 0 ){
		std::cout << "	Fill Lev 0 " << std::endl;
		return true;
	}

	assert( m_ConditionFlag != -1 );
	if(  m_ConditionFlag == -1 ) return false;
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

    ////////////////

	if( m_ConditionFlag == 2 ){
		if( m_valCrs_Blk    != 0 ){ 	// 値配列が確保されている場合は確保し直す
            std::cout << "deallocate value list " << std::endl;
            delete[] m_valCrs_Blk; 
            m_valCrs_Blk = 0; 
        }

		m_pRowLev = new std::vector<CRowLev>;
		m_pRowLev->resize( m_ncrs_Blk );
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			(*m_pRowLev)[icrs].row = m_rowPtr_Blk[icrs];
			(*m_pRowLev)[icrs].lev = 0;
		}
        if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; }
		m_rowPtr_Blk = 0;

		m_ConditionFlag = 1;
	}

    ////////////////
	const unsigned int* ColInd_pre = this->m_colInd_Blk;
	const std::vector<CRowLev>* pRowLev_pre = this->m_pRowLev;

	const unsigned int nblk_col = NBlkMatCol();
	const unsigned int nblk_row = NBlkMatRow();
	const unsigned int ncrs_pre = NCrs();

	m_colInd_Blk = new unsigned int [nblk_col+1];
	m_colInd_Blk[0] = 0;
	assert( m_rowPtr_Blk == 0 );
	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*2);

    CListRLN nonzero( this->NBlkMatRow() );
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        assert( ColInd_pre[iblk+1] - ColInd_pre[iblk] <= nblk_row );
        nonzero.Initialize(*pRowLev_pre,ColInd_pre[iblk],ColInd_pre[iblk+1]);
		if( mat_dia.m_ConditionFlag == 1 ){
            for(unsigned int ikcrs_dia=mat_dia.m_colInd_Blk[iblk];ikcrs_dia<mat_dia.m_DiaInd[iblk];ikcrs_dia++){ assert( ikcrs_dia<mat_dia.NCrs() );
				const unsigned int kblk0 = (*mat_dia.m_pRowLev)[ikcrs_dia].row; assert( kblk0<mat_dia.NBlkMatRow() );
				const int ik_lev0 = (*mat_dia.m_pRowLev)[ikcrs_dia].lev;
				if( ik_lev0+1>lev_fill && lev_fill!=-1 ) continue;
                int jnz_cur = nonzero.ipos_entry;
                for(unsigned int kjcrs=m_colInd_Blk[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<NCrs() );
					const unsigned int jblk0 = (*m_pRowLev)[kjcrs].row; assert( jblk0<nblk_row );
					const int kj_lev0 = (*m_pRowLev)[kjcrs].lev;
					if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
					const int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
				}		
			}
		}
		else if( mat_dia.m_ConditionFlag == 2 ){
            for(unsigned int ikcrs_dia=mat_dia.m_colInd_Blk[iblk];ikcrs_dia<mat_dia.m_DiaInd[iblk];ikcrs_dia++){ assert( ikcrs_dia<mat_dia.NCrs() );
				const unsigned int kblk0 = mat_dia.m_rowPtr_Blk[ikcrs_dia]; assert( kblk0<iblk );
				const int ik_lev0 = 0;
				if( ik_lev0+1>lev_fill && lev_fill!=-1 ) continue;
                int jnz_cur = nonzero.ipos_entry;
				for(unsigned int kjcrs=m_colInd_Blk[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){
					assert( kjcrs>=0 && kjcrs<NCrs() );
					const unsigned int jblk0 = (*m_pRowLev)[kjcrs].row; assert( jblk0<nblk_row );
					const int kj_lev0 = (*m_pRowLev)[kjcrs].lev;
					if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
					const int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
				}		
			}
		}
		else{
			assert(0);
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk);
	}
/*    
	for(unsigned int inode=0;inode<nblk_col;inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = (*m_pRowLev)[ijcrs].row;
			const int lev0 = (*m_pRowLev)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
        std::cout << std::endl;
	}
*/
	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete pRowLev_pre;
	return true;
}

// 左下の行列に対してパターンを追加する
bool CMatFrac_BlkCrs::AddFracPtnLow(const CMatDiaFrac_BlkCrs& mat_dia, const int lev_fill)
{
	assert( mat_dia.NBlkMatRow() == NBlkMatRow() );
	assert( mat_dia.LenBlkRow() == LenBlkRow() );

	assert( m_ConditionFlag != -1 );
	if(  m_ConditionFlag == -1 ) return false;
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

	if( m_ncrs_Blk == NBlkMatCol() * NBlkMatRow() ){
		std::cout << "   Full Matrix " << std::endl;
		return true;
	}

	if( lev_fill == 0 ){
		std::cout << "	Fill Lev 0 " << std::endl;
		return true;
	}

    ////////////////

	if( m_ConditionFlag == 2 ){
		if( m_valCrs_Blk    != 0 ){ 	// 値配列が確保されている場合は確保し直す
            std::cout << "deallocate value list " << std::endl;
            delete[] m_valCrs_Blk; 
            m_valCrs_Blk = 0; 
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
	const unsigned int* ColInd_pre = this->m_colInd_Blk;
	const std::vector<CRowLev>* pRowLev_pre = this->m_pRowLev;

	const unsigned int nblk_col = NBlkMatCol();
	const unsigned int nblk_row = NBlkMatRow();
	const unsigned int ncrs_pre = NCrs();

	m_colInd_Blk = new unsigned int [nblk_col+1];
	m_colInd_Blk[0] = 0;
	assert( m_rowPtr_Blk == 0 );
	m_pRowLev = new std::vector<CRowLev>;
	m_pRowLev->reserve(ncrs_pre*2);

    CListRLN nonzero( this->NBlkMatRow() );
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        nonzero.Initialize(*pRowLev_pre, ColInd_pre[iblk], ColInd_pre[iblk+1] );
        int knz_cur = nonzero.ipos_entry;
        if( knz_cur == -1 ){
            nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk);
            continue;
        }
		if( mat_dia.m_ConditionFlag == 1 ){
			for(;;){
                assert( knz_cur >= 0 && knz_cur < (int)nonzero.list.size() );
				const unsigned int kblk0 = nonzero.list[knz_cur].row; assert( kblk0<NBlkMatRow() );
				const int ik_lev0 = nonzero.list[knz_cur].lev;
			    if( ik_lev0+1>lev_fill && lev_fill!=-1 ){
				    knz_cur = nonzero.list[knz_cur].next;
				    if( knz_cur == -1 ) break;
				    continue;
			    }
                int jnz_cur = knz_cur;
                for(unsigned int kjcrs_dia=mat_dia.m_DiaInd[kblk0];kjcrs_dia<mat_dia.m_colInd_Blk[kblk0+1];kjcrs_dia++){ assert( kjcrs_dia<mat_dia.NCrs() );
					const unsigned int jblk0 = (*mat_dia.m_pRowLev)[kjcrs_dia].row; assert( jblk0<nblk_row && jblk0>kblk0 );
					const int kj_lev0 = (*mat_dia.m_pRowLev)[kjcrs_dia].lev;
					if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				    const int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jblk0,max_lev0+1,jnz_cur);
				}		
			    knz_cur = nonzero.list[knz_cur].next;
			    assert( (knz_cur>=0&&knz_cur<(int)NBlkMatRow()) || knz_cur==-1 );
	    		if( knz_cur == -1 ) break;
    		}
		}
		else if( mat_dia.m_ConditionFlag == 2 ){
			for(;;){
                assert( knz_cur >= 0 && knz_cur < (int)nonzero.list.size() );
				const unsigned int kblk0 = nonzero.list[knz_cur].row; assert( kblk0<NBlkMatRow() );
				const int ik_lev0 = nonzero.list[knz_cur].lev;
			    if( ik_lev0+1>lev_fill && lev_fill!=-1 ){
				    knz_cur = nonzero.list[knz_cur].next;
				    if( knz_cur == -1 ) break;
				    continue;
			    }
                int jnz_cur = knz_cur;
                for(unsigned int kjcrs_dia=mat_dia.m_DiaInd[kblk0];kjcrs_dia<mat_dia.m_colInd_Blk[kblk0+1];kjcrs_dia++){ assert( kjcrs_dia<mat_dia.NCrs() );
					const unsigned int jnode0 = mat_dia.m_rowPtr_Blk[kjcrs_dia]; assert( jnode0<nblk_row && jnode0>kblk0 );
					const int kj_lev0 = 0;
					if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				    const int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
                    nonzero.Insert(jnode0,max_lev0+1,jnz_cur);
				}		
			    knz_cur = nonzero.list[knz_cur].next;
			    assert( (knz_cur>=0&&knz_cur<(int)NBlkMatRow()) || knz_cur==-1 );
	    		if( knz_cur == -1 ) break;
    		}
		}
		else{
			assert(0);
		}
        nonzero.Finalize(iblk,m_colInd_Blk,*m_pRowLev,m_ncrs_Blk);
	}
/*
	for(unsigned int inode=0;inode<nblk_col;inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = (*m_pRowLev)[ijcrs].row;
			const int lev0 = (*m_pRowLev)[ijcrs].lev;
			std::cout << jnode0 << "-" << lev0 << "   ";
		}
        std::cout << std::endl;
	}
*/
	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete pRowLev_pre;
	return true;


	return true;
}

bool CMatFrac_BlkCrs::MakePatternFinalize()
{
	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ){
		return false;
	}

	if( m_ConditionFlag == 2 ){
		assert( m_pRowLev == 0 );
		assert( m_rowPtr_Blk != 0 );
		return true;
	}

	m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
	assert( m_pRowLev->size() == m_ncrs_Blk );
	for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
		m_rowPtr_Blk[icrs] = (*m_pRowLev)[icrs].row;
	}

	delete m_pRowLev;
	m_pRowLev = 0;

    assert( m_ValPtr == 0 );
	assert( m_valCrs_Blk == 0 );

	m_ConditionFlag = 2;
/*
    for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
        std::cout << m_colInd_Blk[iblk] << std::endl;
    }
	for(unsigned int inode=0;inode<this->NBlkMatCol();inode++){
		std::cout << inode << "-->";
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			const int jnode0 = m_rowPtr_Blk[ijcrs];
			std::cout << jnode0 << "  ";
		}
        std::cout << std::endl;
	}
*/
	return true;
}


bool CMatFrac_BlkCrs::DoILUDecompLowUp( 
	const CMatFrac_BlkCrs& mat_low, const CMatFrac_BlkCrs& mat_up )
{
    assert( mat_up.NBlkMatCol() == mat_low.NBlkMatRow() );
	assert( mat_low.NBlkMatCol() == NBlkMatCol() );
	assert( mat_up.NBlkMatRow() == NBlkMatRow() );

	const unsigned int nblk_col = NBlkMatCol();
	const unsigned int nblk_row = NBlkMatRow();
	const unsigned int nblk_dia = mat_up.NBlkMatCol();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
        assert(0);
        return true;
    }

//	assert( mat_low.LenBlkCol() == LenBlkCol() );
//  assert( mat_up.LenBlkRow() == LenBlkRow() );

	const unsigned int len_col = LenBlkCol();
	const unsigned int len_row = LenBlkRow();

	int* row2crs_f = new int [nblk_row];
	for(unsigned  jblk=0;jblk<nblk_row;jblk++){ row2crs_f[jblk] = -1; }
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ 
            assert( ijcrs<this->m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			row2crs_f[jblk0] = ijcrs;
		}
        assert( mat_low.LenBlkCol(iblk) == len_col );
        for(unsigned int ikcrs=mat_low.m_colInd_Blk[iblk];ikcrs<mat_low.m_colInd_Blk[iblk+1];ikcrs++){ assert( ikcrs<mat_low.m_ncrs_Blk );
			const double* ikvalue = mat_low.GetValCrsPtr(ikcrs);
			const unsigned int kblk0 = mat_low.m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk_dia );
            const unsigned int len_dia = mat_low.LenBlkRow(kblk0);
            assert( mat_up.LenBlkCol( kblk0) == len_dia );
            for(unsigned int kjcrs=mat_up.m_colInd_Blk[kblk0];kjcrs<mat_up.m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<mat_up.m_ncrs_Blk ); 
				const double* kjvalue = mat_up.GetValCrsPtr(kjcrs);
				const unsigned int jblk0 = mat_up.m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk_row );
                assert( mat_up.LenBlkRow(jblk0) == len_row );
				const int ijcrs0 = row2crs_f[jblk0];
                if( ijcrs0 == -1 ){ continue; }
                assert( ijcrs0>=0 && ijcrs0<(int)m_ncrs_Blk );
				double* ijvalue = &m_valCrs_Blk[ijcrs0*len_col*len_row];
				CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_col,len_dia,len_row);	// Aij = Aij-Aik*Akk*Akj
			}		
		}
        for(unsigned ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			row2crs_f[jblk0] = -1;
		}
		////////////////
	}	// inode
	delete[] row2crs_f;

	return true;
}



bool CMatFrac_BlkCrs::DoILUDecompUp(const CMatDiaFrac_BlkCrs& mat_dia)
{
//    std::cout << "DoILUDecompUp" << std::endl;
	assert( NBlkMatCol() == mat_dia.NBlkMatCol() );
	assert( LenBlkCol() == mat_dia.LenBlkCol() );

    const unsigned int nblk_col = NBlkMatCol();
    const unsigned int nblk_row = NBlkMatRow();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        double* tmp_buffer;
        {
            unsigned int max_col=0, max_row=0;
            for(unsigned int iblk=0;iblk<nblk_col;iblk++){ max_col = ( LenBlkCol(iblk) > max_col ) ? LenBlkCol(iblk) : max_col; }
            for(unsigned int iblk=0;iblk<nblk_row;iblk++){ max_row = ( LenBlkRow(iblk) > max_row ) ? LenBlkRow(iblk) : max_row; }
	        tmp_buffer = new double [max_col*max_row];	// 小さなバッファ
        }
	    int* row2crs_f = new int [nblk_row];
	    for(unsigned int jblk=0;jblk<nblk_row;jblk++){ row2crs_f[jblk] = -1; }
	    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<this->m_ncrs_Blk );
			    const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			    row2crs_f[jblk0] = ijcrs;
		    }
            const unsigned int len_i = this->LenBlkCol(iblk);
            assert( mat_dia.LenBlkCol(iblk) == len_i );
            for(unsigned int ikcrs=mat_dia.m_colInd_Blk[iblk];ikcrs<mat_dia.m_DiaInd[iblk];ikcrs++){ assert( ikcrs<mat_dia.m_ncrs_Blk );
			    const double* ikvalue = mat_dia.GetValCrsPtr(ikcrs);
			    const unsigned int kblk0 = mat_dia.m_rowPtr_Blk[ikcrs]; assert( kblk0<iblk );
                const unsigned int len_k = mat_dia.LenBlkRow(kblk0);
                for(unsigned int kjcrs=m_colInd_Blk[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<m_ncrs_Blk );
				    const double* kjvalue = &m_valCrs_Blk[ m_ValPtr[kjcrs] ];
				    const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk_row );
                    const unsigned int len_j = this->LenBlkRow(jblk0);
				    const int ijcrs0 = row2crs_f[jblk0];
                    if( ijcrs0 == -1 ){ continue; }
                    assert( ijcrs0>=0 && ijcrs0<(int)m_ncrs_Blk );
				    double* ijvalue = &m_valCrs_Blk[ m_ValPtr[ijcrs0] ];
				    CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_i,len_k,len_j);
			    }
		    }
		    const double* iivalue = mat_dia.GetValDiaPtr(iblk);
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
                const unsigned int len_j = this->LenBlkRow(jblk0);
			    double* ijvalue = &m_valCrs_Blk[ m_ValPtr[ijcrs] ];
			    CalcMatPr(ijvalue,iivalue,tmp_buffer,len_i,len_j);
		    }
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
                const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0>=0 && jblk0<nblk_row );
			    row2crs_f[jblk0] = -1;
		    }
	    }	// end inode

	    delete[] tmp_buffer;
	    delete[] row2crs_f;
        return true;
    }


	const int len_col = LenBlkCol();
	const int len_row = LenBlkRow();
//	const int ncrs_pre = NCrs();

	double* tmp_buffer = new double [len_col*len_row];	// 小さなバッファ
	int* row2crs_f = new int [nblk_row];
	for(unsigned int jblk=0;jblk<nblk_row;jblk++){ row2crs_f[jblk] = -1; }
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<this->m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			row2crs_f[jblk0] = ijcrs;
		}
        for(unsigned int ikcrs=mat_dia.m_colInd_Blk[iblk];ikcrs<mat_dia.m_DiaInd[iblk];ikcrs++){ assert( ikcrs<mat_dia.m_ncrs_Blk );
			const double* ikvalue = &mat_dia.m_valCrs_Blk[ikcrs*len_col*len_col];
			const unsigned int kblk0 = mat_dia.m_rowPtr_Blk[ikcrs]; assert( kblk0<iblk );
            for(unsigned int kjcrs=m_colInd_Blk[kblk0];kjcrs<m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<m_ncrs_Blk );
				const double* kjvalue = &m_valCrs_Blk[kjcrs*len_col*len_row];
				const unsigned int jblk0 = m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk_row );
				const int ijcrs0 = row2crs_f[jblk0];
                if( ijcrs0 == -1 ){ continue; }
                assert( ijcrs0>=0 && ijcrs0<(int)m_ncrs_Blk );
				double* ijvalue = &m_valCrs_Blk[ijcrs0*len_col*len_row];
				CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_col,len_col,len_row);
			}
		}
		const double* iivalue = &mat_dia.m_valDia_Blk[iblk*len_col*len_col];
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
			double* ijvalue = &m_valCrs_Blk[ijcrs*len_col*len_row];
			CalcMatPr(ijvalue,iivalue,tmp_buffer,len_col,len_row);
		}
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<m_ncrs_Blk );
            const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0>=0 && jblk0<nblk_row );
			row2crs_f[jblk0] = -1;
		}
	}	// end inode

	delete[] tmp_buffer;
	delete[] row2crs_f;

	return true;
}

bool CMatFrac_BlkCrs::DoILUDecompLow(const CMatDiaFrac_BlkCrs& mat_dia)
{
//    std::cout << "DoILUDecompLow" << std::endl;
	assert( NBlkMatRow() == mat_dia.NBlkMatRow() );
	assert( LenBlkRow() == mat_dia.LenBlkRow() );

    const unsigned int nblk_col = NBlkMatCol();
    const unsigned int nblk_row = NBlkMatRow();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 )
    {
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
	    int* row2crs_f = new int [nblk_row];
        for(unsigned int jblk=0;jblk<nblk_row;jblk++){ row2crs_f[jblk] = -1; }
	    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            const unsigned int len_i = this->LenBlkCol(iblk);
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<NCrs() );
			    const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			    row2crs_f[jblk0] = ijcrs;
		    }
            for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_colInd_Blk[iblk+1];ikcrs++){ assert( ikcrs<NCrs() );
			    const double* ikvalue = &m_valCrs_Blk[ m_ValPtr[ikcrs] ];
			    const unsigned int kblk0 = m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk_row );
                const unsigned int len_k = this->LenBlkRow(kblk0);
                assert( mat_dia.LenBlkCol(kblk0) == len_k );
                for(unsigned int kjcrs=mat_dia.m_DiaInd[kblk0];kjcrs<mat_dia.m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<mat_dia.NCrs() );
				    const double* kjvalue = mat_dia.GetValCrsPtr(kjcrs);
				    const unsigned int jblk0 = mat_dia.m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk_row && jblk0>kblk0 );
				    const int ijcrs0 = row2crs_f[jblk0];
                    if( ijcrs0 == -1 ){ continue; }
                    assert( ijcrs0>=0 && ijcrs0<(int)NCrs() );
                    const unsigned int len_j = this->LenBlkRow(jblk0);
                    assert( mat_dia.LenBlkRow(jblk0) );
				    double* ijvalue = &m_valCrs_Blk[ m_ValPtr[ijcrs0] ];
				    CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_i,len_k,len_j);
			    }
		    }
            for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<NCrs() );
			    const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0>=0 && jblk0<nblk_row );
			    row2crs_f[jblk0] = -1;
		    }
	    }	// end inode
	    delete[] row2crs_f;
        return true;
    }

	const int len_col = LenBlkCol();
	const int len_row = LenBlkRow();
    assert( len_col >= 0 && len_row >= 0 );
//    const unsigned int ncrs_pre = NCrs();

	int* row2crs_f = new int [nblk_row];
    for(unsigned int jblk=0;jblk<nblk_row;jblk++){ row2crs_f[jblk] = -1; }
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<NCrs() );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0<nblk_row );
			row2crs_f[jblk0] = ijcrs;
		}
        for(unsigned int ikcrs=m_colInd_Blk[iblk];ikcrs<m_colInd_Blk[iblk+1];ikcrs++){ assert( ikcrs<NCrs() );
			const double* ikvalue = &m_valCrs_Blk[ikcrs*len_col*len_row];
			const unsigned int kblk0 = m_rowPtr_Blk[ikcrs]; assert( kblk0<nblk_row );
            for(unsigned int kjcrs=mat_dia.m_DiaInd[kblk0];kjcrs<mat_dia.m_colInd_Blk[kblk0+1];kjcrs++){ assert( kjcrs<mat_dia.NCrs() );
				const double* kjvalue = mat_dia.GetValCrsPtr(kjcrs);
				const unsigned int jblk0 = mat_dia.m_rowPtr_Blk[kjcrs]; assert( jblk0<nblk_row && jblk0>kblk0 );
				const int ijcrs0 = row2crs_f[jblk0];
                if( ijcrs0 == -1 ){ continue; }
                assert( ijcrs0>=0 && ijcrs0<(int)NCrs() );
				double* ijvalue = &m_valCrs_Blk[ijcrs0*len_col*len_row];
				CalcSubMatPr(ijvalue,ikvalue,kjvalue,len_col,len_row,len_row);
			}
		}
        for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){ assert( ijcrs<NCrs() );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs]; assert( jblk0>=0 && jblk0<nblk_row );
			row2crs_f[jblk0] = -1;
		}
		////////////////
	}	// end inode
	delete[] row2crs_f;

	return true;
}
