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
// MatDia_Crs.cpp: Implementation of class "CMatDia_BlkCrs"
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0); else for

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>
#include <algorithm>

#include "delfem/indexed_array.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/matfrac_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/ordering_blk.h"
#include "delfem/matvec/bcflag_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CMatDia_BlkCrs::CMatDia_BlkCrs() : m_valDia_Blk(0), m_DiaValPtr(0)
{
}

CMatDia_BlkCrs::CMatDia_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow)
:CMat_BlkCrs(nblk_colrow, len_colrow, nblk_colrow, len_colrow)
{
//	std::cout << "Construct : CMatDia_BlkCrs(nblk_colrow, len_colrow) " << nblk_colrow << " " << len_colrow << std::endl;
	m_valDia_Blk = new  double [NBlkMatCol() * LenBlkCol()*LenBlkRow()];
    m_DiaValPtr = 0;
}

CMatDia_BlkCrs::~CMatDia_BlkCrs()
{
	if( m_valDia_Blk != 0 ){ delete[] m_valDia_Blk; m_valDia_Blk = 0; }
	if( m_DiaValPtr  != 0 ){ delete[] m_DiaValPtr;  m_DiaValPtr  = 0; }
}


    
bool CMatDia_BlkCrs::Initialize(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
                                unsigned int nblk_row, const std::vector<unsigned int>& alen_row)
{
    // 入力のチェック
    if( nblk_col != nblk_row ){ assert(0); return false; }
    const unsigned int nblk = nblk_col;
    if( alen_col.size() != nblk ){ assert(0); return false; }
    if( alen_row.size() != nblk ){ assert(0); return false; }
    for(unsigned int iblk=0;iblk<nblk;iblk++){
        if( alen_col[iblk] != alen_row[iblk] ){ assert(0); return false; }
    }
    
    // 親クラスのInitialize
    if( !CMat_BlkCrs::Initialize(nblk_col, alen_col, nblk_row,alen_row) ){
        assert(0);  // 何が起きるか分からんのでとりあえずassertion
        return false;
    }
    assert( this->LenBlkCol() == -1 );
    assert( this->LenBlkRow() == -1 );
    
    {
        const unsigned int nblk = this->NBlkMatCol();
        if( m_DiaValPtr != 0 ){ delete[] m_DiaValPtr; }
        m_DiaValPtr = new unsigned int [nblk+1];
        m_DiaValPtr[0] = 0;
        assert( this->m_DofPtrCol != 0 );
        unsigned int ni = 0;
        for(unsigned iblk=0;iblk<nblk;iblk++){
            const unsigned int n = this->LenBlkCol(iblk);
            assert( m_DofPtrCol[iblk] == m_DofPtrRow[iblk] );
            ni += n*n;
            m_DiaValPtr[iblk+1] = ni;
        }
        if( m_valDia_Blk != 0 ){ delete[] m_valDia_Blk; }
        m_valDia_Blk = new double [ni];
    }
    return true;
}


bool CMatDia_BlkCrs::Initialize(unsigned int nblk_col, unsigned int len_col,
                                unsigned int nblk_row, unsigned int len_row )
{
    // 入力のチェック
    if( nblk_col != nblk_row ){ assert(0); return false; }
//    const unsigned int nblk = nblk_col;
    if( len_col != len_row ){ assert(0); return false; }
//    const unsigned int len = len_col;

    
    // 親クラスのInitialize
    if( !CMat_BlkCrs::Initialize(nblk_col, len_col, nblk_row,len_row) ){
        assert(0);  // 何が起きるか分からんのでとりあえずassertion
        return false;
    }
    
    assert( this->m_DofPtrCol == 0 );
    assert( this->m_DofPtrRow == 0 );
    this->m_DiaValPtr = 0;
    if( m_valDia_Blk != 0 ){ delete[] m_valDia_Blk; }
    m_valDia_Blk = new double [nblk_col*len_col*len_row];
    
    return true;
}


bool CMatDia_BlkCrs::DeletePattern(){
	CMat_BlkCrs::DeletePattern();
	return true;
}

bool CMatDia_BlkCrs::SetZero()
{
    assert( this->NBlkMatCol() == this->NBlkMatRow() );
	CMat_BlkCrs::SetZero();
    unsigned int ni = 0;
    if( LenBlkCol() >= 0 && LenBlkRow() >= 0 ){
    	ni = NBlkMatCol()*LenBlkCol()*LenBlkCol();
    }
    else{
        if( LenBlkCol() >=0 || LenBlkRow() >= 0 ){
            std::cout << "Error! -> Not Implemented" << std::endl;
            assert(0);
            return false;
        }
        assert( this->m_ValPtr !=  0 );
        const unsigned int nblk = this->NBlkMatCol();
        assert( this->m_DiaValPtr != 0 );
        assert( this->m_DofPtrCol != 0 );
        assert( this->m_DofPtrRow != 0 );
        ni = m_DiaValPtr[nblk];
    }
    assert( m_valDia_Blk != 0 );    // コンストラクタやInitailizeの時に何かは確保されているはず
	for(unsigned int i=0;i<ni;i++){ m_valDia_Blk[i] = 0.0; }
	return true;
}

bool CMatDia_BlkCrs::Mearge(unsigned int nblkel_col, const unsigned int* blkel_col,
						    unsigned int nblkel_row, const unsigned int* blkel_row,
						    unsigned int blksize, const double* emat)
{
    assert( m_colInd_Blk != 0 );
	assert( m_valCrs_Blk != 0 );
	assert( m_valDia_Blk != 0 );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
        assert( nblkel_col == 1 );
        assert( nblkel_row == 1 );
        const unsigned int iblk0 = blkel_col[0]; assert( iblk0 < this->NBlkMatCol() );
        const unsigned int jblk0 = blkel_row[0]; assert( jblk0 < this->NBlkMatRow() );
        if( iblk0 == jblk0 ){
            const unsigned int len = this->LenBlkCol(iblk0);
            assert( blksize == len*len );
            const unsigned int ivalptr0 = this->m_DiaValPtr[iblk0];
            for(unsigned int i=0;i<blksize;i++){
                this->m_valDia_Blk[ivalptr0+i] += emat[i];
            }
            return true;
        }
        const unsigned int lencol = this->LenBlkCol(iblk0);
        const unsigned int lenrow = this->LenBlkRow(jblk0);
        assert( blksize == lencol*lenrow );
        int icrs0 = -1;
        for(unsigned int icrs=m_colInd_Blk[iblk0];icrs<m_colInd_Blk[iblk0+1];icrs++){
            const unsigned int jblk1 = m_rowPtr_Blk[icrs];
            if( jblk1 == jblk0 ){
                icrs0 = icrs;
                break;
            }
        }
        assert( icrs0 != -1 );
        const unsigned int ivalptr0 = this->m_ValPtr[icrs0];
        for(unsigned int i=0;i<blksize;i++){
            this->m_valCrs_Blk[ivalptr0+i] += emat[i];
        }
        return true;
    }

    if( m_marge_tmp_buffer == 0 ){
        const unsigned int nblkrow = NBlkMatRow();
        m_marge_tmp_buffer = new int [nblkrow];
        for(unsigned int iblk=0;iblk<nblkrow;iblk++){
            m_marge_tmp_buffer[iblk] = -1;
        }
    }

	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
//	const unsigned int BlkLen = LenBlkCol();

	assert( nblkel_col == nblkel_row );
	assert( blksize == BlkSize );

	const unsigned int* colind = m_colInd_Blk;
	const unsigned int* rowptr = m_rowPtr_Blk;
	double* matval_nd = m_valCrs_Blk;
	double* matval_dia = m_valDia_Blk;

	for(unsigned int iblkel=0;iblkel<nblkel_col;iblkel++){
		const int iblk1 = blkel_col[iblkel];
		for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs_Blk );
			const unsigned int jblk1 = rowptr[jpsup];
			m_marge_tmp_buffer[jblk1] = jpsup;
		}
		for(unsigned int jblkel=0;jblkel<nblkel_row;jblkel++){
			if( iblkel == jblkel ){	// Marge Diagonal
				const double* pval_in = &emat[(iblkel*nblkel_row+iblkel)*BlkSize];
				double* pval_out = &matval_dia[iblk1*BlkSize];
				for(unsigned int idof=0;idof<BlkSize;idof++){ pval_out[idof] += pval_in[idof]; }
			}
			else{	// Marge Non-Diagonal
				const unsigned int jblk1 = blkel_row[jblkel];
				assert( jblk1 < this->NBlkMatRow() );
				if( m_marge_tmp_buffer[jblk1] == -1 ) continue;
                assert( m_marge_tmp_buffer[jblk1] >= 0 && m_marge_tmp_buffer[jblk1] < (int)m_ncrs_Blk );
				const unsigned int jpsup1 = m_marge_tmp_buffer[jblk1];
				assert( m_rowPtr_Blk[jpsup1] == jblk1 );
				const double* pval_in = &emat[(iblkel*nblkel_row+jblkel)*BlkSize];
				double* pval_out = &matval_nd[jpsup1*BlkSize];
				for(unsigned int idof=0;idof<BlkSize;idof++){ pval_out[idof] += pval_in[idof]; }
			}
		}
		for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs_Blk );
			const unsigned int jblk1 = rowptr[jpsup];
			m_marge_tmp_buffer[jblk1] = -1;
		}
	}
	return true;
}

// bc_flagが０でない自由度に固定境界条件をセットする
bool CMatDia_BlkCrs::SetBoundaryCondition(const MatVec::CBCFlag& bc_flag)
{
	assert( LenBlkCol() == LenBlkRow() );
	assert( NBlkMatCol() == NBlkMatRow() );
	assert( bc_flag.NBlk() == NBlkMatRow() );
	assert( bc_flag.LenBlk() == LenBlkRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
	    // bc_flagが０でない行列の行を０にする。
	    // 但し対角成分は１にしておく
	    for(unsigned int iblk=0;iblk<NBlkMatRow();iblk++){
            const unsigned int len_i = this->LenBlkCol(iblk);
            assert( bc_flag.LenBlk(iblk) == len_i );
		    for(unsigned int idof=0;idof<len_i;idof++){
			    if( bc_flag.GetBCFlag(iblk,idof) == 0 ) continue;
                {   // 対角ブロックに値を設定
                    assert( m_DiaValPtr != 0 );
                    const unsigned int ipos0 = this->m_DiaValPtr[iblk];
			        for(unsigned int kdof=0;kdof<len_i;kdof++){ m_valDia_Blk[ipos0+kdof*len_i+idof] = 0.0; }
			        for(unsigned int jdof=0;jdof<len_i;jdof++){ m_valDia_Blk[ipos0+idof*len_i+jdof] = 0.0; }
			        m_valDia_Blk[ipos0+idof*len_i+idof] = 1.0;
                }
			    for(unsigned int ipsup=m_colInd_Blk[iblk];ipsup<m_colInd_Blk[iblk+1];ipsup++){
                    assert( ipsup < NCrs() );
                    const unsigned int jblk0 = m_rowPtr_Blk[ipsup];
                    const unsigned int len_j = this->LenBlkRow(jblk0);
                    const unsigned int ipos0 = this->m_ValPtr[ipsup];
			        for(unsigned int jdof=0;jdof<len_j;jdof++){
				        m_valCrs_Blk[ipos0+idof*len_j+jdof] = 0.0; 
			        }
			    }
		    }
	    }
	    // bc_flagが０でない行列の列を０にする。
	    for(unsigned int iblk=0;iblk<NBlkMatRow();iblk++){
            const unsigned int len_i = this->LenBlkCol(iblk);
            assert( bc_flag.LenBlk(iblk) == len_i );
	        for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
		        const unsigned int jblk0 = m_rowPtr_Blk[icrs];
                const unsigned int len_j = this->LenBlkRow(jblk0);
                assert( bc_flag.LenBlk(jblk0) == len_j );
                assert( m_ValPtr != 0 );
                const unsigned int ipos0 = this->m_ValPtr[icrs];
		        for(unsigned int jdof=0;jdof<len_j;jdof++){
			        if( bc_flag.GetBCFlag(jblk0,jdof) == 0 ) continue;
			        for(unsigned int idof=0;idof<len_i;idof++){
				        m_valCrs_Blk[ipos0+idof*len_j+jdof] = 0.0;
			        }
		        }
            }
	    }
        return true;
    }
	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int BlkLen = LenBlkCol();
	
	// bc_flagが０でない行列の行を０にする。
	// 但し対角成分は１にしておく
	for(unsigned int iblk=0;iblk<NBlkMatRow();iblk++){
		for(unsigned int idof=0;idof<BlkLen;idof++){
			if( bc_flag.GetBCFlag(iblk,idof) == 0 ) continue;
			for(unsigned int kdof=0;kdof<BlkLen;kdof++){ m_valDia_Blk[iblk*BlkSize+kdof*BlkLen+idof] = 0.0; }
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){ m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] = 0.0; }
			m_valDia_Blk[iblk*BlkSize+idof*BlkLen+idof] = 1.0;
			for(unsigned int ipsup=m_colInd_Blk[iblk];ipsup<m_colInd_Blk[iblk+1];ipsup++){
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){ 
				m_valCrs_Blk[ipsup*BlkSize+idof*BlkLen+jdof] = 0.0; 
			}
			}
		}
	}
	// bc_flagが０でない行列の列を０にする。
	for(unsigned int ipsup=0;ipsup<this->NCrs();ipsup++){
		const int jblk1 = m_rowPtr_Blk[ipsup];
		for(unsigned int jdof=0;jdof<BlkLen;jdof++){
			if( bc_flag.GetBCFlag(jblk1,jdof) == 0 ) continue;
			for(unsigned int idof=0;idof<BlkLen;idof++){
				m_valCrs_Blk[ipsup*BlkSize+idof*BlkLen+jdof] = 0.0;
			}
		}
	}
	return true;
}

// 非ゼロパターンを加える
bool CMatDia_BlkCrs::AddPattern(const Com::CIndexedArray& crs)
{
	// 入力チェック
	assert( crs.CheckValid() );
	if( !crs.CheckValid() ) return false;
	if( crs.Size() > NBlkMatCol() ) return false;
	{   // サイズが超えてないか
		unsigned int max_val = 0;
		for(unsigned int i=0;i<crs.array.size();i++){
			max_val = (max_val>crs.array[i])?max_val:crs.array[i];
		}
		if( max_val>NBlkMatRow() ) return false;
	}
  {   // 対角に要素があっちゃいけないことにする．
    const unsigned int nblk = crs.Size();
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      for(unsigned int icrs=crs.index[iblk];icrs<crs.index[iblk+1];icrs++){
        const unsigned int jblk0 = crs.array[icrs];
        assert( iblk != jblk0 );
      }
    }
  }
  

	// 入力行列が０だったら何もしなくてよい
	if( crs.array.size() == 0 ) return true;

  assert( m_valDia_Blk != 0); // 行列のサイズが与えられているとする

  ////////////////

	if( m_ncrs_Blk == 0 ){	// 何も要素が入っていない行列にパターンを追加する
		CMatDia_BlkCrs::DeletePattern();
		assert( m_colInd_Blk != 0 );
		const unsigned int nblk = NBlkMatCol();
		{
			const unsigned int nblk_crs = crs.Size();
			for(unsigned int iblk=0;iblk<nblk_crs+1;iblk++){
				m_colInd_Blk[iblk] = crs.index[iblk];
			}
			for(unsigned int iblk=nblk_crs+1;iblk<nblk+1;iblk++){
				m_colInd_Blk[iblk] = crs.index[nblk_crs];
			}
		}
		const unsigned int ncrs = m_colInd_Blk[nblk];
		assert( m_rowPtr_Blk == 0 );
		m_rowPtr_Blk = new unsigned int [ncrs];
		for(unsigned int icrs=0;icrs<ncrs;icrs++){
			m_rowPtr_Blk[icrs] = crs.array[icrs];
		}
		m_ncrs_Blk = ncrs;
		return true;
	}

    ////////////////////////////////////////////////////////////////
  
	bool is_included = true; 
	{	// この要素が行列に含まれているかどうか調べる
		int* tmp_buffer = new int [this->NBlkMatRow()];
		for(unsigned int iblk=0;iblk<this->NBlkMatRow();iblk++){
			tmp_buffer[iblk] = 0;
		}
    //		const unsigned int ncrs_size = crs.size;
    for(unsigned int iblk=0;iblk<crs.Size();iblk++){
			// 既に要素が入っている場所にフラグを立てる
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				unsigned int jnode = m_rowPtr_Blk[icrs];
				tmp_buffer[jnode] = 1;
			}
			// パターンが適合しているかどうか調べる
			for(unsigned int icrs=crs.index[iblk];icrs<crs.index[iblk+1];icrs++){
				unsigned int jnode = crs.array[icrs];
				if( tmp_buffer[jnode] != 1 ){
					is_included = false;
					break;
				}
			}
			if( !is_included ) break;
			// フラグを元に戻す
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				unsigned int jnode = m_rowPtr_Blk[icrs];
				tmp_buffer[jnode] = 0;
			}
			if( !is_included ) break;
		}
		delete[] tmp_buffer;
	}
	if( is_included ) return true;

	// パターンを追加する
	std::vector<unsigned int> tmp_row_ptr;
	tmp_row_ptr.reserve(crs.Size()+this->NCrs());
	{	// とりあえずこのtmp_row_ptrを作ってから、コピーする
		int* tmp_buffer = new int [this->NBlkMatRow()];
		for(unsigned int iblk=0;iblk<NBlkMatRow();iblk++){ tmp_buffer[iblk] = 0; }
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				unsigned int jblk = m_rowPtr_Blk[icrs];
				assert( jblk < this->NBlkMatRow() );
				if( tmp_buffer[jblk] != (int)iblk+1 ){
					tmp_row_ptr.push_back(jblk);
					tmp_buffer[jblk] = iblk+1;
					assert( iblk != jblk );
				}
			}
			if( iblk < crs.Size() ){	// ここは忘れがちで注意が必要
				for(unsigned int icrs=crs.index[iblk];icrs<crs.index[iblk+1];icrs++){
					assert( icrs < crs.array.size() );
					unsigned int jblk = crs.array[icrs];
					assert( jblk < this->NBlkMatRow() );
          if( tmp_buffer[jblk] != (int)iblk+1 ){
						tmp_row_ptr.push_back(jblk);
						tmp_buffer[jblk] = iblk+1;
					}
				}
			}
			m_colInd_Blk[iblk] = tmp_row_ptr.size();
		}
		delete[] tmp_buffer;
	}
	for(unsigned int iblk=NBlkMatCol();iblk>0;iblk--){
		m_colInd_Blk[iblk] = m_colInd_Blk[iblk-1];
	}
	m_colInd_Blk[0] = 0;
	delete[] m_rowPtr_Blk;
	m_ncrs_Blk = m_colInd_Blk[NBlkMatCol()];
	m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
	for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
		m_rowPtr_Blk[icrs] = tmp_row_ptr[icrs];
	}
	// row_ptrを昇順にバブルソートする
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
		const unsigned int is_crs = m_colInd_Blk[iblk  ];
		const unsigned int ie_crs = m_colInd_Blk[iblk+1];
		if( is_crs == ie_crs ) continue;
		assert( is_crs < ie_crs );
		int itmp;
		for(unsigned int icrs=is_crs;icrs<ie_crs-1;icrs++){
			for(unsigned int jcrs=ie_crs-1;jcrs>icrs;jcrs--){
				if( m_rowPtr_Blk[jcrs] < m_rowPtr_Blk[jcrs-1] ){
					itmp = m_rowPtr_Blk[jcrs];
					m_rowPtr_Blk[jcrs] = m_rowPtr_Blk[jcrs-1];
					m_rowPtr_Blk[jcrs-1] = itmp;
				}
			}
		}
	}
	/*
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
		std::cout << iblk << "-->";
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			unsigned int jblk = m_rowPtr_Blk[icrs];
			std::cout << jblk << " ";
		}
		std::cout << std::endl;
	}
	*/

	return true;
}

// 非ゼロパターンを加える
bool CMatDia_BlkCrs::AddPattern(const CMatDia_BlkCrs& rhs, const bool isnt_trans){
	if( isnt_trans ){
		assert( NBlkMatCol() == rhs.NBlkMatRow() );
		assert( NBlkMatRow() == rhs.NBlkMatCol() );
		if( m_ncrs_Blk == 0 ){
			DeletePattern();
			CMat_BlkCrs::AddPattern(rhs,isnt_trans);
		}
		else{
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
			abort();
		}
	}
	else{
		std::cout << "Error!-->Not Implimented" <<  std::endl;
		assert(0);
		abort();
	}
	return true;
}

bool CMatDia_BlkCrs::AddPattern(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order)
{
	assert( rhs.NBlkMatCol() == rhs.NBlkMatRow() );
	assert( rhs.NBlkMatCol() == order.NBlk() );
	if( this->NBlkMatCol() == 0 ){
		assert( this->NBlkMatRow()==0 );
		assert( this->LenBlkCol()==0 );
		assert( this->LenBlkRow()==0 );
		const unsigned int nblk = order.NBlk();
		const unsigned int len = rhs.LenBlkCol();
        CMat_BlkCrs::Initialize(nblk,len,nblk,len);
		assert( this->m_valDia_Blk == 0 );
		m_valDia_Blk = new double [nblk*len];
	}
	assert( NBlkMatCol() == order.NBlk() );
	assert( NBlkMatRow() == NBlkMatCol() );
	if( m_ncrs_Blk == 0 ){
		DeletePattern();
		CMat_BlkCrs::AddPattern(rhs,order,order);
	}
	else{
		std::cout << "Error!-->Not Implimented" << std::endl;
		assert(0);
		abort();
	}
	return true;
}

// M1*M2*M3のパターンを加える（マルチグリッド用)
bool CMatDia_BlkCrs::AddPattern(const CMat_BlkCrs& m1, const CMatDia_BlkCrs& m2, const CMat_BlkCrs& m3)
{
	assert( NBlkMatCol()    == m1.NBlkMatCol() );
	assert( m1.NBlkMatRow() == m2.NBlkMatCol() );
	assert( m2.NBlkMatRow() == m3.NBlkMatCol() );
	assert( m3.NBlkMatRow() == NBlkMatRow() );

	if( m_ncrs_Blk == 0 ){

		DeletePattern();

		const unsigned int mat_len_mid = m2.NBlkMatCol();

		m_colInd_Blk[0] = 0;

		std::vector<unsigned int> row_ptr;
		row_ptr.reserve( m2.m_ncrs_Blk*4 );

		////////////////
		unsigned int* is_l_flag = new unsigned int [mat_len_mid];
		for(unsigned int lblk=0;lblk<mat_len_mid;lblk++){ is_l_flag[lblk] = 0; }
		unsigned int l_size = 0;
		unsigned int* l_buffer = new unsigned int [mat_len_mid];
		////////////////
		unsigned int* is_j_flag = new unsigned int [NBlkMatRow()];
		for(unsigned int jblk=0;jblk<NBlkMatRow();jblk++){ is_j_flag[jblk] = 0; }
		unsigned int j_size = 0;
		unsigned int* j_buffer = new unsigned int [NBlkMatRow()];
		////////////////
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			{	// Make l_buffer, l_size
				l_size = 0;
				unsigned int npsuk;
				const unsigned int* psuk = m1.GetPtrIndPSuP(iblk,npsuk);
				assert( psuk != 0 );
				for(unsigned int k=0;k<npsuk;k++){
					const unsigned int kblk0 = psuk[k];
					assert( kblk0 < mat_len_mid );
					for(unsigned int klcrs=m2.m_colInd_Blk[kblk0];klcrs<m2.m_colInd_Blk[kblk0+1];klcrs++){
						assert( klcrs < m2.m_ncrs_Blk );
						const unsigned int lblk0 = m2.m_rowPtr_Blk[klcrs];
						assert( lblk0 < mat_len_mid );
						if( is_l_flag[lblk0] == 0 ){
							l_buffer[l_size] = lblk0;
							is_l_flag[lblk0] = l_size+1;
							l_size++;
						}
					}
					// Set Diagonal
					if( is_l_flag[kblk0] == 0 ){
						l_buffer[l_size] = kblk0;
						is_l_flag[kblk0] = l_size+1;
						l_size++;
					}
				}
				for(unsigned int l=0;l<l_size;l++){
					unsigned int lblk0 = l_buffer[l];
					is_l_flag[lblk0] = 0;
				}
			}
			{	// Make j_buffer, j_size
				j_size = 0;
				for(unsigned int l=0;l<l_size;l++){
					unsigned int lblk0 = l_buffer[l];
					unsigned int npsul;
					const unsigned int* psul = m3.GetPtrIndPSuP(lblk0,npsul);
					assert( psul != 0 );
					for(unsigned int j=0;j<npsul;j++){
						const unsigned int jblk0 = psul[j];
						assert( jblk0 < NBlkMatRow() );
						if( jblk0 == iblk ) continue;
						if( is_j_flag[jblk0] == 0 ){
							j_buffer[j_size] = jblk0;
							is_j_flag[jblk0] = j_size+1;
							j_size++;
						}
					}
				}
				for(unsigned int j=0;j<j_size;j++){
					unsigned int jblk0 = j_buffer[j];
					is_j_flag[jblk0] = 0;
				}
			}
			{
				std::sort(&j_buffer[0],&j_buffer[j_size]);
				m_colInd_Blk[iblk+1] = m_colInd_Blk[iblk] + j_size;
				unsigned int ijcrs0 = m_colInd_Blk[iblk];
				row_ptr.resize( row_ptr.size() + j_size );
				for(unsigned int j=0;j<j_size;j++){
					row_ptr[ijcrs0] = j_buffer[j];
					ijcrs0++;
				}
			}
		}
		delete[] l_buffer;	delete[] is_l_flag;	
		delete[] j_buffer;	delete[] is_j_flag; 

		m_ncrs_Blk = row_ptr.size();
		m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
		for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
			m_rowPtr_Blk[icrs] = row_ptr[icrs];
		}
	}
	else{
		std::cout << "Error!-->Not Implimented " << std::endl;
		assert(0);
		abort();
	}
 
//	std::cout << "Mat Size : " << m_ncrs_Blk << std::endl;

	return true;
}

bool CMatDia_BlkCrs::SetValue(const CMatDia_BlkCrs& rhs, const bool isnt_trans)
{
	assert( NBlkMatRow() == rhs.NBlkMatRow() );
	assert( NBlkMatCol() == rhs.NBlkMatCol() );
	assert( NBlkMatCol() == NBlkMatRow() );
	
	assert( LenBlkCol() == rhs.LenBlkCol() );
	assert( LenBlkRow() == rhs.LenBlkRow() );
	assert( LenBlkCol() == LenBlkRow() );

    assert( m_valDia_Blk != 0 );
    assert( rhs.m_valDia_Blk != 0 );

	CMat_BlkCrs::SetValue(rhs,isnt_trans);
    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        if( LenBlkCol() >= 0 || LenBlkRow() >= 0 || !isnt_trans ){
            std::cout << "Error!-->Not Implimented!" << std::endl;
            assert(0);
            return false;
        }
        assert( rhs.m_DiaValPtr  != 0 );
        assert( m_DiaValPtr  != 0 );
        const unsigned int nblk = this->NBlkMatCol();
		const unsigned int ndof= this->m_DiaValPtr[nblk];
        assert( rhs.m_DiaValPtr[nblk] == ndof );
		for(unsigned int idof=0;idof<ndof;idof++){
			m_valDia_Blk[idof] = rhs.m_valDia_Blk[idof];
		}
        return true;
    }
    ////////////////
	if( isnt_trans ){
		const unsigned int ndof=NBlkMatCol()*LenBlkCol()*LenBlkRow();
		for(unsigned int idof=0;idof<ndof;idof++){
			m_valDia_Blk[idof] = rhs.m_valDia_Blk[idof];
		}
	}
	else{
		const unsigned int BlkLen = LenBlkCol();
		const unsigned int BlkSize = BlkLen*BlkLen;
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
			for(unsigned int idof=0;idof<BlkLen;idof++){
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){
				m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] 
					= rhs.m_valDia_Blk[iblk*BlkSize+jdof*BlkLen*idof];
			}
			}
		}
	}
	return true;
}

bool CMatDia_BlkCrs::SetValue(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order)
{
	assert( rhs.NBlkMatCol() == order.NBlk() );
	assert( rhs.NBlkMatRow() == order.NBlk() );
	assert( rhs.NBlkMatCol() == rhs.NBlkMatRow() );

	assert( NBlkMatRow() == order.NBlk() );
	assert( NBlkMatCol() == order.NBlk() );
	assert( NBlkMatCol() == NBlkMatRow() );
	
	assert( LenBlkCol() == rhs.LenBlkCol() );
	assert( LenBlkRow() == rhs.LenBlkRow() );
	assert( LenBlkCol() == LenBlkRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implimented!" << std::endl;
        assert(0);
        return false;
    }
	if( m_valCrs_Blk == 0 ){	
		m_valCrs_Blk = new double [m_ncrs_Blk];
	}
	
	CMat_BlkCrs::SetValue(rhs,order,order);
	const unsigned int BlkLen = LenBlkCol();
	const unsigned int BlkSize = BlkLen*BlkLen;
	const unsigned int nblk = NBlkMatCol();
	for(unsigned int i0blk=0;i0blk<nblk;i0blk++){
		const unsigned int i1blk = (unsigned int)order.NewToOld(i0blk);
		assert( i1blk < rhs.NBlkMatCol() );
		for(unsigned int idof=0;idof<BlkSize;idof++){
			m_valDia_Blk[i0blk*BlkSize+idof] 
				= rhs.m_valDia_Blk[i1blk*BlkSize+idof];
		}
	}
	return true;
}

bool CMatDia_BlkCrs::SetValue(const CMat_BlkCrs& m1, const CMatDia_BlkCrs& m2, const CMat_BlkCrs& m3)
{
	assert( NBlkMatRow() == NBlkMatRow() );
	assert( NBlkMatCol() == m1.NBlkMatCol() );
	assert( m1.NBlkMatRow() == m2.NBlkMatCol() );
	assert( m2.NBlkMatRow() == m3.NBlkMatCol() );
	assert( m3.NBlkMatRow() == NBlkMatRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implimented!" << std::endl;
        assert(0);
        return false;
    }

	if( LenBlkCol() != 1 || m2.LenBlkCol() != 1 || LenBlkRow() != 1 ){
		std::cout << "Error!-->Not Implimented " << std::endl;
		assert(0);
		abort();
	}

	const unsigned int mat_len_mid = m2.NBlkMatCol();

	if( m_valCrs_Blk == 0 ){	
		m_valCrs_Blk = new double [m_ncrs_Blk];
	}
	
	////////////////
	unsigned int* is_l_flag = new unsigned int [mat_len_mid];
	for(unsigned int lblk=0;lblk<mat_len_mid;lblk++){ is_l_flag[lblk] = 0; }
	unsigned int l_size = 0;
	unsigned int* l_buffer = new unsigned int [mat_len_mid];
	double* l_val = new double [mat_len_mid];
	////////////////
	unsigned int* is_j_flag = new unsigned int [NBlkMatRow()];
	for(unsigned int jblk=0;jblk<NBlkMatRow();jblk++){ is_j_flag[jblk] = 0; }
	unsigned int j_size = 0;
	unsigned int* j_buffer = new unsigned int [NBlkMatRow()];
	double* j_val = new double [NBlkMatRow()];
	////////////////
	int* row2crs = new int [NBlkMatRow()];
	for(unsigned int jblk=0;jblk<NBlkMatRow();jblk++){ row2crs[jblk] = -1; }
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
		{	// Make l_buffer, l_size
			l_size = 0;
			unsigned int npsuk = 0;
			const unsigned int* ind_psuk = m1.GetPtrIndPSuP(iblk,npsuk);
			assert( ind_psuk != 0 );
			const double* val_psuk = m1.GetPtrValPSuP(iblk,npsuk);
			assert( val_psuk != 0 );
			for(unsigned int k=0;k<npsuk;k++){
				const unsigned int kblk0 = ind_psuk[k];
				assert( kblk0 < mat_len_mid );
				const double& ik_val = val_psuk[k];
				for(unsigned int klcrs=m2.m_colInd_Blk[kblk0];klcrs<m2.m_colInd_Blk[kblk0+1];klcrs++){
					assert( klcrs < m2.m_ncrs_Blk );
					const unsigned int lblk0 = m2.m_rowPtr_Blk[klcrs];
					assert( lblk0 < mat_len_mid );
					if( is_l_flag[lblk0] == 0 ){
						l_buffer[l_size] = lblk0;
						l_val[l_size] = ik_val*m2.m_valCrs_Blk[klcrs];
						is_l_flag[lblk0] = l_size+1;
						l_size++;
					}
					else{
						const unsigned int l = is_l_flag[lblk0]-1;
						assert( l < l_size );
						assert( l_buffer[l] == lblk0 );
						l_val[l] += ik_val*m2.m_valCrs_Blk[klcrs];
					}
				}
				// Set Diagonal
				if( is_l_flag[kblk0] == 0 ){
					l_buffer[l_size] = kblk0;
					l_val[l_size] = ik_val*m2.m_valDia_Blk[kblk0];
					is_l_flag[kblk0] = l_size+1;
					l_size++;
				}
				else{
					const unsigned int l = is_l_flag[kblk0]-1;
					assert(  l_buffer[l] == kblk0 );
					l_val[l] += ik_val*m2.m_valDia_Blk[kblk0];
				}
			}
			for(unsigned int l=0;l<l_size;l++){
				unsigned int lblk0 = l_buffer[l];
				assert( is_l_flag[lblk0] -1 == l );
				is_l_flag[lblk0] = 0;
			}
		}
		{	// Make j_buffer, j_size
			j_size = 0;
			for(unsigned int l=0;l<l_size;l++){
				unsigned int lblk0 = l_buffer[l];
				assert( lblk0  < mat_len_mid );
				const double& lj_val = l_val[l];
				unsigned int npsul = 0;
				const unsigned int* ind_psul = m3.GetPtrIndPSuP(lblk0,npsul);
				assert( ind_psul != 0 );
				const double* val_psul = m3.GetPtrValPSuP(lblk0,npsul);
				assert( val_psul != 0 );
				for(unsigned int l=0;l<npsul;l++){
					const unsigned int jblk0 = ind_psul[l];
					assert( jblk0 < NBlkMatRow() );
					if( is_j_flag[jblk0] == 0 ){
						j_buffer[j_size] = jblk0;
						j_val[j_size] = lj_val*val_psul[l];
						is_j_flag[jblk0] = j_size+1;
						j_size++;
					}
					else{
						const unsigned int j = is_j_flag[jblk0]-1;
						assert( j < j_size );
						assert( j_buffer[j] == jblk0 );
						j_val[j] += lj_val*val_psul[l];
					}
				}
			}
			for(unsigned int j=0;j<j_size;j++){
				unsigned int jblk0 = j_buffer[j];
				assert( is_j_flag[jblk0]-1 == j );
				is_j_flag[jblk0] = 0;
			}
		}
		////////////////
		for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
			assert( ijcrs<m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
			assert( jblk0 < NBlkMatRow() );
			row2crs[jblk0] = ijcrs;
		}
		for(unsigned int j=0;j<j_size;j++){
			const unsigned int jblk0 = j_buffer[j];
			assert( jblk0<NBlkMatRow() );
			if( jblk0 != iblk ){
				const int ijcrs0 = row2crs[jblk0];
				if( ijcrs0 == -1 ){	continue; }
				m_valCrs_Blk[ijcrs0] = j_val[j];
			}
			else{
				m_valDia_Blk[iblk] = j_val[j];
			}
		}
		for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
			const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
			assert( jblk0 < NBlkMatRow() );
			row2crs[jblk0] = -1;
		}
	}
	delete[] l_buffer;	delete[] is_l_flag;	delete[] l_val;
	delete[] j_buffer;	delete[] is_j_flag; delete[] j_val;
	delete[] row2crs;

	return true;
}

// Calc Matrix Vector Product
// {y} = alpha * [A]{x} + beta * {y}
bool CMatDia_BlkCrs::MatVec(double alpha, const CVector_Blk& x, double beta, CVector_Blk& y) const
{
	assert( NBlkMatCol() == NBlkMatRow() );

	assert( x.NBlk() == NBlkMatRow() );
	assert( x.Len()  == LenBlkRow() );

	assert( y.NBlk() == NBlkMatCol() );
	assert( y.Len()  == LenBlkCol() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double* iyval = y.m_Value+y.m_DofPtr[iblk];
            const unsigned int nleni = y.Len(iblk);
            assert( this->LenBlkCol(iblk) == nleni );
			for(unsigned int idof=0;idof<nleni;idof++){ iyval[idof] *= beta; }
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < nblk );
                if( iblk == jblk0 ) continue;
				const double* jxval = x.m_Value + x.m_DofPtr[jblk0];
                const unsigned int nlenj = x.Len(jblk0);
                assert( this->LenBlkRow(jblk0) == nlenj );
                const unsigned int valptr0 = this->m_ValPtr[icrs];
				for(unsigned int idof=0;idof<nleni;idof++){
				for(unsigned int jdof=0;jdof<nlenj;jdof++){
					iyval[idof] += alpha * m_valCrs_Blk[valptr0+idof*nlenj+jdof] * jxval[jdof];
				}
				}
			}
			const double* jxval = x.m_Value + x.m_DofPtr[iblk];
            assert( nleni == x.Len(iblk) );
            const unsigned int valptr0 = this->m_DiaValPtr[iblk];
			for(unsigned int idof=0;idof<nleni;idof++){
			for(unsigned int jdof=0;jdof<nleni;jdof++){
				iyval[idof] += alpha * m_valDia_Blk[valptr0+idof*nleni+jdof] * jxval[jdof];
			}
			}
		}
        return true;
    }

    ////////////////////////////////

	const unsigned int BlkLen = LenBlkCol();
	const unsigned int BlkSize = BlkLen*BlkLen;

	if( BlkLen == 1 ){
		// コンパイラ高速化のためにスコープ内にメンバ変数を取り出しておく
		const double* matval_nd  = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const double* xval = x.m_Value;
		double* yval = y.m_Value;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double* iyval = yval+iblk;
			(*iyval) *= beta;
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk );
				const double* jxval = xval+jblk0;
				(*iyval) += alpha * matval_nd[icrs] * (*jxval);
			}
			const double* ixval = xval + iblk;
			(*iyval) += alpha * matval_dia[iblk] * (*ixval);
		}
	}
	else if( BlkLen == 2 ){
		// コンパイラ高速化のためにスコープ内にメンバ変数を取り出しておく
		const double* matval_nd  = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const double* xval = x.m_Value;
		double* yval = y.m_Value;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double* iyval = yval+iblk*2;
			iyval[0] *= beta; 
			iyval[1] *= beta;
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk );
				const double* jxval = xval+jblk0*2;
				iyval[0] += alpha * ( matval_nd[icrs*4  ]*jxval[0] + matval_nd[icrs*4+1]*jxval[1] );
				iyval[1] += alpha * ( matval_nd[icrs*4+2]*jxval[0] + matval_nd[icrs*4+3]*jxval[1] );
			}
			const double* ixval = xval + iblk*2;
			iyval[0] += alpha * ( matval_dia[iblk*4  ]*ixval[0] + matval_dia[iblk*4+1]*ixval[1] );
			iyval[1] += alpha * ( matval_dia[iblk*4+2]*ixval[0] + matval_dia[iblk*4+3]*ixval[1] );
		}
	}
	else if( BlkLen == 3 ){
		// コンパイラ高速化のためにスコープ内にメンバ変数を取り出しておく
		const double* matval_nd  = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const double* xval = x.m_Value;
		double* yval = y.m_Value;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double* iyval = yval+iblk*3;
			iyval[0] *= beta; 
			iyval[1] *= beta;
			iyval[2] *= beta;
			const unsigned int icrs0 = colind[iblk];
			const unsigned int icrs1 = colind[iblk+1];
			for(unsigned int icrs=icrs0;icrs<icrs1;icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk );
				const double* jxval = xval+jblk0*3;
				iyval[0] += alpha * ( matval_nd[icrs*9  ]*jxval[0] + matval_nd[icrs*9+1]*jxval[1] + matval_nd[icrs*9+2]*jxval[2] );
				iyval[1] += alpha * ( matval_nd[icrs*9+3]*jxval[0] + matval_nd[icrs*9+4]*jxval[1] + matval_nd[icrs*9+5]*jxval[2] );
				iyval[2] += alpha * ( matval_nd[icrs*9+6]*jxval[0] + matval_nd[icrs*9+7]*jxval[1] + matval_nd[icrs*9+8]*jxval[2] );
			}
			const double* ixval = xval + iblk*3;
			iyval[0] += alpha * ( matval_dia[iblk*9  ]*ixval[0] + matval_dia[iblk*9+1]*ixval[1] + matval_dia[iblk*9+2]*ixval[2] );
			iyval[1] += alpha * ( matval_dia[iblk*9+3]*ixval[0] + matval_dia[iblk*9+4]*ixval[1] + matval_dia[iblk*9+5]*ixval[2] );
			iyval[2] += alpha * ( matval_dia[iblk*9+6]*ixval[0] + matval_dia[iblk*9+7]*ixval[1] + matval_dia[iblk*9+8]*ixval[2] );
		}
	}
	else{
		// コンパイラ高速化のためにスコープ内にメンバ変数を取り出しておく
		const double* matval_nd  = m_valCrs_Blk;
		const double* matval_dia = m_valDia_Blk;
		const unsigned int* colind = m_colInd_Blk;
		const unsigned int* rowptr = m_rowPtr_Blk;
		const double* xval = x.m_Value;
		double* yval = y.m_Value;
		const unsigned int nblk = this->NBlkMatCol();
		////////////////
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double* iyval = yval+iblk*BlkLen;
			for(unsigned int idof=0;idof<BlkLen;idof++){ iyval[idof] *= beta; }
			const unsigned int colind0 = colind[iblk];
			const unsigned int colind1 = colind[iblk+1];
			for(unsigned int icrs=colind0;icrs<colind1;icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = rowptr[icrs];
				assert( jblk0 < nblk );
				const double* jxval = xval+jblk0*BlkLen;
				for(unsigned int idof=0;idof<BlkLen;idof++){
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){
					iyval[idof] += alpha * matval_nd[icrs*BlkSize+idof*BlkLen+jdof] * jxval[jdof];
				}
				}
			}
			const double* ixval = xval + iblk*BlkLen;
			for(unsigned int idof=0;idof<BlkLen;idof++){
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){
				iyval[idof] += alpha * matval_dia[iblk*BlkSize+idof*BlkLen+jdof] * ixval[jdof];
			}
			}
		}
	}
	return true;
}
