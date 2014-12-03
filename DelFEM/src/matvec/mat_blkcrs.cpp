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
// Mat_BlkCrs.h : interface of blk crs matrix class (CMat_BlkCrs)
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#ifndef for
#define for if(0); else for
#endif

#include <cassert>
#include <iostream>
#include <algorithm>

#include "delfem/indexed_array.h"
#include "delfem/matvec/mat_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/ordering_blk.h"
#include "delfem/matvec/bcflag_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// construction and destruction
//////////////////////////////////////////////////////////////////////

CMat_BlkCrs::CMat_BlkCrs()
: m_nblk_MatCol(0), m_nblk_MatRow(0), m_len_BlkCol(0), m_len_BlkRow(0)
{
//	std::cout << "Construct : CMat_BlkCrs (Default Constructor) " << std::endl;
	m_colInd_Blk = 0;
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0;
    
    m_DofPtrCol = 0;
    m_DofPtrRow = 0;
    m_ValPtr = 0;

    m_marge_tmp_buffer = 0;
}

CMat_BlkCrs::CMat_BlkCrs(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
                         unsigned int nblk_row, const std::vector<unsigned int>& alen_row )
: m_nblk_MatCol(nblk_col), m_nblk_MatRow(nblk_row), m_len_BlkCol(-1), m_len_BlkRow(-1)
{
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 

    ////////////////
    m_DofPtrCol = new unsigned int [nblk_col+1]; assert( m_DofPtrCol != 0 );
    m_DofPtrCol[0] = 0;
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        m_DofPtrCol[iblk+1] = m_DofPtrCol[iblk] + alen_col[iblk];
    }
    m_len_BlkCol = -1;

    ////////////////
    m_DofPtrRow = new unsigned int [nblk_row+1]; assert( m_DofPtrRow != 0 );
    m_DofPtrRow[0] = 0;
    for(unsigned int iblk=0;iblk<nblk_row;iblk++){
        m_DofPtrRow[iblk+1] = m_DofPtrRow[iblk] + alen_row[iblk];
    }
    m_len_BlkRow = -1;

    ////////////////
    m_ValPtr = 0;
    m_marge_tmp_buffer = 0;
}


CMat_BlkCrs::CMat_BlkCrs(const unsigned int nblk_col, const unsigned int len_col, 
				         const unsigned int nblk_row, const unsigned int len_row )
: m_nblk_MatCol(nblk_col), m_nblk_MatRow(nblk_row), m_len_BlkCol(len_col), m_len_BlkRow(len_row)
{
//	std::cout << "Construct : CMat_BlkCrs " << std::endl;
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 

    m_DofPtrCol = 0;
    m_DofPtrRow = 0;
    m_ValPtr = 0;
    
    m_marge_tmp_buffer = 0;
}

CMat_BlkCrs::CMat_BlkCrs(const CMat_BlkCrs& rhs, bool is_value, bool isnt_trans)
: m_nblk_MatCol( isnt_trans ? rhs.m_nblk_MatCol : rhs.m_nblk_MatRow ),
  m_nblk_MatRow( isnt_trans ? rhs.m_nblk_MatRow : rhs.m_nblk_MatCol ),
  m_len_BlkCol(  isnt_trans ? rhs.m_len_BlkCol  : rhs.m_len_BlkRow ),
  m_len_BlkRow(  isnt_trans ? rhs.m_len_BlkRow  : rhs.m_len_BlkCol )
{
//	std::cout << "Construct : CMat_BlkCrs" << std::endl;
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0;
    
    m_DofPtrCol = 0;
    m_DofPtrRow = 0;
    m_ValPtr = 0;
    
    m_marge_tmp_buffer = 0;

	AddPattern(rhs,isnt_trans);
	if( is_value ){ this->SetValue(rhs,isnt_trans); }
}

bool CMat_BlkCrs::Initialize(const unsigned int nblk_col, const unsigned int len_col, 
						     const unsigned int nblk_row, const unsigned int len_row ){

	if( m_colInd_Blk != 0 ){ delete[] m_colInd_Blk; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; }
    
    if( m_DofPtrCol != 0 ){ delete[] m_DofPtrCol; m_DofPtrCol = 0; }
    if( m_DofPtrRow != 0 ){ delete[] m_DofPtrRow; m_DofPtrRow = 0; }
    if( m_ValPtr    != 0 ){ delete[] m_ValPtr;    m_ValPtr    = 0; }

	m_nblk_MatCol = nblk_col;
	m_nblk_MatRow = nblk_row;
	m_len_BlkCol = len_col;
	m_len_BlkRow = len_row;

	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 

    m_DofPtrCol = 0;
    m_DofPtrRow = 0;
    m_ValPtr    = 0;
    return true;
}

bool CMat_BlkCrs::Initialize(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
                             unsigned int nblk_row, const std::vector<unsigned int>& alen_row )
{
    if( nblk_col != alen_col.size() ){ assert(0); return false; }
    if( nblk_row != alen_row.size() ){ assert(0); return false; }

	if( m_colInd_Blk != 0 ){ delete[] m_colInd_Blk; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; }
    
    if( m_DofPtrCol != 0 ){ delete[] m_DofPtrCol; m_DofPtrCol = 0; }
    if( m_DofPtrRow != 0 ){ delete[] m_DofPtrRow; m_DofPtrRow = 0; }
    if( m_ValPtr    != 0 ){ delete[] m_ValPtr;    m_ValPtr    = 0; }

	m_nblk_MatCol = nblk_col;
	m_nblk_MatRow = nblk_row;

    ////////////////
    m_DofPtrCol = new unsigned int [nblk_col+1];
    m_DofPtrCol[0] = 0;
    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
        m_DofPtrCol[iblk+1] = m_DofPtrCol[iblk] + alen_col[iblk];
    }
    m_len_BlkCol = -1;

    ////////////////
    m_DofPtrRow = new unsigned int [nblk_row+1];
    m_DofPtrRow[0] = 0;
    for(unsigned int iblk=0;iblk<nblk_row;iblk++){
        m_DofPtrRow[iblk+1] = m_DofPtrRow[iblk] + alen_row[iblk];
//        std::cout << iblk << " " << alen_row[iblk] << std::endl;
    }
    m_len_BlkRow = -1;
    
    ////////////////
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 
    return true;
}

CMat_BlkCrs::~CMat_BlkCrs()
{
	if( m_colInd_Blk != 0 ){ delete[] m_colInd_Blk;	m_colInd_Blk = 0; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk;	m_rowPtr_Blk = 0; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk;	m_valCrs_Blk = 0; }
    
    if( m_DofPtrCol != 0 ){ delete[] m_DofPtrCol; m_DofPtrCol = 0; }
    if( m_DofPtrRow != 0 ){ delete[] m_DofPtrRow; m_DofPtrRow = 0; }
    if( m_ValPtr    != 0 ){ delete[] m_ValPtr;    m_ValPtr    = 0; }
}

// パターンを全て消去　RowPtr,Valはメモリ解放
bool CMat_BlkCrs::DeletePattern(){
	m_ncrs_Blk = 0;
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; m_rowPtr_Blk = 0; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; m_valCrs_Blk = 0; }
	return true;
}

bool CMat_BlkCrs::SetZero()
{
    unsigned int ntotdof = 0;
    if( m_len_BlkCol >= 0 && m_len_BlkRow >= 0 ){
	    ntotdof = m_ncrs_Blk*m_len_BlkCol*m_len_BlkRow;
    }
    else{   // 可変ブロックサイズの場合
        assert( m_len_BlkCol == -1 && m_len_BlkRow == -1 );
        assert( m_DofPtrCol != 0 );
        assert( m_DofPtrRow != 0 );
        if( m_ValPtr == 0 ){    // m_ValPtrが作成されて無い場合は作る．
            assert( m_valCrs_Blk == 0 );    // valueが無いならそれを指す領域も確保できない
            m_ValPtr = new unsigned int [m_ncrs_Blk+1];
            for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
                const unsigned int lencol = this->LenBlkCol(iblk);
                for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
                    const unsigned int jblk0 = m_rowPtr_Blk[icrs]; assert( jblk0 < this->NBlkMatRow() );
                    const unsigned int lenrow = this->LenBlkRow(jblk0);
                    m_ValPtr[icrs] = ntotdof;
                    ntotdof += lencol*lenrow;
                }
            }
            m_ValPtr[m_ncrs_Blk] = ntotdof;
        }
        else{ assert( m_valCrs_Blk ); }
	    ntotdof = m_ValPtr[m_ncrs_Blk];
    }
	if( m_valCrs_Blk == 0 ){	// まだ値の領域をを確保していない場合
	    m_valCrs_Blk = new double [ntotdof];
	}
    for(unsigned int idof=0;idof<ntotdof;idof++){ m_valCrs_Blk[idof] = 0.0;	}
	return true;
}


void CMat_BlkCrs::FillPattern()
{
	const unsigned int nblkcol = NBlkMatCol();
	const unsigned int nblkrow = NBlkMatRow();
	m_colInd_Blk[0] = 0;
	for(unsigned int iblk=0;iblk<nblkcol;iblk++){
		m_colInd_Blk[iblk+1] = (iblk+1)*nblkrow;
	}
	m_ncrs_Blk = m_colInd_Blk[nblkcol];
	assert( m_ncrs_Blk == nblkcol*nblkrow );
	m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
	for(unsigned int iblk=0;iblk<nblkcol;iblk++){ 
		for(unsigned int jblk=0;jblk<nblkrow;jblk++){
			m_rowPtr_Blk[iblk*nblkrow+jblk] = jblk; 
		}
	}
}

bool MatVec::CMat_BlkCrs::AddPattern(const Com::CIndexedArray& crs)
{
	// 入力チェック
	if( !crs.CheckValid() ) return false;
	if( crs.Size() > NBlkMatCol() ) return false;
	{
		unsigned int max_val = 0;
		for(unsigned int i=0;i<crs.array.size();i++){
			max_val = (max_val>crs.array[i])?max_val:crs.array[i];
		}
		if( max_val>NBlkMatRow() ) return false;
	}

	if( crs.array.size() == 0 ) return true;

	if( m_ncrs_Blk == 0 ){
		this->DeletePattern();
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
//		for(unsigned int iblk=0;iblk<nblk+1;iblk++){ m_colInd_Blk[iblk] = crs.index[iblk]; }
		const unsigned int ncrs = m_colInd_Blk[nblk];
		assert( m_rowPtr_Blk == 0 );
		assert( crs.array.size() >= ncrs );
		m_rowPtr_Blk = new unsigned int [ncrs];
		for(unsigned int icrs=0;icrs<ncrs;icrs++){ m_rowPtr_Blk[icrs] = crs.array[icrs]; }
		m_ncrs_Blk = ncrs;
//		for(unsigned int iblk=0;iblk<nblk;iblk++){
//			std::cout << iblk  << "-->";
//			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
//				std::cout << m_rowPtr_Blk[icrs] << " ";
//			}
//			std::cout << std::endl;
//		}
		return true;
	}

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
		for(unsigned int iblk=0;iblk>NBlkMatRow();iblk++){ tmp_buffer[iblk] = 0; }
		for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
//			unsigned int icrs_i = tmp_row_ptr.size();
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				unsigned int jblk = m_rowPtr_Blk[icrs];
                if( tmp_buffer[jblk] != (int)iblk+1 ){
					tmp_row_ptr.push_back(jblk);
					tmp_buffer[jblk] = iblk+1;
				}
			}
			if( iblk < crs.Size() ){	// ここは忘れがちで注意が必要
				for(unsigned int icrs=crs.index[iblk];icrs<crs.index[iblk+1];icrs++){
					unsigned int jblk = crs.array[icrs];
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
	return true;
}

bool CMat_BlkCrs::AddPattern(const CMat_BlkCrs& rhs, const bool isnt_trans)
{
	if( rhs.m_ncrs_Blk == 0 ) return true;
	if( isnt_trans ){	// Add Not Transpose Pattern of rhs
		if( this->m_ncrs_Blk == 0 ){
			if( this->m_rowPtr_Blk != 0 ){ delete[] this->m_rowPtr_Blk; }
			assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
			assert( m_nblk_MatRow == rhs.m_nblk_MatRow );
			for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){
				this->m_colInd_Blk[iblk] = rhs.m_colInd_Blk[iblk];
			}
			m_ncrs_Blk = m_colInd_Blk[m_nblk_MatCol];
			assert( m_ncrs_Blk == rhs.m_ncrs_Blk );
			
			if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; m_rowPtr_Blk = 0; }
			assert( this->m_ncrs_Blk = rhs.m_ncrs_Blk );
			this->m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
			for(unsigned int icrs=0;icrs<m_ncrs_Blk;icrs++){
				this->m_rowPtr_Blk[icrs] = rhs.m_rowPtr_Blk[icrs];
			}
		}
		else{
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
			abort();
			/* include されるか調べてから、includeされないならパターンを付け足す */
		}
	}
	else{	// Add Transpose Pattern of rhs
		if( this->m_ncrs_Blk == 0 ){
			if( this->m_rowPtr_Blk != 0 ){ delete[] this->m_rowPtr_Blk; }
			assert( m_nblk_MatCol == rhs.m_nblk_MatRow );
			assert( m_nblk_MatRow == rhs.m_nblk_MatCol );
			for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){
				m_colInd_Blk[iblk] = 0;
			}
			for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){
				for(unsigned int jicrs=rhs.m_colInd_Blk[jblk];jicrs<rhs.m_colInd_Blk[jblk+1];jicrs++){
					assert( jicrs < rhs.m_ncrs_Blk );
					const unsigned int iblk0 = rhs.m_rowPtr_Blk[jicrs];
					assert( iblk0 < m_nblk_MatCol);
					m_colInd_Blk[iblk0+1]++;
				}
			}
			for(unsigned int iblk=1;iblk<m_nblk_MatCol+1;iblk++){
				m_colInd_Blk[iblk] += m_colInd_Blk[iblk-1];
			}
			m_ncrs_Blk = m_colInd_Blk[m_nblk_MatCol];
			assert( m_ncrs_Blk == rhs.m_ncrs_Blk );
			if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; m_rowPtr_Blk = 0; }
			m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
			for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){
				for(unsigned int jicrs=rhs.m_colInd_Blk[jblk];jicrs<rhs.m_colInd_Blk[jblk+1];jicrs++){
					assert( jicrs < rhs.m_ncrs_Blk );
					const unsigned int iblk0 = rhs.m_rowPtr_Blk[jicrs];
					assert( iblk0 < m_nblk_MatCol );
					const unsigned int ijcrs = m_colInd_Blk[iblk0];
					assert( ijcrs < m_ncrs_Blk );
					m_rowPtr_Blk[ijcrs] = jblk;
					m_colInd_Blk[iblk0]++;
				}
			}
			for(unsigned int iblk=m_nblk_MatCol;iblk>=1;iblk--){
				m_colInd_Blk[iblk] = m_colInd_Blk[iblk-1];
			}
			m_colInd_Blk[0] = 0;
		}
		else{
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
			abort();
			/* include されるか調べてから、includeされないならパターンを付け足す */
		}
	}
	return true;
}

bool CMat_BlkCrs::AddPattern(const CMat_BlkCrs& rhs, 
		const COrdering_Blk& order_col, const COrdering_Blk& order_row)
{
	assert( rhs.NBlkMatCol() == order_col.NBlk() );
	assert( rhs.NBlkMatRow() == order_row.NBlk() );
	assert( this->NBlkMatCol() == order_col.NBlk() );
	assert( this->NBlkMatRow() == order_row.NBlk() );

	if( rhs.m_ncrs_Blk == 0 ) return true;	// ???????????????

	if( this->m_ncrs_Blk == 0 ){	// ?????????????
		if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; m_rowPtr_Blk = 0; }
		for(unsigned int i0blk=0;i0blk<m_nblk_MatCol+1;i0blk++){
			m_colInd_Blk[i0blk] = 0;
		}
		for(unsigned int i0blk=0;i0blk<m_nblk_MatCol;i0blk++){
			const unsigned int i1blk0 = order_col.NewToOld(i0blk);
			for(unsigned int ij1crs=rhs.m_colInd_Blk[i1blk0];ij1crs<rhs.m_colInd_Blk[i1blk0+1];ij1crs++){
				assert( ij1crs < rhs.m_ncrs_Blk );
				const unsigned int j1blk0 = rhs.m_rowPtr_Blk[ij1crs];
				assert( j1blk0 < rhs.NBlkMatRow() );
				m_colInd_Blk[i0blk+1]++;
			}
		}
		for(unsigned int i0blk=1;i0blk<m_nblk_MatCol+1;i0blk++){
			m_colInd_Blk[i0blk] += m_colInd_Blk[i0blk-1];
		}
		m_ncrs_Blk = m_colInd_Blk[m_nblk_MatCol];
		m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
		for(unsigned int i0blk=0;i0blk<m_nblk_MatCol;i0blk++){
			const unsigned int i1blk0 = order_col.NewToOld(i0blk);
			for(unsigned int ij1crs=rhs.m_colInd_Blk[i1blk0];ij1crs<rhs.m_colInd_Blk[i1blk0+1];ij1crs++){
				assert( ij1crs < rhs.m_ncrs_Blk );
				const unsigned int j1blk0 = rhs.m_rowPtr_Blk[ij1crs];
				assert( j1blk0 < rhs.NBlkMatRow() );
				const unsigned int j0blk0 = (unsigned int)order_row.OldToNew(j1blk0);
				assert( j0blk0 < m_nblk_MatRow );
				const unsigned int ij0crs0 = m_colInd_Blk[i0blk];
				assert( ij0crs0 < m_ncrs_Blk );
				m_rowPtr_Blk[ij0crs0] = j0blk0;
				m_colInd_Blk[i0blk]++;
			}
		}
		for(unsigned int iblk=m_nblk_MatCol;iblk>=1;iblk--){
			m_colInd_Blk[iblk] = m_colInd_Blk[iblk-1];
		}
		m_colInd_Blk[0] = 0;
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			const unsigned int ijcrs0 = m_colInd_Blk[iblk  ];
			const unsigned int ijcrs1 = m_colInd_Blk[iblk+1];
			std::sort(&m_rowPtr_Blk[ijcrs0],&m_rowPtr_Blk[ijcrs1]);
		}
	}
	else{
		std::cout << "Error!-->Not Implimented" << std::endl;
		assert(0);
		abort();
		/* include ??????????include??????????????? */
	}
	return true;
}



bool CMat_BlkCrs::SetValue(const CMat_BlkCrs& rhs, const bool isnt_trans)
{
	assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
	assert( m_nblk_MatRow == rhs.m_nblk_MatRow );
	assert( m_len_BlkCol == rhs.m_len_BlkCol );
	assert( m_len_BlkRow == rhs.m_len_BlkRow );

    this->SetZero();
    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        if( !isnt_trans ){
            std::cout << "Error!-->Not Implemented" << std::endl;
            assert(0);
            return false;
        }
        if( LenBlkCol() >= 0 || LenBlkRow() >= 0 ){
            std::cout << "Error!-->Not Implemented" << std::endl;
            assert(0);
            return false;
        }
		assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
		assert( m_nblk_MatRow == rhs.m_nblk_MatRow );
		assert( rhs.LenBlkCol() == -1 );
		assert( rhs.LenBlkRow() == -1 );
        assert( this->m_ValPtr != 0 );
        assert( rhs.m_ValPtr != 0 );
		int* row2crs = new int [m_nblk_MatRow];
		for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ row2crs[jblk] = -1; }	
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
            const unsigned int lencol = this->LenBlkCol(iblk);
            assert( rhs.LenBlkCol(iblk) == lencol );
			// iblk行のパターンがある列からcrs番号へのハッシュを作る
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				row2crs[jblk0] = ijcrs;
			}
			for(unsigned int ijcrs=rhs.m_colInd_Blk[iblk];ijcrs<rhs.m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<rhs.m_ncrs_Blk );
				const unsigned int jblk0 = rhs.m_rowPtr_Blk[ijcrs];
                const unsigned int lenrow = this->LenBlkRow(jblk0);
                assert( rhs.LenBlkRow(jblk0) == lenrow );
				assert( jblk0<m_nblk_MatRow );
				const int ijcrs0 = row2crs[jblk0];
				if( ijcrs0 == -1 ) continue;
                assert( ijcrs0 < (int)this->NCrs() );
				////////////////
				const double* pval_in = &rhs.m_valCrs_Blk[ rhs.m_ValPtr[ijcrs] ];
				double* pval_out = &m_valCrs_Blk[ this->m_ValPtr[ijcrs0] ];
                ////////////////
				for(unsigned int idof=0;idof<lencol*lenrow;idof++){ *(pval_out+idof) = *(pval_in+idof); }
			}
			// ハッシュの初期化
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				row2crs[jblk0] = -1;
			}
		}
		delete[] row2crs;
        return true;
    }
    ////////////////
    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemeted" << std::endl;
        assert(0);
        return false;
    }
    const unsigned int len_col = (unsigned int)this->LenBlkCol();
    const unsigned int len_row = (unsigned int)this->LenBlkRow();
    const unsigned int BlkSize = len_col*len_row;
	if( isnt_trans ){	// Set Value -- Not Transpose Pattern
		int* row2crs = new int [m_nblk_MatRow];
		for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ row2crs[jblk] = -1; }	
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			// iblk行のパターンがある列からcrs番号へのハッシュを作る
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				row2crs[jblk0] = ijcrs;
			}
			for(unsigned int ijcrs=rhs.m_colInd_Blk[iblk];ijcrs<rhs.m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<rhs.m_ncrs_Blk );
				const unsigned int jblk0 = rhs.m_rowPtr_Blk[ijcrs];
				assert( jblk0<m_nblk_MatRow );
				const int ijcrs0 = row2crs[jblk0];
				if( ijcrs0 == -1 ) continue;
				////////////////
				const double* pval_in = &rhs.m_valCrs_Blk[ijcrs*BlkSize];
				double* pval_out = &m_valCrs_Blk[ijcrs0*BlkSize];
				for(unsigned int idof=0;idof<BlkSize;idof++){ *(pval_out+idof) = *(pval_in+idof); }
                ////////////////
			}
			// ハッシュの初期化
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				row2crs[jblk0] = -1;
			}
		}
		delete[] row2crs;
	}
	else{	// Set Value -- Transpose Pattern
		assert( m_nblk_MatCol == rhs.m_nblk_MatRow );
		assert( m_nblk_MatRow == rhs.m_nblk_MatCol );
		assert( m_len_BlkCol == rhs.m_len_BlkRow );
        assert( m_len_BlkRow == rhs.m_len_BlkCol );
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				const unsigned int jicrs_begin = rhs.m_colInd_Blk[jblk0];
				const unsigned int jicrs_end = rhs.m_colInd_Blk[jblk0+1];
				if( jicrs_begin == jicrs_end ) continue;
				if( iblk<rhs.m_rowPtr_Blk[jicrs_begin] ) continue;
				if( rhs.m_rowPtr_Blk[ jicrs_end-1 ] < iblk  ) continue;
				for(unsigned int jicrs=rhs.m_colInd_Blk[jblk0];jicrs<rhs.m_colInd_Blk[jblk0+1];jicrs++){
					assert( jicrs<rhs.m_ncrs_Blk );
					const unsigned int iblk0 = rhs.m_rowPtr_Blk[jicrs];
					if( iblk != iblk0 ) continue;
					////////////////
					const double* pval_in = &rhs.m_valCrs_Blk[jicrs*BlkSize];
					double* pval_out = &m_valCrs_Blk[ijcrs*BlkSize];
                    for(unsigned int idof=0;idof<len_col;idof++){
                    for(unsigned int jdof=0;jdof<len_row;jdof++){
                        *(pval_out+idof*m_len_BlkRow+jdof) = *(pval_in+jdof*m_len_BlkCol+idof);
                    }
					}
				}
			}
		}
	}
	return true;
}


bool CMat_BlkCrs::SetValue(const CMat_BlkCrs& rhs, 
		const COrdering_Blk& order_col, const COrdering_Blk& order_row)
{
	assert( rhs.NBlkMatCol() == order_col.NBlk() );
	assert( rhs.NBlkMatRow() == order_row.NBlk() );
	assert( this->NBlkMatCol() == order_col.NBlk() );
	assert( this->NBlkMatRow() == order_row.NBlk() );
	assert( this->LenBlkCol() == rhs.LenBlkCol() );
	assert( this->LenBlkRow() == rhs.LenBlkRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        return false;
    }

	const unsigned int BlkSize = m_len_BlkCol*m_len_BlkRow;
	this->SetZero();

	int* row2crs = new int [m_nblk_MatRow];
	for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ row2crs[jblk] = -1; }
	for(unsigned int i0blk=0;i0blk<m_nblk_MatCol;i0blk++){
		for(unsigned int ij0crs=m_colInd_Blk[i0blk];ij0crs<m_colInd_Blk[i0blk+1];ij0crs++){
			assert( ij0crs<m_ncrs_Blk );
			const unsigned int j0blk0 = m_rowPtr_Blk[ij0crs];
			assert( j0blk0 < m_nblk_MatRow );
			row2crs[j0blk0] = ij0crs;
		}
		const unsigned int i1blk0 = (unsigned int)order_col.NewToOld(i0blk);
		for(unsigned int ij1crs=rhs.m_colInd_Blk[i1blk0];ij1crs<rhs.m_colInd_Blk[i1blk0+1];ij1crs++){
			assert( ij1crs<rhs.m_ncrs_Blk );
			const unsigned int j1blk0 = rhs.m_rowPtr_Blk[ij1crs];
			assert( j1blk0<rhs.NBlkMatRow() );
			////////////////
			const unsigned int j0blk0 = order_row.OldToNew(j1blk0);
			assert( j0blk0<m_nblk_MatRow );
			const int ij0crs0 = row2crs[j0blk0];
			if( ij0crs0 == -1 ) continue;
			////////////////
			const double* pval_in = &rhs.m_valCrs_Blk[ij1crs*BlkSize];
			double* pval_out = &m_valCrs_Blk[ij0crs0*BlkSize];
			for(unsigned int idof=0;idof<BlkSize;idof++){
				*(pval_out+idof) = *(pval_in+idof);
			}
		}
		for(unsigned int ij0crs=m_colInd_Blk[i0blk];ij0crs<m_colInd_Blk[i0blk+1];ij0crs++){
			assert( ij0crs<m_ncrs_Blk );
			const unsigned int j0blk0 = m_rowPtr_Blk[ij0crs];
			assert( j0blk0 < m_nblk_MatRow );
			row2crs[j0blk0] = -1;
		}
	}
	delete[] row2crs;
	return true;
}


bool CMat_BlkCrs::Mearge(
	const unsigned int nblkel_col, const unsigned int* blkel_col,
	const unsigned int nblkel_row, const unsigned int* blkel_row,
	const unsigned int blksize, const double* emat)
{
	assert( m_valCrs_Blk != 0 );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
        assert( nblkel_col == 1 && nblkel_row == 1 );
        const unsigned int iblk0 = blkel_col[0];
        const unsigned int jblk0 = blkel_row[0];
        unsigned int icrs = NCrs();
		for(icrs=m_colInd_Blk[iblk0];icrs<m_colInd_Blk[iblk0+1];icrs++){
			assert( icrs < m_ncrs_Blk );
			const unsigned int jblk1 = m_rowPtr_Blk[icrs];
            if( jblk1 == jblk0 ) break;
        }
        assert( icrs != NCrs() );
        const unsigned int ipos = m_ValPtr[icrs];
        const unsigned int leni = this->LenBlkCol(iblk0);
        const unsigned int lenj = this->LenBlkRow(jblk0);
        assert( blksize == leni*lenj );
        for(unsigned int i=0;i<leni*lenj;i++){
            m_valCrs_Blk[ipos+i] = emat[i];
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
	assert( blksize == BlkSize );

	const unsigned int* colind = m_colInd_Blk;
	const unsigned int* rowptr = m_rowPtr_Blk;
	double* matval_nd = m_valCrs_Blk;

	for(unsigned int iblkel=0;iblkel<nblkel_col;iblkel++){
		const unsigned int iblk1 = blkel_col[iblkel];
		assert( iblk1 < NBlkMatCol() );
		for(unsigned int jpsup=colind[iblk1];jpsup<colind[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs_Blk );
			const unsigned int jblk1 = rowptr[jpsup];
			m_marge_tmp_buffer[jblk1] = jpsup;
		}
		for(unsigned int jblkel=0;jblkel<nblkel_row;jblkel++){
			const unsigned int jblk1 = blkel_row[jblkel];
			assert( jblk1 < NBlkMatRow() );
			if( m_marge_tmp_buffer[jblk1] == -1 ) continue;
            assert( m_marge_tmp_buffer[jblk1] >= 0 && m_marge_tmp_buffer[jblk1] < (int)m_ncrs_Blk );
			const unsigned int jpsup1 = m_marge_tmp_buffer[jblk1];
			assert( jpsup1 < m_ncrs_Blk );
			assert( rowptr[jpsup1] == jblk1 );
			const double* pval_in = &emat[(iblkel*nblkel_row+jblkel)*BlkSize];
			double* pval_out = &matval_nd[jpsup1*BlkSize];
			for(unsigned int idof=0;idof<BlkSize;idof++){ pval_out[idof] += pval_in[idof]; }
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
// bc_flagが０でない行列の列を０にする。
bool CMat_BlkCrs::SetBoundaryCondition_Colum(const CBCFlag& bc_flag)
{    
	assert( bc_flag.NBlk() == NBlkMatRow() );

    const unsigned int nblk_col = this->NBlkMatCol();
//    const unsigned int nblk_row = this->NBlkMatRow();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
        for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            const unsigned int len_col = this->LenBlkCol(iblk);
	        for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
                assert( icrs < this->NCrs() );
		        const unsigned int jblk1 = m_rowPtr_Blk[icrs];
		        assert( jblk1 < NBlkMatRow() );
                const unsigned int len_row = this->LenBlkRow(jblk1);
                assert( len_row == bc_flag.LenBlk(jblk1) );
                const unsigned int ipos0 = this->m_ValPtr[icrs];
		        for(unsigned int jdof=0;jdof<len_row;jdof++){
			        if( bc_flag.GetBCFlag(jblk1,jdof) == 0 ) continue;
			        for(unsigned int idof=0;idof<len_col;idof++){
				        m_valCrs_Blk[ipos0+idof*len_row+jdof] = 0.0;
			        }
		        }
            }
	    }
        return true;
    }
    ////////////////////////////////
	assert( bc_flag.LenBlk() == LenBlkRow() );
	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int lenBlkCol = LenBlkCol();
	const unsigned int lenBlkRow = LenBlkRow();
	for(unsigned int icrs=0;icrs<this->NCrs();icrs++){
		const unsigned int jblk1 = m_rowPtr_Blk[icrs];
		assert( jblk1 < NBlkMatRow() );
		for(unsigned int jdof=0;jdof<lenBlkRow;jdof++){
			if( bc_flag.GetBCFlag(jblk1,jdof) == 0 ) continue;
			for(unsigned int idof=0;idof<lenBlkCol;idof++){
				m_valCrs_Blk[icrs*BlkSize+idof*lenBlkRow+jdof] = 0.0;
			}
		}
	}
	return true;
}

// bc_flagが０でない自由度に固定境界条件をセットする
bool CMat_BlkCrs::SetBoundaryConditionInverse_Colum(const CBCFlag& bc_flag)
{
	assert( bc_flag.NBlk() == NBlkMatRow() );
	assert( bc_flag.LenBlk() == LenBlkRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        return false;
    }

	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int lenBlkCol = LenBlkCol();
	const unsigned int lenBlkRow = LenBlkRow();
	
	// bc_flagが０でない行列の列を０にする。
	const unsigned int ncrs = this->NCrs();
	for(unsigned int icrs=0;icrs<ncrs;icrs++){
		const unsigned int jblk1 = m_rowPtr_Blk[icrs];
		assert( jblk1 < NBlkMatRow() );
		for(unsigned int jdof=0;jdof<lenBlkRow;jdof++){
			if( bc_flag.GetBCFlag(jblk1,jdof) != 0 ) continue;
			for(unsigned int idof=0;idof<lenBlkCol;idof++){
				m_valCrs_Blk[icrs*BlkSize+idof*lenBlkRow+jdof] = 0.0;
			}
		}
	}
	return true;
}



// 行列の行に固定境界条件を設定（非対角サブ行列用）
// bc_flagが０でない自由度の行を０にする
bool CMat_BlkCrs::SetBoundaryCondition_Row(const CBCFlag& bc_flag)
{
	assert( bc_flag.NBlk() == NBlkMatCol() );

	// bc_flagが０でない行列の行を０にする。
	const unsigned int nblk_col = this->NBlkMatCol();
//	const unsigned int nblk_row = this->NBlkMatRow();

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
	    for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            const unsigned int len_i = this->LenBlkCol(iblk); assert( bc_flag.LenBlk(iblk) == len_i );
	        for(unsigned int idof=0;idof<len_i;idof++){
		        if( bc_flag.GetBCFlag(iblk,idof) == 0 ) continue;
		        for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
                    const unsigned int jblk0 = m_rowPtr_Blk[icrs]; assert( jblk0 < this->NBlkMatRow() );
                    const unsigned int len_j = this->LenBlkRow(jblk0);
                    const unsigned int ipos0 = this->m_ValPtr[icrs];
			        for(unsigned int jdof=0;jdof<len_j;jdof++){
				        m_valCrs_Blk[ipos0+idof*len_j+jdof] = 0.0; 
			        }
		        }
	        }
	    }
        return true;
    }

	assert( bc_flag.LenBlk() == LenBlkCol() );
	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int lenBlkCol = this->LenBlkCol();
	const unsigned int lenBlkRow = this->LenBlkRow();
	
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
	for(unsigned int idof=0;idof<lenBlkCol;idof++){
		if( bc_flag.GetBCFlag(iblk,idof) == 0 ) continue;
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			for(unsigned int jdof=0;jdof<lenBlkRow;jdof++){
				m_valCrs_Blk[icrs*BlkSize+idof*lenBlkRow+jdof] = 0.0; 
			}
		}
	}
	}
	return true;
}


bool CMat_BlkCrs::MatVec(double alpha, const CVector_Blk& x, double beta, CVector_Blk& b, const bool isnt_trans) const
{
	assert( x.NBlk() == m_nblk_MatRow );
	assert( b.NBlk() == m_nblk_MatCol );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
		const unsigned int nblk_row = this->NBlkMatRow();
		const unsigned int nblk_col = this->NBlkMatCol();
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
            const unsigned int len_col = this->LenBlkCol(iblk); assert( b.Len(iblk) == len_col );
			double* iyval = b.GetValuePtr(iblk);
			for(unsigned int idof=0;idof<len_col;idof++){ iyval[idof] *= beta; }
            for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){ assert( icrs < NCrs() );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs]; assert( jblk0 < nblk_row );
                const unsigned int len_row = this->LenBlkRow(jblk0); assert( x.Len(jblk0) == len_row );
                const unsigned int ipos0 = m_ValPtr[icrs];
				const double* jxval = x.GetValuePtr(jblk0);
				for(unsigned int idof=0;idof<len_col;idof++){
				for(unsigned int jdof=0;jdof<len_row;jdof++){
					iyval[idof] += alpha*m_valCrs_Blk[ipos0+idof*len_row+jdof]*jxval[jdof];
				}
				}
			}
		}
        return true;
    }
    
    assert( LenBlkCol() >= 0 && LenBlkRow() >= 0 );
    assert( b.Len() >= 0 );
    assert( x.Len() >= 0 );

	assert( x.Len() == m_len_BlkRow );
	assert( b.Len() == m_len_BlkCol );
	if( !isnt_trans ){
		if( LenBlkCol()*LenBlkRow()==1 ) {
			for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ b.m_Value[jblk] *= beta; }
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < m_nblk_MatRow );
					b.m_Value[jblk0] += alpha * m_valCrs_Blk[icrs]*x.m_Value[iblk];
				}
			}
		}
		else{
			std::cout << "Error!-->NotImplimented" << std::endl;
			assert(0);
		}
		return true;
	}

	if( LenBlkCol()*LenBlkRow()==1 ) {
		if( beta != 0.0 ){
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				double* lval = &(b.m_Value[iblk]);
				*lval *= beta;
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < m_nblk_MatRow );
					*lval += alpha * m_valCrs_Blk[icrs]*x.m_Value[jblk0];
				}
			}
		}
		else{
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				double* lval = &(b.m_Value[iblk]);
				*lval = 0.0;
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < m_nblk_MatRow );
					*lval += alpha * m_valCrs_Blk[icrs]*x.m_Value[jblk0];
				}
			}
		}
	}
	////////////////////////////////////////////////////////////////
	else if( this->LenBlkCol() == 1 ){
		if( this->LenBlkRow() == 2 ){
			const unsigned int nblkcol = this->NBlkMatCol();
			const unsigned int nblkrow = this->NBlkMatRow();
			const unsigned int* colind = m_colInd_Blk;
			const unsigned int* rowptr = m_rowPtr_Blk;
			const double* matval = m_valCrs_Blk;
			for(unsigned int iblk=0;iblk<nblkcol;iblk++){
				double* iyval = b.m_Value+iblk;
				*iyval *= beta;
				for(unsigned int icrs=colind[iblk];icrs<colind[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = rowptr[icrs];
					assert( jblk0 < nblkrow );
					const double* jxval = x.m_Value+jblk0*2;
					*iyval += alpha*(matval[icrs*2]*jxval[0]+matval[icrs*2+1]*jxval[1]);
				}
			}
		}
		else if( this->LenBlkRow() == 3 ){
			const unsigned int nblkcol = this->NBlkMatCol();
			const unsigned int nblkrow = this->NBlkMatRow();
			const unsigned int* colind = m_colInd_Blk;
			const unsigned int* rowptr = m_rowPtr_Blk;
			const double* matval = m_valCrs_Blk;
			for(unsigned int iblk=0;iblk<nblkcol;iblk++){
				double* iyval = b.m_Value+iblk;
				*iyval *= beta;
				for(unsigned int icrs=colind[iblk];icrs<colind[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = rowptr[icrs];
					assert( jblk0 < nblkrow );
					const double* jxval = x.m_Value+jblk0*3;
					*iyval += alpha*(matval[icrs*3  ]*jxval[0]+
                                     matval[icrs*3+1]*jxval[1]+
                                     matval[icrs*3+2]*jxval[2]);
				}
			}
		}
		else{
			const unsigned int nlen_row = this->LenBlkRow();
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				double* iyval = b.m_Value+iblk;
				*iyval *= beta;
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < m_nblk_MatRow );
					const double* jxval = x.m_Value+jblk0*nlen_row;
					for(unsigned int jdim=0;jdim<nlen_row;jdim++){
						*iyval += alpha*m_valCrs_Blk[icrs*nlen_row+jdim]*jxval[jdim];
					}
				}
			}
		}
	}
	else if( this->LenBlkRow() == 1 ){
		if( this->LenBlkCol() == 2 ){
			const unsigned int nblkcol = this->NBlkMatCol();
			const unsigned int nblkrow = this->NBlkMatRow();
			const unsigned int* colind = m_colInd_Blk;
			const unsigned int* rowptr = m_rowPtr_Blk;
			const double* matval = m_valCrs_Blk;
			for(unsigned int iblk=0;iblk<nblkcol;iblk++){
				double* iyval = b.m_Value+iblk*2;
				iyval[0] *= beta; 
				iyval[1] *= beta;
				for(unsigned int icrs=colind[iblk];icrs<colind[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = rowptr[icrs];
					assert( jblk0 < nblkrow );
					const double ajxval = *(x.m_Value+jblk0)*alpha;
					iyval[0] += matval[icrs*2  ]*ajxval;
					iyval[1] += matval[icrs*2+1]*ajxval;
				}
			}
		}
		else if( this->LenBlkCol() == 3 ){
			const unsigned int nblkcol = this->NBlkMatCol();
			const unsigned int nblkrow = this->NBlkMatRow();
			const unsigned int* colind = m_colInd_Blk;
			const unsigned int* rowptr = m_rowPtr_Blk;
			const double* matval = m_valCrs_Blk;
			for(unsigned int iblk=0;iblk<nblkcol;iblk++){
				double* iyval = b.m_Value+iblk*3;
				iyval[0] *= beta; 
				iyval[1] *= beta;
				iyval[2] *= beta;
				for(unsigned int icrs=colind[iblk];icrs<colind[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = rowptr[icrs];
					assert( jblk0 < nblkrow );
					const double ajxval = *(x.m_Value+jblk0)*alpha;
					iyval[0] += matval[icrs*3  ]*ajxval;
					iyval[1] += matval[icrs*3+1]*ajxval;
					iyval[2] += matval[icrs*3+2]*ajxval;
				}
			}
		}
		else{
			const unsigned int nlen_col = this->LenBlkCol();
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				double* iyval = b.m_Value+iblk*nlen_col;
				for(unsigned int idof=0;idof<nlen_col;idof++){ iyval[idof] *= beta; }
				for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
					assert( icrs < NCrs() );
					const unsigned int jblk0 = m_rowPtr_Blk[icrs];
					assert( jblk0 < m_nblk_MatRow );
					const double ajxval = *(x.m_Value+jblk0)*alpha;
					for(unsigned int idof=0;idof<nlen_col;idof++){
						iyval[idof] += m_valCrs_Blk[icrs*nlen_col+idof]*ajxval;
					}
				}
			}
		}
	}
	else{
		const unsigned int nblk_row = this->NBlkMatRow();
		const unsigned int nblk_col = this->NBlkMatCol();
		const unsigned int nlen_col = this->LenBlkCol();
		const unsigned int nlen_row = this->LenBlkRow();
		const unsigned int blk_size = nlen_row*nlen_col;
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			double* iyval = b.m_Value+iblk*nlen_col;
			for(unsigned int idof=0;idof<nlen_col;idof++){ iyval[idof] *= beta; }
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < NCrs() );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < nblk_row );
				const double* jxval = x.m_Value+jblk0*nlen_row;
				for(unsigned int idof=0;idof<nlen_col;idof++){
				for(unsigned int jdof=0;jdof<nlen_row;jdof++){
					iyval[idof] += alpha*m_valCrs_Blk[icrs*blk_size+idof*nlen_row+jdof]*jxval[jdof];
				}
				}
			}
		}
	}
	return true;
}


bool CMat_BlkCrs::SetPatternBoundary(const CMat_BlkCrs& rhs, 
		const CBCFlag& bc_flag_col, const CBCFlag& bc_flag_row)
{
	assert( this->NBlkMatCol() == rhs.NBlkMatCol() );
	assert( this->NBlkMatRow() == rhs.NBlkMatRow() );
	assert( this->LenBlkCol() == rhs.LenBlkCol() );
	assert( this->LenBlkRow() == rhs.LenBlkRow() );

	assert( this->NBlkMatCol() == bc_flag_col.NBlk() );
	assert( this->NBlkMatRow() == bc_flag_row.NBlk() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        return false;
    }
	const unsigned int nblk_col = this->NBlkMatCol();
	const unsigned int nblk_row = this->NBlkMatRow();
	const unsigned int nlen_col = this->LenBlkCol();
	const unsigned int nlen_row = this->LenBlkRow();

	this->DeletePattern();

	std::vector<unsigned int> flg_blk_nz_col;
	{	// 行の中一つでも自由境界条件があれば１
		flg_blk_nz_col.resize(nblk_col,0);
		for(unsigned int iblk=0;iblk<nblk_col;iblk++){
			unsigned int ilen;
			for(ilen=0;ilen<nlen_col;ilen++){
				if( bc_flag_col.GetBCFlag(iblk,ilen)==0 ) break;
			}
			if( ilen != nlen_col ) flg_blk_nz_col[iblk] = 1;
		}
	}
	std::vector<unsigned int> flg_blk_nz_row;
	{	// 列の中一つでも固定境界条件だったら１
		flg_blk_nz_row.resize(nblk_row,0);
		for(unsigned int iblk=0;iblk<nblk_row;iblk++){
			unsigned int ilen;
			for(ilen=0;ilen<nlen_row;ilen++){
				if( bc_flag_row.GetBCFlag(iblk,ilen)!=0 ) break;
			}
			if( ilen != nlen_row ) flg_blk_nz_row[iblk] = 1;
		}
	}
	////////////////
	unsigned int ncrs = 0;
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
		if( flg_blk_nz_col[iblk] == 0 ) continue;
		if( flg_blk_nz_row[iblk] == 1 ) ncrs++;
		for(unsigned int icrs=rhs.m_colInd_Blk[iblk];icrs<rhs.m_colInd_Blk[iblk+1];icrs++){
			assert( icrs < rhs.NCrs() );
			unsigned int jblk0 = rhs.m_rowPtr_Blk[icrs];
			assert( jblk0 < nblk_row );
			if( flg_blk_nz_row[jblk0] == 1 ) ncrs++;
		}
	}

	this->m_ncrs_Blk = ncrs;
	this->m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];

	ncrs = 0;
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
		if( flg_blk_nz_col[iblk] == 0 ){
			m_colInd_Blk[iblk+1] = m_colInd_Blk[iblk];
			continue;
		}
		const unsigned int ncrs_begin = m_colInd_Blk[iblk];
		unsigned int inz = 0;
		if( flg_blk_nz_row[iblk] == 1 ){
			assert( ncrs_begin+inz < m_ncrs_Blk );
			m_rowPtr_Blk[ncrs_begin+inz] = iblk;
			inz++;
		}
		for(unsigned int icrs=rhs.m_colInd_Blk[iblk];icrs<rhs.m_colInd_Blk[iblk+1];icrs++){
			assert( icrs < rhs.NCrs() );
			unsigned int jblk0 = rhs.m_rowPtr_Blk[icrs];
			assert( jblk0 < nblk_row );
			if( flg_blk_nz_row[jblk0] == 0 ) continue;
			assert( ncrs_begin+inz < m_ncrs_Blk );
			m_rowPtr_Blk[ncrs_begin+inz] = jblk0;
			inz++;
		}
		m_colInd_Blk[iblk+1] = m_colInd_Blk[iblk] + inz;
		ncrs += inz;
	}
	assert( m_ncrs_Blk == ncrs );
	return true;
}

bool CMat_BlkCrs::SetPatternDia(const CMat_BlkCrs& rhs)
{
	assert( rhs.NBlkMatCol() == rhs.NBlkMatRow() );
	assert( this->NBlkMatCol() == rhs.NBlkMatCol() );
	assert( this->NBlkMatRow() == rhs.NBlkMatRow() );

	assert( rhs.LenBlkCol() == rhs.LenBlkRow() );
	assert( this->LenBlkCol() == rhs.LenBlkCol() );
	assert( this->LenBlkRow() == rhs.LenBlkRow() );

    if( LenBlkCol() == -1 || LenBlkRow() == -1 ){
        std::cout << "Error!-->Not Implemented" << std::endl;
        assert(0);
        return false;
    }

	const unsigned int nblk = this->NBlkMatCol();
//	const unsigned int nlen = this->LenBlkCol();

	this->DeletePattern();

	this->m_ncrs_Blk = rhs.NCrs() + nblk;
	this->m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];

	unsigned int ncrs = 0;
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		const unsigned int ncrs_begin = m_colInd_Blk[iblk];
		unsigned int inz = 0;
		bool flg_dia = false;
		for(unsigned int icrs=rhs.m_colInd_Blk[iblk];icrs<rhs.m_colInd_Blk[iblk+1];icrs++){
			assert( icrs < rhs.NCrs() );
			unsigned int jblk0 = rhs.m_rowPtr_Blk[icrs];
			assert( jblk0 < nblk );
			assert( ncrs_begin+inz < m_ncrs_Blk );
			if( iblk > jblk0 && !flg_dia ){
				m_rowPtr_Blk[ncrs_begin+inz] = iblk;
				inz++;
				flg_dia = true;
			}
			m_rowPtr_Blk[ncrs_begin+inz] = jblk0;
			inz++;
		}
		if( !flg_dia ){
			m_rowPtr_Blk[ncrs_begin+inz] = iblk;
			inz++;
		}
		m_colInd_Blk[iblk+1] = m_colInd_Blk[iblk] + inz;
		ncrs += inz;
	}
	assert( m_ncrs_Blk == ncrs );

	return true;
}

