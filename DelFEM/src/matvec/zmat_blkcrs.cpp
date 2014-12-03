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

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <cstdlib> //(abort)

#include "delfem/matvec/zmat_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/bcflag_blk.h"
#include "delfem/indexed_array.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CZMat_BlkCrs::CZMat_BlkCrs()
: m_nblk_MatCol(0), m_nblk_MatRow(0), m_len_BlkCol(0), m_len_BlkRow(0)
{
//	std::cout << "Construct : CZMat_BlkCrs (Default Constructor) " << std::endl;
	m_colInd_Blk = 0;
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0;
}

CZMat_BlkCrs::CZMat_BlkCrs(const unsigned int nblk_col, const unsigned int len_col, 
				   const unsigned int nblk_row, const unsigned int len_row )
: m_nblk_MatCol(nblk_col), m_nblk_MatRow(nblk_row), m_len_BlkCol(len_col), m_len_BlkRow(len_row)
{
//	std::cout << "Construct : CZMat_BlkCrs " << std::endl;
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 
}

CZMat_BlkCrs::CZMat_BlkCrs(const CZMat_BlkCrs& rhs, const bool is_value, const bool isnt_trans, const bool isnt_conj)
: m_nblk_MatCol( isnt_trans ? rhs.m_nblk_MatCol : rhs.m_nblk_MatRow ),
m_nblk_MatRow( isnt_trans ? rhs.m_nblk_MatRow : rhs.m_nblk_MatCol ),
m_len_BlkCol( isnt_trans ? rhs.m_len_BlkCol : rhs.m_len_BlkRow ),
m_len_BlkRow( isnt_trans ? rhs.m_len_BlkRow : rhs.m_len_BlkCol )
{
//	std::cout << "Construct : CZMat_BlkCrs" << std::endl;
	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 
	AddPattern(rhs,isnt_trans);
	if( is_value ){ this->SetValue(rhs,isnt_trans, isnt_conj); }
}

void CZMat_BlkCrs::Initialize(const unsigned int nblk_col, const unsigned int len_col, 
						 const unsigned int nblk_row, const unsigned int len_row ){

	if( m_colInd_Blk != 0 ){ delete[] m_colInd_Blk; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; }

	m_nblk_MatCol = nblk_col;
	m_nblk_MatRow = nblk_row;
	m_len_BlkCol = len_col;
	m_len_BlkRow = len_row;

	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	m_ncrs_Blk = 0;
	m_rowPtr_Blk = 0;
	m_valCrs_Blk = 0; 
}


CZMat_BlkCrs::~CZMat_BlkCrs()
{
	if( m_colInd_Blk != 0 ){ delete[] m_colInd_Blk;	m_colInd_Blk = 0; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk;	m_rowPtr_Blk = 0; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk;	m_valCrs_Blk = 0; }
}

bool CZMat_BlkCrs::DeletePattern(){
	m_ncrs_Blk = 0;
	for(unsigned int iblk=0;iblk<m_nblk_MatCol+1;iblk++){ m_colInd_Blk[iblk] = 0; }
	if( m_rowPtr_Blk != 0 ){ delete[] m_rowPtr_Blk; m_rowPtr_Blk = 0; }
	if( m_valCrs_Blk != 0 ){ delete[] m_valCrs_Blk; m_valCrs_Blk = 0; }
	return true;
}

bool CZMat_BlkCrs::SetZero(){
	const unsigned int ntotdof = m_ncrs_Blk*m_len_BlkCol*m_len_BlkRow;
	if( m_valCrs_Blk == 0 ){	// まだ値の領域をを確保していない場合
		m_valCrs_Blk = new Com::Complex [ntotdof];
	}
	for(unsigned int idof=0;idof<ntotdof;idof++){ m_valCrs_Blk[idof] = 0.0;	}
	return true;
}

bool CZMat_BlkCrs::AddPattern(const CZMat_BlkCrs& rhs, const bool isnt_trans)
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
			this->m_ncrs_Blk = rhs.m_ncrs_Blk;
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
			for(int iblk=m_nblk_MatCol;iblk>=1;iblk--){
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

bool CZMat_BlkCrs::SetValue(const CZMat_BlkCrs& rhs, const bool isnt_trans, const bool isnt_conj){
	const unsigned int BlkSize = m_len_BlkCol*m_len_BlkRow;
	this->SetZero();

	if( isnt_trans ){	// Set Value -- Not Transpose Pattern
		assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
		assert( m_nblk_MatRow == rhs.m_nblk_MatRow );
		assert( m_len_BlkCol == rhs.m_len_BlkCol );
		assert( m_len_BlkRow == rhs.m_len_BlkRow );
		int* row2crs = new int [m_nblk_MatRow];
		for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ 
			row2crs[jblk] = -1; 
		}
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			for(unsigned int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[ijcrs];
				assert( jblk0 < m_nblk_MatRow );
				row2crs[jblk0] = ijcrs;
			}
			////////////////
			for(unsigned int ijcrs=rhs.m_colInd_Blk[iblk];ijcrs<rhs.m_colInd_Blk[iblk+1];ijcrs++){
				assert( ijcrs<rhs.m_ncrs_Blk );
				const unsigned int jblk0 = rhs.m_rowPtr_Blk[ijcrs];
				assert( jblk0<m_nblk_MatRow );
				const int ijcrs0 = row2crs[jblk0];
				if( ijcrs0 == -1 ) continue;
				const Com::Complex* pval_in = &rhs.m_valCrs_Blk[ijcrs*BlkSize];
				Com::Complex* pval_out = &m_valCrs_Blk[ijcrs0*BlkSize];
				if( isnt_conj ){
					for(unsigned int idof=0;idof<BlkSize;idof++){ *(pval_out+idof) = *(pval_in+idof); }
				}
				else{
					for(unsigned int idof=0;idof<BlkSize;idof++){ *(pval_out+idof) = Com::Conjugate(*(pval_in+idof)); }
				}
			}
			////////////////
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
				////////////////
				for(unsigned int jicrs=rhs.m_colInd_Blk[jblk0];jicrs<rhs.m_colInd_Blk[jblk0+1];jicrs++){
					assert( jicrs<rhs.m_ncrs_Blk );
					const unsigned int iblk0 = rhs.m_rowPtr_Blk[jicrs];
					if( iblk != iblk0 ) continue;
					const Com::Complex* pval_in = &rhs.m_valCrs_Blk[jicrs*BlkSize];
					Com::Complex* pval_out = &m_valCrs_Blk[ijcrs*BlkSize];
					if( isnt_conj ){
						for(unsigned int idof=0;idof<m_len_BlkCol;idof++){
							for(unsigned int jdof=0;jdof<m_len_BlkRow;jdof++){
								*(pval_out+idof*m_len_BlkRow+jdof) = *(pval_in+jdof*m_len_BlkCol+idof);
							}
						}
					}
					else{
						for(unsigned int idof=0;idof<m_len_BlkCol;idof++){
							for(unsigned int jdof=0;jdof<m_len_BlkRow;jdof++){
								*(pval_out+idof*m_len_BlkRow+jdof) = Com::Conjugate(*(pval_in+jdof*m_len_BlkCol+idof));
							}
						}
					}
				}
				////////////////
			}
		}
	}
	return true;
}

bool CZMat_BlkCrs::MatVec(double alpha, const CZVector_Blk& x, double beta, CZVector_Blk& b, const bool isnt_trans) const
{
	
	////////////////
	if( this->m_len_BlkCol != 1 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}
	////////////////

	assert( x.BlkVecLen() == m_nblk_MatRow );
	assert( x.BlkLen() == m_len_BlkRow );
	assert( b.BlkVecLen() == m_nblk_MatCol );
	assert( b.BlkLen() == m_len_BlkCol );

	if( isnt_trans ){
		if( ((beta>=0.0)?beta:-beta) > 1.0e-30 ){
			for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
				Com::Complex* lval = &(b.m_Value[iblk]);
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
				Com::Complex* lval = &(b.m_Value[iblk]);
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
	else{
		for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){
			b.m_Value[jblk] *= beta;
		}
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < m_nblk_MatRow );
				b.m_Value[jblk0] += alpha * m_valCrs_Blk[icrs]*x.m_Value[iblk];
			}
		}
	}
	return true;
}


// bc_flagが０でない自由度に固定境界条件をセットする
bool CZMat_BlkCrs::SetBoundaryCondition_Colum(const CBCFlag& bc_flag)
{
    assert( bc_flag.LenBlk() >= 0 );
	assert( bc_flag.NBlk() == NBlkMatRow() );
    assert( bc_flag.LenBlk() == (int)LenBlkRow() );

	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int lenBlkCol = LenBlkCol();
	const unsigned int lenBlkRow = LenBlkRow();

	const unsigned int ncrs = this->NCrs();
	
	// bc_flagが０でない行列の列を０にする。
	for(unsigned int icrs=0;icrs<ncrs;icrs++){
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
bool CZMat_BlkCrs::SetBoundaryCondition_Row(const CBCFlag& bc_flag)
{
    assert( bc_flag.LenBlk() >= 0 );
	assert( bc_flag.NBlk() == NBlkMatCol() );
    assert( bc_flag.LenBlk() == (int)LenBlkCol() );

	const unsigned int BlkSize = LenBlkCol()*LenBlkRow();
	const unsigned int lenBlkCol = LenBlkCol();
	const unsigned int lenBlkRow = LenBlkRow();
	
	// bc_flagが０でない行列の行を０にする。
	// 但し対角成分は１にしておく
	for(unsigned int iblk=0;iblk<NBlkMatCol();iblk++){
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


bool CZMat_BlkCrs::AddPattern(const Com::CIndexedArray& crs)
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
