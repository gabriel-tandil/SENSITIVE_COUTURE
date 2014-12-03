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

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdio.h>

#include "delfem/complex.h"
#include "delfem/indexed_array.h"

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CZMatDia_BlkCrs::CZMatDia_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow)
    :CZMat_BlkCrs(nblk_colrow, len_colrow, nblk_colrow, len_colrow)
{
//	std::cout << "Construct : CZMatDia_BlkCrs(nblk_colrow, len_colrow) " << nblk_colrow << " " << len_colrow << std::endl;
	m_valDia_Blk = new Com::Complex [m_nblk_MatCol*m_len_BlkCol*m_len_BlkRow];
}

CZMatDia_BlkCrs::CZMatDia_BlkCrs(const CZMatDia_BlkCrs& rhs, const bool is_value, const bool isnt_trans, const bool isnt_conj)
    : CZMat_BlkCrs(rhs,is_value,isnt_trans,isnt_conj)
{
//	std::cout << "Construct : CZMat_BlkCrs" << std::endl;

	assert( m_nblk_MatRow == rhs.NBlkMatRow() );
	assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
	assert( m_nblk_MatCol == m_nblk_MatRow );
	
	assert( m_len_BlkCol == rhs.LenBlkCol() );
	assert( m_len_BlkRow == rhs.LenBlkRow() );
	assert( m_len_BlkCol == m_len_BlkRow );
	
	m_valDia_Blk = new Com::Complex [m_nblk_MatCol*m_len_BlkCol*m_len_BlkRow];

	if( is_value ){ 
		if( isnt_trans ){
			const unsigned int ndof=m_nblk_MatCol*m_len_BlkCol*m_len_BlkRow;
			if( isnt_conj ){
				for(unsigned int idof=0;idof<ndof;idof++){
					m_valDia_Blk[idof] = rhs.m_valDia_Blk[idof];
				}
			}
			else{
				for(unsigned int idof=0;idof<ndof;idof++){
					m_valDia_Blk[idof] = Com::Conjugate(rhs.m_valDia_Blk[idof]);
				}
			}
		}
		else{
			const unsigned int BlkLen = m_len_BlkCol;
			const unsigned int BlkSize = BlkLen*BlkLen;
			if( isnt_conj ){
				for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
					for(unsigned int idof=0;idof<BlkLen;idof++){
					for(unsigned int jdof=0;jdof<BlkLen;jdof++){
						m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] 
							= rhs.m_valDia_Blk[iblk*BlkSize+jdof*BlkLen*idof];
					}
					}
				}
			}
			else{
				for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
					for(unsigned int idof=0;idof<BlkLen;idof++){
					for(unsigned int jdof=0;jdof<BlkLen;jdof++){
						m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] 
//							= rhs.m_valDia_Blk[iblk*BlkSize+jdof*BlkLen*idof];
							= Com::Conjugate(rhs.m_valDia_Blk[iblk*BlkSize+jdof*BlkLen*idof]);
					}
					}
				}
			}
		}
	}
}

CZMatDia_BlkCrs::CZMatDia_BlkCrs(const std::string& file_path)
    :CZMat_BlkCrs()
{
	FILE* fp;
	
	if ((fp = fopen(file_path.c_str(), "r")) == NULL) {
		printf("file open error!!\n");
		exit(1);
	}

	int nblk_col, nnz;
	fscanf(fp,"%d %d",&nblk_col, &nnz);
//	std::cout << nblk_col << " " << nnz << std::endl;
	assert( nblk_col >= 0 );
	assert( nnz >= 0 );

	CZMat_BlkCrs::Initialize(nblk_col,1, nblk_col,1);

	this->m_rowPtr_Blk = new unsigned int [nnz];	
	this->m_valCrs_Blk = new Com::Complex [nnz];
	this->m_valDia_Blk = new Com::Complex [nblk_col];

	this->m_colInd_Blk[0] = 0;
    unsigned int icol_cur = 0;
    unsigned int icrs_cur = 0;
	bool iflag0 = false;
    for(unsigned int inz=0;inz<(unsigned int)nnz;inz++){
		int icol, irow;
		double real, imag;
		fscanf(fp,"%d %d %lf %lf", &icol, &irow, &real, &imag);
		assert( icol >= 0 && icol < nblk_col );
		assert( irow >= 0 && irow < nblk_col );
        if( icol_cur != (unsigned int)icol ){
            assert( icol_cur < (unsigned int)icol );
            for(unsigned int i=icol_cur+1;i<=(unsigned int)icol;i++){ this->m_colInd_Blk[i] = icrs_cur; }
			if( !iflag0 ){ this->m_valDia_Blk[icol_cur] = 0.0; }
			iflag0 = false;
            icol_cur = (unsigned int)icol;
		}
		if( icol == irow ){
			iflag0 = true;
			this->m_valDia_Blk[icol] = Com::Complex(real,imag);
		}
		else{
			this->m_rowPtr_Blk[icrs_cur] = irow;
			this->m_valCrs_Blk[icrs_cur] = Com::Complex(real,imag);
			icrs_cur++;
		}
	}
	this->m_ncrs_Blk = icrs_cur;
    for(unsigned int i=icol_cur+1;i<=(unsigned int)nblk_col;i++){
        this->m_colInd_Blk[nblk_col] = this->m_ncrs_Blk;
    }
	fclose(fp);

/*
	for(unsigned int iblk=0;iblk<nblk_col;iblk++){
		std::cout << iblk << "-->";
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			unsigned int jblk = m_rowPtr_Blk[icrs];
			std::cout << jblk << " ";
		}
		std::cout << std::endl;
	}
*/
}

CZMatDia_BlkCrs::~CZMatDia_BlkCrs()
{
	if( m_valDia_Blk != 0 ){ delete[] m_valDia_Blk; m_valDia_Blk = 0; }
}

bool CZMatDia_BlkCrs::DeletePattern(){
	CZMat_BlkCrs::DeletePattern();
	return true;
}

bool CZMatDia_BlkCrs::SetZero(){
	CZMat_BlkCrs::SetZero();
	const unsigned int ndof = m_nblk_MatCol*m_len_BlkCol*m_len_BlkCol;
	for(unsigned int idof=0;idof<ndof;idof++){ m_valDia_Blk[idof] = 0.0; }
	return true;
}

void CZMatDia_BlkCrs::AddUnitMatrix(const Com::Complex& epsilon){
	assert( m_len_BlkCol == 1 );
	for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){ m_valDia_Blk[iblk] += epsilon; }
	return;
}

bool CZMatDia_BlkCrs::Mearge(
			unsigned int nblkel_col, const unsigned int* blkel_col,
			unsigned int nblkel_row, const unsigned int* blkel_row,
			unsigned int blksize, const Com::Complex* emat, int* tmp_buffer)
{
	assert( m_valCrs_Blk != 0 );
	assert( m_valDia_Blk != 0 );

	const unsigned int BlkSize = m_len_BlkCol*m_len_BlkRow;
//	const unsigned int BlkLen = m_len_BlkCol;

	assert( nblkel_col == nblkel_row );
	assert( blksize == BlkSize );

    for(unsigned int iblkel=0;iblkel<nblkel_col;iblkel++){
		const int iblk1 = blkel_col[iblkel];
		for(unsigned int jpsup=m_colInd_Blk[iblk1];jpsup<m_colInd_Blk[iblk1+1];jpsup++){
			assert( jpsup < m_ncrs_Blk );
			const unsigned int jblk1 = m_rowPtr_Blk[jpsup];
			tmp_buffer[jblk1] = jpsup;
		}
		for(unsigned int jblkel=0;jblkel<nblkel_row;jblkel++){
			if( iblkel == jblkel ){	// Marge Diagonal
				const Com::Complex* pval_in = &emat[(iblkel*nblkel_row+iblkel)*BlkSize];
				Com::Complex* pval_out = &m_valDia_Blk[iblk1*BlkSize];
				for(unsigned int idof=0;idof<BlkSize;idof++){ pval_out[idof] += pval_in[idof]; }
			}
			else{	// Marge Non-Diagonal
				const unsigned int jblk1 = blkel_row[jblkel];
				const unsigned int jpsup1 = tmp_buffer[jblk1];
				assert( jpsup1 < m_ncrs_Blk );
				assert( m_rowPtr_Blk[jpsup1] == jblk1 );
				const Com::Complex* pval_in = &emat[(iblkel*nblkel_row+jblkel)*BlkSize];
				Com::Complex* pval_out = &m_valCrs_Blk[jpsup1*BlkSize];
				for(unsigned int idof=0;idof<BlkSize;idof++){ pval_out[idof] += pval_in[idof]; }
			}
		}
	}
	return true;
}

bool CZMatDia_BlkCrs::SetBoundaryCondition(const MatVec::CBCFlag& bc_flag){

	assert( m_len_BlkCol == m_len_BlkRow );
	const unsigned int BlkSize = m_len_BlkCol*m_len_BlkRow;
	const unsigned int BlkLen = m_len_BlkCol;
	
	for(unsigned int iblk=0;iblk<m_nblk_MatRow;iblk++){
		for(unsigned int idof=0;idof<BlkLen;idof++){
			if( bc_flag.GetBCFlag(iblk,idof) == 0 ) continue;
			{
				for(unsigned int kdof=0;kdof<BlkLen;kdof++){ m_valDia_Blk[iblk*BlkSize+kdof*BlkLen+idof] = 0.0; }
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){ m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] = 0.0; }
				m_valDia_Blk[iblk*BlkSize+idof*BlkLen+idof] = 1.0;
			}
			for(unsigned int ipsup=m_colInd_Blk[iblk];ipsup<m_colInd_Blk[iblk+1];ipsup++){
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){ m_valCrs_Blk[ipsup*BlkSize+idof*BlkLen+jdof] = 0.0; }
			}
		}
	}
	for(unsigned int ipsup=0;ipsup<m_colInd_Blk[m_nblk_MatRow];ipsup++){
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

bool CZMatDia_BlkCrs::AddPattern(const CZMatDia_BlkCrs& rhs, const bool isnt_trans){
	if(  isnt_trans ){
		assert( m_nblk_MatCol == rhs.m_nblk_MatRow );
		assert( m_nblk_MatRow == rhs.m_nblk_MatCol );
		if( m_ncrs_Blk == 0 ){
			DeletePattern();
			CZMat_BlkCrs::AddPattern(rhs,isnt_trans);
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

bool CZMatDia_BlkCrs::AddPattern(const CZMat_BlkCrs& m1, const CZMatDia_BlkCrs& m2, const CZMat_BlkCrs& m3){
	assert( m_nblk_MatCol  == m1.NBlkMatCol() );
	assert( m1.NBlkMatRow() == m2.NBlkMatCol() );
	assert( m2.NBlkMatRow() == m3.NBlkMatCol() );
	assert( m3.NBlkMatRow() == m_nblk_MatRow );

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
		unsigned int* is_j_flag = new unsigned int [m_nblk_MatRow];
		for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ is_j_flag[jblk] = 0; }
		unsigned int j_size = 0;
		unsigned int* j_buffer = new unsigned int [m_nblk_MatRow];
		////////////////
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
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
						assert( jblk0 < m_nblk_MatRow );
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

// 非ゼロパターンを加える
bool CZMatDia_BlkCrs::AddPattern(const Com::CIndexedArray& crs)
{
	// 入力チェック
	assert( crs.CheckValid() );
	if( !crs.CheckValid() ) return false;
	if( crs.Size() > NBlkMatCol() ) return false;
	{
		unsigned int max_val = 0;
		for(unsigned int i=0;i<crs.array.size();i++){
			max_val = (max_val>crs.array[i])?max_val:crs.array[i];
		}
		if( max_val>NBlkMatRow() ) return false;
	}

	// 入力行列が０だったら何もしなくてよい
	if( crs.array.size() == 0 ) return true;

	if( m_ncrs_Blk == 0 ){	// 何も要素が入っていない行列にパターンを追加する
		CZMatDia_BlkCrs::DeletePattern();
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

bool CZMatDia_BlkCrs::SetValue(const CZMatDia_BlkCrs& rhs, const bool isnt_trans){
	assert( m_nblk_MatRow == rhs.m_nblk_MatRow );
	assert( m_nblk_MatCol == rhs.m_nblk_MatCol );
	assert( m_nblk_MatCol == m_nblk_MatRow );
	
	assert( m_len_BlkCol == rhs.LenBlkCol() );
	assert( m_len_BlkRow == rhs.LenBlkRow() );
	assert( m_len_BlkCol == m_len_BlkRow );

	CZMat_BlkCrs::SetValue(rhs,isnt_trans,true);
	if( isnt_trans ){
		const unsigned int ndof=m_nblk_MatCol*m_len_BlkCol*m_len_BlkRow;
		for(unsigned int idof=0;idof<ndof;idof++){
			m_valDia_Blk[idof] = rhs.m_valDia_Blk[idof];
		}
	}
	else{
		const unsigned int BlkLen = m_len_BlkCol;
		const unsigned int BlkSize = BlkLen*BlkLen;
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
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

bool CZMatDia_BlkCrs::SetValue(const CZMat_BlkCrs& m1, const CZMatDia_BlkCrs& m2, const CZMat_BlkCrs& m3){
	assert( m_nblk_MatRow == m_nblk_MatRow );
	assert( m_nblk_MatCol  == m1.NBlkMatCol() );
	assert( m1.NBlkMatRow() == m2.NBlkMatCol() );
	assert( m2.NBlkMatRow() == m3.NBlkMatCol() );
	assert( m3.NBlkMatRow() ==    m_nblk_MatRow );

	if( m_len_BlkCol != 1 || m2.LenBlkCol() != 1 || m_len_BlkRow != 1 ){
		std::cout << "Error!-->Not Implimented " << std::endl;
		assert(0);
		abort();
	}

	const unsigned int mat_len_mid = m2.NBlkMatCol();

	if( m_valCrs_Blk == 0 ){	
		m_valCrs_Blk = new Com::Complex [m_ncrs_Blk];
	}
	
	////////////////
	unsigned int* is_l_flag = new unsigned int [mat_len_mid];
	for(unsigned int lblk=0;lblk<mat_len_mid;lblk++){ is_l_flag[lblk] = 0; }
	unsigned int l_size = 0;
	unsigned int* l_buffer = new unsigned int [mat_len_mid];
	Com::Complex* l_val = new Com::Complex [mat_len_mid];
	////////////////
	unsigned int* is_j_flag = new unsigned int [m_nblk_MatRow];
	for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ is_j_flag[jblk] = 0; }
	unsigned int j_size = 0;
	unsigned int* j_buffer = new unsigned int [m_nblk_MatRow];
	Com::Complex* j_val = new Com::Complex [m_nblk_MatRow];
	////////////////
    int* row2crs = new int [m_nblk_MatRow];
	for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){ 
		row2crs[jblk] = -1; 
	}
	for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
		{	// Make l_buffer, l_size
			l_size = 0;
			unsigned int npsuk = 0;
			const unsigned int* ind_psuk = m1.GetPtrIndPSuP(iblk,npsuk);
			assert( ind_psuk != 0 );
			const Com::Complex* val_psuk = m1.GetPtrValPSuP(iblk,npsuk);
			assert( val_psuk != 0 );
			for(unsigned int k=0;k<npsuk;k++){
				const unsigned int kblk0 = ind_psuk[k];
				assert( kblk0 < mat_len_mid );
				const Com::Complex& ik_val = val_psuk[k];
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
				const Com::Complex& lj_val = l_val[l];
				unsigned int npsul = 0;
				const unsigned int* ind_psul = m3.GetPtrIndPSuP(lblk0,npsul);
				assert( ind_psul != 0 );
				const Com::Complex* val_psul = m3.GetPtrValPSuP(lblk0,npsul);
				assert( val_psul != 0 );
				for(unsigned int l=0;l<npsul;l++){
					const unsigned int jblk0 = ind_psul[l];
					assert( jblk0 < m_nblk_MatRow );
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
			assert( jblk0 < m_nblk_MatRow );
			row2crs[jblk0] = ijcrs;
		}
		for(unsigned int j=0;j<j_size;j++){
			const unsigned int jblk0 = j_buffer[j];
			assert( jblk0<m_nblk_MatRow );
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
			assert( jblk0 < m_nblk_MatRow );
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
bool CZMatDia_BlkCrs::MatVec(double alpha, const CZVector_Blk& x, double beta, CZVector_Blk& y) const
{
	assert( m_nblk_MatCol == m_nblk_MatRow );

	assert( x.BlkVecLen() == m_nblk_MatRow );
    assert( x.BlkLen()  == m_len_BlkRow );

	assert( y.BlkVecLen() == m_nblk_MatCol );
    assert( y.BlkLen()  == m_len_BlkCol );

	const unsigned int BlkLen = m_len_BlkCol;
	const unsigned int BlkSize = BlkLen*BlkLen;

	if( BlkSize == 1 ){
		const Com::Complex* xval = x.m_Value;
		Com::Complex* yval = y.m_Value;
		const unsigned int nblk = m_nblk_MatCol;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			double real0 = 0;
			double imag0 = 0;
			const unsigned int crs0 = m_colInd_Blk[iblk];
			const unsigned int crs1 = m_colInd_Blk[iblk+1];
			for(unsigned int icrs=crs0;icrs<crs1;icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < m_nblk_MatRow );
				real0 += alpha * (m_valCrs_Blk[icrs].Real()*xval[jblk0].Real() - m_valCrs_Blk[icrs].Imag()*xval[jblk0].Imag());
				imag0 += alpha * (m_valCrs_Blk[icrs].Real()*xval[jblk0].Imag() + m_valCrs_Blk[icrs].Imag()*xval[jblk0].Real());
			}
			Com::Complex* iyval = yval+iblk;
			(*iyval) *= beta;
			(*iyval) += Com::Complex(real0,imag0);
			(*iyval) += alpha*m_valDia_Blk[iblk]*xval[iblk];
		}
	}
	else{
		const Com::Complex* xval = x.m_Value;
		Com::Complex* yval = y.m_Value;
		for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
			Com::Complex* iyval = yval+iblk*BlkLen;
			for(unsigned int idof=0;idof<BlkLen;idof++){
				iyval[idof] *= beta;
			}
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jblk0 = m_rowPtr_Blk[icrs];
				assert( jblk0 < m_nblk_MatRow );
				const Com::Complex* jxval = xval+jblk0*BlkLen;
				for(unsigned int idof=0;idof<BlkLen;idof++){
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){
					iyval[idof] += alpha * m_valCrs_Blk[icrs*BlkSize+idof*BlkLen+jdof] * jxval[jdof];
				}
				}
			}
			const Com::Complex* ixval = xval + iblk*BlkLen;
			for(unsigned int idof=0;idof<BlkLen;idof++){
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){
				iyval[idof] += alpha * m_valDia_Blk[iblk*BlkSize+idof*BlkLen+jdof] * ixval[jdof];
			}
			}
		}
	}

	return true;
}

bool CZMatDia_BlkCrs::MatVec_Hermitian(double alpha, const CZVector_Blk& x, double beta, CZVector_Blk& y) const{
	assert( m_nblk_MatCol == m_nblk_MatRow );

	assert( x.BlkVecLen() == m_nblk_MatRow );
    assert( x.BlkLen()  == m_len_BlkRow );

	assert( y.BlkVecLen() == m_nblk_MatCol );
    assert( y.BlkLen()  == m_len_BlkCol );

	const unsigned int BlkLen = m_len_BlkCol;
	const unsigned int BlkSize = BlkLen*BlkLen;

	const Com::Complex* xval = x.m_Value;
	Com::Complex* yval = y.m_Value;

	for(unsigned int iblk=0;iblk<m_nblk_MatCol;iblk++){
		Com::Complex* iyval = yval+iblk*BlkLen;
		{
			for(unsigned int idof=0;idof<BlkLen;idof++){
				iyval[idof] *= beta;
			}
		}
	}
	for(unsigned int jblk=0;jblk<m_nblk_MatRow;jblk++){
		const Com::Complex* jxval = xval+jblk*BlkLen;
		for(unsigned int jicrs=m_colInd_Blk[jblk];jicrs<m_colInd_Blk[jblk+1];jicrs++){
			assert( jicrs < m_ncrs_Blk );
			const unsigned int iblk0 = m_rowPtr_Blk[jicrs];
			assert( iblk0 < m_nblk_MatCol );
			{
				Com::Complex* iyval = yval+iblk0*BlkLen;
				for(unsigned int jdof=0;jdof<BlkLen;jdof++){
				for(unsigned int idof=0;idof<BlkLen;idof++){
					iyval[idof] += alpha * Com::InnerProduct(jxval[jdof],m_valCrs_Blk[jicrs*BlkSize+jdof*BlkLen+idof]);
//					iyval[idof] += alpha * m_valCrs_Blk[jicrs*BlkSize+jdof*BlkLen+idof]*jxval[jdof];
				}
				}
			}
		}
		{
			Com::Complex* jyval = yval+jblk*BlkLen;
			for(unsigned int jdof=0;jdof<BlkLen;jdof++){
			for(unsigned int idof=0;idof<BlkLen;idof++){
				jyval[idof] += alpha * Com::InnerProduct(jxval[jdof],m_valDia_Blk[jblk*BlkSize+jdof*BlkLen+idof]);
//				jyval[idof] += alpha * m_valDia_Blk[jblk*BlkSize+jdof*BlkLen+idof]*jxval[jdof];
			}
			}
		}
	}
	return true;
}
