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

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )   // C4786なんて表示すんな( ﾟДﾟ)ｺﾞﾙｧ
#endif
#define for if(0); else for


#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/bcflag_blk.h"

#include "delfem/ls/linearsystem.h"

void LsSol::CLinearSystem::Clear()
{
	for(unsigned int i=0;i<m_Matrix_NonDia.size();i++){
		for(unsigned int j=0;j<m_Matrix_NonDia[i].size();j++){
			delete m_Matrix_NonDia[i][j];
		}
		m_Matrix_NonDia[i].clear();
	}
	m_Matrix_NonDia.clear();
	////////////////
	for(unsigned int i=0;i<m_Matrix_Dia.size();i++){
		delete m_Matrix_Dia[i];
	}
	m_Matrix_Dia.clear();
	////////////////
	for(unsigned int i=0;i<m_Residual.size();i++){
		delete m_Residual[i];
	}
	m_Residual.clear();
	////////////////
	for(unsigned int i=0;i<m_Update.size();i++){
		delete m_Update[i];
	}
	m_Update.clear();
	////////////////
	for(unsigned int i=0;i<m_BCFlag.size();i++){
		delete m_BCFlag[i];
	}
	m_BCFlag.clear();
	////////////////
	for(unsigned int i=0;i<m_TmpVectorArray.size();i++){
		for(unsigned int j=0;j<m_TmpVectorArray[i].size();j++){
			delete m_TmpVectorArray[i][j];
		}
		m_TmpVectorArray[i].clear();
	}
	m_TmpVectorArray.clear();
    ////////////////
    m_aSeg.clear();
}

// マージ前の初期化
void LsSol::CLinearSystem::InitializeMarge()
{	
	unsigned int nls = m_aSeg.size();
	for(unsigned int ils=0;ils<nls;ils++){
		for(unsigned int jls=0;jls<nls;jls++){
            if( m_Matrix_NonDia[ils][jls] == 0 ){ continue; }
		    m_Matrix_NonDia[ils][jls]->SetZero();
		}
		if( m_Matrix_Dia[ils] != 0 ){ m_Matrix_Dia[ils]->SetZero(); }
	}
	for(unsigned int ils=0;ils<nls;ils++){
		m_Residual[ils]->SetVectorZero();
		m_Update[ils]->SetVectorZero();
	}
}


// マージ後の処理（残差ノルムを返す)
double LsSol::CLinearSystem::FinalizeMarge()
{
	unsigned int nseg = m_aSeg.size();
	{	// 境界条件の行列へのセット
		for(unsigned int iseg=0;iseg<nseg;iseg++){
			if( m_Matrix_Dia[iseg] != 0 ){
				m_Matrix_Dia[iseg]->SetBoundaryCondition(*m_BCFlag[iseg]);
			}
		}
		for(unsigned int iseg=0;iseg<nseg;iseg++){
			for(unsigned int jseg=0;jseg<nseg;jseg++){
				if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
				assert( iseg != jseg );
				m_Matrix_NonDia[iseg][jseg]->SetBoundaryCondition_Row(*m_BCFlag[iseg]);
			}
		}
		for(unsigned int jseg=0;jseg<nseg;jseg++){
			for(unsigned int iseg=0;iseg<nseg;iseg++){
				if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
				assert( iseg != jseg );
				m_Matrix_NonDia[iseg][jseg]->SetBoundaryCondition_Colum(*m_BCFlag[jseg]);
			}
		}
	}

	// 残差への境界条件のセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		m_BCFlag[iseg]->SetZeroToBCDof(*m_Residual[iseg]);
	}

	// 残差ノルムの計算
	double sq_norm_res = 0.0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		sq_norm_res += m_Residual[iseg]->GetSquaredVectorNorm();
	}
	return sqrt(sq_norm_res);
}


void LsSol::CLinearSystem::ClearFixedBoundaryCondition(){
	for(unsigned int ibcflag=0;ibcflag<m_BCFlag.size();ibcflag++){
		m_BCFlag[ibcflag]->SetAllFlagZero();
	}
}

// 連立一次方程式セグメントを一番最後に加える。（行列、ベクトルなどをリサイズ)
int LsSol::CLinearSystem::AddLinSysSeg( unsigned int nnode, unsigned int len )
{
	const unsigned int size_old = m_aSeg.size();
	const unsigned int size_new = size_old+1;
    ////////////////////////////////
    {   // NDへの値の設定
        std::vector< std::vector< MatVec::CMat_BlkCrs* > > old_matrix_nd;
	    old_matrix_nd = m_Matrix_NonDia;
	    m_Matrix_NonDia.resize(0);
	    m_Matrix_NonDia.resize( size_new );
	    for(unsigned int ils=0;ils<size_new;ils++){
		    m_Matrix_NonDia[ils].resize( size_new, 0 );
	    }
	    for(unsigned int ils=0;ils<size_old;ils++){
		    for(unsigned int jls=0;jls<size_old;jls++){
			    m_Matrix_NonDia[ils][jls] = old_matrix_nd[ils][jls];
		    }
	    }
    }
    // その他のサイズ変更
	m_Matrix_Dia.resize( size_new, 0 );
	m_Residual.resize( size_new, 0 );
	m_Update.resize( size_new, 0 );
	m_BCFlag.resize( size_new, 0 );
    m_Residual[ size_old ] = new MatVec::CVector_Blk( nnode, len );
    m_Update[ size_old ] = new MatVec::CVector_Blk( nnode, len );
    m_BCFlag[ size_old ] = new MatVec::CBCFlag( nnode, len );
	m_aSeg.push_back( CLinSysSeg(nnode,len) );
	return size_old;
}


int LsSol::CLinearSystem::AddLinSysSeg( unsigned int nnode, const std::vector<unsigned int>& aLen )
{
	const unsigned int size_old = m_aSeg.size();
	const unsigned int size_new = size_old+1;

    ////////////////////////////////
    {   // NonDiaMatrixへの値の代入
        std::vector< std::vector< MatVec::CMat_BlkCrs* > > old_matrix_nd;
	    old_matrix_nd = m_Matrix_NonDia;
	    m_Matrix_NonDia.resize(0);
	    m_Matrix_NonDia.resize( size_new );
	    for(unsigned int ils=0;ils<size_new;ils++){
		    m_Matrix_NonDia[ils].resize( size_new, 0 );
	    }
	    for(unsigned int ils=0;ils<size_old;ils++){
		for(unsigned int jls=0;jls<size_old;jls++){
		    m_Matrix_NonDia[ils][jls] = old_matrix_nd[ils][jls];
		}
	    }
    }
    ////////////////////////////////
    // その他のサイズの変更
	m_Matrix_Dia.resize( size_new, 0 );
	m_Residual.resize(   size_new, 0 );
	m_Update.resize(     size_new, 0 );
	m_BCFlag.resize(     size_new, 0 );
    m_Residual[ size_old ] = new MatVec::CVector_Blk( nnode, aLen );
    m_Update[   size_old ] = new MatVec::CVector_Blk( nnode, aLen );
    m_BCFlag[   size_old ] = new MatVec::CBCFlag(     nnode, aLen );
	m_aSeg.push_back( CLinSysSeg(nnode,aLen) );
	return size_old;
}


////////////////////////////////////////////////////////////////
// ソルバ用のユーティリティ関数
////////////////////////////////////////////////////////////////

bool LsSol::CLinearSystem::ReSizeTmpVecSolver(unsigned int ntmp_new)
{
	const unsigned int nseg = m_aSeg.size();
	const unsigned int ntmp_old = this->GetTmpVectorArySize();

	if( ntmp_old == ntmp_new ){ return true; }
	else if( ntmp_old < ntmp_new ){
		m_TmpVectorArray.resize(ntmp_new);
		for(unsigned int ivec=ntmp_old;ivec<ntmp_new;ivec++){
			m_TmpVectorArray[ivec].resize(nseg);
			for(unsigned int iseg=0;iseg<nseg;iseg++){
                const unsigned int nnode = m_aSeg[iseg].nnode;
				const int len = m_aSeg[iseg].len;
                if( len == -1 ){
                    std::vector<unsigned int>& aLen = m_aSeg[iseg].aLen;
                    assert( aLen.size() == nnode );
                    m_TmpVectorArray[ivec][iseg] = new MatVec::CVector_Blk(nnode,aLen);
                }
                else{
                    m_TmpVectorArray[ivec][iseg] = new MatVec::CVector_Blk(nnode,len);
				    m_TmpVectorArray[ivec][iseg]->SetVectorZero();
                }
			}
		}
	}
	else{
		assert( ntmp_old > ntmp_new );
		for(unsigned int ivec=ntmp_new-1;ivec>=ntmp_old;ivec--){
			for(unsigned int iseg=0;iseg<nseg;iseg++){
				delete m_TmpVectorArray[ivec][iseg];
			}
			m_TmpVectorArray[ivec].clear();
		}
		m_TmpVectorArray.resize(ntmp_old);
	}
	return true;
}

////////////////////////////////
// 行列ベクトル積
// {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
// v1 != v2     v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
bool LsSol::CLinearSystem::MATVEC(double alpha, int iv1, double beta, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

    std::vector< MatVec::CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

    std::vector< MatVec::CVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	if( alpha == 0.0 ){
		this->SCAL(beta,iv2);
		return true;
	}

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	if( iv1 == iv2 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		if( this->m_Matrix_Dia[iseg] != 0 ){
			m_Matrix_Dia[iseg]->MatVec( alpha, (*(*p_vec1)[iseg]), beta, (*(*p_vec2)[iseg]) );
		}
		else{ (*(*p_vec2)[iseg]) *= beta; }
		for(unsigned int jseg=0;jseg<nseg;jseg++){
			if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
			assert( iseg != jseg );
			m_Matrix_NonDia[iseg][jseg]->MatVec( alpha, (*(*p_vec1)[jseg]), 1.0, (*(*p_vec2)[iseg]), true );
		}
	}
	return true;
}

////////////////////////////////
// ベクトル同士の足し算
// {v2} := alpha*{v1} + {v2}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
bool LsSol::CLinearSystem::AXPY(double alpha, int iv1, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;
	if( alpha == 0 ) return true;

    std::vector< MatVec::CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

    std::vector< MatVec::CVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	if( iv1 == iv2 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec2)[iseg]).AXPY( alpha, (*(*p_vec1)[iseg]) );
	}

	return true;
}

////////////////////////////////
// 内積を求める関数
// return {v1} * {v2}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
double LsSol::CLinearSystem::DOT(int iv1, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

    std::vector< MatVec::CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

    std::vector< MatVec::CVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	double norm = 0.0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		norm += (*(*p_vec1)[iseg]) * (*(*p_vec2)[iseg]);
	}
	return norm;
}

////////////////////////////////
// ベクトルのコピー
// return {v2} := {v1}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル 
bool LsSol::CLinearSystem::COPY(int iv1, int iv2){

	if( iv1 == iv2 ) return true;

	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

    std::vector< MatVec::CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

    std::vector< MatVec::CVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec2)[iseg]) = (*(*p_vec1)[iseg]);
	}
	return true;
}

////////////////////////////////
// ベクトルのスカラー倍
// return {v1} := alpha * {v1}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル 
bool LsSol::CLinearSystem::SCAL(double alpha,int iv1)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

    std::vector< MatVec::CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else{ assert(0); }

	assert( p_vec1->size() == nseg );

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec1)[iseg]) *= alpha;
	}
	return true;
}



bool LsSol::CLinearSystem::AddMat_Dia(unsigned int ils, const Com::CIndexedArray& crs)
{
	assert( ils < m_aSeg.size() );
	assert( crs.CheckValid() );
    
	{
		const unsigned int nnode = m_aSeg[ils].nnode;
		const int len = m_aSeg[ils].len;
    const std::vector<unsigned int>& aLen = m_aSeg[ils].aLen;
		if( m_Matrix_Dia[ils] == 0 ){
      if( len == -1 ){
        assert( aLen.size() == nnode );
        m_Matrix_Dia[ils] = new MatVec::CMatDia_BlkCrs();
        m_Matrix_Dia[ils]->Initialize(nnode,aLen);
      }
      else{
        assert( len > 0 );
        m_Matrix_Dia[ils] = new MatVec::CMatDia_BlkCrs(nnode,len);
      }
		}
    else{
			assert( m_Matrix_Dia[ils]->NBlkMatCol() ==  nnode );
			assert( m_Matrix_Dia[ils]->LenBlkCol()  ==  len );
      assert( m_Matrix_Dia[ils]->NBlkMatRow() ==  nnode );
			assert( m_Matrix_Dia[ils]->LenBlkRow()  ==  len );
    }
  }
	m_Matrix_Dia[ils]->AddPattern(crs);
	return true;
}

bool LsSol::CLinearSystem::AddMat_NonDia(
    unsigned int ils_col, 
    unsigned int ils_row, 
    const Com::CIndexedArray& crs )
{
	assert( ils_col < m_aSeg.size() );
	assert( ils_row < m_aSeg.size() );
	{
		const unsigned int nblk_col = m_aSeg[ils_col].nnode;
		const unsigned int nblk_row = m_aSeg[ils_row].nnode;
		const int len_col = m_aSeg[ils_col].len;
		const int len_row = m_aSeg[ils_row].len;
		if( m_Matrix_NonDia[ils_col][ils_row] == 0 ){
            if( len_col >= 0 && len_row >= 0 ){
                m_Matrix_NonDia[ils_col][ils_row] = new MatVec::CMat_BlkCrs(nblk_col,len_col, nblk_row,len_row);
            }
            else if( len_col == -1 && len_row >=  0 ){
                const std::vector<unsigned int>& aLen_col = m_aSeg[ils_col].aLen;
                std::vector<unsigned int> aLen_row;
                aLen_row.resize( nblk_row );
                for(unsigned int iblk=0;iblk<nblk_row;iblk++){
                    aLen_row[iblk] = len_row;
                }
                m_Matrix_NonDia[ils_col][ils_row] = new MatVec::CMat_BlkCrs(nblk_col,aLen_col, nblk_row,aLen_row);
            }
            else if( len_col >= 0  && len_row == -1 ){
                const std::vector<unsigned int>& aLen_row = m_aSeg[ils_row].aLen;
                std::vector<unsigned int> aLen_col;
                aLen_col.resize( nblk_col );
                for(unsigned int iblk=0;iblk<nblk_col;iblk++){
                    aLen_col[iblk] = len_col;
                }
                m_Matrix_NonDia[ils_col][ils_row] = new MatVec::CMat_BlkCrs(nblk_col,aLen_col, nblk_row,aLen_row);
            }
            else if( len_col == -1 && len_row == -1 ){
                const std::vector<unsigned int>& aLen_col = m_aSeg[ils_col].aLen;
                const std::vector<unsigned int>& aLen_row = m_aSeg[ils_row].aLen;
                m_Matrix_NonDia[ils_col][ils_row] = new MatVec::CMat_BlkCrs(nblk_col,aLen_col, nblk_row,aLen_row);
            }
		}
		else{
			assert( m_Matrix_NonDia[ils_col][ils_row]->NBlkMatCol() ==  nblk_col );
			assert( m_Matrix_NonDia[ils_col][ils_row]->LenBlkCol()  ==  len_col  );
			assert( m_Matrix_NonDia[ils_col][ils_row]->NBlkMatRow() ==  nblk_row );
			assert( m_Matrix_NonDia[ils_col][ils_row]->LenBlkRow()  ==  len_row  );
		}
	}
	assert( crs.CheckValid() );
	m_Matrix_NonDia[ils_col][ils_row]->AddPattern(crs);
	return true;
}
