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
#pragma warning( disable : 4786 )   // C4786Ç»ÇÒÇƒï\é¶Ç∑ÇÒÇ»( ﬂÑDﬂ)∫ﬁŸß
#endif
#define for if(0); else for

#include <cassert>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdio.h>

#include "delfem/matvec/zmatdiafrac_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"

using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// ç\íz/è¡ñ≈
//////////////////////////////////////////////////////////////////////

CZMatDiaFrac_BlkCrs::CZMatDiaFrac_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow)
	:CZMatDia_BlkCrs(nblk_colrow,len_colrow)
{
	m_ConditionFlag = 0;
	m_DiaInd = 0;
	m_pRowLev = 0;
}

CZMatDiaFrac_BlkCrs::CZMatDiaFrac_BlkCrs(const CZMatDia_BlkCrs& rhs)
	:CZMatDia_BlkCrs(rhs.NBlkMatCol(),rhs.LenBlkCol())
{
	m_ConditionFlag = 0;
	m_DiaInd = 0;
	m_pRowLev = 0;

	AddPattern(rhs,true);
	assert( m_nblk_MatRow == m_nblk_MatCol );
	assert( m_len_BlkRow == m_len_BlkCol );
	{	// Make m_DiaInd
		m_DiaInd = new unsigned int [m_nblk_MatCol];
		for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
			m_DiaInd[inode] = m_colInd_Blk[inode+1];
			for(unsigned int icrs=m_colInd_Blk[inode];icrs<m_colInd_Blk[inode+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jnode0 = m_rowPtr_Blk[icrs];
				assert( jnode0 < m_nblk_MatRow );
				if( jnode0 > inode ){
					m_DiaInd[inode] = icrs;
					break;
				}
			}
		}
	}
	m_valCrs_Blk = new Com::Complex [m_ncrs_Blk];
	m_ConditionFlag = 2;
}

CZMatDiaFrac_BlkCrs::CZMatDiaFrac_BlkCrs(const unsigned int lev_fill, const CZMatDia_BlkCrs& rhs)
:CZMatDia_BlkCrs(rhs.NBlkMatCol(),rhs.LenBlkCol())
{	
	m_ConditionFlag = 0;
	m_DiaInd = 0;
	m_pRowLev = 0;
	if( lev_fill == 0 ){
		AddPattern(rhs,true);
		{	// Make m_DiaInd
			m_DiaInd = new unsigned int [m_nblk_MatCol];
			for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
				m_DiaInd[inode] = m_colInd_Blk[inode+1];
				for(unsigned int icrs=m_colInd_Blk[inode];icrs<m_colInd_Blk[inode+1];icrs++){
					assert( icrs < m_ncrs_Blk );
					const unsigned int jnode0 = m_rowPtr_Blk[icrs];
					assert( jnode0 < m_nblk_MatRow );
					if( jnode0 > inode ){
						m_DiaInd[inode] = icrs;
						break;
					}
				}
			}
		}
		m_valCrs_Blk = new Com::Complex [m_ncrs_Blk];
		m_ConditionFlag = 2;
	}
	else{
		FracInitialize(rhs);
		AddFracPtn(lev_fill);
		FracFinalize();
	}
}

CZMatDiaFrac_BlkCrs::~CZMatDiaFrac_BlkCrs()
{
	if( m_DiaInd != 0 ){ delete[] m_DiaInd; }
	if( m_pRowLev != 0 ){ delete m_pRowLev; }
}

bool CZMatDiaFrac_BlkCrs::ForwardSubstitution(CZVector_Blk& vec) const
{
	assert( m_nblk_MatRow == m_nblk_MatCol );
	if( this->m_len_BlkCol == 1 ){
		for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
			Com::Complex lvec_i = vec.GetValue(inode,0);
			const unsigned int ijcrs0 = m_colInd_Blk[inode];
			const unsigned int ijcrs1 = m_DiaInd[inode];
			for(unsigned int ijcrs=ijcrs0;ijcrs<ijcrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
				assert( jnode0<inode );
				lvec_i -= m_valCrs_Blk[ijcrs]*vec.GetValue(jnode0,0);
			}
			vec.SetValue(inode,0,m_valDia_Blk[inode]*lvec_i);
		}
	}
	else{
		std::cout << "Error!-->Not Implimented!" << std::endl;
		assert(0);
	}
	return true;
}

bool CZMatDiaFrac_BlkCrs::BackwardSubstitution(CZVector_Blk& vec) const
{
	if( this->m_len_BlkCol == 1 ){
		for(int inode=m_nblk_MatCol-1;inode>=0;inode--){
            assert( inode < (int)m_nblk_MatCol );
			Com::Complex lvec_i = vec.GetValue(inode,0);
			const unsigned int ijcrs0 = m_DiaInd[inode];
			const unsigned int ijcrs1 = m_colInd_Blk[inode+1];
			for(unsigned int ijcrs=ijcrs0;ijcrs<ijcrs1;ijcrs++){
				assert( ijcrs<m_ncrs_Blk );
				const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
                assert( (int)jnode0>inode && (int)jnode0<(int)m_nblk_MatRow );
				lvec_i -=  m_valCrs_Blk[ijcrs]*vec.GetValue(jnode0,0);
			}
			vec.SetValue(inode,0,lvec_i);
		}
	}
	else{
		std::cout << "Error!-->Not Implimented!" << std::endl;
		assert(0);
	}
	return true;
}
/*
bool CZMatDiaFrac_BlkCrs::Solve(CZVector_Blk& vec) const 
{
	assert( m_nblk_MatRow == m_nblk_MatCol );

	for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
		Com::Complex lvec_i = vec.GetValue(inode,0);
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_DiaInd[inode];ijcrs++){
			assert( ijcrs<m_ncrs_Blk );
			const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
			assert( jnode0<inode );
			lvec_i -= m_valCrs_Blk[ijcrs]*vec.GetValue(jnode0,0);
		}
		vec.SetValue(inode,0,m_valDia_Blk[inode]*lvec_i);
	}
	for(int inode=m_nblk_MatCol-1;inode>=0;inode--){
		assert( inode < m_nblk_MatCol );
		Com::Complex lvec_i = vec.GetValue(inode,0);
		for(unsigned int ijcrs=m_DiaInd[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			assert( ijcrs<m_ncrs_Blk );
			const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
			assert( jnode0>inode && jnode0<m_nblk_MatRow );
			lvec_i -=  m_valCrs_Blk[ijcrs]*vec.GetValue(jnode0,0);
		}
		vec.SetValue(inode,0,lvec_i);
	}
	return true;
}
*/
bool CZMatDiaFrac_BlkCrs::DoILUDecomp()
{
	assert( m_nblk_MatRow == m_nblk_MatCol );

	
	if( this->m_len_BlkCol != 1 ){
		std::cout << "Error!-->Not Implimented!" << std::endl;
		assert(0);
	}

//	int info= 0;

	int* row2crs_f = new int [m_nblk_MatRow];
	for(unsigned int jnode=0;jnode<m_nblk_MatRow;jnode++){
		row2crs_f[jnode] = -1;
	}
	for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
		////////////////
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			assert( ijcrs>=0 && ijcrs<m_ncrs_Blk );
			const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
			assert( jnode0<m_nblk_MatRow );
			row2crs_f[jnode0] = ijcrs;
		}
		////////////////
		for(unsigned int ikcrs=m_colInd_Blk[inode];ikcrs<m_colInd_Blk[inode+1];ikcrs++){
			const unsigned int knode = m_rowPtr_Blk[ikcrs];
			assert( knode<m_nblk_MatCol );
			if( knode >= inode ) break;
			Com::Complex ikvalue = m_valCrs_Blk[ikcrs];
			for(unsigned int kjcrs=m_DiaInd[knode];kjcrs<m_colInd_Blk[knode+1];kjcrs++){
				const unsigned int jnode0 = m_rowPtr_Blk[kjcrs];
				assert( jnode0<m_nblk_MatRow );
				Com::Complex kjvalue = m_valCrs_Blk[kjcrs];
				if( jnode0 != inode ){
					const int ijcrs0 = row2crs_f[jnode0];
					if( ijcrs0 == -1 ) continue;
					m_valCrs_Blk[ijcrs0] -= ikvalue*kjvalue;
				}
				else{
					m_valDia_Blk[inode] -= ikvalue*kjvalue;
				}
			}
		}
		////////////////
		Com::Complex iivalue = m_valDia_Blk[inode];
		if( Com::SquaredNorm(iivalue) < 1.0e-20 ){
			std::cout << "frac false" << inode << std::endl;
			getchar();
		}
		m_valDia_Blk[inode] = 1.0 / iivalue;
		////////////////
		for(unsigned int ijcrs=m_DiaInd[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			assert( ijcrs>=0 && ijcrs<m_ncrs_Blk );
			m_valCrs_Blk[ijcrs] = m_valCrs_Blk[ijcrs] * m_valDia_Blk[inode];
		}
		////////////////
		for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			assert( ijcrs>=0 && ijcrs<m_ncrs_Blk );
			const unsigned int jnode0 = m_rowPtr_Blk[ijcrs];
			assert( jnode0<m_nblk_MatRow );
			row2crs_f[jnode0] = -1;
		}
		////////////////
	}	// end inode

	delete[] row2crs_f;

	return true;
}

bool CZMatDiaFrac_BlkCrs::AddFracPtn(const int lev_fill)
{
	assert( m_nblk_MatCol == m_nblk_MatRow );

	assert( m_ConditionFlag != 0 );
	if(  m_ConditionFlag == 0 ) return false;

//	std::cout << "	ADDFracPtn  1" << std::endl;

	if( m_ncrs_Blk == 0 ){
		std::cout << "	Dia Center " << std::endl;
		return true;
	}

	if( m_ncrs_Blk == m_nblk_MatCol * ( m_nblk_MatRow - 1 ) ){
		std::cout << "   Full Matrix " << std::endl;
		return true;
	}

	if( lev_fill == 0 ){
		std::cout << "    Fill Lev 0 " << std::endl;
		return true;
	}

	if( m_ConditionFlag == 2 ){

		if( m_valCrs_Blk    != 0 ){ delete[] m_valCrs_Blk;    }

		for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
			m_DiaInd[inode] = m_colInd_Blk[inode+1];
			for(unsigned int icrs=m_colInd_Blk[inode];icrs<m_colInd_Blk[inode+1];icrs++){
				assert( icrs < m_ncrs_Blk );
				const unsigned int jnode0 = m_rowPtr_Blk[icrs];
				assert( jnode0 < m_nblk_MatRow );
				if( jnode0 > inode ){
					m_DiaInd[inode] = icrs;
					break;
				}
			}
		}

		m_pRowLev = new std::vector<SRowLev>;
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
	const unsigned int* ColInd_pre = m_colInd_Blk;
	const unsigned int* DiaInd_pre = m_DiaInd;
	const std::vector<SRowLev>* pRowLev_pre = m_pRowLev;

	const unsigned int ncrs_pre = m_ncrs_Blk;

	m_colInd_Blk = new unsigned int [m_nblk_MatCol+1];
	m_colInd_Blk[0] = 0;
	m_DiaInd = new unsigned int [m_nblk_MatCol];

	m_pRowLev = new std::vector<SRowLev>;
	m_pRowLev->reserve(ncrs_pre*4);

	std::vector<SRowLevNext> nonzero;
	nonzero.reserve(m_nblk_MatCol);
	for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
		{	// copy row matrix pattern
			nonzero.resize(ColInd_pre[inode+1]-ColInd_pre[inode]);
			unsigned int inz = 0;
			for(unsigned int ijcrs=ColInd_pre[inode];ijcrs<ColInd_pre[inode+1];ijcrs++){
				assert( ijcrs<ncrs_pre );
				const unsigned int jnode0 = (*pRowLev_pre)[ijcrs].row;
				assert( jnode0<m_nblk_MatRow );
				const unsigned int ij_lev0 = (*pRowLev_pre)[ijcrs].lev;
				nonzero[inz].row = jnode0;
				nonzero[inz].lev = ij_lev0;
				nonzero[inz].next = inz+1;
				inz++;
			}
			nonzero[inz-1].next = -1;
		}
		int knz_cur = 0;
		for(;;){
			const unsigned int knode0 = nonzero[knz_cur].row;
			assert( knode0<m_nblk_MatCol );
			const int ik_lev0 = nonzero[knz_cur].lev;
			if( ik_lev0+1>lev_fill && lev_fill!=-1 ){
				knz_cur = nonzero[knz_cur].next;
				if( knz_cur == -1 ) break;
				continue;
			}
			if( knode0 >= inode ) break;

			unsigned int jnz_cur = knz_cur;
			for(unsigned int kjcrs=m_DiaInd[knode0];kjcrs<m_colInd_Blk[knode0+1];kjcrs++){
				const int kj_lev0 = (*m_pRowLev)[kjcrs].lev;
				if( kj_lev0+1>lev_fill && lev_fill!=-1 ) continue;
				const int unsigned jnode0 = (*m_pRowLev)[kjcrs].row;
				assert( jnode0>knode0 && jnode0<m_nblk_MatRow );
				assert( nonzero[jnz_cur].row < jnode0 );
				if( jnode0 == inode ) continue;

				// check if this is fill in
				bool is_fill_in = false;
				for(;;){
					const int jnz_nex = nonzero[jnz_cur].next;
                    assert( (jnz_nex>=0&&jnz_nex<(int)m_nblk_MatRow) || jnz_nex==-1 );
					if( jnz_nex == -1 ){ is_fill_in = true; break; }
					if( nonzero[jnz_nex].row  > jnode0 ){ is_fill_in = true; break; }
					if( nonzero[jnz_nex].row == jnode0 ){ break; }
					assert( nonzero[jnz_nex].row < jnode0 );
					jnz_cur = jnz_nex;
				}
				if( !is_fill_in ){ continue; }

				// pick up fill in
				const unsigned int max_lev0 = ( ik_lev0 > kj_lev0 ) ? ik_lev0 : kj_lev0;
				const unsigned  int inz_last = nonzero.size();
				nonzero.resize( nonzero.size()+1 );
				nonzero[inz_last].row = jnode0;
				nonzero[inz_last].lev = max_lev0 + 1;
				nonzero[inz_last].next = nonzero[jnz_cur].next;
				nonzero[jnz_cur].next = inz_last;
				jnz_cur = inz_last;
			}
			knz_cur = nonzero[knz_cur].next;
            assert( (knz_cur>=0&&knz_cur<(int)m_nblk_MatRow) || knz_cur==-1 );
			if( knz_cur == -1 ) break;
		}
		////////////////
		{
			if( ColInd_pre[inode] + nonzero.size() > m_pRowLev->max_size() ){
				std::cout << "		overflow and memory reallocate in ilu frac" << std::endl;
			}
			m_pRowLev->resize( m_colInd_Blk[inode] + nonzero.size() );
			int icrs0 = m_colInd_Blk[inode];
			for(int inz=0;inz!=-1;inz=nonzero[inz].next){
				const unsigned int jnode = nonzero[inz].row;
				const unsigned int jlev = nonzero[inz].lev;
				assert( jnode<m_nblk_MatRow );
				assert( jnode != inode );
				(*m_pRowLev)[icrs0].row = jnode;
				(*m_pRowLev)[icrs0].lev = jlev;
				icrs0++;
			}

			m_colInd_Blk[inode+1] = icrs0;
			m_ncrs_Blk = icrs0;
			m_DiaInd[inode] = icrs0;
			for(unsigned  int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
				const unsigned int jnode0 = (*m_pRowLev)[ijcrs].row;
				if( jnode0 > inode ){
					m_DiaInd[inode] = ijcrs;
					break;
				}
			}
		}
		////////////////
	}	// end inode

	delete[] const_cast<unsigned int*>(ColInd_pre);
	delete[] const_cast<unsigned int*>(DiaInd_pre);
	delete const_cast< std::vector<SRowLev>* >(pRowLev_pre);
	
	return true;
}

bool CZMatDiaFrac_BlkCrs::FracFinalize(){
	assert( m_nblk_MatRow == m_nblk_MatCol );

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

	assert( m_valCrs_Blk == 0 );
	m_valCrs_Blk = new Com::Complex [m_ncrs_Blk];

	m_ConditionFlag = 2;	// ILUï™âäÆóπ

	return true;
}


bool CZMatDiaFrac_BlkCrs::SetValue(const CZMatDia_BlkCrs& rhs)
{
	assert( m_nblk_MatCol == rhs.NBlkMatCol() );
	assert( m_nblk_MatRow == rhs.NBlkMatRow() );
	assert( m_nblk_MatRow == m_nblk_MatCol );

	if( m_ConditionFlag == 0 ){
		// make lev_fill 0 pattern
		CZMatDia_BlkCrs::AddPattern(rhs,true);
		m_valCrs_Blk = new Com::Complex [m_ncrs_Blk];
		m_ConditionFlag = 2;
	}
	else if( m_ConditionFlag == 1 ){
		FracFinalize();
	}
	CZMatDia_BlkCrs::SetValue(rhs,true);
	DoILUDecomp();
	return true;
}


bool CZMatDiaFrac_BlkCrs::FracInitialize(const CZMatDia_BlkCrs& rhs)
{

	CZMatDiaFrac_BlkCrs::AddPattern(rhs,true);

	assert( m_nblk_MatCol == rhs.NBlkMatCol() );
	assert( m_nblk_MatRow == rhs.NBlkMatRow() );
	assert( m_nblk_MatCol == m_nblk_MatRow );

	m_DiaInd = new unsigned int [m_nblk_MatCol];
	for(unsigned int inode=0;inode<m_nblk_MatCol;inode++){
		m_DiaInd[inode] = m_colInd_Blk[inode+1];
		for(unsigned int icrs=m_colInd_Blk[inode];icrs<m_colInd_Blk[inode+1];icrs++){
			assert( icrs < m_ncrs_Blk );
			const unsigned int jnode0 = m_rowPtr_Blk[icrs];
			assert( jnode0 < m_nblk_MatRow );
			if( jnode0 > inode ){
				m_DiaInd[inode] = icrs;
				break;
			}
		}
	}

	m_pRowLev = new std::vector<SRowLev>;
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
