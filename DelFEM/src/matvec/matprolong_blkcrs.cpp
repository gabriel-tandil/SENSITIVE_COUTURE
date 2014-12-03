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
// Prolongation_Crs.cpp: CMatProlong_BlkCrs クラスのインプリメンテーション
////////////////////////////////////////////////////////////////

#ifndef for 
#define for if(0); else for
#endif

#include <iostream>
#include <cassert>
#include <math.h>
#include <vector>

#include "delfem/matvec/matprolong_blkcrs.h"
#include "delfem/matvec/matdia_blkcrs.h"

using namespace MatVec;


//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CMatProlong_BlkCrs::CMatProlong_BlkCrs(const CMatDia_BlkCrs& mat,const double& theta){
	const unsigned int nblk = mat.NBlkMatCol();

	////////////////
	// Make Crs_S
	////////////////
	unsigned int  ncrs_s = 0;
	unsigned int* crs_s = 0;
	unsigned int* crs_ind_s = 0;
	{
		ncrs_s = 0;
		crs_s = new unsigned int [mat.NCrs()];
		crs_ind_s = new unsigned int [nblk+1];
		crs_ind_s[0] = 0;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int ncrs_i;
			const unsigned int* crs_i = mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i = 0;
			const double* val_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );
			double val_max = 0.0;
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
//				if( val_crs_i[icrs_i]<0.0 ){
//					const double dtmp1 = -val_crs_i[icrs_i];
					const double dtmp1 = fabs(val_crs_i[icrs_i]);
					val_max = (dtmp1>val_max) ? dtmp1 : val_max;
//				}
			}
			const double* ptr_val_dia_i = mat.GetPtrValDia(iblk);
			const double val_dia = ptr_val_dia_i[0];
			if( val_max > 1.0e-20*val_dia ){
				for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
					const unsigned int jblk0 = crs_i[icrs_i];
//					if( -val_crs_i[icrs_i] >= theta*val_max+1.0e-20 ){ 
					if( fabs(val_crs_i[icrs_i]) >= theta*val_max+1.0e-20 ){ 
						crs_s[ncrs_s] = jblk0;
						ncrs_s++;
					}
				}
			}
			crs_ind_s[iblk+1] = ncrs_s;
		}
	}

	////////////////
	// Make Crs_ST
	////////////////
	unsigned int  ncrs_st = 0;
	unsigned int* crs_st = 0;
	unsigned int* crs_ind_st = 0;
	{
		crs_ind_st = new unsigned int [nblk+1];
		for(unsigned int iblk=0;iblk<nblk+1;iblk++){
			crs_ind_st[iblk] = 0;
		}
		for(unsigned int jblk=0;jblk<nblk;jblk++){
			for(unsigned int jicrs=crs_ind_s[jblk];jicrs<crs_ind_s[jblk+1];jicrs++){
				assert( jicrs < ncrs_s );
				const unsigned int iblk0 = crs_s[jicrs];
				assert( iblk0 < nblk );
				crs_ind_st[iblk0+1]++;
			}
		}
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			crs_ind_st[iblk+1] += crs_ind_st[iblk];
		}
		ncrs_st = crs_ind_st[nblk];
		assert( ncrs_st == ncrs_s );
		crs_st = new unsigned int [ncrs_st];
		for(unsigned int jblk=0;jblk<nblk;jblk++){
			for(unsigned int jicrs=crs_ind_s[jblk];jicrs<crs_ind_s[jblk+1];jicrs++){
				assert( jicrs < ncrs_s );
				const unsigned int iblk0 = crs_s[jicrs];
				assert( iblk0 < nblk );
				const unsigned int ijcrs = crs_ind_st[iblk0];
				assert( ijcrs < ncrs_st );
				crs_st[ijcrs] = jblk;
				crs_ind_st[iblk0]++;
			}
		}
		assert( crs_ind_st[nblk] == ncrs_st );
		for(unsigned int iblk=nblk;iblk>0;iblk--){
			crs_ind_st[iblk] = crs_ind_st[iblk-1];
		}
		crs_ind_st[0] = 0;
	}

	// Initialize lambda 
	int* lambda_flag = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		lambda_flag[iblk] = crs_ind_st[iblk+1]-crs_ind_st[iblk];
//		std::cout << "initial " << iblk << " Lv:" << lambda_flag[iblk] << std::endl;
	}
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] == 0 ){	// No neighberhood
			lambda_flag[iblk] = -3;	// set Fine point
			for(unsigned int icrs_s=crs_ind_s[iblk];icrs_s<crs_ind_s[iblk+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				unsigned int jblk0 = crs_s[icrs_s];
				assert( jblk0 < nblk );
				if( lambda_flag[jblk0] > 0 ){
					lambda_flag[jblk0]++;
				}
			}
		}
	}

	// make max_lambda
	unsigned int max_lambda=0;
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		max_lambda = ( lambda_flag[iblk] > (int)max_lambda ) ? lambda_flag[iblk] : max_lambda;
	}
//	std::cout << max_lambda << std::endl;

	// Initialize Bucket
	std::vector<SBucket>  aBucket;
	aBucket.reserve(40);
	aBucket.resize(max_lambda+1);
	for(unsigned int ilambda=0;ilambda<max_lambda+1;ilambda++){
		aBucket[ilambda].data = new unsigned int [nblk];
		aBucket[ilambda].size = 0;
	}
	unsigned int* blk_location = new unsigned int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		const int ilambda0 = lambda_flag[iblk];
		if( ilambda0 > 0 ){
			assert( (unsigned int)ilambda0 < aBucket.size() );
			assert( ilambda0 != 0 );
			const unsigned int isize0 = aBucket[ilambda0].size;
			assert( isize0 <= nblk  );
			blk_location[iblk] = isize0;
			aBucket[ilambda0].data[isize0] = iblk;
			aBucket[ilambda0].size++;
		}
		else{
			assert( ilambda0 == -3 );
		}
	}

	// F/C Splitting
	for(;;){
		// pick up 'iblk0' from 'max_lambda' bucket
		unsigned int iblk0;
		{
			if( aBucket[max_lambda].size == 0 ){
				for(int ilambda=max_lambda-1;ilambda>=0;ilambda--){
					max_lambda = ilambda;
					if( aBucket[ilambda].size != 0 ){ break; }
				}
				if( max_lambda == 0 && aBucket[0].size == 0 ) break;
			}
//			std::cout << max_lambda << " " << aBucket[max_lambda].size << std::endl;
			assert( aBucket[max_lambda].size > 0 );
			const unsigned int isize0 = aBucket[max_lambda].size;
			iblk0 = aBucket[max_lambda].data[ isize0-1 ];
            assert( lambda_flag[iblk0] == (int)max_lambda );
			assert( blk_location[iblk0] == isize0-1 );
			aBucket[max_lambda].size--;
			lambda_flag[iblk0] = -1;	// C Point
		}

//		std::cout << " Pick Up iblk0 : " << iblk0 << " Lv: " << max_lambda << std::endl;

		for(unsigned int icrs_st=crs_ind_st[iblk0];icrs_st<crs_ind_st[iblk0+1];icrs_st++){
			assert( icrs_st < ncrs_st );
			unsigned int jblk0 = crs_st[icrs_st];
			assert( jblk0 < nblk );
			if( lambda_flag[jblk0] < 0 ) continue;	// C/F unknown
			{	// make 'jblk0' st_around 'iblk0' to F 
				unsigned int ilambda0 = lambda_flag[jblk0];
				assert( ilambda0 < aBucket.size() );
				unsigned int iloc0 = blk_location[jblk0];
				assert( aBucket[ilambda0].data[iloc0] == jblk0 );
				lambda_flag[jblk0] = -3;	// F Point
				if( iloc0 == aBucket[ilambda0].size-1 ){
					aBucket[ilambda0].size--;
				}
				else{
					unsigned int iloc_last = aBucket[ilambda0].size-1;
					unsigned int iblk_last = aBucket[ilambda0].data[iloc_last];
					assert( blk_location[iblk_last] == iloc_last );
                    assert( lambda_flag[iblk_last] == (int)ilambda0 );
					blk_location[iblk_last] = iloc0;
					aBucket[ilambda0].data[iloc0] = iblk_last;
					aBucket[ilambda0].size--;
				}
			}

			// update lambda of 'kblk0' s_around 'jblk0'
			for(unsigned int icrs_s=crs_ind_s[jblk0];icrs_s<crs_ind_s[jblk0+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				unsigned int kblk0 = crs_s[icrs_s];
				assert( kblk0 < nblk );
				if( lambda_flag[kblk0] >= 0 ){
					unsigned int ilambda0 = lambda_flag[kblk0];
					assert( ilambda0 < aBucket.size() );
					{	// delete kblk from ilambda0
						unsigned int iloc0 = blk_location[kblk0];
						assert( aBucket[ilambda0].data[iloc0] == kblk0 );
						if( iloc0 == aBucket[ilambda0].size-1 ){
							aBucket[ilambda0].size--;
						}
						else{
							unsigned int iloc_last = aBucket[ilambda0].size-1;
							unsigned int iblk_last = aBucket[ilambda0].data[iloc_last];
							assert( blk_location[iblk_last] == iloc_last );
                            assert( lambda_flag[iblk_last] == (int)ilambda0 );
							blk_location[iblk_last] = iloc0;
							aBucket[ilambda0].data[iloc0] = iblk_last;
							aBucket[ilambda0].size--;
						}
					}
					{	// add kblk to ilambda0+1
						if(  aBucket.size() <= ilambda0+1 ){
							aBucket.resize( ilambda0+1+1 );
							aBucket[ilambda0+1].size = 0;
							aBucket[ilambda0+1].data = new unsigned int [nblk];
						}
						if( max_lambda < ilambda0+1 ){
							max_lambda = ilambda0+1;
						}
						assert( max_lambda >= ilambda0+1 );
						unsigned int iloc_last = aBucket[ilambda0+1].size;
						assert( iloc_last < nblk );
						aBucket[ilambda0+1].data[iloc_last] = kblk0;
						aBucket[ilambda0+1].size++;
						lambda_flag[kblk0] = ilambda0+1;
						blk_location[kblk0] = iloc_last;
					}
				}
			} // end loop for (icrs_s)
		} // end loop for (icrs_st)
	}

	delete[] blk_location;
	for(unsigned int ilambda=0;ilambda<aBucket.size();ilambda++){
		delete[] aBucket[ilambda].data;
	}

	////////////////////////////////////////
	// Corrct C/F Splitting by Second Pass
	////////////////////////////////////////
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] == -1 ) continue;
		unsigned int iflag=0;
		for(;;){
			for(unsigned int icrs_s=crs_ind_s[iblk];icrs_s<crs_ind_s[iblk+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				const unsigned int jblk0 = crs_s[icrs_s];
				assert( jblk0 < nblk );
				if( lambda_flag[jblk0] == -3 ){
					lambda_flag[jblk0] = -4;	// Mark Si
				}
				if( lambda_flag[jblk0] == -1 ){
					lambda_flag[jblk0] = -2;	// Mark Ci
				}
			}
			unsigned int ncrs_i;
			const unsigned int* crs_i = mat.GetPtrIndPSuP(iblk,ncrs_i);
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int lblk0 = crs_i[icrs_i];
				assert( lblk0 < nblk );
				if( lambda_flag[lblk0] != -2 ) continue;	
				// lbkl0 is Ci
				for(unsigned int icrs_st=crs_ind_st[lblk0];icrs_st<crs_ind_st[lblk0+1];icrs_st++){
					assert( icrs_st < ncrs_st );
					const unsigned int kblk0 = crs_st[icrs_st];
					assert( kblk0 < nblk );
					if( lambda_flag[kblk0] == -4 ){
						lambda_flag[kblk0] = -3;
					}
				}
			}
			for(unsigned int icrs_s=crs_ind_s[iblk];icrs_s<crs_ind_s[iblk+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				const unsigned int jblk0 = crs_s[icrs_s];
				assert( jblk0 < nblk );
				if( lambda_flag[jblk0] == -2 ){ // jblk0 is Point Ci
					lambda_flag[jblk0] = -1;
				}
				else if( lambda_flag[jblk0] == -4 ){	// jblk0 is Poit Si
					if( iflag == 0 ){
						lambda_flag[jblk0] = -1;
						iflag = 1;
						continue;
					}
					else if( iflag == 2 ){
						lambda_flag[iblk] = -1;
						lambda_flag[jblk0] = -3;
						iflag = 3;
						continue;
					}
					lambda_flag[jblk0] = -3;
				}
			}
			if( iflag != 1 ) break;
			iflag = 2;
		}
	}
	
	delete[] crs_st;
	delete[] crs_ind_st;

	unsigned int nblk_c=0;
	unsigned int nblk_f=0;
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] == -3 ){ 	// fine point
			nblk_f++; 
			continue; 
		}	
		else if( lambda_flag[iblk] == -1 ){ // coase point
			lambda_flag[iblk] = nblk_c;
			nblk_c++; 
			continue; 
		}
		else{ 
			std::cout << "Unknown : " << iblk << " " << lambda_flag[iblk] << std::endl; 
			assert(0);
		}
	}

//	std::cout << " map " << nblk << " to " << nblk_c << std::endl;
	assert( nblk_c + nblk_f == nblk );

	////////////////////////////////
	// Make Prolongation Matrix Pattern
	////////////////////////////////
	this->Initialize(nblk,1,nblk_c,1);
	assert( m_colInd_Blk != 0 );
	for(unsigned int iblk=0;iblk<nblk+1;iblk++){
		m_colInd_Blk[iblk] = 0;
	}
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] < 0 ){ // iblk is F point 
			assert( lambda_flag[iblk] == -3 );
			for(unsigned int icrs_s=crs_ind_s[iblk];icrs_s<crs_ind_s[iblk+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				unsigned int jblk0 = crs_s[icrs_s];
				assert( jblk0 < nblk );
				if( lambda_flag[jblk0] >= 0 ){ // jblk0 is Ci
					m_colInd_Blk[iblk+1]++;
				}
			}
		}
		else{	// iblk is C point
			m_colInd_Blk[iblk+1]++;
		}
	}
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		m_colInd_Blk[iblk+1] += m_colInd_Blk[iblk];
	}
	m_ncrs_Blk = m_colInd_Blk[nblk];
	m_rowPtr_Blk = new unsigned int [m_ncrs_Blk];
//	std::cout << " Prolongation Crs Size : " << m_ncrs_Blk << std::endl;
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] < 0 ){	// F point
			assert( lambda_flag[iblk] == -3 );
			for(unsigned int icrs_s=crs_ind_s[iblk];icrs_s<crs_ind_s[iblk+1];icrs_s++){
				assert( icrs_s < ncrs_s );
				unsigned int jblk0 = crs_s[icrs_s];
				assert( jblk0 < nblk );
				if( lambda_flag[jblk0] >= 0 ){
					unsigned int icrs0 = m_colInd_Blk[iblk];
					assert( icrs0 < m_ncrs_Blk );
					m_rowPtr_Blk[icrs0] = lambda_flag[jblk0];
					m_colInd_Blk[iblk]++;
				}
			}
		}
		else{	// C Point
			assert( (unsigned int)lambda_flag[iblk] < nblk_c );
			unsigned int icrs0 = m_colInd_Blk[iblk];
			assert( icrs0 < m_ncrs_Blk );
			m_rowPtr_Blk[icrs0] = lambda_flag[iblk];
			m_colInd_Blk[iblk]++;
		}
	}
	for(unsigned int iblk=nblk;iblk>0;iblk--){
		m_colInd_Blk[iblk] = m_colInd_Blk[iblk-1];
	}
	m_colInd_Blk[0] = 0;
	assert( m_ncrs_Blk == m_colInd_Blk[nblk] );

	delete[] crs_s;
	delete[] crs_ind_s;

	////////////////////////////////
	// Make C/F Mapping
	////////////////////////////////
	this->m_MapperCF = new unsigned int [nblk_c];
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( lambda_flag[iblk] >= 0 ){
			unsigned int iblk_c = lambda_flag[iblk];
			m_MapperCF[iblk_c] = iblk;
//			std::cout << " coarse " << iblk << std::endl;
		}
		else{
			assert( lambda_flag[iblk] == -3 );
//			std::cout << " fine   " << iblk << std::endl;
		}
	}

	delete[] lambda_flag;

	////////////////
	// Fill Value
	////////////////
	this->MakeValueAMG(mat,theta);
}

CMatProlong_BlkCrs::~CMatProlong_BlkCrs()
{

}

bool CMatProlong_BlkCrs::MakeValueAMG(const CMatDia_BlkCrs& mat,const double& theta)
{
	assert( mat.NBlkMatCol() == this->NBlkMatCol() );
	assert( mat.LenBlkCol() == 1 );

	if( m_valCrs_Blk == 0 ){ 
		m_valCrs_Blk = new double [NCrs() * LenBlkCol()*LenBlkRow()];
	}

	const unsigned int nblk = NBlkMatCol();
	const unsigned int nblk_c = NBlkMatRow();

	int* row_type = 0;
	{
		row_type = new int [nblk];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			row_type[iblk] = -1; 
		}
	}

	int* row2crs = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ row2crs[iblk] = -1; }

	double* nume_k = new double [nblk];

	for(unsigned int iblk=0;iblk<nblk;iblk++){

		// iblk is Originally Coarse point
		if( m_colInd_Blk[iblk+1]-m_colInd_Blk[iblk] == 1 ){
			const unsigned int icrs0 = m_colInd_Blk[iblk];
			assert( icrs0 < m_ncrs_Blk );
			const unsigned int jblk_c = m_rowPtr_Blk[icrs0];
			assert( jblk_c < nblk_c );
			const unsigned int jblk0 = m_MapperCF[jblk_c];
			assert( jblk0 < nblk );
			if( jblk0 == iblk ){
				m_valCrs_Blk[icrs0] = 1.0;
				continue;
			}
		}

		// iblk is No Cuppling Point
		if( m_colInd_Blk[iblk+1] == m_colInd_Blk[iblk] ) continue;
		
		unsigned int ncrs_i;
		const unsigned int* crs_i = mat.GetPtrIndPSuP(iblk,ncrs_i);

		unsigned int tmp_ncrs_i = 0;
		const double* val_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
		assert( tmp_ncrs_i == ncrs_i );

		{	// Set Row-type Strong/Week Cuppled to iblk
			
			double val_max_row_i = 0.0;
			if( ncrs_i > 0 ){
				for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
//					if(  val_crs_i[icrs_i] < 0.0 ){
//						const double dtmp1 = -val_crs_i[icrs_i];
						const double dtmp1 = fabs(val_crs_i[icrs_i]);
						val_max_row_i = ( dtmp1 > val_max_row_i ) ? dtmp1 : val_max_row_i;
//					}
				}
			}
			else{ val_max_row_i = 0.0; }

			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk = crs_i[icrs_i];
				assert( jblk < nblk  );
				assert( row_type[jblk] == -1 );
//				assert( -val_crs_i[icrs_i] <= val_max_row_i+1.0e-20 );
//				if( -val_crs_i[icrs_i] >= theta*val_max_row_i+1.0e-20 ){ row_type[jblk] = -2; } // Strong Cuppled
				if( fabs(val_crs_i[icrs_i]) >= theta*val_max_row_i+1.0e-20 ){ row_type[jblk] = -2; } // Strong Cuppled
				else{ row_type[jblk] = -3; } // Week Cuppled
			}
		}

		// Set Row-type of Ci point to index
		const unsigned int ncoarse = m_colInd_Blk[iblk+1]-m_colInd_Blk[iblk];
		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			const unsigned int jblk_c = m_rowPtr_Blk[icrs];
			assert( jblk_c < nblk_c );
			const unsigned int jblk0 = m_MapperCF[jblk_c];
			assert( jblk0 < nblk );
			assert( jblk0 != iblk );
			assert( row_type[jblk0] == -2 );
			row_type[jblk0] = icrs-m_colInd_Blk[iblk];
			assert( row_type[jblk0] >= 0 );
		}

		{	// set a_ik to nume_k
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk = crs_i[icrs_i];
				assert( jblk < nblk );
				row2crs[jblk] = icrs_i;
			}
			for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
				const unsigned int jblk_c = m_rowPtr_Blk[icrs];
				assert( jblk_c < nblk_c );
				const unsigned int jblk0 = m_MapperCF[jblk_c];
				assert( jblk0 < nblk );
				assert( jblk0 != iblk );
                assert( row_type[jblk0] == (int)icrs-(int)m_colInd_Blk[iblk] );

				const unsigned int jcrs_i0 = row2crs[jblk0];
				assert( jcrs_i0 < ncrs_i );
				assert( crs_i[jcrs_i0] == jblk0 );
				const double val_ik = val_crs_i[jcrs_i0];

				const unsigned int kcoarse = icrs-m_colInd_Blk[iblk];
				assert( kcoarse < ncoarse );
				nume_k[kcoarse] = val_ik;
			}
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk = crs_i[icrs_i];
				assert( jblk < nblk );
				row2crs[jblk]  = -1;
			}
		}	// end set a_ik to nume_k

		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){	// loop for Si
			const unsigned int jblk = crs_i[icrs_i];
			assert( jblk  < nblk );
			if( row_type[jblk] != -2 ) continue;
			// jblk is Si and not Ci
			const double val_ij = val_crs_i[icrs_i];

			unsigned int ncrs_j;
			const unsigned int* crs_j = mat.GetPtrIndPSuP(jblk,ncrs_j); 

			unsigned int tmp_ncrs_j = 0;
			const double* val_crs_j = mat.GetPtrValPSuP(jblk,tmp_ncrs_j);
			assert( tmp_ncrs_j == ncrs_j );

			double inv_a_jl = 0.0;
			{	// add val1
				double a_jl = 0.0;
				int iflag = 0;
				for(unsigned int jcrs_j=0;jcrs_j<ncrs_j;jcrs_j++){
					unsigned int kblk = crs_j[jcrs_j];
					assert( kblk < nblk );
					if( row_type[kblk] < 0  ) continue;
					// kblk is Ci point
					const double val_jk = val_crs_j[jcrs_j];
					a_jl += val_jk;
					iflag = 1;
				}
				{	// calc inv_a_jl
					if( iflag == 0){ 
						std::cout<< "uncuppled fc " << iblk << " " << jblk << std::endl;
						inv_a_jl = 0.0; 
					}
					else {
						if( fabs(a_jl) > 1.0e-20 ){ 
							assert( iflag == 1 );
							inv_a_jl = 1.0 / a_jl; 
						}
						else{ 
							std::cout << iflag << std::endl;
							std::cout << " uncuppled " << iblk << " " << jblk << std::endl;
							inv_a_jl = 0.0; 
						}
					}
				}
			}
			// add nume_k
			double dtmp1 = 0.0;
			for(unsigned int jcrs_j=0;jcrs_j<ncrs_j;jcrs_j++){
				unsigned int kblk = crs_j[jcrs_j];
				assert( kblk < nblk );
				if( row_type[kblk] < 0  ) continue;
				// kblk is Ci point
				const double val_jk = val_crs_j[jcrs_j];
				const unsigned int kcoarse0 = row_type[kblk];
				assert( kcoarse0 < ncoarse );
				nume_k[kcoarse0] += (val_ij*val_jk*inv_a_jl);
				dtmp1 += val_jk*inv_a_jl;
			}
			assert( fabs(dtmp1-1.0) < 1.0e-5 );
		}

		double val_deno;
		{
			// make val2
			double val_dwi = 0.0;
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				assert( row_type[jblk0] != -1 );
				if( row_type[jblk0] != -3 ) continue; 
				// jblk is Dwi point
				val_dwi += val_crs_i[icrs_i];
			}
			double val_ii;
			{
				const double* val_dia = mat.GetPtrValDia(iblk);
				val_ii = val_dia[0];
			}
			val_deno = val_ii+val_dwi;
		}

		for(unsigned int icrs=m_colInd_Blk[iblk];icrs<m_colInd_Blk[iblk+1];icrs++){
			const unsigned int kcoarse = icrs-m_colInd_Blk[iblk];
			assert( kcoarse < ncoarse );
			const double val = -nume_k[kcoarse]/val_deno;
			m_valCrs_Blk[icrs] = val;
		}
		
		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
			const unsigned int jblk0 = crs_i[icrs_i];
			row_type[jblk0] = -1;
		}
	}

	delete[] nume_k;
	delete[] row_type;
	delete[] row2crs;

	return true;
}

bool CMatProlong_BlkCrs::ScatterCoaseValue(unsigned int* distination, const unsigned int* coarse_value)
{
	const unsigned int nblk = NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		distination[iblk] = 0;
	}
	const unsigned int nblk_c = NBlkMatRow();
	for(unsigned int iblk_c=0;iblk_c<nblk_c;iblk_c++){
		const unsigned int ival0 = coarse_value[iblk_c];
		const unsigned int idest0 = m_MapperCF[iblk_c];
		distination[idest0] = ival0;
	}
	return true;
}


bool CMatProlong_BlkCrs::MaskCoarse(unsigned int* distination)
{
	const unsigned int nblk = NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		distination[iblk] = 0;
	}
	const unsigned int nblk_c = NBlkMatRow();
	for(unsigned int iblk_c=0;iblk_c<nblk_c;iblk_c++){
		const unsigned int idest0 = m_MapperCF[iblk_c];
		distination[idest0] = 1;
	}
	return true;
}
