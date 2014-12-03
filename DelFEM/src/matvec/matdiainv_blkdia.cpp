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
// MatDiaCrsFrac.cpp: implementtion of ILU factorization class CMatDiaCrsFrac
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )   // C4786‚È‚ñ‚Ä•\Ž¦‚·‚ñ‚È( ß„Dß)ºÞÙ§
#endif
#define for if(0); else for

#include <math.h>
#include <assert.h>
#include <iostream>
#include <cstdlib> //(abort)

#include "delfem/matvec/matdiainv_blkdia.h"
#include "delfem/matvec/vector_blk.h"

using namespace MatVec;


CMatDiaInv_BlkDia::CMatDiaInv_BlkDia(unsigned int nblk, unsigned int blklen) 
:CMatDia_BlkCrs(nblk,blklen)
{
	m_valDia_Blk = 0;
}


CMatDiaInv_BlkDia::CMatDiaInv_BlkDia(const CMatDia_BlkCrs& mat) 
:CMatDia_BlkCrs(mat.NBlkMatCol(),mat.LenBlkCol())
{
	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Impliment" << std::endl;
		assert(0);
		abort();
	}
	////////////////
	for(unsigned int iblk=0;iblk<this->NBlkMatCol();iblk++){
		const double* ptr_dia_val =  mat.GetPtrValDia(iblk);
		const double dia_val = *ptr_dia_val;
		assert( fabs(dia_val) > 1.0e-30 );
		this->m_valDia_Blk[iblk] = 1.0/dia_val;
	}
}

CMatDiaInv_BlkDia::~CMatDiaInv_BlkDia()
{

}


int CMatDiaInv_BlkDia::SolveUpdate_GaussSidel(const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int ncrs_i;
			const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];
			double* pval_update_i = &update.m_Value[iblk];

			pval_update_i[0] = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_update_j = &update.m_Value[jblk0];
				pval_update_i[0] -= pval_crs_i[icrs_i]*pval_update_j[0];
			}
			pval_update_i[0] *= this->m_valDia_Blk[iblk];
		}
	}
	return 0;
}

int CMatDiaInv_BlkDia::SolveUpdate_GaussSidel_BackWard(   
		const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update ) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		for(int iblk=nblk-1;iblk>=0;iblk--){
			unsigned int ncrs_i;
			const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];
			double* pval_update_i = &update.m_Value[iblk];

			pval_update_i[0] = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_update_j = &update.m_Value[jblk0];
				pval_update_i[0] -= pval_crs_i[icrs_i]*pval_update_j[0];
			}
			pval_update_i[0] *= this->m_valDia_Blk[iblk];
		}
	}
	return 0;
	return true;
}

int CMatDiaInv_BlkDia::SolveUpdate_GaussSidel_Mask(  
		const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		unsigned int* mask, bool is_mask_later ) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			assert( mask[iblk]==0 || mask[iblk]==1 );
			if( (mask[iblk]==1)==is_mask_later ) continue;
			unsigned int ncrs_i;
			const unsigned int* crs_i = mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];
			double* pval_update_i = &update.m_Value[iblk];

			pval_update_i[0] = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_update_j = &update.m_Value[jblk0];
				pval_update_i[0] -= pval_crs_i[icrs_i]*pval_update_j[0];
			}
			pval_update_i[0] *= this->m_valDia_Blk[iblk];
		}
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			assert( mask[iblk]==0 || mask[iblk]==1 );
			if( (mask[iblk]==0)==is_mask_later ) continue;
			unsigned int ncrs_i;
			const unsigned int* crs_i = mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];
			double* pval_update_i = &update.m_Value[iblk];

			pval_update_i[0] = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_update_j = &update.m_Value[jblk0];
				pval_update_i[0] -= pval_crs_i[icrs_i]*pval_update_j[0];
			}
			pval_update_i[0] *= this->m_valDia_Blk[iblk];
		}
	}
	return 0;
}

int CMatDiaInv_BlkDia::SolveUpdate_SOR(const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		double omega) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int ncrs_i;
			const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];

			double tmp_val0 = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_update_j = &update.m_Value[jblk0];
				tmp_val0 -= pval_crs_i[icrs_i]*pval_update_j[0];
			}
			tmp_val0 *= this->m_valDia_Blk[iblk];

			double* pval_update_i = &update.m_Value[iblk];
			pval_update_i[0] += (tmp_val0-pval_update_i[0])*omega;
		}
	}
	return 0;
}

int CMatDiaInv_BlkDia::SolveUpdate_Jacobi(const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		CVector_Blk& tmp_vec) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		tmp_vec = update;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int ncrs_i;
			const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];
			double* pval_update_i = &update.m_Value[iblk];

			pval_update_i[0] = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_tmp_vec_j = &tmp_vec.m_Value[jblk0];
				pval_update_i[0] -= pval_crs_i[icrs_i]*pval_tmp_vec_j[0];
			}
			pval_update_i[0] *= this->m_valDia_Blk[iblk];
		}
	}
	return 0;
}


int CMatDiaInv_BlkDia::SolveUpdate_OmegaJacobi(const unsigned int max_iter, 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		CVector_Blk& tmp_vec, double omega) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iitr=0;iitr<max_iter;iitr++){
		tmp_vec = update;
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int ncrs_i;
			const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

			unsigned int tmp_ncrs_i;
			const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
			assert( tmp_ncrs_i == ncrs_i );

			const double* pval_res_i = &residual.m_Value[iblk];

			double tmp_val0 = pval_res_i[0];
			for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
				const unsigned int jblk0 = crs_i[icrs_i];
				assert( jblk0 < nblk );
				double* pval_tmp_vec_j = &tmp_vec.m_Value[jblk0];
				tmp_val0 -= pval_crs_i[icrs_i]*pval_tmp_vec_j[0];
			}
			tmp_val0 *= this->m_valDia_Blk[iblk];
			
			double* pval_update_i = &update.m_Value[iblk];
			pval_update_i[0] += (tmp_val0-pval_update_i[0])*omega;
		}
	}
	return 0;
}


int CMatDiaInv_BlkDia::Solve_IniGaussSidel( 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		unsigned int ncrs_i;
		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

		unsigned int tmp_ncrs_i;
		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];
		double tmp_val0 = pval_res_i[0];

		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
			const unsigned int jblk0 = crs_i[icrs_i];
			assert( jblk0 < nblk );
			if( jblk0 > iblk ) break;
			double* pval_update_j = &update.m_Value[jblk0];
			tmp_val0 -= pval_crs_i[icrs_i]*pval_update_j[0];
		}
		tmp_val0 *= this->m_valDia_Blk[iblk];

		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0;
	}
	return 0;
}

int CMatDiaInv_BlkDia::Solve_IniGaussSidel_Mask( 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		unsigned int* mask, bool is_mask_later ) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( (mask[iblk]==1)==is_mask_later ) continue;
		unsigned int ncrs_i;
		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

		unsigned int tmp_ncrs_i;
		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];
		double tmp_val0 = pval_res_i[0];

		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
			const unsigned int jblk0 = crs_i[icrs_i];
			assert( jblk0 < nblk );
			if( (mask[jblk0]==1)==is_mask_later ) continue;
			if( jblk0 > iblk ) break;
			double* pval_update_j = &update.m_Value[jblk0];
			tmp_val0 -= pval_crs_i[icrs_i]*pval_update_j[0];
		}
		tmp_val0 *= this->m_valDia_Blk[iblk];

		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0;
	}
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		if( (mask[iblk]==0)==is_mask_later ) continue;
		unsigned int ncrs_i;
		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

		unsigned int tmp_ncrs_i;
		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];
		double tmp_val0 = pval_res_i[0];

		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
			const unsigned int jblk0 = crs_i[icrs_i];
			assert( jblk0 < nblk );
			if( jblk0 > iblk ) break;
			double* pval_update_j = &update.m_Value[jblk0];
			tmp_val0 -= pval_crs_i[icrs_i]*pval_update_j[0];
		}
		tmp_val0 *= this->m_valDia_Blk[iblk];

		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0;
	}
	return 0;
}

int CMatDiaInv_BlkDia::Solve_IniSOR( 
	const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
	double omega) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		unsigned int ncrs_i;
		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

		unsigned int tmp_ncrs_i;
		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];

		double tmp_val0 = pval_res_i[0];
		for(unsigned int icrs_i=0;icrs_i<ncrs_i;icrs_i++){
			const unsigned int jblk0 = crs_i[icrs_i];
			assert( jblk0 < nblk );
			if( jblk0 > iblk ) break;
			double* pval_update_j = &update.m_Value[jblk0];
			tmp_val0 -= pval_crs_i[icrs_i]*pval_update_j[0];
		}
		tmp_val0 *= this->m_valDia_Blk[iblk];

		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0*omega;
	}
	return 0;
}

int CMatDiaInv_BlkDia::Solve_IniJacobi( 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update ) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
//        unsigned int ncrs_i;
//		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

//		unsigned int tmp_ncrs_i;
//		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
//		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];

		double tmp_val0 = pval_res_i[0];
		tmp_val0 *= this->m_valDia_Blk[iblk];
		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0;
	}
	return 0;
}

int CMatDiaInv_BlkDia::Solve_IniOmegaJacobi( 
		const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, 
		double omega ) const
{
	assert( mat.LenBlkCol() == this->LenBlkCol() );

	////////////////
	if( mat.LenBlkCol() != 1 ){
		std::cout << "Not Implimented" << std::endl;
		assert(0);
	}
	////////////////
	const unsigned int nblk = mat.NBlkMatCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
//        unsigned int ncrs_i;
//		const unsigned int* crs_i =mat.GetPtrIndPSuP(iblk,ncrs_i);

//        unsigned int tmp_ncrs_i;
//		const double* pval_crs_i = mat.GetPtrValPSuP(iblk,tmp_ncrs_i);
//		assert( tmp_ncrs_i == ncrs_i );

		const double* pval_res_i = &residual.m_Value[iblk];

		double tmp_val0 = pval_res_i[0];
		tmp_val0 *= this->m_valDia_Blk[iblk];
			
		double* pval_update_i = &update.m_Value[iblk];
		pval_update_i[0] = tmp_val0*omega;
	}
	return 0;
}
