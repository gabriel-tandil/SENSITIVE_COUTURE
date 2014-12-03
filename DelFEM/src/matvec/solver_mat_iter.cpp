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

#include <iostream>
#include <math.h>

#include "delfem/matvec/solver_mat_iter.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/matprecond_blk.h"

using namespace MatVec;

////////////////////////////////////////////////////////////////
// Solve Matrix with CG Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_CG
(double& conv_ratio,
 unsigned int& iteration,
 const MatVec::CMatDia_BlkCrs& mat, MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec)
{
	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.NBlk() == nblk );
	assert( u_vec.NBlk() == nblk );

  assert( mat.LenBlkCol() >= 0 );
  const unsigned int blk_len = mat.LenBlkCol();
  assert( r_vec.Len() == (int)blk_len );
  assert( u_vec.Len() == (int)blk_len );

 	u_vec.SetVectorZero();

	double sq_norm_res;
	double sq_inv_norm_res;
	{
		sq_norm_res = r_vec*r_vec;
		std::cout << sqrt(sq_norm_res) << std::endl;
		if( sq_norm_res < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / sq_norm_res;
	}

	CVector_Blk  p_vec(nblk,blk_len);
	CVector_Blk Ap_vec(nblk,blk_len);

	// Set Initial Serch Direction 
	// {p} = {r}
	p_vec = r_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		double alpha;
		{	// alpha = (r,r) / (p,Ap)
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			const double pAp = p_vec*Ap_vec;
			alpha = sq_norm_res / pAp;
		}

		// {r} = {r} - alpha * [A]{p}
		r_vec.AXPY(-alpha,Ap_vec);

		double sq_norm_res1 = r_vec*r_vec;
		{	// Converge Judgement
//			std::cout << iitr << " " << sqrt(sq_norm_res * sq_inv_norm_res) << "  " << sqrt(sq_norm_res) << std::endl;
			if( sq_norm_res * sq_inv_norm_res < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}

		// {u} = {u} + alpha * {p}
		u_vec.AXPY(alpha,p_vec);

		// beta = (r1,r1) / (r0,r0)
		const double beta = sq_norm_res1 / sq_norm_res;
		sq_norm_res = sq_norm_res1;


		// {p} = {r} + beta*{p}
		p_vec *= beta;
		p_vec += r_vec;
	}

	return true;
}

bool MatVec::Sol::Solve_PCG
(double& conv_ratio,
 unsigned int& iteration,
 const CMatDia_BlkCrs& mat, CVector_Blk& r_vec, CVector_Blk& u_vec,
 const CPrecond_Blk& precond, const CMatDia_BlkCrs& mat_p)
{

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nnode = mat.NBlkMatCol();
	const unsigned int nlen = mat.LenBlkCol();

	u_vec.SetVectorZero();	
	
	double sq_inv_norm_res0;
	{
		const double sq_norm_res0 = r_vec*r_vec;
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res0 = 1.0 / sq_norm_res0;
	}

	CVector_Blk p_vec(nnode,nlen);
	CVector_Blk z_vec(nnode,nlen);

	z_vec = r_vec;
	precond.SolvePrecond(mat_p,z_vec);

	p_vec = z_vec;

	double product_rz;
	for(unsigned int iitr=0;iitr<mx_iter;iitr++){
		double alpha;
		{
			product_rz = r_vec * z_vec;
			mat.MatVec(1.0,p_vec,0.0,z_vec);
			double tmp2 = p_vec*z_vec;
			alpha = product_rz / tmp2;
		}

		r_vec.AXPY(-alpha,z_vec);
		u_vec.AXPY(alpha, p_vec);

		{	// Converge Judgement
			double sq_norm_res = r_vec*r_vec;
//			std::cout << iitr << " " << sqrt(sq_norm_res * sq_inv_norm_res0) << std::endl;
			if( sq_norm_res * sq_inv_norm_res0 < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res0 );
				iteration = iitr;
				return true;
			}
		}

		double beta;
		{
			z_vec = r_vec;
			precond.SolvePrecond(mat_p,z_vec);
			double tmp2 = r_vec*z_vec;
			beta = tmp2/product_rz;
		}

		p_vec *= beta;
		p_vec += z_vec;
	}
	return false;
}



////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_BiCGSTAB(double& conv_ratio,
				unsigned int& iteration,
                const MatVec::CMatDia_BlkCrs& mat, 
                MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec)
{

	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.NBlk() == nblk );
	assert( u_vec.NBlk() == nblk );

    assert( mat.LenBlkCol() >= 0 );
	const unsigned int blk_len = mat.LenBlkCol();
    assert( r_vec.Len() == (int)blk_len );
    assert( u_vec.Len() == (int)blk_len );

	CVector_Blk  s_vec(nblk,blk_len);
	CVector_Blk As_vec(nblk,blk_len);
	CVector_Blk  p_vec(nblk,blk_len);
	CVector_Blk Ap_vec(nblk,blk_len);

	CVector_Blk r0_conjugate_vec(nblk,blk_len);

	u_vec.SetVectorZero();

	double sq_inv_norm_res;
	{
		double dtmp1 = r_vec*r_vec;
		std::cout << sqrt(dtmp1) << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

	r0_conjugate_vec = r_vec;

	// {p} = {r}
	p_vec = r_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// calc (r,r0*)
		const double r_r0conj = r_vec * r0_conjugate_vec;

		// calc {Ap} = [A]*{p}
		mat.MatVec(1.0,p_vec,0.0,Ap_vec);

		// calc alpha
		double alpha;
		{
			const double denominator = Ap_vec * r0_conjugate_vec;
			alpha = r_r0conj / denominator;
		}

		// calc s_vector
		s_vec = r_vec;
		s_vec.AXPY(-alpha,Ap_vec);

		// calc {As} = [A]*{s}
		mat.MatVec(1.0,s_vec,0.0,As_vec);

		// calc omega
		double omega;
		{
			const double denominator = As_vec * As_vec;
			const double numerator = As_vec * s_vec;
			omega = numerator / denominator;
		}

		// update solution
		u_vec.AXPY(alpha,p_vec);
		u_vec.AXPY(omega,s_vec);

		// update residual
		r_vec = s_vec;
		r_vec.AXPY(-omega,As_vec);

		{
			const double sq_norm_res = r_vec*r_vec;
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}

		// calc beta
		double beta;
		{
			const double tmp1 = r_vec * r0_conjugate_vec;
			beta = tmp1 * alpha / (r_r0conj*omega);
		}

		// update p_vector
		p_vec *= beta;
		p_vec += r_vec;
		p_vec.AXPY(-beta*omega,Ap_vec);
	}

	return true;
}


////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_PBiCGSTAB(double& conv_ratio,
				unsigned int& iteration,
				const CMatDia_BlkCrs& mat, 
				CVector_Blk& r_vec, CVector_Blk& u_vec,
				const CPrecond_Blk& precond, const MatVec::CMatDia_BlkCrs& mat_p)
{
	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.NBlk() == nblk );
	assert( u_vec.NBlk() == nblk );

    assert( mat.LenBlkCol() >= 0 );
	const unsigned int blk_len = mat.LenBlkCol();
    assert( r_vec.Len() == (int)blk_len );
    assert( u_vec.Len() == (int)blk_len );

	CVector_Blk   s_vec(nblk,blk_len);
	CVector_Blk  Ms_vec(nblk,blk_len);
	CVector_Blk AMs_vec(nblk,blk_len);
	CVector_Blk   p_vec(nblk,blk_len);
	CVector_Blk  Mp_vec(nblk,blk_len);
	CVector_Blk AMp_vec(nblk,blk_len);

	CVector_Blk r0_conjugate_vec(nblk,blk_len);

	u_vec.SetVectorZero();

	double sq_inv_norm_res;
	{
		double dtmp1 = r_vec*r_vec;
		std::cout << sqrt(dtmp1) << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

	r0_conjugate_vec = r_vec;

	// {p} = {r}
	p_vec = r_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// {Mp_vec} = [M^-1]*{p}
		Mp_vec = p_vec;
		precond.SolvePrecond(mat,Mp_vec);

		// calc (r,r0*)
		const double r_r0conj = r_vec * r0_conjugate_vec;

		// calc {AMp_vec} = [A]*{Mp_vec}
		mat.MatVec(1.0,Mp_vec,0.0,AMp_vec);

		// calc alpha
		double alpha;
		{
			const double denominator = AMp_vec * r0_conjugate_vec;
			alpha = r_r0conj / denominator;
		}

		// calc s_vector
		s_vec = r_vec;
		s_vec.AXPY(-alpha,AMp_vec);

		// {Ms_vec} = [M^-1]*{s}
		Ms_vec = s_vec;
		precond.SolvePrecond(mat,Ms_vec);

		// calc {AMs_vec} = [A]*{Ms_vec}
		mat.MatVec(1.0,Ms_vec,0.0,AMs_vec);

		// calc omega
		double omega;
		{
			const double denominator = AMs_vec * AMs_vec;
			const double numerator = s_vec * AMs_vec;
			omega = numerator / denominator;
		}

		// update solution
		u_vec.AXPY(alpha,Mp_vec);
		u_vec.AXPY(omega,Ms_vec);

		// update residual
		r_vec = s_vec;
		r_vec.AXPY(-omega,AMs_vec);

		{
			const double sq_norm_res = r_vec*r_vec;
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}

		// calc beta
		double beta;
		{
			const double tmp1 = r_vec * r0_conjugate_vec;
			beta = tmp1 * alpha / (r_r0conj*omega);
		}

		// update p_vector
		p_vec *= beta;
		p_vec += r_vec;
		p_vec.AXPY(-beta*omega,AMp_vec);
	}

	return true;
}
