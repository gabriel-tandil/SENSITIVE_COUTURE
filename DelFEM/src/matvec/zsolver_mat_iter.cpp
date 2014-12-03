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

#include <iostream>

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/zsolver_mat_iter.h"
#include "delfem/matvec/zmatprecond_blk.h"

using namespace MatVec;

////////////////////////////////////////////////////////////////
// Solve Matrix with CG Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_CG(double& conv_ratio,
				unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec){

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
    assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	u_vec.SetVectorZero();

	double sq_norm_res0;
	double sq_inv_norm_res_ini;
	{
		sq_norm_res0 = r_vec.GetSquaredVectorNorm();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res0;
	}

	CZVector_Blk  p_vec(nblk,blk_len);
	CZVector_Blk Ap_vec(nblk,blk_len);

	////////////////////////////////
	// Set Initial Serch Direction 
	// {p} = {r}
	////////////////////////////////
	p_vec = r_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		mat.MatVec(1.0,p_vec,0.0,Ap_vec);
		const double val_pAp = InnerProduct(p_vec,Ap_vec).Real();
		const double alpha = sq_norm_res0 / val_pAp;

		////////////////////////////////
		// Update update and residual
		// {u} = {u} + alpha * {p}
		// {r} = {r} - alpha * [A]{p}
		////////////////////////////////
		u_vec.AXPY( alpha, p_vec);
		r_vec.AXPY(-alpha,Ap_vec);

		////////////////////////////////
		// Converge Judgement
		////////////////////////////////
		double sq_norm_res1;
		{
			sq_norm_res1 = r_vec.GetSquaredVectorNorm();
			std::cout << iitr << " " << sqrt(sq_norm_res1 * sq_inv_norm_res_ini) << std::endl;
			if( sq_norm_res1 * sq_inv_norm_res_ini < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res1 * sq_inv_norm_res_ini );
				iteration = iitr;
				return true;
			}
		}
		////////////////////////////////
		// Calc Beta
		// beta = - (r,Ap) / (p,Ap)
		////////////////////////////////
		const double beta = sq_norm_res1 / sq_norm_res0;
		sq_norm_res0 = sq_norm_res1;

		////////////////////////////////
		// Update Search Direction
		// {p} = {r} + beta*{p}
		////////////////////////////////
		p_vec *= beta;
		p_vec.AXPY(1.0,r_vec);
	}

	return true;
}



////////////////////////////////////////////////////////////////
// Solve Matrix with Preconditioned-Conjugate Gradient Methods	
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_PCG(double& conv_ratio,
				unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec,
				const CZMatPrecond_Blk& precond, const CZMatDia_BlkCrs& mat_p){

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nnode = mat.NBlkMatCol();
	const unsigned int nlen = mat.LenBlkCol();

	u_vec.SetVectorZero();	
	double sq_inv_norm_res;
	{
		double dtmp1 = r_vec.GetSquaredVectorNorm();
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		std::cout << "Initial Norm Res " << sqrt(dtmp1) << std::endl;
		sq_inv_norm_res = 1.0 / dtmp1;
	}

	CZVector_Blk  p_vec(nnode,nlen);
	CZVector_Blk Ap_vec(nnode,nlen);
	CZVector_Blk  z_vec(nnode,nlen);

	z_vec = r_vec;
	
	precond.SolvePrecond(mat_p,z_vec);

	p_vec = z_vec;

	double val_rz0 = InnerProduct(r_vec, z_vec).Real();

	for(unsigned int iitr=0;iitr<mx_iter;iitr++){

		////////////////////////////////
		// Calc Alpha
		////////////////////////////////
		double alpha;
		{
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			double dtmp0 = InnerProduct(p_vec,Ap_vec).Real();
			alpha = val_rz0 / dtmp0;
		}

		////////////////////////////////
		// Update Solution
		////////////////////////////////
		u_vec.AXPY( alpha, p_vec);

		////////////////////////////////
		// Update residual
		////////////////////////////////
		r_vec.AXPY(-alpha,Ap_vec);

		////////////////////////////////
		// Converge Judgement
		////////////////////////////////
		{
			double sq_norm_res = r_vec.GetSquaredVectorNorm();
			std::cout << iitr << " " << sqrt(sq_norm_res * sq_inv_norm_res) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_norm_res * sq_inv_norm_res < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}

		////////////////////////////////
		// Update Z vector
		////////////////////////////////
		z_vec = r_vec;
		precond.SolvePrecond(mat_p,z_vec);

		////////////////////////////////
		// Calc Beta, RZ
		////////////////////////////////
		double beta;
		{
			double val_rz1 = InnerProduct(r_vec, z_vec).Real();
			beta = val_rz1 / val_rz0;
			val_rz0 = val_rz1;
		}

		////////////////////////////////
		// Update Direction
		////////////////////////////////
		p_vec *= beta;
		p_vec += z_vec;
	}
	return false;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_COCG(double& conv_ratio,
				unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec){

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	u_vec.SetVectorZero();

	double sq_inv_norm_res_ini;
	{
		const double sq_norm_res0 = r_vec.GetSquaredVectorNorm();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res0;
	}

	CZVector_Blk  p_vec(nblk,blk_len);
	CZVector_Blk Ap_vec(nblk,blk_len);

	// {p} = {r}
	p_vec = r_vec;

	Com::Complex res_res0 = r_vec*r_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		Com::Complex alpha;
		{
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			const Com::Complex val_pAp = Ap_vec*p_vec;
//			std::cout << val_pAp.Imag() << std::endl;
			alpha = res_res0 / val_pAp;
//			std::cout << " alpha " << alpha.Real() << " " << alpha.Imag() << std::endl;
		}

		// {r} = {r} - alpha * [A]{p}
		r_vec.AXPY(-alpha,Ap_vec);

		{	// Converge Judgement
			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double ratio = sqrt(sq_norm_res*sq_inv_norm_res_ini);
			std::cout << iitr << " " << ratio << " " << Com::Norm(res_res0) << std::endl;
			if( ratio < conv_ratio_tol ){
				conv_ratio = ratio;
				iteration = iitr;
				return true;
			}
		}

		// {u} = {u} + alpha * {p}
		u_vec.AXPY( alpha, p_vec);

		Com::Complex beta;
		{
			const Com::Complex res_res1 = r_vec*r_vec;
			beta = res_res1/res_res0;
			res_res0 = res_res1;
		}

		// {p} = {r} + beta*{p}
		p_vec *= beta;
		p_vec += r_vec;
	}

	return true;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_PCOCG(double& conv_ratio,
				unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec,
				const CZMatPrecond_Blk& precond, const CZMatDia_BlkCrs& mat_p){

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	u_vec.SetVectorZero();

	double sq_inv_norm_res_ini;
	{
		const double sq_norm_res0 = r_vec.GetSquaredVectorNorm();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res0;
	}

	CZVector_Blk  p_vec(nblk,blk_len);
	CZVector_Blk Ap_vec(nblk,blk_len);
	CZVector_Blk  w_vec(nblk,blk_len);

	w_vec = r_vec;
	precond.SolvePrecond(mat_p,w_vec);

	// Set Initial Serch Direction 
	// {p} = {r}
	p_vec = w_vec;

	Com::Complex val_rw = r_vec*w_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		Com::Complex alpha;
		{
			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
			const Com::Complex val_pAp = p_vec*Ap_vec;
			alpha = val_rw / val_pAp;
		}

		// {r} = {r} - alpha * [A]{p}
		r_vec.AXPY(-alpha,Ap_vec);

		{	// Converge Judgement
			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double ratio = sqrt(sq_norm_res*sq_inv_norm_res_ini);
			std::cout << iitr << " " << ratio << " " << std::endl;
			if( ratio < conv_ratio_tol ){
				conv_ratio = ratio;
				iteration = iitr;
				return true;
			}
		}

		// {u} = {u} + alpha * {p}
		u_vec.AXPY( alpha, p_vec);

		w_vec = r_vec;
		precond.SolvePrecond(mat_p,w_vec);

		// Calc Beta
		Com::Complex beta;
		{
			Com::Complex val_rw2 = r_vec*w_vec;
			beta = val_rw2/val_rw;
			val_rw = val_rw2;
		}

		// {p} = {r} + beta*{p}
		p_vec *= beta;
		p_vec += w_vec;
	}

	return true;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_CGNR(double& conv_ratio, unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec){

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	u_vec.SetVectorZero();

	double sq_inv_norm_res_ini;
	{
		const double sq_norm_res0 = r_vec.GetSquaredVectorNorm();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res0;
	}

	CZVector_Blk  p_vec(nblk,blk_len);
	CZVector_Blk  z_vec(nblk,blk_len);
	CZVector_Blk Ap_vec(nblk,blk_len);

	mat.MatVec_Hermitian(1.0,r_vec,0.0,z_vec);
	double sq_norm_z0 = z_vec.GetSquaredVectorNorm();

	////////////////////////////////
	// Set Initial Serch Direction 
	// {p} = {r}
	////////////////////////////////
	p_vec = z_vec;

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		mat.MatVec(1.0,p_vec,0.0,Ap_vec);
		const double sq_norm_Ap = Ap_vec.GetSquaredVectorNorm();
	
		const double alpha = sq_norm_z0 / sq_norm_Ap;

		////////////////////////////////
		// Update update and residual
		// {u} = {u} + alpha * {p}
		// {r} = {r} - alpha * [A]{p}
		////////////////////////////////
		u_vec.AXPY( alpha, p_vec);
		r_vec.AXPY(-alpha,Ap_vec);

		////////////////////////////////
		// Converge Judgement
		////////////////////////////////
		{
			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double ratio = sqrt(sq_norm_res*sq_inv_norm_res_ini);
			std::cout << iitr << " " << ratio << " " << Com::Norm(sq_norm_z0) << std::endl;
			if( ratio < conv_ratio_tol ){
				conv_ratio = ratio;
				iteration = iitr;
				return true;
			}
		}

		////////////////////////////////
		// Calc Beta
		// beta = - (r,Ap) / (p,Ap)
		////////////////////////////////
		mat.MatVec_Hermitian(1.0,r_vec,0.0,z_vec);
		double sq_norm_z1 = z_vec.GetSquaredVectorNorm();
		const double beta = sq_norm_z1 / sq_norm_z0;
		sq_norm_z0 = sq_norm_z1;

		////////////////////////////////
		// Update Search Direction
		// {p} = {r} + beta*{p}
		////////////////////////////////
		p_vec *= beta;
		p_vec += z_vec;
	}

	return true;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool MatVec::Sol::Solve_BiCGSTAB(double& conv_ratio, unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, 
				CZVector_Blk& r_vec, CZVector_Blk& u_vec){

//	std::cout.precision(18);

	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	CZVector_Blk  s_vec(nblk,blk_len);
	CZVector_Blk As_vec(nblk,blk_len);
	CZVector_Blk  p_vec(nblk,blk_len);
	CZVector_Blk Ap_vec(nblk,blk_len);

	CZVector_Blk r0_conjugate_vec(nblk,blk_len);

	u_vec.SetVectorZero();
/*
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		u_vec.SetValue(iblk,0,(iblk+1)*3.14);
	}
	mat.MatVec(-1.0,u_vec,1.0,r_vec);
*/
/*
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		const Complex& c = r_vec.GetValue(iblk,0);
		std::cout << iblk << " " << c.Real() << " " << c.Imag() << std::endl;
	}
*/

	double sq_inv_norm_res;
	{
		double dtmp1 = r_vec.GetSquaredVectorNorm();
		std::cout << "Initial Residual: " << sqrt(dtmp1) << " " << dtmp1 << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

	r0_conjugate_vec = r_vec;
//	r_vec.SetVectorConjugate();

	// {p} = {r}
	p_vec = r_vec;

	// calc (r,r0*)
	Com::Complex r_r0conj = InnerProduct(r_vec,r0_conjugate_vec);

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// calc {Ap} = [A]*{p}
		mat.MatVec(1.0,p_vec,0.0,Ap_vec);

		// calc alpha
		Com::Complex alpha;
		{
			const Com::Complex denominator = InnerProduct(Ap_vec,r0_conjugate_vec);
			alpha = r_r0conj / denominator;
		}

//		std::cout << "Alpha    " << alpha.Real() << " " << alpha.Imag() << std::endl;

		// calc s_vector
		s_vec = r_vec;
		s_vec.AXPY(-alpha,Ap_vec);

		// calc {As} = [A]*{s}
		mat.MatVec(1.0,s_vec,0.0,As_vec);

		// calc omega
		Com::Complex omega;
		{
			const double denominator = As_vec.GetSquaredVectorNorm();
			const Com::Complex numerator = InnerProduct(s_vec,As_vec);
			omega = numerator / denominator;
		}

//		std::cout << "Omega    " << omega.Real() << " " << omega.Imag() << std::endl;

		// update solution
		u_vec.AXPY(alpha,p_vec);
		u_vec.AXPY(omega,s_vec);

		// update residual
		r_vec = s_vec;
		r_vec.AXPY(-omega,As_vec);

		{
			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}
//		getchar();

		// calc beta
		Com::Complex beta;
		{
			const Com::Complex tmp1 = InnerProduct(r_vec,r0_conjugate_vec);
//			beta = (alpha*tmp1)/(omega*r_r0conj);
			beta = (tmp1/r_r0conj) * (alpha/omega);
//			beta = alpha*tmp1/omega/r_r0conj;
			r_r0conj = tmp1;
		}

//		std::cout << "Beta    " << beta.Real() << " " << beta.Imag() << std::endl;

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
bool MatVec::Sol::Solve_BiCGStabP(double& conv_ratio, unsigned int& iteration,
				const CZMatDia_BlkCrs& mat, CZVector_Blk& r_vec, CZVector_Blk& u_vec,
				const CZMatPrecond_Blk& precond, CZMatDia_BlkCrs& mat_p)
{
	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	const unsigned int nblk = mat.NBlkMatCol();
	assert( r_vec.BlkVecLen() == nblk );
	assert( u_vec.BlkVecLen() == nblk );

	const unsigned int blk_len = mat.LenBlkCol();
	assert( r_vec.BlkLen() == blk_len );
	assert( u_vec.BlkLen() == blk_len );

	CZVector_Blk   s_vec(nblk,blk_len);
	CZVector_Blk  Ms_vec(nblk,blk_len);
	CZVector_Blk AMs_vec(nblk,blk_len);
	CZVector_Blk   p_vec(nblk,blk_len);
	CZVector_Blk  Mp_vec(nblk,blk_len);
	CZVector_Blk AMp_vec(nblk,blk_len);

	CZVector_Blk r0_conjugate_vec(nblk,blk_len);

	u_vec.SetVectorZero();

	double sq_inv_norm_res;
	{
		double dtmp1 = r_vec.GetSquaredVectorNorm();
		std::cout << "Initial Residual : " << sqrt(dtmp1) << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

	r0_conjugate_vec = r_vec;
	r_vec.SetVectorConjugate();

	// {p} = {r}
	p_vec = r_vec;

	// calc (r,r0*)
	Com::Complex r_r0conj = InnerProduct(r_vec,r0_conjugate_vec);

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// {Mp_vec} = [M^-1]*{p}
		Mp_vec = p_vec;
		precond.SolvePrecond(mat_p,Mp_vec);

		// calc {AMp_vec} = [A]*{Mp_vec}
		mat.MatVec(1.0,Mp_vec,0.0,AMp_vec);

		// calc alpha
		Com::Complex alpha;
		{
			const Com::Complex denominator = InnerProduct(AMp_vec,r0_conjugate_vec);
			alpha = r_r0conj / denominator;
		}

		// calc s_vector
		s_vec = r_vec;
		s_vec.AXPY(-alpha,AMp_vec);

		// {Ms_vec} = [M^-1]*{s}
		Ms_vec = s_vec;
		precond.SolvePrecond(mat_p,Ms_vec);

		// calc {AMs} = [A]*{Ms_vec}
		mat.MatVec(1.0,Ms_vec,0.0,AMs_vec);

		// calc omega
		Com::Complex omega;
		{
			const double denominator = AMs_vec.GetSquaredVectorNorm();
			const Com::Complex numerator = InnerProduct(s_vec,AMs_vec);
			omega = numerator / denominator;
		}

		// update solution
		u_vec.AXPY(alpha,Mp_vec);
		u_vec.AXPY(omega,Ms_vec);

		// update residual
		r_vec = s_vec;
		r_vec.AXPY(-omega,AMs_vec);

		{
			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				iteration = iitr;
				return true;
			}
		}

		// calc beta
		Com::Complex beta;
		{
			const Com::Complex tmp1 = InnerProduct(r_vec,r0_conjugate_vec);
			beta = (alpha/omega) * (tmp1/r_r0conj);
			r_r0conj = tmp1;
		}

		// update p_vector
		p_vec *= beta;
		p_vec += r_vec;
		p_vec.AXPY(-beta*omega,AMp_vec);
	}

	return true;
}
