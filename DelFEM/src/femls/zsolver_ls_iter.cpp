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

#include "delfem/femls/zsolver_ls_iter.h"
#include "delfem/femls/zpreconditioner.h"

//Å@Solve Hermetizn matrix "ls" with Conjugate Gradient Method
bool Fem::Ls::Solve_CG(double& conv_ratio, unsigned int& num_iter, 
			  Fem::Ls::CZLinearSystem& ls)
{	
	const unsigned int max_iter = num_iter;
	const double tolerance = conv_ratio;

	if( ls.GetTmpVectorArySize() < 2 ){ ls.ReSizeTmpVecSolver(2); }

	const int ix = -2;
	const int ir = -1;
	const int ip = 0;
	const int iAp = 1;

	// x = 0.0
	ls.SCAL(0.0,ix);

	double sq_norm_res;
	{
		sq_norm_res = ls.DOT(ir,ir).Real();
		if( sq_norm_res < 1.0e-60 ){
			num_iter = 0;
			conv_ratio = sqrt(sq_norm_res);
			return true;
		}
	}
	const double sq_inv_norm_res_ini = 1.0 / sq_norm_res;

	// p = r
	ls.COPY(ir,ip);

	for(num_iter=0;num_iter<max_iter;num_iter++){

		double alpha;
		{	// alpha = (r,r) / (p,Ap)
			ls.MatVec(1.0,ip,0.0,iAp);
			const double pAp = ls.DOT(ip,iAp).Real();
			alpha = sq_norm_res / pAp;
		}

		// update x 
		// x = x + alpha*p
		ls.AXPY(alpha,ip,ix);

		// update residual
		// r = r - alpha*Ap
		ls.AXPY(-alpha,iAp,ir);

		// calc residual 
		const double sq_norm_res_new = ls.DOT(ir,ir).Real();
		std::cout << num_iter << " " << sqrt( sq_norm_res_new*sq_inv_norm_res_ini ) << std::endl;
		if( sq_norm_res_new*sq_inv_norm_res_ini < tolerance*tolerance ){
			conv_ratio = sqrt( sq_norm_res_new*sq_inv_norm_res_ini );
			break;
		}

		double beta;
		{	// calc beta
			beta = sq_norm_res_new / sq_norm_res;
			sq_norm_res = sq_norm_res_new;
		}

		// update direction
		// {p} = {r}+beta*{p}
		ls.SCAL(beta,ip);
		ls.AXPY(1.0,ir,ip);
	}
	return true;
}





bool Fem::Ls::Solve_PCG(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond )
{

	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	if( ls.GetTmpVectorArySize() < 2 ){ ls.ReSizeTmpVecSolver(2); }

	const int ix  = -2;
	const int ir  = -1;
	const int ip  =  0;
	const int iz  =  1;

	// x = 0.0
	ls.SCAL(0.0,ix);
	
	double sq_inv_norm_res0;
	{
		const double sq_norm_res0 = ls.DOT(ir,ir).Real();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res0 = 1.0 / sq_norm_res0;
	}

	ls.COPY(ir,iz);
	precond.SolvePrecond(ls,iz);

	ls.COPY(iz,ip);

	double inpro_rz = ls.DOT(ir,iz).Real();
	for(unsigned int iitr=0;iitr<mx_iter;iitr++){

		double alpha;
		{	// calc alpha
			ls.MatVec(1.0,ip,0.0,iz);
			const double val_pAp = ls.DOT(ip,iz).Real();
			alpha = inpro_rz / val_pAp;
		}

		ls.AXPY(-alpha,iz,ir);
		ls.AXPY(alpha,ip,ix);	// Converge JudgementÇÃëOÇ…ì¸ÇÍÇÈ

		{	// Converge Judgement
			double sq_norm_res = ls.DOT(ir,ir).Real();
			std::cout << iitr << " " << sqrt(sq_norm_res * sq_inv_norm_res0) << std::endl;
			if( sq_norm_res * sq_inv_norm_res0 < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res0 );
				iteration = iitr;
				return true;
			}
		}

		double beta;
		{	// calc beta
			ls.COPY(ir,iz);
			precond.SolvePrecond(ls,iz);
			const double inpro_rz_new = ls.DOT(ir,iz).Real();
			beta = inpro_rz_new/inpro_rz;
			inpro_rz = inpro_rz_new;
		}

		ls.SCAL(beta,ip);
		ls.AXPY(1.0,iz,ip);
	}
	return false;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods
////////////////////////////////////////////////////////////////
bool Fem::Ls::Solve_PCOCG(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond )
{
	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = iteration;

	if( ls.GetTmpVectorArySize() < 3 ){ ls.ReSizeTmpVecSolver(3); }

	const int ix = -2;
	const int ir = -1;
	const int ip = 0;
	const int iAp = 1;
	const int iw = 2;

//	u_vec.SetVectorZero();
	ls.SCAL(0,ix);

	double sq_inv_norm_res_ini;
	{
//		const double sq_norm_res0 = r_vec.GetSquaredVectorNorm();
		const double sq_norm_res0 = ls.INPROCT(ir,ir).Real();
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res0;
	}

//	w_vec = r_vec;
//	precond.SolvePrecond(mat_p,w_vec);
	ls.COPY(ir,iw);
	precond.SolvePrecond(ls,iw);

	// Set Initial Serch Direction 
	// {p} = {r}
//	p_vec = w_vec;
	ls.COPY(iw,ip);

//	Com::Complex val_rw = r_vec*w_vec;
	Com::Complex val_rw = ls.DOT(ir,iw);

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++)
	{
		Com::Complex alpha;
		{
//			mat.MatVec(1.0,p_vec,0.0,Ap_vec);
//			const Com::Complex val_pAp = p_vec*Ap_vec;
			ls.MatVec(1,ip,0,iAp);
			const Com::Complex val_pAp = ls.DOT(ip,iAp);
			alpha = val_rw / val_pAp;
		}

		// {u} = {u} + alpha * {p}
//		u_vec.AXPY( alpha, p_vec);
		ls.AXPY(alpha,ip,ix);

		// {r} = {r} - alpha * [A]{p}
//		r_vec.AXPY(-alpha,Ap_vec);
		ls.AXPY(-alpha,iAp,ir);

		{	// Converge Judgement
//			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double sq_norm_res = ls.INPROCT(ir,ir).Real();
			const double ratio = sqrt(sq_norm_res*sq_inv_norm_res_ini);
//			std::cout << iitr << " " << ratio << " " << std::endl;
			if( ratio < conv_ratio_tol ){
				conv_ratio = ratio;
				iteration = iitr;
				return true;
			}
		}

//		w_vec = r_vec;
//		precond.SolvePrecond(mat_p,w_vec);
		ls.COPY(ir,iw);
		precond.SolvePrecond(ls,iw);

		// Calc Beta
		Com::Complex beta;
		{
//			Com::Complex val_rw2 = r_vec*w_vec;	
			Com::Complex val_rw2 = ls.DOT(ir,iw);
			beta = val_rw2/val_rw;
			val_rw = val_rw2;
		}

		// {p} = {r} + beta*{p}
//		p_vec *= beta;
//		p_vec += w_vec;
		ls.SCAL(beta,ip);
		ls.AXPY(1,iw,ip);
	}

	return true;
}


////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods
////////////////////////////////////////////////////////////////
bool Fem::Ls::Solve_CGNR(double& conv_ratio, unsigned int& num_iter,
				CZLinearSystem& ls)
{
	const double conv_ratio_tol = conv_ratio;
	const unsigned int mx_iter = num_iter;

	if( ls.GetTmpVectorArySize() < 3 ){ ls.ReSizeTmpVecSolver(3); }

	const int ix = -2;
	const int ir = -1;
	const int ip = 0;
	const int iAp = 1;
	const int iz = 2;

	// x = 0.0
	ls.SCAL(0.0,ix);

	double sq_norm_res;
	{
		sq_norm_res = ls.DOT(ir,ir).Real();
		if( sq_norm_res < 1.0e-60 ){
			num_iter = 0;
			conv_ratio = sqrt(sq_norm_res);
			return true;
		}
	}
	const double sq_inv_norm_res_ini = 1.0 / sq_norm_res;

	ls.MatVec_Hermitian(1,ir,0,iz);

	// p = z
	ls.COPY(iz,ip);

	double sq_norm_z0 = ls.INPROCT(iz,iz).Real();

	num_iter = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		////////////////////////////////
		// {Ap} = [A]{p}
		ls.MatVec(1.0,ip,0.0,iAp);
	
		////////////////////////////////
		// alpha = ({z}*{z})/({Ap}*{Ap})
		const double sq_norm_Ap = ls.INPROCT(iAp,iAp).Real();
		const double alpha = sq_norm_z0 / sq_norm_Ap;

		////////////////////////////////
		// Update update
		// {u} = {u} + alpha * {p}
		////////////////////////////////
		ls.AXPY(alpha,ip,ix);

		////////////////////////////////
		// Update residual
		// {r} = {r} - alpha * [A]{p}
		////////////////////////////////
		ls.AXPY(-alpha,iAp,ir);

		////////////////////////////////
		// Converge Judgement
		////////////////////////////////
		const double sq_norm_res_new = ls.INPROCT(ir,ir).Real();
		std::cout << iitr << " " << sqrt( sq_norm_res_new*sq_inv_norm_res_ini ) << std::endl;
		if( sq_norm_res_new*sq_inv_norm_res_ini < conv_ratio_tol*conv_ratio_tol ){
			conv_ratio = sqrt( sq_norm_res_new*sq_inv_norm_res_ini );
			num_iter = iitr;
			break;
		}

		////////////////////////////////
		// Calc Beta
		// beta = - (r,Ap) / (p,Ap)
		////////////////////////////////
		ls.MatVec_Hermitian(1.0,ir,0.0,iz);
		double sq_norm_z1 = ls.INPROCT(iz,iz).Real();
		const double beta = sq_norm_z1 / sq_norm_z0;
		sq_norm_z0 = sq_norm_z1;

		////////////////////////////////
		// Update Search Direction
		// {p} = {r} + beta*{p}
		////////////////////////////////
		ls.SCAL(beta,ip);
		ls.AXPY(1.0,iz,ip);
	}

	return true;
}


////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool Fem::Ls::Solve_BiCGSTAB(double& conv_ratio, unsigned int& iteration,
					CZLinearSystem& ls)
{

//	std::cout.precision(18);

	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	if( ls.GetTmpVectorArySize() < 5 ){ ls.ReSizeTmpVecSolver(5); }

	const int ix = -2;
	const int ir = -1;
	const int is = 0;
	const int iAs = 1;
	const int ip = 2;
	const int iAp = 3;
	const int ir0 = 4;

	ls.SCAL(0,ix);
//	u_vec.SetVectorZero();

	double sq_inv_norm_res;
	{
//		double dtmp1 = r_vec.GetSquaredVectorNorm();
		double dtmp1 = ls.INPROCT(ir,ir).Real();
		std::cout << "Initial Residual: " << sqrt(dtmp1) << " " << dtmp1 << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

//	r0_conjugate_vec = r_vec;
	ls.COPY(ir,ir0);

	// {p} = {r}
//	p_vec = r_vec;
	ls.COPY(ir,ip);

	// calc (r,r0*)
//	Com::Complex r_r0conj = InnerProduct(r_vec,r0_conjugate_vec);
	Com::Complex r_r0conj = ls.INPROCT(ir,ir0);

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// calc {Ap} = [A]*{p}
//		mat.MatVec(1.0,p_vec,0.0,Ap_vec);
		ls.MatVec(1,ip,0,iAp);

		// calc alpha
		Com::Complex alpha;
		{
//			const Com::Complex denominator = InnerProduct(Ap_vec,r0_conjugate_vec);
			const Com::Complex denominator = ls.INPROCT(iAp,ir0);
			alpha = r_r0conj / denominator;
		}

//		std::cout << "Alpha    " << alpha.Real() << " " << alpha.Imag() << std::endl;

		// calc s_vector
//		s_vec = r_vec;
//		s_vec.AXPY(-alpha,Ap_vec);
		ls.COPY(ir,is);	
		ls.AXPY(-alpha,iAp,is);

		// calc {As} = [A]*{s}
//		mat.MatVec(1.0,s_vec,0.0,As_vec);
		ls.MatVec(1,is,0,iAs);

		// calc omega
		Com::Complex omega;
		{
//			const double denominator = As_vec.GetSquaredVectorNorm();
//			const Com::Complex numerator = InnerProduct(s_vec,As_vec);
			const double denominator = ls.INPROCT(iAs,iAs).Real();
			const Com::Complex numerator = ls.INPROCT(is,iAs);
			omega = numerator / denominator;
		}

//		std::cout << "Omega    " << omega.Real() << " " << omega.Imag() << std::endl;

		// update solution
//		u_vec.AXPY(alpha,p_vec);
//		u_vec.AXPY(omega,s_vec);
		ls.AXPY(alpha,ip,ix);
		ls.AXPY(omega,is,ix);

		// update residual
//		r_vec = s_vec;
//		r_vec.AXPY(-omega,As_vec);
		ls.COPY(is,ir);
		ls.AXPY(-omega,iAs,ir);

		{
//			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double sq_norm_res = ls.INPROCT(ir,ir).Real();
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res;
			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res );
				iteration = iitr;
				return true;
			}
		}

		// calc beta
		Com::Complex beta;
		{
//			const Com::Complex tmp1 = InnerProduct(r_vec,r0_conjugate_vec);
			const Com::Complex tmp1 = ls.INPROCT(ir,ir0);
			beta = (tmp1/r_r0conj) * (alpha/omega);
			r_r0conj = tmp1;
		}

//		std::cout << "Beta    " << beta.Real() << " " << beta.Imag() << std::endl;

		// update p_vector
//		p_vec *= beta;
//		p_vec += r_vec;
//		p_vec.AXPY(-beta*omega,Ap_vec);
		ls.SCAL(beta,ip);
		ls.AXPY(1,ir,ip);
		ls.AXPY(-beta*omega,iAp,ip);
	}

	return true;
}


////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool Fem::Ls::Solve_BiCGStabP(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond )
{
	const double conv_ratio_tol = conv_ratio;	
	const unsigned int mx_iter = iteration;

	if( ls.GetTmpVectorArySize() < 7 ){ ls.ReSizeTmpVecSolver(7); }

	const int ix = -2;
	const int ir = -1;
	const int is = 0;
	const int iMs = 1;
	const int iAMs = 2;
	const int ip = 3;
	const int iMp = 4;
	const int iAMp = 5;
	const int ir0 = 6;

//	u_vec.SetVectorZero();
	ls.SCAL(0,ix);

	double sq_inv_norm_res;
	{
//		double dtmp1 = r_vec.GetSquaredVectorNorm();
		double dtmp1 = ls.INPROCT(ir,ir).Real();
		std::cout << "Initial Residual : " << sqrt(dtmp1) << std::endl;
		if( dtmp1 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res = 1.0 / dtmp1;
	}

//	r0_conjugate_vec = r_vec;
//	r_vec.SetVectorConjugate();
	ls.COPY(ir,ir0);
//	ls.Conjugate(ir0);

	// {p} = {r}
//	p_vec = r_vec;
	ls.COPY(ir,ip);

	// calc (r,r0*)
//	Com::Complex r_r0conj = InnerProduct(r_vec,r0_conjugate_vec);
	Com::Complex r_r0conj = ls.INPROCT(ir,ir0);

	iteration = mx_iter;
	for(unsigned int iitr=1;iitr<mx_iter;iitr++){

		// {Mp_vec} = [M^-1]*{p}
//		Mp_vec = p_vec;
//		precond.SolvePrecond(mat_p,Mp_vec);
		ls.COPY(ip,iMp);
		precond.SolvePrecond(ls,iMp);

		// calc {AMp_vec} = [A]*{Mp_vec}
//		mat.MatVec(1.0,Mp_vec,0.0,AMp_vec);
		ls.MatVec(1,iMp,0,iAMp);

		// calc alpha
		Com::Complex alpha;
		{
//			const Com::Complex denominator = InnerProduct(AMp_vec,r0_conjugate_vec);
			const Com::Complex denominator = ls.INPROCT(iAMp,ir0);
			alpha = r_r0conj / denominator;
		}

		// calc s_vector
//		s_vec = r_vec;
//		s_vec.AXPY(-alpha,AMp_vec);
		ls.COPY(ir,is);
		ls.AXPY(-alpha,iAMp,is);

		// {Ms_vec} = [M^-1]*{s}
//		Ms_vec = s_vec;
//		precond.SolvePrecond(mat_p,Ms_vec);
		ls.COPY(is,iMs);
		precond.SolvePrecond(ls,iMs);

		// calc {AMs} = [A]*{Ms_vec}
//		mat.MatVec(1.0,Ms_vec,0.0,AMs_vec);
		ls.MatVec(1,iMs,0,iAMs);

		// calc omega
		Com::Complex omega;
		{
//			const double denominator = AMs_vec.GetSquaredVectorNorm();
			const double denominator = ls.INPROCT(iAMs,iAMs).Real();
//			const Com::Complex numerator = InnerProduct(s_vec,AMs_vec);
			const Com::Complex numerator = ls.INPROCT(is,iAMs);
			omega = numerator / denominator;
		}

		// update solution
//		u_vec.AXPY(alpha,Mp_vec);
//		u_vec.AXPY(omega,Ms_vec);
		ls.AXPY(alpha,iMp,ix);
		ls.AXPY(omega,iMs,ix);

		// update residual
//		r_vec = s_vec;
//		r_vec.AXPY(-omega,AMs_vec);
		ls.COPY(is,ir);
		ls.AXPY(-omega,iAMs,ir);

		{
//			const double sq_norm_res = r_vec.GetSquaredVectorNorm();
			const double sq_norm_res = ls.INPROCT(ir,ir).Real();
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
//			const Com::Complex tmp1 = InnerProduct(r_vec,r0_conjugate_vec);
			const Com::Complex tmp1 = ls.INPROCT(ir,ir0);
			beta = (alpha/omega) * (tmp1/r_r0conj);
			r_r0conj = tmp1;
		}

		// update p_vector
//		p_vec *= beta;
//		p_vec += r_vec;
//		p_vec.AXPY(-beta*omega,AMp_vec);
		ls.SCAL(beta,ip);
		ls.AXPY(1,ir,ip);
		ls.AXPY(-beta*omega,iAMp,ip);
	}

	return true;
}
