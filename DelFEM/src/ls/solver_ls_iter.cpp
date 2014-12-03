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
#include <math.h>

#include "delfem/ls/solver_ls_iter.h"
#include "delfem/ls/linearsystem_interface_solver.h"

using namespace LsSol;

////////////////////////////////////////////////////////////////
// Solve Matrix with CG Methods
////////////////////////////////////////////////////////////////
bool LsSol::Solve_CG(double& conv_ratio, unsigned int& num_iter, LsSol::ILinearSystem_Sol& ls)
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
		sq_norm_res = ls.DOT(ir,ir);
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
			ls.MATVEC(1.0,ip,0.0,iAp);
			const double pAp = ls.DOT(ip,iAp);
			alpha = sq_norm_res / pAp;
		}

		// update x 
		// x = x + alpha*p
		ls.AXPY(alpha,ip,ix);

		// update residual
		// r = r - alpha*Ap
		ls.AXPY(-alpha,iAp,ir);

		// calc residual 
		const double sq_norm_res_new = ls.DOT(ir,ir);
//		std::cout << num_iter << " " << sqrt( sq_norm_res_new*sq_inv_norm_res_ini ) << std::endl;
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


bool LsSol::Solve_PCG(double& conv_ratio, unsigned int& iteration,
                LsSol::ILinearSystemPreconditioner_Sol& ls)
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
		const double sq_norm_res0 = ls.DOT(ir,ir);
		if( sq_norm_res0 < 1.0e-30 ){
			conv_ratio = 0.0;
			iteration = 0;
			return true;			
		}
		sq_inv_norm_res0 = 1.0 / sq_norm_res0;
	}

	ls.COPY(ir,iz);
	ls.SolvePrecond(iz);

	ls.COPY(iz,ip);

	double inpro_rz = ls.DOT(ir,iz);
	for(unsigned int iitr=0;iitr<mx_iter;iitr++){

		double alpha;
		{	// calc alpha
			ls.MATVEC(1.0,ip,0.0,iz);
			const double val_pAp = ls.DOT(ip,iz);
			alpha = inpro_rz / val_pAp;
		}

		ls.AXPY(-alpha,iz,ir);
		ls.AXPY(alpha,ip,ix);	// Converge Judgement‚Ì‘O‚É“ü‚ê‚é

		{	// Converge Judgement
			double sq_norm_res = ls.DOT(ir,ir);
//			std::cout << iitr << " " << sqrt(sq_norm_res * sq_inv_norm_res0) << std::endl;
			if( sq_norm_res * sq_inv_norm_res0 < conv_ratio_tol*conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res0 );
				iteration = iitr;
				return true;
			}
		}

		double beta;
		{	// calc beta
			ls.COPY(ir,iz);
			ls.SolvePrecond(iz);
			const double inpro_rz_new = ls.DOT(ir,iz);
			beta = inpro_rz_new/inpro_rz;
			inpro_rz = inpro_rz_new;
		}

		ls.SCAL(beta,ip);
		ls.AXPY(1.0,iz,ip);
	}
	// Converge Judgement
  double sq_norm_res = ls.DOT(ir,ir);
  conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res0 );
  return false;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with CG Methods
////////////////////////////////////////////////////////////////
bool LsSol::Solve_BiCGSTAB(double& conv_ratio, unsigned int& num_iter, 
						   LsSol::ILinearSystem_Sol& ls)
{
	const unsigned int max_iter = num_iter;
	const double tolerance = conv_ratio;
	
	if( ls.GetTmpVectorArySize()<5 ){ ls.ReSizeTmpVecSolver(5); }

	const int ix = -2;
	const int ir = -1;
	const int is  = 0;	const int iAs = 1;
	const int ip  = 2;	const int iAp = 3;
	const int ir2 = 4;
    
    std::cout << "ini norm : " << sqrt(ls.DOT(ir,ir)) << std::endl;

	ls.SCAL(0.0,ix);

	double sq_inv_norm_res_ini;
	{
		const double sq_norm_res_ini = ls.DOT(ir,ir);
		if( sq_norm_res_ini < 1.0e-60 ){
			conv_ratio = sqrt( sq_norm_res_ini );
			num_iter = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
	}


	// {r2} = {r}
	ls.COPY(ir,ir2);

	// {p} = {r}
	ls.COPY(ir,ip);

	// calc ({r},{r2})
	double r_r2 = ls.DOT(ir,ir2);

	num_iter = max_iter;
	for(unsigned int iitr=1;iitr<max_iter;iitr++){

		// calc {Ap} = [A]*{p}
		ls.MATVEC(1.0,ip,0.0,iAp);
        std::cout << " sq_norm iAp : " << ls.DOT(iAp,iAp) << std::endl;

		// calc alpha
		// alhpa = ({r},{r2}) / ({Ap},{r2})
		double alpha;
		{
			const double denominator = ls.DOT(iAp,ir2);
            std::cout << " alpha deno : " << denominator << std::endl;
			alpha = r_r2 / denominator;
		}

		// {s} = {r} - alpha*{Ap}
		ls.COPY(ir,is);
		ls.AXPY(-alpha,iAp,is);

		// calc {As} = [A]*{s}
		ls.MATVEC(1.0,is,0.0,iAs);

		// calc omega
		// omega = ({As},{s}) / ({As},{As})
		double omega;
		{
			const double denominator = ls.DOT(iAs,iAs);
			const double numerator = ls.DOT(iAs,is);
			omega = numerator / denominator;
		}

		// update solution
		// ix += alpha*{p} + omega*{s}
		ls.AXPY(alpha,ip,ix);
		ls.AXPY(omega,is,ix);

		// update residual
		// {r} = {s} - omega*{As}
		ls.COPY(is,ir);
		ls.AXPY(-omega,iAs,ir);

		{
			const double sq_norm_res = ls.DOT(ir,ir);
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res_ini;
			std::cout << iitr << " " << sq_norm_res << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < tolerance*tolerance ){
				conv_ratio = sqrt( sq_conv_ratio );
				num_iter = iitr;
				break;
			}
		}

		// calc beta
		// beta = ({r},{r2})^new/({r},{r2})^old * alpha / omega
		double beta;
		{
			const double tmp1 = ls.DOT(ir,ir2);
			beta = (tmp1*alpha) / (r_r2*omega);
			r_r2 = tmp1;
		}

		// update p_vector
		// {p} = {r} + beta*({p}-omega*[A]*{p})
		ls.SCAL(beta,ip);
		ls.AXPY(1.0,ir,ip);
		ls.AXPY(-beta*omega,iAp,ip);
	}

	return true;
}

////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods
////////////////////////////////////////////////////////////////
bool LsSol::Solve_PBiCGSTAB(double& conv_ratio, unsigned int& num_iter, 
                          LsSol::ILinearSystemPreconditioner_Sol& ls)
{
	const double conv_ratio_tol = conv_ratio;	
	const unsigned int max_iter = num_iter;
	    
	if( ls.GetTmpVectorArySize()<7 ){ ls.ReSizeTmpVecSolver(7); }

	const int ix = -2;
	const int ir = -1;
	const int is  = 0;	const int iMs = 1;	const int iAMs = 2;
	const int ip  = 3;	const int iMp = 4;	const int iAMp = 5;
	const int ir2 = 6;

	double sq_inv_norm_res_ini;
	{
		const double sq_norm_res_ini = ls.DOT(ir,ir);
		if( sq_norm_res_ini < 1.0e-60 ){
			conv_ratio = sqrt( sq_norm_res_ini );
			num_iter = 0;
			return true;			
		}
		sq_inv_norm_res_ini = 1.0 / sq_norm_res_ini;
	}
	
//    std::cout << "SqIniRes : " << ls.DOT(ir,ir) << std::endl;

	// {u} = 0
	ls.SCAL(0.0,ix);
	
	// {r2} = {r}
	ls.COPY(ir,ir2);

	// {p} = {r}
	ls.COPY(ir,ip);

	num_iter = max_iter;
	for(unsigned int iitr=1;iitr<max_iter;iitr++)
	{
		// {Mp_vec} = [M^-1]*{p}
		ls.COPY(ip,iMp);
		ls.SolvePrecond(iMp);

		// calc (r,r0*)
		const double r_r2 = ls.DOT(ir,ir2);

//        std::cout << "r_r2 : " << r_r2 << std::endl;   

		// calc {AMp_vec} = [A]*{Mp_vec}
		ls.MATVEC(1.0,iMp,0.0,iAMp);

		// calc alpha
		double alpha;
		{
			const double denominator = ls.DOT(iAMp,ir2);
			alpha = r_r2 / denominator;
		}

//        std::cout << "Alpha : " << alpha << std::endl;

		// calc s_vector
		ls.COPY(ir,is);
		ls.AXPY(-alpha,iAMp,is);

//        std::cout << "ir iAMp is " << ls.DOT(ir,ir) << " " << ls.DOT(iAMp,iAMp) << " " << ls.DOT(is,is) << std::endl;

		// {Ms_vec} = [M^-1]*{s}
		ls.COPY(is,iMs);
		ls.SolvePrecond(iMs);

//        std::cout << "Is iMs " << ls.DOT(is,is) << " " << ls.DOT(iMs,iMs) << std::endl;

		// calc {AMs_vec} = [A]*{Ms_vec}
		ls.MATVEC(1.0,iMs,0.0,iAMs);

		double omega;
		{	// calc omega
			const double denominator = ls.DOT(iAMs,iAMs);
			const double numerator = ls.DOT(is,iAMs);
//            std::cout << "Omega0 : " << denominator << " " << numerator << std::endl;
			omega = numerator / denominator;
		}

//        std::cout << "Omega : " << omega << std::endl;

		// update solution
		ls.AXPY(alpha,iMp,ix);
		ls.AXPY(omega,iMs,ix);

		// update residual
		ls.COPY(is,ir);
		ls.AXPY(-omega,iAMs,ir);

		{
			const double sq_norm_res = ls.DOT(ir,ir);
			const double sq_conv_ratio = sq_norm_res * sq_inv_norm_res_ini;
//			std::cout << iitr << " " << sqrt(sq_conv_ratio) << " " << sqrt(sq_norm_res) << std::endl;
			if( sq_conv_ratio < conv_ratio_tol * conv_ratio_tol ){
				conv_ratio = sqrt( sq_norm_res * sq_inv_norm_res_ini );
				num_iter = iitr;
				return true;
			}
		}

		double beta;
		{	// calc beta
			const double tmp1 = ls.DOT(ir,ir2);
			beta = tmp1 * alpha / (r_r2*omega);
		}

		// update p_vector
		ls.SCAL(beta,ip);
		ls.AXPY(1.0,ir,ip);
		ls.AXPY(-beta*omega,iAMp,ip);
	}

	return true;
}

