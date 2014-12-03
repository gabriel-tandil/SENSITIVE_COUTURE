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

/*! @file
@brief 行列のための連立一次方程式ソルバ
@author Nobuyuki Umetani
*/

#if !defined(SOLVER_MAT_ITER_H)
#define SOLVER_MAT_ITER_H

#include <vector>

namespace MatVec{
	class CMatDia_BlkCrs;
	class CVector_Blk;
	class CPrecond_Blk;

namespace Sol
{

/*!
@addtogroup MatVec
*/
//@{	

//! 前処理無しのCG法
bool Solve_CG(double& conv_ratio, unsigned int& iteration,
			  const MatVec::CMatDia_BlkCrs& mat, MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec);

//! 前処理付きのCG法
bool Solve_PCG(double& conv_ratio, unsigned int& iteration,
			   const MatVec::CMatDia_BlkCrs& mat, MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec,
			   const MatVec::CPrecond_Blk& precond, const MatVec::CMatDia_BlkCrs& mat_p);


//! 前処理無しのBiCGSTAB法
bool Solve_BiCGSTAB(double& conv_ratio, unsigned int& iteration,
			  const MatVec::CMatDia_BlkCrs& mat, MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec);
//! 前処理付きのBiCGSTAB法			  
bool Solve_PBiCGSTAB(double& conv_ratio, unsigned int& iteration,
			   const MatVec::CMatDia_BlkCrs& mat, MatVec::CVector_Blk& r_vec, MatVec::CVector_Blk& u_vec,
			   const MatVec::CPrecond_Blk& precond, const MatVec::CMatDia_BlkCrs& mat_p);
//@}			   
}

}

#endif
