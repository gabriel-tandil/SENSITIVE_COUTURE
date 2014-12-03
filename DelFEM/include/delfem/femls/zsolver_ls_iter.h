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
@brief 複素数反復法関数群
@author Nobuyuki Umetani
*/

#if !defined(_Z_SOLVER_CG_H)
#define Z_SOLVER_CG_H

#include "delfem/femls/zpreconditioner.h"

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/zmatprecond_blk.h"

namespace Fem{
namespace Ls{

//! @ingroup FemLs
//@{
////////////////////////////////////////////////////////////////
// Solve Hermetizn matrix "ls" with Conjugate Gradient Method

//! CG法でエルミート行列を解く
bool Solve_CG(double& conv_ratio, unsigned int& num_iter, 
			  Fem::Ls::CZLinearSystem& ls);
			  
//! 前処理付きCG法でエルミート行列を解く
bool Solve_PCG(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond );
				
////////////////////////////////////////////////////////////////
// Solve Matrix with COCG Methods

//! COCG法で行列を解く
bool Solve_PCOCG(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond );

////////////////////////////////////////////////////////////////
// Solve Matrix with Conjugate Gradient NR Methods

//! CGNR法で行列を解く
bool Solve_CGNR(double& conv_ratio, unsigned int& num_iter,
				CZLinearSystem& ls);


////////////////////////////////////////////////////////////////
// Solve Matrix with BiCGSTAB Methods

//! BiCGSTAB法で行列を解く
bool Solve_BiCGSTAB(double& conv_ratio, unsigned int& iteration,
					CZLinearSystem& ls);
					
//! 前処理付きBiCGSTAB法で行列を解く
bool Solve_BiCGStabP(double& conv_ratio, unsigned int& iteration,
				CZLinearSystem& ls, CZPreconditioner& precond );

//@}
}
}




#endif


