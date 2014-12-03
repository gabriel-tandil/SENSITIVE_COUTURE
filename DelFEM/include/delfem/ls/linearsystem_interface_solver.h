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
@brief 連立一次方程式クラス(Sol::CLinearSystem_SolInterface, Sol::CLinearSystemPreconditioner_SolInterface, Sol::CLinearSystemPreconditioner_SolInterface)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(LINEARSYSTEM_INTERFACE_SOLVER_H)
#define LINEARSYSTEM_INTERFACE_SOLVER_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <assert.h>

namespace LsSol{

/*! 
@brief interface class of linear system for solver
@ingroup LsSol
*/
class ILinearSystem_Sol
{
public:
	////////////////////////////////
	// function for linear solver
	// v=-1:residual    v=-2:update

	//! ソルバに必要な作業ベクトルの数を得る
	virtual unsigned int GetTmpVectorArySize() const = 0;
	//! ソルバに必要な作業ベクトルの数を設定
	virtual bool ReSizeTmpVecSolver(unsigned int size_new) = 0;

	virtual double DOT(int iv1, int iv2) = 0; //!< ベクトルの内積 (return {v1} * {v2})
	virtual bool COPY(int iv1, int iv2) = 0; //!< ベクトルのコピー ({v2} := {v1})
	virtual bool SCAL(double alpha, int iv1) = 0; //!< ベクトルのスカラー倍 ({v1} := alpha * {v1})
	virtual bool AXPY(double alpha, int iv1, int iv2) = 0; //!< ベクトルの足し算({v2} := alpha*{v1} +　{v2})	
	virtual bool MATVEC(double alpha, int iv1, double beta, int iv2) = 0; //!< 行列ベクトル積 ({v2} := alpha*[MATRIX]*{v1} + beta*{v2})
};


/*! 
@brief Interface class of linear system and preconditioner for solver
@ingroup LsSol
*/
class ILinearSystemPreconditioner_Sol
{
public:
	////////////////////////////////
	// function for linear solver
	// v=-1:residual    v=-2:update

	//! ソルバに必要な作業ベクトルの数を得る
	virtual unsigned int GetTmpVectorArySize() const = 0;
	//! ソルバに必要な作業ベクトルの数を設定
	virtual bool ReSizeTmpVecSolver(unsigned int size_new) = 0;

	virtual double DOT(int iv1, int iv2) = 0; //!< ベクトルの内積 (return {v1} * {v2})
	virtual bool COPY(int iv1, int iv2) = 0; //!< ベクトルのコピー ({v2} := {v1})
	virtual bool SCAL(double alpha, int iv1) = 0; //!< ベクトルのスカラー倍 ({v1} := alpha * {v1})
	virtual bool AXPY(double alpha, int iv1, int iv2) = 0; //!< ベクトルの足し算({v2} := alpha*{v1} +　{v2})	
	virtual bool MATVEC(double alpha, int iv1, double beta, int iv2) = 0; //!< 行列ベクトル積 ({v2} := alpha*[MATRIX]*{v1} + beta*{v2})

	virtual bool SolvePrecond(int iv) = 0;
};



}

#endif

