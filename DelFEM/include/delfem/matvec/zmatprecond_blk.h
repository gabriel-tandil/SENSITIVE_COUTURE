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
@brief 複素行列に関する抽象前処理クラス(MatVec::CZPrecond_Blk)
@author Nobuyuki Umetani
*/

#if !defined(ZMAT_PRECOND_BLK_H)
#define ZMAT_PRECOND_BLK_H

#include <assert.h>

namespace MatVec{

class CZVector_Blk;
class CZMatDia_BlkCrs;
//! 複素行列に対する前処理行列クラス
class CZMatPrecond_Blk{
public:
	//! 前処理行列を解く
	virtual bool SolvePrecond(const CZMatDia_BlkCrs& mat, CZVector_Blk& vec) const = 0;
};

}

#endif // PRECONDITIONER_H


