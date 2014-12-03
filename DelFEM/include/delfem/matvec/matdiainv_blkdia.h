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
@brief 対角逆行列クラス(MatVec::CMatDiaInv_BlkDia)のインターフェイス
@author Nobuyuki Umetani
*/



#if !defined(MAT_DIA_INV_H)
#define MAT_DIA_INV_H

#include "delfem/matvec/matdia_blkcrs.h"

namespace MatVec{

/*! 
@brief 対角逆行列クラス
@ingroup MatVec
*/
class CMatDiaInv_BlkDia : public CMatDia_BlkCrs{
public:
	CMatDiaInv_BlkDia(unsigned int nblk, unsigned int blklen);
	CMatDiaInv_BlkDia(const CMatDia_BlkCrs& mat);
	virtual  ~CMatDiaInv_BlkDia();

	int SetValue(const CMatDia_BlkCrs& mat);

	int SolveUpdate_GaussSidel(  const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update                              ) const;
	int SolveUpdate_SOR(         const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update,                double omega ) const;
	int SolveUpdate_Jacobi(      const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, CVector_Blk& tmp_vec               ) const;
	int SolveUpdate_OmegaJacobi( const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, CVector_Blk& tmp_vec, double omega ) const;
	int SolveUpdate_GaussSidel_BackWard(   const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update ) const;
	int SolveUpdate_GaussSidel_Mask(  const unsigned int max_iter, const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, unsigned int* mask, bool is_mask_later ) const;

	int Solve_IniGaussSidel(  const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update               ) const;
	int Solve_IniSOR(         const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, double omega ) const;
	int Solve_IniJacobi(      const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update               ) const;
	int Solve_IniOmegaJacobi( const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, double omega ) const;
	int Solve_IniGaussSidel_Mask(  const CMatDia_BlkCrs& mat, const CVector_Blk& residual, CVector_Blk& update, unsigned int* mask, bool is_mask_later ) const;
private:
};

}	// end of namespace Ls





#endif
