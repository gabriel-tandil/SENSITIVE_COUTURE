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

#if !defined(SOLVERMG_H)
#define SOLVERMG_H

#include <vector>

#include "delfem/matvec/matprecond_blk.h"

namespace field{
	class CElemAry_SquarePow2;
}

namespace MatVec{

class CMatDia_BlkCrs;
class CVector_Blk;
class CPrecond_Blk;
class CMat_BlkCrs;
class CMatDia_BlkCrs;
class CMatDiaFrac_BlkCrs;
class CMatDiaInv_BlkDia;
class CMatProlong_BlkCrs;
class CVector_Blk;

enum MG_SMOOTHER{ GS, JAC, O_JAC, ILU, SOR, GS_FB };

class CSolverMG : public CPrecond_Blk{
public:
	CSolverMG(const CMatDia_BlkCrs& mat,MG_SMOOTHER smoother_type);
	CSolverMG(const field::CElemAry_SquarePow2& ea, const CMatDia_BlkCrs& mat,MG_SMOOTHER smoother_type );
	CSolverMG(const CMatDia_BlkCrs& mat);
	~CSolverMG();

	bool SolveCycleV(const CMatDia_BlkCrs& mat,CVector_Blk& vec) const;
	bool SetSmoothingNumberOfTimes(const unsigned int& npre, const unsigned int& npos){
		this->m_niter_pre = npre;
		this->m_niter_pos = npos;
		if( this->m_LayerMG_Coarse != 0 ){
			this->m_LayerMG_Coarse->SetSmoothingNumberOfTimes(npre,npos);
		}
		return true;
	}

	virtual bool SolvePrecond(const CMatDia_BlkCrs& mat, CVector_Blk& vec) const{
		return this->SolveCycleV(mat,vec);
	}
	unsigned int SumOfCrsMatSize(const CMatDia_BlkCrs& mat) const;

private:
	CMatProlong_BlkCrs* m_pPro;
	CMat_BlkCrs* m_pRes;
	CMatDia_BlkCrs* m_pMat;

	CVector_Blk* m_pV_c;
	CVector_Blk* m_pVf1;
	CVector_Blk* m_pVf2;
	CVector_Blk* m_pVf3;	// Use only ILU_CF

	CSolverMG* m_LayerMG_Coarse;

	CMatDiaFrac_BlkCrs* m_pFrac;
	CMatDiaInv_BlkDia* m_pDiaInv;
	unsigned int* m_mask_coarse_GS_CF;
	MG_SMOOTHER m_smoother_type;
	double m_Omega;
	unsigned int m_niter_pre;
	unsigned int m_niter_pos;
};

}
#endif // !defined(SOLVERMG_H)
