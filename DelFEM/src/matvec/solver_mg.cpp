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
#include <cstdlib> //(abort)
#include <math.h>

#include "delfem/matvec/solver_mg.h"
#include "delfem/matvec/matdia_blkcrs.h"

#include "delfem/matvec/matdiainv_blkdia.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"
#include "delfem/matvec/matprolong_blkcrs.h"
#include "delfem/matvec/vector_blk.h"

using namespace MatVec;

CSolverMG::CSolverMG(const CMatDia_BlkCrs& mat){
	m_pPro = 0; m_pRes = 0; m_pMat = 0;
	m_pDiaInv = 0; m_pFrac = 0;
	m_mask_coarse_GS_CF = 0;

	m_pVf1=0; m_pVf2=0;
	m_pV_c=0;

	m_LayerMG_Coarse = 0;

	m_Omega = 0.9;
	m_niter_pre = 1;
	m_niter_pos = 1;

	if( mat.NBlkMatCol() < 50 ){
        m_pFrac = new CMatDiaFrac_BlkCrs(-1,mat);
		m_pFrac->SetValue(mat);
		return;
	}

	// Make Prolongation
	m_pPro = new CMatProlong_BlkCrs(mat,0.25);
	if( m_pPro->NBlkMatRow() == 0 ){
		delete m_pPro; m_pPro = 0;
		m_pFrac = new CMatDiaFrac_BlkCrs(-1,mat);
		m_pFrac->SetValue(mat);
		return;
	}

	// Make Restriction
	m_pRes = new CMat_BlkCrs(*m_pPro,true,false);
			
	{	// Make Value of Matrix
		const unsigned int nblk_c = m_pPro->NBlkMatRow();
		m_pMat = new CMatDia_BlkCrs(nblk_c,1); 
		m_pMat->AddPattern( *m_pRes, mat, *m_pPro );
		m_pMat->SetValue(	*m_pRes, mat, *m_pPro );
	}

	{	// Allocate Vector
		const unsigned int nblk_c = m_pPro->NBlkMatRow();
		m_pV_c = new CVector_Blk(nblk_c,1);
		const unsigned int nblk_f = mat.NBlkMatCol();
		this->m_pVf1 = new CVector_Blk(nblk_f,1);
		this->m_pVf2 = new CVector_Blk(nblk_f,1);
	}

	{	// Smoother
		m_smoother_type = GS;
		m_pDiaInv = new CMatDiaInv_BlkDia(mat);
	}

	// Make Coarse Grid
	m_LayerMG_Coarse = new CSolverMG(*m_pMat);
}

CSolverMG::CSolverMG(const CMatDia_BlkCrs& mat, MG_SMOOTHER smoother_type)
{
	m_pPro = 0; m_pRes = 0; m_pMat = 0;
	m_pDiaInv = 0; m_pFrac = 0;
	m_mask_coarse_GS_CF = 0;

	this->m_pVf1 = 0; this->m_pVf2 = 0;
	this->m_pV_c = 0;

	m_LayerMG_Coarse = 0;

	this->m_Omega = 0.9;

	if( mat.NBlkMatCol() < 2000 ){
		m_pFrac = new CMatDiaFrac_BlkCrs(-1,mat);
		m_pFrac->SetValue(mat);
		return;
	}

	// Make Prorongation Matrix
	m_pPro = new CMatProlong_BlkCrs(mat,0.25);
	if( m_pPro->NBlkMatRow() == 0 ){
		delete m_pPro; m_pPro = 0;
		m_pFrac = new CMatDiaFrac_BlkCrs(-1,mat);
		m_pFrac->SetValue(mat);
		return;
	}

	// Make Restriction Matrix
	m_pRes = new CMat_BlkCrs(*m_pPro,true,false);		

	{	// Make Value of Matrix
		const unsigned int nblk_c = m_pPro->NBlkMatRow();
		m_pMat = new CMatDia_BlkCrs(nblk_c,1); 
		m_pMat->AddPattern( *m_pRes, mat, *m_pPro );
		m_pMat->SetValue(	*m_pRes, mat, *m_pPro );
	}

	{	// Allocate Vector
		const unsigned int nblk_c = m_pPro->NBlkMatRow();
		m_pV_c = new CVector_Blk(nblk_c,1);
		const unsigned int nblk_f = mat.NBlkMatCol();
		this->m_pVf1 = new CVector_Blk(nblk_f,1);
		this->m_pVf2 = new CVector_Blk(nblk_f,1);
	}

	{	// Smoother
		m_niter_pre = 1;
		m_niter_pos = 1;
		m_smoother_type = smoother_type;
		if( smoother_type == GS  || smoother_type == GS_FB || 
			smoother_type == SOR || 
			smoother_type == JAC || smoother_type == O_JAC ){
			m_pDiaInv = new CMatDiaInv_BlkDia(mat);
		}
		else if( smoother_type == ILU ){
			m_pFrac = new CMatDiaFrac_BlkCrs(0,mat);
			m_pFrac->SetValue(mat);
		}
		else{
			std::cout << "Error!-->Unknown Smoother type" << std::endl;
			assert(0);
			abort();
		}
	}

	// Make Coarse Grid
	m_LayerMG_Coarse = new CSolverMG(*m_pMat,smoother_type);

	return;
}

CSolverMG::~CSolverMG(){
	if( m_pMat    != 0 ) delete m_pMat;
	if( m_pFrac   != 0 ) delete m_pFrac;
	if( m_pDiaInv != 0 ) delete m_pDiaInv;
	if( m_mask_coarse_GS_CF != 0 ) delete m_mask_coarse_GS_CF;
	////////////////
	if( m_pV_c != 0 ) delete m_pV_c;
	if( m_pVf1 != 0 ) delete m_pVf1;
	if( m_pVf2 != 0 ) delete m_pVf2;
	////////////////
	if( m_LayerMG_Coarse != 0 ) delete m_LayerMG_Coarse;
}

unsigned int CSolverMG::SumOfCrsMatSize(const CMatDia_BlkCrs& mat) const
{
	if( m_pMat == 0 ){
		assert( m_LayerMG_Coarse == 0 );
		return m_pFrac->NCrs() + m_pFrac->NBlkMatCol();
	}
	assert( m_LayerMG_Coarse != 0 );
	return m_LayerMG_Coarse->SumOfCrsMatSize(*m_pMat) + mat.NCrs() + mat.NBlkMatCol();
}

bool CSolverMG::SolveCycleV(const CMatDia_BlkCrs& mat,CVector_Blk& v_f) const 
{
	if( m_LayerMG_Coarse == 0 ){
		// Calc Exact Solution
		assert( m_pFrac != 0 );
		m_pFrac->Solve(v_f);
		return true;
	}

	if( m_smoother_type == ILU ){

		{	// Calc Update of Pre Smoothing	
			if( m_niter_pre > 0 ){
				*m_pVf1 = v_f;
				assert( m_pFrac != 0 );
				m_pFrac->Solve(*m_pVf1);
				mat.MatVec(-1.0, *m_pVf1, 1.0, v_f);
				for(unsigned int iiter=1;iiter<m_niter_pre;iiter++){
					*m_pVf2 = v_f;
					m_pFrac->Solve(*m_pVf2);
					mat.MatVec(-1.0, *m_pVf2, 1.0, v_f);
					*m_pVf1 += *m_pVf2;
				}
			}
			else{
				m_pVf1->SetVectorZero();
			}
		}
		
		{	// Coarse Grid Corrction
			m_pRes->MatVec(1.0, v_f,0.0,*m_pV_c,true);			// Restriction
			assert( m_LayerMG_Coarse != 0 );
			m_LayerMG_Coarse->SolveCycleV(*m_pMat,*m_pV_c);	// Solve Coarse Grid
			m_pPro->MatVec(1.0,*m_pV_c,0.0,*m_pVf2,true);			// Prolongation
		}
		
		{	// Post Smoothing
			if( m_niter_pos == 0 ){
				v_f  = *m_pVf1;
				v_f += *m_pVf2;
				return true;
			}
			*m_pVf1 += *m_pVf2;
			mat.MatVec(-1.0,*m_pVf2,1.0,v_f);
			if( m_niter_pos == 1 ){
				assert( m_pFrac != 0 );
				m_pFrac->Solve(v_f);
				v_f += *m_pVf1;
				return true;
			}
			for(unsigned int iiter=0;iiter<m_niter_pos;iiter++){
				*m_pVf2 = v_f;
				assert( m_pFrac != 0 );
				m_pFrac->Solve(*m_pVf2);
				mat.MatVec(-1.0,*m_pVf2,1.0,v_f);
				*m_pVf1 += *m_pVf2;
			}
			v_f = *m_pVf1;
		}
		return true;
	}

	////////////////////////////////

	{	// Pre Smooting
		if( m_niter_pre > 0 ){
			assert( m_pDiaInv != 0 );
			if(     m_smoother_type==JAC  ){m_pDiaInv->Solve_IniJacobi(     mat,v_f,*m_pVf1        );}
			else if(m_smoother_type==GS || m_smoother_type==GS_FB ){
				m_pDiaInv->Solve_IniGaussSidel( mat,v_f,*m_pVf1        );
			}
			else if(m_smoother_type==O_JAC){m_pDiaInv->Solve_IniOmegaJacobi(mat,v_f,*m_pVf1,m_Omega);}
			else if(m_smoother_type==SOR  ){m_pDiaInv->Solve_IniSOR(        mat,v_f,*m_pVf1,m_Omega);}
			else{ assert(0); }
			if( m_niter_pre > 1 ){
				if(     m_smoother_type==JAC   ){m_pDiaInv->SolveUpdate_Jacobi(     m_niter_pre-1, mat,v_f,*m_pVf1,*m_pVf2        );}
				else if(m_smoother_type==GS || m_smoother_type==GS_FB ){
					m_pDiaInv->SolveUpdate_GaussSidel( m_niter_pre-1, mat,v_f,*m_pVf1                );
				}
				else if(m_smoother_type==O_JAC ){m_pDiaInv->SolveUpdate_OmegaJacobi(m_niter_pre-1, mat,v_f,*m_pVf1,*m_pVf2,m_Omega);}
				else if(m_smoother_type==SOR   ){m_pDiaInv->SolveUpdate_SOR(        m_niter_pre-1, mat,v_f,*m_pVf1,        m_Omega);}
				else{ assert(0); }
			}
			*m_pVf2 = v_f;
			mat.MatVec(-1.0, *m_pVf1, 1.0, *m_pVf2);
		}
		else{
			*m_pVf2 = v_f;
			m_pVf1->SetVectorZero();
		}
	}

	{	// Coarse Grid Corrction		
		m_pRes->MatVec(1.0, *m_pVf2,0.0,*m_pV_c,true);		// Restriction
		assert( m_LayerMG_Coarse != 0 );
		m_LayerMG_Coarse->SolveCycleV(*m_pMat,*m_pV_c);	// Solve Coarse Grid
		m_pPro->MatVec(1.0,*m_pV_c,1.0,*m_pVf1,true);			// Prolongation
	}			
		
	{	// Post Smoothing
		if( m_niter_pos > 0 ){
			assert( m_pDiaInv != 0 );
			if(     m_smoother_type==JAC  ){m_pDiaInv->SolveUpdate_Jacobi(     m_niter_pos, mat, v_f,*m_pVf1,*m_pVf2        );}
			else if(m_smoother_type==GS   ){m_pDiaInv->SolveUpdate_GaussSidel( m_niter_pos, mat, v_f,*m_pVf1                );}
			else if(m_smoother_type==O_JAC){m_pDiaInv->SolveUpdate_OmegaJacobi(m_niter_pos, mat, v_f,*m_pVf1,*m_pVf2,m_Omega);}
			else if(m_smoother_type==SOR  ){m_pDiaInv->SolveUpdate_SOR(        m_niter_pos, mat, v_f,*m_pVf1,        m_Omega);}
			else if(m_smoother_type==GS_FB){m_pDiaInv->SolveUpdate_GaussSidel_BackWard(  m_niter_pos, mat, v_f,*m_pVf1                );}
			else{ assert(0); }
		}
		v_f = *m_pVf1;
	}
	return true;

}


