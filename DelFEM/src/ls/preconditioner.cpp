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

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )   // C4786なんて表示すんな( ﾟДﾟ)ｺﾞﾙｧ
#endif
#define for if(0); else for

#include <assert.h>
#include <time.h>
#include <stdio.h>

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"
#include "delfem/matvec/matfrac_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/solver_mg.h"
#include "delfem/matvec/ordering_blk.h"

#include "delfem/ls/preconditioner.h"

void LsSol::CPreconditioner_ILU::Clear()
{
	for(unsigned int i=0;i<m_Matrix_NonDia.size();i++){
		for(unsigned int j=0;j<m_Matrix_NonDia[i].size();j++){
			if( m_Matrix_NonDia[i][j] != 0 ) delete m_Matrix_NonDia[i][j];
		}
		m_Matrix_NonDia[i].clear();
	}
	m_Matrix_NonDia.clear();
	////////////////
	for(unsigned int i=0;i<m_Matrix_Dia.size();i++){
		if( m_Matrix_Dia[i] != 0 ) delete m_Matrix_Dia[i];
	}
	m_Matrix_Dia.clear();
  ////////////////
  m_alev_input.clear();  
  m_is_ordering = false;
}
	
// symbolic factorization
void LsSol::CPreconditioner_ILU::SetLinearSystem(const CLinearSystem& ls)
{
//    std::cout << "0 prec : set linsys " << std::endl;
	if( m_is_ordering ){
    /*
		clock_t start,mid,end;
		start = clock();
		const unsigned int nlss = ls.GetNLinSysSeg();
		assert( nlss == 1 );
    assert( m_alev_input.size() == 0 || (m_alev_input.size() == 1 && m_alev_input[0].second == -1 ) );
    unsigned int lev = 0;
    if( m_alev_input.size() == 1 ){
      assert( m_alev_input[0].second == -1 );
      lev = m_alev_input[0].first;
    }
		m_order.MakeOrdering_AMD( ls.GetMatrix(0) );
//		m_order.MakeOrdering_RCM2( ls.GetMatrix(0) );
//		m_order.MakeOrdering_RCM( ls.GetMatrix(0) );
		mid = clock();
		end = clock();
     printf("Ordering:%.4f  Pattern:%.4f\n",(double)(mid-start)/CLOCKS_PER_SEC,(double)(end-mid)/CLOCKS_PER_SEC);     
     */
		m_Matrix_NonDia.resize(1);
		m_Matrix_NonDia[0].push_back(0);
    m_Matrix_Dia.push_back( new MatVec::CMatDiaFrac_BlkCrs(0,ls.GetMatrix(0),m_order) );    
		const unsigned int nblk = ls.GetMatrix(0).NBlkMatCol();
    const unsigned int nlen = ls.GetMatrix(0).LenBlkCol();
		m_vec.Initialize(nblk,nlen);
//		std::cout << ls.GetMatrix(0).NBlkMatCol() << " " << ls.GetMatrix(0).NCrs() << " " << m_Matrix_Dia[0]->NCrs() << std::endl;
		return;
	}

    ////////////////
//    std::cout << "1 prec : set linsys " << std::endl;

	const unsigned int nlss = ls.GetNLinSysSeg();
	m_Matrix_NonDia.resize(nlss);
	for(unsigned int ilss=0;ilss<nlss;ilss++){
        m_Matrix_Dia.push_back( new MatVec::CMatDiaFrac_BlkCrs(0,ls.GetMatrix(ilss)) );
		for(unsigned int jlss=0;jlss<nlss;jlss++){
		    if( ilss == jlss || !ls.IsMatrix(ilss,jlss) ){ 
			    m_Matrix_NonDia[ilss].push_back(0);
    			continue; 
	    	}
            m_Matrix_NonDia[ilss].push_back( new MatVec::CMatFrac_BlkCrs(ls.GetMatrix(ilss,jlss)) );
        }
    }

    if( m_alev_input.size() == 0 && m_afill_blk.empty() ) return;  // 指定が何も無ければ0 レベルのフィルイン
    
//    std::cout << "2 prec : set linsys " << std::endl;

    std::vector<int> alev;
    alev.resize(nlss*nlss,0);
    for(unsigned int iin=0;iin<m_alev_input.size();iin++){
        const int lev  = m_alev_input[iin].first;
        const int ilss = m_alev_input[iin].second;
        if( ilss == -1 ){
            for(unsigned int i=0;i<nlss*nlss;i++){ alev[i] = lev; }
        }
        else{
            for(unsigned int i=0;i<nlss;i++){ 
                alev[i*nlss+ilss] = lev; 
                alev[ilss*nlss+i] = lev; 
            }
        }
    }

    {
        bool iflag = false;
        for(unsigned int i=0;i<nlss*nlss;i++){ 
            if( alev[i] != 0 ) iflag = true;
        }
        if( !iflag && m_afill_blk.empty() ) return;
    }
    
//    std::cout << "3 prec : set linsys " << std::endl;

    for(unsigned int ilss=0;ilss<nlss;ilss++){
    for(unsigned int jlss=0;jlss<nlss;jlss++){
//        std::cout << "add frac ptn" << ilss << " " << jlss << std::endl;
        const int lev_ij = alev[ilss*nlss+jlss];
        if( ilss == jlss ){
            assert( m_Matrix_Dia[ilss] != 0 );
            for(unsigned int klss=0;klss<ilss;klss++){
                assert( klss < ilss );
                if( m_Matrix_NonDia[ilss][klss] == 0 || m_Matrix_NonDia[klss][jlss] == 0) continue;
//                std::cout << "   add frac ptn low up " << ilss << " " << klss << " " << jlss << std::endl;
                m_Matrix_Dia[ilss]->AddFracPtnLowUp(lev_ij, 
                    *m_Matrix_NonDia[ilss][klss],  *m_Matrix_NonDia[klss][jlss] );
            }
            std::cout << "   add frac ptn dia " << ilss << " " << ilss << " " << jlss << std::endl;
            if( ilss == 0 && m_afill_blk.size() > 0 ){
                m_Matrix_Dia[ilss]->AddFracPtn(lev_ij,m_afill_blk);
            }
            else{
                m_Matrix_Dia[ilss]->AddFracPtn(lev_ij);
            }
        }
        else{
            const unsigned int kmax = ( ilss < jlss ) ? ilss : jlss;
            for(unsigned int klss=0;klss<kmax+1;klss++){
                if( klss == ilss && klss < jlss ){
                    if( m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
//                    std::cout << "   add frac ptn up " << ilss << " " << klss << " " << jlss << std::endl;
                    m_Matrix_NonDia[ilss][jlss]->AddFracPtnUp( *m_Matrix_Dia[klss], lev_ij);
                    continue;
                }
                if( klss == jlss && klss < ilss ){
                    if( m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
//                    std::cout << "   add frac ptn low " << ilss << " " << klss << " " << jlss << std::endl;
                    m_Matrix_NonDia[ilss][jlss]->AddFracPtnLow(*m_Matrix_Dia[klss], lev_ij);
                    continue;
                }
                if( m_Matrix_NonDia[ilss][klss] == 0 || m_Matrix_NonDia[klss][jlss] == 0 ){ continue; }
                if( m_Matrix_NonDia[ilss][jlss] == 0 ){
                    const MatVec::CVector_Blk& res_i = ls.GetVector(-1,ilss);
                    const MatVec::CVector_Blk& res_j = ls.GetVector(-1,jlss);
                    const unsigned int nblk_i = res_i.NBlk();
                    const unsigned int nblk_j = res_j.NBlk();
                    if( res_i.Len() >= 0 && res_i.Len() >= 0 ){
//                        std::cout << "add matrix " << ilss << " " << jlss << std::endl;
                        m_Matrix_NonDia[ilss][jlss] 
                            = new MatVec::CMatFrac_BlkCrs(nblk_i,res_i.Len(), nblk_j, res_j.Len() );
                    }
                    else{
                        std::cout << "Error!-->Not Implemented" << std::endl;
                        assert(0);
                    }
                }
//                std::cout << "   add frac ptn low up" << ilss << " " << klss << " " << jlss << std::endl;
                m_Matrix_NonDia[ilss][jlss]
                    ->AddFracPtnLowUp(*m_Matrix_NonDia[ilss][klss], *m_Matrix_NonDia[klss][jlss], lev_ij);
            }
        }
    }
    }

//    std::cout << "4 prec : set linsys " << std::endl;

    for(unsigned int ilss=0;ilss<nlss;ilss++){
        assert( m_Matrix_Dia[ilss] != 0 );
        m_Matrix_Dia[ilss]->MakePatternFinalize();
        for(unsigned int jlss=0;jlss<nlss;jlss++){
            if( ilss == jlss ){ continue; }
            if( m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
//            std::cout << "Pattern finalize " << ilss << " " << jlss << std::endl;
            m_Matrix_NonDia[ilss][jlss]->MakePatternFinalize();
        }
    }

//    std::cout << "5 prec : set linsys " << std::endl;
    
//    std::cout << "End Finalize Pattern " << std::endl;
}


// numerical factorization
// 値を設定してILU分解を行う関数
// ILU分解が成功しているかどうかはもっと詳細なデータを返したい
bool LsSol::CPreconditioner_ILU::SetValue(const LsSol::CLinearSystem& ls)
{
  
  //    std::cout << "0 prec : set linsys " << std::endl;
  //    std::cout << "SetValue and LU decompose " << std::endl;
  
	const unsigned int nlss = ls.GetNLinSysSeg();
	if( m_is_ordering ){
		assert( nlss == 1 );
		m_Matrix_Dia[0]->SetValue_Initialize( ls.GetMatrix(0), m_order );
		if( !m_Matrix_Dia[0]->DoILUDecomp() ){ return false; }
		return true;
	}
  
	////////////////
	// 値をコピー
	for(unsigned int ilss=0;ilss<nlss;ilss++){
    //        std::cout << "   0SetValue Dia : " << ilss << std::endl;
    assert( m_Matrix_Dia[ilss] != 0 );
		m_Matrix_Dia[ilss]->SetValue_Initialize( ls.GetMatrix(ilss) );
    //        std::cout << "   1SetValue Dia : " << ilss << std::endl;
		for(unsigned int jlss=0;jlss<nlss;jlss++){
      if( ilss == jlss ){ continue; }
      if( !ls.IsMatrix(ilss,jlss) && m_Matrix_NonDia[ilss][jlss] != 0 ){ 
        m_Matrix_NonDia[ilss][jlss] -> SetZero();
        continue;
      }
      if( !ls.IsMatrix(ilss,jlss) ) continue;
      //            std::cout << "   0Set Value : " << ilss << " " << jlss << std::endl;
      m_Matrix_NonDia[ilss][jlss]->SetValue_Initialize( ls.GetMatrix(ilss,jlss) );
      //            std::cout << "   1Set Value : " << ilss << " " << jlss << std::endl;
    }
	}
  
  //    std::cout << "1 prec : set linsys " << std::endl;
	////////////////
	// ILU分解
  for(unsigned int ilss=0;ilss<nlss;ilss++){
    for(unsigned int jlss=0;jlss<nlss;jlss++){
      if( ilss == jlss ){
        for(unsigned int klss=0;klss<ilss;klss++){
          if( m_Matrix_NonDia[ilss][klss] && m_Matrix_NonDia[klss][ilss] ){
            //					std::cout << "ILU Frac Dia LowUp:" << ilss << " " << klss << " " << ilss << std::endl;
            if( !m_Matrix_Dia[ilss]->DoILUDecompLowUp( 
                                                      *m_Matrix_NonDia[ilss][klss],
                                                      *m_Matrix_NonDia[klss][ilss] ) ){ return false; }
          }
        }
        //		    std::cout << "ILU Frac Dia:" << ilss << std::endl;
        if( !m_Matrix_Dia[ilss]->DoILUDecomp() ){ 
          std::cout << "ilu frac false 33 matrix non-ordered" << std::endl;
          return false; 
        }
      }
      else{
        if( m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
        const unsigned int kmax = ( ilss < jlss ) ? ilss : jlss;
        for(unsigned int klss=0;klss<kmax+1;klss++){
          if( klss == ilss && klss < jlss ){
            //				    std::cout << "ILU Frac Up:" << ilss << " " << klss << " " << jlss << std::endl;
            if( !m_Matrix_NonDia[ilss][jlss]->DoILUDecompUp(  *m_Matrix_Dia[ilss] ) ){ return false; }
            continue;
          }
          if( klss == jlss && klss < ilss ){
            //				    std::cout << "ILU Frac Low:" << ilss << " " << klss << " " << jlss << std::endl;
            if( !m_Matrix_NonDia[ilss][jlss]->DoILUDecompLow( *m_Matrix_Dia[jlss] ) ){ return false; }
            continue;
          }
          if( m_Matrix_NonDia[ilss][klss] == 0 || m_Matrix_NonDia[klss][jlss] == 0 ){ continue; }
          assert( klss < ilss );
          assert( klss < jlss );
          //				std::cout << "ILU Frac LowUp:" << ilss << " " << klss << " " << jlss << std::endl;
          if( !m_Matrix_NonDia[ilss][jlss]->DoILUDecompLowUp( 
                                                             *m_Matrix_NonDia[ilss][klss],
                                                             *m_Matrix_NonDia[klss][jlss] ) ){ return false; }
        }
      }
    }
  }
	return true;  
}


// Solve Preconditioning System
bool LsSol::CPreconditioner_ILU::SolvePrecond(LsSol::CLinearSystem& ls, unsigned int iv)
{
	const unsigned int nlss = ls.GetNLinSysSeg();

	if( m_is_ordering ){
		assert( nlss == 1 );
		m_order.OrderingVector_OldToNew(m_vec,ls.GetVector(iv,0));
		m_Matrix_Dia[0]->ForwardSubstitution( m_vec);
		m_Matrix_Dia[0]->BackwardSubstitution(m_vec);
		m_order.OrderingVector_NewToOld(ls.GetVector(iv,0),m_vec);
		return true;
	}

    ////////////////

    // Forward Substitution
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		for(unsigned int jlss=0;jlss<ilss;jlss++){
            if( !m_Matrix_NonDia[ilss][jlss] ){ continue; }
			m_Matrix_NonDia[ilss][jlss]->MatVec(
				-1.0,ls.GetVector(iv,jlss),1.0,ls.GetVector(iv,ilss),true);
		}
		m_Matrix_Dia[ilss]->ForwardSubstitution(ls.GetVector(iv,ilss));
	}
    // Backward Substitution
	for(int ilss=(int)nlss-1;ilss>=0;ilss--){
		for(unsigned int jlss=ilss+1;jlss<nlss;jlss++){
            if( !m_Matrix_NonDia[ilss][jlss] ){ continue; }
			m_Matrix_NonDia[ilss][jlss]->MatVec(
				-1.0,ls.GetVector(iv,jlss),1.0,ls.GetVector(iv,ilss),true);
		}
		m_Matrix_Dia[ilss]->BackwardSubstitution(ls.GetVector(iv,ilss));
    }
	return true;
}

