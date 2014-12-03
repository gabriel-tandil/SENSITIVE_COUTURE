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

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 ) 
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif  /* _WIN32 */
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

//#include <GL/glut.h>



#include <iostream>
#include "delfem/rigid/linearsystem_rigidfield.h"
#include "delfem/rigid/rigidbody.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/field_world.h"
#include "delfem/field.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
void Ls::CLinearSystem_RigidField::Solve(Ls::CPreconditioner_RigidField& prec)
{
    Ls::CPreconditioner_RigidBody_CRS2& prec_rigid = prec.GetPreconditionerRigid();
    {
        ls_rigid.COPY(-1,-2);
        prec_rigid.Solve( ls_rigid.GetVector(-2) );
    }

//    ls_field.GetVector(-2,0).SetVectorZero();
//    ls_field.GetVector(-2,1).SetVectorZero();
    {	// 行列を解く
        double conv_ratio = 1.0e-9;
        unsigned int max_iter = 10000;
		// Solve with Preconditioned Conjugate Gradient
		Sol::Solve_CG(conv_ratio,max_iter,ls_field);
		std::cout << max_iter << " " << conv_ratio << std::endl;
	}
}

*/

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// add pattern to the diagonal
bool Ls::CLinearSystem_RigidField2::AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
  const Fem::Field::CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}

	const unsigned int nlen_value = field.GetNLenValue();

	int ilss_add;
	{	// Cornerブロックを作る
    unsigned int id_na_val = field.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
      const Fem::Field::CNodeAry& na = world.GetNA(id_na_val);
      assert( m_ls.GetNLinSysSeg() == this->m_aSegRF.size() );
      ilss_add = m_ls.GetNLinSysSeg();
      this->m_aSegRF.push_back( CLinSysSegRF(id_field_parent,Fem::Field::CORNER) );
      int ires = m_ls.AddLinSysSeg( na.Size(), field.GetNLenValue() );
      assert( ires == ilss_add );
		}
		else{ ilss_add = -1; }
	}
	////////////////////////////////
	const std::vector<unsigned int> aIdEA = field.GetAryIdEA();
  assert( aIdEA.size() > 0 );
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
    const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
		// CORNER節点について
    if( field.GetIdElemSeg(id_ea,Fem::Field::CORNER,true,world) != 0 ){
			assert( world.IsIdEA(id_ea) );
      const unsigned int id_es_c = field.GetIdElemSeg(id_ea,Fem::Field::CORNER,true,world);
			assert( ea.IsSegID(id_es_c) );
      {
        Com::CIndexedArray crs;
        ea.MakePattern_FEM(id_es_c,crs);
        m_ls.AddMat_Dia(ilss_add, crs );			// cc行列を作る
      }
      if( field.GetIdElemSeg(id_ea,Fem::Field::BUBBLE,true,world) != 0 ){	// CORNER-BUBBLE
        assert(0);
			}
      if( field.GetIdElemSeg(id_ea,Fem::Field::EDGE,true,world) != 0 ){	// CONRER-EDGE
        assert(0);
			}
		}
		// EDGE節点について
    if( field.GetIdElemSeg(id_ea,Fem::Field::EDGE,true,world) != 0 ){
      assert(0);
		}
		// BUBBLE節点について
    if( field.GetIdElemSeg(id_ea,Fem::Field::BUBBLE,true,world) != 0 ){
      assert(0);
		}
	}
	return true;
}

// fieldで初期化する、fieldの中の非ゼロパターンで行列を作る
bool Ls::CLinearSystem_RigidField2::AddPattern_Field(
        unsigned int id_field1, 
        unsigned int id_field2, 
        const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field1) ) return false;
    const Fem::Field::CField& field1 = world.GetField(id_field1);
	unsigned int id_field_parent;
	{
		if( field1.GetIDFieldParent() == 0 ){ id_field_parent = id_field1; }
		else{ id_field_parent = field1.GetIDFieldParent(); }
	}
	const unsigned int nlen_value = field1.GetNLenValue();

	int ilss_add_c;
    int ilss_add_b;
	{	// Cornerセグメントを作る
        unsigned int id_na_val = field1.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
            const Fem::Field::CNodeAry& na = world.GetNA(id_na_val);
            assert( m_ls.GetNLinSysSeg() == this->m_aSegRF.size() );
            ilss_add_c = m_ls.GetNLinSysSeg();
            this->m_aSegRF.push_back( CLinSysSegRF(id_field_parent,Fem::Field::CORNER) );
            int ires = m_ls.AddLinSysSeg( na.Size(), field1.GetNLenValue() );
            assert( ires == ilss_add_c );
		}
		else{ ilss_add_c = -1; }
	}
	{	// Bubbleセグメントを作る
        unsigned int id_na_val = field1.GetNodeSegInNodeAry(Fem::Field::BUBBLE).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
            const Fem::Field::CNodeAry& na = world.GetNA(id_na_val);
            assert( m_ls.GetNLinSysSeg() == this->m_aSegRF.size() );
            ilss_add_b = m_ls.GetNLinSysSeg();
            this->m_aSegRF.push_back( CLinSysSegRF(id_field_parent,Fem::Field::BUBBLE) );
            int ires = m_ls.AddLinSysSeg( na.Size(), field1.GetNLenValue() );
            assert( ires == ilss_add_b );
        }
	}

    const Fem::Field::CField& field2 = world.GetField(id_field2);
	const std::vector<unsigned int>& aIdEA1 = field1.GetAryIdEA();
	const std::vector<unsigned int>& aIdEA2 = field2.GetAryIdEA();

	for(;;){	// ダミーのfor文を使ってbreakで抜けられるようにする
		// Corner-Corner関係を作る
		if( ilss_add_c == -1 ) break;
        const unsigned int id_na_va1 = field1.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
        const unsigned int id_na_co1 = field1.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_co;
        const unsigned int id_na_va2 = field2.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
        const unsigned int id_na_co2 = field2.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_co;
        assert( aIdEA1.size() > 0 );
        if( id_na_co1 != id_na_co2 ) break;

		for(unsigned int iiea1=0;iiea1<aIdEA1.size();iiea1++){
			const unsigned int id_ea1 = aIdEA1[iiea1];
            if( field1.GetIdElemSeg(id_ea1,Fem::Field::CORNER,true,world) == 0 ) continue;
			assert( world.IsIdEA(id_ea1) );
            const Fem::Field::CElemAry& ea1 = world.GetEA(id_ea1);
            const unsigned int id_es_c1 = field1.GetIdElemSeg(id_ea1,Fem::Field::CORNER,true,world);
			assert( ea1.IsSegID(id_es_c1) );
			// CORNER1-CORNER1
            {
                Com::CIndexedArray crs;
                ea1.MakePattern_FEM(id_es_c1,crs);
		        m_ls.AddMat_Dia(ilss_add_c, crs );			// cc行列を作る
            }
			////////////////
			for(unsigned int iiea2=0;iiea2<aIdEA2.size();iiea2++)
			{
				const unsigned int id_ea2 = aIdEA2[iiea2];
				if( id_ea1 == id_ea2 ){
					// CORNER1-CORNER2
                    const int ils2_c = this->FindIndexArray_Seg(id_field2,Fem::Field::CORNER,world);
					if( ils2_c != -1 ){
						assert( ils2_c >= 0 && ils2_c < this->GetNLinSysSeg() );
                        const unsigned int id_es_c2 = field2.GetIdElemSeg(id_ea2,Fem::Field::CORNER,true,world);
						assert( id_es_c1 == id_es_c2 );
						Com::CIndexedArray crs;
						ea1.MakePattern_FEM(id_es_c1,id_es_c2,crs);
						assert( crs.CheckValid() );
						m_ls.AddMat_NonDia(ilss_add_c,ils2_c, crs);		// c1c2行列を足す
						const unsigned int nnode2 = m_ls.GetBlkSizeLinSeg(ils2_c);
						Com::CIndexedArray crs_inv;
						crs_inv.SetTranspose( nnode2, crs );
						m_ls.AddMat_NonDia(ils2_c,ilss_add_c, crs_inv);	// c2c1行列を足す
					}
					// CORNER1-BUBBLE2
          const int ils2_b = this->FindIndexArray_Seg(id_field2,Fem::Field::BUBBLE,world);
					if( ils2_b != -1 ){
            assert(0);
					}
				}
				else{
          const Fem::Field::CNodeAry& na1 = world.GetNA(id_na_co1);
          const unsigned int id_es_c_co1 = field1.GetIdElemSeg(id_ea1,Fem::Field::CORNER,false,world);
          const unsigned int id_es_c_co2 = field2.GetIdElemSeg(id_ea2,Fem::Field::CORNER,false,world);
					if( na1.IsIncludeEaEs_InEaEs( 
						std::make_pair(id_ea1,id_es_c_co1),
						std::make_pair(id_ea2,id_es_c_co2) ) )
					{
						std::cout << "ea : " << id_ea1 << " is included in " << id_ea2 << std::endl;
						assert( ea1.IsSegID(id_es_c_co1) );
            Com::CIndexedArray crs;
						ea1.MakePattern_FEM(id_es_c1,id_es_c1,crs);	// 自分も含む
						assert( crs.CheckValid() );
						if( field2.IsPartial() ){
							std::cout << "Error!-->未実装" << std::endl;
							assert(0);
						}
						for(unsigned int icrs=0;icrs<crs.array.size();icrs++){
							unsigned int jno_va = crs.array[icrs];
							unsigned int jno_co = field1.GetMapVal2Co(jno_va);
//							std::cout << jno_van << " " << jno_co << std::endl;
							// 本当はfield2.GetMapVal2Coの逆写像を使って求める．
							unsigned int jno_va2 = jno_co;
							crs.array[icrs] = jno_va2;
						}
            int ils2 = this->FindIndexArray_Seg(id_field2,Fem::Field::CORNER,world);
						assert( ils2 >= 0 && ils2 < this->GetNLinSysSeg() );
//						std::cout << "ils_seg : " << ils0 << " " << ils2 << std::endl;
						m_ls.AddMat_NonDia(ilss_add_c,ils2, crs);
						const unsigned int nnode2 = m_ls.GetBlkSizeLinSeg(ils2);
						Com::CIndexedArray crs_inv;
						crs_inv.SetTranspose( nnode2, crs );
						m_ls.AddMat_NonDia(ils2,ilss_add_c, crs_inv);
					}
				}
			}
		}
		break;
	}
			
	for(unsigned int iiea=0;iiea<aIdEA1.size();iiea++){
		const unsigned int id_ea1 = aIdEA1[iiea];
        const unsigned int id_es_b1 = field1.GetIdElemSeg(id_ea1,Fem::Field::BUBBLE,true,world);
		if( id_es_b1 == 0 ) continue;
        const Fem::Field::CElemAry& ea1 = world.GetEA(id_ea1);
		assert( world.IsIdEA(id_ea1) );
		assert( ea1.IsSegID(id_es_b1) );
		const unsigned int ils0 = ilss_add_b;
		// BUBLE1-BUBBLE2
        {
			Com::CIndexedArray crs;
			ea1.MakePattern_FEM(id_es_b1,crs);
		    m_ls.AddMat_Dia(ilss_add_b, crs );
        }
		const unsigned int id_ea2 = aIdEA2[iiea];
		assert( id_ea1 == id_ea2 );
        const unsigned int id_es_c2 = field2.GetIdElemSeg(id_ea2,Fem::Field::CORNER,true,world);
		assert( ea1.IsSegID(id_es_c2) );
        int ils2 = this->FindIndexArray_Seg(id_field2,Fem::Field::CORNER,world);
		assert( ils2 >= 0 && ils2 < this->GetNLinSysSeg() );
		assert( ils2 >= 0 && ils2 < m_ls.GetNLinSysSeg() );
		{
			Com::CIndexedArray crs;
			ea1.MakePattern_FEM(id_es_b1,id_es_c2,crs);
			assert( crs.CheckValid() );
			m_ls.AddMat_NonDia(ils0,ils2, crs);		// b1c2行列を作る
			Com::CIndexedArray crs_inv;
			const unsigned nnode2 = m_ls.GetBlkSizeLinSeg(ils2);
			crs_inv.SetTranspose( nnode2, crs );
			m_ls.AddMat_NonDia(ils2,ils0, crs_inv);	// c2b1行列を作る
		}
	}

	// いろいろと追加が必要
	// まずは足りない部分を要求するようになったらエラーを出す関数を実装しよう。

	return true;
}


////////////////////////////////
// 剛体弾性体連成インターフェース
bool Ls::CLinearSystem_RigidField2::AddPattern_RigidField(
    const unsigned int id_field_l, const unsigned int id_field_u, const Fem::Field::CFieldWorld& world, 
    unsigned int irb, std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst)
{
  this->AddPattern_Field(id_field_l,id_field_u,world);
  const int ils_l = this->FindIndexArray_Seg(id_field_l,Fem::Field::CORNER,world);
  const int ils_u = this->FindIndexArray_Seg(id_field_u,Fem::Field::CORNER,world);
  assert( ils_l >= 0 && ils_l < this->GetNLinSysSeg() );
  assert( ils_u >= 0 && ils_u < this->GetNLinSysSeg() );
  assert( ilss_rigid >= 0 && ilss_rigid < this->GetNLinSysSeg() );
  Com::CIndexedArray crs;
	const unsigned int nblk_r = m_ls.GetBlkSizeLinSeg(ilss_rigid);
	const unsigned int nblk_l = m_ls.GetBlkSizeLinSeg(ils_l);
  {
//    crs.size = nblk_r;
    crs.index.clear();
    crs.index.resize( nblk_r+1,0);
    crs.index[0] = 0;
    for(unsigned int iblk=0;iblk<nblk_r;iblk++){
      if( iblk == irb ){
        crs.index[iblk+1] = crs.index[iblk] + nblk_l;
      }
      else{
        crs.index[iblk+1] = crs.index[iblk];
      }
    }
    crs.array.resize( nblk_l );
    for(unsigned int iblk=0;iblk<nblk_l;iblk++){
      crs.array[iblk] = iblk;
    }
    /*        for(unsigned int iblk=0;iblk<nblk_r;iblk++){
     for(unsigned int icrs=crs.index[iblk];icrs<crs.index[iblk+1];icrs++){
     const unsigned int jblk = crs.array[icrs];
     std::cout << iblk << " " << jblk << std::endl;
     }
     }*/
  }
	m_ls.AddMat_NonDia(ilss_rigid,ils_l, crs);
	Com::CIndexedArray crs_inv;
	crs_inv.SetTranspose( nblk_l, crs );
	m_ls.AddMat_NonDia(ils_l,ilss_rigid, crs_inv);
  return true;  
}



int Ls::CLinearSystem_RigidField2::FindIndexArray_Seg( 
        unsigned int id_field, 
        Fem::Field::ELSEG_TYPE type, 
        const Fem::Field::CFieldWorld& world )
{	
	if( !world.IsIdField(id_field) ) return -1;
    const Fem::Field::CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}
	for(unsigned int ils=0;ils<m_aSegRF.size();ils++){
		if( m_aSegRF[ils].id_field==id_field_parent && m_aSegRF[ils].node_config==type ){
			return ils;
		}
	}
	return -1;
}

bool Ls::CLinearSystem_RigidField2::SetFixedBoundaryCondition_Field
(unsigned int id_field, const Fem::Field::CFieldWorld& world )
{
  
	if( !world.IsIdField(id_field) ) return false;
  const Fem::Field::CField& field = world.GetField(id_field);
	unsigned int id_field_parent = field.GetIDFieldParent();
	if( id_field_parent == 0 ) id_field_parent = id_field;
  
	{
    int ils0 = this->FindIndexArray_Seg(id_field_parent,Fem::Field::CORNER,world);
		if( ils0 >= 0 && ils0 < m_aSegRF.size() ){
			const CLinSysSegRF& ls0 = this->m_aSegRF[ils0];
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
			if( ls0.id_field== id_field_parent ){
        Fem::Ls::BoundaryCondition(id_field,Fem::Field::CORNER,bc_flag,world);
			}
		}
	}
	{
    int ils0 = this->FindIndexArray_Seg(id_field_parent,Fem::Field::EDGE,world);
		if( ils0 >= 0 && ils0 < m_aSegRF.size() ){
      assert(0);
		}
	}
	{
    int ils0 = this->FindIndexArray_Seg(id_field_parent,Fem::Field::BUBBLE,world);
		if( ils0 >= 0 && ils0 < m_aSegRF.size() ){
      assert(0);
		}
	}
	return true;
}

bool Ls::CLinearSystem_RigidField2::UpdateValueOfField_NewmarkBeta(
		double gamma, double beta, double dt, 
		unsigned int id_field_val, Fem::Field::CFieldWorld& world, bool IsInitial )
{

	if( !world.IsIdField(id_field_val) ) return false;
    const Fem::Field::CField& field_val = world.GetField(id_field_val);

    if( field_val.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va != 0 ){	// 角節点を更新
		unsigned int id_na_val=0, id_ns_a,id_ns_v,id_ns_u;
		{
            const Fem::Field::CField::CNodeSegInNodeAry& nsna = field_val.GetNodeSegInNodeAry(Fem::Field::CORNER);
			id_na_val = nsna.id_na_va;
			id_ns_u = nsna.id_ns_va;
			id_ns_v = nsna.id_ns_ve;
			id_ns_a = nsna.id_ns_ac;
		}
        const int ilss = this->FindIndexArray_Seg(id_field_val,Fem::Field::CORNER,world);
		assert( ilss != -1 );
		const CLinSysSegRF& lss = this->m_aSegRF[ilss];
        const MatVec::CVector_Blk& upd = m_ls.GetVector(-2,ilss);
        Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val);
		// 値を更新する
		if( IsInitial ){
			na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v, dt);
			na_val.AddValueToNodeSegment(id_ns_u ,id_ns_a, 0.5*dt*dt);
		}
		na_val.AddValueToNodeSegment(id_ns_u, upd, beta*dt*dt);
		// 速度を更新する．
		if( IsInitial ){ 
		    na_val.AddValueToNodeSegment(id_ns_v ,id_ns_a, dt); 
		}
		na_val.AddValueToNodeSegment(id_ns_v, upd, gamma*dt);
		// 加速度を更新する
		na_val.AddValueToNodeSegment(id_ns_a, upd, 1.0);
	}
    if( field_val.GetNodeSegInNodeAry(Fem::Field::BUBBLE).id_na_va != 0 ){	// バブル節点を更新
		unsigned int id_na_val=0, id_ns_a,id_ns_v,id_ns_u;
		{
            const Fem::Field::CField::CNodeSegInNodeAry& nsna = field_val.GetNodeSegInNodeAry(Fem::Field::BUBBLE);
			id_na_val = nsna.id_na_va;
			id_ns_u = nsna.id_ns_va;
			id_ns_v = nsna.id_ns_ve;
			id_ns_a = nsna.id_ns_ac;
		}
        const int ilss = this->FindIndexArray_Seg(id_field_val,Fem::Field::BUBBLE,world);
		assert( ilss != -1 );
		const CLinSysSegRF& lss = this->m_aSegRF[ilss];
        const MatVec::CVector_Blk& upd = m_ls.GetVector(-2,ilss);
        Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val);
		// 値を更新する
		if( IsInitial ){
			na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v, dt);
			na_val.AddValueToNodeSegment(id_ns_u ,id_ns_a, 0.5*dt*dt);
		}
		na_val.AddValueToNodeSegment(id_ns_u, upd, beta*dt*dt);
		// 速度を更新する．
		if( IsInitial ){ 
		    na_val.AddValueToNodeSegment(id_ns_v ,id_ns_a, dt); 
		}
		na_val.AddValueToNodeSegment(id_ns_v, upd, gamma*dt);
		// 加速度を更新する
		na_val.AddValueToNodeSegment(id_ns_a, upd, 1.0);
    }
    if( field_val.GetNodeSegInNodeAry(Fem::Field::EDGE).id_na_va != 0 ){	// 角節点を更新
        assert(0);
    }
	return true;
}

////////////////////////////////////////////////////////////////



void Ls::CLinearSystem_RigidField2::SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB, 
                                                   const std::vector<Rigid::CConstraint*>& aConst)
{   
  //    std::cout << "0ls : set rigid system" << std::endl;
  assert( m_ls.GetNLinSysSeg() == this->m_aSegRF.size() );
  ilss_rigid = m_ls.GetNLinSysSeg();
  const unsigned int nRB = aRB.size();
  const unsigned int nConst = aConst.size();
  m_aSegRF.push_back( CLinSysSegRF(nRB,nConst) );
  
  const unsigned int nblk = nRB+nConst;
  {
    std::vector<unsigned int> aBlkSize;
    aBlkSize.resize(nblk);
    for(unsigned int irb=0;irb<nRB;irb++){
      aBlkSize[irb] = aRB[irb].GetDOF();
    }
    for(unsigned int icst=0;icst<nConst;icst++){
      aBlkSize[nRB+icst] = aConst[icst]->GetDOF();
    }
    const int ires = m_ls.AddLinSysSeg(nblk,aBlkSize);
    assert( ires == ilss_rigid );
  }
  Com::CIndexedArray crs;
  {
    crs.index.clear();
    crs.index.resize(nblk+1,0);
    crs.index[0] = 0;
    for(unsigned int i=0;i<nRB;i++){ crs.index[i+1]=0; }
    for(unsigned int i=nRB;i<nRB+nConst;i++){ crs.index[i+1]=0; }
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        crs.index[irb0+1] += 1;
      }
      crs.index[icst+nRB+1] += aIndRB.size();
    }
    for(unsigned int i=0;i<nblk;i++){ 
      crs.index[i+1] = crs.index[i+1] + crs.index[i];
    }
    const unsigned int ncrs = crs.index[nblk];
    crs.array.resize(ncrs);
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        const unsigned int icrs0 = crs.index[icst+nRB];
        crs.array[icrs0] = irb0;
        crs.index[icst+nRB] += 1;
        const unsigned int icrs1 = crs.index[irb0];
        crs.array[icrs1] = icst+nRB;
        crs.index[irb0] += 1;
      }
    }
    for(int i=nblk;i>0;i--){ 
      crs.index[i] = crs.index[i-1];
    }
    crs.index[0] = 0;
    crs.Sort();
    assert( crs.index[nRB+nConst] == ncrs );
  }
  m_ls.AddMat_Dia(ilss_rigid,crs);  
}


bool Ls::CLinearSystem_RigidField2::UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta,         bool is_first) const
{
    assert( ilss_rigid < m_ls.GetNLinSysSeg() );
    const MatVec::CVector_Blk& upd = m_ls.GetVector(-2,ilss_rigid);
    const unsigned int nRB    = this->GetSizeRigidBody();
    const unsigned int nConst = this->GetSizeConstraint();
    assert( nRB    == aRB.size()    );
    assert( nConst == aConst.size() );
    assert( upd.NBlk() == nRB + nConst );
    double tmp[6];
    for(unsigned int irb=0;irb<aRB.size();irb++){
        const unsigned int nlen = upd.Len(irb);
        assert( nlen <= 6 );
        for(unsigned int ilen=0;ilen<nlen;ilen++){
            tmp[ilen] = upd.GetValue(irb,ilen);
        }
        aRB[irb].UpdateSolution( tmp, 
            dt, newmark_gamma, newmark_beta, 
            is_first);
    }
    for(unsigned int icst=0;icst<aConst.size();icst++){
        const unsigned int nlen = upd.Len(icst+nRB);
        assert( nlen <= 6 );
        for(unsigned int ilen=0;ilen<nlen;ilen++){
            tmp[ilen] = upd.GetValue(icst+nRB,ilen);
        }
        aConst[icst]->UpdateSolution( tmp, 
            dt, newmark_gamma, newmark_beta );
    }
    return true;
    return true;
}

void Ls::CLinearSystem_RigidField2::AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
                                                const Com::CVector3D& vres, double d )
{
    MatVec::CVector_Blk& res = m_ls.GetVector(-1,ilss_rigid);
    const unsigned int nRB    = this->GetSizeRigidBody();
    unsigned int iblk = (is_rb) ? ind : ind+nRB;
    assert( offset+3 <= res.Len(iblk) );
    res.AddValue(iblk,offset+0,vres.x*d);
    res.AddValue(iblk,offset+1,vres.y*d);
    res.AddValue(iblk,offset+2,vres.z*d);
}

void Ls::CLinearSystem_RigidField2::AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
                                                const double* eres, double d )
{
    MatVec::CVector_Blk& res = m_ls.GetVector(-1,ilss_rigid);
    const unsigned int nRB    = this->GetSizeRigidBody();
    unsigned int iblk = (is_rb) ? ind : ind+nRB;
    assert( offset+size <= res.Len(iblk) );
    for(unsigned int i=0;i<size;i++){
        res.AddValue(iblk,offset+i,eres[i]*d);
    }
}

void Ls::CLinearSystem_RigidField2::SubResidual(const unsigned int ind, bool is_rb, const double* vec)
{
    MatVec::CVector_Blk& res = m_ls.GetVector(-1,ilss_rigid);
    const unsigned int nRB    = this->GetSizeRigidBody();
    unsigned int iblk = (is_rb) ? ind : ind+nRB;
    unsigned int len = res.Len(iblk);
    for(unsigned int i=0;i<len;i++){
        res.AddValue(iblk,i,vec[i]*-1);
    }
}

void Ls::CLinearSystem_RigidField2::AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans )
{
    assert( m_ls.IsMatrix(ilss_rigid,ilss_rigid) );
    MatVec::CMatDia_BlkCrs& mat = m_ls.GetMatrix(ilss_rigid);

    const unsigned int nRB    = this->GetSizeRigidBody();
    const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
    const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
    const unsigned int lencol = mat.LenBlkCol(iblk0);
    const unsigned int lenrow = mat.LenBlkRow(iblk1);

    assert( lencol*lenrow <= 36 );
    double tmp[36];
    for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
    if( isnt_trans ){
        for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[i*3+j]*d;
        }
        }
    }
    else{
        for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[j*3+i]*d;
        }
        }
    }
    mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
}

void Ls::CLinearSystem_RigidField2::AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_column)
{
    assert( m_ls.IsMatrix(ilss_rigid,ilss_rigid) );
    MatVec::CMatDia_BlkCrs& mat = m_ls.GetMatrix(ilss_rigid);

    const unsigned int nRB    = this->GetSizeRigidBody();
    const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
    const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
    const unsigned int lencol = mat.LenBlkCol(iblk0);
    const unsigned int lenrow = mat.LenBlkRow(iblk1);
    assert( lencol*lenrow <= 36 );

    double tmp[36];
    for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
    if( is_column ){
        tmp[ (0+offsetr)*lenrow + (offsetl)] += vec.x*d;
        tmp[ (1+offsetr)*lenrow + (offsetl)] += vec.y*d;
        tmp[ (2+offsetr)*lenrow + (offsetl)] += vec.z*d;
    }
    else{
        tmp[ (offsetr)*lenrow + (0+offsetl)] += vec.x*d;
        tmp[ (offsetr)*lenrow + (1+offsetl)] += vec.y*d;
        tmp[ (offsetr)*lenrow + (2+offsetl)] += vec.z*d;
    }
    mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
}

void Ls::CLinearSystem_RigidField2::AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
               unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
               const double* emat, double d)
{
  assert( m_ls.IsMatrix(ilss_rigid,ilss_rigid) );
  MatVec::CMatDia_BlkCrs& mat = m_ls.GetMatrix(ilss_rigid);
  const unsigned int nRB    = this->GetSizeRigidBody();
  const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
  const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
  const unsigned int lencol = mat.LenBlkCol(iblk0);
  const unsigned int lenrow = mat.LenBlkRow(iblk1);
  assert( lencol*lenrow <= 36 );
  double tmp[36];
  for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
  for(unsigned int i=0;i<sizer;i++){
  for(unsigned int j=0;j<sizel;j++){
    tmp[(i+offsetr)*lenrow + (j+offsetl)] += emat[i*sizel+j]*d;
  }
  }  
  mat.Mearge(1,&iblk0, 1,&iblk1, lencol*lenrow, tmp );
}
    