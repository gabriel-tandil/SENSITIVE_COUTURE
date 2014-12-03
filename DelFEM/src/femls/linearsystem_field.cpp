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

////////////////////////////////////////////////////////////////
// LinearSystem.cpp : 連立一次方程式クラス(CLinearSystem.h)の実装
////////////////////////////////////////////////////////////////


#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#ifndef for 
#define for if(0); else for
#endif

#include "math.h"

#include "delfem/field_world.h"
#include "delfem/femls/linearsystem_field.h"

#include "delfem/indexed_array.h"
#include "delfem/vector3d.h"
#include "delfem/quaternion.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/bcflag_blk.h"

using namespace MatVec;
using namespace Fem::Ls;
using namespace Fem::Field;



// elem_aryに含まれる節点全てidofblk番目の自由度を固定する
static bool SetBCFlagToES
(unsigned int id_field,
 MatVec::CBCFlag& bc_flag, 
 const Fem::Field::CElemAry& ea, unsigned int id_es, unsigned int idofblk)
{
  assert( (int)idofblk < bc_flag.LenBlk() );
  if( (int)idofblk >= bc_flag.LenBlk() ) return false;
  assert( ea.IsSegID(id_es) );
  const Fem::Field::CElemAry::CElemSeg& es = ea.GetSeg(id_es);
  unsigned int noes[256];
  unsigned int nnoes = es.Length();
  for(unsigned int ielem=0;ielem<ea.Size();ielem++){
    es.GetNodes(ielem,noes);
    for(unsigned int inoes=0;inoes<nnoes;inoes++){
      const unsigned int inode0 = noes[inoes];
      bc_flag.SetBC(inode0,idofblk);
      //			m_Flag[inode0*m_lenBlk+idofblk] = 1;
    }
  }
  return true;
}


void Fem::Ls::BoundaryCondition
(unsigned int id_field, const ELSEG_TYPE& elseg_type, unsigned int idofns, 
 MatVec::CBCFlag& bc_flag, const CFieldWorld& world)
{
  if( !world.IsIdField(id_field) ) return;
  const Fem::Field::CField& field = world.GetField(id_field);
  assert( (int)idofns < bc_flag.LenBlk()  );
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
  for(unsigned int iea=0;iea<aIdEA.size();iea++){
    const unsigned int id_ea = aIdEA[iea];
    const CElemAry& ea = world.GetEA(id_ea);
    if( elseg_type == CORNER && field.GetIdElemSeg(id_ea,CORNER,true,world) != 0 ){
      SetBCFlagToES(id_field,bc_flag, ea, field.GetIdElemSeg(id_ea,CORNER,true,world), idofns);
    }
    if( elseg_type == BUBBLE && field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){
      SetBCFlagToES(id_field,bc_flag, ea, field.GetIdElemSeg(id_ea,BUBBLE,true,world), idofns);
    }
    if( elseg_type == EDGE   && field.GetIdElemSeg(id_ea,EDGE,  true,world) != 0 ){
      SetBCFlagToES(id_field,bc_flag, ea, field.GetIdElemSeg(id_ea,EDGE,  true,world), idofns);
    }
  }
}


void Fem::Ls::BoundaryCondition
(unsigned int id_field, const ELSEG_TYPE& elseg_type,  
 MatVec::CBCFlag& bc_flag, const CFieldWorld& world,  
 unsigned int ioffset)
{
  if( !world.IsIdField(id_field) ) return;
  const Fem::Field::CField& field = world.GetField(id_field);  
  {	// Assert
    const unsigned int len = bc_flag.LenBlk();
    assert( ioffset < len );
  }
  const unsigned int nlen = field.GetNLenValue();
  
  const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
  for(unsigned int iea=0;iea<aIdEA.size();iea++){
    unsigned int id_ea = aIdEA[iea];
    const CElemAry& ea = world.GetEA(id_ea);
    unsigned int noes[256];
    if( elseg_type == CORNER && field.GetIdElemSeg(id_ea,CORNER,true,world) != 0 ){
      const Fem::Field::CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,CORNER,true,world);
      unsigned int nnoes = es.Length();
      for(unsigned int ielem=0;ielem<ea.Size();ielem++){
        es.GetNodes(ielem,noes);
        for(unsigned int inoes=0;inoes<nnoes;inoes++){
          for(unsigned int ilen=0;ilen<nlen;ilen++){
            bc_flag.SetBC(noes[inoes],ilen+ioffset);
          }
        }
      }
    }
    if( elseg_type == BUBBLE && field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){
      const Fem::Field::CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,BUBBLE,true,world);
      unsigned int nnoes = es.Length();
      for(unsigned int ielem=0;ielem<ea.Size();ielem++){
        es.GetNodes(ielem,noes);
        for(unsigned int inoes=0;inoes<nnoes;inoes++){
          for(unsigned int ilen=0;ilen<nlen;ilen++){
            bc_flag.SetBC(noes[inoes],ilen+ioffset);
          }
        }
      }
    }
    if( elseg_type == EDGE   && field.GetIdElemSeg(id_ea,EDGE,true,world) != 0 ){
      const Fem::Field::CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,EDGE,true,world);
      unsigned int nnoes = es.Length();
      for(unsigned int ielem=0;ielem<ea.Size();ielem++){
        es.GetNodes(ielem,noes);
        for(unsigned int inoes=0;inoes<nnoes;inoes++){
          for(unsigned int ilen=0;ilen<nlen;ilen++){
            bc_flag.SetBC(noes[inoes],ilen+ioffset);
          }
        }
      }
    }
  }
}



//////////////////////////////////////


int CLinearSystem_Field::AddLinSysSeg_Field(
		const unsigned int id_field, ELSEG_TYPE es_type, const CFieldWorld& world)
{
	{
		int iseg = this->FindIndexArray_Seg(id_field,es_type,world);
		if( iseg != -1 ) return iseg;
	}
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}
	const unsigned int nlen_value = field.GetNLenValue();
	////////////////
	unsigned int id_na_val = field.GetNodeSegInNodeAry(es_type).id_na_va;
	if( id_na_val == 0 ) return -1;	
	////////////////
	assert( world.IsIdNA(id_na_val) );
	const CNodeAry& na = world.GetNA(id_na_val);
	CLinSysSeg_Field seg;
	{
		seg.id_field = id_field_parent; 
		seg.id_field2 = 0;
		seg.node_config = es_type;
		seg.len=field.GetNLenValue(); 
		seg.nnode=na.Size();
	}
	return this->AddLinSysSeg(seg);
}

// add pattern into diagonal sub matrix
bool CLinearSystem_Field::AddPattern_Field(const unsigned int id_field, const CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	const int ils_b = AddLinSysSeg_Field(id_field,BUBBLE,world);
	const int ils_e = AddLinSysSeg_Field(id_field,EDGE,  world);
	const int ils_c = AddLinSysSeg_Field(id_field,CORNER,world);
//	std::cout << "AddPattern_Field : " << ils_c << " " << ils_b << " " << ils_e << std::endl;
	////////////////////////////////
	const std::vector<unsigned int> aIdEA = field.GetAryIdEA();
  assert( aIdEA.size() > 0 );
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
	{
		const unsigned int id_ea = aIdEA[iiea];
		const CElemAry& ea = world.GetEA(id_ea);
		// CORNER節点について
		if( field.GetIdElemSeg(id_ea,CORNER,true,world) != 0 ){
			assert( world.IsIdEA(id_ea) );
			const unsigned int id_es_c = field.GetIdElemSeg(id_ea,CORNER,true,world);
			assert( ea.IsSegID(id_es_c) );
			const unsigned int ils0 = ils_c;
      {
        Com::CIndexedArray crs;
        ea.MakePattern_FEM(id_es_c,crs);
        this->AddMat_Dia(ils0, crs );			// cc行列を作る        
      }
			if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){	// CORNER-BUBBLE
				const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);	assert( ea.IsSegID(id_es_b) );
				const unsigned int ils1 = ils_b;
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_c,id_es_b,crs);	assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// cb行列を作る
				const unsigned int nnode1 = m_aSegField[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// bc行列を作る
			}
			if( field.GetIdElemSeg(id_ea,EDGE,true,world) != 0 ){	// CONRER-EDGE
				const unsigned int id_es_e = field.GetIdElemSeg(id_ea,EDGE,true,world);		assert( ea.IsSegID(id_es_e) );
				const unsigned int ils1 = ils_e;
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_c,id_es_e,crs);	assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// ce行列を作る
				const unsigned int nnode1 = m_aSegField[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// ec行列を作る
			}
		}
		// EDGE節点について
		if( field.GetIdElemSeg(id_ea,EDGE,true,world) != 0 ){
			const unsigned int id_es_e = field.GetIdElemSeg(id_ea,EDGE,true,world);
			assert( ea.IsSegID(id_es_e) );
			const unsigned int ils0 = ils_e;
      {
        Com::CIndexedArray crs;
        ea.MakePattern_FEM(id_es_e,crs);
        this->AddMat_Dia(ils0, crs );			// cc行列を作る
      }
			if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){	// EDGE-BUBBLE
				const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);	assert( ea.IsSegID(id_es_b) );
				const unsigned int ils1 = ils_b;
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_e,id_es_b,crs);	assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// eb行列を作る
				const unsigned int nnode1 = m_aSegField[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// be行列を作る
			}
		}
		// BUBBLE節点について
		if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){
			const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
			assert( ea.IsSegID(id_es_b) );
			const unsigned int ils0 = ils_b;			
      {
        Com::CIndexedArray crs;
        ea.MakePattern_FEM(id_es_b,crs);
        this->AddMat_Dia(ils0, crs );			// bb行列を作る
      }
		}
	}
	return true;
}

// fieldとfield2がパターンが同じだとして，ブロックが結合された一つの行列を作る
bool CLinearSystem_Field::AddPattern_CombinedField(unsigned id_field1, unsigned int id_field2, const Fem::Field::CFieldWorld& world)
{
	std::cout << "AddPattern Combined    (CLinearSystem)" << std::endl;
	if( !world.IsIdField(id_field1) ) return false;
	const CField& field1 = world.GetField(id_field1);
	unsigned int id_field1_parent;
	{
		if( field1.GetIDFieldParent() == 0 ){ id_field1_parent = id_field1; }
		else{ id_field1_parent = field1.GetIDFieldParent(); }
	}
	////////////////
	if( !world.IsIdField(id_field2) ) return false;
	const CField& field2 = world.GetField(id_field2);
	unsigned int id_field2_parent;
	{
		if( field2.GetIDFieldParent() == 0 ){ id_field2_parent = id_field2; }
		else{ id_field2_parent = field2.GetIDFieldParent(); }
	}

	unsigned int nnode;
	{
		const Fem::Field::CField::CNodeSegInNodeAry& nsna1 = field1.GetNodeSegInNodeAry(CORNER);
		const Fem::Field::CField::CNodeSegInNodeAry& nsna2 = field2.GetNodeSegInNodeAry(CORNER);
		if( nsna1.id_na_va != nsna2.id_na_va ){
			assert(0);
			return false;
		}
		const unsigned int id_na = nsna1.id_na_va;
		assert( world.IsIdNA(id_na) );
		const CNodeAry& na = world.GetNA(id_na);
		nnode = na.Size();
	}
	if(   !field1.IsNodeSeg(CORNER,true,world) 
		|| field1.IsNodeSeg(BUBBLE,true,world) 
		|| field1.IsNodeSeg(EDGE,  true,world) ){
		std::cout << "corner " << field1.GetNodeSegInNodeAry(CORNER).id_na_va << " " << field1.IsNodeSeg(CORNER,true,world) << std::endl;
		std::cout << "bubble " << field1.GetNodeSegInNodeAry(BUBBLE).id_na_va << std::endl;
		std::cout << "edge   " << field1.GetNodeSegInNodeAry(EDGE  ).id_na_va << std::endl;
		std::cout << "Error!-->Not Implimented!" << std::endl;
		assert(0);
		return false;
	}
	////////////////////////////////////////////////////////////////
	const std::vector<unsigned int>& aIdEA1 = field1.GetAryIdEA();
	const std::vector<unsigned int>& aIdEA2 = field2.GetAryIdEA();
	if( aIdEA1.size() != aIdEA2.size() ){
		assert(0);
		return false;
	}
	////////////////
	const unsigned int nIdEA = aIdEA1.size();
	{
		for(unsigned int iIdEA=0;iIdEA<nIdEA;iIdEA++){
			if( aIdEA1[iIdEA] != aIdEA2[iIdEA] ){	
				assert(0); return false;
			}
			const unsigned int id_ea = aIdEA1[iIdEA];
			if( field1.GetInterpolationType(id_ea,world) != field2.GetInterpolationType(id_ea,world) ){
				assert(0); return false;
			}
			const unsigned int id_es1 = field1.GetIdElemSeg(id_ea,CORNER,true,world);
			const unsigned int id_es2 = field2.GetIdElemSeg(id_ea,CORNER,true,world);
			if( id_es1 != id_es2 ){
				assert(0); return false;
			}
		}
	}

	CLinSysSeg_Field lss;
	{
		lss.id_field = id_field1; lss.id_field2 = id_field2;
		lss.len = field1.GetNLenValue() + field2.GetNLenValue();
		lss.nnode = nnode;
		lss.node_config = CORNER;
	}
	const int ilss0 = this->AddLinSysSeg(lss);
	assert( ilss0 >= 0 );
	for(unsigned int iIdEA=0;iIdEA<nIdEA;iIdEA++){
		const unsigned int id_ea = aIdEA1[iIdEA];
		const unsigned int id_es = field1.GetIdElemSeg(id_ea,CORNER,true,world);
		const CElemAry& ea = world.GetEA(id_ea);
		assert( ea.IsSegID(id_es) );
        {
            Com::CIndexedArray crs;
            ea.MakePattern_FEM(id_es,crs);
		    this->AddMat_Dia(ilss0, crs );			// cc行列を作る
        }
	}

	return true;
}

// fieldで初期化する、fieldの中の非ゼロパターンで行列を作る
bool CLinearSystem_Field::AddPattern_Field(
        unsigned int id_field1, 
        unsigned int id_field2, 
        const CFieldWorld& world)
{
	if( !world.IsIdField(id_field1) ) return false;
	const CField& field1 = world.GetField(id_field1);
	unsigned int id_field_parent;
	{
		if( field1.GetIDFieldParent() == 0 ){ id_field_parent = id_field1; }
		else{ id_field_parent = field1.GetIDFieldParent(); }
	}

	const unsigned int nlen_value = field1.GetNLenValue();

	const int ils_c = this->AddLinSysSeg_Field(id_field1,CORNER,world);
	const int ils_b = this->AddLinSysSeg_Field(id_field1,BUBBLE,world);

	const CField& field2 = world.GetField(id_field2);
	const std::vector<unsigned int>& aIdEA1 = field1.GetAryIdEA();

	if( field1.GetNodeSegInNodeAry(CORNER).id_na_va != 0 &&
		field1.GetNodeSegInNodeAry(CORNER).id_na_co == 0 &&
		field1.GetNodeSegInNodeAry(BUBBLE).id_na_co != 0 &&
		field1.GetNodeSegInNodeAry(BUBBLE).id_na_va == 0 )
	{
		// Lagrange未定乗数のための接続（これはどこかに移す可能性が大きい）
		const unsigned int id_na_co1 = field1.GetNodeSegInNodeAry(BUBBLE).id_na_co;
		const unsigned int id_na_va2 = field1.GetNodeSegInNodeAry(CORNER).id_na_va;
			
		const unsigned int id_ea1 = aIdEA1[0];
		assert( field1.GetIdElemSeg(id_ea1,CORNER,true,world) != 0 );
		assert( ils_c != -1 );

		assert( world.IsIdEA(id_ea1) );
		const CElemAry& ea1 = world.GetEA(id_ea1);
		const unsigned int id_es_c1 = field1.GetIdElemSeg(id_ea1,CORNER,true,world);
		assert( ea1.IsSegID(id_es_c1) );
		// CORNER1-CORNER1
		{
			Com::CIndexedArray crs;
      ea1.MakePattern_FEM(id_es_c1,crs);
			this->AddMat_Dia(ils_c, crs );			// cc行列を作る
		}
		const int ils2_c = this->FindIndexArray_Seg(id_field2,CORNER,world);
		assert( ils2_c >= 0 && ils2_c < (int)this->GetNLynSysSeg() );
		const unsigned int id_es_b1 = field1.GetIdElemSeg(id_ea1,BUBBLE,false,world);
		Com::CIndexedArray crs;
		ea1.MakePattern_FEM(id_es_c1,id_es_b1,crs);
		assert( crs.CheckValid() );
		this->AddMat_NonDia(ils_c,ils2_c, crs);		// c1c2行列を足す
		unsigned int nnode2 = m_aSegField[ils2_c].nnode;
		Com::CIndexedArray crs_inv;
		crs_inv.SetTranspose( nnode2, crs );
		this->AddMat_NonDia(ils2_c,ils_c, crs_inv);	// c2c1行列を足す
		return true;
	}

	const std::vector<unsigned int>& aIdEA2 = field2.GetAryIdEA();

	for(;;){	// ダミーのfor文を使ってbreakで抜けられるようにする
		// Corner-Corner関係を作る
		if( ils_c == -1 ) break;
//		const unsigned int id_na_va1 = field1.GetNodeSegInNodeAry(CORNER).id_na_va;
		const unsigned int id_na_co1 = field1.GetNodeSegInNodeAry(CORNER).id_na_co;
//		const unsigned int id_na_va2 = field2.GetNodeSegInNodeAry(CORNER).id_na_va;
		const unsigned int id_na_co2 = field2.GetNodeSegInNodeAry(CORNER).id_na_co;
        assert( aIdEA1.size() > 0 );
        if( id_na_co1 != id_na_co2 ) break;

		for(unsigned int iiea1=0;iiea1<aIdEA1.size();iiea1++){
			const unsigned int id_ea1 = aIdEA1[iiea1];
			if( field1.GetIdElemSeg(id_ea1,CORNER,true,world) == 0 ) continue;
			assert( ils_c != -1 );

			assert( world.IsIdEA(id_ea1) );
			const CElemAry& ea1 = world.GetEA(id_ea1);
			const unsigned int id_es_c1 = field1.GetIdElemSeg(id_ea1,CORNER,true,world);
			assert( ea1.IsSegID(id_es_c1) );
			// CORNER1-CORNER1
            {
                Com::CIndexedArray crs;
                ea1.MakePattern_FEM(id_es_c1,crs);
		        this->AddMat_Dia(ils_c, crs );			// cc行列を作る
            }
			////////////////
			for(unsigned int iiea2=0;iiea2<aIdEA2.size();iiea2++)
			{
				const unsigned int id_ea2 = aIdEA2[iiea2];
				if( id_ea1 == id_ea2 ){
					// CORNER1-CORNER2
					const int ils2_c = this->FindIndexArray_Seg(id_field2,CORNER,world);
					if( ils2_c != -1 ){
                        assert( ils2_c >= 0 && ils2_c < (int)this->GetNLynSysSeg() );
						const unsigned int id_es_c2 = field2.GetIdElemSeg(id_ea2,CORNER,true,world);
						Com::CIndexedArray crs;
						ea1.MakePattern_FEM(id_es_c1,id_es_c2,crs);
						assert( crs.CheckValid() );
						this->AddMat_NonDia(ils_c,ils2_c, crs);		// c1c2行列を足す
						unsigned int nnode2 = m_aSegField[ils2_c].nnode;
						Com::CIndexedArray crs_inv;
						crs_inv.SetTranspose( nnode2, crs );
						this->AddMat_NonDia(ils2_c,ils_c, crs_inv);	// c2c1行列を足す
					}
					// CORNER1-BUBBLE2
					const int ils2_b = this->FindIndexArray_Seg(id_field2,BUBBLE,world);
					if( ils2_b != -1 ){
                        assert( ils2_b >= 0 && ils2_b < (int)this->GetNLynSysSeg() );
						const unsigned int id_es_b2 = field2.GetIdElemSeg(id_ea2,BUBBLE,true,world);
						assert( id_es_c1 != id_es_b2 );
						Com::CIndexedArray crs;
						ea1.MakePattern_FEM(id_es_c1,id_es_b2,crs);
						assert( crs.CheckValid() );
						this->AddMat_NonDia(ils_c,ils2_b, crs);		// c1b2行列を足す
						unsigned int nnode2 = m_aSegField[ils2_b].nnode;
						Com::CIndexedArray crs_inv;
						crs_inv.SetTranspose( nnode2, crs );
						this->AddMat_NonDia(ils2_b,ils_c, crs_inv);	// b2c1行列を足す
					}
				}
				else{	// 含まれる場合
					const CNodeAry& na1 = world.GetNA(id_na_co1);
					const unsigned int id_es_c_co1 = field1.GetIdElemSeg(id_ea1,CORNER,false,world);
					const unsigned int id_es_c_co2 = field2.GetIdElemSeg(id_ea2,CORNER,false,world);
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
						int ils2 = this->FindIndexArray_Seg(id_field2,CORNER,world);
            assert( ils2 >= 0 && ils2 < (int)this->GetNLynSysSeg() );
//						std::cout << "ils_seg : " << ils0 << " " << ils2 << std::endl;
						this->AddMat_NonDia(ils_c,ils2, crs);
						unsigned int nnode2 = m_aSegField[ils2].nnode;
						Com::CIndexedArray crs_inv;
						crs_inv.SetTranspose( nnode2, crs );
						this->AddMat_NonDia(ils2,ils_c, crs_inv);
					}
				}
			}
		}
		break;
	}
			
	for(unsigned int iiea=0;iiea<aIdEA1.size();iiea++){
		const unsigned int id_ea1 = aIdEA1[iiea];
		const unsigned int id_es_b1 = field1.GetIdElemSeg(id_ea1,BUBBLE,true,world);
		if( id_es_b1 == 0 ) continue;
		const CElemAry& ea1 = world.GetEA(id_ea1);
		assert( world.IsIdEA(id_ea1) );
		assert( ea1.IsSegID(id_es_b1) );
		// BUBLE1-BUBBLE2
        {
			Com::CIndexedArray crs;
			ea1.MakePattern_FEM(id_es_b1,crs);
		    this->AddMat_Dia(ils_b, crs );
        }
		const unsigned int id_ea2 = aIdEA2[iiea];
		assert( id_ea1 == id_ea2 );
		const unsigned int id_es_c2 = field2.GetIdElemSeg(id_ea2,CORNER,true,world);
		assert( ea1.IsSegID(id_es_c2) );
		int ils2 = this->FindIndexArray_Seg(id_field2,CORNER,world);
    assert( ils2 >= 0 && ils2 < (int)this->GetNLynSysSeg() );
		{
			Com::CIndexedArray crs;
			ea1.MakePattern_FEM(id_es_b1,id_es_c2,crs);
			assert( crs.CheckValid() );
			this->AddMat_NonDia(ils_b,ils2, crs);		// b1c2行列を作る
			Com::CIndexedArray crs_inv;
			const unsigned nnode2 = m_aSegField[ils2].nnode;
			crs_inv.SetTranspose( nnode2, crs );
			this->AddMat_NonDia(ils2,ils_b, crs_inv);	// c2b1行列を作る
		}
	}



	// いろいろと追加が必要
	// まずは足りない部分を要求するようになったらエラーを出す関数を実装しよう。

	return true;
}


// 残差ベクトルをゲットする
MatVec::CVector_Blk& CLinearSystem_Field::GetResidual
(unsigned int id_field, const ELSEG_TYPE elseg_type, 
 const CFieldWorld& world)
{
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
    if( ils0 < 0 || ils0 >= (int)this->GetNLynSysSeg() ){ assert(0); throw 0; }
	return m_ls.GetVector(-1,ils0);
}

// 更新ベクトルをゲットする
MatVec::CVector_Blk& CLinearSystem_Field::GetUpdate
(unsigned int id_field, const ELSEG_TYPE elseg_type, 
 const CFieldWorld& world)
{
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
  if( ils0 < 0 || ils0 >= (int)this->GetNLynSysSeg() ){ assert(0); throw 0; }
	return m_ls.GetVector(-2,ils0);
}

CMatDia_BlkCrs& CLinearSystem_Field::GetMatrix
(unsigned int id_field,
 const ELSEG_TYPE elseg_type, 
 const CFieldWorld& world)
{
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
  if( ils0 < 0 || ils0 >= (int)this->GetNLynSysSeg() ){ assert(0); throw 0; }
	return m_ls.GetMatrix(ils0);
}

CMat_BlkCrs& CLinearSystem_Field::GetMatrix
(unsigned int id_field_col,const ELSEG_TYPE elseg_type_col,
 unsigned int id_field_row,const ELSEG_TYPE elseg_type_row,
 const CFieldWorld& world)
{
	int ils_col = FindIndexArray_Seg(id_field_col,elseg_type_col,world);
    if( ils_col < 0 || ils_col >= (int)this->GetNLynSysSeg() ){ assert(0); throw 0; }
	int ils_row = FindIndexArray_Seg(id_field_row,elseg_type_row,world);
    if( ils_row < 0 || ils_row >= (int)this->GetNLynSysSeg() ){ assert(0); throw 0; }
	return m_ls.GetMatrix(ils_col,ils_row);
}

bool CLinearSystem_Field::SetFixedBoundaryCondition_Field( 
        unsigned int id_field, 
        unsigned int idofns, 
        const CFieldWorld& world )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent = field.GetIDFieldParent();
	if( id_field_parent == 0 ){ id_field_parent = id_field; }
	{
    const int ils0 = this->FindIndexArray_Seg(id_field_parent,CORNER,world);
    if( ils0 >=0 && ils0 < (int)this->GetNLynSysSeg() ){
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_BCFlag[ils0];
			const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
			if( ls0.id_field == id_field_parent ){
        Fem::Ls::BoundaryCondition(id_field,CORNER,idofns,bc_flag,world);
			}
			else{   // 結合された場
				assert( ls0.id_field2 == id_field_parent );
				unsigned int id_field1 = ls0.id_field;
				assert( world.IsIdField(id_field1) );
				const CField& field1 = world.GetField(id_field1);
				const unsigned int len1 = field1.GetNLenValue();
        Fem::Ls::BoundaryCondition(id_field,CORNER,idofns+len1,bc_flag,world);
			}
		}
	}
	{
    const int ils0 = this->FindIndexArray_Seg(id_field_parent,EDGE,world);
    if( ils0 >= 0 && ils0 < (int)m_aSegField.size() ){
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
			const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
			assert( ls0.id_field== id_field_parent );
      Fem::Ls::BoundaryCondition(id_field,EDGE,idofns,bc_flag,world);
		}
	}
	{
		const int ils0 = this->FindIndexArray_Seg(id_field_parent,BUBBLE,world);
    if( ils0 >= 0 && ils0 < (int)m_aSegField.size() ){
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
			const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
			assert( ls0.id_field== id_field_parent );
      Fem::Ls::BoundaryCondition(id_field,BUBBLE,idofns,bc_flag,world);
		}
	}
	return true;
}

bool CLinearSystem_Field::SetFixedBoundaryCondition_Field( unsigned int id_field, const CFieldWorld& world )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent = field.GetIDFieldParent();
  if( id_field_parent == 0 ) id_field_parent = id_field;
	{
    int ils0 = this->FindIndexArray_Seg(id_field,CORNER,world);
    if( ils0 >= 0 && ils0 < (int)m_aSegField.size() ){
			const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
			if( ls0.id_field == id_field_parent ){
				BoundaryCondition(id_field,CORNER,bc_flag,world);
			}
			else{   // combined field
				assert( ls0.id_field2 == id_field_parent );
				unsigned int id_field1 = ls0.id_field;
				assert( world.IsIdField(id_field1) );
				const CField& field1 = world.GetField(id_field1);
				const unsigned int len1 = field1.GetNLenValue();
        assert( bc_flag.LenBlk() == (int)len1 + (int)field.GetNLenValue() );
				BoundaryCondition(id_field,CORNER,bc_flag,world,len1);
			}
		}
	}
	{
		int ils0 = this->FindIndexArray_Seg(id_field_parent,EDGE,world);
    if( ils0 >= 0 && ils0 < (int)m_aSegField.size() ){            
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
			const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
			assert( ls0.id_field== id_field_parent );
			BoundaryCondition(id_field,EDGE,bc_flag,world);
		}
	}
	{
		int ils0 = this->FindIndexArray_Seg(id_field_parent,BUBBLE,world);
    if( ils0 >= 0 && ils0 < (int)m_aSegField.size() ){
      MatVec::CBCFlag& bc_flag = m_ls.GetBCFlag(ils0);//*m_ls.m_BCFlag[ils0];
      const CLinSysSeg_Field& ls0 = this->m_aSegField[ils0];
			assert( ls0.id_field== id_field_parent );
			BoundaryCondition(id_field,BUBBLE,bc_flag,world);
		}
	}
	return true;
}

void CLinearSystem_Field::ClearFixedBoundaryCondition(){
    m_ls.ClearFixedBoundaryCondition();
}

bool CLinearSystem_Field::UpdateValueOfField( 
	unsigned int id_field, Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// CORNERについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(CORNER);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
			int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
            assert( ilss0 >= 0 && ilss0 < (int)this->GetNLynSysSeg() );
            const CLinSysSeg_Field& lss0 = this->m_aSegField[ilss0];
			const CVector_Blk* m_Update = m_ls.m_Update[ilss0];
			assert( m_Update != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			if( lss0.id_field == id_field ){
				na.AddValueToNodeSegment(id_ns_val,*m_Update,1.0);
			}
			else{	// Combined Segmentに対するアップデート
				assert( lss0.id_field2 == id_field );
				const Fem::Field::CField& field1 = world.GetField( lss0.id_field );
				const unsigned int nlen1 = field1.GetNLenValue();
				na.AddValueToNodeSegment(id_ns_val,*m_Update,1.0,nlen1);
				{	// assertルーティン
					const Fem::Field::CField& field2 = world.GetField( lss0.id_field2 );
					const unsigned int nlen2 = field2.GetNLenValue();
                    assert( (int)(nlen1 + nlen2) == m_Update->Len() );
				}
			}
		}
	}
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// BUBBLEについて更新
			const CField::CNodeSegInNodeAry& nsna_b = field.GetNodeSegInNodeAry(BUBBLE);
			id_na_val = nsna_b.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_b.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_b.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_b.id_ns_ac;
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,BUBBLE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			na.AddValueToNodeSegment(id_ns_val,upd,1.0);
		} 
    }
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// EDGEについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(EDGE);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,EDGE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			na.AddValueToNodeSegment(id_ns_val,upd,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
//			ns.is_updated = true;
		}
	}
	return true;
}


bool CLinearSystem_Field::UpdateValueOfField_RotCRV( 
	unsigned int id_field, Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// CORNERについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(CORNER);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
			int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
            assert( ilss0 >= 0 && ilss0 < (int)this->GetNLynSysSeg() );
			const CLinSysSeg_Field& lss0 = this->m_aSegField[ilss0];
            MatVec::CVector_Blk* m_Update = m_ls.m_Update[ilss0];
			assert( m_Update != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
			const unsigned int nblk = m_Update->NBlk();
			assert( ns.Size() == nblk );
			assert( ns.Length() == 3 );
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				Com::CMatrix3 rot0;
				{
					double rot0a[3];
					ns.GetValue(iblk,rot0a);
                    rot0.SetRotMatrix_CRV(rot0a);
				}
				Com::CMatrix3 rot1;
				if( lss0.id_field == id_field ){
					double rot1a[3];
					rot1a[0] = m_Update->GetValue(iblk,0);
					rot1a[1] = m_Update->GetValue(iblk,1);
					rot1a[2] = m_Update->GetValue(iblk,2);
                    rot1.SetRotMatrix_Rodrigues(rot1a);
				}
				else{
					assert( lss0.id_field2 == id_field );
					assert( m_Update->Len() == 6 );
					double rot1a[3];
					rot1a[0] = m_Update->GetValue(iblk,3);
					rot1a[1] = m_Update->GetValue(iblk,4);
					rot1a[2] = m_Update->GetValue(iblk,5);
                    rot1.SetRotMatrix_Rodrigues(rot1a);
				}
                Com::CMatrix3 rot10 = rot1.MatMat(rot0);
                double crv10[3];
                rot10.GetCRV_RotMatrix(crv10);
				ns.SetValue(iblk,0,crv10[0]);
				ns.SetValue(iblk,1,crv10[1]);
				ns.SetValue(iblk,2,crv10[2]);
			}
		}
	}
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// EDGEについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(EDGE);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		assert( id_na_val == 0 );
	}
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// EDGEについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(BUBBLE);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		assert( id_na_val == 0 );
	}
	return true;
}


bool CLinearSystem_Field::UpdateValueOfField_Newmark(
		double gamma, double dt, 
		unsigned int id_field_val, Fem::Field::CFieldWorld& world,
		Fem::Field::FIELD_DERIVATION_TYPE fdt,
		bool IsInitial )
{

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	if( field_val.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){	// 角節点を更新
		unsigned int id_na_val=0, id_ns_v,id_ns_u;
		{
			const CField::CNodeSegInNodeAry& nsna = field_val.GetNodeSegInNodeAry(CORNER);
			id_na_val = nsna.id_na_va;
			if( fdt == VELOCITY ){
				id_ns_v = nsna.id_ns_ve; 
				id_ns_u = nsna.id_ns_va;
			}
			else if( fdt == ACCELERATION ){
				id_ns_v = nsna.id_ns_ac; 
				id_ns_u = nsna.id_ns_ve;
			}
			else{ assert(0); }
		}
		const int ilss = this->FindIndexArray_Seg(id_field_val,CORNER,world);
		assert( ilss != -1 );
		const CLinSysSeg_Field& lss = this->m_aSegField[ilss];
		const CVector_Blk* m_Update = m_ls.m_Update[ilss]; assert( m_Update != 0 );
		CNodeAry& na_val = world.GetNA(id_na_val);
		if( lss.id_field == id_field_val ){
			// 値を更新する．
			if( IsInitial ){ na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v,dt); }
			na_val.AddValueToNodeSegment(id_ns_u ,*m_Update, gamma*dt);
			// 速度を更新する
			na_val.AddValueToNodeSegment(id_ns_v,*m_Update,1.0);
		}
		else{
			assert( lss.id_field2 == id_field_val );
			unsigned int len0 = 0;
			{
				const CField& field = world.GetField(lss.id_field);
				len0 = field.GetNLenValue();
			}
			// 値を更新する．
			if( IsInitial ){ na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v,dt); }
			na_val.AddValueToNodeSegment(id_ns_u, *m_Update, gamma*dt, len0);
			// 速度を更新する
			na_val.AddValueToNodeSegment(id_ns_v, *m_Update, 1.0,      len0);
		}
	}
	{	// 気泡節点を更新
		unsigned int id_na_val=0,  id_ns_u,id_ns_v;
		{
			const CField::CNodeSegInNodeAry& nsna = field_val.GetNodeSegInNodeAry(BUBBLE);
			id_na_val = nsna.id_na_va;
			if( fdt == VELOCITY ){
				id_ns_v = nsna.id_ns_ve; 
				id_ns_u = nsna.id_ns_va;
			}
			else if( fdt == ACCELERATION ){
				id_ns_v = nsna.id_ns_ac; 
				id_ns_u = nsna.id_ns_ve;
			}
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field_val,BUBBLE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na_val = world.GetNA(id_na_val);	
			// 値を更新する
			if( IsInitial ){ na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v,dt);	}
			na_val.AddValueToNodeSegment(id_ns_u, upd, gamma*dt);
			// 速度を更新する
			na_val.AddValueToNodeSegment(id_ns_v, upd, 1.0);
		}
	}
	return true;
}


bool CLinearSystem_Field::UpdateValueOfField_NewmarkBeta
(
 double gamma, double beta, double dt, 
 unsigned int id_field, Fem::Field::CFieldWorld& world, bool IsInitial )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}

	if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){	// 角節点を更新
		unsigned int id_na_val=0, id_ns_a,id_ns_v,id_ns_u;
		{
			const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(CORNER);
			id_na_val = nsna.id_na_va;
			id_ns_u = nsna.id_ns_va;
			id_ns_v = nsna.id_ns_ve;
			id_ns_a = nsna.id_ns_ac;
		}
		const int ilss = this->FindIndexArray_Seg(id_field,CORNER,world);
		assert( ilss != -1 );
		const CLinSysSeg_Field& lss = this->m_aSegField[ilss];
		const CVector_Blk* m_Update = m_ls.m_Update[ilss]; assert( m_Update != 0 );
		CNodeAry& na_val = world.GetNA(id_na_val);
		if( lss.id_field == id_field_parent ){
			// 値を更新する
			if( IsInitial ){
				na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v, dt);
				na_val.AddValueToNodeSegment(id_ns_u ,id_ns_a, 0.5*dt*dt);
			}
			na_val.AddValueToNodeSegment(id_ns_u, *m_Update, beta*dt*dt);
			// 速度を更新する．
			if( IsInitial ){ 
				na_val.AddValueToNodeSegment(id_ns_v ,id_ns_a, dt); 
			}
			na_val.AddValueToNodeSegment(id_ns_v, *m_Update, gamma*dt);
			// 加速度を更新する
			na_val.AddValueToNodeSegment(id_ns_a, *m_Update, 1.0);
		}
		else{   // CombinedFieldの場合
			assert( lss.id_field2 == id_field_parent );
			unsigned int len0 = 0;
			{
				const CField& field = world.GetField(lss.id_field);
				len0 = field.GetNLenValue();
			}
			// 値を更新する
			if( IsInitial ){
				na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v, dt);
				na_val.AddValueToNodeSegment(id_ns_u ,id_ns_a, 0.5*dt*dt);
			}
			na_val.AddValueToNodeSegment(id_ns_u, *m_Update, beta*dt*dt, len0);
			// 速度を更新する．
			if( IsInitial ){ 
				na_val.AddValueToNodeSegment(id_ns_v ,id_ns_a, dt); 
			}
			na_val.AddValueToNodeSegment(id_ns_v, *m_Update, gamma*dt, len0);
			// 加速度を更新する
			na_val.AddValueToNodeSegment(id_ns_a, *m_Update, 1.0, len0);
		}
	}
	{	// バブル節点を更新
		unsigned int id_na_val=0, id_ns_a,id_ns_v,id_ns_u;
		{
			const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(BUBBLE);
			id_na_val = nsna.id_na_va;
			id_ns_u = nsna.id_ns_va;
			id_ns_v = nsna.id_ns_ve;
			id_ns_a = nsna.id_ns_ac;
		}
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,BUBBLE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na_val = world.GetNA(id_na_val);
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
	}
	return true;
}

bool CLinearSystem_Field::UpdateValueOfField_BackwardEular
(double dt, 
 unsigned int id_field, Fem::Field::CFieldWorld& world, 
 bool IsInitial )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}
	
	if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){	// update value for corner node
		unsigned int id_na_val=0, id_ns_v,id_ns_u;
		{
			const CField::CNodeSegInNodeAry& nsna = field.GetNodeSegInNodeAry(CORNER);
			id_na_val = nsna.id_na_va;
			id_ns_u = nsna.id_ns_va;
			id_ns_v = nsna.id_ns_ve;
		}
		const int ilss = this->FindIndexArray_Seg(id_field,CORNER,world);
		assert( ilss != -1 );
		const CLinSysSeg_Field& lss = this->m_aSegField[ilss];
		const CVector_Blk* m_Update = m_ls.m_Update[ilss]; assert( m_Update != 0 );
		CNodeAry& na_val = world.GetNA(id_na_val);
		if( lss.id_field == id_field_parent ){
			na_val.AddValueToNodeSegment(id_ns_v, *m_Update, 1.0);
			if( IsInitial ){
				na_val.AddValueToNodeSegment(id_ns_u ,id_ns_v, dt);
			}
			else{				
				na_val.AddValueToNodeSegment(id_ns_u, *m_Update, dt);
			}
		}
		else{   // in case of combined field
			assert(0);
		}
	}
	
}


////////////////////////////////////////////////////////////////
// Private 関数
////////////////////////////////////////////////////////////////
	
int CLinearSystem_Field::FindIndexArray_Seg( unsigned int id_field, const ELSEG_TYPE& type, const CFieldWorld& world )
{	
	if( !world.IsIdField(id_field) ) return -1;
	const CField& field = world.GetField(id_field);
  unsigned int id_field_parent = field.GetIDFieldParent();
  if( id_field_parent == 0 ){ id_field_parent = id_field; }  
  ////
  for(unsigned int ils=0;ils<m_aSegField.size();ils++){
    if( m_aSegField[ils].id_field ==id_field_parent && m_aSegField[ils].node_config==type ){
			return ils;
		}
		if( m_aSegField[ils].id_field2==id_field_parent && m_aSegField[ils].node_config==type ){
			return ils;
		}
	}
	return -1;
}

// 連立一次方程式セグメントを一番最後に加える。（行列、ベクトルなどをリサイズ)
int CLinearSystem_Field::AddLinSysSeg( const CLinSysSeg_Field& seg )
{
    const unsigned int nseg_pre = this->m_aSegField.size();
    assert( m_ls.GetNLinSysSeg() == nseg_pre );
    int res = m_ls.AddLinSysSeg( seg.nnode, seg.len );
    assert( res != - 1 );
    assert( res == (int)nseg_pre );
    this->m_aSegField.push_back( seg );
    return res;
}


