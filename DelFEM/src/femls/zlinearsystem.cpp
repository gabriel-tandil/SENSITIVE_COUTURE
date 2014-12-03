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
#pragma warning( disable : 4786 )
#endif

#include <math.h>
#include <cstdlib> //(abort)

#include "delfem/field_world.h"

#include "delfem/indexed_array.h"
#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/zmat_blkcrs.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/diamat_blk.h"

#include "delfem/femls/zlinearsystem.h"

#ifndef for 
#define for if(0); else for
#endif

using namespace MatVec;
using namespace Fem::Ls;
using namespace Fem::Field;

CZLinearSystem::CZLinearSystem(){}

CZLinearSystem::~CZLinearSystem(){
	this->Clear();
}

// 対角にパターンを追加
bool CZLinearSystem::AddPattern_Field(const unsigned int id_field, const CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}

	unsigned int nlen_value;
	{	// 複素数の場合は２で割る
		const unsigned int nlen = field.GetNLenValue();
		assert( nlen % 2 == 0 );
		nlen_value = nlen / 2;
	}

	int ESType2iLS[3];
	{	// Bubbleブロックを作る
		unsigned int id_na_val = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			CLinSysSeg seg;
			{
				seg.id_field = id_field_parent; seg.node_config = BUBBLE;
				seg.len=nlen_value; seg.nnode=na.Size();
			}
			ESType2iLS[2] = this->AddLinSysSeg(seg);
		}
		else{ ESType2iLS[2] = -1; }
	}
	{	// Edgeブロックを作る
		unsigned int id_na_val = field.GetNodeSegInNodeAry(EDGE).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			CLinSysSeg seg;
			{
				seg.id_field = id_field_parent; seg.node_config = EDGE;
				seg.len=nlen_value; seg.nnode=na.Size();
			}
			ESType2iLS[1] = this->AddLinSysSeg(seg);
		}
		else{ ESType2iLS[1] = -1; }
	}
	{	// Cornerブロックを作る
		unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		if( id_na_val != 0 ){
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			CLinSysSeg seg;
			{
				seg.id_field = id_field_parent; seg.node_config = CORNER;
				seg.len=nlen_value; seg.nnode=na.Size();
			}
			ESType2iLS[0] = this->AddLinSysSeg(seg);
		}
		else{ ESType2iLS[0] = -1; }
	}
	////////////////////////////////
	const std::vector<unsigned int> aIdEA = field.GetAryIdEA();
	if( aIdEA.size() == 0 ){ // 剛体モードのための行列
		unsigned int ils0 = ESType2iLS[0];
		if( m_Matrix_Dia[ils0] == 0 ){
			m_Matrix_Dia[ils0] = new CZMatDia_BlkCrs(1, nlen_value);
		}
		return true;
	}

	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
	{
		const unsigned int id_ea = aIdEA[iiea];
		const CElemAry& ea = world.GetEA(id_ea);
		// CORNER節点について
		if( field.GetIdElemSeg(id_ea,CORNER,true,world) != 0 ){
			assert( world.IsIdEA(id_ea) );
			const unsigned int id_es_c = field.GetIdElemSeg(id_ea,CORNER,true,world);
			assert( ea.IsSegID(id_es_c) );
			const unsigned int ils0 = ESType2iLS[0];
			this->AddMat_Dia(ils0, ea, id_es_c );			// cc行列を作る
			if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){	// CORNER-BUBBLE
				const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
				assert( ea.IsSegID(id_es_b) );
				const unsigned int ils1 = ESType2iLS[2];
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_c,id_es_b,crs);
				assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// cb行列を作る
				const unsigned int nnode1 = m_aSeg[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// bc行列を作る
			}
			if( field.GetIdElemSeg(id_ea,EDGE,true,world) != 0 ){	// CONRER-EDGE
				const unsigned int id_es_e = field.GetIdElemSeg(id_ea,EDGE,true,world);
				assert( ea.IsSegID(id_es_e) );
				const unsigned int ils1 = ESType2iLS[1];
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_c,id_es_e,crs);
				assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// ce行列を作る
				const unsigned int nnode1 = m_aSeg[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// ec行列を作る
			}
		}
		// EDGE節点について
		if( field.GetIdElemSeg(id_ea,EDGE,true,world) != 0 ){
			const unsigned int id_es_e = field.GetIdElemSeg(id_ea,EDGE,true,world);
			assert( ea.IsSegID(id_es_e) );
			const unsigned int ils0 = ESType2iLS[1];
			this->AddMat_Dia(ils0, ea, id_es_e);			// ee行列を作る
			if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){	// EDGE-BUBBLE
				const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
				assert( ea.IsSegID(id_es_b) );
				const unsigned int ils1 = ESType2iLS[2];
				Com::CIndexedArray crs;
				ea.MakePattern_FEM(id_es_e,id_es_b,crs);
				assert( crs.CheckValid() );
				this->AddMat_NonDia(ils0,ils1, crs);		// eb行列を作る
				const unsigned int nnode1 = m_aSeg[ils1].nnode;
				Com::CIndexedArray crs_inv;
				crs_inv.SetTranspose(nnode1,crs);
				this->AddMat_NonDia(ils1,ils0, crs_inv);	// be行列を作る
			}
		}
		// BUBBLE節点について
		if( field.GetIdElemSeg(id_ea,BUBBLE,true,world) != 0 ){
			const unsigned int id_es_b = field.GetIdElemSeg(id_ea,BUBBLE,true,world);
			assert( ea.IsSegID(id_es_b) );
			const unsigned int ils0 = ESType2iLS[2];
			this->AddMat_Dia(ils0, ea, id_es_b);
		}
	}
	return true;
}

void CZLinearSystem::Clear()
{
	m_aSeg.clear();
	////////////////
	for(unsigned int i=0;i<m_Matrix_NonDia.size();i++){
		for(unsigned int j=0;j<m_Matrix_NonDia[i].size();j++){
			delete m_Matrix_NonDia[i][j];
		}
		m_Matrix_NonDia[i].clear();
	}
	m_Matrix_NonDia.clear();
	////////////////
	for(unsigned int i=0;i<m_Matrix_Dia.size();i++){
		delete m_Matrix_Dia[i];
	}
	m_Matrix_Dia.clear();
	////////////////
	for(unsigned int i=0;i<m_Residual.size();i++){
		delete m_Residual[i];
	}
	m_Residual.clear();
	////////////////
	for(unsigned int i=0;i<m_Update.size();i++){
		delete m_Update[i];
	}
	m_Update.clear();
	////////////////
	for(unsigned int i=0;i<m_BCFlag.size();i++){
		delete m_BCFlag[i];
	}
	m_BCFlag.clear();
	////////////////
	for(unsigned int i=0;i<m_TmpVectorArray.size();i++){
		for(unsigned int j=0;j<m_TmpVectorArray[i].size();j++){
			delete m_TmpVectorArray[i][j];
		}
		m_TmpVectorArray[i].clear();
	}
	m_TmpVectorArray.clear();
}

// 残差ベクトルをゲットする
CZVector_Blk* CZLinearSystem::GetResidualPtr(unsigned int id_field, const ELSEG_TYPE& elseg_type, 
										const CFieldWorld& world){
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
	if( ils0 < 0 || ils0 >= (int)m_aSeg.size() ) return 0;
	return m_Residual[ils0];
}

// 更新ベクトルをゲットする
CZVector_Blk* CZLinearSystem::GetUpdatePtr(unsigned int id_field, const ELSEG_TYPE& elseg_type, 
									  const CFieldWorld& world){	
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
	if( ils0 < 0 || ils0 >= (int)m_aSeg.size() ) return 0;
	return m_Update[ils0];
}

CZMatDia_BlkCrs* CZLinearSystem::GetMatrixPtr(unsigned int id_field,const ELSEG_TYPE& elseg_type, 
										 const CFieldWorld& world){
	int ils0 = this->FindIndexArray_Seg(id_field,elseg_type,world);
	if( ils0 < 0 || ils0 >= (int)m_aSeg.size() ) return 0;
	return m_Matrix_Dia[ils0];
}

CZMat_BlkCrs* CZLinearSystem::GetMatrixPtr(unsigned int id_field_col,const ELSEG_TYPE& elseg_type_col,
	unsigned int id_field_row,const ELSEG_TYPE& elseg_type_row,
	const CFieldWorld& world){

	int ils_col = FindIndexArray_Seg(id_field_col,elseg_type_col,world);
	if( ils_col < 0 || ils_col >= (int)m_aSeg.size() ) return 0;
	int ils_row = FindIndexArray_Seg(id_field_row,elseg_type_row,world);
	if( ils_row < 0 || ils_row >= (int)m_aSeg.size() ) return 0;

	return m_Matrix_NonDia[ils_col][ils_row];
}

unsigned int CZLinearSystem::GetTmpBufferSize(){
	unsigned int size=0;
	for(unsigned int ils=0;ils<m_aSeg.size();ils++){
		size = ( size > m_aSeg[ils].nnode ) ? size : m_aSeg[ils].nnode;
	}
	return size;
}

// マージ前の初期化
void CZLinearSystem::InitializeMarge()
{	
	unsigned int nls = m_aSeg.size();
	for(unsigned int ils=0;ils<nls;ils++){
		for(unsigned int jls=0;jls<nls;jls++){
			if( m_Matrix_NonDia[ils][jls] != 0 ){
				m_Matrix_NonDia[ils][jls]->SetZero();
			}
		}
		if( m_Matrix_Dia[ils] != 0 ){ m_Matrix_Dia[ils]->SetZero(); }
	}
	for(unsigned int ils=0;ils<nls;ils++){
		m_Residual[ils]->SetVectorZero();
		m_Update[ils]->SetVectorZero();
	}
}


// マージ後の処理（残差ノルムを返す)
double CZLinearSystem::FinalizeMarge()
{
	unsigned int nseg = m_aSeg.size();
	{	// 境界条件の行列へのセット
		for(unsigned int iseg=0;iseg<nseg;iseg++){
			if( m_Matrix_Dia[iseg] != 0 ){
				m_Matrix_Dia[iseg]->SetBoundaryCondition(*m_BCFlag[iseg]);
			}
		}
		for(unsigned int iseg=0;iseg<nseg;iseg++){
			for(unsigned int jseg=0;jseg<nseg;jseg++){
				if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
				assert( iseg != jseg );
				m_Matrix_NonDia[iseg][jseg]->SetBoundaryCondition_Row(*m_BCFlag[iseg]);
			}
		}
		for(unsigned int jseg=0;jseg<nseg;jseg++){
			for(unsigned int iseg=0;iseg<nseg;iseg++){
				if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
				assert( iseg != jseg );
				m_Matrix_NonDia[iseg][jseg]->SetBoundaryCondition_Colum(*m_BCFlag[jseg]);
			}
		}
	}

	// 残差への境界条件のセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		m_BCFlag[iseg]->SetZeroToBCDof(*m_Residual[iseg]);
	}

	// 残差ノルムの計算
	double sq_norm_res = 0.0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		sq_norm_res += m_Residual[iseg]->GetSquaredVectorNorm();
	}
	return sqrt(sq_norm_res);
}



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


static void BoundaryCondition
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

static void BoundaryCondition
(unsigned int id_field, const ELSEG_TYPE& elseg_type,  
 MatVec::CBCFlag& bc_flag, const CFieldWorld& world,  
 unsigned int ioffset=0)
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




bool CZLinearSystem::SetFixedBoundaryCondition_Field( unsigned int id_field, unsigned int idofns, const CFieldWorld& world ){
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent = field.GetIDFieldParent();
	if( id_field_parent == 0 ){ id_field_parent = id_field; }
	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,CORNER,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,CORNER,idofns,*m_BCFlag[ils0],world);
		}
	}
	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,EDGE,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,EDGE,idofns,*m_BCFlag[ils0],world);
		}
	}
	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,BUBBLE,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,BUBBLE,idofns,*m_BCFlag[ils0],world);
		}
	}
	return true;
}

bool CZLinearSystem::SetFixedBoundaryCondition_Field( unsigned int id_field, const CFieldWorld& world )
{
	if( !world.IsIdField(id_field) ) return false;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent = field.GetIDFieldParent();
	if( id_field_parent == 0 ) id_field_parent = id_field;

	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,CORNER,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,CORNER,*m_BCFlag[ils0],world);
		}
	}
	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,EDGE,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,EDGE,*m_BCFlag[ils0],world);
		}
	}
	{
		unsigned int ils0 = this->FindIndexArray_Seg(id_field_parent,BUBBLE,world);
		if( ils0 < m_aSeg.size() ){
			BoundaryCondition(id_field,BUBBLE,*m_BCFlag[ils0],world);
		}
	}
	return true;
}

void CZLinearSystem::ClearFixedBoundaryCondition(){
	for(unsigned int ibcflag=0;ibcflag<m_BCFlag.size();ibcflag++){
		m_BCFlag[ibcflag]->SetAllFlagZero();
	}
}

bool CZLinearSystem::UpdateValueOfField( 
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
			CZVector_Blk* m_Update = this->GetUpdatePtr(id_field,CORNER,world);
			assert( m_Update != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			na.AddValueToNodeSegment(id_ns_val,*m_Update,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
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
			CZVector_Blk* m_Update = this->GetUpdatePtr(id_field,EDGE,world);
			assert( m_Update != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			na.AddValueToNodeSegment(id_ns_val,*m_Update,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
//			ns.is_updated = true;
		}
	}
	{
		unsigned int id_na_val=0, id_ns_val=0;
		{	// BUBBLEについて更新
			const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(BUBBLE);
			id_na_val = nsna_c.id_na_va;
			if(      fdt == VALUE )        id_ns_val = nsna_c.id_ns_va;
			else if( fdt == VELOCITY )     id_ns_val = nsna_c.id_ns_ve;
			else if( fdt == ACCELERATION ) id_ns_val = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		if( id_na_val != 0 ){
			CZVector_Blk* m_Update = this->GetUpdatePtr(id_field,BUBBLE,world);
			assert( m_Update != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			assert( na.IsSegID(id_ns_val) );
			na.AddValueToNodeSegment(id_ns_val,*m_Update,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
//			ns.is_updated = true;
		}
	}
	return true;
}


////////////////////////////////////////////////////////////////
// ソルバ用のユーティリティ関数
////////////////////////////////////////////////////////////////
bool CZLinearSystem::ReSizeTmpVecSolver(unsigned int ntmp_new)
{
	const unsigned int nseg = m_aSeg.size();
	const unsigned int ntmp_old = this->GetTmpVectorArySize();

	if( ntmp_old == ntmp_new ){ return true; }
	else if( ntmp_old < ntmp_new ){
		m_TmpVectorArray.resize(ntmp_new);
		for(unsigned int ivec=ntmp_old;ivec<ntmp_new;ivec++){
			m_TmpVectorArray[ivec].resize(nseg);
			for(unsigned int iseg=0;iseg<nseg;iseg++){
				const int nnode = m_aSeg[iseg].nnode;
				const int len = m_aSeg[iseg].len;
				m_TmpVectorArray[ivec][iseg] = new CZVector_Blk(nnode,len);
				m_TmpVectorArray[ivec][iseg]->SetVectorZero();
			}
		}
	}
	else{
		assert( ntmp_old > ntmp_new );
		for(unsigned int ivec=ntmp_new-1;ivec>=ntmp_old;ivec--){
			for(unsigned int iseg=0;iseg<nseg;iseg++){
				delete m_TmpVectorArray[ivec][iseg];
			}
			m_TmpVectorArray[ivec].clear();
		}
		m_TmpVectorArray.resize(ntmp_old);
	}
	return true;
}

// ベクトルの共役を取る
bool CZLinearSystem::Conjugate(int iv1){
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec1)[iseg]).SetVectorConjugate();
	}

	return true;
}

////////////////////////////////
// 行列ベクトル積
// {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
// v1 != v2     v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
bool CZLinearSystem::MatVec(double alpha, int iv1, double beta, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	if( alpha == 0.0 ){
		this->SCAL(beta,iv2);
		return true;
	}

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	if( iv1 == iv2 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		if( this->m_Matrix_Dia[iseg] != 0 ){
			m_Matrix_Dia[iseg]->MatVec( alpha, (*(*p_vec1)[iseg]), beta, (*(*p_vec2)[iseg]) );
		}
		else{ (*(*p_vec2)[iseg]) *= beta; }
		for(unsigned int jseg=0;jseg<nseg;jseg++){
			if( m_Matrix_NonDia[iseg][jseg] == 0 ) continue;
			assert( iseg != jseg );
			m_Matrix_NonDia[iseg][jseg]->MatVec( alpha, (*(*p_vec1)[jseg]), 1.0, (*(*p_vec2)[iseg]), true );
		}
	}
	return true;
}

////////////////////////////////
// 行列ベクトル積
// {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
// v1 != v2     v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
bool CZLinearSystem::MatVec_Hermitian(double alpha, int iv1, double beta, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	if( alpha == 0.0 ){
		this->SCAL(beta,iv2);
		return true;
	}

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	if( iv1 == iv2 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}

	assert( nseg == 1 );

	if( this->m_Matrix_Dia[0] != 0 ){
		m_Matrix_Dia[0]->MatVec_Hermitian( alpha, (*(*p_vec1)[0]), beta, (*(*p_vec2)[0]) );
	}
	else{ (*(*p_vec2)[0]) *= beta; }
	return true;
}

////////////////////////////////
// ベクトル同士の足し算
// {v2} := alpha*{v1} + {v2}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
bool CZLinearSystem::AXPY(Com::Complex alpha, int iv1, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;
//	if( alpha == 0 ) return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	if( iv1 == iv2 ){
		std::cout << "Error!-->未実装" << std::endl;
		assert(0);
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec2)[iseg]).AXPY( alpha, (*(*p_vec1)[iseg]) );
	}

	return true;
}

////////////////////////////////
// 内積を求める関数
// return {v1} * {v2}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル
Com::Complex CZLinearSystem::DOT(int iv1, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	Com::Complex norm = 0.0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		norm += (*(*p_vec1)[iseg]) * (*(*p_vec2)[iseg]);
	}
	return norm;
}


// ベクトルの内積
// return {v1} * {v2}^H
Com::Complex CZLinearSystem::INPROCT(int iv1, int iv2)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	Com::Complex norm = 0.0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		norm += InnerProduct( (*(*p_vec1)[iseg]), (*(*p_vec2)[iseg]) );
	}
	return norm;
}

////////////////////////////////
// ベクトルのコピー
// return {v2} := {v1}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル 
bool CZLinearSystem::COPY(int iv1, int iv2){

	if( iv1 == iv2 ) return true;

	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else assert(0);

	std::vector< CZVector_Blk* >* p_vec2 = 0;
	if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = &m_TmpVectorArray[iv2];
	else if( iv2 == -1 ) p_vec2 = &this->m_Residual;
	else if( iv2 == -2 ) p_vec2 = &this->m_Update;
	else assert(0);

	assert( p_vec1->size() == nseg );
	assert( p_vec2->size() == nseg );

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec2)[iseg]) = (*(*p_vec1)[iseg]);
	}
	return true;
}

////////////////////////////////
// ベクトルのスカラー倍
// return {v1} := alpha * {v1}
// v?=-1 : v?が右辺ベクトル    v?=-2 : v?が左辺ベクトル 
bool CZLinearSystem::SCAL(Com::Complex alpha,int iv1){

	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else{ assert(0); }

	assert( p_vec1->size() == nseg );

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*(*p_vec1)[iseg]) *= alpha;
	}
	return true;
}

////////////////////////////////////////////////////////////////
// Private 関数
////////////////////////////////////////////////////////////////
	
int CZLinearSystem::FindIndexArray_Seg( unsigned int id_field, const ELSEG_TYPE& type, const CFieldWorld& world )
{	
	if( !world.IsIdField(id_field) ) return -1;
	const CField& field = world.GetField(id_field);
	unsigned int id_field_parent;
	{
		if( field.GetIDFieldParent() == 0 ){ id_field_parent = id_field; }
		else{ id_field_parent = field.GetIDFieldParent(); }
	}
	for(unsigned int ils=0;ils<m_aSeg.size();ils++){
		if( m_aSeg[ils].id_field==id_field_parent && m_aSeg[ils].node_config==type ){
			return ils;
		}
	}
	return -1;
}

// 連立一次方程式セグメントを一番最後に加える。（行列、ベクトルなどをリサイズ)
int CZLinearSystem::AddLinSysSeg( const CLinSysSeg& seg ){
	std::vector< std::vector< CZMat_BlkCrs* > > old_matrix_nd;
	old_matrix_nd = m_Matrix_NonDia;
	m_Matrix_NonDia.resize(0);
	const unsigned int size_old = m_aSeg.size();
	const unsigned int size_new = size_old+1;
	m_Matrix_NonDia.resize( size_new );
	for(unsigned int ils=0;ils<size_new;ils++){
		m_Matrix_NonDia[ils].resize( size_new, 0 );
	}
	for(unsigned int ils=0;ils<size_old;ils++){
		for(unsigned int jls=0;jls<size_old;jls++){
			m_Matrix_NonDia[ils][jls] = old_matrix_nd[ils][jls];
		}
	}
	m_Matrix_Dia.resize( size_new, 0 );
	m_Residual.resize( size_new, 0 );
	m_Update.resize( size_new, 0 );
	m_BCFlag.resize( size_new, 0 );
	m_Residual[ size_old ] = new CZVector_Blk( seg.nnode, seg.len );
	m_Update[ size_old ] = new CZVector_Blk( seg.nnode, seg.len );
	m_BCFlag[ size_old ] = new CBCFlag( seg.nnode, seg.len );
	m_aSeg.push_back(seg);
	return size_old;
}

bool CZLinearSystem::AddMat_Dia(unsigned int ils, const CElemAry& ea, unsigned int id_es)
{
	assert( ils < m_aSeg.size() );
	assert( ea.IsSegID(id_es) );
	{	// 行列がなければ作るよ
		const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
		const unsigned int nnode = m_aSeg[ils].nnode;
		const unsigned int len = m_aSeg[ils].len;
		assert( es.GetMaxNoes() < nnode );
		if( m_Matrix_Dia[ils] == 0 ){
			m_Matrix_Dia[ils] = new CZMatDia_BlkCrs(nnode,len);
		}
		else{	// assertする
			assert( m_Matrix_Dia[ils]->NBlkMatCol() ==  nnode );
			assert( m_Matrix_Dia[ils]->LenBlkCol()  ==  len );
			assert( m_Matrix_Dia[ils]->NBlkMatRow() ==  nnode );
			assert( m_Matrix_Dia[ils]->LenBlkRow()  ==  len );
		}
	}
    Com::CIndexedArray crs;
	ea.MakePattern_FEM(id_es,crs);
	assert( crs.CheckValid() );
	m_Matrix_Dia[ils]->AddPattern(crs);
	return true;
}

bool CZLinearSystem::AddMat_NonDia(unsigned int ils_col, unsigned int ils_row, const Com::CIndexedArray& crs )
{
	assert( ils_col < m_aSeg.size() );
	assert( ils_row < m_aSeg.size() );
	{
		const unsigned int nnode_col = m_aSeg[ils_col].nnode;
		const unsigned int len_col = m_aSeg[ils_col].len;
		const unsigned int nnode_row = m_aSeg[ils_row].nnode;
		const unsigned int len_row = m_aSeg[ils_row].len;
		if( m_Matrix_NonDia[ils_col][ils_row] == 0 ){
			m_Matrix_NonDia[ils_col][ils_row] = new CZMat_BlkCrs(nnode_col,len_col, nnode_row,len_row);
		}
		else{
			assert( m_Matrix_NonDia[ils_col][ils_row]->NBlkMatCol() ==  nnode_col );
			assert( m_Matrix_NonDia[ils_col][ils_row]->LenBlkCol()  ==  len_col );
			assert( m_Matrix_NonDia[ils_col][ils_row]->NBlkMatRow() ==  nnode_row );
			assert( m_Matrix_NonDia[ils_col][ils_row]->LenBlkRow()  ==  len_row );
		}
	}
	assert( crs.CheckValid() );
	m_Matrix_NonDia[ils_col][ils_row]->AddPattern(crs);
	return true;
}




bool CZLinearSystem::NormalizeVector(int iv1){
		
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* paVec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) paVec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) paVec1 = &this->m_Residual;
	else if( iv1 == -2 ) paVec1 = &this->m_Update;
	else assert(0);

	assert( (*paVec1).size() == nseg );

	double sq_norm = 0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		const CZVector_Blk* pVec = (*paVec1)[iseg];
		sq_norm += pVec->GetSquaredVectorNorm();
	}
	std::cout << "Norm : " << sqrt(sq_norm) << std::endl;
	const double inv_norm = 1.0 / sqrt(sq_norm);
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		CZVector_Blk* pVec = (*paVec1)[iseg];
		(*pVec) *= inv_norm;
	}
	return true;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



bool CZLinearSystem_GeneralEigen::SetVector_fromField(int iv1,
	unsigned int id_field, const Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
	if( !world.IsIdField(id_field) ){ assert(0); return false; }
	const CField& field = world.GetField(id_field);
	if( field.GetFieldType() != ZSCALAR ) return false;
	
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return true;

	std::vector< CZVector_Blk* >* paVec1 = 0;
    if(      iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) paVec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) paVec1 = &this->m_Residual;
	else if( iv1 == -2 ) paVec1 = &this->m_Update;
	else assert(0);


	// コーナー節点について値を更新
	{
/*		const int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
		assert( ilss0 != -1 );
		std::cout << ilss0 << std::endl;*/
        unsigned int ilss0 = 0;
		assert( ilss0 < (*paVec1).size() );
		CZVector_Blk* pVec = (*paVec1)[ilss0];
		const CField& field = world.GetField(id_field);
		unsigned int id_na = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		const CNodeAry& na = world.GetNA(id_na);
		const CNodeAry::CNodeSeg& ns = field.GetNodeSeg(CORNER,true,world,fdt);
		const unsigned int nblk = na.Size();
		const unsigned int nlen = ns.Length();
		assert( nlen == 2 );
		double* val = new double [nlen];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			ns.GetValue(iblk,val);
			pVec->SetValue(iblk,0, Com::Complex(val[0],val[1]) );
		}
		delete[] val;
	}
	return true;
}

void CZLinearSystem_GeneralEigen::Clear(){
	CZLinearSystem::Clear();
	const unsigned int nlss = this->GetNLynSysSeg();

	assert( m_DiaMassMatrix.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		delete m_DiaMassMatrix[ilss];
	}
	m_DiaMassMatrix.clear();
}

bool CZLinearSystem_GeneralEigen::AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world){
	CZLinearSystem::AddPattern_Field(id_field,world);

	int ilss = this->FindIndexArray_Seg(id_field,CORNER,world);
    assert( ilss != -1 && ilss < (int)this->GetNLynSysSeg() );
    assert( ilss == (int)this->m_DiaMassMatrix.size() );
	const unsigned int len = this->m_Update[ilss]->BlkLen();
	const unsigned int nblk = this->m_Update[ilss]->BlkVecLen();
	m_DiaMassMatrix.resize(ilss+1,0);
	m_DiaMassMatrix[ilss] = new CDiaMat_Blk(nblk,len);
	return true;
}

void CZLinearSystem_GeneralEigen::InitializeMarge(){
	CZLinearSystem::InitializeMarge();
	assert( this->GetNLynSysSeg() == 1 );
	m_DiaMassMatrix[0]->SetZero();
}


CDiaMat_Blk* CZLinearSystem_GeneralEigen::GetDiaMassMatrixPtr(
	unsigned int id_field, const Fem::Field::ELSEG_TYPE& elseg_type, const Fem::Field::CFieldWorld& world)
{
	int ilss = this->FindIndexArray_Seg(id_field,elseg_type,world);
	if( ilss == -1 ) return 0;
	return m_DiaMassMatrix[ilss];
}


bool CZLinearSystem_GeneralEigen::MultVecMassDecomp(int ivec){
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	CZVector_Blk* pVec = 0;
	if( ivec >= 0 ){
        assert( ivec < (int)this->GetTmpVectorArySize() );
		pVec = this->m_TmpVectorArray[ivec][0];
	}
	else if( ivec == -2 ){	pVec = this->m_Update[0];   }
	else{					pVec = this->m_Residual[0]; }
	CDiaMat_Blk* dmat = this->m_DiaMassMatrix[0];
	const unsigned int nblk = dmat->NBlk();
	const unsigned int nlen = dmat->LenBlk();
	assert( nlen == 1 );
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		const double* pValMi = dmat->GetPtrValDia(iblk);
		Com::Complex val0 = pVec->GetValue(iblk,0);
		const Com::Complex val1 = val0/pValMi[0];
		pVec->SetValue(iblk,0,val1);
	}
	return true;
}

bool CZLinearSystem_GeneralEigen::MultUpdateInvMassDecomp()
{
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	{
		CZVector_Blk* pUpdate = this->m_Update[0];
		CDiaMat_Blk* dmat = this->m_DiaMassMatrix[0];
		const unsigned int nblk = dmat->NBlk();
		const unsigned int nlen = dmat->LenBlk();
		assert( nlen == 1 );
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			const double* pValMi = dmat->GetPtrValDia(iblk);
			Com::Complex val = pUpdate->GetValue(iblk,0);
			val *= (*pValMi);
			pUpdate->SetValue(iblk,0,val);
		}
	}
	return true;
}

void CZLinearSystem_GeneralEigen::OffsetDiagonal(double lambda){
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	CZMatDia_BlkCrs* dmat = this->m_Matrix_Dia[0];
	const unsigned int nblk = dmat->NBlkMatCol();
	const unsigned int nlen = dmat->LenBlkCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		Com::Complex* pVal = dmat->GetPtrValDia(iblk);
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			pVal[ilen*nlen+ilen] -= lambda;
		}
	}
}

bool CZLinearSystem_GeneralEigen::DecompMultMassMatrix()
{
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	m_DiaMassMatrix[0]->CholeskyDecomp();
	{
		CZMatDia_BlkCrs* dmat = this->m_Matrix_Dia[0];
		const unsigned int nblk = dmat->NBlkMatCol();
		const unsigned int nlen = dmat->LenBlkCol();
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			const double* pVal_Mi = m_DiaMassMatrix[0]->GetPtrValDia(iblk);
			unsigned int npsup;
			const unsigned int* psup = dmat->GetPtrIndPSuP(iblk,npsup);
			for(unsigned int ipsup=0;ipsup<npsup;ipsup++){
				const unsigned int jblk = psup[ipsup];
				assert( iblk != jblk );
				unsigned int npsup0;
				Com::Complex* pVal_Cij = dmat->GetPtrValPSuP(iblk,npsup0);
				assert( npsup0 == npsup );
				const double* pVal_Mj = m_DiaMassMatrix[0]->GetPtrValDia(jblk);
				assert( nlen == 1 );
				pVal_Cij[ipsup] = (*pVal_Mi)*pVal_Cij[ipsup]*(*pVal_Mj);
			}
			{
				Com::Complex* pVal_C = dmat->GetPtrValDia(iblk);
				assert( nlen == 1 );
				(*pVal_C) = (*pVal_Mi)*(*pVal_C)*(*pVal_Mi);
			}
		}
	}
	return true;
}

