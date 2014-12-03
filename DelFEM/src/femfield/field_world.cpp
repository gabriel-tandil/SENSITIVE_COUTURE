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
// implementation of field administration class (CFieldWorld)
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
    #pragma warning ( disable : 4996 )
#endif

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <set>

#include "delfem/matvec/vector_blk.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/elem_ary.h"
#include "delfem/mesh_interface.h"


using namespace Fem::Field;

////////////////////////////////////////////////////////////////
// 生成/消滅
////////////////////////////////////////////////////////////////

CFieldWorld::CFieldWorld(){
//	std::cout << "CFieldWorld::CFieldWorld" << std::endl;
}

CFieldWorld::CFieldWorld(const CFieldWorld& world)
{
  std::cout << " Copy Constructor World" << std::endl;
  m_map_field_conv = world.m_map_field_conv;
  {
    const std::vector<unsigned int>& aIdEA = world.GetAry_IdEA();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
      const unsigned int id_ea = aIdEA[iiea];
      const CElemAry& ea = world.GetEA(id_ea);
      CElemAry* pEA = new CElemAry(ea);
      this->m_apEA.AddObj(std::make_pair(id_ea,pEA));
    }
  }
  {
    const std::vector<unsigned int>& aIdNA = world.GetAry_IdNA();
    for(unsigned int iina=0;iina<aIdNA.size();iina++){
      const unsigned int id_na = aIdNA[iina];
      const CNodeAry& na = world.GetNA(id_na);
      CNodeAry* pNA = new CNodeAry(na);
      this->m_apNA.AddObj(std::make_pair(id_na,pNA));
    }        
  }
  {
    const std::vector<unsigned int>& aIdField = world.GetAry_IdField();
    for(unsigned int iifd=0;iifd<aIdField.size();iifd++){
      const unsigned int id_field = aIdField[iifd];
      const CField& fd = world.GetField(id_field);
      CField* pField = new CField(fd);
//      std::cout << "copy" << pField->IsNodeSeg(CORNER,false,world) << std::endl;
      this->m_apField.AddObj(std::make_pair(id_field,pField));
    }        
  }  
//	m_apEA;		//!< set of element array
//	Com::CObjSet<CNodeAry*> m_apNA;		//!< set of node array
//	Com::CObjSet<CField*> m_apField;	//!< set of field
}

CFieldWorld& CFieldWorld::operator = (const CFieldWorld& world)
{  
  std::cout << " Copy World" << std::endl;
  this->Clear();
  m_map_field_conv = world.m_map_field_conv;  
  {
    const std::vector<unsigned int>& aIdEA = world.GetAry_IdEA();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
      const unsigned int id_ea = aIdEA[iiea];
      const CElemAry& ea = world.GetEA(id_ea);
      CElemAry* pEA = new CElemAry(ea);
      this->m_apEA.AddObj(std::make_pair(id_ea,pEA));
    }
  }
  {
    const std::vector<unsigned int>& aIdNA = world.GetAry_IdNA();
    for(unsigned int iina=0;iina<aIdNA.size();iina++){
      const unsigned int id_na = aIdNA[iina];
      const CNodeAry& na = world.GetNA(id_na);
      CNodeAry* pNA = new CNodeAry(na);
      this->m_apNA.AddObj(std::make_pair(id_na,pNA));
    }        
  }
  {
    const std::vector<unsigned int>& aIdField = world.GetAry_IdField();
    for(unsigned int iifd=0;iifd<aIdField.size();iifd++){
      const unsigned int id_field = aIdField[iifd];
      const CField& fd = world.GetField(id_field);
      CField* pField = new CField(fd);
      this->m_apField.AddObj(std::make_pair(id_field,pField));
    }        
  }  
  return *this;
}


CFieldWorld::~CFieldWorld(){
//	std::cout << "CFieldWorld::~CFieldWorld" << std::endl;
	this->Clear();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void CFieldWorld::Clear()
{
	// 要素配列の削除
	std::vector<unsigned int> id_ary_pea = m_apEA.GetAry_ObjID();
	unsigned int iid_pea;
	for(iid_pea=0;iid_pea<id_ary_pea.size();iid_pea++){ 
		unsigned int id_pEA = id_ary_pea[iid_pea];
		assert( m_apEA.IsObjID(id_pEA) );
		CElemAry* pEA = m_apEA.GetObj(id_pEA);
		delete pEA;
	}
	m_apEA.Clear();

	// 節点配列の削除
	std::vector<unsigned int> id_ary_pna = m_apNA.GetAry_ObjID();
	unsigned int iid_pna;
	for(iid_pna=0;iid_pna<id_ary_pna.size();iid_pna++){ 
		unsigned int id_pNA = id_ary_pna[iid_pna];
		assert( m_apNA.IsObjID(id_pNA) );
		CNodeAry* pna = m_apNA.GetObj(id_pNA);
		delete pna;
	}
	m_apNA.Clear();

	// 場配列の削除
	std::vector<unsigned int> id_ary_pField = m_apField.GetAry_ObjID();
	unsigned int iid_pField;
	for(iid_pField=0;iid_pField<id_ary_pField.size();iid_pField++){ 
		unsigned int id_pField = id_ary_pField[iid_pField];
		assert( m_apField.IsObjID(id_pField) );
		CField* pField = m_apField.GetObj(id_pField);
		delete pField;
	}
	m_apField.Clear();

	m_map_field_conv.clear();
}

bool CFieldWorld::UpdateMeshCoord(const unsigned int id_base, const Msh::IMesh& mesh)
{
  assert( this->IsIdField(id_base) );
  if( !this->IsIdField(id_base) ) return false;
	CField& field_base = this->GetField(id_base);
	const unsigned int ndim = field_base.GetNDimCoord();
  assert( ndim == mesh.GetDimention() );
  if( ndim != mesh.GetDimention() ) return false;

	unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	Field::CNodeAry& na = this->GetNA( id_na_co );
	Field::CNodeAry::CNodeSeg& ns_coord = field_base.GetNodeSeg(CORNER,false,*this,VALUE);

	std::vector<double> coord;
  mesh.GetCoord(coord);

  assert( coord.size() == ndim*na.Size() );
  if( coord.size() != ndim*na.Size() ) return false;
		
	if( ndim == 2 ){
		for(unsigned int inode=0;inode<na.Size();inode++){
			double x_value = coord[inode*2+0];
			ns_coord.SetValue(inode, 0, x_value );
			double y_value = coord[inode*2+1];
			ns_coord.SetValue(inode, 1, y_value );
		}
	}
	else if( ndim == 3 ){ 
		for(unsigned int inode=0;inode<na.Size();inode++){
			double x_value = coord[inode*3+0];
			ns_coord.SetValue(inode, 0, x_value );
			double y_value = coord[inode*3+1];
			ns_coord.SetValue(inode, 1, y_value );
			double z_value = coord[inode*3+2];
			ns_coord.SetValue(inode, 2, z_value );
		}
	}
	
	return true;
}

bool CFieldWorld::UpdateMeshCoord(const unsigned int id_base, const unsigned int id_field_disp, 
                                  const Msh::IMesh& mesh)
{
  assert( this->IsIdField(id_field_disp) );
  if( !this->IsIdField(id_field_disp) ) return false;
	CField& disp = this->GetField(id_field_disp);
	const unsigned int ndim = disp.GetNDimCoord();
  assert( ndim == mesh.GetDimention() );
  if( ndim != mesh.GetDimention() ) return false;
  
	Field::CNodeAry::CNodeSeg& ns_c = disp.GetNodeSeg(CORNER,false,*this,VALUE);
	Field::CNodeAry::CNodeSeg& ns_u = disp.GetNodeSeg(CORNER,true, *this,VALUE);
  assert( ns_c.Size() == ns_u.Size() );  
  
	std::vector<double> coord;
  mesh.GetCoord(coord);
  
  assert( coord.size() == ndim*ns_c.Size() );
  if( coord.size() != ndim*ns_c.Size() ) return false;  
  
	if( ndim == 2 ){
		for(unsigned int inode=0;inode<ns_c.Size();inode++){
      double co0[2]; ns_c.GetValue(inode,co0);
			double co1_x = coord[inode*2+0];
			ns_c.SetValue(inode, 0, co1_x );
      ns_u.AddValue(inode, 0, co0[0]-co1_x);
			double co1_y = coord[inode*2+1];
			ns_c.SetValue(inode, 1, co1_y );
      ns_u.AddValue(inode, 1, co0[1]-co1_y);
		}
	}
	else if( ndim == 3 ){ 
		for(unsigned int inode=0;inode<ns_c.Size();inode++){
      double co0[3]; ns_c.GetValue(inode,co0);  
			double co1_x = coord[inode*3+0];
			ns_c.SetValue(inode, 0, co1_x );
      ns_u.AddValue(inode, 0, co0[0]-co1_x);      
			double co1_y = coord[inode*3+1];
			ns_c.SetValue(inode, 1, co1_y );
      ns_u.AddValue(inode, 1, co0[1]-co1_y);            
			double co1_z = coord[inode*3+2];
			ns_c.SetValue(inode, 2, co1_z );
      ns_u.AddValue(inode, 2, co0[2]-co1_z);            
		}
	}  
}


bool CFieldWorld::UpdateConnectivity( const unsigned int id_base, const Msh::IMesh& mesh )
{
	assert( this->IsIdField(id_base) );
  if( !this->IsIdField(id_base) ) return false;
	CField& field_base = this->GetField(id_base);
	unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	const CNodeAry& na = this->GetNA(id_na_co);
	const std::vector< std::pair<unsigned int, unsigned int> >& aEaEs = na.GetAryEaEs();
  const std::vector<unsigned int>& aIdMsh = mesh.GetAry_ID();
  CIDConvEAMshCad conv = this->GetIDConverter(id_base);
  std::vector<int> lnods;
  
	for(unsigned int iid_msh=0;iid_msh<aIdMsh.size();iid_msh++){
    const unsigned int id_msh = aIdMsh[iid_msh];
		const unsigned int id_ea  = conv.GetIdEA_fromMsh(id_msh);
		unsigned int id_es = 0;
		for(unsigned int iEaEs=0;iEaEs<aEaEs.size();iEaEs++){
			if( aEaEs[iEaEs].first == id_ea ){
				id_es = aEaEs[iEaEs].second;
				break;
			}
		}
    assert( this->IsIdEA(id_ea) );
		CElemAry& ea = this->GetEA(id_ea);
    assert( ea.IsSegID(id_es) );
		CElemAry::CElemSeg& es = ea.GetSeg(id_es);
    if(      ea.ElemType() == TRI ){
      std::cout << "Updated connectivity :" << id_ea << " " << id_es << std::endl;
      Msh::MSH_TYPE type = mesh.GetConnectivity(id_msh,lnods);
      assert( type == Msh::TRI );
      assert( lnods.size() == ea.Size()*3 );
      for(unsigned int itri=0;itri<ea.Size();itri++){
        for(unsigned int inotri=0;inotri<3;inotri++){
          int ino0 = lnods[itri*3+inotri];
          es.SetNodes(itri,inotri,ino0);
        }
      }
    }
    else if( ea.ElemType() == TET ){
      Msh::MSH_TYPE type = mesh.GetConnectivity(id_msh,lnods);
      assert( type == Msh::TET );
      for(unsigned int itet=0;itet<ea.Size();itet++){
        for(unsigned int inotet=0;inotet<4;inotet++){
          int ino0 = lnods[itet*4+inotet];
          es.SetNodes(itet,inotet,ino0);
        }
      }
    }
	}
	return true;  
}

unsigned int CFieldWorld::AddMesh(const Msh::IMesh& mesh)
{
  std::cout << "CFieldWorld::AddMesh" << std::endl;
	unsigned int id_na, id_ns_co;
	{
		std::vector<double> coord;
		mesh.GetCoord(coord);
		const unsigned int ndim = mesh.GetDimention();
		assert( coord.size() % ndim == 0 );
		const unsigned int nnode = coord.size() / ndim;
		CNodeAry* pna = new CNodeAry(nnode);
		id_ns_co = pna->GetFreeSegID();
		CNodeAry::CNodeSeg ns_co(ndim,"COORD");
		std::vector< std::pair<unsigned int,CNodeAry::CNodeSeg> > ns_input_ary;
		ns_input_ary.push_back( std::make_pair(id_ns_co,ns_co) );
		pna->AddSegment( ns_input_ary, coord );
		id_na = this->m_apNA.GetFreeObjID();
		const unsigned int tmp_id = this->m_apNA.AddObj( std::make_pair(id_na,pna) );
		assert( tmp_id == id_na );
	}
  std::cout << id_na << " " << id_ns_co << std::endl;
  
	unsigned int max_id_msh;
	{	// 要素IDの最大値を求める
		max_id_msh = 0;
		const std::vector<unsigned int>& aID = mesh.GetAry_ID();
		for(unsigned int iid=0;iid<aID.size();iid++){
			if( max_id_msh < aID[iid] ) max_id_msh = aID[iid];
		}
	}
  
  CIDConvEAMshCad conv;
	conv.m_aIdAry.reserve(max_id_msh+1);
	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	std::vector< std::pair<unsigned int, unsigned int> > aEaEs;
	{
		const std::vector<unsigned int>& aID = mesh.GetAry_ID();
		for(unsigned int iid=0;iid<aID.size();iid++){
			const unsigned int id_msh = aID[iid];
			std::vector<int> lnods;
			Msh::MSH_TYPE msh_type = mesh.GetConnectivity(id_msh,lnods);
			unsigned int id_cad_part, id_msh_before_ext, inum_ext;
			int ilayer;
      mesh.GetInfo(id_msh, id_cad_part, id_msh_before_ext, inum_ext, ilayer);
			Cad::CAD_ELEM_TYPE itype_cad_part;
			unsigned int nnoel = 0;
			ELEM_TYPE elem_type;
			if(      msh_type == Msh::HEX     ){ 
        nnoel = 8; elem_type = HEX;   
        if( inum_ext == 0 ){ itype_cad_part = Cad::SOLID; }
        else{ assert( inum_ext%2 == 0 ); itype_cad_part = Cad::LOOP; }
      }
			else if( msh_type == Msh::TET     ){ 
        nnoel = 4; elem_type = TET;
        if( inum_ext == 0 ){ itype_cad_part = Cad::SOLID; }
        else{ assert( inum_ext%2 == 0 ); itype_cad_part = Cad::LOOP; }
      }
			else if( msh_type == Msh::QUAD    ){ 
        nnoel = 4; elem_type = QUAD;  
        if( inum_ext == 0 ){        itype_cad_part = Cad::LOOP; }
        else if( inum_ext%2 == 0 ){ itype_cad_part = Cad::EDGE; }   
        else{                       itype_cad_part = Cad::LOOP; }
      }
			else if( msh_type == Msh::TRI     ){ 
        nnoel = 3; elem_type = TRI;   
        if( inum_ext == 0 ){        itype_cad_part = Cad::LOOP; }
        else if( inum_ext%2 == 0 ){ itype_cad_part = Cad::EDGE; } 
        else{                       itype_cad_part = Cad::LOOP; }
      }
			else if( msh_type == Msh::BAR     ){ 
        nnoel = 2; elem_type = LINE;  
        if( inum_ext == 0 ){        itype_cad_part = Cad::EDGE; }
        else if( inum_ext%2 == 0 ){ itype_cad_part = Cad::VERTEX; }   
        else{                       itype_cad_part = Cad::EDGE; }
      }
			else if( msh_type == Msh::VERTEX  ){ 
        nnoel = 1; elem_type = POINT; 
        if( inum_ext != 0 ){ assert( inum_ext % 2 == 1 ); }
        itype_cad_part = Cad::VERTEX; 
      }
			else{ assert(0); }
			assert( lnods.size() % nnoel == 0 );
			const unsigned int nelem = lnods.size()/nnoel;
			CElemAry* pea = new CElemAry(nelem,elem_type);
			const unsigned int id_es = pea->AddSegment(pea->GetFreeSegID(),CElemAry::CElemSeg(id_na,CORNER),lnods);
			const unsigned int id_ea = this->m_apEA.GetFreeObjID();
			const unsigned int tmp_id = this->m_apEA.AddObj( std::make_pair(id_ea,pea)  );
			assert( tmp_id == id_ea );
			aEaEs.push_back( std::make_pair(id_ea,id_es) );
			CIDConvEAMshCad::CInfoCadMshEA info;
			{
				info.id_ea = id_ea;
				info.id_part_msh = id_msh;
				info.id_part_cad = id_cad_part;
				info.itype_part_cad = itype_cad_part;
        info.id_part_msh_before_extrude = id_msh_before_ext;
        info.inum_extrude = inum_ext;
			}
			conv.m_aIdAry.push_back( info );
			////////////////
			CField::CElemInterpolation ei( id_ea,  0,id_es,  0,0,  0,0 );
			ei.ilayer = ilayer;
			aElemIntp.push_back(ei);
		}
	}
  
	{	// 包含関係(include relation)を作る
		CNodeAry& na = this->GetNA(id_na);
		for(unsigned int ieaes=0;ieaes<aEaEs.size();ieaes++){
			na.AddEaEs( aEaEs[ieaes] );
		}
		const std::vector<unsigned int>& id_ary_msh = mesh.GetAry_ID();
		for(unsigned int iid_ary_msh=0;iid_ary_msh<id_ary_msh.size();iid_ary_msh++){
			const unsigned int id_msh = id_ary_msh[iid_ary_msh];
			const unsigned int id_ea = conv.GetIdEA_fromMsh(id_msh); assert( this->IsIdEA(id_ea) );
			unsigned int iEaEs = 0;
			for(;iEaEs<aEaEs.size();iEaEs++){
				if( id_ea == aEaEs[iEaEs].first ) break;
			}
			assert( iEaEs != aEaEs.size() );
			std::vector<unsigned int> id_ary_msh_inc = mesh.GetIncludeElemIDAry(id_msh);
			for(unsigned int iid_ary_msh_inc=0;iid_ary_msh_inc<id_ary_msh_inc.size();iid_ary_msh_inc++){
				unsigned int id_msh_inc = id_ary_msh_inc[iid_ary_msh_inc];
				unsigned int id_ea_inc = conv.GetIdEA_fromMsh(id_msh_inc);
				assert( this->IsIdEA(id_ea_inc) );
				unsigned int iEaEs_inc;
				for(iEaEs_inc=0;iEaEs_inc<aEaEs.size();iEaEs_inc++){
					if( id_ea_inc == aEaEs[iEaEs_inc].first ) break;
				}
				assert( iEaEs != aEaEs.size() );
				na.SetIncludeEaEs_InEaEs( aEaEs[iEaEs_inc], aEaEs[iEaEs] );
			}
		}
	}
  
	unsigned int id_field_base = 0;
	{
		CField* pField = new CField( 0,	// 親フィールド
                                aElemIntp,	// 要素Index
                                CField::CNodeSegInNodeAry(id_na,false,id_ns_co,  0,false,0,0,0), // CORNER節点
                                CField::CNodeSegInNodeAry(),	// BUBBLE節点
                                *this );
		id_field_base = this->m_apField.AddObj( std::make_pair(0,pField) );
	}
	if( id_field_base == 0 ){ return 0; }
  m_map_field_conv.insert( std::make_pair(id_field_base, conv) );
	return id_field_base;  
}

unsigned int CFieldWorld::SetCustomBaseField(unsigned int id_base,
	std::vector<unsigned int> aIdEA_Inc,
	std::vector< std::vector<int> >& aLnods,
	std::vector<unsigned int>& mapVal2Co)
{
	Fem::Field::CField& field_base = this->GetField(id_base);
	const unsigned int nno_v = mapVal2Co.size();
	unsigned int id_na_v = this->m_apNA.GetFreeObjID();
	id_na_v = m_apNA.AddObj( std::make_pair(id_na_v,new CNodeAry(nno_v)) );
	unsigned int id_na_c = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_ns_c = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	for(unsigned int iidea=0;iidea<aIdEA_Inc.size();iidea++){
		const unsigned int id_ea0 = aIdEA_Inc[iidea];
		unsigned int id_es_c = field_base.GetIdElemSeg(id_ea0,CORNER,false,*this);
		Fem::Field::CElemAry& ea = this->GetEA(id_ea0);
		unsigned int id_es_v = ea.AddSegment(ea.GetFreeSegID(),CElemAry::CElemSeg(id_na_v,CORNER), aLnods[iidea] );
		aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea0, id_es_v,id_es_c, 0,0, 0,0) );
	}
	CField* pField = new CField( 0,	// 親フィールド
		aElemIntp,	// 要素Index
		CField::CNodeSegInNodeAry(id_na_c,false,id_ns_c,  id_na_v,false,0,0,0),	// CORNER補間節点
		CField::CNodeSegInNodeAry(), // BUBBLE補間節点
		*this );
	const unsigned int id_base_new = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_base_new;
}


bool CFieldWorld::UpdateConnectivity_CustomBaseField(const unsigned int id_base,
	const std::vector<unsigned int>& aIdEA_Inc, 
	const std::vector< std::vector<int> >& aLnods,
	const std::vector<unsigned int>& mapVal2Co)
{
	Fem::Field::CField& fb = this->GetField(id_base);
	const unsigned int id_na = fb.GetNodeSegInNodeAry(Fem::Field::CORNER).id_na_va;
	const std::vector<unsigned int>& aIdEA = fb.GetAryIdEA();
	assert( aIdEA_Inc.size() == aIdEA.size() );
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		unsigned int id_ea = aIdEA[iiea];
		assert( aIdEA_Inc[iiea] == id_ea );
    assert( this->IsIdEA(id_ea) );
		CElemAry& ea = this->GetEA(id_ea);
		assert( ea.ElemType() == Fem::Field::TRI );
    unsigned int id_es = 0;
		{
			const std::vector<unsigned int>& aIdES = ea.GetAry_SegID();
			for(unsigned int iies=0;iies<aIdES.size();iies++){
				unsigned int id_es0 = aIdES[iies];
				const Fem::Field::CElemAry::CElemSeg& es = ea.GetSeg(id_es0);
				if( id_na == es.GetIdNA() ){ 
					id_es = id_es0; 
					break;
				}
			}
		}
    assert( ea.IsSegID(id_es) );
		CElemAry::CElemSeg& es = ea.GetSeg(id_es);
		const std::vector<int>& lnods = aLnods[iiea];
		assert( lnods.size() == ea.Size() * 3 );
		for(unsigned int itri=0;itri<ea.Size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				int ino0 = lnods[itri*3+inotri];
				es.SetNodes(itri,inotri,ino0);
			}
    }
	}
	return true;  
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 節点配列に関係する関数群

bool CFieldWorld::IsIdNA( unsigned int id_na ) const
{
	if( this->m_apNA.IsObjID(id_na) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdNA() const 
{ 
	return this->m_apNA.GetAry_ObjID(); 
}

const CNodeAry& CFieldWorld::GetNA(unsigned int id_na) const
{
	assert( this->m_apNA.IsObjID(id_na) );
	if( !this->m_apNA.IsObjID(id_na) ) throw;
	return *m_apNA.GetObj(id_na);
}

CNodeAry& CFieldWorld::GetNA(unsigned int id_na)
{
	assert( this->m_apNA.IsObjID(id_na) );
	if( !this->m_apNA.IsObjID(id_na) ) throw;
	return *m_apNA.GetObj(id_na);
}

unsigned int CFieldWorld::AddNodeAry(unsigned int size)
{
	CNodeAry* pNA = new CNodeAry( size );
	return this->m_apNA.AddObj( std::make_pair(0,pNA) );
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 要素配列に関係する関数群

bool CFieldWorld::IsIdEA( unsigned int id_ea ) const{
	if( this->m_apEA.IsObjID(id_ea) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdEA() const { return this->m_apEA.GetAry_ObjID(); }

const CElemAry& CFieldWorld::GetEA(unsigned int id_ea) const{
	assert( this->m_apEA.IsObjID(id_ea) );
	if( !this->m_apEA.IsObjID(id_ea) ) throw;
	return *m_apEA.GetObj(id_ea);
}

CElemAry& CFieldWorld::GetEA(unsigned int id_ea){
	assert( this->m_apEA.IsObjID(id_ea) );
	if( !this->m_apEA.IsObjID(id_ea) ) throw;
	return *m_apEA.GetObj(id_ea);
}

unsigned int CFieldWorld::AddElemAry(unsigned int size, ELEM_TYPE elem_type){
	CElemAry* pEA = new CElemAry( size, elem_type );
	unsigned int id_ea = this->m_apEA.AddObj( std::make_pair(0,pEA) );
	return id_ea;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// Fieldに関係する関数群

bool CFieldWorld::IsIdField( unsigned int id_field ) const{
	if( this->m_apField.IsObjID(id_field) ){ return true; }
	return false;
}

const std::vector<unsigned int>& CFieldWorld::GetAry_IdField() const { return this->m_apField.GetAry_ObjID(); }

// Fieldの参照を得る関数（ポインタじゃなくて参照なのは、勝手にdeleteされないため）
const CField& CFieldWorld::GetField(unsigned int id_field) const{
	assert( this->m_apField.IsObjID(id_field) );
	if( !this->m_apField.IsObjID(id_field) ) throw "invalid field id";
	return *m_apField.GetObj(id_field);
}

CField& CFieldWorld::GetField(unsigned int id_field){
	assert( this->m_apField.IsObjID(id_field) );
	if( !this->m_apField.IsObjID(id_field) ) throw "invalid field id";
	return *m_apField.GetObj(id_field);
}

unsigned int CFieldWorld::AddField
(unsigned int id_field_parent,	// parent field
 const std::vector<CField::CElemInterpolation>& aEI, 
 const CField::CNodeSegInNodeAry& nsna_c, const CField::CNodeSegInNodeAry& nsna_b,
 unsigned int id_field_candidate)
{
	CField* pField = new CField
  (id_field_parent,
   aEI,
   nsna_c, nsna_b,
   *this );
  
	const unsigned int id_field = this->m_apField.AddObj( std::make_pair(id_field_candidate,pField) );
	assert( id_field != 0 );  
  return id_field;
}


unsigned int CFieldWorld::MakeField_FieldElemAry
(unsigned int id_field_base,                                               
 unsigned int id_ea, 
 Fem::Field::FIELD_TYPE field_type, const int fdt, const int nct )
{
	if( !this->IsIdEA(id_ea) ) return 0;
	assert( nct == CORNER );

	unsigned int id_es_val = 0;
	unsigned int id_na_val = 0;
	{   // 既存のNAで対応できるか調べる
		CElemAry& ea = this->GetEA(id_ea);
		const std::vector<unsigned int> aIdEs = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
			unsigned int id_es0 = aIdEs[iid_es];
			assert( ea.IsSegID(id_es0) );
			const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
			if( es0.GetElSegType() != CORNER ){ continue; }
			const unsigned int id_na = es0.GetIdNA();
			assert( this->IsIdNA(id_na) );
			const CNodeAry& na = this->GetNA(id_na);
			assert( es0.GetMaxNoes() < na.Size() );
            if( es0.GetMaxNoes()+1 != na.Size() ){ continue; }
            ////////////////
			const unsigned int nnode = na.Size();
			std::vector<unsigned int> flg_vec;
			flg_vec.resize(nnode,0);
			unsigned int noes[10];
			unsigned int nnoes = es0.Length();
			for(unsigned int ielem=0;ielem<ea.Size();ielem++){
				es0.GetNodes(ielem,noes);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){ flg_vec[ noes[inoes] ] = 1; }
			}
			bool iflg = true;
			for(unsigned int inode=0;inode<nnode;inode++){
				if( flg_vec[inode] == 0 ){ iflg = false; break; }
			}
			if( iflg == true ){	// 既存のセグメントで対応できる場合
//				std::cout << " Use Not New " << std::endl;
				id_es_val = id_es0;
				id_na_val = id_na;
				break;
			}
		}
		if( id_na_val == 0 && id_es_val == 0 ){	// 新しい値セグメントを追加
			unsigned int id_es1 = 0;
			for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
				unsigned int id_es0 = aIdEs[iid_es];
				assert( ea.IsSegID(id_es0) );
				const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
				if( es0.GetElSegType() == CORNER ){ 
					id_es1 = id_es0;
					break; 
				}
			}
			const CElemAry::CElemSeg es1 = ea.GetSeg(id_es1);
			unsigned int nnode;
			std::vector<int> flg_vec;
			{	// coからvaのテーブルであるflg_vecを作る
				const unsigned int mx_noes = es1.GetMaxNoes();
				flg_vec.resize(mx_noes+1,-2);
				unsigned int noes[10];
				unsigned int nnoes = es1.Length();
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es1.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						flg_vec[ noes[inoes] ] = -1;
					}
				}
				nnode = 0;
				for(unsigned int i=0;i<mx_noes+1;i++){
					if( flg_vec[i] == -1 ){
						flg_vec[i] = nnode;
						nnode++;
					}
				}
			}
//			std::cout << "New Nnode for Val : " << nnode << std::endl;
			id_na_val = this->AddNodeAry(nnode);
			std::vector<int> lnods;
			{	// lnodsを作る
				const unsigned int nnoes = es1.Length();
				const unsigned int nelem = ea.Size();
				lnods.resize( nnoes*nelem );;
				unsigned int noes[16];
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es1.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						int inode0 = noes[inoes];
						int inode1 = flg_vec[inode0];
						assert( inode1 >= 0 );
						lnods[ielem*nnoes+inoes] = inode1;
					}
				}
			}
			id_es_val = ea.AddSegment(ea.GetFreeSegID(),CElemAry::CElemSeg(id_na_val,CORNER),lnods);
			{	// 包含関係を入れる
				CNodeAry& na_val = this->GetNA(id_na_val);
				na_val.AddEaEs( std::make_pair(id_ea,id_es_val) );
			}
		}	// 新しい要素/節点セグメントを追加
	}

//	std::cout << "VAL   EA:" << id_ea << "  ES:" << id_es_val << "   NA:" << id_na_val << std::endl;

    assert( id_field_base != 0 );

	unsigned int ndim_coord;
	unsigned int id_na_co = 0;
	unsigned int id_ns_co = 0;
	bool is_part_co;
	unsigned int id_es_co = 0;
	{
		assert( this->IsIdField(id_field_base) );
		const CField& field_base = this->GetField(id_field_base);
		ndim_coord = field_base.GetNDimCoord();
		id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
		id_ns_co = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;
		CElemAry& ea = this->GetEA(id_ea);
		const std::vector<unsigned int> aIdEs = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<aIdEs.size();iid_es++){
			unsigned int id_es0 = aIdEs[iid_es];
			assert( ea.IsSegID(id_es0) );
			const CElemAry::CElemSeg es0 = ea.GetSeg(id_es0);
			if( es0.GetElSegType() != CORNER ){ continue; }
			if( id_na_co == es0.GetIdNA() ){
				id_es_co = id_es0;
			}
		}
		assert( id_es_co != 0 );
		////////////////
		assert( this->IsIdNA(id_na_co) );
		const CNodeAry& na = this->GetNA(id_na_co);
		const std::vector< std::pair<unsigned int,unsigned int> >& aEaEs = na.GetAry_EaEs_Min();
		if( aEaEs.size() == 1 && aEaEs[0].first == id_ea ){ is_part_co = false; }
		else{ is_part_co = true; }
	}

	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c.id_na_co = id_na_co;
		na_c.id_ns_co = id_ns_co;
		na_c.id_na_va = id_na_val;
		na_c.is_part_va = false;
		na_c.is_part_co = is_part_co;
	}

	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	{
		aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea, id_es_val,id_es_co, 0,0, 0,0) );
	}
  
  const unsigned int id_field = this->AddField(0,aElemIntp,na_c,Field::CField::CNodeSegInNodeAry());  

  assert( this->IsIdField(id_field) );  
  CField& field = this->GetField(id_field);
	field.SetValueType(field_type,fdt,*this);
  assert( field.AssertValid(*this) );
	return id_field;
}

unsigned int CFieldWorld::MakeField_FieldElemDim
(unsigned int id_field_base, 
 int idim_elem,  // 要素の次元（－１なら全ての要素）
 Fem::Field::FIELD_TYPE field_type, 
 const int fdt, const int nct )
{
  //	std::cout << "CFieldWorld::MakeField_Field" << std::endl;
	if( !(fdt&VALUE) && !(fdt&VELOCITY) && !(fdt&ACCELERATION) ) return 0;
	if( this->m_apNA.GetAry_ObjID().size() == 0 ) return 0;
	if( this->m_apEA.GetAry_ObjID().size() == 0 ) return 0;
  
	assert( this->IsIdField(id_field_base) );
	const CField& field_base = this->GetField(id_field_base);
  
	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
	Fem::Field::CField::CNodeSegInNodeAry na_c, na_b;
  
	na_c = field_base.GetNodeSegInNodeAry(CORNER);
	na_c.id_na_va = 0;
  
  {
    const std::vector<unsigned int>& aIdEA = field_base.GetAryIdEA();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
      const unsigned int id_ea = aIdEA[iiea];
      if( idim_elem != -1 ){
        const Fem::Field::CElemAry& ea = this->GetEA(id_ea);
        if( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) != idim_elem ) continue;
      }
      const unsigned int id_es_c_co = field_base.GetIdElemSeg(id_ea,CORNER,false,*this);
      aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea, 0,id_es_c_co, 0,0, 0,0 ) );
    }
  }
  
	if( nct & CORNER ){
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
      {
        Fem::Field::CElemAry& ea = this->GetEA(id_ea);
        assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
      }
      if( field_base.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
				aElemIntp[iei].id_es_c_va = field_base.GetIdElemSeg(id_ea,CORNER,true, *this);
			}
      else{
				aElemIntp[iei].id_es_c_va = field_base.GetIdElemSeg(id_ea,CORNER,false,*this);
      }
		}
		na_c.id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
		na_c.id_ns_co = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;
    na_c.is_part_va = false;
    if( field_base.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
      na_c.id_na_va = field_base.GetNodeSegInNodeAry(CORNER).id_na_va;
    }
    else{
      na_c.id_na_va = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
    }
	}
	if( nct &  BUBBLE ){
    // 節点配列の追加
    unsigned int nbubble = 0;
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
			Fem::Field::CElemAry& ea = this->GetEA(id_ea);
      assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
			nbubble += ea.Size();
		}
		unsigned int id_na = this->AddNodeAry(nbubble);
		assert( this->IsIdNA(id_na) );
		na_b.id_na_va = id_na;
    // 要素配列の追加
		unsigned int ibubble = 0;
		std::vector<int> lnods;
		for(unsigned int iei=0;iei<aElemIntp.size();iei++){
			const unsigned int id_ea = aElemIntp[iei].id_ea;
			CElemAry& ea = this->GetEA(id_ea);
      assert( (int)CElemAry::CElemSeg::GetElemDim( ea.ElemType() ) == idim_elem );
			{	// コネクティビティ(lnods)を作る
				lnods.resize(ea.Size());
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					lnods[ielem] = ibubble;
					ibubble++;
				}
			}
			std::vector<CElemAry::CElemSeg> es_ary;
//			const unsigned int id_es_b_va = ea.GetFreeSegID();
//			es_ary.push_back( CElemAry::CElemSeg(id_es_b_va,id_na,BUBBLE) );
      const unsigned int id_es_b_va = ea.AddSegment(ea.GetFreeSegID(),CElemAry::CElemSeg(id_na,BUBBLE),lnods);
//      assert( ares.size() == 1 && ares[0] == (int)id_es_b_va );
      aElemIntp[iei].id_es_b_va = id_es_b_va;
		}
		assert( ibubble == nbubble );
	}
	for(unsigned int iei=0;iei<aElemIntp.size();iei++){
		unsigned int id_ea = aElemIntp[iei].id_ea;
		aElemIntp[iei].ilayer = field_base.GetLayer(id_ea);
    //		std::cout << "layer " << aElemIntp[iei].ilayer << std::endl;
	}
  
  const unsigned int id_field= this->AddField(0,aElemIntp,na_c,na_b); 
	assert( id_field != 0 );
  assert( this->IsIdField(id_field) );
  CField& field = this->GetField(id_field);  
	field.SetValueType(field_type,fdt,*this);
  assert( field.AssertValid(*this) );
	return id_field;  
}

bool CFieldWorld::UpdateConnectivity_EdgeField_Tri(unsigned int id_field, unsigned int id_field_base)
{
	assert( this->IsIdField(id_field_base) );
	const CField& field_base = this->GetField(id_field_base);
	
//	unsigned int id_na = field_base.GetNodeSegInNodeAry(CORNER).id_na_va;

	std::vector<unsigned int> edge_ary;
	unsigned int nedge;
	{
		unsigned int id_ea0 = field_base.GetAryIdEA()[0];
		unsigned int id_es_co = field_base.GetIdElemSeg(id_ea0,CORNER,false,*this);
		const CElemAry& ea = this->GetEA(id_ea0);
		ea.MakeEdge(id_es_co,nedge,edge_ary);
		assert( edge_ary.size() == nedge*2 );
	}
	assert( this->IsIdField(id_field) );
	const CField& field = this->GetField(id_field);
	{
		unsigned int id_ea0 = field.GetAryIdEA()[0];
		CElemAry& ea0 = this->GetEA(id_ea0);
		unsigned int id_es0 = ea0.GetAry_SegID()[0];
		CElemAry::CElemSeg& es0 = ea0.GetSeg(id_es0);
		assert( es0.Size() == nedge );
		assert( es0.Length() == 2 );
		for(unsigned int iedge=0;iedge<nedge;iedge++){
			es0.SetNodes(iedge,0, edge_ary[iedge*2+0]);
			es0.SetNodes(iedge,1, edge_ary[iedge*2+1]);
		}
	}	
	return true;
}

bool CFieldWorld::UpdateConnectivity_HingeField_Tri( unsigned int id_field_edge, unsigned int id_field_base)
{
	assert( this->IsIdField(id_field_base) );
	const CField& field_base = this->GetField(id_field_base);
  
  const CField& field_edge = this->GetField(id_field_edge);


//	unsigned int id_na = field_base.GetNodeSegInNodeAry(CORNER).id_na_va;
  
  const std::vector<unsigned int>& aIdEA_base = field_base.GetAryIdEA();  
  
  for(unsigned int iiea=0;iiea<aIdEA_base.size();iiea++){
    unsigned int id_ea0 = aIdEA_base[iiea];
    std::vector<int> lnods;
    unsigned int nedge = 0;
    {
      unsigned int id_es_co = field_base.GetIdElemSeg(id_ea0,CORNER,false,*this);
      const CElemAry& ea = this->GetEA(id_ea0);	
      int* elsuel = new int [ea.Size()*3];
      ea.MakeElemSurElem(id_es_co,elsuel);
      
      {
        const CElemAry::CElemSeg& es = ea.GetSeg(id_es_co);
        for(unsigned int ielem=0;ielem<ea.Size();ielem++){
          unsigned int no[3];
          es.GetNodes(ielem,no);
          if( elsuel[ielem*3+0] >= 0 && no[1] < no[2] ){ nedge++; }
          if( elsuel[ielem*3+1] >= 0 && no[2] < no[0] ){ nedge++; }
          if( elsuel[ielem*3+2] >= 0 && no[0] < no[1] ){ nedge++; }
        }
      }
      {
        lnods.resize(nedge*4);
        unsigned int iedge=0;
        const CElemAry::CElemSeg& es = ea.GetSeg(id_es_co);
        for(unsigned int ielem=0;ielem<ea.Size();ielem++){
          unsigned int no[3];	es.GetNodes(ielem,no);
          if( elsuel[ielem*3+0] >= 0 && no[1] < no[2] ){
            lnods[iedge*4+0]=no[0];	lnods[iedge*4+2]=no[1];	lnods[iedge*4+3]=no[2];
            const unsigned int ielem0 = elsuel[ielem*3+0];
            unsigned int no0[3];	es.GetNodes(ielem0,no0);
            if(      no0[1]==no[2] && no0[2]==no[1] ){ lnods[iedge*4+1] = no0[0]; }
            else if( no0[2]==no[2] && no0[0]==no[1] ){ lnods[iedge*4+1] = no0[1]; }
            else if( no0[0]==no[2] && no0[1]==no[1] ){ lnods[iedge*4+1] = no0[2]; }
            else{ assert(0); }
            iedge++; 
          }
          if( elsuel[ielem*3+1] >= 0 && no[2] < no[0] ){
            lnods[iedge*4+0]=no[1];	lnods[iedge*4+2]=no[2];	lnods[iedge*4+3]=no[0];
            const unsigned int ielem0 = elsuel[ielem*3+1];
            unsigned int no0[3];	es.GetNodes(ielem0,no0);
            if(      no0[1]==no[0] && no0[2]==no[2] ){ lnods[iedge*4+1] = no0[0]; }
            else if( no0[2]==no[0] && no0[0]==no[2] ){ lnods[iedge*4+1] = no0[1]; }
            else if( no0[0]==no[0] && no0[1]==no[2] ){ lnods[iedge*4+1] = no0[2]; }
            else{ assert(0); }
            iedge++; 
          }
          if( elsuel[ielem*3+2] >= 0 && no[0] < no[1] ){
            lnods[iedge*4+0]=no[2];	lnods[iedge*4+2]=no[0];	lnods[iedge*4+3]=no[1];
            const unsigned int ielem0 = elsuel[ielem*3+2];
            unsigned int no0[3];	es.GetNodes(ielem0,no0);
            if(      no0[1]==no[1] && no0[2]==no[0] ){ lnods[iedge*4+1] = no0[0]; }
            else if( no0[2]==no[1] && no0[0]==no[0] ){ lnods[iedge*4+1] = no0[1]; }
            else if( no0[0]==no[1] && no0[1]==no[0] ){ lnods[iedge*4+1] = no0[2]; }
            else{ assert(0); }
            iedge++; 
          }
        }
        assert( iedge == nedge );
      }
//      std::cout << "nedge : " << nedge << std::endl;
    }
    
    {
      unsigned int id_ea1 = field_edge.GetAryIdEA()[iiea];
      CElemAry& ea1 = this->GetEA(id_ea1);
      unsigned int id_es_co = field_edge.GetIdElemSeg(id_ea1,CORNER,false,*this);
      assert( ea1.IsSegID(id_es_co) );
      CElemAry::CElemSeg& es1 = ea1.GetSeg(id_es_co);
      assert( es1.Size() == nedge );
      assert( es1.Length() == 4 );
      for(unsigned int iedge=0;iedge<nedge;iedge++){
        es1.SetNodes(iedge,0, lnods[iedge*4+0]);
        es1.SetNodes(iedge,1, lnods[iedge*4+1]);
        es1.SetNodes(iedge,2, lnods[iedge*4+2]);
        es1.SetNodes(iedge,3, lnods[iedge*4+3]);
      }
    }    
  }
    
	return true;
}


/*
unsigned int CFieldWorld::GetPartialField(unsigned int id_field_val, unsigned int id_ea)
{
//	std::cout << "GetPartialField" << std::endl;
	if( !this->IsIdEA(id_ea) ){ return 0; }

	if( !this->IsIdField(id_field_val) ){ return 0; }
	const CField& field_val = this->GetField(id_field_val);
	if( !field_val.IsValid() ){ return 0; }
	if( !(field_val.GetIDFieldParent()==0) ){ return 0; }	// 親フィールドでなければならない．

	// 既に作られているPartialFieldで適合するものを探す
	{
		const std::vector<unsigned int> aIdField = this->m_apField.GetAry_ObjID();
		for(unsigned int iid_f=0;iid_f<aIdField.size();iid_f++){
			unsigned int id_f = aIdField[iid_f];
			const Fem::Field::CField& field = *m_apField.GetObj(id_f);
			if( field.GetIDFieldParent() != id_field_val ) continue;
			const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
			if( aIdEA.size() != 1 ) continue;
			if( aIdEA[0] == id_ea ){ return id_f; }
		}
	}

	{
		bool iflag = false;
		const unsigned int id_na_va_c = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		assert( this->IsIdEA(id_na_va_c) );
		const CNodeAry& na_va = this->GetNA(id_na_va_c);
		const unsigned int id_na_co_c = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
		assert( this->IsIdEA(id_na_co_c) );
		const CNodeAry& na_co = this->GetNA(id_na_co_c);
		const std::vector<unsigned int>& aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea0 = aIdEA[iiea];
			unsigned int id_es0 = field_val.GetIdElemSeg(id_ea0,CORNER,true,*this);
//			std::cout << "Try : " << id_ea0 << " " << id_es0 << " " << id_ea << std::endl;
			const unsigned int id_es1 = na_va.IsContainEa_InEaEs( std::make_pair(id_ea0,id_es0), id_ea );
			if( id_es1 != 0 ){
				iflag = true;
				break;
			}
			unsigned int id_es2 = field_val.GetIdElemSeg(id_ea0,CORNER,false,*this);
//			std::cout << "Try : " << id_ea0 << " " << id_es2 << " " << id_ea << std::endl;
			const unsigned int id_es3 = na_co.IsContainEa_InEaEs( std::make_pair(id_ea0,id_es2), id_ea );
			if( id_es3 != 0 ){
				iflag = true;
				break;
			}
		}
		if( !iflag ){ 
            std::cout << " not partial field " << std::endl;
            return 0; 
        }
	}

	Field::CField::CElemInterpolation ei(id_ea, 0,0, 0,0, 0,0);
	{	// corner節点の追加
		unsigned int id_na_c_co  = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
		unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		assert( this->IsIdEA(id_ea) );
		const CElemAry& ea =this->GetEA(id_ea);
		const std::vector<unsigned int>& id_es_ary = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<id_es_ary.size();iid_es++){
			const unsigned int id_es = id_es_ary[iid_es];
			assert( ea.IsSegID(id_es) );
			const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
			const unsigned int id_na0 = es.GetIdNA();
			if( id_na0 == id_na_c_val ){ ei.id_es_c_va = id_es; }
			if( id_na0 == id_na_c_co  ){ ei.id_es_c_co = id_es; }
		}
	}
	assert( ei.id_es_c_co != 0 );
	if( ei.id_es_c_va == 0 ){	// この値場は新しく作られたものなので、部分場に新しい要素配列を追加する
		unsigned int id_es_co = ei.id_es_c_co;
		assert( this->IsIdEA(id_ea) );
		CElemAry& ea =this->GetEA(id_ea);
		assert( ea.IsSegID(id_es_co) );
		const CElemAry::CElemSeg& es_co = ea.GetSeg(id_es_co);
		std::vector<int> lnods;
		{
			std::vector<int> map_co2va;
			{
				const CNodeAry& na_co = this->GetNA( field_val.GetNodeSegInNodeAry(CORNER).id_na_co );
				unsigned int nnode_co = na_co.Size();
				map_co2va.resize( nnode_co, -1);
				const CNodeAry& na_va = this->GetNA( field_val.GetNodeSegInNodeAry(CORNER).id_na_va );
				unsigned int nnode_va = na_va.Size();
				for(unsigned int inode_va=0;inode_va<nnode_va;inode_va++){
					unsigned int inode_co = field_val.GetMapVal2Co(inode_va);
					map_co2va[inode_co] = inode_va;
				}
			}
			unsigned int nelem = ea.Size();
			unsigned int nnoes = es_co.GetSizeNoes();
			lnods.resize( nelem*nnoes );
			unsigned int noes[16];
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				es_co.GetNodes(ielem,noes);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					unsigned int inode_co = noes[inoes];
					int inode_va = map_co2va[inode_co];
					assert( inode_va != -1 );   // これが失敗した場合は，このEAは場の一部ではない
					lnods[ielem*nnoes+inoes] = inode_va;
				}
			}
		}
		unsigned int id_es_va = ea.GetFreeSegID();
		const unsigned int id_na_va = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		ea.AddSegment(CElemAry::CElemSeg(id_es_va,id_na_va,CORNER),lnods);
		ei.id_es_c_va = id_es_va;
	}

	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	{
		aElemIntp.push_back(ei);
	}

	
	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c = field_val.GetNodeSegInNodeAry(CORNER);
		na_c.is_part_co = true;
		na_c.is_part_va = true;
	}

	Field::CField::CNodeSegInNodeAry na_b;
	// TODO:部分場のBUBBLE節点が含まれるかどうかを判断しなければならない．
	// とりあえず今は境界条件のためのPartialFieldと思って保留
//	{
//		na_b = field_val.GetNodeSegInNodeAry(BUBBLE);
//		na_b.is_part_co = true;
//		na_b.is_part_va = true;
//	}

	// ここに辺節点の追加に関する物を書く
	CField* pField = new CField(
		id_field_val,
		aElemIntp,
		na_c, na_b,
		*this );

	unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_field;
}
*/


unsigned int CFieldWorld::GetPartialField(unsigned int id_field_val, unsigned int id_ea)
{
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->GetPartialField(id_field_val,aIdEA);
}

unsigned int CFieldWorld::GetPartialField(unsigned int id_field_val, 
										  std::vector<unsigned int> id_ea_ary)
{
    if( id_ea_ary.size() == 0 ) return 0;
//	std::cout << "GetPartialField" << std::endl;
	for(unsigned int iid_ea=0;iid_ea<id_ea_ary.size();iid_ea++){
		unsigned int id_ea = id_ea_ary[iid_ea];
		if( !this->IsIdEA(id_ea) ){	return 0; }
	}
	if( !this->IsIdField(id_field_val) ) return 0;
	const Fem::Field::CField& field_val = this->GetField(id_field_val);
	if( !field_val.IsValid() ){ return 0; }
	assert( field_val.AssertValid(*this) );

	std::vector< Field::CField::CElemInterpolation > aElemIntp;
	unsigned int id_na_c_co = field_val.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;

	std::vector<int> map_co2val;
	if( id_na_c_co != id_na_c_val ){
		const CNodeAry& na_co = this->GetNA(id_na_c_co);
		map_co2val.resize( na_co.Size(), -1 );
		const CNodeAry& na_va = this->GetNA(id_na_c_val);
		const unsigned int npoin_val = na_va.Size();
		for(unsigned int ipoin=0;ipoin<npoin_val;ipoin++){
			unsigned int ipoin_co0 = field_val.GetMapVal2Co(ipoin);
            assert( ipoin_co0 < map_co2val.size() );
			map_co2val[ipoin_co0] = ipoin;
		}
	}
	for(unsigned int iid_ea=0;iid_ea<id_ea_ary.size();iid_ea++){
		unsigned int id_ea = id_ea_ary[iid_ea];
		Field::CField::CElemInterpolation ei(id_ea, 0,0, 0,0, 0,0);
		assert( this->IsIdEA(id_ea) );
		CElemAry& ea =this->GetEA(id_ea);
		const std::vector<unsigned int>& id_es_ary = ea.GetAry_SegID();
		for(unsigned int iid_es=0;iid_es<id_es_ary.size();iid_es++){
			const unsigned int id_es = id_es_ary[iid_es];
			assert( ea.IsSegID(id_es) );
			const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
			if( es.GetElSegType() != CORNER ) continue;
			const unsigned int id_na0 = es.GetIdNA();
			if( id_na0 == id_na_c_val ){ ei.id_es_c_va = id_es; }
			if( id_na0 == id_na_c_co  ){ ei.id_es_c_co = id_es; }
		}
		assert( ei.id_es_c_co != 0 );
		if( ei.id_es_c_va == 0 ){
			assert( map_co2val.size() != 0 );
			const CElemAry::CElemSeg& es_co = ea.GetSeg( ei.id_es_c_co );
			std::vector<int> lnods;
			const unsigned int nelem = ea.Size();
			const unsigned int nnoes = es_co.Length();
			lnods.resize(nelem*nnoes);
			unsigned int noes_co[32];
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				es_co.GetNodes(ielem,noes_co);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					unsigned int inode_co = noes_co[inoes];
                    assert( inode_co < map_co2val.size() );
					int inode_va = map_co2val[inode_co];
					assert( inode_va >= 0 );
					lnods[ielem*nnoes+inoes] = inode_va;
				}
			}
			ei.id_es_c_va = ea.AddSegment(ea.GetFreeSegID(),CElemAry::CElemSeg(id_na_c_val,CORNER), lnods);
		}
		aElemIntp.push_back( ei );
	}

	Field::CField::CNodeSegInNodeAry na_c;
	{
		na_c = field_val.GetNodeSegInNodeAry(CORNER);
		na_c.is_part_co = true;
		na_c.is_part_va = true;
	}

	Field::CField::CNodeSegInNodeAry na_b;
	// TODO:部分場のBUBBLE節点が含まれるかどうかを判断しなければならない．
	// とりあえず今は境界条件のためのPartialFieldと思って保留
/*	{
		na_b = field_val.GetNodeSegInNodeAry(BUBBLE);
		na_b.
		na_b.is_part_co = true;
		na_b.is_part_va = true;
	}*/

	// ここに辺節点の追加に関する物を書く
	CField* pField = new CField(
		id_field_val,
		aElemIntp,
		na_c,  na_b,
		*this
	);

	unsigned int id_field = this->m_apField.AddObj( std::make_pair(0,pField) );
	return id_field;
}

/*
void CFieldWorld::FieldValueExec(double time){
	const std::vector<unsigned int>& id_field_ary = this->m_apField.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_field_ary.size();iid++){
		unsigned int id_field = id_field_ary[iid];
		CField* pField = m_apField.GetObj(id_field);
		if( !pField->IsDepend() ){
			pField->ExecuteValue(time,*this);
		}
	}
}
*/
/*
void CFieldWorld::FieldValueDependExec(){
	const std::vector<unsigned int>& id_field_ary = this->m_apField.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_field_ary.size();iid++){
		unsigned int id_field = id_field_ary[iid];
		CField* pField = m_apField.GetObj(id_field);
		if( pField->IsDepend() ){
			pField->ExecuteValue(0.0,*this);
		}
	}
}
 */

void InsertIdToSetSet(std::map<unsigned int, std::set<unsigned int> >& mapset, unsigned int iv1, unsigned int iv2){
  if( iv1 == 0 || iv2  == 0 ) return;
  std::map<unsigned int, std::set<unsigned int> >::iterator itrmap = mapset.find(iv1);
  if( itrmap == mapset.end() ){ 
    std::set<unsigned int> set1;
    set1.insert(iv2);
    mapset.insert( std::make_pair(iv1,set1) ); 
  }
  else{
    itrmap->second.insert(iv2);
  }
}


// Delete Field and referenced EA and NA. the EAs and NAs that is referenced from Field in use are not deleted
void CFieldWorld::DeleteField( const std::vector<unsigned int>& aIdFieldDel )
{
//  std::cout << "delete field" << std::endl;
  std::vector<unsigned int> aIdField_InUse;
  {
    const std::vector<unsigned int>& aIdField = this->m_apField.GetAry_ObjID();  
    for(unsigned int iidf=0;iidf<aIdField.size();iidf++){
      const unsigned int id_field = aIdField[iidf];
      assert(this->IsIdField(id_field));
      bool is_use_flg = true;
      for(unsigned int iidf_del=0;iidf_del<aIdFieldDel.size();iidf_del++){
        if( aIdFieldDel[iidf_del] != id_field ) continue;
        is_use_flg = false;
        break;
      }
      if( is_use_flg ){ aIdField_InUse.push_back(id_field); }
    }
  }
  
  typedef std::map<unsigned int, std::set<unsigned int> > mapsetID;

  mapsetID EaEs_InUse, NaNs_InUse;
  for(unsigned int iidf_use=0;iidf_use<aIdField_InUse.size();iidf_use++){
    unsigned int id_field_use = aIdField_InUse[iidf_use];
//    std::cout << "Id Field In Use : " << id_field_use << std::endl;
    const CField& field = this->GetField(id_field_use);
    {
      const Fem::Field::CField::CNodeSegInNodeAry& nans = field.GetNodeSegInNodeAry(CORNER);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_co, nans.id_ns_co);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_va);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_ve);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_ac);          
    }
    {
      const Fem::Field::CField::CNodeSegInNodeAry& nans = field.GetNodeSegInNodeAry(BUBBLE);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_co, nans.id_ns_co);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_va);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_ve);
      InsertIdToSetSet(NaNs_InUse, nans.id_na_va, nans.id_ns_ac);
    }
    const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
    for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
      const unsigned int id_ea = aIdEA[iid_ea];
      InsertIdToSetSet(EaEs_InUse, id_ea, field.GetIdElemSeg(id_ea,CORNER,false,*this) );
      InsertIdToSetSet(EaEs_InUse, id_ea, field.GetIdElemSeg(id_ea,CORNER,true, *this) ); 
      InsertIdToSetSet(EaEs_InUse, id_ea, field.GetIdElemSeg(id_ea,BUBBLE,true, *this) );
      InsertIdToSetSet(EaEs_InUse, id_ea, field.GetIdElemSeg(id_ea,EDGE,  true, *this) );
    }
  }

  for(mapsetID::iterator itr=NaNs_InUse.begin();itr!=NaNs_InUse.end();itr++){
    unsigned int id_na = itr->first;
//    std::cout << "Id NA in use" << id_na << "  NS : ";
    const std::set<unsigned int>& setID = itr->second;
    for(std::set<unsigned int>::iterator itrID=setID.begin();itrID!=setID.end();itrID++){
//      std::cout << *itrID << " ";
    }
//    std::cout << std::endl;
  }
  for(mapsetID::iterator itr=EaEs_InUse.begin();itr!=EaEs_InUse.end();itr++){
    unsigned int id_ea = itr->first;
//    std::cout << "Id EA in use" << id_ea << "  ES : ";
    const std::set<unsigned int>& setID = itr->second;
    for(std::set<unsigned int>::iterator itrID=setID.begin();itrID!=setID.end();itrID++){
//      std::cout << *itrID << " ";
    }
//    std::cout << std::endl;
  }  
  {
    std::vector<unsigned int> aIdNA = this->m_apNA.GetAry_ObjID();
    for(unsigned int iina=0;iina<aIdNA.size();iina++){
      const unsigned int id_na = aIdNA[iina];
      if( NaNs_InUse.find(id_na) != NaNs_InUse.end() ) continue;
      CNodeAry* pNA = m_apNA.GetObj(id_na);
      delete pNA;
      m_apNA.DeleteObj(id_na);
    }
  }
  {
    std::vector<unsigned int> aIdEA = this->m_apEA.GetAry_ObjID();
    for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
      const unsigned int id_ea = aIdEA[iiea];
      if( EaEs_InUse.find(id_ea) != EaEs_InUse.end() ) continue;
      CElemAry* pEA = m_apEA.GetObj(id_ea);
      delete pEA;
      m_apEA.DeleteObj(id_ea);
    }
  }
  for(unsigned int iidf_del=0;iidf_del<aIdFieldDel.size();iidf_del++){
    const unsigned int id_field_del = aIdFieldDel[iidf_del];
    if( m_map_field_conv.find(id_field_del) != m_map_field_conv.end() ){
      std::map<unsigned int,CIDConvEAMshCad>::iterator itr = m_map_field_conv.find(id_field_del);
      m_map_field_conv.erase(itr);
    }
    CField* pField = m_apField.GetObj(id_field_del);
    delete pField;
    m_apField.DeleteObj(id_field_del);
  }  
}


