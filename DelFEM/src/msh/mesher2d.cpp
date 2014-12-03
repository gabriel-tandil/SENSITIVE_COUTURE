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
#pragma warning(disable: 4786)
#pragma warning(disable: 4996)
#endif

#define for if(0);else for

#include <stdio.h>
#include <set>
#include <vector>
#include <stack>
#include <cassert>
#include <math.h>
#include <cstdlib> //(abort)

#include <iostream>
#include <fstream>

//#include "../Cad/CadObj2D.h"

#include "delfem/msh/meshkernel2d.h"
#include "delfem/mesher2d.h"

using namespace Msh;
using namespace Com;

static const char indexRot3[3][3] = { 
	{ 0, 1, 2 },
	{ 1, 2, 0 },
	{ 2, 0, 1 },
};

void CMesher2D::Clear()
{
	setIdLCad_CutMesh.clear();
	this->m_imode_meshing = 1;
	this->m_elen = 0.1;
	this->m_esize = 1000;
	this->ClearMeshData();
}

CMesher2D::CMesher2D(const CMesher2D& rhs)
{
  this->Clear();
	setIdLCad_CutMesh = rhs.setIdLCad_CutMesh;
	m_imode_meshing = rhs.m_imode_meshing;
	m_elen = rhs.m_elen;
	m_esize = rhs.m_esize;
  
	m_ElemType = rhs.m_ElemType;
	m_ElemLoc = rhs.m_ElemLoc;
	m_include_relation = rhs.m_include_relation;
  
	m_aVertex = rhs.m_aVertex;
	m_aBarAry = rhs.m_aBarAry;
	m_aTriAry = rhs.m_aTriAry;
	m_aQuadAry = rhs.m_aQuadAry;
  
	aVec2D = rhs.aVec2D;
}


unsigned int CMesher2D::FindMaxID() const
{
	unsigned int max_id = 0;
	{	// get the maximum id
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			if( max_id < m_aVertex[iver].id ) max_id = m_aVertex[iver].id;
		}
		for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
			if( max_id < m_aBarAry[ibarary].id ) max_id = m_aBarAry[ibarary].id;
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			if( max_id < m_aTriAry[itriary].id ) max_id = m_aTriAry[itriary].id;
		}
		for(unsigned int iquadary=0;iquadary<m_aQuadAry.size();iquadary++){
			if( max_id < m_aQuadAry[iquadary].id ) max_id = m_aQuadAry[iquadary].id;
		}
	}
	return max_id;
}

unsigned int CMesher2D::GetFreeObjID(){
	unsigned int max_id = this->FindMaxID();
	std::vector<unsigned int> is_used_flg_ary;
	{	// このIDが使われているかどうかを示すハッシュを作る
		is_used_flg_ary.resize( max_id+1, 0 );
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			assert( is_used_flg_ary[ m_aVertex[iver].id ] == 0 );
			assert( m_aVertex[iver].id >= 1 && m_aVertex[iver].id <= max_id );
			is_used_flg_ary[ m_aVertex[iver].id ] = 1;
		}
		for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
			assert( is_used_flg_ary[ m_aBarAry[ibarary].id ] == 0 );
			assert( m_aBarAry[ibarary].id >= 1 && m_aBarAry[ibarary].id <= max_id );
			is_used_flg_ary[ m_aBarAry[ibarary].id ] = 1;
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			assert( is_used_flg_ary[ m_aTriAry[itriary].id ] == 0 );
			assert( m_aTriAry[itriary].id >= 1 && m_aTriAry[itriary].id <= max_id );
			is_used_flg_ary[ m_aTriAry[itriary].id ] = 1;
		}
		for(unsigned int iquadary=0;iquadary<m_aQuadAry.size();iquadary++){
			assert( is_used_flg_ary[ m_aQuadAry[iquadary].id ] == 0 );
			assert( m_aQuadAry[iquadary].id >= 1 && m_aQuadAry[iquadary].id <= max_id );
			is_used_flg_ary[ m_aQuadAry[iquadary].id ] = 1;
		}
	}
	assert( is_used_flg_ary[0] == 0 );
	for(unsigned int i=1;i<is_used_flg_ary.size();i++){
		if( is_used_flg_ary[i] == 0 ){ return i; }
	}
	return max_id+1;
}

bool CMesher2D::IsID(unsigned int id) const
{
//	assert( this->m_ElemLoc.size() > id );
	if( this->m_ElemLoc.size() <= id ) return false;
	const int iloc = this->m_ElemLoc[id];
	if( iloc == -1 ) return false;

	// 以下assertルーティン
	assert( this->m_ElemType.size() > id );
	const int itype = this->m_ElemType[id];
	assert( itype >= 0 );
	assert( iloc >= 0 );
	if( itype == 0 ){
		assert( m_aVertex.size() > (unsigned int)iloc );
		assert( m_aVertex[iloc].id == id );
	}
	else if( itype == 1 ){
		assert( m_aBarAry.size() > (unsigned int)iloc );
		assert( m_aBarAry[iloc].id == id );
	}
	else if( itype == 2 ){
		assert( m_aTriAry.size() > (unsigned int)iloc );
		assert( m_aTriAry[iloc].id == id );
	}
	else if( itype == 3 ){
		assert( m_aQuadAry.size() > (unsigned int)iloc );
		assert( m_aQuadAry[iloc].id == id );
	}
	return true;
}


bool CMesher2D::GetMshInfo(unsigned int id, unsigned int& nelem, MSH_TYPE& msh_type, unsigned int& iloc, unsigned int& id_cad ) const
{
	if( !this->IsID(id) ) return false;
	assert( id < this->m_ElemLoc.size() );
	assert( id < this->m_ElemType.size() );
	const int itype = this->m_ElemType[id];
	iloc = this->m_ElemLoc[id];
	assert( itype >= 0 );
	switch(itype){
	case 0:
		{
            assert( iloc < m_aVertex.size() );
			const SVertex& vertex = this->m_aVertex[iloc];
			id_cad = vertex.id_v_cad;
			nelem = 1;
			msh_type = VERTEX;
			break;
		}
	case 1:
		{
            assert( iloc < m_aBarAry.size() );
			const CBarAry& aBar = this->m_aBarAry[iloc];
			id_cad = aBar.id_e_cad;
			nelem = aBar.m_aBar.size();
			msh_type = BAR;
			break;
		}
	case 2:
		{
            assert( iloc < m_aTriAry.size() );
			const CTriAry2D& aTri = this->m_aTriAry[iloc];
			id_cad = aTri.id_l_cad;
			nelem = aTri.m_aTri.size();
			msh_type = TRI;
			break;
		}
	case 3:
		{
            assert( iloc < m_aQuadAry.size() );
			const CQuadAry2D& aQuad = this->m_aQuadAry[iloc];
			id_cad = aQuad.id_l_cad;
			nelem = aQuad.m_aQuad.size();
			msh_type = QUAD;
			break;
		}
	default:
		assert(0);
		abort();
	}
	return true;
}


MSH_TYPE CMesher2D::GetConnectivity(unsigned int id_msh,std::vector<int>& lnods) const
{
	assert( this->IsID(id_msh) );
	MSH_TYPE msh_type;
	unsigned int nnoel, nelem;
	const int itype = m_ElemType[id_msh];
	const int iloc = m_ElemLoc[id_msh];
	assert( itype != -1 && iloc != -1 );
	if(      itype == 0 ){ 
		msh_type = VERTEX; nnoel = 1; nelem = 1;
		lnods.resize(nnoel*nelem);
		lnods[0] = m_aVertex[iloc].v;
	}
	else if( itype == 1 ){ 
		msh_type = BAR;    nnoel = 2; 
		const std::vector<SBar>& aBar = m_aBarAry[iloc].m_aBar;
		nelem = aBar.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*2  ] = aBar[ielem].v[0];
			lnods[ielem*2+1] = aBar[ielem].v[1];
		}
	}
	else if( itype == 2 ){ 
		msh_type = TRI;    nnoel = 3;
		const std::vector<STri2D>& aTri = m_aTriAry[iloc].m_aTri;
		nelem = aTri.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*3  ] = aTri[ielem].v[0];
			lnods[ielem*3+1] = aTri[ielem].v[1];
			lnods[ielem*3+2] = aTri[ielem].v[2];
		}
	}
	else if( itype == 3 ){ 
		msh_type = QUAD;   nnoel = 4;
		const std::vector<SQuad2D>& aQuad = m_aQuadAry[iloc].m_aQuad;
		nelem = aQuad.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*4  ] = aQuad[ielem].v[0];
			lnods[ielem*4+1] = aQuad[ielem].v[1];
			lnods[ielem*4+2] = aQuad[ielem].v[2];
			lnods[ielem*4+3] = aQuad[ielem].v[3];
		}
	}
	else{ assert(0); }
	return msh_type;
}


// CADのIDからMshの要素IDを見つける関数
// puclic関数
unsigned int CMesher2D::GetElemID_FromCadID(unsigned int id_cad, Cad::CAD_ELEM_TYPE itype_cad)	const
{
	switch(itype_cad){
	case Cad::VERTEX:
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			if( m_aVertex[iver].id_v_cad == id_cad ){
				unsigned int id_msh = m_aVertex[iver].id;
                assert( m_ElemLoc[id_msh] == (int)iver );
				assert( m_ElemType[id_msh] == 0 );
				return id_msh;
			}
		}
		break;
	case Cad::EDGE:
		for(unsigned int ibar=0;ibar<m_aBarAry.size();ibar++){
			if( m_aBarAry[ibar].id_e_cad == id_cad ){
				unsigned int id_msh = m_aBarAry[ibar].id;
                assert( m_ElemLoc[id_msh] == (int)ibar );
				assert( m_ElemType[id_msh] == 1 );
				return id_msh;
			}
		}
		break;
	case Cad::LOOP:
		for(unsigned int itri=0;itri<m_aTriAry.size();itri++){
			if( m_aTriAry[itri].id_l_cad == id_cad ){ 
				unsigned int id_msh = m_aTriAry[itri].id;
                assert( m_ElemLoc[id_msh] == (int)itri );
				assert( m_ElemType[id_msh] == 2 );
				return id_msh;
			}
		}
		for(unsigned int iquad=0;iquad<m_aQuadAry.size();iquad++){
			if( m_aQuadAry[iquad].id_l_cad == id_cad ){ 
				unsigned int id_msh = m_aTriAry[iquad].id;
                assert( m_ElemLoc[id_msh] == (int)iquad );
				assert( m_ElemType[id_msh] == 3 );
				return id_msh;
			}
		}
		break;
	default:
		assert(0);
		return 0;
	}
	return 0;
}

// 要素のIndexと種類をCADのIDと種類から見つける
// private関数
bool CMesher2D::FindElemLocType_CadIDType
(unsigned int& iloc, unsigned int& itype, 
 unsigned int id_cad, Cad::CAD_ELEM_TYPE itype_cad )
{
	switch(itype_cad){
	case Cad::VERTEX:
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			if( m_aVertex[iver].id_v_cad == id_cad ){
				iloc = iver;
				itype = 0;
				unsigned int id_msh = m_aVertex[iver].id;
        assert( m_ElemLoc[id_msh]  == (int)iloc );
        assert( m_ElemType[id_msh] == (int)itype );
				return true;
			}
		}
		break;
	case Cad::EDGE:
		for(unsigned int ibar=0;ibar<m_aBarAry.size();ibar++){
			if( m_aBarAry[ibar].id_e_cad == id_cad ){
				iloc = ibar;
				itype = 1; 
				unsigned int id_msh = m_aBarAry[ibar].id;
        assert( m_ElemLoc[id_msh]  == (int)iloc );
                assert( m_ElemType[id_msh] == (int)itype );
				return true;
			}
		}
		break;
	case Cad::LOOP:
		for(unsigned int itri=0;itri<m_aTriAry.size();itri++){
			if( m_aTriAry[itri].id_l_cad == id_cad ){ 
				iloc = itri;
				itype = 2;
				unsigned int id_msh = m_aTriAry[itri].id;
        assert( m_ElemLoc[id_msh]  == (int)iloc );
        assert( m_ElemType[id_msh] == (int)itype );
				return true;
			}
		}
		for(unsigned int iquad=0;iquad<m_aQuadAry.size();iquad++){
			if( m_aQuadAry[iquad].id_l_cad == id_cad ){ 
				iloc = iquad;
				itype = 3;
				unsigned int id_msh = m_aQuadAry[iquad].id;
        assert( m_ElemLoc[id_msh]  == (int)iloc );
        assert( m_ElemType[id_msh] == (int)itype );
				return true;
			}
		}
		break;
	default:
		assert(0);
		break;
	}
	return false;
}

// 要素の場所と種類をハッシュ（IDが引数）している配列を初期化
void CMesher2D::MakeElemLocationType()
{
	unsigned int max_id = this->FindMaxID();
	////////////////
	this->m_ElemLoc.clear();
	this->m_ElemLoc.resize(max_id+1,-1);
	this->m_ElemType.clear();
	this->m_ElemType.resize(max_id+1,-1);
	////////////////
	for(unsigned int iver=0;iver<m_aVertex.size();iver++){
		unsigned int id0 = m_aVertex[iver].id;
		m_ElemLoc[id0] = iver;
		m_ElemType[id0] = 0;
	}
	for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
		unsigned int id0 = m_aBarAry[ibarary].id;
		m_ElemLoc[id0] = ibarary;
		m_ElemType[id0] = 1;
	}
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		unsigned int id0 = m_aTriAry[itriary].id;
		m_ElemLoc[id0] = itriary;
		m_ElemType[id0] = 2;
	}
	for(unsigned int iquadary=0;iquadary<m_aQuadAry.size();iquadary++){
		unsigned int id0 = m_aQuadAry[iquadary].id;
		m_ElemLoc[id0] = iquadary;
		m_ElemType[id0] = 3;
	}
}


// include_relationを作る
void CMesher2D::MakeIncludeRelation(const Cad::ICad2D_Msh& cad)
{
  this->m_include_relation.clear();
  
	if( m_ElemLoc.size() == 0 || m_ElemType.size() == 0 ) return;

	unsigned int max_id = this->FindMaxID();
	this->m_include_relation.resize( max_id+1 );
	////////////////
	std::vector<unsigned int> id_l_ary = cad.GetAryElemID(Cad::LOOP);
	for(unsigned int iid_l=0;iid_l<id_l_ary.size();iid_l++){
		unsigned int id_l = id_l_ary[iid_l];
		unsigned int id_atri = this->GetElemID_FromCadID(id_l,Cad::LOOP);
		if( !this->IsID(id_atri) ) continue;
    std::auto_ptr<Cad::IItrLoop> pItr_l = cad.GetPtrItrLoop(id_l);
		for(;;){
			for(;!pItr_l->IsEnd();(*pItr_l)++){
				unsigned int id_v_cad = pItr_l->GetIdVertex();
				unsigned int id_v_msh = this->GetElemID_FromCadID(id_v_cad,Cad::VERTEX);
				this->m_include_relation[id_atri].push_back(id_v_msh);
				assert( this->IsID(id_v_msh) );
				////////////////
				unsigned int id_e;
				bool is_same_dir;
				if( !pItr_l->GetIdEdge(id_e,is_same_dir) ) continue;
				unsigned int id_abar = this->GetElemID_FromCadID(id_e,Cad::EDGE);
				assert( this->IsID(id_abar) );
				this->m_include_relation[id_atri].push_back(id_abar);
			}
			if( !pItr_l->ShiftChildLoop() ) break;
		}
	}
	////////////////
  const std::vector<unsigned int>& id_e_ary = cad.GetAryElemID(Cad::EDGE);
	for(unsigned int iid_e=0;iid_e<id_e_ary.size();iid_e++){
		unsigned int id_e = id_e_ary[iid_e];
    assert( cad.IsElemID(Cad::EDGE,id_e) );
		unsigned int id_abar = this->GetElemID_FromCadID(id_e,Cad::EDGE);
    if( !this->IsID(id_abar) ){ continue; } // 浮いている辺があって，辺メッシュが切られなかった場合
		unsigned int id_v_s, id_v_e;
		if( !cad.GetIdVertex_Edge(id_v_s,id_v_e,id_e) ) assert(0);
		unsigned int id_v_s_msh = this->GetElemID_FromCadID(id_v_s,Cad::VERTEX);
		unsigned int id_v_e_msh = this->GetElemID_FromCadID(id_v_e,Cad::VERTEX);
		assert( this->IsID(id_v_s_msh) );
		assert( this->IsID(id_v_e_msh) );
		this->m_include_relation[id_abar].push_back(id_v_s_msh);
		this->m_include_relation[id_abar].push_back(id_v_e_msh);
	}
}

int CMesher2D::CheckMesh()
{
	{	// Check ID Overlap
		unsigned int max_id = 0;
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			const unsigned int id0 = m_aVertex[iver].id;
			if( max_id < id0 ) max_id = id0;
			assert( this->m_ElemType.size() > id0 );
			assert( this->m_ElemType[id0] == 0 );
			assert( this->m_ElemLoc.size() > id0 );
			const int iloc0 = this->m_ElemLoc[id0];
      assert( iloc0 == (int)iver );
			assert( this->m_aVertex.size() > (unsigned int)iloc0 );
			const SVertex& ver0 = this->m_aVertex[iloc0];
			assert( ver0.id == id0 );
		}
		for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
			const unsigned int id0 = m_aBarAry[ibarary].id;
			if( max_id < id0 ) max_id = id0;
			assert( this->m_ElemType.size() > id0 );
			assert( this->m_ElemType[id0] == 1 );
			assert( this->m_ElemLoc.size() > id0 );
			const int iloc0 = this->m_ElemLoc[id0];
      assert( iloc0 == (int)ibarary );
			assert( this->m_aBarAry.size() > (unsigned int)iloc0 );
			const CBarAry& ver0 = this->m_aBarAry[iloc0];
			assert( ver0.id == id0 );
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			const unsigned int id0 = m_aTriAry[itriary].id;
			if( max_id < id0 ) max_id = id0;
			assert( this->m_ElemType.size() > id0 );
			assert( this->m_ElemType[id0] == 2 );
			assert( this->m_ElemLoc.size() > id0 );
			const int iloc0 = this->m_ElemLoc[id0];
      assert( iloc0 == (int)itriary );
			assert( this->m_aTriAry.size() > (unsigned int)iloc0 );
			const CTriAry2D& ver0 = this->m_aTriAry[iloc0];
			assert( ver0.id == id0 );
		}
		for(unsigned int iquadary=0;iquadary<m_aQuadAry.size();iquadary++){
			const unsigned int id0 = m_aQuadAry[iquadary].id;
			if( max_id < id0 ) max_id = id0;
			assert( this->m_ElemType.size() > id0 );
			assert( this->m_ElemType[id0] == 3 );
			assert( this->m_ElemLoc.size() > id0 );
			const int iloc0 = this->m_ElemLoc[id0];
      assert( iloc0 == (int)iquadary );
			assert( this->m_aQuadAry.size() > (unsigned int)iloc0 );
			const CQuadAry2D& ver0 = this->m_aQuadAry[iloc0];
			assert( ver0.id == id0 );
		}
		assert( max_id == this->FindMaxID() );
		std::vector<unsigned int> is_used_flg_ary;
		is_used_flg_ary.resize( max_id+1, 0 );
		for(unsigned  int iver=0;iver<m_aVertex.size();iver++){
			assert( is_used_flg_ary[ m_aVertex[iver].id ] == 0 );
			is_used_flg_ary[ m_aVertex[iver].id ] = 1;
		}
		for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
			assert( is_used_flg_ary[ m_aBarAry[ibarary].id ] == 0 );
			is_used_flg_ary[ m_aBarAry[ibarary].id ] = 1;
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			assert( is_used_flg_ary[ m_aTriAry[itriary].id ] == 0 );
			is_used_flg_ary[ m_aTriAry[itriary].id ] = 1;
		}
		assert( is_used_flg_ary[0] == 0 );
	}
	for(unsigned int iver=0;iver<m_aVertex.size();iver++){
		assert( IsID( m_aVertex[iver].id ) );
	}
	for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
		assert( IsID( m_aBarAry[ibarary].id ) );
	}
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		assert( IsID( m_aTriAry[itriary].id ) );
	}
	{
		for(unsigned int ind=0;ind<m_ElemLoc.size();ind++){
			if( m_ElemLoc[ind] == -1 ) continue;
			assert( m_ElemLoc[ind] >= 0 );
			assert( IsID(ind) );
		}
	}
	////////////////////////////////	
	for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
		const unsigned int id_msh_bar = m_aBarAry[ibarary].id;
		const int iloc_bar = this->m_ElemLoc[id_msh_bar];
		const std::vector<SBar>& aBar = this->m_aBarAry[iloc_bar].m_aBar;
		for(unsigned int isidebar=0;isidebar<2;isidebar++){
			int id_msh_adj = this->m_aBarAry[iloc_bar].id_lr[isidebar];
			if( id_msh_adj <= 0 ) continue;	// 外部と接している場合
//			std::cout << id_msh_adj << " " << m_ElemLoc.size() << std::endl;
			assert(id_msh_adj<(int)m_ElemLoc.size());
			const int iloc_adj = this->m_ElemLoc[id_msh_adj];
			if( this->m_ElemType[id_msh_adj] == 2 ){	// 三角形と接している場合
				const std::vector<STri2D>& aTri = this->m_aTriAry[iloc_adj].m_aTri;
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					unsigned int itri = aBar[ibar].s2[isidebar];
					unsigned int inotri = aBar[ibar].r2[isidebar];
                    assert( aTri[itri].g2[inotri] == (int)id_msh_bar );
					assert( aTri[itri].s2[inotri] == ibar );
					assert( aTri[itri].r2[inotri] == isidebar );
				}
			}
		}
	}
	return 0;
}

void CMesher2D::SmoothingMesh_Delaunay(unsigned int& num_flip)
{
	num_flip = 0;
	std::vector<int> flg_vec(aVec2D.size(),0);
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;
		for(unsigned int itri=0;itri<aTri.size();itri++){	// 点周りの点を探索して調べる。
		for(unsigned int inotri=0;inotri<3;inotri++){
			const unsigned int ipoin = aTri[itri].v[inotri];
			assert( ipoin < aVec2D.size() );
			if( flg_vec[ipoin] != 0 ) continue;
			flg_vec[ipoin] = 1;
			DelaunayAroundPoint(itri,inotri,  aVec2D,aTri,  num_flip);
		}	// inotri
		}	// itri
	}	// itriary
	if( num_flip == 0 ) return;
	////////////////
	// 境界における要素との整合性をとる
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;		
		unsigned int id_this_loop = m_aTriAry[itriary].id;
		for(unsigned int itri=0;itri<aTri.size();itri++){
		for(unsigned int ifatri=0;ifatri<3;ifatri++){
			if( aTri[itri].g2[ifatri] < 0 ) continue;
			const unsigned int id0 = aTri[itri].g2[ifatri];
			const unsigned int iele0 = aTri[itri].s2[ifatri];
			assert( id0 < this->m_ElemType.size() );
			const unsigned int itype0 = m_ElemType[id0];
			assert( id0 < this->m_ElemLoc.size() );
			const int iloc0 = m_ElemLoc[id0];
			if( itype0 == 1 ){
				assert( (unsigned int)iloc0 < this->m_aBarAry.size() );
				CBarAry& bar_ary = m_aBarAry[iloc0];
				assert( bar_ary.id == id0 );
				assert( iele0 < bar_ary.m_aBar.size() );
				SBar& bar = bar_ary.m_aBar[iele0];
				const unsigned int iver0 = aTri[itri].v[ noelTriEdge[ifatri][0] ];
				const unsigned int iver1 = aTri[itri].v[ noelTriEdge[ifatri][1] ];
				if( iver0 == bar.v[0] && iver1 == bar.v[1] ){
                    assert( bar_ary.id_lr[0] == id_this_loop );
					bar.s2[0] = itri;
					bar.r2[0] = ifatri;
					aTri[itri].r2[ifatri] = 0;
				}
				else{
					assert( iver0 == bar.v[1] && iver1 == bar.v[0] );
                    assert( bar_ary.id_lr[1] == id_this_loop );
					bar.s2[1] = itri;
					bar.r2[1] = ifatri;
					aTri[itri].r2[ifatri] = 1;
				}
			}
			else{
				std::cout << "Error!-->not defined type" << itype0 << std::endl;
				assert(0);
			}
		}
		}
	}
}

void CMesher2D::SmoothingMesh_Laplace(unsigned int num_iter)
{
	std::vector<int> flg_vec(aVec2D.size(),0);
	for(unsigned int iiter=0;iiter<num_iter;iiter++){
		if( iiter != 0 ){
			for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){
				flg_vec[ivec] = 0;
			}
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			const std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;
			for(unsigned int itri=0;itri<aTri.size();itri++){	// 点周りの点を探索して調べる。
			for(unsigned int inotri=0;inotri<3;inotri++){
				const unsigned int ipoin = aTri[itri].v[inotri];
				assert( ipoin < aVec2D.size() );
				if( flg_vec[ipoin] != 0 ) continue;
				flg_vec[ipoin] = 1;
				const unsigned int itri_ini = itri;
				const unsigned int inoel_c_ini = inotri;

				unsigned int itri0= itri_ini;
				unsigned int inoel_c0 = inoel_c_ini;
				unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
				bool is_bound_flg = false;
				Com::CVector2D vec_delta = aVec2D[ipoin];
				unsigned int ntri_around = 1;
				for(;;){
					assert( itri0 < aTri.size() );
					assert( inoel_c0 < 3 );
					assert( aTri[itri0].v[inoel_c0] == ipoin );
					{
						vec_delta.x += aVec2D[ aTri[itri0].v[inoel_b0] ].x;
						vec_delta.y += aVec2D[ aTri[itri0].v[inoel_b0] ].y;
						ntri_around++;
					}
					if( aTri[itri0].g2[inoel_b0] == -2 ){
						unsigned int itri1 = aTri[itri0].s2[inoel_b0];
            const unsigned int rel01 = (int)aTri[itri0].r2[inoel_b0];
            unsigned int inoel_c1 = (unsigned int)relTriTri[rel01][inoel_c0];
            unsigned int inoel_b1 = (unsigned int)relTriTri[rel01][ (int)noelTriEdge[inoel_c0][1] ];
						assert( itri1 < aTri.size() );
            assert( aTri[itri1].s2[ (unsigned int)relTriTri[rel01][inoel_b0] ] == itri0 );
						assert( aTri[itri1].v[inoel_c1] == ipoin );
						if( itri1 == itri_ini ) break;
						itri0 = itri1;
						inoel_c0 = inoel_c1;
						inoel_b0 = inoel_b1;
					}
					else{	// この点は境界上の点だから動かしてはならない。
						is_bound_flg = true;
						break;
					}
				}
				if( is_bound_flg ) continue;		
				aVec2D[ipoin].x = vec_delta.x / ntri_around;
				aVec2D[ipoin].y = vec_delta.y / ntri_around;
			}	// inotri
			}	// itri
		}	// itriary
	}	// iiter
}


	// 負の体積や歪んだ要素がないか調べる
void CMesher2D::CheckMeshQuality(bool& is_inverted, double& max_aspect, const double ave_edge_len )
{
	is_inverted = false;
	max_aspect = 0.0;
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		const std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;
		for(unsigned int itri=0;itri<aTri.size();itri++){	// 点周りの点を探索して調べる。
			const double area = TriArea( aVec2D[ aTri[itri].v[0] ], aVec2D[ aTri[itri].v[1] ], aVec2D[ aTri[itri].v[2] ] );
			const double len0 = sqrt( SquareLength(aVec2D[aTri[itri].v[1]], aVec2D[aTri[itri].v[2]] ) );
			const double len1 = sqrt( SquareLength(aVec2D[aTri[itri].v[0]], aVec2D[aTri[itri].v[2]] ) );
			const double len2 = sqrt( SquareLength(aVec2D[aTri[itri].v[0]], aVec2D[aTri[itri].v[1]] ) );
			double max_len = len0;
			if( len1 > max_len ) max_len = len1;
			if( len2 > max_len ) max_len = len2;
			double aspect = max_len * (len0+len1+len2) * 0.5 / area;
			if( aspect > max_aspect ) max_aspect = aspect;
/*			if( max_len > max_edge_len_ratio * ave_edge_len ){
				max_edge_len_ratio = max_len / ave_edge_len;
			}*/
            if( area < ave_edge_len*ave_edge_len*1.0e-5 ){
//				std::cout << "Nega : " << itri << std::endl;
				is_inverted = true;
				break;
			}
		}
		if( is_inverted ) break;
	}
//	std::cout << max_edge_len_ratio << " " << max_aspect << std::endl;
}



// 平均メッシュ長の計算
double CMesher2D::GetAverageEdgeLength
(const Cad::ICad2D_Msh& cad_2d, 
 const std::set<unsigned int>& setIdL)
{
	double area = 0; // この点周りのループの面積
	unsigned int ntri = 0;	// この点周りのループの３角形分割数
  for(std::set<unsigned int>::const_iterator itrl=setIdL.begin();itrl!=setIdL.end();itrl++){
    unsigned int id_l_cad = (*itrl);
    if( !cad_2d.IsElemID(Cad::LOOP,id_l_cad) ) continue;
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_l_cad,Cad::LOOP) ) continue;	// このループにメッシュが切られてなかったら無視
		ntri += this->m_aTriAry[iloc].m_aTri.size();
		area += cad_2d.GetArea_Loop(id_l_cad);
	}
	const double ave_tri_area = area/ntri;
	return sqrt( ave_tri_area )*1.4;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

bool CMesher2D::Tessalate_Edge( const Cad::ICad2D_Msh& cad_2d, const unsigned int id_e )
{
	/*
	このEdgeに含まれるVertexが引数に存在するかどうか
	調べて存在しなかったら付け足すルーティンをそのうち追加
	*/

	unsigned int id_v_s, id_v_e;
    assert( cad_2d.IsElemID(Cad::EDGE,id_e) );
    if( !cad_2d.GetIdVertex_Edge(id_v_s,id_v_e,id_e) ){
        std::cout << "error edge : " << id_e << std::endl;
        assert(0);
    }

	unsigned int ip_s, ip_e;
	unsigned int id_msh_s, id_msh_e;
	{
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_v_s,Cad::VERTEX) ) assert(0);
		assert(itype==0 && iloc<m_aVertex.size());
		const Msh::SVertex& pvs = m_aVertex[iloc];
		ip_s = pvs.v;
		id_msh_s = pvs.id;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_v_e,Cad::VERTEX) ) assert(0);
		assert(itype==0 && iloc<m_aVertex.size());
		const Msh::SVertex& pve = m_aVertex[iloc];
		ip_e = pve.v;
		id_msh_e = pve.id;
	}

//	unsigned int curve_type = cad_2d.GetEdgeCurveType(id_e);
	unsigned int id_new_elem_ary = this->GetFreeObjID();
	const unsigned int ibarary0 = m_aBarAry.size();
	{
		this->m_ElemLoc.resize(id_new_elem_ary+1,-1);
		this->m_ElemType.resize(id_new_elem_ary+1);
		this->m_ElemLoc[id_new_elem_ary] = ibarary0;
		this->m_ElemType[id_new_elem_ary] = 1;
	}
	m_aBarAry.resize( m_aBarAry.size()+1 );
	CBarAry& ref_bar_ary = m_aBarAry[ibarary0];
    ////////////////
    std::vector<Com::CVector2D> aCo;
    cad_2d.GetCurveAsPolyline(id_e,aCo,-1);
    ////////////////
    const unsigned int ndiv = aCo.size()+1;
    std::vector<unsigned int> ipo_ary;
    {
      ipo_ary.resize(ndiv+1);
      ipo_ary[0] = ip_s;
      for(unsigned int i=1;i<ndiv;i++){
        ipo_ary[i] = aVec2D.size();
        aVec2D.push_back(aCo[i-1]);
      }
      ipo_ary[ndiv] = ip_e;
    }
    {
      ref_bar_ary.id = id_new_elem_ary;
      ref_bar_ary.id_e_cad = id_e;
      ref_bar_ary.ilayer = cad_2d.GetLayer(Cad::EDGE,id_e);
      ref_bar_ary.m_aBar.resize(ndiv);
      ref_bar_ary.id_se[0] = id_msh_s;
      ref_bar_ary.id_se[1] = id_msh_e;
      ref_bar_ary.id_lr[0] = 0;
      ref_bar_ary.id_lr[1] = 0;
      for(unsigned int ibar=0;ibar<ndiv;ibar++){
        ref_bar_ary.m_aBar[ibar].v[0] = ipo_ary[ibar];
        ref_bar_ary.m_aBar[ibar].v[1] = ipo_ary[ibar+1];
        ref_bar_ary.m_aBar[ibar].s2[0] = 0;
        ref_bar_ary.m_aBar[ibar].s2[1] = 0;
        ref_bar_ary.m_aBar[ibar].r2[0] = 0;
        ref_bar_ary.m_aBar[ibar].r2[1] = 0;
      }
    }
	assert( this->CheckMesh()==0 );

	return true;
}

bool CMesher2D::MakeMesh_Edge
(const Cad::ICad2D_Msh& cad_2d, unsigned int id_e , const double len)
{
	{
		assert( this->GetElemID_FromCadID(id_e,Cad::EDGE)==0 );
		/*
		このEdgeに含まれるVertexが引数に存在するかどうか
		調べて存在しなかったら付け足すルーティンをそのうち追加
		*/
	}

	unsigned int id_v_s, id_v_e;
	if( !cad_2d.GetIdVertex_Edge(id_v_s,id_v_e,id_e) ) assert(0);

	// 始点、終点のメッシュ点番号をip_s,ip_eに代入
	unsigned int ip_s, ip_e;
	unsigned int id_msh_s, id_msh_e;
	{		
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_v_s,Cad::VERTEX) ) assert(0);
		assert(itype==0 && iloc<m_aVertex.size());
		const Msh::SVertex& pvs = m_aVertex[iloc];
		ip_s = pvs.v;
		id_msh_s = pvs.id;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_v_e,Cad::VERTEX) ) assert(0);
		assert(itype==0 && iloc<m_aVertex.size());
		const Msh::SVertex& pve = m_aVertex[iloc];
		ip_e = pve.v;
		id_msh_e = pve.id;
	}

	unsigned int id_new_elem_ary = this->GetFreeObjID();
	const unsigned int ibarary0 = m_aBarAry.size();
	m_aBarAry.resize( m_aBarAry.size()+1 );
	{
		this->m_ElemLoc.resize(id_new_elem_ary+1,-1);
		this->m_ElemType.resize(id_new_elem_ary+1);
		this->m_ElemLoc[id_new_elem_ary] = ibarary0;
		this->m_ElemType[id_new_elem_ary] = 1;
	}
	CBarAry& ref_bar_ary = m_aBarAry[ibarary0];
  std::vector<Com::CVector2D> aCo;
  cad_2d.GetCurveAsPolyline(id_e,aCo,len);
  ////////////////
  const unsigned int ndiv = aCo.size()+1;
  std::vector<unsigned int> ipo_ary;
  {
    ipo_ary.resize(ndiv+1);
    ipo_ary[0] = ip_s;
    for(unsigned int i=1;i<ndiv;i++){
      ipo_ary[i] = aVec2D.size();
      aVec2D.push_back(aCo[i-1]);
    }
    ipo_ary[ndiv] = ip_e;
  }
  {
    ref_bar_ary.id = id_new_elem_ary;
    ref_bar_ary.id_e_cad = id_e;
    ref_bar_ary.ilayer = cad_2d.GetLayer(Cad::EDGE,id_e);
    ref_bar_ary.id_se[0] = id_msh_s;
    ref_bar_ary.id_se[1] = id_msh_e;
    ref_bar_ary.id_lr[0] = 0;
    ref_bar_ary.id_lr[1] = 0;
    ref_bar_ary.m_aBar.resize(ndiv);
    for(unsigned int ibar=0;ibar<ndiv;ibar++){
      ref_bar_ary.m_aBar[ibar].v[0] = ipo_ary[ibar];
      ref_bar_ary.m_aBar[ibar].v[1] = ipo_ary[ibar+1];
      ref_bar_ary.m_aBar[ibar].s2[0] = 0;
      ref_bar_ary.m_aBar[ibar].s2[1] = 0;
      ref_bar_ary.m_aBar[ibar].r2[0] = 0;
      ref_bar_ary.m_aBar[ibar].r2[1] = 0;
    }
  }  
	assert( this->CheckMesh()==0 );    
	return true;
}

bool CMesher2D::Tesselate_Loop
( const Cad::ICad2D_Msh& cad_2d, const unsigned int id_l )
{
	/*
	このLoop中のEdgeに含まれるVertexが引数に存在するかどうか
	調べて存在しなかったら付け足すルーティンをそのうち追加
	*/

	/*
	このLoop中のEdgが引数に存在するかどうか
	調べて存在しなかったら付け足すルーティンをそのうち追加
	*/

	std::vector<CPoint2D> aPo2D;
	std::vector<int> vec2po;
	{
		// 要素分割する領域の節点　aPo2Dを作成
		// 辺に属する節点の全体番号から、要素分割する領域のローカル番号への対応(vec2po)を作成

		////////////////////////////////
		vec2po.resize( aVec2D.size(), -1 );
		{	// このループで使用される節点のフラグを立てる
      std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
			for(;;){	// ループをめぐる
				for(;!pItrEdgeLoop->IsEnd();(*pItrEdgeLoop)++){	// このループの中のエッジをめぐる
					{
						const unsigned int id_v = pItrEdgeLoop->GetIdVertex();
						unsigned int iloc, itype;
						if( !this->FindElemLocType_CadIDType(iloc,itype,id_v,Cad::VERTEX) ) assert(0);
						assert( itype == 0 && iloc < m_aVertex.size() );
						const SVertex& ver = m_aVertex[iloc];
						vec2po[ ver.v ] = 1;
					}
					unsigned int id_e; bool is_same_dir;
					if( !pItrEdgeLoop->GetIdEdge(id_e,is_same_dir) ) continue; // 浮遊点は飛ばす
					unsigned int itype, iloc;
					if( !this->FindElemLocType_CadIDType(iloc,itype,id_e,Cad::EDGE) ) assert(0);
					assert( iloc < m_aBarAry.size() );
					assert( itype == 1 );
					const Msh::CBarAry& BarAry = m_aBarAry[iloc];
					assert( BarAry.id_e_cad == id_e );
					const std::vector<SBar>& aBar = BarAry.m_aBar;
					for(unsigned int ibar=0;ibar<aBar.size();ibar++){ 
						vec2po[ aBar[ibar].v[0] ] = 1;
						vec2po[ aBar[ibar].v[1] ] = 1;
					}
				}
				if( !pItrEdgeLoop->ShiftChildLoop() ) break;
			}
		}
		{	// vec2poを作る、aPo2Dを確保する
			int ipo=0;
			for(unsigned int ivec=0;ivec<vec2po.size();ivec++){
				if( vec2po[ivec] != -1 ){
					vec2po[ivec] = ipo;
					ipo++;
				}
			}
//			std::cout << "Loop : " << id_l << "  Npoint " << ipo << std::endl;
			aPo2D.resize(ipo);
		}
		for(unsigned int ivec=0;ivec<vec2po.size();ivec++){
			if( vec2po[ivec] != -1 ){
				const unsigned int ipoin = vec2po[ivec];
				assert( ipoin < aPo2D.size() );
				aPo2D[ipoin].p.x = aVec2D[ivec].x;
				aPo2D[ipoin].p.y = aVec2D[ivec].y;
				aPo2D[ipoin].e = -1;
        aPo2D[ipoin].d = 0;
			}
		}
	}

	std::vector<STri2D> aTri;
	{	// 与えられた点群を内部に持つ、大きな三角形を作る
		assert( aVec2D.size() >= 3 );
		double max_len;
		double center[2];
		{
			double bound_2d[4];
			bound_2d[0] = aPo2D[0].p.x;
			bound_2d[1] = aPo2D[0].p.x;
			bound_2d[2] = aPo2D[0].p.y;
			bound_2d[3] = aPo2D[0].p.y;
			for(unsigned int ipoin=1;ipoin<aPo2D.size();ipoin++){
				if( aPo2D[ipoin].p.x < bound_2d[0] ){ bound_2d[0] = aPo2D[ipoin].p.x; }
				if( aPo2D[ipoin].p.x > bound_2d[1] ){ bound_2d[1] = aPo2D[ipoin].p.x; }
				if( aPo2D[ipoin].p.y < bound_2d[2] ){ bound_2d[2] = aPo2D[ipoin].p.y; }
				if( aPo2D[ipoin].p.y > bound_2d[3] ){ bound_2d[3] = aPo2D[ipoin].p.y; }
			}
			max_len = (bound_2d[1]-bound_2d[0]>bound_2d[3]-bound_2d[2]) ? bound_2d[1]-bound_2d[0] : bound_2d[3]-bound_2d[2];
			center[0] = (bound_2d[1]+bound_2d[0])*0.5;
			center[1] = (bound_2d[3]+bound_2d[2])*0.5;
		}

//		std::cout << "center " << center[0] << " " << center[1] << std::endl;
//		std::cout << "max_len " << max_len << std::endl;

		const double tri_len = max_len * 4.0;
		const double tmp_len = tri_len * sqrt(3.0) / 6.0;

		const int npo = aPo2D.size();
		aPo2D.resize(npo+3);
		aPo2D[npo+0].p.x = center[0];				aPo2D[npo+0].p.y = center[1]+2.0*tmp_len;	aPo2D[npo+0].e = 0;	aPo2D[npo+0].d = 0;
		aPo2D[npo+1].p.x = center[0]-0.5*tri_len;	aPo2D[npo+1].p.y = center[1]-tmp_len;		aPo2D[npo+1].e = 0;	aPo2D[npo+1].d = 1;
		aPo2D[npo+2].p.x = center[0]+0.5*tri_len;	aPo2D[npo+2].p.y = center[1]-tmp_len;		aPo2D[npo+2].e = 0;	aPo2D[npo+2].d = 2;

		aTri.resize(1);
		aTri[0].v[0] = npo+0;	aTri[0].v[1] = npo+1;	aTri[0].v[2] = npo+2;
		aTri[0].g2[0] = -1; aTri[0].g2[1] = -1; aTri[0].g2[2] = -1;
		aTri[0].s2[0] =  0; aTri[0].s2[1] =  0; aTri[0].s2[2] =  0;
		aTri[0].r2[0] =  0; aTri[0].r2[1] =  0; aTri[0].r2[2] =  0;
//		assert( CheckTri(aPo2D,aTri) );
	}
//	OutInp("hoge1.inp",aPo2D,aTri);

	// Make Delaunay Division
	for(unsigned int ipoin=0;ipoin<aPo2D.size();ipoin++){
		if( aPo2D[ipoin].e >= 0 ) continue;	// 既にメッシュの一部である。
		const CVector2D& po_add = aPo2D[ipoin].p;
		int itri_in = -1;
		int iedge = -1;
		unsigned int iflg1 = 0, iflg2 = 0;
		for(unsigned int itri=0;itri<aTri.size();itri++){
			iflg1 = 0; iflg2 = 0;
			const STri2D& ref_tri = aTri[itri];
			if( TriArea(po_add, aPo2D[ref_tri.v[1]].p, aPo2D[ref_tri.v[2]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 0;
			}
			if( TriArea(po_add, aPo2D[ref_tri.v[2]].p, aPo2D[ref_tri.v[0]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 1;
			}
			if( TriArea(po_add, aPo2D[ref_tri.v[0]].p, aPo2D[ref_tri.v[1]].p ) > MIN_TRI_AREA ){
				iflg1++; iflg2 += 2;
			}
			if( iflg1 == 3 ){
				itri_in = itri;
				break;
			}
			else if( iflg1 == 2 ){
				const unsigned int ied0 = 3-iflg2;
				const unsigned int ipo_e0 = ref_tri.v[ noelTriEdge[ied0][0] ];
				const unsigned int ipo_e1 = ref_tri.v[ noelTriEdge[ied0][1] ];
        const unsigned int* rel = relTriTri[ ref_tri.r2[ied0] ];
				const unsigned int itri_s = ref_tri.s2[ied0];
        assert( aTri[itri_s].v[ rel[ noelTriEdge[ied0][0] ] ] == ipo_e0 );
        assert( aTri[itri_s].v[ rel[ noelTriEdge[ied0][1] ] ] == ipo_e1 );
				const unsigned int inoel_d = rel[ied0];
				assert( aTri[itri_s].s2[inoel_d] == itri );
				const unsigned int ipo_d = aTri[itri_s].v[inoel_d];
				assert( TriArea( po_add, aPo2D[ipo_e1].p, aPo2D[ aTri[itri].v[ied0] ].p ) > MIN_TRI_AREA );
				assert( TriArea( po_add, aPo2D[ aTri[itri].v[ied0] ].p, aPo2D[ipo_e0].p ) > MIN_TRI_AREA );
				if( TriArea( po_add, aPo2D[ipo_e0].p, aPo2D[ipo_d ].p ) < MIN_TRI_AREA ){ continue;	}
				if( TriArea( po_add, aPo2D[ipo_d ].p, aPo2D[ipo_e1].p ) < MIN_TRI_AREA ){ continue; }
				const unsigned int det_d =  DetDelaunay(po_add,aPo2D[ipo_e0].p,aPo2D[ipo_e1].p,aPo2D[ipo_d].p);
				if( det_d == 2 || det_d == 1 ) continue;
				itri_in = itri;
				iedge = ied0;
				break;
			}
		}
		if( itri_in == -1 ){
      std::cout << "Super Triangle Failure " << ipoin << " " << po_add.x << " " << po_add.y << std::endl;
			std::cout << aTri.size() << std::endl;
      return false;
			assert(0);
		}
		if( iedge == -1 ){
			InsertPoint_Elem(ipoin,itri_in,aPo2D,aTri);
		}
		else{
			InsertPoint_ElemEdge(ipoin,itri_in,iedge,aPo2D,aTri);
		}
		DelaunayAroundPoint(ipoin,aPo2D,aTri);
	}
	assert( CheckTri(aPo2D,aTri) );
//	OutInp("hoge2.inp",aPo2D,aTri);

	const unsigned int id_new_tri_ary = this->GetFreeObjID();

	{	// エッジを回復する
    std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
		for(;;){	// 子ループのためのループ
    for(;!pItrEdgeLoop->IsEnd();(*pItrEdgeLoop)++){
      unsigned int id_e;
			bool is_same_dir;
			if( !pItrEdgeLoop->GetIdEdge(id_e,is_same_dir) ) continue;	// ループの中の点
			unsigned int iloc, itype;
			if( !FindElemLocType_CadIDType(iloc,itype, id_e,Cad::EDGE) ) assert(0);
			assert( iloc < m_aBarAry.size() );
			assert( itype == 1 );
			const Msh::CBarAry& BarAry = m_aBarAry[iloc];
			assert( id_e == BarAry.id_e_cad );
			const std::vector<SBar>& aBar = BarAry.m_aBar;
			const unsigned int id_elem_bar = BarAry.id;
			assert( id_elem_bar != id_new_tri_ary );
			for(unsigned int ibar=0;ibar<aBar.size();ibar++){
				for(;;){ // EdgeをFlipしたら同じ辺について繰り返す				
				unsigned int ipoi0,ipoi1; // ipoi0は左周りのbarの始点、ipoi1は終点
				if( is_same_dir ){	ipoi0 = vec2po[ aBar[ibar].v[0] ]; ipoi1 = vec2po[ aBar[ibar].v[1] ]; }
				else{				ipoi0 = vec2po[ aBar[ibar].v[1] ]; ipoi1 = vec2po[ aBar[ibar].v[0] ]; }
				assert( ipoi0 < aPo2D.size() ); assert( ipoi1 < aPo2D.size() );
				////////////////
				unsigned int itri0;
				unsigned int inotri0,inotri1;
				if( FindEdge(ipoi0,ipoi1,itri0,inotri0,inotri1,aPo2D,aTri) ){	// ループの内側に接する要素を見つける
					// Split Triangle
					assert( inotri0 != inotri1 );
					assert( inotri0 < 3 );
					assert( inotri1 < 3 );
					assert( aTri[itri0].v[ inotri0 ] == ipoi0 );
					assert( aTri[itri0].v[ inotri1 ] == ipoi1 );
					const unsigned int ied0 = 3 - inotri0 - inotri1;
					{
						const unsigned int itri1 = aTri[itri0].s2[ied0];
            const unsigned int ied1 = (unsigned int)relTriTri[ (int)aTri[itri0].r2[ied0] ][ied0];
						assert( aTri[itri1].s2[ied1] == itri0 );
						aTri[itri1].g2[ied1] = -3;	
						aTri[itri0].g2[ied0] = -3;
					}
					break;	// 次のBarへ　for(;;)を抜ける
				}
				else{
					double ratio;
					if( !FindEdgePoint_AcrossEdge(ipoi0,ipoi1,   itri0,inotri0,inotri1,ratio,   aPo2D,aTri) ){ 
            std::cout << "歪んだメッシュ" << std::endl;
						return false;
						assert(0); 
					}
					// return false if degeneration
					if( ratio < -1.0e-20 || ratio > 1.0+1.0e-20 ){ return false; }
					if( TriArea( aPo2D[ipoi0].p, aPo2D[ aTri[itri0].v[inotri0] ].p, aPo2D[ipoi1].p ) < 1.0e-20 ){ return false; }
					if( TriArea( aPo2D[ipoi0].p, aPo2D[ipoi1].p, aPo2D[ aTri[itri0].v[inotri1] ].p ) < 1.0e-20 ){ return false; }
/*					assert( ratio > -1.0e-20 && ratio < 1.0+1.0e-20 );
					assert( TriArea( aPo2D[ipoi0].p, aPo2D[ aTri[itri0].v[inotri0] ].p, aPo2D[ipoi1].p ) > 1.0e-20 );
					assert( TriArea( aPo2D[ipoi0].p, aPo2D[ipoi1].p, aPo2D[ aTri[itri0].v[inotri1] ].p ) > 1.0e-20 );*/
					if( ratio < 1.0e-20 ){
//						std::cout << "未実装 辺上に点がある場合" << std::endl;
						return false;
						assert(0);
					}
					else if( ratio > 1.0-1.0e-10 ){
//						std::cout << "未実装 辺上に点がある場合" << std::endl;
						return false;
						assert(0);
					}
					else{
						const unsigned int ied0 = 3 - inotri0 - inotri1;
						if( aTri[itri0].g2[ied0] != -2 ) return false;
						assert( aTri[itri0].g2[ied0] == -2 );
						const unsigned int itri1 = aTri[itri0].s2[ied0];
                        const unsigned int ied1 = (unsigned int)relTriTri[ (int)aTri[itri0].r2[ied0] ][ied0];
						assert( aTri[itri1].s2[ied1] == itri0 );
						assert( aTri[itri1].g2[ied1] == -2 );
						FlipEdge(itri0,ied0,aPo2D,aTri);
						continue;
					}
				}
			}
			}
		}
		if( !pItrEdgeLoop->ShiftChildLoop() ) break;
		}
		assert( CheckTri(aPo2D,aTri) );
	}

	////////////////////////////////////////////////
	// ここからはクラスの内容を変更する
	// エラーを出して戻るなら、ここ以前にすること
	////////////////////////////////////////////////

	// ここから辺要素の隣接関係を変更する．３角形についてはそのまま

	{	// 辺要素から３角形要素への隣接情報を作成
    std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
		for(;;){	// 子ループのためのループ
		for(;!pItrEdgeLoop->IsEnd();(*pItrEdgeLoop)++){
			unsigned int id_e;
			bool is_same_dir;
			if( !pItrEdgeLoop->GetIdEdge(id_e,is_same_dir) ) continue;	// ループの中の点
			unsigned int iloc, itype;
			if( !FindElemLocType_CadIDType(iloc,itype, id_e,Cad::EDGE) ) assert(0);
			assert( iloc < m_aBarAry.size() );
			assert( itype == 1 );
			Msh::CBarAry& BarAry = m_aBarAry[iloc];
			assert( id_e == BarAry.id_e_cad );
			std::vector<SBar>& aBar = BarAry.m_aBar;
			const unsigned int id_elem_bar = BarAry.id;
			assert( id_elem_bar != id_new_tri_ary );
			for(unsigned int ibar=0;ibar<aBar.size();ibar++){					
				unsigned int ipoi0,ipoi1; // ipoi0は左周りのbarの始点、ipoi1は終点
				if( is_same_dir ){	ipoi0 = vec2po[ aBar[ibar].v[0] ]; ipoi1 = vec2po[ aBar[ibar].v[1] ]; }
				else{				ipoi0 = vec2po[ aBar[ibar].v[1] ]; ipoi1 = vec2po[ aBar[ibar].v[0] ]; }
				assert( ipoi0 < aPo2D.size() ); assert( ipoi1 < aPo2D.size() );
				////////////////
				unsigned int itri0;
				unsigned int inotri0,inotri1;
				if( !FindEdge(ipoi0,ipoi1,itri0,inotri0,inotri1,aPo2D,aTri) ){	assert(0); }// ループの内側に接する要素を見つける
				assert( inotri0 != inotri1 );
				assert( inotri0 < 3 );
				assert( inotri1 < 3 );
				assert( aTri[itri0].v[ inotri0 ] == ipoi0 );
				assert( aTri[itri0].v[ inotri1 ] == ipoi1 );
				const unsigned int ied0 = 3 - inotri0 - inotri1;
				// 辺要素の隣接情報を作る
				if( is_same_dir ){	
          assert( BarAry.id_lr[0] == id_new_tri_ary || BarAry.id_lr[0] == 0 );
					BarAry.id_lr[0] = id_new_tri_ary;
					aBar[ibar].s2[0]=itri0; aBar[ibar].r2[0]=ied0; 
				}
				else{				
          assert( BarAry.id_lr[1] == id_new_tri_ary || BarAry.id_lr[1] == 0 );
					BarAry.id_lr[1] = id_new_tri_ary;
					aBar[ibar].s2[1]=itri0; aBar[ibar].r2[1]=ied0; 
				}
			}
		}
		if( !pItrEdgeLoop->ShiftChildLoop() ) break;
		}
		assert( CheckTri(aPo2D,aTri) );
	}
	
	// 今後は辺要素を変更するのは，TriAryの番号付けを変化させるとき

//	OutInp("hoge2.inp",aPo2D,aTri);

	{	// 辺との隣接番号の整合性をとる
    std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
		for(;;){	// 子ループのためのループ
		for(;!pItrEdgeLoop->IsEnd();(*pItrEdgeLoop)++){
			unsigned int id_e;
			bool is_same_dir;
			if( !pItrEdgeLoop->GetIdEdge(id_e,is_same_dir) ) continue;	// 子ループが点の場合
			unsigned int iloc, itype;
			if( !FindElemLocType_CadIDType(iloc,itype, id_e,Cad::EDGE) ) assert(0);
			assert( iloc < m_aBarAry.size() );
			assert( itype == 1 );
			const Msh::CBarAry& BarAry = m_aBarAry[iloc];
			assert( id_e == BarAry.id_e_cad );
			const std::vector<SBar>& aBar = BarAry.m_aBar;
			const unsigned int id_elem_bar = BarAry.id;
			assert( id_elem_bar != id_new_tri_ary );
      if( BarAry.id_lr[0] == id_new_tri_ary ){	// 左側を切り離す
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					const SBar& bar = aBar[ibar];
					const unsigned int itri0 = bar.s2[0];
					const unsigned int ied0  = bar.r2[0];
					assert( itri0 < aTri.size() );
					assert( ied0 < 3 );
          assert( (int)aTri[itri0].v[noelTriEdge[ied0][0]]==vec2po[bar.v[0]] || (int)aTri[itri0].v[noelTriEdge[ied0][0]]==vec2po[bar.v[1]]);
          assert( (int)aTri[itri0].v[noelTriEdge[ied0][1]]==vec2po[bar.v[0]] || (int)aTri[itri0].v[noelTriEdge[ied0][1]]==vec2po[bar.v[1]]);
          if( aTri[itri0].g2[ied0] == (int)id_elem_bar ) continue;	// すでに切り離されてる
					{	// 向かい側の要素の処理
						const unsigned int itri1 = aTri[itri0].s2[ied0];
            const unsigned int ied1 = (unsigned int)relTriTri[ (int)aTri[itri0].r2[ied0] ][ied0];
						assert( aTri[itri1].s2[ied1] == itri0 );
            if( BarAry.id_lr[1] != id_new_tri_ary ){	// 外側の要素を切り離す
							assert( aTri[itri1].s2[ied1] == itri0 );
							aTri[itri1].g2[ied1] = -1;
						}
						else{ // 辺をはさんで向かい側の要素も内側だから辺にくっつける
							aTri[itri1].g2[ied1] = id_elem_bar;
							aTri[itri1].s2[ied1] = ibar;
							aTri[itri1].r2[ied1] = 1;
						}
					}
					{	// 内側の要素を辺にくっつける
						aTri[itri0].g2[ied0] = id_elem_bar;
						aTri[itri0].s2[ied0] = ibar;
						aTri[itri0].r2[ied0] = 0;
					}
				}
			}			
      if( BarAry.id_lr[1] == id_new_tri_ary ){	// 辺の右側を切り離す
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					const SBar& bar = aBar[ibar];
					const unsigned int itri0 = bar.s2[1];
					const unsigned int ied0= bar.r2[1];
          if( aTri[itri0].g2[ied0] == (int)id_elem_bar ) continue;	// すでに切り離されてる
					{	// 外側の要素を切り離す
						const unsigned int itri1 = aTri[itri0].s2[ied0];
            const unsigned int ied1 = relTriTri[ (int)aTri[itri0].r2[ied0] ][ied0];
						assert( itri1 < aTri.size() );
						assert( ied1 < 3 );
            assert( (int)aTri[itri0].v[noelTriEdge[ied0][0]]==vec2po[bar.v[1]] );
            assert( (int)aTri[itri0].v[noelTriEdge[ied0][1]]==vec2po[bar.v[0]] );
            if( BarAry.id_lr[0] != id_new_tri_ary ){	// 外側の要素を切り離す
							assert( aTri[itri1].s2[ied1] == itri0 );
							aTri[itri1].g2[ied1] = -1;
						}
						else{ // 辺をはさんで向かい側の要素も内側だから辺にくっつける
							aTri[itri1].g2[ied1] = id_elem_bar;
							aTri[itri1].s2[ied1] = ibar;
							aTri[itri1].r2[ied1] = 0;
						}
					}
					{	// 内側の要素を辺にくっつける
						aTri[itri0].g2[ied0] = id_elem_bar;
						aTri[itri0].s2[ied0] = ibar;
						aTri[itri0].r2[ied0] = 1;
					}
				}
			}
		}
		if( !pItrEdgeLoop->ShiftChildLoop() ) break;
		}	// ループのfor文終わり
		// ここから先はFlip禁止フラグ(隣接要素配列番号-3)はないはず
		assert( CheckTri(aPo2D,aTri) );	
	}

	// 外側の３角形の消去
	////////////////////////////////////////////////

	std::vector<STri2D> aTri_in;	// 内側の三角形
	{	// 外側の三角形の除去
		// 内側にある三角形をひとつ(itri0_ker)見つける
		unsigned int itri0_ker = aTri.size();
		{
      std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
			Cad::IItrLoop& itrEdgeLoop = *pItrEdgeLoop;
			for(;!itrEdgeLoop.IsEnd();itrEdgeLoop++){
				unsigned int id_e; bool is_same_dir;
				itrEdgeLoop.GetIdEdge(id_e,is_same_dir);
				unsigned int iloc, itype;
				if( !FindElemLocType_CadIDType(iloc,itype, id_e,Cad::EDGE) ) assert(0);
				assert( itype == 1 );
				assert( iloc < m_aBarAry.size() );
				const Msh::CBarAry& BarAry = m_aBarAry[iloc];
				assert( id_e == BarAry.id_e_cad );
				const std::vector<SBar>& aBar = BarAry.m_aBar;
        if(      BarAry.id_lr[0] == id_new_tri_ary ){
					for(unsigned int ibar=0;ibar<aBar.size();ibar++){
						itri0_ker = aBar[ibar].s2[0];
						break;
					}
				}
        else if( BarAry.id_lr[1] == id_new_tri_ary ){
					for(unsigned int ibar=0;ibar<aBar.size();ibar++){
						itri0_ker = aBar[ibar].s2[1];
						break;
					}
				}
			}
		}
		assert( itri0_ker < aTri.size() );

		// 領域の外の要素ならフラグが-1、そうでなければフラグは昇順の要素番号が入った配列inout_flgを作る
		unsigned int ntri_in;
		std::vector<int> inout_flg;	// フラグ配列
		{	// 上で見つけた内側の三角形を核として内側の三角形を周囲に拡大していく
			inout_flg.resize(aTri.size(),-1);
			inout_flg[itri0_ker] = 0;
			ntri_in = 1;
			std::stack<unsigned int> ind_stack;	// 周囲が探索されていない三角形
			ind_stack.push(itri0_ker);
			for(;;){
				if( ind_stack.empty() ) break;
				const unsigned int itri_cur = ind_stack.top();
				ind_stack.pop();
				for(unsigned int inotri=0;inotri<3;inotri++){
					if( aTri[itri_cur].g2[inotri] != -2 ) continue;
					const unsigned int itri_s = aTri[itri_cur].s2[inotri];
					if( inout_flg[itri_s] == -1 ){
						inout_flg[itri_s] = ntri_in;
						ntri_in++;
						ind_stack.push(itri_s);
					}
				}
			}
		}

		// フラグ配列に沿って内側の三角形を集めた配列aTri_inを作る
		aTri_in.resize( ntri_in );
		for(unsigned int itri=0;itri<aTri.size();itri++){
			if( inout_flg[itri] != -1 ){
				int itri_in = inout_flg[itri];
				assert( itri_in >= 0 && (unsigned int)itri_in < ntri_in );
				aTri_in[itri_in] = aTri[itri];
			}
		}
		// 内側の三角形配列のの隣接情報を作る
		for(unsigned int itri=0;itri<aTri_in.size();itri++){
			for(unsigned int ifatri=0;ifatri<3;ifatri++){
				if( aTri_in[itri].g2[ifatri] != -2 ) continue;
				int itri_s0 = aTri_in[itri].s2[ifatri];
				assert( itri_s0 >= 0 && (unsigned int)itri_s0 < aTri.size() );
				int itri_in_s0 = inout_flg[itri_s0];
				assert( itri_in_s0 >= 0 && (unsigned int)itri_in_s0 < aTri_in.size() );
				aTri_in[itri].s2[ifatri] = itri_in_s0;
			}
		}
		{	// 辺の隣接情報を更新
			std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
			for(;;){	// 子ループのためのループ
			for(;!pItrEdgeLoop->IsEnd();(*pItrEdgeLoop)++){
				unsigned int id_e;
				bool is_same_dir;
				if( !pItrEdgeLoop->GetIdEdge(id_e,is_same_dir) ) continue;	// ループの中の点
				unsigned int iloc, itype;
				if( !FindElemLocType_CadIDType(iloc,itype, id_e,Cad::EDGE) ) assert(0);
				assert( iloc < m_aBarAry.size() );
				assert( itype == 1 );
				Msh::CBarAry& BarAry = m_aBarAry[iloc];
				assert( id_e == BarAry.id_e_cad );
				std::vector<SBar>& aBar = BarAry.m_aBar;
				const unsigned int id_elem_bar = BarAry.id;
				assert( id_elem_bar != id_new_tri_ary );
				unsigned int iside = ( is_same_dir ) ? 0 : 1;
				assert( BarAry.id_lr[iside] == id_new_tri_ary );
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					SBar& bar = aBar[ibar];
					int itri_s0 = bar.s2[iside];
					assert( itri_s0 >= 0 && (unsigned int)itri_s0 < aTri.size() );
					int itri_in_s0 = inout_flg[itri_s0];
					assert( itri_in_s0 >= 0 && (unsigned int)itri_in_s0 < aTri_in.size() );
					bar.s2[iside] = itri_in_s0;
				}
			}
			if( !pItrEdgeLoop->ShiftChildLoop() ) break;
			}
		}
		inout_flg.clear();
		for(unsigned int ipo=0;ipo<aPo2D.size();ipo++){ aPo2D[ipo].e = -1; }
		assert( CheckTri(aPo2D,aTri_in) );
	}
	{	// Remove not used point
		std::vector<int> po2vec;
		po2vec.resize(aPo2D.size(),-2);
		for(unsigned int itri=0;itri<aTri_in.size();itri++){
			po2vec[ aTri_in[itri].v[0] ] = -1;
			po2vec[ aTri_in[itri].v[1] ] = -1;
			po2vec[ aTri_in[itri].v[2] ] = -1;
		}
		for(unsigned int ivec=0;ivec<vec2po.size();ivec++){
			if( vec2po[ivec] != -1 ){
				const unsigned int ipo0 = vec2po[ivec];
        if( po2vec[ipo0] != -1 ){ std::cout << "対応しない点" << std::endl; return false; }
				assert( po2vec[ipo0] == -1 );
				po2vec[ipo0] = ivec;
			}
		}
		for(unsigned  int ipo=0;ipo<po2vec.size();ipo++){
			if( po2vec[ipo] == -1 ){
				std::cout << ipo << " " << aPo2D[ipo].p.x << " " << aPo2D[ipo].p.y << std::endl;
//				std::cout << "未実装  Ｌｏｏｐに新しい節点の追加したときの処理" << std::endl;
				return false;
				assert(0);
			}
		}
		for(unsigned int itri=0;itri<aTri_in.size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				const int ipo0 = aTri_in[itri].v[inotri];
				assert( ipo0 >= 0 && (unsigned int)ipo0 < aPo2D.size() );
				const int ivec0 = po2vec[ipo0];
				assert( ivec0 >=0 && (unsigned int)ivec0 < aVec2D.size() );
				aTri_in[itri].v[inotri] = ivec0;
			}
		}
	}
//	OutInp("hoge3.inp",aVec2D, aTri_in);

	{
		unsigned int itriary = m_aTriAry.size();
		m_aTriAry.resize( m_aTriAry.size()+1 );
		m_aTriAry[ itriary ].m_aTri = aTri_in;
		m_aTriAry[ itriary ].id_l_cad = id_l;
		m_aTriAry[ itriary ].id = id_new_tri_ary;
		m_aTriAry[ itriary ].ilayer = cad_2d.GetLayer(Cad::LOOP,id_l);
		this->m_ElemType.resize(id_new_tri_ary+1,-1);
		this->m_ElemLoc.resize(id_new_tri_ary+1);
		this->m_ElemType[id_new_tri_ary] = 2;	// TRI
		this->m_ElemLoc[id_new_tri_ary] = itriary;
	}

	assert( this->CheckMesh() == 0 );
	return true;
}

bool CMesher2D::MakeMesh_Loop
( const Cad::ICad2D_Msh& cad_2d, unsigned int id_cad_l, const double len )
{
	if( !Tesselate_Loop(cad_2d,id_cad_l) ){
		std::cout << "Tesselation_Loop Fail" << std::endl;
		assert(0);
		return false;
	}

	// されたLoopを編集するため
	// vec2poを作成
	// aPo2Dを作成
	// aTriをコピーして作成
	std::vector<CPoint2D> aPo2D;
	std::vector<STri2D> aTri;
	std::vector<int> vec2po; // MSH節点番号vecからローカル節点番号poへのフラグ、対応してない場合は-2が入る
	{
		vec2po.resize( aVec2D.size(), -2 );
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_cad_l,Cad::LOOP) ) assert(0);
		assert( itype == 2 );
		assert( iloc < m_aTriAry.size() );
		const CTriAry2D& TriAry = m_aTriAry[iloc];
		const std::vector<STri2D>& aTri_ini = TriAry.m_aTri;
		for(unsigned int itri=0;itri<aTri_ini.size();itri++){	// ３角形に使われている全ての節点をマーク
			for(unsigned int inotri=0;inotri<3;inotri++){
				const unsigned int ivec0 = aTri_ini[itri].v[inotri];
				assert( ivec0 < aVec2D.size() );
				vec2po[ivec0] = -1;
			}
		}
		unsigned int npo = 0;
		for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){
			if( vec2po[ivec] == -1 ){
				vec2po[ivec] = npo;
				npo++;
			}
			else{ assert( vec2po[ivec] == -2 );	}
		}
		aPo2D.resize(npo);
		for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){
			if( vec2po[ivec] >= 0 ){
				const int ipo0 = vec2po[ivec];
				assert( ipo0 >=0 && ipo0 < (int)aPo2D.size() );
				aPo2D[ipo0].p.x = aVec2D[ivec].x;
				aPo2D[ipo0].p.y = aVec2D[ivec].y;
			}
			else{
				assert( vec2po[ivec] == -2 );
			}
		}
		aTri = aTri_ini;
		for(unsigned int itri=0;itri<aTri_ini.size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				const unsigned int ivec0 = aTri_ini[itri].v[inotri];
				assert( ivec0 < aVec2D.size() );
				const int ipo0 = vec2po[ivec0];
				assert( ipo0 >=0 && ipo0< (int)aPo2D.size() );
				assert( aTri[itri].v[inotri] == aTri_ini[itri].v[inotri] );
				aTri[itri].v[inotri] = ipo0;
			}
		}
		for(unsigned int itri=0;itri<aTri.size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				const unsigned int ipo0 = aTri[itri].v[inotri];
                assert( ipo0< aPo2D.size() );
				aPo2D[ipo0].e = itri;
				aPo2D[ipo0].d = inotri;
			}
		}
		assert( CheckTri(aPo2D,aTri) );
	}

	std::vector<unsigned int> aflag_isnt_move;	// フラグが１なら動かさない
	{
		aflag_isnt_move.resize( aPo2D.size(), 0 );
		for(unsigned int iver=0;iver<m_aVertex.size();iver++){
			const unsigned int ivec = this->m_aVertex[iver].v;
			if( ivec < vec2po.size() ){
				if( vec2po[ivec] == -2 ) continue;
				assert( vec2po[ivec] >= 0 );
				const unsigned int ipo = vec2po[ivec];
				if( ipo < aPo2D.size() ){
					aflag_isnt_move[ ipo ] = 1;
				}
			}
		}
	}

	{	// aTriに節点を追加
		double ratio = 3.0;
		for(;;){
			unsigned int nadd = 0;
			for(unsigned int itri=0;itri<aTri.size();itri++){
				const double area = TriArea(
					aPo2D[aTri[itri].v[0]].p, 
					aPo2D[aTri[itri].v[1]].p, 
					aPo2D[aTri[itri].v[2]].p);
				if( area > len * len * ratio ){
					// itriの重心に新しい節点を追加
					const unsigned int ipo0 = aPo2D.size();	// ipo0は新しい節点番号
					aPo2D.resize( aPo2D.size()+1 );
					aPo2D[ipo0].p.x = (aPo2D[aTri[itri].v[0]].p.x+aPo2D[aTri[itri].v[1]].p.x+aPo2D[aTri[itri].v[2]].p.x)/3.0;
					aPo2D[ipo0].p.y = (aPo2D[aTri[itri].v[0]].p.y+aPo2D[aTri[itri].v[1]].p.y+aPo2D[aTri[itri].v[2]].p.y)/3.0;
					InsertPoint_Elem(ipo0,itri,aPo2D,aTri);
					DelaunayAroundPoint(ipo0,aPo2D,aTri);
					nadd++;
				}
			}
			LaplacianSmoothing(aPo2D,aTri,aflag_isnt_move);
//			LaplaceDelaunaySmoothing(aPo2D,aTri);
			if( nadd != 0 ){ ratio *= 0.8; }
			else{ ratio *= 0.5; }
			if( ratio < 0.65 ) break;
		}
	}

	LaplaceDelaunaySmoothing(aPo2D,aTri,aflag_isnt_move);

	for(unsigned int i=0;i<20;i++){
//		PliantBossenHeckbertSmoothing(len,aPo2D,aTri);
//		LaplaceDelaunaySmoothing(aPo2D,aTri);
//		LaplacianSmoothing(aPo2D,aTri);
	}

	// 全体節点番号へ直す
	std::vector<int> po2vec;
	po2vec.resize(aPo2D.size(),-2);
	for(unsigned int itri=0;itri<aTri.size();itri++){
		po2vec[ aTri[itri].v[0] ] = -1;
		po2vec[ aTri[itri].v[1] ] = -1;
		po2vec[ aTri[itri].v[2] ] = -1;
	}
	for(unsigned int ivec=0;ivec<vec2po.size();ivec++){
		if( vec2po[ivec] >= 0 ){
			const unsigned int ipo0 = vec2po[ivec];
			assert( po2vec[ipo0] == -1 );
			po2vec[ipo0] = ivec;
		}
		else{
			assert( vec2po[ivec] == -2 );
		}
	}
	{	// 全体節点を追加
		unsigned int npo_add = 0;
		for(unsigned  int ipo=0;ipo<po2vec.size();ipo++){
			if( po2vec[ipo] == -1 ){
				npo_add++;
			}
		}
		aVec2D.reserve( aVec2D.size()+npo_add );
		for(unsigned  int ipo=0;ipo<po2vec.size();ipo++){
			if( po2vec[ipo] == -1 ){
				CVector2D vec0;
				vec0.x = aPo2D[ipo].p.x;
				vec0.y = aPo2D[ipo].p.y;
				const unsigned int ivec0 = aVec2D.size();
				aVec2D.push_back(vec0);
				po2vec[ipo] = ivec0;
			}
		}
	}
	{	// ローカル節点番号から全体節点番号への並び替え
		for(unsigned int itri=0;itri<aTri.size();itri++){
			for(unsigned int inotri=0;inotri<3;inotri++){
				const int ipo0 = aTri[itri].v[inotri];
				assert( ipo0 >=0 && (unsigned int)ipo0 < aPo2D.size() );
				const unsigned int ivec0 = po2vec[ipo0];
				assert( ivec0 < aVec2D.size() );
				aTri[itri].v[inotri] = ivec0;
			}
		}
	}

	unsigned int id_this_loop;
	{
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype,id_cad_l,Cad::LOOP) ) assert(0);
		assert( itype == 2 );
		assert( iloc < m_aTriAry.size() );
		const CTriAry2D& TriAry = m_aTriAry[iloc];
		id_this_loop = TriAry.id;
	}

	{	// 境界における要素との整合性をとる
		for(unsigned int itri=0;itri<aTri.size();itri++){
			for(unsigned int ifatri=0;ifatri<3;ifatri++){
				if( aTri[itri].g2[ifatri] < 0 ) continue;
				const unsigned int id0 = aTri[itri].g2[ifatri];
				const unsigned int iele0 = aTri[itri].s2[ifatri];
				assert( id0 < this->m_ElemType.size() );
				const unsigned int itype0 = m_ElemType[id0];
				assert( id0 < this->m_ElemLoc.size() );
				const int iloc0 = m_ElemLoc[id0];
				if( itype0 == 1 ){
					assert( (unsigned int)iloc0 < this->m_aBarAry.size() );
					CBarAry& bar_ary = m_aBarAry[iloc0];
					assert( bar_ary.id == id0 );
					assert( iele0 < bar_ary.m_aBar.size() );
					SBar& bar = bar_ary.m_aBar[iele0];
					const unsigned int iver0 = aTri[itri].v[ noelTriEdge[ifatri][0] ];
					const unsigned int iver1 = aTri[itri].v[ noelTriEdge[ifatri][1] ];
					if( iver0 == bar.v[0] && iver1 == bar.v[1] ){
                        assert( bar_ary.id_lr[0] == id_this_loop );
						bar.s2[0] = itri;
						bar.r2[0] = ifatri;
						aTri[itri].r2[ifatri] = 0;
					}
					else{
						assert( iver0 == bar.v[1] && iver1 == bar.v[0] );
                        assert( bar_ary.id_lr[1] == id_this_loop );
						bar.s2[1] = itri;
						bar.r2[1] = ifatri;
						aTri[itri].r2[ifatri] = 1;
					}
				}
				else{
					std::cout << "Error!-->Not defined type" << itype0 << std::endl;
					assert(0);
				}
			}
		}
	}

	{
		unsigned int iloc,itype;
		if( !this->FindElemLocType_CadIDType(iloc,itype, id_cad_l,Cad::LOOP) ) assert(0);
		assert( itype == 2 );
		assert( iloc < m_aTriAry.size() );
		assert( m_aTriAry[iloc].id_l_cad == id_cad_l );
		m_aTriAry[iloc].m_aTri = aTri;
	}

	return true;
}

bool CMesher2D::Tesselate_LoopAround
(const Cad::ICad2D_Msh& cad_2d,const unsigned int id_l)
{
	this->aVec2D.reserve(256);
	aVec2D.reserve(256);

	// Make Points belong to Loop
	{
    std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
		Cad::IItrLoop& itrEdgeLoop = *pItrEdgeLoop;
		for(;!itrEdgeLoop.IsEnd();itrEdgeLoop++){
			const unsigned int id_v = itrEdgeLoop.GetIdVertex();
			assert( this->GetElemID_FromCadID(id_v,Cad::VERTEX) == 0 );
			const unsigned int id_add = this->GetFreeObjID();
			const CVector2D& vec2d = cad_2d.GetVertexCoord(id_v);
			aVec2D.push_back(vec2d);
			{
				Msh::SVertex v;
				v.id_v_cad = id_v;
				v.id = id_add;
				v.v = aVec2D.size()-1;
				this->m_aVertex.push_back(v);
			}
			{
				this->m_ElemLoc.resize( id_add+1,-1 );
				this->m_ElemType.resize( id_add+1 );
				this->m_ElemLoc[id_add] = m_aVertex.size()-1;
				this->m_ElemType[id_add] = 0;
			}
			assert( this->CheckMesh() == 0 );
		}
	}

	// Tessalation Edge
	{
    std::auto_ptr<Cad::IItrLoop> pItrEdgeLoop = cad_2d.GetPtrItrLoop(id_l);
		Cad::IItrLoop& itrEdgeLoop = *pItrEdgeLoop;
		for(;!itrEdgeLoop.IsEnd();itrEdgeLoop++){
			unsigned int id_e;
			bool is_same_dir;
			if( !itrEdgeLoop.GetIdEdge(id_e,is_same_dir) ){
				assert(0);
			}
			Tessalate_Edge(cad_2d,id_e);
		}
	}

	// Tessalation Loop
	std::vector<Msh::CTriAry2D> aTriAry;
	aTriAry.reserve(16);
	this->Tesselate_Loop(cad_2d,id_l);

	return true;
}

bool CMesher2D::Tesselation(const Cad::ICad2D_Msh& cad_2d, const std::vector<unsigned int>& aIdLoop)
{	
	this->aVec2D.reserve(cad_2d.GetAryElemID(Cad::VERTEX).size()*10);

	{	// Make Points belong to Loop
		const std::vector<unsigned int>& aIdV = cad_2d.GetAryElemID(Cad::VERTEX);
		for(unsigned int iid_v=0;iid_v<aIdV.size();iid_v++){
			const unsigned int id_v = aIdV[iid_v];
			assert( this->GetElemID_FromCadID(id_v,Cad::VERTEX)==0 );      
			const unsigned int id_add = this->GetFreeObjID();
			const CVector2D& vec2d = cad_2d.GetVertexCoord(id_v);
			aVec2D.push_back(vec2d);
			{
				Msh::SVertex tmp_ver;
				tmp_ver.id = id_add;
				tmp_ver.id_v_cad = id_v;
				tmp_ver.v = aVec2D.size()-1;
				this->m_aVertex.push_back( tmp_ver );
			}
			{
				this->m_ElemLoc.resize( id_add+1,-1 );
				this->m_ElemType.resize( id_add+1 );
				this->m_ElemLoc[id_add] = m_aVertex.size()-1;
				this->m_ElemType[id_add] = 0;
			}
			assert( this->CheckMesh() == 0 );
		}
	}
	{	// Tessalation Edge
		const std::vector<unsigned int>& id_e_ary = cad_2d.GetAryElemID(Cad::EDGE);
		for(unsigned int iid_e=0;iid_e<id_e_ary.size();iid_e++){
			const unsigned int id_e = id_e_ary[iid_e];
			this->Tessalate_Edge(cad_2d, id_e);
			assert( this->CheckMesh() == 0 );
		}
	}
	{	// Tessalation Loop
		for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++){
			const unsigned int id_l = aIdLoop[iid_l];
			this->Tesselate_Loop(cad_2d, id_l);
			assert( this->CheckMesh() == 0 );
		}
	}

	this->MakeIncludeRelation(cad_2d);

	return true;
}

bool CMesher2D::Meshing_ElemLength
(const Cad::ICad2D_Msh& cad_2d, const double len, const std::vector<unsigned int>& aIdLoop)
{
  this->ClearMeshData();
	assert( len > 0.0 ); 	
	this->aVec2D.reserve(256);

  {   // ループに使われているVtxをメッシュに生成してフラグを0から1に立てる
		std::vector<unsigned int> aFlgVtx;
    for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++){
			const unsigned int id_l = aIdLoop[iid_l];
      std::auto_ptr<Cad::IItrLoop> pItr = cad_2d.GetPtrItrLoop(id_l);
			Cad::IItrLoop& itr = *pItr;
			for(;;){
				for(;!itr.IsEnd();itr++){
					unsigned int id_v = itr.GetIdVertex();
					if( aFlgVtx.size() <= id_v ){
						aFlgVtx.resize(id_v+1,0);
					}
					aFlgVtx[id_v] = 1;
				}
				if( !itr.ShiftChildLoop() ) break;
			}
		}
		for(unsigned int id_v=0;id_v<aFlgVtx.size();id_v++){
			if( aFlgVtx[id_v] == 0 ) continue;
			const unsigned int id_add = this->GetFreeObjID();
			const CVector2D& vec2d = cad_2d.GetVertexCoord(id_v);
			aVec2D.push_back(vec2d);
			{
				Msh::SVertex tmp_ver;
				tmp_ver.id = id_add;
				tmp_ver.id_v_cad = id_v;
				tmp_ver.ilayer = cad_2d.GetLayer(Cad::VERTEX,id_v);
				tmp_ver.v = aVec2D.size()-1;
				this->m_aVertex.push_back( tmp_ver );
			}
			{
				this->m_ElemLoc.resize( id_add+1,-1 );
				this->m_ElemType.resize( id_add+1 );
				this->m_ElemLoc[id_add] = m_aVertex.size()-1;
				this->m_ElemType[id_add] = 0;
			}
		}
		assert( this->CheckMesh() == 0 );
	}

  for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++)
  {	// ループに必要な辺を作る
		const unsigned int id_l = aIdLoop[iid_l];
		for(std::auto_ptr<Cad::IItrLoop> pItr=cad_2d.GetPtrItrLoop(id_l);!pItr->IsEndChild();pItr->ShiftChildLoop()){
			for(pItr->Begin();!pItr->IsEnd();(*pItr)++){				
				unsigned int id_e;   bool is_same_dir;
				if( !pItr->GetIdEdge(id_e,is_same_dir) ) continue;
				if( this->GetElemID_FromCadID(id_e,Cad::EDGE) != 0 ){ continue; }	// 既にこの辺はMeshに存在
        this->MakeMesh_Edge(cad_2d, id_e, len);
				assert( this->CheckMesh() == 0 );
			}
		}
	}
	
  for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++)
  {	// ループを作る
		const unsigned int id_l = aIdLoop[iid_l];
		this->MakeMesh_Loop(cad_2d, id_l, len);
		assert( this->CheckMesh() == 0 );
	}

	this->MakeIncludeRelation(cad_2d);
	return true;
}

bool CMesher2D::Meshing_ElemSize
(const Cad::ICad2D_Msh& cad_2d, unsigned int esize, const std::vector<unsigned int>& aIdLoop)
{
	assert( esize != 0 );
	if( esize == 0 ) return false;

	double area = 0;
	for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++){
		const unsigned int id_l = aIdLoop[iid_l];
		area += cad_2d.GetArea_Loop(id_l);
	}
	const double elen = sqrt( area / (double)esize ) * 1.4;

  return CMesher2D::Meshing_ElemLength(cad_2d,elen,aIdLoop);
}

bool CMesher2D::GetTriMesh_Around(
	std::vector< CTriAround >& aTriAround,
	////////////////
	unsigned int id_msh_bar0, 
	unsigned int ibar0, unsigned int inobar0, bool is_left,
	const std::vector<unsigned int>& aFlgMshCut ) const
{
	unsigned int ivec_cent;
	int id_msh_tri0;
	unsigned int itri0, ied0;
	{
		assert( aFlgMshCut[id_msh_bar0] == 1 );
		int ibarary0 = m_ElemLoc[id_msh_bar0];
        assert( ibarary0 >= 0 && ibarary0 < (int)m_aBarAry.size() );
		const std::vector<SBar>& aBar = m_aBarAry[ibarary0].m_aBar;
		assert( id_msh_bar0 == m_aBarAry[ibarary0].id );
		////////////////////////////////
		ivec_cent = aBar[ibar0].v[inobar0];
		const unsigned int iside0 = ( is_left ) ? 0 : 1;
		id_msh_tri0 = m_aBarAry[ibarary0].id_lr[iside0];
		itri0 = aBar[ibar0].s2[iside0];
		ied0  = aBar[ibar0].r2[iside0];
		int itriary0 = m_ElemLoc[id_msh_tri0];
        assert( itriary0 >= 0 && itriary0 < (int)m_aTriAry.size() );
		const std::vector<STri2D>& aTri0 = m_aTriAry[itriary0].m_aTri;
        assert( aTri0[itri0].g2[ied0] == (int)id_msh_bar0 );
		assert( aTri0[itri0].s2[ied0] == ibar0 );
	}

	unsigned int inoTriEdge_d, inoTriEdge_c;
	if( is_left == (inobar0==0) ){	// 反時計回り
			inoTriEdge_d = 1;	inoTriEdge_c = 0; }
	else{	inoTriEdge_d = 0;	inoTriEdge_c = 1; }
	
	////////////////
	aTriAround.clear();
	for(;;){
		unsigned int inotri_d = noelTriEdge[ied0][inoTriEdge_d];
		unsigned int inotri_c = noelTriEdge[ied0][inoTriEdge_c];	
		int itriary0 = m_ElemLoc[id_msh_tri0];
        assert( itriary0 >= 0 && itriary0 < (int)m_aTriAry.size() );
		const std::vector<STri2D>& aTri0 = m_aTriAry[itriary0].m_aTri;
		assert( aTri0[itri0].v[inotri_c] == ivec_cent );
		aTriAround.push_back( CTriAround(id_msh_tri0,itri0,inotri_c) );
		////////////////
		int id_msh_tri1;
		unsigned int itri1, ied1;
		if( aTri0[itri0].g2[inotri_d] != -2 ){
			const int id_msh_adj = aTri0[itri0].g2[inotri_d];
            if( id_msh_adj == -1 ){ return false; }	// 一周しない(外側にぶつかった)
            if( id_msh_adj == (int)id_msh_bar0 && aTri0[itri0].s2[inotri_d] == ibar0 ){
				return true;	// 一周
			}
            assert( id_msh_adj >= 0 && id_msh_adj < (int)aFlgMshCut.size() );
            if( aFlgMshCut[id_msh_adj] == 1 ){ return false; }	// 一周しない(CutMeshにぶつかった)
			assert( m_ElemType[id_msh_adj] == 1 );
			int iloc_adj = m_ElemLoc[id_msh_adj];	assert( iloc_adj != - 1);
			const std::vector<SBar>& aBar_adj = m_aBarAry[(unsigned int)iloc_adj].m_aBar;
			const unsigned int ibar_adj  = aTri0[itri0].s2[inotri_d];
			const unsigned int iside_adj = aTri0[itri0].r2[inotri_d];
            assert( (int)m_aBarAry[(unsigned int)iloc_adj].id_lr[iside_adj] == id_msh_tri0 );
			assert( aBar_adj[ibar_adj].s2[iside_adj] == itri0 );
			assert( aBar_adj[ibar_adj].r2[iside_adj] == inotri_d );
			const unsigned int jside_adj = 1-iside_adj;
            if( m_aBarAry[(unsigned int)iloc_adj].id_lr[jside_adj] == 0 ){    // 一周しない(外側にぶつかった)
                return false;
            }
			id_msh_tri1 = m_aBarAry[(unsigned int)iloc_adj].id_lr[jside_adj];
            assert( id_msh_tri1 >= 0 && id_msh_tri1 < (int)m_ElemType.size() );
			assert( m_ElemType[id_msh_tri1] == 2 );
			itri1 = aBar_adj[ibar_adj].s2[jside_adj];
			ied1 = aBar_adj[ibar_adj].r2[jside_adj];
			int itriary1 = m_ElemLoc[id_msh_tri1];
            assert( itriary1 >= 0 && itriary1 < (int)m_aTriAry.size() );
			const std::vector<STri2D>& aTri1 = m_aTriAry[itriary1].m_aTri;
			assert( aTri1[itri1].g2[ied1] == id_msh_adj );
			assert( aTri1[itri1].s2[ied1] == ibar_adj );
			assert( aTri1[itri1].r2[ied1] == jside_adj );
		}
		else{
			id_msh_tri1 = id_msh_tri0;
			itri1 = aTri0[itri0].s2[inotri_d];
			ied1 = relTriTri[ aTri0[itri0].r2[inotri_d] ][inotri_d];
			assert( aTri0[itri1].s2[ied1] == itri0 );
		}
		id_msh_tri0 = id_msh_tri1;
		itri0 = itri1;
		ied0 = ied1;
	}
	assert(0);
	return true;
}

// 辺で分離したメッシュを生成する
// 辺の両隣の３角形メッシュのIDは同じであるとする
// 辺の端点には両隣の３角形メッシュ以外
bool CMesher2D::GetClipedMesh(
		std::vector< std::vector<int> >& aLnods,
		std::vector<unsigned int>& mapVal2Co,
		const std::vector<unsigned int>& aIdMsh_Inc,
		const std::vector<unsigned int>& aIdMshBar_Cut ) const
{
	if( aIdMsh_Inc.empty() ) return true;
	const unsigned int nTriAryInc =aIdMsh_Inc.size();

	////////////////
	std::vector<int> aTriId2Inc;	// 三角形のIdからIncのmapを作る．対応がなければ-1が返る
	// aLnodsを作る
	{
		unsigned int max_id = this->FindMaxID();
		aTriId2Inc.resize( max_id+1,-1 );
		aLnods.resize(nTriAryInc);
		for(unsigned int iidmsh_inc=0;iidmsh_inc<aIdMsh_Inc.size();iidmsh_inc++){
			unsigned int id_msh = aIdMsh_Inc[iidmsh_inc];
//			std::cout << "IdMshInc : " << id_msh << std::endl;
			if( id_msh >= m_ElemLoc.size() ) return false;
			int iloc = m_ElemLoc[id_msh];
            if( iloc == -1 ) return false;  // aIdMsh_Incは存在するメッシュIDでなければならない
            if( m_ElemType[id_msh] != 2 ) return false; // aIdMsh_Incは三角形のメッシュIDでなければならない
			aTriId2Inc[id_msh] = iidmsh_inc;
            if( iloc<0 || iloc>=(int)m_aTriAry.size() ){ return false; }
			unsigned int ntri = m_aTriAry[iloc].m_aTri.size();
			aLnods[iidmsh_inc].resize(ntri*3,-1);
		}
	}

	////////////////
    // カットする辺メッシュが1，カットしないメッシュが０である配列が入っている
	std::vector<unsigned int> aFlgMshCut;
	{	
		unsigned int max_id = this->FindMaxID();
		aFlgMshCut.resize( max_id+1,0 );
		for(unsigned int iid_msh_cut=0;iid_msh_cut<aIdMshBar_Cut.size();iid_msh_cut++){
			unsigned int id_msh_cut = aIdMshBar_Cut[iid_msh_cut];
            int iloc = m_ElemLoc[id_msh_cut];
            if( iloc == -1 ) return false;  // aIdMsh_Incは存在するメッシュIDでなければならない
            if( m_ElemType[id_msh_cut] != 1 ) return false; // aIdMsh_Incは線メッシュIDでなければならない
            if( iloc<0 || iloc>=(int)m_aBarAry.size() ){ return false; }
//			std::cout << " id_msh_cut " << id_msh_cut << std::endl;
			aFlgMshCut[id_msh_cut] = 1;
		}
	}	

	////////////////
	mapVal2Co.resize(aVec2D.size());
	for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){ mapVal2Co[ivec] = ivec; }

    std::vector<unsigned int> aDupli(aVec2D.size(),0);
	////////////////////////////////
	for(unsigned int iid_msh_cut=0;iid_msh_cut<aIdMshBar_Cut.size();iid_msh_cut++)
	{
		unsigned int id_msh_cut = aIdMshBar_Cut[iid_msh_cut];
		int ibarary = m_ElemLoc[id_msh_cut];
        if( ibarary<0 || ibarary>=(int)m_aBarAry.size() ){ return false; }
        const CBarAry& BarAry = m_aBarAry[ibarary];
        {
            unsigned int id_msh_l = BarAry.id_lr[0];
            unsigned int id_msh_r = BarAry.id_lr[1];
            if( aTriId2Inc[id_msh_l]==-1 ) continue;
            if( aTriId2Inc[id_msh_r]==-1 ) continue;
        }
        const std::vector<SBar>& aBar = BarAry.m_aBar;
        ////////////////////////////////
        std::vector< CTriAround > aTriAround0, aTriAround1;
		for(unsigned int ivecbar=0;ivecbar<aBar.size()+1;ivecbar++)
		{
			unsigned int ibar, inobar;
			if( ivecbar == 0 ){ ibar=0; inobar=0; }
			else{ ibar=ivecbar-1; inobar=1; }
			const unsigned int ivec0 = aBar[ibar].v[inobar];
			if( !this->GetTriMesh_Around(aTriAround0,  id_msh_cut,ibar,inobar,true,  aFlgMshCut) ){
				// この点はこの辺要素周りに分離されている
				// 左側分離
				if( aDupli[ivec0] == 0 ){
					for(unsigned int itria=0;itria<aTriAround0.size();itria++){
						const int id_msh_tri0 = aTriAround0[itria].id_msh;
						const unsigned int itri0 = aTriAround0[itria].ielem;
						const unsigned int inotri0 = aTriAround0[itria].inoel;
                        assert( id_msh_tri0 >= 0 && id_msh_tri0 < (int)m_ElemLoc.size() );
						assert( m_ElemType[id_msh_tri0] == 2 );
						int itriary0 = m_ElemLoc[id_msh_tri0];
                        assert( itriary0 >= 0 && itriary0 < (int)m_aTriAry.size() );
						const std::vector<STri2D>& aTri0 = m_aTriAry[itriary0].m_aTri;
						assert( aTri0[itri0].v[inotri0] == ivec0 );
						int iidmsh_inc = aTriId2Inc[id_msh_tri0];
                        assert( iidmsh_inc >=0 && iidmsh_inc < (int)aLnods.size() );
						assert( aLnods[iidmsh_inc][itri0*3+inotri0] == -1 );
						aLnods[iidmsh_inc][itri0*3+inotri0] = aTri0[itri0].v[inotri0];	// この三角形頂点は今後考慮しませんフラグを立てる
					}
					aDupli[ivec0] += 1;
				}
				else{
					const unsigned int ivec_add = mapVal2Co.size();
					bool iflag = false;
					for(unsigned int itria=0;itria<aTriAround0.size();itria++){
						const int id_msh_tri0 = aTriAround0[itria].id_msh;
						const unsigned int itri0 = aTriAround0[itria].ielem;
						const unsigned int inotri0 = aTriAround0[itria].inoel;
                        assert( id_msh_tri0 >= 0 && id_msh_tri0 < (int)m_ElemLoc.size() );
						assert( m_ElemType[id_msh_tri0] == 2 );
						int itriary0 = m_ElemLoc[id_msh_tri0];
                        assert( itriary0 >= 0 && itriary0 < (int)m_aTriAry.size() );
						const std::vector<STri2D>& aTri0 = m_aTriAry[itriary0].m_aTri;
						assert( aTri0[itri0].v[inotri0] == ivec0 );
						int iidmsh_inc = aTriId2Inc[id_msh_tri0];
                        assert( iidmsh_inc >=0 && iidmsh_inc < (int)aLnods.size() );
						if( aLnods[iidmsh_inc][itri0*3+inotri0] != -1 ){ assert(itria==0); break; }
						iflag = true;
						aLnods[iidmsh_inc][itri0*3+inotri0] = ivec_add;	// この三角形頂点は今後考慮しませんフラグを立てる
					}
					if( iflag ){ 
						mapVal2Co.resize( mapVal2Co.size()+1, ivec0 ); 
						aDupli[ivec0] += 1;
					}
				}
				{	// 右側分離
					bool res = this->GetTriMesh_Around(aTriAround1,  id_msh_cut,ibar,inobar,false,  aFlgMshCut);
					assert( !res );	// 左側で一周回らなければ右側で一周回らないはず
					const unsigned int ivec_add = mapVal2Co.size();
					bool iflag = false;
					for(unsigned int itria=0;itria<aTriAround1.size();itria++){
						const int id_msh_tri0 = aTriAround1[itria].id_msh;
						const unsigned int itri0 = aTriAround1[itria].ielem;
						const unsigned int inotri0 = aTriAround1[itria].inoel;
                        assert( id_msh_tri0 >= 0 && id_msh_tri0 < (int)m_ElemLoc.size() );
						assert( m_ElemType[id_msh_tri0] == 2 );
						int itriary0 = m_ElemLoc[id_msh_tri0];
                        assert( itriary0 >= 0 && itriary0 < (int)m_aTriAry.size() );
						const std::vector<STri2D>& aTri0 = m_aTriAry[itriary0].m_aTri;
						assert( aTri0[itri0].v[inotri0] == ivec0 );
						int iidmsh_inc = aTriId2Inc[id_msh_tri0];
                        assert( iidmsh_inc >=0 && iidmsh_inc < (int)aLnods.size() );
						if( aLnods[iidmsh_inc][itri0*3+inotri0] != -1 ){ assert(itria==0); break; }
						iflag = true;
						aLnods[iidmsh_inc][itri0*3+inotri0] = ivec_add;	// この三角形頂点は今後考慮しませんフラグを立てる
					}
					if( iflag ){ 
						mapVal2Co.resize( mapVal2Co.size()+1, ivec0 ); 
						aDupli[ivec0] += 1;
					}
				}
			}
		}
	}	
	////////////////////////////////
	
	{
		unsigned int max_id = this->FindMaxID();
		aTriId2Inc.resize( max_id+1,-1 );
		aLnods.resize(nTriAryInc);
		for(unsigned int iidmsh_inc=0;iidmsh_inc<aIdMsh_Inc.size();iidmsh_inc++){
			unsigned int id_msh = aIdMsh_Inc[iidmsh_inc];
			if( id_msh >= m_ElemLoc.size() ) return false;
			int iloc = m_ElemLoc[id_msh];
			if( iloc == -1 ) return false;
			if( m_ElemType[id_msh] != 2 ) return false;
			aTriId2Inc[id_msh] = iidmsh_inc;
			const std::vector<STri2D>& aTri = m_aTriAry[iloc].m_aTri;
			unsigned int ntri = aTri.size();
			for(unsigned int itri=0;itri<ntri;itri++){
				for(unsigned int inotri=0;inotri<3;inotri++){
					if( aLnods[iidmsh_inc][itri*3+inotri] != -1 ) continue;
					aLnods[iidmsh_inc][itri*3+inotri] = aTri[itri].v[inotri];
				}
			}
		}
	}

//	std::cout << aVec2D.size() << " " << mapVal2Co.size() << std::endl;

	return true;
}


bool CMesher2D::Serialize( Com::CSerializer& arch, bool is_only_cadmshlink )
{
	if( arch.IsLoading() ){	// 読み込み時の処理
        CMesher2D::Clear();
		const unsigned int buff_size = 256;
		char class_name[buff_size];
		arch.ReadDepthClassName(class_name,buff_size);
		assert( strncmp(class_name,"CMesher2D",9) == 0 );
		arch.ShiftDepth(true);

		////////////////////////////////
		// CADとの接続関係のロード

        {
            arch.ReadDepthClassName(class_name,buff_size); 
			assert( strncmp(class_name,"setIdLCad_CutMesh",17)==0);
            int nl;
            arch.Get("%d",&nl);   assert( nl >= 0 );
            unsigned int ind, id_l;
            setIdLCad_CutMesh.clear();
            for(unsigned int il=0;il<(unsigned int)nl;il++){
                arch.Get("%d%d",&ind,&id_l);   assert( ind == il );
                setIdLCad_CutMesh.insert(id_l);
            }
        }
		{	// メッシュ生成モードのロード
			int ntmp0, ntmp1;
			double dtmp0;
			arch.Get("%d%d%lf",&ntmp0,&ntmp1,&dtmp0);
			assert( ntmp0 >= 0 && ntmp0 < 3 );
			m_imode_meshing = ntmp0; 
			assert( ntmp1 > 0 );
			m_esize = ntmp1;
			assert( dtmp0 > 0 );
			m_elen = dtmp0;
		}
		if( is_only_cadmshlink ){
			arch.ShiftDepth(false);
			return true;
		}

		////////////////////////////////
		// メッシュ情報のロード

		{	// 座標をロード
			int nvec, ndim;	arch.Get("%d%d",&nvec,&ndim);	assert(nvec>0 && (ndim>0&&ndim<4) );
			this->aVec2D.resize(nvec);
            for(unsigned int ivec=0;ivec<(unsigned int)nvec;ivec++){
				int itmp0;
				double x,y;
                arch.Get("%d%lf%lf",&itmp0,&x,&y);	assert( itmp0 == (int)ivec );
				aVec2D[ivec].x = x;
				aVec2D[ivec].y = y;
			}
		}
		int naVer, naBar, naTri, naQuad;
		arch.Get("%d%d%d%d",&naVer,&naBar,&naTri,&naQuad);
        assert( naVer>=0 && naBar>=0 && naTri>=0 && naQuad>=0 );
        for(unsigned int iaVer=0;iaVer<(unsigned int)naVer;iaVer++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"SVertex",7)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int iv;		arch.Get("%d",&iv);		assert(iv>=0);
			const unsigned int iver = m_aVertex.size();
			m_aVertex.resize( m_aVertex.size()+1 );
			SVertex& Ver = m_aVertex[iver];
			Ver.id = id;
			Ver.id_v_cad = id_cad;
			Ver.v = iv;
		}
        for(unsigned int iaBar=0;iaBar<(unsigned int)naBar;iaBar++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CBarAry",7)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int id_s, id_e; arch.Get("%d %d",&id_s,&id_e);
			int id_l, id_r; arch.Get("%d %d",&id_l,&id_r);			
			int nbar;	arch.Get("%d",&nbar);	assert(nbar>0);
			const unsigned int ibar_ary = m_aBarAry.size();
			m_aBarAry.resize( m_aBarAry.size()+1 );
			CBarAry& aBar = m_aBarAry[ibar_ary];
			aBar.id = id;
			aBar.id_e_cad = id_cad;
			aBar.id_se[0] = id_s;	aBar.id_se[1] = id_e;
			aBar.id_lr[0] = id_l;	aBar.id_lr[1] = id_r;
			aBar.m_aBar.resize(nbar);
			for(int ibar=0;ibar<nbar;ibar++){
				int tmp_ibar,iv0,iv1, s0,s1, r0,r1;
				arch.Get("%d%d%d%d%d%d%d",&tmp_ibar, &iv0,&iv1, &s0,&s1, &r0,&r1);
				assert( tmp_ibar == ibar );
				assert( iv0>=0 && iv1>=0 );
				aBar.m_aBar[ibar].v[ 0] = iv0;
				aBar.m_aBar[ibar].v[ 1] = iv1;
				aBar.m_aBar[ibar].s2[0] = s0;
				aBar.m_aBar[ibar].s2[1] = s1;
				aBar.m_aBar[ibar].r2[0] = r0;
				aBar.m_aBar[ibar].r2[1] = r1;
			}
		}
        for(unsigned int iaTri=0;iaTri<(unsigned int)naTri;iaTri++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CTriAry2D",9)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int ntri;	arch.Get("%d",&ntri);	assert(ntri>0);
			const unsigned int itri_ary = m_aTriAry.size();
			m_aTriAry.resize( m_aTriAry.size()+1 );
			CTriAry2D& aTri = m_aTriAry[itri_ary];
			aTri.id = id;
			aTri.id_l_cad = id_cad;
			aTri.m_aTri.resize(ntri);
			for(int itri=0;itri<ntri;itri++){
				int tmp_itri,iv0,iv1,iv2;
				arch.Get("%d %d %d %d",&tmp_itri,&iv0,&iv1,&iv2);
				assert( tmp_itri == itri );
				assert( iv0>=0 && iv1>=0 && iv2>=0 );
				aTri.m_aTri[itri].v[0] = iv0;
				aTri.m_aTri[itri].v[1] = iv1;
				aTri.m_aTri[itri].v[2] = iv2;
			}
			{
				unsigned int npoin = aVec2D.size();
				unsigned int* elsup_ind = new unsigned int [npoin+1];
				unsigned int nelsup;
				unsigned int* elsup;
				Msh::MakePointSurTri(aTri.m_aTri,npoin,elsup_ind,nelsup,elsup);
				Msh::MakeInnerRelationTri(aTri.m_aTri,npoin,elsup_ind,nelsup,elsup);
				delete[] elsup_ind;
				delete[] elsup;
			}
		}
        for(unsigned int iaQuad=0;iaQuad<(unsigned int)naQuad;iaQuad++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CQuadAry2D",10)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int nquad;	arch.Get("%d",&nquad);	assert(nquad>0);
			const unsigned int iquad_ary = m_aQuadAry.size();
			m_aQuadAry.resize( m_aQuadAry.size()+1 );
			CQuadAry2D& aQuad = m_aQuadAry[iquad_ary];
			aQuad.id = id;
			aQuad.id_l_cad = id_cad;
			aQuad.m_aQuad.resize(nquad);
			for(int iquad=0;iquad<nquad;iquad++){
				int tmp_iquad,iv0,iv1,iv2,iv3;
				arch.Get("%d %d %d %d %d",&tmp_iquad,&iv0,&iv1,&iv2,&iv3);
				assert( tmp_iquad == iquad );
				assert( iv0>=0 && iv1>=0 && iv2>=0 && iv3>=0 );
				aQuad.m_aQuad[iquad].v[0] = iv0;
				aQuad.m_aQuad[iquad].v[1] = iv1;
				aQuad.m_aQuad[iquad].v[2] = iv2;
				aQuad.m_aQuad[iquad].v[3] = iv3;
			}
		}
		this->MakeElemLocationType();
		////////////////////////////////
		// build the internal relationship
		for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
			const unsigned int id_msh_bar = m_aBarAry[ibarary].id;
			const int iloc_bar = this->m_ElemLoc[id_msh_bar];
			const std::vector<SBar>& aBar = this->m_aBarAry[iloc_bar].m_aBar;
			for(unsigned int isidebar=0;isidebar<2;isidebar++){
				int id_msh_adj = this->m_aBarAry[iloc_bar].id_lr[isidebar];
				if( id_msh_adj <= 0 ) continue;
//				std::cout << id_msh_adj << " " << m_ElemLoc.size() << std::endl;
				assert(id_msh_adj<(int)m_ElemLoc.size());
				const int iloc_adj = this->m_ElemLoc[id_msh_adj];
				if( this->m_ElemType[id_msh_adj] == 2 ){	// tri mesh
					std::vector<STri2D>& aTri = this->m_aTriAry[iloc_adj].m_aTri;
					for(unsigned int ibar=0;ibar<aBar.size();ibar++){
						unsigned int itri   = aBar[ibar].s2[isidebar];
						unsigned int inotri = aBar[ibar].r2[isidebar];
						aTri[itri].g2[inotri] = (int)id_msh_bar;
						aTri[itri].s2[inotri] = ibar;
						aTri[itri].r2[inotri] = isidebar;
					}
				}
			}
		}
		assert( this->CheckMesh() == 0 );
		arch.ShiftDepth(false);
		return true;
	}
	else{ // write file
		arch.WriteDepthClassName("CMesher2D");	
		arch.ShiftDepth(true);
		
		////////////////////////////////
		// write relation to CAD

		{
			arch.WriteDepthClassName("setIdLCad_CutMesh");
			arch.Out("%d\n",setIdLCad_CutMesh.size());
			std::set<unsigned int>::iterator itr = setIdLCad_CutMesh.begin();
			unsigned int icnt = 0;
			for(;itr!=setIdLCad_CutMesh.end();itr++){
				unsigned int id_l = *itr;
				arch.Out("%d %d\n",icnt,id_l);
				icnt++;
			}
		}
		arch.Out("%d %d %lf\n",m_imode_meshing,m_esize,m_elen);
		if( is_only_cadmshlink ){ // 
			arch.ShiftDepth(false); 
			return true; 
		}

		////////////////////////////////
		// Write information of Msh

		arch.Out("%d %d\n",aVec2D.size(),2);
		for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){
			arch.Out("%d %lf %lf\n",ivec,aVec2D[ivec].x,aVec2D[ivec].y);
		}
		arch.Out("%d %d %d %d\n",m_aVertex.size(), m_aBarAry.size(), m_aTriAry.size(), m_aQuadAry.size() );
		{	// Vertexの出力
			for(unsigned int iver=0;iver<m_aVertex.size();iver++){
				arch.WriteDepthClassName("SVertex");
				arch.Out("%d\n",m_aVertex[iver].id);
				arch.Out("%d\n",m_aVertex[iver].id_v_cad);
				arch.Out("%d\n",m_aVertex[iver].v);
			}
		}
		{	// Barの出力
			for(unsigned int ibar_ary=0;ibar_ary<m_aBarAry.size();ibar_ary++){
				arch.WriteDepthClassName("CBarAry");
				arch.Out("%d\n",m_aBarAry[ibar_ary].id);
				arch.Out("%d\n",m_aBarAry[ibar_ary].id_e_cad);
				arch.Out("%d %d\n",m_aBarAry[ibar_ary].id_se[0],m_aBarAry[ibar_ary].id_se[1]);
				arch.Out("%d %d\n",m_aBarAry[ibar_ary].id_lr[0],m_aBarAry[ibar_ary].id_lr[1]);				
				arch.Out("%d\n",m_aBarAry[ibar_ary].m_aBar.size());
				const std::vector<SBar>& aBar = m_aBarAry[ibar_ary].m_aBar;
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					arch.Out("%d %d %d  %d %d  %d %d\n",ibar,
						aBar[ibar].v[0], aBar[ibar].v[1],
						aBar[ibar].s2[0],aBar[ibar].s2[1],
						aBar[ibar].r2[0],aBar[ibar].r2[1]);
				}
			}
		}
		{	// Triの出力
			for(unsigned int itri_ary=0;itri_ary<m_aTriAry.size();itri_ary++){
				arch.WriteDepthClassName("CTriAry2D");
				arch.Out("%d\n",m_aTriAry[itri_ary].id);
				arch.Out("%d\n",m_aTriAry[itri_ary].id_l_cad);
				arch.Out("%d\n",m_aTriAry[itri_ary].m_aTri.size());
				const std::vector<STri2D>& aTri = m_aTriAry[itri_ary].m_aTri;
				for(unsigned int itri=0;itri<aTri.size();itri++){
					arch.Out("%d %d %d %d\n",itri,aTri[itri].v[0],aTri[itri].v[1],aTri[itri].v[2]);
				}
			}
		}
		{	// Quadの出力
			for(unsigned int iquad_ary=0;iquad_ary<m_aQuadAry.size();iquad_ary++){
				arch.WriteDepthClassName("CQuadAry2D");
				arch.Out("%d\n",m_aQuadAry[iquad_ary].id);
				arch.Out("%d\n",m_aQuadAry[iquad_ary].id_l_cad);
				arch.Out("%d\n",m_aQuadAry[iquad_ary].m_aQuad.size());
				const std::vector<SQuad2D>& aQuad = m_aQuadAry[iquad_ary].m_aQuad;
				for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
					arch.Out("%d %d %d %d %d\n",iquad,aQuad[iquad].v[0],aQuad[iquad].v[1],aQuad[iquad].v[2],aQuad[iquad].v[3]);
				}
			}
		}
		arch.ShiftDepth(false);
	}
	return true;
}


bool CMesher2D::ReadFromFile_GiDMsh(const std::string& file_name)
{
	FILE *fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];
	char stmp2[64], stmp3[64], stmp4[64], stmp6[64];

	if( (fp = ::fopen(file_name.c_str(),"r"))== NULL ){
		std::cout << "Cannot Open Gid File : " << std::endl;//<< file_name << std::endl;
		return false;
	}

	char str_elem_type[64];
	int ndim, nnoel;
	fgets(stmp1,buff_size,fp);
	sscanf(stmp1,"%s %s %d %s %s %s %d", 
		stmp2,stmp3,&ndim,stmp4,str_elem_type,stmp6,&nnoel);

	if( strcmp(str_elem_type,"Triangle")==0){
		if( nnoel != 3 ) return false;
	}
	else if( strcmp(str_elem_type,"Quadrilateral")==0 ){
		if( nnoel != 4 ) return false;
	}
	else{
		std::cout << "Error!-->Not defined elem type" << std::endl;
		return false;
	}

	////////////////////////////////
	// Change Data from here

    CMesher2D::Clear();	// 現在のメッシュ情報を全てクリア

	{	// 節点の読み込み
		fgets(stmp1,buff_size,fp);	// Coorinates
		std::vector<double> tmp_buffer;
		tmp_buffer.reserve(16384);
		unsigned int counter = 0;
		if( ndim == 3 ){	// 次元が３の場合で点群が同一平面上にある場合は２次元に縮退させる
			double v_x=0.0, v_y=0.0, v_z=0.0;
			double a_x=0.0, a_y=0.0, a_z=0.0;
			for(;;){
				fgets(stmp1,buff_size,fp);
				if( strncmp(stmp1 ,"end coordinates",15)==0 ) break;
				int inode;
				double co_x, co_y, co_z;
				sscanf(stmp1,"%d%lf%lf%lf",&inode, &co_x, &co_y, &co_z);
                assert( (int)counter+1 == inode );
				counter++;
				unsigned int i0 = tmp_buffer.size();
				tmp_buffer.resize( i0+3 );
				tmp_buffer[i0  ] = co_x;
				tmp_buffer[i0+1] = co_y;
				tmp_buffer[i0+2] = co_z;
				a_x+=co_x; a_y+=co_y; a_z+=co_z;
				v_x+=co_x*co_x; v_y+=co_y*co_y; v_z+=co_z*co_z;
			}
			const unsigned int nnode = counter;
			assert( nnode != 0 );
			v_x/=nnode; v_y/=nnode; v_z/=nnode;
			a_x/=nnode; a_y/=nnode; a_z/=nnode;
			v_x-=a_x*a_x; v_y-=a_y*a_y; v_z-=a_z*a_z;
			double v_min = v_x;
			unsigned int idim_min = 0;
			if( v_y < v_min ){ v_min=v_y; idim_min = 1; }
			if( v_z < v_min ){ v_min=v_z; idim_min = 2; }
			if( v_min > (v_x+v_y+v_z)*0.0001 ) return false;
            int idim0=0, idim1=0;
            if(      idim_min==0 ){ idim0=1; idim1=2; }
			else if( idim_min==1 ){ idim0=2; idim1=0; }
			else if( idim_min==2 ){ idim0=0; idim1=1; }
			aVec2D.resize(nnode);
			for(unsigned int inode=0;inode<nnode;inode++){
				aVec2D[inode].x = tmp_buffer[inode*3+idim0];
				aVec2D[inode].y = tmp_buffer[inode*3+idim1];
			}
		}
	}	// 節点の読み込み終了

	fgets(stmp1,buff_size,fp);	// 改行

	std::vector<SBar> aBar;	// 境界要素
	{	// 要素の読み込み
		fgets(stmp1,buff_size,fp);	// Elements
		std::vector<int> tmp_buffer;
		tmp_buffer.reserve(16384);
		if( strcmp(str_elem_type,"Triangle")==0 ){
			unsigned int counter = 0;
			for(;;){
				fgets(stmp1,buff_size,fp);
				if( strncmp(stmp1 ,"end elements", 12)==0 ) break;
				int ielem;
				int v1, v2, v3;
				sscanf(stmp1,"%d%d%d%d",&ielem, &v1, &v2, &v3);
                assert( (int)counter+1 == ielem );
				assert( v1>0 && v2>0 && v3>0 );
				counter++;
				unsigned int i0 = tmp_buffer.size();
				tmp_buffer.resize( i0+3 );
				tmp_buffer[i0  ] = v1-1;
				tmp_buffer[i0+1] = v2-1;
				tmp_buffer[i0+2] = v3-1;
			}
			const unsigned int ntri = counter;
			m_ElemLoc.resize(2);
			m_ElemLoc[1] = 0;
			m_ElemType.resize(2);
			m_ElemType[1] = 2;
			m_aTriAry.resize(1);
			m_aTriAry[0].id_l_cad = 0;
			m_aTriAry[0].id = 1;
			m_aTriAry[0].m_aTri.resize(ntri);
			for(unsigned int itri=0;itri<ntri;itri++){
				m_aTriAry[0].m_aTri[itri].v[0] = tmp_buffer[itri*3  ];
				m_aTriAry[0].m_aTri[itri].v[1] = tmp_buffer[itri*3+1];
				m_aTriAry[0].m_aTri[itri].v[2] = tmp_buffer[itri*3+2];
			}
			m_include_relation.resize(2);
			unsigned int npoin = aVec2D.size();
			{
				unsigned int* elsup_ind = new unsigned int [npoin+1];
				unsigned int nelsup;
				unsigned int* elsup;
				Msh::MakePointSurTri(m_aTriAry[0].m_aTri,npoin,elsup_ind,nelsup,elsup);
				Msh::MakeInnerRelationTri(m_aTriAry[0].m_aTri,npoin,elsup_ind,nelsup,elsup);
				delete[] elsup_ind;
				delete[] elsup;
			}
			Msh::CheckTri(m_aTriAry[0].m_aTri);
			Msh::MakeOuterBoundTri(m_aTriAry[0].m_aTri,aBar);
		}
		else if(  strcmp(str_elem_type,"Quadrilateral")==0 )
		{
			unsigned int counter = 0;
			for(;;){
				fgets(stmp1,buff_size,fp);
				if( strncmp(stmp1 ,"end elements", 12)==0 )break;
				int ielem;
				int v1, v2, v3, v4;
				sscanf(stmp1,"%d%d%d%d%d",&ielem, &v1, &v2, &v3, &v4);
                assert( (int)counter+1 == ielem );
				assert( v1>0 && v2>0 && v3>0 && v4>0 );
				counter++;
				unsigned int i0 = tmp_buffer.size();
				tmp_buffer.resize( i0+4 );
				tmp_buffer[i0  ] = v1-1;
				tmp_buffer[i0+1] = v2-1;
				tmp_buffer[i0+2] = v3-1;
				tmp_buffer[i0+3] = v4-1;
			}
			const unsigned int nquad = counter;
			m_ElemLoc.resize(2);
			m_ElemLoc[1] = 0;
			m_ElemType.resize(2);
			m_ElemType[1] = 3;
			m_aQuadAry.resize(1);
			m_aQuadAry[0].id_l_cad = 0;
			m_aQuadAry[0].id = 1;
			m_aQuadAry[0].m_aQuad.resize(nquad);
			for(unsigned int iquad=0;iquad<nquad;iquad++){
				m_aQuadAry[0].m_aQuad[iquad].v[0] = tmp_buffer[iquad*4  ];
				m_aQuadAry[0].m_aQuad[iquad].v[1] = tmp_buffer[iquad*4+1];
				m_aQuadAry[0].m_aQuad[iquad].v[2] = tmp_buffer[iquad*4+2];
				m_aQuadAry[0].m_aQuad[iquad].v[3] = tmp_buffer[iquad*4+3];
			}
			m_include_relation.resize(2);
			unsigned int npoin = aVec2D.size();
			{
				unsigned int* elsup_ind = new unsigned int [npoin+1];
				unsigned int nelsup;
				unsigned int* elsup;
				Msh::MakePointSurQuad(m_aQuadAry[0].m_aQuad,npoin,elsup_ind,nelsup,elsup);
				Msh::MakeInnerRelationQuad(m_aQuadAry[0].m_aQuad,npoin,elsup_ind,nelsup,elsup);
				delete[] elsup_ind;
				delete[] elsup;
			}
//			Msh::CheckQuad(m_aQuadAry[0].m_aQuad);	// 未実装
			Msh::MakeOuterBoundQuad(m_aQuadAry[0].m_aQuad,aBar);
		}
	}
	fclose(fp);
	////////////////////////////////
	// 境界辺の抽出
	{
		std::vector< std::vector<unsigned  int> > aIndBarAry;
		Msh::ColorCodeBarAry(aBar,aVec2D, aIndBarAry);
		for(unsigned int iibarary=0;iibarary<aIndBarAry.size();iibarary++){		
			unsigned int id_bar_ary_new = this->GetFreeObjID();
			m_aBarAry.resize( m_aBarAry.size()+1 );
			CBarAry& refBarAry = m_aBarAry[ m_aBarAry.size()-1 ];
			refBarAry.id_lr[0] = 0;
			refBarAry.id_lr[1] = 0;
			refBarAry.id_se[0] = 0;
			refBarAry.id_se[1] = 0;
			refBarAry.id = id_bar_ary_new;
			refBarAry.id_e_cad = 0;
			const unsigned int nbar = aIndBarAry[iibarary].size();
			refBarAry.m_aBar.reserve(nbar);
			for(unsigned int iibar=0;iibar<nbar;iibar++){
				const unsigned int ibar0 = aIndBarAry[iibarary][iibar];
				assert( ibar0 < aBar.size() );
				refBarAry.m_aBar.push_back( aBar[ibar0] );
			}
			{
				m_ElemLoc.resize(id_bar_ary_new+1);
				m_ElemType.resize(id_bar_ary_new+1);
				m_ElemLoc[id_bar_ary_new] = m_aBarAry.size()-1;
				m_ElemType[id_bar_ary_new] = 1;
			}
			m_include_relation[1].push_back( id_bar_ary_new );
		}
	}
	////////////////////////////////
	// 境界節点を抽出
	{
		std::map<unsigned int,unsigned int> aIndVec2IdVer;
		for(unsigned int iibarary=0;iibarary<m_aBarAry.size();iibarary++){
			const CBarAry& barary = m_aBarAry[iibarary];
			const unsigned int iv0 = barary.m_aBar[0].v[0];
			const unsigned int nbar = barary.m_aBar.size();
			const unsigned int iv1 = barary.m_aBar[nbar-1].v[1];
			aIndVec2IdVer.insert( std::make_pair(iv0,0) );
			aIndVec2IdVer.insert( std::make_pair(iv1,0) );
		}
		for(std::map<unsigned int,unsigned int>::iterator itr=aIndVec2IdVer.begin();itr!=aIndVec2IdVer.end();itr++){
			unsigned int ivec0 = itr->first;
			unsigned int id_vtx_new = this->GetFreeObjID();
			m_aVertex.resize( m_aVertex.size()+1 );
			SVertex& refVertex = m_aVertex[ m_aVertex.size()-1 ];
			refVertex.id = id_vtx_new;
			refVertex.id_v_cad = 0;
			refVertex.v = ivec0;
			{
				m_ElemLoc.resize(id_vtx_new+1);
				m_ElemType.resize(id_vtx_new+1);
				m_ElemLoc[id_vtx_new] = m_aVertex.size()-1;
				m_ElemType[id_vtx_new] = 0;
			}
			////////////////
			itr->second = id_vtx_new;
		}
		for(unsigned int iibarary=0;iibarary<m_aBarAry.size();iibarary++){
			CBarAry& barary = m_aBarAry[iibarary];
			const unsigned int id_bar = barary.id;
			const unsigned int iv0 = barary.m_aBar[0].v[0];
			const unsigned int nbar = barary.m_aBar.size();
			const unsigned int iv1 = barary.m_aBar[nbar-1].v[1];
			std::map<unsigned int,unsigned int>::iterator itr0 = aIndVec2IdVer.find(iv0);
			std::map<unsigned int,unsigned int>::iterator itr1 = aIndVec2IdVer.find(iv1);
			assert( itr0 != aIndVec2IdVer.end() );
			assert( itr1 != aIndVec2IdVer.end() );
			barary.id_se[0] = itr0->second;
			barary.id_se[1] = itr1->second;
			if( m_include_relation.size() <= id_bar ){ m_include_relation.resize(id_bar+1); }
			m_include_relation[id_bar].push_back( barary.id_se[0] );
			m_include_relation[id_bar].push_back( barary.id_se[1] );
		}
	}
	return true;
}
