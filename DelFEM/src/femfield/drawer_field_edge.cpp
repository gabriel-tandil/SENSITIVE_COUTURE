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
// CDrawerEdge.cpp : implementation of edge drawing class(CDrawerEdge)
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#  pragma warning ( disable : 4786 )
#endif

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")     /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <memory>


#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "delfem/drawer_field_edge.h"
#include "delfem/elem_ary.h"
#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

CDrawerEdge::CDrawerEdge(){
  this->line_width_ = 1;
	m_paVer = 0;
	m_nline = 0;
}

CDrawerEdge::CDrawerEdge(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world)
{
  this->line_width_ = 1;  
	m_paVer = 0;
	m_nline = 0;
	this->Set(id_field, isnt_value_disp, world);
}

CDrawerEdge::~CDrawerEdge(){
	if( m_paVer != 0 ){ delete m_paVer; }
}

Com::CBoundingBox3D CDrawerEdge::GetBoundingBox( double rot[] ) const
{
	return m_paVer->GetBoundingBox(rot);
}

bool CDrawerEdge::Update(const Fem::Field::CFieldWorld& world)
{
	assert( world.IsIdField(m_IdField) );
	const Fem::Field::CField& field = world.GetField(m_IdField);

	{	// 頂点配列をセット
		unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
		unsigned int id_ns_c_co = field.GetNodeSegInNodeAry(CORNER).id_ns_co;
		assert( id_na_c_co != 0 );
		const Fem::Field::CNodeAry& na_c_co = world.GetNA(id_na_c_co);
		const Fem::Field::CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);
		const unsigned int ndim_field = ns_c_co.Length();
		unsigned int ndim_draw;
		if( isnt_value_disp ){ 
			ndim_draw = ndim_field;
		}
		else{
			if( ndim_field == 2 && (field.GetFieldType()==SCALAR||field.GetFieldType()==ZSCALAR) ){ 
				ndim_draw = 3;
			}
			else{ ndim_draw = ndim_field; }
		}
		const unsigned int npoin = m_nline*2;
		if( m_paVer == 0 ){	m_paVer = new Com::View::CVertexArray(npoin,ndim_draw); }
		else if( m_paVer->NDim() != ndim_draw || m_paVer->NPoin() != npoin ){ 
			assert(0);
			return false; 
		}

		if( ndim_draw == 2 ){ sutable_rot_mode = 1; }
		else if( ndim_field == 3 ){ sutable_rot_mode = 3; }
		else{ sutable_rot_mode = 2; }

		const unsigned int nedge = m_nline;
		if( isnt_value_disp ){
			for(unsigned int iedge=0;iedge<nedge;iedge++){
				const unsigned int ipoin_va0 = m_EdgeAry[iedge*2  ];	// 点０
				const unsigned int ipoin_va1 = m_EdgeAry[iedge*2+1];	// 点１
				double* pval0 = &m_paVer->pVertexArray[iedge*2*ndim_draw];
				double* pval1 = &m_paVer->pVertexArray[iedge*2*ndim_draw+ndim_draw];
				const unsigned int ipoin_co0 = field.GetMapVal2Co(ipoin_va0);
				const unsigned int ipoin_co1 = field.GetMapVal2Co(ipoin_va1);
				ns_c_co.GetValue(ipoin_co0,pval0);
				ns_c_co.GetValue(ipoin_co1,pval1);
			}
		}
		else{	// 変位を伴う場合
			const unsigned int id_na_c_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_c_val) );
			const Fem::Field::CNodeAry& na_c_va = world.GetNA(id_na_c_val);
			const unsigned int id_ns_c_val = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			assert( na_c_va.IsSegID(id_ns_c_val) );
			const Fem::Field::CNodeAry::CNodeSeg& ns_c_va = na_c_va.GetSeg(id_ns_c_val);
			if( ndim_field == 2 && ndim_draw == 3 ){
				double coord[3], value[3];
				for(unsigned int iedge=0;iedge<nedge;iedge++){
					const unsigned int ipoin_va0 = m_EdgeAry[iedge*2  ];	// 点０
					const unsigned int ipoin_co0 = field.GetMapVal2Co(ipoin_va0);
					ns_c_co.GetValue(ipoin_co0,coord);
					ns_c_va.GetValue(ipoin_va0,value);
//					double val0=0.0, val1=0.0;
					m_paVer->pVertexArray[iedge*2*ndim_draw+0] = coord[0];
					m_paVer->pVertexArray[iedge*2*ndim_draw+1] = coord[1];
					m_paVer->pVertexArray[iedge*2*ndim_draw+2] = value[0];
					const unsigned int ipoin_va1 = m_EdgeAry[iedge*2+1];	// 点１
					const unsigned int ipoin_co1 = field.GetMapVal2Co(ipoin_va1);
					ns_c_co.GetValue(ipoin_co1,coord);
					ns_c_va.GetValue(ipoin_va1,value);
					m_paVer->pVertexArray[iedge*2*ndim_draw+ndim_draw+0] = coord[0];
					m_paVer->pVertexArray[iedge*2*ndim_draw+ndim_draw+1] = coord[1];
					m_paVer->pVertexArray[iedge*2*ndim_draw+ndim_draw+2] = value[0];
				}
			}
			else{
				assert( ndim_field == ns_c_va.Length() );
				double coord[3], value[3];
				for(unsigned int iedge=0;iedge<nedge;iedge++){
					const unsigned int ipoin_va0 = m_EdgeAry[iedge*2  ];	// 点０
					const unsigned int ipoin_co0 = field.GetMapVal2Co(ipoin_va0);
					ns_c_co.GetValue(ipoin_co0,coord);
					ns_c_va.GetValue(ipoin_va0,value);
//					double val0=0.0, val1=0.0;
					for(unsigned int idim=0;idim<ndim_draw;idim++){	// 点０の座標をセット
						m_paVer->pVertexArray[iedge*2*ndim_draw+idim] = coord[idim]+value[idim];
					}
					const unsigned int ipoin_va1 = m_EdgeAry[iedge*2+1];	// 点１
					const unsigned int ipoin_co1 = field.GetMapVal2Co(ipoin_va1);
					ns_c_co.GetValue(ipoin_co1,coord);
					ns_c_va.GetValue(ipoin_va1,value);
					for(unsigned int idim=0;idim<ndim_draw;idim++){ // 点１の座標をセット
						m_paVer->pVertexArray[iedge*2*ndim_draw+ndim_draw+idim] = coord[idim]+value[idim];
					}
				}
			}
		}
	}
	return true;
}

bool CDrawerEdge::Set(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	this->m_IdField = id_field;
	this->isnt_value_disp = isnt_value_disp;
	if( m_paVer != 0 ){ delete m_paVer; m_paVer=0; }
	const Fem::Field::CField& field = world.GetField(id_field);

	unsigned int nedge = 0;
	m_EdgeAry.resize(0);
	{	// 辺をゲットする関数
		std::vector<unsigned int> edge_ary_tmp;
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
			unsigned int id_es_c_va = field.GetIdElemSeg(id_ea,CORNER,true, world);
			unsigned int id_es_c_co = field.GetIdElemSeg(id_ea,CORNER,false,world);
			unsigned int nedge_tmp = 0;
			edge_ary_tmp.resize(0);
			Fem::Field::ELEM_TYPE elemtype = ea.ElemType();
//			std::cout << elemtype << std::endl;
			if( elemtype == TRI || elemtype == QUAD ){
				if(      ea.IsSegID(id_es_c_va) ){
					ea.MakeEdge(id_es_c_va,nedge_tmp,edge_ary_tmp);
				}
				else if( ea.IsSegID(id_es_c_co) ){
					ea.MakeEdge(id_es_c_co,nedge_tmp,edge_ary_tmp);
				}
			}
			else if( elemtype == TET || elemtype == HEX ){	// 境界要素を作ってから辺を抽出
				unsigned int id_es_add;
				std::vector<unsigned int> aIndElemFace;
				assert( ea.IsSegID(id_es_c_co) );
				CElemAry* pEA = ea.MakeBoundElemAry(id_es_c_co,id_es_add,aIndElemFace);	
				assert( pEA != 0 );
				assert( pEA->IsSegID(id_es_add) );
				pEA->MakeEdge(id_es_add,nedge_tmp,edge_ary_tmp);
				delete pEA;
			}
			else{
//				std::cout << "Not Supporting ElemType : " << elemtype << std::endl;
//				assert(0);
			}
//			std::cout << nedge_tmp << std::endl;
			nedge += nedge_tmp;
			m_EdgeAry.reserve( m_EdgeAry.size()+edge_ary_tmp.size() );
			for(unsigned int i=0;i<edge_ary_tmp.size();i++){
				m_EdgeAry.push_back( edge_ary_tmp[i] );
			}
		}
//		std::cout << nedge << std::endl;
	}
	m_nline = nedge;

	this->Update(world);
	return true;
}

void CDrawerEdge::Draw() const
{
	if( m_nline == 0 ) return;
	assert( m_paVer != 0 );
  const bool is_texture = glIsEnabled(GL_TEXTURE_2D); ::glDisable(GL_TEXTURE_2D);
  const bool is_lighting = glIsEnabled(GL_LIGHTING);  ::glDisable(GL_LIGHTING);
  
  ::glColor3d(0.0,0.0,0.0);
	::glLineWidth(line_width_);
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(m_paVer->NDim(),GL_DOUBLE,0,m_paVer->pVertexArray);
   if( m_paVer->NDim() == 2 ){ ::glTranslated(0,0,+0.01); }
	::glDrawArrays(GL_LINES,0,m_nline*2);
   if( m_paVer->NDim() == 2 ){ ::glTranslated(0,0,-0.01); }
	::glDisableClientState(GL_VERTEX_ARRAY);
  if( is_texture ){ ::glEnable(GL_TEXTURE_2D); }
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }
}
