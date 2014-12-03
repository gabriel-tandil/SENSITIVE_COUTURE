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
// DrawerField.cpp : implementation of the field visualization class (DrawerField)
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
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "delfem/drawer_field_face.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

CDrawerFace::CDrawerFace()
{
  pNormalArray = 0;  
  pUVArray = 0;
  tex_cent_x = 0;
  tex_cent_y = 0;
  tex_scale = 1;
  ////
	pColorArray = 0;
	is_draw_color_legend = false;
	color_map = std::auto_ptr<CColorMap>(new CColorMap());
}

CDrawerFace::CDrawerFace
(const unsigned int id_field, bool isnt_value_disp, 
 const Fem::Field::CFieldWorld& world, unsigned int id_field_color)
{
  pNormalArray = 0;
  pUVArray = 0;  
  tex_cent_x = 0;
  tex_cent_y = 0;  
  tex_scale = 1;  
  ////
	pColorArray = 0;
	if( world.IsIdField(id_field_color) ){ is_draw_color_legend = true;  }
	else{                                  is_draw_color_legend = false; }
	color_map = std::auto_ptr<CColorMap>(new CColorMap());
	this->Set( id_field, world, isnt_value_disp, id_field_color);
}

CDrawerFace::CDrawerFace
(const unsigned int id_field, bool isnt_value_disp, 
 const Fem::Field::CFieldWorld& world, unsigned int id_field_color,
 double min, double max)
{
  pUVArray = 0;    
  tex_cent_x = 0;
  tex_cent_y = 0;  
  pNormalArray = 0;
  ////
	pColorArray = 0;
	if( world.IsIdField(id_field_color) ){ is_draw_color_legend = true;  }
	else{                                  is_draw_color_legend = false; }
	color_map = std::auto_ptr<CColorMap>(new CColorMap(min,max));
	this->Set( id_field, world, isnt_value_disp, id_field_color);
}

CDrawerFace::CDrawerFace
(const unsigned int id_field, bool isnt_value_disp, 
 const Fem::Field::CFieldWorld& world, unsigned int id_field_color, 
 std::auto_ptr<CColorMap> color_map )
{
  pUVArray = 0;    
  tex_cent_x = 0;
  tex_cent_y = 0;    
  pNormalArray = 0;
  ////
	pColorArray = 0;
	is_draw_color_legend = false;
	this->color_map = color_map;
	this->Set( id_field, world, isnt_value_disp, id_field_color);
}

CDrawerFace::~CDrawerFace()
{
	if( pColorArray  != 0 ){ delete[] pColorArray; }
  if( pNormalArray != 0 ){ delete[] pNormalArray; }
  if( pUVArray     != 0 ){ delete[] pUVArray; }
	for(unsigned int i=0;i<this->m_apIndexArrayElem.size();i++){
		delete this->m_apIndexArrayElem[i];
	}
}

void CDrawerFace::EnableNormal(bool is_lighting){
  if( (pNormalArray != 0 ) == is_lighting ){ return; }
  if( is_lighting ){
    const unsigned int nnode = this->m_vertex_ary.NPoin();
    pNormalArray = new double [nnode*3];
    MakeNormal();
  }
  else{
    delete[] pNormalArray;
    pNormalArray = 0;
  }
}

void CDrawerFace::EnableUVMap(bool is_uv_map, const Fem::Field::CFieldWorld& world)
{
  if( (pUVArray != 0 ) == is_uv_map ){ return; }
  if( is_uv_map ){
    const unsigned int nnode = this->m_vertex_ary.NPoin();
    const unsigned int ndim  = this->m_vertex_ary.NDim();
    delete[] pUVArray;
    pUVArray = new double [nnode*2];
    if( !world.IsIdField(m_id_field) ) return;
    const Fem::Field::CField& field = world.GetField(m_id_field);
//    unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
//    assert( world.IsIdNA(id_na_c_co) );
//    const Fem::Field::CNodeAry& na_c_co = world.GetNA(id_na_c_co);
//    assert( field.IsNodeSeg(CORNER,false,world) );
    const Fem::Field::CNodeAry::CNodeSeg& ns_c_co = field.GetNodeSeg(CORNER,false,world);
//    const unsigned int ndim = ns_c_co.Length();      
//    assert( ndim >= 2 );
    for(unsigned int ino=0;ino<ns_c_co.Size();ino++){
      double c[3]; ns_c_co.GetValue(ino,c);
      pUVArray[ino*2+0] = c[0]*tex_scale;
      pUVArray[ino*2+1] = c[1]*tex_scale;
    }
  }
  else{
    delete[] pUVArray;
    pUVArray = 0;
  }  
}

void CDrawerFace::SetTexScale(double scale, const Fem::Field::CFieldWorld& world){
  tex_scale = scale;  
  if( pUVArray != 0 ){
    if( !world.IsIdField(m_id_field) ) return;
    const Fem::Field::CField& field = world.GetField(m_id_field);
    // set the vertex array
    unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
    assert( world.IsIdNA(id_na_c_co) );
    const Fem::Field::CNodeAry& na_c_co = world.GetNA(id_na_c_co);
    assert( field.IsNodeSeg(CORNER,false,world) );
    const Fem::Field::CNodeAry::CNodeSeg& ns_c_co = field.GetNodeSeg(CORNER,false,world);      
    for(unsigned int ino=0;ino<ns_c_co.Size();ino++){
      double c[3]; ns_c_co.GetValue(ino,c);
      pUVArray[ino*2+0] = c[0]*tex_scale;
      pUVArray[ino*2+1] = c[1]*tex_scale;
    }
  }   
}


void CDrawerFace::Draw() const 
{
	if( m_vertex_ary.NDim() == 2 ){	// cannot see the opposite side
		::glEnable(GL_CULL_FACE);
		::glCullFace(GL_BACK);
	}
	else{ ::glDisable(GL_CULL_FACE); }
	
	int ilayer_min, ilayer_max;
	{
		if( m_apIndexArrayElem.size() > 0 ){
			ilayer_min = this->m_apIndexArrayElem[0]->ilayer;
			ilayer_max = ilayer_min;
		}
		else{
			ilayer_min=0;	ilayer_max=0;
		}
		for(unsigned int idp=1;idp<this->m_apIndexArrayElem.size();idp++){
			const int ilayer = this->m_apIndexArrayElem[idp]->ilayer;
			ilayer_min = (ilayer<ilayer_min) ? ilayer : ilayer_min;
			ilayer_max = (ilayer>ilayer_max) ? ilayer : ilayer_max;
		}
//		std::cout << ilayer_min << " " << ilayer_max << std::endl;
	}
	const double layer_height = 1.0/(ilayer_max-ilayer_min+1);

	if( this->pColorArray == 0 ){ // color is assigned to element face
		::glLineWidth(3);
		::glEnableClientState(GL_VERTEX_ARRAY);
		::glVertexPointer(m_vertex_ary.NDim(),GL_DOUBLE,0,m_vertex_ary.pVertexArray);
    if( this->pNormalArray != 0 && glIsEnabled(GL_LIGHTING) ){
      ::glEnableClientState(GL_NORMAL_ARRAY);
      ::glNormalPointer(GL_DOUBLE,0,pNormalArray);
      float shine[4] = {0,0,0,0};  
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
      ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.0);      
    }
    if( this->pUVArray != 0 && glIsEnabled( GL_TEXTURE_2D ) ){
      ::glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      ::glTexCoordPointer(2,GL_DOUBLE,0,pUVArray);
      ::glMatrixMode(GL_TEXTURE);
      ::glLoadIdentity();
      ::glTranslated(-tex_cent_x, -tex_cent_y, 0.0);      
      ::glMatrixMode(GL_MODELVIEW);
    }
		for(unsigned int idp=0;idp<this->m_apIndexArrayElem.size();idp++){ 
			View::CIndexArrayElem* pIndexArray = this->m_apIndexArrayElem[idp];
			if( pIndexArray->GetElemDim() == 1 ){	// draw line
				::glColor3d(0.0,0.0,0.0);
				::glLineWidth(3);
			}
			if( pIndexArray->GetElemDim() == 2 ||  pIndexArray->GetElemDim() == 3 ){	// draw face
        ::glColor3fv(m_apIndexArrayElem[idp]->color);
			}      
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, m_apIndexArrayElem[idp]->color);              
//      ::glColor3d(0.8,0.8,0.8);      
			this->m_apIndexArrayElem[idp]->DrawElements(); 
		}
		::glDisableClientState(GL_VERTEX_ARRAY);
    ::glDisableClientState(GL_NORMAL_ARRAY);
    ::glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	}
	else{ // color is assigned to vertex
    ::glShadeModel(GL_SMOOTH);
		::glEnableClientState(GL_VERTEX_ARRAY);
		::glVertexPointer(m_vertex_ary.NDim(),GL_DOUBLE,0,m_vertex_ary.pVertexArray);
    if( this->pNormalArray != 0 && glIsEnabled(GL_LIGHTING) ){
      ::glEnableClientState(GL_NORMAL_ARRAY);
      ::glNormalPointer(GL_DOUBLE,0,pNormalArray);
      ::glEnable(GL_COLOR_MATERIAL);
      ::glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
      float shine[4] = {1,1,1,1};
      ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
      ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
    }    
    if( this->pUVArray != 0 && glIsEnabled( GL_TEXTURE_2D ) ){
      ::glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      ::glTexCoordPointer(2,GL_DOUBLE,0,pUVArray);
    }    
		::glEnableClientState(GL_COLOR_ARRAY);
		::glColorPointer(4,GL_FLOAT,0,pColorArray);
		for(unsigned int idp=0;idp<this->m_apIndexArrayElem.size();idp++){ 
			const unsigned int ilayer = m_apIndexArrayElem[idp]->ilayer;
			const double height = (ilayer-ilayer_min)*layer_height;
      ::glTranslated(0,0,+height);
			this->m_apIndexArrayElem[idp]->DrawElements(); 
      ::glTranslated(0,0,-height);
		}
		::glDisableClientState(GL_COLOR_ARRAY);
		::glDisableClientState(GL_VERTEX_ARRAY);
    ::glDisableClientState(GL_NORMAL_ARRAY);
    ::glDisableClientState(GL_TEXTURE_COORD_ARRAY);    
	}

	if( this->is_draw_color_legend )	
	{ // draw legend
    View::DrawColorLegend(*color_map);
	}
}
bool CDrawerFace::Update
(const Fem::Field::CFieldWorld& world)
{
	const Fem::Field::CField& field = world.GetField(m_id_field);
	// set the vertex array
	unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
	assert( world.IsIdNA(id_na_c_co) );
	const Fem::Field::CNodeAry& na_c_co = world.GetNA(id_na_c_co);
  assert( field.IsNodeSeg(CORNER,false,world) );
	const Fem::Field::CNodeAry::CNodeSeg& ns_c_co = field.GetNodeSeg(CORNER,false,world);
	const unsigned int ndim = ns_c_co.Length();
	const unsigned int npoin_co = na_c_co.Size();
	unsigned int npoin;	
	{ // calc the number of point
		if( m_is_draw_nsv ){
			unsigned int id_na_c_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_c_val) );
			const Fem::Field::CNodeAry& na_c_val = world.GetNA(id_na_c_val);
			npoin = na_c_val.Size();
		}
		else{ npoin = npoin_co; }
	}
//	std::cout << m_vertex_ary.NPoin() << " " << npoin << std::endl;
	assert( m_vertex_ary.NPoin() == npoin );

	if( !m_isnt_value_disp ){	// including displacement
		unsigned int id_na_c_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		assert( world.IsIdNA(id_na_c_val) );
		const Fem::Field::CNodeAry& na_c_val = world.GetNA(id_na_c_val);
		const Fem::Field::CNodeAry::CNodeSeg& ns_c_val = field.GetNodeSeg(CORNER,true,world,VALUE|VELOCITY|ACCELERATION);
		if( ndim == 2 && (field.GetFieldType()==SCALAR||field.GetFieldType()==ZSCALAR) )	// 垂直方向の変位として捉える
		{
			assert( m_vertex_ary.NDim() == 3 );
			double coord[2], value[2];
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
				assert( ipoin_co < npoin_co );
				ns_c_val.GetValue(ipoin,value);
				ns_c_co.GetValue(ipoin_co,coord);
				this->m_vertex_ary.pVertexArray[ipoin*3+0] = coord[0];
				this->m_vertex_ary.pVertexArray[ipoin*3+1] = coord[1];
				this->m_vertex_ary.pVertexArray[ipoin*3+2] = value[0];
			}
		}
		else{
			assert( m_vertex_ary.NDim() == ndim );
			assert( ndim == ns_c_val.Length() ); // Coordの次元とValueの次元が合ってなければならない
			assert( npoin_co == na_c_co.Size() );
			assert( na_c_val.Size() == npoin );
			assert( ndim <= 3 );
			double coord[3], value[3];
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
				assert( ipoin_co < npoin_co );
				ns_c_val.GetValue(ipoin,value);
				ns_c_co.GetValue(ipoin_co,coord);
				for(unsigned int idim=0;idim<ndim;idim++){
					this->m_vertex_ary.pVertexArray[ipoin*ndim+idim] = coord[idim]+value[idim];
				}
			}
		}
	}
	else{		
		assert( m_vertex_ary.NDim() == ndim );
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
			assert( ipoin_co < npoin_co );
			double* pval = &this->m_vertex_ary.pVertexArray[ipoin*ndim];
			ns_c_co.GetValue(ipoin_co,pval);
		}
	}

	////////////////////////////////////////////////
	// make color
  
	if( world.IsIdField(id_field_val) )
	{
		const Fem::Field::CField& field_val = world.GetField(id_field_val);
		Fem::Field::FIELD_DERIVATION_TYPE fdt;
		{	// 描画する値の微分タイプ(できるだけVALUEを選択)
			unsigned int fdt_all = field_val.GetFieldDerivativeType();
			if(      fdt_all & VALUE        ){ fdt = VALUE; }
			else if( fdt_all & VELOCITY     ){ fdt = VELOCITY; }
			else if( fdt_all & ACCELERATION ){ fdt = ACCELERATION; }
			else{ assert(0); }
		}
		if( !color_map->IsMinMaxFix() ){	// 値の最大値最小値を求める
			double min_val, max_val;
			field_val.GetMinMaxValue(min_val,max_val,world,0,fdt);
			color_map->SetMinMax(min_val,max_val);
		}
		unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		unsigned int id_na_b_val = field_val.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		if(      world.IsIdNA(id_na_c_val) ){
			const CNodeAry::CNodeSeg& ns_v = field_val.GetNodeSeg(CORNER,true,world,fdt);
			double val[10];
			if( pColorArray == 0 ){ pColorArray = new float [npoin*4]; }
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				ns_v.GetValue(ipoin,val);
        color_map->GetColor(&pColorArray[ipoin*4],val[0]);
				pColorArray[ipoin*4+3] = 0.0f;
			}
		}
		else if( world.IsIdNA(id_na_b_val) ){
			unsigned int id_ns_v = 0;
			if(      fdt == VALUE        ){ id_ns_v = field_val.GetNodeSegInNodeAry(BUBBLE).id_ns_va; }
			else if( fdt == VELOCITY     ){ id_ns_v = field_val.GetNodeSegInNodeAry(BUBBLE).id_ns_ve; }
			else if( fdt == ACCELERATION ){ id_ns_v = field_val.GetNodeSegInNodeAry(BUBBLE).id_ns_ac; }
			for(unsigned int idp=0;idp<this->m_apIndexArrayElem.size();idp++){ 
				View::CIndexArrayElem* pIA = this->m_apIndexArrayElem[idp];
				unsigned int id_ea = pIA->GetIdEA();
				unsigned int id_es_v = field_val.GetIdElemSeg(id_ea,BUBBLE,true,world);
				pIA->SetColor(id_es_v,id_ns_v,world,color_map);
			}
		}
	}
  
  /////////////////////
  // make normal
  if( pNormalArray != 0 ){ MakeNormal(); }  
  
  
  /////////////////////
  // make uv map  
  if( pUVArray != 0 ){
    for(unsigned int ino=0;ino<ns_c_co.Size();ino++){
      double c[3]; ns_c_co.GetValue(ino,c);
      pUVArray[ino*2+0] = c[0]*tex_scale;
      pUVArray[ino*2+1] = c[1]*tex_scale;
    }
  }
	return true;
}

bool CDrawerFace::Set
(unsigned int id_field, const Fem::Field::CFieldWorld& world, bool isnt_value_disp,
 unsigned int id_field_val)
{
	if( !world.IsIdField(id_field) ){ return false; }
	////////////////
	this->m_id_field = id_field;
	this->id_field_val = id_field_val;
	this->m_isnt_value_disp = isnt_value_disp;

	const Fem::Field::CField& field = world.GetField(id_field);

	// setting of vertex array
	unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_na_c_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
	////////////////////////////////
	// decide whether draw ns of value or coord
	if( id_na_c_val == 0 ){ 
		this->m_is_draw_nsv = false;	// draw NS of Coord
		this->m_isnt_value_disp = true;	// don't include displacement
	}
	else{ this->m_is_draw_nsv = true;}
	////////////////
	assert( field.IsNodeSeg(CORNER,false,world,VALUE) );
	unsigned int ndim_field = field.GetNDimCoord();
	////////////////
	unsigned int npoin;
	if( m_is_draw_nsv ){
		assert( world.IsIdNA(id_na_c_val) );
		const Fem::Field::CNodeAry& na_c_val = world.GetNA(id_na_c_val);
		npoin = na_c_val.Size();
	}
	else{ 
		assert( id_na_c_co != 0 );
		const Fem::Field::CNodeAry& na_c_co = world.GetNA(id_na_c_co);
		npoin = na_c_co.Size(); 
	}
	////////////////
  // set size to vertex array
  
	unsigned int ndim_draw;
	if( this->m_isnt_value_disp == false 
     && ndim_field == 2 
     && (field.GetFieldType()==SCALAR||field.GetFieldType()==ZSCALAR) ){
		ndim_draw = 3;
	}
	else{ ndim_draw = ndim_field; }
  this->m_vertex_ary.SetSize(npoin,ndim_draw);
  
  { // normal
    const bool is_normal = ( pNormalArray != 0 );
    if( pNormalArray != 0 ){ delete pNormalArray; pNormalArray = 0; }
    if( is_normal ){ pNormalArray = new double [npoin*3]; }
  }
  { // uv map
    const bool is_uv = ( pUVArray != 0 );
    if( pUVArray != 0 ){ delete pUVArray; pUVArray = 0; }
    if( is_uv ){ pUVArray = new double [npoin*2]; }  
  }
    
	////////////////
	if(      ndim_draw  == 2 ){ sutable_rot_mode = 1; }
	else if( ndim_field == 3 ){ sutable_rot_mode = 3; }
	else                      { sutable_rot_mode = 2; }
	CDrawerFace::Update(world);

	////////////////////////////////
	{	// setting of element array        
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			unsigned int id_es_c;
			if( m_is_draw_nsv ){ id_es_c = field.GetIdElemSeg(id_ea,CORNER,true, world); }
			else{                id_es_c = field.GetIdElemSeg(id_ea,CORNER,false,world); }
			assert( id_es_c != 0 );
			CIndexArrayElem* pIAE = new CIndexArrayElem(id_ea,id_es_c,world);
			pIAE->ilayer = field.GetLayer(id_ea);
			this->m_apIndexArrayElem.push_back( pIAE );
		}
	}

	////////////////////////////////
	// color setting

	if( world.IsIdField(id_field_val) ){
		const Fem::Field::CField& field_val = world.GetField(id_field_val);
		unsigned int id_na_c_val = field_val.GetNodeSegInNodeAry(CORNER).id_na_va;
		unsigned int id_na_b_val = field_val.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		if(      world.IsIdNA(id_na_c_val) ){
			if( pColorArray != 0 ){ delete[] pColorArray;  pColorArray=0; }
		}
		else if( world.IsIdNA(id_na_b_val) ){
		}
	}
	this->Update(world);
	return true;
}



inline void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
	n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
	const double invlen = 0.5/a;
	n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}


void CDrawerFace::MakeNormal()
{  
  const unsigned int npoin = m_vertex_ary.NPoin();
  if( pNormalArray == 0 ){ pNormalArray = new double [npoin*3]; }
  for(unsigned int i=0;i<3*npoin;i++){ pNormalArray[i] = 0; }   
  for(unsigned int idp=0;idp<this->m_apIndexArrayElem.size();idp++){ 
    View::CIndexArrayElem* pIA = this->m_apIndexArrayElem[idp];
    const unsigned int nelem = pIA->GetSize();
    for(unsigned int ielem=0;ielem<nelem;ielem++){
      unsigned int no[8];
      pIA->GetNoes(ielem, no);
      const double* c0 = &this->m_vertex_ary.pVertexArray[no[0]*3];
      const double* c1 = &this->m_vertex_ary.pVertexArray[no[1]*3];
      const double* c2 = &this->m_vertex_ary.pVertexArray[no[2]*3];
      double n[3], area;
      UnitNormalAreaTri3D(n,area,c0,c1,c2);
      pNormalArray[no[0]*3+0] += n[0];
      pNormalArray[no[0]*3+1] += n[1];
      pNormalArray[no[0]*3+2] += n[2];
      pNormalArray[no[1]*3+0] += n[0];
      pNormalArray[no[1]*3+1] += n[1];
      pNormalArray[no[1]*3+2] += n[2];
      pNormalArray[no[2]*3+0] += n[0];
      pNormalArray[no[2]*3+1] += n[1];
      pNormalArray[no[2]*3+2] += n[2];        
    }
  }
  for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
    double* p = &pNormalArray[ipoin*3];
    const double invlen = 1.0/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    p[0] *= invlen;
    p[1] *= invlen;
    p[2] *= invlen;    
    //      std::cout << ipoin << "   " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }
}


