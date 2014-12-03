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
// drawer_field_image_based_flow_vis.cpp : 場可視化クラス(DrawerField)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
#endif

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <stdio.h>

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "delfem/drawer_field_image_based_flow_vis.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

// 三角形の面積を求める関数
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

bool CEdgeTextureColor::Update(const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(this->id_field_velo) ){ return false; }
	if( !world.IsIdField(this->id_part_field_velo) ){ return false; }
	if( !world.IsIdEA(this->id_ea) ){ return false; }
	const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == Fem::Field::LINE );
	if( nelem != ea.Size() ){
		nelem  = ea.Size();
		if( m_aXYVeloElem != 0 ){ delete[] m_aXYVeloElem; }
		if( m_aXYElem     != 0 ){ delete[] m_aXYElem; }
		m_aXYVeloElem = new double [nelem*4];
		m_aXYElem     = new double [nelem*4];
	}
	const Fem::Field::CField& fv = world.GetField(this->id_part_field_velo);
//	const Fem::Field::CField::CNodeSegInNodeAry& nans = fv.GetNodeSegInNodeAry(Fem::Field::CORNER);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea,CORNER,true, world);
	const Fem::Field::CElemAry::CElemSeg& es_c = fv.GetElemSeg(id_ea,CORNER,false,world);
	assert( es_v.Length() == 2 );
	assert( es_c.Length() == 2 );
	for(unsigned int ielem=0;ielem<nelem;ielem++)
	{
		unsigned int noes[2];
		es_c.GetNodes(ielem,noes);
		double co[2][2];
		ns_c.GetValue(noes[0],co[0]);
		ns_c.GetValue(noes[1],co[1]);
		////////////////
		es_v.GetNodes(ielem,noes);
		double velo[2][2];
		ns_v.GetValue(noes[0],velo[0]);	
		ns_v.GetValue(noes[1],velo[1]);
		////////////////
		m_aXYElem[ielem*4+0] = co[0][0];
		m_aXYElem[ielem*4+1] = co[0][1];
		m_aXYElem[ielem*4+2] = co[1][0];
		m_aXYElem[ielem*4+3] = co[1][1];
		////////////////
		const double r = velo_scale;
		m_aXYVeloElem[ielem*4+0] = co[0][0] + r*velo[0][0];	
		m_aXYVeloElem[ielem*4+1] = co[0][1] + r*velo[0][1];
		m_aXYVeloElem[ielem*4+2] = co[1][0] + r*velo[1][0];
		m_aXYVeloElem[ielem*4+3] = co[1][1] + r*velo[1][1];
	}
	return true;
}

void CEdgeTextureColor::Draw() const
{
	::glColor3d(color[0], color[1], color[2]);
	::glBegin(GL_QUADS);
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		::glVertex2dv(m_aXYElem+ielem*4+2);
		::glVertex2dv(m_aXYElem+ielem*4+0);
		::glVertex2dv(m_aXYVeloElem+ielem*4+0);
		::glVertex2dv(m_aXYVeloElem+ielem*4+2);
	}
	::glEnd();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


CDrawerImageBasedFlowVis::CDrawerImageBasedFlowVis()
{
	m_IdFieldVelo = 0;
	m_IdFieldColor = 0;
	////////////////
	nnode = 0;
	aVelo = 0;
	aCoord = 0;
	aValColor = 0;
	////////////////
	iPtn = 0;
    alpha = 0;
	imode = 1;
    ////////////////
    velo_scale = 0.1;
	this->MakePattern();
}

CDrawerImageBasedFlowVis::CDrawerImageBasedFlowVis(const unsigned int id_field, 
												   const Fem::Field::CFieldWorld& world,
												   unsigned int imode)
{
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);
	////////////////
	nnode = 0;
	aVelo = 0;
	aCoord = 0;
	aValColor = 0;
	////////////////
    this->imode = imode;
    if( imode == 0 ){ alpha = 0; }
	else{ alpha = 40; }
	iPtn = 0;
	m_nameDisplayList = 0;
    ////////////////
    velo_scale = 0.02;
    ////////////////
	this->MakePattern();
	m_IdFieldVelo = 0;
	m_IdFieldColor = 0;
	this->Set( id_field, world );
}

void CDrawerImageBasedFlowVis::ClearDisplayList()
{
	if( m_nameDisplayList == 0 ) return;
	::glDeleteLists(m_nameDisplayList, m_nPattern);
	m_nameDisplayList = 0;
}

void CDrawerImageBasedFlowVis::MakePattern() 
{ 
	this->ClearDisplayList();
	const unsigned int NPN = 32;	// テクスチャのサイズ
	if( imode == 0 ){
		m_nPattern = 0;
		m_nameDisplayList = 0;
	}
	else if( imode == 1 )
	{
		m_nPattern = 32;
		m_nameDisplayList = glGenLists(m_nPattern);
		int lut[256];
//		for(unsigned int i=0; i<256; i++){ lut[i] =  (i < 127) ? 0 : 255; }
		for(unsigned int i=0; i<256; i++){ lut[i] =  (i < 230) ? 0 : 150; }
		int phase[NPN][NPN];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){ 
			phase[i][j] = rand() % 256; 
		}
		}
		GLubyte pat[NPN][NPN][4];
		for(unsigned int ipat=0; ipat<m_nPattern; ipat++){
			const unsigned int t = ipat*256/m_nPattern;
			for(unsigned int i=0; i<NPN; i++){
			for(unsigned int j=0; j<NPN; j++){
				GLubyte val = lut[(t + phase[i][j]) % 255];
				pat[i][j][0] = val;
				pat[i][j][1] = val;
				pat[i][j][2] = val;
				pat[i][j][3] = alpha;
			}
			}
			glNewList(m_nameDisplayList + ipat, GL_COMPILE);
			glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
			glEndList();
		}
	}
	else if( imode == 2 )
	{
		m_nPattern = 1;
		m_nameDisplayList = glGenLists(m_nPattern);
		GLubyte pat[NPN][NPN][4];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){
			const int phase = ( abs((int)(i % 16)) < 3 || abs((int)(j % 16)) < 3 ) ?  0 : 255;
			pat[i][j][0] = phase;
			pat[i][j][1] = phase;
			pat[i][j][2] = phase;
            pat[i][j][3] = alpha;
		}
		}
		glNewList(m_nameDisplayList, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();
	}
	else if( imode == 3 )
	{
		m_nPattern = 1;
		m_nameDisplayList = glGenLists(m_nPattern);
		GLubyte pat[NPN][NPN][4];
		for(unsigned int i=0; i<NPN; i++){
		for(unsigned int j=0; j<NPN; j++){			
			const int phase = ( abs((int)(i % 16)) < 3 && abs((int)(j % 16)) < 3 ) ?  0 : 255;
			pat[i][j][0] = phase;
			pat[i][j][1] = phase;
			pat[i][j][2] = phase;
			pat[i][j][3] = alpha;
		}
		}
		glNewList(m_nameDisplayList, GL_COMPILE);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0,   GL_RGBA, GL_UNSIGNED_BYTE, pat);
		glEndList();		
   }
}

void CDrawerImageBasedFlowVis::SetColorField(unsigned int id_field_color, const CFieldWorld& world,
											 std::auto_ptr<CColorMap> color_map)
{
	m_IdFieldColor = id_field_color;
//	if( aColorElem == 0 ){ delete[] aColorElem; }
	if( !world.IsIdField(this->m_IdFieldVelo) ) return;
	const Fem::Field::CField& fv = world.GetField(m_IdFieldVelo);
	const Fem::Field::CField& fc = world.GetField(id_field_color);
	const std::vector<unsigned int>& aIdEA_v = fv.GetAryIdEA();
	const std::vector<unsigned int>& aIdEA_c = fc.GetAryIdEA();
	assert( aIdEA_v.size() == aIdEA_c.size() );
	if( aIdEA_v.size() != aIdEA_c.size() ) return;
	for(unsigned int iiea=0;iiea<aIdEA_v.size();iiea++){
		unsigned int id_ea_v = aIdEA_v[iiea];
		unsigned int id_ea_c = aIdEA_c[iiea];
		assert( id_ea_v == id_ea_c );
		unsigned int id_es_v = fv.GetIdElemSeg(id_ea_v,CORNER,true,world);
		unsigned int id_es_c = fc.GetIdElemSeg(id_ea_c,CORNER,true,world);
		assert( id_es_v == id_es_c );
	}
//	aColorElem = new double [nelem*3];
	this->color_map = color_map;
	this->Update(world);
}

bool CDrawerImageBasedFlowVis::Set(unsigned int id_field_velo, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field_velo) ){ return false; }
	this->m_IdFieldVelo = id_field_velo;
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
	
	////////////////////////////////
	{	// 要素配列の設定
		const std::vector<unsigned int>& aIdEA = fv.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const unsigned int id_es_c = fv.GetIdElemSeg(id_ea,CORNER,true, world);
			assert( id_es_c != 0 );
			CIndexArrayElem* pIAE = new CIndexArrayElem(id_ea,id_es_c,world);
			pIAE->ilayer = fv.GetLayer(id_ea);
			this->m_apIndexArrayElem.push_back( pIAE );
		}
	}
	
//	if( aVelo ==0 ){ aVelo  = new double [nnode*2]; }
//	if( aCoord == 0 ){ aCoord = new double [nnode*2]; }
/*
	const std::vector<unsigned int>& aIdEA = fv.GetAryIdEA();
	unsigned int nelem0 = 0;
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		assert( world.IsIdEA(id_ea) );
		const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
		nelem0 += ea.Size();
	}
	if( nelem != nelem0 ){
		nelem = nelem0;
		if( aXYVeloElem != 0 ){ delete[] aXYVeloElem; }
        if( aXYElem     != 0 ){ delete[] aXYElem; }
		aXYVeloElem = new double [nelem*6];
		aXYElem = new double [nelem*6];
		if( world.IsIdField(m_IdFieldColor) ){
			if( aColorElem != 0 ){ delete[] aColorElem; }
			aColorElem = new double [nelem*3];
		}
	}
	*/
	return this->Update(world);
}
		
bool CDrawerImageBasedFlowVis::Update(const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(this->m_IdFieldVelo) ){ return false; }
	const Fem::Field::CField& fv = world.GetField(this->m_IdFieldVelo);
	////////////////////////////////
	{	// 要素配列の設定
		const std::vector<unsigned int>& aIdEA = fv.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const unsigned int id_es_v_c = fv.GetIdElemSeg(id_ea,CORNER,true, world);
			assert( id_es_v_c != 0 );
			CIndexArrayElem* pIAE = new CIndexArrayElem(id_ea,id_es_v_c,world);
			pIAE->ilayer = fv.GetLayer(id_ea);
			this->m_apIndexArrayElem.push_back( pIAE );
		}
	}
/*	const std::vector<unsigned int>& aIdEA = fv.GetAryIdEA();
	{
		unsigned int nelem0 = 0;
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			assert( world.IsIdEA(id_ea) );
			const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
			nelem0 += ea.Size();
		}
		if( nelem != nelem0 ){
			nelem = nelem0;
			if( aXYVeloElem != 0 ){ delete[] aXYVeloElem; }
			if( aXYElem     != 0 ){ delete[] aXYElem; }
			aXYVeloElem = new double [nelem*6];
			aXYElem     = new double [nelem*6];
			if( world.IsIdField(m_IdFieldColor) ){
				if( aColorElem != 0 ){ delete[] aColorElem; }
				aColorElem = new double [nelem*3];
			}
		}
	}*/
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	nnode = ns_v.Size();
	if( aVelo ==0 ){ aVelo  = new double [nnode*2]; }
	if( aCoord == 0 ){ aCoord = new double [nnode*2]; }
	for(unsigned int inode=0;inode<nnode;inode++){
		ns_v.GetValue(inode,aVelo+inode*2);
		unsigned int inoc = fv.GetMapVal2Co(inode);
		ns_c.GetValue(inoc,aCoord+inode*2);
	}
	/*
	unsigned int ielem_cur = 0;
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
	{
		const unsigned int id_ea = aIdEA[iiea];
		assert( world.IsIdEA(id_ea) );
		const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea,CORNER,true,world);
		const Fem::Field::CElemAry::CElemSeg& es_c = fv.GetElemSeg(id_ea,CORNER,false,world);
		for(unsigned int ielem=0;ielem<es_v.GetSizeElem();ielem++){
			unsigned int noes[3];
			es_c.GetNodes(ielem,noes);
			double co[3][2];
			ns_c.GetValue(noes[0],co[0]);
			ns_c.GetValue(noes[1],co[1]);
			ns_c.GetValue(noes[2],co[2]);
			////////////////
			es_v.GetNodes(ielem,noes);
			double velo[3][2];
			ns_v.GetValue(noes[0],velo[0]);
			ns_v.GetValue(noes[1],velo[1]);
			ns_v.GetValue(noes[2],velo[2]);
			////////////////
			aXYElem[ielem_cur*6+0] = co[0][0];
			aXYElem[ielem_cur*6+1] = co[0][1];
			aXYElem[ielem_cur*6+2] = co[1][0];
			aXYElem[ielem_cur*6+3] = co[1][1];
			aXYElem[ielem_cur*6+4] = co[2][0];
			aXYElem[ielem_cur*6+5] = co[2][1];
			////////////////
            const double r = velo_scale;
			aXYVeloElem[ielem_cur*6+0] = co[0][0] + r*velo[0][0];
			aXYVeloElem[ielem_cur*6+1] = co[0][1] + r*velo[0][1];
			aXYVeloElem[ielem_cur*6+2] = co[1][0] + r*velo[1][0];
			aXYVeloElem[ielem_cur*6+3] = co[1][1] + r*velo[1][1];
			aXYVeloElem[ielem_cur*6+4] = co[2][0] + r*velo[2][0];
			aXYVeloElem[ielem_cur*6+5] = co[2][1] + r*velo[2][1];
			ielem_cur++;
		}
	}
    assert( ielem_cur == nelem );
	*/
	if( world.IsIdField(m_IdFieldColor) ){
		if( aValColor==0 ){ aValColor = new double [nnode]; }
		const Fem::Field::CField& fc = world.GetField(this->m_IdFieldColor);
		assert( fc.GetNodeSegInNodeAry(CORNER).id_na_va == fv.GetNodeSegInNodeAry(CORNER).id_na_va );
		const Fem::Field::CNodeAry::CNodeSeg& nsv_c = fc.GetNodeSeg(CORNER,true, world,VALUE);
		for(unsigned int inode=0;inode<nnode;inode++){
			nsv_c.GetValue(inode,aValColor+inode);
		}
	}
	else{
		if( aValColor != 0 ){ delete[] aValColor; aValColor = 0; }
	}
	for(unsigned int iec=0;iec<aEdgeColor.size();iec++){
		aEdgeColor[iec].Update(world);
	}
	return true;
}

Com::CBoundingBox3D CDrawerImageBasedFlowVis::GetBoundingBox( double rot[] ) const
{
	if( nnode == 0 ){ return Com::CBoundingBox3D(-1,1, -1,1, 0,0); }
	double x_min,x_max, y_min,y_max;
	x_min = x_max = aCoord[0];
	y_min = y_max = aCoord[1];
	for(unsigned int inode=0;inode<nnode*3;inode++){
		const double x0 = aCoord[inode*2+0];
		const double y0 = aCoord[inode*2+1];
		x_min = ( x0 < x_min ) ? x0 : x_min;
		x_max = ( x0 > x_max ) ? x0 : x_max;
		y_min = ( y0 < y_min ) ? y0 : y_min;
		y_max = ( y0 > y_max ) ? y0 : y_max;
	}
	return Com::CBoundingBox3D(x_min,x_max, y_min,y_max, -1,1);
}

void CDrawerImageBasedFlowVis::Draw() const
{
	for(unsigned int iec=0;iec<aEdgeColor.size();iec++){ aEdgeColor[iec].Draw(); }

	////////////////////////////////
	// ModelView-Projection行列を取得
	GLint view[4];
	::glGetIntegerv(GL_VIEWPORT,view);
	const int win_w = view[2];
	const int win_h = view[3];
	const double invww = 1.0/win_w;
	const double invwh = 1.0/win_h;
	GLdouble model[16], prj[16];
	::glGetDoublev(GL_MODELVIEW_MATRIX,  model);	// ModelVeiw行列の取得
	::glGetDoublev(GL_PROJECTION_MATRIX, prj);		// Projection行列の取得
	////////////////
	::glEnable(GL_TEXTURE_2D);
	::glDisable(GL_DEPTH_TEST);
	if( imode != 0 ){	// 前画像のテクスチャを流す
		glEnable(GL_BLEND); 
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		::glColor4f(1.0f, 1.0f, 1.0f, 0.2f);
/*		::glDisable(GL_BLEND);
		::glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		::glBegin(GL_TRIANGLES);*/
		// 要素配列の設定
		for(unsigned int iIAE=0;iIAE<m_apIndexArrayElem.size();iIAE++){
			unsigned int nelem = m_apIndexArrayElem[iIAE]->GetSize();
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				unsigned int  noes[3];
				m_apIndexArrayElem[iIAE]->GetNoes(ielem,noes);
				const double co[3][2] = {
					{ aCoord[noes[0]*2+0], aCoord[noes[0]*2+1] },
					{ aCoord[noes[1]*2+0], aCoord[noes[1]*2+1] },
					{ aCoord[noes[2]*2+0], aCoord[noes[2]*2+1] } };
				double co_win[3][2], dtmp1;
				::gluProject(co[0][0],co[0][1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
				::gluProject(co[1][0],co[1][1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
				::gluProject(co[2][0],co[2][1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);				
				const double velo[3][2] = {
					{ aVelo[noes[0]*2+0], aVelo[noes[0]*2+1] },
					{ aVelo[noes[1]*2+0], aVelo[noes[1]*2+1] },
					{ aVelo[noes[2]*2+0], aVelo[noes[2]*2+1] } };
				const double vs = this->velo_scale*10000;
				::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(co[0][0]+vs*velo[0][0], co[0][1]+vs*velo[0][1]);
				::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(co[1][0]+vs*velo[1][0], co[1][1]+vs*velo[1][1]);
				::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(co[2][0]+vs*velo[2][0], co[2][1]+vs*velo[2][1]);
			}
		}
		::glEnd();
			/*
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			double co_obj[3][2];
			aXYElem[ielem*6];
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			const double* velo_co = &aXYVeloElem[ielem*6];
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2dv(velo_co+0);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2dv(velo_co+2);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2dv(velo_co+4);
		}
		*/
	}
	/*
	else{ // imode == 0 : 純粋移流拡散
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		::glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			const double* co_obj = &aXYElem[ielem*6];
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			const double* velo_co = &aXYVeloElem[ielem*6];
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(velo_co[0], velo_co[1]);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(velo_co[2], velo_co[3]);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(velo_co[4], velo_co[5]);
		}
		::glEnd();
		////////////////
		::glDisable(GL_DEPTH_TEST);
		::glEnable(GL_BLEND);		
		::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		////////////////
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		::glColor4f(1.0f, 1.0f, 1.0f, 0.8f);
		::glBegin(GL_TRIANGLES);
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			const double* co_obj = &aXYElem[ielem*6];
			const double area = TriArea2D(co_obj,co_obj+2,co_obj+4);
			double co_win[3][2], dtmp1;
			::gluProject(co_obj[0*2+0],co_obj[0*2+1],0, model,prj,view, &co_win[0][0],&co_win[0][1],&dtmp1);
			::gluProject(co_obj[1*2+0],co_obj[1*2+1],0, model,prj,view, &co_win[1][0],&co_win[1][1],&dtmp1);
			::gluProject(co_obj[2*2+0],co_obj[2*2+1],0, model,prj,view, &co_win[2][0],&co_win[2][1],&dtmp1);
			::glTexCoord2d(co_win[0][0]*invww, co_win[0][1]*invwh);  ::glVertex2d(co_obj[0], co_obj[1]);
			::glTexCoord2d(co_win[1][0]*invww, co_win[1][1]*invwh);  ::glVertex2d(co_obj[2], co_obj[3]);
			::glTexCoord2d(co_win[2][0]*invww, co_win[2][1]*invwh);  ::glVertex2d(co_obj[4], co_obj[5]);
		}
		::glEnd();
		::glDisable(GL_BLEND);
	}
	*/
	iPtn = iPtn + 1;

	if( m_nPattern != 0 ){
		::glShadeModel(GL_SMOOTH);
		::glDisable(GL_DEPTH_TEST);
		::glEnable(GL_BLEND); 
//		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
//		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
//		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
		::glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_ADD);
		::glCallList(iPtn % m_nPattern + m_nameDisplayList);	//　パターンを呼ぶ
		if( aValColor == 0 ){
			::glBegin(GL_TRIANGLES);
			for(unsigned int iIAE=0;iIAE<m_apIndexArrayElem.size();iIAE++){
				const unsigned int nelem = m_apIndexArrayElem[iIAE]->GetSize();
				for(unsigned int ielem=0;ielem<nelem;ielem++){
					unsigned  int noes[3];
					m_apIndexArrayElem[iIAE]->GetNoes(ielem,noes);					
					double co[3][2] = {
						{ aCoord[noes[0]*2+0], aCoord[noes[0]*2+1] },
						{ aCoord[noes[1]*2+0], aCoord[noes[1]*2+1] },
						{ aCoord[noes[2]*2+0], aCoord[noes[2]*2+1] } };
					double r = 10.0;	// 流す画像の周期をきめる(大きければ細かい)	
					::glTexCoord2d(r*co[0][0],r*co[0][1]);		::glVertex2dv(co[0]);
					::glTexCoord2d(r*co[1][0],r*co[1][1]);		::glVertex2dv(co[1]);
					::glTexCoord2d(r*co[2][0],r*co[2][1]);		::glVertex2dv(co[2]);
				}
			}
			::glEnd();
			/*
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				const double co[3][2] = {
					{ aXYElem[ielem*6+0], aXYElem[ielem*6+1] },
					{ aXYElem[ielem*6+2], aXYElem[ielem*6+3] },
					{ aXYElem[ielem*6+4], aXYElem[ielem*6+5] } };
				double r = 10.0;	// 流す画像の周期をきめる(大きければ細かい)	
				::glTexCoord2d(r*co[0][0],r*co[0][1]);		::glVertex2dv(co[0]);
				::glTexCoord2d(r*co[1][0],r*co[1][1]);		::glVertex2dv(co[1]);
				::glTexCoord2d(r*co[2][0],r*co[2][1]);		::glVertex2dv(co[2]);
			}
			*/
		}
		else{
			::glBegin(GL_TRIANGLES);
			for(unsigned int iIAE=0;iIAE<m_apIndexArrayElem.size();iIAE++){
				const unsigned int nelem = m_apIndexArrayElem[iIAE]->GetSize();
				for(unsigned int ielem=0;ielem<nelem;ielem++){
					unsigned int  noes[3];
					m_apIndexArrayElem[iIAE]->GetNoes(ielem,noes);				
					double co[3][2] = {
						{ aCoord[noes[0]*2+0], aCoord[noes[0]*2+1] },
						{ aCoord[noes[1]*2+0], aCoord[noes[1]*2+1] },
						{ aCoord[noes[2]*2+0], aCoord[noes[2]*2+1] } };
					double r = 10.0;	// 流す画像の周期をきめる(大きければ細かい)	
					for(unsigned int ino=0;ino<3;ino++){
						const double d = aValColor[noes[ino]];
						float color[3];
						color_map->GetColor(color,d);
						::glColor4d(color[0],color[1],color[2],1.0);
						::glTexCoord2d(r*co[ino][0],r*co[ino][1]);		::glVertex2dv(co[ino]);
					}
				}
			}
			::glEnd();
			/*
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				const double co[3][2] = {
					{ aXYElem[ielem*6+0], aXYElem[ielem*6+1] },
					{ aXYElem[ielem*6+2], aXYElem[ielem*6+3] },
					{ aXYElem[ielem*6+4], aXYElem[ielem*6+5] } };
				double r = 10.0;	// 流す画像の周期をきめる(大きければ細かい)	
	//			std::cout<< ielem << " " << aColorElem[ielem*3+0] << " " << aColorElem[ielem*3+1] << " " << aColorElem[ielem*3+2] << std::endl;
				for(unsigned int ino=0;ino<3;ino++){
					const double d = aColorElem[ielem*3+ino];
					float color[3];
					color_map->GetColor(color,d);
					::glColor4d(color[0],color[1],color[2],1.0);
					::glTexCoord2d(r*co[ino][0],r*co[ino][1]);		::glVertex2dv(co[ino]);
				}
			}
			::glEnd();
			*/
		}
		::glDisable(GL_BLEND);
	}
	::glEnable(GL_DEPTH_TEST);
	
	////////////////
	::glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, win_w, win_h, 0);	// テクスチャを作る
	::glDisable(GL_BLEND);
}
