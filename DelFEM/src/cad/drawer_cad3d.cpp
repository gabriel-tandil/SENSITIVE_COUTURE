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
// DrawerCad.cpp : implementation of the class CDrawer_Cad3D which visualize class CCadObj3D
// Please make it clash-less because this is just visualization
// Hence, Don't put assertion in the program 
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#include <windows.h>
#if defined(__VISUALC__)
#pragma comment (lib, "winmm.lib")     /* link with Windows MultiMedia lib */
#pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */


#if defined(__APPLE__) && defined(__MACH__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <assert.h>
#include <iostream>

#include "delfem/drawer_cad3d.h"
#include "delfem/cad2d_interface.h"
#include "delfem/mesher2d.h"

using namespace Cad::View;
using namespace Com;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



void CDrawer_Cad3D::CDrawPart::DrawElements() const
{
//	if( !this->is_show ) return;    // セレクションに必要
	if( npoel == 2 ){ 
		::glDrawElements(GL_LINES    ,nelem*npoel,GL_UNSIGNED_INT,pIndexArray); 
		return;
	}
	else if( npoel == 3 ){ 
		::glDrawElements(GL_TRIANGLES,nelem*npoel,GL_UNSIGNED_INT,pIndexArray); 
		return;
	}
	// assert( npoel == 0 );
}

////////////////////////////////////////////////////////////////

class CLoopPxy2D : public Cad::ICad2D_Msh
{
public:
  CLoopPxy2D(const Cad::CCadObj3D& cad, unsigned int id_l) : cad_(cad), id_l_(id_l) {
    if( !cad_.IsElemID(Cad::LOOP, id_l_) ) return;
    const Cad::CLoop3D& l = cad_.GetLoop(id_l_);
    vec_n = l.normal;    
    vec_x = l.dirx;
    vec_o = l.org;
    vec_y = Cross(vec_n,vec_x);  
  }
	virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const{
    return cad_.GetIdVertex_Edge(id_v_s,id_v_e,id_e);
  }  
	//! idが使われているかどうかを調べる関数
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE type, unsigned int id) const {
    if( type == Cad::LOOP ){ return id == id_l_; }
    return cad_.IsElemID(type,id);
  }
	//! すべてのIDを配列にして返す関数
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE type) const{
    if( type == Cad::LOOP ){ return std::vector<unsigned int>(1,id_l_); }
    return cad_.GetAryElemID(type);
  }	
	// レイヤ関係の関数
  virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const { return 0; }
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const { layer_min=0; layer_max=0; }
  
	////////////////////////////////
	// member fucntion for loop
  
  //! ID:id_lのループの色を返す(本来このクラスは位相と幾何情報以外を持つべきではないかもしれないので暫定的)
  virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const { 
    if( id_l != id_l_ ) return false;
    color[0] = 0.8; color[1] = 0.8; color[2] = 0.8;    
    return true;
  }
	//! get area loop (ID:id_l)
	virtual double GetArea_Loop(unsigned int id_l) const{
    return 1;
  }
	//! ID:id_lのループを構成する頂点や辺をめぐるイテレータを返す関数
  virtual std::auto_ptr<Cad::IItrLoop> GetPtrItrLoop(unsigned int id_l) const{
    return std::auto_ptr<Cad::IItrLoop>( new Cad::CBRepSurface::CItrLoop(cad_.GetItrLoop(id_l)) );	// instance
	}
  
  //! ID:id_eの辺のメッシュ分割を得る(elen<=0ならできるだけ詳細にメッシュを切ろうとする)
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const{
    aCo.clear();
    return true;
  }  
	virtual Com::CVector2D GetVertexCoord(unsigned int id_v) const {
    const CVector3D& v = cad_.GetVertexCoord(id_v);
    const CVector3D& dir = v-vec_o;
    const double x = Dot(dir,vec_x);
    const double y = Dot(dir,vec_y);
    return Com::CVector2D(x,y);
  }
private:    
  const Cad::CCadObj3D& cad_;
  Com::CVector3D vec_o, vec_x, vec_y, vec_n;
  unsigned int id_l_;
};

bool CDrawer_Cad3D::UpdateCAD_TopologyGeometry(const Cad::CCadObj3D& cad_3d)
{
	this->sutable_rot_mode = 3;	
  //! 今までのDrawerPartの配列を一端バッファにコピーして，必要な物だけを戻す  
  for(unsigned int idp=0;idp<m_apIndexAry.size();idp++){ delete m_apIndexAry[idp]; }  
  m_apIndexAry.clear();
  m_aIndexVertex.clear();  
  
  std::vector<CVector3D> aVec3D;
  aVec3D.reserve(256);
  const std::vector<unsigned int>& aIdV = cad_3d.GetAryElemID(Cad::VERTEX);  
  for(unsigned int i=0;i<aIdV.size();i++){
    aVec3D.push_back( cad_3d.GetVertexCoord(aIdV[i]) );
  }
  
  {
    const std::vector<unsigned int>& aIdE = cad_3d.GetAryElemID(Cad::EDGE);
    for(unsigned int iie=0;iie<aIdE.size();iie++){
      unsigned int id_e = aIdE[iie];
      assert( cad_3d.IsElemID(Cad::EDGE,id_e) );
      unsigned int id_vs, id_ve;
      cad_3d.GetIdVertex_Edge(id_vs,id_ve,id_e);
      CDrawPart* pcp = new CDrawPart();
      CDrawPart& cp = *pcp;
      cp.is_show = true;
      cp.itype = Cad::EDGE;
      cp.id_cad = id_e;
      cp.nelem = 1;
      cp.npoel = 2;
      cp.pIndexArray = new unsigned int [2];
      cp.pIndexArray[0] = id_vs-1;
      cp.pIndexArray[1] = id_ve-1; 
      m_apIndexAry.push_back(pcp);
    }
  }
  
  {
    const std::vector<unsigned int>& aIdL = cad_3d.GetAryElemID(Cad::LOOP);
    for(unsigned int iil=0;iil<aIdL.size();iil++){
      unsigned int id_l = aIdL[iil];
      assert( cad_3d.IsElemID(Cad::LOOP,id_l) );
      Msh::CMesher2D msh2d( CLoopPxy2D(cad_3d,id_l) );
      const std::vector<Msh::CTriAry2D>& aTriAry = msh2d.GetTriArySet();
      const std::vector<CVector2D>& aVec2D = msh2d.GetVectorAry();
      if( aTriAry.size() != 1 ) continue;
      const unsigned int ntri = aTriAry[0].m_aTri.size();
//      std::cout << aTriAry[0].m_aTri.size() << std::endl;
//      std::cout << aVec2D.size() << std::endl;
      CDrawPart* pcp = new CDrawPart();
      CDrawPart& cp = *pcp;
      cp.is_show = true;
      cp.itype = Cad::LOOP;
      cp.id_cad = id_l;
      {
        const Cad::CLoop3D& l = cad_3d.GetLoop(id_l);
        const Com::CVector3D& n = l.normal;
        cp.is_const_normal = true;
        cp.normal[0] = n.x;
        cp.normal[1] = n.y;
        cp.normal[2] = n.z;        
      }      
      cp.nelem = ntri;
      cp.npoel = 3;
      cp.pIndexArray = new unsigned int [ntri*3];
      for(unsigned int itri=0;itri<ntri;itri++){
        cp.pIndexArray[itri*3+0] = aTriAry[0].m_aTri[itri].v[0];
        cp.pIndexArray[itri*3+1] = aTriAry[0].m_aTri[itri].v[1];
        cp.pIndexArray[itri*3+2] = aTriAry[0].m_aTri[itri].v[2];        
      }
      m_apIndexAry.push_back(pcp);      
    }
  }
      
  

  {
    for(unsigned int i=0;i<aIdV.size();i++){
      CDrawPart_CadVertex dpv;
      dpv.id_cad = aIdV[i];
      dpv.is_selected = false;
      dpv.iv = i;
      dpv.is_show = true;
      m_aIndexVertex.push_back(dpv);    
    }    
  }
	{	// 座標をセット
		const unsigned int npoin = aVec3D.size();
		const unsigned int ndim = 3;
		m_vertex_ary.SetSize(npoin,ndim);
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			m_vertex_ary.pVertexArray[ipoin*ndim  ] = aVec3D[ipoin].x;
			m_vertex_ary.pVertexArray[ipoin*ndim+1] = aVec3D[ipoin].y;
			m_vertex_ary.pVertexArray[ipoin*ndim+2] = aVec3D[ipoin].z;
		}
	}
	return true;
}

void CDrawer_Cad3D::AddSelected(const int selec_flag[])
{
  Cad::CAD_ELEM_TYPE type = Cad::NOT_SET;
  if(      selec_flag[1] == 1 ){ type = Cad::VERTEX; }
  else if( selec_flag[1] == 2 ){ type = Cad::EDGE;   }    
  else if( selec_flag[1] == 3 ){ type = Cad::LOOP;   }
  unsigned int id = selec_flag[2];
  if( type == Cad::LOOP || type == Cad::EDGE ){
    for(unsigned int idp=0;idp<m_apIndexAry.size();idp++){  
      CDrawPart* part = m_apIndexAry[idp];
      if( part->itype == type && part->id_cad == id ){ part->is_selected = true; }
    }
  }
  else if( type == Cad::VERTEX ){    
    for(unsigned int iver=0;iver<this->m_aIndexVertex.size();iver++){
      if( !this->m_aIndexVertex[iver].is_show ) continue;
      unsigned int id_cad = this->m_aIndexVertex[iver].id_cad;
      if( id_cad == id ){ this->m_aIndexVertex[iver].is_selected = true; }
    }
  }
}

void CDrawer_Cad3D::ClearSelected()
{
  for(unsigned int idp=0;idp<m_apIndexAry.size();idp++){  
    CDrawPart* part = m_apIndexAry[idp];
    part->is_selected = false;
  }  
  for(unsigned int iver=0;iver<this->m_aIndexVertex.size();iver++){
    this->m_aIndexVertex[iver].is_selected = false;    
  }
}

void CDrawer_Cad3D::DrawSelection(unsigned int idraw) const
{
  const bool is_blend       = ::glIsEnabled(GL_BLEND);
  const bool is_line_smooth = ::glIsEnabled(GL_LINE_SMOOTH);
  const bool is_texture     = ::glIsEnabled(GL_TEXTURE_2D);
  const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_TEXTURE_2D);
  ::glDisable(GL_LIGHTING);
  
  ::glPushName(idraw);
  
	const unsigned int ndim = this->m_vertex_ary.NDim();
  
	////////////////////////////////////////////////////////////////
	// モデルの描画
  { // draw vertecies
    ::glColor3d(0.0,0.0,0.0);      
    ////////////////
    ::glPointSize(m_pointsize);
    for(unsigned int iver=0;iver<this->m_aIndexVertex.size();iver++){
      if( !this->m_aIndexVertex[iver].is_show ) continue;
      unsigned int ipo0 = this->m_aIndexVertex[iver].iv;
      ::glPushName(1);
      ::glPushName(m_aIndexVertex[iver].id_cad);
      ::glBegin(GL_POINTS);      
      ::glVertex3d(m_vertex_ary.pVertexArray[ipo0*ndim+0], 
                   m_vertex_ary.pVertexArray[ipo0*ndim+1], 
                   m_vertex_ary.pVertexArray[ipo0*ndim+2]);
      ::glEnd();            
      ::glPopName();
      ::glPopName();
    }
    if( is_texture ){ glEnable(GL_TEXTURE_2D); }
  }  
  ////
  // regist vertex array
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);
	for(unsigned int idp=0;idp<m_apIndexAry.size();idp++){    
    const CDrawPart* part = m_apIndexAry[idp];
    if(      part->itype == Cad::EDGE ){ ::glPushName(2); }
    else if( part->itype == Cad::LOOP ){ ::glPushName(3); }
    ::glPushName( part->id_cad );
    part->DrawElements();
    ::glPopName();
    ::glPopName();
	}
  ::glDisableClientState(GL_VERTEX_ARRAY);
  ::glDisableClientState(GL_TEXTURE_COORD_ARRAY);      
  if( is_lighting    ){ ::glEnable(GL_LIGHTING);    } else{ ::glDisable(GL_LIGHTING);    }
  if( is_blend       ){ ::glEnable(GL_BLEND);       } else{ ::glDisable(GL_BLEND);       }
  if( is_line_smooth ){ ::glEnable(GL_LINE_SMOOTH); } else{ ::glDisable(GL_LINE_SMOOTH); }
  if( is_texture     ){ ::glEnable(GL_TEXTURE_2D);  } else{ ::glDisable(GL_TEXTURE_2D);  }
  ::glPopName();
}

void CDrawer_Cad3D::Draw() const
{
	::glEnable(GL_DEPTH_TEST);
  ::glDisable(GL_CULL_FACE);    // to make the program simple...
	const bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  const bool is_texture  = ::glIsEnabled(GL_TEXTURE_2D);  
  const bool is_blend    = ::glIsEnabled(GL_BLEND);
	const unsigned int ndim = this->m_vertex_ary.NDim();

	////////////////////////////////////////////////////////////////
	// モデルの描画
  { // draw vertecies
    ::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glColor3d(0.0,0.0,0.0);      
    ////////////////
    ::glPointSize(m_pointsize);
    ::glBegin(GL_POINTS);
    for(unsigned int iver=0;iver<this->m_aIndexVertex.size();iver++){
      if( !this->m_aIndexVertex[iver].is_show ) continue;
      if( this->m_aIndexVertex[iver].is_selected ){ ::glColor3d(1.0,1.0,0.0); }
      else{ ::glColor3d(0.0,0.0,0.0);	}
      unsigned int ipo0 = this->m_aIndexVertex[iver].iv;      
      ::glVertex3d(m_vertex_ary.pVertexArray[ipo0*ndim+0], 
                   m_vertex_ary.pVertexArray[ipo0*ndim+1], 
                   m_vertex_ary.pVertexArray[ipo0*ndim+2]);
    }
    ::glEnd();      
    if( is_texture ){ glEnable(GL_TEXTURE_2D); }
  }  
  ////
  // regist vertex array
  ::glLineWidth(m_linewidth);
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);
	for(unsigned int idp=0;idp<m_apIndexAry.size();idp++){
    const CDrawPart* part = m_apIndexAry[idp];
    if(      part->itype == Cad::LOOP ){
      if( is_lighting ){ 
        ::glEnable(GL_LIGHTING); 
        if( part->is_const_normal ){
          ::glNormal3d(part->normal[0],part->normal[1],part->normal[2]);
        }
      }
      if( part->is_selected ){
        ::glDisable(GL_LIGHTING);
				::glEnable(GL_POLYGON_STIPPLE);
				::glPolygonStipple((const GLubyte*)m_mask);
        ::glColor3d(1.0,1.0,0.0);
				part->DrawElements();
				::glDisable(GL_POLYGON_STIPPLE);
			}      
      ::glColor3d(0.8,0.8,0.8);
      part->DrawElements();
    }
    else if( part->itype == Cad::EDGE ){
      ::glDisable(GL_LIGHTING);
      if( part->is_selected ){ ::glColor3d(1.0,1.0,0.0);  }
			else{                    ::glColor3d(0.0,0.0,0.0);  }
      part->DrawElements();
    }
	}
  ::glDisableClientState(GL_VERTEX_ARRAY);
  if( is_lighting ){ ::glEnable(GL_LIGHTING);   } else{ ::glDisable(GL_LIGHTING);    }
  if( is_blend    ){ ::glEnable(GL_BLEND);      } else{ ::glDisable(GL_BLEND);      }
  if( is_texture  ){ ::glEnable(GL_TEXTURE_2D); } else{ ::glDisable(GL_TEXTURE_2D); }  
	return;
}

