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

/*! @file
@brief CADを可視化するためのクラス(Cad::View::CDrawer_Cad2D)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_CAD_3D_H)
#define DRAWER_CAD_3D_H

#include <vector>

#include "cad_obj3d.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

namespace Msh{
	class CBarAry;
	class CTriAry2D;
}

namespace Cad{
namespace View{

/*! 
@brief CadをOpenGLで可視化するクラス
@ingroup CAD
*/
class CDrawer_Cad3D : public Com::View::CDrawer
{
public : 
  CDrawer_Cad3D(const Cad::CCadObj3D& cad)
  {    
		this->m_is_anti_aliasing = false;
    m_linewidth = 3;
    m_pointsize = 5;
		// 2ptのチェック柄を作る
		for(unsigned int j=0;j<4;j++){
			unsigned int i;
			for(i=0; i<8 ;i++) m_mask[j*32   +i] = 0x33;
			for(i=0; i<8 ;i++) m_mask[j*32+ 8+i] = 0xcc;
			for(i=0; i<8 ;i++) m_mask[j*32+16+i] = 0x33;
			for(i=0; i<8 ;i++) m_mask[j*32+24+i] = 0xcc;
		}	
    this->UpdateCAD_TopologyGeometry(cad);
	}  
	virtual ~CDrawer_Cad3D(){
		for(unsigned int i=0;i<m_apIndexAry.size();i++){ delete m_apIndexAry[i]; }
	}
	// virtual関数
  
  //! トポロジーと幾何を更新する
  bool UpdateCAD_TopologyGeometry(const Cad::CCadObj3D&);

	//! 描画
	virtual void Draw() const;
  virtual void DrawSelection(unsigned int) const;
	virtual void AddSelected(const int selec_flag[]);
	virtual void ClearSelected();
  //! 線の太さを設定
  void SetLineWidth(unsigned int linewidth){ m_linewidth = linewidth; }
  //! 点の大きさを設定
  void SetPointSize(unsigned int pointsize){ m_pointsize = pointsize; }

	//! バウンディング・ボックスを得る
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const {
		return m_vertex_ary.GetBoundingBox(	rot );
	}

private:	
	////////////////
	class CDrawPart{
	public : 
		CDrawPart(){
			nelem = 0; npoel = 0; pIndexArray = 0; 
			is_selected = false; itype = Cad::VERTEX; is_show=true;
      color[0]=0; color[1]=0; color[2]=0;
//			xdisp=0; ydisp=0;
		}
		virtual ~CDrawPart(){
			this->Clear();
		}
		void Clear(){
			delete[] pIndexArray;
			pIndexArray = 0;
			nelem = 0;
			npoel = 0;
		}
		void DrawElements() const;
//		bool Set(const Msh::CBarAry& BarAry);
//		bool Set(const Msh::CTriAry2D& TriAry);
//		void SetHeight(double h){ height = h; }
	public:
		bool is_show;
		bool is_selected;
    float color[3];
    bool is_const_normal;
    double normal[3];
    ////////////////
		Cad::CAD_ELEM_TYPE itype;
		unsigned int id_cad;
//    unsigned int id_msh;
    ////////////////
		unsigned int nelem;
		unsigned int npoel;
    unsigned int* pIndexArray;
		////////////////
//		double height;
//		double xdisp, ydisp;
	};
	////////////////
	class CDrawPart_CadVertex{
	public:
		unsigned int iv;
		unsigned int id_cad;
		unsigned int id_msh;
    ////////////////
		bool is_selected;
		bool is_show;
		////////////////
		double height;
	};	
private:
  unsigned char m_mask[128];	//!< mask showing selected region with dots
  std::vector<CDrawPart*> m_apIndexAry;				//! edge and loops
  std::vector<CDrawPart_CadVertex> m_aIndexVertex;	//! vertex
  Com::View::CVertexArray m_vertex_ary;	//! vertex array
  unsigned int m_linewidth;   //!< the width of the line
  unsigned int m_pointsize;   //!< the size of vertices
};

} // end namespace View
} // end namespace Cad

#endif
