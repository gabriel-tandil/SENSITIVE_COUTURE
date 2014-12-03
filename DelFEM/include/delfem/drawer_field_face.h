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
@brief Interfacde for field visualization class of element surface(Fem::Field::View::CDrawerFace)
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_FIELD_FACE_H)
#define DRAWER_FIELD_FACE_H

#include <memory>

#include "delfem/drawer_field.h"

namespace Fem{
namespace Field{
namespace View{

//! drawing face class
class CDrawerFace : public CDrawerField
{
public:
	CDrawerFace();
	CDrawerFace(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world, unsigned int id_field_color=0);
	CDrawerFace(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world, unsigned int id_field_color, double min, double max);
	CDrawerFace(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world, unsigned int id_field_color, std::auto_ptr<CColorMap> color_map );
	virtual ~CDrawerFace();
  
	////////////////////////////////
	// declaration of virtual functions
  
	virtual void Draw() const;
	virtual void DrawSelection(unsigned int idraw) const{};
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const{
		return m_vertex_ary.GetBoundingBox(	rot );
	}
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual bool Update(const Fem::Field::CFieldWorld& world);
  
	////////////////////////////////
	// declaration of non-virtual functions
  
  void SetColor(double r, double g, double b, unsigned int id_ea = 0){
    const unsigned int niea = m_apIndexArrayElem.size();
    if( id_ea == 0 ){
      for(unsigned int i=0;i<niea;i++){ m_apIndexArrayElem[i]->SetColor(r,g,b); }
      return;
    }
    for(unsigned int i=0;i<niea;i++){
      if( m_apIndexArrayElem[i]->GetIdEA() != id_ea ){ continue; }
      m_apIndexArrayElem[i]->SetColor(r,g,b);
    }
  }
	void SetColor(unsigned int id_es_v, unsigned int id_ns_v, const Fem::Field::CFieldWorld& world,
		const std::auto_ptr<CColorMap>& color_map);
  
  void EnableNormal(bool is_lighting);
  void EnableUVMap(bool is_uv_map, const Fem::Field::CFieldWorld& world);
  void SetTexCenter(double cent_x, double cent_y){
    this->tex_cent_x = cent_x;
    this->tex_cent_y = cent_y;    
  }
  void GetTexCenter(double& cent_x, double& cent_y){
    cent_x = this->tex_cent_x;
    cent_y = this->tex_cent_y;
  }  
  void SetTexScale(double scale, const Fem::Field::CFieldWorld& world);
protected:
	bool Set(unsigned int id_field, const Fem::Field::CFieldWorld& world, bool isnt_value_disp, unsigned int id_field_color);
  void MakeNormal();
protected:
	std::vector<CIndexArrayElem*> m_apIndexArrayElem;
	Com::View::CVertexArray m_vertex_ary;

	unsigned int m_id_field;
	bool m_isnt_value_disp;
	bool m_is_draw_nsv;	//!< value‚ÌNS‚ð•`‰æ‚·‚é‚©coord‚ÌNS‚ð•`‰æ‚·‚é‚©DValue‚Ìns_c‚ª–³‚¯‚ê‚ÎCoord‚Ìns_c‚ð•`‰æ‚·‚é

	////////////////
	// color
  float* pColorArray;	//!< array of color ( rgb for each node )
	unsigned int id_field_val;
	std::auto_ptr<CColorMap> color_map; // use this function to calc rgb from field_val 
	bool is_draw_color_legend;	//!< draw legent if true

  
  ////////////////
  // Lighting
  double* pNormalArray;
  
  ////////////////
  // uv map
  double* pUVArray;  
  double tex_cent_x, tex_cent_y;
  double tex_scale;
};

}	// View
}	// Field
}	// Fem


#endif
