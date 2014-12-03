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
@brief 辺で場を可視化するクラス(Fem::Field::View::CDrawerEdge)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_FIELD_EDGE_H)
#define DRAWER_FIELD_EDGE_H

#include <memory>

#include "delfem/drawer_field.h"

namespace Fem{
namespace Field{
namespace View{

//! 辺の描画クラス
class CDrawerEdge : public CDrawerField{
public:
	CDrawerEdge();
	CDrawerEdge(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world );
	virtual ~CDrawerEdge();
	Com::CBoundingBox3D GetBoundingBox( double rot[] ) const;
	virtual void DrawSelection(unsigned int idraw) const{};
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual void Draw() const;
	virtual bool Update(const Fem::Field::CFieldWorld& world);
  void SetLineWidth(unsigned int iw){ this->line_width_ = iw; }
private:
	bool Set(unsigned int id_field, bool isnt_value_disp, const Fem::Field::CFieldWorld& world );
private:
	unsigned int m_nline;	// 辺の数
  unsigned int line_width_;
	Com::View::CVertexArray* m_paVer;	// 頂点配列
	unsigned int m_IdField;
	bool isnt_value_disp;
	std::vector<unsigned int> m_EdgeAry;
};


}
}
}

#endif
