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
@brief メッシュ描画クラス(CDrawerMsh2D,CDrawerMsh3D)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_MSH_H)
#define DRAWER_MSH_H

#include <vector>

#include "delfem/drawer.h"
#include "delfem/vector3d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"

namespace Msh{
	class CMesher2D;
	class CMesh3D;

namespace View
{
class CDrawPart;
class CDrawPart_MshVertex{
public:
	unsigned int id_v;
	unsigned int id_cad;
	unsigned int id_msh;
	bool is_selected;
	double height;
};

/*! 
@brief ２次元メッシュ描画クラス
@ingroup Msh2D
*/
class CDrawerMsh2D : public Com::View::CDrawer{
public : 
	CDrawerMsh2D(const Msh::CMesher2D& msh, bool is_draw_face = true){ 
        this->is_draw_face = is_draw_face;
        this->Set(msh); 
        SetLineWidth(1);
    }
	virtual ~CDrawerMsh2D();

	virtual void Draw() const;
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const{
		return m_vertex_ary.GetBoundingBox( rot );
	}
    void SetLineWidth(unsigned int iwidth);

	//! 座標を更新するためのルーティン
	void UpdateCoord(const Msh::CMesher2D& msh );

	/*! 
	@brief セレクションバッファへの書き出し
	@param[in] idraw 名前付けの最初に付けられる数
	*/
	virtual void DrawSelection(unsigned int idraw) const;
	virtual void AddSelected(const int selec_flag[]);
	virtual void ClearSelected();
	void GetMshPartID(int selec_flag[], unsigned int& msh_part_id);
private:
	bool Set(const Msh::CMesher2D& msh);
private:
	bool is_front_and_back;
	std::vector<CDrawPart*> m_pIndexArySet;
	std::vector<CDrawPart_MshVertex> m_DrawPartCadVertex;
	Com::View::CVertexArray m_vertex_ary;
    bool is_draw_face;
    unsigned int line_width;
};

////////////////////////////////
////////////////////////////////

/*! 
@brief ３次元メッシュ描画クラス
@ingroup Msh3D
*/
class CDrawerMsh3D : public Com::View::CDrawer
{
public : 
	CDrawerMsh3D(const Msh::CMesh3D& msh){ this->Set(msh); }
	virtual ~CDrawerMsh3D();
	// virtual関数
	virtual void Draw() const;
	virtual void DrawSelection(unsigned int idraw) const;
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const {
		return m_vertex_ary.GetBoundingBox(	rot );
	}
	virtual void AddSelected(const int selec_flag[]);
	virtual void ClearSelected();
	virtual void Hide(unsigned int id_msh);
    void SetColor(unsigned int id_msh, double r, double g, double b);
private:
	bool Set(const Msh::CMesh3D& msh);
private:
	bool is_front_and_back;
	std::vector<CDrawPart*> m_pIndexArySet;
	std::vector<CDrawPart_MshVertex> m_DrawPartCadVertex;
	Com::View::CVertexArray m_vertex_ary;
};

}	// end namespace View
}	// end namespace Msh


#endif
