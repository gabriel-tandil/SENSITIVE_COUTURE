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
@brief 面で場を可視化するクラス(Fem::Field::View::CDrawerFaceContour)のインターフェース
@author Nobuyuki Umetani
*/


#if !defined(DRAWER_FIELD_IMAGE_BASED_FLOW_VIS_H)
#define DRAWER_FIELD_IMAGE_BASED_FLOW_VIS_H

#include "delfem/drawer_field.h"
#include "delfem/field_world.h"

namespace Fem{
namespace Field{
namespace View{
/*
class CRect
{
public:
	CRect(double min_x, double max_x, double min_y, double max_y,    double r, double g, double b)
		: min_x(min_x), max_x(max_x), min_y(min_y), max_y(max_y){
		color[0] = r;	color[1] = g;	color[2] = b;
		is_active = true;
	}
	double min_x, max_x, min_y, max_y;
	double color[3];
	bool is_active;
};
*/
class CEdgeTextureColor
{
public:
	CEdgeTextureColor(unsigned int id_field_velo, unsigned int id_ea, Fem::Field::CFieldWorld& world, 
		double r, double g, double b)
	{
		color[0]=r;	color[1]=g;	color[2]=b;
		nelem = 0;
		m_aXYVeloElem = 0;
		m_aXYElem = 0;
        velo_scale = 0.32;
		this->Set(id_field_velo,id_ea,world);
	}
	CEdgeTextureColor( const CEdgeTextureColor& rhs )
	{
		color[0] = rhs.color[0];
		color[1] = rhs.color[1];
		color[2] = rhs.color[2];
		nelem = rhs.nelem;
		m_aXYVeloElem = rhs.m_aXYVeloElem;
		m_aXYElem = rhs.m_aXYElem;
		velo_scale = rhs.velo_scale;
		id_ea = rhs.id_ea;
		id_field_velo = rhs.id_field_velo;
		id_part_field_velo = rhs.id_part_field_velo;
	}
	virtual ~CEdgeTextureColor(){
		if( m_aXYElem     != 0 ){ delete[] m_aXYElem; }
		if( m_aXYVeloElem != 0 ){ delete[] m_aXYVeloElem; }
	}
	bool Set(unsigned int id_field_velo, unsigned int id_ea, Fem::Field::CFieldWorld& world){
		this->id_ea = id_ea;
		this->id_field_velo = id_field_velo;
		assert( world.IsIdField(id_field_velo) );
		assert( world.IsIdEA(id_ea) );
		if( !world.IsIdField(id_field_velo) ) return false;
		if( !world.IsIdEA(id_ea) ) return false;
		const Fem::Field::CField& fv = world.GetField(id_field_velo);
		const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
		assert( ea.ElemType() == Fem::Field::LINE );
		if( ea.ElemType() != Fem::Field::LINE ) return false;
		id_part_field_velo = world.GetPartialField(this->id_field_velo,id_ea);
		std::cout << "part field id : " << id_part_field_velo << std::endl;
		return true;
	}
	bool Update(const Fem::Field::CFieldWorld& world);
	void Draw() const;
private: 
	////////////////
	unsigned int id_ea;
	unsigned int id_field_velo;
	unsigned int id_part_field_velo;
	////////////////
	double velo_scale;
	////////////////
	unsigned int nelem;
	double* m_aXYVeloElem;
	double* m_aXYElem;
	////////////////
	double color[3];
};

//! 辺の描画クラス
class CDrawerImageBasedFlowVis : public CDrawerField
{
public:
	CDrawerImageBasedFlowVis();
	// imode (0:拡散なし) (1:ランダムノイズ) (2:格子) (3:ドット)
	CDrawerImageBasedFlowVis(unsigned int id_field, const Fem::Field::CFieldWorld& world, unsigned int imode = 1 );
	virtual ~CDrawerImageBasedFlowVis(){
		for(unsigned int i=0;i<this->m_apIndexArrayElem.size();i++){
			delete this->m_apIndexArrayElem[i];
		}
		if( aCoord    != 0 ){ delete[] aCoord; }
		if( aVelo     != 0 ){ delete[] aVelo; }
		if( aValColor != 0 ){ delete[] aValColor; }
		this->ClearDisplayList();
	}
	////////////////	
	virtual Com::CBoundingBox3D GetBoundingBox( double rot[] ) const;
    virtual void DrawSelection(unsigned int idraw) const{}
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual void Draw() const;
	virtual bool Update(const Fem::Field::CFieldWorld& world);
	void AddFlowInOutEdgeColor(unsigned int id_e, Fem::Field::CFieldWorld& world, 
		double r, double g, double b)
    {
		aEdgeColor.push_back( CEdgeTextureColor(m_IdFieldVelo,id_e, world, r,g,b) );
	}
	void SetColorField(unsigned int id_field_color, const CFieldWorld& world, std::auto_ptr<CColorMap> color_map);
private:
	bool Set(unsigned int id_field, const Fem::Field::CFieldWorld& world);
	void ClearDisplayList();
	void MakePattern();
public:
	std::vector<CEdgeTextureColor> aEdgeColor;
private:
    ////////////////
	unsigned int m_IdFieldVelo;
	unsigned int m_IdFieldColor;
	unsigned int m_nPattern;
    mutable unsigned int iPtn;    // フレームの番号　(const関数のDraw()でカウントアップする必要があるのでmutable)
	unsigned int m_nameDisplayList;
	////////////////
	// 0:パターンなし(移流拡散)
	// 1:ランダムノイズ
	// 2:格子
	// 3:ドット
    unsigned int imode; // パターンの種類
	////////////////
	int alpha;	// 0-256
    double velo_scale;
    ////////////////
	std::vector<CIndexArrayElem*> m_apIndexArrayElem;
	unsigned int nnode;
	double* aVelo;
	double* aCoord;
	double* aValColor;
//	unsigned int nelem;
//	double* aXYVeloElem;
//	double* aXYElem;
//	double* aColorElem;

	std::auto_ptr<CColorMap> color_map;
};


}
}
}

#endif
