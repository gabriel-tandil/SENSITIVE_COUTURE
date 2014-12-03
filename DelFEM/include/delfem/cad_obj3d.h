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
@brief interface of 3D cad class (Cad::CCadObj2Dm)
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_3D_H)
#define CAD_OBJ_3D_H

#include <vector>

#include "delfem/vector3d.h"
#include "delfem/vector2d.h"
#include "delfem/objset.h"
#include "delfem/cad_com.h"
#include "delfem/cad/brep2d.h"
#include "delfem/cad/cad_elem2d.h"
#include "delfem/cad/cad_elem3d.h"

namespace Cad{

class CVertex3D;
class CLoop3D;
class CEdge3D;
class CTopology;

/*! 
@brief 3 dimentional cad model class
@ingroup CAD
*/
class CCadObj3D
{
public:  
	//! iterator go around vertex
	class CItrVertex
	{
	public:		
		CItrVertex(const CCadObj3D* pCadObj3D, unsigned int id_v)
    : itrv(pCadObj3D->m_BRep.GetItrVertex(id_v)){}
		void operator++(){ itrv++; }	//!< go around (cc-wise) loop around vertex 
		void operator++(int n){ itrv++; }	//!< dummy operator (for ++)		
		//! cc-wise ahead  edge-id and its direction(true root of edge is this vertex)
		bool GetIdEdge_Ahead( unsigned int& id_e, bool& is_same_dir) const { return itrv.GetIdEdge_Ahead( id_e,is_same_dir); }
		//! cc-wise behind edge-id and its direction(true root of edge is this vertex)
		bool GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const { return itrv.GetIdEdge_Behind(id_e,is_same_dir); }
		unsigned int GetIdLoop() const{ return itrv.GetIdLoop(); } //!< get loop-id		
		bool IsEnd() const { return itrv.IsEnd(); }	//!< return true if iterator go around
	private:
		CBRepSurface::CItrVertex itrv;
	};
  friend class CItrVertex;  
public:
  ////////////////////////////////
	// constructor & destructor
  
	//! default constructor
	CCadObj3D(){}
  
  void Clear(){
    m_LoopSet.Clear();
    m_EdgeSet.Clear();
    m_VertexSet.Clear();
    m_BRep.Clear();
  }
  
	//! function gives iterator which travel edge and loop around vertex (ID:id_v)
  
  bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
	const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;  

//	unsigned int AddSolid(const std::vector<Com::CVector3D>& vec_ary, const Com::CVector3D& dir);
  unsigned int AddCuboid(double len_x, double len_y, double len_z);
  unsigned int AddPolygon( const std::vector<Com::CVector3D>& aVec, unsigned int id_l);
  unsigned int AddRectLoop(unsigned int id_l, const Com::CVector2D& p0, const Com::CVector2D& p1);
  void LiftLoop(unsigned int id_l, Com::CVector3D dir);
  unsigned int AddPoint( Cad::CAD_ELEM_TYPE type, unsigned int id_elem, const Com::CVector3D& po );
//  unsigned int AddVertex(const Com::CVector3D& dir);
  CBRepSurface::CResConnectVertex ConnectVertex(unsigned int id_v1, unsigned int id_v2);
  
  const Com::CVector3D& GetVertexCoord(unsigned int id_v) const;
	bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;  
  unsigned int GetIdLoop_Edge(unsigned int id_e, bool is_left){ return m_BRep.GetIdLoop_Edge(id_e,is_left); }
  
  CBRepSurface::CItrLoop GetItrLoop(unsigned int id_l) const{ 
    assert(m_BRep.IsElemID(Cad::LOOP,id_l)); 
    return CBRepSurface::CItrLoop(m_BRep,id_l); 
  }
  const Cad::CLoop3D& GetLoop(unsigned int id_l) const;
  unsigned int AssertValid();
private:
  bool CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& point) const;
  bool FindIntersectionEdge(unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe,
                            unsigned int& id_e_nearest, bool& is_same_dir_e_near,
                            Com::CVector2D& p_nearest);
  Cad::CEdge2D GetEdge2D(unsigned int id_e, unsigned int id_l) const;
protected:
	////////////////
	Com::CObjSet<CLoop3D>   m_LoopSet;
	Com::CObjSet<CEdge3D>   m_EdgeSet;
	Com::CObjSet<CVertex3D> m_VertexSet;
	////////////////
	CBRepSurface m_BRep;	// class which have topology
};

}	// end namespace CAD

#endif
