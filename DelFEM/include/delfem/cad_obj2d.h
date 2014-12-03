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
@brief interface of 2D cad class (Cad::CCadObj2Dm)
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_2D_H)
#define CAD_OBJ_2D_H

#include <vector>

#include "delfem/vector2d.h"
#include "delfem/serialize.h"
#include "delfem/cad2d_interface.h"
#include "delfem/objset.h"
#include "delfem/cad/brep2d.h"
#include "delfem/cad/cad_elem2d.h"

namespace Cad{

class CVertex2D;
class CLoop2D;
class CEdge2D;
class CTopology;

/*! 
@brief 2 dimentional cad model class
@ingroup CAD
*/
class CCadObj2D : public Cad::ICad2D_Msh
{
public:     
  class CResAddVertex{
  public:
    CResAddVertex(){
      id_v_add = 0;
      id_e_add = 0;
    }
  public:
    unsigned int id_v_add;
    unsigned int id_e_add;
  };
  class CResAddPolygon{
  public:
    CResAddPolygon(){
      id_l_add = 0;
    }
    CResAddPolygon(const CResAddPolygon& lhs){
      this->id_l_add = lhs.id_l_add;
      this->aIdV = lhs.aIdV;
      this->aIdE = lhs.aIdE;
    }
  public:
    unsigned int id_l_add;
    std::vector<unsigned int> aIdV;
    std::vector<unsigned int> aIdE;
  };
  
  
  ////////////////////////////////
	// constructor & destructor
  
	//! default constructor
	CCadObj2D();
  //! copy constructor
  CCadObj2D(const CCadObj2D& cad);
	//! destructor
	virtual ~CCadObj2D();
	//! initialization clear all element
	void Clear();

	////////////////////////////////
	// Get method

	//! function gives iterator which travel vtx and edge inside the loop (ID:id_l)
  virtual std::auto_ptr<Cad::IItrLoop> GetPtrItrLoop(unsigned int id_l) const{
    return std::auto_ptr<Cad::IItrLoop>( new CBRepSurface::CItrLoop(m_BRep,id_l) );	// instance
	}
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;  
  virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const;
  unsigned int GetIdVertex_Edge(unsigned int id_e, bool is_s) const;
  virtual bool GetIdLoop_Edge(unsigned int &id_l_l, unsigned int& id_l_r, unsigned int id_e) const;    
  CBRepSurface::CItrVertex GetItrVertex(unsigned int id_v) const { return CBRepSurface::CItrVertex(m_BRep,id_v); }
  CBRepSurface::CItrLoop GetItrLoop(  unsigned int id_l) const { return CBRepSurface::CItrLoop(m_BRep,id_l); }
    
	// functions related to layer
  virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const;
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const;
  bool ShiftLayer_Loop(unsigned int id_l, bool is_up);
  
  double GetMinClearance() const { return min_clearance; }
  
	// loop functions
	//! @{
	bool CheckIsPointInsideLoop(unsigned int id_l1, const Com::CVector2D& point) const;    
  double SignedDistPointLoop(unsigned int id_l1, const Com::CVector2D& point, unsigned int id_v_ignore=0) const;  
  //! get color(double[3]) of loop(ID:id_l), return false if there is no loop(ID:id_l)
  virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const;
  //! ID:id_l set color of loop
  virtual bool SetColor_Loop(unsigned int id_l, const double color[3] );
	//! ID:id_l return are of the loop
	virtual double GetArea_Loop(unsigned int id_l) const;
	//! @}

  ////////////////////////////////
	// Edge member functions

	const CEdge2D& GetEdge(unsigned int id_e) const;  

  // Get Geometric Type of the Curve  (0:line, 1:arc, 2:polyline)
	virtual int GetEdgeCurveType(const unsigned int& id_e) const;
  
	//! get information of edge(ID:id_e)
  
  // Get Arc property if this curve is not arc, this returns false
  // if is_left_side is true, this arc lies left side from the line connect start and end point of this edge (ID:id_e).
  // The dist means how far is the center of the circle from the line connect start and end point of edge (ID:id_e).
	virtual bool GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const;
	virtual bool GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const;
  
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen = -1) const;    
  // Get edge (ID:id_e) geometry as polyoine wich have ndiv divisions. The start and end points is ps and pe each
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
                                 unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const;

  bool GetPointOnCurve_OnCircle(unsigned int id_e,
                                const Com::CVector2D& v0, double len, bool is_front,
                                bool& is_exceed, Com::CVector2D& out) const;
  //! Get point on the edge (ID:id_e) that is nearest from point (po_in)
  Com::CVector2D GetNearestPoint(unsigned int id_e, const Com::CVector2D& po_in) const;

	////////////////////////////////
	// Vertex

  // get position of the vertex
	virtual Com::CVector2D GetVertexCoord(unsigned int id_v) const;

	////////////////////////////////////////////////
	// Toplogy affecting shape edit functions
	
  // Add Polyton to loop
  // The vec_ary is a array of vertex points (both clockwise and anti-clockwise is possible.)
  // If id_l is 0 or omitted, the loop will be added outside.
  // This returns CResAddPolygon which contains IDs of all vertices and edges.
	CResAddPolygon AddPolygon(const std::vector<Com::CVector2D>& vec_ary, unsigned int id_l = 0);

  // Add vertex to Cad element
  // add vertex (with the position vec) to elemnet (type:itype,id:ielem)
  // if itype is Cad::NOT_SET, the vertex is added outside of cad shape
  CResAddVertex AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id_elem, const Com::CVector2D& vec);
	
  // Remove CAD element
  // return true if this operation was sucessfull. 
  // in case it returns false, the cad shape is intact.
	bool RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id_elem);

	CBRepSurface::CResConnectVertex ConnectVertex(CEdge2D edge);
  CBRepSurface::CResConnectVertex ConnectVertex_Line(unsigned int id_v1, unsigned int id_v2){
		Cad::CEdge2D e(id_v1,id_v2);
		return this->ConnectVertex(e);    
	}
	//! @}

	////////////////////////////////////////////////
	// Geometry editing functions (topoloty intact)

  bool SetCurve_Polyline(const unsigned int id_e);
	//! set edge (ID:id_e) mesh
	bool SetCurve_Polyline(unsigned int id_e, const std::vector<Com::CVector2D>& aVec);
	//! set edge (ID:id_e) arc 
  bool SetCurve_Arc(const unsigned int id_e, bool is_left_side=false, double rdist=10.0);
  //! set edge (ID:id_e) straight line
	bool SetCurve_Line(const unsigned int id_e);

	////////////////////////////////////////////////
	// IO routines
	
	virtual bool WriteToFile_dxf(const std::string& file_name, double scale) const;
	bool Serialize( Com::CSerializer& serialize );		
protected:
	// return edge with vertex id(id_vs, id_ve) and vertex coord  
	CEdge2D& GetEdgeRef(unsigned int id_e);  
	int AssertValid() const;
	bool CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& p1) const;
  double DistPointItrLoop(CBRepSurface::CItrLoop& itrl, const Com::CVector2D& point) const;
  
  
  // ret:0 all the points in id_ul1 are inside id_ul2
  // ret:1 there are both inside/outside id_ul1 points in id_ul1 or ambiguous
  // ret:2 all the points in id_ul1 are outside id_ul2
	unsigned int CheckInOut_ItrLoopPoint_ItrLoop(CBRepSurface::CItrLoop& itrl1, CBRepSurface::CItrLoop& itrl2) const;
	
	// assert loop : return 0 if loop is OK
  int CheckLoop(unsigned int id_l) const;

	// ret iterator is setted to the half-edge that meet from half line between id_v and point anti-clockwise for the first time
  // the loop that this half-edge belongs to is overlapped with this half line
	CBRepSurface::CItrVertex FindCorner_HalfLine(unsigned int id_v, const Com::CVector2D& point) const;
	bool CheckIntersection_Loop(unsigned int id_l=0) const;
  bool CheckIntersection_EdgeAgainstLoop(const CEdge2D& edge,unsigned int id_l=0) const;
	double GetArea_ItrLoop(CBRepSurface::CItrLoop& itrl) const;
protected:
	////////////////
	Com::CObjSet<CLoop2D>   m_LoopSet;
	Com::CObjSet<CEdge2D>   m_EdgeSet;
	Com::CObjSet<CVertex2D> m_VertexSet;
	////////////////
	CBRepSurface m_BRep;	// class which have topology
  double min_clearance;
};

}	// end namespace CAD

#endif
