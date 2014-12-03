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
@brief Interfaces define the geometry of 2d cad elements
@author Nobuyuki Umetani
*/

#if !defined(CAD_ELEM_3D_H)
#define CAD_ELEM_3D_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <assert.h>
#include <iostream> // needed only in debug

#include "delfem/vector3d.h"
#include "delfem/vector2d.h"
#include "delfem/cad/cad_elem2d.h"

////////////////////////////////////////////////////////////////

namespace Cad{

/*!
@addtogroup CAD
*/
//!@{

//! 2dim loop class
class CLoop3D{
public:
  CLoop3D() : normal(0,0,1), dirx(1,0,0), org(0,0,0){}    
	CLoop3D(const CLoop3D& rhs){ 
    this->normal   = rhs.normal; 
    this->org      = rhs.org; 
    this->dirx     = rhs.dirx; 
    this->aEdge    = rhs.aEdge;
    this->aIndEdge = rhs.aIndEdge;
    this->bb_      = rhs.bb_;
  }
	CLoop3D(const Com::CVector3D& o, const Com::CVector3D& n, const Com::CVector3D& x ) : org(o), normal(n), dirx(x){
    normal.Normalize();
    dirx.Normalize();
  }
public:
  Com::CVector2D Project(const Com::CVector3D& p) const{
    double x = Dot((p-org),dirx);
    double y = Dot((p-org),Cross(normal,dirx));
    return Com::CVector2D(x,y);
  }
  Com::CVector3D UnProject(const Com::CVector2D& p) const{
    return org+dirx*p.x+Cross(normal,dirx)*p.y;
  }
  
  ////
  Com::CVector3D GetNearestPoint(const Com::CVector3D& p) const;
  int NumIntersecRay(const Com::CVector3D& org, const Com::CVector3D& dir) const;
  const Com::CBoundingBox3D& GetBoundingBox() const;
private:
  bool CheckPointInside2D(const Com::CVector2D& p) const;
  Com::CVector2D GetNearestPointInEdge2D(const Com::CVector2D& p) const;
public:
  Com::CVector3D org, normal, dirx; // this is when this loop is planer
  mutable std::vector< std::pair<Cad::CEdge2D,bool> > aEdge;
  mutable std::vector<unsigned int> aIndEdge;
  mutable Com::CBoundingBox3D bb_;
};
  

//! 2dim edge
class CEdge3D{
public:
	CEdge3D()
	{
    id_v_s = 0;
    id_v_e = 0;
		po_s = Com::CVector3D(0,0,0);
		po_e = Com::CVector3D(0,0,0);
	}
public:
  mutable unsigned int id_v_s, id_v_e;	//!< start vertex
  mutable Com::CVector3D po_s, po_e;
};

//! ２次元幾何頂点クラス
class CVertex3D{
public:
	CVertex3D(const Com::CVector3D& point) : point(point){}	
	CVertex3D(const CVertex3D& rhs)
		: point(rhs.point){}	
public:
  Com::CVector3D point;   //!< coordinate
};


//! @}
}

#endif
