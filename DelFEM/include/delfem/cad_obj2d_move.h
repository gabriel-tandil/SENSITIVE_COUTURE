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
@brief interface of class (Cad::CCadObj2D_Move)
@author Nobuyuki Umetani
*/

#if !defined(CAD_OBJ_2D_MOVE_H)
#define CAD_OBJ_2D_MOVE_H

#include <set>

#include "delfem/vector2d.h"
#include "delfem/cad_obj2d.h"

#include "delfem/cad/cad_edge2d_polyline.h"

namespace Cad{

/*! 
@brief 2D cad class that element can be moved
@ingroup CAD
*/
class CCadObj2D_Move : public Cad::CCadObj2D
{
public:
	//! move vertex (ID:id_v) to vec_dist
	bool MoveVertex( unsigned int id_v, const Com::CVector2D& vec);
  //! move vertices to each distinations 
	bool MoveVertex( const std::vector< std::pair<unsigned int,Com::CVector2D> >& vec );
	//! move veretex (ID:id_e) with direction vec_delta
	bool MoveEdge(unsigned int id_e, const Com::CVector2D& vec_delta);
	//! move loop (ID:id_l) with direction vec_delta
	bool MoveLoop(unsigned int id_l, const Com::CVector2D& vec_delta);

  //! if edge (ID:id_e) is an arc make it go through point(dist)
  bool DragArc(unsigned int id_e, const Com::CVector2D& dist);
  //! if edge (ID:id_e) is an polyline make it go through point(dist)
  bool DragPolyline(unsigned int id_e, const Com::CVector2D& dist);
  
  bool PreCompDragPolyline(unsigned int id_e, const Com::CVector2D& pick_pos);
  // smoothing edge (id_e) if radius is negative smooth whole edge 
  bool SmoothingPolylineEdge(unsigned int id_e, unsigned int niter,
                             const Com::CVector2D& pos, double radius);
  
private:
  CCadEdge2DPolyline polyline;
};

}	// end namespace CAD

#endif
