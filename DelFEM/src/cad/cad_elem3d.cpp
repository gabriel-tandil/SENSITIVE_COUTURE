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

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#include <iostream>
#include <cassert>
#include <cstdlib>	// abort

#include "delfem/cad/cad_elem3d.h"

const Com::CBoundingBox3D& Cad::CLoop3D::GetBoundingBox() const
{
  if( bb_.isnt_empty ){ return bb_; }
  Com::CBoundingBox2D bb2;
  for(unsigned int iie=0;iie<aEdge.size();iie++){
    bb2 += aEdge[iie].first.GetBoundingBox();
  }
  const double ep0 = (bb2.x_max-bb2.x_min+bb2.y_max-bb2.y_min)*1.0e-10;
  bb_.AddPoint( this->UnProject( Com::CVector2D(bb2.x_min,bb2.y_min) ), ep0 );
  bb_.AddPoint( this->UnProject( Com::CVector2D(bb2.x_max,bb2.y_min) ), ep0 );
  bb_.AddPoint( this->UnProject( Com::CVector2D(bb2.x_min,bb2.y_max) ), ep0 );
  bb_.AddPoint( this->UnProject( Com::CVector2D(bb2.x_max,bb2.y_max) ), ep0 );
  return bb_;  
}

Com::CVector3D Cad::CLoop3D::GetNearestPoint(const Com::CVector3D& p) const
{
  Com::CVector2D vf2 = this->Project(p);
  bool is_inside = this->CheckPointInside2D(vf2);
  if( is_inside ){ return this->UnProject(vf2); }
  Com::CVector2D mfe = this->GetNearestPointInEdge2D(vf2);
  return this->UnProject(mfe);
}

int Cad::CLoop3D::NumIntersecRay(const Com::CVector3D& org0, const Com::CVector3D& dir0) const
{
  double t0 = Dot(normal,org-org0);
  double t1 = Dot(normal,dir0);
  if( fabs(t0) < 1.0e-10*fabs(t1) ){
    Com::CVector2D vf2 = this->Project(org0);
    Com::CVector3D vf3 = this->UnProject(vf2);
    const double sqheight = SquareDistance(org0,vf3);    
    if( sqheight < 1.0e-10 ) return -1;
    return 0; // pararell
  }
  const double r = t0/t1;
  if( r < -1.0e-10 ) return 0;
  Com::CVector3D vi3 = org0+r*dir0;
  Com::CVector2D vi2 = this->Project(vi3);  
  bool is_inside = this->CheckPointInside2D(vi2);
  if( is_inside ){
    if( r < 1.0e-10 ) return -1;
    return 1;
  }
  return 0;
}

bool Cad::CLoop3D::CheckPointInside2D(const Com::CVector2D& p) const
{  
	for(unsigned int i=1;i<29;i++){	// 29 is handy prim number
		unsigned int cross_counter = 0;
		bool iflg = true;
    Com::CVector2D dir(sin(6.28*i/29.0),cos(6.28*i/29.0));
		for(unsigned int iie=0;iie<aEdge.size();iie++){
			const CEdge2D& e = aEdge[iie].first;
			const int ires = e.NumIntersect_AgainstHalfLine( p,dir );
			// -1 is vague so let's try again!
			if( ires == -1 ){ iflg = false; break; }
			cross_counter += ires;
		}
		if( iflg == true ){ 
			if( cross_counter % 2 == 0 ) return false;
			return true;
		}
	}
	assert(0);	// I hope process don't come here !!!!
	return false;  
}
                 

Com::CVector2D Cad::CLoop3D::GetNearestPointInEdge2D(const Com::CVector2D& p) const
{
  assert( aEdge.size() > 0 );
  double min_dist = -1; 
  Com::CVector2D p_nearest;
  for(unsigned int iie=0;iie<aEdge.size();iie++){
    const CEdge2D& e = aEdge[iie].first;
    const Com::CVector2D& p_near = e.GetNearestPoint(p);
    const double dist = Distance(p,p_near);
    if( min_dist > 0 && dist > min_dist ){ continue; }
    min_dist = dist; 
    p_nearest = p_near;
  }  
  return p_nearest;
}
                                    
                                    


