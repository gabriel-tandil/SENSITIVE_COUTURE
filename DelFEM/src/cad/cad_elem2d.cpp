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
#include <cstdlib>	// abort

#include "delfem/cad/cad_elem2d.h"


// 辺のバウンディングボックスを得る
// po_s, po_eに値がセットされていないと正常に動かないので注意
void Cad::CEdge2D::GetBoundingBox( double& x_min, double& x_max, double& y_min, double& y_max ) const
{
	x_min = ( po_s.x < po_e.x ) ? po_s.x : po_e.x;
	x_max = ( po_s.x > po_e.x ) ? po_s.x : po_e.x;
	y_min = ( po_s.y < po_e.y ) ? po_s.y : po_e.y;
	y_max = ( po_s.y > po_e.y ) ? po_s.y : po_e.y;
	if( itype == 0 ){ return; }	// 直線の場合はそのまま
	else if( itype == 1 ){	// 円弧の場合
		Com::CVector2D po_c;
		double radius;
		this->GetCenterRadius(po_c,radius);
		Com::CVector2D tmp_v;

		tmp_v.x = po_c.x+radius;
		tmp_v.y = po_c.y;
		if( this->IsDirectionArc(tmp_v) == 1 ){
			x_max = ( tmp_v.x > x_max ) ? tmp_v.x : x_max; 
		}
		tmp_v.x = po_c.x-radius;
		tmp_v.y = po_c.y;
		if( this->IsDirectionArc(tmp_v) == 1 ){
			x_min = ( tmp_v.x < x_min ) ? tmp_v.x : x_min; 
		}

		tmp_v.x = po_c.x;
		tmp_v.y = po_c.y+radius;
		if( this->IsDirectionArc(tmp_v) == 1 ){
			y_max = ( tmp_v.y > y_max ) ? tmp_v.y : y_max; 
		}

		tmp_v.x = po_c.x;
		tmp_v.y = po_c.y-radius;
		if( this->IsDirectionArc(tmp_v) == 1 ){
			y_min = ( tmp_v.y < y_min ) ? tmp_v.y : y_min; 
		}
	}
	else if( itype == 2 ){
		const unsigned int n = aRelCoMesh.size()/2;
		Com::CVector2D v0 = po_e-po_s;
		Com::CVector2D v1(-v0.y,v0.x);
		for(unsigned int i=0;i<n;i++){
			Com::CVector2D po0 = po_s + v0*aRelCoMesh[i*2+0] + v1*aRelCoMesh[i*2+1];
			x_min = ( po0.x < x_min ) ? po0.x : x_min; 
			x_max = ( po0.x > x_max ) ? po0.x : x_max; 
			y_min = ( po0.y < y_min ) ? po0.y : y_min; 
			y_max = ( po0.y > y_max ) ? po0.y : y_max; 
		}
	}
	else{ assert(0); abort(); }
}

// 辺と頂点を結ぶ直線で囲まれる面積を計算(右側にあれば＋)
double Cad::CEdge2D::AreaEdge() const
{
	if( this->itype == 0 ){
		return 0;
	}
	else if( this->itype == 1 ){
		Com::CVector2D pc;
		double radius;
		this->GetCenterRadius(pc,radius);
		Com::CVector2D v_sc = po_s - pc;
		Com::CVector2D ax = v_sc * (1/radius);
		Com::CVector2D ay;
		if( is_left_side ){  ay.x =  ax.y;  ay.y = -ax.x;  }
		else{                ay.x = -ax.y;  ay.y =  ax.x;  }
		Com::CVector2D v_ec = po_e - pc;
		double x=Com::Dot(v_ec,ax);
		double y=Com::Dot(v_ec,ay);
		double theta = atan2(y,x);
		const double PI = 3.14159265;
		if( theta < 0.0 ) theta += 2.0*PI;
		// 円弧と弦の間の領域の面積を求める
		double segment_area = theta*radius*radius*0.5;
		segment_area -= fabs( Com::TriArea(po_s,pc,po_e) );
		if( this->is_left_side ){ segment_area *= -1; }
		return segment_area;
	}
	else if( this->itype == 2 ){
		const unsigned int n = aRelCoMesh.size()/2;
		double area_ratio = 0;
		for(unsigned int i=0;i<n+1;i++){
			double divx, h0, h1;
			if( i== 0 ){
				divx = aRelCoMesh[0];
				h0 = 0;
				h1 = aRelCoMesh[1];
			}
			else if( i== n ){
				divx = 1-aRelCoMesh[n*2-2];
				h0 = aRelCoMesh[n*2-1];
				h1 = 0;
			}
			else{
				divx = aRelCoMesh[i*2] - aRelCoMesh[i*2-2];
				h0 = aRelCoMesh[i*2-1];
				h1 = aRelCoMesh[i*2+1];
			}
			area_ratio += 0.5*divx*(h0+h1);
		}
		const double sqlen = Com::SquareLength(po_e-po_s);
		double segment_area = -sqlen * area_ratio;
		return segment_area;
	}
	else{ assert(0); }
	return 0;
}

Com::CVector2D Cad::CEdge2D::GetTangentEdge(bool is_s) const {
	if( this->itype == 0 ){
		Com::CVector2D d = (is_s ) ? po_e-po_s : po_s-po_e;
		d.Normalize();
		return d;
	}
	else if( this->itype == 1 ){
		double radius;
		Com::CVector2D po_c;
		this->GetCenterRadius(po_c,radius);
		Com::CVector2D h = (is_s) ? po_s-po_c : po_e-po_c;
		Com::CVector2D v;
		if( (is_s && !is_left_side) || (!is_s && is_left_side ) ){ v.x = -h.y; v.y = h.x; }
		else{ v.x = h.y; v.y = -h.x; }
		v.Normalize();
		return v;
	}
	else if( this->itype == 2 ){
		const std::vector<double>& relcomsh = aRelCoMesh;
		const unsigned int ndiv = relcomsh.size()/2+1;
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		if( is_s ){
			Com::CVector2D d = h0*relcomsh[0] + v0*relcomsh[1];
			d.Normalize();
			return d;
		}
		else{
			Com::CVector2D d = po_s + h0*relcomsh[(ndiv-2)*2+0] + v0*relcomsh[(ndiv-2)*2+1] - po_e;
			d.Normalize();
			return d;
		}
	}
	assert(0);
	Com::CVector2D v(0,0);
	return v;
}

//! return the point on this edge that is nearest from (po_in)
Com::CVector2D Cad::CEdge2D::GetNearestPoint(const Com::CVector2D& po_in) const
{
	if( itype == 0 ){
		const double t = FindNearestPointParameter_Line_Point(po_in,po_s,po_e);
		if( t < 0 ){ return po_s; }
		if( t > 1 ){ return po_e; }
		else{ return po_s + (po_e-po_s)*t; }
	}
	else if( itype ==1 ){
		double radius;
		Com::CVector2D po_c;
		this->GetCenterRadius(po_c,radius);
		double len = Com::Distance(po_c,po_in);
		if( len < 1.0e-10 ){ return po_s; }
//		double dist = fabs(len-radius);
		const Com::CVector2D& po_p = po_c*(1-radius/len)+po_in*(radius/len);	// projected point
		if( this->IsDirectionArc(po_p) ){ return po_p; }
		const double dist_s = Com::SquareLength(po_in,po_s);
		const double dist_e = Com::SquareLength(po_in,po_e);
		if( dist_s < dist_e ){ return po_s; }
		else{ return po_e; }
	}
	else if( itype ==2 ){
		const std::vector<double>& relcomsh = aRelCoMesh;
		const unsigned int ndiv = relcomsh.size()/2+1;
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		double min_dist = Com::Distance(po_s,po_in);
		Com::CVector2D cand = po_s;
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
			Com::CVector2D poi0;
			if( idiv == 0 ){ poi0 = po_s; }
			else{ poi0 = po_s + h0*relcomsh[(idiv-1)*2+0] + v0*relcomsh[(idiv-1)*2+1]; }
			Com::CVector2D poi1;
			if( idiv == ndiv-1 ){ poi1 = po_e; }
			else{ poi1 = po_s + h0*relcomsh[idiv*2+0]     + v0*relcomsh[idiv*2+1]; }
			if( Com::Distance(po_in,poi1) < min_dist ){ 
				min_dist = Com::Distance(po_in,poi1);	
				cand = poi1; 
			}
			const double t = FindNearestPointParameter_Line_Point(po_in,poi0,poi1);
			if( t < 0.01 || t > 0.99 ) continue;
			Com::CVector2D po_mid = poi0 + (poi1-poi0)*t;
			if( Com::Distance(po_in,po_mid) < min_dist ){ 
				min_dist = Com::Distance(po_in,po_mid);	
				cand = po_mid;
			}
		}
		return cand;
	}
	assert(0);
	Com::CVector2D v(0,0);
	return v;
}

//! check self interaction inside edge
bool Cad::CEdge2D::IsCrossEdgeSelf() const
{
	if( this->itype == 0 || this->itype == 1 ){ return false; } // line segment and arc never intersect in itself
	if( this->itype == 2 ){	// Mesh
		const std::vector<double>& relcomsh = aRelCoMesh;
		const unsigned int ndiv = relcomsh.size()/2+1;
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
			Com::CVector2D poi0;
			if( idiv == 0 ){ poi0 = po_s; }
			else{ poi0 = po_s + h0*relcomsh[(idiv-1)*2+0] + v0*relcomsh[(idiv-1)*2+1]; }
			Com::CVector2D poi1;
			if( idiv == ndiv-1 ){ poi1 = po_e; }
			else{ poi1 = po_s + h0*relcomsh[idiv*2+0]     + v0*relcomsh[idiv*2+1]; }
			for(unsigned int jdiv=idiv+2;jdiv<ndiv;jdiv++){
				Com::CVector2D poj0;
				if( jdiv == 0 ){ poj0 = po_s; }
				else{ poj0 = po_s + h0*relcomsh[(jdiv-1)*2+0] + v0*relcomsh[(jdiv-1)*2+1]; }
				Com::CVector2D poj1;
				if( jdiv == ndiv-1 ){ poj1 = po_e; }
				else{ poj1 = po_s + h0*relcomsh[jdiv*2+0]     + v0*relcomsh[jdiv*2+1]; }
				if( IsCross_LineSeg_LineSeg(poi0,poi1, poj0,poj1) != 0 ){ return true; }
			}
		}
		return false;
	}
	else{ 
		assert(0);
	}
	return false;
}

double Cad::GetDist_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                                    const Com::CVector2D& po_s1, const Com::CVector2D& po_e1)
{
  if( IsCross_LineSeg_LineSeg(po_s0,po_e0, po_s1,po_e1) ) return -1;
  const double ds1 = GetDist_LineSeg_Point(po_s0,po_s1,po_e1);
  const double de1 = GetDist_LineSeg_Point(po_e0,po_s1,po_e1);
  const double ds0 = GetDist_LineSeg_Point(po_s1,po_s0,po_e0);
  const double de0 = GetDist_LineSeg_Point(po_e1,po_s0,po_e0);
  double min_dist = ds1;
  min_dist = ( de1 < min_dist ) ? de1 : min_dist;
  min_dist = ( ds0 < min_dist ) ? ds0 : min_dist;
  min_dist = ( de0 < min_dist ) ? de0 : min_dist;    
  return min_dist;  
}

Com::CVector2D Cad::GetProjectedPointOnCircle(const Com::CVector2D& c, double r, 
                                              const Com::CVector2D& v)
{
  Com::CVector2D cv = v-c;
  return (r/cv.Length())*cv + c;
}

// 円弧の中心からみて，点poと円弧が同じ方向に重なっているか？
static bool IsDirectionArc(const Com::CVector2D& po, 
                           const Com::CVector2D& po_s, const Com::CVector2D& po_e,
                           const Com::CVector2D& po_c, bool is_left_side)
{
	if( is_left_side ){
		if( Com::TriArea(po_s,po_c,po_e) > 0.0 ){
			if(    Com::TriArea(po_s,po_c,po) > 0.0 && Com::TriArea(po,po_c,po_e) > 0.0 ){ return true; }
			else{ return false; }
		}
		else{
			if(    Com::TriArea(po_s,po_c,po) > 0.0 || Com::TriArea(po,po_c,po_e) > 0.0 ){ return true; }
			else{ return false; }
		}
	}
	else{
		if( Com::TriArea(po_e,po_c,po_s) > 0.0 ){
			if(    Com::TriArea(po_e,po_c,po) > 0.0 && Com::TriArea(po,po_c,po_s) > 0.0 ){ return true; }
			else{ return false; }
		}
		else{
			if(    Com::TriArea(po_e,po_c,po) > 0.0 || Com::TriArea(po,po_c,po_s) > 0.0 ){ return true; }
			else{ return false; }
		}
	}
	return true;
}

static double GetDist_Point_Arc(const Com::CVector2D& po,
                                const Com::CVector2D& po_s1, const Com::CVector2D& po_e1,
                                const Com::CVector2D& po_c1, double radius1, bool is_left_side1)
{  
  double min_dist = Distance(po,po_s1);
  {
    double d0 = Distance(po,po_e1);
    min_dist = ( min_dist < d0 ) ? min_dist : d0;
  }
  if( IsDirectionArc( Cad::GetProjectedPointOnCircle(po_c1,radius1,po), po_s1,po_e1,po_c1,is_left_side1) ){
    double d0 = fabs(Distance(po,po_c1)-radius1);      
    min_dist = ( d0 < min_dist ) ? d0 : min_dist;
  }
  return min_dist;
}

static double GetDist_LineSeg_Arc(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,   
                                  const Com::CVector2D& po_s1, const Com::CVector2D& po_e1,
                                  const Com::CVector2D& po_c1, double radius1, bool is_left_side1)
{
  std::cout << "GetDist_Linseg_Arc " << std::endl;
  {
    double t0,t1;
    if( Cad::IsCross_Line_Circle(po_c1,radius1,  po_s0,po_e0, t0,t1) ){ 
      if( 0 < t0 && t0 < 1 && IsDirectionArc(po_s0 + (po_e0-po_s0)*t0, po_s1,po_e1,po_c1,is_left_side1) ){ return -1; }
      if( 0 < t1 && t1 < 1 && IsDirectionArc(po_s0 + (po_e0-po_s0)*t1, po_s1,po_e1,po_c1,is_left_side1) ){ return -1; }
    }
  }
  const double min_dist_s0 = GetDist_Point_Arc(po_s0, po_s1,po_e1,po_c1,radius1,is_left_side1);
  const double min_dist_e0 = GetDist_Point_Arc(po_e0, po_s1,po_e1,po_c1,radius1,is_left_side1);
  std::cout << "arc seg dist s0 e0 " << min_dist_s0 << " " << min_dist_e0 << std::endl;
  double min_dist = ( min_dist_s0 < min_dist_e0 ) ? min_dist_s0 : min_dist_e0;
  {
    const double t = Cad::FindNearestPointParameter_Line_Point(po_c1,po_s0,po_e0);
    if( t > 0 && t < 1){
      const Com::CVector2D& v = po_s0 + (po_e0-po_s0)*t;
      double d0 = Distance(v,po_c1)-radius1;
      std::cout << "arc seg shortese cand dist : " << d0 << std::endl;
      if( d0 > 0 ){
        if( IsDirectionArc( Cad::GetProjectedPointOnCircle(po_c1,radius1,v), po_s1,po_e1,po_c1,is_left_side1 ) ){
          min_dist = ( d0 < min_dist ) ? d0 : min_dist;
        }
      }        
    }      
  }    
  return min_dist;
}


// return minimum distance between this edge and e1
// if the edge obviouly intersects, this function return 0 or -1
double Cad::CEdge2D::Distance(const CEdge2D& e1) const
{
	const Com::CVector2D& po_s1 = e1.po_s;
	const Com::CVector2D& po_e1 = e1.po_e;
	if( this->itype == 0 && e1.itype == 0 ){	// intersection between lines
    return GetDist_LineSeg_LineSeg(po_s,po_e,po_s1,po_e1);
	}
	else if( this->itype == 0 && e1.itype == 1 ){
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
    return GetDist_LineSeg_Arc(po_s,po_e, e1.po_s,e1.po_e,po_c1,radius1,e1.is_left_side);
	}
	else if( this->itype == 1 && e1.itype == 0 ){
		return e1.Distance(*this);
	}
	else if( this->itype == 1 && e1.itype == 1 ){
		Com::CVector2D po_c0; double radius0;  this->GetCenterRadius(po_c0,radius0);
		Com::CVector2D po_c1; double radius1;  e1.GetCenterRadius(po_c1,radius1);
		////////////////
		Com::CVector2D po0,po1;
    bool is_cross_circle01 = false;
		if( IsCross_Circle_Circle(po_c0,radius0,  po_c1,radius1, po0,po1) ){ 
      is_cross_circle01 = true;
      if( this->IsDirectionArc(po0) && e1.IsDirectionArc(po0) ){ return -1; }
      if( this->IsDirectionArc(po1) && e1.IsDirectionArc(po1) ){ return -1; }      
    }
    const double min_dist_s0 = GetDist_Point_Arc(this->po_s, e1.po_s,e1.po_e,po_c1,radius1,e1.is_left_side);
    const double min_dist_e0 = GetDist_Point_Arc(this->po_e, e1.po_s,e1.po_e,po_c1,radius1,e1.is_left_side);
    const double min_dist_s1 = GetDist_Point_Arc(e1.po_s, this->po_s,this->po_e,po_c0,radius0,this->is_left_side);
    const double min_dist_e1 = GetDist_Point_Arc(e1.po_e, this->po_s,this->po_e,po_c0,radius0,this->is_left_side);
    const double min_dist0 = ( min_dist_s0 < min_dist_e0 ) ? min_dist_s0 : min_dist_e0;
    const double min_dist1 = ( min_dist_s1 < min_dist_e1 ) ? min_dist_s1 : min_dist_e1;    
    double min_dist = ( min_dist0 < min_dist1 ) ? min_dist0 : min_dist1;
    if( !is_cross_circle01 ){
      const bool is_c0_inside_c1 = Com::Distance(po_c0,po_c1) < radius1;
      const bool is_c1_inside_c0 = Com::Distance(po_c1,po_c0) < radius0;
      if( !is_c0_inside_c1 && !is_c1_inside_c0 ){
        const Com::CVector2D& v1 = Cad::GetProjectedPointOnCircle(po_c1,radius1, po_c0);
        const Com::CVector2D& v0 = Cad::GetProjectedPointOnCircle(po_c0,radius0, po_c1);      
        if( e1.IsDirectionArc( v1 ) && this->IsDirectionArc( v0 ) ){        
          double d0 = Com::Distance(v0,v1);
          min_dist = ( d0 < min_dist ) ? d0 : min_dist;
        }
      }
      else{
        if( radius0 < radius1 ){
          const Com::CVector2D& v1 = Cad::GetProjectedPointOnCircle(po_c1,radius1, po_c0);        
          const Com::CVector2D& v0 = Cad::GetProjectedPointOnCircle(po_c0,radius0, v1);               
          if( e1.IsDirectionArc( v1 ) && this->IsDirectionArc( v0 ) ){        
            double d0 = Com::Distance(v0,v1);
            min_dist = ( d0 < min_dist ) ? d0 : min_dist;
          }        
        }
        else{
          assert( is_c1_inside_c0 );
          const Com::CVector2D& v0 = Cad::GetProjectedPointOnCircle(po_c0,radius0, po_c1);              
          const Com::CVector2D& v1 = Cad::GetProjectedPointOnCircle(po_c1,radius1, v0);
          if( e1.IsDirectionArc( v1 ) && this->IsDirectionArc( v0 ) ){        
            double d0 = Com::Distance(v0,v1);
            min_dist = ( d0 < min_dist ) ? d0 : min_dist;
          }
        }
      }
    }            
		return min_dist;
	}
	else if( this->itype == 0 && e1.itype == 2 ){
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
    double min_dist = -1;
		for(unsigned int idiv=0;idiv<ndiv1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
      const double dist = GetDist_LineSeg_LineSeg(po_s,po_e,po0,po1);
      if( dist < -0.5 ) return dist;
      if( dist < min_dist || min_dist < -0.5 ){ min_dist = dist; }
		}
		return min_dist;
	}
	else if( this->itype == 2 && e1.itype == 0 ){
		return e1.Distance(*this);
	}
	else if( this->itype == 2 && e1.itype == 1 ){
		return e1.Distance(*this);
	}
	else if( this->itype == 1 && e1.itype == 2 ){
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
    double min_dist = -1;    
		for(unsigned int idiv=0;idiv<ndiv1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
      const double dist = GetDist_LineSeg_Arc(po0,po1, po_s,po_e,po_c0,radius0,is_left_side);
      if( dist < -0.5 ) return dist;
      if( dist < min_dist || min_dist < -0.5 ){ min_dist = dist; }      
		}
		return min_dist;
	}
	else if( this->itype == 2 && e1.itype == 2 ){
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh0 = this->aRelCoMesh;
		const unsigned int ndiv0 = relcomsh0.size()/2+1;
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		////////////////
    double min_dist = -1;GetDist_LineSeg_LineSeg(po_s,po_e,po_s1,po_e1);
		for(unsigned int idiv=0;idiv<ndiv0;idiv++){
			Com::CVector2D po0_i;
			if( idiv == 0 ){ po0_i = po_s; }
			else{ po0_i = po_s + h0*aRelCoMesh[(idiv-1)*2+0] + v0*aRelCoMesh[(idiv-1)*2+1]; }
			Com::CVector2D po1_i;
			if( idiv == ndiv0-1 ){ po1_i = po_e; }
			else{ po1_i = po_s + h0*aRelCoMesh[idiv*2+0]     + v0*aRelCoMesh[idiv*2+1]; }
			for(unsigned int jdiv=0;jdiv<ndiv1;jdiv++){
				Com::CVector2D po0_j;
				if( jdiv == 0 ){ po0_j = e1.po_s; }
				else{ po0_j = e1.po_s + h1*relcomsh1[(jdiv-1)*2+0] + v1*relcomsh1[(jdiv-1)*2+1]; }
				Com::CVector2D po1_j;
				if( jdiv == ndiv1-1 ){ po1_j = e1.po_e; }
				else{ po1_j = e1.po_s + h1*relcomsh1[jdiv*2+0]     + v1*relcomsh1[jdiv*2+1]; }
        const double dist = GetDist_LineSeg_LineSeg(po0_i,po1_i,po0_j,po1_j);
        if( dist < -0.5 ) return -1;
        if( min_dist < -0.5 || dist < min_dist ){ min_dist = dist; }        
			}
		}
		return min_dist;
	}
	return 1;
}

bool Cad::CEdge2D::IsCrossEdge(const CEdge2D& e1) const
{
	const Com::CVector2D& po_s1 = e1.po_s;
	const Com::CVector2D& po_e1 = e1.po_e;
	if( this->itype == 0 && e1.itype == 0 ){	// intersection between lines
		return ( IsCross_LineSeg_LineSeg(po_s,po_e, po_s1,po_e1) != 0 );
	}
	else if( this->itype == 0 && e1.itype == 1 ){
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
		double t0,t1;
		if( !IsCross_Line_Circle(po_c1,radius1,  po_s,po_e, t0,t1) ){ return false; }
		if( 0 < t0 && t0 < 1 && e1.IsDirectionArc(po_s + (po_e-po_s)*t0) == 1 ){ return true; }
		if( 0 < t1 && t1 < 1 && e1.IsDirectionArc(po_s + (po_e-po_s)*t1) == 1 ){ return true; }
		return false;
	}
	else if( this->itype == 1 && e1.itype == 0 ){
		return e1.IsCrossEdge(*this);
	}
	else if( this->itype == 1 && e1.itype == 1 ){
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
		////////////////
		Com::CVector2D po0,po1;
		if( !IsCross_Circle_Circle(po_c0,radius0,  po_c1,radius1, po0,po1) ){ return false; }
		if( this->IsDirectionArc(po0) && e1.IsDirectionArc(po0) ){ return true; }
		if( this->IsDirectionArc(po1) && e1.IsDirectionArc(po1) ){ return true; }
		return false;
	}
	else if( this->itype == 0 && e1.itype == 2 ){
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		for(unsigned int idiv=0;idiv<ndiv1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
			if( IsCross_LineSeg_LineSeg(po_s,po_e, po0,po1) != 0 ){ return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 0 ){
		return e1.IsCrossEdge(*this);
	}
	else if( this->itype == 2 && e1.itype == 1 ){
		return e1.IsCrossEdge(*this);
	}
	else if( this->itype == 1 && e1.itype == 2 ){
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		for(unsigned int idiv=0;idiv<ndiv1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
			double t0,t1;
			if( !IsCross_Line_Circle(po_c0,radius0,  po0,po1, t0,t1) ){ continue; }
			if( 0 < t0 && t0 < 1 && this->IsDirectionArc(po0 + (po1-po0)*t0) == 1 ){ return true; }
			if( 0 < t1 && t1 < 1 && this->IsDirectionArc(po0 + (po1-po0)*t1) == 1 ){ return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 2 ){
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh0 = this->aRelCoMesh;
		const unsigned int ndiv0 = relcomsh0.size()/2+1;
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		////////////////
		for(unsigned int idiv=0;idiv<ndiv0;idiv++){
			Com::CVector2D po0_i;
			if( idiv == 0 ){ po0_i = po_s; }
			else{ po0_i = po_s + h0*aRelCoMesh[(idiv-1)*2+0] + v0*aRelCoMesh[(idiv-1)*2+1]; }
			Com::CVector2D po1_i;
			if( idiv == ndiv0-1 ){ po1_i = po_e; }
			else{ po1_i = po_s + h0*aRelCoMesh[idiv*2+0]     + v0*aRelCoMesh[idiv*2+1]; }
			for(unsigned int jdiv=0;jdiv<ndiv1;jdiv++){
				Com::CVector2D po0_j;
				if( jdiv == 0 ){ po0_j = e1.po_s; }
				else{ po0_j = e1.po_s + h1*relcomsh1[(jdiv-1)*2+0] + v1*relcomsh1[(jdiv-1)*2+1]; }
				Com::CVector2D po1_j;
				if( jdiv == ndiv1-1 ){ po1_j = e1.po_e; }
				else{ po1_j = e1.po_s + h1*relcomsh1[jdiv*2+0]     + v1*relcomsh1[jdiv*2+1]; }
				if( IsCross_LineSeg_LineSeg(po0_i,po1_i, po0_j,po1_j) != 0 ){ return true; }
			}
		}
		return false;
	}
	return true;
}

bool Cad::CEdge2D::IsCrossEdge_ShareOnePoint(const CEdge2D& e1, bool is_share_s0, bool is_share_s1) const
{
	const Com::CVector2D& po_s1 = e1.po_s;
	const Com::CVector2D& po_e1 = e1.po_e;
  if(  is_share_s0 &&  is_share_s1 ){ assert( Com::SquareLength(po_s,po_s1) < 1.0e-20 ); }
	if(  is_share_s0 && !is_share_s1 ){ assert( Com::SquareLength(po_s,po_e1) < 1.0e-20 ); }
	if( !is_share_s0 &&  is_share_s1 ){ assert( Com::SquareLength(po_e,po_s1) < 1.0e-20 ); }
	if( !is_share_s0 && !is_share_s1 ){ assert( Com::SquareLength(po_e,po_e1) < 1.0e-20 ); }
	////////////////////////////////
	if( this->itype == 0 && e1.itype == 0 ){	// line-line intersection
		return false;
	}
	else if( this->itype == 0 && e1.itype == 1 ){	// line-arc intersection
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
		////////////////
		double t0,t1;
		if( !IsCross_Line_Circle(po_c1,radius1,  po_s,po_e, t0,t1) ) return false;
		const Com::CVector2D p0 = po_s + (po_e-po_s)*t0;
		const Com::CVector2D p1 = po_s + (po_e-po_s)*t1;
		// この後t0,t1は共有点との距離計算に使われる
		if( !is_share_s0 ){ t0=1-t0; t1=1-t1; }
		if( fabs(t0) < fabs(t1) && 0 < t1 && t1 < 1 && e1.IsDirectionArc(p1) ){ assert(fabs(t0)<1.0e-5); return true; }
		if( fabs(t0) > fabs(t1) && 0 < t0 && t0 < 1 && e1.IsDirectionArc(p0) ){ assert(fabs(t1)<1.0e-5); return true; }
		return false;
	}
	else if( this->itype == 1 && e1.itype == 0 ){
		return e1.IsCrossEdge_ShareOnePoint(*this,is_share_s1,is_share_s0);
	}
	else if( this->itype == 1 && e1.itype == 1 ){	// 円弧同士の交差
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
		////////////////
		Com::CVector2D po0,po1;
		bool is_cross = IsCross_Circle_Circle(po_c0,radius0,  po_c1,radius1, po0,po1);
		if( !is_cross ) return false;
		////////////////
		Com::CVector2D po_share;
		if( is_share_s0 ){ po_share = po_s; }
		else{ po_share = po_e; }
		const double sqdist0 = Com::SquareLength(po_share,po0);
		const double sqdist1 = Com::SquareLength(po_share,po1);
		if( sqdist0 < sqdist1 && this->IsDirectionArc(po1) && e1.IsDirectionArc(po1) ){ assert(sqdist0<1.0e-20); return true; }
		if( sqdist0 > sqdist1 && this->IsDirectionArc(po0) && e1.IsDirectionArc(po0) ){ assert(sqdist1<1.0e-20); return true; }
		return false;
	}
	else if( this->itype == 0 && e1.itype == 2 ){	// 直線と折れ線の交差
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		Com::CVector2D v0 = e1.po_e-e1.po_s;
		Com::CVector2D v1(-v0.y,v0.x);
		const unsigned int idiv_s = ( is_share_s1 ) ? 1 : 0;
		const unsigned int idiv_e = ( is_share_s1 ) ? ndiv1 : ndiv1-1;
		for(unsigned int idiv=idiv_s;idiv<idiv_e;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + v0*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + v0*relcomsh1[idiv*2+0] + v1*relcomsh1[idiv*2+1]; }
			if( IsCross_LineSeg_LineSeg(po_s,po_e, po0,po1) != 0 ){ return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 0 ){
		return e1.IsCrossEdge_ShareOnePoint(*this,is_share_s1,is_share_s0);
	}
	else if( this->itype == 1 && e1.itype == 2 ){	// 円弧と折れ線の交差
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		const unsigned int idiv_s = ( is_share_s1 ) ? 1 : 0;
		const unsigned int idiv_e = ( is_share_s1 ) ? ndiv1 : ndiv1-1;
		for(unsigned int idiv=idiv_s;idiv<idiv_e;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
			double t0,t1;
			if( !IsCross_Line_Circle(po_c0,radius0,  po0,po1, t0,t1) ){ continue; }
			if( 0 < t0 && t0 < 1 && this->IsDirectionArc( po0+(po1-po0)*t0 ) == 1 ){ return true; }
			if( 0 < t1 && t1 < 1 && this->IsDirectionArc( po0+(po1-po0)*t1 ) == 1 ){ return true; }
		}
		{
			Com::CVector2D po0, po1;
			if( is_share_s1 ){	po0 = e1.po_s;	po1 = e1.po_s + h1*relcomsh1[0        ] + v1*relcomsh1[1        ]; }
			else{				po1 = e1.po_e;	po0 = e1.po_s + h1*relcomsh1[ndiv1*2-4] + v1*relcomsh1[ndiv1*2-3]; }	
			////////////////
			double t0,t1;
			if( !IsCross_Line_Circle(po_c0,radius0,  po0,po1, t0,t1) ){ return false; }
			const Com::CVector2D r0 = po0 + (po1-po0)*t0;
			const Com::CVector2D r1 = po0 + (po1-po0)*t1;
			// この後t0,t1は共有点との距離計算に使われる
			if( !is_share_s1 ){ t0=1-t0; t1=1-t1; }
			if( fabs(t0) < fabs(t1) && 0 < t1 && t1 < 1 && this->IsDirectionArc(r1) ){ assert(fabs(t0)<1.0e-5); return true; }
			if( fabs(t0) > fabs(t1) && 0 < t0 && t0 < 1 && this->IsDirectionArc(r0) ){ assert(fabs(t1)<1.0e-5); return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 1 ){
		return e1.IsCrossEdge_ShareOnePoint(*this,is_share_s1,is_share_s0);
	}
	else if( this->itype == 2 && e1.itype == 2 ){	// 折れ線同士の交差
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh0 = this->aRelCoMesh;
		const unsigned int ndiv0 = relcomsh0.size()/2+1;
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		const unsigned int idiv0_exc = ( is_share_s0 ) ? 0 : ndiv0-1;
		const unsigned int idiv1_exc = ( is_share_s1 ) ? 0 : ndiv1-1;
		////////////////
		for(unsigned int idiv=0;idiv<ndiv0;idiv++){
			Com::CVector2D po0_i;
			if( idiv == 0 ){ po0_i = po_s; }
			else{ po0_i = po_s + h0*aRelCoMesh[(idiv-1)*2+0] + v0*aRelCoMesh[(idiv-1)*2+1]; }
			Com::CVector2D po1_i;
			if( idiv == ndiv0-1 ){ po1_i = po_e; }
			else{ po1_i = po_s + h0*aRelCoMesh[idiv*2+0]     + v0*aRelCoMesh[idiv*2+1]; }
			for(unsigned int jdiv=0;jdiv<ndiv1;jdiv++){
				if( idiv == idiv0_exc && jdiv == idiv1_exc ) continue;
				Com::CVector2D po0_j;
				if( jdiv == 0 ){ po0_j = e1.po_s; }
				else{ po0_j = e1.po_s + h1*relcomsh1[(jdiv-1)*2+0] + v1*relcomsh1[(jdiv-1)*2+1]; }
				Com::CVector2D po1_j;
				if( jdiv == ndiv1-1 ){ po1_j = e1.po_e; }
				else{ po1_j = e1.po_s + h1*relcomsh1[jdiv*2+0]     + v1*relcomsh1[jdiv*2+1]; }
				if( IsCross_LineSeg_LineSeg(po0_i,po1_i, po0_j,po1_j) != 0 ){ return true; }
			}
		}
		return false;
	}
	else{
		assert(0);
	}
	return false;
}

//! 一端が共有された辺同士の交差判定
bool Cad::CEdge2D::IsCrossEdge_ShareBothPoints(const CEdge2D& e1, bool is_share_s0s1) const
{
	const Com::CVector2D& po_s1 = e1.po_s;
	const Com::CVector2D& po_e1 = e1.po_e;
	if(  is_share_s0s1 ){ assert( Com::SquareLength(po_s,po_s1) < 1.0e-20 && Com::SquareLength(po_e,po_e1) < 1.0e-20); }
	else{                 assert( Com::SquareLength(po_s,po_e1) < 1.0e-20 && Com::SquareLength(po_e,po_s1) < 1.0e-20); }
	////////////////////////////////
	if( this->itype == 0 && e1.itype == 0 ){	// intersection between lines
		return true;
	}
	else if( this->itype == 0 && e1.itype == 1 ){	// intersection between line and arc
		return false;
	}
	else if( this->itype == 1 && e1.itype == 0 ){
		return e1.IsCrossEdge_ShareBothPoints(*this,is_share_s0s1);
	}
	else if( this->itype == 1 && e1.itype == 1 ){	// 円弧同士の交差
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		Com::CVector2D po_c1;
		double radius1;
		e1.GetCenterRadius(po_c1,radius1);
		if( Com::SquareLength(po_c0-po_c1) < 1.0e-10 ) return true;
		return false;
	}
	else if( this->itype == 0 && e1.itype == 2 ){	// 直線と折れ線の交差
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		Com::CVector2D v0 = e1.po_e-e1.po_s;
		Com::CVector2D v1(-v0.y,v0.x);
		for(unsigned int idiv=1;idiv<ndiv1-1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + v0*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + v0*relcomsh1[idiv*2+0] + v1*relcomsh1[idiv*2+1]; }
			if( IsCross_LineSeg_LineSeg(po_s,po_e, po0,po1) != 0 ){ return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 0 ){
		return e1.IsCrossEdge_ShareBothPoints(*this,is_share_s0s1);
	}
	else if( this->itype == 1 && e1.itype == 2 ){	// 円弧と折れ線の交差
		Com::CVector2D po_c0;
		double radius0;
		this->GetCenterRadius(po_c0,radius0);
		////////////////
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		for(unsigned int idiv=1;idiv<ndiv1-1;idiv++){
			Com::CVector2D po0;
			if( idiv == 0 ){ po0 = e1.po_s; }
			else{ po0 = e1.po_s + h1*relcomsh1[(idiv-1)*2+0] + v1*relcomsh1[(idiv-1)*2+1]; }
			Com::CVector2D po1;
			if( idiv == ndiv1-1 ){ po1 = e1.po_e; }
			else{ po1 = e1.po_s + h1*relcomsh1[idiv*2+0]     + v1*relcomsh1[idiv*2+1]; }
			////////////////
			double t0,t1;
			if( !IsCross_Line_Circle(po_c0,radius0,  po0,po1, t0,t1) ){ continue; }
			if( 0 < t0 && t0 < 1 && this->IsDirectionArc( po0+(po1-po0)*t0 ) == 1 ){ return true; }
			if( 0 < t1 && t1 < 1 && this->IsDirectionArc( po0+(po1-po0)*t1 ) == 1 ){ return true; }
		}
		return false;
	}
	else if( this->itype == 2 && e1.itype == 1 ){
		return e1.IsCrossEdge_ShareBothPoints(*this,is_share_s0s1);
	}
	else if( this->itype == 2 && e1.itype == 2 ){	// 折れ線同士の交差
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		const Com::CVector2D& h1 = e1.po_e-e1.po_s;
		const Com::CVector2D v1(-h1.y,h1.x);
		////////////////
		const std::vector<double>& relcomsh0 = this->aRelCoMesh;
		const unsigned int ndiv0 = relcomsh0.size()/2+1;
		const std::vector<double>& relcomsh1 = e1.aRelCoMesh;
		const unsigned int ndiv1 = relcomsh1.size()/2+1;
		////////////////
		for(unsigned int idiv=1;idiv<ndiv0-1;idiv++){
			Com::CVector2D po0_i;
			if( idiv == 0 ){ po0_i = po_s; }
			else{ po0_i = po_s + h0*aRelCoMesh[(idiv-1)*2+0] + v0*aRelCoMesh[(idiv-1)*2+1]; }
			Com::CVector2D po1_i;
			if( idiv == ndiv0-1 ){ po1_i = po_e; }
			else{ po1_i = po_s + h0*aRelCoMesh[idiv*2+0]     + v0*aRelCoMesh[idiv*2+1]; }
			for(unsigned int jdiv=1;jdiv<ndiv1-1;jdiv++){
				Com::CVector2D po0_j;
				if( jdiv == 0 ){ po0_j = e1.po_s; }
				else{ po0_j = e1.po_s + h1*relcomsh1[(jdiv-1)*2+0] + v1*relcomsh1[(jdiv-1)*2+1]; }
				Com::CVector2D po1_j;
				if( jdiv == ndiv1-1 ){ po1_j = e1.po_e; }
				else{ po1_j = e1.po_s + h1*relcomsh1[jdiv*2+0]     + v1*relcomsh1[jdiv*2+1]; }
				if( IsCross_LineSeg_LineSeg(po0_i,po1_i, po0_j,po1_j) != 0 ){ return true; }
			}
		}
		return false;
	}
	else{
		assert(0);
	}
	return false;
}

bool Cad::CEdge2D::GetNearestIntersectionPoint_AgainstHalfLine(Com::CVector2D& sec, const Com::CVector2D& org, const Com::CVector2D& dir) const
{
	const Com::CVector2D dir1 = dir*( 1.0/sqrt( Com::SquareLength(dir) ) );
	double lenlong = 0;
	{	// BoundingBoxとの干渉チェック
		double min_x,max_x, min_y,max_y;
		this->GetBoundingBox(min_x,max_x, min_y,max_y);
		Com::CVector2D po_d = org + dir1;
		double area1 = Com::TriArea(org,po_d, Com::CVector2D(min_x,min_y) );
		double area2 = Com::TriArea(org,po_d, Com::CVector2D(min_x,max_y) );
		double area3 = Com::TriArea(org,po_d, Com::CVector2D(max_x,min_y) );
		double area4 = Com::TriArea(org,po_d, Com::CVector2D(max_x,max_y) );
		if( area1<0 && area2<0 && area3<0 && area4<0 ) return false;
		if( area1>0 && area2>0 && area3>0 && area4>0 ) return false;
		const double len0 = Com::Distance(org, Com::CVector2D(0.5*(min_x+max_x),0.5*(min_y+max_y)) );
		const double len1 = max_x - min_x;
		const double len2 = max_y - min_y;
		lenlong = 2*(len0+len1+len2);
	}
	const Com::CVector2D end = org + dir1*lenlong;
  if(      itype == 0 ){
    const double area1 = Com::TriArea(po_s,po_e,org);
    const double area2 = Com::TriArea(po_s,po_e,end);    
    if( (area1>0) == (area2>0) ) return false;
    const double area3 = Com::TriArea(org,end,po_s);
    const double area4 = Com::TriArea(org,end,po_e); 
    if( (area3>0) == (area4>0) ) return false;
    double t = area1/(area1-area2);
    sec = org + (end-org)*t;
    return true;
	}
	else if( itype == 1 ){    
		Com::CVector2D po_c1;
		double radius1;
		this->GetCenterRadius(po_c1,radius1);
		////////////////
		double t0,t1;
		if( !IsCross_Line_Circle(po_c1,radius1,  org,end, t0,t1) ) return false;
		const Com::CVector2D p0 = org + (end-org)*t0;
		const Com::CVector2D p1 = org + (end-org)*t1;
    const bool is_sec0 = 0 < t0 && t0 < 1 && this->IsDirectionArc(p0);
    const bool is_sec1 = 0 < t1 && t1 < 1 && this->IsDirectionArc(p1);
    if(  is_sec0 && !is_sec1 ){ sec=p0; return true; }
    if( !is_sec0 &&  is_sec1 ){ sec=p1; return true; }
    if(  is_sec0 &&  is_sec1 ){ sec = ( t0 < t1 ) ? p0 : p1; return true; }
    return false;
	}
	else if( itype == 2 ){
		const std::vector<double>& relcomsh = aRelCoMesh;
		const unsigned int ndiv = relcomsh.size()/2+1;
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		unsigned int icnt = 0;
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
			Com::CVector2D poi0;
			if( idiv == 0 ){ poi0 = po_s; }
			else{ poi0 = po_s + h0*relcomsh[(idiv-1)*2+0] + v0*relcomsh[(idiv-1)*2+1]; }
			Com::CVector2D poi1;
			if( idiv == ndiv-1 ){ poi1 = po_e; }
			else{ poi1 = po_s + h0*relcomsh[idiv*2+0]     + v0*relcomsh[idiv*2+1]; }
      const double area1 = Com::TriArea(poi0,poi1,org);
      const double area2 = Com::TriArea(poi0,poi1,end);    
      if( (area1>0) == (area2>0) ) continue;
      const double area3 = Com::TriArea(org,end,poi0);
      const double area4 = Com::TriArea(org,end,poi1); 
      if( (area3>0) == (area4>0) ) continue;
      double t = area1/(area1-area2);
      sec = org + (end-org)*t;      
      return false;
		}
		return true;
	}
	assert(0);
	return false;  
}


// 辺と半直線の交差回数を得る
// 領域の内側or外側判定に使われる
int Cad::CEdge2D::NumIntersect_AgainstHalfLine(const Com::CVector2D& po_b, const Com::CVector2D& dir0) const 
{
	const Com::CVector2D dir1 = dir0*( 1.0/sqrt( Com::SquareLength(dir0) ) );
	double lenlong = 0;
	{	// BoundingBoxとの干渉チェック
		double min_x,max_x, min_y,max_y;
		this->GetBoundingBox(min_x,max_x, min_y,max_y);
		Com::CVector2D po_d = po_b + dir1;
		double area1 = Com::TriArea(po_b,po_d, Com::CVector2D(min_x,min_y) );
		double area2 = Com::TriArea(po_b,po_d, Com::CVector2D(min_x,max_y) );
		double area3 = Com::TriArea(po_b,po_d, Com::CVector2D(max_x,min_y) );
		double area4 = Com::TriArea(po_b,po_d, Com::CVector2D(max_x,max_y) );
		if( area1<0 && area2<0 && area3<0 && area4<0 ) return 0;
		if( area1>0 && area2>0 && area3>0 && area4>0 ) return 0;
		const double len0 = Com::Distance(po_b, Com::CVector2D(0.5*(min_x+max_x),0.5*(min_y+max_y)) );
		const double len1 = max_x - min_x;
		const double len2 = max_y - min_y;
		lenlong = 2*(len0+len1+len2);
	}
	const Com::CVector2D po_d = po_b + dir1*lenlong;
  if(      itype == 0 ){
    return IsCross_LineSeg_LineSeg(po_s,po_e,po_b,po_d);
	}
	else if( itype == 1 ){
		return this->NumCross_Arc_LineSeg(po_b,po_d);
	}
	else if( itype == 2 ){
		const std::vector<double>& relcomsh = aRelCoMesh;
		const unsigned int ndiv = relcomsh.size()/2+1;
		const Com::CVector2D& h0 = po_e-po_s;
		const Com::CVector2D v0(-h0.y,h0.x);
		unsigned int icnt = 0;
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
			Com::CVector2D poi0;
			if( idiv == 0 ){ poi0 = po_s; }
			else{ poi0 = po_s + h0*relcomsh[(idiv-1)*2+0] + v0*relcomsh[(idiv-1)*2+1]; }
			Com::CVector2D poi1;
			if( idiv == ndiv-1 ){ poi1 = po_e; }
			else{ poi1 = po_s + h0*relcomsh[idiv*2+0]     + v0*relcomsh[idiv*2+1]; }
			int res = IsCross_LineSeg_LineSeg(po_b,po_d,poi0,poi1);
			if( res == -1 ) return -1;
			icnt += res;
		}
		return icnt;
	}
	assert(0);
	return 0;
}

double Cad::CEdge2D::GetCurveLength() const
{
    if(      this->itype == 0 )
    {
        return sqrt( (po_s.x-po_e.x)*(po_s.x-po_e.x)+(po_s.y-po_e.y)*(po_s.y-po_e.y) );
    }
    else if( this->itype == 1 )
    {
        double radius, theta;
        Com::CVector2D pc, lx,ly;
        this->GetCenterRadiusThetaLXY(pc,radius, theta,lx,ly);
        return radius*theta;
    }
    else if( this->itype == 2 )
    {
        const unsigned int npo_cad = this->aRelCoMesh.size()/2;
        double lenrel_tot = 0;
        for(unsigned int idiv=0;idiv<npo_cad+1;idiv++){
            double x0,y0;
            if( idiv == 0 ){ x0=0; y0=0; }
            else{ x0=aRelCoMesh[idiv*2-2]; y0=aRelCoMesh[idiv*2-1]; }
            double x1,y1;
            if( idiv == npo_cad ){ x1=1; y1=0; }
            else{ x1=aRelCoMesh[idiv*2+0]; y1=aRelCoMesh[idiv*2+1]; }
            lenrel_tot += sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
        }
        const double edge_len = sqrt( (po_s.x-po_e.x)*(po_s.x-po_e.x)+(po_s.y-po_e.y)*(po_s.y-po_e.y) );
        return lenrel_tot*edge_len;
    }
    return 0;
}

// 円の中心と半径を計算する
bool Cad::CEdge2D::GetCenterRadius(Com::CVector2D& po_c, double& radius) const
{
	if( itype != 1 ) return false;

	const double len_edge = Com::Distance(po_s,po_e);
	{
		Com::CVector2D ph((po_s.x+po_e.x)*0.5, (po_s.y+po_e.y)*0.5);
		Com::CVector2D vv(po_s.y-ph.y, ph.x-po_s.x);
		double len_vv = sqrt( vv.x*vv.x + vv.y*vv.y );
		vv.x /= len_vv;
		vv.y /= len_vv;
		po_c.x = ph.x + vv.x*dist*len_edge;
		po_c.y = ph.y + vv.y*dist*len_edge;
	}
	{
		Com::CVector2D vcs(po_s.x-po_c.x,po_s.y-po_c.y);
		radius = sqrt( vcs.x*vcs.x + vcs.y*vcs.y );
	}
	return true;
}

bool Cad::CEdge2D::GetCenterRadiusThetaLXY(Com::CVector2D& pc, double& radius,
        double& theta, Com::CVector2D& lx, Com::CVector2D& ly) const
{
    if( itype != 1 ) return false;
    this->GetCenterRadius(pc,radius);
    Com::CVector2D prs(po_s.x-pc.x,po_s.y-pc.y);
    lx.x = prs.x / radius;
    lx.y = prs.y / radius;
    if( is_left_side ){  ly.x =  lx.y;  ly.y = -lx.x;  }
    else{                ly.x = -lx.y;  ly.y =  lx.x;  }
    Com::CVector2D pre(po_e.x-pc.x,po_e.y-pc.y);
    assert( fabs( Com::SquareLength(prs)-Com::SquareLength(pre) ) < 1.0e-10 );
    double x=Com::Dot(pre,lx);
    double y=Com::Dot(pre,ly);
    theta = atan2(y,x);
    const double PI = 3.14159265;
    if( theta < 0.0 ) theta += 2.0*PI;
    return true;
}

bool Cad::CEdge2D::GetCurve_Mesh(std::vector<Com::CVector2D>& aCo, int ndiv) const
{
    aCo.clear();
    if( ndiv <= 0 ){    // 出来るだけ点多くして忠実な曲線の表現
        if(      this->itype == 0 ) // 直線の場合
        {
            return true;
        }
        else if( this->itype == 1 ) // 円弧の場合
        {
            double radius, theta;
            Com::CVector2D pc, lx,ly;
            this->GetCenterRadiusThetaLXY(pc,radius, theta,lx,ly);
            const unsigned int ndiv = (unsigned int)( theta*360/(6.0*6.28) ) + 1; // 6度ごとの分割
            // 再帰呼び出しのときはndiv=0は絶対入れちゃだめ(無限ループ)
            return this->GetCurve_Mesh(aCo,ndiv);
        }
        else if( this->itype == 2 ) // メッシュの場合
        {
            const Com::CVector2D& v0 = po_e - po_s;
            const Com::CVector2D v1(-v0.y, v0.x);
            ////////////////
            const unsigned int npo = aRelCoMesh.size()/2;
            for(unsigned int i=0;i<npo;i++){
               Com::CVector2D vec0 = po_s + v0*aRelCoMesh[i*2+0] + v1*aRelCoMesh[i*2+1];
               aCo.push_back(vec0);
            }
            return true;
        }
    }
    else{    // elenの長さでメッシュを切る
        aCo.reserve(ndiv);
        if(      this->itype == 0 ) // 直線の場合
        {
            Com::CVector2D div_t = ( po_e - po_s )*(1.0/ndiv);
            for(unsigned int idiv=1;idiv<(unsigned int)ndiv;idiv++){
                Com::CVector2D vec0 = po_s + div_t*idiv;
                aCo.push_back(vec0);
            }
            return true;
        }
        else if( this->itype == 1 ) // 円弧の場合
        {
            double radius, theta;
            Com::CVector2D pc, lx,ly;
            this->GetCenterRadiusThetaLXY(pc,radius, theta,lx,ly);
            const double div_theta = theta/ndiv;
//            std::cout << ndiv << std::endl;
            for(unsigned int i=1;i<(unsigned int)ndiv;i++){
                const double cur_theta = i*div_theta;
                Com::CVector2D vec0 = ( cos(cur_theta)*lx + sin(cur_theta)*ly )*radius + pc;
                aCo.push_back(vec0);
            }
            return true;
        }
        else if( this->itype == 2 ) // メッシュの場合
        {
            const unsigned int npo_cad = aRelCoMesh.size()/2;
            double lenrel_tot = 0;
            for(unsigned int idiv=0;idiv<npo_cad+1;idiv++){
                double x0,y0;
                if( idiv == 0 ){ x0=0; y0=0; }
                else{ x0=aRelCoMesh[idiv*2-2]; y0=aRelCoMesh[idiv*2-1]; }
                double x1,y1;
                if( idiv == npo_cad ){ x1=1; y1=0; }
                else{ x1=aRelCoMesh[idiv*2+0]; y1=aRelCoMesh[idiv*2+1]; }
                lenrel_tot += sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
            }
            const unsigned int ndiv1 = ndiv;
            const unsigned int npo1 = ndiv-1;
            const double divrel = lenrel_tot / ndiv1;
            const Com::CVector2D& gh = po_e - po_s;
            const Com::CVector2D gv(-gh.y, gh.x);
            aCo.reserve(npo1);
            unsigned int idiv0_cur = 0;
            double ratio0_cur = 0;
            for(unsigned int ipo1=0;ipo1<npo1;ipo1++){
                double restlen0_cur = divrel;
                for(;;){
                    assert( idiv0_cur < npo_cad+1 );
                    assert( ratio0_cur >= 0 && ratio0_cur <= 1);
                    double x0,y0;
                    if( idiv0_cur == 0 ){ x0=0; y0=0; }
                    else{ x0=aRelCoMesh[idiv0_cur*2-2]; y0=aRelCoMesh[idiv0_cur*2-1]; }
                    double x1,y1;
                    if( idiv0_cur == npo_cad ){ x1=1; y1=0; }
                    else{ x1=aRelCoMesh[idiv0_cur*2+0]; y1=aRelCoMesh[idiv0_cur*2+1]; }
                    const double len_idiv0 = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
                    if( len_idiv0*(1-ratio0_cur) > restlen0_cur ){
                        ratio0_cur += restlen0_cur/len_idiv0;
                        const double xintp = x0*(1.0-ratio0_cur)+x1*ratio0_cur;
                        const double yintp = y0*(1.0-ratio0_cur)+y1*ratio0_cur;
                        aCo.push_back( po_s + gh*xintp + gv*yintp );
                        break;
                    }
                    restlen0_cur -= len_idiv0*(1-ratio0_cur);
                    idiv0_cur += 1;
                    ratio0_cur = 0;
                }
            }
            return true;
        }
    }
    return false;
}

////////////////////////////////

bool Cad::IsCross_LineSeg_LineSeg(const Com::CVector2D& po_s0, const Com::CVector2D& po_e0,
                                  const Com::CVector2D& po_s1, const Com::CVector2D& po_e1 )
{
	{
		const double min0x = ( po_s0.x < po_e0.x ) ? po_s0.x : po_e0.x;
		const double max0x = ( po_s0.x > po_e0.x ) ? po_s0.x : po_e0.x;
		const double max1x = ( po_s1.x > po_e1.x ) ? po_s1.x : po_e1.x;
		const double min1x = ( po_s1.x < po_e1.x ) ? po_s1.x : po_e1.x;
		const double min0y = ( po_s0.y < po_e0.y ) ? po_s0.y : po_e0.y;
		const double max0y = ( po_s0.y > po_e0.y ) ? po_s0.y : po_e0.y;
		const double max1y = ( po_s1.y > po_e1.y ) ? po_s1.y : po_e1.y;
		const double min1y = ( po_s1.y < po_e1.y ) ? po_s1.y : po_e1.y;
		const double len = ((max0x-min0x)+(max0y-min0y)+(max1x-min1x)+(max1y-min1y))*0.0001;
//		std::cout << len << std::endl;
		if( max1x+len < min0x ) return false;
		if( max0x+len < min1x ) return false;
		if( max1y+len < min0y ) return false;
		if( max0y+len < min1y ) return false;
	}
	const double area1 = Com::TriArea(po_s0,po_e0,po_s1);
	const double area2 = Com::TriArea(po_s0,po_e0,po_e1);
	const double area3 = Com::TriArea(po_s1,po_e1,po_s0);
	const double area4 = Com::TriArea(po_s1,po_e1,po_e0);  
//	std::cout << area1 << " " << area2 << " " << area3 << " " << area4 << std::endl;
  const double a12 = area1*area2; if( a12 > 0 ) return false;
  const double a34 = area3*area4; if( a34 > 0 ) return false;
  return true;
  /*
  if( fabs(a12) > fabs(a34) ){
    if( a12 > 0 ){
      if( fabs(area1) > fabs(area2)*0.000000001 && fabs(area1) > fabs(area2)*0.000000001 ){ return 0; }
      return -1;
    }
    if( a34 > 0 ){
      if( fabs(area3) > fabs(area4)*0.000000001 && fabs(area3) > fabs(area4)*0.000000001 ){ return 0; }
      return -1;
    }
  }
  else{
    if( a34 > 0 ){
      if( fabs(area3) > fabs(area4)*0.000000001 && fabs(area3) > fabs(area4)*0.000000001 ){ return 0; }
      return -1;
    }
    else if( a12 > 0 ){
      if( fabs(area1) > fabs(area2)*0.000000001 && fabs(area1) > fabs(area2)*0.000000001 ){ return 0; }
      return -1;
    }    
  }
  return 1;
   */
}

bool Cad::IsCross_Circle_Circle(
		const Com::CVector2D& po_c0, double radius0,
		const Com::CVector2D& po_c1, double radius1,
		Com::CVector2D& po0, Com::CVector2D& po1 )
{
	const double sq_dist = Com::SquareLength(po_c0,po_c1);
	const double dist = sqrt( sq_dist );
	if( radius0 + radius1 < dist ) return false;
	if( fabs(radius0 - radius1) > dist ) return false;		
	if( dist < 1.0e-30 ) return false;
	const double ct = 0.5*(sq_dist+radius0*radius0-radius1*radius1)/(radius0*dist);
	assert( ct >= -1 && ct <= 1 );
	const double st = sqrt( 1-ct*ct );
	Com::CVector2D e0 = (po_c1-po_c0)*(1/dist);
	Com::CVector2D e1;
	e1.x = e0.y;	e1.y = -e0.x;
	po0 = po_c0 + e0*(radius0*ct) + e1*(radius0*st);
	po1 = po_c0 + e0*(radius0*ct) - e1*(radius0*st);
	return true;
}

// 線分と円弧の交錯を判定する
int Cad::CEdge2D::NumCross_Arc_LineSeg(const Com::CVector2D& po_s1, const Com::CVector2D& po_e1) const
{
	if( itype != 1 ) return -1;

	Com::CVector2D po_c0;
	double radius0;
	this->GetCenterRadius(po_c0,radius0);

	bool is_out0 = Com::SquareLength(po_c0,po_s1) > radius0*radius0;
	bool is_out1 = Com::SquareLength(po_c0,po_e1) > radius0*radius0;
	// 頂点が両方とも円の中なら交錯しない
	if( !is_out0 && !is_out1 ){	return 0; }

	// どちらか一つの点が辺から見て円弧の側になければ交錯しない
  bool is_arc_side0 = (Com::TriArea(po_e,po_s1,po_s) > 0.0) == is_left_side;
  bool is_arc_side1 = (Com::TriArea(po_e,po_e1,po_s) > 0.0) == is_left_side;
	if( !is_arc_side0 && !is_arc_side1 ){ return 0; }
		
	// 片方が円の内側，片方が円の外側の場合
	if( (!is_out0 && is_out1 ) || (is_out0 && !is_out1 ) ){
		if( is_arc_side0 && is_arc_side1 ){ return 1; }
		// 弧と直線が交錯するか調べる
		bool is_cross;
		if( IsCross_LineSeg_LineSeg(po_s,po_e, po_s1,po_e1) != 0 ){
			if( is_out0 ){ is_cross = is_arc_side0; }
			else{ is_cross = is_arc_side1; }
		}
		else{
			if( is_out0 ){ is_cross = !is_arc_side0; }
			else{ is_cross = !is_arc_side1; }
		}
		if( is_cross ) return 1;
		return 0;
	}

	// 両方が円の外側の場合
	if( is_out0 && is_out1 ){
		// 弧と直線が交錯するか調べる．交錯なら円弧とも交錯
		if( IsCross_LineSeg_LineSeg(po_s,po_e, po_s1,po_e1) != 0 ){ return 1; }
		Com::CVector2D foot;
		{	// c0を線分s1e1に下ろした垂線の足を計算する
			Com::CVector2D v01(po_e1.x-po_s1.x,  po_e1.y-po_s1.y );
			Com::CVector2D vc0(po_s1.x-po_c0.x, po_s1.y-po_c0.y);
			Com::CVector2D vc1(po_e1.x-po_c0.x, po_e1.y-po_c0.y);
			const double dc0 = Com::Dot(vc0,v01);
			const double dc1 = Com::Dot(vc1,v01);
			// 垂線の足が辺上に無い場合は交錯しない
			if( dc0*dc1 > 0.0 ){ return 0; }
			const double r0 = -dc0 / (-dc0+dc1);
			const double r1 =  dc1 / (-dc0+dc1);
			foot = po_s1*r1+po_e1*r0;
		}
		 // 垂線の足が円の外なら交錯しない
		if( Com::SquareLength(po_c0,foot) > radius0*radius0 ){ return 0; }
		// 垂線の足が弧の側になければ交錯しない
        if( (Com::TriArea( po_s, foot, po_e ) > 0.0) == is_left_side ){ return 0; }
		return 2;
	}
	assert(0);
	return 0;
}

// 円弧と辺(直線でなければならない)の交点を求める
// ２つある場合は交点のposからpoeへのパラメータがt1,t2に入り、返り値が１となる
// 交点が無い場合は０が返り値となる
bool Cad::IsCross_Line_Circle(
		const Com::CVector2D& poc, const double radius, 
		const Com::CVector2D& pos, const Com::CVector2D& poe,
		double& t0, double& t1)
{
	{
		double min_x = ( pos.x < poe.x ) ? pos.x : poe.x; if( poc.x+radius < min_x ) return false;
		double max_x = ( pos.x > poe.x ) ? pos.x : poe.x; if( poc.x-radius > max_x ) return false;
		double min_y = ( pos.y < poe.y ) ? pos.y : poe.y; if( poc.y+radius < min_y ) return false;
		double max_y = ( pos.y > poe.y ) ? pos.y : poe.y; if( poc.y-radius > max_y ) return false;
	}

	const Com::CVector2D& es = poe-pos;
	const Com::CVector2D& cs = poc-pos;
	const double a = es.SqLength();
	const double b = Com::Dot(es,cs);
	const double c = cs.SqLength() - radius*radius;
	double det = b*b-a*c;
	if( det < 0 ) return false;
	t0 = (b-sqrt(det))/a;
	t1 = (b+sqrt(det))/a;
	return true;
}


// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
double Cad::FindNearestPointParameter_Line_Point(const Com::CVector2D& po_c,
	const Com::CVector2D& po_s, const Com::CVector2D& po_e)
{
	const Com::CVector2D& es = po_e-po_s;
	const Com::CVector2D& sc = po_s-po_c;
	const double a = Com::SquareLength(es);
	const double b = Com::Dot(es,sc);
	return - b/a;
}

// get parameter 't' of the line against point. t=0 is po_s, t=1 is po_e
double Cad::GetDist_LineSeg_Point(const Com::CVector2D& po_c,
                                  const Com::CVector2D& po_s, const Com::CVector2D& po_e)
{
	const Com::CVector2D& es = po_e-po_s;
	const Com::CVector2D& sc = po_s-po_c;
	const double a = Com::SquareLength(es);
	const double b = Com::Dot(es,sc);
	const double t =  - b/a;
  if( t < 0 ){ return Distance(po_s,po_c); }
  if( t > 1 ){ return Distance(po_e,po_c); }
  Com::CVector2D p = po_s + t*(po_e-po_s);
  return Distance(p,po_c);
}


// 弦と弧で張られる領域内部に点poが入っているかを調べる
int Cad::CEdge2D::IsInsideArcSegment(const Com::CVector2D& po) const
{
	if( itype != 1 ) return -1;
	Com::CVector2D po_c0;
	double radius;
	this->GetCenterRadius(po_c0,radius);
	if( (po_c0.x-po.x)*(po_c0.x-po.x)+(po_c0.y-po.y)*(po_c0.y-po.y) > radius*radius ){ return 0; }
    if( (Com::TriArea( po_e, po, po_s ) > 0.0) == is_left_side ){ return 1; }
	return 0;
}

// 円弧の中心からみて，点poと円弧が同じ方向に重なっているか？
int Cad::CEdge2D::IsDirectionArc(const Com::CVector2D& po) const{
	if( itype != 1 ) return -1;
	Com::CVector2D po_c;
	double radius;
	this->GetCenterRadius(po_c,radius);
	if( is_left_side ){
		if( Com::TriArea(po_s,po_c,po_e) > 0.0 ){
			if(    Com::TriArea(po_s,po_c,po) > 0.0 && Com::TriArea(po,po_c,po_e) > 0.0 ){ return 1; }
			else{ return 0; }
		}
		else{
			if(    Com::TriArea(po_s,po_c,po) > 0.0 || Com::TriArea(po,po_c,po_e) > 0.0 ){ return 1; }
			else{ return 0; }
		}
	}
	else{
		if( Com::TriArea(po_e,po_c,po_s) > 0.0 ){
			if(    Com::TriArea(po_e,po_c,po) > 0.0 && Com::TriArea(po,po_c,po_s) > 0.0 ){ return 1; }
			else{ return 0; }
		}
		else{
			if(    Com::TriArea(po_e,po_c,po) > 0.0 || Com::TriArea(po,po_c,po_s) > 0.0 ){ return 1; }
			else{ return 0; }
		}
	}
	return -1;
}

bool Cad::CEdge2D::ConnectEdge(const Cad::CEdge2D& e1, bool is_add_ahead, bool is_same_dir)
{
	if(       is_add_ahead &&  is_same_dir ){ assert( id_v_e == e1.id_v_s ); }
	else if(  is_add_ahead && !is_same_dir ){ assert( id_v_e == e1.id_v_e ); }
	else if( !is_add_ahead &&  is_same_dir ){ assert( id_v_s == e1.id_v_e ); }
	else if( !is_add_ahead && !is_same_dir ){ assert( id_v_s == e1.id_v_s ); }
	if( this->itype == 2 ){		
		std::vector<Com::CVector2D> aPo0;
		{
			const Com::CVector2D& ps = this->po_s;
			const Com::CVector2D& pe = this->po_e;
			const Com::CVector2D& h0 = pe-ps;
			const Com::CVector2D v0(-h0.y,h0.x);
			const std::vector<double>& relcomsh = this->aRelCoMesh;
			const unsigned int npo = relcomsh.size()/2;
			aPo0.resize(npo);
			for(unsigned ipo=0;ipo<npo;ipo++){
				aPo0[ipo] = ps + h0*relcomsh[ipo*2+0] + v0*relcomsh[ipo*2+1];
			}
		}
		const double ave_elen = this->GetCurveLength()/(aPo0.size()+1.0);
		std::vector<Com::CVector2D> aPo1;
		{
			const unsigned int ndiv1 = e1.GetCurveLength()/ave_elen;
			e1.GetCurve_Mesh(aPo1,ndiv1);
		}
		std::vector<Com::CVector2D> aPo2;
		if( is_add_ahead ){
			aPo2 = aPo0;
			aPo2.push_back(this->po_e);
			if( is_same_dir ){
				for(unsigned int i=0;i<aPo1.size();i++){ aPo2.push_back(aPo1[i]); }
				this->id_v_e = e1.id_v_e;
				this->po_e = e1.po_e;
			}
			else{
				for(int i=aPo1.size()-1;i>=0;i--){ aPo2.push_back(aPo1[i]); }
				this->id_v_e = e1.id_v_s;
				this->po_e = e1.po_s;
			}
		}
		else{
			if( is_same_dir ){
				aPo2 = aPo1;
				aPo2.push_back(this->po_s);
				for(unsigned int i=0;i<aPo0.size();i++){ aPo2.push_back(aPo0[i]); }
				this->id_v_s = e1.id_v_s;
				this->po_s = e1.po_s;
			}
			else{
				for(int i=aPo1.size()-1;i>=0;i--){ aPo2.push_back(aPo1[i]); }
				aPo2.push_back(this->po_s);
				for(unsigned int i=0;i<aPo0.size();i++){ aPo2.push_back(aPo0[i]); }
				this->id_v_s = e1.id_v_e;
				this->po_s = e1.po_e;
			}
		}
/*		std::cout << po_s.x << " " << po_s.y << std::endl;
		for(unsigned int i=0;i<aPo2.size();i++){
			std::cout << i << " " << aPo2[i].x << " " << aPo2[i].y << std::endl;
		}
		std::cout << po_e.x << " " << po_e.y << std::endl;*/
		{
			const Com::CVector2D& ps = this->po_s;
			const Com::CVector2D& pe = this->po_e;
			const unsigned int npo = aPo2.size();
			this->aRelCoMesh.resize(npo*2);
			const double sqlen = Com::SquareLength(pe-ps);
			const Com::CVector2D& eh = (pe-ps)*(1/sqlen);
			const Com::CVector2D ev(-eh.y,eh.x);
			for(unsigned int ipo=0;ipo<npo;ipo++){
				double x1 = Com::Dot(aPo2[ipo]-ps,eh);
				double y1 = Com::Dot(aPo2[ipo]-ps,ev);
				this->aRelCoMesh[ipo*2+0] = x1;
				this->aRelCoMesh[ipo*2+1] = y1;
			}
		}
	}
	return true;
}

bool Cad::CEdge2D::Split(Cad::CEdge2D& edge_a, const Com::CVector2D& pa)
{
	if(      this->itype == 0 ){}
	else if( this->itype == 1 ){
		const Com::CVector2D& ps = this->po_s;
		const Com::CVector2D& pe = this->po_e;
		Com::CVector2D pc; double r;
		this->GetCenterRadius(pc,r);
		this->dist = Com::TriHeight(pc,ps,pa)/Com::Distance(ps,pa);
		edge_a.itype = 1;
		edge_a.is_left_side = this->is_left_side;
		edge_a.dist = Com::TriHeight(pc,pa,pe)/Com::Distance(pa,pe);
	}
	else if( this->itype == 2 ){
		const Com::CVector2D& ps = this->po_s;
		const Com::CVector2D& pe = this->po_e;
		std::vector<Com::CVector2D> aPo;
		{
			const Com::CVector2D& h0 = pe-ps;
			const Com::CVector2D v0(-h0.y,h0.x);
			const std::vector<double>& relcomsh = this->aRelCoMesh;
			const unsigned int npo = relcomsh.size()/2;
			aPo.resize(npo);
			for(unsigned ipo=0;ipo<npo;ipo++){
				aPo[ipo] = ps + h0*relcomsh[ipo*2+0] + v0*relcomsh[ipo*2+1];
			}
		}
		int ipo0_e;
		int ipo1_s;
		{
	 		bool is_segment = false;
			int ind = -1;
			double min_dist = Com::Distance(ps,pa);
			for(unsigned int idiv=0;idiv<aPo.size()+1;idiv++){
				Com::CVector2D poi0 = ( idiv==0          ) ? ps : aPo[idiv-1];
				Com::CVector2D poi1 = ( idiv==aPo.size() ) ? pe : aPo[idiv];
				if( Com::Distance(pa,poi1) < min_dist ){ 
					is_segment = false;
					ind = idiv;
					min_dist = Com::Distance(pa,poi1);	
				}
				const double t = FindNearestPointParameter_Line_Point(pa,poi0,poi1);
				if( t < 0.01 || t > 0.99 ) continue;
				Com::CVector2D po_mid = poi0 + (poi1-poi0)*t;
				if( Com::Distance(pa,po_mid) < min_dist ){ 
					is_segment = true;
					min_dist = Com::Distance(pa,po_mid);	
					ind = idiv;
				}
			}
			////////////////
			if( is_segment ){	ipo0_e = ind-2;	ipo1_s = ind+1;	}
			else{				ipo0_e = ind-1;	ipo1_s = ind+1; }
		}
		////////////////////////////////
		if( ipo0_e > 0 ){
			this->itype = 2;
			const unsigned int npo = ipo0_e+1;
			this->aRelCoMesh.resize(npo*2);
			const double sqlen = Com::SquareLength(pa-ps);
			const Com::CVector2D& eh = (pa-ps)*(1/sqlen);
			const Com::CVector2D ev(-eh.y,eh.x);
	//		double x0=0, y0=0;
	//		const double elen = 0.1;
			for(unsigned int ipo=0;ipo<npo;ipo++){
				double x1 = Com::Dot(aPo[ipo]-ps,eh);
				double y1 = Com::Dot(aPo[ipo]-ps,ev);
				this->aRelCoMesh[ipo*2+0] = x1;
				this->aRelCoMesh[ipo*2+1] = y1;
			}
		}
		else{ this->itype = 0; }
		////////////////
		if( ipo1_s < (int)aPo.size() ){
			edge_a.itype = 2;
			const unsigned int npo = aPo.size()-ipo1_s;
			edge_a.aRelCoMesh.resize(npo*2);
			const double sqlen = Com::SquareLength(pe-pa);
			const Com::CVector2D& eh = (pe-pa)*(1/sqlen);
			const Com::CVector2D ev(-eh.y,eh.x);
	//		double x0=0, y0=0;
	//		const double elen = 0.1;
			for(unsigned int ipo=0;ipo<npo;ipo++){
				double x1 = Com::Dot(aPo[ipo+ipo1_s]-pa,eh);
				double y1 = Com::Dot(aPo[ipo+ipo1_s]-pa,ev);
				edge_a.aRelCoMesh[ipo*2+0] = x1;
				edge_a.aRelCoMesh[ipo*2+1] = y1;
			}
		}
		else{ edge_a.itype = 0; }
	}
	return true;
}


// v0からエッジに沿った距離でlenの長さにある点を得る．
// is_front==true 沿ならエッジに沿って，is_front==falseならエッジに沿わない
bool Cad::CEdge2D::GetPointOnCurve_OnCircle
(const Com::CVector2D& v0, double len, bool is_front,
 bool& is_exceed, Com::CVector2D& out) const
{
  if( len <= 0 ) return false;
  if( this->itype == 0 ) // 直線の場合
  {
    const Com::CVector2D& v1 = this->GetNearestPoint(v0);
    if( is_front ){
      double len1e = Com::Distance(v1,this->po_e);
      if( len1e < len ){ is_exceed=true; out=po_e; return true; }
      out = (len/len1e)*po_e + (1-len/len1e)*v1;
      is_exceed = false;
      return true;
    }
    else{
      double len1s = Com::Distance(v1,this->po_s);
      if( len1s < len ){ is_exceed=true; out=po_s; return true; }
      out = (len/len1s)*po_s + (1-len/len1s)*v1;
      is_exceed = false;
      return true;
    }
  }
  return false;
}


//! そのうち交錯位置の情報も返したい
int Cad::CheckEdgeIntersection(const std::vector<CEdge2D>& aEdge)
{
	const unsigned int nedge = aEdge.size();
	for(unsigned int iedge=0;iedge<nedge;iedge++){
		const CEdge2D& e_i = aEdge[iedge];
		if( e_i.IsCrossEdgeSelf() ){ return 1; }
		const unsigned int ipo0 = e_i.id_v_s; 
		const unsigned int ipo1 = e_i.id_v_e;
		// edge_iのバウンディングボックスを取得
    const Com::CBoundingBox2D& bb_i = e_i.GetBoundingBox();
    //		double x_min_i, x_max_i, y_min_i, y_max_i;
    //		e_i.GetBoundingBox(x_min_i,x_max_i, y_min_i,y_max_i);
		////////////////
		for(unsigned int jedge=iedge+1;jedge<nedge;jedge++){
			const CEdge2D& e_j = aEdge[jedge];
			const unsigned int jpo0 = e_j.id_v_s;
			const unsigned int jpo1 = e_j.id_v_e;
			// 共有点が無い場合
			if( (ipo0-jpo0)*(ipo0-jpo1)*(ipo1-jpo0)*(ipo1-jpo1) != 0 ){
				// BoundingBoxを用いて，交錯しないパターンを除外
				// edge_jのバウンディングボックスを取得
        const Com::CBoundingBox2D& bb_j = e_j.GetBoundingBox();        
        //				double x_min_j, x_max_j, y_min_j, y_max_j;
        //				e_j.GetBoundingBox(x_min_j,x_max_j, y_min_j,y_max_j);
        if( bb_j.x_min > bb_i.x_max || bb_j.x_max < bb_i.x_min ) continue;
        if( bb_j.y_min > bb_i.y_max || bb_j.y_max < bb_i.y_min ) continue;
				////////////////
        //				if( x_min_j > x_max_i || x_max_j < x_min_i ) continue;	// 交錯がありえないパターンを除外
        //				if( y_min_j > y_max_i || y_max_j < y_min_i ) continue;	// 上に同じ
				// 交点が無いか判断する
				if( !e_i.IsCrossEdge(e_j) ){ continue; }
				return 1;
			}			
			if(      ipo0 == jpo0 && ipo1 == jpo1 ){
				if( e_i.IsCrossEdge_ShareBothPoints(e_j,true) == 1 ){ return 1; }
			}
			else if( ipo0 == jpo1 && ipo1 == jpo0 ){
				if( e_i.IsCrossEdge_ShareBothPoints(e_j,false) == 1 ){ return 1; }
			}
			else if( ipo0 == jpo0 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j, true, true)==1 ){ return 1; }
			}
			else if( ipo0 == jpo1 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j, true,false)==1 ){ return 1; }
			}
			else if( ipo1 == jpo0 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j,false, true)==1 ){ return 1; }
			}
			else if( ipo1 == jpo1 ){ 
				if( e_i.IsCrossEdge_ShareOnePoint(e_j,false,false)==1 ){ return 1; }
			}
			continue;
		}
	}
  //	std::cout << "Intersect Dosen't Occur" << std::endl;
	return 0;
}
