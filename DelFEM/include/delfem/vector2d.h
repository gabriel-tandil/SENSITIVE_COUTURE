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
@brief two-dimensional vector class (Com::CVector2D)
@author Nobuyuki Umetani
*/

#if !defined(VECTOR2D_H)
#define VECTOR2D_H

#if defined(__VISUALC__)
	#pragma warning(disable:4786)
#endif

#include <assert.h>
#include <vector>
#include "math.h"

namespace Com
{

class CVector2D;

CVector2D operator*(double, const CVector2D&);
CVector2D operator*(const CVector2D&, double);

//! ２次元ベクトルクラス
class CVector2D{
public:
	CVector2D(){}
	CVector2D( const CVector2D& rhs ){
		this->x = rhs.x;
		this->y = rhs.y;
	}
	CVector2D(double x, double y){
		this->x = x;
		this->y = y;
	}
	
	friend double Dot(const CVector2D&, const CVector2D&);

	// オペレータ定義
	inline CVector2D& operator+=(const CVector2D& rhs){
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	inline CVector2D& operator-=(const CVector2D& rhs){
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}
	inline CVector2D& operator*=(double scale){
		x *= scale;
		y *= scale;
		return *this;
    }
	inline CVector2D operator+(const CVector2D& rhs) const {
		CVector2D v = *this;
		return v += rhs;
	}
	inline CVector2D operator-(const CVector2D& rhs) const {
		CVector2D v = *this;
		return v -= rhs;
	}
	//! 長さを正規化する
	inline void Normalize(){
		const double mag = Length();
		x /= mag;
		y /= mag;
	}
	//! 座標値に０を代入する
	inline void SetZero(){
		x = 0.0;
		y = 0.0;
	}
	//! ベクトルの長さを計算する
	double Length() const{
		return sqrt( x*x+y*y );
	}
	//! ベクトルの長さの２乗を計算する
	double SqLength() const{
		return x*x+y*y;
	}
public:
	double x;	//!< ｘ座標値
	double y;	//!< ｙ座標値
};
  
  
//! 2D bounding box class
class CBoundingBox2D
{
public:
  CBoundingBox2D(){
    x_min=0;	x_max=0;
    y_min=0;	y_max=0;
    isnt_empty = false;
  }
  CBoundingBox2D(double x_min0,double x_max0,  double y_min0,double y_max0)
  : x_min(x_min0),x_max(x_max0),  y_min(y_min0),y_max(y_max0)
  {
    assert( x_min <= x_max );
    assert( y_min <= y_max );
    isnt_empty = true;
  }
  CBoundingBox2D( const CBoundingBox2D& bb )
  : x_min(bb.x_min),x_max(bb.x_max), y_min(bb.y_min),y_max(bb.y_max),  
  isnt_empty(bb.isnt_empty){}
  
  CBoundingBox2D& operator+=(const CBoundingBox2D& bb)
  {
    if( !bb.isnt_empty ) return *this;
    if( !isnt_empty ){
      x_max = bb.x_max;	x_min = bb.x_min;
      y_max = bb.y_max;	y_min = bb.y_min;
      this->isnt_empty = bb.isnt_empty;
      return *this;
    }
    x_max = ( x_max > bb.x_max ) ? x_max : bb.x_max;
    x_min = ( x_min < bb.x_min ) ? x_min : bb.x_min;
    y_max = ( y_max > bb.y_max ) ? y_max : bb.y_max;
    y_min = ( y_min < bb.y_min ) ? y_min : bb.y_min;
    return *this;
  }
  bool IsInside(const CVector2D& vec)
  {
    if( !isnt_empty ) return false; // 何もない場合は常に偽
    if(   vec.x >= x_min && vec.x <= x_max
       && vec.y >= y_min && vec.y <= y_max ) return true;
    return false;
  }
  bool IsIntersectSphere(const CVector2D& vec, const double radius ) const
  {
    if( !isnt_empty ) return false; // 何もない場合は常に偽
    if( vec.x < x_min-radius || vec.x > x_max+radius ||
       vec.y < y_min-radius || vec.y > y_max+radius ) return false;
    // 干渉しないやつが含まれているが、アバウトなままでよい．
    return true;
  }
  bool IsIntersect(const CBoundingBox2D& bb_j, double clearance) const
  {
    if( bb_j.x_min > x_max+clearance || bb_j.x_max < x_min-clearance ) return false;
    if( bb_j.y_min > y_max+clearance || bb_j.y_max < y_min-clearance ) return false;
    return true;
  }
public:
  double x_min,x_max,  y_min,y_max;
  bool isnt_empty;	//!< false if there is nothing inside
};  
  


inline CVector2D operator*(double c, const CVector2D& v0)
{
    return Com::CVector2D(v0.x*c,v0.y*c);
}

inline CVector2D operator*(const CVector2D& v0, double c)
{
    return Com::CVector2D(v0.x*c,v0.y*c);
}


////////////////////////////////////////////////////////////////

//! Area of the Triangle
inline double TriArea(const CVector2D& v1, const CVector2D& v2, const CVector2D& v3){
	return 0.5*( (v2.x-v1.x)*(v3.y-v1.y) - (v3.x-v1.x)*(v2.y-v1.y) );
}
//! Area of the Triangle (3 indexes and vertex array)
inline double TriArea(const int iv1, const int iv2, const int iv3, 
				   const std::vector<CVector2D>& point ){
	return TriArea(point[iv1],point[iv2],point[iv3]);
}


////////////////

//! 長さの２乗
inline double SquareLength(const CVector2D& ipo0, const CVector2D& ipo1){
	return	( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y );
}

//! 長さの２乗
inline double SquareLength(const CVector2D& point){
	return	point.x*point.x + point.y*point.y;
}

////////////////

//! 長さ
inline double Distance(const CVector2D& ipo0, const CVector2D& ipo1){
	return	sqrt( ( ipo1.x - ipo0.x )*( ipo1.x - ipo0.x ) + ( ipo1.y - ipo0.y )*( ipo1.y - ipo0.y ) );
}

////////////////

//! ３角形の高さ
inline double TriHeight(const CVector2D& v1, const CVector2D& v2, const CVector2D& v3){
	const double area = TriArea(v1,v2,v3);
	const double len = sqrt( SquareLength(v2,v3) );
	return area*2.0/len;
}

////////////////

//! 内積の計算
inline double Dot(const CVector2D& ipo0, const CVector2D& ipo1){
	return	ipo0.x*ipo1.x + ipo0.y*ipo1.y;
}

//! 外接円の半径の２乗
inline double SquareCircumradius(
		const CVector2D& p0, 
		const CVector2D& p1, 
		const CVector2D& p2 )
{
	const double area = TriArea(p0,p1,p2);

	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);

	return dtmp0*dtmp1*dtmp2/(16.0*area*area);
}

//! 外接円の中心
inline bool CenterCircumcircle(
      const CVector2D& p0, 
      const CVector2D& p1, 
      const CVector2D& p2, 
      CVector2D& center){

   const double area = TriArea(p0,p1,p2);
   if( fabs(area) < 1.0e-10 ){ return false; }
   const double tmp_val = 1.0/(area*area*16.0);

   const double dtmp0 = SquareLength(p1,p2);
   const double dtmp1 = SquareLength(p0,p2);
   const double dtmp2 = SquareLength(p0,p1);

   const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
   const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
   const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);

   center.x = etmp0*p0.x + etmp1*p1.x + etmp2*p2.x;
   center.y = etmp0*p0.y + etmp1*p1.y + etmp2*p2.y;
   return true;
}


////////////////////////////////

//! ドロネー条件を満たすかどうか調べる
inline int DetDelaunay(
		const CVector2D& p0, 
		const CVector2D& p1, 
		const CVector2D& p2, 
		const CVector2D& p3)
{
	const double area = TriArea(p0,p1,p2);
	if( fabs(area) < 1.0e-10 ){
		return 3;
	}
	const double tmp_val = 1.0/(area*area*16.0);

	const double dtmp0 = SquareLength(p1,p2);
	const double dtmp1 = SquareLength(p0,p2);
	const double dtmp2 = SquareLength(p0,p1);

	const double etmp0 = tmp_val*dtmp0*(dtmp1+dtmp2-dtmp0);
	const double etmp1 = tmp_val*dtmp1*(dtmp0+dtmp2-dtmp1);
	const double etmp2 = tmp_val*dtmp2*(dtmp0+dtmp1-dtmp2);

	const CVector2D out_center(
		etmp0*p0.x + etmp1*p1.x + etmp2*p2.x,
		etmp0*p0.y + etmp1*p1.y + etmp2*p2.y );

	const double qradius = SquareLength(out_center,p0);
	const double qdistance = SquareLength(out_center,p3);

//	assert( fabs( qradius - SquareLength(out_center,p1) ) < 1.0e-10*qradius );
//	assert( fabs( qradius - SquareLength(out_center,p2) ) < 1.0e-10*qradius );

	const double tol = 1.0e-20;
	if( qdistance > qradius*(1.0+tol) ){ return 2; }	// 外接円の外
	else{
		if( qdistance < qradius*(1.0-tol) ){ return 0; }	// 外接円の中
		else{ return 1;	}	// 外接円上
	}
	return 0;
}

////////////////////////////////////////////////
  
  
static inline double TriArea2D(const double v1[2], const double v2[2], const double v3[2]){
  double z = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
  return z*0.5;
}
  
inline double SqDistance2D(const double v1[2], const double v2[2]){
    return (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]);
}

inline double Distance2D(const double v1[2], const double v2[2]){
  return sqrt( (v1[0]-v2[0])*(v1[0]-v2[0]) + (v1[1]-v2[1])*(v1[1]-v2[1]) );
}
  
  


} // end namespace Com

#endif // VEC_KER_H


