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

////////////////////////////////////////////////////////////////
// Quaternion.cpp: CQuaternion クラスのインプリメンテーション
////////////////////////////////////////////////////////////////

#include "delfem/quaternion.h"

#define PAI 3.141592

using namespace Com;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

Com::CQuaternion::CQuaternion(const CVector3D& axis ){
	AxisToQuat( axis );
}

Com::CQuaternion::CQuaternion(double real, const CVector3D& vector){
	this->real = real;
	this->vector = vector;
}

Com::CQuaternion::CQuaternion(const CVector3D& a_vector, const CVector3D& b_vector){
	VectorTrans(a_vector, b_vector);
}

Com::CQuaternion::~CQuaternion(){
}

//////////////////////////////////////////////////////////////////////
// 非メンバ関数の非フレンド関数
//////////////////////////////////////////////////////////////////////

namespace Com
{

CQuaternion operator+(const CQuaternion& lhs, const CQuaternion& rhs){  /* 四元数同士の加算 */
    CQuaternion temp = lhs;
	temp += rhs;
    return temp;
}

CQuaternion operator-(const CQuaternion& lhs, const CQuaternion& rhs){  /* 四元数同士の減算 */
    CQuaternion temp = lhs;
	temp -= rhs;
    return temp;
}

CQuaternion operator*(const CQuaternion& quat, double d){  /* 実数との乗算 */
    CQuaternion temp = quat;
	temp *= d;
    return temp;
}

CQuaternion operator*(double d, const CQuaternion& quat){/* 実数が左の乗算 */
    CQuaternion temp = quat;
    temp *= d;
	return temp;
}

Com::CQuaternion operator*(const CQuaternion& lhs, const CQuaternion& rhs){ /* 四元数同士の積 */
    CQuaternion temp = lhs;
    temp *= rhs;
	return temp;
}

CVector3D Rotate(const CQuaternion& quat, const CVector3D& vec){
	CQuaternion tmp(0, vec);
	tmp = quat.GetConjugate() *  tmp * quat ;
	return tmp.GetVector();
}

CVector3D UnRotate(const CQuaternion& quat, const CVector3D& vec){
	CQuaternion tmp(0, vec);
	tmp = quat * tmp * quat.GetConjugate();
	return tmp.GetVector();
}

}

//////////////////////////////////////////////////////////////////////
// メンバ関数のフレンド関数
//////////////////////////////////////////////////////////////////////

namespace Com
{

bool operator==(const CQuaternion& lhs, const CQuaternion& rhs){
	if( fabs(lhs.real - rhs.real) < NEARLY_ZERO && lhs == rhs )	return true;
	else return false;
}

bool operator!=(const CQuaternion& lhs, const CQuaternion& rhs){
	if( lhs == rhs ) return false;
	else return true;
}

}

//////////////////////////////////////////////////////////////////////
// メンバ関数の非フレンド関数　
//////////////////////////////////////////////////////////////////////

CQuaternion& CQuaternion::operator = (const CQuaternion& rhs){
	if( this != &rhs ){
		real = rhs.real;
		vector = rhs.vector;
	}
	return *this;
}

CQuaternion& CQuaternion::operator *= (const CQuaternion& rhs){
	const double last_real = real;
	real = real * rhs.real - Dot( vector, rhs.vector );
	vector = ( (rhs.vector * last_real) + (vector * rhs.real) ) + Cross( vector, rhs.vector );
	return *this;
}

CQuaternion& CQuaternion::operator *= (double d){
	real *= d;
	vector *= d;
	return *this;
}

CQuaternion& CQuaternion::operator += (const CQuaternion& rhs){
	real += rhs.real;
	vector += rhs.vector;
	return *this;
}

CQuaternion& CQuaternion::operator -= (const CQuaternion& rhs){
	real -= rhs.real;
	vector -= rhs.vector;
	return *this;
}

CQuaternion CQuaternion::GetConjugate() const{
	CQuaternion tmp(real, vector*(-1.0));
	return tmp;
}

void CQuaternion::RotMatrix33(double* m) const
{
	double vx=vector.x, vy=vector.y, vz=vector.z;

	m[0] = 1.0 - 2.0 * ( vy * vy + vz * vz );
	m[1] = 2.0 * ( vx * vy - vz * real );
	m[2] = 2.0 * ( vz * vx + vy * real );

	m[3] = 2.0 * ( vx * vy + vz * real );
	m[4] = 1.0 - 2.0 * ( vz * vz + vx * vx );
	m[5] = 2.0 * ( vy * vz - vx * real );

	m[6] = 2.0 * ( vz * vx - vy * real );
	m[7] = 2.0 * ( vy * vz + vx * real );
	m[8] = 1.0 - 2.0 * ( vy * vy + vx * vx );
}

void CQuaternion::Normalize()
{	
	const double len = ( real * real + vector.DLength() );
	real /= len;
	vector /= len;
}

void CQuaternion::AxisToQuat(const CVector3D &axis )
{
	const double phi = axis.Length();
	if( phi < 1.0e-30 ){
		vector = CVector3D(0,0,0);
		real = 1.0;
		return;
	}
	vector = axis;
	vector.Normalize();
	vector *= sin( phi * 0.5 );
	real = cos( phi * 0.5 );
}

void CQuaternion::VectorTrans(const CVector3D& a_vector, const CVector3D& b_vector){
	vector = Cross(a_vector, b_vector);
	if( vector.DLength() < NEARLY_ZERO ){
		real = 1.0;
		vector.SetZero();
		return;
	}
	vector.Normalize();
	double cos_theta = Dot(a_vector, b_vector) / ( a_vector.Length() * b_vector.Length() );
	real = sqrt( 0.5*(1+cos_theta) );
	vector *= sqrt( 0.5*(1-cos_theta) );
}


double CQuaternion::Length()
{
	return sqrt( real*real + vector.DLength() );
}
