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
@brief 複素数クラス(Com::Complex)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(_COMPLEX_H)
#define _COMPLEX_H

#include <math.h>

namespace Com{

//! 複素数クラス
class Complex
{
	inline friend Complex operator +(const Complex& rhs, const Complex& lhs);
	inline friend Complex operator -(const Complex& rhs, const Complex& lhs);
	inline friend Complex operator *(const Complex& rhs, const Complex& lhs);
	inline friend Complex operator *(const double&, const Complex&);
	inline friend Complex operator *(const Complex&, const double&);
	inline friend Complex operator /(const Complex& rhs, const Complex& lhs);
	inline friend Complex operator /(const double&, const Complex&);
	inline friend Complex operator /(const Complex&, const double&);
	inline friend Complex operator -(const Complex&);
	inline friend Complex Conjugate(const Complex&);

	inline friend double SquaredNorm(const Complex& c);
	inline friend double Norm(const Complex& c);
	inline friend void SetZero(Complex& c);

	//! 内積
	inline friend Complex InnerProduct(const Complex& rhs, const Complex& lhs);

public:
	inline Complex& operator=(const double& lhs);	//!< 代入

	inline Complex& operator *= ( const double& );	//!< スカラー倍
	inline Complex& operator *= ( const Complex& );	//!< 複素数との積
	inline Complex& operator += ( const Complex& );	//!< 加える
	inline Complex& operator -= ( const Complex& );	//!< 引く

	Complex(const double& real, const double& imag );	//!< 実部，虚部で初期化
	Complex(const double& real);
	Complex(){};

	inline double Real() const { return real; }	//!< Get Real Part
	inline double Imag() const { return imag; }	//!< Get Imaginary Part

private:
	double real;	//!< 実部
	double imag;	//!< 虚部
};

inline Complex::Complex(const double& real, const double& imag ) : real(real), imag(imag){}

inline Complex::Complex(const double& real) : real(real), imag(0.0){}

inline Complex& Complex::operator = ( const double& d ){
	real=d; imag=0.0;
	return *this;
}

inline Complex& Complex::operator *= ( const double& d ){
	real *= d; imag *= d;
	return *this;
}

inline Complex& Complex::operator *= ( const Complex& c ){
	double rel0 = real;
	real = real*c.real-imag*c.imag;
	imag = rel0*c.imag+imag*c.real;
	return *this;
}

inline Complex& Complex::operator += ( const Complex& c ){
	real += c.real; imag += c.imag;
	return *this;
}

inline Complex& Complex::operator -= ( const Complex& c ){
	real -= c.real; imag -= c.imag;
	return *this;
}

inline Complex operator -(const Complex& rhs){
	return Complex(-rhs.real,-rhs.imag);
}

inline Complex Conjugate(const Complex& rhs){
	return Complex(rhs.real,-rhs.imag);
}

inline Complex operator +(const Complex& rhs, const Complex& lhs){
	return Complex(rhs.real+lhs.real,rhs.imag+lhs.imag);
}

inline Complex operator -(const Complex& rhs, const Complex& lhs){
	return Complex(rhs.real-lhs.real,rhs.imag-lhs.imag);
}

inline Complex operator *( const Complex& rhs, const Complex& lhs ){
	return Complex(rhs.real*lhs.real-rhs.imag*lhs.imag, 
		rhs.real*lhs.imag+rhs.imag*lhs.real);
}

inline Complex operator *(const double& d, const Complex& lhs){
	return Complex(d*lhs.real,d*lhs.imag);
}

inline Complex operator *(const Complex& lhs, const double& d){
	return Complex(d*lhs.real,d*lhs.imag);
}

inline Complex operator /( const Complex& rhs, const Complex& lhs ){
	const double fact = 1.0 / (lhs.real*lhs.real+lhs.imag*lhs.imag);
	return Complex((rhs.real*lhs.real+rhs.imag*lhs.imag)*fact,
		(-rhs.real*lhs.imag+rhs.imag*lhs.real)*fact );
}

inline Complex operator /(const double& d, const Complex& c){
	const double fact = d / (c.real*c.real+c.imag*c.imag);
	return Complex(fact*c.real,-fact*c.imag);
}

inline Complex operator /(const Complex& c, const double& d){
	return Complex( c.real/d, c.imag/d );
}

inline double SquaredNorm(const Complex& c){ return c.real*c.real+c.imag*c.imag; }
inline double Norm(const Complex& c){ return sqrt(c.real*c.real+c.imag*c.imag); }

inline void SetZero(Complex& c){ c.real = 0.0; c.imag = 0.0; }

inline Complex InnerProduct(const Complex& lhs, const Complex& rhs){
	return Complex(lhs.real*rhs.real+lhs.imag*rhs.imag, 
		+lhs.imag*rhs.real-lhs.real*rhs.imag);
}

}

#endif
