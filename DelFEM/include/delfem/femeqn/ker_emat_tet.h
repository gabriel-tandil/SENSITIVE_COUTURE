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


/*! @file:
@brief 四面体要素で要素剛性行列を作る際のUtility関数の集まり
@author Nobuyuki Umetani
*/


#if !defined(KER_EMAT_TET_H)
#define KER_EMAT_TET_H

#include <cassert>

static inline double TetVolume(const double p0[], const double p1[], const double p2[], const double p3[])
{
	double vol = (p1[0]-p0[0])*( (p2[1]-p0[1])*(p3[2]-p0[2]) - (p3[1]-p0[1])*(p2[2]-p0[2]) )
	            -(p1[1]-p0[1])*( (p2[0]-p0[0])*(p3[2]-p0[2]) - (p3[0]-p0[0])*(p2[2]-p0[2]) )
	            +(p1[2]-p0[2])*( (p2[0]-p0[0])*(p3[1]-p0[1]) - (p3[0]-p0[0])*(p2[1]-p0[1]) );
	vol *= 0.1666666666666;
	return vol;
}

// caluculate Derivative of Area Coord
static inline void TetDlDx(double dldx[][3], double a[],
			 const double p0[], const double p1[], const double p2[], const double p3[])
{
	const double vol = TetVolume(p0,p1,p2,p3);
	const double dtmp1 = 1.0 / ( vol * 6.0 );

	a[0] = +dtmp1*( p1[0]*(p2[1]*p3[2]-p3[1]*p2[2])-p1[1]*(p2[0]*p3[2]-p3[0]*p2[2])+p1[2]*(p2[0]*p3[1]-p3[0]*p2[1]) );
	a[1] = -dtmp1*( p2[0]*(p3[1]*p0[2]-p0[1]*p3[2])-p2[1]*(p3[0]*p0[2]-p0[0]*p3[2])+p2[2]*(p3[0]*p0[1]-p0[0]*p3[1]) );
	a[2] = +dtmp1*( p3[0]*(p0[1]*p1[2]-p1[1]*p0[2])-p3[1]*(p0[0]*p1[2]-p1[0]*p0[2])+p3[2]*(p0[0]*p1[1]-p1[0]*p0[1]) );
	a[3] = -dtmp1*( p0[0]*(p1[1]*p2[2]-p2[1]*p1[2])-p0[1]*(p1[0]*p2[2]-p2[0]*p1[2])+p0[2]*(p1[0]*p2[1]-p2[0]*p1[1]) );

	dldx[0][0] = -dtmp1*( (p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2]) );
	dldx[0][1] = +dtmp1*( (p2[0]-p1[0])*(p3[2]-p1[2])-(p3[0]-p1[0])*(p2[2]-p1[2]) );
	dldx[0][2] = -dtmp1*( (p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]) );
	
	dldx[1][0] = +dtmp1*( (p3[1]-p2[1])*(p0[2]-p2[2])-(p0[1]-p2[1])*(p3[2]-p2[2]) );
	dldx[1][1] = -dtmp1*( (p3[0]-p2[0])*(p0[2]-p2[2])-(p0[0]-p2[0])*(p3[2]-p2[2]) );
	dldx[1][2] = +dtmp1*( (p3[0]-p2[0])*(p0[1]-p2[1])-(p0[0]-p2[0])*(p3[1]-p2[1]) );
	
	dldx[2][0] = -dtmp1*( (p0[1]-p3[1])*(p1[2]-p3[2])-(p1[1]-p3[1])*(p0[2]-p3[2]) );
	dldx[2][1] = +dtmp1*( (p0[0]-p3[0])*(p1[2]-p3[2])-(p1[0]-p3[0])*(p0[2]-p3[2]) );
	dldx[2][2] = -dtmp1*( (p0[0]-p3[0])*(p1[1]-p3[1])-(p1[0]-p3[0])*(p0[1]-p3[1]) );
	
	dldx[3][0] = +dtmp1*( (p1[1]-p0[1])*(p2[2]-p0[2])-(p2[1]-p0[1])*(p1[2]-p0[2]) );
	dldx[3][1] = -dtmp1*( (p1[0]-p0[0])*(p2[2]-p0[2])-(p2[0]-p0[0])*(p1[2]-p0[2]) );
	dldx[3][2] = +dtmp1*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );

//	std::cout << dldx[0][0]+dldx[1][0]+dldx[2][0]+dldx[3][0] << std::endl;
//	std::cout << dldx[0][1]+dldx[1][1]+dldx[2][1]+dldx[3][1] << std::endl;
//	std::cout << dldx[0][2]+dldx[1][2]+dldx[2][2]+dldx[3][2] << std::endl;

//	std::cout << a[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1]+dldx[0][2]*p0[2] << std::endl;
//	std::cout << a[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1]+dldx[1][2]*p1[2] << std::endl;
//	std::cout << a[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1]+dldx[2][2]*p2[2] << std::endl;
//	std::cout << a[3]+dldx[3][0]*p3[0]+dldx[3][1]*p3[1]+dldx[3][2]*p3[2] << std::endl;
}


const static unsigned int NIntTetGauss[4] = { // 積分点の数
	1, 4, 5, 16
};
const static double TetGauss[4][16][4] =		// 積分点の位置(r1,r2)と重みの配列
{
	{	// order-1    1point
		{ 0.25, 0.25, 0.25, 1.0 },
	},
	{	// order-2    4point
		{ 0.585410196624968, 0.138196601125015, 0.138196601125015, 0.25 },
		{ 0.138196601125015, 0.585410196624968, 0.138196601125015, 0.25 },
		{ 0.138196601125015, 0.138196601125015, 0.585410196624968, 0.25 },
		{ 0.138196601125015, 0.138196601125015, 0.138196601125015, 0.25 },
	},
	{	// order-3    5point
		{ 0.25, 0.25, 0.25, -0.8 },
		{ 0.5               , 0.1666666666666667, 0.1666666666666667, 0.45 },
		{ 0.1666666666666667, 0.5,                0.1666666666666667, 0.45 },
		{ 0.1666666666666667, 0.1666666666666667, 0.5,                0.45 },
		{ 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.45 },
	},
	{	// order-4    16point
		{ 0.7716429020672371 , 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },
		{ 0.07611903264425430, 0.7716429020672371,  0.07611903264425430, 0.05037379410012282 },
		{ 0.07611903264425430, 0.07611903264425430, 0.7716429020672371,  0.05037379410012282 },
		{ 0.07611903264425430, 0.07611903264425430, 0.07611903264425430, 0.05037379410012282 },

		{ 0.1197005277978019,  0.4042339134672644,  0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.1197005277978019,  0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.4042339134672644,  0.1197005277978019,  0.06654206863329239 },

		{ 0.07183164526766925, 0.4042339134672644,  0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.07183164526766925, 0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.4042339134672644,  0.07183164526766925, 0.06654206863329239 },

		{ 0.1197005277978019,  0.07183164526766925, 0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.1197005277978019,  0.07183164526766925, 0.06654206863329239 },
		{ 0.07183164526766925, 0.4042339134672644,  0.1197005277978019,  0.06654206863329239 },

		{ 0.07183164526766925, 0.1197005277978019,  0.4042339134672644,  0.06654206863329239 },
		{ 0.4042339134672644,  0.07183164526766925, 0.1197005277978019,  0.06654206863329239 },
		{ 0.1197005277978019,  0.4042339134672644,  0.07183164526766925, 0.06654206863329239 },
	}
};

#endif
