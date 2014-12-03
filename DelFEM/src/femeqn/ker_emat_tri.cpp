
#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif
#define for if(0);else for

#include <cassert>
#include <math.h>

#include "delfem/femeqn/ker_emat_tri.h"

// 三角形の面積を求める関数
double TriArea(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

// ある点の面積座標を求める関数
void TriAreaCoord(double vc_p[],
				  const double p0[], const double p1[], const double p2[], const double pb[] ){

	vc_p[0] = TriArea( pb, p1, p2 );
	vc_p[1] = TriArea( p0, pb, p2 );
	vc_p[2] = TriArea( p0, p1, pb );

	const double area = TriArea(p0,p1,p2);
	const double inv_area = 1.0 / area;

	vc_p[0] = vc_p[0]*inv_area;
	vc_p[1] = vc_p[1]*inv_area;
	vc_p[2] = vc_p[2]*inv_area;

	assert( fabs( vc_p[0]+vc_p[1]+vc_p[2] - 1.0 ) < 1.0e-15 );
}

// 三角形の一次補間関数の微分とその定数成分
void TriDlDx(double dldx[][2], double const_term[],
			 const double p0[], const double p1[], const double p2[]){

	const double area = TriArea(p0,p1,p2);
	const double tmp1 = 0.5 / area;

	const_term[0] = tmp1*(p1[0]*p2[1]-p2[0]*p1[1]);
	const_term[1] = tmp1*(p2[0]*p0[1]-p0[0]*p2[1]);
	const_term[2] = tmp1*(p0[0]*p1[1]-p1[0]*p0[1]);

	dldx[0][0] = tmp1*(p1[1]-p2[1]);
	dldx[1][0] = tmp1*(p2[1]-p0[1]);
	dldx[2][0] = tmp1*(p0[1]-p1[1]);

	dldx[0][1] = tmp1*(p2[0]-p1[0]);
	dldx[1][1] = tmp1*(p0[0]-p2[0]);
	dldx[2][1] = tmp1*(p1[0]-p0[0]);
/*
	assert( fabs( dldx[0][0]+dldx[1][0]+dldx[2][0] ) < 1.0e-15 );
	assert( fabs( dldx[0][1]+dldx[1][1]+dldx[2][1] ) < 1.0e-15 );

	assert( fabs( const_term[0]+dldx[0][0]*p0[0]+dldx[0][1]*p0[1] - 1.0 ) < 1.0e-10 );
	assert( fabs( const_term[0]+dldx[0][0]*p1[0]+dldx[0][1]*p1[1] ) < 1.0e-10 );
	assert( fabs( const_term[0]+dldx[0][0]*p2[0]+dldx[0][1]*p2[1] ) < 1.0e-10 );

	assert( fabs( const_term[1]+dldx[1][0]*p0[0]+dldx[1][1]*p0[1] ) < 1.0e-10 );
	assert( fabs( const_term[1]+dldx[1][0]*p1[0]+dldx[1][1]*p1[1] - 1.0 ) < 1.0e-10 );
	assert( fabs( const_term[1]+dldx[1][0]*p2[0]+dldx[1][1]*p2[1] ) < 1.0e-10 );

	assert( fabs( const_term[2]+dldx[2][0]*p0[0]+dldx[2][1]*p0[1] ) < 1.0e-10 );
	assert( fabs( const_term[2]+dldx[2][0]*p1[0]+dldx[2][1]*p1[1] ) < 1.0e-10 );
	assert( fabs( const_term[2]+dldx[2][0]*p2[0]+dldx[2][1]*p2[1] - 1.0 ) < 1.0e-10 );
*/
}
