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
@brief 六面体要素で要素剛性行列を作る際のUtility関数の集まり
@author Nobuyuki Umetani
*/

#if !defined(KER_EMAT_HEX_H)
#define KER_EMAT_HEX_H

#include <cassert>

#include "ker_emat_bar.h"

/*!
@brief 積分点における形状関数の値とその微分を作る関数
*/
static inline void ShapeFunc_Hex8(
	const double& r1, const double& r2,	const double& r3,	// (入力)積分点の自然座標における位置
	const double coords[][3],			// (入力)節点座標
	double& detjac,		// 積分点におけるヤコビアンの値
	double dndx[][3],	// 積分点における形状関数の微分値
	double an[] )		// 積分点における形状関数の値
{
	unsigned int inode;

	an[0] = 0.125*(1.0-r1)*(1.0-r2)*(1.0-r3);
	an[1] = 0.125*(1.0+r1)*(1.0-r2)*(1.0-r3);
	an[2] = 0.125*(1.0+r1)*(1.0+r2)*(1.0-r3);
	an[3] = 0.125*(1.0-r1)*(1.0+r2)*(1.0-r3);
	an[4] = 0.125*(1.0-r1)*(1.0-r2)*(1.0+r3);
	an[5] = 0.125*(1.0+r1)*(1.0-r2)*(1.0+r3);
	an[6] = 0.125*(1.0+r1)*(1.0+r2)*(1.0+r3);
	an[7] = 0.125*(1.0-r1)*(1.0+r2)*(1.0+r3);

	double dndr[8][3];
	dndr[0][0] = -0.125*(1.0-r2)*(1.0-r3);
	dndr[1][0] = -dndr[0][0];
	dndr[2][0] =  0.125*(1.0+r2)*(1.0-r3);
	dndr[3][0] = -dndr[2][0];
	dndr[4][0] = -0.125*(1.0-r2)*(1.0+r3);
	dndr[5][0] = -dndr[4][0];
	dndr[6][0] =  0.125*(1.0+r2)*(1.0+r3);
	dndr[7][0] = -dndr[6][0];

	dndr[0][1] = -0.125*(1.0-r1)*(1.0-r3);
	dndr[1][1] = -0.125*(1.0+r1)*(1.0-r3);
	dndr[2][1] = -dndr[1][1];
	dndr[3][1] = -dndr[0][1];
	dndr[4][1] = -0.125*(1.0-r1)*(1.0+r3);
	dndr[5][1] = -0.125*(1.0+r1)*(1.0+r3);
	dndr[6][1] = -dndr[5][1];
	dndr[7][1] = -dndr[4][1];

	dndr[0][2] = -0.125*(1.0-r1)*(1.0-r2);
	dndr[1][2] = -0.125*(1.0+r1)*(1.0-r2);
	dndr[2][2] = -0.125*(1.0+r1)*(1.0+r2);
	dndr[3][2] = -0.125*(1.0-r1)*(1.0+r2);
	dndr[4][2] = -dndr[0][2];
	dndr[5][2] = -dndr[1][2];
	dndr[6][2] = -dndr[2][2];
	dndr[7][2] = -dndr[3][2];

	double dxdr[3][3]  = {
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
		{ 0.0, 0.0, 0.0 },
	};

	for(inode=0;inode<8;inode++){
		dxdr[0][0] += coords[inode][0]*dndr[inode][0];
		dxdr[0][1] += coords[inode][0]*dndr[inode][1];
		dxdr[0][2] += coords[inode][0]*dndr[inode][2];
		dxdr[1][0] += coords[inode][1]*dndr[inode][0];
		dxdr[1][1] += coords[inode][1]*dndr[inode][1];
		dxdr[1][2] += coords[inode][1]*dndr[inode][2];
		dxdr[2][0] += coords[inode][2]*dndr[inode][0];
		dxdr[2][1] += coords[inode][2]*dndr[inode][1];
		dxdr[2][2] += coords[inode][2]*dndr[inode][2];
	}

	detjac = dxdr[0][0]*dxdr[1][1]*dxdr[2][2] 
		   + dxdr[1][0]*dxdr[2][1]*dxdr[0][2]
		   + dxdr[2][0]*dxdr[0][1]*dxdr[1][2]
		   - dxdr[0][0]*dxdr[2][1]*dxdr[1][2]
		   - dxdr[1][0]*dxdr[0][1]*dxdr[2][2]
		   - dxdr[2][0]*dxdr[1][1]*dxdr[0][2];

	const double inv_jac = 1.0 / detjac;

	double drdx[3][3];
	drdx[0][0] = inv_jac*( dxdr[1][1]*dxdr[2][2]-dxdr[1][2]*dxdr[2][1] );
	drdx[0][1] = inv_jac*( dxdr[0][2]*dxdr[2][1]-dxdr[0][1]*dxdr[2][2] );
	drdx[0][2] = inv_jac*( dxdr[0][1]*dxdr[1][2]-dxdr[0][2]*dxdr[1][1] );
	drdx[1][0] = inv_jac*( dxdr[1][2]*dxdr[2][0]-dxdr[1][0]*dxdr[2][2] );
	drdx[1][1] = inv_jac*( dxdr[0][0]*dxdr[2][2]-dxdr[0][2]*dxdr[2][0] );
	drdx[1][2] = inv_jac*( dxdr[0][2]*dxdr[1][0]-dxdr[0][0]*dxdr[1][2] );
	drdx[2][0] = inv_jac*( dxdr[1][0]*dxdr[2][1]-dxdr[1][1]*dxdr[2][0] );
	drdx[2][1] = inv_jac*( dxdr[0][1]*dxdr[2][0]-dxdr[0][0]*dxdr[2][1] );
	drdx[2][2] = inv_jac*( dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0] );

	for(inode=0;inode<8;inode++){
		dndx[inode][0] = dndr[inode][0]*drdx[0][0] + dndr[inode][1]*drdx[1][0] + dndr[inode][2]*drdx[2][0];
		dndx[inode][1] = dndr[inode][0]*drdx[0][1] + dndr[inode][1]*drdx[1][1] + dndr[inode][2]*drdx[2][1];
		dndx[inode][2] = dndr[inode][0]*drdx[0][2] + dndr[inode][1]*drdx[1][2] + dndr[inode][2]*drdx[2][2];
	}
}

#endif
