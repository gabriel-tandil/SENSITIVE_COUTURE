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
@brief 線要素で要素剛性行列を作る際のUtility関数の集まり
@author Nobuyuki Umetani
*/


#if !defined(KER_EMAT_BAR_H)
#define KER_EMAT_BAR_H

#include <cassert>

const static unsigned int NIntLineGauss[4] = {
	1, 2, 3, 4
};
const static double LineGauss[4][4][2] = 
{
	{
		{ 0.0, 2.0 },
		{ 0.0, 0.0 },
		{ 0.0, 0.0 },
		{ 0.0, 0.0 },
	},
	{
		{ -0.577350269189626, 1.0 },
		{  0.577350269189626, 1.0 },
		{  0.0,               0.0 },
		{  0.0,               0.0 },
	},
	{
		{ -0.774596669241483, 0.555555555555556 },
		{  0.0,               0.888888888888889 },
		{  0.774596669241483, 0.555555555555556 },
		{  0.0,               0.0               },
	},
	{
		{ -0.861136311594053, 0.347854845137454 },
		{ -0.339981043584856, 0.652145154862546 },
		{  0.339981043584856, 0.652145154862546 },
		{  0.861136311594053, 0.347854845137454 },
	}
};



#endif
