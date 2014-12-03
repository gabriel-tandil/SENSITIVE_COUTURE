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
@brief interface of vector drawer vector class (Fem::Field::View::CDrawerVector)
@author Nobuyuki Umetani
*/


#if !defined(DRAWER_FIELD_VECTOR_H)
#define DRAWER_FIELD_VECTOR_H

#include <memory>

#include "delfem/drawer_field.h"

namespace Fem{
namespace Field{
namespace View{
	
//! visualization class using vector
class CDrawerVector : public CDrawerField
{
public:
	CDrawerVector();
	CDrawerVector(unsigned int id_field, const Fem::Field::CFieldWorld& world );
	virtual ~CDrawerVector();
	Com::CBoundingBox3D GetBoundingBox( double rot[] ) const;
	virtual void DrawSelection(unsigned int idraw) const{};
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual void Draw() const;
	virtual bool Update(const Fem::Field::CFieldWorld& world);
private:
	bool Set(unsigned int id_field, const Fem::Field::CFieldWorld& world );
	virtual bool Update_VECTOR(const Fem::Field::CFieldWorld& world);	// update vector
	virtual bool Update_SSTR2(const Fem::Field::CFieldWorld& world);	// update principal stress
private:
	unsigned int id_field;
	unsigned int npo;	// number of base point
	unsigned int ndim_co;	// dimention of space
	unsigned int ndim_va;	// dimention of value
	unsigned int itype;	// 0:vector 1:sstr2
	
	/*
	 data : size = (ndim_co+ndim_va)*npo
	 ndim=2,itype=0	: cx,cy, vx,vy
	 ndim=3,itype=0	: cx,cy,cz, vx,vy,vz
	 ndim=2,itype=1	: cx,cy, pvsx,pvsy,flgs, pvlx,pvly,flgl
	*/
	double* pData;
};

}
}
}

#endif
