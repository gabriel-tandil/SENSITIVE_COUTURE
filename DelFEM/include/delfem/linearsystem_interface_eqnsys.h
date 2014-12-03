/*
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
@brief 連立一次方程式の方程式へのインターフェスクラス(Fem::Eqn::CLinearSystem_EqnInterface)
@author Nobuyuki Umetani
*/

#if !defined(EQUATION_INTERFACE_H)
#define EQUATION_INTERFACE_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <assert.h>

#include "delfem/elem_ary.h"   // CORNER, EDDGE, BUBBLEの名前のために追加（何とかしてこのincludeを削除したい）

namespace MatVec{
    class CMatDia_BlkCrs;
    class CMat_BlkCrs;
    class CVector_Blk;
}

namespace Fem{
namespace Field{
    class CFieldWorld;
}
namespace Eqn{

//! interface class of linear system for equation
class ILinearSystem_Eqn
{
public:
    virtual MatVec::CMat_BlkCrs& GetMatrix(
        unsigned int id_field1,
        Fem::Field::ELSEG_TYPE node_config1,
        unsigned int id_field2,
        Fem::Field::ELSEG_TYPE node_config2,
        const Fem::Field::CFieldWorld& world) = 0;

    virtual MatVec::CMatDia_BlkCrs& GetMatrix(
        unsigned int id_field_disp,
        Fem::Field::ELSEG_TYPE node_config,
        const Fem::Field::CFieldWorld& world) = 0;

    virtual MatVec::CVector_Blk& GetResidual(
        unsigned int id_field_disp,
        Fem::Field::ELSEG_TYPE node_config,
        const Fem::Field::CFieldWorld& world) = 0;
public:
};
	
}	// Eqn
}	// Fem

#endif
