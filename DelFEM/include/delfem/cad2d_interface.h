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
@brief Interface for Msh::CMesher2D
@remarks A class derieve this class can be cut mesh using class (Msh::CMesher2D)
@author Nobuyuki Umetani
*/

#if !defined(CAD_2D_INTERFACE_H)
#define CAD_2D_INTERFACE_H

#include <memory>   // autoptr

#include "delfem/vector2d.h"
#include "delfem/cad_com.h"

namespace Cad{

/*! 
@ingroup CAD
@brief 2D CAD model class (Model class of 2D CAD for mesher)
*/
class ICad2D_Msh
{
public:	
  ICad2D_Msh(){}	//!< needs defalut constructor
  virtual ~ICad2D_Msh(){}	//!< virtual destructor is must for interface class
	virtual bool GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const = 0;  
	virtual bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const = 0;
	virtual const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE) const = 0;
  virtual int GetLayer(Cad::CAD_ELEM_TYPE, unsigned int id) const = 0;
	virtual void GetLayerMinMax(int& layer_min, int& layer_max) const = 0;
  virtual bool GetColor_Loop(unsigned int id_l, double color[3] ) const = 0;
	//! get area loop (ID:id_l)
	virtual double GetArea_Loop(unsigned int id_l) const = 0;
  virtual std::auto_ptr<IItrLoop> GetPtrItrLoop(unsigned int id_l) const = 0;
  virtual bool GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const = 0;
	virtual Com::CVector2D GetVertexCoord(unsigned int id_v) const = 0;
};

}	// end namespace CAD

#endif
