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
@brief interface of 2D editable mesh class (CMsher2D)
@author Nobuyuki Umetani
*/

#if !defined(MSHER_2D_EDIT_H)
#define MSHER_2D_EDIT_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <map>

#include "delfem/mesh_interface.h"
#include "delfem/cad_obj2d.h"
#include "delfem/vector2d.h"
#include "delfem/serialize.h"

#include "delfem/mesher2d.h"

////////////////////////////////////////////////

namespace Msh{

/*!
@ingroup Msh
*/
class CMesher2D_Edit : public CMesher2D
{
public:
	virtual bool Meshing(const Cad::CCadObj2D& cad_2d){
		return CMesher2D::Meshing(cad_2d);
	}
	virtual bool FitMeshToCad_UsingPrecomp(const Cad::CCadObj2D& cad_2d, 
                                         Cad::CAD_ELEM_TYPE itype_elem, unsigned int id_elem, 
                                         unsigned int& itype_ope ){
		if(      itype_elem == Cad::VERTEX ){ return this->FitMeshToCad_Vertex(cad_2d,id_elem,itype_ope); }
		else if( itype_elem == Cad::EDGE   ){ return this->FitMeshToCad_Edge(  cad_2d,id_elem,itype_ope); }
		else if( itype_elem == Cad::LOOP   ){ return this->FitMeshToCad_Loop(  cad_2d,id_elem,itype_ope); }
		return false;
	}
  virtual bool FitMeshToCad_Slider_UsingPrecomp(const Cad::CCadObj2D& cad_2d, unsigned int& itype_ope, double delta );
    
//	Cad::CAD_ELEM_TYPE GetMoveCadElemType() const { return move_cad_elem_type; }
//	unsigned int GetMoveCadElemID() const {   return move_cad_elem_id; }
	void Precomp_FitMeshToCad(const Cad::CCadObj2D& cad_2d, 
                            //                            Cad::CAD_ELEM_TYPE itype_elem, unsigned int id_elem, 
                            const std::vector<double>& aLamTnsr_Vtx);
  void GetXYHarmonicFunction( std::vector<double>& har) const { har = aLamTnsr; }
	////////////////
	virtual bool FitMeshToCad_Vertex(const Cad::CCadObj2D& cad_2d, unsigned int id_v_cad,  unsigned int& itype_operation );
	virtual bool FitMeshToCad_Edge(  const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad,  unsigned int& itype_operation );
	virtual bool FitMeshToCad_Loop(  const Cad::CCadObj2D& cad_2d, unsigned int id_l_cad,  unsigned int& itype_operation );
	//! Try fitting all the point (few arg but taking time)
	virtual bool FitMeshToCad_All(   const Cad::CCadObj2D& cad_2d, unsigned int& itype_operation );
private:
	// for precomputation
//	Cad::CAD_ELEM_TYPE move_cad_elem_type;
//	unsigned int move_cad_elem_id;
	std::vector<double> aLamTnsr;
  std::vector<double> aLamTnsr_Vtx;
protected:
	void SmoothingMesh_Laplace(unsigned int num_iter, double elen);

	// locate a end point (id_v_cad_mov) to the distination (dist_mov)
	bool MovePointsOnEdge(
		const Cad::CCadObj2D& cad_2d,unsigned int id_e_cad_mov,  
		std::vector<int>& aflg_ismov,
		unsigned int id_v_cad_mov, const Com::CVector2D& dist_mov );

	// set end poins of the edge to distination (dist_mov)
	bool MovePointsOnEdge(
		const Cad::CCadObj2D& cad_2d,unsigned int id_e_cad_mov,  
		std::vector<int>& aflg_ismov,
		const Com::CVector2D& dist_mov_s, const Com::CVector2D& dist_mov_e );

	// fit end points of the edge to cad
	bool MovePointsOnEdge(
		const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad_mov, 
		std::vector<int>& aflg_ismov );
	
	bool LambdaEdge
  (const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad,
   const std::vector<double>& aValVtx,
   std::vector<int>& aflg_ismov);
  
	bool LambdaLoop_MVC
  (const Cad::CCadObj2D& cad_2d, unsigned int id_l_cad,
   const std::vector<double>& aValVtx,
   std::vector<int>& aflg_ismov);
};


}   // end namespace

#endif
