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
@brief interface of field administration class (Fem::Field::CFieldWorld)
@author Nobuyuki Umetani
*/

#if !defined(FIELD_WORLD_H)
#define FIELD_WORLD_H

#include <map>

#include "delfem/elem_ary.h"	// need because reference of class "CElemAry" used
#include "delfem/node_ary.h"	// need because reference of class "CNodeAry" used
#include "delfem/field.h"
#include "delfem/objset.h"		// template for container with ID
#include "delfem/cad_com.h"		// need for enum CAD_ELEM_TYPE

namespace Msh{
	class IMesh;	
}

namespace Fem{
namespace Field{
  
/*! 
@brief ID converter between element array, mesh, CAD
@ingroup Fem
*/
class CIDConvEAMshCad
{
	friend class CFieldWorld;
public:
	void Clear(){
		m_aIdAry.clear();
	}
    bool IsIdEA(unsigned int id_ea) const{
        for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
            if( m_aIdAry[iid].id_ea == id_ea ) return true;
        }
        return false;
    }
	unsigned int GetIdEA_fromMsh(unsigned int id_part_msh) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if( m_aIdAry[iid].id_part_msh == id_part_msh ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	unsigned int GetIdEA_fromMshExtrude(unsigned int id_part_msh, unsigned int inum_ext) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
			if(    m_aIdAry[iid].id_part_msh_before_extrude == id_part_msh
                && m_aIdAry[iid].inum_extrude == inum_ext ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	// itype_cad_part : Vertex(0), Edge(1), Loop(2)
  unsigned int GetIdEA_fromCad(unsigned int id_part_cad, Cad::CAD_ELEM_TYPE itype_cad_part, unsigned int inum_ext = 0) const {
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){      
			if(   m_aIdAry[iid].id_part_cad    == id_part_cad 
         && m_aIdAry[iid].itype_part_cad == itype_cad_part 
         && m_aIdAry[iid].inum_extrude   == inum_ext ){
				return m_aIdAry[iid].id_ea;
			}
		}
		return 0;
	}
	void GetIdCad_fromIdEA(unsigned int id_ea, unsigned int& id_part_cad, Cad::CAD_ELEM_TYPE& itype_part_cad ) const{
		for(unsigned int iid=0;iid<m_aIdAry.size();iid++){
//      std::cout << iid << " " << id_ea << " " << m_aIdAry[iid].id_ea << " " << m_aIdAry[iid].id_part_cad << " " << m_aIdAry[iid].itype_part_cad << std::endl;
			if( m_aIdAry[iid].id_ea == id_ea ){        
				id_part_cad    = m_aIdAry[iid].id_part_cad;
				itype_part_cad = m_aIdAry[iid].itype_part_cad;
				return;
			}
		}
		id_part_cad    = 0;
		itype_part_cad = Cad::NOT_SET;
	}
private:
	class CInfoCadMshEA{
	public:
		unsigned int id_ea;
		unsigned int id_part_msh;
		unsigned int id_part_cad;
		Cad::CAD_ELEM_TYPE itype_part_cad;
		unsigned int id_part_msh_before_extrude;    // mesh id before extrude
		unsigned int inum_extrude;  // notextruded(0), ground(1), middle-face(2), top(3)
	};
	std::vector<CInfoCadMshEA> m_aIdAry;
};

/*! 
@brief field administration class
@ingroup Fem
*/
class CFieldWorld{
public:
	CFieldWorld();
  CFieldWorld(const CFieldWorld& world);
	~CFieldWorld();
  
  CFieldWorld& operator = (const CFieldWorld& world);
	
	// import mesh into FEM world
	unsigned int AddMesh(const Msh::IMesh& mesh);

	unsigned int SetCustomBaseField(unsigned int id_base,
		std::vector<unsigned int> aIdEA_Inc,
		std::vector< std::vector<int> >& aLnods,
		std::vector<unsigned int>& mapVal2Co);

	CIDConvEAMshCad GetIDConverter(unsigned int id_field_base) const {
    std::map<unsigned int,CIDConvEAMshCad>::const_iterator itr = m_map_field_conv.find(id_field_base);
    if( itr == m_map_field_conv.end() ){
      std::cout << "There is no base field Empty Conv " << std::endl;
      CIDConvEAMshCad conv;
      return conv;
    }
    return itr->second;
  }

	void Clear();	// Delate all field, elem_ary, node_ary
  
	////////////////////////////////////////////////////////////////
	// functionas for element array
	
	bool IsIdEA( unsigned int id_ea ) const;	// check if id_ea is a ID of element array
	const std::vector<unsigned int>& GetAry_IdEA() const;	// get all ID of element array
	const CElemAry& GetEA(unsigned int id_ea) const;	// get element array ( const )
	CElemAry& GetEA(unsigned int id_ea);	// get element array ( without-const )
	// add element array ( return id>0, return 0 if fail )
	unsigned int AddElemAry(unsigned int size, ELEM_TYPE elem_type);
	bool AddIncludeRelation(unsigned int id_ea, unsigned int id_ea_inc);	 // AddElemAryÇ∆ìùàÍÇµÇΩÇ¢ÅD

	////////////////////////////////////////////////////////////////
	// functions for node arary

	bool IsIdNA( unsigned int id_na ) const;	// check if id_na is a ID of node array
	const std::vector<unsigned int>& GetAry_IdNA() const;	// get all ID of node array
	const CNodeAry& GetNA(unsigned int id_na) const;	// Get Node Array ( const )
	CNodeAry& GetNA(unsigned int id_na);	// Get Node Array ( without-const )
	// add node array ( return id>0, return 0 if fail )
	unsigned int AddNodeAry(unsigned int size);
	
	////////////////////////////////////////////////////////////////
	// functions for field
	
	bool IsIdField( unsigned int id_field ) const;	// check if id_field is a ID of field
	const std::vector<unsigned int>& GetAry_IdField() const;	// get all ID of field
	const CField& GetField(unsigned int id_field) const;	// get field ( const )
	CField& GetField(unsigned int id_field);	// get field ( without field )

	////////////////////////////////////////////////////////////////
	// functions to add field

	unsigned int MakeField_FieldElemAry(unsigned int id_field, unsigned int id_ea, Field::FIELD_TYPE field_type = NO_VALUE, const int derivative_type = 1, const int node_configuration_type = 1 );
	unsigned int MakeField_FieldElemDim(unsigned int id_field, int idim_elem,      Field::FIELD_TYPE field_type = NO_VALUE, const int derivative_type = 1, const int node_configuration_type = 1 );
  unsigned int AddField(unsigned int id_field_parent,	// parent field
                        const std::vector<CField::CElemInterpolation>& aEI, 
                        const CField::CNodeSegInNodeAry& nsna_c, const CField::CNodeSegInNodeAry& nsna_b,
                        unsigned int id_field_candidate = 0);  
  
	//! Get partial field consists of element array with ID:id_ea
	unsigned int GetPartialField(unsigned int id_field, unsigned int IdEA );
	//! Get partial field consists of IDs of element array:aIdEA
	unsigned int GetPartialField(unsigned int id_field, std::vector<unsigned int> aIdEA);
  
  // update field
	bool UpdateMeshCoord(    const unsigned int id_base, const Msh::IMesh& mesh);
	bool UpdateMeshCoord(    const unsigned int id_base, const unsigned int id_field_disp, const Msh::IMesh& mesh);  
  ////
	bool UpdateConnectivity( const unsigned int id_base, const Msh::IMesh& mesh );
	bool UpdateConnectivity_CustomBaseField(const unsigned int id_base,
                                          const std::vector<unsigned int>& aIdEA_Inc, 
                                          const std::vector< std::vector<int> >& aLnods,
                                          const std::vector<unsigned int>& mapVal2Co);
	bool UpdateConnectivity_HingeField_Tri(unsigned int id_field, unsigned int id_field_base);
	bool UpdateConnectivity_EdgeField_Tri( unsigned int id_field, unsigned int id_field_base);
    
  // set value to field
//	void FieldValueExec(double time);
//	void FieldValueDependExec();
  
  // Delete Field and referenced EA and NA. the EAs and NAs that is referenced from Field in use are not deleted
  void DeleteField( const std::vector<unsigned int>& aIdFieldDel );

private:
	Com::CObjSet<CElemAry*> m_apEA;		//!< set of element array
	Com::CObjSet<CNodeAry*> m_apNA;		//!< set of node array
	Com::CObjSet<CField*> m_apField;	//!< set of field

  std::map<unsigned int,CIDConvEAMshCad> m_map_field_conv;
};

}	// end namespace Field
}	// end namespace Fem

#endif
