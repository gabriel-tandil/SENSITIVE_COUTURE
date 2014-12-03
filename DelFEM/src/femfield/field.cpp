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

////////////////////////////////////////////////////////////////
// Field.cpp：場クラス(CField)の実装
////////////////////////////////////////////////////////////////


#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#define for if(0);else for

#include <math.h>
#include <fstream>
#include <stdio.h>

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/eval.h"
#include "delfem/femeqn/ker_emat_hex.h"
#include "delfem/femeqn/ker_emat_tet.h"

using namespace Fem::Field;
using namespace MatVec;

static double TriArea(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

CField::CField(const CField& rhs)
{
	m_is_valid = rhs.m_is_valid;
  m_id_field_parent = rhs.m_id_field_parent;
	m_ndim_coord = rhs.m_ndim_coord;
	m_aElemIntp = rhs.m_aElemIntp;
  
	m_na_c = rhs.m_na_c;
	m_na_e = rhs.m_na_e;
	m_na_b = rhs.m_na_b;
  
	m_field_type = rhs.m_field_type;
	m_field_derivative_type = rhs.m_field_derivative_type;
	m_map_val2co = rhs.m_map_val2co;
  
	////////////////
	m_DofSize = rhs.m_DofSize;
//	m_aValueFieldDof = rhs.m_aValueFieldDof;  
//	m_id_field_dep = rhs.m_id_field_dep;
//	m_is_gradient = rhs.m_is_gradient;
}

CField::CField
(unsigned int id_field_parent,	// parent field
 const std::vector<CElemInterpolation>& aEI, 
 const CNodeSegInNodeAry& nsna_c, const CNodeSegInNodeAry& nsna_b, 
 CFieldWorld& world)
{
	m_aElemIntp = aEI;

	m_na_c = nsna_c;
	m_na_b = nsna_b;

	m_ndim_coord = 0;
	if( m_na_c.id_na_co ){
		const CNodeAry& na = world.GetNA( m_na_c.id_na_co );
    assert( na.IsSegID( m_na_c.id_ns_co ) );
		const CNodeAry::CNodeSeg& ns_co = na.GetSeg( m_na_c.id_ns_co );
		m_ndim_coord = ns_co.Length();
	}
	else if( m_na_b.id_na_co ){
		const CNodeAry& na = world.GetNA( m_na_b.id_na_co );
    assert( na.IsSegID( m_na_b.id_ns_co ) );
		const CNodeAry::CNodeSeg& ns_co = na.GetSeg( m_na_b.id_ns_co );
		m_ndim_coord = ns_co.Length();
	}
  
//  std::cout << "field : " << m_na_b.id_na_co << " " << m_na_b.id_ns_co << std::endl;

	m_field_type = NO_VALUE;
	m_field_derivative_type = 0;
	m_DofSize = 0;

	m_id_field_parent = id_field_parent;
	if( m_id_field_parent == 0 ){	
		assert( m_na_c.is_part_va == false );
	}
	else{
		assert( world.IsIdField(m_id_field_parent) );
		const CField& field_parent = world.GetField(m_id_field_parent);
		this->m_DofSize = field_parent.GetNLenValue();
		this->m_field_type = field_parent.GetFieldType();
		this->m_field_derivative_type = field_parent.GetFieldDerivativeType();
	}

	// map_val2co
	if(    m_na_c.id_na_va != 0 
		&& m_na_c.id_na_va != m_na_c.id_na_co
		&& m_na_c.id_na_co != 0 )
	{
		const CNodeAry& na_va = world.GetNA(m_na_c.id_na_va);
		const unsigned int nnode_va = na_va.Size();
		m_map_val2co.resize(nnode_va);
		for(unsigned int iei=0;iei<m_aElemIntp.size();iei++){
			const CElemInterpolation& ei = m_aElemIntp[iei];
			const CElemAry& ea = world.GetEA(ei.id_ea);
			const CElemAry::CElemSeg& es_co = ea.GetSeg(ei.id_es_c_co);
			const CElemAry::CElemSeg& es_val = ea.GetSeg(ei.id_es_c_va);
			unsigned int noes_co[10];
			unsigned int noes_val[10];
			unsigned int nnoes = es_co.Length();
			assert( es_co.Size() == es_val.Size() );
			for(unsigned int ielem=0;ielem<ea.Size();ielem++){
				es_co.GetNodes(ielem,noes_co);
				es_val.GetNodes(ielem,noes_val);
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					m_map_val2co[ noes_val[inoes] ] = noes_co[inoes];
				}
			}
		}
	}

	m_is_valid = this->AssertValid(world);
	assert( m_is_valid );
}


bool CField::AssertValid(CFieldWorld& world) const
{
	{	// fieldが持っているeaが正しくnaを参照してるか調べる
		for(unsigned int iei=0;iei<m_aElemIntp.size();iei++){
			const Field::CField::CElemInterpolation& ei = m_aElemIntp[iei];
			unsigned int id_ea = ei.id_ea;
			assert( world.IsIdEA(id_ea) );
			const CElemAry& ea = world.GetEA(id_ea);
			if( ei.id_es_c_co != 0 ){
				unsigned int id_es_c_coord = ei.id_es_c_co;
				assert( ea.IsSegID(id_es_c_coord) );
				const CElemAry::CElemSeg& es = ea.GetSeg(id_es_c_coord);
				const unsigned int id_na_co = es.GetIdNA();
				assert( this->m_na_c.id_na_co == id_na_co );
				assert( world.IsIdNA(id_na_co) );
				const CNodeAry& na_co = world.GetNA(id_na_co);
				assert( es.GetMaxNoes() < na_co.Size() );
			}
			if( ei.id_es_c_va != 0 ){
				unsigned int id_es_c_va = ei.id_es_c_va;
				assert( ea.IsSegID(id_es_c_va) );
				const CElemAry::CElemSeg& es_va = ea.GetSeg(id_es_c_va);
				unsigned int id_na_va = es_va.GetIdNA();
				assert( this->m_na_c.id_na_va == id_na_va );
				assert( world.IsIdNA(id_na_va) );
				const CNodeAry& na_va = world.GetNA(id_na_va);
				assert( es_va.GetMaxNoes() < na_va.Size() );
			}
		}
	}
	{	// m_dim_coordを調べる
		if( world.IsIdNA(m_na_c.id_na_co) ){
			assert( world.IsIdNA(m_na_c.id_na_co) );
			const CNodeAry& na_c_coord = world.GetNA(m_na_c.id_na_co);
			assert( na_c_coord.IsSegID(m_na_c.id_ns_co) );
			const CNodeAry::CNodeSeg& ns_c_coord = na_c_coord.GetSeg(m_na_c.id_ns_co);
			assert( this->m_ndim_coord == ns_c_coord.Length() );
		}
		// TODO : bubbleやedgeについても調べたい
	}
	{	// cornerのvalueが部分場かどうか示すフラグが合ってるかどうか調べる
		if( m_na_c.id_na_va != 0 ){
			assert( world.IsIdNA(m_na_c.id_na_va) );
			CNodeAry& na_val = world.GetNA(m_na_c.id_na_va);
			unsigned int nnode_val = na_val.Size();
			std::vector<unsigned int> flag_vec( nnode_val, 0 );
			////////////////
			for(unsigned int iei=0;iei<m_aElemIntp.size();iei++){
				const unsigned int id_ea = m_aElemIntp[iei].id_ea;
				const CElemAry& ea = world.GetEA(id_ea);
				unsigned int id_es_val = m_aElemIntp[iei].id_es_c_va;
				assert( ea.IsSegID(id_es_val) );
				const CElemAry::CElemSeg& es_val = ea.GetSeg(id_es_val);
				assert( es_val.GetMaxNoes() < nnode_val );
				unsigned int nnoes = es_val.Length();
				unsigned int nno_c[64];
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es_val.GetNodes(ielem,nno_c);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						assert( nno_c[inoes] < nnode_val );
						flag_vec[ nno_c[inoes] ] = 1;
					}
				}
			}
			bool iflag1 = false;
			for(unsigned int ino=0;ino<nnode_val;ino++){
				if( flag_vec[ino] == 0 ){ iflag1 = true; break; }
			}
//			assert( m_na_c.is_part_va == iflag1 );
		}
	}
	{	// cornerのcoordが部分場かどうか示すフラグが合ってるかどうか調べる
		if( m_na_c.id_na_co != 0 ){
			assert( world.IsIdNA(m_na_c.id_na_co) );
			CNodeAry& na = world.GetNA(m_na_c.id_na_co);
			unsigned int nnode_co = na.Size();
			std::vector<unsigned int> flag_vec( nnode_co, 0 );
			////////////////
			for(unsigned int iei=0;iei<m_aElemIntp.size();iei++){
				const unsigned int id_ea = m_aElemIntp[iei].id_ea;
				assert( world.IsIdEA(id_ea) );
				const CElemAry& ea = world.GetEA(id_ea);
				unsigned int id_es_co = m_aElemIntp[iei].id_es_c_co;
				assert( ea.IsSegID(id_es_co) );
				const CElemAry::CElemSeg& es_co = ea.GetSeg(id_es_co);
				assert( es_co.GetMaxNoes() < nnode_co );
				unsigned int nnoes = es_co.Length();
				unsigned int no_c[64];
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					es_co.GetNodes(ielem,no_c);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						assert( no_c[inoes] < nnode_co );
						flag_vec[ no_c[inoes] ] = 1;
					}
				}
			}
			bool iflag1 = false;
			for(unsigned int ino=0;ino<nnode_co;ino++){
				if( flag_vec[ino] == 0 ){ iflag1 = true; break; }
			}
//			assert( m_na_c.is_part_co == iflag1 );
		}
	}
	return true;
}

const CNodeAry::CNodeSeg& CField::GetNodeSeg
(ELSEG_TYPE elseg_type, 
 bool is_value, const CFieldWorld& world, unsigned int fdt ) const
{
	unsigned int id_na = 0;
	unsigned int id_ns = 0;
	if(      elseg_type == CORNER ){
		if( is_value ){
			id_na = this->m_na_c.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_c.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ac; }
		}
		else{
			id_na = this->m_na_c.id_na_co;
			id_ns = this->m_na_c.id_ns_co;
		}
	}
	else if( elseg_type == BUBBLE ){
		if( is_value ){
			id_na = this->m_na_b.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_b.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ac; }
		}
		else{
			id_na = this->m_na_b.id_na_co;
			id_ns = this->m_na_b.id_ns_co;
		}
	}
//  std::cout << "IsNodeSeg : " << id_na << " " << id_ns << std::endl;
	assert( world.IsIdNA(id_na) );
	const CNodeAry& na = world.GetNA(id_na);
	assert( na.IsSegID(id_ns) );
	return na.GetSeg(id_ns);

}

bool CField::IsNodeSeg(ELSEG_TYPE elseg_type,
					   bool is_value, const CFieldWorld& world, unsigned int fdt) const
{
	unsigned int id_na = 0;
	unsigned int id_ns = 0;
	if(      elseg_type == CORNER ){
		if( is_value ){
			id_na = this->m_na_c.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_c.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ac; }
		}
		else{
			id_na = this->m_na_c.id_na_co;
			id_ns = this->m_na_c.id_ns_co;
		}
	}
	else if( elseg_type == BUBBLE ){
		if( is_value ){
			id_na = this->m_na_b.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_b.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ac; }
		}
		else{
			id_na = this->m_na_b.id_na_co;
			id_ns = this->m_na_b.id_ns_co;
		}
	}
	if( !world.IsIdNA(id_na) ) return false;
	const CNodeAry& na = world.GetNA(id_na);
	if( !na.IsSegID(id_ns) ) return false;
	return true;
}


CNodeAry::CNodeSeg& CField::GetNodeSeg(
	ELSEG_TYPE elseg_type, 
	bool is_value, CFieldWorld& world, unsigned int fdt )
{
	unsigned int id_na = 0;
	unsigned int id_ns = 0;
	if(      elseg_type == CORNER ){
		if( is_value ){
			id_na = this->m_na_c.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_c.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_c.id_ns_ac; }
		}
		else{
			id_na = this->m_na_c.id_na_co;
			id_ns = this->m_na_c.id_ns_co;
		}
	}
	else if( elseg_type == BUBBLE ){
		if( is_value ){
			id_na = this->m_na_b.id_na_va;
			if( fdt & 1               ){ id_ns = this->m_na_b.id_ns_va; }
			if( fdt & 2 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ve; }
			if( fdt & 4 && id_ns == 0 ){ id_ns = this->m_na_b.id_ns_ac; }
		}
		else{
			id_na = this->m_na_b.id_na_co;
			id_ns = this->m_na_b.id_ns_co;
		}
	}
	assert( world.IsIdNA(id_na) );
	CNodeAry& na = world.GetNA(id_na);
	assert( na.IsSegID(id_ns) );
	return na.GetSeg(id_ns);

}

const CElemAry::CElemSeg& CField::GetElemSeg(unsigned int id_ea, ELSEG_TYPE elseg_type, bool is_value, const CFieldWorld& world ) const
{
	unsigned int id_es = this->GetIdElemSeg(id_ea,elseg_type,is_value,world);
	if( id_es == 0 ){// try partiel element
    const CElemAry& ea = world.GetEA(id_ea);
    unsigned int id_na_co = this->GetNodeSegInNodeAry(CORNER).id_na_co;
    unsigned int id_na_va = this->GetNodeSegInNodeAry(CORNER).id_na_va;    
    const std::vector<unsigned int>& aIdES = ea.GetAry_SegID();
    for(unsigned int iies=0;iies<aIdES.size();iies++){
      unsigned int id_es0 = aIdES[iies];
      const CElemAry::CElemSeg& es = ea.GetSeg(id_es0);
      if( is_value ){ if( es.GetIdNA() == id_na_va ){ id_es = id_es0; break; } }
      else{           if( es.GetIdNA() == id_na_co ){ id_es = id_es0; break; } }
    }
    if( id_es == 0 ){
      assert(0);
      throw 0;
    }
	}
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.IsSegID(id_es) );
	return ea.GetSeg(id_es);
}

unsigned int CField::GetIdElemSeg(unsigned int id_ea, ELSEG_TYPE elseg_type, bool is_value, const CFieldWorld& world ) const
{
	if( !world.IsIdEA(id_ea) ){ return false; }
	unsigned int iei;
	for(iei=0;iei<m_aElemIntp.size();iei++){
		if( m_aElemIntp[iei].id_ea == id_ea ){ break; }
	}
	if( iei == m_aElemIntp.size() ){ return 0; }
	const CElemInterpolation& ei = m_aElemIntp[iei];
	unsigned int id_es = 0;
	if( is_value ){
		if(      elseg_type == CORNER ){ id_es = ei.id_es_c_va; }
		else if( elseg_type == BUBBLE ){ id_es = ei.id_es_b_va; }
		else if( elseg_type == EDGE   ){ id_es = ei.id_es_e_va; }
		else{ assert(0); }
	}
	else{
		if(      elseg_type == CORNER ){ id_es = ei.id_es_c_co; }
		else if( elseg_type == BUBBLE ){ id_es = ei.id_es_b_co; }
		else if( elseg_type == EDGE   ){ id_es = ei.id_es_e_co; }
		else{ assert(0); }
	}
	const CElemAry& ea = world.GetEA(id_ea);
	if( !ea.IsSegID(id_es) ) return 0;
	return id_es;
}


int CField::GetLayer(unsigned int id_ea) const
{
	unsigned int iei;
	for(iei=0;iei<m_aElemIntp.size();iei++){
		if( m_aElemIntp[iei].id_ea == id_ea ){ break; }
	}
	if( iei == m_aElemIntp.size() ){ return 0; }
	const CElemInterpolation& ei = m_aElemIntp[iei];
	return ei.ilayer;
}

// 最大値最小値を取得する関数
void CField::GetMinMaxValue(double& min, double& max, const CFieldWorld& world, 
							unsigned int idof, int fdt) const
{
	min = 0; max = 0;
	const unsigned int ndof = this->GetNLenValue();
	assert( idof < ndof );
	if( idof >= ndof ) return;
	double *pVal = new double [ndof];
	if( this->m_na_c.id_na_va != 0 ){
		assert( this->IsNodeSeg(CORNER,true,world,fdt) );
		const CNodeAry::CNodeSeg& ns_c = this->GetNodeSeg(CORNER,true,world,fdt);
		{
			ns_c.GetValue(0,pVal);
			const double val = pVal[idof];
			min = val; max = val;
		}
		const unsigned int nnode = ns_c.Size();
		for(unsigned int inode=0;inode<nnode;inode++){
			ns_c.GetValue(inode,pVal);
			double val = pVal[idof];
			min = ( val < min ) ? val : min;
			max = ( val > max ) ? val : max;
		}
	}
	if( this->m_na_b.id_na_va != 0 ){
		assert( this->IsNodeSeg(BUBBLE,true,world,fdt) );
		const CNodeAry::CNodeSeg& ns_b = this->GetNodeSeg(BUBBLE,true,world,fdt);
		if( this->m_na_c.id_na_va == 0 ){
			ns_b.GetValue(0,pVal);
			const double val = pVal[idof];
			min = val; max = val;
		}
		for(unsigned int inode=0;inode<ns_b.Size();inode++){
			ns_b.GetValue(inode,pVal);
			const double val = pVal[idof];
			min = ( val < min ) ? val : min;
			max = ( val > max ) ? val : max;
		}
	}
/*	if( this->GetFieldDerivativeType() & fdt ){
		const CNodeAry::CNodeSeg& ns_c = this->GetNodeSeg(CORNER,true,world,fdt);
		ns_c.GetValue(0,val);
		min = val; max = val;
		for(unsigned int inode=0;inode<ns_c.GetNnode();inode++){
			ns_c.GetValue(0,val);
			min = ( val < min ) ? val : min;
			max = ( val > max ) ? val : max;
		}
	}*/
	delete pVal;
}


// 場の追加
// 古い場を見て足りないNSだけを追加するようにしたい．
bool CField::SetValueType( FIELD_TYPE field_type, const int fdt, CFieldWorld& world)
{
	if( this->GetIDFieldParent() != 0 ){
		assert(0);
		return false;	// 親フィールドのみに有効
	}

	// CORNER節点についてNAやNSを追加する
	if( this->m_na_c.id_na_va != 0 ){
		Fem::Field::CField::CNodeSegInNodeAry& nsna = this->m_na_c;
		unsigned int id_na = nsna.id_na_va;
		assert( world.IsIdNA(id_na) );
		CNodeAry& na = world.GetNA(id_na);
		std::vector<unsigned int> tmp_id_ns_ary = na.GetFreeSegID(3);
		std::vector< std::pair<unsigned int,CNodeAry::CNodeSeg> > add_ns;
		if( fdt&VALUE ){
			nsna.id_ns_va = tmp_id_ns_ary[0];
			if(      field_type == SCALAR  ){
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(1,"SCAL_VAL")) ); 
			}
			else if( field_type == VECTOR2 ){
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(2,"VEC2_VAL")) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"VEC3_VAL")) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_VAL")) );
			}
			else if( field_type == ZSCALAR ){	
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(2,"ZSCAL_VAL")) ); 
			}
			else{
				assert(0);
			}
		}
		if( fdt&VELOCITY ){
			nsna.id_ns_ve = tmp_id_ns_ary[1];
			if(      field_type == SCALAR  ){	
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(1,"SCAL_VELO")) ); 
			}
			else if( field_type == VECTOR2 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(2,"VEC2_VELO")) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(3,"VEC3_VELO")) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_VELO")) );
			}
		}
		if( fdt&ACCELERATION ){
			nsna.id_ns_ac = tmp_id_ns_ary[2];
			if(      field_type == SCALAR  ){	
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(1,"SCAL_ACC")) ); 
			}
			else if( field_type == VECTOR2 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(2,"VEC2_ACC")) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(3,"VEC3_ACC")) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_ACC")) );
			}
		}
		std::vector<int> id_ary = na.AddSegment( add_ns, 0.0 );
	}
	// BUBBLE節点についてNAとNSを追加する
	if( this->m_na_b.id_na_va != 0 ){
		Fem::Field::CField::CNodeSegInNodeAry& nsna = this->m_na_b;
//		unsigned int id_na = nsna.id_na_va;
		assert( world.IsIdNA(nsna.id_na_va) );
		CNodeAry& na = world.GetNA(nsna.id_na_va);
		std::vector<unsigned int> tmp_id_ns_ary = na.GetFreeSegID(3);
		std::vector< std::pair<unsigned int,CNodeAry::CNodeSeg> > add_ns;
		if( fdt & VALUE ){
			nsna.id_ns_va = tmp_id_ns_ary[0];
			if( field_type == SCALAR ){	
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(1,"SCAL_VAL")) ); 
			}
			else if( field_type == VECTOR2 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(2,"VEC2_VAL")) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"VEC3_VAL")) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_VELO")) );
			}
		}
		if( fdt & VELOCITY ){
			nsna.id_ns_ve = tmp_id_ns_ary[1];
			if(      field_type == SCALAR  ){
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(1,"SCALAR_VELO")) ); 
			}
			else if( field_type == VECTOR2 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(2,"VEC2_VELO"  )) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ve, CNodeAry::CNodeSeg(3,"VEC3_VELO"  )) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_VELO" )) );
			}
		}
		if( fdt & ACCELERATION ){
			nsna.id_ns_ac = tmp_id_ns_ary[2];
			if( field_type == SCALAR ){	
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(1,"SCAL_ACC" )) ); 
			}
			else if( field_type == VECTOR2 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(2,"VEC2_ACC" )) );	
			}
			else if( field_type == VECTOR3 ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_ac, CNodeAry::CNodeSeg(3,"VEC3_ACC" )) );	
			}
			else if( field_type == STSR2   ){ 
				add_ns.push_back( std::make_pair(nsna.id_ns_va, CNodeAry::CNodeSeg(3,"STSR2_ACC")) );
			}
		}
		std::vector<int> id_ary = na.AddSegment( add_ns, 0.0 );
	}

	this->m_field_type = field_type;
	this->m_field_derivative_type = fdt;

	unsigned int ndofsize = 0;
	{
		if(      field_type == SCALAR  ) ndofsize = 1;
		else if( field_type == VECTOR2 ) ndofsize = 2;
		else if( field_type == VECTOR3 ) ndofsize = 3;
		else if( field_type == STSR2   ) ndofsize = 3;
		else if( field_type == ZSCALAR ) ndofsize = 2;
		else{ assert(0); }
	}
	this->m_DofSize = ndofsize;
	return true;
}

// 要素補間のタイプを取得する
INTERPOLATION_TYPE CField::GetInterpolationType(unsigned int id_ea,const CFieldWorld& world) const{
	unsigned int iei;
	for(iei=0;iei<m_aElemIntp.size();iei++){
		if( m_aElemIntp[iei].id_ea == id_ea ) break;
	}
	assert( iei != m_aElemIntp.size() );
	if( iei == m_aElemIntp.size() ) throw;
	CElemInterpolation ei = m_aElemIntp[iei];
	ELEM_TYPE elem_type;
	{
		const unsigned int id_ea = ei.id_ea;
		assert( world.IsIdEA(id_ea) );
		const CElemAry& ea = world.GetEA(id_ea);
		elem_type = ea.ElemType();
	}
	if( elem_type == LINE ){
		if( ei.id_es_c_va != 0 && ei.id_es_b_va == 0 && ei.id_es_e_va==0 ){	     return LINE11; }	// 直線一次補間
		else{
			std::cout << "Error!-->Interpolation Not Defined(Tri)" << std::endl;
			std::cout << ei.id_es_c_va << " " << ei.id_es_b_va << std::endl;
			assert(0);
		}
	}
	else if( elem_type == TRI ){
		if(      ei.id_es_c_va != 0 && ei.id_es_b_va == 0 && ei.id_es_e_va==0 ){ return TRI11;   }	// 三角形一次補間
		else if( ei.id_es_c_va != 0 && ei.id_es_b_va != 0 && ei.id_es_e_va==0 ){ return TRI1011; }	// 三角形バブル補間(サブパラメトリック)
		else if( ei.id_es_c_va == 0 && ei.id_es_b_va != 0 && ei.id_es_e_va==0 ){ return TRI1001; }	// 三角形要素一定補間
		else{
			std::cout << "Error!-->Interpolation Not Defined(Tri)" << std::endl;
			std::cout << ei.id_es_c_va << " " << ei.id_es_b_va << std::endl;
			assert(0);
		}
	}
	else if( elem_type == TET ){
		if(      ei.id_es_c_va!=0 && ei.id_es_b_va == 0 && ei.id_es_e_va==0 ){ return TET11;   }	// 四面体一次補間
		else if( ei.id_es_c_va==0 && ei.id_es_b_va != 0 && ei.id_es_e_va==0 ){ return TET1001; }	// 四面体要素一定補間
	}
	else if( elem_type == HEX ){
		if(      ei.id_es_c_va!=0 && ei.id_es_b_va==0 && ei.id_es_e_va==0 ){ return HEX11;   }	// 六面体一次補間
		else if( ei.id_es_c_va==0 && ei.id_es_b_va!=0 && ei.id_es_e_va==0 ){ return HEX1001; }	// 六面体要素一定補間
	}
	std::cout << iei << "   " << ei.id_es_c_va << " " << ei.id_es_b_va << " " << ei.id_es_e_va << std::endl;
	std::cout << "Error!-->Not Implimented : " << elem_type << std::endl;
	assert(0);
	return TRI11;
}

// 三角形の面積を求める関数
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}

// TODO: 一度この関数を呼んだら，次に呼んだ時は高速に処理されるように，ハッシュを構築する
// TODO: 座標やコネクティビティの変更があった場合は，ハッシュを削除する
bool CField::FindVelocityAtPoint(double velo[],  
	unsigned int& id_ea_stat, unsigned int& ielem_stat, double& r1, double& r2,
	const double co[], const Fem::Field::CFieldWorld& world) const 
{
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = this->GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = this->GetNodeSeg(CORNER,false,world,VELOCITY);
	assert( ns_v.Length() == 2 );
	assert( ns_c.Length() == 2 );
	for(unsigned int iiea=0;iiea<m_aElemIntp.size();iiea++)
	{
		unsigned int id_ea = m_aElemIntp[iiea].id_ea;
		const Fem::Field::CElemAry::CElemSeg& es_v = this->GetElemSeg(id_ea,CORNER,true, world);
		const Fem::Field::CElemAry::CElemSeg& es_c = this->GetElemSeg(id_ea,CORNER,false,world);
		const unsigned int nnoes = 3;		
		assert( es_v.Length() == nnoes );
		assert( es_c.Length() == nnoes );
		assert( es_v.Size() == es_c.Size() );
		const unsigned int nelem = es_v.Size();
		for(unsigned int ielem=0;ielem<nelem;ielem++)
		{
			unsigned int noes_c[nnoes];
			es_c.GetNodes(ielem,noes_c);
			double ec[nnoes][2];
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				const unsigned int ino = noes_c[inoes];
				ns_c.GetValue(ino,ec[inoes]);
			}
			////////////////
			const double at = TriArea2D(ec[0],ec[1],ec[2]);
            const double a0 = TriArea2D(co,ec[1],ec[2]);    if( a0 < -at*1.0e-3 ) continue;
            const double a1 = TriArea2D(co,ec[2],ec[0]);    if( a1 < -at*1.0e-3 ) continue;
            const double a2 = TriArea2D(co,ec[0],ec[1]);    if( a2 < -at*1.0e-3 ) continue;
			////////////////
			unsigned int noes_v[nnoes];
			es_v.GetNodes(ielem,noes_v);
			double ev[nnoes][2];
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				const unsigned int ino = noes_v[inoes];
				ns_v.GetValue(ino,ev[inoes]);
			}
			////////////////
			velo[0] = (a0*ev[0][0] + a1*ev[1][0] + a2*ev[2][0])/at;
			velo[1] = (a0*ev[0][1] + a1*ev[1][1] + a2*ev[2][1])/at;
			id_ea_stat = id_ea;
			ielem_stat = ielem; 
			r1 = a1/at;
			r2 = a2/at;
			return true;
		}
	}
	id_ea_stat = 0;
	ielem_stat = 0;
	r1 = 0;
	r2 = 0;
	return false;
}

// MicroAVS inpファイルへの書き出し
bool CField::ExportFile_Inp(const std::string& file_name, const CFieldWorld& world)
{
	if( m_aElemIntp.size() != 1 ){
		std::cout << "未実装" << std::endl;
		assert(0);
	}
	const unsigned int id_ea = m_aElemIntp[0].id_ea;
	if( this->GetInterpolationType(id_ea,world) != HEX1001 && 
		this->GetInterpolationType(id_ea,world) != TET1001  )
	{
		std::cout << "未実装" << std::endl;
		assert(0);
	}

	std::ofstream fout(file_name.c_str());
	fout << 1 << "\n";
	fout << "data\n";
	fout << "step1\n";
	{
		unsigned int id_na_c_co = m_na_c.id_na_co;
		const CNodeAry& na = world.GetNA(id_na_c_co);
		const unsigned int nnode = na.Size();
		unsigned int id_ea = m_aElemIntp[0].id_ea;
		const CElemAry& ea = world.GetEA(id_ea);
		const unsigned int nelem = ea.Size();
		fout << nnode << " " << nelem << "\n";
	}
	{
		unsigned int id_na_c_co = m_na_c.id_na_co;
		unsigned int id_ns_c_co = m_na_c.id_ns_co;
		const CNodeAry& na = world.GetNA(id_na_c_co);
		const CNodeAry::CNodeSeg& ns_c_co = na.GetSeg(id_ns_c_co);
		double coord[3];
		for(unsigned int inode=0;inode<na.Size();inode++){
			fout << inode+1 << " ";
			ns_c_co.GetValue(inode,coord);
			for(unsigned int idim=0;idim<3;idim++){
				fout << coord[idim] << " ";
			}
			fout << "\n";
		}
	}

	const CElemAry& ea = world.GetEA(id_ea);
	std::string str_elem_type;
	{
		if( ea.ElemType() == TET ){
			str_elem_type = "tet";
		}
		else if( ea.ElemType() == HEX ){
			str_elem_type = "hex";
		}
		else{ assert(0); }
	}
	{
		const CElemAry::CElemSeg& es_c_co = this->GetElemSeg(id_ea,CORNER,false,world);
		unsigned int id_na_c_co = m_na_c.id_na_co;
		assert( id_na_c_co == es_c_co.GetIdNA() );
		unsigned int noes_c[16];
		unsigned int nnoes_c = es_c_co.Length();
		for(unsigned int ielem=0;ielem<ea.Size();ielem++){
			es_c_co.GetNodes(ielem,noes_c);
			fout << ielem+1 << " 0 " << str_elem_type << " ";
			for(unsigned int inoes=0;inoes<nnoes_c;inoes++){
				fout << noes_c[inoes]+1 << " ";
			}
			fout << "\n";
		}
	}
	fout << "0 3\n";
	fout << "3 1 1 1\n";
	fout << "hoge_x,\n";
	fout << "hoge_y,\n";
	fout << "hoge_z,\n";
	{
		const CElemAry::CElemSeg& es_b_va = this->GetElemSeg(id_ea,BUBBLE,true,world);
		unsigned int id_na_b_va = m_na_b.id_na_va;
		assert( id_na_b_va == es_b_va.GetIdNA() );
		const CNodeAry& na_b_va = world.GetNA(id_na_b_va);
		unsigned int id_ns_b_va = m_na_b.id_ns_va;
		const CNodeAry::CNodeSeg& ns_b_va = na_b_va.GetSeg(id_ns_b_va);
		unsigned int inoes_b;
		double value[3];
		const unsigned int nnoes = es_b_va.Length();
		assert( nnoes == 1 );
		for(unsigned int ielem=0;ielem<ea.Size();ielem++){
			es_b_va.GetNodes(ielem,&inoes_b);
			ns_b_va.GetValue(inoes_b,value);
			fout << ielem+1 << " " << value[0] << " " << value[1] << " " << value[2] << "\n";
		}
	}
	
	return true;
}
