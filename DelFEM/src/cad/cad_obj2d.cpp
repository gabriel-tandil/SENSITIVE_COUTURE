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
// CadObj2D.cpp : ２次元ＣＡＤモデルクラス(CCadObj2D)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#define for if(0);else for

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <cassert>	
#include <math.h>	
#include <cstring>	// strlen

#include "delfem/cad_obj2d.h"
#include "delfem/cad/cad_elem2d.h"

using namespace Cad;
using namespace Com;

CCadObj2D::CCadObj2D()
{
  this->min_clearance = 1.0e-3;
//  this->min_clearance = 0.1;  
}

CCadObj2D::CCadObj2D(const CCadObj2D& rhs)
{
  this->Clear();
  this->m_LoopSet = rhs.m_LoopSet;
  this->m_EdgeSet = rhs.m_EdgeSet;
  this->m_VertexSet = rhs.m_VertexSet;
  this->m_BRep = rhs.m_BRep;  
  this->min_clearance = rhs.min_clearance;
}

CCadObj2D::~CCadObj2D()
{
	this->Clear();
}

void CCadObj2D::Clear()
{
	this->m_LoopSet.Clear();
	this->m_EdgeSet.Clear();
	this->m_VertexSet.Clear();
	this->m_BRep.Clear();
}

////////////////////////////////////////////////////////////////
 
bool CCadObj2D::IsElemID(Cad::CAD_ELEM_TYPE itype,unsigned int id) const
{
  if(      itype == Cad::NOT_SET ){ return false; }
  else if( itype == Cad::VERTEX  ){ return m_VertexSet.IsObjID(id); }
  else if( itype == Cad::EDGE    ){ return m_EdgeSet.IsObjID(id);   }
  else if( itype == Cad::LOOP    ){ return m_LoopSet.IsObjID(id);   }  
	else{ assert(0); }
	return false;
}

const std::vector<unsigned int> CCadObj2D::GetAryElemID(Cad::CAD_ELEM_TYPE itype) const
{
	if(      itype == Cad::VERTEX ){ return m_VertexSet.GetAry_ObjID(); }
	else if( itype == Cad::EDGE   ){ return m_EdgeSet.GetAry_ObjID();   }
	else if( itype == Cad::LOOP   ){ return m_LoopSet.GetAry_ObjID();   }
	assert(0);
	std::vector<unsigned int> null_vec;
	return null_vec;
}

CVector2D CCadObj2D::GetVertexCoord(unsigned int id_v) const
{
	assert( m_VertexSet.IsObjID(id_v) );
	const CVertex2D& v = m_VertexSet.GetObj(id_v);
	return v.point;
}

const CEdge2D& CCadObj2D::GetEdge(unsigned int id_e) const
{
	assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	assert( this->m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	unsigned int id_vs, id_ve;
	m_BRep.GetIdVertex_Edge(id_e, id_vs, id_ve);
	e.id_v_s = id_vs;
	e.id_v_e = id_ve;
	assert( m_BRep.IsElemID(Cad::VERTEX,id_vs) );
	assert( m_BRep.IsElemID(Cad::VERTEX,id_ve) );
	assert( m_VertexSet.IsObjID(id_vs) );
	assert( m_VertexSet.IsObjID(id_ve) );
	e.po_s = this->GetVertexCoord(id_vs);
	e.po_e = this->GetVertexCoord(id_ve);
  e.bb_.isnt_empty = false;
	return e;
}

CEdge2D& CCadObj2D::GetEdgeRef(unsigned int id_e)
{
	assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	assert( this->m_EdgeSet.IsObjID(id_e) );
	CEdge2D& e = m_EdgeSet.GetObj(id_e);
	unsigned int id_vs, id_ve;
	m_BRep.GetIdVertex_Edge(id_e, id_vs, id_ve);
	e.id_v_s = id_vs;
	e.id_v_e = id_ve;
	assert( m_BRep.IsElemID(Cad::VERTEX,id_vs) );
	assert( m_BRep.IsElemID(Cad::VERTEX,id_ve) );
	assert( m_VertexSet.IsObjID(id_vs) );
	assert( m_VertexSet.IsObjID(id_ve) );
	e.po_s = this->GetVertexCoord(id_vs);
	e.po_e = this->GetVertexCoord(id_ve);
  e.bb_.isnt_empty = false;  
	return e;
}

bool CCadObj2D::GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const
{
  assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	return m_BRep.GetIdVertex_Edge(id_e,id_v_s,id_v_e);
}

unsigned int CCadObj2D::GetIdVertex_Edge(unsigned int id_e, bool is_s) const
{
  assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	return m_BRep.GetIdVertex_Edge(id_e,is_s);
}



bool CCadObj2D::GetIdLoop_Edge(unsigned int& id_l_l, unsigned int& id_l_r, unsigned int id_e) const
{
	return this->m_BRep.GetIdLoop_Edge(id_e,id_l_l,id_l_r);
}

int CCadObj2D::GetEdgeCurveType(const unsigned int& id_e) const
{
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	return e.itype;
}

bool CCadObj2D::GetCurve_Arc(unsigned int id_e, bool& is_left_side, double& dist) const
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	is_left_side = e.is_left_side;
	dist = e.dist;
	return true;
}

bool CCadObj2D::GetCurve_Polyline(unsigned int id_e, std::vector<double>& aRelCoMesh) const
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	const CEdge2D& e = m_EdgeSet.GetObj(id_e);
	aRelCoMesh = e.aRelCoMesh;
	return true;
}

bool CCadObj2D::GetCurveAsPolyline(
        unsigned int id_e, std::vector<Com::CVector2D>& aCo, double elen) const
{
  if( !m_EdgeSet.IsObjID(id_e) ){
    //		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
    return false;
  }
  const CEdge2D& e = this->GetEdge(id_e);
  double len = e.GetCurveLength();
  if( elen > 0 ){
    const unsigned int ndiv = (unsigned int)(len/elen)+1;
    return e.GetCurve_Mesh(aCo, ndiv);
  }
  return e.GetCurve_Mesh(aCo,-1);  
}

//! get division points of edge(ID:id_e) into ndiv
//! the begining and end of devision are ps,pe.
bool CCadObj2D::GetCurveAsPolyline(unsigned int id_e, std::vector<Com::CVector2D>& aCo,
        unsigned int ndiv, const Com::CVector2D& ps, const Com::CVector2D& pe) const
{
  if( !m_EdgeSet.IsObjID(id_e) ){
    //		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
    return false;
  }
  const CEdge2D& e = m_EdgeSet.GetObj(id_e);
  e.po_s = ps;    e.po_e = pe;
  return e.GetCurve_Mesh(aCo,ndiv);  
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

double CCadObj2D::GetArea_ItrLoop(CBRepSurface::CItrLoop& itrl) const
{
	double area = 0.0;
	for(itrl.Begin();!itrl.IsEnd();itrl++){
		unsigned int id_e;   bool is_same_dir;
    itrl.GetIdEdge(id_e,is_same_dir);
    if( !this->IsElemID(Cad::EDGE,id_e) ) return 0;
    assert( this->IsElemID(Cad::EDGE,id_e) );
		const CEdge2D& e = this->GetEdge(id_e);
		assert( ((is_same_dir) ? e.id_v_s : e.id_v_e) == itrl.GetIdVertex() );
		assert( ((is_same_dir) ? e.id_v_e : e.id_v_s) == itrl.GetIdVertex_Ahead() );
		// add area here
		const double earea = TriArea(e.po_s, e.po_e, CVector2D(0.0,0.0)) + e.AreaEdge();
		if( is_same_dir ){ area += earea; }
		else{              area -= earea; }
	}
	return area;
}

// return the area of loop
double CCadObj2D::GetArea_Loop(unsigned int id_l) const
{
	assert( m_LoopSet.IsObjID(id_l) );
	double area = 0.0;
    for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
		area += this->GetArea_ItrLoop(itrl);
	}
	return area;
}

// return the color of the loop (ID:id_l)
bool CCadObj2D::GetColor_Loop(unsigned int id_l, double color[3] ) const
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    const CLoop2D& l = m_LoopSet.GetObj(id_l);
    color[0] = l.m_color[0];
    color[1] = l.m_color[1];
    color[2] = l.m_color[2];
    return true;
}

int CCadObj2D::GetLayer(Cad::CAD_ELEM_TYPE type, unsigned int id) const
{
	if(      type == Cad::LOOP ){
		if( !m_LoopSet.IsObjID(id) ) return 0;
		const CLoop2D& l = m_LoopSet.GetObj(id);
		return l.ilayer;
	}
	else if( type == Cad::EDGE ){
		unsigned int id_l_l, id_l_r;
		this->GetIdLoop_Edge(id_l_l,id_l_r,id);
		const bool bl = this->IsElemID(Cad::LOOP,id_l_l);
		const bool br = this->IsElemID(Cad::LOOP,id_l_r);
		if( !bl && !br ){ return 0; }
		if(  bl && !br ){ return this->GetLayer(Cad::LOOP,id_l_l); }
		if( !bl &&  br ){ return this->GetLayer(Cad::LOOP,id_l_r); }
		const int ilayer_l = this->GetLayer(Cad::LOOP,id_l_l);
		const int ilayer_r = this->GetLayer(Cad::LOOP,id_l_r);
		return ( ilayer_l > ilayer_r )	? ilayer_l : ilayer_r;
	}
	else if( type == Cad::VERTEX ){		
		int ilayer = 0;
		bool iflg = true;
		for(CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id);!itrv.IsEnd();itrv++){
			unsigned int id_l0 = itrv.GetIdLoop();
			if( !this->IsElemID(Cad::LOOP,id_l0) ) continue;
			const int ilayer0 = this->GetLayer(Cad::LOOP,id_l0);
      if( iflg == true ){ ilayer = ilayer0; iflg = false; }
			else{ ilayer =  ( ilayer0 > ilayer )	? ilayer0 : ilayer; }
		}
		return ilayer;
	}
	return 0;
}

void CCadObj2D::GetLayerMinMax(int& ilayer_min, int& ilayer_max) const
{
	const std::vector<unsigned int>& aIdL = this->GetAryElemID(Cad::LOOP);
  if( aIdL.size() == 0 ){
    ilayer_min = 0;
    ilayer_max = 0;    
    return;
  }
	{
		assert( aIdL.size() > 0 );
		unsigned int id_l0 = aIdL[0];
		ilayer_min = this->GetLayer(Cad::LOOP,id_l0);
		ilayer_max = ilayer_min;
	}
	for(unsigned int il=0;il<aIdL.size();il++){
		unsigned int id_l = aIdL[il];
		int ilayer = this->GetLayer(Cad::LOOP,id_l);
		ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
		ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
	}
}

// ID:id_lのループの色を設定する
bool CCadObj2D::SetColor_Loop(unsigned int id_l, const double color[3] )
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    CLoop2D& l = m_LoopSet.GetObj(id_l);
    l.m_color[0] = color[0];
    l.m_color[1] = color[1];
    l.m_color[2] = color[2];
    return true;
}

bool CCadObj2D::ShiftLayer_Loop(unsigned int id_l, bool is_up)
{
    if( !m_LoopSet.IsObjID(id_l) ) return false;
    CLoop2D& l = m_LoopSet.GetObj(id_l);
	if( is_up ){
		l.ilayer++;
	}
	else{
		l.ilayer--;
	}
	return true;
}

bool CCadObj2D::CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const CVector2D& point) const
{
//	std::cout << "CCadObj2D::CheckIsPointInside_ItrLoop" << std::endl;
	for(unsigned int i=1;i<29;i++){	// 29 is handy prim number
		unsigned int cross_counter = 0;
		bool iflg = true;
		CVector2D dir(sin(6.28*i/29.0),cos(6.28*i/29.0));
		for(itrl.Begin();!itrl.IsEnd();itrl++){
			unsigned int id_e;   bool is_same_dir;
			itrl.GetIdEdge(id_e,is_same_dir);
			if( id_e == 0 ) return false;
			assert( m_EdgeSet.IsObjID(id_e) );
			const CEdge2D& e = this->GetEdge(id_e);
			const int ires = e.NumIntersect_AgainstHalfLine( point,dir );
			// -1 is vague so let's try again!
			if( ires == -1 ){ iflg = false; break; }
			cross_counter += ires;
		}
		if( iflg == true ){ 
			if( cross_counter % 2 == 0 ) return false;
			return true;
		}
	}
	assert(0);	// I hope process don't come here !!!!
	return false;
}

// ret:0 all the points in id_ul1 are inside id_ul2
// ret:1 there are both inside/outside id_ul1 points in id_ul1 or ambiguous
// ret:2 all the points in id_ul1 are outside id_ul2
unsigned int CCadObj2D::CheckInOut_ItrLoopPoint_ItrLoop
(CBRepSurface::CItrLoop& itrl1, CBRepSurface::CItrLoop& itrl2) const
{
//	std::cout << "CheckInOut_ItrLoopPoint_ItrLoop " << std::endl;
	unsigned int count_out=0, count_in=0;
	for(itrl1.Begin();!itrl1.IsEnd();itrl1++){
		const unsigned int id_v = itrl1.GetIdVertex();
		const CVertex2D& v = this->m_VertexSet.GetObj(id_v);
    const double dist = this->DistPointItrLoop(itrl2,v.point);
    if( fabs(dist) < min_clearance ){ return 1; } // ambiguous
    if( this->CheckIsPointInside_ItrLoop(itrl2, v.point) ){    
			if( count_out!=0 ) return 1;
			count_in++;
		}
		else{ // outside
			if( count_in!=0 ) return 1;
			count_out++;
		}
	}
	if( count_in == 0 ){
		assert( count_out != 0 );
		return 2;
	}	
	assert( count_out == 0 );
	assert( count_in != 0 );
	return 0;
} 

bool CCadObj2D::CheckIsPointInsideLoop(unsigned int id_l1, const CVector2D& point) const
{
	assert( m_LoopSet.IsObjID(id_l1) );
	for(CBRepSurface::CItrLoop itrl = m_BRep.GetItrLoop(id_l1);!itrl.IsEndChild();itrl.ShiftChildLoop()){    
    if( itrl.IsParent() ){
      if( !CheckIsPointInside_ItrLoop(itrl,point) ) return false;   // should be inside parent use loop
    }
    else{
      if(  CheckIsPointInside_ItrLoop(itrl,point) ) return false;  // should be outside child use loop
    }
	} 
	return true;
}


double CCadObj2D::SignedDistPointLoop
(unsigned int id_l1, const Com::CVector2D& point, 
 unsigned int id_v_ignore) const
{
  double min_sd = 0;
	assert( m_LoopSet.IsObjID(id_l1) );
	for(CBRepSurface::CItrLoop itrl = m_BRep.GetItrLoop(id_l1);!itrl.IsEndChild();itrl.ShiftChildLoop()){    
    if( itrl.IsParent() ){
      min_sd = +DistPointItrLoop(itrl,point);
      assert( min_sd >= 0 );
      if( !CheckIsPointInside_ItrLoop(itrl,point) ){ min_sd = -min_sd; }
    }
    else{
      if( itrl.GetIdVertex() == itrl.GetIdVertex_Ahead() ){
        unsigned int id_v = itrl.GetIdVertex();
        if( id_v == id_v_ignore ) continue;
      }
      double sd0 = DistPointItrLoop(itrl,point);  
      if( sd0 < 0 ) continue;      
      if( CheckIsPointInside_ItrLoop(itrl,point) ){ sd0 = -sd0; }
      if( fabs(sd0) < fabs(min_sd) ){ min_sd = sd0; } 
    }
	} 
	return min_sd;
}


double CCadObj2D::DistPointItrLoop(CBRepSurface::CItrLoop& itrl, const CVector2D& point) const
{
  double min_dist = -1;
  for(itrl.Begin();!itrl.IsEnd();itrl++){
    unsigned int id_e;   bool is_same_dir;
    itrl.GetIdEdge(id_e,is_same_dir);
    if( id_e == 0 ){
      unsigned int id_v0 = itrl.GetIdVertex();
      assert( this->IsElemID(Cad::VERTEX, id_v0) );
      const CVector2D& p1 = this->GetVertexCoord(id_v0);
      return Distance(point,p1);
    }
    assert( m_EdgeSet.IsObjID(id_e) );
    const CEdge2D& e = this->GetEdge(id_e);
    const Com::CVector2D& v = e.GetNearestPoint(point);
    const double d0 = Distance(v,point);
    if( min_dist < 0 || d0 < min_dist ){ min_dist = d0; }
	}
  return min_dist;
}


int CCadObj2D::CheckLoop(unsigned int id_l) const
{
//  std::cout << "Check Loop " << id_l << std::endl;
  {	// check self interaction of the loop
		if( this->CheckIntersection_Loop(id_l) ){ return 1; }
	}
	{	// 親ループの面積が正で子ループの面積が負であることを調べる
		for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			if( itrl.IsParent() ){
				if( itrl.GetType() != 2 ) return 2;
				if( this->GetArea_ItrLoop(itrl) < 0 ) return 2;
			}
			else if( itrl.GetType() == 2 ){
				if( this->GetArea_ItrLoop(itrl) > 0 ) return 2;
			}
		}
	}
	{	// 子ループ内の点が親ループの中に入っているかどうか調べる
    CBRepSurface::CItrLoop itrl_p = m_BRep.GetItrLoop(id_l);
		for(CBRepSurface::CItrLoop itrl_c=m_BRep.GetItrLoop(id_l);!itrl_c.IsEndChild();itrl_c.ShiftChildLoop()){
			if( itrl_c.IsParent() ) continue;
      if( this->CheckInOut_ItrLoopPoint_ItrLoop(itrl_c,itrl_p)!=0 ) return 3;
		}
	}	
	{   // 子ループ同士が互いにお互いの頂点を含まないかを調べる
		for(CBRepSurface::CItrLoop itrl1=m_BRep.GetItrLoop(id_l);!itrl1.IsEndChild();itrl1.ShiftChildLoop()){		
			if( itrl1.IsParent() ){ continue; }
			for(CBRepSurface::CItrLoop itrl2=m_BRep.GetItrLoop(id_l);!itrl2.IsEndChild();itrl2.ShiftChildLoop()){
				if( itrl2.IsParent() ){ continue; }
				if( itrl1.IsSameUseLoop(itrl2) ){ continue; }
				if( this->CheckInOut_ItrLoopPoint_ItrLoop(itrl1,itrl2)!=2 ) return 4;
			}
		}
	}
	{	// ループのの中の各角についてその角の中に他の点が入っていないかを確認する
		Com::CVector2D vec_zero(0,0);
		for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(itrl.Begin();!itrl.IsEnd();itrl++){	// ループの中の点をめぐる
				unsigned int id_vm, id_vf;
				Com::CVector2D dir;
				{
					id_vm = itrl.GetIdVertex();
					id_vf = itrl.GetIdVertex_Ahead();
					unsigned int id_e0;   bool is_same_dir0;
					itrl.GetIdEdge(id_e0,is_same_dir0);
					if( !this->m_EdgeSet.IsObjID(id_e0) ) continue;
					const CEdge2D& e0 = this->GetEdge(id_e0);
					dir = e0.GetTangentEdge(is_same_dir0);
					assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_vm );
					assert( ((is_same_dir0) ? e0.id_v_e : e0.id_v_s) == id_vf );
				}
				for(CBRepSurface::CItrVertex itrv=m_BRep.GetItrVertex(id_vm);!itrv.IsEnd();itrv++){	// 点周りの辺をめぐる
					unsigned int id_e0;   bool is_same_dir0;
					itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
					if( !this->m_EdgeSet.IsObjID(id_e0) ){ continue; }
					const Cad::CEdge2D& e0 = this->GetEdge(id_e0);
					assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_vm );
					if( ((is_same_dir0) ? e0.id_v_e : e0.id_v_s) == id_vf ){ continue; }
					const CVector2D& tan0 = e0.GetTangentEdge(is_same_dir0);
					////////////////
					unsigned int id_e1;   bool is_same_dir1;
					itrv.GetIdEdge_Ahead(id_e1,is_same_dir1);
					if( !this->m_EdgeSet.IsObjID(id_e1) ){ continue; }
					const Cad::CEdge2D& e1 = this->GetEdge(id_e1);
					assert( ((is_same_dir1) ? e1.id_v_s : e1.id_v_e) == id_vm );
					if( ((is_same_dir1) ? e1.id_v_e : e1.id_v_s) == id_vf ){ continue; }
					const CVector2D& tan1 = e1.GetTangentEdge(is_same_dir1);
					////////////////
					const double area0 = TriArea(tan1,vec_zero,tan0);
					const double area1 = TriArea(tan1,vec_zero,dir);	
					const double area2 = TriArea(dir, vec_zero,tan0);
					if( (area0 > 0.0 &&  area1 > 0.0 && area2 > 0.0) || 
   						(area0 < 0.0 && (area1 > 0.0 || area2 > 0.0) ) ){ return 5; }
				}
			}
		}
	}
  return 0;
}

int CCadObj2D::AssertValid() const 
{
	{	// Check each loops
		std::vector<unsigned int> id_l_ary = m_LoopSet.GetAry_ObjID();
		for(unsigned iid_l=0;iid_l<id_l_ary.size();iid_l++){
			const unsigned int id_l = id_l_ary[iid_l];
      int res = CheckLoop(id_l);
      if( res!=0 ){
        if(      res == 1 ){ std::cout << "Intersectoin in the loop" << std::endl; }
				else if( res == 2 ){ std::cout << "Check area parent plus, childe minus" << std::endl; }
				else if( res == 3 ){ std::cout << "Check whether childe loop included in parent loop" << std::endl; }
				else if( res == 4 ){ std::cout << "Check childe loop excluded from other child loop" << std::endl; }
				else if( res == 5 ){ std::cout << "Check positive angle around vertex on the loop" << std::endl; }
				return res;
      }
		}
	}
  if( !m_BRep.AssertValid() ){ return 6; }
  {
    std::vector<unsigned int> aIdL = m_BRep.GetAryElemID(Cad::LOOP);
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){ 
      if( !this->m_LoopSet.IsObjID(aIdL[iidl]  ) ){
        std::cout << aIdL[iidl] << std::endl;
        return 7; 
      }
    }
  }
  {
    std::vector<unsigned int> aIdE = m_BRep.GetAryElemID(Cad::EDGE);
    for(unsigned int iide=0;iide<aIdE.size();iide++){ if( !this->m_EdgeSet.IsObjID(aIdE[iide]  ) ) return 7; }
  }
  {
    std::vector<unsigned int> aIdV = m_BRep.GetAryElemID(Cad::VERTEX);
    for(unsigned int iidv=0;iidv<aIdV.size();iidv++){ if( !this->m_VertexSet.IsObjID(aIdV[iidv]) ) return 7; }
  }                                       
	return 0;
}


Com::CVector2D CCadObj2D::GetNearestPoint(unsigned int id_e, const Com::CVector2D& po_in) const
{
  if( !this->IsElemID(Cad::EDGE,id_e) ) return Com::CVector2D();
  const CEdge2D& e1 = this->GetEdge(id_e);
  return e1.GetNearestPoint(po_in);
}

bool CCadObj2D::GetPointOnCurve_OnCircle
(unsigned int id_e,
 const Com::CVector2D& v0, double len, bool is_front,
 bool& is_exceed, Com::CVector2D& out) const
{
  if( !this->IsElemID(Cad::EDGE,id_e) ) return false;
  const CEdge2D& e1 = this->GetEdge(id_e);
  return e1.GetPointOnCurve_OnCircle(v0,len,is_front,is_exceed,out);
}


// Private関数
// 返り値は半辺のIDでid_vを基点としたdir1の方向に伸びる半直線からid_vを中心に時計周りに最初に出会う半辺
// この半辺の属するループはこの半直線と重なっている。
// id_vが浮遊点の場合は浮遊点周りのループが帰る
CBRepSurface::CItrVertex CCadObj2D::FindCorner_HalfLine(unsigned int id_v, const CVector2D& dir1) const
{
	assert( m_VertexSet.IsObjID(id_v) );
	Com::CVector2D dir = dir1;
	dir.Normalize();
	const Com::CVector2D vec_zero(0,0);
	Cad::CBRepSurface::CItrVertex itrv = this->m_BRep.GetItrVertex(id_v);
	if( itrv.CountEdge() < 2 ){ return itrv; }
	for(;!itrv.IsEnd();itrv++){
		unsigned int id_e0;   bool is_same_dir0;
		itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
		assert( this->m_EdgeSet.IsObjID(id_e0) );	// id_vが浮遊点の場合は省かれている
		const Cad::CEdge2D& e0 = this->GetEdge(id_e0);
		assert( ((is_same_dir0) ? e0.id_v_s : e0.id_v_e) == id_v );
		const CVector2D& tan0 = e0.GetTangentEdge(is_same_dir0);
		////////////////
		unsigned int id_e1;   bool is_same_dir1;
		itrv.GetIdEdge_Ahead(id_e1,is_same_dir1);
		assert( this->m_EdgeSet.IsObjID(id_e1) );	// id_vが浮遊点の場合は省かれている
		const Cad::CEdge2D& e1 = this->GetEdge(id_e1);
		assert( ((is_same_dir1) ? e1.id_v_s : e1.id_v_e) == id_v );
		const CVector2D& tan1 = e1.GetTangentEdge(is_same_dir1);
		////////////////
		assert( id_e0 != id_e1 );	// id_vが端点の場合は省かれているはず
		const double area0 = TriArea(tan1,vec_zero,tan0);
		const double area1 = TriArea(tan1,vec_zero,dir);	
		const double area2 = TriArea(dir, vec_zero,tan0);
		if( area0 > 0.0 ){	if( area1 > 0.0 && area2 > 0.0 ){ return itrv; } }
		else{				if( area1 > 0.0 || area2 > 0.0 ){ return itrv; } }
	}
	return itrv;
}

CBRepSurface::CResConnectVertex CCadObj2D::ConnectVertex(CEdge2D edge)
{
	const unsigned int id_v1 = edge.id_v_s;
	const unsigned int id_v2 = edge.id_v_e;
  CBRepSurface::CResConnectVertex res;
  res.id_v1 = id_v1; res.id_v2 = id_v2;
  ////  
	if( !m_VertexSet.IsObjID(id_v1) ){ return res; }	
	if( !m_VertexSet.IsObjID(id_v2) ){ return res; }
	if( id_v1 == id_v2 ){ return res; }

	if( edge.itype == 0 ){	// check if there is line edge with end points (v1,v2)
		const std::vector<unsigned int>& id_ary_e = m_EdgeSet.GetAry_ObjID();
		for(unsigned int iid_e=0;iid_e<id_ary_e.size();iid_e++){
			const unsigned int id_e = id_ary_e[iid_e]; 	assert( m_EdgeSet.IsObjID(id_e) );
			const CEdge2D& e = this->GetEdge(id_e);
			if( e.itype != 0 ) continue;
			const unsigned int id_v_e = e.id_v_e;
			const unsigned int id_v_s = e.id_v_s;
			if( (id_v_s-id_v1)*(id_v_s-id_v2) != 0 ) continue;
			if( (id_v_e-id_v1)*(id_v_e-id_v2) != 0 ) continue;
			return res;
		}
	}
	edge.po_s = m_VertexSet.GetObj(id_v1).point;
	edge.po_e = m_VertexSet.GetObj(id_v2).point;
	if( edge.IsCrossEdgeSelf() ){ return res; }	// if this edge have self interesection, it is impossible to create edge
  
	////////////////
	{
		const CBRepSurface::CItrVertex& itrv1 = this->FindCorner_HalfLine(id_v1, edge.GetTangentEdge(true)  );
		const CBRepSurface::CItrVertex& itrv2 = this->FindCorner_HalfLine(id_v2, edge.GetTangentEdge(false) );
		if( itrv1.GetIdLoop() != itrv2.GetIdLoop() ){ return res; }
		const unsigned int id_l = itrv1.GetIdLoop();
    /////////////
    if( CheckIntersection_EdgeAgainstLoop(edge, id_l) ){ return res; }
    /////////////
    bool is_left_ladd = false;    
		if( itrv1.IsSameUseLoop(itrv2) && ( !itrv1.IsParent() || id_l == 0 ) ){
      // array of edge that begin id_v2 and end id_v1 that is part of loop generated if v1 and v2 are connected
			const std::vector< std::pair<unsigned int,bool> >& aIdEDir = m_BRep.GetItrLoop_ConnectVertex(itrv1,itrv2);
			assert( !aIdEDir.empty() );
			double area = TriArea(edge.po_s, edge.po_e, CVector2D(0,0)) + edge.AreaEdge();
			for(unsigned int ie=0;ie<aIdEDir.size();ie++){
				unsigned int id_e = aIdEDir[ie].first;
				assert( this->IsElemID(Cad::EDGE,id_e) );
				const CEdge2D& e = this->GetEdge(id_e);
				const double earea = e.AreaEdge() + TriArea(e.po_s, e.po_e, CVector2D(0,0));
				if( aIdEDir[ie].second ){ area += earea; }
				else{                     area -= earea; }
			}
			is_left_ladd = ( area > 0 );
		}
    ////////////
    // Leave input check section    
		res = m_BRep.ConnectVertex(itrv1,itrv2, is_left_ladd);
	}
	{ // register the edge
		int tmp_id = this->m_EdgeSet.AddObj( std::make_pair(res.id_e_add,edge) );
    assert( tmp_id == (int)res.id_e_add );
	}

	{
		unsigned int id_l_l, id_l_r;
		m_BRep.GetIdLoop_Edge(res.id_e_add, id_l_l,id_l_r);
		if( id_l_l == id_l_r ){ // the loop was not divided
			assert( this->AssertValid()==0 );
			return res;
		}	
		res.id_l_add = (res.is_left_l_add) ? id_l_l : id_l_r; 
		assert( res.id_l_add != res.id_l || res.id_l == 0 );
		assert( ((res.is_left_l_add) ? id_l_l : id_l_r) == res.id_l_add );
	}
  
  if( !m_BRep.IsElemID(Cad::LOOP,res.id_l) ){	// 外の領域に切った場合に，新しいループが自己干渉している場合はこのループは無し
    if( this->CheckIntersection_Loop(res.id_l_add) ){
			m_BRep.MakeHole_fromLoop(res.id_l_add);
			assert( this->AssertValid()==0 );
			return res;
		}
	}  
     
	if( m_BRep.IsElemID(Cad::LOOP,res.id_l)  && m_BRep.IsElemID(Cad::LOOP,res.id_l_add)  ){
    // Here (id_l != id_l_add) and we sort child loop of id_l into id_l or ld_l_add    
		for(;;){
			bool iflg = true;
			CBRepSurface::CItrLoop itrl_add_inner = m_BRep.GetItrLoop_SideEdge(res.id_e_add, res.is_left_l_add);
			CBRepSurface::CItrLoop itrl_add_outer = m_BRep.GetItrLoop_SideEdge(res.id_e_add,!res.is_left_l_add);
			for(CBRepSurface::CItrLoop itrl_c=m_BRep.GetItrLoop(res.id_l);!itrl_c.IsEndChild();itrl_c.ShiftChildLoop()){
				if( itrl_c.IsParent() ) continue;
				if( itrl_c.IsSameUseLoop(itrl_add_outer) ) continue;
				unsigned int ires = this->CheckInOut_ItrLoopPoint_ItrLoop(itrl_c,itrl_add_inner);
				assert( ires == 0 || ires == 2 );
				if( ires == 0 ){ 
					m_BRep.SwapItrLoop(itrl_c,res.id_l_add); 
					iflg = false;
					break;
				}
			}
			if( iflg ) break;
		}
	}
  
  if( m_LoopSet.IsObjID(res.id_l) ){  // created new loop (id_l_add) by splitting existing loop (id_l)
		CLoop2D loop_add = this->m_LoopSet.GetObj(res.id_l);
		this->m_LoopSet.AddObj( std::make_pair(res.id_l_add,loop_add) );
	}
	else{ // created entirely new loop
		this->m_LoopSet.AddObj( std::make_pair(res.id_l_add,CLoop2D()) );
	}

	assert( this->AssertValid()==0 );
	return res;
}

bool CCadObj2D::RemoveElement(Cad::CAD_ELEM_TYPE itype, unsigned int id)
{
	if( !this->IsElemID(itype,id) ) return false;
	if(      itype == Cad::EDGE   ){
		CBRepSurface::CItrLoop itrl_l = m_BRep.GetItrLoop_SideEdge(id,true );
		CBRepSurface::CItrLoop itrl_r = m_BRep.GetItrLoop_SideEdge(id,false);
		unsigned int id_l_l = itrl_l.GetIdLoop();
		unsigned int id_l_r = itrl_r.GetIdLoop();
		unsigned int id_v1, id_v2;
		m_BRep.GetIdVertex_Edge(id, id_v1,id_v2);
		CBRepSurface::CItrVertex itrv1 = m_BRep.GetItrVertex(id_v1);
		CBRepSurface::CItrVertex itrv2 = m_BRep.GetItrVertex(id_v2);
		bool is_del_cp = false;
		if( itrl_l.IsSameUseLoop(itrl_r) && itrl_l.IsParent() && itrl_l.GetIdLoop()!=0  
			&& itrv1.CountEdge() > 1 && itrv2.CountEdge() > 1 )	// 親ループと子ループを繋ぐ辺を削除した場合
		{
			const std::vector< std::pair<unsigned int,bool> >& aIdEDir = m_BRep.GetItrLoop_RemoveEdge(id);
			assert( !aIdEDir.empty() );
			{	// 明らかに面積が無いかどうか調べる
				unsigned int ie=0;
				for(;ie<aIdEDir.size();ie++){
					unsigned int je = 0;
					for(;je<aIdEDir.size();je++){
						if( ie == je ) continue;
						if( aIdEDir[ie].first == aIdEDir[je].first ){ 
							assert( aIdEDir[ie].second != aIdEDir[je].second );
							break;
						}
					}
					if( je == aIdEDir.size() ){ break; }
				}
				is_del_cp = ( ie == aIdEDir.size() );
			}
			if( !is_del_cp ){	// 面積を実際に調べてみる
				double area = 0.0;
				for(unsigned int ie=0;ie<aIdEDir.size();ie++){
					unsigned int id_e = aIdEDir[ie].first;
					assert( this->IsElemID(Cad::EDGE,id_e) );
					const CEdge2D& e = this->GetEdge(id_e);
					const double earea = e.AreaEdge() + TriArea(e.po_s, e.po_e, CVector2D(0,0));
					if( aIdEDir[ie].second ){ area += earea; }
					else{                     area -= earea; }
				}
				if( area < 0 ){ is_del_cp = true; }
			}
		}
		if( !m_BRep.RemoveEdge(id,is_del_cp) ){ 
      std::cout << "Remove Edge B-Rep unsuccessfull : " << id << std::endl;
      return false; 
    }
		m_EdgeSet.DeleteObj(id);
		if( !m_BRep.IsElemID(Cad::LOOP,id_l_l) ){ m_LoopSet.DeleteObj(id_l_l); }
		if( !m_BRep.IsElemID(Cad::LOOP,id_l_r) ){ m_LoopSet.DeleteObj(id_l_r); }
		assert( this->AssertValid()==0 );
		return true;
	}
	else if( itype == Cad::VERTEX )
	{ 
		CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id);
		if( itrv.CountEdge() == 2 ){
			unsigned int id_e1,id_e2;   bool btmp1, btmp2;
			itrv.GetIdEdge_Ahead( id_e1,btmp1); 
			itrv.GetIdEdge_Behind(id_e2,btmp2);
			CEdge2D e_tmp = this->GetEdge(id_e1);
      {
        //                const unsigned int id_v1 = m_BRep.GetIdVertex_Edge(id_e1,!btmp1);
        const unsigned int id_v2 = m_BRep.GetIdVertex_Edge(id_e2,!btmp2);
        assert( m_BRep.GetIdVertex_Edge(id_e1,btmp1) == id );
        assert( m_BRep.GetIdVertex_Edge(id_e2,btmp2) == id );
        e_tmp.ConnectEdge(this->GetEdge(id_e2),!btmp1,btmp1!=btmp2);
        if( btmp1 ){ e_tmp.id_v_s = id_v2; e_tmp.po_s = this->GetVertexCoord(id_v2); }
        else{        e_tmp.id_v_e = id_v2; e_tmp.po_e = this->GetVertexCoord(id_v2); }
      }
      ////////////////
      {   // 辺がとe_tmpが交差していないかを全ての辺について調べる
        const unsigned int ipo0 = e_tmp.id_v_s;
        const unsigned int ipo1 = e_tmp.id_v_e;
        const Com::CBoundingBox2D& bb_i = e_tmp.GetBoundingBox(); // get bounding box of e_tmp
				const std::vector<unsigned int>& aIdE = m_BRep.GetAryElemID(Cad::EDGE);
				for(unsigned int ije=0;ije<aIdE.size();ije++){
					const unsigned int id_je = aIdE[ije];
          if( id_je == id_e2 || id_je == id_e1 ) continue;
					const CEdge2D& e_j = this->GetEdge(id_je);
					const unsigned int jpo0 = e_j.id_v_s;
					const unsigned int jpo1 = e_j.id_v_e;			
					if( (ipo0-jpo0)*(ipo0-jpo1)*(ipo1-jpo0)*(ipo1-jpo1) != 0 ){ // no-sharing point
            const Com::CBoundingBox2D& bb_j = e_j.GetBoundingBox(); // get bounding box of e_j
            if( !bb_i.IsIntersect(bb_j,min_clearance) ) continue;                    
            const double dist = e_tmp.Distance(e_j);
            if( dist > min_clearance ) continue;
            return true;
					}
          else if( ipo0 == jpo0 && ipo1 == jpo1 ){ if( e_tmp.IsCrossEdge_ShareBothPoints(e_j,true ) == 1 ){ return false; } }
          else if( ipo0 == jpo1 && ipo1 == jpo0 ){ if( e_tmp.IsCrossEdge_ShareBothPoints(e_j,false) == 1 ){ return false; } }
          else if( ipo0 == jpo0 ){ if( e_tmp.IsCrossEdge_ShareOnePoint(e_j, true, true)==1 ){ return false; } }
          else if( ipo0 == jpo1 ){ if( e_tmp.IsCrossEdge_ShareOnePoint(e_j, true,false)==1 ){ return false; } }
          else if( ipo1 == jpo0 ){ if( e_tmp.IsCrossEdge_ShareOnePoint(e_j,false, true)==1 ){ return false; } }
          else if( ipo1 == jpo1 ){ if( e_tmp.IsCrossEdge_ShareOnePoint(e_j,false,false)==1 ){ return false; } }
				}      
      }
			if( !m_BRep.RemoveVertex(id) ){ return false; }
      assert(  m_BRep.IsElemID(Cad::EDGE,id_e1) );
			assert( !m_BRep.IsElemID(Cad::EDGE,id_e2) );
			m_EdgeSet.DeleteObj(id_e2);
			this->GetEdgeRef(id_e1) = e_tmp;
			m_VertexSet.DeleteObj(id);
			assert( this->AssertValid()==0 );
			return true;
		}
		else if( itrv.CountEdge() == 0 ){
			if( !m_BRep.RemoveVertex(id) ){ return false; }
			m_VertexSet.DeleteObj(id);
			assert( this->AssertValid()==0 );
			return true;
		}
	}
	return false;
}


// this function returns loop id
// in case of failure retun 0
// (return variavle should be class contains more informations)
CCadObj2D::CResAddPolygon  Cad::CCadObj2D::AddPolygon(const std::vector<Com::CVector2D>& aPoint_input, unsigned int id_l )
{
  CResAddPolygon res;
	const unsigned int npoint = aPoint_input.size();
	if( npoint < 3 ) return res;

  std::vector<Com::CVector2D> aPoint = aPoint_input;
	{	// 交線判断
		const unsigned int n = aPoint_input.size();
		std::vector<CEdge2D> aEdge;
		for(unsigned int i=0;i<n-1;i++){
			CEdge2D e(i,i+1);
			e.po_s = aPoint_input[i  ];
			e.po_e = aPoint_input[i+1];
			aEdge.push_back( e );
		}
		{
			CEdge2D e(n-1,0);
			e.po_s = aPoint_input[n-1];
			e.po_e = aPoint_input[0];
			aEdge.push_back( e );
		}
		if( Cad::CheckEdgeIntersection(aEdge) != 0 ) return res;
	}
	// (注)時計回りになってる場合でも以下は動きます
//	std::vector<unsigned int> aIdV;
//	std::vector<unsigned int> aIdE;
	// add vertex
  res.aIdV.clear(); res.aIdE.clear();
	res.aIdV.reserve(npoint);
	for(unsigned int ipoint=0;ipoint<npoint;ipoint++){
		unsigned int id_v0 = this->AddVertex(Cad::LOOP, id_l, aPoint[ipoint] ).id_v_add;
		if( id_v0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		res.aIdV.push_back(id_v0);
	}
	// add edge
	res.aIdE.reserve(npoint);
	for(unsigned int iedge=0;iedge<npoint-1;iedge++){
		unsigned int id_e0 = this->ConnectVertex_Line( res.aIdV[iedge], res.aIdV[iedge+1] ).id_e_add;
		if( id_e0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		res.aIdE.push_back(id_e0);
	}
	{
		unsigned int id_e0 = this->ConnectVertex_Line( res.aIdV[npoint-1], res.aIdV[0] ).id_e_add;
		if( id_e0 == 0 ) goto FAIL_ADD_POLYGON_INSIDE_LOOP;
		res.aIdE.push_back(id_e0);
	}
	assert( this->AssertValid() == 0 );
	// 新しく出来たループのIDを取得
  
	{	// 辺の両側のループを調べる
		unsigned int id_e0 = res.aIdE[ npoint-1 ];
		unsigned int id_l0, id_l1;
		m_BRep.GetIdLoop_Edge(id_e0, id_l0,id_l1);
		res.id_l_add = ( id_l0 == id_l ) ? id_l1 : id_l0;
	}
  return res;
	////////////////////////////////
	// failure procedure
FAIL_ADD_POLYGON_INSIDE_LOOP :
	for(unsigned int iie=0;iie<res.aIdE.size();iie++){
		unsigned int id_e0 = res.aIdE[iie];
		this->RemoveElement(Cad::EDGE,id_e0);
	}
	for(unsigned int iiv=0;iiv<res.aIdV.size();iiv++){
		unsigned int id_v0 = res.aIdV[iiv];
		this->RemoveElement(Cad::VERTEX,id_v0);
	}
	assert( this->AssertValid()==0 );
	return CResAddPolygon();	
}

CCadObj2D::CResAddVertex CCadObj2D::AddVertex(Cad::CAD_ELEM_TYPE itype, unsigned int id, const  Com::CVector2D& vec)
{
//  std::cout << "CadObj2D::AddVertex" << std::endl;
  CResAddVertex res;
	if(      itype == Cad::NOT_SET || id == 0 )
	{
		unsigned int id_v_add = m_BRep.AddVertex_Loop(0);
		const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec)) );
		assert( tmp_id ==(int)id_v_add );
    res.id_v_add = id_v_add;
		return res;
	}
	else if( itype == Cad::LOOP )
	{
		unsigned int id_l = id;
		assert( m_LoopSet.IsObjID(id_l) );    
		if( !m_LoopSet.IsObjID(id_l) ) return res;    
    { // check this point is inside the loop with clearance
      const double dist = this->SignedDistPointLoop(id_l,vec);
      if( dist < this->min_clearance ){ return res; }
    }
		unsigned int id_v_add = m_BRep.AddVertex_Loop(id_l);
		const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec)) );
		assert( tmp_id == (int)id_v_add );
		assert( this->AssertValid()==0 );
    res.id_v_add = id_v_add;
		return res;
	}
	else if( itype == Cad::EDGE )
	{
		unsigned int id_e = id;
		if( !m_EdgeSet.IsObjID(id_e) ) return res;
		CEdge2D edge_old = this->GetEdge(id_e);
		// check if the projection on edge is possible
		Com::CVector2D vec_add = edge_old.GetNearestPoint(vec);
		if(   SquareLength(vec_add-edge_old.po_e) < 1.0e-20 
       || SquareLength(vec_add-edge_old.po_s) < 1.0e-20 ){ return res; }
		////////////////////////////////
		// Leave Input Check Section

		unsigned int id_v_add = m_BRep.AddVertex_Edge(id_e);
		unsigned int id_e_add;
		{
			CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id_v_add);
			bool is_same_dir0;
			unsigned int id_e_b,id_e_a;
			itrv.GetIdEdge_Behind(id_e_b,is_same_dir0);
			itrv.GetIdEdge_Ahead( id_e_a,is_same_dir0);
			id_e_add = ( id_e_b == id_e ) ? id_e_a : id_e_b;
		}		
		{	// add vertex
			const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex2D(vec_add)) );
			assert( tmp_id == (int)id_v_add );
		}
		{	// add edge
			const int tmp_id = m_EdgeSet.AddObj( std::make_pair(id_e_add,CEdge2D(id_v_add,edge_old.id_v_e)) );
			assert( tmp_id == (int)id_e_add );
		}
		{
			CEdge2D& edge_a = this->GetEdgeRef(id_e_add);
			edge_old.Split(edge_a, vec_add);
			CEdge2D& edge = this->GetEdgeRef(id_e);
			assert( edge.id_v_e == id_v_add );
			edge = edge_old;
		}
		assert( this->AssertValid() == 0 );
    res.id_v_add = id_v_add;
    res.id_e_add = id_e_add;
		return res;
	}
	return res;
}

bool CCadObj2D::SetCurve_Line(const unsigned int id_e)
{
//    std::cout << "CCadObj2D::SetCurve_Line" << id_e << std::endl;
	if( !m_EdgeSet.IsObjID(id_e) ){
		assert(0);
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
  CEdge2D& e = m_EdgeSet.GetObj(id_e);
  CEdge2D e_old = e;
	////////////////
	e.itype = 0;
	////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
    if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

bool CCadObj2D::SetCurve_Arc(const unsigned int id_e, bool is_left_side, double rdist)
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		assert(0);
		return false;
  }
	assert( m_EdgeSet.IsObjID(id_e) );
  CEdge2D& e = this->GetEdgeRef(id_e);
	CEdge2D e_old = e;
  ////////////////////////////////
  // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
	e.itype = 1;
	e.is_left_side = is_left_side;
	e.dist = rdist;
  // ここまで
  ////////////////////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
    if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

bool CCadObj2D::SetCurve_Polyline(const unsigned int id_e)
{
  if( !m_EdgeSet.IsObjID(id_e) ){
//    std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
    assert(0);
    return false;
  }
  assert( m_EdgeSet.IsObjID(id_e) );
  CEdge2D& e = this->GetEdgeRef(id_e);
  CEdge2D e_old = e;
  ////////////////////////////////
  // ここからを現在のCurveTypeによって決める,次の設定は直線の場合
  const CVector2D& pos = e_old.po_s;
  const CVector2D& poe = e_old.po_e;
  std::vector<Com::CVector2D> aCo;
  e_old.GetCurve_Mesh(aCo,20);
  const double sqlen = Com::SquareLength(poe-pos);
  const Com::CVector2D& eh = (poe-pos)*(1/sqlen);
  const Com::CVector2D ev(-eh.y,eh.x);
  e.itype = 2;
  e.aRelCoMesh.clear();
  for(unsigned int ico=0;ico<aCo.size();ico++){
    double x1 = Com::Dot(aCo[ico]-pos,eh);
    double y1 = Com::Dot(aCo[ico]-pos,ev);
    e.aRelCoMesh.push_back(x1);
    e.aRelCoMesh.push_back(y1);
  }
  // ここまで
  ////////////////////////////////
  const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
  for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
    unsigned int id_l = aID_Loop[iid_l];
    if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
  }
  return true;
FAILURE:
  e = e_old;
  return false;  
}

bool Cad::CCadObj2D::SetCurve_Polyline(unsigned int id_e, const std::vector<Com::CVector2D>& aCo)
{
	if( !m_EdgeSet.IsObjID(id_e) ){
//		std::cout << "SetCurve_Arc Failure(this Edge is not exist)" << std::endl;
		assert(0);
		return false;
	}
	assert( m_EdgeSet.IsObjID(id_e) );
	CEdge2D& e = this->GetEdgeRef(id_e);
	CEdge2D e_old = e;
	////////////////
	e.itype = 2;
	{	// 相対座標を作る
		const unsigned int n = aCo.size();
		e.aRelCoMesh.resize(0);
		e.aRelCoMesh.reserve(n*2);
		const Com::CVector2D& pos = e.po_s;
		const Com::CVector2D& poe = e.po_e;
		const double sqlen = Com::SquareLength(poe-pos);
		const Com::CVector2D& eh = (poe-pos)*(1/sqlen);
		const Com::CVector2D ev(-eh.y,eh.x);
		for(unsigned int i=0;i<n;i++){
			double x0 = Com::Dot(aCo[i]-pos,eh);
			double y0 = Com::Dot(aCo[i]-pos,ev);
			e.aRelCoMesh.push_back(x0);
			e.aRelCoMesh.push_back(y0);
		}
	}
	////////////////
	const std::vector<unsigned int>& aID_Loop = m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aID_Loop.size();iid_l++){
		unsigned int id_l = aID_Loop[iid_l];
    if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
	}
	return true;
FAILURE:
	e = e_old;
	return false;
}

// id_l is not exist(e.g. 0), check intersection for all edges
bool CCadObj2D::CheckIntersection_Loop(unsigned int id_l) const
{
	std::vector<unsigned int> aIdEdge;
  if( m_BRep.IsElemID(Cad::LOOP,id_l) ){
		for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(itrl.Begin();!itrl.IsEnd();itrl++){
				unsigned int id_e;   bool is_same_dir;
				if( !itrl.GetIdEdge(id_e,is_same_dir) ) continue;	// skip floating points
				if( itrl.IsEdge_BothSideSameLoop() && !is_same_dir ){ continue; }	// if both side is same loop add list onece
				aIdEdge.push_back(id_e);
			}
		}
	}
	else{ aIdEdge = this->GetAryElemID(Cad::EDGE); }  

	const unsigned int ne = aIdEdge.size();
	for(unsigned int ie=0;ie<ne;ie++){
    const CEdge2D& edge = this->GetEdge( aIdEdge[ie] );
		if( edge.IsCrossEdgeSelf() ){ return true; }
		const unsigned int ipo0 = edge.id_v_s; 
		const unsigned int ipo1 = edge.id_v_e;
    const Com::CBoundingBox2D& bb_i = edge.GetBoundingBox();		// get bounding box of edge_i
		for(unsigned int je=ie+1;je<ne;je++){
      const CEdge2D& e_j = this->GetEdge( aIdEdge[je] );
			const unsigned int jpo0 = e_j.id_v_s;
			const unsigned int jpo1 = e_j.id_v_e;			
			if( (ipo0-jpo0)*(ipo0-jpo1)*(ipo1-jpo0)*(ipo1-jpo1) != 0 ){ // there is no shared point
        const Com::CBoundingBox2D& bb_j = e_j.GetBoundingBox();				// get bounding box of edge_j
        if( !bb_i.IsIntersect(bb_j,min_clearance) ) continue;        
        const double dist = edge.Distance(e_j);
        if( dist < min_clearance ) return true;
        continue;
			}
      else if( ipo0 == jpo0 && ipo1 == jpo1){ if( edge.IsCrossEdge_ShareBothPoints(e_j,true ) == 1 ){ return true; } }
      else if( ipo0 == jpo1 && ipo1 == jpo0){ if( edge.IsCrossEdge_ShareBothPoints(e_j,false) == 1 ){ return true; } }
      else if( ipo0 == jpo0 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j, true, true)==1 ){ return true; } }
      else if( ipo0 == jpo1 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j, true,false)==1 ){ return true; } }
      else if( ipo1 == jpo0 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j,false, true)==1 ){ return true; } }
      else if( ipo1 == jpo1 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j,false,false)==1 ){ return true; } }
		}
	}  
	return false;
}

// id_l is not exist(e.g. 0), check intersection for all edges
bool CCadObj2D::CheckIntersection_EdgeAgainstLoop(const CEdge2D& edge,unsigned int id_l) const
{
	std::vector<unsigned int> aIdEdge;
  if( m_BRep.IsElemID(Cad::LOOP,id_l) ){
		for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(itrl.Begin();!itrl.IsEnd();itrl++){
				unsigned int id_e;   bool is_same_dir;
				if( !itrl.GetIdEdge(id_e,is_same_dir) ) continue;	// skip floating points
				if( itrl.IsEdge_BothSideSameLoop() && !is_same_dir ){ continue; }	// if both side is same loop add list onece
				aIdEdge.push_back(id_e);
			}
		}
	}
	else{ aIdEdge = this->GetAryElemID(Cad::EDGE); }  
  
  ////
	const unsigned int ne = aIdEdge.size();
  if( edge.IsCrossEdgeSelf() ){ return true; }
  const unsigned int ipo0 = edge.id_v_s; 
  const unsigned int ipo1 = edge.id_v_e;
  const Com::CBoundingBox2D& bb_i = edge.GetBoundingBox(); // get bounding box of edge_i  
  for(unsigned int je=0;je<ne;je++){
    const CEdge2D& e_j = this->GetEdge( aIdEdge[je] );
    const unsigned int jpo0 = e_j.id_v_s;
    const unsigned int jpo1 = e_j.id_v_e;			
    if( (ipo0-jpo0)*(ipo0-jpo1)*(ipo1-jpo0)*(ipo1-jpo1) != 0 ){ // there is no shared point
      const Com::CBoundingBox2D& bb_j = e_j.GetBoundingBox(); // get bounding box of e_j 
      if( !bb_i.IsIntersect(bb_j,min_clearance) ) continue;
      double dist = edge.Distance(e_j);
      if( dist < min_clearance ) return true;
      continue;      
    }
    else if( ipo0 == jpo0 && ipo1 == jpo1){ if( edge.IsCrossEdge_ShareBothPoints(e_j,true ) == 1 ){ return true; } }
    else if( ipo0 == jpo1 && ipo1 == jpo0){ if( edge.IsCrossEdge_ShareBothPoints(e_j,false) == 1 ){ return true; } }
    else if( ipo0 == jpo0 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j, true, true)==1 ){ return true; } }
    else if( ipo0 == jpo1 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j, true,false)==1 ){ return true; } }
    else if( ipo1 == jpo0 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j,false, true)==1 ){ return true; } }
    else if( ipo1 == jpo1 ){ if( edge.IsCrossEdge_ShareOnePoint(e_j,false,false)==1 ){ return true; } }
  }  
	return false;
}


// write to DXF file
bool CCadObj2D::WriteToFile_dxf(const std::string& file_name, double scale) const
{
	FILE *fp;
	if( (fp = ::fopen(file_name.c_str(),"w"))== NULL ){
		fclose(fp);
		assert(0);
		return false;
	}
  CBoundingBox2D bb;
	{	// Get Bounding Box of this Object
		const std::vector<unsigned int>& aIdEdge = this->m_EdgeSet.GetAry_ObjID();
		assert( aIdEdge.size() > 0 );
		for(unsigned int iid_e=0;iid_e<aIdEdge.size();iid_e++){
			const unsigned int id_e = aIdEdge[iid_e];
			const CEdge2D& edge = this->GetEdge(id_e);
      bb += edge.GetBoundingBox();
		}
	}
	// header section
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nHEADER\n");
	fprintf(fp, "  9\n$ACADVER\n  1\nAC1009\n");
	fprintf(fp, "  9\n$EXTMIN\n  10\n%lf\n  20\n%lf\n",bb.x_min*scale,bb.y_min*scale);
	fprintf(fp, "  9\n$EXTMAX\n  10\n%lf\n  20\n%lf\n",bb.x_max*scale,bb.y_max*scale);
	fprintf(fp, "  0\nENDSEC\n");
	// table section
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nTABLES\n");
	fprintf(fp, "  0\nENDSEC\n");
	// block section
	fprintf(fp, "  0\nSECTION\n");
	fprintf(fp, "  2\nBLOCKS\n");
	fprintf(fp, "  0\nENDSEC\n");
	// entity section
	fprintf(fp,"  0\nSECTION\n");
	fprintf(fp,"  2\nENTITIES\n");
	const std::vector<unsigned int>& aIdLoop = this->m_LoopSet.GetAry_ObjID();
	for(unsigned int iid_l=0;iid_l<aIdLoop.size();iid_l++){
		const unsigned int id_l = aIdLoop[iid_l];
		assert( m_LoopSet.IsObjID(id_l) );
    CBRepSurface::CItrLoop pItr = m_BRep.GetItrLoop(id_l);
		for(;;){
      for(;!pItr.IsEnd();pItr++){
				bool is_same_dir;   unsigned int id_e;
				if( !pItr.GetIdEdge(id_e,is_same_dir) ){ assert(0); fclose(fp); return false; }
				unsigned int id_vs = this->GetIdVertex_Edge(id_e, true );	assert( this->IsElemID(Cad::VERTEX,id_vs) );
        unsigned int id_ve = this->GetIdVertex_Edge(id_e, false); assert( this->IsElemID(Cad::VERTEX,id_ve) );
				const CVector2D& ps = this->GetVertexCoord(id_vs);
				const CVector2D& pe = this->GetVertexCoord(id_ve);
				if( this->GetEdgeCurveType(id_e) == 0 ){
					fprintf(fp,"  0\nLINE\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n",id_l);
					fprintf(fp,"  10\n%lf\n",ps.x*scale);
					fprintf(fp,"  20\n%lf\n",ps.y*scale);
					fprintf(fp,"  11\n%lf\n",pe.x*scale);
					fprintf(fp,"  21\n%lf\n",pe.y*scale);
				}
				else if( this->GetEdgeCurveType(id_e) == 1 ){ // Arc
					const CEdge2D& edge = this->GetEdge(id_e);
					CVector2D pc;	double r;
					edge.GetCenterRadius(pc,r);
					double d1, d2;
					{
						CVector2D vs = ps - pc;
						CVector2D ve = pe - pc;
						double ds = atan2(vs.y,vs.x); ds = ds * 180.0 / 3.14159265; if( ds < 0.0 ) ds += 360;
						double de = atan2(ve.y,ve.x); de = de * 180.0 / 3.14159265; if( de < 0.0 ) de += 360;
						if( edge.is_left_side ){ d1 = de; d2 = ds; }
						else{                    d1 = ds; d2 = de; }
					}
					fprintf(fp,"  0\nARC\n  8\n%d\n  6\nCONTINUOUS\n  62\n7\n  100\nAcDbCircle\n",id_l);
					fprintf(fp,"  10\n%lf\n",pc.x*scale);	// x coord
					fprintf(fp,"  20\n%lf\n",pc.y*scale);	// y coord
					fprintf(fp,"  40\n%lf\n",r*scale);	// radius
					fprintf(fp,"  100\nAcDbArc\n");
					fprintf(fp,"  50\n%lf\n",d1);
					fprintf(fp,"  51\n%lf\n",d2);
				}
        else if( this->GetEdgeCurveType(id_e) == 2 ){ // polyline
					const CEdge2D& edge = this->GetEdge(id_e);
					fprintf(fp,"  0\nPOLYLINE\n  8\n%d\n  6\nCONTINUOUS\n",id_l);
					fprintf(fp,"  10\n0.0\n");
					fprintf(fp,"  20\n0.0\n");
					fprintf(fp,"  30\n0.0\n");
					fprintf(fp,"  70\n8\n");
					fprintf(fp,"  66\n1\n");          
          ////
          const std::vector<double>& axys = edge.aRelCoMesh;
          assert( axys.size() % 2 == 0 );
          const unsigned int nno = axys.size()/2;
          const Com::CVector2D& po_s = this->GetVertexCoord( edge.id_v_s );
          const Com::CVector2D& po_e = this->GetVertexCoord( edge.id_v_e );
          Com::CVector2D v0 = po_e-po_s;
          Com::CVector2D v1(-v0.y,v0.x);          
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_s.x*scale, po_s.y*scale, 0.0);
          for(unsigned int ino=0;ino<nno;ino++){
            const Com::CVector2D& p = po_s + v0*axys[ino*2+0] + v1*axys[ino*2+1];
            fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", p.x*scale, p.y*scale, 0.0);
          }
          fprintf(fp,"  0\nVERTEX\n 8\n0\n 10\n%lf\n 20\n%lf\n 30\n%lf\n", po_e.x*scale, po_e.y*scale, 0.0);          
          fprintf(fp,"  0\nSEQEND\n");
				}        
			}
			if( !pItr.ShiftChildLoop() ) break;
		}
	}
	fprintf(fp, "  0\nENDSEC\n  0\nEOF\n");
	fclose(fp);
	return true;
}


bool CCadObj2D::Serialize( Com::CSerializer& arch )
{
	if( arch.IsLoading() ){	// load
		this->Clear();
		////////////////
		const unsigned int buff_size = 256; 
		char class_name[buff_size];
		arch.ReadDepthClassName(class_name,buff_size);
		if( strncmp(class_name,"CadObj2D",8) != 0 ) return false;
		int nv, ne, nl;
		{
			arch.Get("%d%d%d",&nv, &ne, &nl);
//			assert( nv>0 ); assert( ne>0 ); assert( nl>0 );
			m_VertexSet.Reserve(nv*2);
			m_EdgeSet.Reserve(ne*2);
			m_LoopSet.Reserve(nl*2);
		}
		arch.ShiftDepth(true);
		////////////////////////////////////////////////
    for(int iv=0;iv<nv;iv++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CVertex2D",9) == 0 );
			int id;		arch.Get("%d",&id);		assert( id>0 );
			double x,y;	arch.Get("%lf%lf",&x,&y);
			int tmp_id = m_VertexSet.AddObj( std::make_pair(id,CVertex2D(CVector2D(x,y) )) );
			assert(tmp_id==id);
		}
    for(int ie=0;ie<ne;ie++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CEdge2D",7) == 0 );
			int id;				arch.Get("%d",&id);					assert( id>0 );
			int id_v_s,id_v_e;	arch.Get("%d%d",&id_v_s,&id_v_e);	assert( id_v_s>0 && id_v_e>0 );
      int itype;			arch.Get("%d",&itype);				assert( itype == 0 || itype == 1 || itype == 2);
			int i_is_left_side;
			double dist;
			arch.Get("%d%lf",&i_is_left_side,&dist);
			assert( i_is_left_side == 0 || i_is_left_side == 1 );
			const bool is_left_side = (i_is_left_side != 0 );
      std::vector<double> aRelCo;
      int npo;
      arch.Get("%d",&npo);
      assert( npo >= 0 );
      for(unsigned int ipo=0;ipo<(unsigned int)npo;ipo++){
        double x,y;
        arch.Get("%lf%lf",&x,&y);
        aRelCo.push_back(x);
        aRelCo.push_back(y);
      }
      Cad::CEdge2D e(id_v_s,id_v_e);
      e.itype = itype;
      e.is_left_side = is_left_side;  
      e.dist = dist;
      e.aRelCoMesh = aRelCo;
      const int tmp_id = m_EdgeSet.AddObj( std::make_pair(id,e) );
			assert( tmp_id == id );            
		}
    for(int il=0;il<nl;il++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CLoop2D",7) == 0 );
      int id;       arch.Get("%d",&id);		assert( id>0 );
			int ilayer;		arch.Get("%d",&ilayer);
      double c[3];  arch.Get("%lf%lf%lf",&c[0],&c[1],&c[2]);
      CLoop2D l;
			l.ilayer = ilayer;
      l.m_color[0]=c[0];  l.m_color[1]=c[1];  l.m_color[2]=c[2];
      const int tmp_id = m_LoopSet.AddObj( std::make_pair(id,l) );
			assert( tmp_id == id );
		}
		m_BRep.Serialize(arch);
		this->AssertValid();
		arch.ShiftDepth(false);
		return true;
	}
	else{ // write
    // class name and size
		arch.WriteDepthClassName("CadObj2D");
		arch.Out("%d %d %d\n",m_VertexSet.GetAry_ObjID().size(), m_EdgeSet.GetAry_ObjID().size(),m_LoopSet.GetAry_ObjID().size());
		arch.ShiftDepth(true);
    { // print Vertex2D
			const std::vector<unsigned int>& id_ary = m_VertexSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_v = id_ary[iid];
				assert( m_VertexSet.IsObjID(id_v) );
				const CVertex2D& v = m_VertexSet.GetObj(id_v);
				arch.WriteDepthClassName("CVertex2D");
				arch.Out("%d\n",id_v);
				arch.Out("%lf %lf\n",v.point.x,v.point.y);
			}
		}
    { // print Edge2D
			const std::vector<unsigned int> id_ary = m_EdgeSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_e = id_ary[iid];
				assert( m_EdgeSet.IsObjID(id_e) );
				const CEdge2D& e = m_EdgeSet.GetObj(id_e);
				arch.WriteDepthClassName("CEdge2D");
				arch.Out("%d\n",id_e);
				arch.Out("%d %d\n",e.id_v_s,e.id_v_e);
				arch.Out("%d\n",e.itype);
        const unsigned int i_is_left_side = (e.is_left_side) ? 1 : 0;
				arch.Out("%d %lf\n",i_is_left_side,e.dist);
        const unsigned int n = e.aRelCoMesh.size()/2;
        arch.Out("%d\n",n);
        for(unsigned int i=0;i<n;i++){
          arch.Out("%lf %lf\n",e.aRelCoMesh[i*2+0], e.aRelCoMesh[i*2+1]);
        }
			}
		}    
    { // print Loop2D
      const std::vector<unsigned int> id_ary= m_LoopSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_l = id_ary[iid];
				assert( m_LoopSet.IsObjID(id_l) );
				const CLoop2D& l = m_LoopSet.GetObj(id_l);
				arch.WriteDepthClassName("CLoop2D");
				arch.Out("%d\n",id_l);
				arch.Out("%d\n",l.ilayer);
        arch.Out("%lf %lf %lf\n",l.m_color[0],l.m_color[1],l.m_color[2]);
			}
		}
		m_BRep.Serialize(arch);
		arch.ShiftDepth(false);
	}
	return true;
}

