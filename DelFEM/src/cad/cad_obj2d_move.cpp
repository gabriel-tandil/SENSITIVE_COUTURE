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
// CadObj2D.cpp : ２次元ＣＡＤモデルクラス(CCadObj2Dm)の実装
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

#include "delfem/cad/cad_elem2d.h"

#include "delfem/cad_obj2d_move.h"

using namespace Cad;
using namespace Com;


// id_eの辺ががメッシュなら，スムージングをかけて滑らかにする
//! id_eがメッシュなら，スムージングをかけて滑らかにする
bool CCadObj2D_Move::SmoothingPolylineEdge(unsigned int id_e, unsigned int nitr,
                                   const Com::CVector2D& po_c, double radius)
{
  if( !this->IsElemID(Cad::EDGE,id_e) ) return false;
  if( this->GetEdgeCurveType(id_e) != 2 ) return true;
  if( nitr == 0 ) return true;
  CEdge2D& edge = this->m_EdgeSet.GetObj(id_e);
  std::vector<double>& axys = edge.aRelCoMesh;
  assert( axys.size() % 2 == 0 );
  const unsigned int nno = axys.size()/2;
  std::vector<unsigned int> aIndNo;
  if( radius > 0 ){
    aIndNo.reserve(nno);
    const Com::CVector2D& po_s = this->GetVertexCoord( edge.id_v_s );
    const Com::CVector2D& po_e = this->GetVertexCoord( edge.id_v_e );
    Com::CVector2D v0 = po_e-po_s;
    Com::CVector2D v1(-v0.y,v0.x);
    for(unsigned int ino=0;ino<nno;ino++){
      const Com::CVector2D& p = po_s + v0*axys[ino*2+0] + v1*axys[ino*2+1];
      if( Com::SquareLength(p-po_c) < radius*radius ){ aIndNo.push_back(ino); }
    }
  }
  else{
    aIndNo.resize(nno);
    for(unsigned int ino=0;ino<nno;ino++){ aIndNo[ino] = ino; }
  }
  const double w = 0.8;   // strength of smoothing(0から1)
  for(unsigned int iitr=0;iitr<nitr;iitr++){
    for(unsigned int iindno=0;iindno<aIndNo.size();iindno++){
      unsigned int ino = aIndNo[iindno];
      double po0[2],po1[2],po2[2];
      if( ino==0 ){ po0[0]=0; po0[1]=0; }
      else{ po0[0]=axys[ino*2-2]; po0[1]=axys[ino*2-1]; }
      po1[0]=axys[ino*2+0]; po1[1]=axys[ino*2+1];
      if( ino==nno-1 ){ po2[0]=1; po2[1]=0; }
      else{ po2[0]=axys[ino*2+2]; po2[1]=axys[ino*2+3]; }
      ////////////////
      const double pom[2] = { (po0[0]+po2[0])*0.5, (po0[1]+po2[1])*0.5 };
      const double pod[2] = { (1.0-w)*po1[0]+w*pom[0], (1.0-w)*po1[1]+w*pom[1] };
      axys[ino*2+0] = pod[0];
      axys[ino*2+1] = pod[1];
    }
  }
  return true;
}

// id_eのカーブをdistに通るようにする
bool CCadObj2D_Move::DragArc(unsigned int id_e, const CVector2D& vec)
{  
	if( !this->IsElemID(Cad::EDGE,id_e) ) return false;
	if( this->GetEdgeCurveType(id_e) != 1 ) return true;
  
	bool is_left_side_old;
	double dist_old;
	this->GetCurve_Arc(id_e,is_left_side_old,dist_old);
  
  CEdge2D& edge = this->GetEdgeRef(id_e);
  assert( edge.itype == 1 );
  const double base_len = sqrt( SquareLength(edge.po_s,edge.po_e) );
  if( fabs(TriHeight(vec,edge.po_s,edge.po_e)) > base_len*0.02 ){
    Com::CVector2D& ps = edge.po_s;
    Com::CVector2D& pe = edge.po_e;
    CVector2D pc;
    Com::CenterCircumcircle(ps,pe,vec,pc);
    double dist = TriHeight(pc,ps,pe);
    edge.is_left_side = ( TriArea(ps,pe,vec) > 0 );
    edge.dist = dist;
  }
  else{
    return true;
  }
  
  /*
	if( tol > 0 ){
		const double dist90 = base_len*0.5;
		if( fabs(edge.dist) < tol ){ edge.dist = 0; }
		else if( edge.dist >  dist90-tol*2 && edge.dist <  dist90+tol*2 ){ edge.dist =  dist90; }
		else if( edge.dist > -dist90-tol*2 && edge.dist < -dist90+tol*2 ){ edge.dist = -dist90; }
	}
   */
  
  std::set<unsigned int> setIdL;
  for(CBRepSurface::CItrVertex itrv=m_BRep.GetItrVertex(edge.id_v_s);!itrv.IsEnd();itrv++){
    setIdL.insert( itrv.GetIdLoop() );
  }
  for(CBRepSurface::CItrVertex itrv=m_BRep.GetItrVertex(edge.id_v_e);!itrv.IsEnd();itrv++){
    setIdL.insert( itrv.GetIdLoop() );
  }
  
  for(std::set<unsigned int>::const_iterator itrl=setIdL.begin();itrl!=setIdL.end();itrl++){
    const unsigned int id_l = *itrl;
    if( !m_BRep.IsElemID(Cad::LOOP,id_l) ) continue;
    if( this->CheckLoop(id_l)!=0 ) goto FAILURE;
	}
	return true;
  
FAILURE:
	edge.dist = dist_old;
	edge.is_left_side = is_left_side_old;
	return true;  
}

bool CCadObj2D_Move::PreCompDragPolyline(unsigned int id_e, const Com::CVector2D& pick_pos)
{
  if( !this->IsElemID(Cad::EDGE,id_e) ) return false;
	if( this->GetEdgeCurveType(id_e) != 2 ) return true;  // exit if cad elem(ID:id_e) is not polyline     
  
  polyline.SetCadEdge(*this,id_e,pick_pos);
  
  return true;
}

//! if edge (ID:id_e) is an polyline make it go through point(dist)
bool CCadObj2D_Move::DragPolyline(unsigned int id_e, const Com::CVector2D& dist)
{
  if( !this->IsElemID(Cad::EDGE,id_e) ) return false;
	if( this->GetEdgeCurveType(id_e) != 2 ) return true;  // exit if cad elem(ID:id_e) is not polyline 
  
  polyline.Drag(*this,dist);
  
  return true;
}


// move all child loop
bool CCadObj2D_Move::MoveLoop(unsigned int id_l, const CVector2D& vec_del)
{
	if( !this->IsElemID(Cad::LOOP,id_l) ) return false;
  std::map<unsigned int, Com::CVector2D> map_vec_old;
  std::set<unsigned int> setIdL;  // check these loop for intersection detection
	for(CBRepSurface::CItrLoop itrl=m_BRep.GetItrLoop(id_l);!itrl.IsEndChild();itrl.ShiftChildLoop()){
		for(itrl.Begin();!itrl.IsEnd();itrl++){
			unsigned int id_v = itrl.GetIdVertex();
			if( map_vec_old.find(id_v) != map_vec_old.end() ) continue; // this point is already moved
			map_vec_old.insert( std::make_pair(id_v,this->GetVertexCoord(id_v)) );
			CVertex2D& ver = m_VertexSet.GetObj(id_v);
			ver.point.x += vec_del.x;
			ver.point.y += vec_del.y;
      unsigned int id_e0;   bool is_same_dir0;
      itrl.GetIdEdge(id_e0,is_same_dir0);
      if( !this->IsElemID(Cad::EDGE,id_e0) ) continue;  // this is point
      unsigned int id_l_l, id_l_r;
      this->GetIdLoop_Edge(id_l_l,id_l_r,id_e0);
      setIdL.insert(id_l_l);
      setIdL.insert(id_l_r);
		}		
	}
  for(std::set<unsigned int>::const_iterator itr=setIdL.begin();itr!=setIdL.end();itr++){
    unsigned int id_l0 = *itr;
    if( !this->IsElemID(Cad::LOOP,id_l0) ) continue;  // id_0 can be 0 if this loop(id_l) have open boundary
    if( this->CheckLoop(id_l0)!=0 ){ goto FAILURE; }
  }
  return true;
FAILURE:
  for(std::map<unsigned int,Com::CVector2D>::const_iterator itr = map_vec_old.begin();itr!=map_vec_old.end();itr++){
    unsigned int id_v = itr->first;
    CVertex2D& ver = m_VertexSet.GetObj(id_v);
    ver.point = itr->second;
	}
	return false;  
}

bool CCadObj2D_Move::MoveEdge(unsigned int id_e, const CVector2D& vec_del)
{
	assert( this->IsElemID(Cad::EDGE,id_e) );
	unsigned int id_v_s = this->GetIdVertex_Edge(id_e,true );
  unsigned int id_v_e = this->GetIdVertex_Edge(id_e,false);
  
	CVector2D vec_pre_s;
	{	// 点を動かす→駄目だったら元に戻す
		CVertex2D& ver = m_VertexSet.GetObj(id_v_s);
		vec_pre_s = ver.point;
		ver.point.x += vec_del.x;
		ver.point.y += vec_del.y;
	}

	CVector2D vec_pre_e;
	{	// 点を動かす→駄目だったら元に戻す
		CVertex2D& ver = m_VertexSet.GetObj(id_v_e);
		vec_pre_e = ver.point;
		ver.point.x += vec_del.x;
		ver.point.y += vec_del.y;
	}
	
	{	// Check Interfarance
    std::set<unsigned int> aIdL;
    for(CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id_v_s);!itrv.IsEnd();itrv++){
      const unsigned int id_l = itrv.GetIdLoop();
      if( this->IsElemID(Cad::LOOP,id_l) ){
        std::pair< std::set<unsigned int>::iterator, bool> res = aIdL.insert(id_l);
        if( !res.second ) continue; // this loop have already examined
        if( this->CheckLoop(id_l) != 0 ){ goto FAILURE; }
      }
    }
    for(CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id_v_e);!itrv.IsEnd();itrv++){
      const unsigned int id_l = itrv.GetIdLoop();
      if( this->IsElemID(Cad::LOOP,id_l) ){
        std::pair< std::set<unsigned int>::iterator, bool> res = aIdL.insert(id_l);
        if( !res.second ) continue; // this loop have already examined
        if( this->CheckLoop(id_l) != 0 ){ goto FAILURE; }
      }
    }    
    /*
		unsigned int id_l_l, id_l_r;
		this->GetIdLoop_Edge(id_l_l,id_l_r, id_e);
		if( this->IsElemID(Cad::LOOP,id_l_l) ){ 
			if( this->CheckLoop(id_l_l)!=0 ) goto FAILURE;
		}
		if( this->IsElemID(Cad::LOOP,id_l_r) ){ 
			if( this->CheckLoop(id_l_r)!=0 ) goto FAILURE;
		}
     */
	}
	return true;
	////////////////////////////////  
FAILURE:	// if the operation fails
	{	// 動かした点を元に戻す
		CVertex2D& ver = m_VertexSet.GetObj(id_v_s);
		ver.point = vec_pre_s;
	}
	{	// 動かした点を元に戻す
		CVertex2D& ver = m_VertexSet.GetObj(id_v_e);
		ver.point = vec_pre_e;
	}
	return false;
}

bool CCadObj2D_Move::MoveVertex(const unsigned int id_v, const CVector2D& vec)
{
	if( !m_VertexSet.IsObjID(id_v) ) return false;
  
  //	std::vector<unsigned int> aLoopID, aEdgeID;
  //	this->GetSurroundingObject_Vertex(id_v,aEdgeID,aLoopID);
  
	CVector2D vec_pre;
	{	// store point to move point back in case it fails
		CVertex2D& ver = m_VertexSet.GetObj(id_v);
		vec_pre = ver.point;
	}
  
	CVector2D dist = vec;
	{	// move point
		CVertex2D& ver = m_VertexSet.GetObj(id_v);
		ver.point = dist;
	}
  
  {
    CBRepSurface::CItrVertex itrv = m_BRep.GetItrVertex(id_v);
    if( itrv.CountEdge() == 0 ){	// move point inside loop
      unsigned int id_l = itrv.GetIdLoop();
      if( this->IsElemID(Cad::LOOP,id_l) ){
        const double dist = this->SignedDistPointLoop(id_l,vec,id_v); // ignore vtx(id_v) in the signd distance computation
        if( dist < min_clearance ){ goto FAILURE; }
      }
    }
    else{	// move point adjacent to loop
      std::set<unsigned int> aIdL;
      for(;!itrv.IsEnd();itrv++){
        const unsigned int id_l = itrv.GetIdLoop();
        if( this->IsElemID(Cad::LOOP,id_l) ){
          std::pair< std::set<unsigned int>::iterator, bool> res = aIdL.insert(id_l);
          if( !res.second ) continue; // this loop have already examined
          if( this->CheckLoop(id_l) != 0 ){ goto FAILURE; }
        }
      }
    }
  }
  return true;
	////////////////////////////////
FAILURE:	
	{	// reset the moved point
		CVertex2D& ver = m_VertexSet.GetObj(id_v);
		ver.point = vec_pre;
	}
	return false;  
}

bool CCadObj2D_Move::MoveVertex( const std::vector< std::pair<unsigned int,CVector2D> >& aIdVec )
{
	std::vector<CVector2D> aVecOld;
	for(unsigned int i=0;i<aIdVec.size();i++){
		unsigned int id_v = aIdVec[i].first;
		if( !m_VertexSet.IsObjID(id_v) ) goto FAILURE;
		CVertex2D& ver = m_VertexSet.GetObj(id_v);
		aVecOld.push_back( ver.point );
		ver.point = aIdVec[i].second;
	}

	{
		std::vector<unsigned int> aLoopID = this->GetAryElemID(Cad::LOOP);
		for(unsigned int iid_l=0;iid_l<aLoopID.size();iid_l++){
			unsigned int id_l = aLoopID[iid_l];
      if( this->CheckLoop(id_l) != 0 ) goto FAILURE;
		}
	}
	return true;
	////////////////////////////////
FAILURE:	
	// 動かした点を元に戻す
	for(unsigned int i=0;i<aVecOld.size();i++){
		unsigned int id_v = aIdVec[i].first;
		CVertex2D& ver = m_VertexSet.GetObj(id_v);
		ver.point = aVecOld[i];
	}
	return false;
}
