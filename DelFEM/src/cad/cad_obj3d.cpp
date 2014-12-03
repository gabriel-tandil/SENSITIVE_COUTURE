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
// CadObj3D.cpp : implementation of 3D cad class
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

#include "delfem/cad_obj3d.h"
#include "delfem/cad/cad_elem3d.h"
#include "delfem/cad/cad_elem2d.h"


using namespace Cad;
using namespace Com;


bool CCadObj3D::IsElemID(Cad::CAD_ELEM_TYPE itype,unsigned int id) const
{
  if(      itype == Cad::NOT_SET ){ return false; }
  else if( itype == Cad::VERTEX  ){ return m_VertexSet.IsObjID(id); }
  else if( itype == Cad::EDGE    ){ return m_EdgeSet.IsObjID(id);   }
  else if( itype == Cad::LOOP    ){ return m_LoopSet.IsObjID(id);   }  
	else{ assert(0); }
	return false;
}

const std::vector<unsigned int> CCadObj3D::GetAryElemID(Cad::CAD_ELEM_TYPE itype) const
{
	if(      itype == Cad::VERTEX ){ return m_VertexSet.GetAry_ObjID(); }
	else if( itype == Cad::EDGE   ){ return m_EdgeSet.GetAry_ObjID();   }
	else if( itype == Cad::LOOP   ){ return m_LoopSet.GetAry_ObjID();   }
	assert(0);
	std::vector<unsigned int> null_vec;
	return null_vec;
}

bool CCadObj3D::GetIdVertex_Edge(unsigned int &id_v_s, unsigned int& id_v_e, unsigned int id_e) const
{
  assert( m_BRep.IsElemID(Cad::EDGE,id_e) );
	return m_BRep.GetIdVertex_Edge(id_e,id_v_s,id_v_e);
}

const CVector3D& CCadObj3D::GetVertexCoord(unsigned int id_v) const
{
	assert( m_VertexSet.IsObjID(id_v) );
	const CVertex3D& v = m_VertexSet.GetObj(id_v);
	return v.point;
}

const Cad::CLoop3D& CCadObj3D::GetLoop(unsigned int id_l) const
{
	assert( m_LoopSet.IsObjID(id_l) );
	const CLoop3D& l = m_LoopSet.GetObj(id_l);  
  l.aEdge.clear();
  l.aIndEdge.clear();
  ////
  l.aIndEdge.push_back(0);
  for(CBRepSurface::CItrLoop itr=m_BRep.GetItrLoop(id_l);!itr.IsEndChild();itr.ShiftChildLoop()){
    for(;!itr.IsEnd();itr++){
      unsigned int id_e0; bool is_same_dir0; itr.GetIdEdge(id_e0,is_same_dir0);
      if( id_e0 == 0 ) continue;
      const Cad::CEdge2D& e2 = this->GetEdge2D(id_e0,id_l);
      l.aEdge.push_back( std::make_pair(e2,is_same_dir0) );
    }
    l.aIndEdge.push_back( l.aEdge.size() );
  }      
  l.bb_.isnt_empty = false;
  return l;
}

unsigned int Cad::CCadObj3D::AddCuboid(double len_x, double len_y, double len_z){  
  CVector3D pos_o(0,0,0), dir_x(len_x,0,0), dir_y(0,len_y,0), dir_z(0,0,len_z);
  unsigned int aIdV[8];
  {
    unsigned int aIdUV[8];
    for(unsigned int i=0;i<8;i++){ aIdUV[i] = m_BRep.AddVertex_Loop(0); }
    aIdV[0] = m_VertexSet.AddObj( std::make_pair(aIdUV[0],CVertex3D(pos_o                  )) );
    aIdV[1] = m_VertexSet.AddObj( std::make_pair(aIdUV[1],CVertex3D(pos_o      +dir_y      )) );
    aIdV[2] = m_VertexSet.AddObj( std::make_pair(aIdUV[2],CVertex3D(pos_o+dir_x+dir_y      )) );
    aIdV[3] = m_VertexSet.AddObj( std::make_pair(aIdUV[3],CVertex3D(pos_o+dir_x            )) );
    aIdV[4] = m_VertexSet.AddObj( std::make_pair(aIdUV[0],CVertex3D(pos_o            +dir_z)) );
    aIdV[5] = m_VertexSet.AddObj( std::make_pair(aIdUV[1],CVertex3D(pos_o      +dir_y+dir_z)) );
    aIdV[6] = m_VertexSet.AddObj( std::make_pair(aIdUV[2],CVertex3D(pos_o+dir_x+dir_y+dir_z)) );
    aIdV[7] = m_VertexSet.AddObj( std::make_pair(aIdUV[3],CVertex3D(pos_o+dir_x      +dir_z)) );    
  }
  ////
  unsigned int aIdE[12];
  {    
    aIdE[0] = m_BRep.ConnectVertex(CBRepSurface::CItrVertex(m_BRep,aIdV[0]),CBRepSurface::CItrVertex(m_BRep,aIdV[1]),true).id_e_add;
    aIdE[1] = m_BRep.ConnectVertex(CBRepSurface::CItrVertex(m_BRep,aIdV[1]),CBRepSurface::CItrVertex(m_BRep,aIdV[2]),true).id_e_add;
    aIdE[2] = m_BRep.ConnectVertex(CBRepSurface::CItrVertex(m_BRep,aIdV[2]),CBRepSurface::CItrVertex(m_BRep,aIdV[3]),true).id_e_add;    
    aIdE[3] = m_BRep.ConnectVertex(CBRepSurface::CItrVertex(m_BRep,aIdV[3]),CBRepSurface::CItrVertex(m_BRep,aIdV[0]),true).id_e_add;
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[0],CEdge3D()) );
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[1],CEdge3D()) );
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[2],CEdge3D()) );
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[3],CEdge3D()) );
  }
  for(unsigned int i=0;i<4;i++){
    CBRepSurface::CItrVertex itrv(m_BRep,aIdV[i]);
    if( itrv.GetIdLoop() != 0 ){ itrv++; }  assert( itrv.GetIdLoop()==0 );
    aIdE[i+4] = m_BRep.ConnectVertex(itrv,CBRepSurface::CItrVertex(m_BRep,aIdV[i+4]),true).id_e_add;
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[i+4],CEdge3D()) );    
  }
  {
    unsigned int id_l0 = m_BRep.GetIdLoop_Edge(aIdE[0],true);
    this->m_LoopSet.AddObj( std::make_pair(id_l0,CLoop3D(pos_o,Com::CVector3D(0,0,-1),Com::CVector3D(0,1,0))) );
  }
  for(unsigned int iiv=0;iiv<4;iiv++){
    CBRepSurface::CItrVertex itrv(m_BRep,aIdV[iiv+4]);
    if( itrv.GetIdLoop() != 0 ){ itrv++; }  assert( itrv.GetIdLoop()==0 );
    unsigned int jiv = ( iiv==3 ) ? 0 : iiv+1 ;
    aIdE[iiv+8] = m_BRep.ConnectVertex(itrv,CBRepSurface::CItrVertex(m_BRep,aIdV[jiv+4]),true).id_e_add;
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[iiv+8],CEdge3D()) );
    unsigned int id_l0 = m_BRep.GetIdLoop_Edge(aIdE[iiv+8],true);
    {
      const Com::CVector3D& v0 = this->GetVertexCoord(aIdV[iiv]);
      const Com::CVector3D& v1 = this->GetVertexCoord(aIdV[jiv]); 
      Com::CVector3D x = (v0-v1);        x.Normalize();
      Com::CVector3D n = Cross(x,dir_z); n.Normalize();
      this->m_LoopSet.AddObj( std::make_pair(id_l0,CLoop3D(v0,n,x)) );
    }
  }
  {
    unsigned int id_l_add = m_BRep.SealHole(aIdE[11],false);
    this->m_LoopSet.AddObj( std::make_pair(id_l_add,CLoop3D(pos_o+dir_z,Com::CVector3D(0,0,1),Com::CVector3D(0,1,0))) );
  }
	assert( this->AssertValid() == 0);
  return 0;
	////////////////////////////////
	// failure procedure
FAIL_ADD_POLYGON_INSIDE_LOOP :
  return 0;
}

unsigned int Cad::CCadObj3D::AddPolygon( const std::vector<Com::CVector3D>& aVec, unsigned int id_l)
{
//  std::cout << "AddPolygon" << std::endl;
  if( !this->IsElemID(Cad::LOOP,id_l) ) return 0;
  if( aVec.size() < 3 ) return 0;
  const unsigned int npo = aVec.size();
  std::vector<unsigned int> aIdV(npo,0);
  for(unsigned int iiv=0;iiv<npo;iiv++){
    unsigned int id_uv = m_BRep.AddVertex_Loop(id_l);
    aIdV[iiv] = m_VertexSet.AddObj( std::make_pair(id_uv,CVertex3D(aVec[iiv])) );
  }
  std::vector<unsigned int> aIdE(npo,0);
  for(unsigned int iiv=0;iiv<npo;iiv++){
    unsigned int jiv= (iiv==npo-1) ? 0 : iiv+1;
    CBRepSurface::CItrVertex itrv0(m_BRep,aIdV[iiv]); assert( itrv0.GetIdLoop()==id_l );
    CBRepSurface::CItrVertex itrv1(m_BRep,aIdV[jiv]); assert( itrv1.GetIdLoop()==id_l );
    aIdE[iiv] = m_BRep.ConnectVertex(itrv0,itrv1,true).id_e_add;
    this->m_EdgeSet.AddObj( std::make_pair(aIdE[iiv],CEdge3D()) );    
  }
  unsigned int id_l_add = m_BRep.GetIdLoop_Edge(aIdE[npo-1],true);
  {
    this->m_LoopSet.AddObj( std::make_pair(id_l_add,m_LoopSet.GetObj(id_l)) );
  }
  assert( this->AssertValid() == 0);
  return id_l_add;
}

class CLiftEdgeProp{
public:
  CLiftEdgeProp(unsigned int id_v, unsigned int id_e, bool is_same_dir, bool is_vert){
    this->id_v_out = id_v;
    this->id_e=id_e;
    this->is_same_dir = is_same_dir;
    this->is_vert = is_vert;
  }
public:
  unsigned int id_v_out;  
  unsigned int id_v_in;
  unsigned int id_e;
  bool is_same_dir;
  bool is_vert;
};

Cad::CEdge2D CCadObj3D::GetEdge2D(unsigned int id_e, unsigned int id_l) const{
  assert( this->IsElemID(Cad::EDGE, id_e) );
  assert( this->IsElemID(Cad::LOOP, id_l) );  
  const CLoop3D& l = m_LoopSet.GetObj(id_l);  // don't use GetLoop(). GetLoop() call this function and cause endless loop
  unsigned int id_vs, id_ve;
  m_BRep.GetIdVertex_Edge(id_e, id_vs,id_ve);
  CEdge2D e(id_vs,id_ve);  
  e.po_s = l.Project( this->GetVertexCoord(id_vs) );
  e.po_e = l.Project( this->GetVertexCoord(id_ve) );  
  return e;
}


bool CCadObj3D::FindIntersectionEdge
(unsigned int id_l,const Com::CVector2D& pos, const Com::CVector2D& poe,
 unsigned int& id_e_nearest, bool& is_same_dir_e_near,
 Com::CVector2D& p_nearest)
{  
  id_e_nearest = 0;
  Com::CVector2D dir = poe-pos; dir.Normalize();
  double sqdist = -1;
  for(Cad::CBRepSurface::CItrLoop pItr = m_BRep.GetItrLoop(id_l);!pItr.IsEndChild();pItr.ShiftChildLoop()){
    for(;!pItr.IsEnd();pItr++){
      unsigned int id_e;  bool is_same_dir;
      pItr.GetIdEdge(id_e,is_same_dir);
      const Cad::CEdge2D& edge = this->GetEdge2D(id_e,id_l);
      Com::CVector2D sec;
      if( !edge.GetNearestIntersectionPoint_AgainstHalfLine(sec,pos,dir) ) continue;
      if( sqdist < 0 || (sec-pos).SqLength() < sqdist ){
        id_e_nearest = id_e;
        is_same_dir_e_near = is_same_dir;
        p_nearest = sec;
        sqdist = (sec-pos).SqLength();
      }
    }
  }        
  if( !this->IsElemID(Cad::EDGE,id_e_nearest) ){ return false; }
  return true;
}


bool CCadObj3D::CheckIsPointInside_ItrLoop(CBRepSurface::CItrLoop& itrl, const CVector2D& point) const
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
			const CEdge2D& e = this->GetEdge2D(id_e,itrl.GetIdLoop());
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



unsigned int Cad::CCadObj3D::AddRectLoop(unsigned int id_l, const Com::CVector2D& qs, const Com::CVector2D& qe)
{
  Com::CVector2D p0(qs), p1(qe.x,qs.y), p2(qe), p3(qs.x,qe.y);
  if( (qs-qe).x * (qs-qe).y < 0 ){ 
    p1 = Com::CVector2D(qs.x,qe.y);
    p3 = Com::CVector2D(qe.x,qs.y);    
  }
  assert( this->IsElemID(Cad::LOOP, id_l) );  
  CBRepSurface::CItrLoop itrl = m_BRep.GetItrLoop(id_l);
  bool in0 = CheckIsPointInside_ItrLoop(itrl,p0);
  bool in1 = CheckIsPointInside_ItrLoop(itrl,p1);  
  bool in2 = CheckIsPointInside_ItrLoop(itrl,p2);    
  bool in3 = CheckIsPointInside_ItrLoop(itrl,p3);
  if( !in0 ) return 0;
  const CLoop3D& l = this->GetLoop(id_l);  
  if( in1 && in2 && in3 ){ 
    std::vector<CVector3D> aPo;
    aPo.push_back( l.UnProject(p0) );
    aPo.push_back( l.UnProject(p1) );        
    aPo.push_back( l.UnProject(p2) );    
    aPo.push_back( l.UnProject(p3) );
    return this->AddPolygon(aPo,id_l);
  }
  if( in1 && !in2 && !in3 ){
    unsigned int id_e7; bool is_same_dir7; Com::CVector2D p7;
    FindIntersectionEdge(id_l,p0,p3,id_e7,is_same_dir7,p7);
    unsigned int id_v7 = this->AddPoint(Cad::EDGE,id_e7,l.UnProject(p7));    
    unsigned int id_e5; bool is_same_dir5; Com::CVector2D p5;
    FindIntersectionEdge(id_l,p1,p2,id_e5,is_same_dir5,p5);    
    unsigned int id_v5 = this->AddPoint(Cad::EDGE,id_e5,l.UnProject(p5));    
    unsigned int id_v0 = this->AddPoint(Cad::LOOP,id_l, l.UnProject(p0));
    unsigned int id_v1 = this->AddPoint(Cad::LOOP,id_l, l.UnProject(p1));    
    this->ConnectVertex(id_v7,id_v0);
    this->ConnectVertex(id_v0,id_v1);
    this->ConnectVertex(id_v1,id_v5);
  }
  if( !in1 && !in2 && in3 ){
    unsigned int id_e6; bool is_same_dir6; Com::CVector2D p6;
    FindIntersectionEdge(id_l,p3,p2,id_e6,is_same_dir6,p6);
    unsigned int id_v6 = this->AddPoint(Cad::EDGE,id_e6,l.UnProject(p6));    
    unsigned int id_e4; bool is_same_dir4; Com::CVector2D p4;
    FindIntersectionEdge(id_l,p0,p1,id_e4,is_same_dir4,p4);    
    unsigned int id_v4 = this->AddPoint(Cad::EDGE,id_e4,l.UnProject(p4));    
    unsigned int id_v3 = this->AddPoint(Cad::LOOP,id_l, l.UnProject(p3));
    unsigned int id_v0 = this->AddPoint(Cad::LOOP,id_l, l.UnProject(p0));    
    this->ConnectVertex(id_v6,id_v3);
    this->ConnectVertex(id_v3,id_v0);
    this->ConnectVertex(id_v0,id_v4);
  }  
  if( !in1 && !in3 ){
    unsigned int id_e7; bool is_same_dir7; Com::CVector2D p7;
    FindIntersectionEdge(id_l,p0,p3,id_e7,is_same_dir7,p7);
    unsigned int id_v7 = this->AddPoint(Cad::EDGE,id_e7,l.UnProject(p7));    
    unsigned int id_e4; bool is_same_dir4; Com::CVector2D p4;
    FindIntersectionEdge(id_l,p0,p1,id_e4,is_same_dir4,p4);    
    unsigned int id_v4 = this->AddPoint(Cad::EDGE,id_e4,l.UnProject(p4));        
    unsigned int id_v0 = this->AddPoint(Cad::LOOP,id_l, l.UnProject(p0));        
    this->ConnectVertex(id_v7,id_v0);
    this->ConnectVertex(id_v0,id_v4);    
  }
  return 0;
}


void Cad::CCadObj3D::LiftLoop(unsigned int id_l, Com::CVector3D dir)
{
  if( !this->m_BRep.IsElemID(Cad::LOOP,id_l) ) return;
  std::vector< std::vector<CLiftEdgeProp> > LoopEV;
  {
    Com::CVector3D udir = dir;
    udir.Normalize();
    for(CBRepSurface::CItrLoop itr=m_BRep.GetItrLoop(id_l);!itr.IsEndChild();itr.ShiftChildLoop()){
      unsigned int iiul = LoopEV.size();
      LoopEV.resize(iiul+1);
      for(itr.Begin();!itr.IsEnd();itr++){
        unsigned int id_v = itr.GetIdVertex();
        unsigned int id_e0; bool is_same_dir0; itr.GetIdEdge(id_e0,is_same_dir0);
        unsigned int id_l0 = m_BRep.GetIdLoop_Edge(id_e0,!is_same_dir0); 
        assert( id_l0 != id_l );
        assert( m_LoopSet.IsObjID(id_l0) );
        const Com::CVector3D& n0 =  m_LoopSet.GetObj(id_l0).normal;
        bool is_vert0 = fabs( Dot(n0,dir) ) < 1.0e-5;
        LoopEV[iiul].push_back( CLiftEdgeProp(id_v,id_e0,is_same_dir0,is_vert0) );
      }
    }
  }
  for(unsigned int iiul=0;iiul<LoopEV.size();iiul++){    
    std::vector<CLiftEdgeProp>& aIdE = LoopEV[iiul];
    const unsigned int npo = aIdE.size();  
    for(unsigned int iiv=0;iiv<npo;iiv++){
      unsigned int jiv = (iiv==0) ? npo-1 : iiv-1;
      unsigned int id_v_out = aIdE[iiv].id_v_out;
      if( !aIdE[iiv].is_vert && !aIdE[jiv].is_vert ){
        unsigned int id_uv = m_BRep.AddVertex_Loop(id_l);
        Com::CVector3D vout = this->GetVertexCoord(id_v_out);
        aIdE[iiv].id_v_in = m_VertexSet.AddObj( std::make_pair(id_uv,CVertex3D(vout+dir)) );
        ////
        CBRepSurface::CItrVertex itrv0(m_BRep,id_v_out); 
        for(;itrv0.GetIdLoop()!=id_l&&!itrv0.IsEnd();itrv0++){} assert( itrv0.GetIdLoop()==id_l );
        CBRepSurface::CItrVertex itrv1(m_BRep,aIdE[iiv].id_v_in); assert( itrv1.GetIdLoop()==id_l );
        unsigned int id_e_add0 = m_BRep.ConnectVertex(itrv0,itrv1,true).id_e_add;
        this->m_EdgeSet.AddObj( std::make_pair(id_e_add0,CEdge3D()) );
      }
      else if( aIdE[iiv].is_vert &&  aIdE[jiv].is_vert ){
        aIdE[iiv].id_v_in = id_v_out;      
        CVertex3D& ver = m_VertexSet.GetObj(aIdE[iiv].id_v_in);
        ver.point += dir;
      }
      else {
        unsigned int id_e_insert;
        {
          CBRepSurface::CItrVertex itrv(m_BRep,aIdE[iiv].id_v_out); 
          for(;itrv.GetIdLoop()!=id_l&&!itrv.IsEnd();itrv++){} assert( itrv.GetIdLoop()==id_l );
          unsigned int id_e0; bool is_same_dir0; itrv.GetIdEdge_Ahead( id_e0, is_same_dir0);
          unsigned int id_e1; bool is_same_dir1; itrv.GetIdEdge_Behind(id_e1, is_same_dir1);
          if(  aIdE[iiv].is_vert ){          id_e_insert = id_e1; }
          else{ assert( aIdE[jiv].is_vert ); id_e_insert = id_e0; }
        }     
        unsigned int id_uv = m_BRep.AddVertex_Edge(id_e_insert);
        Com::CVector3D vout = this->GetVertexCoord(id_v_out);
        aIdE[iiv].id_v_in = m_VertexSet.AddObj( std::make_pair(id_uv,CVertex3D(vout+dir)) );
        /////
        unsigned int id_e_add;
        {
          CBRepSurface::CItrVertex itrv(m_BRep,aIdE[iiv].id_v_in);
          unsigned int id_e0; bool is_same_dir0; itrv.GetIdEdge_Ahead( id_e0, is_same_dir0);
          unsigned int id_e1; bool is_same_dir1; itrv.GetIdEdge_Behind(id_e1, is_same_dir1);
          id_e_add = (id_e0 == id_e_insert) ? id_e1 : id_e0;
          assert( id_e_add != id_e_insert );          
        }
        this->m_EdgeSet.AddObj( std::make_pair(id_e_add,CEdge3D()) );
      }
    }
  }
      
  for(unsigned int iiul=0;iiul<LoopEV.size();iiul++){    
    std::vector<CLiftEdgeProp>& aIdE = LoopEV[iiul];
    const unsigned int npo = aIdE.size();      
    for(unsigned int iiv=0;iiv<npo;iiv++){
      if( aIdE[iiv].is_vert ) continue;
      unsigned int jiv = (iiv==npo-1) ? 0 : iiv+1;
      CBRepSurface::CItrVertex itrv0(m_BRep,aIdE[iiv].id_v_in); 
      for(;itrv0.GetIdLoop()!=id_l&&!itrv0.IsEnd();itrv0++){} assert( itrv0.GetIdLoop()==id_l );
      CBRepSurface::CItrVertex itrv1(m_BRep,aIdE[jiv].id_v_in); 
      for(;itrv1.GetIdLoop()!=id_l&&!itrv1.IsEnd();itrv1++){} assert( itrv1.GetIdLoop()==id_l );
      CBRepSurface::CResConnectVertex res_cv = m_BRep.ConnectVertex(itrv0,itrv1,iiul==0);
      const unsigned int id_e_add0 = res_cv.id_e_add;
      this->m_EdgeSet.AddObj( std::make_pair(id_e_add0,CEdge3D()) );
      unsigned int id_l_add = m_BRep.GetIdLoop_Edge(id_e_add0,res_cv.is_left_l_add);
      assert( id_l_add != id_l );
      {
        const Com::CVector3D& v0 = this->GetVertexCoord(aIdE[iiv].id_v_out);
        const Com::CVector3D& v1 = this->GetVertexCoord(aIdE[jiv].id_v_out); 
        Com::CVector3D x = v1-v0;        x.Normalize();
        Com::CVector3D n = Cross(x,dir); n.Normalize();
        this->m_LoopSet.AddObj( std::make_pair(id_l_add,CLoop3D(v0,n,x)) );
      }
    }
  }
  {
    Cad::CLoop3D& l = this->m_LoopSet.GetObj(id_l);
    l.org += dir;
  }
  assert( AssertValid() == 0 );
}

unsigned int Cad::CCadObj3D::AddPoint( Cad::CAD_ELEM_TYPE type, unsigned int id_elem, const Com::CVector3D& vec_add )
{
  if( type == Cad::EDGE ){
    const unsigned int id_e = id_elem;
    if( !this->IsElemID(Cad::EDGE,id_e) ) return 0;
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
			const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex3D(vec_add)) );
			assert( tmp_id == (int)id_v_add );
		}
		{	// add edge
			const int tmp_id = m_EdgeSet.AddObj( std::make_pair(id_e_add,CEdge3D()) );
			assert( tmp_id == (int)id_e_add );
		}
		assert( this->AssertValid() == 0 );    
    return id_v_add;
  }
  if( type == Cad::LOOP ){
    const unsigned int id_l = id_elem;    
    if( !this->IsElemID(Cad::EDGE,id_l) ) return 0;
    /*
    { // check this point is inside the loop with clearance
      const double dist = this->SignedDistPointLoop(id_l,vec);
      if( dist < this->min_clearance ){ return res; }
    }
     */
		unsigned int id_v_add = m_BRep.AddVertex_Loop(id_l);
		const int tmp_id = m_VertexSet.AddObj( std::make_pair(id_v_add,CVertex3D(vec_add)) );
		assert( tmp_id == (int)id_v_add );
		assert( this->AssertValid()==0 );
    return id_v_add;    
  }
  return 0;
}

CBRepSurface::CResConnectVertex Cad::CCadObj3D::ConnectVertex(unsigned int id_v1, unsigned int id_v2)
{
  CBRepSurface::CResConnectVertex res;
  res.id_v1 = id_v1; res.id_v2 = id_v2;
  unsigned int max_id=0;
  {
    const std::vector<unsigned int>& aIdL = m_BRep.GetAryElemID(Cad::LOOP);
    for(unsigned int iil=0;iil<aIdL.size();iil++){ 
      max_id = ( aIdL[iil] > max_id ) ? aIdL[iil] : max_id;
    }
  }
  std::vector<int> aflg(max_id+1,0);
  CBRepSurface::CItrVertex itrv0 = m_BRep.GetItrVertex(id_v1);
  for(;!itrv0.IsEnd();itrv0++){ 
    unsigned int id_l = itrv0.GetIdLoop();
    assert( m_BRep.IsElemID(Cad::LOOP,id_l) );
    assert( aflg[id_l] == 0 );
    aflg[id_l]=1; 
  }
  unsigned int id_l =0;
  CBRepSurface::CItrVertex itrv1 = m_BRep.GetItrVertex(id_v2);
  for(;!itrv1.IsEnd();itrv1++){ 
    unsigned int id_l0 = itrv1.GetIdLoop();
    assert( m_BRep.IsElemID(Cad::LOOP,id_l0) );
    if( aflg[id_l0] == 1 ){
      id_l = id_l0;
      break;
    }
  }  
  res.id_l = id_l;
  if( id_l == 0 ) return res;
  for(itrv0.Begin();!itrv0.IsEnd();itrv0++){ 
    unsigned int id_l0 = itrv0.GetIdLoop();
    assert( m_BRep.IsElemID(Cad::LOOP,id_l0) );
    if( id_l0 == id_l ) break;
  }  
  assert( itrv0.GetIdLoop() == id_l );
  assert( itrv1.GetIdLoop() == id_l );
  res = m_BRep.ConnectVertex(itrv0,itrv1, true);
  { // register the edge
		this->m_EdgeSet.AddObj( std::make_pair(res.id_e_add,CEdge3D()) );
	}  
  if( res.id_l_add != id_l && res.id_l_add != 0 ){
		CLoop3D loop_add = this->m_LoopSet.GetObj(id_l);
		this->m_LoopSet.AddObj( std::make_pair(res.id_l_add,loop_add) );
	}  
  assert( this->AssertValid() == 0 );
  return res;
}


unsigned int Cad::CCadObj3D::AssertValid()
{
  if( !m_BRep.AssertValid() ){ return 6; }
  {
    std::vector<unsigned int> aIdL = m_BRep.GetAryElemID(Cad::LOOP);
    for(unsigned int iidl=0;iidl<aIdL.size();iidl++){ if( !this->m_LoopSet.IsObjID(aIdL[iidl]  ) ) return 7; }
  }
  {
    std::vector<unsigned int> aIdE = m_BRep.GetAryElemID(Cad::EDGE);
    for(unsigned int iide=0;iide<aIdE.size();iide++){ if( !this->m_EdgeSet.IsObjID(aIdE[iide]  ) ) return 8; }
  }
  {
    std::vector<unsigned int> aIdV = m_BRep.GetAryElemID(Cad::VERTEX);
    for(unsigned int iidv=0;iidv<aIdV.size();iidv++){ if( !this->m_VertexSet.IsObjID(aIdV[iidv]) ) return 9; }
  }                                         
  return 0;
}

