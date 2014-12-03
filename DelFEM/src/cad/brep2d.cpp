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
// brep2d.cpp : ２次元位相クラス(Cad::CBrep2D)の実装
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

//#include "delfem/cad_obj2d.h"
#include "delfem/cad/brep2d.h"
#include "delfem/cad/brep.h"

using namespace Cad;
//using namespace Com;

CBRepSurface::CItrLoop::CItrLoop(const CBRepSurface& ptr_cad_2d, unsigned int id_l)
	: m_pBRep2D(ptr_cad_2d)
{
	is_valid = false;
	m_IdUL = 0;
	{
		const std::map<unsigned int, unsigned int>::const_iterator itr =  m_pBRep2D.map_l2ul.find(id_l);
		assert( itr != m_pBRep2D.map_l2ul.end() );
		m_IdUL = itr->second;
	}
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul.id == m_IdUL );
	m_IdHE = ul.id_he;
	is_initial = true;
	is_valid = true;
	is_end_child = false;
}

CBRepSurface::CItrLoop::CItrLoop(const CBRepSurface& ptr_cad_2d, unsigned int id_he, unsigned int id_ul)
	: m_pBRep2D(ptr_cad_2d)
{
	is_valid = false;
	m_IdUL = id_ul;
	m_IdHE = id_he;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	assert( he.id_ul == m_IdUL );
	is_initial = false;
	is_valid = true;
	is_end_child = false;
}

void CBRepSurface::CItrLoop::Begin(){
	if( !is_valid ) return;
	is_initial = true;
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul.id == m_IdUL );
	m_IdHE = ul.id_he;
}

void CBRepSurface::CItrLoop::operator++(){
	if( !is_valid ) return;
	is_initial = false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	m_IdHE = he.id_he_f;
	return;
}

void CBRepSurface::CItrLoop::operator++(int n){
	if( !is_valid ) return;
	is_initial = false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	m_IdHE = he.id_he_f;
	return;
}

bool CBRepSurface::CItrLoop::IsEnd() const {
	if( is_initial ) return false;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul.id == m_IdUL );
	if( m_IdHE == ul.id_he ) return true;
	return false;
}

bool CBRepSurface::CItrLoop::GetIdEdge(unsigned int& id_e, bool& is_same_dir) const {
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	id_e = he.id_e;
	is_same_dir = he.is_same_dir;
	if( id_e == 0 ){ return false; }
	return true;
}

unsigned int CBRepSurface::CItrLoop::GetIdVertex() const {
	if( !is_valid ) return 0;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	assert( m_pBRep2D.m_BRep.IsID_UseVertex(he.id_uv ) );
	const CUseVertex& uv = m_pBRep2D.m_BRep.GetUseVertex( he.id_uv );
	const unsigned int id_v = uv.id_v;
	m_pBRep2D.IsElemID(Cad::VERTEX,id_v);
	return id_v;
}

unsigned int CBRepSurface::CItrLoop::GetIdVertex_Ahead() const {
	if( !is_valid ) return 0;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_f = he.id_he_f;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_f) );
	const CHalfEdge& he_f = m_pBRep2D.m_BRep.GetHalfEdge(id_he_f);
	const unsigned int id_uv_f = he_f.id_uv;
	const CUseVertex& uv_f = m_pBRep2D.m_BRep.GetUseVertex( id_uv_f );
	const unsigned int id_v = uv_f.id_v;
	m_pBRep2D.IsElemID(Cad::VERTEX,id_v);
	return id_v;
}

unsigned int CBRepSurface::CItrLoop::GetIdVertex_Behind() const {
	if( !is_valid ) return 0;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_b = he.id_he_b;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_b) );
	const CHalfEdge& he_b = m_pBRep2D.m_BRep.GetHalfEdge(id_he_b);
	const unsigned int id_uv_b = he_b.id_uv;
	const CUseVertex& uv_b = m_pBRep2D.m_BRep.GetUseVertex( id_uv_b );
	const unsigned int id_v = uv_b.id_v;
	m_pBRep2D.IsElemID(Cad::VERTEX,id_v);
	return id_v;
}

bool CBRepSurface::CItrLoop::ShiftChildLoop()
{
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul_p = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul_p.id == m_IdUL );
	if( ul_p.id_ul_c == 0 ){ 
		is_end_child = true; 
		return false;
	}
	m_IdUL = ul_p.id_ul_c;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul_c = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul_c.id == m_IdUL );
	m_IdHE = ul_c.id_he;
	is_initial = true;
	is_valid = true;
	return true;
}

unsigned int CBRepSurface::CItrLoop::GetType() const // 0:浮遊点 1:浮遊辺 2:面積がある
{
	return m_pBRep2D.TypeUseLoop(this->m_IdUL);
}

unsigned int CBRepSurface::CItrLoop::GetIdLoop() const
{
	if( !is_valid ) return 0;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul.id == m_IdUL );
	return ul.id_l;
}

bool CBRepSurface::CItrLoop::IsParent() const
{
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
	assert( ul.id == m_IdUL );
	return (ul.id_ul_p == m_IdUL);
}

bool CBRepSurface::CItrLoop::IsEdge_BothSideSameLoop() const
{
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_o = he.id_he_o;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_o) );
	const CHalfEdge& he_o = m_pBRep2D.m_BRep.GetHalfEdge(id_he_o);
	return he.id_ul == he_o.id_ul;
}

unsigned int CBRepSurface::CItrLoop::CountVertex_UseLoop() const
{
	if( !is_valid ) return 0;
	unsigned int id_he_ini;
	{
		assert( m_pBRep2D.m_BRep.IsID_UseLoop(m_IdUL) );
		const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(m_IdUL);
		id_he_ini = ul.id_he;
	}
	unsigned int id_he = id_he_ini;
	unsigned int icnt = 0;
	for(;;){
		icnt++;
		const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(id_he);
		const unsigned int id_he_f = he.id_he_f;
		if( id_he_f == id_he_ini ) break;
		{
			assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_f) );
			const CHalfEdge& he_f = m_pBRep2D.m_BRep.GetHalfEdge(id_he_f);
			assert( he_f.id_ul == m_IdUL );
		}
		id_he = id_he_f;
	}
	return icnt;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// 頂点周りのループを巡ることができるイテレータ
CBRepSurface::CItrVertex::CItrVertex(const CBRepSurface& ptr_cad_2d, unsigned int id_v)
	: m_pBRep2D(ptr_cad_2d)
{
	is_valid = false;
	m_IdUV = id_v;
	const CUseVertex& uv = m_pBRep2D.m_BRep.GetUseVertex(m_IdUV);
	assert( uv.id == m_IdUV );
	m_IdHE = uv.id_he;
	is_initial = true;
	is_valid = true;
}

// 反時計周りに頂点まわりをめぐる
void CBRepSurface::CItrVertex::operator++(){
	if( !is_valid ) return;
	is_initial = false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_b = he.id_he_b;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_b) );
	const CHalfEdge& he_b = m_pBRep2D.m_BRep.GetHalfEdge(id_he_b);
	m_IdHE = he_b.id_he_o;
	return;
}

// ダミーのオペレータ(++と同じ働き)
void CBRepSurface::CItrVertex::operator++(int n){
	if( !is_valid ) return;
	is_initial = false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_b = he.id_he_b;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_b) );
	const CHalfEdge& he_b = m_pBRep2D.m_BRep.GetHalfEdge(id_he_b);
	m_IdHE = he_b.id_he_o;
	return;
}

bool CBRepSurface::CItrVertex::GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const {	
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	id_e = he.id_e;
	is_same_dir = he.is_same_dir;
	return true;
}

// 頂点周りの辺のIDと、その辺の始点がid_vと一致しているかどうか
bool CBRepSurface::CItrVertex::GetIdEdge_Ahead(unsigned int& id_e, bool& is_same_dir) const {
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	const unsigned int id_he_b = he.id_he_b;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_b) );
	const CHalfEdge& he_b = m_pBRep2D.m_BRep.GetHalfEdge(id_he_b);
	id_e = he_b.id_e;
	is_same_dir = !he_b.is_same_dir;
	return true;
}

// ループのIDを得る
unsigned int CBRepSurface::CItrVertex::GetIdLoop() const 
{
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	unsigned int id_ul = he.id_ul;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(id_ul) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(id_ul);
	return ul.id_l;
}

// 面周りの辺が一周したらtrueを返す
bool CBRepSurface::CItrVertex::IsEnd() const
{
	if( is_initial ) return false;
	assert( m_pBRep2D.m_BRep.IsID_UseVertex(m_IdUV) );
	// いちいちUVを取ってくる必要が無いように高速化が必要
	const CUseVertex& uv = m_pBRep2D.m_BRep.GetUseVertex(m_IdUV); 
	if( m_IdHE == uv.id_he ) return true;
	return false;
}


void CBRepSurface::CItrVertex::Begin(){
	if( !is_valid ) return;
	is_initial = true;
	const CUseVertex& uv = m_pBRep2D.m_BRep.GetUseVertex(m_IdUV);
	assert( uv.id == m_IdUV );
	m_IdHE = uv.id_he;
}

unsigned int CBRepSurface::CItrVertex::CountEdge() const
{
	if( !is_valid ) return 0;
	assert( m_pBRep2D.m_BRep.IsID_UseVertex(m_IdUV) );
	const CUseVertex& uv = m_pBRep2D.m_BRep.GetUseVertex(m_IdUV);
	const unsigned int id_he0 = uv.id_he;
	unsigned int id_he = id_he0;
	unsigned int icnt = 0;
	for(;;)
	{
		assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he) );
		const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(id_he);		
		if( he.id_e == 0 ){ return 0; } // 浮遊点
		icnt++;
		const unsigned int id_he_b = he.id_he_b;
		assert( m_pBRep2D.m_BRep.IsID_HalfEdge(id_he_b) );
		const CHalfEdge& he_b = m_pBRep2D.m_BRep.GetHalfEdge(id_he_b);
		id_he = he_b.id_he_o;
		if( id_he == id_he0 ) break;
	}
	return icnt;
}

bool CBRepSurface::CItrVertex::IsParent() const
{
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);
	unsigned int id_ul = he.id_ul;
	assert( m_pBRep2D.m_BRep.IsID_UseLoop(id_ul) );
	const CUseLoop& ul = m_pBRep2D.m_BRep.GetUseLoop(id_ul);
	return (ul.id_ul_p == id_ul);
}

bool CBRepSurface::CItrVertex::IsSameUseLoop(const CItrVertex& itrl) const
{
	if( !is_valid ) return false;
	assert( m_pBRep2D.m_BRep.IsID_HalfEdge(m_IdHE) );
	const CHalfEdge& he0 = m_pBRep2D.m_BRep.GetHalfEdge(m_IdHE);

	if( !itrl.is_valid ) return false;
	assert( itrl.m_pBRep2D.m_BRep.IsID_HalfEdge(itrl.m_IdHE) );
	const CHalfEdge& he1 = itrl.m_pBRep2D.m_BRep.GetHalfEdge(itrl.m_IdHE);
	return (he0.id_ul == he1.id_ul);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

bool CBRepSurface::IsElemID(Cad::CAD_ELEM_TYPE type, unsigned int id) const
{
	if(      type == Cad::VERTEX ){
		return m_BRep.IsID_UseVertex(id);
	}
	else if( type == Cad::EDGE ){
		return map_e2he.find(id) != map_e2he.end();
	}
	else if( type == Cad::LOOP ){
		return map_l2ul.find(id) != map_l2ul.end();
	}
	return false;
}

const std::vector<unsigned int> CBRepSurface::GetAryElemID(Cad::CAD_ELEM_TYPE type) const
{
	if( type == Cad::VERTEX ){
		return m_BRep.GetAry_UseVertexID();
	}
	std::vector<unsigned int> res;
	if(      type == Cad::EDGE ){
		std::map<unsigned int,unsigned int>::const_iterator itr = map_e2he.begin();
		for(;itr!=map_e2he.end();itr++){ res.push_back(itr->first); }
	}
	else if( type == Cad::LOOP ){
		std::map<unsigned int,unsigned int>::const_iterator itr = map_l2ul.begin();
		for(;itr!=map_l2ul.end();itr++){ res.push_back(itr->first); }
	}
	return res;
}

bool CBRepSurface::GetIdLoop_Edge(unsigned int id_e, unsigned int& id_l_l, unsigned int& id_l_r) const
{
	unsigned int id_he0;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ) return false;
		id_he0 = itr->second;
	}
	assert( m_BRep.IsID_HalfEdge(id_he0) );
  const Cad::CHalfEdge& he0 = m_BRep.GetHalfEdge(id_he0);
	const unsigned int id_he1 = he0.id_he_o;
  assert( m_BRep.IsID_HalfEdge(id_he1));
  const Cad::CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
  unsigned int id_ul_l, id_ul_r;
  if( he0.is_same_dir ){
    assert( !he1.is_same_dir );
    id_ul_l = he0.id_ul;
    id_ul_r = he1.id_ul;
  }
  else{
    assert(  he1.is_same_dir );
    id_ul_l = he1.id_ul;
    id_ul_r = he0.id_ul;
  }
  const Cad::CUseLoop& ul_l = m_BRep.GetUseLoop(id_ul_l); assert( m_BRep.IsID_UseLoop(id_ul_l) );
  const Cad::CUseLoop& ul_r = m_BRep.GetUseLoop(id_ul_r); assert( m_BRep.IsID_UseLoop(id_ul_r) );
  id_l_l = ul_l.id_l;
  id_l_r = ul_r.id_l;
  return true;
}


unsigned int CBRepSurface::GetIdLoop_Edge(unsigned int id_e, bool is_left) const
{
	unsigned int id_he0;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ) return 0;
		id_he0 = itr->second;
	}
	assert( m_BRep.IsID_HalfEdge(id_he0) );
  const Cad::CHalfEdge& he0 = m_BRep.GetHalfEdge(id_he0);
	const unsigned int id_he1 = he0.id_he_o;
  assert( m_BRep.IsID_HalfEdge(id_he1));
  const Cad::CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
  unsigned int id_ul = 0;
  assert( he0.is_same_dir != he1.is_same_dir );
  if( he0.is_same_dir == is_left ){ id_ul = he0.id_ul; } 
  else{                             id_ul = he1.id_ul; }
  assert( m_BRep.IsID_UseLoop(id_ul) );
  const Cad::CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
  return ul.id_l;
}

CBRepSurface::CItrLoop CBRepSurface::GetItrLoop_SideEdge(unsigned int id_e, bool is_left) const
{
	std::map<unsigned int,unsigned int>::const_iterator itr = this->map_e2he.find(id_e);
	assert( itr != map_e2he.end() );
	const unsigned int id_he1 = itr->second;
	assert( m_BRep.IsID_HalfEdge(id_he1) );
	const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
	if( is_left ){
		unsigned int id_ul1 = he1.id_ul;
		return CItrLoop(*this,id_he1,id_ul1);
	}
	unsigned int id_he2 =  he1.id_he_o;
	assert( m_BRep.IsID_HalfEdge(id_he2) );
	const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
	unsigned int id_ul2 = he2.id_ul;
	return CItrLoop(*this,id_he2,id_ul2);
}


void CBRepSurface::Clear()
{
	this->m_BRep.Clear();
	this->map_e2he.clear();
	this->map_l2ul.clear();
}

bool CBRepSurface::AssertValid() const
{
	// Check Loop
	for(std::map<unsigned int,unsigned int>::const_iterator itr=map_l2ul.begin();itr!=map_l2ul.end();itr++){
    const unsigned int id_l = itr->first;
		unsigned int id_ul_p = itr->second;
		assert( m_BRep.IsID_UseLoop(id_ul_p) );
		{
			const CUseLoop& ul = m_BRep.GetUseLoop(id_ul_p);
			assert( ul.id == id_ul_p ); 
			assert( ul.id_l == id_l );
		}
		assert( this->TypeUseLoop(id_ul_p) == 2 );
		unsigned int id_ul = id_ul_p;
		for(;;){
      assert( m_BRep.IsID_UseLoop(id_ul) );
			const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
			assert( ul.id == id_ul ); 
			assert( ul.id_l == id_l );  // stops here      
      assert( ul.id_ul_p == id_ul_p );
			// goto next loop
			id_ul = ul.id_ul_c;
			if( id_ul == 0 ) break;      
		}
	}

	// Check Edge
	for(std::map<unsigned int,unsigned int>::const_iterator itr=map_e2he.begin();itr!=map_e2he.end();itr++){
//        const unsigned int id_e = itr->first;
		unsigned int id_he = itr->second;
		assert( m_BRep.IsID_HalfEdge(id_he) );
	}

	// Check UseLoop
	const std::vector<unsigned int>& id_ul_ary = this->m_BRep.m_UseLoopSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_ul_ary.size();iid++){
    const unsigned int id_ul = id_ul_ary[iid];
		assert( m_BRep.IsID_UseLoop(id_ul) );
		const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
		assert( ul.id == id_ul );
		const unsigned int id_l = ul.id_l;
		if( id_l == 0 ){ 
      assert( ul.id_ul_p == 0 );
      assert( ul.id_ul_c == 0 );
      continue;
    }
		std::map<unsigned int,unsigned int>::const_iterator itr = map_l2ul.find(id_l);
		assert( itr != map_l2ul.end() );
	}

	// Check UseEdge
	const std::vector<unsigned int>& id_he_ary = this->m_BRep.m_HalfEdgeSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_he_ary.size();iid++){
        const unsigned int id_he = id_he_ary[iid];
		assert( m_BRep.IsID_HalfEdge(id_he) );
		const CHalfEdge& hedge = m_BRep.GetHalfEdge(id_he);
		assert( hedge.id == id_he );

        unsigned int id_v1;
		{
            const unsigned int id_uv1 = hedge.id_uv;
			assert( m_BRep.IsID_UseVertex(id_uv1) );
			const CUseVertex& uv = m_BRep.GetUseVertex(id_uv1);
			assert( uv.id == id_uv1 );
			id_v1 = uv.id_v;
		}

        unsigned int id_v2;
		{
			const unsigned int id_he_f = hedge.id_he_f;
			assert( m_BRep.IsID_HalfEdge(id_he_f) );
			const CHalfEdge& edge_cw = m_BRep.GetHalfEdge(id_he_f);
			assert( edge_cw.id == id_he_f );
			assert( edge_cw.id_he_b == id_he );
			assert( edge_cw.id_ul == hedge.id_ul );
			const unsigned int id_uv2 = edge_cw.id_uv;
			assert( m_BRep.IsID_UseVertex(id_uv2) );
			const CUseVertex& uv = m_BRep.GetUseVertex(id_uv2);
			assert( uv.id == id_uv2 );
			id_v2 = uv.id_v;
		}

		const bool is_same_dir = hedge.is_same_dir;
		const unsigned int id_e = hedge.id_e;
		if( id_e == 0 ){// UVを囲むHEの場合（面に点が一点だけ)
			assert( hedge.id_he_o == id_he );
			assert( hedge.id_he_b == id_he );
			assert( hedge.id_he_f == id_he );
			continue;	
		}
		std::map<unsigned int,unsigned int>::const_iterator itr = map_e2he.find(id_e);
		assert( itr != map_e2he.end() );
//		unsigned int id_he0 = itr->second;
		unsigned int id_vs, id_ve;
		this->GetIdVertex_Edge(id_e,id_vs,id_ve);

		if( is_same_dir ){
			assert( id_v1 == id_vs );
			assert( id_v2 == id_ve );
		}
		else{
			assert( id_v1 == id_ve );
			assert( id_v2 == id_vs );
		}
	}

	// Check UseVertex
	const std::vector<unsigned int>& id_uv_ary = this->m_BRep.m_UseVertexSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_uv_ary.size();iid++){
        const unsigned int id_uv = id_uv_ary[iid];
		assert( m_BRep.IsID_UseVertex(id_uv) );
		const CUseVertex& uv = m_BRep.GetUseVertex(id_uv);
		assert( uv.id == id_uv );
		const unsigned int id_v = uv.id_v;
		assert( id_uv == id_v );
	}
	return true;
}

unsigned int CBRepSurface::AddVertex_Loop(unsigned int id_l)
{
	unsigned int id_ul = 0;
	{
		const std::map<unsigned int, unsigned int>::const_iterator itr =  this->map_l2ul.find(id_l);
		if( itr != this->map_l2ul.end() ){ id_ul = itr->second; }
	}
	/* UseVertex追加 */
	unsigned int id_uv_add,id_he_add,id_ul_add;
	m_BRep.MVEL(id_uv_add,id_he_add,id_ul_add,id_ul);
	this->m_BRep.AssertValid_Use();
	m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l);
	/* Vertex追加 */
	unsigned int id_v_add = id_uv_add;
	m_BRep.SetVertexIDtoUseVertex(id_uv_add,id_v_add);
	assert( this->AssertValid() );
	return id_v_add;
}

bool CBRepSurface::GetIdVertex_Edge(unsigned int id_e, unsigned int& id_v1, unsigned int& id_v2) const
{
	unsigned int id_he;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = map_e2he.find(id_e);
        if( itr==map_e2he.end() ) return false;
		id_he = itr->second;
	}
    assert( this->m_BRep.IsID_HalfEdge(id_he) );
	const Cad::CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
	id_v1 = he.id_uv;
	unsigned int id_he_f = he.id_he_f;
    assert( this->m_BRep.IsID_HalfEdge(id_he_f) );
	const Cad::CHalfEdge& he_f = m_BRep.GetHalfEdge(id_he_f);
	id_v2 = he_f.id_uv;
	return true;
}

unsigned int CBRepSurface::GetIdVertex_Edge(unsigned int id_e, bool is_root ) const
{
	unsigned int id_he;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = map_e2he.find(id_e);
		id_he = itr->second;
	}
	if( !this->m_BRep.IsID_HalfEdge(id_he) ) return 0;
	const Cad::CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
	const unsigned int id_v1 = he.id_uv;
	unsigned int id_he_f = he.id_he_f;
	if( !this->m_BRep.IsID_HalfEdge(id_he_f) ) return 0;
	const Cad::CHalfEdge& he_f = m_BRep.GetHalfEdge(id_he_f);
	const unsigned int id_v2 = he_f.id_uv;
	if( is_root ) return id_v1;
	return id_v2;
}

unsigned int CBRepSurface::GetFreeKey(const std::map<unsigned int,unsigned int>& map)
{
	if( map.empty() ) return 1;
	////////////////
	{
		std::map<unsigned int,unsigned int>::const_reverse_iterator ritr = map.rbegin();
		const unsigned int max_id = ritr->first;
		if( max_id < map.size() * 2 ){ return max_id+1; }
	}
	////////////////
	unsigned int cand = 1;
	std::map<unsigned int,unsigned int>::const_iterator itr = map.begin();
	for(;itr!=map.end();itr++){
		const unsigned int key = itr->first;
		if( key > cand ) return cand;
		cand++;
	}
	return cand;
}

unsigned int CBRepSurface::AddVertex_Edge(unsigned int id_e)
{
	unsigned int id_he;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ){ return 0; }
		id_he = itr->second;
	}
	unsigned int id_v1,id_v2;
	if( !this->GetIdVertex_Edge(id_e,id_v1,id_v2) ) return 0;

	////////////////////////////////
	// Leave Input Check Section

	unsigned int id_uv_add;
//	const std::vector<unsigned int> id_ary_he = this->m_BRep.FindHalfEdge_Edge(id_e);
//	assert( id_ary_he.size() == 2 );
//	const unsigned int id_he = id_ary_he[0];
	const CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
	assert( he.id_e == id_e );
	////////////////
	const unsigned int id_he_o = he.id_he_o;
//	assert( id_he_o == id_ary_he[1] );
	assert( m_BRep.IsID_HalfEdge(id_he_o) );
	unsigned int id_he1, id_he2;
	if( he.is_same_dir ){	id_he1 = id_he;		id_he2 = id_he_o; }
	else{					id_he1 = id_he_o; 	id_he2 = id_he;	  }
	unsigned int id_he_add1, id_he_add2;
	////////////////
	if( !this->m_BRep.MVE(id_he_add1,id_he_add2,id_uv_add,id_he1) ) assert(0);
	this->m_BRep.AssertValid_Use();
	const unsigned int id_e_add = GetFreeKey(this->map_e2he);
	m_BRep.SetEdgeIDtoHalfEdge(id_he1,    id_e,    true);
	m_BRep.SetEdgeIDtoHalfEdge(id_he_add1,id_e_add,true);
	m_BRep.SetEdgeIDtoHalfEdge(id_he2,    id_e_add,false);
	m_BRep.SetEdgeIDtoHalfEdge(id_he_add2,id_e,    false);
	unsigned int id_v_add = id_uv_add;
	m_BRep.SetVertexIDtoUseVertex(id_uv_add,id_v_add);
	map_e2he.insert( std::make_pair(id_e_add,id_he_add1) );
	assert( this->AssertValid() );
	return id_v_add;
}

std::vector< std::pair< unsigned int, bool > > CBRepSurface::GetItrLoop_RemoveEdge(unsigned int id_e) const
{
	std::vector< std::pair<unsigned int,bool> > aIdEDir;
	unsigned int id_he1;
	{
		std::map<unsigned int,unsigned int>::const_iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ){ return aIdEDir; }
		id_he1 = itr->second;
	}
	assert( m_BRep.IsID_HalfEdge(id_he1) );
	const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
	unsigned int id_he2 = he1.id_he_o;
	if( id_he2 == id_he1 ){ return aIdEDir; }	// 浮遊点
	assert( m_BRep.IsID_HalfEdge(id_he2) );
	const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
	if( he1.id_ul != he2.id_ul ){ return aIdEDir; }	// 領域を２分する
	////////////////
	unsigned int id_he = he2.id_he_f;
	for(;;){
		assert( m_BRep.IsID_HalfEdge(id_he) );
		const CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
		aIdEDir.push_back( std::make_pair(he.id_e,he.is_same_dir) );
		const unsigned int id_he_f =he.id_he_f;
		assert( id_he_f != id_he2 );
		if( id_he_f == id_he1 ) break;
		id_he = id_he_f;
	}
	return aIdEDir;
}

std::vector< std::pair<unsigned int,bool> > CBRepSurface::GetItrLoop_ConnectVertex(
		const CItrVertex& itrv1, const CItrVertex& itrv2) const
{
	std::vector< std::pair<unsigned int,bool> > aIdEDir;

	const unsigned int id_uv1 = itrv1.GetIdUseVertex();
	const unsigned int id_uv2 = itrv2.GetIdUseVertex();
	if( !m_BRep.IsID_UseVertex(id_uv1) ){ return aIdEDir; }
	if( !m_BRep.IsID_UseVertex(id_uv2) ){ return aIdEDir; }
	if( id_uv1 == id_uv2 ){ return aIdEDir; }
	if( itrv1.GetIdLoop() != itrv2.GetIdLoop() ){ return aIdEDir; }	// 違うループに属する点を結ぼうとしている
	////////////////
	unsigned int id_he1 = itrv1.GetIdHalfEdge();
	unsigned int id_ul1;
	{
		assert( m_BRep.IsID_HalfEdge(id_he1) );
		const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
		id_ul1 = he1.id_ul; 
		assert( he1.id_e != 0 );
	}
	////////////////
	unsigned int id_he2 = itrv2.GetIdHalfEdge();
	unsigned int  id_ul2;
	{
		assert( m_BRep.IsID_HalfEdge(id_he2) );
		const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
		id_ul2 = he2.id_ul;	
		assert( he2.id_e != 0 );
	}
	////////////////
	if( id_ul1 != id_ul2 ) return aIdEDir;	// UseLoopが違う
	unsigned int id_he = id_he2;
	for(;;){
		const CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
		aIdEDir.push_back( std::make_pair(he.id_e,he.is_same_dir) );
		id_he = he.id_he_f;
		if( id_he == id_he1 ){ break; }
	}
	return aIdEDir;
}

CBRepSurface::CResConnectVertex CBRepSurface::ConnectVertex
(const CItrVertex& itrv1, const CItrVertex& itrv2, bool is_left_ladd)
{
  CBRepSurface::CResConnectVertex res;
//  std::cout << "CBRep2D::ConnectVertex" << std::endl;
	const unsigned int id_uv1 = itrv1.GetIdUseVertex();
	const unsigned int id_uv2 = itrv2.GetIdUseVertex();
	if( !m_BRep.IsID_UseVertex(id_uv1) ){ return res; }
	if( !m_BRep.IsID_UseVertex(id_uv2) ){ return res; }
	if( id_uv1 == id_uv2 ){ return res; }
  res.id_v1 = m_BRep.GetUseVertex(id_uv1).id_v;
  res.id_v2 = m_BRep.GetUseVertex(id_uv2).id_v;  

	if( itrv1.GetIdLoop() != itrv2.GetIdLoop() ){ return res; }
	unsigned int id_l = itrv1.GetIdLoop();
  res.id_l = id_l;

	unsigned int id_he1 = itrv1.GetIdHalfEdge();
	unsigned int id_he2 = itrv2.GetIdHalfEdge();

	unsigned int id_ul1; bool is_float1;
	{
		assert( m_BRep.IsID_HalfEdge(id_he1) );
		const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
		id_ul1 = he1.id_ul; assert( m_BRep.IsID_UseLoop(id_ul1) );
		is_float1 = (he1.id_e==0);
	}
	////////////////
	unsigned int id_ul2; bool is_float2;
	{
		assert( m_BRep.IsID_HalfEdge(id_he2) );
		const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
		id_ul2 = he2.id_ul;	assert( m_BRep.IsID_UseLoop(id_ul2) );
		is_float2 = (he2.id_e==0);
	}

	// Leave Input Check Section
	////////////////////////////////

	if( id_ul1 != id_ul2 ){	// class of MEKL, connect [parent and child loop], [two child loop] or [two boundray loop]
    if( is_float1 && !is_float2 ){	// id_uv1 is floating vertex
      //				std::cout << "Add Edge In Loop One Floating Vertex(1)" << std::endl;
      unsigned int id_he_add1;
      if( !m_BRep.MEKL_OneFloatingVertex(id_he_add1, id_he2,id_he1) ){ assert(0); }
      const unsigned int id_e_add = this->GetFreeKey(map_e2he);
      m_BRep.SetEdgeIDtoHalfEdge(id_he_add1,id_e_add,false);
      m_BRep.SetEdgeIDtoHalfEdge(id_he1,    id_e_add,true);
      map_e2he.insert( std::make_pair(id_e_add,id_he1) );
      assert( this->AssertValid() );
      res.id_e_add = id_e_add;
      return res;
    }
    ////////////////
    else if( is_float2 && !is_float1 ){	// id_uv2 is floating vertex
      //				std::cout << "Add Edge In Loop One Floating Vertex(2)" << std::endl;
      unsigned int id_he_add1;
      if( !m_BRep.MEKL_OneFloatingVertex(id_he_add1, id_he1,id_he2) ){ assert(0); }
      const unsigned int id_e_add = this->GetFreeKey(map_e2he);
      m_BRep.SetEdgeIDtoHalfEdge(id_he_add1,id_e_add,true );
      m_BRep.SetEdgeIDtoHalfEdge(id_he2,    id_e_add,false);
      this->map_e2he.insert( std::make_pair(id_e_add,id_he_add1) );
      assert( this->AssertValid() );
      res.id_e_add = id_e_add;
      return res;
    }
    ////////////////
    else if( is_float1 && is_float2 ){	// id_uv1, id_uv2 are floating vertex
      //				std::cout << "Add Edge In Loop Two Floating Vertex" << std::endl;
      if( !m_BRep.MEKL_TwoFloatingVertex(id_he1,id_he2) ){ assert(0); }
      const unsigned int id_e_add = this->GetFreeKey(map_e2he);
      m_BRep.SetEdgeIDtoHalfEdge(id_he1,id_e_add,true );
      m_BRep.SetEdgeIDtoHalfEdge(id_he2,id_e_add,false);
      this->map_e2he.insert( std::make_pair(id_e_add,id_he1) );				
      assert( this->AssertValid() );
      res.id_e_add = id_e_add;
      return res;
    }
    ////////////////
    else{   // connect [parent and its child loops], [two child loop] or [two boundary loop]
      //				std::cout << "Add Edge with Child Loop " << std::endl;
      unsigned int id_he_add1, id_he_add2;
      if( !m_BRep.MEKL(id_he_add1, id_he_add2, id_he1,id_he2) ){ assert(0); }
      const unsigned int id_e_add = this->GetFreeKey(map_e2he);
      m_BRep.SetEdgeIDtoHalfEdge(id_he_add1,id_e_add,true );
      m_BRep.SetEdgeIDtoHalfEdge(id_he_add2,id_e_add,false);
      this->map_e2he.insert( std::make_pair(id_e_add,id_he_add1) );
      if( id_l != 0 ){ // connect [parent and its child loop] or [two child loop]
        unsigned int id_ul1p = m_BRep.GetUseLoop(id_ul1).id_ul_p; assert(id_ul1p!=0);
        std::map<unsigned int,unsigned int>::iterator itr = this->map_l2ul.find(id_l); assert( itr != this->map_l2ul.end() );
        itr->second = id_ul1p;
      }
      else{} // connect boundary loops
      assert( this->AssertValid() );
      res.id_e_add = id_e_add;
      return res;
    }
	}
	else if( id_ul1 == id_ul2 ){	// class of MEL, split loop into two
		// 辺v1-v2の左側が古いループ、右側が新しいループとなる
		unsigned int id_he_add1,id_he_add2,id_ul_add;
    if( !m_BRep.MEL(id_he_add1,id_he_add2,id_ul_add,  id_he1,id_he2) ){ assert(0); }
		assert( m_BRep.AssertValid_Use()==0 );
		const unsigned int id_e_add = this->GetFreeKey(map_e2he);
		m_BRep.SetEdgeIDtoHalfEdge(id_he_add1,id_e_add,true);
		m_BRep.SetEdgeIDtoHalfEdge(id_he_add2,id_e_add,false);
		this->map_e2he.insert( std::make_pair(id_e_add,id_he_add1) );
		////////////////
		bool iflag = true;	// trueなら新しく加わるループはid_ul_addを親ループとする,falseならid_ul1が親
		if( id_l == 0 ){	// ループの外側にループを追加する場合
			iflag = !is_left_ladd;
		}
		else{
			unsigned int id_ul_p = map_l2ul.find(id_l)->second;
//      std::cout << "split l1 ladd lp : " << id_ul1 << " " << id_ul_add << " " << id_ul_p << std::endl;
			if( id_ul1 != id_ul_p ){	// 子ループの外側にループを追加する場合
				if( is_left_ladd ){	// id_ul_addが子ループの外側のループになった場合
					m_BRep.SwapUseLoop(id_ul_add,id_ul1);
					m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l);  
					m_BRep.AssertValid_Use();
					iflag = false;
				}				
				else{}	// id_ul1が子ループの外側になった場合
			}
			else{}	// ループを二分する場合
		}
		////////////////
		const unsigned int id_l_add = this->GetFreeKey(map_l2ul);
		if( iflag ){	// 新しいループはid_ul_addのループ
			this->map_l2ul.insert( std::make_pair(id_l_add,id_ul_add) );
			m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l_add);
		}
		else{	// 新しいループはid_ul1のループ
			this->map_l2ul.insert( std::make_pair(id_l_add,id_ul1) );
			m_BRep.SetLoopIDtoUseLoop(id_ul1   ,id_l_add);
		}
		assert( this->AssertValid() );
    res.id_e_add = id_e_add;
    res.id_l_add = id_l_add;
    res.is_left_l_add = !iflag;
    return res;
	}
	assert(0);
	return res;
}

bool CBRepSurface::SwapItrLoop(const CItrLoop& itrl, unsigned int id_l_to )
{
	unsigned int id_ul_p_to;
	{	// 新しくできるループの親半ループ
		id_ul_p_to = map_l2ul.find(id_l_to)->second;
		const CUseLoop& ul = m_BRep.GetUseLoop(id_ul_p_to);
		assert( ul.id_ul_p == id_ul_p_to );
	}
	const unsigned int id_ul_fr = itrl.GetIdUseLoop();
	{
		assert( this->m_BRep.IsID_UseLoop(id_ul_fr) );
		const CUseLoop& ul = m_BRep.GetUseLoop(id_ul_fr);
		assert( ul.id_ul_p != id_ul_fr );
	}
	m_BRep.MoveUseLoop(id_ul_fr,id_ul_p_to);
	m_BRep.SetLoopIDtoUseLoop(id_ul_fr,id_l_to);
	return true;
}

bool CBRepSurface::RemoveEdge(unsigned int id_e, bool is_del_cp)
{
	unsigned int id_he1 = 0;
	{
		std::map<unsigned int,unsigned int>::iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ) return false;
		id_he1 = itr->second;
	}
	unsigned int id_he2;
	{
		assert( m_BRep.IsID_HalfEdge(id_he1) );
		const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
		id_he2 = he1.id_he_o;
	}
	////////////////////////////////
	unsigned int id_ul1;
	unsigned int id_l1;
	bool is_edge_vertex2=false;
	{
		const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
		assert( he1.id_he_o == id_he2 );
		id_ul1 = he1.id_ul;
		const CUseLoop& ul1 = m_BRep.GetUseLoop(id_ul1);
		id_l1 = ul1.id_l;
		if( he1.id_he_f == id_he2 ){
//			assert( he2.id_he_b == id_he1 );
			is_edge_vertex2 = true;
		}
	}
	////////////////////////////////
	unsigned int id_ul2;
	unsigned int id_l2;
	bool is_edge_vertex1=false;
	{
		const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
		assert( he2.id_he_o == id_he1 );
		id_ul2 = he2.id_ul;
		const CUseLoop& ul2 = m_BRep.GetUseLoop(id_ul2);
		id_l2 = ul2.id_l;
		if( he2.id_he_f == id_he1 ){
//			assert( he1.id_he_b == id_he2 );
			is_edge_vertex1 = true;
		}
	}
	////////////////////////////////
	if( id_l1 != id_l2 ){	// ２つのループを隔てていた辺を消す
//		std::cout << "Remove Splited Edge" << std::endl;
		assert( id_ul1 != id_ul2 );
		assert( !is_edge_vertex1 && !is_edge_vertex2 );
    unsigned int id_lp_new;	// id_ul1の新しい幾何親ループ              
    if( id_l1 != 0 ){ 
      if( id_l2 == 0 ){ id_lp_new = 0; }
      else{             id_lp_new = id_l1; }
    }
    else{               id_lp_new = 0; }
    ////////
    {
      unsigned int id_ul = id_ul1;
      const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);         
      id_ul = ul.id_ul_p;
      for(;;){	// 全ての子ループに新しい幾何親ループを設定
        if( id_ul == 0 ) break;          
        m_BRep.SetLoopIDtoUseLoop(id_ul,id_lp_new);
        assert( m_BRep.IsID_UseLoop(id_ul) );
        const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);        
        id_ul = ul.id_ul_c;
      }
    }
    {
      unsigned int id_ul = id_ul2;
      const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);         
      id_ul = ul.id_ul_p;
      for(;;){	// 全ての子ループに新しい幾何親ループを設定
        if( id_ul == 0 ) break;          
        m_BRep.SetLoopIDtoUseLoop(id_ul,id_lp_new);
        assert( m_BRep.IsID_UseLoop(id_ul) );
        const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);        
        id_ul = ul.id_ul_c;
      }
    }  
    /////
    if( id_l1 != 0 ){
      if( !m_BRep.KEL(id_he1) ) assert(0);      
      if( id_l2 == 0 ){ map_l2ul.erase(id_l1); }
      else{  // id_l1!=0 && id_l2!=0 
        map_l2ul.erase(id_l1);
        map_l2ul.erase(id_l2); 
        // regist the parent use loop for id_l1
        assert( m_BRep.IsID_UseLoop(id_ul1) );
        const CUseLoop& ul = m_BRep.GetUseLoop(id_ul1);
        assert( m_BRep.IsID_UseLoop(ul.id_ul_p) );
        map_l2ul.insert( std::make_pair(id_l1,ul.id_ul_p) );
      }
    }
    else{ // id_l1 == 0 
      assert( id_l2 != 0 );
      if( !m_BRep.KEL(id_he2) ) assert(0);              
      map_l2ul.erase(id_l2);
    }    
		map_e2he.erase(id_e);
		assert( this->AssertValid() );
		return true;
	}
	if( id_ul1 == id_ul2 ){	// MEKL系、面を分割していない辺を消す
		assert( id_l1 == id_l2 );
		if( !is_edge_vertex1 && !is_edge_vertex2 ){
			if( id_l1 != 0 ){	// 両方どちらも端点でなく、ループの内側にある辺
				unsigned int id_ul_add;
				if( !m_BRep.KEML(id_ul_add,id_he1) ) assert(0);
				if( is_del_cp ){
					assert( TypeUseLoop(id_ul_add)==2 );
					if( !m_BRep.SwapUseLoop(id_ul_add,id_ul1) ) assert(0); 
					this->map_l2ul.find(id_l1)->second = id_ul_add;
				}
				m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l1);
			}
			else{	// 二つのループを外側でつないでる辺
        unsigned int id_ul_add;
				if( !m_BRep.KEML(id_ul_add,id_he1) ) assert(0);        
				m_BRep.SetLoopIDtoUseLoop(id_ul_add,0);        
			}
		}
		else if( is_edge_vertex1 && is_edge_vertex2 ){	// 浮遊辺を消す
//			std::cout << "Reverse MEKL_TwoFloatingVertex" << std::endl;
			unsigned int id_ul_add;
			if( !m_BRep.KEML_TwoFloatingVertex(id_ul_add,id_he1) ) assert(0);
			m_BRep.SetEdgeIDtoHalfEdge(id_he1,0,true);
			m_BRep.SetEdgeIDtoHalfEdge(id_he2,0,true);
			m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l1);
		}
		else if( is_edge_vertex1 && !is_edge_vertex2 ){	// uv1が端点な辺を消す
//			std::cout << "Reverse MEKL_OneFloatingVertex(1)" << std::endl;
			unsigned int id_ul_add;
			if( !m_BRep.KEML_OneFloatingVertex(id_ul_add,id_he1) ) assert(0);
			m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l1);
		}
		else if( is_edge_vertex2 && !is_edge_vertex1){	// uv2が端点な辺を消す
//			std::cout << "Reverse MEKL_OneFloatingVertex(2)" << std::endl;
			unsigned int id_ul_add;
			if( !m_BRep.KEML_OneFloatingVertex(id_ul_add,id_he2) ) assert(0);
			m_BRep.SetLoopIDtoUseLoop(id_ul_add,id_l1);
		}
		map_e2he.erase(id_e);
		assert( this->AssertValid() );
		return true;
	}
	return false;
}

unsigned int CBRepSurface::TypeUseLoop(unsigned int id_ul) const
{
	assert( m_BRep.IsID_UseLoop(id_ul) );
	unsigned int id_he_ini;
	{	// 浮遊点かどうかチェックする
		const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
		id_he_ini = ul.id_he;
		assert( m_BRep.IsID_HalfEdge(id_he_ini) );
		const CHalfEdge& he_ini = m_BRep.GetHalfEdge(id_he_ini);
		if( he_ini.id_he_f == id_he_ini ){
			assert( he_ini.id_he_b == id_he_ini );
			return 0;
		}
	}
	unsigned int id_he = id_he_ini;
	for(;;){	// 全ての半辺を回ってみて、一つでも浮遊辺じゃないところがあったら面
		assert( m_BRep.IsID_HalfEdge(id_he) );
		const CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
		const unsigned int id_he_f = he.id_he_f;
		{
			const unsigned int id_he_o = he.id_he_o;
			assert( m_BRep.IsID_HalfEdge(id_he_o) );
			const CHalfEdge& he_o = m_BRep.GetHalfEdge(id_he_o);
			if( he_o.id_ul != id_ul ) return 2;
		}
		if( id_he_f == id_he_ini ) break;
		id_he = id_he_f;
	}
	return 1;
}


bool CBRepSurface::RemoveVertex(unsigned int id_v)
{
	unsigned int id_uv = id_v;
	if( !this->m_BRep.IsID_UseVertex(id_uv) ) return false;

	unsigned int nedge_around_vtx = 0;
	{
		CItrVertex itrv(*this,id_v);
		nedge_around_vtx = itrv.CountEdge();
	}

	if(      nedge_around_vtx == 0 ){	// 面上の浮遊点を消す
    assert( m_BRep.IsID_UseVertex(id_uv) );
		if( !m_BRep.KVEL(id_uv) ) assert(0);
	}
	else if( nedge_around_vtx == 2 ){	// 辺上の点を消す
		const CUseVertex& uv = this->m_BRep.GetUseVertex(id_uv);
		const unsigned int id_he1 = uv.id_he;
		const CHalfEdge& he1 = this->m_BRep.GetHalfEdge(id_he1);
		const unsigned int id_e1 = he1.id_e;
		const unsigned int id_he2 = he1.id_he_o;
		const CHalfEdge& he2 = this->m_BRep.GetHalfEdge(id_he2);
		{	// ループ上の点が２だったら省く
			unsigned int id_ul1 = he1.id_ul;
			unsigned int id_ul2 = he2.id_ul;
			const unsigned int nvl1 = CBRepSurface::CItrLoop(*this,id_he1,id_ul1).CountVertex_UseLoop();
			const unsigned int nvl2 = CBRepSurface::CItrLoop(*this,id_he2,id_ul2).CountVertex_UseLoop();
			assert( nvl1 > 1 && nvl2 > 1);
			if( nvl1 == 2 || nvl2 == 2 ) return false;
		}
		const unsigned int id_uv2 = he2.id_uv;  
		assert( m_BRep.IsID_UseVertex(id_uv2) );
		unsigned int id_he_remove1 = id_he1;
//		unsigned int id_he_replace1 = he1.id_he_b;
    if( !m_BRep.KVE(id_he_remove1) ) assert(0);
		assert( m_BRep.AssertValid_Use()==0 );
		map_e2he.erase(id_e1);
  }
	assert( this->AssertValid() );
	return true;
}

bool CBRepSurface::MakeHole_fromLoop(unsigned int id_l)
{
	unsigned int id_ul1 = 0;
	{
		std::map<unsigned int,unsigned int>::iterator itr = this->map_l2ul.find(id_l);
		if( itr == this->map_l2ul.end() ) return false;
		id_ul1 = itr->second;
	}
	assert( this->m_BRep.IsID_UseLoop(id_ul1) );
	m_BRep.SetLoopIDtoUseLoop(id_ul1,0);
	map_l2ul.erase(id_l);
	return true;
}

unsigned int CBRepSurface::SealHole(unsigned int id_e, bool is_left)
{  
  unsigned int id_ul = 0;
	{
    unsigned int id_he1 = 0;    
		std::map<unsigned int,unsigned int>::iterator itr = this->map_e2he.find(id_e);
		if( itr == this->map_e2he.end() ) return 0;
		id_he1 = itr->second;
		assert( m_BRep.IsID_HalfEdge(id_he1) );
		const CHalfEdge& he1 = m_BRep.GetHalfEdge(id_he1);
    if( he1.is_same_dir == is_left ){ id_ul = he1.id_ul; }
    else{
      const unsigned int id_he2 = he1.id_he_o;
      if( id_he1 == id_he2 ) return 0;  // uv is floating vertex
      assert( m_BRep.IsID_HalfEdge(id_he2) );
      const CHalfEdge& he2 = m_BRep.GetHalfEdge(id_he2);
      id_ul = he2.id_ul;
    }
	}
  assert( m_BRep.IsID_UseLoop(id_ul) );
  {
    const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
    if( ul.id_l != 0 ){ return 0; }
  }
  assert( this->TypeUseLoop(id_ul) == 2 );
  unsigned int id_l_new = GetFreeKey(map_l2ul);
  m_BRep.SetLoopIDtoUseLoop(id_ul,id_l_new);
	map_l2ul.insert( std::make_pair(id_l_new,id_ul) );
  return id_l_new;
}


bool CBRepSurface::Serialize( Com::CSerializer& arch )
{
	if( arch.IsLoading() ){	// 読み込み時の処理
		this->Clear();
		const unsigned int buff_size = 256;
		char class_name[buff_size];
		arch.ReadDepthClassName(class_name,buff_size);
		if( strncmp(class_name,"BRep2D",6)!=0 ) return true;
		{
			int ne;
			arch.Get("%d",&ne);
			assert(ne>=0);
      for(unsigned int ie=0;ie<(unsigned int)ne;ie++){
				int id_e, id_he, itmp;
				arch.Get("%d %d %d", &itmp,&id_e,&id_he);
				map_e2he.insert( std::make_pair(id_e,id_he) );
			}
		}
		{
			int nl;
			arch.Get("%d",&nl);
			assert(nl>=0);
      for(unsigned int il=0;il<(unsigned int)nl;il++){
				int id_l, id_ul, itmp;
				arch.Get("%d %d %d", &itmp,&id_l,&id_ul);
				map_l2ul.insert( std::make_pair(id_l,id_ul) );
			}
		}
		////////////////
		int nuv, nhe, nul;
		{
			arch.Get("%d%d%d",&nuv, &nhe, &nul);
			assert( nuv>0 ); assert( nhe>0 ); assert( nul>0 );
			this->m_BRep.m_UseVertexSet.Reserve(nuv*2);
			this->m_BRep.m_HalfEdgeSet.Reserve(nhe*2);
			this->m_BRep.m_UseLoopSet.Reserve(nul*2);
		}
		////////////////////////////////////////////////
		arch.ShiftDepth(true);
    for(int iuv=0;iuv<nuv;iuv++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CUseVertex",10)  ==0 );
			int id;		arch.Get("%d",&id);		assert( id>0 );
			int id_v;	arch.Get("%d",&id_v);	assert( id_v>0 );
			int id_he;	arch.Get("%d",&id_he);	assert( id_he>0 );
			int tmp_id = this->m_BRep.m_UseVertexSet.AddObj( CUseVertex(id,id_he) );
			assert( tmp_id == id );
			m_BRep.SetVertexIDtoUseVertex(id,id_v);
		}
    for(int ihe=0;ihe<nhe;ihe++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CHalfEdge",9)==0 );
			int id;				arch.Get("%d",&id);				assert( id>0 );
      int id_e;			arch.Get("%d",&id_e);			assert( id_e>=0 );
			int i_is_same_dir;	arch.Get("%d",&i_is_same_dir);	assert( i_is_same_dir>=0 );
			int id_uv;			arch.Get("%d",&id_uv);			assert( id_uv>0 );
			int id_he_f, id_he_ccw, id_he_o;
			arch.Get("%d%d%d",&id_he_f,&id_he_ccw,&id_he_o);
			assert( id_he_f>0 && id_he_ccw>0 && id_he_o>0 );
			int id_ul;			arch.Get("%d",&id_ul);			assert( id_ul>0 );
			const bool is_same_dir = (i_is_same_dir!=0);
			const int tmp_id = this->m_BRep.m_HalfEdgeSet.AddObj( CHalfEdge(id, id_uv,  id_he_f,id_he_ccw,id_he_o,  id_ul) );
			assert( tmp_id == id );
			m_BRep.SetEdgeIDtoHalfEdge(id,id_e,is_same_dir);
		}
    for(int iul=0;iul<nul;iul++){
			arch.ReadDepthClassName(class_name,buff_size);
			assert( strncmp(class_name,"CUseLoop",8)==0 );
			int id;			arch.Get("%d",&id);			assert( id>0 );
			int id_l;		arch.Get("%d",&id_l);		assert( id_l>=0 );
			int id_he;		arch.Get("%d",&id_he);		assert( id_he>0 );
			int id_ul_c;	arch.Get("%d",&id_ul_c);	assert( id_ul_c>=0 );
			int id_ul_p;	arch.Get("%d",&id_ul_p);	assert( id_ul_p>=0 );
			const int tmp_id = this->m_BRep.m_UseLoopSet.AddObj( CUseLoop(id,id_he,id_ul_c,id_ul_p) );
			assert( tmp_id == id );
			m_BRep.SetLoopIDtoUseLoop(id,id_l);
		}
		arch.ShiftDepth(false);
		this->AssertValid();
		return true;
	}
	else{ // 書き込み時の処理
    // クラスの名前の指定，サイズの指定
		arch.WriteDepthClassName("BRep2D");
		{
			arch.Out("%d\n",map_e2he.size());
			std::map<unsigned int,unsigned int>::const_iterator itr = map_e2he.begin();
			unsigned icnt = 0;
			for(;itr!=map_e2he.end();itr++){
				arch.Out("%d %d %d\n",icnt,itr->first,itr->second);
				icnt++;
			}
		}
		{
			arch.Out("%d\n",map_l2ul.size());
			std::map<unsigned int,unsigned int>::const_iterator itr = map_l2ul.begin();
			unsigned icnt = 0;
			for(;itr!=map_l2ul.end();itr++){
				arch.Out("%d %d %d\n",icnt,itr->first,itr->second);
				icnt++;
			}
		}
		arch.Out("%d %d %d\n",this->m_BRep.m_UseVertexSet.GetAry_ObjID().size(), this->m_BRep.m_HalfEdgeSet.GetAry_ObjID().size(),this->m_BRep.m_UseLoopSet.GetAry_ObjID().size());
		arch.ShiftDepth(true);
		////////////////
        // UseVertexの出力
        {
			const std::vector<unsigned int> id_ary = this->m_BRep.m_UseVertexSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_uv = id_ary[iid];
				assert( m_BRep.IsID_UseVertex(id_uv) );
				const CUseVertex& uv = m_BRep.GetUseVertex(id_uv);
				assert( uv.id == id_uv );
				arch.WriteDepthClassName("CUseVertex");
				arch.Out("%d\n",id_uv);
				arch.Out("%d\n",uv.id_v);
				arch.Out("%d\n",uv.id_he);
			}
		}
        // HalfEdgeの出力
        {
			const std::vector<unsigned int> id_ary = this->m_BRep.m_HalfEdgeSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_he = id_ary[iid];
				assert( m_BRep.IsID_HalfEdge(id_he) );
				const CHalfEdge& he = m_BRep.GetHalfEdge(id_he);
				assert( he.id == id_he );
				arch.WriteDepthClassName("CHalfEdge");
				arch.Out("%d\n",id_he);
				arch.Out("%d\n",he.id_e);
				arch.Out("%d\n",(int)he.is_same_dir);
				arch.Out("%d\n",he.id_uv);
				arch.Out("%d %d %d\n",he.id_he_f, he.id_he_b, he.id_he_o);
				arch.Out("%d\n",he.id_ul);
			}
		}
    // UseLoopの出力
    {
			const std::vector<unsigned int> id_ary = this->m_BRep.m_UseLoopSet.GetAry_ObjID();
			for(unsigned int iid=0;iid<id_ary.size();iid++){
				const unsigned int id_ul = id_ary[iid];
				assert( m_BRep.IsID_UseLoop(id_ul) );
				const CUseLoop& ul = m_BRep.GetUseLoop(id_ul);
				assert( ul.id == id_ul );
				arch.WriteDepthClassName("CUseLoop");
				arch.Out("%d\n",id_ul);
				arch.Out("%d\n",ul.id_l);
				arch.Out("%d\n",ul.id_he);
				arch.Out("%d\n",ul.id_ul_c);
				arch.Out("%d\n",ul.id_ul_p);
			}
		}
		arch.ShiftDepth(false);
	}
	return true;
}
