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

#define for if(0);else for

#ifdef __VISUALC__
	#pragma warning(disable:4786)
	#pragma warning(disable:4996)
#endif

#include <set>
#include <iostream>
#include <assert.h>

#include "delfem/cad/brep.h"

using namespace Cad;

void CBRep::Clear()
{
	this->m_HalfEdgeSet.clear();
	this->m_UseLoopSet.clear();
	this->m_UseVertexSet.clear();
}

const CUseLoop& CBRep::GetUseLoop(unsigned int id_ul) const 
{ 
	assert( m_UseLoopSet.IsObjID(id_ul) );
	if( !m_UseLoopSet.IsObjID(id_ul) ) throw;
	return m_UseLoopSet.GetObj(id_ul); 
}


const CUseVertex& CBRep::GetUseVertex(unsigned int id_uv) const
{ 
	assert( m_UseVertexSet.IsObjID(id_uv) );
	if( !m_UseVertexSet.IsObjID(id_uv) ) throw;
	return m_UseVertexSet.GetObj(id_uv); 
}

const CHalfEdge& CBRep::GetHalfEdge(unsigned int id_he) const
{ 
	assert( m_HalfEdgeSet.IsObjID(id_he) );
	if( !m_HalfEdgeSet.IsObjID(id_he) ) throw;
	return m_HalfEdgeSet.GetObj(id_he); 
}


bool CBRep::SetLoopIDtoUseLoop(unsigned int id_ul, unsigned int id_l)
{
	if( !m_UseLoopSet.IsObjID(id_ul) ) return false;
	CUseLoop& ul = this->m_UseLoopSet.GetObj(id_ul);
	ul.id_l = id_l;
  if( ul.id_ul_p == 0 || ul.id_ul_p == id_ul ){
    if( id_l == 0 ){ ul.id_ul_p == 0; }
    else{            ul.id_ul_p = id_ul; }
  }
	return true;
}

bool CBRep::SetVertexIDtoUseVertex(unsigned int id_uv, unsigned int id_v)
{
	if( !m_UseVertexSet.IsObjID(id_uv) ) return false;
	CUseVertex& uv = this->m_UseVertexSet.GetObj(id_uv);
	uv.id_v = id_v;
	return true;
}

bool CBRep::SetEdgeIDtoHalfEdge(unsigned int id_he, unsigned int id_e, bool is_same_dir)
{
	if( !m_HalfEdgeSet.IsObjID(id_he) ) return false;
	CHalfEdge& he = this->m_HalfEdgeSet.GetObj(id_he);
	he.id_e = id_e;
	he.is_same_dir = is_same_dir;
	return true;
}

int CBRep::AssertValid_Use() const 
{
	// Check UseVertex
	const std::vector<unsigned int>& id_uv_ary = m_UseVertexSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_uv_ary.size();iid++){
    const unsigned int id_uv = id_uv_ary[iid];
		assert( m_UseVertexSet.IsObjID(id_uv) );
		const CUseVertex& uv = m_UseVertexSet.GetObj(id_uv);
		assert( uv.id == id_uv );
		const unsigned int id_he = uv.id_he;
		assert( m_HalfEdgeSet.IsObjID(id_he) );
		const CHalfEdge& hedge = m_HalfEdgeSet.GetObj(id_he);
		assert( hedge.id == id_he );
		assert( hedge.id_uv == id_uv );
	}

	// Check UseEdge
	const std::vector<unsigned int>& id_hedge_ary = m_HalfEdgeSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_hedge_ary.size();iid++){
    const unsigned int id_he = id_hedge_ary[iid];
		assert( m_HalfEdgeSet.IsObjID(id_he) );
		const CHalfEdge& hedge = m_HalfEdgeSet.GetObj(id_he);
		assert( hedge.id == id_he );

		{
      const unsigned int id_uv1 = hedge.id_uv;
			assert( m_UseVertexSet.IsObjID(id_uv1) );
			const CUseVertex& uv = m_UseVertexSet.GetObj(id_uv1);
			assert( uv.id == id_uv1 );
		}

		unsigned int id_uv2;
		{
			const unsigned int id_he_f = hedge.id_he_f;
			assert( m_HalfEdgeSet.IsObjID(id_he_f) );
			const CHalfEdge& edge_cw = m_HalfEdgeSet.GetObj(id_he_f);
			assert( edge_cw.id == id_he_f );
			assert( edge_cw.id_he_b == id_he );
			assert( edge_cw.id_ul == hedge.id_ul );
			id_uv2 = edge_cw.id_uv;
			assert( m_UseVertexSet.IsObjID(id_uv2) );
			const CUseVertex& uv = m_UseVertexSet.GetObj(id_uv2);
			assert( uv.id == id_uv2 );
		}

		{
      const unsigned int id_he_ccw = hedge.id_he_b;
			assert( m_HalfEdgeSet.IsObjID(id_he_ccw) );
			const CHalfEdge& edge_ccw = m_HalfEdgeSet.GetObj(id_he_ccw);
			assert( edge_ccw.id == id_he_ccw );
			assert( edge_ccw.id_he_f == id_he );
			assert( edge_ccw.id_ul == hedge.id_ul );
		}

		{
      const unsigned int id_he_o = hedge.id_he_o;
			assert( m_HalfEdgeSet.IsObjID(id_he_o) );
			const CHalfEdge& edge_o = m_HalfEdgeSet.GetObj(id_he_o);
			assert( edge_o.id == id_he_o );
			assert( edge_o.id_he_o == id_he );
			assert( edge_o.id_uv == id_uv2 );
		}
	}
	
	// Check UseLoop
	const std::vector<unsigned int>& id_ul_ary = m_UseLoopSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_ul_ary.size();iid++){
    const unsigned int id_ul = id_ul_ary[iid];
		assert( m_UseLoopSet.IsObjID(id_ul) );
		const CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		assert( ul.id == id_ul );
    {
      std::set<unsigned int> passed_hedge;
      unsigned int id_he = ul.id_he;
      for(;;){
        assert( passed_hedge.find(id_he) == passed_hedge.end() );
        passed_hedge.insert(id_he);
        assert( m_HalfEdgeSet.IsObjID(id_he) );
        const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
        assert( he.id == id_he );
        assert( he.id_ul == id_ul );
        const unsigned int id_he_next = he.id_he_f;
        if( id_he_next == ul.id_he ) break;
        id_he = id_he_next;
      }      
    }
		if( ul.id_ul_p != id_ul && ul.id_ul_p != 0 ){	// 自分は子ループ
			const unsigned int id_ul_p = ul.id_ul_p;
			unsigned int id_ul2 = id_ul_p;
			bool iflag = false;
			for(;;){
				assert( m_UseLoopSet.IsObjID(id_ul2) );
				const CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
        assert( ul2.id_ul_p == id_ul_p );
        if( id_ul2 == id_ul ) iflag = true;
				id_ul2 = ul2.id_ul_c;
				if( id_ul2 == 0 ) break;
			}
			assert( iflag == true );
		}
	}
	return 0;
}

std::vector<unsigned int> CBRep::FindHalfEdge_Edge(const unsigned int& id_e) const
{
	std::vector<unsigned int> res;
	std::vector<unsigned int> id_ary_he = m_HalfEdgeSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_ary_he.size();iid++){
		unsigned int id_he = id_ary_he[iid];
		const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
		if( he.id_e != id_e ) continue;
		res.push_back(id_he);
	}
	return res;
}

std::vector<unsigned int> CBRep::FindHalfEdge_Vertex(const unsigned int& id_v) const
{
	std::vector<unsigned int> res;
	std::vector<unsigned int> id_ary_uv = m_UseVertexSet.GetAry_ObjID();
	for(unsigned int iid=0;iid<id_ary_uv.size();iid++){
		unsigned int id_uv = id_ary_uv[iid];
		const CUseVertex& uv = m_UseVertexSet.GetObj(id_uv);
		if( uv.id_v != id_v ) continue;
		unsigned int id_he = uv.id_he;
		const unsigned int id_he_ini = id_he;
		res.push_back(id_he);
		for(;;){
			const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
			const unsigned int id_he_o = he.id_he_o;
			const CHalfEdge& he_o = m_HalfEdgeSet.GetObj(id_he_o);
			const unsigned int id_he_next = he_o.id_he_f;
			if( id_he_next == id_he_ini ) break;
			id_he = id_he_next;
			res.push_back(id_he);
		}
	}

	return res;
}


// 浮遊点の削除
bool CBRep::KVEL(const unsigned int id_uv_rem)
{
	unsigned int id_he_rem, id_ul_rem;
	{
		const CUseVertex& uv = m_UseVertexSet.GetObj(id_uv_rem);
		id_he_rem = uv.id_he;
		const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he_rem);
		id_ul_rem = he.id_ul;
		assert( he.id_he_b == id_he_rem );
		assert( he.id_he_f == id_he_rem );
		assert( he.id_he_o == id_he_rem );
	}

	unsigned int id_ul_p;
	{
		const CUseLoop& ul_rem = m_UseLoopSet.GetObj(id_ul_rem);
		id_ul_p = ul_rem.id_ul_p;
		if( id_ul_p != id_ul_rem && id_ul_p != 0 ){ // assertion
			assert( m_UseLoopSet.IsObjID(id_ul_p) );
			const CUseLoop& ul_p = m_UseLoopSet.GetObj(id_ul_p);
			assert( ul_p.id_ul_p == id_ul_p );
			assert( ul_p.id_ul_c != 0 );
		}
	}

	// Leave Imput Check Section
	////////////////////////////////

	if( id_ul_p != id_ul_rem && id_ul_p != 0 ){	// 子ループリストからid_ul_remを削除
		unsigned int id_ul = id_ul_p;
		for(;;){
			CUseLoop& ul =m_UseLoopSet.GetObj(id_ul);
			const unsigned int id_ul_c = ul.id_ul_c;
			assert( id_ul_c != 0 );
			if( id_ul_c == id_ul_rem ){
				CUseLoop& ul_c =m_UseLoopSet.GetObj(id_ul_c);
				const unsigned int id_ul_cc = ul_c.id_ul_c;
				ul.id_ul_c = id_ul_cc;
				break;
			}
			id_ul = id_ul_c;
		}
	}
  
	m_UseLoopSet.DeleteObj(id_ul_rem);
	m_UseVertexSet.DeleteObj(id_uv_rem);
	m_HalfEdgeSet.DeleteObj(id_he_rem);

	return true;
}

// 浮遊点を作る
// id_ulが属する位相ループに位相頂点を付け加える
bool CBRep::MVEL(unsigned int& id_uv_add, unsigned int& id_he_add, unsigned int& id_ul_add, 
					 const unsigned int id_ul1)	
{
	id_uv_add = m_UseVertexSet.GetFreeObjID();
	id_he_add = m_HalfEdgeSet.GetFreeObjID();
	id_ul_add = m_UseLoopSet.GetFreeObjID();

	{	// UseVertexを追加
		CUseVertex uv_add(id_uv_add,id_he_add);
		const unsigned int tmp_id = m_UseVertexSet.AddObj(uv_add);
		assert( tmp_id == id_uv_add );
	}
	{	// HalfEdgeの追加
		CHalfEdge he_add(id_he_add, id_uv_add,  id_he_add,id_he_add,id_he_add,  id_ul_add);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(he_add);
		assert( tmp_id == id_he_add );
	}
	{	// UseLoopの追加
		CUseLoop ul_add(id_ul_add, id_he_add,0,id_ul1);
		const unsigned int tmp_id = m_UseLoopSet.AddObj(ul_add);
		assert( tmp_id == id_ul_add );
	}
	if( id_ul1 != 0 ){
		unsigned int id_ul = id_ul1;
		// P/Cループリストの最後を探す
		for(;;){
			assert( m_UseLoopSet.IsObjID(id_ul) );
			CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
			assert( ul.id == id_ul );
			if( ul.id_ul_c == 0 ){
				ul.id_ul_c = id_ul_add;
				break;
			}
			id_ul = ul.id_ul_c;
		}
	}
	return true;
}


bool CBRep::MEVVL(unsigned int& id_he_add1, unsigned int& id_he_add2,
					   unsigned int& id_uv_add1, unsigned int& id_uv_add2, unsigned int& id_ul_add )
{
	id_he_add1 = 0; id_he_add2 = 0;
	id_uv_add1 = 0; id_uv_add2 = 0;
	id_ul_add = 0;

	////////////////////////////////
	// leave input check section

	id_ul_add = m_UseLoopSet.GetFreeObjID();

	{
		const std::vector<unsigned int>& free_id_ary = m_HalfEdgeSet.GetFreeObjID(2);
		assert( free_id_ary.size() == 2 );
		id_he_add1 = free_id_ary[0];
		id_he_add2 = free_id_ary[1];
		assert( !m_HalfEdgeSet.IsObjID(id_he_add1) );
		assert( !m_HalfEdgeSet.IsObjID(id_he_add2) );
		assert( id_he_add1 != id_he_add2 );
	}

	{
		const std::vector<unsigned int>& free_id_ary = m_UseVertexSet.GetFreeObjID(2);
		assert( free_id_ary.size() == 2 );
		id_uv_add1 = free_id_ary[0];
		id_uv_add2 = free_id_ary[1];
		assert( !m_UseVertexSet.IsObjID(id_uv_add1) );
		assert( !m_UseVertexSet.IsObjID(id_uv_add2) );
		assert( id_uv_add1 != id_uv_add2 );
	}

	////////////////////////////////
	// enter topological change section

	{	// Add New Use Loop
		const unsigned int tmp_id = m_UseLoopSet.AddObj( CUseLoop(id_ul_add,id_he_add1,0,0) );
		assert( tmp_id == id_ul_add );
		assert( m_UseLoopSet.IsObjID(id_ul_add) );
	}
	////////////////
	{	// Add New Half Edge 1
		CHalfEdge tmp_he(id_he_add1, id_uv_add1,  id_he_add2,id_he_add2,id_he_add2,  id_ul_add);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add1 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add1) );
	}
	{	// Add New Half Edge 2
		CHalfEdge tmp_he(id_he_add2, id_uv_add2,  id_he_add1,id_he_add1,id_he_add1,  id_ul_add);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add2 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add2) );
	}
	////////////////
	{	// Add New Vertex 1 
		const unsigned int tmp_id = m_UseVertexSet.AddObj( CUseVertex(id_uv_add1,id_he_add1) );
		assert( tmp_id == id_uv_add1 );
		assert( m_UseVertexSet.IsObjID(id_uv_add1) );
	}
	{	// Add New Vertex 2
		const unsigned int tmp_id = m_UseVertexSet.AddObj( CUseVertex(id_uv_add2,id_he_add2) );
		assert( tmp_id == id_uv_add2 );
		assert( m_UseVertexSet.IsObjID(id_uv_add2) );
	}
	////////////////////////////////
	// leave topologycal change section

	assert( AssertValid_Use()==0 );
	return true;
}


// disconnect loop
bool CBRep::KEML(unsigned int& id_ul_add, const unsigned int& id_he1 )
{
	unsigned int id_he1f, id_he1b, id_uv1, id_ul1, id_he2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id = id_he1 );
		id_he1b = he.id_he_b;
		id_he1f = he.id_he_f;
		id_uv1 = he.id_uv;
		id_ul1 = he.id_ul;
		id_he2 = he.id_he_o;
	}

	unsigned int id_he2f, id_he2b, id_uv2, id_ul2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		if( !m_HalfEdgeSet.IsObjID(id_he2) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id = id_he2 );
		id_he2b = he.id_he_b;
		id_he2f = he.id_he_f;
		id_uv2 = he.id_uv;
		id_ul2 = he.id_ul;
	}
	assert( id_ul1 == id_ul2 );
	{
		assert( id_he1f != id_he1 ); // not floating vertex
		assert( id_he1b != id_he1 ); // not floating vertex
		assert( m_HalfEdgeSet.IsObjID(id_he1f) );
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he_b.id = id_he1b );
		assert( he_b.id_he_f == id_he1 );
		assert( he_b.id_ul == id_ul1 );
	}
	{
		assert( id_he2f != id_he2 ); // not floating vertex
		assert( id_he2b != id_he2 ); // not floating vertex
		assert( m_HalfEdgeSet.IsObjID(id_he2f) );
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he2b);
		assert( he_b.id = id_he2b );
		assert( he_b.id_he_f == id_he2 );
		assert( he_b.id_ul == id_ul2 );
	}
	
	////////////////////////////////
	// the data will be changed from here

	id_ul_add = m_UseLoopSet.GetFreeObjID();
	unsigned int id_ul_p;
	{
		assert( m_UseLoopSet.IsObjID(id_ul1) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
		ul.id_he = id_he1b;
		id_ul_p = ul.id_ul_p;
		if( id_ul_p != 0 ){ // assertion
			assert( m_UseLoopSet.IsObjID(id_ul_p) );
			CUseLoop& ul_p = m_UseLoopSet.GetObj(id_ul_p);
			assert( ul_p.id_ul_p == id_ul_p );
      unsigned int id_ul = id_ul1;
      for(;;){ // 子ループリストの一番最後に追加したループを付け加える
        CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
        if( ul.id_ul_c == 0 ){
          ul.id_ul_c = id_ul_add;
          break;
        }
        id_ul = ul.id_ul_c;
      }
    }      
	}
	{
		const unsigned int tmp_id = m_UseLoopSet.AddObj( CUseLoop(id_ul_add,id_he1f,0,id_ul_p) );
		assert( tmp_id == id_ul_add  );
	}
	////////////////
	{
		assert( m_UseVertexSet.IsObjID(id_uv1) );
		CUseVertex& uv = m_UseVertexSet.GetObj(id_uv1);
		uv.id_he = id_he2f;
	}
	{
		assert( m_UseVertexSet.IsObjID(id_uv2) );
		CUseVertex& uv = m_UseVertexSet.GetObj(id_uv2);
		uv.id_he = id_he1f;
	}
	////////////////
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1b);
		he.id_he_f = id_he2f;
	}
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2f) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2f);
		he.id_he_b = id_he1b;
	}
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1f) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1f);
		he.id_he_b = id_he2b;
	}
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2b);
		he.id_he_f = id_he1f;
	}
	{	// 新しくできたループ上の半辺に全て新しいループ番号をセット
		unsigned int id_he = id_he2b;
		for(;;){
			assert( m_HalfEdgeSet.IsObjID(id_he) );
			CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
			assert( he.id_ul == id_ul1 );
			he.id_ul = id_ul_add;
			id_he = he.id_he_f;
			if( id_he == id_he2b ) break;
		}
	}
	m_HalfEdgeSet.DeleteObj(id_he1);
	m_HalfEdgeSet.DeleteObj(id_he2);

	return true;
}

// ループと浮遊点をつなげる,he1,he2が浮遊点Edge
// he1は[uv1-uv2],he2は[uv2-uv1]のHalfEdgeとなる
bool CBRep::MEKL_TwoFloatingVertex(const unsigned int id_he1, const unsigned int id_he2 )
{
	unsigned int id_uv1, id_ul1, id_ul1p, id_ul1c;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id = id_he1 );
		assert( he.id_he_f == id_he1 );
		assert( he.id_he_b == id_he1 );
		assert( he.id_he_o == id_he1 );
		id_uv1 = he.id_uv;
		id_ul1 = he.id_ul;
		const CUseLoop ul = m_UseLoopSet.GetObj(id_ul1);
//		assert( ul.id_ul_p != 0 );
		id_ul1p = ul.id_ul_p;
		id_ul1c = ul.id_ul_c;
	}
	unsigned int id_uv2, id_ul2, id_ul2p, id_ul2c;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		if( !m_HalfEdgeSet.IsObjID(id_he2) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id = id_he2 );
		assert( he.id_he_f == id_he2 );
		assert( he.id_he_b == id_he2 );
		assert( he.id_he_o == id_he2 );
		id_uv2 = he.id_uv;
		id_ul2 = he.id_ul;
		const CUseLoop ul = m_UseLoopSet.GetObj(id_ul2);
//		assert( ul.id_ul_p != 0 );
		id_ul2p = ul.id_ul_p;
		id_ul2c = ul.id_ul_c;
	}
	assert( id_ul1p == id_ul2p );

	{
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		he.id_he_b = id_he2;
		he.id_he_f = id_he2;
		he.id_he_o = id_he2;
		he.id_ul = id_ul1;
	}
	{
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		he.id_he_b = id_he1;
		he.id_he_f = id_he1;
		he.id_he_o = id_he1;
		he.id_ul = id_ul1;
	}
	////////////////
	m_UseLoopSet.DeleteObj(id_ul2);
	if( id_ul1p != 0 ){
		unsigned int id_ul = id_ul1p;
		for(;;){
			CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
			if( ul.id_ul_c == id_ul2 ){
				ul.id_ul_c = id_ul2c;
				break;
			}
			id_ul = ul.id_ul_c;
			assert( id_ul != 0 );
		}
	}

	return true;
}


// 両方が端点であるEdgeを削除する。
bool CBRep::KEML_TwoFloatingVertex
(unsigned int& id_ul_add,
 const unsigned int id_he1)
{
	unsigned int id_he2;
	unsigned int id_uv1;
	unsigned int id_ul;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id == id_he1 );
		assert( he.id_he_o == he.id_he_f );
		assert( he.id_he_o == he.id_he_b );
		id_uv1 = he.id_uv;
		id_he2 = he.id_he_o;
		assert( he.id_he_b == id_he2 );
		assert( he.id_he_f == id_he2 );
		assert( id_he1 != id_he2 );
		id_ul = he.id_ul;
	}
	unsigned int id_ul_p, id_ul_c;
	{
		const CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		id_ul_p = ul.id_ul_p;
		id_ul_c = ul.id_ul_c;
	}
	unsigned int id_uv2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		const CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id == id_he2 );
		assert( he.id_he_b == id_he1 );
		assert( he.id_he_f == id_he1 );
		assert( he.id_he_o == id_he1 );
		id_uv2 = he.id_uv;
		assert( he.id_ul == id_ul );
	}

	////////////////

	id_ul_add = m_UseLoopSet.GetFreeObjID();
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		he.id_he_b = id_he1;
		he.id_he_f = id_he1;
		he.id_he_o = id_he1;
		assert( he.id_ul == id_ul );
	}
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		he.id_he_b = id_he2;
		he.id_he_f = id_he2;
		he.id_he_o = id_he2;
		he.id_ul = id_ul_add;
	}
	////////////////
	{
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		ul.id_ul_c = id_ul_add;
		assert( ul.id_ul_p == id_ul_p );
		assert( ul.id_he = id_he1 );
	}
	{
		unsigned int tmp_id 
			= m_UseLoopSet.AddObj( CUseLoop(id_ul_add,id_he2,id_ul_c,id_ul_p) );
		assert( tmp_id == id_ul_add );
	}

	return true;
}

// 片方が端点である、HalfEdgeを削除する。
// he1の起点uv1は他の辺につながっていない端点である
bool CBRep::KEML_OneFloatingVertex(unsigned int& id_ul_add, const unsigned int id_he1)
{
	id_ul_add = 0;

	unsigned int id_he1f, id_he2;
	unsigned int id_uv1;
	unsigned int id_ul1, id_ul_p;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id = id_he1 );
		id_he1f = he.id_he_f;
		id_uv1 = he.id_uv;
		id_he2 = he.id_he_o;
		assert( he.id_he_b == id_he2 );
		assert( he.id_he_f != id_he2 );
		assert( id_he1 != id_he2 );
		id_ul1 = he.id_ul;
		const CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
		id_ul_p = ul.id_ul_p;
    if( id_ul_p != 0 ){ // if this edge is outside, id_ul_p is 0
      const CUseLoop ul_p = m_UseLoopSet.GetObj(id_ul_p);
      assert( ul_p.id_ul_p == id_ul_p );
    }
	}
	unsigned int id_uv2;
	unsigned int id_he2b;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id == id_he2 );
		id_he2b = he.id_he_b;
		id_uv2 = he.id_uv;
		assert( id_he2b != id_he1f );
		assert( he.id_he_f == id_he1 );
		assert( he.id_ul == id_ul1 );
	}
	
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1f) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1f);
		assert( he.id == id_he1f );
		assert( he.id_uv == id_uv2 );
		assert( he.id_ul == id_ul1 );
	}
	////////////////

	id_ul_add = m_UseLoopSet.GetFreeObjID();
	{
		const unsigned int tmp_id = m_UseLoopSet.AddObj( CUseLoop(id_ul_add,id_he1,0,id_ul_p) );
		assert( tmp_id == id_ul_add );
	}
	if( id_ul_p !=0 ){	// 子ループリストにid_ul_addを登録
		unsigned int id_ul = id_ul1;
		for(;;){
			CUseLoop&  ul = m_UseLoopSet.GetObj(id_ul);
			assert( ul.id_ul_p==id_ul_p );
			if( ul.id_ul_c == 0 ){
				ul.id_ul_c = id_ul_add;
				break;
			}
			id_ul = ul.id_ul_c;
		}
	}
	{
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
		ul.id_he = id_he1f;
	}
	////////////////
	{
		CUseVertex& uv = m_UseVertexSet.GetObj(id_uv2);
		uv.id_he = id_he1f;
	}
	////////////////
	m_HalfEdgeSet.DeleteObj(id_he2);
	{
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		he.id_he_b = id_he1;
		he.id_he_f = id_he1;
		he.id_he_o = id_he1;
		he.id_ul = id_ul_add;
		he.id_uv = id_uv1;
		he.id_e = 0;
		he.is_same_dir = true;
	}
	{
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1f);
		he.id_he_b = id_he2b;
	}
	{
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2b);
		he.id_he_f = id_he1f;
	}
	return true;
}

// ループと浮遊点をつなげる,he1がLoop上のEdgeでhe2が浮遊点Edge
// he2は[uv2-uv1], he_add1は[uv1-uv2]のHalfEdgeとなる
bool CBRep::MEKL_OneFloatingVertex(unsigned int& id_he_add1,
	const unsigned int id_he1, const unsigned int id_he2 )
{
	unsigned int id_he1b, id_uv1, id_ul1;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id = id_he1 );
		id_he1b = he.id_he_b;
		id_uv1 = he.id_uv;
		id_ul1 = he.id_ul;
		assert( id_he1b != id_he1 ); // not floating vertex
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he_b.id = id_he1b );
		assert( he_b.id_he_f == id_he1 );
		assert( he_b.id_ul == id_ul1 );
	}
	unsigned int id_ul1p;
	{
		assert( m_UseLoopSet.IsObjID(id_ul1) );
		const CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul.id_ul_p;
	}

	unsigned int id_uv2, id_ul2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		if( !m_HalfEdgeSet.IsObjID(id_he2) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id = id_he2 );
		id_uv2 = he.id_uv;
		id_ul2 = he.id_ul;
		assert( he.id_he_f==id_he2 ); // floating vertex
		assert( he.id_he_b==id_he2 ); // floating vertex
	}

	assert( id_ul1 != id_ul2 );
	if( id_ul1 == id_ul2 ) return false;

	////////////////////////////////
	// leave input check section

	{
		id_he_add1 = m_HalfEdgeSet.GetFreeObjID();
		assert( !m_HalfEdgeSet.IsObjID(id_he_add1) );
	}

	////////////////////////////////
	// enter topological change section

	{	// Add New Half Edge 1
		CHalfEdge tmp_he(id_he_add1,  id_uv1,  id_he2,id_he1b,id_he2,  id_ul1);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add1 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add1) );
	}
	{	// Modefy Half Edge 1b
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he.id == id_he1b );
		he.id_he_f = id_he_add1;
	}
	{	// Modefy Half Edge 1
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id == id_he1 );
		he.id_he_b = id_he2;
	}
	{	// Modefy Half Edge 2 
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id == id_he2 );
		he.id_he_b = id_he_add1;
		he.id_he_f  = id_he1;
		he.id_he_o = id_he_add1;
		he.id_ul = id_ul1;
		he.id_uv = id_uv2;
	}
	{	// Delete Childe Loop List
		unsigned int id_ul_cc;
		{
			assert( m_UseLoopSet.IsObjID(id_ul2) );
			const CUseLoop& ul = m_UseLoopSet.GetObj(id_ul2);
			id_ul_cc = ul.id_ul_c;
		}
    if( id_ul1p != 0 ){
      unsigned int id_ul = id_ul1p;
      for(;;){
        assert( m_UseLoopSet.IsObjID(id_ul) );
        CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
        const unsigned int id_ul_c = ul.id_ul_c;
        if( id_ul_c==id_ul2 ){
          ul.id_ul_c = id_ul_cc;
          break;
        }
        if( id_ul_c == 0 ){ break; }	// if floating vertex is outside
        assert( id_ul_c != id_ul_cc );
        id_ul = id_ul_c;
      }      
    }
	}
	m_UseLoopSet.DeleteObj(id_ul2);
	assert( !m_UseLoopSet.IsObjID(id_ul2) );
	assert( this->AssertValid_Use()==0 );

	return true;
}
	
// he1の起点uv1とhe2の起点uv2を結んで、２つのループをつなげる
// ul2は削除されてul1に統一される。
// he1は[uv1-uv2]、he2は[uv2-uv1]
bool CBRep::MEKL(unsigned int& id_he_add1, unsigned int& id_he_add2,  
	const unsigned int id_he1, const unsigned int  id_he2 )
{
	id_he_add1 = 0; id_he_add2 = 0;

	unsigned int id_he1b, id_uv1, id_ul1;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id = id_he1 );
		id_he1b = he.id_he_b;
		id_uv1 = he.id_uv;
		id_ul1 = he.id_ul;
		assert( id_he1b != id_he1 ); // not floating vertex
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he_b.id = id_he1b );
		assert( he_b.id_he_f == id_he1 );
		assert( he_b.id_ul == id_ul1 );
	}

	unsigned int id_he2b, id_uv2, id_ul2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		if( !m_HalfEdgeSet.IsObjID(id_he2) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id = id_he2 );
		id_he2b = he.id_he_b;
		id_uv2 = he.id_uv;
		id_ul2 = he.id_ul;
		assert( id_he2b != id_he2 ); // not floating vertex
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he2b);
		assert( he_b.id = id_he2b );
		assert( he_b.id_he_f == id_he2 );
		assert( he_b.id_ul == id_ul2 );
	}

	assert( id_ul1 != id_ul2 );
	if( id_ul1 == id_ul2 ) return false;

	unsigned int id_ul1p, id_ul1c;
	{
		assert( m_UseLoopSet.IsObjID(id_ul1) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul.id_ul_p;
		id_ul1c = ul.id_ul_c;
	}

	unsigned int id_ul2p, id_ul2c;
	{
		assert( m_UseLoopSet.IsObjID(id_ul2) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul2);
		id_ul2p = ul.id_ul_p;
		id_ul2c = ul.id_ul_c;
	}

	////////////////////////////////
	// leave input check section

	{
		const std::vector<unsigned int>& free_edge_id_ary = m_HalfEdgeSet.GetFreeObjID(2);
		assert( free_edge_id_ary.size() == 2 );
		id_he_add1 = free_edge_id_ary[0];
		id_he_add2 = free_edge_id_ary[1];
		assert( !m_HalfEdgeSet.IsObjID(id_he_add1) );
		assert( !m_HalfEdgeSet.IsObjID(id_he_add2) );
		assert( id_he_add1 != id_he_add2 );
	}

	////////////////////////////////
	// enter topological change section

	{	// Add New Half Edge 1
		CHalfEdge tmp_he(id_he_add1,  id_uv1,  id_he2,id_he1b,id_he_add2,  id_ul1);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add1 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add1) );
	}
	{	// Add New Half Edge 2
		CHalfEdge tmp_he(id_he_add2,  id_uv2,  id_he1,id_he2b,id_he_add1,  id_ul1);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add2 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add2) );
	}
	{	// Modefy Half Edge 1b
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he.id == id_he1b );
		he.id_he_f = id_he_add1;
	}
	{	// Modefy Half Edge 2b
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2b);
		assert( he.id == id_he2b );
		he.id_he_f = id_he_add2;
	}
	{	// Modefy Half Edge 1
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he1);
		assert( he.id == id_he1 );
		he.id_he_b = id_he_add2;
	}
	{	// Modefy Half Edge 2 
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he2);
		assert( he.id == id_he2 );
		he.id_he_b = id_he_add1;
	}
	{	// Modefy Old Half Edges ( Connecting Loop )
    int id_he = id_he_add2; // 途中で-1になるかもしれないからint型
		for(;;){
			assert( m_HalfEdgeSet.IsObjID(id_he) );
			CHalfEdge& edge = m_HalfEdgeSet.GetObj(id_he);
      assert( (int)edge.id == id_he );
			edge.id_ul = id_ul1;
			id_he = edge.id_he_f;
      if( id_he == (int)id_he_add2 ) break;
		}
	}
	////////////////
	m_UseLoopSet.DeleteObj(id_ul2);
	assert( !m_UseLoopSet.IsObjID(id_ul2) );
	if( id_ul2p != id_ul2 && id_ul2p != 0 ){	// ul2は子ループ
		unsigned int id_ul;
		if( id_ul1p != id_ul1 ){	// ２つのループの外側に辺を追加した場合
			assert( id_ul1p == id_ul2p );
			id_ul = id_ul1p;
		}
		else{ // 親ループul1の内側の子ループul2と辺を結んだ場合
			id_ul = id_ul1;
		}
		for(;;){
			assert( m_UseLoopSet.IsObjID(id_ul) );
			CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
			if( ul.id_ul_c == id_ul2 ){
				ul.id_ul_c = id_ul2c;
				break;
			}
			id_ul = ul.id_ul_c;
			assert( id_ul != 0 );
		}
	}
	// ul2は親ループ
	else if( id_ul1p != 0 ){ // 親ループul2の内側の子ループul1と辺を結んだ場合
		// ul1を親にして、ul1を削除
		if( id_ul2c == id_ul1 ){
			{
				assert( m_UseLoopSet.IsObjID(id_ul1) );
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul1);
				ul.id_ul_p = id_ul1;
			}
			unsigned int id_ul = id_ul1c;
			for(;;){
				if( id_ul == 0 ) break;
				assert( m_UseLoopSet.IsObjID(id_ul) );
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul2 );
				ul.id_ul_p = id_ul1;
				id_ul = ul.id_ul_c;
			}
		}
		else{
			assert( id_ul2c != 0 );
			unsigned int id_ul = id_ul2c;
			for(;;){
				assert( m_UseLoopSet.IsObjID(id_ul) );
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul2 );
				const unsigned int id_ul_c = ul.id_ul_c;
				assert( ul.id_ul_p == id_ul2 );
				if( id_ul == id_ul1 ){
					ul.id_ul_p = id_ul;
					ul.id_ul_c = id_ul2c;
				}
				else{
					ul.id_ul_p = id_ul1;
					if( ul.id_ul_c == id_ul1 ){
						ul.id_ul_c = id_ul1c;
					}
				}
				id_ul = id_ul_c;
				if( id_ul == 0 ) break;
			}
		}
	}

	assert( this->AssertValid_Use()==0 );

	return true;
}

// ul2側のループが消去される
bool CBRep::KEL(const unsigned int id_he1){

	unsigned int id_uv1;
	unsigned int id_he2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id = id_he1 );
		id_he2 = he1.id_he_o;
		id_uv1 = he1.id_uv;
	}
	unsigned int id_uv2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id = id_he2 );
		assert( he2.id_he_o == id_he1 );
		id_uv2 = he2.id_uv;
	}

	unsigned int id_he1b, id_he1f, id_ul1;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id = id_he1 );
		assert( he1.id_uv == id_uv1 );
		id_ul1 =  he1.id_ul;
		id_he1b = he1.id_he_b;
		id_he1f = he1.id_he_f;
		{	// Assert ルーティン
			assert( m_HalfEdgeSet.IsObjID(id_he1b) );
			CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
			assert( he1b.id == id_he1b );
			assert( he1b.id_he_f == id_he1 );
			assert( he1b.id_ul == id_ul1 );
		}
		{	// Assert ルーティン
			assert( m_HalfEdgeSet.IsObjID(id_he1f) );
			CHalfEdge& he1f = m_HalfEdgeSet.GetObj(id_he1f);
			assert( he1f.id == id_he1f );
			assert( he1f.id_he_b == id_he1 );
			assert( he1f.id_ul == id_ul1 );
			assert( he1f.id_uv == id_uv2 );
		}
		{
			assert( m_UseLoopSet.IsObjID(id_ul1) );
		}
	}

	unsigned int id_he2b, id_he2f, id_ul2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id = id_he2 );
		assert( he2.id_uv == id_uv2 );
		id_ul2 = he2.id_ul;
		id_he2b = he2.id_he_b;
		id_he2f = he2.id_he_f;
		{	// Assert ルーティン
			assert( m_HalfEdgeSet.IsObjID(id_he2b) );
			const CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
			assert( he2b.id == id_he2b );
			assert( he2b.id_he_f == id_he2 );
			assert( he2b.id_ul == id_ul2 );
		}
		{	// Assert ルーティン
			assert( m_HalfEdgeSet.IsObjID(id_he2f) );
			const CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
			assert( he2f.id == id_he2f );
			assert( he2f.id_he_b == id_he2 );
			assert( he2f.id_ul == id_ul2 );
			assert( he2f.id_uv == id_uv1 );
		}
		{
			assert( m_UseLoopSet.IsObjID(id_ul2) );
		}
	}

	assert( id_he2b != id_he1f );
	assert( id_he1b != id_he2f );
	assert( id_ul1 != id_ul2 );

	// Reave Input Check Section
	////////////////////////////////


	////////////////
	m_HalfEdgeSet.DeleteObj(id_he1);
	m_HalfEdgeSet.DeleteObj(id_he2);
	{
		CHalfEdge& he1f = m_HalfEdgeSet.GetObj(id_he1f);
		he1f.id_he_b = id_he2b;
	}
	{
		CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
		he1b.id_he_f = id_he2f;
	}
	{
		CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
		he2f.id_he_b = id_he1b;
	}
	{
		CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
		he2b.id_he_f = id_he1f;
	}
	{	// Modefy Old Half Edges ( Connecting Loop )
		int id_he = id_he1f;
		for(;;){
			assert( m_HalfEdgeSet.IsObjID(id_he) );
			CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
            assert( (int)he.id == id_he );
			assert( he.id_ul == id_ul1 || he.id_ul == id_ul2 );
			////////////////
			he.id_ul = id_ul1;
			////////////////
			id_he = he.id_he_f;
            if( id_he == (int)id_he1f ) break;
		}
	}
	////////////////
	{
		CUseVertex& uv1 = m_UseVertexSet.GetObj(id_uv1);
		uv1.id_he = id_he2f;
	}
	{
		CUseVertex& uv2 = m_UseVertexSet.GetObj(id_uv2);
		uv2.id_he = id_he1f;
	}
	////////////////
	{
		CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		ul1.id_he = id_he1f;
	}
	unsigned int id_ul1p,id_ul1c;
	{
		CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul1.id_ul_p;
		id_ul1c = ul1.id_ul_c;
	}
	unsigned int id_ul2p,id_ul2c;
	{
		CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		id_ul2p = ul2.id_ul_p;
		id_ul2c = ul2.id_ul_c;
	}
	m_UseLoopSet.DeleteObj(id_ul2);
	if( id_ul2p == id_ul2 ){	// ul1,ul2の全ての子ループの親はid_ul1pになる
//		std::cout << "parent 2 " << std::endl;
		{	// ul2に属していた子ループの親ul1pに移し変える
			unsigned int id_ul = id_ul2c;
			for(;;){
				if( id_ul == 0 ) break;
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul2 );
				////////////////
				ul.id_ul_p = id_ul1p;
				////////////////
				id_ul = ul.id_ul_c;
			}
		}
		{
			CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
			ul1.id_he = id_he1f;
			unsigned int id_ul = id_ul1;
			for(;;){	// ul2に属していた子ループをul1に移し変える
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul1p || (ul.id_ul_p == 0&&id_ul==id_ul1) );
				id_ul = ul.id_ul_c;
				if( id_ul == 0 ){	// 最後に付け足す
					ul.id_ul_c = id_ul2c;
					break;
				}
			}
		}
	}
	else{	// ul1と、ul1,ul2の全ての子ループの親はid_ul2pになる
//		std::cout << "parent 1 " << std::endl;
		assert( id_ul1p == id_ul1 );
		{	// ul1に属していた子ループの親ul2pに移し変える
			unsigned int id_ul = id_ul1;
			for(;;){
				assert( m_UseLoopSet.IsObjID(id_ul) );
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul1 || ( ul.id_ul_p==0&&id_ul==id_ul1) );
				////////////////
				ul.id_ul_p = id_ul2p; // 親をid_ul2pにする
				////////////////
				id_ul = ul.id_ul_c;
				if( id_ul == 0 ) break;
			}
		}
		if( id_ul2p != 0 ){
			unsigned int id_ul = id_ul2p;
			for(;;){
				assert( m_UseLoopSet.IsObjID(id_ul) );
				CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
				assert( ul.id_ul_p == id_ul2p || (ul.id_ul_p==0&&id_ul==id_ul2p) );
				if( ul.id_ul_c == id_ul2 ){
					ul.id_ul_c = id_ul2c;
				}
				if( ul.id_ul_c == 0 ){
					ul.id_ul_c = id_ul1;
					break;
				}
				id_ul = ul.id_ul_c;
			}
		}
	}
	return true;
}

// ループを２つに分ける
bool CBRep::MEL(unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_ul_add,
		 const unsigned int id_he1, const unsigned int id_he2)
{  
	id_he_add1 = 0; id_he_add2 = 0;
	id_ul_add = 0;

	unsigned int id_ul;
	unsigned int id_he1b, id_uv1;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id = id_he1 );
		id_ul =  he1.id_ul;
		id_uv1 = he1.id_uv;
		id_he1b = he1.id_he_b;
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he1b.id == id_he1b );
		assert( he1b.id_he_f == id_he1 );
		assert( he1b.id_ul == id_ul );
	}

	unsigned int id_he2b, id_uv2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		if( !m_HalfEdgeSet.IsObjID(id_he2) ) return false;
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id = id_he2 );
		assert( he2.id_ul == id_ul );
		id_uv2 = he2.id_uv;
		id_he2b = he2.id_he_b;
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
		assert( he2b.id == id_he2b );
		assert( he2b.id_he_f == id_he2 );
		assert( he2b.id_ul == id_ul );
	}

	////////////////////////////////
	// leave input check section

  id_ul_add = m_UseLoopSet.GetFreeObjID();

	{
		const std::vector<unsigned int>& free_edge_id_ary = m_HalfEdgeSet.GetFreeObjID(2);
		assert( free_edge_id_ary.size() == 2 );
		id_he_add1 = free_edge_id_ary[0];
		id_he_add2 = free_edge_id_ary[1];
		assert( !m_HalfEdgeSet.IsObjID(id_he_add1) );
		assert( !m_HalfEdgeSet.IsObjID(id_he_add2) );
		assert( id_he_add1 != id_he_add2 );
	}

	////////////////////////////////
	// enter topological change section

	{	// Add New Use Loops
    const int tmp_id = m_UseLoopSet.AddObj( CUseLoop(id_ul_add,id_he_add2,0,0) );
    assert( (int)id_ul_add == tmp_id );
		assert( m_UseLoopSet.IsObjID(id_ul_add) );
	}
	{	// Modefy Old Use Loop ( Connecting Edge )
		assert( m_UseLoopSet.IsObjID(id_ul) );
		CUseLoop& loop = m_UseLoopSet.GetObj(id_ul);    
		assert( loop.id == id_ul );
		loop.id_he = id_he_add1;
	}
	////////////////
	{	// Add New Half Edge 1
		CHalfEdge tmp_he(id_he_add1,  id_uv1,  id_he2,id_he1b,id_he_add2,  id_ul);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add1 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add1) );
	}
	{	// Add New Half Edge 2
		CHalfEdge tmp_he(id_he_add2,  id_uv2,  id_he1,id_he2b,id_he_add1,  id_ul_add);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add2 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add2) );
	}
	{	// Modefy Half Edge 1
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
		assert( he1b.id == id_he1b );
		he1b.id_he_f = id_he_add1;
	}
	{	// Modefy Half Edge 4 
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id == id_he2 );
		he2.id_he_b = id_he_add1;
	}
	{	// Modefy Half Edge 3
		assert( m_HalfEdgeSet.IsObjID(id_he2b) );
		CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
		assert( he2b.id == id_he2b );
		he2b.id_he_f = id_he_add2;
	}
	{	// Modefy Half Edge 2
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id == id_he1 );
		he1.id_he_b = id_he_add2;
	}
	{	// Modefy Old Half Edges ( Connecting Loop )
		int id_he = id_he_add2;
		for(;;){
			assert( m_HalfEdgeSet.IsObjID(id_he) );
			CHalfEdge& edge = m_HalfEdgeSet.GetObj(id_he);
      assert( (int)edge.id == id_he );
			edge.id_ul = id_ul_add;
			id_he = edge.id_he_f;
      if( id_he == (int)id_he_add2 ) break;
		}
	}

	////////////////////////////////
	// leave topological change section
	return true;
}

// id_heの起点から、線分を伸ばす
bool CBRep::MEV(unsigned int& id_he_add1,unsigned int& id_he_add2, unsigned int& id_uv_add, 
		 const unsigned int id_he )
{
	id_he_add1 = 0;
	id_he_add2 = 0;
	id_uv_add = 0;

	unsigned int id_ul;
	unsigned int id_he_b, id_uv;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he) );
		if( !m_HalfEdgeSet.IsObjID(id_he) ) return false;
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
		assert( he.id = id_he );
		id_ul = he.id_ul;
		id_uv = he.id_uv;
		id_he_b = he.id_he_b;
		assert( m_HalfEdgeSet.IsObjID(id_he_b) );
		CHalfEdge& he_b = m_HalfEdgeSet.GetObj(id_he_b);
		assert( he_b.id = id_he_b );
		assert( he_b.id_he_f == id_he );
		assert( he_b.id_ul == id_ul );
	}

	////////////////////////////////
	// leave input check section

	{
		const std::vector<unsigned int>& free_edge_id_ary = m_HalfEdgeSet.GetFreeObjID(2);
		assert( free_edge_id_ary.size() == 2 );
		id_he_add1 = free_edge_id_ary[0];
		id_he_add2 = free_edge_id_ary[1];
		assert( !m_HalfEdgeSet.IsObjID(id_he_add1) );
		assert( !m_HalfEdgeSet.IsObjID(id_he_add2) );
		assert( id_he_add1 != id_he_add2 );
	}

	id_uv_add = m_UseVertexSet.GetFreeObjID();

	////////////////////////////////
	// enter topological change section

	{	// Add New Half Edge 1
		CHalfEdge tmp_he(id_he_add1,  id_uv,  id_he_add2,id_he_b,id_he_add2,  id_ul);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add1 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add1) );
	}
	{	// Add New Half Edge 2
		CHalfEdge tmp_he(id_he_add2,  id_uv_add,  id_he,id_he_add1,id_he_add1,  id_ul);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(tmp_he);
		assert( tmp_id == id_he_add2 );
		assert( m_HalfEdgeSet.IsObjID(id_he_add2) );
	}
	{	// Modefy Half Edge1
		assert( m_HalfEdgeSet.IsObjID(id_he_b) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he_b);
		assert( he.id == id_he_b );
		he.id_he_f = id_he_add1;
	}
	{	// Modefy Half Edge2
		assert( m_HalfEdgeSet.IsObjID(id_he) );
		CHalfEdge& he = m_HalfEdgeSet.GetObj(id_he);
		assert( he.id == id_he );
		he.id_he_b = id_he_add2;
	}
	////////////////
	{	// Add New Vertex
		const unsigned int tmp_id = m_UseVertexSet.AddObj( CUseVertex(id_uv_add,id_he_add2) );
		assert( tmp_id == id_uv_add );
		assert( m_UseVertexSet.IsObjID(id_uv_add) );
	}
	{	// Modefy Vertex
		assert( m_UseVertexSet.IsObjID(id_uv) );
		CUseVertex& uv = m_UseVertexSet.GetObj(id_uv);
		uv.id_he = id_he_add1;
	}
	////////////////////////////////
	// leave topological change section

	return true;
}
	

// id_heの起点を消去して２つの辺を１つにする。
bool CBRep::KVE( unsigned int& id_he1 )
{
	unsigned int id_uv1;
	unsigned int id_he2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		if( !m_HalfEdgeSet.IsObjID(id_he1) ) return false;
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id = id_he1 );
		id_he2 = he1.id_he_o;
		id_uv1 = he1.id_uv;
	}
	unsigned int id_uv2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id = id_he2 );
		assert( he2.id_he_o == id_he1 );
		id_uv2 = he2.id_uv;
	}

	unsigned int id_he1b, id_he1f, id_ul1;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he1) );
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		assert( he1.id = id_he1 );
		assert( he1.id_uv == id_uv1 );
		id_ul1 =  he1.id_ul;
		id_he1b = he1.id_he_b;
		id_he1f = he1.id_he_f;
		{	// Assert ルーティン he1b
			assert( m_HalfEdgeSet.IsObjID(id_he1b) );
			const CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
			assert( he1b.id == id_he1b );
			assert( he1b.id_he_f == id_he1 );
			assert( he1b.id_ul == id_ul1 );
		}
		{	// Assert ルーティン he1f
			assert( m_HalfEdgeSet.IsObjID(id_he1f) );
			const CHalfEdge& he1f = m_HalfEdgeSet.GetObj(id_he1f);
			assert( he1f.id == id_he1f );
			assert( he1f.id_he_b == id_he1 );
			assert( he1f.id_ul == id_ul1 );
			assert( he1f.id_uv == id_uv2 );
		}
		assert( m_UseLoopSet.IsObjID(id_ul1) );
	}

	unsigned int id_he2b, id_he2f, id_ul2;
	{
		assert( m_HalfEdgeSet.IsObjID(id_he2) );
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		assert( he2.id = id_he2 );
		assert( he2.id_uv == id_uv2 );
		id_ul2 = he2.id_ul;
		id_he2b = he2.id_he_b;
		id_he2f = he2.id_he_f;
		{	// Assert ルーティン he2b
			assert( m_HalfEdgeSet.IsObjID(id_he2b) );
			const CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
			assert( he2b.id == id_he2b );
			assert( he2b.id_he_f == id_he2 );
			assert( he2b.id_ul == id_ul2 );
		}
		{	// Assert ルーティン he2f
			assert( m_HalfEdgeSet.IsObjID(id_he2f) );
			const CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
			assert( he2f.id == id_he2f );
			assert( he2f.id_he_b == id_he2 );
			assert( he2f.id_ul == id_ul2 );
			assert( he2f.id_uv == id_uv1 );
		}
		assert( m_UseLoopSet.IsObjID(id_ul2) );
	}

	{	// he1bとhe2fが向かい合っているかをチェック
		assert( m_HalfEdgeSet.IsObjID(id_he1b) );
		const CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
		if( he1b.id_he_o != id_he2f ) return false;
		{
			assert( m_HalfEdgeSet.IsObjID(id_he2f) );
			const CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
			assert( he2f.id_he_o == id_he1b );
		}
	}

	assert( id_he1b != id_he2f );

	// Leave Input Check Section
	////////////////////////////////

	if( id_he1 != id_he2b ){	// uv2が端点でなかった場合
        {
            CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
            he1b.id_he_f = id_he1f;
        }
        {
            CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
            he2f.id_uv = id_uv2;
            he2f.id_he_b = id_he2b;
        }
        {
            CHalfEdge& he1f = m_HalfEdgeSet.GetObj(id_he1f);
            he1f.id_he_b = id_he1b;
        }
        {
            CHalfEdge& he2b = m_HalfEdgeSet.GetObj(id_he2b);
            he2b.id_he_f = id_he2f;
        }
	}
	else{	// uv2が端点の場合
		assert( id_he2 == id_he1f );
		{
			CHalfEdge& he1b = m_HalfEdgeSet.GetObj(id_he1b);
			he1b.id_he_f = id_he2f;
		}
		{
			CHalfEdge& he2f = m_HalfEdgeSet.GetObj(id_he2f);
			he2f.id_uv = id_uv2;
			he2f.id_he_b = id_he1b;
		}
	}
	m_HalfEdgeSet.DeleteObj(id_he1);	// id_he1 == id_he2bの可能性があるので最後に消す
	m_HalfEdgeSet.DeleteObj(id_he2);	// id_he2 == id_he1fの可能性があるので最後に消す

	////////////////
	m_UseVertexSet.DeleteObj(id_uv1);
	{
		CUseVertex& uv2 = m_UseVertexSet.GetObj(id_uv2);
		uv2.id_he = id_he2f;
	}

	////////////////
	{
		CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		ul1.id_he = id_he1b;
	}
	{
		CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		ul2.id_he = id_he2f;
	}

	return true;
}

// 辺を２つに分ける
// 入力ではhe2はhe1に向かいあう半辺
// 出力ではhe1の前にhe_add1,he1の向かいにhe_add2,he_add2の後ろにhe2
bool CBRep::MVE(unsigned int& id_he_add1,unsigned int& id_he_add2, unsigned int& id_uv_add, 
					const unsigned int id_he1)
{
	const CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
	assert( he1.id == id_he1 );
	const unsigned int id_he2 = he1.id_he_o;
	const CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
	assert( he2.id == id_he2 );

	const unsigned int id_ul1 = he1.id_ul;
	const unsigned int id_ul2 = he2.id_ul;

	const unsigned int id_he1_cw = he1.id_he_f;
	const unsigned int id_he2_cw = he2.id_he_f;

	id_uv_add = m_UseVertexSet.GetFreeObjID();
	const std::vector<unsigned int> free_id_he_ary = m_HalfEdgeSet.GetFreeObjID(2);
	assert( free_id_he_ary.size() == 2 );
	id_he_add1 = free_id_he_ary[0];
	id_he_add2 = free_id_he_ary[1];
	{	// UseVertexの追加
		CUseVertex uv_add(id_uv_add,id_he_add1);
		const unsigned int tmp_id = m_UseVertexSet.AddObj(uv_add);
		assert( tmp_id == id_uv_add );
	}
	{	// HalfEdge1の追加
		CHalfEdge he_add1(id_he_add1,  id_uv_add,  id_he1_cw,id_he1,id_he2,  id_ul1);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(he_add1);
		assert( tmp_id == id_he_add1 );
	}
	{	// HalfEdge2の追加
		CHalfEdge he_add2(id_he_add2,  id_uv_add,  id_he2_cw,id_he2,id_he1,  id_ul2);
		const unsigned int tmp_id = m_HalfEdgeSet.AddObj(he_add2);
		assert( tmp_id == id_he_add2 );
	}
	{
		CHalfEdge& he1 = m_HalfEdgeSet.GetObj(id_he1);
		he1.id_he_f = id_he_add1;
		he1.id_he_o = id_he_add2;
	}
	{
		CHalfEdge& he2 = m_HalfEdgeSet.GetObj(id_he2);
		he2.id_he_f = id_he_add2;
		he2.id_he_o = id_he_add1;
		he2.id_e = 0;
	}
	{
		assert( m_HalfEdgeSet.IsObjID( id_he1_cw ) );
		CHalfEdge& he1_cw  = m_HalfEdgeSet.GetObj(id_he1_cw);
		assert( he1_cw.id == id_he1_cw  );
		he1_cw.id_he_b = id_he_add1;
	}
	{
		assert( m_HalfEdgeSet.IsObjID( id_he2_cw ) );
		CHalfEdge& he2_cw  = m_HalfEdgeSet.GetObj(id_he2_cw);
		assert( he2_cw.id == id_he2_cw  );
		he2_cw.id_he_b = id_he_add2;
	}
	return true;
}

// 位相ループを動かす、id_ul1のループをid_ul2の親の子にする。
// id_ul1は親ループであってはならない。
// id_ul2は親でなくてもいい。
bool CBRep::MoveUseLoop(unsigned int id_ul1, unsigned int id_ul2){
	assert( m_UseLoopSet.IsObjID(id_ul1) );
	assert( m_UseLoopSet.IsObjID(id_ul2) );
	unsigned int id_ul1p, id_ul1c;
	{
		const CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul1.id_ul_p;
		assert( id_ul1p != id_ul1 );	// id_ul1は親であってはならない
		id_ul1c = ul1.id_ul_c;
	}
	unsigned int id_ul2c, id_ul2p;
	{
		const CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		id_ul2c = ul2.id_ul_c;
		id_ul2p = ul2.id_ul_p;
	}

	////////////////
	CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
	ul1.id_ul_p = id_ul2p;
	ul1.id_ul_c = 0;
	////////////////
	{
		unsigned int id_ul = id_ul2;
		for(;;){
			CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
			assert( ul.id_ul_p == id_ul2p || (id_ul==id_ul2p&&ul.id_ul_p==0) );
			if( ul.id_ul_c == 0 ){
				ul.id_ul_c = id_ul1;
				break;
			}
			id_ul = ul.id_ul_c;
		}
	}
	////////////////
	{
		unsigned int id_ul = id_ul1p;
		for(;;){
			CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
			assert( ul.id_ul_p == id_ul1p || (id_ul==id_ul1p&&ul.id_ul_p==0) );
			if( ul.id_ul_c == id_ul1 ){
				ul.id_ul_c = id_ul1c;
				break;
			}
			assert( ul.id_ul_c != 0 );
			id_ul = ul.id_ul_c;
		}
	}

	return true;
}

bool CBRep::SwapUseLoop(unsigned int id_ul1, unsigned int id_ul2){
	assert( m_UseLoopSet.IsObjID(id_ul1) );
	assert( m_UseLoopSet.IsObjID(id_ul2) );
	unsigned int id_ul1p;
	{
		const CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul1.id_ul_p;
	}
	unsigned int id_ul2p;
	{
		const CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		id_ul2p = ul2.id_ul_p;
	}
	if(      id_ul1p == id_ul1 && id_ul2p == id_ul1 ){ 
		return SwapUseLoop_CP_SameLoop(id_ul2,id_ul1); 
	}
	else if( id_ul2p == id_ul2 && id_ul1p == id_ul2 ){ 
		return SwapUseLoop_CP_SameLoop(id_ul1,id_ul2); 
	}
	////////////////
	else if( (id_ul1p == id_ul1 || id_ul1p == 0) && id_ul2p != id_ul1 ){
		return SwapUseLoop_CP_DifferentLoop(id_ul2,id_ul1); 
	}
	else if( (id_ul2p == id_ul2 || id_ul2p == 0) && id_ul1p != id_ul2 ){
		return SwapUseLoop_CP_DifferentLoop(id_ul1,id_ul2); 
	}
	////////////////
	assert(0);
	return true;
}


bool CBRep::SwapUseLoop_CP_DifferentLoop(unsigned int id_ul1, unsigned int id_ul2){
	assert( m_UseLoopSet.IsObjID(id_ul1) );
	assert( m_UseLoopSet.IsObjID(id_ul2) );
	unsigned int id_ul1p;
	unsigned int id_ul1c;
	{
		const CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		id_ul1p = ul1.id_ul_p;
		assert( id_ul1p != id_ul1 );
		assert( id_ul2 != id_ul1p );
		id_ul1c = ul1.id_ul_c;
	}
	unsigned int id_ul2c;
	{
		assert( m_UseLoopSet.IsObjID(id_ul2) );
		const CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		assert( ul2.id_ul_p == id_ul2 || ul2.id_ul_p == 0 );
		id_ul2c = ul2.id_ul_c;
	}
	unsigned int id_ul = id_ul1p;
	for(;;){
		assert( m_UseLoopSet.IsObjID(id_ul) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		assert( ul.id == id_ul );
		unsigned int id_ul_c = ul.id_ul_c;
//		unsigned int id_ul_p = ul.id_ul_p;
		////////////////
		if( id_ul == id_ul1 ){
			ul.id_ul_c = id_ul2c;
			ul.id_ul_p = id_ul;
		}
		else{
			if( ul.id_ul_c == id_ul1 ){ ul.id_ul_c = id_ul2; }
			assert( ul.id_ul_p == id_ul1p || ( id_ul==id_ul1p&&ul.id_ul_p==0) );
		}
		////////////////
		if( id_ul_c == 0 ) break;
		id_ul = id_ul_c;
	}
	id_ul = id_ul2;
	for(;;){
		assert( m_UseLoopSet.IsObjID(id_ul) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		assert( ul.id == id_ul );
		unsigned int id_ul_c = ul.id_ul_c;
//		unsigned int id_ul_p = ul.id_ul_p;
		////////////////
		if( id_ul == id_ul2 ){
			ul.id_ul_c = id_ul1c;
			ul.id_ul_p = id_ul1p;
		}
		else{
			ul.id_ul_p = id_ul1;
		}
		////////////////
		if( id_ul_c == 0 ) break;
		id_ul = id_ul_c;
	}

	return true;
}


bool CBRep::SwapUseLoop_CP_SameLoop(unsigned int id_ul1, unsigned int id_ul2){
	assert( m_UseLoopSet.IsObjID(id_ul1) );
	assert( m_UseLoopSet.IsObjID(id_ul2) );
	unsigned int id_ul1c;
	{
		const CUseLoop& ul1 = m_UseLoopSet.GetObj(id_ul1);
		assert( id_ul2 == ul1.id_ul_p );
		id_ul1c = ul1.id_ul_c;
	}
	unsigned int id_ul2c;
	{
		assert( m_UseLoopSet.IsObjID(id_ul2) );
		const CUseLoop& ul2 = m_UseLoopSet.GetObj(id_ul2);
		assert( ul2.id_ul_p == id_ul2 );
		id_ul2c = ul2.id_ul_c;
		assert( id_ul2c != 0 );
	}

	unsigned int id_ul = id_ul2;
	for(;;){
		assert( m_UseLoopSet.IsObjID(id_ul) );
		CUseLoop& ul = m_UseLoopSet.GetObj(id_ul);
		assert( ul.id == id_ul );
		unsigned int id_ul_c = ul.id_ul_c;
//		unsigned int id_ul_p = ul.id_ul_p;
		////////////////
		if( id_ul == id_ul2 ){ 
			assert( id_ul1c != id_ul2 );
			ul.id_ul_c = id_ul1c;
			ul.id_ul_p = id_ul1;
		}
		else if( id_ul == id_ul1 ){
			if( id_ul2c != id_ul1 ){ ul.id_ul_c = id_ul2c; }
			else{ ul.id_ul_c = id_ul2; }
			ul.id_ul_p = id_ul;
		}
		else{
			assert( ul.id_ul_c != id_ul2 );
			if( ul.id_ul_c == id_ul1 ){ ul.id_ul_c = id_ul2; }
			ul.id_ul_p = id_ul1;
		}
		////////////////
		if( id_ul_c == 0 ) break;
		id_ul = id_ul_c;
	}

	return true;
}
