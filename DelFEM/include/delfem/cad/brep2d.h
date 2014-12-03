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
@brief B-repを用いた位相情報格納クラス(Cad::CBrep)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(B_REP_2D_H)
#define B_REP_2D_H

#ifdef __VISUALC__
	#pragma warning( disable : 4786 )
#endif

#include <map>
#include "delfem/cad/brep.h"
#include "delfem/cad_com.h"
#include "delfem/serialize.h"

namespace Cad
{

class CBRepSurface
{
public:
	//! iterator goes edge and vetex around one specific loop
	class CItrLoop : public Cad::IItrLoop
	{
	public:
		CItrLoop(const CBRepSurface& pBRep2D, unsigned int id_l);
		CItrLoop(const CBRepSurface& pBRep2D, unsigned int id_he, unsigned int id_ul);
		void Begin();		
		bool IsEnd() const;
		void operator++(); //!< move to next edge
		void operator++(int n);	//!< dummy operator(to implement ++)
		bool GetIdEdge(unsigned int& id_e, bool& is_same_dir) const;
		////////////////
		bool ShiftChildLoop();	
		bool IsEndChild() const { return is_end_child; }
		unsigned int GetIdVertex() const;
		unsigned int GetIdVertex_Ahead()  const;
		unsigned int GetIdVertex_Behind() const;
		////////////////
		unsigned int GetIdHalfEdge() const { return m_IdHE; }
		unsigned int GetIdUseLoop() const { return m_IdUL; }
		unsigned int GetIdLoop() const;
		unsigned int GetType() const; // 0:浮遊点 1:浮遊辺 2:面積がある
		unsigned int CountVertex_UseLoop() const;
		bool IsParent() const;
		bool IsSameUseLoop(const CItrLoop& itrl) const{
			return ( itrl.m_IdUL == this->m_IdUL );
		}
		bool IsEdge_BothSideSameLoop() const;
	private:
		bool is_valid;
		bool is_initial;
		bool is_end_child;
		unsigned int m_IdHE;
		unsigned int m_IdUL;
		const CBRepSurface& m_pBRep2D;
	};
  
	//! iterator goes through loops and edges around one vertex
	class CItrVertex : public Cad::IItrVertex
	{
	public:
		CItrVertex(const CBRepSurface& pBRep2D, unsigned int id_v);
		//! 反時計周りに頂点まわりをめぐる
		void operator++();
		void operator++(int n);	//!< ダミーのオペレータ(++と同じ働き)
		//! 頂点周りの辺のIDと、その辺の始点がid_vと一致しているかどうか
		bool GetIdEdge_Ahead( unsigned int& id_e, bool& is_same_dir) const;
		bool GetIdEdge_Behind(unsigned int& id_e, bool& is_same_dir) const;
		//! Get ID of the loop
		unsigned int GetIdLoop() const;
		void Begin();		    
		bool IsEnd() const; //! return true if iterator goes around vertex
		////////////////
    // non virtual hrom here
		unsigned int GetIdHalfEdge()  const { return m_IdHE; }
		unsigned int GetIdUseVertex() const { return m_IdUV; }
		unsigned int CountEdge() const;
		bool IsParent() const;
		bool IsSameUseLoop(const CItrVertex& itrl) const;
	private:
		bool is_valid;
		bool is_initial;
		unsigned int m_IdUV;
		unsigned int m_IdHE;
		const CBRepSurface& m_pBRep2D;
	};
	friend class CItrLoop;
	friend class CItrVertex;
  ////
  class CResConnectVertex{
  public:
    CResConnectVertex() : id_v1(0), id_v2(0), id_l(0), 
    id_e_add(0), id_l_add(0), is_left_l_add(true){}
  public:
    unsigned int id_v1, id_v2, id_l;
    unsigned int id_e_add, id_l_add;
    bool is_left_l_add;
  };  
public:
	bool AssertValid() const;
	bool IsElemID(Cad::CAD_ELEM_TYPE,unsigned int id) const;
	const std::vector<unsigned int> GetAryElemID(Cad::CAD_ELEM_TYPE itype) const;
	bool GetIdLoop_Edge(unsigned int id_e, unsigned int& id_l_l, unsigned int& id_l_r) const;
  unsigned int GetIdLoop_Edge(unsigned int id_e, bool is_left) const;
	bool GetIdVertex_Edge(unsigned int id_e, unsigned int& id_v1, unsigned int& id_v2) const;
	unsigned int GetIdVertex_Edge(unsigned int id_e, bool is_root ) const;
	
	CItrLoop GetItrLoop(unsigned int id_l) const { return CItrLoop(*this,id_l); }
	CItrLoop GetItrLoop_SideEdge(unsigned int id_e, bool is_left) const;
	CItrVertex GetItrVertex(unsigned int id_v) const { return CItrVertex(*this,id_v); }

	////////////////

	void Clear();

	// Add vertex to edge; ret val is id of vertex; ret 0 if it fails
	unsigned int AddVertex_Edge(unsigned int id_e); 
	// 面に頂点を加える関数(失敗したら０を返す)
	unsigned int AddVertex_Loop(unsigned int id_l);
  
	bool RemoveEdge(unsigned int id_e, bool is_del_cp);
	bool RemoveVertex(unsigned int id_v);
  
	bool MakeHole_fromLoop(unsigned int id_l);
  unsigned int SealHole(unsigned int id_e, bool is_left);
  
	// function to make edge from 2 vertices (ID:id_v1,id_v2)
	CResConnectVertex ConnectVertex(const CItrVertex& itrv1, const CItrVertex& itrv2, bool is_id_l_add_left);
	std::vector< std::pair<unsigned int,bool> > GetItrLoop_ConnectVertex(const CItrVertex& itrv1, const CItrVertex& itrv2) const;
	std::vector< std::pair<unsigned int,bool> > GetItrLoop_RemoveEdge(unsigned int id_e) const;

	bool SwapItrLoop(const CItrLoop& itrl, unsigned int id_l_to );
	
	//! save and load file
	bool Serialize( Com::CSerializer& serialize );	
private:
	static unsigned int GetFreeKey(const std::map<unsigned int,unsigned int>& map);
	unsigned int TypeUseLoop(unsigned int id_ul) const;
private:
	CBRep m_BRep;
	std::map<unsigned int, unsigned int> map_l2ul;
	std::map<unsigned int, unsigned int> map_e2he;
};



}


#endif
