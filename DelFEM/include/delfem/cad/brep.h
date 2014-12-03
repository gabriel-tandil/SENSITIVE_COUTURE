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
@brief interface of the class (Cad::CBRep) wich represents topology with B-Rep data strcutre
@author Nobuyuki Umetani
*/

#if !defined(B_REP_H)
#define B_REP_H

#ifdef __VISUALC__
	#pragma warning( disable : 4786 )
#endif

#include "delfem/cad/objset_cad.h"

namespace Cad{

//! @addtogroup CAD
//! @{

//! topology loop class
class CUseLoop{
public:
	CUseLoop(const CUseLoop& rhs)
  : id(rhs.id), id_l(rhs.id_l), id_he(rhs.id_he), id_ul_c(rhs.id_ul_c), id_ul_p(rhs.id_ul_p){}
  CUseLoop(const unsigned int id,
           const unsigned int id_he, const unsigned int id_ul_c, const unsigned int id_ul_p)
  : id(id), id_l(0), id_he(id_he), id_ul_c(id_ul_c), id_ul_p(id_ul_p){}
public:
  unsigned int id;    //!< ID
  
	// geometry element ID (brep.cppやbrep2d.cppでは参照されないはず？)
  unsigned int id_l;	//!< 外側のループの場合は０
  
	// topology elemnet ID
  unsigned int id_he;   //!< HalfEdgeのID
  unsigned int id_ul_c;	//!< 子ループ、id_ul_c=0の場合はリストの終わり
  unsigned int id_ul_p;	//!< 親ループ。id_lはid_ul_pを持っている．id_ul_p==idの場合は自分が親，id_ul_p==0の場合は外側のループ
};

//! topoloty HalfEdge class
class CHalfEdge{
public:
	CHalfEdge(const CHalfEdge& rhs)
		: id(rhs.id),  
		id_uv(rhs.id_uv),
		id_he_f(rhs.id_he_f), id_he_b(rhs.id_he_b), id_he_o(rhs.id_he_o),
		id_ul(rhs.id_ul),
		id_e(rhs.id_e), is_same_dir(rhs.is_same_dir){}
  CHalfEdge(const unsigned int id,
            const unsigned int id_uv,
            const unsigned int id_he_f, const unsigned int id_he_b, const unsigned int id_he_o,
            const unsigned int id_ul )
		: id(id), 
		id_uv(id_uv), 
		id_he_f(id_he_f), id_he_b(id_he_b), id_he_o(id_he_o), 
		id_ul(id_ul), 
		id_e(0), is_same_dir(true){}
public:
  unsigned int id;        //!< ID
  unsigned int id_uv;     //!< ID of UseVertex
  unsigned int id_he_f;   //!< ID of previous HalfEdge
  unsigned int id_he_b;   //!< ID of later HalfEdge
  unsigned int id_he_o;	  //!< ID of opposite side HalfEdge (if UV is floating point, this become self )
  unsigned int id_ul;     //!< ID of UseLoop
  
  ////
  unsigned int id_e;      //!< id of edge geometry. 0 if uv is a floating vertex
  bool is_same_dir;       //!< is the geometry edge goes same direction as topology edge
};

//! 位相頂点クラス
class CUseVertex{
public:
  CUseVertex(const unsigned int id, const unsigned int id_he)
		: id(id), id_he(id_he), id_v(0){}	
	CUseVertex(const CUseVertex& rhs)
		: id(rhs.id), id_he(rhs.id_he), id_v(rhs.id_v){}	
public:
  unsigned int id;    //!< ID
  unsigned int id_he; //!< HalfEdgeのID
	//ここからは幾何要素ID
  unsigned int id_v;  //!< 幾何頂点のID
};

//! B-repを用いた位相情報格納クラス
class CBRep{
public:
	//! 全ての要素を削除して，初期状態に戻る
	void Clear();

	bool IsID_UseLoop(unsigned int id_ul) const { return m_UseLoopSet.IsObjID(id_ul); }
	bool IsID_HalfEdge(unsigned int id_he) const { return m_HalfEdgeSet.IsObjID(id_he); }
	bool IsID_UseVertex(unsigned int id_uv) const { return m_UseVertexSet.IsObjID(id_uv); }
	std::vector<unsigned int> GetAry_UseVertexID() const { return m_UseVertexSet.GetAry_ObjID(); }

	const CUseLoop& GetUseLoop(unsigned int id_ul) const;
	const CUseVertex& GetUseVertex(unsigned int id_uv) const;
	const CHalfEdge& GetHalfEdge(unsigned int id_he) const;

	bool SetLoopIDtoUseLoop(unsigned int id_ul, unsigned int id_l);
	bool SetVertexIDtoUseVertex(unsigned int id_uv, unsigned int id_v);
	bool SetEdgeIDtoHalfEdge(unsigned int id_he, unsigned int id_e, bool is_same_dir);
	
	int AssertValid_Use() const;
	std::vector<unsigned int> FindHalfEdge_Edge(const unsigned int& id_e) const;
	std::vector<unsigned int> FindHalfEdge_Vertex(const unsigned int& id_v) const;
	
	////////////////////////////////
	// オイラー操作

  //! 点と点を繋げてループとエッジを作る
	bool MEVVL(unsigned int& id_he_add1, unsigned int& id_he_add2,
             unsigned int& id_uv_add1, unsigned int& id_uv_add2, unsigned int& id_ul_add );
	//! ループを２つに分ける
	bool MEL(unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_ul_add,
           const unsigned int id_he1, const unsigned int id_he2);
	/*!
	@brief 辺を消去して２つのループを１つにする 半辺he_remの反対側の半辺が接する半ループを消去する
	*/
	bool KEL(const unsigned int id_he_rem1);
	//! id_heの起点から、線分を伸ばす
	bool MEV(unsigned int& id_he_add1,unsigned int& id_he_add2, unsigned int& id_uv_add, 
           const unsigned int id_he);
	//! id_heの起点を消去して２つの辺を１つにする。
	bool KVE( unsigned int& id_he_rem1 );
	/*!
	@brief id_heを２つに分ける
	入力ではhe2はheに向かいあう半辺とすると
	出力ではheの前にhe_add1,heの向かいにhe_add2,he_add2の後ろにhe2
	*/
	bool MVE( unsigned int& id_he_add1, unsigned int& id_he_add2, unsigned int& id_uv_add, 
           const unsigned int id_he);
	/*!
	@brief he1の起点uv1とhe2の起点uv2を結んで、２つのループをつなげる
	he1は[uv1-uv2]、he2は[uv2-uv1]
	*/
	bool MEKL(unsigned int& id_he_add1, unsigned int& id_he_add2,  
            const unsigned int id_he1, const unsigned int id_he2 );
	//! ループをつなげる
  bool KEML(unsigned int& id_ul_add1,
            const unsigned int& id_he1 );
	/*!
	@brief ループと浮遊点をつなげる,he1がLoop上のEdgeでhe2が浮遊点Edge
	he2は[uv2-uv1], he_add1は[uv1-uv2]のHalfEdgeとなる
	*/
	bool MEKL_OneFloatingVertex(unsigned int& id_he_add1,
                              const unsigned int id_he1, const unsigned int id_he2 );
	/*!
	@brief ループと浮遊点をつなげる,he1,he2が浮遊点Edge
	he1は[uv1-uv2],he2は[uv2-uv1]のHalfEdgeとなる
	*/
	bool MEKL_TwoFloatingVertex( const unsigned int id_he1, const unsigned int id_he2 );
	/*! 
	@brief 片方が端点であるEdgeを削除する。
	he1の起点uv1は他の辺につながっていない端点である
	*/
	bool KEML_OneFloatingVertex(unsigned int& id_ul_add, 
                              const unsigned int id_he1);
	//! 両方が端点であるEdgeを削除する。
	bool KEML_TwoFloatingVertex(unsigned int& id_ul_add,
                              const unsigned int id_he1);
	//! 浮遊点を作る
	bool MVEL(unsigned int& id_uv_add, unsigned int& id_he_add, unsigned int& id_ul_add, 
            const unsigned int id_ul_p);
	//! 浮遊点を消す
	bool KVEL(const unsigned int id_uv_rem);

	//! 位相ループを入れ替える
	bool SwapUseLoop(unsigned int id_ul1, unsigned int id_ul2);
	//! 位相ループを動かす、id_ul1のループをid_ul2の親の子にする。
	bool MoveUseLoop(unsigned int id_ul1, unsigned int id_ul2);
private:
	bool SwapUseLoop_CP_SameLoop(unsigned int id_ul1,unsigned int id_ul2);
	bool SwapUseLoop_CP_DifferentLoop(unsigned int id_ul1,unsigned int id_ul2);

//private:	// そのうちprivateにする予定(残りはSerializ)
public:
  CObjSetCad<CUseLoop>   m_UseLoopSet;   //!< UseLoopのIDで管理された集合
  CObjSetCad<CHalfEdge>  m_HalfEdgeSet;  //!< HalfEdgeのIDで管理された集合
  CObjSetCad<CUseVertex> m_UseVertexSet; //!< UseVertexのIDで管理された集合
};

//! @}

}

#endif

