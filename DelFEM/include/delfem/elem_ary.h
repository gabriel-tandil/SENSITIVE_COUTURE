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
@brief interface of element array class (Fem::Field::CElemAry)
@author Nobuyuki Umetani
*/

#if defined(__VISUALC__)
    #pragma warning( disable : 4786 )
#endif

#if !defined(ELEM_ARY_H)
#define ELEM_ARY_H

#include <vector>
#include <assert.h>
#include <cstdlib> //(abs)
#include <cstring> //(strspn, strlen, strncmp, strtok)

#include "delfem/objset.h"
#include "delfem/indexed_array.h"

namespace Fem{
namespace Field{

//! type of element segment
enum ELSEG_TYPE{ 
	CORNER=1,	//!< corner node
	EDGE=2,		//!< edge node
	BUBBLE=4	//!< in the element node
};

// this integet index will be used in "elem_ary.cpp". Please don't change without consideration.
// ここだけは別のヘッダファイルに移したほうがいいかも．このIndexだけ欲しいクラス(drawer_field.hとか)のために
//! type of element
enum ELEM_TYPE{	
	ELEM_TYPE_NOT_SET=0,
	POINT=1,	//!< point element
	LINE=2,		//!< line element
	TRI=3, 		//!< triangle element
	QUAD=4, 	//!< quadratic element
	TET=5, 		//!< tetrahedra element
	HEX=6 		//!< hexagonal element
};

/*! 
@brief class of Element Array
@ingroup Fem
*/
class CElemAry
{
	
public:
	//! 要素セグメント．要素の場所(Corner,Bubble..etc)ごとの節点番号を取得するインターフェース関数
	class CElemSeg{
		friend class CElemAry;
	public:
		CElemSeg(unsigned int id_na, ELSEG_TYPE elseg_type)
			: m_id_na(id_na), m_elseg_type(elseg_type){}

		unsigned int GetMaxNoes() const { return max_noes; }	//!< ノード番号の一番大きなものを得る（このnoesを格納するためには一つ大きな配列が必要なので注意）
		unsigned int Length() const { return m_nnoes; }	//!< return the node size per elem seg  ( will be renamed to Length() );
		unsigned int Size() const { return nelem; }	  //!< get the number of elements ( will be renamed to Size() )
		unsigned int GetIdNA() const { return m_id_na; }	    //!< get ID of node array this element segment refers
		ELSEG_TYPE GetElSegType() const { return m_elseg_type; }	//!< 要素セグメントタイプ(Fem::Field::CORNER,Fem::Field::BUBBLE,Fem::Field::EDGE)を得る

		//! get node indexes of element (ielem)
		void GetNodes(const unsigned int& ielem, unsigned int* noes ) const {
			for(unsigned int inoes=0;inoes<m_nnoes;inoes++){ 
				assert( ielem*npoel+begin+inoes < nelem*npoel );
				noes[inoes] = abs(pLnods[ielem*npoel+begin+inoes]); 
			}
		}
		//! 節点番号を設定
		void SetNodes(unsigned int ielem, unsigned int idofes, int ino ){
			pLnods[ielem*npoel+begin+idofes] = ino;
		}
		//! 各場所(Corner,Bubble)に定義されている要素節点の数を出す
		static unsigned GetLength(ELSEG_TYPE elseg_type, ELEM_TYPE elem_type)
		{	// スタティックな関数なので実体は持たなくて良い
			if( elem_type == POINT ){
				if( elseg_type == CORNER ){ return 1; }
				if( elseg_type == EDGE   ){ return 2; }
				if( elseg_type == BUBBLE ){ return 2; }	// 点と点を結ぶLagrange未定乗数を定義するときに使うから０じゃダメ
				assert(0);
			}
			else if( elem_type == LINE ){
				if( elseg_type == CORNER ){ return 2; }
				if( elseg_type == EDGE   ){ return 1; }
				if( elseg_type == BUBBLE ){ return 0; }
				assert(0);
			}
			else if( elem_type == TRI ){
				if( elseg_type == CORNER ){ return 3; }
				if( elseg_type == EDGE   ){ return 3; }
				if( elseg_type == BUBBLE ){ return 1; }
				assert(0);
			}
			else if( elem_type == QUAD ){
				if( elseg_type == CORNER ){ return 4; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 4; }
				assert(0);
			}
			else if( elem_type == TET ){
				if( elseg_type == CORNER ){ return 4; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 6; }
				assert(0);
			}
			else if( elem_type == HEX ){
				if( elseg_type == CORNER ){ return 8; }
				if( elseg_type == BUBBLE ){ return 1; }
				if( elseg_type == EDGE   ){ return 12; }
				assert(0);
			}
			else{ assert(0); }
			return 0;
		}
		//! 要素の次元を取得する
    static unsigned GetElemDim(ELEM_TYPE elem_type){
      if(      elem_type == POINT ){ return 0; }
      else if( elem_type == LINE  ){ return 1; }
      else if( elem_type == TRI   ){ return 2; }
      else if( elem_type == QUAD  ){ return 2; }
      else if( elem_type == TET   ){ return 3; }
      else if( elem_type == HEX   ){ return 3; }
			else{ assert(0); }      
			return 0;
    }
	private: // variable given at construction
		unsigned int m_id_na; // <- 復活、ないと大変なことになった。
		enum ELSEG_TYPE m_elseg_type;
		// int id_es_corner <- そのうち追加予定
	private: // 
		unsigned int max_noes;	// このElementSegmentに属する節点番号の最大のもの
		unsigned int begin;		// m_aLnodsのどこから始まるか (<npoel)
		unsigned int m_nnoes;		// ElementSegmentの長さ
	private: // variable given by CElemAry (CELemAry is friend class)
		mutable unsigned int* pLnods;
		mutable unsigned int npoel;
		mutable unsigned int nelem;
	};
public:
	//! default constructor
	CElemAry(){
		m_nElem = 0; npoel = 0; m_pLnods = 0;
	}
  CElemAry(const CElemAry& ea);
  
  // Constructor with elem number (nelem) and elem type (elemtype)
	CElemAry(unsigned int nelem, ELEM_TYPE elem_type) : m_nElem(nelem), m_ElemType(elem_type){
		npoel = 0; m_pLnods = 0;
	}
	//! destructor
	virtual ~CElemAry(){
		if( this->m_pLnods != 0 ) delete[] m_pLnods;
	}

	bool IsSegID( unsigned int id_es ) const { return m_aSeg.IsObjID(id_es); }	//!< Check if there are elem segment with ID:id_es
	const std::vector<unsigned int>& GetAry_SegID() const { return this->m_aSeg.GetAry_ObjID(); } //!< 要素セグメントIDの配列を得る関数
	unsigned int GetFreeSegID() const{ return m_aSeg.GetFreeObjID(); }	//!< 使われていない要素セグメントIDを取得する関数
  
	const CElemSeg& GetSeg(unsigned int id_es) const{			
		assert( this->m_aSeg.IsObjID(id_es) );
		if( !m_aSeg.IsObjID(id_es) ) throw;
		const CElemSeg& es = m_aSeg.GetObj(id_es);
		es.pLnods = this->m_pLnods;
		es.npoel = this->npoel;
		es.nelem = this->m_nElem;
		return es;
	}
	CElemSeg& GetSeg(unsigned int id_es){
		assert( this->m_aSeg.IsObjID(id_es) );
		if( !m_aSeg.IsObjID(id_es) ) throw;
		CElemSeg& es = m_aSeg.GetObj(id_es);
		es.pLnods = this->m_pLnods;
		es.npoel = this->npoel;
		es.nelem = this->m_nElem;
		return es;
	}

	virtual ELEM_TYPE ElemType() const { return m_ElemType; }	//!< get type of element
	virtual unsigned int Size() const { return m_nElem; }	//!< get number of elements

	//! make CRS data (for the off-diaglnal block)
	virtual bool MakePattern_FEM(	
		const unsigned int& id_es0, const unsigned int& id_es1, 
        Com::CIndexedArray& crs ) const;

	//! make CRS data (for the diaglnal)
	virtual bool MakePattern_FEM(	
		const unsigned int& id_es0, 
        Com::CIndexedArray& crs ) const;

	// make edge ( for 2nd order interporation )
	virtual bool MakeEdge(const unsigned int& id_es_co, unsigned int& nedge, std::vector<unsigned int>& edge_ary) const;
	
	// make map from elem to edge
	virtual bool MakeElemToEdge(const unsigned int& id_es_corner, 
		const unsigned int& nedge, const std::vector<unsigned int>& edge_ary,
		std::vector<int>& el2ed ) const;
	/*!
	@brief 境界要素を作る(可視化のための関数)
	*/
	virtual CElemAry* MakeBoundElemAry(unsigned int id_es_corner, unsigned int& id_es_add, std::vector<unsigned int>& aIndElemFace) const;
	
	//! IO functions
	int InitializeFromFile(const std::string& file_name, long& offset);
	int WriteToFile(       const std::string& file_name, long& offset, unsigned int id) const;

	// lnods should be unsigned int?
	//! Add element segment
	std::vector<int> AddSegment(const std::vector< std::pair<unsigned int,CElemSeg> >& es_ary, const std::vector<int>& lnods );
	int              AddSegment(unsigned int id, const CElemSeg& es,					 const std::vector<int>& lnods );

	// 要素を囲む要素を作る．内部でelsuelのメモリ領域確保はしないので，最初から確保しておく
	bool MakeElemSurElem( const unsigned int& id_es_corner, int* elsuel) const;
private:

	bool MakePointSurElem( const unsigned int id_es, Com::CIndexedArray& elsup ) const;

	bool MakePointSurPoint
  (const unsigned int id_es, const Com::CIndexedArray& elsup, bool isnt_self,
   Com::CIndexedArray& psup ) const;

  bool MakeElemSurElem( const unsigned int& id_es_corner, const Com::CIndexedArray& elsup, int* elsuel) const;
private:


protected:
  unsigned int m_nElem;	//!< number of elements
  ELEM_TYPE m_ElemType;	//!< type of element
	unsigned int npoel;		//!< the total number of nodes including in one elemement
	unsigned int * m_pLnods;			//!< connectivity
	Com::CObjSet<CElemSeg> m_aSeg;	//!< the set of elem segment
};

}
}

#endif
