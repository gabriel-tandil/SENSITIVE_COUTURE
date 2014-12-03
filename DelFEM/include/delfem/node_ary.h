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
@brief the interface of node array class (Fem::Field::CNodeAry)
@author Nobuyuki Umetani
*/

#if !defined(NODE_ARY_H)
#define NODE_ARY_H

#if defined(__VISUALC__)
    #pragma warning( disable : 4786 )
#endif

#ifndef for 
#define for if(0); else for
#endif

#include <iostream>
#include <cassert>

#include "delfem/complex.h"
#include "delfem/objset.h"
#include "delfem/elem_ary.h"	// remove this dependency in future

// the temp class declearation
namespace MatVec{
	class CVector_Blk;
	class CZVector_Blk;
  class CBCFlag;
}

namespace Fem{
namespace Field
{

/*! 
@brief class which contains nodes value (coordinte,displacement,temparature....etc )
@ingroup Fem
*/
class CNodeAry
{
public:
	//! Class for Node Segment (Holding Data of Node)
	class CNodeSeg{
		friend class CNodeAry;
	public:
		CNodeSeg(const unsigned int& len, const std::string& name)
			: len(len), name(name){}
		unsigned int Length() const { return len; }	//!< The length of value
		unsigned int Size() const { return nnode; }	//!< The number of nodes
		inline void GetValue(unsigned int inode, double* aVal ) const	//!< get value from node
		{
			for(unsigned int i=0;i<len;i++){
				aVal[i] = paValue[inode*DofSize+idofval_begin+i];
			}
		}
		void GetValue(unsigned int inode, Com::Complex* aVal ) const	//!< get complex value from node
		{
			const unsigned int n = len/2;
			for(unsigned int i=0;i<n;i++){
				double dr = paValue[inode*DofSize+idofval_begin+i*2];
				double di = paValue[inode*DofSize+idofval_begin+i*2+1];
				aVal[i] = Com::Complex(dr,di);
			}
		}
		inline void SetValue(unsigned int inode, unsigned int idofns, double val )	//!< set value to node 
		{
			paValue[inode*DofSize+idofval_begin+idofns] = val;
		}
		inline void AddValue(unsigned int inode, unsigned int idofns, double val )	//!< add value to node
		{
			paValue[inode*DofSize+idofval_begin+idofns] += val;
		}
		void SetZero()	//!< set zero to all value
		{
			for(unsigned int ino=0;ino<nnode;ino++){
			for(unsigned int ilen=0;ilen<len;ilen++){
				paValue[ino*DofSize+idofval_begin+ilen] = 0;
			}
			}
		}
	private:
    unsigned int len;	//!< the size of value
    std::string name;	//!< name
	private: // not need when initialize 
		unsigned int idofval_begin;	//!< offset of value
	private: // the variables given by CNodeAry
		mutable double* paValue;	//!< value list
		mutable unsigned int DofSize;	//!< the width of node array
		mutable unsigned int nnode;	//!< number of nodes
	};

private:
	class CEaEsInc
	{
	public:
		unsigned int id_ea;
		unsigned int id_es;
		std::vector<unsigned int> aIndEaEs_Include;	//!< このEsが含んでるiidEaEs
	};

public:
	CNodeAry(const unsigned int size);
	CNodeAry();
  CNodeAry(const CNodeAry& na);
	virtual ~CNodeAry();

	bool ClearSegment();

	////////////////////////////////
	// Get Methods

//	const std::string& Name() const { return m_str_name; }  //!< name
	unsigned int Size()	const { return m_Size; }  //!< number of nodes

	//! Get a not-used ID of Node Segment
	unsigned int GetFreeSegID() const { return m_aSeg.GetFreeObjID(); }
	//! Get not-used num IDs of Node Segment
	std::vector<unsigned int> GetFreeSegID(unsigned int num) const { return m_aSeg.GetFreeObjID(num); }
	//! Check if id_ns is whether ID of Node Segment or not
	bool IsSegID( unsigned int id_ns ) const{ return m_aSeg.IsObjID(id_ns); }
	const std::vector<unsigned int>& GetAry_SegID() const { return this->m_aSeg.GetAry_ObjID(); }
	//! Get Node Segment (const)
	const CNodeSeg& GetSeg(unsigned int id_ns) const{		
		assert( this->m_aSeg.IsObjID(id_ns) );
		if( !m_aSeg.IsObjID(id_ns) ) throw;
		const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
		assert( m_paValue != 0 );
		ns.paValue = m_paValue;
		ns.DofSize = m_DofSize;
		ns.nnode = m_Size;
		return ns;
	}
	//! Get Node Segment 
	CNodeSeg& GetSeg(unsigned int id_ns){			
		assert( this->m_aSeg.IsObjID(id_ns) );
		if( !m_aSeg.IsObjID(id_ns) ) throw;
		CNodeSeg& ns = m_aSeg.GetObj(id_ns);
		assert( m_paValue != 0 );
		ns.paValue = m_paValue;
		ns.DofSize = m_DofSize;
		ns.nnode = m_Size;
		return ns;
	}

	//! 節点セグメントの値をvecに代入する
	bool GetValueFromNodeSegment(unsigned int id_ns, MatVec::CVector_Blk& vec, unsigned int ioffset=0) const; 
	bool AddValueFromNodeSegment(double alpha, unsigned int id_ns, MatVec::CVector_Blk& vec, unsigned int ioffset=0) const;

	////////////////
	// 参照要素セグメント追加メソッド

	void AddEaEs( std::pair<unsigned int, unsigned int> eaes );
	std::vector< std::pair<unsigned int, unsigned int> > GetAryEaEs() const;
	void SetIncludeEaEs_InEaEs( std::pair<unsigned int, unsigned int> eaes_included,
                             std::pair<unsigned int, unsigned int> eaes_container );
	bool IsIncludeEaEs_InEaEs( std::pair<unsigned int, unsigned int> eaes_inc,
                            std::pair<unsigned int, unsigned int> eaes_in ) const;
	std::vector< std::pair<unsigned int, unsigned int> > GetAry_EaEs_Min() const;
	unsigned int IsContainEa_InEaEs(std::pair<unsigned int, unsigned int>eaes, unsigned int id_ea) const;


	////////////////
	// 変更・削除・追加メソッド

	const std::vector<int> AddSegment( const std::vector< std::pair<unsigned int,CNodeSeg> >& seg_vec ); //!< 初期化しない
	const std::vector<int> AddSegment( const std::vector< std::pair<unsigned int,CNodeSeg> >& seg_vec, const double& val);  //!< 値valで初期化
	const std::vector<int> AddSegment( const std::vector< std::pair<unsigned int,CNodeSeg> >& seg_vec, const std::vector<double>& val_vec);  //!< ベクトルvalで初期化

	////////////////
	// 値の変更メソッド

	//! set value of vector to the node segment
	bool SetValueToNodeSegment(unsigned int id_ns, const MatVec::CVector_Blk& vec, unsigned int ioffset=0); // Segmentにvecをセットする

	//! 境界条件設定に使われる。
	bool SetValueToNodeSegment(const Field::CElemAry& ea, const unsigned int id_es, const unsigned int id_ns, const unsigned int idofns, const double val){
		assert( m_aSeg.IsObjID(id_ns) );
		if( !m_aSeg.IsObjID(id_ns) ) return false;
		const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
		assert( idofns < ns.len );
		assert( ea.IsSegID(id_es) );
		const CElemAry::CElemSeg& es = ea.GetSeg(id_es);
		unsigned int noes[256]; // ARHHHHH!!!!
		unsigned int nnoes = es.Length();
		for(unsigned int ielem=0;ielem<ea.Size();ielem++){
			es.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				const unsigned int inode0 = noes[inoes];
				m_paValue[inode0*m_DofSize+ns.idofval_begin+idofns] = val;
			}
		}
		return true;
	}

	
	bool AddValueToNodeSegment(unsigned int id_ns, const MatVec::CVector_Blk& vec, double alpha,  unsigned int ioffset = 0); //!< Segmentにalpha倍されたvecを加える
	bool AddValueToNodeSegment(unsigned int id_ns, const MatVec::CZVector_Blk& vec, double alpha); //!< Segmentにalpha倍されたvecを加える
	bool AddValueToNodeSegment(unsigned int id_ns_to, unsigned int id_ns_from, double alpha );	//!< ns_toへalpha倍されたns_fromを加える


	//! load from file
	int InitializeFromFile(const std::string& file_name, long& offset);
	//! write to file
	int WriteToFile(const std::string& file_name, long& offset, unsigned int id ) const;
	int DumpToFile_UpdatedValue(const std::string& file_name, long& offset, unsigned int id ) const;

private:
	unsigned int GetIndEaEs( std::pair<unsigned int, unsigned int> eaes ) const
	{
		unsigned int ieaes;
		for(ieaes=0;ieaes<m_aEaEs.size();ieaes++){
			if( m_aEaEs[ieaes].id_ea == eaes.first && m_aEaEs[ieaes].id_es == eaes.second ) break;
		}
		return ieaes;
	}
private:
//	std::string m_str_name;	//!< name
	unsigned int m_Size;	//!< number of nodes
	unsigned int m_DofSize;		//!< the size of DOF in node  	
	double* m_paValue;		//!< the values in nodes  
	Com::CObjSet<CNodeSeg> m_aSeg;	//!< the array of node segment
	std::vector< CEaEsInc > m_aEaEs;	//!< whitch element segments this node is included
};

}	// end namespace field
}	// end namespace Fem

#endif // !defined(NODE_ARY_H)
