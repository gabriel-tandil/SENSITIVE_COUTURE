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
@brief CADのためのID管理テンプレートクラス(Cad::CObjSet)の宣言・実装
@author Nobuyuki Umetani
*/

#if !defined(OBJ_SET_CAD_H)
#define OBJ_SET_CAD_H

#ifdef __VISUALC__
	#pragma warning( disable : 4786 )
#endif

#include <vector>

namespace Cad{

//! ID管理テンプレート
template <class T>
class CObjSetCad
{		
public:
	void clear(){
		m_aIndex2ID.clear();
		m_aID2Index.clear();
		m_aObj.clear();
	}
    //! 要素の実体参照(const)を得る
    const T& GetObj(unsigned int id_e) const {
		const int ie = this->GetAryInd(id_e);
		assert( ie >= 0 && ie < (int)m_aObj.size() );
		return m_aObj[ie];
	}
    //! 要素の実体参照を得る
    T& GetObj(unsigned int id_e){
		const int ie = this->GetAryInd(id_e);
		assert( ie >= 0 && ie < (int)m_aObj.size() );
		return m_aObj[ie];
	}
	bool DeleteObj(unsigned int id_e){
		const int ie = this->GetAryInd(id_e);
		assert( ie >= 0 && ie < (int)m_aObj.size() );
		m_aObj.erase( m_aObj.begin()+ie );
		assert( ie >= 0 && ie < (int)m_aIndex2ID.size() );
		m_aIndex2ID.erase( m_aIndex2ID.begin()+ie );
		for(unsigned int i=0;i<m_aID2Index.size();i++){ m_aID2Index[i] = -1; }
		for(unsigned int iie=0;iie<m_aIndex2ID.size();iie++){
			unsigned int id_e = m_aIndex2ID[iie];
			assert( id_e<m_aID2Index.size() );
			m_aID2Index[id_e] = iie;
		}
		return true;
	}
	int AddObj(const T& in_obj){
		unsigned int id1 = AddID(in_obj.id);
		assert( IsObjID(id1) );
		T tmp_obj = in_obj;
		tmp_obj.id = id1;
		this->m_aObj.push_back( tmp_obj );
		assert( m_aIndex2ID.size() == m_aObj.size() );
		return id1;
	}
	bool IsObjID( unsigned int id_e ) const {
		const int ie1 = this->GetAryInd(id_e);
		if( ie1 != -1 ) return true;
		return false;
	}
	const std::vector<unsigned int>& GetAry_ObjID() const { return this->m_aIndex2ID; }
	unsigned int GetFreeObjID() const {
		if( m_aID2Index.size() == 0 )	return 1;
		for(unsigned int iid=1;iid<m_aID2Index.size();iid++){
			if( m_aID2Index[iid] == -1 )	return iid;
		}
		return m_aID2Index.size();
	}
	std::vector<unsigned int> GetFreeObjID(const unsigned int size) const {
		std::vector<unsigned int> res;
		res.resize(size);
		unsigned int isize = 0;
		assert( m_aID2Index.size() != 1 );
		if( m_aID2Index.size() == 0 ){
			for(;;){
				res[isize] = isize+1;
				isize++;
				if(isize==size) return res;
			}
		}
		for(unsigned int iid=1;iid<m_aID2Index.size();iid++){
			if( m_aID2Index[iid] == -1 ){
				res[isize] = iid;
				isize++;
				if(isize==size) return res;
			}
		}
		for(unsigned int i=0;;i++){
			res[isize] = m_aID2Index.size()+i;
			isize++;
			if(isize==size) return res;
		}
		return res;
	}
	void Reserve(unsigned int size){
		m_aIndex2ID.reserve(size);
		m_aID2Index.reserve(size+1);
		m_aObj.reserve(size);
	}
protected:
	unsigned int AddID(const int& tmp_id){
		unsigned int id1;
		if( !IsObjID( tmp_id ) && tmp_id>0 && tmp_id<255 ){
			id1 = tmp_id;
			if( m_aID2Index.size() <= id1 ){
				m_aID2Index.resize( id1+1, -1 );
				m_aID2Index[id1] = m_aIndex2ID.size();
			}
			m_aIndex2ID.push_back(id1);
			m_aID2Index[id1] = m_aIndex2ID.size()-1;
		}
		else{
			id1 = this->GetFreeObjID();
			assert( !IsObjID(id1) );
			if( m_aID2Index.size() <= id1 ){
				m_aID2Index.resize( id1+1, -1 );
				m_aID2Index[id1] = m_aIndex2ID.size();
			}
			m_aIndex2ID.push_back(id1);
			m_aID2Index[id1] = m_aIndex2ID.size()-1;
		}
		assert( IsObjID(id1) );
		return id1;
	}
    //! IDからIndexを得る．もしこのIDを持つ要素が存在しなければ-1を返す
	int GetAryInd(unsigned int id_e) const {
		int ie1;
		if( m_aID2Index.size() <= id_e ) ie1 = -1;
		else ie1 = m_aID2Index[id_e];
		return ie1;
	}
protected:
	std::vector<unsigned int> m_aIndex2ID;
    std::vector<int> m_aID2Index;   // なければ-1
	std::vector<T> m_aObj;
};

} // end namespace Cad

#endif
