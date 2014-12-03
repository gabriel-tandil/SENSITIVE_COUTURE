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
// Ordering_Blk.cpp : オーダリングクラス(COrdering_Blk)の実装
////////////////////////////////////////////////////////////////


#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#ifndef for 
#define for if(0); else for
#endif

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <iostream>

#include "delfem/matvec/ordering_blk.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/matdia_blkcrs.h"

using namespace MatVec;

void COrdering_Blk::SetOrdering(const std::vector<int>& ord)
{
  m_nblk = ord.size();
  if( this->m_pOrder != 0 ) delete[] m_pOrder;
  if( this->m_pInvOrder != 0 ) delete[] m_pInvOrder;  
  this->m_pOrder = new int [m_nblk];
  this->m_pInvOrder = new int [m_nblk];
  for(unsigned int iblk=0;iblk<m_nblk;iblk++){ m_pOrder[iblk] = ord[iblk]; }
  for(unsigned int iblk=0;iblk<m_nblk;iblk++){ m_pInvOrder[iblk] = -1;  }    
  for(unsigned int iblk=0;iblk<m_nblk;iblk++){
    unsigned int i1 = m_pOrder[iblk];
    assert( m_pInvOrder[i1] == -1 );
    m_pInvOrder[i1] = iblk;
  }
}

void COrdering_Blk::MakeOrdering_RCM(const MatVec::CMatDia_BlkCrs& mat)
{
	const unsigned int nblk = mat.NBlkMatCol();
	this->m_nblk = nblk;
	if( this->m_pOrder != 0 ) delete[] m_pOrder;
	this->m_pOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ m_pOrder[iblk] = -1; }

	std::vector<int> aLevel;
	aLevel.resize(nblk,-1);

	std::vector<unsigned int> aNext;
	std::multimap<unsigned int,unsigned int> nxt_map;
	aNext.reserve(nblk);
	unsigned int cur_level = 0;
	aLevel[0] = 0;
	aNext.push_back(0);
	unsigned int iblk_cur = 0;
	for(;;){
		if( aNext.empty() ){
			if( iblk_cur != nblk ){
				std::cout << "Error!-->未実装(グラフに島がある場合)" << std::endl;
				assert(0);
			}
			break;
		}
		nxt_map.clear();
		for(unsigned int inxt=0;inxt<aNext.size();inxt++){
			unsigned int iblk0 = aNext[inxt];
			unsigned int ndeg;
			mat.GetPtrIndPSuP(iblk0,ndeg);
			nxt_map.insert( std::make_pair(ndeg,iblk0) );
      //			nxt_map.insert( std::make_pair(100-ndeg,iblk0) );
      //			nxt_map.insert( std::make_pair(inxt,iblk0) );
		}
		aNext.clear();
		std::multimap<unsigned int,unsigned int>::iterator itr;
		for(itr=nxt_map.begin();itr!=nxt_map.end();itr++){
			unsigned int iblk0 = itr->second;
			{	// iblk0を登録
        assert( aLevel[iblk0] == (int)cur_level );
				m_pOrder[iblk0] = iblk_cur;
				iblk_cur++;
			}
			{	// iblk0の周りの次のlevelを検索
				unsigned int npsup;
				const unsigned int* psup = mat.GetPtrIndPSuP(iblk0,npsup);
				assert( npsup == itr->first );
				for(unsigned int ipsup=0;ipsup<npsup;ipsup++){
					const unsigned int jblk0 = psup[ipsup];
					if( aLevel[jblk0] != -1 ){
            assert( aLevel[jblk0] == (int)cur_level
                   || aLevel[jblk0] == (int)cur_level+1
                   || aLevel[jblk0] == (int)cur_level-1 );
						continue;
					}
					aLevel[jblk0] = cur_level+1;
					aNext.push_back(jblk0);
				}
			}
		}
		cur_level++;    
	}
	for(unsigned int iblk=0;iblk<nblk;iblk++){ 
		m_pOrder[iblk] = nblk-1 - m_pOrder[iblk]; 
	}
	m_pInvOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ 
		const unsigned int iblk0 = m_pOrder[iblk];
		m_pInvOrder[iblk0] = iblk;
	}
}


void MatVec::COrdering_Blk::OrderingVector_NewToOld(CVector_Blk& vec_to, const CVector_Blk& vec_from){
	const unsigned int nblk = this->m_nblk;
	assert( vec_to.NBlk() == nblk );
	assert( vec_from.NBlk() == nblk );
  assert( vec_to.Len() >= 0 );
	const unsigned int nlen = vec_to.Len();
  assert( vec_from.Len() == vec_to.Len() );
	for(unsigned int iblk=0;iblk<nblk;iblk++){	// old
		const unsigned int jblk0 = this->OldToNew(iblk);	// new 
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			const double val = vec_from.GetValue(jblk0,ilen);
			vec_to.SetValue(iblk,ilen,val);
		}
	}
}

void COrdering_Blk::OrderingVector_OldToNew(CVector_Blk& vec_to, const CVector_Blk& vec_from){
	const unsigned int nblk = this->m_nblk;
	assert( vec_to.NBlk() == nblk );
	assert( vec_from.NBlk() == nblk );
  assert( vec_to.Len() >= 0 );
	const unsigned int nlen = vec_to.Len();
  assert( vec_from.Len() == vec_to.Len() );
	for(unsigned int iblk=0;iblk<nblk;iblk++){	// new
		const unsigned int jblk0 = this->NewToOld(iblk);	// old
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			const double val = vec_from.GetValue(jblk0,ilen);
			vec_to.SetValue(iblk,ilen,val);
		}
	}
}


void MatVec::COrdering_Blk::MakeOrdering_RCM2(const CMatDia_BlkCrs& mat)
{
	const unsigned int nblk = mat.NBlkMatCol();
	this->m_nblk = nblk;
	if( this->m_pOrder != 0 ) delete[] m_pOrder;
	this->m_pOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ m_pOrder[iblk] = -1; }

	std::vector<int> aLevel;
	aLevel.resize(nblk,-1);

	std::vector<unsigned int> aNextCand;
	std::vector<unsigned int> aNext;
	aNext.reserve(nblk);
	aNextCand.reserve(nblk);
	unsigned int cur_level = 0;
	unsigned int iblk_cur = 0;
	// level0のノードを見つける
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		unsigned int npsup;
		const unsigned int* psup = mat.GetPtrIndPSuP(iblk,npsup);
		unsigned int ipsup;
		for(ipsup=0;ipsup<npsup;ipsup++){
			const unsigned int jblk0 = psup[ipsup];
			if( jblk0 < iblk ) break;
		}
		if( ipsup < npsup ) continue;
		aLevel[iblk] = cur_level;
		m_pOrder[iblk_cur] = iblk;
		iblk_cur++;
		for(ipsup=0;ipsup<npsup;ipsup++){
			const unsigned int jblk0 = psup[ipsup];
      if( aLevel[jblk0] == (int)cur_level || aLevel[jblk0] == (int)cur_level+1 ) continue;
			assert( jblk0 > iblk );
			aNextCand.push_back(jblk0);
			aLevel[jblk0] = cur_level+1;
		}
	}
	cur_level++;
	for(;;){
		if( aNextCand.empty() ){
			if( iblk_cur != nblk ){
				std::cout << "Error!-->not implemented(there is isolated island in the graph)" << std::endl;
				assert(0);
			}
			break;
		}
		aNext.clear();
		for(unsigned int inxt=0;inxt<aNextCand.size();inxt++){
			unsigned int iblk0 = aNextCand[inxt];
      assert( aLevel[iblk0] == (int)cur_level );
			unsigned int npsup;
			const unsigned int* psup = mat.GetPtrIndPSuP(iblk0,npsup);
			unsigned int ipsup;
			for(ipsup=0;ipsup<npsup;ipsup++){
				const unsigned int jblk0 = psup[ipsup];
				if( jblk0 < iblk0 ){
					const int jlev0 = aLevel[jblk0];
          if( jlev0 == (int)cur_level || jlev0 == -1 ){ break; }
				}
			}
			if( ipsup != npsup ){
				aLevel[iblk0] = -1;
				continue;
			}
			aNext.push_back(iblk0);
		}
		aNextCand.clear();
		for(unsigned int inxt=0;inxt<aNext.size();inxt++){
			unsigned int iblk0 = aNext[inxt];
			{	// iblk0を登録
        assert( aLevel[iblk0] == (int)cur_level );
				m_pOrder[iblk_cur] = iblk0;
				iblk_cur++;
			}
			{	// iblk0の周りの次のlevelを検索
				unsigned int npsup;
				const unsigned int* psup = mat.GetPtrIndPSuP(iblk0,npsup);
				for(unsigned int ipsup=0;ipsup<npsup;ipsup++){
					const unsigned int jblk0 = psup[ipsup];
					if( aLevel[jblk0] != -1 || jblk0 < iblk0 ) continue;
					aLevel[jblk0] = cur_level+1;
					aNextCand.push_back(jblk0);
				}
			}
		}
		cur_level++;
	}
	m_pInvOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ 
		const unsigned int iblk0 = m_pOrder[iblk];
		m_pInvOrder[iblk0] = iblk;
	}
}

unsigned AroundEnode(unsigned int iblk_begin, 
	const unsigned int nblk, const std::vector<unsigned int>& ColInd, const std::vector<unsigned int>& RowPtr,
	const std::vector<int>& aParent,
	std::vector<unsigned int>& aStackBlkAroundEnode, std::vector<unsigned int>& aStackEnode )
{	
	aStackBlkAroundEnode.clear();
	aStackEnode.clear();
	aStackBlkAroundEnode.reserve(nblk);
	aStackEnode.reserve(nblk);
	std::vector<int> aTmpFlg;
	aTmpFlg.resize(nblk,-1);
	if( aParent[iblk_begin] == -2 ) return 0;	// Not Enode
    assert( aParent[iblk_begin] >= 0 );
    const unsigned int iblk_per = (unsigned int)aParent[iblk_begin];
    assert( aParent[iblk_per] == (int)iblk_per );

	std::vector<unsigned int> aNext;
	aNext.reserve(nblk);
	aNext.push_back(iblk_begin);
	aTmpFlg[iblk_begin] = iblk_per;

	for(;;){
		if( aNext.empty() ) break;
		const unsigned int iblk_cur = aNext[ aNext.size()-1 ];
		aNext.resize( aNext.size()-1 );
		for(unsigned int ipsup=ColInd[iblk_cur];ipsup<ColInd[iblk_cur+1];ipsup++){
			const unsigned int jblk0 = RowPtr[ipsup];
			if( aParent[jblk0] == -2 ){
                if( aTmpFlg[jblk0] == (int)iblk_per ) continue;
				aStackBlkAroundEnode.push_back(jblk0);
				aTmpFlg[jblk0] = iblk_per;
			}
			else{
                if( aTmpFlg[jblk0] == (int)iblk_per ) continue;
				aTmpFlg[jblk0] = iblk_per;
				aNext.push_back(jblk0);
				aStackEnode.push_back(jblk0);	
			}
		}
	}
	/*
	std::cout << "node around " << iblk_begin << "     parent " << iblk_per << std::endl;
	for(unsigned int i=0;i<aStackBlkAroundEnode.size();i++){
		std::cout << aStackBlkAroundEnode[i] << " ";
	}
	std::cout << std::endl;
	for(unsigned int i=0;i<aStackEnode.size();i++){
		std::cout << aStackEnode[i] << " ";
	}
	std::cout << std::endl;
	getchar();
	*/
	return 0;
}


void MatVec::COrdering_Blk::MakeOrdering_AMD(const MatVec::CMatDia_BlkCrs& mat)
{
	const unsigned int nblk = mat.NBlkMatCol();
	this->m_nblk = nblk;
	if( this->m_pOrder != 0 ) delete[] m_pOrder;
	this->m_pOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ m_pOrder[iblk] = -1; }

	std::vector<unsigned int> aDegree;
	aDegree.resize(nblk);
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		unsigned int npsup;
		mat.GetPtrIndPSuP(iblk,npsup);
		aDegree[iblk] = npsup;
	}

	std::vector<unsigned int> ColInd;
	std::vector<unsigned int> RowPtr;
	{
		ColInd.resize(nblk+1);
		RowPtr.reserve( mat.NCrs() );
		ColInd.push_back(0);
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			unsigned int npsup;
			const unsigned int* psup = mat.GetPtrIndPSuP(iblk,npsup);
			ColInd[iblk+1] += ColInd[iblk]+npsup;
//			std::cout << iblk << " --> ";
			for(unsigned int ipsup=0;ipsup<npsup;ipsup++){
//				std::cout << psup[ipsup] << " ";
				RowPtr.push_back( psup[ipsup] );
			}
//			std::cout << std::endl;
		}
	}

	std::vector<int> aBlk2BlkNext(nblk,-1);
	std::vector<int> aDeg2Blk(nblk,-1);
  for(unsigned int iblk=0;iblk<nblk;iblk++){
		const unsigned int d = aDegree[iblk];
		aBlk2BlkNext[iblk] = aDeg2Blk[d];
		aDeg2Blk[d] = iblk;
	}

	std::vector<int> aEnodeParentBlk(nblk,-2);	// -2はEnodeでないことを表す
/*
	// Degreeの順番に表示
	for(unsigned int ideg=0;ideg<nblk;ideg++){
		if( aDeg2Blk[ideg] == -1 ) continue;
		int iblk0 = aDeg2Blk[ideg];
		std::cout << "Degree : " << ideg << std::endl;
		for(;iblk0!=-1;iblk0=aBlk2BlkNext[iblk0]){
			std::cout << iblk0 << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
*/
	unsigned int iblk_cur = 0;
	for(;;){
		////////////////////////////////
		// Choose Minimum Degree
		unsigned int min_deg, iblk_del;
/*		{
			for(min_deg=0;min_deg<nblk;min_deg++){
				if( aDeg2Blk[min_deg] != -1 ) break;
			}
			if( min_deg == nblk ) break;
			iblk_del = aDeg2Blk[min_deg];
		}*/
		{
			min_deg = nblk;
			iblk_del = nblk;
			for(unsigned int iblk=0;iblk<nblk;iblk++){
//				std::cout << iblk << " " << aDegree[iblk] << std::endl;
				if( aEnodeParentBlk[iblk] == -2 && aDegree[iblk] < min_deg ){
					min_deg = aDegree[iblk];
					iblk_del = iblk;
				}
			}
		}
		if( iblk_del == nblk ) break;
//		std::cout << "iblk_del:" << iblk_del << "   Deg:" << min_deg << std::endl;
		aDeg2Blk[min_deg] = aBlk2BlkNext[iblk_del];
		m_pOrder[iblk_cur] = iblk_del;
		iblk_cur++;
		////////////////////////////////
		// Degree Update 
		int iblk_parent = iblk_del;
		for(unsigned int ipsup=ColInd[iblk_del];ipsup<ColInd[iblk_del+1];ipsup++){
			const unsigned int jblk0 = RowPtr[ipsup];
			if( aEnodeParentBlk[jblk0] != -2 ){	// これはEnode
				iblk_parent = aEnodeParentBlk[jblk0];
				break;
			}
		}
		aEnodeParentBlk[iblk_del] = iblk_parent;
		std::vector<unsigned int> aAroundEnode, aEnode;
/*		unsigned int nno = AroundEnode(iblk_del,
			nblk,ColInd,RowPtr,
			aEnodeParentBlk,
            aAroundEnode, aEnode);*/
//		std::cout << "Del:" << iblk_del << "   Enode:" << aEnode.size() << "   Around:" << aAroundEnode.size() << std::endl;
		for(unsigned int i=0;i<aEnode.size();i++){
			unsigned int iblk0 = aEnode[i];
			aEnodeParentBlk[iblk0] = iblk_parent;
			aDegree[iblk0] = aAroundEnode.size();
		}
		std::vector<int> aFlg(nblk,-1);
		for(unsigned int iar=0;iar<aAroundEnode.size();iar++){
			unsigned int iblk0 = aAroundEnode[iar];
			unsigned int deg = 0;
			for(unsigned int ipsup=ColInd[iblk0];ipsup<ColInd[iblk0+1];ipsup++){
				const unsigned int jblk0 = RowPtr[ipsup];
				if( aEnodeParentBlk[jblk0] == -2 ){
                    if( aFlg[jblk0] == (int)iblk0 ) continue;
					deg++; 
					aFlg[jblk0] = iblk0;
//					std::cout << "a " << iblk0 << " " << jblk0 << std::endl;
				}
				else{
					if( aEnodeParentBlk[jblk0] == iblk_parent ){
						for(unsigned int i=0;i<aAroundEnode.size();i++){
							const unsigned int kblk0 = aAroundEnode[i];
							if( kblk0 == iblk0 ) continue;
                            if( aFlg[kblk0] == (int)iblk0 ) continue;
//							std::cout << "b " << iblk0 << " " << kblk0 << std::endl;
							deg++;
							aFlg[kblk0] = iblk0;
						}
					}
					else{
//						unsigned int iblkper0 = aEnodeParentBlk[jblk0];
						std::vector<unsigned int> aAroundEnode1, aEnode1;
/*						unsigned int nno = AroundEnode(jblk0,
							nblk,ColInd,RowPtr,
							aEnodeParentBlk,
                            aAroundEnode1, aEnode1);*/
						for(unsigned int jar=0;jar<aAroundEnode1.size();jar++){
							const unsigned int kblk0 = aAroundEnode1[jar];
							if( kblk0 == iblk0 ) continue;
                            if( aFlg[kblk0] == (int)iblk0 ) continue;
//							std::cout << "c " << iblk0 << " " << kblk0 << std::endl;
							deg++;
							aFlg[kblk0] = iblk0;
						}
					}
				}
			}
			aDegree[iblk0] = deg;
		}
		////////////////
	}

//	std::cout << "hoge" << std::endl;

/*	for(unsigned int iblk=0;iblk<nblk;iblk++){
		m_pOrder[iblk] = iblk;
	}*/
	m_pInvOrder = new int [nblk];
	for(unsigned int iblk=0;iblk<nblk;iblk++){ 
		int iblk0 = m_pOrder[iblk];
        assert( iblk0 >= 0 && iblk0 < (int)nblk );
		m_pInvOrder[iblk0] = iblk;
	}
}

