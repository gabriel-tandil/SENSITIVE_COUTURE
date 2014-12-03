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
// ElemAry.cpp: interface of element array class (CElemAry)
////////////////////////////////////////////////////////////////


#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
    #pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm>
#include <iostream>

#include "delfem/elem_ary.h"

using namespace Fem::Field;

////////////////
// Point
const unsigned int nlpoelPoint = 1;


////////////////
// LINE
const unsigned int nlpoelLine = 2;
const unsigned int nfaelLine = 2;
const unsigned int mxlpofaLine = 1;
const unsigned int nlpofaLine[nfaelLine] = {
	1, 1
};
const unsigned int lpofaLine[nfaelLine*mxlpofaLine] = {
	0, 
	1
};

////////////////
// Triangle
const unsigned int nlpoelTri = 3;
const unsigned int nfaelTri = 3;
const unsigned int mxlpofaTri = 2;
const unsigned int nlpofaTri[nfaelTri] = {
	2, 2, 2
};
const unsigned int lpofaTri[nfaelTri*mxlpofaTri] = {
	1, 2,
	2, 0,
	0, 1 
};


////////////////
// Quadric
const unsigned int nlpoelQuad = 4;
const unsigned int nfaelQuad = 4;
const unsigned int mxlpofaQuad = 2;
const unsigned int nlpofaQuad[nfaelQuad] = {
	2, 2, 2, 2
};
const unsigned int lpofaQuad[nfaelQuad*mxlpofaQuad] = {
	0, 1, 
	1, 2, 
	2, 3, 
	3, 0
};


////////////////
// Tetrahedra
const unsigned int nlpoelTet = 4;
const unsigned int nfaelTet = 4;
const unsigned int mxlpofaTet = 3;
const unsigned int nlpofaTet[nfaelTet] = {
	3, 3, 3, 3
};
const unsigned int lpofaTet[nfaelTet*mxlpofaTet] = {
	1, 2, 3,
	0, 3, 2,
	0, 1, 3,
	0, 2, 1
};


////////////////
// Hexahedra
const unsigned int nlpoelHex = 8;
const unsigned int nfaelHex = 6;
const unsigned int mxlpofaHex = 4;
const unsigned int nlpofaHex[nfaelHex] = {
	4, 4, 4, 4, 4, 4
};
const unsigned int lpofaHex[nfaelHex*mxlpofaHex] = {
	0, 3, 2, 1,
	0, 1, 5, 4,
	1, 2, 6, 5, 
	2, 3, 7, 6,
	3, 0, 4, 7,
	4, 5, 6, 7,
};
const unsigned int mxlp2lpedHex = 3;
const unsigned int nlp2lpedHex[nlpoelHex] = {
	3, 3, 3, 3, 3, 3, 3, 3
};
const unsigned int lp2lpedHex[nlpoelHex*mxlp2lpedHex] = {
	1, 3, 4,
	0, 2, 5,
	1, 3, 6,
	0, 2, 7, 
	0, 5, 7,
	1, 4, 6, 
	2, 5, 7, 
	3, 4, 6
};

namespace Fem{
namespace Field{

struct SElemInfo{
	unsigned int nlpoel;

	unsigned int nfael;
	unsigned int mxlpofa;
	const unsigned int* nlpofa;
	const unsigned int* lpofa;

	unsigned int  mxlp2lped;
	const unsigned int* nlp2lped;
	const unsigned int* lp2lped;
};

}
}

const SElemInfo ElemInfoAry[7] = {
	{	// 0
		0
	},
	{	// 1
		nlpoelPoint
	},
	{	// 2
		nlpoelLine,
		nfaelLine,  mxlpofaLine,  nlpofaLine,  lpofaLine
	},
	{	// 3
		nlpoelTri,
		nfaelTri,  mxlpofaTri,  nlpofaTri,  lpofaTri
	},
	{	// 4
		nlpoelQuad,
		nfaelQuad, mxlpofaQuad, nlpofaQuad, lpofaQuad
	},
	{	// 5
		nlpoelTet,
		nfaelTet,  mxlpofaTet,  nlpofaTet,  lpofaTet
	},
	{	// 6
		nlpoelHex,
		nfaelHex,  mxlpofaHex,  nlpofaHex,  lpofaHex,
		mxlp2lpedHex, nlp2lpedHex, lp2lpedHex
	}
};

CElemAry::CElemAry(const CElemAry& ea){
  m_nElem = ea.m_nElem;
  m_ElemType = ea.m_ElemType;
  npoel = ea.npoel;
  {
    const unsigned int n = m_nElem*npoel;
    m_pLnods = new unsigned int [n];
    for(unsigned int i=0;i<n;i++){ m_pLnods[i] = ea.m_pLnods[i]; }
  }
  m_aSeg = ea.m_aSeg;
}


// FEM用の行列パターンを作る（非対角ブロック行列用）
bool CElemAry::MakePattern_FEM(
	const unsigned int& id_es0, const unsigned int& id_es1, Com::CIndexedArray& crs ) const
{
	assert( m_aSeg.IsObjID(id_es0) );
	assert( m_aSeg.IsObjID(id_es1) );
	if( !m_aSeg.IsObjID(id_es0) ) return false;
	if( !m_aSeg.IsObjID(id_es1) ) return false;
  Com::CIndexedArray elsup;
	if( !this->MakePointSurElem( id_es0, elsup) ) goto FAILURE;
	{	// npoinの長さチェック
		const CElemSeg& es = m_aSeg.GetObj(id_es0);
		if( es.max_noes >= crs.Size() ){ crs.InitializeSize(es.max_noes+1); }
	}
	if( !this->MakePointSurPoint(id_es1, elsup, false, crs) ) goto FAILURE;
	crs.Sort();
	return true;

	// 失敗時の処理
FAILURE:
	assert(0);
	return false;
}

// FEM用の行列パターンを作る（対角ブロック行列用）
bool CElemAry::MakePattern_FEM
(const unsigned int& id_es, Com::CIndexedArray& crs ) const
{	
	assert( m_aSeg.IsObjID(id_es) );
	if( !m_aSeg.IsObjID(id_es) ) return false;
  Com::CIndexedArray elsup;
	if( !this->MakePointSurElem( id_es, elsup) ) goto FAILURE;
	{	// npoinの長さチェック
		const CElemSeg& es = m_aSeg.GetObj(id_es);
		if( es.max_noes >= crs.Size() ){	crs.InitializeSize(es.max_noes+1); }
	}  
	if( !this->MakePointSurPoint(id_es, elsup, true, crs) ) goto FAILURE;
	crs.Sort();
	assert( crs.CheckValid() );
	return true;

	// 失敗時の処理
FAILURE:
	assert(0);
	return false;
}

CElemAry* CElemAry::MakeBoundElemAry(
		unsigned int id_es_corner, unsigned int& id_es_add, 
		std::vector<unsigned int>& aIndElemFace) const
{
	if( !this->IsSegID(id_es_corner) ) return 0;
	
	unsigned int nfael, npoel_b;
	{
		nfael = ElemInfoAry[ m_ElemType ].nfael;
		npoel_b = ElemInfoAry[ m_ElemType ].mxlpofa;
	}

	int* elsuel = new int [m_nElem*nfael];
	this->MakeElemSurElem(id_es_corner,elsuel);

	std::string str_elem_type_b;
    ELEM_TYPE elem_type_b = Fem::Field::POINT;
	if( m_ElemType == TRI || m_ElemType == QUAD ){
		str_elem_type_b = "LINE";
		elem_type_b = LINE;
	}
	else if( m_ElemType == TET ){
		str_elem_type_b = "TRI";
		elem_type_b = TRI;
	}
	else if( m_ElemType == HEX ){
		str_elem_type_b = "QUAD";
		elem_type_b = QUAD;
	}

	unsigned int nelem_b;
	{
		unsigned int icoun0 = 0; 
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			for(unsigned int ifael=0;ifael<nfael;ifael++){
				if( elsuel[ielem*nfael+ifael] == -1 ) icoun0++;
			}
		}
		nelem_b = icoun0;
	}
	aIndElemFace.resize(nelem_b);
	std::vector<int> lnods_b;
	{
		const CElemSeg& es = this->GetSeg(id_es_corner);
		unsigned int inoes_s = es.begin;
		assert( inoes_s+es.m_nnoes <= npoel );
		lnods_b.resize(nelem_b*npoel_b,-1);
		const unsigned int* nlpofa = ElemInfoAry[ m_ElemType ].nlpofa;
		const unsigned int* lpofa = ElemInfoAry[ m_ElemType ].lpofa;
		unsigned int icoun0 = 0;
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
		for(unsigned int ifael=0;ifael<nfael;ifael++){
			if( elsuel[ielem*nfael+ifael] == -1 ){
				const unsigned int nilpofa = nlpofa[ifael];
				assert( nilpofa == npoel_b );
				for(unsigned int ilpofa=0;ilpofa<nilpofa;ilpofa++){
					lnods_b[icoun0*npoel_b+ilpofa] 
						= m_pLnods[ielem*npoel+ lpofa[ifael*npoel_b+ilpofa] + inoes_s];
				}
				aIndElemFace[icoun0] = ielem;
				icoun0++;
			}
		}
		}
	}
	delete[] elsuel;

  unsigned int id_na = this->GetSeg(id_es_corner).GetIdNA();
	CElemAry* ea_b = new CElemAry(nelem_b,elem_type_b);
	id_es_add = ea_b->AddSegment(ea_b->GetFreeSegID(),CElemSeg(id_na,CORNER),lnods_b);
	return ea_b;
}

// 辺を作る関数
bool CElemAry::MakeEdge(
		const unsigned int& id_es_corner, 
		unsigned int& nedge, 
		std::vector<unsigned int>& edge_ary) const
{
	if( !this->IsSegID(id_es_corner) ) return false;
	const CElemSeg& es = this->GetSeg(id_es_corner);
	unsigned int ibegin_es = es.begin;
	unsigned int len = es.m_nnoes;
	/*
	if( elem_type == TET ){
		const unsigned int npoin = end_ipoin;

		unsigned int* elsup_ind = new unsigned int [npoin+1];
		unsigned int* elsup = 0;
		unsigned int nelsup;
		this->MakePointSurElem(npoin,  elsup_ind,nelsup,elsup);

		unsigned int* psup_ind = new unsigned int [npoin+1];
		unsigned int npsup;
		unsigned int* psup;
		this->MakePointSurPoint(npoin,elsup_ind,nelsup,elsup,  psup_ind,npsup,psup );

		delete[] elsup_ind;	delete[] elsup;

		unsigned int icoun0 = 0;
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			for(unsigned ipsup=psup_ind[ipoin];ipsup<psup_ind[ipoin+1];ipsup++){
				unsigned int jpoin0 = psup[ipsup];
				if( ipoin < jpoin0 ) icoun0++;
			}
		}

		nedge = icoun0;
		edge_ary = new unsigned int [nedge*2];

		icoun0 = 0;
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			for(unsigned ipsup=psup_ind[ipoin];ipsup<psup_ind[ipoin+1];ipsup++){
				unsigned int jpoin0 = psup[ipsup];
				if( ipoin < jpoin0 ){
					edge_ary[icoun0*2  ] = ipoin;
					edge_ary[icoun0*2+1] = jpoin0;
					icoun0++;
				}
			}
		}

		delete[] psup_ind;	delete[] psup;

		return true;
	}
	else */
	if( this->ElemType() == TRI ){
		if( len != 3 ) return false;
		const SElemInfo& elem_info = ElemInfoAry[ this->ElemType() ];
		const unsigned int nfael = elem_info.nfael;
		assert( nfael == 3 );
		int* elsuel = new int [m_nElem*nfael];
		this->MakeElemSurElem(id_es_corner,elsuel);

		unsigned int icoun0 = 0;
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			if( m_pLnods[ielem*npoel+ibegin_es  ] > m_pLnods[ielem*npoel+ibegin_es+1] || elsuel[ielem*nfael+2] == -1 ){ icoun0++; }
			if( m_pLnods[ielem*npoel+ibegin_es+1] > m_pLnods[ielem*npoel+ibegin_es+2] || elsuel[ielem*nfael  ] == -1 ){ icoun0++; }
			if( m_pLnods[ielem*npoel+ibegin_es+2] > m_pLnods[ielem*npoel+ibegin_es  ] || elsuel[ielem*nfael+1] == -1 ){ icoun0++; }
		}
		nedge = icoun0;
		edge_ary.resize(nedge*2);

		icoun0 = 0;
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			if( m_pLnods[ielem*npoel+ibegin_es  ] > m_pLnods[ielem*npoel+ibegin_es+1] || elsuel[ielem*nfael+2] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+ibegin_es  ]; 
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+ibegin_es+1];
				icoun0++;
			}
			if( m_pLnods[ielem*npoel+ibegin_es+1] > m_pLnods[ielem*npoel+ibegin_es+2] || elsuel[ielem*nfael  ] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+ibegin_es+1];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+ibegin_es+2];
				icoun0++;
			}
			if( m_pLnods[ielem*npoel+ibegin_es+2] > m_pLnods[ielem*npoel+ibegin_es  ] || elsuel[ielem*nfael+1] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+ibegin_es+2];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+ibegin_es  ];
				icoun0++;
			}
		}
		assert( nedge == icoun0 );

		delete[] elsuel;
		return true;
	}
	else if( this->ElemType() == LINE ){
		if( len != 2 ) return false;
		nedge = m_nElem;
		edge_ary.resize(nedge*2);
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			edge_ary[ielem*2  ] = m_pLnods[ielem*npoel+ibegin_es  ];
			edge_ary[ielem*2+1] = m_pLnods[ielem*npoel+ibegin_es+1];
		}
	}
	else if( this->ElemType() == QUAD ){
		const SElemInfo& elem_info = ElemInfoAry[ this->ElemType() ];
		const unsigned int nfael = elem_info.nfael;
		assert( nfael == 4 );
		int* elsuel = new int [m_nElem*nfael];
		this->MakeElemSurElem(id_es_corner,elsuel);

		unsigned int icoun0 = 0;

		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			if( m_pLnods[ielem*npoel  ] > m_pLnods[ielem*npoel+1] || elsuel[ielem*nfael  ] == -1 ){ icoun0++; }
			if( m_pLnods[ielem*npoel+1] > m_pLnods[ielem*npoel+2] || elsuel[ielem*nfael+1] == -1 ){ icoun0++; }
			if( m_pLnods[ielem*npoel+2] > m_pLnods[ielem*npoel+3] || elsuel[ielem*nfael+2] == -1 ){ icoun0++; }
			if( m_pLnods[ielem*npoel+3] > m_pLnods[ielem*npoel  ] || elsuel[ielem*nfael+3] == -1 ){ icoun0++; }
		}
		nedge = icoun0;
		edge_ary.resize(nedge*2);

		icoun0 = 0;
		for(unsigned int ielem=0;ielem<m_nElem;ielem++){
			if( m_pLnods[ielem*npoel  ] > m_pLnods[ielem*npoel+1] || elsuel[ielem*nfael  ] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel  ];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+1];
				icoun0++;
			}
			if( m_pLnods[ielem*npoel+1] > m_pLnods[ielem*npoel+2] || elsuel[ielem*nfael+1] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+1];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+2];
				icoun0++;
			}
			if( m_pLnods[ielem*npoel+2] > m_pLnods[ielem*npoel+3] || elsuel[ielem*nfael+2] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+2];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel+3];
				icoun0++;
			}
			if( m_pLnods[ielem*npoel+3] > m_pLnods[ielem*npoel  ] || elsuel[ielem*nfael+3] == -1 ){
				edge_ary[icoun0*2  ] = m_pLnods[ielem*npoel+3];
				edge_ary[icoun0*2+1] = m_pLnods[ielem*npoel  ];
				icoun0++;
			}
		}
		assert( nedge == icoun0 );

		delete[] elsuel;
		return true;
	}/*
	else if( elem_type == HEX ){
		const unsigned int npoin = end_ipoin;

		unsigned int* elsup_ind = 0;
		unsigned int nelsup;
		unsigned int* elsup2 = 0;
		{
			elsup_ind = new unsigned int [npoin+1];
			for(unsigned int ipoin=0;ipoin<npoin+1;ipoin++){
				elsup_ind[ipoin] = 0;
			}
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				for(unsigned int inoel=0;inoel<npoel;inoel++){
					const unsigned int ipoin0 = m_pLnods[ielem*npoel+inoel];
					assert( ipoin0 < npoin );
					elsup_ind[ipoin0+1]++;
				}
			}
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				elsup_ind[ipoin+1] += elsup_ind[ipoin];
			}
			nelsup = elsup_ind[npoin];
			elsup2 = new unsigned int [nelsup*2];
			for(unsigned int ielem=0;ielem<nelem;ielem++){
				for(unsigned int inoel=0;inoel<npoel;inoel++){
					const unsigned int ipoin0 = m_pLnods[ielem*npoel+inoel];
					assert( ipoin0 < npoin );
					const unsigned int ipsup0 = elsup_ind[ipoin0];
					assert( ipsup0 < nelsup );
					elsup2[ipsup0*2  ] = ielem;
					elsup2[ipsup0*2+1] = inoel;
					elsup_ind[ipoin0]++;
				}
			}
			for(unsigned int ipoin=npoin;ipoin>0;ipoin--){
				elsup_ind[ipoin] = elsup_ind[ipoin-1];
			}
			elsup_ind[0] = 0;
		}

		unsigned int* psup_ind = new unsigned int [npoin+1];
		unsigned int npsup;
		unsigned int* psup = 0;
		{
			const unsigned int mxlp2lped = ElemInfoAry[ elem_type ].mxlp2lped;
			const unsigned int* nlp2lped = ElemInfoAry[ elem_type ].nlp2lped;
			const unsigned int* lp2lped = ElemInfoAry[ elem_type ].lp2lped;

			int* tmp_poin = new int [npoin];
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				tmp_poin[ipoin] = -1;
			}
			int icoun0 = 0;
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				tmp_poin[ipoin] = ipoin;
				for(unsigned int ielsup=elsup_ind[ipoin];ielsup<elsup_ind[ipoin+1];ielsup++){
					const unsigned int ielem0 = elsup2[ielsup*2  ];
					const unsigned int ipoel0 = elsup2[ielsup*2+1];
					const unsigned int nilp2lped = nlp2lped[ipoel0];
					for(unsigned int ilp2lped=0;ilp2lped<nilp2lped;ilp2lped++){
						const unsigned int jpoel0 = lp2lped[ipoel0*mxlp2lped+ilp2lped];
						const unsigned int jpoin0 = m_pLnods[ielem0*npoel+jpoel0];
						if( tmp_poin[jpoin0] != ipoin ){
							tmp_poin[jpoin0] = ipoin;
							icoun0++;
						}
					}
				}
			}

			npsup = icoun0;
			psup = new unsigned int [npsup];

			icoun0 = 0;
			psup_ind[0] = 0;
			for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
				tmp_poin[ipoin] = ipoin;
				for(unsigned int ielsup=elsup_ind[ipoin];ielsup<elsup_ind[ipoin+1];ielsup++){
					const unsigned int ielem0 = elsup2[ielsup*2  ];
					const unsigned int ipoel0 = elsup2[ielsup*2+1];
					const unsigned int nilp2lped = nlp2lped[ipoel0];
					for(unsigned int ilp2lped=0;ilp2lped<nilp2lped;ilp2lped++){
						const unsigned int jpoel0 = lp2lped[ipoel0*mxlp2lped+ilp2lped];
						const unsigned int jpoin0 = m_pLnods[ielem0*npoel+jpoel0];
						if( tmp_poin[jpoin0] != ipoin ){
							tmp_poin[jpoin0] = ipoin;
							psup[icoun0] = jpoin0;
							icoun0++;
						}
					}
				}
				psup_ind[ipoin+1] = icoun0;
			}
			delete[] tmp_poin;
		}

		delete[] elsup_ind;	delete[] elsup2;

		unsigned int icoun0 = 0;
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			for(unsigned ipsup=psup_ind[ipoin];ipsup<psup_ind[ipoin+1];ipsup++){
				unsigned int jpoin0 = psup[ipsup];
				if( ipoin < jpoin0 ) icoun0++;
			}
		}

		nedge = icoun0;
		edge_ary = new unsigned int [nedge*2];

		icoun0 = 0;
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			for(unsigned ipsup=psup_ind[ipoin];ipsup<psup_ind[ipoin+1];ipsup++){
				unsigned int jpoin0 = psup[ipsup];
				if( ipoin < jpoin0 ){
					edge_ary[icoun0*2  ] = ipoin;
					edge_ary[icoun0*2+1] = jpoin0;
					icoun0++;
				}
			}
		}

		delete[] psup_ind;	delete[] psup;
	}*/
	else{
		std::cout << "Error!-->Not Implimented" << std::endl;
		std::cout << "Elem Type : " << this->ElemType() << std::endl;
		assert(0);
		abort();
		nedge = 0;
	}
	return true;
}

bool CElemAry::MakeElemSurElem(const unsigned int& id_es_corner,int* elsuel) const
{
    Com::CIndexedArray elsup;
	this->MakePointSurElem(id_es_corner,elsup);
	this->MakeElemSurElem(id_es_corner,elsup,elsuel);
	return true;
}

bool CElemAry::MakeElemSurElem( const unsigned int& id_es, 
		const Com::CIndexedArray& elsup, int* elsuel ) const
{
	assert( m_aSeg.IsObjID(id_es) );
	if( !m_aSeg.IsObjID(id_es) ){ return false; }
	const CElemSeg& es = m_aSeg.GetObj(id_es);
	const unsigned int inoel_s = es.begin;
	const unsigned int inoel_e = inoel_s + es.m_nnoes;
	const unsigned int max_noes = es.max_noes;
	assert( max_noes < elsup.Size() );
	assert( inoel_e <= npoel );

	SElemInfo elem_info = ElemInfoAry[ this->ElemType() ];
	const unsigned int nfael = elem_info.nfael;
	const unsigned int mxlpofa = elem_info.mxlpofa;
	const unsigned int* nlpofa = elem_info.nlpofa;
	const unsigned int* lpofa = elem_info.lpofa;

	const unsigned int npoin = elsup.Size();

	unsigned int* tmp_poin = new unsigned int [npoin];
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){ tmp_poin[ipoin] = 0; }
	unsigned int inpofa[4];

	for(unsigned int ielem=0;ielem<m_nElem;ielem++){
	for(unsigned int ifael=0;ifael<nfael;ifael++){
		const unsigned int nilpofa = nlpofa[ifael];
		for(unsigned int ipofa=0;ipofa<nilpofa;ipofa++){
			inpofa[ipofa] = m_pLnods[ielem*npoel + lpofa[ifael*mxlpofa+ipofa] + inoel_s];
			tmp_poin[ inpofa[ipofa] ] = 1;
		}
		const unsigned int ipoin0= inpofa[0];
		bool iflg = false;
		for(unsigned int ielsup=elsup.index[ipoin0];ielsup<elsup.index[ipoin0+1];ielsup++){
			const unsigned int jelem0 = elsup.array[ielsup];
			if( jelem0 == ielem ) continue;
			for(unsigned int jfael=0;jfael<nfael;jfael++){
				iflg = true;
				const unsigned int njlpofa = nlpofa[jfael];
				for(unsigned int jpofa=0;jpofa<njlpofa;jpofa++){
					const unsigned int jpoin0 = m_pLnods[jelem0*npoel + lpofa[jfael*mxlpofa+jpofa] + inoel_s];
					assert( jpoin0 < npoin );
					if( tmp_poin[ jpoin0 ] == 0 ){
						iflg = false;
						break;
					}
				}
				if( iflg ){
					elsuel[ielem*nfael+ifael] = jelem0;
					break;
				}
			}
			if( iflg ) break;
		}
		if( !iflg ){ elsuel[ielem*nfael+ifael] = -1; }
		for(unsigned int ipofa=0;ipofa<nilpofa;ipofa++){
			tmp_poin[ inpofa[ipofa] ] = 0;
		}
	}	
	}

	delete[] tmp_poin;
	return true;
}

bool CElemAry::MakePointSurElem( const unsigned int id_es, Com::CIndexedArray& elsup ) const
{
	assert( m_aSeg.IsObjID(id_es) );
	if( !m_aSeg.IsObjID(id_es) ){ return false; }
	const CElemSeg& es = m_aSeg.GetObj(id_es);
	const unsigned int inoel_s = es.begin;
	const unsigned int inoel_e = inoel_s + es.m_nnoes;
	if( elsup.Size() <= es.GetMaxNoes() ){ elsup.InitializeSize(es.GetMaxNoes()+1); }
	const unsigned int npoin = elsup.Size();

	elsup.index.resize(npoin+1);
	for(unsigned int ipoin=0;ipoin<npoin+1;ipoin++){ elsup.index[ipoin] = 0; }
	for(unsigned int ielem=0;ielem<m_nElem;ielem++){
		for(unsigned int inoel=inoel_s;inoel<inoel_e;inoel++){
			const unsigned int ipoin0 = m_pLnods[ielem*npoel+inoel];
			assert( ipoin0 < npoin );
			elsup.index[ipoin0+1]++;
		}
	}
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		elsup.index[ipoin+1] += elsup.index[ipoin];
	}
	const unsigned int nelsup = elsup.index[npoin];
	elsup.array.resize(nelsup);
	for(unsigned int ielem=0;ielem<m_nElem;ielem++){
		for(unsigned int inoel=inoel_s;inoel<inoel_e;inoel++){
			const unsigned int ipoin0 = m_pLnods[ielem*npoel+inoel];
			assert( ipoin0 < npoin );
			const unsigned int ipsup0 = elsup.index[ipoin0];
			assert( ipsup0 < nelsup );
			elsup.array[ipsup0] = ielem;
			elsup.index[ipoin0]++;
		}
	}
	for(unsigned int ipoin=npoin;ipoin>0;ipoin--){
		elsup.index[ipoin] = elsup.index[ipoin-1];
	}
	elsup.index[0] = 0;

//	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
//		std::cout << ipoin << "--> ";
//		for(unsigned int ielsup=elsup.index[ipoin];ielsup<elsup.index[ipoin+1];ielsup++){
//			std::cout << elsup.array[ielsup] << " ";
//		}
//		std::cout << std::endl;
//	}


	return true;
}

/*
bool CElemAry_UNS::MakePointSurPoint(
	const unsigned int& npoin, unsigned int* psup_ind, unsigned int& npsup, unsigned int*& psup ) const {

	assert( npoin >= end_ipoin );
	unsigned int* elsup_ind = new unsigned int [npoin+1];
	unsigned int* elsup;
	unsigned int nelsup;
	this->MakePointSurElem(npoin, elsup_ind,nelsup,elsup);
	this->MakePointSurPoint(npoin,elsup_ind,nelsup,elsup,  psup_ind,npsup,psup);
	delete[] elsup_ind;
	delete[] elsup;
	return true;
}
*/

bool CElemAry::MakeElemToEdge
(const unsigned int& id_es_corner, 
 const unsigned int& nedge, const std::vector<unsigned int>& edge_ary,
 std::vector<int>& el2ed ) const
{
	assert( nedge*2 == edge_ary.size () );
	assert( m_aSeg.IsObjID(id_es_corner) );
	if( !m_aSeg.IsObjID(id_es_corner) ) return false;
	const CElemSeg& es = m_aSeg.GetObj(id_es_corner);
	assert( es.GetElSegType() == CORNER );
	const unsigned int neled = CElemSeg::GetLength( EDGE, this->ElemType() );
	el2ed.resize( neled * this->Size(), -1 );

  Com::CIndexedArray elsup;
	this->MakePointSurElem(id_es_corner,elsup);

	if( this->ElemType() == TRI ){
		const unsigned int nEdTri = 3;
		for(unsigned int iedge=0;iedge<nedge;iedge++){
			const unsigned int ipo0 = edge_ary[iedge*2  ];
			const unsigned int ipo1 = edge_ary[iedge*2+1];
			assert( ipo0 < elsup.Size() );
			for(unsigned int ielsup=elsup.index[ipo0];ielsup<elsup.index[ipo0+1];ielsup++){
				if( ielsup >= elsup.array.size() ) continue;
				const unsigned int ielem0 = elsup.array[ielsup];
				for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
					unsigned int jpo0 = m_pLnods[ielem0*npoel+es.begin+lpofaTri[iedtri*2+0]];
					unsigned int jpo1 = m_pLnods[ielem0*npoel+es.begin+lpofaTri[iedtri*2+1]];
					if( (ipo0-jpo0)*(ipo1-jpo0)==0 && (ipo0-jpo1)*(ipo1-jpo1)==0 ){
						el2ed[ielem0*nEdTri+iedtri] = iedge;
						break;
					}
				}
			}
		}
	}
	else if( this->ElemType() == QUAD ){
		const unsigned int nEdQuad = 4;
		for(unsigned int iedge=0;iedge<nedge;iedge++){
			const unsigned int ipo0 = edge_ary[iedge*2  ];
			const unsigned int ipo1 = edge_ary[iedge*2+1];
			assert( ipo0 < elsup.Size() );
			for(unsigned int ielsup=elsup.index[ipo0];ielsup<elsup.index[ipo0+1];ielsup++){
				if( ielsup >= elsup.array.size() ) continue;
				const unsigned int ielem0 = elsup.array[ielsup];
				for(unsigned int iedquad=0;iedquad<nEdQuad;iedquad++){
					unsigned int jpo0 = m_pLnods[ielem0*npoel+es.begin+lpofaQuad[iedquad*2+0]];
					unsigned int jpo1 = m_pLnods[ielem0*npoel+es.begin+lpofaQuad[iedquad*2+1]];
					if( (ipo0-jpo0)*(ipo1-jpo0)==0 && (ipo0-jpo1)*(ipo1-jpo1)==0 ){
						el2ed[ielem0*nEdQuad+iedquad] = iedge;
						break;
					}
				}
			}
		}
	}
	else if( this->ElemType() == LINE ){
		// elemに対してedgeが多い場合、効率が悪い
		for(unsigned int iedge=0;iedge<nedge;iedge++){
			const unsigned int ipo0 = edge_ary[iedge*2  ];
			const unsigned int ipo1 = edge_ary[iedge*2+1];
			if( ipo0+1 >= elsup.index.size() ) continue;
			for(unsigned int ielsup=elsup.index[ipo0];ielsup<elsup.index[ipo0+1];ielsup++){
				if( ielsup >= elsup.array.size() ) break;
				const unsigned int ielem0 = elsup.array[ielsup];
				unsigned int jpo0 = m_pLnods[ielem0*npoel+es.begin+0];
				unsigned int jpo1 = m_pLnods[ielem0*npoel+es.begin+1];
				if( (ipo0-jpo0)*(ipo1-jpo0)==0 && (ipo0-jpo1)*(ipo1-jpo1)==0 ){
					el2ed[ielem0] = iedge;
					break;
				}
			}
		}
	}
	else{ assert(0); }
	return true;
}


bool CElemAry::MakePointSurPoint( 
	const unsigned int id_es, const Com::CIndexedArray& elsup, bool isnt_self, 
	Com::CIndexedArray& psup ) const 
{
	assert( m_aSeg.IsObjID(id_es) );
	const CElemSeg& es = m_aSeg.GetObj(id_es);
	unsigned int jnoel_s = es.begin;
	unsigned int jnoel_e = jnoel_s + es.m_nnoes;
	unsigned int npoin_j = es.GetMaxNoes()+1;

	if( psup.Size() < elsup.Size() ){ psup.InitializeSize(elsup.Size()); }
	const unsigned int npoin_i = psup.Size();

	int* lpoin = new int [npoin_j];
	for(unsigned int jpoin=0;jpoin<npoin_j;jpoin++){ lpoin[jpoin] = -1; }
	int icoun0 = 0;
	for(unsigned int ipoin=0;ipoin<elsup.Size();ipoin++){
		if( isnt_self ){ lpoin[ipoin] = ipoin; }
		for(unsigned int ielsup=elsup.index[ipoin];ielsup<elsup.index[ipoin+1];ielsup++){
			const unsigned int jelem0 = elsup.array[ielsup];
			assert( jelem0 < m_nElem );
			for(unsigned int jnoel=jnoel_s;jnoel<jnoel_e;jnoel++){
				const unsigned int jpoin0 = m_pLnods[jelem0*npoel+jnoel];
				assert( jpoin0 < npoin_j );
        if( lpoin[jpoin0] != (int)ipoin ){
					lpoin[jpoin0] = ipoin;
					icoun0++;
				}
			}
		}
	}

	const unsigned int npsup = icoun0;
	psup.index.resize(npoin_i+1);
	psup.array.resize(npsup);

	for(unsigned int jpoin=0;jpoin<npoin_j;jpoin++){ lpoin[jpoin] = -1; }
	icoun0 = 0;
	psup.index[0] = 0;
	for(unsigned int ipoin=0;ipoin<elsup.Size();ipoin++){
		if( isnt_self ){ lpoin[ipoin] = ipoin; }
		for(unsigned int ielsup=elsup.index[ipoin];ielsup<elsup.index[ipoin+1];ielsup++){
			const unsigned int jelem0 = elsup.array[ielsup];
			assert( jelem0 < m_nElem );
			for(unsigned int jnoel=jnoel_s;jnoel<jnoel_e;jnoel++){
				const unsigned int jpoin0 = m_pLnods[jelem0*npoel+jnoel];
				assert( jpoin0 < npoin_j );
        if( lpoin[jpoin0] != (int)ipoin ){
					lpoin[jpoin0] = ipoin;
					psup.array[icoun0] = jpoin0;
					icoun0++;
				}
			}
		}
		psup.index[ipoin+1] = icoun0;
	}
	for(unsigned int ipoin=elsup.Size();ipoin<npoin_i;ipoin++){
		psup.index[ipoin+1] = icoun0;
	}
	delete[] lpoin;
  
/*
	for(unsigned int ipoin=0;ipoin<npoin_i;ipoin++){
		std::cout << ipoin << "--> ";
		for(unsigned int ipsup=psup.index[ipoin];ipsup<psup.index[ipoin+1];ipsup++){
			std::cout << psup.array[ipsup] << " ";
		}
		std::cout << std::endl;
	}
*/
	return true;
}


int CElemAry::AddSegment
(unsigned int id, 
 const CElemAry::CElemSeg& es, 
 const std::vector<int>& lnods )
{
	std::vector< std::pair<unsigned int,CElemSeg> > aEs;
	aEs.push_back( std::make_pair(id,es) );
	const std::vector<int>& res = this->AddSegment( aEs, lnods );
	assert( res.size() == 1 );
	return res[0];
}

std::vector<int> CElemAry::AddSegment
(const std::vector< std::pair<unsigned int,CElemSeg> >& es_ary_in, 
 const std::vector<int>& lnods )
{
	std::vector< std::pair<unsigned int,CElemSeg> > es_ary = es_ary_in;

	unsigned int npoel_add=0;
	for(unsigned int ies=0;ies<es_ary.size();ies++){
		es_ary[ies].second.m_nnoes = CElemSeg::GetLength(es_ary[ies].second.GetElSegType(),this->ElemType());
		npoel_add += es_ary[ies].second.m_nnoes;
	}
	assert( lnods.size() == npoel_add*Size() );
	const unsigned int npoel_new = npoel+npoel_add;
	unsigned int* lnods_new = new unsigned int [npoel_new*Size()];
	for(unsigned int ielem=0;ielem<m_nElem;ielem++){
		for(unsigned int inoel=0;inoel<npoel;inoel++){
			lnods_new[ielem*npoel_new+inoel] = m_pLnods[ielem*npoel+inoel];
		}
		for(unsigned int inoel=0;inoel<npoel_add;inoel++){
			lnods_new[ielem*npoel_new+npoel+inoel] = lnods[ielem*npoel_add+inoel];
		}
	}
	delete[] m_pLnods;
	m_pLnods = lnods_new;

	std::vector<int> add_es_id_ary;

	{	// 要素セグメント情報の追加
		// beginを作る
		for(unsigned int ies=0;ies<es_ary.size();ies++){
			es_ary[ies].second.begin = npoel;
			npoel += es_ary[ies].second.m_nnoes;
		}
		assert( npoel == npoel_new );
		// 要素セグメントの最大の節点番号を作る
		for(unsigned int ies=0;ies<es_ary.size();ies++){
			unsigned int inoel_s = es_ary[ies].second.begin;
			unsigned int inoel_e = inoel_s + es_ary[ies].second.m_nnoes;
			unsigned int max_noes = 0;
			for(unsigned int ielem=0;ielem<m_nElem;ielem++){
				for(unsigned int inoel=inoel_s;inoel<inoel_e;inoel++){
					max_noes = ( (unsigned int)abs(m_pLnods[ielem*npoel+inoel]) > max_noes ) 
						? (unsigned int)abs(m_pLnods[ielem*npoel+inoel]) : max_noes;
				}
			}
			es_ary[ies].second.max_noes = max_noes;
		}
		// 要素セグメントの追加
		for(unsigned int ies=0;ies<es_ary.size();ies++){
			const int add_es_id = m_aSeg.AddObj( std::make_pair(es_ary[ies].first,es_ary[ies].second) );
			add_es_id_ary.push_back(add_es_id);
		}
	}
	return add_es_id_ary;
}

int CElemAry::WriteToFile(const std::string& file_name, long& offset, unsigned int id) const {

	FILE *fp;
	if( (fp = fopen(file_name.c_str(),"a"))== NULL ){
//		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return false;
	}

	fprintf(fp,"$$$$$$\n");
	fprintf(fp,"ELEM_ARY\n");

	fprintf(fp,"%d\n",id);

	switch( this->ElemType() ){
	case HEX:
		fprintf(fp,"HEX\n");	break;
	case TRI:
		fprintf(fp,"TRI\n");	break;
	case QUAD:
		fprintf(fp,"QUAD\n");	break;
	case LINE:
		fprintf(fp,"LINE\n");	break;
	case POINT:
		fprintf(fp,"POINT\n");	break;
	default:
		assert(0);
	}

	{	// 要素セグメントの出力
		std::vector<unsigned int> id_ary_es = this->m_aSeg.GetAry_ObjID();
		fprintf(fp,"%d\n",(int)id_ary_es.size());
		for(unsigned int iid_es=0;iid_es<id_ary_es.size();iid_es++){
			unsigned int id_es = id_ary_es[iid_es];
			assert( m_aSeg.IsObjID(id_es) );
			const CElemSeg& es = m_aSeg.GetObj(id_es);
			std::string str_es_type;
			{
				if(      es.GetElSegType() == CORNER ){ str_es_type = "CORNER"; }
				else if( es.GetElSegType() == EDGE   ){ str_es_type = "EDGE";   }
				else if( es.GetElSegType() == BUBBLE ){ str_es_type = "BUBBLE"; }
				else{ assert(0); }
			}
			fprintf(fp,"%d %d %d %s\n", iid_es+1, id_es, es.GetIdNA(), str_es_type.c_str() );
		}
	}
	{	// 節点番号の出力
		fprintf(fp,"%d %d\n",this->Size(),npoel);
		for(unsigned int ielem=0;ielem<this->Size();ielem++){
			fprintf(fp,"%d ",ielem+1);
			for(unsigned int inoel=0;inoel<npoel;inoel++){
				fprintf(fp,"%d ",this->m_pLnods[ielem*npoel+inoel]+1);
			}
			fprintf(fp,"\n");
		}
	}

	fclose(fp);

	return 0;
}

int CElemAry::InitializeFromFile(const std::string& file_name, long& offset)
{

	FILE *fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];

	if( (fp = fopen(file_name.c_str(),"r"))== NULL ){
//		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return 1;
	}

	fseek(fp,offset,SEEK_SET);

	////////////////////////////////
	
	std::string str_elem_type;
	int nseg;
	std::vector< std::pair<unsigned int,CElemSeg> > tmp_es_ary;			
	int tmp_nelem, tmp_nnoel;
	std::vector<unsigned int> tmp_lnods;
	{
		while( 1 ){
			if( fgets(stmp1,buff_size,fp) == NULL ){ return 1; }
			if( stmp1[0] == '#' ) continue;
			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
		}
		assert( strncmp(stmp1,"$$$$$$$$",6)==0 );

		while( 1 ){
			if( fgets(stmp1,buff_size,fp) == NULL ){ return -1; }
			if( stmp1[0] == '#' ) continue;
			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
		}
		assert( strncmp(stmp1,"ELEM_ARY",8)==0 );

		{	// IDの読み込み（使用しない）
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			int tmp_id;
			sscanf(stmp1,"%d",&tmp_id);
			assert( tmp_id > 0 );
//			unsigned int id = tmp_id;
		}

		{	// 要素タイプの読み込み
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			const unsigned int len = strlen(stmp1);
			stmp1[len-1] = '\0';
			str_elem_type = stmp1;
		}

		{	// セグメント数の読み込み
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			int tmp_nseg;
			sscanf(stmp1,"%d",&tmp_nseg);
			assert( tmp_nseg > 0 );
			nseg = tmp_nseg;
		}

		{	// 要素セグメントの読み込み
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			int tmp_iseg, tmp_id, tmp_id_na;
			char tmp_type_name[256];
			for(int iseg=0;;iseg++){
				sscanf(stmp1,"%d%d%d%s",&tmp_iseg, &tmp_id, &tmp_id_na, tmp_type_name);
				assert( tmp_iseg == iseg+1 );
				assert( tmp_id > 0 );
				assert( tmp_id_na > 0 );
				ELSEG_TYPE tmp_type;
				if(      strncmp(tmp_type_name,"CORNER",6)==0 ){ tmp_type = CORNER; }
				else if( strncmp(tmp_type_name,"EDGE"  ,4)==0 ){ tmp_type = EDGE;   }
				else if( strncmp(tmp_type_name,"BUBBLE",6)==0 ){ tmp_type = BUBBLE; }
				else{ assert(0); }
				CElemSeg tmp_es(tmp_id_na,tmp_type);
				tmp_es_ary.push_back( std::make_pair(tmp_id,tmp_es) );
				if( iseg >= nseg-1 ) break;
				fgets(stmp1,buff_size,fp);
			}
		}

		{	// 要素数、要素あたりの節点数の読み込み
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			sscanf(stmp1,"%d%d",&tmp_nelem, &tmp_nnoel);
			assert( tmp_nelem > 0 && tmp_nnoel > 0 );
		}

		{	// 要素節点番号の読み込み
			tmp_lnods.resize(tmp_nelem*tmp_nnoel);
			while( fgets(stmp1,buff_size,fp) != NULL ){
				if( stmp1[0] == '#' ) continue;
				if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			}
			char *tp;
			int inode0;
			int tmp_ielem;
			for(unsigned int ielem=0;;ielem++){
				tp = strtok(stmp1," ");
				tmp_ielem = atoi(tp);
                assert( tmp_ielem-1 == (int)ielem );
				for(unsigned int inoel=0;inoel<(unsigned int)tmp_nnoel;inoel++){
					tp = strtok(NULL," ");
					inode0 = atoi(tp)-1; 
					assert( inode0 >= 0 );
					tmp_lnods[ielem*tmp_nnoel+inoel] = inode0;
				}
				if( ielem >= (unsigned int)tmp_nelem-1 ) break;
				fgets(stmp1,buff_size,fp);
			}
		}
	}

	offset = ftell(fp);
	fclose(fp);

	{
		this->m_nElem = tmp_nelem;
		this->npoel = tmp_nnoel;
		if( this->m_pLnods != 0 ){ delete[] this->m_pLnods; }
		this->m_pLnods = new unsigned int [m_nElem*npoel];
		for(unsigned int i=0;i<m_nElem*npoel;i++){ this->m_pLnods[i] = tmp_lnods[i]; }
		if(      str_elem_type == "LINE" ){ m_ElemType = LINE; }
		else if( str_elem_type == "TRI"  ){ m_ElemType = TRI;  }
		else if( str_elem_type == "QUAD" ){ m_ElemType = QUAD; }
		else if( str_elem_type == "TET"  ){ m_ElemType = TET;  }
		else if( str_elem_type == "HEX"  ){ m_ElemType = HEX;  }
		else if( str_elem_type == "POINT"){ m_ElemType = POINT;}
		else{
//			std::cout << "Error!-->Not Implimented (Elem Type:" << str_elem_type << ")" << std::endl;
			assert(0);
			abort();
		}
		{	// 要素セグメントの追加

			// len と beignを作る
			unsigned int ipoel_cnt = 0;
			for(unsigned int ies=0;ies<tmp_es_ary.size();ies++){
				tmp_es_ary[ies].second.begin = ipoel_cnt;
				tmp_es_ary[ies].second.m_nnoes = CElemSeg::GetLength(tmp_es_ary[ies].second.GetElSegType(),this->ElemType());
				ipoel_cnt += tmp_es_ary[ies].second.m_nnoes;
			}
			// 要素セグメントの最大の節点番号を作る
			assert( ipoel_cnt == npoel );
			for(unsigned int ies=0;ies<tmp_es_ary.size();ies++){
				unsigned int inoel_s = tmp_es_ary[ies].second.begin;
				unsigned int inoel_e = inoel_s + tmp_es_ary[ies].second.m_nnoes;
				unsigned int max_noes = 0;
				for(unsigned int ielem=0;ielem<m_nElem;ielem++){
					for(unsigned int inoel=inoel_s;inoel<inoel_e;inoel++){
						max_noes = ( (unsigned int)abs(m_pLnods[ielem*npoel+inoel]) > max_noes ) 
							? (unsigned int)abs(m_pLnods[ielem*npoel+inoel]) : max_noes;
					}
				}
				tmp_es_ary[ies].second.max_noes = max_noes;
			}
			// 要素セグメントの追加
			for(unsigned int ies=0;ies<tmp_es_ary.size();ies++){
				m_aSeg.AddObj( std::make_pair(tmp_es_ary[ies].first, tmp_es_ary[ies].second) );
			}
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


/*
CElemAry_Rect::CElemAry_Rect(double len_x, double len_y, unsigned int div_x, unsigned int div_y)
: m_len_x(len_x), m_len_y(len_y), m_div_x(div_x), m_div_y(div_y)
{
	this->elem_type = QUAD;
	this->nelem = m_div_x*m_div_y;
	this->npoel = 4;
	this->m_pLnods = new unsigned int [nelem*npoel];
	for(unsigned int jdiv=0;jdiv<m_div_y;jdiv++){
	for(unsigned int idiv=0;idiv<m_div_x;idiv++){
		const unsigned int ielem = idiv+jdiv*m_div_x;
		m_pLnods[ielem*npoel  ] = idiv  +(m_div_y+1)*jdiv;
		m_pLnods[ielem*npoel+1] = idiv+1+(m_div_y+1)*jdiv;
		m_pLnods[ielem*npoel+2] = idiv+1+(m_div_y+1)*(jdiv+1);
		m_pLnods[ielem*npoel+3] = idiv  +(m_div_y+1)*(jdiv+1);
	}
	}

	CElemSeg es;
	es.begin = 0;
	es.len = 4;
	es.max_noes = (m_div_x+1)*(m_div_y+1)-1;
	es.elseg_type = CORNER;
	es.id_na = 0;
	m_id_es = this->m_aSeg.AddObj( std::make_pair(1,es) );
	assert( this->IsSegID(m_id_es) );
}


void CElemAry_Rect::MakePattern_FEM(Fem::CIndexedArray& crs) const
{
	CElemAry::MakePattern_FEM(m_id_es,crs);
}
*/
