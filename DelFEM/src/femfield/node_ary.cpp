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
// NodeAry.cpp : 節点配列クラス(NodeAry)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
    #pragma warning ( disable : 4996 )
#endif

#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

#include "delfem/node_ary.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/zvector_blk.h"

using namespace Fem::Field;
using namespace MatVec;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CNodeAry::CNodeAry(const unsigned int size) : m_Size(size)
{
//	std::cout << "Construction of class CNodeAry : from size" << std::endl;
	m_DofSize = 0;
	m_paValue = 0;
}

CNodeAry::CNodeAry() : m_Size(0)
{
//	std::cout << "Construction of class CNodeAry : default" << std::endl;
	m_DofSize = 0;
	m_paValue = 0;
}

CNodeAry::CNodeAry(const CNodeAry& na){
//	m_str_name = na.m_str_name;
	m_Size = na.m_Size;
	m_DofSize = na.m_DofSize;
  {
    const unsigned int n = m_Size*m_DofSize;
    m_paValue = new double [n];
    for(unsigned int i=0;i<n;i++){ m_paValue[i] = na.m_paValue[i]; }
  }
	m_aSeg = na.m_aSeg;
	m_aEaEs = na.m_aEaEs;
}


CNodeAry::~CNodeAry()
{
//	std::cout << "Destruction of class CNodeAry" << std::endl;

	ClearSegment();
}

bool CNodeAry::ClearSegment(){
	if( m_paValue != 0 ){ 
		delete[] m_paValue; 
		m_paValue=0; 
	}
	m_aSeg.Clear();
	return true;
}

//////////////////////////////////////////////////////////////////////
// publicなメンバ関数
//////////////////////////////////////////////////////////////////////

const std::vector<int> CNodeAry::AddSegment( const std::vector< std::pair<unsigned int,CNodeSeg> >& id_seg_vec, const double& value)
{
	const unsigned int ndofn_pre = m_DofSize;
	const std::vector<int>& add_id_ary = this->AddSegment(id_seg_vec);
	const unsigned int ndofn_pos = m_DofSize;
	const unsigned int ndofn_add = ndofn_pos - ndofn_pre;

	const unsigned int nnode = m_Size;
	for(unsigned int inode=0;inode<nnode;inode++){
		for(unsigned int idof=0;idof<ndofn_add;idof++){
			m_paValue[inode*ndofn_pos+ndofn_pre+idof] = value;
		}
	}
	return add_id_ary;
}

const std::vector<int> CNodeAry::AddSegment( 
		const std::vector< std::pair<unsigned int,CNodeSeg> >& id_seg_vec, 
		const std::vector<double>& val_vec)
{
	unsigned int ndofn_add = 0;
	for(unsigned int iseg=0;iseg<id_seg_vec.size();iseg++){
		ndofn_add += id_seg_vec[iseg].second.len;
	}
	if( val_vec.size() != m_Size*ndofn_add ){
		std::cout << val_vec.size() << " " << m_Size << " " << ndofn_add << std::endl;
		std::cout << "Error!-->Not Mattching" << std::endl;
		assert(0);
		std::vector<int> tmp_id_ary;
		return tmp_id_ary;
	}

	const unsigned int ndofn_pre = m_DofSize;
	const std::vector<int>& add_id_ary = this->AddSegment(id_seg_vec);
	const unsigned int ndofn_pos = m_DofSize;
	assert( ndofn_pos - ndofn_pre == ndofn_add );

	const unsigned int nnode = m_Size;
	for(unsigned int inode=0;inode<nnode;inode++){
		for(unsigned int idof=0;idof<ndofn_add;idof++){
			m_paValue[inode*ndofn_pos+ndofn_pre+idof] = val_vec[inode*ndofn_add+idof];
		}
	}
	return add_id_ary;
}

const std::vector<int> CNodeAry::AddSegment( 
		const std::vector< std::pair<unsigned int,CNodeSeg> >& id_seg_vec )
{
	std::vector<int> add_id_ary;
	add_id_ary.resize(id_seg_vec.size());

	const unsigned int ndofval_begin = m_DofSize;
	unsigned int idofval0 = ndofval_begin;
	for(unsigned int iseg=0;iseg<id_seg_vec.size();iseg++){
		add_id_ary[iseg] = m_aSeg.AddObj(id_seg_vec[iseg]);
		assert( m_aSeg.IsObjID( add_id_ary[iseg] ) );
		CNodeSeg& ns = m_aSeg.GetObj( add_id_ary[iseg] );
		ns.idofval_begin = idofval0;
		idofval0 += id_seg_vec[iseg].second.len;
	}
	const unsigned int ndofval_end = idofval0;

	{
		double* paVal_old = m_paValue;
		double* paVal_new = new double [m_Size*ndofval_end];
		assert( paVal_new != 0 );
		for(unsigned int inode=0;inode<m_Size;inode++){
			for(unsigned int idofval=0;idofval<ndofval_begin;idofval++){
				paVal_new[inode*ndofval_end+idofval] = paVal_old[inode*ndofval_begin+idofval];
			}
		}
		m_DofSize = ndofval_end;
		delete[] m_paValue;
		m_paValue = paVal_new;
	}
	return add_id_ary;
}

//! ioffsetはvecの一つのblkにの何番目の自由度から始めるかを決める(e.g., シェルにおいて変位回転を一体化する場合)
bool CNodeAry::SetValueToNodeSegment(unsigned int id_ns,
                                     const MatVec::CVector_Blk& vec, unsigned int ioffset )
{
	if( !m_aSeg.IsObjID(id_ns) ) return false;
	const CNodeSeg& ns = m_aSeg.GetObj(id_ns);

	const unsigned int idofval_begin = ns.idofval_begin;
	const unsigned int ndofns = ns.len;
	const unsigned int nnode = this->Size();
	const unsigned int ndofval = this->m_DofSize;

    assert( vec.Len() >= 0 );   // 可変長ではない
    assert( vec.Len() >= (int)(ndofns+ioffset) );
	assert( vec.NBlk() == nnode );

	for(unsigned int inode=0;inode<nnode;inode++){
		for(unsigned int idofns=0;idofns<ndofns;idofns++){
            m_paValue[inode*ndofval+idofval_begin+idofns] = vec.GetValue(inode,idofns+ioffset);
		}
	}
	return true;
}

// ns0へalpha倍されたns1を加える
bool CNodeAry::AddValueToNodeSegment(unsigned int id_ns0, unsigned int id_ns1, double alpha)
{
	if( !this->IsSegID(id_ns0) ) return false;
	if( !this->IsSegID(id_ns1) ) return false;
	if( id_ns0 == id_ns1 ) return false;

	const CNodeSeg& ns0 = m_aSeg.GetObj(id_ns0);
	const unsigned int idofval_begin0 = ns0.idofval_begin;
	const CNodeSeg& ns1 = m_aSeg.GetObj(id_ns1);
	const unsigned int idofval_begin1 = ns1.idofval_begin;
	assert( idofval_begin0 != idofval_begin1 );
	assert( ns0.len == ns1.len );
	const unsigned int ndofns = ns0.len;

	const unsigned int nnode = this->Size();
	const unsigned int ndofn = m_DofSize;

	for(unsigned int inode=0;inode<nnode;inode++){
	for(unsigned int idof=0;idof<ndofns;idof++){
		m_paValue[inode*ndofn+idofval_begin0+idof] 
			+= alpha*m_paValue[inode*ndofn+idofval_begin1+idof];
	}
	}
	return true;
}

bool CNodeAry::AddValueToNodeSegment(unsigned int id_ns, const MatVec::CVector_Blk& vec, double alpha, unsigned int ioffset)
{
	if( !this->IsSegID(id_ns) ) return false;

	const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
	const unsigned int idofval_begin = ns.idofval_begin;
	const unsigned int ndofns = ns.len;

	const unsigned int nnode = this->Size();
	const unsigned int ndofn = m_DofSize;

    assert( vec.Len() >= 0 );   // 可変長配列は取り扱わない
    if( (int)(ndofns+ioffset) > vec.Len() ){
		std::cout << ndofns << " " << ioffset << " " << vec.Len() << std::endl;
		assert(0);
		return false;
	}
	if( vec.NBlk() != nnode ){
		assert(0);
		return false;
	}

	for(unsigned int inode=0;inode<nnode;inode++){
	for(unsigned int idof=0;idof<ndofns;idof++){
		m_paValue[inode*ndofn+idofval_begin+idof] += alpha*vec.GetValue(inode,idof+ioffset);
	}
	}
	return true;
}

bool CNodeAry::AddValueToNodeSegment(unsigned int id_ns, const MatVec::CZVector_Blk& vec, double alpha)
{
	if( !this->IsSegID(id_ns) ) return false;

	const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
	const unsigned int idofval_begin = ns.idofval_begin;
	const unsigned int ndofns = ns.len;
	const unsigned int nval = ndofns / 2;

	const unsigned int nnode = this->Size();
	const unsigned int ndofn = m_DofSize;

    assert( vec.BlkLen() >= 0 );     // 可変長配列は取り扱わない
    if( vec.BlkLen() != nval ) return false;
    if( vec.BlkVecLen() != nnode ) return false;

	for(unsigned int inode=0;inode<nnode;inode++){
	for(unsigned int ival=0;ival<nval;ival++){
		m_paValue[inode*ndofn+idofval_begin+ival*2  ] += alpha*vec.GetValue(inode,ival).Real();
		m_paValue[inode*ndofn+idofval_begin+ival*2+1] += alpha*vec.GetValue(inode,ival).Imag();
	}
	}
	return true;
}

bool CNodeAry::GetValueFromNodeSegment(unsigned int id_ns, MatVec::CVector_Blk& vec, unsigned int ioffset) const
{
	if( !this->IsSegID(id_ns) ) return false;

	const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
	const unsigned int idofval_begin = ns.idofval_begin;
	const unsigned int ndofns = ns.len;

	const unsigned int nnode = this->Size();
	const unsigned int ndofn = m_DofSize;

    assert( vec.Len() >= 0 );
    assert( vec.Len() >= (int)(ndofns+ioffset) );
	assert( vec.NBlk() == nnode );

	for(unsigned int inode=0;inode<nnode;inode++){
	for(unsigned int idof=0;idof<ndofns;idof++){
		 double val = m_paValue[inode*ndofn+idofval_begin+idof];
		 vec.SetValue(inode,idof+ioffset,val);
	}
	}
	return true;
}

bool CNodeAry::AddValueFromNodeSegment(double alpha, unsigned int id_ns, 
	MatVec::CVector_Blk& vec, unsigned int ioffset) const
{
	if( !this->IsSegID(id_ns) ) return false;

	const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
	const unsigned int idofval_begin = ns.idofval_begin;
	const unsigned int ndofns = ns.len;

	const unsigned int nnode = this->Size();
	const unsigned int ndofn = m_DofSize;

    assert( vec.Len() >= 0 );
    assert( vec.Len() >= (int)(ndofns+ioffset) );
	assert( vec.NBlk() == nnode );

	for(unsigned int inode=0;inode<nnode;inode++){
	for(unsigned int idof=0;idof<ndofns;idof++){
		 const double val = m_paValue[inode*ndofn+idofval_begin+idof];
		 vec.AddValue(inode,idof+ioffset,val*alpha);
	}
	}
	return true;
}

int CNodeAry::WriteToFile(const std::string& file_name, long& offset, unsigned int id ) const{

	FILE *fp;
	if( (fp = fopen(file_name.c_str(),"a"))== NULL ){
		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return 1;
	}

	fprintf(fp,"$$$$$$\n");
	fprintf(fp,"NODE\n");
	fprintf(fp,"%d\n",id);
	fprintf(fp,"%d %d\n",this->m_Size,this->m_DofSize);

	std::vector<unsigned int> id_ary_ns = this->m_aSeg.GetAry_ObjID();
	fprintf(fp,"%d\n",(int)id_ary_ns.size());
	for(unsigned int iid_ns=0;iid_ns<id_ary_ns.size();iid_ns++){
		unsigned int id_ns = id_ary_ns[iid_ns];
		const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
		fprintf(fp,"%d %d %d %s\n", iid_ns+1, id_ns, ns.len, ns.name.c_str());
	}

	const unsigned int ndofn = this->m_DofSize;
	for(unsigned int inode=0;inode<this->Size();inode++){
		fprintf(fp,"%d ",inode+1);
		for(unsigned int idofn=0;idofn<ndofn;idofn++){
			fprintf(fp,"%lf ",this->m_paValue[inode*ndofn+idofn]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return 0;
}

int CNodeAry::InitializeFromFile(const std::string& file_name, long& offset){

	FILE *fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];

	if( (fp = fopen(file_name.c_str(),"r"))== NULL ){
		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return false;
	}

	fseek(fp,offset,SEEK_SET);

	////////////////////////////////

	// read "$$$$$$"
	while( fgets(stmp1,buff_size,fp) != NULL ){
		if( stmp1[0] == '#' ) continue;
		if( strspn(stmp1," \n") != strlen(stmp1) ) break;
	}
	assert( strncmp(stmp1,"$$$$$$$$",6)==0 );

	// read "NODE"
	while( fgets(stmp1,buff_size,fp) != NULL ){
		if( stmp1[0] == '#' ) continue;
//		if( strspn(stmp1," \n") != strlen(stmp1) ) break;
		break;
	}
	assert( strncmp(stmp1,"NODE",4)==0 );

	// read id_na
	unsigned int id;
	{
		while( fgets(stmp1,buff_size,fp) != NULL ){
			if( stmp1[0] == '#' ) continue;
//			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			break;
		}
		int tmp_id;
		sscanf(stmp1,"%d",&tmp_id);
		assert( tmp_id > 0 );
		id = tmp_id;
	}
		
	// read nnode, ndofna
	unsigned int nnode, ndofn;
	{
		while( fgets(stmp1,buff_size,fp) != NULL ){
			if( stmp1[0] == '#' ) continue;
//			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			break;
		}
		int tmp_nnode, tmp_ndofn;
		sscanf(stmp1,"%d%d",&tmp_nnode, &tmp_ndofn);
		assert( tmp_nnode > 0 && tmp_ndofn > 0 );
		nnode = tmp_nnode;
		ndofn = tmp_ndofn;
	}

	// read nseg_na
	unsigned int nseg;
	{
		while( fgets(stmp1,buff_size,fp) != NULL ){
			if( stmp1[0] == '#' ) continue;
//			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			break;
		}
		int tmp_nseg;
		sscanf(stmp1,"%d",&tmp_nseg);
		assert( tmp_nseg > 0 );
		nseg = tmp_nseg;
	}

	std::vector< std::pair<unsigned int,CNodeSeg> > pair_id_ns_ary;
	{
		while( fgets(stmp1,buff_size,fp) != NULL ){
			if( stmp1[0] == '#' ) continue;
//			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			break;
		}
		int tmp_iseg, tmp_id, tmp_len;
		char tmp_name[256];
		for(unsigned int iseg=0;;iseg++){
			sscanf(stmp1,"%d%d%d%s",&tmp_iseg, &tmp_id, &tmp_len, tmp_name);
            assert( tmp_iseg-1 == (int)iseg );
			assert( tmp_id > 0 );
			assert( tmp_len > 0 );
			pair_id_ns_ary.push_back( std::make_pair(tmp_id, CNodeSeg(tmp_len, tmp_name)) );
			if( iseg >= nseg-1 ) break;
			fgets(stmp1,buff_size,fp);
		}
	}
		
	std::vector<double> val_vec;
	//CVector_Blk val_vec(nnode,ndofn);
	{
		while( fgets(stmp1,buff_size,fp) != NULL ){
			if( stmp1[0] == '#' ) continue;
//			if( strspn(stmp1," \n") != strlen(stmp1) ) break;
			break;
		}
		char *tp, *e;
		int tmp_inode;
		double dtmp1;
		val_vec.resize(nnode*ndofn);
		for(unsigned int inode=0;;inode++){
			tp = strtok(stmp1," ");
			tmp_inode = atoi(tp);
            assert( tmp_inode-1 == (int)inode );
			for(unsigned int idofn=0;idofn<ndofn;idofn++){
				tp = strtok(NULL," ");
				dtmp1 = strtod(tp,&e);
				assert( *e == '\0' || *e == '\n' );
				val_vec[inode*ndofn+idofn] = dtmp1;
//				val_vec.SetValue(inode,idofn,dtmp1);
			}
			if( inode >= nnode-1 ) break;
			fgets(stmp1,buff_size,fp);
		}
	}

	offset = ftell(fp);
	fclose(fp);

	this->m_aSeg.Clear();
	if( m_paValue != 0 ) delete m_paValue;
	this->m_DofSize = 0;

	this->m_Size = nnode;
	std::vector<int> add_id_ary = this->AddSegment(pair_id_ns_ary,val_vec);

	assert( add_id_ary.size() == pair_id_ns_ary.size() );
	for(unsigned int ins=0;ins<pair_id_ns_ary.size();ins++){
        if( add_id_ary[ins] != (int)pair_id_ns_ary[ins].first ){
			std::cout << add_id_ary[ins] <<  " " << pair_id_ns_ary[ins].first << std::endl;
			assert(0);
		}
	}

	return 0;
}

/*
int CNodeAry::DumpToFile_UpdatedValue(const std::string& file_name, long& offset,unsigned int id) const
{
	
	bool is_update = false;
	std::vector<int> is_update_flag;
	{	// UpdeteされたDOFのフラグを立てる
		is_update_flag.resize(m_DofSize,0);
		std::vector<unsigned int> id_ns_ary = this->m_aSeg.GetAry_ObjID();
		for(unsigned int iid=0;iid<id_ns_ary.size();iid++){
			unsigned int id_ns = id_ns_ary[iid];
			const CNodeSeg& ns = m_aSeg.GetObj(id_ns);
			if( !ns.is_updated ) continue;
			assert( ns.idofval_begin+ns.len <= m_DofSize );
			for(unsigned int idof=0;idof<ns.len;idof++){
				is_update_flag[ns.idofval_begin+idof] = 1;
				is_update = true;
			}
		}
	}

	// Updateされていないなら戻る
	if( !is_update ) return false;
	
	FILE *fp;
	if( (fp = fopen(file_name.c_str(),"a"))== NULL ){
		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return 1;
	}

	fprintf(fp,"NODE\n");
	for(unsigned int idof=0;idof<m_DofSize;idof++){
		if( is_update_flag[idof] == 1 ){ fprintf(fp,"%d ",idof); }
	}
	fprintf(fp,"\n");
	for(unsigned int inode=0;inode<this->m_Size;inode++){
		for(unsigned int idof=0;idof<m_DofSize;idof++){
			if( is_update_flag[idof] == 1 ){
				double val = this->m_paValue[inode*m_DofSize+idof];
				fprintf(fp,"%lf ",val);
			}
		}
		fprintf(fp,"\n");
	}

	offset = ftell(fp);
	fclose(fp);

	return 0;
}
*/


void CNodeAry::AddEaEs( std::pair<unsigned int, unsigned int> eaes ){
  unsigned int ieaes = this->GetIndEaEs( eaes );
  if( ieaes != m_aEaEs.size() ) return;
  {
    CEaEsInc eaesinc;
    eaesinc.id_ea = eaes.first;
    eaesinc.id_es = eaes.second;
    m_aEaEs.push_back( eaesinc );
  }
}

std::vector< std::pair<unsigned int, unsigned int> > CNodeAry::GetAryEaEs() const {
  std::vector< std::pair<unsigned int, unsigned int> > aEaEs;
  for(unsigned int ieaes=0;ieaes<m_aEaEs.size();ieaes++){
    aEaEs.push_back( std::make_pair(m_aEaEs[ieaes].id_ea,m_aEaEs[ieaes].id_es) );
  }
  return aEaEs;
}

void CNodeAry::SetIncludeEaEs_InEaEs( std::pair<unsigned int, unsigned int> eaes_included,
                           std::pair<unsigned int, unsigned int> eaes_container )
{
  unsigned int ieaes_included = this->GetIndEaEs( eaes_included );
  if( ieaes_included == m_aEaEs.size() ) return;
  unsigned int ieaes_container = this->GetIndEaEs( eaes_container );
  if( ieaes_container == m_aEaEs.size() ) return;
  //        std::cout << "SetInclude " << eaes_included.first << " " << eaes_included.second << "    in    ";
  //        std::cout << eaes_container.first << " " << eaes_container.second << std::endl;
  m_aEaEs[ieaes_container].aIndEaEs_Include.push_back(ieaes_included);
}


bool CNodeAry::IsIncludeEaEs_InEaEs( std::pair<unsigned int, unsigned int> eaes_inc,
                          std::pair<unsigned int, unsigned int> eaes_in ) const
{
  unsigned int ieaes_inc = this->GetIndEaEs( eaes_inc );
  if( ieaes_inc == m_aEaEs.size() ) return false;
  unsigned int ieaes_in = this->GetIndEaEs( eaes_in );
  if( ieaes_in == m_aEaEs.size() ) return false;
  const std::vector<unsigned int>& inc = m_aEaEs[ieaes_in].aIndEaEs_Include;
  for(unsigned int iieaes=0;iieaes<inc.size();iieaes++){
    if( inc[iieaes] == ieaes_inc ){
      return true;
    }
  }
  return false;
}

std::vector< std::pair<unsigned int, unsigned int> > CNodeAry::GetAry_EaEs_Min() const
{
  std::vector< std::pair<unsigned int, unsigned int> > aEaEs;
  std::vector<unsigned int> aflg;
  aflg.resize(m_aEaEs.size(),1);
  for(unsigned int ieaes=0;ieaes<m_aEaEs.size();ieaes++){
    const std::vector<unsigned int>& inc = m_aEaEs[ieaes].aIndEaEs_Include;
    for(unsigned int jjeaes=0;jjeaes<inc.size();jjeaes++){
      unsigned int jeaes = inc[jjeaes];
      assert( jeaes < m_aEaEs.size() );
      aflg[jeaes] = 0;
    }
  }
  for(unsigned int ieaes=0;ieaes<m_aEaEs.size();ieaes++){
    if( aflg[ieaes] == 1 ){
      aEaEs.push_back( std::make_pair(m_aEaEs[ieaes].id_ea,m_aEaEs[ieaes].id_es) );
    }
  }
  return aEaEs;
}

unsigned int CNodeAry::IsContainEa_InEaEs(std::pair<unsigned int, unsigned int>eaes, unsigned int id_ea) const
{
  const unsigned int ieaes = this->GetIndEaEs( eaes );
  if( ieaes == m_aEaEs.size() ) return 0;
  if( m_aEaEs[ieaes].id_ea == id_ea ){
    return m_aEaEs[ieaes].id_es;
  }
  const std::vector<unsigned int>& inc = m_aEaEs[ieaes].aIndEaEs_Include;
  for(unsigned int jjeaes=0;jjeaes<inc.size();jjeaes++){
    unsigned int jeaes = inc[jjeaes];
    assert( jeaes < m_aEaEs.size() );
    if( m_aEaEs[jeaes].id_ea == id_ea ){
      return m_aEaEs[jeaes].id_es;
    }
  }
  return 0;
}


