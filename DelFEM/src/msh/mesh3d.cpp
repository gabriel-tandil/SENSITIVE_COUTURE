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

#if defined(__VISUALC__)
#pragma warning(disable: 4786)
#pragma warning(disable: 4996)
#endif
#define for if(0);else for

#include <stdio.h>
#include <set>
#include <vector>
#include <queue>
#include <cassert>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "delfem/msh/meshkernel2d.h"
#include "delfem/msh/meshkernel3d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"

using namespace Com;

void Msh::CMesh3D::Clear(){
	this->aVec.clear();

	this->m_aVertex.clear();
	this->m_aBarAry.clear();
	this->m_aTriAry.clear();
	this->m_aQuadAry.clear();
	this->m_aTetAry.clear();
	this->m_aHexAry.clear();

	this->m_ElemLoc.clear();
	this->m_ElemType.clear();
	this->m_include_relation.clear();
}


unsigned int Msh::CMesh3D::FindMaxID() const
{
	unsigned int max_id = 0;
	{	// 要素IDの最大値を求める
		for(unsigned int ia=0;ia<m_aVertex.size();ia++){
			if( max_id < m_aVertex[ia].id ) max_id = m_aVertex[ia].id;
		}
		for(unsigned int ia=0;ia<m_aBarAry.size();ia++){
			if( max_id < m_aBarAry[ia].id ) max_id = m_aBarAry[ia].id;
		}
		for(unsigned int ia=0;ia<m_aTriAry.size();ia++){
			if( max_id < m_aTriAry[ia].id ) max_id = m_aTriAry[ia].id;
		}
		for(unsigned int ia=0;ia<m_aQuadAry.size();ia++){
			if( max_id < m_aQuadAry[ia].id ) max_id = m_aQuadAry[ia].id;
		}
		for(unsigned int ia=0;ia<m_aTetAry.size();ia++){
			if( max_id < m_aTetAry[ia].id ) max_id = m_aTetAry[ia].id;
		}
		for(unsigned int ia=0;ia<m_aHexAry.size();ia++){
			if( max_id < m_aHexAry[ia].id ) max_id = m_aHexAry[ia].id;
		}
	}
	return max_id;
}

unsigned int Msh::CMesh3D::GetFreeObjID()
{
	unsigned int max_id = this->FindMaxID();
	std::vector<unsigned int> is_used_flg_ary;
	{	// このIDが使われているかどうかを示すハッシュを作る
		is_used_flg_ary.resize( max_id+1, 0 );
		for(unsigned int ia=0;ia<m_aVertex.size();ia++){
			if( m_aVertex[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aVertex[ia].id ] == 0 );
			assert( m_aVertex[ia].id <= max_id );
			is_used_flg_ary[ m_aVertex[ia].id ] = 1;
		}
		for(unsigned int ia=0;ia<m_aBarAry.size();ia++){
			if( m_aBarAry[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aBarAry[ia].id ] == 0 );
			assert( m_aBarAry[ia].id <= max_id );
			is_used_flg_ary[ m_aBarAry[ia].id ] = 1;
		}
		for(unsigned int ia=0;ia<m_aTriAry.size();ia++){
			if( m_aTriAry[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aTriAry[ia].id ] == 0 );
			assert( m_aTriAry[ia].id <= max_id );
			is_used_flg_ary[ m_aTriAry[ia].id ] = 1;
		}
		for(unsigned int ia=0;ia<m_aQuadAry.size();ia++){
			if( m_aQuadAry[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aQuadAry[ia].id ] == 0 );
			assert( m_aQuadAry[ia].id <= max_id );
			is_used_flg_ary[ m_aQuadAry[ia].id ] = 1;
		}
		for(unsigned int ia=0;ia<m_aTetAry.size();ia++){
			if( m_aTetAry[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aTetAry[ia].id ] == 0 );
			assert( m_aTetAry[ia].id <= max_id );
			is_used_flg_ary[ m_aTetAry[ia].id ] = 1;
		}
		for(unsigned int ia=0;ia<m_aHexAry.size();ia++){
			if( m_aHexAry[ia].id == 0 ) continue;
			assert( is_used_flg_ary[ m_aHexAry[ia].id ] == 0 );
			assert( m_aHexAry[ia].id <= max_id );
			is_used_flg_ary[ m_aHexAry[ia].id ] = 1;
		}
	}
	for(unsigned int i=1;i<is_used_flg_ary.size();i++){
		if( is_used_flg_ary[i] == 0 ){ return i; }
	}
	return max_id+1;
}

bool Msh::CMesh3D::IsID(unsigned int id) const
{
	assert( this->m_ElemLoc.size() > id );
	if( this->m_ElemLoc.size() <= id ) return false;
	const int iloc = this->m_ElemLoc[id];
	if( iloc == -1 ) return false;

	// 以下assertルーティン
	assert( this->m_ElemType.size() > id );
	const int itype = this->m_ElemType[id];
	assert( itype >= 0 );
	assert( iloc >= 0 );
	if( itype == 0 ){
		assert( m_aVertex.size() > (unsigned int)iloc );
		assert( m_aVertex[iloc].id == id );
	}
	else if( itype == 1 ){
		assert( m_aBarAry.size() > (unsigned int)iloc );
		assert( m_aBarAry[iloc].id == id );
	}
	else if( itype == 2 ){
		assert( m_aTriAry.size() > (unsigned int)iloc );
		assert( m_aTriAry[iloc].id == id );
	}
	else if( itype == 3 ){
		assert( m_aQuadAry.size() > (unsigned int)iloc );
		assert( m_aQuadAry[iloc].id == id );
	}
	else if( itype == 4 ){
		assert( m_aTetAry.size() > (unsigned int)iloc );
		assert( m_aTetAry[iloc].id == id );
	}
	else if( itype == 5 ){
		assert( m_aHexAry.size() > (unsigned int)iloc );
		assert( m_aHexAry[iloc].id == id );
	}
	return true;
}


// 要素の場所と種類をハッシュ（IDが引数）している配列を初期化
void Msh::CMesh3D::MakeElemLocationType()
{
	unsigned int max_id = this->FindMaxID();
	////////////////
	this->m_ElemLoc.clear();
	this->m_ElemLoc.resize(max_id+1,-1);
	this->m_ElemType.clear();
	this->m_ElemType.resize(max_id+1,-1);
	////////////////
	for(unsigned int iver=0;iver<m_aVertex.size();iver++){
		unsigned int id0 = m_aVertex[iver].id;
		m_ElemLoc[id0] = iver;
		m_ElemType[id0] = 0;
	}
	for(unsigned int ibarary=0;ibarary<m_aBarAry.size();ibarary++){
		unsigned int id0 = m_aBarAry[ibarary].id;
		m_ElemLoc[id0] = ibarary;
		m_ElemType[id0] = 1;
	}
	for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
		unsigned int id0 = m_aTriAry[itriary].id;
		m_ElemLoc[id0] = itriary;
		m_ElemType[id0] = 2;
	}
	for(unsigned int ia=0;ia<m_aQuadAry.size();ia++){
		unsigned int id0 = m_aQuadAry[ia].id;
		m_ElemLoc[id0] = ia;
		m_ElemType[id0] = 3;
	}
	for(unsigned int itetary=0;itetary<m_aTetAry.size();itetary++){
		unsigned int id0 = m_aTetAry[itetary].id;
		m_ElemLoc[id0] = itetary;
		m_ElemType[id0] = 4;
	}
	for(unsigned int ihexary=0;ihexary<m_aHexAry.size();ihexary++){
		unsigned int id0 = m_aHexAry[ihexary].id;
		m_ElemLoc[id0] = ihexary;
		m_ElemType[id0] = 5;
	}
}


Msh::MSH_TYPE Msh::CMesh3D::GetConnectivity(unsigned int id_msh, std::vector<int>& lnods) const
{
	assert( this->IsID(id_msh) );
	MSH_TYPE msh_type;
	unsigned int nnoel, nelem;
	const unsigned int itype = m_ElemType[id_msh];
	const unsigned int iloc = m_ElemLoc[id_msh];
	if(      itype == 0 ){ 
		msh_type = VERTEX; nnoel = 1; nelem = 1;
		lnods.resize(nnoel*nelem);
		lnods[0] = m_aVertex[iloc].v;
	}
	else if( itype == 1 ){ 
		msh_type = BAR;    nnoel = 2; 
		const std::vector<SBar>& aBar = m_aBarAry[iloc].m_aBar;
		nelem = aBar.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*2  ] = aBar[ielem].v[0];
			lnods[ielem*2+1] = aBar[ielem].v[1];
		}
	}
	else if( itype == 2 ){ 
		msh_type = TRI;    nnoel = 3;
		const std::vector<STri3D>& aTri = m_aTriAry[iloc].m_aTri;
		nelem = aTri.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*3  ] = aTri[ielem].v[0];
			lnods[ielem*3+1] = aTri[ielem].v[1];
			lnods[ielem*3+2] = aTri[ielem].v[2];
		}
	}
	else if( itype == 3 ){ 
		msh_type = QUAD;   nnoel = 4;
		const std::vector<SQuad3D>& aQuad = m_aQuadAry[iloc].m_aQuad;
		nelem = aQuad.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*4  ] = aQuad[ielem].v[0];
			lnods[ielem*4+1] = aQuad[ielem].v[1];
			lnods[ielem*4+2] = aQuad[ielem].v[2];
			lnods[ielem*4+3] = aQuad[ielem].v[3];
		}
	}
	else if( itype == 4 ){ 
		msh_type = TET;   nnoel = 4;
		const std::vector<STet>& aTet = m_aTetAry[iloc].m_aTet;
		nelem = aTet.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*4  ] = aTet[ielem].v[0];
			lnods[ielem*4+1] = aTet[ielem].v[1];
			lnods[ielem*4+2] = aTet[ielem].v[2];
			lnods[ielem*4+3] = aTet[ielem].v[3];
		}
	}
	else if( itype == 5 ){ 
		msh_type = HEX;   nnoel = 8;
		const std::vector<SHex>& aHex = m_aHexAry[iloc].m_aHex;
		nelem = aHex.size();
		lnods.resize(nnoel*nelem);
		for(unsigned int ielem=0;ielem<nelem;ielem++){
			lnods[ielem*8  ] = aHex[ielem].v[0];
			lnods[ielem*8+1] = aHex[ielem].v[1];
			lnods[ielem*8+2] = aHex[ielem].v[2];
			lnods[ielem*8+3] = aHex[ielem].v[3];
			lnods[ielem*8+4] = aHex[ielem].v[4];
			lnods[ielem*8+5] = aHex[ielem].v[5];
			lnods[ielem*8+6] = aHex[ielem].v[6];
			lnods[ielem*8+7] = aHex[ielem].v[7];
		}
	}
	else{ assert(0); }
	return msh_type;
}


bool Msh::CMesh3D::Serialize( Com::CSerializer& arch )
{
	if( arch.IsLoading() ){	// 読み込み時の処理
		this->Clear();
		char stmp1[256];
		arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$$$$$",8)==0 );
		arch.Get("%s",stmp1);	assert( strncmp(stmp1,"CMesh3D",   5)==0 );
		int naVer, naBar, naTri, naQuad, naTet, naHex;
		arch.Get("%d%d%d%d%d%d",&naVer,&naBar,&naTri,&naQuad,&naTet,&naHex);
		{
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );	
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"VEC3D",5)==0 );
			int nvec, ndim;	arch.Get("%d%d",&nvec,&ndim);	assert(nvec>0 && (ndim>0&&ndim<4) );
			this->aVec.resize(nvec);
			for(unsigned int ivec=0;ivec<(unsigned int)nvec;ivec++){
				int itmp0;
				double x,y,z;
				arch.Get("%d%lf%lf%lf",&itmp0,&x,&y,&z);// assert( itmp0 == ivec );
				aVec[ivec].x = x;
				aVec[ivec].y = y;
				aVec[ivec].z = z;
			}
		}
		for(unsigned int iaVer=0;iaVer<(unsigned int)naVer;iaVer++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"VER",  3)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int iv;		arch.Get("%d",&iv);		assert(iv>=0);
			const unsigned int iver = m_aVertex.size();
			m_aVertex.resize( m_aVertex.size()+1 );
			SVertex3D& Ver = m_aVertex[iver];
			Ver.id = id;
			Ver.id_cad = id_cad;
			Ver.v = iv;
		}
		for(unsigned int iaBar=0;iaBar<(unsigned int)naBar;iaBar++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"BAR",  3)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int nbar;	arch.Get("%d",&nbar);	assert(nbar>0);
			const unsigned int ibar_ary = m_aBarAry.size();
			m_aBarAry.resize( m_aBarAry.size()+1 );
			CBarAry3D& aBar = m_aBarAry[ibar_ary];
			aBar.id = id;
			aBar.id_cad = id_cad;
			aBar.m_aBar.resize(nbar);
			for(int ibar=0;ibar<nbar;ibar++){
				int tmp_ibar,iv0,iv1;
				arch.Get("%d%d%d",&tmp_ibar,&iv0,&iv1);
				assert( tmp_ibar == ibar );
				assert( iv0>=0 && iv1>=0 );
				aBar.m_aBar[ibar].v[0] = iv0;
				aBar.m_aBar[ibar].v[1] = iv1;
			}
		}
		for(unsigned int iaTri=0;iaTri<(unsigned int)naTri;iaTri++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"TRI",  3)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int ntri;	arch.Get("%d",&ntri);	assert(ntri>0);
			const unsigned int itri_ary = m_aTriAry.size();
			m_aTriAry.resize( m_aTriAry.size()+1 );
			CTriAry3D& aTri = m_aTriAry[itri_ary];
			aTri.id = id;
			aTri.m_aTri.resize(ntri);
			for(int itri=0;itri<ntri;itri++){
				int tmp_itri,iv0,iv1,iv2;
				arch.Get("%d %d %d %d",&tmp_itri,&iv0,&iv1,&iv2);
				assert( tmp_itri == itri );
				assert( iv0>=0 && iv1>=0 && iv2>=0 );
				aTri.m_aTri[itri].v[0] = iv0;
				aTri.m_aTri[itri].v[1] = iv1;
				aTri.m_aTri[itri].v[2] = iv2;
			}
		}
		for(unsigned int iaQuad=0;iaQuad<(unsigned int)naQuad;iaQuad++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"QUAD", 4)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int nquad;	arch.Get("%d",&nquad);	assert(nquad>0);
			const unsigned int iquad_ary = m_aQuadAry.size();
			m_aQuadAry.resize( m_aQuadAry.size()+1 );
			CQuadAry3D& aQuad = m_aQuadAry[iquad_ary];
			aQuad.id = id;
			aQuad.m_aQuad.resize(nquad);
			for(int iquad=0;iquad<nquad;iquad++){
				int tmp_iquad,iv0,iv1,iv2,iv3;
				arch.Get("%d %d %d %d %d",&tmp_iquad,&iv0,&iv1,&iv2,&iv3);
//				assert( tmp_iquad == iquad );
				assert( iv0>=0 && iv1>=0 && iv2>=0 && iv3>=0 );
				aQuad.m_aQuad[iquad].v[0] = iv0;
				aQuad.m_aQuad[iquad].v[1] = iv1;
				aQuad.m_aQuad[iquad].v[2] = iv2;
				aQuad.m_aQuad[iquad].v[3] = iv3;
			}
		}
		for(unsigned int iaTet=0;iaTet<(unsigned int)naTet;iaTet++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"TET",  3)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int ntet;	arch.Get("%d",&ntet);	assert(ntet>0);
			const unsigned int itet_ary = m_aQuadAry.size();
			m_aTetAry.resize( m_aTetAry.size()+1 );
			CTetAry& aTet = m_aTetAry[itet_ary];
			aTet.id = id;
			aTet.m_aTet.resize(ntet);
			for(int itet=0;itet<ntet;itet++){
				int tmp_itet,iv0,iv1,iv2,iv3;
				arch.Get("%d %d %d %d %d",&tmp_itet,&iv0,&iv1,&iv2,&iv3);
				assert( tmp_itet == itet );
				assert( iv0>=0 && iv1>=0 && iv2>=0 && iv3>=0 );
				aTet.m_aTet[itet].v[0] = iv0;
				aTet.m_aTet[itet].v[1] = iv1;
				aTet.m_aTet[itet].v[2] = iv2;
				aTet.m_aTet[itet].v[3] = iv3;
			}
		}
		for(unsigned int iaHex=0;iaHex<(unsigned int)naHex;iaHex++){
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
			arch.Get("%s",stmp1);	assert( strncmp(stmp1,"HEX",  4)==0 );
			int id;		arch.Get("%d",&id);		assert(id>0);
			int id_cad;	arch.Get("%d",&id_cad);	assert(id_cad>=0);
			int nhex;	arch.Get("%d",&nhex);	assert(nhex>0);
			const unsigned int ihex_ary = m_aHexAry.size();
			m_aHexAry.resize( m_aHexAry.size()+1 );
			CHexAry& aHex = m_aHexAry[ihex_ary];
			aHex.id = id;
			aHex.m_aHex.resize(nhex);
			for(int ihex=0;ihex<nhex;ihex++){
				int tmp_ihex,iv0,iv1,iv2,iv3,iv4,iv5,iv6,iv7;
				arch.Get("%d%d%d%d%d%d%d%d%d",&tmp_ihex,&iv0,&iv1,&iv2,&iv3,&iv4,&iv5,&iv6,&iv7);
//				assert( tmp_ihex == ihex );
				assert( iv0>=0 && iv1>=0 && iv2>=0 && iv3>=0 && iv4>=0 && iv5>=0 && iv6>=0 && iv7>=0 );
				aHex.m_aHex[ihex].v[0] = iv0;
				aHex.m_aHex[ihex].v[1] = iv1;
				aHex.m_aHex[ihex].v[2] = iv2;
				aHex.m_aHex[ihex].v[3] = iv3;
				aHex.m_aHex[ihex].v[4] = iv4;
				aHex.m_aHex[ihex].v[5] = iv5;
				aHex.m_aHex[ihex].v[6] = iv6;
				aHex.m_aHex[ihex].v[7] = iv7;
			}
		}
		this->MakeElemLocationType();
        ////////////////////////////////
        m_include_relation.resize( this->FindMaxID()+1 );
	    arch.Get("%s",stmp1);	assert( strncmp(stmp1,"$$$$", 4)==0 );
	    arch.Get("%s",stmp1);	assert( strncmp(stmp1,"INC_REL", 7)==0 );
		int nrel;		arch.Get("%d",&nrel);		assert(nrel>=0);
        for(unsigned int irel=0;irel<(unsigned int)nrel;irel++){
            char stmp1[256];
            int itmp1, id_msh=-1;
            arch.GetLine(stmp1,256);
//            std::cout << "irel : " << stmp1 << std::endl;
            // スペースを区切りに文字列を抽出
            std::vector<unsigned int> id_msh_inc_ary;
            unsigned int icnt = 0;
            char *tp;
            tp = strtok( stmp1, " \n" );
            while ( tp != NULL ) {
                sscanf(tp,"%d",&itmp1);
                if( icnt == 0 ){ assert( itmp1 == (int)irel ); }
                else if( icnt == 1 ){ id_msh = itmp1; }
                else{ id_msh_inc_ary.push_back(itmp1); }
                icnt++;
                tp = strtok( NULL," \n" );
                if( tp == NULL ) break;
            }
			assert( id_msh != -1 );
/*			{
				std::cout << irel << "  " << id_msh << " : ";
				for(unsigned int i=0;i<id_msh_inc_ary.size();i++){
					std::cout << id_msh_inc_ary[i] << " ";
				}
				std::cout << std::endl;
			}*/
            m_include_relation[ id_msh ] = id_msh_inc_ary;
        }
		return true;
	}
	else{ // 書き込み時の処理	
		{	// ヘッダを書き込む
			arch.Out("$$$$$$$$\n");
			arch.Out("CMesh2D\n");
			arch.Out("%d %d %d %d %d %d\n",m_aVertex.size(), m_aBarAry.size(), m_aTriAry.size(), m_aQuadAry.size(), m_aTetAry.size(), m_aHexAry.size());
		}
		{	// Vector2Dの出力
			arch.Out("$$$$\n");
			arch.Out("VEC3D\n");
			arch.Out("%d %d\n",aVec.size(),3);
			for(unsigned int ivec=0;ivec<aVec.size();ivec++){
				arch.Out("%d %lf %lf %lf\n",ivec,aVec[ivec].x,aVec[ivec].y,aVec[ivec].z);
			}
		}
		{	// Vertexの出力
			for(unsigned int iver=0;iver<m_aVertex.size();iver++){
				arch.Out("$$$$\n");
				arch.Out("VER\n");
				arch.Out("%d\n",m_aVertex[iver].id);
				arch.Out("%d\n",m_aVertex[iver].id_cad);
				arch.Out("%d\n",m_aVertex[iver].v);
			}
		}
		{	// Barの出力
			for(unsigned int ibar_ary=0;ibar_ary<m_aBarAry.size();ibar_ary++){
				arch.Out("$$$$\n");
				arch.Out("BAR\n");
				arch.Out("%d\n",m_aBarAry[ibar_ary].id);
				arch.Out("%d\n",m_aBarAry[ibar_ary].id_cad);
				arch.Out("%d\n",m_aBarAry[ibar_ary].m_aBar.size());
				const std::vector<SBar>& aBar = m_aBarAry[ibar_ary].m_aBar;
				for(unsigned int ibar=0;ibar<aBar.size();ibar++){
					arch.Out("%d %d %d\n",ibar,aBar[ibar].v[0],aBar[ibar].v[1]);
				}
			}
		}
		{	// Triの出力
			for(unsigned int itri_ary=0;itri_ary<m_aTriAry.size();itri_ary++){
				arch.Out("$$$$\n");
				arch.Out("TRI\n");
				arch.Out("%d\n",m_aTriAry[itri_ary].id);
				arch.Out("%d\n",0);
				arch.Out("%d\n",m_aTriAry[itri_ary].m_aTri.size());
				const std::vector<STri3D>& aTri = m_aTriAry[itri_ary].m_aTri;
				for(unsigned int itri=0;itri<aTri.size();itri++){
					arch.Out("%d %d %d %d\n",itri,aTri[itri].v[0],aTri[itri].v[1],aTri[itri].v[2]);
				}
			}
		}
		{	// Quadの出力
			for(unsigned int iquad_ary=0;iquad_ary<m_aQuadAry.size();iquad_ary++){
				arch.Out("$$$$\n");
				arch.Out("QUAD\n");
				arch.Out("%d\n",m_aQuadAry[iquad_ary].id);
				arch.Out("%d\n",0);
				arch.Out("%d\n",m_aQuadAry[iquad_ary].m_aQuad.size());
				const std::vector<SQuad3D>& aQuad = m_aQuadAry[iquad_ary].m_aQuad;
				for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
					arch.Out("%d %d %d %d %d\n",iquad,aQuad[iquad].v[0],aQuad[iquad].v[1],aQuad[iquad].v[2],aQuad[iquad].v[3]);
				}
			}
		}
		{	// Tetの出力
			for(unsigned int itet_ary=0;itet_ary<m_aTetAry.size();itet_ary++){
				arch.Out("$$$$\n");
				arch.Out("TET\n");
				arch.Out("%d\n",m_aTetAry[itet_ary].id);
				arch.Out("%d\n",0);
				arch.Out("%d\n",m_aTetAry[itet_ary].m_aTet.size());
				const std::vector<STet>& aTet = m_aTetAry[itet_ary].m_aTet;
				for(unsigned int itet=0;itet<aTet.size();itet++){
					arch.Out("%d %d %d %d %d\n",itet,aTet[itet].v[0],aTet[itet].v[1],aTet[itet].v[2],aTet[itet].v[3]);
				}
			}
		}
		{	// Hexの出力
			for(unsigned int ihex_ary=0;ihex_ary<m_aHexAry.size();ihex_ary++){
				arch.Out("$$$$\n");
				arch.Out("HEX\n");
				arch.Out("%d\n",m_aHexAry[ihex_ary].id);
				arch.Out("%d\n",0);
				arch.Out("%d\n",m_aHexAry[ihex_ary].m_aHex.size());
				const std::vector<SHex>& aHex = m_aHexAry[ihex_ary].m_aHex;
				for(unsigned int ihex=0;ihex<aHex.size();ihex++){
					arch.Out("%d %d %d %d %d %d %d %d %d\n",ihex,
						aHex[ihex].v[0],aHex[ihex].v[1],aHex[ihex].v[2],aHex[ihex].v[3],
						aHex[ihex].v[4],aHex[ihex].v[5],aHex[ihex].v[6],aHex[ihex].v[7] );
				}
			}
		}
		{	// 包含関係の出力
			arch.Out("$$$$\n");
			arch.Out("INC_REL\n");
			unsigned int nrel = 0;
			for(unsigned int id_msh=0;id_msh<m_include_relation.size();id_msh++){
				if( !m_include_relation[id_msh].empty() ) nrel++;
			}
			arch.Out("%d\n",nrel);
			char stmp[256], stmp1[256];
			unsigned int irel = 0;
			for(unsigned int id_msh=0;id_msh<m_include_relation.size();id_msh++){
				if( m_include_relation[id_msh].empty() ){ continue; }
				sprintf(stmp,"%d",irel); irel++;
				sprintf(stmp1," %d",id_msh);
				strcat(stmp,stmp1);
				for(unsigned int i=0;i<m_include_relation[irel].size();i++){
					sprintf(stmp1," %d",m_include_relation[irel][i]);
					strcat(stmp,stmp1);
				}
				strcat(stmp,"\n");
//				std::cout << stmp << std::endl;
				arch.Out(stmp);
			}
			assert( irel == nrel );
		}
	}
	return true;
}


bool Msh::CMesh3D::ReadFromFile_GiDMsh(const std::string& file_name)
{
	FILE *fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];
	char stmp2[64], stmp3[64], stmp4[64], stmp6[64];

	if( (fp = ::fopen(file_name.c_str(),"r"))== NULL ){
		std::cout << "Cannot Open Gid File : " << file_name << std::endl;
		return false;
	}

	char str_elem_type[64];
	int ndim, nnoel;
	fgets(stmp1,buff_size,fp);
	sscanf(stmp1,"%s %s %d %s %s %s %d", 
		stmp2,stmp3,&ndim,stmp4,str_elem_type,stmp6,&nnoel);

	if( strcmp(str_elem_type,"Hexahedra")==0 ){
		if( nnoel != 8 ) return false;
	}
	else if( strcmp(str_elem_type,"Tetrahedra")==0 ){
		if( nnoel != 4 ) return false;
	}
	else if( strcmp(str_elem_type,"Triangle")==0 ){
		if( nnoel != 3 ) return false;
	}
	else{
		std::cout << "四面体と六面体と三角形以外の要素は未対応：" << str_elem_type << std::endl;
		return false;
	}

	////////////////////////////////
	// ここからは内容を変更する

	this->aVec.clear();
	this->m_aBarAry.clear();
	this->m_aTriAry.clear();
	this->m_aQuadAry.clear();
	this->m_aTetAry.clear();
	this->m_aHexAry.clear();
	this->m_aVertex.clear();
	this->m_ElemLoc.clear();
	this->m_ElemType.clear();


	{	// 節点の読み込み
		fgets(stmp1,buff_size,fp);	// Coorinates
		std::vector<double> tmp_buffer;
		tmp_buffer.reserve(16384);	// 2^14
		unsigned int counter = 0;
		assert( ndim == 3 );
		for(;;){
			fgets(stmp1,buff_size,fp);
			if( strncmp(stmp1 ,"end coordinates", 15)==0 )break;
			int inode;
			double co_x, co_y, co_z;
			sscanf(stmp1,"%d%lf%lf%lf",&inode, &co_x, &co_y, &co_z);
			assert( (int)counter+1 == inode );
			counter++;
			unsigned int i0 = tmp_buffer.size();
			tmp_buffer.resize( i0+3 );
			tmp_buffer[i0  ] = co_x;
			tmp_buffer[i0+1] = co_y;
			tmp_buffer[i0+2] = co_z;
		}
		const unsigned int nnode = counter;
		assert( nnode != 0 );
		aVec.resize(nnode);
		for(unsigned int inode=0;inode<nnode;inode++){
			aVec[inode].x = tmp_buffer[inode*3  ];
			aVec[inode].y = tmp_buffer[inode*3+1];
			aVec[inode].z = tmp_buffer[inode*3+2];
		}
	}	// 節点の読み込み終了

	fgets(stmp1,buff_size,fp);	// 改行

	// 要素の読み込み
	fgets(stmp1,buff_size,fp);	// Elements
	std::vector<int> tmp_buffer;
	tmp_buffer.reserve(16384);
	if( strcmp(str_elem_type,"Hexahedra")==0 ){
		unsigned int counter = 0;
		for(;;){
			fgets(stmp1,buff_size,fp);
			if( strncmp(stmp1 ,"end elements", 12)==0 ) break;
			int ielem;
			int v1, v2, v3, v4, v5, v6, v7, v8;
			sscanf(stmp1,"%d%d%d%d%d%d%d%d%d",&ielem, &v1,&v2,&v3,&v4,&v5,&v6,&v7,&v8);
			assert( (int)counter+1 == ielem );
			assert( v1>0 && v2>0 && v3>0 && v4>0 && v5>0 && v6>0 && v7>0 && v8>0);
			counter++;
			unsigned int i0 = tmp_buffer.size();
			tmp_buffer.resize( i0+8 );
			tmp_buffer[i0  ] = v1-1;
			tmp_buffer[i0+1] = v2-1;
			tmp_buffer[i0+2] = v3-1;
			tmp_buffer[i0+3] = v4-1;
			tmp_buffer[i0+4] = v5-1;
			tmp_buffer[i0+5] = v6-1;
			tmp_buffer[i0+6] = v7-1;
			tmp_buffer[i0+7] = v8-1;
		}
		const unsigned int nhex = counter;
		m_aHexAry.resize(1);
//		m_aHexAry[0].m_CadLoopID = 0;
		m_aHexAry[0].id = 1;
		m_aHexAry[0].m_aHex.resize(nhex);
		for(unsigned int ihex=0;ihex<nhex;ihex++){
			m_aHexAry[0].m_aHex[ihex].v[0] = tmp_buffer[ihex*8  ];
			m_aHexAry[0].m_aHex[ihex].v[1] = tmp_buffer[ihex*8+1];
			m_aHexAry[0].m_aHex[ihex].v[2] = tmp_buffer[ihex*8+2];
			m_aHexAry[0].m_aHex[ihex].v[3] = tmp_buffer[ihex*8+3];
			m_aHexAry[0].m_aHex[ihex].v[4] = tmp_buffer[ihex*8+4];
			m_aHexAry[0].m_aHex[ihex].v[5] = tmp_buffer[ihex*8+5];
			m_aHexAry[0].m_aHex[ihex].v[6] = tmp_buffer[ihex*8+6];
			m_aHexAry[0].m_aHex[ihex].v[7] = tmp_buffer[ihex*8+7];
		}
		Msh::MakeHexSurHex(m_aHexAry[0].m_aHex);
		std::vector<SQuad3D> aQuadBound;
//		Msh::MakeOuterBoundHex(m_aHexAry[0].m_aHex,aQuadBound);
//		Msh::MakeQuadSurQuad(aQuadBound);
	}
	else if( strcmp(str_elem_type,"Tetrahedra")==0 ){
		unsigned int counter = 0;
		for(;;){
			fgets(stmp1,buff_size,fp);
			if( strncmp(stmp1 ,"end elements", 12)==0 ) break;
			int ielem;
			int v1, v2, v3, v4;
			sscanf(stmp1,"%d %d%d%d%d",&ielem, &v1,&v2,&v3,&v4);
			assert( (int)counter+1 == ielem );
			assert( v1>0 && v2>0 && v3>0 && v4>0);
			counter++;
			unsigned int i0 = tmp_buffer.size();
			tmp_buffer.resize( i0+4 );
			tmp_buffer[i0  ] = v1-1;
			tmp_buffer[i0+1] = v2-1;
			tmp_buffer[i0+2] = v3-1;
			tmp_buffer[i0+3] = v4-1;
		}
		const unsigned int ntet = counter;
		m_aTetAry.resize(1);
		m_aTetAry[0].id = 1;
		m_aTetAry[0].m_aTet.resize(ntet);
		for(unsigned int itet=0;itet<ntet;itet++){
			m_aTetAry[0].m_aTet[itet].v[0] = tmp_buffer[itet*4  ];
			m_aTetAry[0].m_aTet[itet].v[1] = tmp_buffer[itet*4+1];
			m_aTetAry[0].m_aTet[itet].v[2] = tmp_buffer[itet*4+2];
			m_aTetAry[0].m_aTet[itet].v[3] = tmp_buffer[itet*4+3];
		}
		m_include_relation.resize(2);
		Msh::MakeTetSurTet(m_aTetAry[0].m_aTet);
		std::vector<STri3D> aTriBound;
		Msh::MakeOuterBoundTet(m_aTetAry[0].m_aTet,aTriBound);
		Msh::MakeTriSurTri(aTriBound);
		std::vector<int> aColorTri;
		unsigned int nColorTri;
		Msh::MakeColorCoding_Tri3D(aTriBound,aVec,aColorTri,nColorTri,0.2);
		m_aTriAry.resize(nColorTri);
		for(unsigned int icolor=0;icolor<nColorTri;icolor++){
			unsigned int ntri_color = 0;
			for(unsigned int itri=0;itri<aTriBound.size();itri++){
				if( aColorTri[itri] == (int)icolor ) ntri_color++;
			}
			m_aTriAry[icolor].m_aTri.resize(ntri_color);
			unsigned int icnt0 = 0;
			for(unsigned int itri=0;itri<aTriBound.size();itri++){
				if( aColorTri[itri] != (int)icolor ) continue;
				m_aTriAry[icolor].m_aTri[icnt0].v[0] = aTriBound[itri].v[0];
				m_aTriAry[icolor].m_aTri[icnt0].v[1] = aTriBound[itri].v[1];
				m_aTriAry[icolor].m_aTri[icnt0].v[2] = aTriBound[itri].v[2];
				icnt0++;
			}
			m_aTriAry[icolor].id = 2+icolor;
			m_include_relation[1].push_back(2+icolor);
		}
		/*
		std::vector<SBar> aBarBound;
		MakeBar_fromColorCodingTri( aTriBound, aColorTri, nColorTri, aBarBound);
		std::vector<int> aColorBar;
		unsigned int nColorBar;		
		MakeColorCodingBar( aBarBound, aColorBar, nColorBar);
		std::cout << nColorBar << std::endl;
		m_include_relation.resize(2+nColorTri);
		m_aBarAry.resize( nColorBar );
		for(unsigned int icolor=0;icolor<nColorBar;icolor++){
			unsigned int nbar_color = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] == icolor ) nbar_color++;
			}
			m_aBarAry[icolor].m_aBar.resize(nbar_color);
			unsigned int icnt0 = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] != icolor ) continue;
				m_aBarAry[icolor].m_aBar[icnt0].v[0] = aBarBound[ibar].v[0];
				m_aBarAry[icolor].m_aBar[icnt0].v[1] = aBarBound[ibar].v[1];
				icnt0++;
			}
			int icolor0=-1, icolor1=-1;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] != icolor ) continue;
				unsigned int ic0 = aColorTri[ aBarBound[ibar].s2[0] ];
				unsigned int ic1 = aColorTri[ aBarBound[ibar].s2[1] ];
//				std::cout << ic0 << " " << ic1 << std::endl;
				if( icolor0 == -1 ){ icolor0 = ic0; }
				else if( icolor1 == -1 ){
					if( icolor0 == ic0 ){ icolor1 = ic1; }
					else if( icolor0 == ic1 ){ icolor1 = ic0; }
				}
			}
//			std::cout << icolor0 << " " << icolor1 << std::endl;
			m_aBarAry[icolor].id = 2+nColorTri+icolor;
			m_include_relation[1].push_back(2+nColorTri+icolor);
			m_include_relation[2+icolor0].push_back(2+nColorTri+icolor);
			m_include_relation[2+icolor1].push_back(2+nColorTri+icolor);
		}
		*/
	}
	else if( strcmp(str_elem_type,"Triangle")==0 ){
		unsigned int counter = 0;
		for(;;){
			fgets(stmp1,buff_size,fp);
			if( strncmp(stmp1 ,"end elements", 12)==0 ) break;
			int ielem;
			int v1, v2, v3;
			sscanf(stmp1,"%d %d%d%d",&ielem, &v1,&v2,&v3);
			assert( (int)counter+1 == ielem );
			assert( v1>0 && v2>0 && v3>0 );
			counter++;
			unsigned int i0 = tmp_buffer.size();
			tmp_buffer.resize( i0+3 );
			tmp_buffer[i0  ] = v1-1;
			tmp_buffer[i0+1] = v2-1;
			tmp_buffer[i0+2] = v3-1;
		}
		const unsigned int ntri = counter;
		m_aTriAry.resize(1);
//		m_aTriAry[0].m_CadLoopID = 0;
		m_aTriAry[0].id = 1;
		m_aTriAry[0].m_aTri.resize(ntri);
		for(unsigned int itri=0;itri<ntri;itri++){
			m_aTriAry[0].m_aTri[itri].v[0] = tmp_buffer[itri*3  ];
			m_aTriAry[0].m_aTri[itri].v[1] = tmp_buffer[itri*3+1];
			m_aTriAry[0].m_aTri[itri].v[2] = tmp_buffer[itri*3+2];
		}
		/*
		m_include_relation.resize(2);
		Msh::MakeTriSurTri(m_aTriAry[0].m_aTri);
		std::vector<SBar> aBarBound;
		Msh::MakeOuterBoundTri(m_aTriAry[0].m_aTri,aBarBound);
		MakeInnerRelationBar( aBarBound );
		std::vector<int> aColorBar;
		unsigned int nColor;
		Msh::MakeColorCodingBar(aBarBound,aVec,aColorBar,nColor);
		m_aBarAry.resize(nColor);
		for(unsigned int icolor=0;icolor<nColor;icolor++){
			unsigned int nbar_color = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] == icolor ) nbar_color++;
			}
			m_aBarAry[icolor].m_aBar.resize(nbar_color);
			unsigned int icnt0 = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] != icolor ) continue;
				m_aBarAry[icolor].m_aBar[icnt0].v[0] = aBarBound[ibar].v[0];
				m_aBarAry[icolor].m_aBar[icnt0].v[1] = aBarBound[ibar].v[1];
				icnt0++;
			}
			m_aBarAry[icolor].id = 2+icolor;
			m_include_relation[1].push_back(2+icolor);
		}*/
	}
	fclose(fp);
	this->MakeElemLocationType();
	return true;
}







bool Msh::CMesh3D::ReadFromFile_TetgenMsh(const std::string& file_name){
	FILE *fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];
	unsigned int itmp0, itmp1;

	////////////////////////////////
	// ここからは内容を変更する

	this->aVec.clear();
	this->m_aBarAry.clear();
	this->m_aTriAry.clear();
	this->m_aQuadAry.clear();
	this->m_aTetAry.clear();
	this->m_aHexAry.clear();
	this->m_aVertex.clear();
	this->m_ElemLoc.clear();
	this->m_ElemType.clear();

    ////////////////

	{	// 節点の読み込み
        std::string fname_node = file_name + ".node";
	    if( (fp = ::fopen(fname_node.c_str(),"r"))== NULL ){
		    std::cout << "Cannot Open Tetgen node File : " << fname_node << std::endl;
		    return false;
	    }
	    int ndim, nnode;
	    fgets(stmp1,buff_size,fp);
	    sscanf(stmp1,"%d %d %d %d", &nnode,&ndim,&itmp0,&itmp1);
        std::cout << nnode << " " << ndim << std::endl;
		aVec.resize(nnode);
		assert( ndim == 3 );
		for(unsigned int inode=0;inode<(unsigned int)nnode;inode++){
			double co_x, co_y, co_z;
			fgets(stmp1,buff_size,fp);
			sscanf(stmp1,"%d%lf%lf%lf",&itmp0, &co_x, &co_y, &co_z);
			assert( itmp0 == inode );
			aVec[inode].x = co_x;
			aVec[inode].y = co_y;
			aVec[inode].z = co_z;
		}
        fclose(fp);
	}	// 節点の読み込み終了

	
    {   // 要素の読み込み  
        std::string fname_ele = file_name + ".ele";
	    if( (fp = ::fopen(fname_ele.c_str(),"r"))== NULL ){
		    std::cout << "Cannot Open Tetgen ele File : " << fname_ele << std::endl;
		    return false;
	    }
	    int nnoel, nelem;
	    fgets(stmp1,buff_size,fp);
	    sscanf(stmp1,"%d %d %d", &nelem,&nnoel,&itmp0);
		const unsigned int ntet = nelem;
		m_aTetAry.resize(1);
		m_aTetAry[0].id = 1;
		m_aTetAry[0].m_aTet.resize(ntet);
		for(unsigned int ielem=0;ielem<(unsigned int)nelem;ielem++){
			fgets(stmp1,buff_size,fp);
			int v1, v2, v3, v4;
			sscanf(stmp1,"%d %d%d%d%d",&itmp0,&v1,&v2,&v3,&v4);
			assert( itmp0 == ielem );
			assert( v1>=0 && v2>=0 && v3>=0 && v4>=0);
			m_aTetAry[0].m_aTet[ielem].v[0] = v1;
			m_aTetAry[0].m_aTet[ielem].v[1] = v2;
			m_aTetAry[0].m_aTet[ielem].v[2] = v3;
			m_aTetAry[0].m_aTet[ielem].v[3] = v4;
		}
		m_include_relation.resize(2);
		Msh::MakeTetSurTet(m_aTetAry[0].m_aTet);
		std::vector<STri3D> aTriBound;
		Msh::MakeOuterBoundTet(m_aTetAry[0].m_aTet,aTriBound);
		Msh::MakeTriSurTri(aTriBound);
		std::vector<int> aColorTri;
		unsigned int nColorTri;
		Msh::MakeColorCoding_Tri3D(aTriBound,aVec,aColorTri,nColorTri,0.8);
		m_aTriAry.resize(nColorTri);
		for(unsigned int icolor=0;icolor<nColorTri;icolor++){
			unsigned int ntri_color = 0;
			for(unsigned int itri=0;itri<aTriBound.size();itri++){
				if( aColorTri[itri] == (int)icolor ) ntri_color++;
			}
			m_aTriAry[icolor].m_aTri.resize(ntri_color);
			unsigned int icnt0 = 0;
			for(unsigned int itri=0;itri<aTriBound.size();itri++){
				if( aColorTri[itri] != (int)icolor ) continue;
				m_aTriAry[icolor].m_aTri[icnt0].v[0] = aTriBound[itri].v[0];
				m_aTriAry[icolor].m_aTri[icnt0].v[1] = aTriBound[itri].v[1];
				m_aTriAry[icolor].m_aTri[icnt0].v[2] = aTriBound[itri].v[2];
				icnt0++;
			}
			m_aTriAry[icolor].id = 2+icolor;
			m_include_relation[1].push_back(2+icolor);
		}
		/*
		std::vector<SBar> aBarBound;
		MakeBar_fromColorCodingTri( aTriBound, aColorTri, nColorTri, aBarBound);
		std::vector<int> aColorBar;
		unsigned int nColorBar;		
		MakeColorCodingBar( aBarBound, aColorBar, nColorBar);
		std::cout << nColorBar << std::endl;
		m_include_relation.resize(2+nColorTri);
		m_aBarAry.resize( nColorBar );
		for(unsigned int icolor=0;icolor<nColorBar;icolor++){
			unsigned int nbar_color = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] == icolor ) nbar_color++;
			}
			m_aBarAry[icolor].m_aBar.resize(nbar_color);
			unsigned int icnt0 = 0;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] != icolor ) continue;
				m_aBarAry[icolor].m_aBar[icnt0].v[0] = aBarBound[ibar].v[0];
				m_aBarAry[icolor].m_aBar[icnt0].v[1] = aBarBound[ibar].v[1];
				icnt0++;
			}
			int icolor0=-1, icolor1=-1;
			for(unsigned int ibar=0;ibar<aBarBound.size();ibar++){
				if( aColorBar[ibar] != icolor ) continue;
				unsigned int ic0 = aColorTri[ aBarBound[ibar].s2[0] ];
				unsigned int ic1 = aColorTri[ aBarBound[ibar].s2[1] ];
//				std::cout << ic0 << " " << ic1 << std::endl;
				if( icolor0 == -1 ){ icolor0 = ic0; }
				else if( icolor1 == -1 ){
					if( icolor0 == ic0 ){ icolor1 = ic1; }
					else if( icolor0 == ic1 ){ icolor1 = ic0; }
				}
			}
//			std::cout << icolor0 << " " << icolor1 << std::endl;
			m_aBarAry[icolor].id = 2+nColorTri+icolor;
			m_include_relation[1].push_back(2+nColorTri+icolor);
			m_include_relation[2+icolor0].push_back(2+nColorTri+icolor);
			m_include_relation[2+icolor1].push_back(2+nColorTri+icolor);
		}*/
	}
	fclose(fp);
	this->MakeElemLocationType();
	return true;
}


