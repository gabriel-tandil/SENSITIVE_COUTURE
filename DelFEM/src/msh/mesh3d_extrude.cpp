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

#include "delfem/msh/meshkernel3d.h"
#include "delfem/msh/meshkernel2d.h"
#include "delfem/mesh3d.h"
#include "delfem/mesher2d.h"

using namespace Com;

bool Msh::CMesh3D_Extrude::UpdateMeshCoord(const CMesher2D& msh_2d)
{
	const unsigned int nvec3 = this->aVec.size();
	const unsigned int nvec2 = msh_2d.GetVectorAry().size();
	if( nvec3 % nvec2 != 0 ) return false;
	const std::vector<CVector2D>& aVec2D = msh_2d.GetVectorAry();
	const unsigned int nlev = nvec3 / nvec2;
	for(unsigned int ilev=0;ilev<nlev;ilev++){
		for(unsigned int ivec=0;ivec<nvec2;ivec++){
			this->aVec[ilev*nvec2+ivec].x = aVec2D[ivec].x;
			this->aVec[ilev*nvec2+ivec].y = aVec2D[ivec].y;
		}
	}
	return true;
}


bool Msh::CMesh3D_Extrude::Extrude(const CMesher2D& msh_2d, const double height, const double elen)
{
	if( elen <= 0 ) return false;
	if( height <= 0 ) return false;

	////////////////////////////////
	const unsigned int naTri = msh_2d.GetTriArySet().size();
	const unsigned int naQuad = msh_2d.GetQuadArySet().size();

	if(      naTri == 1 && naQuad == 0 ){ return this->Extrude_Tri( msh_2d,height,elen); }
	else if( naTri == 0 && naQuad == 1 ){ return this->Extrude_Quad(msh_2d,height,elen); }
	assert(0);
	return false;
}

bool Msh::CMesh3D_Extrude::Extrude_Quad(const CMesher2D& msh_2d, const double height, const double elen)
{
	const unsigned int nNo2d = msh_2d.GetVectorAry().size();
	unsigned int ndiv = (unsigned int)(height / elen);

	const unsigned int nVer2d    = msh_2d.GetVertexAry().size();
	const unsigned int nQuadAry2d = msh_2d.GetQuadArySet().size();
	const unsigned int nBarAry2d = msh_2d.GetBarArySet().size();

	std::vector<unsigned int> id_ary_3d;
	{
		const unsigned int elem_size = 
			 nVer2d*2+				// 点 : 下＋上
			 nBarAry2d+nVer2d*2+	// 辺 : 垂直辺、水平辺(下＋上)
			 nQuadAry2d*2+nBarAry2d+	// 面 : 下＋上＋垂直面
			 nQuadAry2d;				// 体
		id_ary_3d.resize( elem_size, 0 );
	}

	{	// 頂点を作る
		const std::vector<SVertex>& aVer2d = msh_2d.GetVertexAry();
		m_aVertex.resize( nVer2d*2 );
		for(unsigned int iver=0;iver<nVer2d;iver++){    // 底面を作る
			m_aVertex[iver].id_cad                = aVer2d[iver].id_v_cad;
			m_aVertex[iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aVertex[iver].inum_extrude = 1;
			const unsigned int id0 = GetFreeObjID(); 
			m_aVertex[iver].id = id0; 
			m_aVertex[iver].v = aVer2d[iver].v;
			id_ary_3d[iver] = id0;
		}
		for(unsigned int iver=0;iver<nVer2d;iver++){    // 上面を作る
			m_aVertex[nVer2d+iver].id_cad                = aVer2d[iver].id_v_cad;
			m_aVertex[nVer2d+iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aVertex[nVer2d+iver].inum_extrude = 3;
			const unsigned int id0 = GetFreeObjID();
			m_aVertex[nVer2d+iver].id = id0;
			m_aVertex[nVer2d+iver].v = aVer2d[iver].v+nNo2d*ndiv;
			id_ary_3d[nVer2d+iver] = id0;
		}
	}

	{ // 辺を作る
		const std::vector<SVertex>& aVer2d = msh_2d.GetVertexAry();
		m_aBarAry.resize( nBarAry2d*2+nVer2d );
		// 下の辺を作る
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			m_aBarAry[ibar_ary].id_cad                = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aBarAry[ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aBarAry[ibar_ary].inum_extrude = 1;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[ibar_ary].id = id0;
            {
			    const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			    std::vector<SBar>& aBar_d = m_aBarAry[ibar_ary*2  ].m_aBar;
			    aBar_d.resize(aBar2d.size());
			    for(unsigned int ibar=0;ibar<aBar2d.size();ibar++){
				    aBar_d[ibar].v[0] = aBar2d[ibar].v[0];
				    aBar_d[ibar].v[1] = aBar2d[ibar].v[1];
			    }
            }
			id_ary_3d[nVer2d*2+ibar_ary] = id0;
		}
		// 上の辺を作る
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			m_aBarAry[nBarAry2d+ibar_ary].id_cad                = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aBarAry[nBarAry2d+ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aBarAry[nBarAry2d+ibar_ary].inum_extrude = 3;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[nBarAry2d+ibar_ary].id = id0;
            {
			    const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			    std::vector<SBar>& aBar_u = m_aBarAry[ibar_ary*2+1].m_aBar;
			    aBar_u.resize(aBar2d.size());
			    for(unsigned int ibar=0;ibar<aBar2d.size();ibar++){
				    aBar_u[ibar].v[0] = aBar2d[ibar].v[0]+nNo2d*ndiv;
				    aBar_u[ibar].v[1] = aBar2d[ibar].v[1]+nNo2d*ndiv;
			    }
            }
			id_ary_3d[nVer2d*2+nBarAry2d+ibar_ary] = id0;
		}
		// 垂直辺を作る
		for(unsigned int iver=0;iver<nVer2d;iver++){
    		m_aBarAry[nBarAry2d*2+iver].id_cad                = aVer2d[iver].id_v_cad;
			m_aBarAry[nBarAry2d*2+iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aBarAry[nBarAry2d*2+iver].inum_extrude = 2;
			const unsigned int iv2d = aVer2d[iver].v;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[nBarAry2d*2+iver].id = id0;
            {
			    std::vector<SBar>& aBar = m_aBarAry[nBarAry2d*2+iver].m_aBar;
			    aBar.resize(ndiv);
			    for(unsigned int idiv=0;idiv<ndiv;idiv++){
				    aBar[idiv].v[0] = iv2d+nNo2d*idiv;
				    aBar[idiv].v[1] = iv2d+nNo2d*(idiv+1);
			    }
            }
			id_ary_3d[nVer2d*2+nBarAry2d*2+iver] = id0;
		}
	}

	{	// ３次元３角形を作る
		m_aQuadAry.resize(nQuadAry2d*2+nBarAry2d);
		{	// 上下の３角形を作る
			const std::vector<CQuadAry2D>& aQuadAry2D = msh_2d.GetQuadArySet();
			assert( aQuadAry2D.size() == 1 );
			const std::vector<SQuad2D>& aQuad2d = aQuadAry2D[0].m_aQuad;
			const unsigned int nNo2d = msh_2d.GetVectorAry().size();
			{	// 下の三角形
				m_aQuadAry[0].id_cad                = aQuadAry2D[0].id_l_cad;
				m_aQuadAry[0].id_msh_before_extrude = aQuadAry2D[0].id;
				m_aQuadAry[0].inum_extrude = 1;
				const unsigned int id0 = GetFreeObjID();
				m_aQuadAry[0].id = id0;
				std::vector<SQuad3D>& aQuad = m_aQuadAry[0].m_aQuad;
				aQuad.resize( aQuad2d.size() );
				for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
					aQuad[iquad].v[0] = aQuad2d[iquad].v[3];
					aQuad[iquad].v[1] = aQuad2d[iquad].v[2];
					aQuad[iquad].v[2] = aQuad2d[iquad].v[1];
					aQuad[iquad].v[3] = aQuad2d[iquad].v[0];
				}
				id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d  ] = id0;
			}
			{	// 上の三角形配列
				m_aQuadAry[1].id_cad                = aQuadAry2D[0].id_l_cad;
				m_aQuadAry[1].id_msh_before_extrude = aQuadAry2D[0].id;
				m_aQuadAry[1].inum_extrude = 3;
				const unsigned int id0 = GetFreeObjID();
				m_aQuadAry[1].id = id0;
				std::vector<SQuad3D>& aQuad = m_aQuadAry[1].m_aQuad;
				aQuad.resize( aQuad2d.size() );
				for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
					aQuad[iquad].v[0] = aQuad2d[iquad].v[0]+nNo2d*ndiv;
					aQuad[iquad].v[1] = aQuad2d[iquad].v[1]+nNo2d*ndiv;
					aQuad[iquad].v[2] = aQuad2d[iquad].v[2]+nNo2d*ndiv;
					aQuad[iquad].v[3] = aQuad2d[iquad].v[3]+nNo2d*ndiv;
				}
				id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d] = id0;
			}
		}
		// 周囲の4角形を作る
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			const unsigned int nBar2d = aBar2d.size();
			std::vector<SQuad3D>& aQuad = m_aQuadAry[nQuadAry2d*2+ibar_ary].m_aQuad;
			aQuad.resize( nBar2d*ndiv );
			m_aQuadAry[nQuadAry2d*2+ibar_ary].id_cad                = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aQuadAry[nQuadAry2d*2+ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aQuadAry[nQuadAry2d*2+ibar_ary].inum_extrude = 2;
			const unsigned int id0 = GetFreeObjID();
			m_aQuadAry[nQuadAry2d*2+ibar_ary].id = id0;
			for(unsigned int idiv=0;idiv<ndiv;idiv++){
			for(unsigned int ibar=0;ibar<aBar2d.size();ibar++){
				aQuad[idiv*nBar2d+ibar].v[0] = aBar2d[ibar].v[0]+nNo2d*idiv;
				aQuad[idiv*nBar2d+ibar].v[1] = aBar2d[ibar].v[1]+nNo2d*idiv;
				aQuad[idiv*nBar2d+ibar].v[2] = aBar2d[ibar].v[1]+nNo2d*(idiv+1);
				aQuad[idiv*nBar2d+ibar].v[3] = aBar2d[ibar].v[0]+nNo2d*(idiv+1);
			}
			}
			id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d*2+ibar_ary] = id0;
		}
	}

	{	// 六面体を作るルーティン
		assert( nQuadAry2d == 1 );
		const std::vector<CQuadAry2D>& aQuadAry2D = msh_2d.GetQuadArySet();
		m_aHexAry.resize(nQuadAry2d);
		const std::vector<SQuad2D>& aQuad     = aQuadAry2D[0].m_aQuad;
		const unsigned int nQuad2d = aQuad.size();
		m_aHexAry[0].id_cad                = aQuadAry2D[0].id_l_cad;
		m_aHexAry[0].id_msh_before_extrude = aQuadAry2D[0].id;
		m_aHexAry[0].inum_extrude = 2;
		const unsigned int id0 = GetFreeObjID();
		m_aHexAry[0].id = id0;
		std::vector<SHex>& aHex = m_aHexAry[0].m_aHex;
		aHex.resize( aQuad.size()*ndiv );
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
		for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
			const unsigned int ihex0 = nQuad2d*idiv+iquad;
			aHex[ihex0].v[0] = nNo2d* idiv   +aQuad[iquad].v[0];
			aHex[ihex0].v[1] = nNo2d* idiv   +aQuad[iquad].v[1];
			aHex[ihex0].v[2] = nNo2d* idiv   +aQuad[iquad].v[2];
			aHex[ihex0].v[3] = nNo2d* idiv   +aQuad[iquad].v[3];
			aHex[ihex0].v[4] = nNo2d*(idiv+1)+aQuad[iquad].v[0];
			aHex[ihex0].v[5] = nNo2d*(idiv+1)+aQuad[iquad].v[1];
			aHex[ihex0].v[6] = nNo2d*(idiv+1)+aQuad[iquad].v[2];
			aHex[ihex0].v[7] = nNo2d*(idiv+1)+aQuad[iquad].v[3];
		}
		}
		MakeHexSurHex(aHex);
//		assert( CheckHex(aHex) );
		id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d*2+nBarAry2d  ] = id0;
	}

	{	// 節点を作るルーティン
		double hlen = height / ndiv;
		const unsigned int nNo2d = msh_2d.GetVectorAry().size();
		this->aVec.resize((ndiv+1)*nNo2d);
		const std::vector<CVector2D>& aVec2D = msh_2d.GetVectorAry();
		for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
		for(unsigned int ino=0;ino<nNo2d;ino++){
			this->aVec[idiv*nNo2d+ino].x = aVec2D[ino].x;
			this->aVec[idiv*nNo2d+ino].y = aVec2D[ino].y;
			this->aVec[idiv*nNo2d+ino].z = hlen*idiv;
		}
		}
	}

	{	// 包含関係を作る
		m_include_relation.resize( this->FindMaxID()+1 );
		for(unsigned int iver2d=0;iver2d<nBarAry2d;iver2d++){
			const unsigned int id_bar_v_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+iver2d];	// 垂直辺
			m_include_relation[id_bar_v_3d].push_back( id_ary_3d[iver2d] );	// 下頂点
			m_include_relation[id_bar_v_3d].push_back( id_ary_3d[nVer2d+iver2d] ); // 上頂点
		}
		for(unsigned int ibar2d=0;ibar2d<nBarAry2d;ibar2d++){
			const unsigned int id_bar_d_3d = id_ary_3d[nVer2d*2+  ibar2d]; // 下辺
			const unsigned int id_bar_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d+ibar2d]; // 上辺
			const unsigned int id_quad_3d  = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d*2+ibar2d]; // 垂直面

			const unsigned int id_bar_2d = msh_2d.GetBarArySet()[ibar2d].id;
			std::vector<unsigned int> id_ary_inc = msh_2d.GetIncludeElemIDAry(id_bar_2d);
			for(unsigned int iinc=0;iinc<id_ary_inc.size();iinc++){
				const unsigned int id_inc = id_ary_inc[iinc];
				unsigned int iver2d;
				{
					unsigned int id_cad;
					unsigned int nelem;
					MSH_TYPE msh_type;
					msh_2d.GetMshInfo(id_inc,nelem,msh_type,iver2d, id_cad);
					assert( msh_type == VERTEX );
					assert( iver2d < nVer2d );
				}
				m_include_relation[id_bar_d_3d].push_back( id_ary_3d[iver2d] );	// 下辺
				m_include_relation[id_bar_u_3d].push_back( id_ary_3d[nVer2d+iver2d] );	// 上辺
				m_include_relation[id_quad_3d].push_back( id_ary_3d[iver2d] );		// 垂直面
				m_include_relation[id_quad_3d].push_back( id_ary_3d[nVer2d+iver2d] );	// 垂直面
				m_include_relation[id_quad_3d].push_back( id_ary_3d[nVer2d*2+nBarAry2d*2+iver2d] );	// 垂直面
			}
			m_include_relation[id_quad_3d].push_back( id_ary_3d[nVer2d*2+ibar2d] );	// 垂直面
			m_include_relation[id_quad_3d].push_back( id_ary_3d[nVer2d*2+nBarAry2d+ibar2d] );	// 垂直面
		}
		for(unsigned int iquad2d=0;iquad2d<nQuadAry2d;iquad2d++){
			const unsigned int id_quad_d_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  iquad2d];
			const unsigned int id_quad_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  iquad2d+nQuadAry2d];
			const unsigned int id_hex_3d    = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d*2       + nBarAry2d + iquad2d];

			const unsigned int id_quad_2d = msh_2d.GetQuadArySet()[iquad2d].id;
			std::vector<unsigned int> id_ary_inc = msh_2d.GetIncludeElemIDAry(id_quad_2d);
			for(unsigned int iinc=0;iinc<id_ary_inc.size();iinc++){
				const unsigned int id_inc = id_ary_inc[iinc];
				unsigned int iloc2d, id_cad;
				unsigned int nelem;
				MSH_TYPE msh_type;
				msh_2d.GetMshInfo(id_inc,nelem,msh_type,iloc2d, id_cad);
				if( msh_type == VERTEX ){
					assert( iloc2d < nVer2d );
					const unsigned int id_ver_d_3d = id_ary_3d[iloc2d];
					const unsigned int id_ver_u_3d = id_ary_3d[nVer2d+iloc2d];
					const unsigned int id_bar_v_3d   = id_ary_3d[nVer2d*2+  nBarAry2d*2+iloc2d];
					m_include_relation[id_quad_d_3d].push_back( id_ver_d_3d );
					m_include_relation[id_quad_u_3d].push_back( id_ver_u_3d );
					m_include_relation[id_hex_3d   ].push_back( id_ver_d_3d );
					m_include_relation[id_hex_3d   ].push_back( id_ver_u_3d );
					m_include_relation[id_hex_3d   ].push_back( id_bar_v_3d );
				}
				else if( msh_type == BAR ){
					assert( iloc2d < nBarAry2d );
					const unsigned int id_bar_d_3d = id_ary_3d[nVer2d*2+  iloc2d];
					const unsigned int id_bar_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d+iloc2d];
					const unsigned int id_quad_v_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d*2+iloc2d];
					m_include_relation[id_quad_d_3d].push_back( id_bar_d_3d );
					m_include_relation[id_quad_u_3d].push_back( id_bar_u_3d );
					m_include_relation[id_hex_3d   ].push_back( id_bar_d_3d );
					m_include_relation[id_hex_3d   ].push_back( id_bar_u_3d );
					m_include_relation[id_hex_3d   ].push_back( id_quad_v_3d );
				}
				else{ assert(0); }
			}
			{
				unsigned int id_quad_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  iquad2d];
				unsigned int id_quad_d_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nQuadAry2d+iquad2d];
				m_include_relation[id_hex_3d].push_back( id_quad_u_3d );
				m_include_relation[id_hex_3d].push_back( id_quad_d_3d );
			}
		}
	}

	this->MakeElemLocationType();
/*
	for(unsigned int id_msh=0;id_msh<m_include_relation.size();id_msh++){
		std::cout << id_msh << " --> ";
		for(unsigned int iid_msh_inc=0;iid_msh_inc<m_include_relation[id_msh].size();iid_msh_inc++){
			unsigned int id_msh_inc = m_include_relation[id_msh][iid_msh_inc];
			std::cout << id_msh_inc << " ";
		}
		std::cout << std::endl;
	}
	getchar();
*/
	return true;
}

bool Msh::CMesh3D_Extrude::MakeDirectionFlagForExtrusion(const CMesher2D& msh_2d,
								   std::vector<std::vector<int> >& tri_flg_ary, 
								   std::vector< std::vector<int> >& bar_flg_ary )
{
	{	// ３角形のフラグを初期化する
		const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
		assert( aTriAry2D.size() == 1 );
		const std::vector<STri2D>& aTri = aTriAry2D[0].m_aTri;
		tri_flg_ary.resize(1);
		tri_flg_ary[0].resize( aTri.size()*3, 0 );
	}

	{	// 辺のフラグを初期化する
		const unsigned int nBarAry = msh_2d.GetBarArySet().size();
		bar_flg_ary.resize( nBarAry );
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry;ibar_ary++){
			const std::vector<SBar>& aBar = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			bar_flg_ary[ibar_ary].resize(aBar.size(),0);
			int idmshl = msh_2d.GetBarArySet()[ibar_ary].id_lr[0];
			int idmshr = msh_2d.GetBarArySet()[ibar_ary].id_lr[1];
			assert( idmshl==0 || idmshl==-1 || msh_2d.IsID(idmshl) );
			assert( idmshr==0 || idmshr==-1 || msh_2d.IsID(idmshr) );
		}
	}

	{
		// 初期設定として隣り合う要素の同士に違うフラグを与える
		const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
		assert( aTriAry2D.size() == 1 );
		const std::vector<STri2D>& aTri = aTriAry2D[0].m_aTri;
		std::vector<int>& tri_flg = tri_flg_ary[0];
		for(unsigned int itri=0;itri<aTri.size();itri++){
		for(unsigned int inotri=0;inotri<3;inotri++){
			int iflg0 = tri_flg[itri*3+inotri]; assert( iflg0==0 || iflg0==1 );
			if( aTri[itri].g2[inotri] == -2 ){
				const unsigned int itri_s = aTri[itri].s2[inotri];
				const unsigned int ied_s = relTriTri[ aTri[itri].r2[inotri] ][inotri];
				assert( aTri[itri_s].s2[ied_s] == itri );
				tri_flg[itri_s*3+ied_s] = 1-iflg0;
			}
		}
		}
		std::vector<int> elem_msk;
		elem_msk.resize( aTri.size(), 0 );
		unsigned int icoun_pre = 0;
		for(unsigned int jloop=0;jloop<10;jloop++){
		for(unsigned int iloop=0;iloop<3;iloop++){
			unsigned int icoun = 0;
			for(unsigned int itri=0;itri<aTri.size();itri++){
				if( elem_msk[itri] == 1 ) continue;
				int iflg0 = tri_flg[itri*3  ];
				int iflg1 = tri_flg[itri*3+1];
				int iflg2 = tri_flg[itri*3+2];
				int iflg_a = -1;
				if( iflg0*iflg1*iflg2 != 0 ){ iflg_a = 1; } 
				if( (iflg0-1)*(iflg1-1)*(iflg2-1) != 0 ){ iflg_a = 0; }
				if( iflg_a == -1 ) continue;
				icoun++;
				// 始めは厳しい条件で判定する。つまり、隣と入れ替えるときに隣が駄目になってはならない。
				for(unsigned int inotri=0;inotri<3;inotri++){
					if( aTri[itri].g2[inotri] == -2 ){
						const unsigned int itri_s = aTri[itri].s2[inotri];
						const unsigned int ied_s = relTriTri[ aTri[itri].r2[inotri] ][inotri];
						assert( aTri[itri_s].s2[ied_s] == itri );
						assert( tri_flg[itri_s*3+ied_s] == 1-iflg_a );
						if( elem_msk[itri_s] == 1 ) continue;
						unsigned int ied_s0 = noelTriEdge[ied_s][0];
						unsigned int ied_s1 = noelTriEdge[ied_s][1];
						if( tri_flg[itri_s*3+ied_s0] == iflg_a &&
							tri_flg[itri_s*3+ied_s1] == iflg_a ){
							continue;
						}
						tri_flg[itri_s*3+ied_s] = iflg_a;
						tri_flg[itri*3+inotri] = 1-iflg_a;
						elem_msk[itri] = 1;
						break;
					}
				}
				if( elem_msk[itri] == 1 ) continue; // 上で成功すれば次の要素へ
				// 緩い条件で判定する。つまり、隣と入れ替えるときに、隣が駄目になっても構わない。
				for(unsigned int inotri=0;inotri<3;inotri++){
					if( aTri[itri].g2[inotri] == -2 ){
						const unsigned int itri_s = aTri[itri].s2[inotri];
						const unsigned int ied_s = relTriTri[ aTri[itri].r2[inotri] ][inotri];
						assert( aTri[itri_s].s2[ied_s] == itri );
						assert( tri_flg[itri_s*3+ied_s] == 1-iflg_a );
						if( elem_msk[itri_s] == 1 ) continue;
						tri_flg[itri_s*3+ied_s] = iflg_a;
						tri_flg[itri*3+inotri] = 1-iflg_a;
						elem_msk[itri] = 1;
						break;
					}
				}
			}
			if( icoun == 0 ){	// 計算終了
				icoun_pre = 0;
				break;
			}
			if( (jloop!=0||iloop!=0) && icoun < icoun_pre*0.5 ){ // 半分になったらフラグをクリアして計算しなおし
				icoun_pre = icoun;
				for(unsigned int itri=0;itri<aTri.size();itri++){
					elem_msk[itri] = 0;
				}
				break;
			}
			icoun_pre = icoun;
		}
			if( icoun_pre == 0 ) break;
		}
	} // elem_flgを作るルーティン終わり
	return true;
}

bool Msh::CMesh3D_Extrude::MakeConnectivitryTriTetForExtrusion(
	const CMesher2D& msh_2d,
	const std::vector<std::vector<int> >& tri_flg_ary, 
	const std::vector< std::vector<int> >& bar_flg_ary, 
	const unsigned int ndiv )
{
	const unsigned int nTriAry2d = tri_flg_ary.size();
	const unsigned int nBarAry2d = bar_flg_ary.size();
	const unsigned int nNo2d = msh_2d.GetVectorAry().size();

	{	// ３次元３角形を作る
		assert( nTriAry2d == tri_flg_ary.size() );
		assert( nBarAry2d == bar_flg_ary.size() );
		{	// 上下の３角形を作る
			const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
			assert( aTriAry2D.size() == 1 );
			const std::vector<STri2D>& aTri2d = aTriAry2D[0].m_aTri;
			const unsigned int nNo2d = msh_2d.GetVectorAry().size();
			{	// 下の三角形
				std::vector<STri3D>& aTri = m_aTriAry[0].m_aTri;
				for(unsigned int itri=0;itri<aTri.size();itri++){
					aTri[itri].v[0] = aTri2d[itri].v[2];
					aTri[itri].v[1] = aTri2d[itri].v[1];
					aTri[itri].v[2] = aTri2d[itri].v[0];
				}
			}
			{	// 上の三角形配列
				std::vector<STri3D>& aTri = m_aTriAry[1].m_aTri;
				for(unsigned int itri=0;itri<aTri.size();itri++){
					aTri[itri].v[0] = aTri2d[itri].v[0]+nNo2d*ndiv;
					aTri[itri].v[1] = aTri2d[itri].v[1]+nNo2d*ndiv;
					aTri[itri].v[2] = aTri2d[itri].v[2]+nNo2d*ndiv;
				}
			}
		}
		// 周囲の３角形を作る
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			const std::vector<int>& bar_flg = bar_flg_ary[ibar_ary];
			std::vector<STri3D>& aTri = m_aTriAry[nTriAry2d*2+ibar_ary].m_aTri;
			const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			const unsigned int nBar2d = aBar2d.size();
			for(unsigned int idiv=0;idiv<ndiv;idiv++){
			for(unsigned int ibar=0;ibar<bar_flg.size();ibar++){
				if( bar_flg[ibar] == 0 ){
					aTri[(idiv*nBar2d+ibar)*2  ].v[0] = aBar2d[ibar].v[0]+nNo2d*idiv;
					aTri[(idiv*nBar2d+ibar)*2  ].v[1] = aBar2d[ibar].v[1]+nNo2d*idiv;
					aTri[(idiv*nBar2d+ibar)*2  ].v[2] = aBar2d[ibar].v[1]+nNo2d*(idiv+1);

					aTri[(idiv*nBar2d+ibar)*2+1].v[0] = aBar2d[ibar].v[0]+nNo2d*idiv;
					aTri[(idiv*nBar2d+ibar)*2+1].v[1] = aBar2d[ibar].v[1]+nNo2d*(idiv+1);
					aTri[(idiv*nBar2d+ibar)*2+1].v[2] = aBar2d[ibar].v[0]+nNo2d*(idiv+1);
				}
				else{
					assert(0);
				}
			}
			}
		}
	}

	{	// 四面体を作るルーティン
		const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
		m_aTetAry.resize(nTriAry2d);
		assert( aTriAry2D.size() == 1 );
		const std::vector<STri2D>& aTri = aTriAry2D[0].m_aTri;
		const unsigned int nTri2d = aTri.size();
		std::vector<STet>& aTet = m_aTetAry[0].m_aTet;
		const std::vector<int>& tri_flg = tri_flg_ary[0];
		for(unsigned int idiv=0;idiv<ndiv;idiv++){
		for(unsigned int itri=0;itri<aTri.size();itri++){
			int iflg0 = tri_flg[itri*3  ];
			int iflg1 = tri_flg[itri*3+1];
			int iflg2 = tri_flg[itri*3+2];
			unsigned int ino0=0, ino1=0, ino2=0, ino3=0, ino4=0, ino5=0;
			if(      iflg0 == 1 && iflg1 == 0 && iflg2 == 0 ){ // m u d
				ino0=1; ino1=2; ino2=0; ino3=1; ino4=0; ino5=2; }
			else if( iflg0 == 0 && iflg1 == 1 && iflg2 == 0 ){ // d m u
				ino0=2; ino1=0; ino2=1; ino3=2; ino4=1; ino5=0; }
			else if( iflg0 == 0 && iflg1 == 0 && iflg2 == 1 ){ // u d m
				ino0=0; ino1=1; ino2=2; ino3=0; ino4=2; ino5=1; }
			else if( iflg0 == 0 && iflg1 == 1 && iflg2 == 1 ){ // m d u
				ino0=2; ino1=0; ino2=1; ino3=2; ino4=0; ino5=1; }
			else if( iflg0 == 1 && iflg1 == 0 && iflg2 == 1 ){ // u m d
				ino0=0; ino1=1; ino2=2; ino3=0; ino4=1; ino5=2; }
			else if( iflg0 == 1 && iflg1 == 1 && iflg2 == 0 ){ // d u m
				ino0=1; ino1=2; ino2=0; ino3=1; ino4=2; ino5=0; }
			else{ assert(0); }
			const unsigned int itri0 = nTri2d*idiv+itri;
			aTet[itri0*3+0].v[0] = nNo2d* idiv   +aTri[itri].v[0];
			aTet[itri0*3+0].v[1] = nNo2d* idiv   +aTri[itri].v[1];
			aTet[itri0*3+0].v[2] = nNo2d* idiv   +aTri[itri].v[2];
			aTet[itri0*3+0].v[3] = nNo2d*(idiv+1)+aTri[itri].v[ino0];

			aTet[itri0*3+1].v[0] = nNo2d* idiv   +aTri[itri].v[ino1];
			aTet[itri0*3+1].v[1] = nNo2d* idiv   +aTri[itri].v[ino2];
			aTet[itri0*3+1].v[2] = nNo2d*(idiv+1)+aTri[itri].v[ino3];
			aTet[itri0*3+1].v[3] = nNo2d*(idiv+1)+aTri[itri].v[ino4];

			aTet[itri0*3+2].v[0] = nNo2d*(idiv+1)+aTri[itri].v[2];
			aTet[itri0*3+2].v[1] = nNo2d*(idiv+1)+aTri[itri].v[1];
			aTet[itri0*3+2].v[2] = nNo2d*(idiv+1)+aTri[itri].v[0];
			aTet[itri0*3+2].v[3] = nNo2d*(idiv  )+aTri[itri].v[ino5];
		}
		}
		MakeTetSurTet(aTet);
		assert( CheckTet(aTet) );
	}

	return true;
}

// メッシュを突き出す
// height : メッシュ高さ
// elen : 高さ方向のメッシュ幅
bool Msh::CMesh3D_Extrude::Extrude_Tri(const CMesher2D& msh_2d, const double height, const double elen)
{
	if( elen <= 0 ) return false;
	if( height <= 0 ) return false;

	////////////////

	std::vector< std::vector<int> > tri_flg_ary;
	std::vector< std::vector<int> > bar_flg_ary;
	this->MakeDirectionFlagForExtrusion(msh_2d,tri_flg_ary,bar_flg_ary);

	const unsigned int nNo2d = msh_2d.GetVectorAry().size();
	unsigned int ndiv = (int)(height / elen);

	const unsigned int nVer2d = msh_2d.GetVertexAry().size();
	const unsigned int nTriAry2d = msh_2d.GetTriArySet().size();
	const unsigned int nBarAry2d = msh_2d.GetBarArySet().size();

	std::vector<unsigned int> id_ary_3d;
	{
		const unsigned int elem_size = 
			 nVer2d*2+				// 点 : 下＋上
			 nBarAry2d+nVer2d*2+	// 辺 : 垂直辺、水平辺(下＋上)
			 nTriAry2d*2+nBarAry2d+	// 面 : 下＋上＋垂直面
			 nTriAry2d;				// 体
		id_ary_3d.resize( elem_size, 0 );
	}

	{ // 頂点を作る
		const std::vector<SVertex>& aVer2d = msh_2d.GetVertexAry();
		m_aVertex.resize( nVer2d*2 );
		for(unsigned int iver=0;iver<nVer2d;iver++){    // 底面を作る
			m_aVertex[iver].id_cad = aVer2d[iver].id_v_cad;
			m_aVertex[iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aVertex[iver].inum_extrude = 1;
			const unsigned int id0 = GetFreeObjID(); 
			m_aVertex[iver].id = id0; 
			m_aVertex[iver].v = aVer2d[iver].v;
			id_ary_3d[iver] = id0;
		}
		for(unsigned int iver=0;iver<nVer2d;iver++){    // 上面を作る
			m_aVertex[nVer2d+iver].id_cad = aVer2d[iver].id_v_cad;
			m_aVertex[nVer2d+iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aVertex[nVer2d+iver].inum_extrude = 3;
			const unsigned int id0 = GetFreeObjID();
			m_aVertex[nVer2d+iver].id = id0;
			m_aVertex[nVer2d+iver].v = aVer2d[iver].v+nNo2d*ndiv;
			id_ary_3d[nVer2d+iver] = id0;
		}
	}

	{	// 辺を作る
		m_aBarAry.resize( nBarAry2d*2+nVer2d );
		// 下の辺の領域を確保
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			m_aBarAry[ibar_ary].id_cad = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aBarAry[ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aBarAry[ibar_ary].inum_extrude = 1;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[ibar_ary].id = id0;
			const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
            {
			    std::vector<SBar>& aBar_d = m_aBarAry[ibar_ary*2  ].m_aBar;
			    aBar_d.resize(aBar2d.size());
			    for(unsigned int ibar=0;ibar<aBar2d.size();ibar++){
				    aBar_d[ibar].v[0] = aBar2d[ibar].v[0];
				    aBar_d[ibar].v[1] = aBar2d[ibar].v[1];
			    }
            }
			id_ary_3d[nVer2d*2+ibar_ary] = id0;
		}
		// 上の辺の領域を確保
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			m_aBarAry[nBarAry2d+ibar_ary].id_cad = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aBarAry[nBarAry2d+ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aBarAry[nBarAry2d+ibar_ary].inum_extrude = 3;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[nBarAry2d+ibar_ary].id = id0;
			const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
            {
			    std::vector<SBar>& aBar_u = m_aBarAry[ibar_ary*2+1].m_aBar;
			    aBar_u.resize(aBar2d.size());
			    for(unsigned int ibar=0;ibar<aBar2d.size();ibar++){
				    aBar_u[ibar].v[0] = aBar2d[ibar].v[0]+nNo2d*ndiv;
				    aBar_u[ibar].v[1] = aBar2d[ibar].v[1]+nNo2d*ndiv;
			    }
            }
			id_ary_3d[nVer2d*2+nBarAry2d+ibar_ary] = id0;
		}
		// 周囲の辺を確保
		const std::vector<SVertex>& aVer2d = msh_2d.GetVertexAry();
		for(unsigned int iver=0;iver<nVer2d;iver++){
    		m_aBarAry[nBarAry2d*2+iver].id_cad = aVer2d[iver].id_v_cad;
			m_aBarAry[nBarAry2d*2+iver].id_msh_before_extrude = aVer2d[iver].id;
			m_aBarAry[nBarAry2d*2+iver].inum_extrude = 2;
			const unsigned int id0 = GetFreeObjID();
			m_aBarAry[nBarAry2d*2+iver].id = id0;
            {
			    std::vector<SBar>& aBar = m_aBarAry[nBarAry2d*2+iver].m_aBar;
			    aBar.resize(ndiv);
			    const unsigned int iv2d = aVer2d[iver].v;
			    for(unsigned int idiv=0;idiv<ndiv;idiv++){
				    aBar[idiv].v[0] = iv2d+nNo2d*idiv;
				    aBar[idiv].v[1] = iv2d+nNo2d*(idiv+1);
			    }
            }
			id_ary_3d[nVer2d*2+nBarAry2d*2+iver] = id0;
		}
	}
	
	{	// ３次元３角形の領域の確保
		assert( nTriAry2d == tri_flg_ary.size() );
		assert( nBarAry2d == bar_flg_ary.size() );
		m_aTriAry.resize(nTriAry2d*2+nBarAry2d);
		{	// 上下の３角形
			const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
			assert( aTriAry2D.size() == 1 );
			const std::vector<STri2D>& aTri2d = aTriAry2D[0].m_aTri;
//			const unsigned int nNo2d = msh_2d.GetVectorAry().size();
			{	// 下の三角形
				m_aTriAry[0].id_cad = aTriAry2D[0].id_l_cad;
				m_aTriAry[0].id_msh_before_extrude = aTriAry2D[0].id;
				m_aTriAry[0].inum_extrude = 1;
				const unsigned int id0 = GetFreeObjID();
				m_aTriAry[0].id = id0;
				std::vector<STri3D>& aTri = m_aTriAry[0].m_aTri;
				aTri.resize( aTri2d.size() );
				id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d  ] = id0;
			}
			{	// 上の三角形配列
				m_aTriAry[1].id_cad = aTriAry2D[0].id_l_cad;
				m_aTriAry[1].id_msh_before_extrude = aTriAry2D[0].id;
				m_aTriAry[1].inum_extrude = 3;
				const unsigned int id0 = GetFreeObjID();
				m_aTriAry[1].id = id0;
				std::vector<STri3D>& aTri = m_aTriAry[1].m_aTri;
				aTri.resize( aTri2d.size() );
				id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d] = id0;
			}
		}
		// 周囲の３角形
		for(unsigned int ibar_ary=0;ibar_ary<nBarAry2d;ibar_ary++){
			const std::vector<SBar>& aBar2d = msh_2d.GetBarArySet()[ibar_ary].m_aBar;
			const unsigned int nBar2d = aBar2d.size();
            {
			    std::vector<STri3D>& aTri = m_aTriAry[nTriAry2d*2+ibar_ary].m_aTri;
			    aTri.resize( nBar2d*ndiv*2 );
            }
			m_aTriAry[nTriAry2d*2+ibar_ary].id_cad = msh_2d.GetBarArySet()[ibar_ary].id_e_cad;
			m_aTriAry[nTriAry2d*2+ibar_ary].id_msh_before_extrude = msh_2d.GetBarArySet()[ibar_ary].id;
			m_aTriAry[nTriAry2d*2+ibar_ary].inum_extrude = 2;
			const unsigned int id0 = GetFreeObjID();
			m_aTriAry[nTriAry2d*2+ibar_ary].id = id0;
			id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d*2+ibar_ary] = id0;
		}
	}

	
	{	// 四面体の領域を確保する
		const std::vector<CTriAry2D>& aTriAry2D = msh_2d.GetTriArySet();
		assert( aTriAry2D.size() == 1 );
		m_aTetAry.resize(nTriAry2d);
		const std::vector<STri2D>& aTri = aTriAry2D[0].m_aTri;
//		const unsigned int nTri2d = aTri.size();
		m_aTetAry[0].id_cad = aTriAry2D[0].id_l_cad;
		m_aTetAry[0].id_msh_before_extrude = aTriAry2D[0].id;
		m_aTetAry[0].inum_extrude = 2;
		const unsigned int id0 = GetFreeObjID();
		m_aTetAry[0].id = id0;
		std::vector<STet>& aTet = m_aTetAry[0].m_aTet;
		aTet.resize( aTri.size()*3*ndiv );
		id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d*2+nBarAry2d  ] = id0;
	}


	////////////////
	// 三角形と四面体のConnectivityを作る
	this->MakeConnectivitryTriTetForExtrusion(msh_2d,tri_flg_ary,bar_flg_ary,ndiv);

	{	// 節点を作るルーティン
		double hlen = height / ndiv;
		const unsigned int nNo2d = msh_2d.GetVectorAry().size();
		this->aVec.resize((ndiv+1)*nNo2d);
		const std::vector<CVector2D>& aVec2D = msh_2d.GetVectorAry();
		for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
		for(unsigned int ino=0;ino<nNo2d;ino++){
			this->aVec[idiv*nNo2d+ino].x = aVec2D[ino].x;
			this->aVec[idiv*nNo2d+ino].y = aVec2D[ino].y;
			this->aVec[idiv*nNo2d+ino].z = hlen*idiv;
		}
		}
	}

	{
		m_include_relation.resize( this->FindMaxID()+1 );
		for(unsigned int iver2d=0;iver2d<nBarAry2d;iver2d++){
			const unsigned int id_bar_v_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+iver2d];	// 垂直辺
			m_include_relation[id_bar_v_3d].push_back( id_ary_3d[iver2d] );
			m_include_relation[id_bar_v_3d].push_back( id_ary_3d[nVer2d+iver2d] );
		}
		for(unsigned int ibar2d=0;ibar2d<nBarAry2d;ibar2d++){
			const unsigned int id_bar_d_3d = id_ary_3d[nVer2d*2+  ibar2d];
			const unsigned int id_bar_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d+ibar2d];
			const unsigned int id_tri_3d   = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d*2+ibar2d];

			const unsigned int id_bar_2d = msh_2d.GetBarArySet()[ibar2d].id;
			std::vector<unsigned int> id_ary_inc = msh_2d.GetIncludeElemIDAry(id_bar_2d);
			for(unsigned int iinc=0;iinc<id_ary_inc.size();iinc++){
				const unsigned int id_inc = id_ary_inc[iinc];
				unsigned int iver2d;
				{
					unsigned int nelem, id_cad;
					MSH_TYPE msh_type;
					msh_2d.GetMshInfo(id_inc,nelem,msh_type,iver2d, id_cad);
					assert( msh_type == VERTEX );
					assert( iver2d < nVer2d );
				}
				m_include_relation[id_bar_d_3d].push_back( id_ary_3d[iver2d] );	// 下の辺
				m_include_relation[id_bar_u_3d].push_back( id_ary_3d[nVer2d+iver2d] );	// 上の辺
				m_include_relation[id_tri_3d  ].push_back( id_ary_3d[iver2d] );		// 垂直面
				m_include_relation[id_tri_3d  ].push_back( id_ary_3d[nVer2d+iver2d] );	// 垂直面
				m_include_relation[id_tri_3d  ].push_back( id_ary_3d[nVer2d*2+nBarAry2d*2+iver2d] );	// 垂直面
			}
			m_include_relation[id_tri_3d].push_back( id_ary_3d[nVer2d*2+ibar2d] );	// 垂直面
			m_include_relation[id_tri_3d].push_back( id_ary_3d[nVer2d*2+nBarAry2d+ibar2d] );	// 垂直面
		}
		for(unsigned int itri2d=0;itri2d<nTriAry2d;itri2d++){
			const unsigned int id_tri_d_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  itri2d];
			const unsigned int id_tri_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d+itri2d];
			const unsigned int id_tet_3d   = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d*2+nBarAry2d+  itri2d];

			const unsigned int id_tri_2d = msh_2d.GetTriArySet()[itri2d].id;
			std::vector<unsigned int> id_ary_inc = msh_2d.GetIncludeElemIDAry(id_tri_2d);
			for(unsigned int iinc=0;iinc<id_ary_inc.size();iinc++){
				const unsigned int id_inc = id_ary_inc[iinc];
				unsigned int iloc2d;
				unsigned int nelem, id_cad;
				MSH_TYPE msh_type;
				msh_2d.GetMshInfo(id_inc,nelem,msh_type,iloc2d, id_cad);
				if( msh_type == VERTEX ){
					assert( iloc2d < nVer2d );
					const unsigned int id_ver_d_3d = id_ary_3d[iloc2d];
					const unsigned int id_ver_u_3d = id_ary_3d[nVer2d+iloc2d];
					const unsigned int id_bar_v_3d   = id_ary_3d[nVer2d*2+  nBarAry2d*2+iloc2d];
					m_include_relation[id_tri_d_3d].push_back( id_ver_d_3d );
					m_include_relation[id_tri_u_3d].push_back( id_ver_u_3d );
					m_include_relation[id_tet_3d  ].push_back( id_ver_d_3d );
					m_include_relation[id_tet_3d  ].push_back( id_ver_u_3d );
					m_include_relation[id_tet_3d  ].push_back( id_bar_v_3d );
				}
				else if( msh_type == BAR ){
					assert( iloc2d < nBarAry2d );
					const unsigned int id_bar_d_3d = id_ary_3d[nVer2d*2+  iloc2d];
					const unsigned int id_bar_u_3d = id_ary_3d[nVer2d*2+  nBarAry2d+iloc2d];
					const unsigned int id_tri_v_3d = id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d*2+iloc2d];
					m_include_relation[id_tri_d_3d].push_back( id_bar_d_3d );
					m_include_relation[id_tri_u_3d].push_back( id_bar_u_3d );
					m_include_relation[id_tet_3d  ].push_back( id_bar_d_3d );
					m_include_relation[id_tet_3d  ].push_back( id_bar_u_3d );
					m_include_relation[id_tet_3d  ].push_back( id_tri_v_3d );
				}
				else{
					assert(0);
				}
			}
			m_include_relation[id_tet_3d  ].push_back( id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  itri2d] );
			m_include_relation[id_tet_3d  ].push_back( id_ary_3d[nVer2d*2+  nBarAry2d*2+nVer2d+  nTriAry2d+itri2d] );
		}
	}

	this->MakeElemLocationType();

	return true;
}


bool Msh::CMesh3D_Extrude::UpdateMeshConnectivity(const CMesher2D& msh_2d){

	std::vector< std::vector<int> > tri_flg_ary;
	std::vector< std::vector<int> > bar_flg_ary;
	this->MakeDirectionFlagForExtrusion(msh_2d,tri_flg_ary,bar_flg_ary);

	const unsigned int nNo2d = msh_2d.GetVectorAry().size();
	const unsigned int nNo3d = this->GetVectorAry().size();
	assert( nNo3d % nNo2d == 0 );
	const unsigned int ndiv= nNo3d / nNo2d - 1;

	// 三角形と四面体のConnectivityを作る
	this->MakeConnectivitryTriTetForExtrusion(msh_2d,tri_flg_ary,bar_flg_ary,ndiv);

	return true;
}

