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
#pragma warning( disable : 4786 )   // C4786なんて表示すんな( ﾟДﾟ)ｺﾞﾙｧ
#endif
#define for if(0);else for

#include <iostream>
#include <set>
#include <vector>
#include <stack>
#include <cassert>
#include <math.h>

#include "delfem/msh/meshkernel2d.h"

using namespace Com;

static const unsigned int invRelTriTri[3] = {
	0, 1, 2
};

// (0に相当するノード番号)*3+(1に相当するノード番号)  →→　関係番号
static const int noel2RelTriTri[9] = {
	-1,	// 0 00
	-1,	// 1 01
	 0,	// 2 02
	 2, // 3 10
	-1, // 4 11
	-1,	// 5 12
	-1,	// 6 20
	 1,	// 7 21
	-1, // 8 22
};

// (こちら側の辺番号)*3+(相手側の辺番号)　→→　関係番号
static const unsigned int ed2RelTriTri[9] = {
	0,	// 0 00
	2,	// 1 01
	1,	// 2 02
	2,	// 3 10
	1,	// 4 11
	0,	// 5 12
	1,	// 6 20
	0,	// 7 21
	2,	// 8 22
}; 

static const unsigned int indexRot3[3][3] = {
	{ 0, 1, 2 },
	{ 1, 2, 0 },
	{ 2, 0, 1 },
};

bool Msh::CheckTri( const std::vector<STri2D>& aTri )
{
	const unsigned int ntri = aTri.size();
	for(unsigned int itri=0;itri<ntri;itri++){
		const STri2D& ref_tri = aTri[itri];
/*		for(int inotri=0;inotri<nNoTri;inotri++){
			assert( ref_tri.v[inotri] >= 0 );
		}*/
		for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
			if( ref_tri.g2[iedtri] == -2 || ref_tri.g2[iedtri] == -3 ){
				const unsigned int itri_s = ref_tri.s2[iedtri];
				const unsigned int irel = ref_tri.r2[iedtri];
				assert( itri_s < ntri );
				assert( irel < 3 );
				// check sorounding
				{
					const unsigned int noel_dia = relTriTri[irel][iedtri];
					assert( aTri[itri_s].s2[noel_dia] == itri );
//					std::cout << itri << " " << itri_s << std::endl;
				}
				// check relation 
				for(unsigned int inoed=0;inoed<nNoEd;inoed++){
					const unsigned int inoel = noelTriEdge[iedtri][inoed];
                    assert( ref_tri.v[inoel] == aTri[itri_s].v[ (int)relTriTri[irel][inoel] ] );
				}
			}
		}
	}

	return true;
}

bool Msh::CheckTri(
		const std::vector<CPoint2D>& po,
		const std::vector<STri2D>& tri )
{
//	std::cout << "Check Tri" << std::endl;

	const unsigned int npo = po.size();
	const unsigned int ntri = tri.size();

    ////////////////////////////////
    // 要素Indexのチェック

	for(unsigned int itri=0;itri<ntri;itri++){
		const STri2D& ref_tri = tri[itri];
        for(unsigned int inotri=0;inotri<nNoTri;inotri++){
			assert( ref_tri.v[inotri] < npo );
		}
		for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
			if( ref_tri.g2[iedtri] == -2 || ref_tri.g2[iedtri] == -3 ){
				const unsigned int itri_s = ref_tri.s2[iedtri];
				const unsigned int irel = ref_tri.r2[iedtri];
				assert( itri_s < ntri );
				assert( irel < 3 );
				// check sorounding
				{
					const unsigned int noel_dia = relTriTri[irel][iedtri];
					assert( noel_dia < 3 );
					if( tri[itri_s].s2[noel_dia] != itri ){
						std::cout << itri << " " << iedtri << std::endl;
					}
					assert( tri[itri_s].s2[noel_dia] == itri );
				}
				// check relation 
				for(unsigned int inoed=0;inoed<nNoEd;inoed++){
					const unsigned int inoel = noelTriEdge[iedtri][inoed];
                    if( ref_tri.v[inoel] != tri[itri_s].v[ (int)relTriTri[irel][inoel] ] ){
						std::cout << itri << " " << iedtri << std::endl;
					}
                    assert( ref_tri.v[inoel] == tri[itri_s].v[ (int)relTriTri[irel][inoel] ] );
				}
            }
		}
/*		{
			if( ref_tri.g2[0]==-1 && ref_tri.g2[1]==-1 && ref_tri.g2[2]==-1 ){
				std::cout << "Isolated Triangle " << itri << std::endl;
			}
		}*/
	}

    ////////////////////////////////
    // 頂点-要素間の一貫性のチェック

    for(unsigned int ipoin=0;ipoin<npo;ipoin++){
        if( po[ipoin].e >= 0 ){
            assert( po[ipoin].d >= 0 && po[ipoin].d < 3 );
            const int itri0 = po[ipoin].e;
            const int inoel0 = po[ipoin].d;
            if( tri[itri0].v[inoel0] != ipoin ){
//				std::cout << itri0 << " " << inoel0 << "   " << tri[itri0].v[inoel0] << " " << ipoin << std::endl;
            }
            assert( tri[itri0].v[inoel0] == ipoin );
        }
    }

    ////////////////////////////////
    // Geometryのチェック

	for(unsigned int itri=0;itri<ntri;itri++){
		const STri2D& ref_tri = tri[itri];
		{
			double area = Com::TriArea( po[ref_tri.v[0]].p, po[ref_tri.v[1]].p, po[ref_tri.v[2]].p);
			if( area < 1.0e-10 ){
//				std::cout << "Negative Volume : " << itri << " " << area << std::endl;
			}
		}
/*		{	
			CVector2D v0 = po[ref_tri.v[0]].p;
			CVector2D v1 = po[ref_tri.v[1]].p;
			CVector2D v2 = po[ref_tri.v[2]].p;

			double area = TriArea( v0, v1, v2 );
			const double tmp1 = 0.5 / area;

			double const_term[3];
			const_term[0] = tmp1*(v1.x*v2.y-v2.x*v1.y);
			const_term[1] = tmp1*(v2.x*v0.y-v0.x*v2.y);
			const_term[2] = tmp1*(v0.x*v1.y-v1.x*v0.y);

			double dldx[2][2];
			dldx[0][0] = tmp1*(v1.y-v2.y);
			dldx[1][0] = tmp1*(v2.y-v0.y);
			dldx[2][0] = tmp1*(v0.y-v1.y);

			dldx[0][1] = tmp1*(v2.x-v1.x);
			dldx[1][1] = tmp1*(v0.x-v2.x);
			dldx[2][1] = tmp1*(v1.x-v0.x);

			assert( fabs( dldx[0][0]+dldx[1][0]+dldx[2][0] ) < 1.0e-15 );
			assert( fabs( dldx[0][1]+dldx[1][1]+dldx[2][1] ) < 1.0e-15 );

			assert( fabs( const_term[0]+dldx[0][0]*v0.x+dldx[0][1]*v0.y - 1.0 ) < 1.0e-10 );
			assert( fabs( const_term[0]+dldx[0][0]*v1.x+dldx[0][1]*v1.y ) < 1.0e-10 );
			assert( fabs( const_term[0]+dldx[0][0]*v2.x+dldx[0][1]*v2.y ) < 1.0e-10 );

			assert( fabs( const_term[1]+dldx[1][0]*v0.x+dldx[1][1]*v0.y ) < 1.0e-10 );
			assert( fabs( const_term[1]+dldx[1][0]*v1.x+dldx[1][1]*v1.y - 1.0 ) < 1.0e-10 );
			assert( fabs( const_term[1]+dldx[1][0]*v2.x+dldx[1][1]*v2.y ) < 1.0e-10 );

			assert( fabs( const_term[2]+dldx[2][0]*v0.x+dldx[2][1]*v0.y ) < 1.0e-10 );
			assert( fabs( const_term[2]+dldx[2][0]*v1.x+dldx[2][1]*v1.y ) < 1.0e-10 );
			assert( fabs( const_term[2]+dldx[2][0]*v2.x+dldx[2][1]*v2.y - 1.0 ) < 1.0e-10 );
		}*/
	}

	return true;
}

bool Msh::InsertPoint_ElemEdge( const unsigned int ipo_ins, 
							   const unsigned int itri_ins, const unsigned int ied_ins,
							   std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	assert( itri_ins < tri.size() );
	assert( ipo_ins < po.size() );

	if( tri[itri_ins].g2[ied_ins] != -2 ){
//		std::cout << "未実装" << std::endl;
		assert(0);
	}

	const unsigned int itri_adj = tri[itri_ins].s2[ied_ins];
    const unsigned int ied_adj  = (unsigned int)relTriTri[ (int)tri[itri_ins].r2[ied_ins] ][ied_ins];
	assert( itri_adj < tri.size() );
	assert( ied_ins < 3 );

    const unsigned int itri0 = itri_ins;
    const unsigned int itri1 = itri_adj;
    const unsigned int itri2 = tri.size();
    const unsigned int itri3 = tri.size()+1;

	tri.resize( tri.size()+2 );

	STri2D old0 = tri[itri_ins];
	STri2D old1 = tri[itri_adj];

	const unsigned int ino0_0 = ied_ins;
	const unsigned int ino1_0 = noelTriEdge[ied_ins][0];
	const unsigned int ino2_0 = noelTriEdge[ied_ins][1];

	const unsigned int ino0_1 = ied_adj;
	const unsigned int ino1_1 = noelTriEdge[ied_adj][0];
	const unsigned int ino2_1 = noelTriEdge[ied_adj][1];
	
	assert( old0.v[ino1_0] == old1.v[ino2_1] );
	assert( old0.v[ino2_0] == old1.v[ino1_1] );
	assert( old0.s2[ino0_0 ] == itri1 );
	assert( old1.s2[ino0_1 ] == itri0 );

	po[ipo_ins].e = itri0;			po[ipo_ins].d = 0;
	po[ old0.v[ino2_0] ].e = itri0;	po[ old0.v[ino2_0] ].d = 1;
	po[ old0.v[ino0_0] ].e = itri1;	po[ old0.v[ino0_0] ].d = 1;
	po[ old1.v[ino2_1] ].e = itri2;	po[ old1.v[ino2_1] ].d = 1;
	po[ old1.v[ino0_1] ].e = itri3;	po[ old1.v[ino0_1] ].d = 1;

	{
		STri2D& ref_tri = tri[itri0];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old0.v[ino2_0];	ref_tri.v[2]  = old0.v[ino0_0];
		ref_tri.g2[0] = old0.g2[ino1_0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old0.s2[ino1_0];	ref_tri.s2[1] = itri1;			ref_tri.s2[2] = itri3;
		////////////////
		if( old0.g2[ino1_0] == -2 || old0.g2[ino1_0] == -3 ){
			assert( old0.r2[ino1_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[ino1_0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino1_0]*3 + rel[ino2_0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old0.s2[ino1_0] < tri.size() );
            tri[ old0.s2[ino1_0] ].s2[ rel[ino1_0] ] = itri0;
            tri[ old0.s2[ino1_0] ].r2[ rel[ino1_0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri1];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old0.v[ino0_0];	ref_tri.v[2]  = old0.v[ino1_0];
		ref_tri.g2[0] = old0.g2[ino2_0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old0.s2[ino2_0];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;
		////////////////
		if( old0.g2[ino2_0] == -2 || old0.g2[ino2_0] == -3 ){
			assert( old0.r2[ino2_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[ino2_0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino2_0]*3 + rel[ino0_0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old0.s2[ino2_0] < tri.size() );
            tri[ old0.s2[ino2_0] ].s2[ rel[ino2_0] ] = itri1;
            tri[ old0.s2[ino2_0] ].r2[ rel[ino2_0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri2];
		////////////////
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old1.v[ino2_1];	ref_tri.v[2]  = old1.v[ino0_1];
		ref_tri.g2[0] = old1.g2[ino1_1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old1.s2[ino1_1];	ref_tri.s2[1] = itri3;			ref_tri.s2[2] = itri1;
		////////////////
		if( old1.g2[ino1_1] == -2 || old0.g2[ino2_0] == -3 ){
			assert( old1.r2[ino1_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[ino1_1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino1_1]*3 + rel[ino2_1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old1.s2[ino1_1] < tri.size() );
            tri[ old1.s2[ino1_1] ].s2[ rel[ino1_1] ] = itri2;
            tri[ old1.s2[ino1_1] ].r2[ rel[ino1_1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri3];
		ref_tri.v[0]  = ipo_ins;			ref_tri.v[1]  = old1.v[ino0_1];	ref_tri.v[2]  = old1.v[ino1_1];
		ref_tri.g2[0] = old1.g2[ino2_1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old1.s2[ino2_1];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri2;
		if( old1.g2[ino2_1] == -2 || old1.g2[ino2_1] == -3 ){
			assert( old1.r2[ino2_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[ino2_1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[ino2_1]*3 + rel[ino0_1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old1.s2[ino2_1] < tri.size() );
			tri[ old1.s2[ino2_1] ].s2[ rel[ino2_1] ] = itri3;
			tri[ old1.s2[ino2_1] ].r2[ rel[ino2_1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	return true;
}

// 辺[ipo0-ipo1]の左側の３角形itri0を探索する
// 三角形がなければ->falseを返す。
// 三角形があれば  ->true を返す。
// 但し、その場合
// tri[itri0].v[inotri0]==ipo0
// tri[itri0].v[inotri1]==ipo1
// を満たす
bool Msh::FindEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1,
	const std::vector<CPoint2D>& po, const std::vector<STri2D>& tri )
{
	const unsigned int itri_ini = po[ipo0].e;
	const unsigned int inotri_ini = po[ipo0].d;
	unsigned int inotri_cur = inotri_ini;
	unsigned int itri_cur = itri_ini;
	for(;;){	//　時計周りに検索する。
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{	// この要素がOKか調べる
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].v[inotri2] == ipo1 ){
				itri0 = itri_cur;
				inotri0 = inotri_cur;
				inotri1 = inotri2;
				assert( tri[itri0].v[ inotri0 ] == ipo0 );
				assert( tri[itri0].v[ inotri1 ] == ipo1 );
				return true;
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){ break; }
			const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
            const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ) return false;
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	inotri_cur = inotri_ini;
	itri_cur = itri_ini;
	for(;;){	//　反時計周りの検索
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 || tri[itri_cur].g2[inotri2] != -3 ){ break; }
			const unsigned int itri_nex = tri[itri_cur].s2[inotri2];
            const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){	// 一周したら終わり
				itri0 = 0;
				inotri0 = 0; inotri1 = 0;
				return false;
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
		{	// 要素の向きを調べる
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].v[inotri2] == ipo1 ){
				itri0 = itri_cur;
				inotri0 = inotri_cur;
				inotri1 = inotri2;
				assert( tri[itri0].v[ inotri0 ] == ipo0 );
				assert( tri[itri0].v[ inotri1 ] == ipo1 );
				return true;
			}
		}
	}

	return false;
}


bool Msh::FlipEdge( unsigned int itri0, unsigned int ied0, std::vector<STri2D>& aTri )
{
	assert( itri0 < aTri.size() );
	assert( ied0 < 3 );
	assert( aTri[itri0].g2[ied0] == -2 );

	const unsigned int itri1 = aTri[itri0].s2[ied0];
	const unsigned int ied1  = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
	assert( itri1 < aTri.size() );
	assert( ied1 < 3 );
	assert( aTri[itri1].g2[ied1] == -2 );

//	std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;

	STri2D old0 = aTri[itri0];
	STri2D old1 = aTri[itri1];

	const unsigned int no0_0 = ied0;
	const unsigned int no1_0 = noelTriEdge[ied0][0];
	const unsigned int no2_0 = noelTriEdge[ied0][1];

	const unsigned int no0_1 = ied1;
	const unsigned int no1_1 = noelTriEdge[ied1][0];
	const unsigned int no2_1 = noelTriEdge[ied1][1];
	
	assert( old0.v[no1_0] == old1.v[no2_1] );
	assert( old0.v[no2_0] == old1.v[no1_1] );

	{
		STri2D& ref_tri = aTri[itri0];
		////////////////
		ref_tri.v[0]  = old0.v[no1_0];	ref_tri.v[1]  = old1.v[no0_1];	ref_tri.v[2]  = old0.v[no0_0];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old0.g2[no2_0];	ref_tri.g2[2] = old1.g2[no1_1];
		ref_tri.s2[0] = itri1;			ref_tri.s2[1] = old0.s2[no2_0];	ref_tri.s2[2] = old1.s2[no1_1];
		////////////////
		ref_tri.r2[0] = 0;
		if( old0.g2[no2_0] == -2 || old0.g2[no2_0] == -3 ){
			assert( old0.r2[no2_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[no2_0] ];
			assert( old0.s2[no2_0] < aTri.size() );
			assert( old0.s2[no2_0] != itri0 );
			assert( old0.s2[no2_0] != itri1 );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_0]*3 + rel[no2_0] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old0.s2[no2_0] ].s2[ rel[no2_0] ] = itri0;
			aTri[ old0.s2[no2_0] ].r2[ rel[no2_0] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old1.g2[no1_1] == -2 || old1.g2[no1_1] == -3 ){
			assert( old1.r2[no1_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[no1_1] ];
			assert( old1.s2[no1_1] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_1]*3 + rel[no0_1] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old1.s2[no1_1] ].s2[ rel[no1_1] ] = itri0;
			aTri[ old1.s2[no1_1] ].r2[ rel[no1_1] ] = invRelTriTri[ ref_tri.r2[2] ];			
		}
	}

	{
		STri2D& ref_tri = aTri[itri1];
		////////////////
		ref_tri.v[0] = old1.v[no1_1];	ref_tri.v[1]  = old0.v[no0_0];	ref_tri.v[2]  = old1.v[no0_1];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old1.g2[no2_1];	ref_tri.g2[2] = old0.g2[no1_0];
		ref_tri.s2[0] = itri0;			ref_tri.s2[1] = old1.s2[no2_1];	ref_tri.s2[2] = old0.s2[no1_0];
		////////////////
		ref_tri.r2[0] = 0;
		if( old1.g2[no2_1] == -2 || old1.g2[no2_1] == -3 ){
			assert( old1.r2[no2_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[no2_1] ];
			assert( old1.s2[no2_1] < aTri.size() );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_1]*3 + rel[no2_1] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old1.s2[no2_1] ].s2[ rel[no2_1] ] = itri1;
			aTri[ old1.s2[no2_1] ].r2[ rel[no2_1] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old0.g2[no1_0] == -2 || old0.g2[no1_0] == -3 ){
			assert( old0.r2[no1_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[no1_0] ];
			assert( old0.s2[no1_0] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_0]*3 + rel[no0_0] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old0.s2[no1_0] ].s2[ rel[no1_0] ] = itri1;
			aTri[ old0.s2[no1_0] ].r2[ rel[no1_0] ] = invRelTriTri[ ref_tri.r2[2] ];
		}
	}
	
	return true;
}

bool Msh::FlipEdge( unsigned int itri0, unsigned int ied0,
			  std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri )
{
	assert( itri0 < aTri.size() );
	assert( ied0 < 3 );
	assert( aTri[itri0].g2[ied0] == -2 );

	const unsigned int itri1 = aTri[itri0].s2[ied0];
	const unsigned int ied1  = relTriTri[ aTri[itri0].r2[ied0] ][ied0];
	assert( itri1 < aTri.size() );
	assert( ied1 < 3 );
	assert( aTri[itri1].g2[ied1] == -2 );

//	std::cout << itri0 << "-" << ied0 << "    " << itri1 << "-" << ied1 << std::endl;

	STri2D old0 = aTri[itri0];
	STri2D old1 = aTri[itri1];

	const unsigned int no0_0 = ied0;
	const unsigned int no1_0 = noelTriEdge[ied0][0];
	const unsigned int no2_0 = noelTriEdge[ied0][1];

	const unsigned int no0_1 = ied1;
	const unsigned int no1_1 = noelTriEdge[ied1][0];
	const unsigned int no2_1 = noelTriEdge[ied1][1];
	
	assert( old0.v[no1_0] == old1.v[no2_1] );
	assert( old0.v[no2_0] == old1.v[no1_1] );

	aPo[ old0.v[no1_0] ].e = itri0;	aPo[ old0.v[no1_0] ].d = 0;
	aPo[ old0.v[no0_0] ].e = itri0;	aPo[ old0.v[no0_0] ].d = 2;
	aPo[ old1.v[no1_1] ].e = itri1;	aPo[ old1.v[no1_1] ].d = 0;
	aPo[ old1.v[no0_1] ].e = itri1;	aPo[ old1.v[no0_1] ].d = 2;

	{
		STri2D& ref_tri = aTri[itri0];
		////////////////
		ref_tri.v[0]  = old0.v[no1_0];	ref_tri.v[1]  = old1.v[no0_1];	ref_tri.v[2]  = old0.v[no0_0];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old0.g2[no2_0];	ref_tri.g2[2] = old1.g2[no1_1];
		ref_tri.s2[0] = itri1;			ref_tri.s2[1] = old0.s2[no2_0];	ref_tri.s2[2] = old1.s2[no1_1];
		////////////////
		ref_tri.r2[0] = 0;
		if( old0.g2[no2_0] == -2 || old0.g2[no2_0] == -3 ){
			assert( old0.r2[no2_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[no2_0] ];
			assert( old0.s2[no2_0] < aTri.size() );
			assert( old0.s2[no2_0] != itri0 );
			assert( old0.s2[no2_0] != itri1 );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_0]*3 + rel[no2_0] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old0.s2[no2_0] ].s2[ rel[no2_0] ] = itri0;
			aTri[ old0.s2[no2_0] ].r2[ rel[no2_0] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old1.g2[no1_1] == -2 || old1.g2[no1_1] == -3 ){
			assert( old1.r2[no1_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[no1_1] ];
			assert( old1.s2[no1_1] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_1]*3 + rel[no0_1] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old1.s2[no1_1] ].s2[ rel[no1_1] ] = itri0;
			aTri[ old1.s2[no1_1] ].r2[ rel[no1_1] ] = invRelTriTri[ ref_tri.r2[2] ];			
		}
	}

	{
		STri2D& ref_tri = aTri[itri1];
		////////////////
		ref_tri.v[0] = old1.v[no1_1];	ref_tri.v[1]  = old0.v[no0_0];	ref_tri.v[2]  = old1.v[no0_1];
		ref_tri.g2[0] = -2;				ref_tri.g2[1] = old1.g2[no2_1];	ref_tri.g2[2] = old0.g2[no1_0];
		ref_tri.s2[0] = itri0;			ref_tri.s2[1] = old1.s2[no2_1];	ref_tri.s2[2] = old0.s2[no1_0];
		////////////////
		ref_tri.r2[0] = 0;
		if( old1.g2[no2_1] == -2 || old1.g2[no2_1] == -3 ){
			assert( old1.r2[no2_1] < 3 );
            const unsigned int* rel = relTriTri[ old1.r2[no2_1] ];
			assert( old1.s2[no2_1] < aTri.size() );
			ref_tri.r2[1] = noel2RelTriTri[ rel[no1_1]*3 + rel[no2_1] ];
			assert( ref_tri.r2[1] >= 0 && ref_tri.r2[1] < 3 );
			aTri[ old1.s2[no2_1] ].s2[ rel[no2_1] ] = itri1;
			aTri[ old1.s2[no2_1] ].r2[ rel[no2_1] ] = invRelTriTri[ ref_tri.r2[1] ];
		}
		if( old0.g2[no1_0] == -2 || old0.g2[no1_0] == -3 ){
			assert( old0.r2[no1_0] < 3 );
            const unsigned int* rel = relTriTri[ old0.r2[no1_0] ];
			assert( old0.s2[no1_0] < aTri.size() );
			ref_tri.r2[2] = noel2RelTriTri[ rel[no2_0]*3 + rel[no0_0] ];
			assert( ref_tri.r2[2] >= 0 && ref_tri.r2[2] < 3 );
			aTri[ old0.s2[no1_0] ].s2[ rel[no1_0] ] = itri1;
			aTri[ old0.s2[no1_0] ].r2[ rel[no1_0] ] = invRelTriTri[ ref_tri.r2[2] ];
		}
	}
	
	return true;
}

bool Msh::InsertPoint_Elem( const unsigned int ipo_ins, const unsigned int itri_ins, 
					  std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	assert( itri_ins < tri.size() );
	assert( ipo_ins < po.size() );

	const int itri0 = itri_ins;
	const int itri1 = tri.size();
	const int itri2 = tri.size()+1;

	tri.resize( tri.size()+2 );

	const STri2D old_tri = tri[itri_ins];

	po[ipo_ins].e = itri0;			po[ipo_ins].d = 0;
	po[ old_tri.v[0] ].e = itri1;	po[ old_tri.v[0] ].d = 2;
	po[ old_tri.v[1] ].e = itri2;	po[ old_tri.v[1] ].d = 2;
	po[ old_tri.v[2] ].e = itri0;	po[ old_tri.v[2] ].d = 2;

	{
		STri2D& ref_tri = tri[itri0];

		ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[1];	ref_tri.v[2] = old_tri.v[2];
		ref_tri.g2[0] = old_tri.g2[0];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[0];	ref_tri.s2[1] = itri1;			ref_tri.s2[2] = itri2;

		if( old_tri.g2[0] == -2 || old_tri.g2[0] == -3 ){
			assert( old_tri.r2[0] < 3 );
            const unsigned int* rel = relTriTri[ old_tri.r2[0] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[0]*3 + rel[1] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[0] < tri.size() );
			tri[ old_tri.s2[0] ].s2[ rel[0] ] = itri0;
			tri[ old_tri.s2[0] ].r2[ rel[0] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri1];

		ref_tri.v[0] = ipo_ins;			ref_tri.v[1] = old_tri.v[2];	ref_tri.v[2] = old_tri.v[0];
		ref_tri.g2[0] = old_tri.g2[1];	ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[1];	ref_tri.s2[1] = itri2;			ref_tri.s2[2] = itri0;

		if( old_tri.g2[1] == -2 || old_tri.g2[1] == -3 ){
			assert( old_tri.r2[1] < 3 );
            const unsigned int* rel = relTriTri[ old_tri.r2[1] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[1]*3 + rel[2] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[1] < tri.size() );
			tri[ old_tri.s2[1] ].s2[ rel[1] ] = itri1;
			tri[ old_tri.s2[1] ].r2[ rel[1] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}
	{
		STri2D& ref_tri = tri[itri2];

		ref_tri.v[0] = ipo_ins; 		ref_tri.v[1] = old_tri.v[0];	ref_tri.v[2] = old_tri.v[1];
		ref_tri.g2[0] = old_tri.g2[2];  ref_tri.g2[1] = -2;				ref_tri.g2[2] = -2;
		ref_tri.s2[0] = old_tri.s2[2];	ref_tri.s2[1] = itri0;			ref_tri.s2[2] = itri1;

		if( old_tri.g2[2] == -2 || old_tri.g2[2] == -3 ){
			assert( old_tri.r2[2] < 3 );
            const unsigned int* rel = relTriTri[ old_tri.r2[2] ];
			ref_tri.r2[0] = noel2RelTriTri[ rel[2]*3 + rel[0] ];
			assert( ref_tri.r2[0] >= 0 && ref_tri.r2[0] < 3 );
			assert( old_tri.s2[2] < tri.size() );
			tri[ old_tri.s2[2] ].s2[ rel[2] ] = itri2;
			tri[ old_tri.s2[2] ].r2[ rel[2] ] = invRelTriTri[ ref_tri.r2[0] ];
		}
		ref_tri.r2[1] = 0;
		ref_tri.r2[2] = 0;
	}

	return true;
}

bool Msh::DelaunayAroundPoint( const unsigned int itri0, unsigned int inotri0, 
	std::vector<Com::CVector2D>& aVec, std::vector<STri2D>& aTri, 
	unsigned int& num_flip)
{
	assert( itri0 < aTri.size() );
	assert( inotri0 < 3 );
	const unsigned int ivec0 = aTri[itri0].v[inotri0];
	assert( ivec0 < aVec.size() );

	unsigned int itri_cur = itri0;
	unsigned int inotri_cur = inotri0;
	bool flag_is_wall = false;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ivec0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
            const unsigned int* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( Com::DetDelaunay(
				aVec[ aTri[itri_cur].v[0] ],
				aVec[ aTri[itri_cur].v[1] ],
				aVec[ aTri[itri_cur].v[2] ],
				aVec[ ipo_dia ] ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aTri);	// 辺を切り替える
				inotri_cur = 2;
				assert( aTri[itri_cur].v[inotri_cur] == ivec0 );
				num_flip++;
				// FlipによってaTri[itri0].v[inotri0] != ivec0 でなくなってしまうのを防ぐため
				if( itri_cur == itri0 ) inotri0 = inotri_cur;
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri1 = indexRot3[1][inotri_cur];
			if( aTri[itri_cur].g2[inotri1] != -2 && aTri[itri_cur].g2[inotri1] != -3 ){
				flag_is_wall = true;
				break;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri1];
            const unsigned int* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri1] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ivec0 );
			if( itri_nex == itri0 ) break;	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}
	if( !flag_is_wall ) return true;

	////////////////////////////////
	// 逆向きへの回転

	itri_cur = itri0;
	inotri_cur = inotri0;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ivec0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
            const unsigned int* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( Com::DetDelaunay(
				aVec[ aTri[itri_cur].v[0] ],
				aVec[ aTri[itri_cur].v[1] ],
				aVec[ aTri[itri_cur].v[2] ],
				aVec[ ipo_dia ] ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aTri);	// 辺を切り替える
				itri_cur = itri_dia;
				inotri_cur = 1;
				num_flip++;
				assert( aTri[itri_cur].v[inotri_cur] == ivec0 );
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( aTri[itri_cur].g2[inotri2] != -2 && aTri[itri_cur].g2[inotri2] != -3 ){
				return true;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri2];
            const unsigned int* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri2] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ivec0 );
			assert( itri_nex != itri0 );	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}
	return true;
}

bool Msh::DelaunayAroundPoint( unsigned int ipo0,
	std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri )
{
	assert( ipo0 < aPo.size() );
	if( aPo[ipo0].e == -1 ) return true;

	assert( aPo[ipo0].e >= 0 && (unsigned int)aPo[ipo0].e < aTri.size() );
	assert( aTri[ aPo[ipo0].e ].v[ aPo[ipo0].d ] == ipo0 );

	const unsigned int itri0 = aPo[ipo0].e;
	unsigned int inotri0 = aPo[ipo0].d;

	unsigned int itri_cur = itri0;
	unsigned int inotri_cur = aPo[ipo0].d;
	bool flag_is_wall = false;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ipo0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
            const unsigned int* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( Com::DetDelaunay(
				aPo[ aTri[itri_cur].v[0] ].p,
				aPo[ aTri[itri_cur].v[1] ].p,
				aPo[ aTri[itri_cur].v[2] ].p,
				aPo[ ipo_dia ].p ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aPo,aTri);	// 辺を切り替える
				// FlipEdgeによってitri_curは時計回り側の３角形に切り替わる
				inotri_cur = 2;
				assert( aTri[itri_cur].v[inotri_cur] == ipo0 );
				// FlipによってaTri[itri0].v[inotri0] != ipo0 でなくなってしまうのを防ぐため
				if( itri_cur == itri0 ) inotri0 = inotri_cur;
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri1 = indexRot3[1][inotri_cur];
			if( aTri[itri_cur].g2[inotri1] != -2 && aTri[itri_cur].g2[inotri1] != -3 ){
				flag_is_wall = true;
				break;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri1];
            const unsigned int* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri1] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ipo0 );
			if( itri_nex == itri0 ) break;	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}
	if( !flag_is_wall ) return true;

	////////////////////////////////
	// 逆向きへの回転

	itri_cur = itri0;
	inotri_cur = inotri0;
	for(;;){
		assert( aTri[itri_cur].v[inotri_cur] == ipo0 );

		if( aTri[itri_cur].g2[inotri_cur] == -2 ){
			// 向かいの要素を調べる
			const unsigned int itri_dia = aTri[itri_cur].s2[inotri_cur];
            const unsigned int* rel_dia = relTriTri[ aTri[itri_cur].r2[inotri_cur] ];
			const unsigned int inotri_dia = rel_dia[inotri_cur];
			assert( aTri[itri_dia].g2[inotri_dia] == -2 );
			assert( aTri[itri_dia].s2[inotri_dia] == itri_cur );
			const unsigned int ipo_dia = aTri[itri_dia].v[inotri_dia];
			if( Com::DetDelaunay(
				aPo[ aTri[itri_cur].v[0] ].p,
				aPo[ aTri[itri_cur].v[1] ].p,
				aPo[ aTri[itri_cur].v[2] ].p,
				aPo[ ipo_dia ].p ) == 0 )	// Delaunay条件が満たされない場合
			{
				FlipEdge(itri_cur,inotri_cur,aPo,aTri);	// 辺を切り替える
				itri_cur = itri_dia;
				inotri_cur = 1;
				assert( aTri[itri_cur].v[inotri_cur] == ipo0 );
				continue;	// ループの始めに戻る
			}
		}

		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( aTri[itri_cur].g2[inotri2] != -2 && aTri[itri_cur].g2[inotri2] != -3 ){
				return true;
			}
			const unsigned int itri_nex = aTri[itri_cur].s2[inotri2];
            const unsigned int* rel_nex = relTriTri[ aTri[itri_cur].r2[inotri2] ];
			const unsigned int inotri_nex = rel_nex[inotri_cur];
			assert( aTri[itri_nex].v[inotri_nex] == ipo0 );
			assert( itri_nex != itri0 );	// 一周したら終わり
			itri_cur = itri_nex;
			inotri_cur = inotri_nex;
		}
	}

	return true;
}

bool Msh::FindEdgePoint_AcrossEdge( const unsigned int& ipo0, const unsigned int& ipo1,
	unsigned int& itri0, unsigned int& inotri0, unsigned int& inotri1, double& ratio,
	std::vector<CPoint2D>& po, std::vector<STri2D>& tri )
{
	const unsigned int itri_ini = po[ipo0].e;
	const unsigned int inotri_ini = po[ipo0].d;
	unsigned int inotri_cur = inotri_ini;
	unsigned int itri_cur = itri_ini;
	for(;;){	//　反時計周りの検索
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			const unsigned int inotri3 = indexRot3[2][inotri_cur];
			double area0 = Com::TriArea( po[ipo0].p, po[ tri[itri_cur].v[inotri2] ].p, po[ipo1].p );
			if( area0 > -1.0e-20 ){
				double area1 =  Com::TriArea( po[ipo0].p, po[ipo1].p, po[ tri[itri_cur].v[inotri3] ].p );
				if( area1 > -1.0e-20 ){
					assert( area0 + area1 > 1.0e-20 );
					ratio = area0 / ( area0 + area1 );
					itri0 = itri_cur;
					inotri0 = inotri2;
					inotri1 = inotri3;
					return true;
				}
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
				break;
			}
			unsigned int itri_nex = tri[itri_cur].s2[inotri2];
            const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){	// 一周したら終わり
				itri0 = 0;
				inotri0 = 0; inotri1 = 0;
				ratio = 0.0;
				return false;
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	// itri_iniを２回計算しているので少し無駄、いつか直す

	inotri_cur = inotri_ini;
	itri_cur = itri_ini;
	for(;;){	//　時計周りに検索する。
		assert( tri[itri_cur].v[inotri_cur] == ipo0 );
		{
			const unsigned int inotri2 = indexRot3[1][inotri_cur];
			const unsigned int inotri3 = indexRot3[2][inotri_cur];
			double area0 = Com::TriArea( po[ipo0].p, po[ tri[itri_cur].v[inotri2] ].p, po[ipo1].p );
			if( area0 > -1.0e-20 ){
				double area1 =  Com::TriArea( po[ipo0].p, po[ipo1].p, po[ tri[itri_cur].v[inotri3] ].p );
				if( area1 > -1.0e-20 ){
					assert( area0 + area1 > 1.0e-20 );
					ratio = area0 / ( area0 + area1 );
					itri0 = itri_cur;
					inotri0 = inotri2;
					inotri1 = inotri3;
					return true;
				}
			}
		}
		{	// 次の要素へ進める
			const unsigned int inotri2 = indexRot3[2][inotri_cur];
			if( tri[itri_cur].g2[inotri2] != -2 && tri[itri_cur].g2[inotri2] != -3 ){
				break;
			}
			unsigned int itri_nex = tri[itri_cur].s2[inotri2];
            const unsigned int* rel = relTriTri[ tri[itri_cur].r2[inotri2] ];
			const unsigned int inotri3 = rel[inotri_cur];
			assert( tri[itri_nex].v[inotri3] == ipo0 );
			if( itri_nex == itri_ini ){
				assert(0);	// 一周しないはず
			}
			itri_cur = itri_nex;
			inotri_cur = inotri3;
		}
	}

	// 失敗したときの値を入れる
	itri0 = 0;
	inotri0 = 0; inotri1 = 0;
	ratio = 0.0;

	return false;
}

bool Msh::MakePointSurBar( const std::vector<SBar>& aBar, const unsigned int npoin, 
					 unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup )
{
	const unsigned int nnobar = 2;

	for(unsigned int ipoin=0;ipoin<npoin+1;ipoin++){
		elsup_ind[ipoin] = 0;
	}
	for(unsigned int ibar=0;ibar<aBar.size();ibar++){
		for(unsigned int inobar=0;inobar<nnobar;inobar++){
			elsup_ind[ aBar[ibar].v[inobar]+1 ]++;
		}
	}
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		elsup_ind[ipoin+1] += elsup_ind[ipoin];
	}
	nelsup = elsup_ind[npoin];
	elsup = new unsigned int [nelsup];
	for(unsigned int ibar=0;ibar<aBar.size();ibar++){
		for(unsigned int inobar=0;inobar<nnobar;inobar++){
			const unsigned int ipoin0 = aBar[ibar].v[inobar];
			const unsigned int ielsup = elsup_ind[ipoin0];
			elsup[ielsup] = ibar;
			elsup_ind[ipoin0]++;
		}
	}
	for(int ipoin=npoin;ipoin>0;ipoin--){
		elsup_ind[ipoin] = elsup_ind[ipoin-1];
	}
	elsup_ind[0] = 0;

	return true;
}

bool Msh::MakePointSurTri( const std::vector<STri2D>& aTri, const unsigned int npoin, 
					  unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup ){

	const unsigned int nnotri = 3;

	for(unsigned int ipoin=0;ipoin<npoin+1;ipoin++){
		elsup_ind[ipoin] = 0;
	}
	for(unsigned int itri=0;itri<aTri.size();itri++){
		for(unsigned int inotri=0;inotri<nnotri;inotri++){
			elsup_ind[ aTri[itri].v[inotri]+1 ]++;
		}
	}
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		elsup_ind[ipoin+1] += elsup_ind[ipoin];
	}
	nelsup = elsup_ind[npoin];
	elsup = new unsigned int [nelsup];
	for(unsigned int itri=0;itri<aTri.size();itri++){
		for(unsigned int inotri=0;inotri<nnotri;inotri++){
			const unsigned int ipoin0 = aTri[itri].v[inotri];
			const unsigned int ielsup = elsup_ind[ipoin0];
			elsup[ielsup] = itri;
			elsup_ind[ipoin0]++;
		}
	}
	for(int ipoin=npoin;ipoin>0;ipoin--){
		elsup_ind[ipoin] = elsup_ind[ipoin-1];
	}
	elsup_ind[0] = 0;
/*
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		std::cout << ipoin << " ";
		for(unsigned int ielsup=elsup_ind[ipoin];ielsup<elsup_ind[ipoin+1];ielsup++){
			std::cout << elsup[ielsup] << " ";
		}
		std::cout << std::endl;
	}
*/
	return true;
}

bool Msh::MakePointSurQuad( const std::vector<SQuad2D>& aQuad, const unsigned int npoin, 
					  unsigned int* const elsup_ind, unsigned int& nelsup, unsigned int*& elsup )
{
	const unsigned int nnoquad = 4;

	for(unsigned int ipoin=0;ipoin<npoin+1;ipoin++){
		elsup_ind[ipoin] = 0;
	}
	for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
		for(unsigned int inoquad=0;inoquad<nnoquad;inoquad++){
			elsup_ind[ aQuad[iquad].v[inoquad]+1 ]++;
		}
	}
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		elsup_ind[ipoin+1] += elsup_ind[ipoin];
	}
	nelsup = elsup_ind[npoin];
	elsup = new unsigned int [nelsup];
	for(unsigned int iquad=0;iquad<aQuad.size();iquad++){
		for(unsigned int inoquad=0;inoquad<nnoquad;inoquad++){
			const unsigned int ipoin0 = aQuad[iquad].v[inoquad];
			const unsigned int ielsup = elsup_ind[ipoin0];
			elsup[ielsup] = iquad;
			elsup_ind[ipoin0]++;
		}
	}
	for(int ipoin=npoin;ipoin>0;ipoin--){
		elsup_ind[ipoin] = elsup_ind[ipoin-1];
	}
	elsup_ind[0] = 0;
/*
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		std::cout << ipoin << " ";
		for(unsigned int ielsup=elsup_ind[ipoin];ielsup<elsup_ind[ipoin+1];ielsup++){
			std::cout << elsup[ielsup] << " ";
		}
		std::cout << std::endl;
	}
*/
	return true;
}
/*
bool Msh::MakeInnerRelationBar( std::vector<SBar>& aBar ) // std::vector<aVertex>の入力をつける
{
	unsigned int npoin;
	{	// 属する点の最大値＋１を求める
		npoin=0;
		for(unsigned int ibar=0;ibar<aBar.size();ibar++){
			npoin = ( aBar[ibar].v[0] > npoin ) ? aBar[ibar].v[0] : npoin;
			npoin = ( aBar[ibar].v[1] > npoin ) ? aBar[ibar].v[1] : npoin;
		}
		npoin = npoin+1;
	}

	unsigned int* elsup_ind = new unsigned int [npoin+1];
	unsigned int nelsup;
	unsigned int* elsup;

	Msh::MakePointSurBar(aBar,npoin,elsup_ind,nelsup,elsup);

	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		const unsigned int nelsup_i = elsup_ind[ipoin+1]-elsup_ind[ipoin];
		if( nelsup_i != 0 && nelsup_i != 2 ) return false;
	}

	for(unsigned int ibar=0;ibar<aBar.size();ibar++){
		for(unsigned int inoed=0;inoed<nNoEd;inoed++){
			const unsigned int ipoin0 = aBar[ibar].v[inoed];
			assert( elsup_ind[ipoin0+1]-elsup_ind[ipoin0] == 2 );
//			aBar[ibar].g1[1-inoed] = -2;
			const unsigned int ibar0 = elsup[ elsup_ind[ipoin0]   ];
			const unsigned int ibar1 = elsup[ elsup_ind[ipoin0]+1 ];
			unsigned int ibar_s;
			if( ibar0 == ibar ){ ibar_s = ibar1; }
			else{
				assert( ibar1 == ibar );
				ibar_s = ibar0;
			}
			aBar[ibar].s1[1-inoed] = ibar_s;
			if( aBar[ibar].v[inoed] == aBar[ibar_s].v[inoed] ){
				aBar[ibar].r1[1-inoed] = 1;
			}
			else{
				assert( aBar[ibar].v[inoed] == aBar[ibar_s].v[1-inoed] );
				aBar[ibar].r1[1-inoed] = 0;
			}
		}
	}

	delete[] elsup_ind;
	delete[] elsup;

	return true;
}
*/

bool Msh::MakeInnerRelationTri( std::vector<STri2D>& aTri, const unsigned int npoin, 
	const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup )
{
	const unsigned int EdEd2Rel[nEdTri][nEdTri] = {
        { 0, 2, 1 },
        { 2, 1, 0 },
        { 1, 0, 2 } };

	unsigned int* tmp_poin = new unsigned int [npoin];
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){ tmp_poin[ipoin] = 0; }
	unsigned int inpofa[2];

	const unsigned int nTri = aTri.size();
	for(unsigned int itri=0;itri<nTri;itri++){
		for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
			for(unsigned int ipoed=0;ipoed<nNoEd;ipoed++){
				inpofa[ipoed] = aTri[itri].v[ noelTriEdge[iedtri][ipoed] ];
				tmp_poin[ inpofa[ipoed] ] = 1;
			}
			const unsigned int ipoin0= inpofa[0];
			bool iflg = false;
			for(unsigned int ielsup=elsup_ind[ipoin0];ielsup<elsup_ind[ipoin0+1];ielsup++){
				const unsigned int jtri0 = elsup[ielsup];
				if( jtri0 == itri ) continue;
				for(unsigned int jedtri=0;jedtri<nEdTri;jedtri++){
					iflg = true;
					for(unsigned int jpoed=0;jpoed<nNoEd;jpoed++){
						const unsigned int jpoin0 =  aTri[jtri0].v[ noelTriEdge[jedtri][jpoed] ];
						if( tmp_poin[ jpoin0 ] == 0 ){ iflg = false; break; }
					}
					if( iflg ){
						aTri[itri].g2[iedtri] = -2;
						aTri[itri].s2[iedtri] = jtri0;
						aTri[itri].r2[iedtri] = EdEd2Rel[iedtri][jedtri];
						break;
					}
				}
				if( iflg ) break;
			}
			if( !iflg ){ 
				aTri[itri].g2[iedtri] = -1;
			}
			for(unsigned int ipofa=0;ipofa<nNoEd;ipofa++){
				tmp_poin[ inpofa[ipofa] ] = 0;
			}
		}
	}

	delete[] tmp_poin;
	return true;
}

bool Msh::MakeInnerRelationQuad( std::vector<SQuad2D>& aQuad, const unsigned int npoin, 
	const unsigned int* elsup_ind, const unsigned int nelsup, const unsigned int* elsup )
{
	// エッジ番号から関係番号へのテーブル
	const unsigned int EdEd2Rel[nEdQuad][nEdQuad] = {
        { 1, 2, 3, 0 },
        { 2, 3, 0, 1 },
        { 3, 0, 1, 2 },
        { 0, 1, 2, 3 } };

	unsigned int* tmp_poin = new unsigned int [npoin];
	for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
		tmp_poin[ipoin] = 0;
	}
	unsigned int inpofa[2];

	const unsigned int nQuad = aQuad.size();
	for(unsigned int iquad=0;iquad<nQuad;iquad++){
		for(unsigned int iedquad=0;iedquad<nEdQuad;iedquad++){
			for(unsigned int ipoed=0;ipoed<nNoEd;ipoed++){
				inpofa[ipoed] = aQuad[iquad].v[ noelQuadEdge[iedquad][ipoed] ];
				tmp_poin[ inpofa[ipoed] ] = 1;
			}
			const unsigned int ipoin0= inpofa[0];
			bool iflg = false;
			for(unsigned int ielsup=elsup_ind[ipoin0];ielsup<elsup_ind[ipoin0+1];ielsup++){
				const unsigned int jquad0 = elsup[ielsup];
				if( jquad0 == iquad ) continue;
				for(unsigned int jedquad=0;jedquad<nEdQuad;jedquad++){
					iflg = true;
					for(unsigned int jpoed=0;jpoed<nNoEd;jpoed++){
						const unsigned int jpoin0 =  aQuad[jquad0].v[ noelQuadEdge[jedquad][jpoed] ];
						if( tmp_poin[ jpoin0 ] == 0 ){
							iflg = false;
							break;
						}
					}
					if( iflg ){
						aQuad[iquad].g2[iedquad] = -2;
						aQuad[iquad].s2[iedquad] = jquad0;
						aQuad[iquad].r2[iedquad] = EdEd2Rel[iedquad][jedquad];
						break;
					}
				}
				if( iflg ) break;
			}
			if( !iflg ){ 
				aQuad[iquad].g2[iedquad] = -1;
			}
			for(unsigned int ipofa=0;ipofa<nNoEd;ipofa++){
				tmp_poin[ inpofa[ipofa] ] = 0;
			}
		}
	}

	delete[] tmp_poin;
	return true;
}

bool Msh::MakeOuterBoundTri( const std::vector<STri2D>& aTri, std::vector<SBar>& aBar )
{
	unsigned int counter;
	{
		counter = 0;
		const unsigned int ntri = aTri.size();
		for(unsigned int itri=0;itri<ntri;itri++){
			for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
				if( aTri[itri].g2[iedtri] != -2 ) counter++;
			}
		}
	}
//	const unsigned int nbar = counter;
	aBar.reserve(counter);
	const unsigned int ntri = aTri.size();
	for(unsigned int itri=0;itri<ntri;itri++){
	for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
		if( aTri[itri].g2[iedtri] != -2 ){
			unsigned int ibar0 = aBar.size();
			aBar.resize( aBar.size()+1 );
			aBar[ibar0].v[0] = aTri[itri].v[ noelTriEdge[iedtri][0] ];
			aBar[ibar0].v[1] = aTri[itri].v[ noelTriEdge[iedtri][1] ];
		}
	}
	}
	return true;
}

bool Msh::MakeOuterBoundQuad( const std::vector<SQuad2D>& aQuad, std::vector<SBar>& aBar )
{
	unsigned int counter;
	{
		counter = 0;
		const unsigned int nquad = aQuad.size();
		for(unsigned int iquad=0;iquad<nquad;iquad++){
			for(unsigned int iedquad=0;iedquad<nEdQuad;iedquad++){
				if( aQuad[iquad].g2[iedquad] != -2 ) counter++;
			}
		}
	}
//	const unsigned int nbar = counter;
	aBar.reserve(counter);
	const unsigned int nquad = aQuad.size();
	for(unsigned int iquad=0;iquad<nquad;iquad++){
	for(unsigned int iedquad=0;iedquad<nEdQuad;iedquad++){
		if( aQuad[iquad].g2[iedquad] != -2 ){
			unsigned int ibar0 = aBar.size();
			aBar.resize( aBar.size()+1 );
			aBar[ibar0].v[0] = aQuad[iquad].v[ noelQuadEdge[iedquad][0] ];
		aBar[ibar0].v[1] = aQuad[iquad].v[ noelQuadEdge[iedquad][1] ];
		}
	}
	}
	return true;
}

void Msh::LaplacianSmoothing( std::vector<CPoint2D>& aPo, const std::vector<STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move)
{	
	for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++){	// 点周りの点を探索して調べる。
		if( ipoin < aflg_isnt_move.size() ){
			if( aflg_isnt_move[ipoin] == 1 ) continue;
		}
		const unsigned int itri_ini = aPo[ipoin].e;
		const unsigned int inoel_c_ini = aPo[ipoin].d;
		assert( itri_ini < aTri.size() );
		assert( inoel_c_ini < 3 );
		assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
		unsigned int itri0= itri_ini;
		unsigned int inoel_c0 = inoel_c_ini;
		unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
		bool is_bound_flg = false;
		Com::CVector2D vec_delta = aPo[ipoin].p;
		unsigned int ntri_around = 1;
		for(;;){
			assert( itri0 < aTri.size() );
			assert( inoel_c0 < 3 );
			assert( aTri[itri0].v[inoel_c0] == ipoin );
			{
				vec_delta.x += aPo[ aTri[itri0].v[inoel_b0] ].p.x;
				vec_delta.y += aPo[ aTri[itri0].v[inoel_b0] ].p.y;
				ntri_around++;
			}
			if( aTri[itri0].g2[inoel_b0] == -2 ){
				unsigned int itri1 = aTri[itri0].s2[inoel_b0];
                const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
                const unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
                const unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
				assert( itri1 < aTri.size() );
				assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
				assert( aTri[itri1].v[inoel_c1] == ipoin );
				if( itri1 == itri_ini ) break;
				itri0 = itri1;
				inoel_c0 = inoel_c1;
				inoel_b0 = inoel_b1;
			}
			else{	// この点は境界上の点だから動かしてはならない。
				is_bound_flg = true;
				break;
			}
		}
		if( is_bound_flg ) continue;		
		aPo[ipoin].p.x = vec_delta.x / ntri_around;
		aPo[ipoin].p.y = vec_delta.y / ntri_around;
	}
}


void Msh::LaplaceDelaunaySmoothing( std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri,
	const std::vector<unsigned int>& aflg_isnt_move )
{
	for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++){	// 点周りの点を探索して調べる。
		if( ipoin < aflg_isnt_move.size() ){
			if( aflg_isnt_move[ipoin] == 1 ) continue;
		}
		const unsigned int itri_ini = aPo[ipoin].e;
		const unsigned int inoel_c_ini = aPo[ipoin].d;
		assert( itri_ini < aTri.size() );
		assert( inoel_c_ini < 3 );
		assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
		unsigned int itri0= itri_ini;
		unsigned int inoel_c0 = inoel_c_ini;
		unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
		bool is_bound_flg = false;
		Com::CVector2D vec_delta = aPo[ipoin].p;
		unsigned int ntri_around = 1;
		for(;;){	// 点の周りの要素を一回りする
			assert( itri0 < aTri.size() );
			assert( inoel_c0 < 3 );
			assert( aTri[itri0].v[inoel_c0] == ipoin );
			{
				vec_delta.x += aPo[ aTri[itri0].v[inoel_b0] ].p.x;
				vec_delta.y += aPo[ aTri[itri0].v[inoel_b0] ].p.y;
				ntri_around++;
			}
			if( aTri[itri0].g2[inoel_b0] == -2 ){
				unsigned int itri1 = aTri[itri0].s2[inoel_b0];
                const unsigned  rel01 = aTri[itri0].r2[inoel_b0];
                const unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
                const unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
				assert( itri1 < aTri.size() );
				assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
				assert( aTri[itri1].v[inoel_c1] == ipoin );
				if( itri1 == itri_ini ) break;
				itri0 = itri1;
				inoel_c0 = inoel_c1;
				inoel_b0 = inoel_b1;
			}
			else{	// この点は境界上の点だから動かしてはならない。
				is_bound_flg = true;
				break;
			}
		}
		if( is_bound_flg ) continue;		
		aPo[ipoin].p.x = vec_delta.x / ntri_around;
		aPo[ipoin].p.y = vec_delta.y / ntri_around;
		DelaunayAroundPoint(ipoin,aPo,aTri);
	}
}



// TODO : aTriEdを消去する
bool Msh::Tesselate_GiftWrapping(
	unsigned int ipo_st, unsigned int ipo_ed,
	std::vector< std::pair<unsigned int,unsigned int> >& aTriEd,
	const std::vector<CVector2D>& aVec, std::vector<STri2D>& aTri)
{
	assert( ipo_st < aVec.size() );
	assert( ipo_ed < aVec.size() );

	if( aTriEd.size() < aVec.size() ){
		aTriEd.resize( aVec.size() );
	}

	const unsigned int npo = aVec.size();

	// edからstの間の節点をを反時計周りで調べる
	unsigned int ipo0 = ipo_ed;
	unsigned int ipo1 = ipo_st;
	for(;;){
		int ipo_min = -1;
    double r_min = 0;
		for(unsigned int ipo=(ipo0==npo-1)?0:ipo0+1;ipo!=ipo1;ipo=(ipo==npo-1)?0:ipo+1){
			if( Com::TriArea(aVec[ipo1],aVec[ipo0],aVec[ipo]) < 1.0e-10 ) continue;
			if( ipo_min == -1 ){
				r_min = SquareCircumradius(aVec[ipo0],aVec[ipo1],aVec[ipo]);
				ipo_min = ipo;
//				std::cout << "    " << ipo << " " << r_min << "  |||  " << ipo0 << " " << ipo1 << " " << npo << std::endl;
			}
			else{
				const double r = SquareCircumradius(aVec[ipo0],aVec[ipo1],aVec[ipo]);
//				std::cout << "    " << ipo << " " << r << "  |||  " << ipo0 << " " << ipo1 << " " << npo << std::endl;
				if( r < r_min ){ 
					r_min = r;
					ipo_min = ipo;
				}
			}
		}
		if( ipo_min == -1 ){
			std::cout << "Triangle Fail " << std::endl;
			return false;
		}

//		std::cout << "next point " << ipo_min << std::endl;

		unsigned int itri0 = aTri.size();
		aTri.resize( aTri.size()+1 );
		{
			aTri[itri0].v[0] = ipo1;
			aTri[itri0].v[1] = ipo0;
			aTri[itri0].v[2] = ipo_min;
			if( aTri.size() == 1 ){
				aTri[itri0].g2[2] = -1;
				aTri[itri0].s2[2] = 0;
				aTri[itri0].r2[2] = 0;
			}
			else{
				int itri1=-1, iedtri1=-1;
				for(unsigned int itri=0;itri<aTri.size()-1;itri++){
				for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
					unsigned int jpo1 = aTri[itri].v[ noelTriEdge[iedtri][0] ];
					unsigned int jpo0 = aTri[itri].v[ noelTriEdge[iedtri][1] ];
					if( jpo1 == ipo0 && jpo0 == ipo1 ){ 
						itri1 = itri;
						iedtri1 = iedtri;
						break;
					}
				}
				}
				assert( itri1 != -1 && iedtri1 != -1 );
				aTri[itri0].g2[2] = -2;
				aTri[itri0].s2[2] = itri1;
                const unsigned int irel0 = noel2RelTriTri[ noelTriEdge[iedtri1][1]*3+noelTriEdge[iedtri1][0] ];
				aTri[itri0].r2[2] = irel0;
				{
					assert( aTri[itri1].g2[iedtri1] == -2 );
					assert( aTri[itri1].s2[iedtri1] == itri0 );
					assert( aTri[itri1].r2[iedtri1] == invRelTriTri[irel0] );
                    const unsigned int* rel1 = relTriTri[ aTri[itri1].r2[iedtri1] ];
					assert( rel1[iedtri1] == 2 );
					assert( rel1[ noelTriEdge[iedtri1][0] ] == 1 );
					assert( rel1[ noelTriEdge[iedtri1][1] ] == 0 );
				}
			}
		}

    if(      ipo_min==(int)ipo0+1 || (ipo_min==0&&ipo0==npo-1) ){
			aTri[itri0].g2[0] = -1;
			aTri[itri0].s2[0] = ipo0;
			aTriEd[ipo0].first = itri0;
			aTriEd[ipo0].second = 0;
//			std::cout << "left  " << ipo_min << " " << ipo1 << std::endl;
      if( ipo_min==(int)ipo1-1 || (ipo_min==(int)npo-1&&ipo1==0) ){
				aTri[itri0].g2[1] = -1;
				aTri[itri0].s2[1] = ipo_min;
				aTriEd[ipo_min].first = itri0;
				aTriEd[ipo_min].second = 1;
				break;
			}
			aTri[itri0].g2[1] = -2;
			aTri[itri0].s2[1] = aTri.size();
			aTri[itri0].r2[1] = noel2RelTriTri[ 0*3+2 ];
			ipo0 = ipo_min;
		}
    else if( ipo_min==(int)ipo1-1 || (ipo_min==(int)npo-1&&ipo1==0) ){
//			std::cout << "right " << ipo0 << " " << ipo_min << std::endl;
			aTri[itri0].g2[1] = -1;
			aTri[itri0].s2[1] = ipo_min;
			aTriEd[ipo_min].first = itri0;
			aTriEd[ipo_min].second = 1;
      if( ipo_min==(int)ipo0+1 || (ipo_min==0&&ipo0==npo-1) ){
				aTri[itri0].g2[0] = -1;
				aTri[itri0].s2[0] = ipo0;
				aTriEd[ipo0].first = itri0;
				aTriEd[ipo0].second = 0;
				break;
			}
			aTri[itri0].g2[0] = -2;
			aTri[itri0].s2[0] = aTri.size();
			aTri[itri0].r2[0] = noel2RelTriTri[ 2*3+1 ];
			ipo1 = ipo_min;
		}
		else{
//			std::cout << "center right : " << ipo0 << " " << ipo_min << std::endl;
			aTri[itri0].g2[0] = -2;
			aTri[itri0].s2[0] = aTri.size();
			aTri[itri0].r2[0] = noel2RelTriTri[ 2*3+1 ];
			////////////////
			aTriEd[ipo0].first = itri0;
			aTriEd[ipo0].second = 0;
			if( !Tesselate_GiftWrapping(ipo_min,ipo0,aTriEd,aVec,aTri) ) return false;
			////////////////////////////////////////////////////////////////
//			std::cout << "center left  : " << ipo_min << " " << ipo1 << std::endl;
			aTri[itri0].g2[1] = -2;
			aTri[itri0].s2[1] = aTri.size();
			aTri[itri0].r2[1] = noel2RelTriTri[ 0*3+2 ];
			////////////////
			aTriEd[ipo_min].first = itri0;
			aTriEd[ipo_min].second = 1;
			if( !Tesselate_GiftWrapping(ipo1,ipo_min,aTriEd,aVec,aTri) ) return false;
			break;
		}
	}
	return true;
}


bool Msh::DeleteTri(unsigned int itri_to, std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri)
{
	if( itri_to >= aTri.size() ) return true;
	{
		assert( aTri[itri_to].g2[0] == -1 );
		assert( aTri[itri_to].g2[1] == -1 );
		assert( aTri[itri_to].g2[2] == -1 );
	}
	const unsigned int itri_from = aTri.size()-1;
	if( itri_to == itri_from ){
		aTri.resize( aTri.size()-1 );
		return true;
	}
	aTri[itri_to] = aTri[itri_from];
	aTri.resize( aTri.size()-1 );
	for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
		if( aTri[itri_to].g2[iedtri] != -2 ) continue;
		const unsigned int itri_adj = aTri[itri_to].s2[iedtri];
		assert( itri_adj < aTri.size() );
    const unsigned int* rel = relTriTri[ (int)aTri[itri_to].r2[iedtri] ];
		const unsigned int iedtri_adj = rel[iedtri];
		assert( aTri[itri_adj].g2[iedtri_adj] == -2 );
		assert( aTri[itri_adj].s2[iedtri_adj] == itri_from );
		aTri[itri_adj].s2[iedtri_adj] = itri_to;
	}
	for(unsigned int inotri=0;inotri<nNoTri;inotri++){
		const unsigned int ipo0 = aTri[itri_to].v[inotri];
		aPo[ipo0].e = itri_to;
		aPo[ipo0].d = inotri;
	}
	return true;
}


bool Msh::DeletePointFromMesh(
        unsigned int ipo_del,
        std::vector<CPoint2D>& aPo,
        std::vector<STri2D>& aTri)
{
	std::vector<CVector2D> aVecAround;
    std::vector< std::pair<unsigned int,unsigned int> > aIndexTriAroundOld;
	std::vector<unsigned int> aIndexPoAround;
	{	// 節点配置を求め、境界にあるかどうか調べる
		const unsigned int itri_ini = aPo[ipo_del].e;
		const unsigned int inoel_c_ini = aPo[ipo_del].d;
		assert( itri_ini < aTri.size() );
		assert( inoel_c_ini < 3 );
		assert( aTri[itri_ini].v[inoel_c_ini] == ipo_del );
		unsigned int itri0= itri_ini;
		unsigned int inoel_c0 = inoel_c_ini;
		unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
		for(;;){
			assert( itri0 < aTri.size() );
			assert( inoel_c0 < 3 );
			assert( aTri[itri0].v[inoel_c0] == ipo_del );
			{
				aVecAround.push_back( aPo[ aTri[itri0].v[ noelTriEdge[inoel_c0][0] ] ].p );
				aIndexPoAround.push_back( aTri[itri0].v[ noelTriEdge[inoel_c0][0] ] );
				aIndexTriAroundOld.push_back( std::make_pair(itri0,inoel_c0) );
			}
			if( aTri[itri0].g2[inoel_b0] == -2 ){
				unsigned int itri1 = aTri[itri0].s2[inoel_b0];
        const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
        const unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
        const unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
				assert( itri1 < aTri.size() );
        assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
				assert( aTri[itri1].v[inoel_c1] == ipo_del );
				if( itri1 == itri_ini ) break;
				itri0 = itri1;
				inoel_c0 = inoel_c1;
				inoel_b0 = inoel_b1;
			}
			else{	// この点は境界上の点だから動かしてはならない。
				return false;
			}
		}
	}

	std::vector<STri2D> aTri_tmp;
	std::vector< std::pair<unsigned int,unsigned int> > aEdAround;
	{	// 点を取り除いた後の多角形の３角形分割を行う
		aTri_tmp.reserve( aVecAround.size()-2 );
		if( !Tesselate_GiftWrapping(0,1,aEdAround,aVecAround,aTri_tmp) ) return false;
		assert( aEdAround.size() == aVecAround.size() );
		aEdAround[0].first = 0;
		aEdAround[0].second = 2;
		aTri_tmp[0].g2[2] = -1;
		aTri_tmp[0].s2[2] = 0;
		for(unsigned int itri=0;itri<aTri_tmp.size();itri++){
		for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
			if( aTri_tmp[itri].g2[iedtri] == -1 ){	// 外側の３角形と接している
        const unsigned int ied1 = aTri_tmp[itri].s2[iedtri];
				assert( ied1 < aEdAround.size() );
				assert( aEdAround[ied1].first == itri );
        assert( aEdAround[ied1].second == iedtri );
			}
			else{	// 内側の３角形と接している．
				assert( aTri_tmp[itri].g2[iedtri] == -2 );
				const unsigned int itri1 = aTri_tmp[itri].s2[iedtri];
				assert( itri1 < aTri_tmp.size() );
        const unsigned int* rel0 = relTriTri[ aTri_tmp[itri].r2[iedtri] ];
        assert( aTri_tmp[itri1].g2[ rel0[iedtri] ] == -2 );
        assert( aTri_tmp[itri1].s2[ rel0[iedtri] ] == itri );
			}
		}
		}
	}

	const unsigned int ntri_new = aTri_tmp.size();
	{	// 新しい三角形をセット
		std::vector<STri2D> aTriAroundOld;
		{	// 古い３角形の情報を保存
			aTriAroundOld.resize( aIndexTriAroundOld.size() );
			for(unsigned int iitri=0;iitri<aIndexTriAroundOld.size();iitri++){
				unsigned int itri0 = aIndexTriAroundOld[iitri].first;
				aTriAroundOld[iitri] = aTri[itri0];
				unsigned int iedtri0 = aIndexTriAroundOld[iitri].second;
				assert( aTriAroundOld[iitri].v[iedtri0] == ipo_del );
			}
		}
		for(unsigned int iitri=0;iitri<ntri_new;iitri++){
			unsigned int itri0 = aIndexTriAroundOld[iitri].first;
      assert( aTriAroundOld[iitri].v[ (int)aIndexTriAroundOld[iitri].second ] == ipo_del );
			aTri[itri0].v[0] = aIndexPoAround[ aTri_tmp[iitri].v[0] ];
			aTri[itri0].v[1] = aIndexPoAround[ aTri_tmp[iitri].v[1] ];
			aTri[itri0].v[2] = aIndexPoAround[ aTri_tmp[iitri].v[2] ];
			for(unsigned int iedtri=0;iedtri<nEdTri;iedtri++){
				if( aTri_tmp[iitri].g2[iedtri] == -2 ){
					unsigned int jitri0 = aTri_tmp[iitri].s2[iedtri];
					assert( jitri0 < aTri_tmp.size() );
					unsigned int jtri0 = aIndexTriAroundOld[jitri0].first;
					aTri[itri0].g2[iedtri] = -2;
					aTri[itri0].s2[iedtri] = jtri0;
					aTri[itri0].r2[iedtri] = aTri_tmp[iitri].r2[iedtri];
				}
				else{
					assert( aTri_tmp[iitri].g2[iedtri] == -1 );
					unsigned int ied0 = aTri_tmp[iitri].s2[iedtri];
					assert( ied0 < aIndexTriAroundOld.size() );
					unsigned int iedtri_old = aIndexTriAroundOld[ied0].second;
					assert( aTriAroundOld[ied0].v[iedtri_old] == ipo_del );
					aTri[itri0].g2[iedtri] = aTriAroundOld[ied0].g2[iedtri_old];
					aTri[itri0].s2[iedtri] = aTriAroundOld[ied0].s2[iedtri_old];
					if( aTri[itri0].g2[iedtri] == -2 ){
						unsigned int itri_dia = aTriAroundOld[ied0].s2[iedtri_old];
						assert( itri_dia < aTri.size() );
            const unsigned int* rel0 = relTriTri[ aTriAroundOld[ied0].r2[iedtri_old] ];
						unsigned int iedtri_dia = rel0[iedtri_old];
						assert( aTri[itri_dia].g2[iedtri_dia] == -2 );
						assert( aTri[itri_dia].s2[iedtri_dia] == aIndexTriAroundOld[ied0].first );
						aTri[itri0].r2[iedtri] = ed2RelTriTri[ iedtri*3+iedtri_dia ];
						aTri[itri_dia].s2[iedtri_dia] = itri0;
            aTri[itri_dia].r2[iedtri_dia] = invRelTriTri[ aTri[itri0].r2[iedtri] ];
					}
				}
			}
		}
	}

	{	// 点を囲む要素一つを登録
		for(unsigned int iitri=0;iitri<ntri_new;iitri++){
			unsigned int itri0 = aIndexTriAroundOld[iitri].first;
			assert( itri0 < aTri.size() );
			for(unsigned int inotri=0;inotri<nNoTri;inotri++){
				unsigned int ipo0 = aTri[itri0].v[inotri];
				aPo[ipo0].e = itri0;
				aPo[ipo0].d = inotri;
			}
		}
	}
	
	{	// 古い三角形二つと頂点を孤立化
		aPo[ipo_del].e = -1;
		unsigned int itri0 = aIndexTriAroundOld[ntri_new  ].first;
		aTri[itri0].g2[0] = -1;  aTri[itri0].g2[1] = -1;  aTri[itri0].g2[2] = -1;
		unsigned int itri1 = aIndexTriAroundOld[ntri_new+1].first;
		aTri[itri1].g2[0] = -1;  aTri[itri1].g2[1] = -1;  aTri[itri1].g2[2] = -1;
		// ３角形を消去するときは番号が大きいものから消さないといけない．
		// 一番後ろにある三角形が、そのうち消す予定の３角形だったら、その３角形は消すことができないから．
		const unsigned int itri_1st = ( itri0 > itri1 ) ? itri0 : itri1;
		const unsigned int itri_2nd = ( itri0 < itri1 ) ? itri0 : itri1;
		DeleteTri(itri_1st,aPo,aTri);
		DeleteTri(itri_2nd,aPo,aTri);
	}

//	assert( CheckTri(aPo,aTri) );
	return true;
}

void Msh::PliantBossenHeckbertSmoothing( double elen, std::vector<CPoint2D>& aPo, std::vector<STri2D>& aTri )
{
	elen *= 0.75;

	for(unsigned int ipoin=0;ipoin<aPo.size();ipoin++)	// 点周りの点を探索して調べる。
	{
		if( aPo[ipoin].e == -1 ) continue;

		bool is_bound_flg = false;
        std::vector< std::pair<unsigned int,unsigned int> > aIndexTriAround;
		{	// 節点の周りの要素と向きを調べ、境界にあるかどうか調べる
			const unsigned int itri_ini = aPo[ipoin].e;
			const unsigned int inoel_c_ini = aPo[ipoin].d;
			assert( itri_ini < aTri.size() );
			assert( inoel_c_ini < 3 );
			assert( aTri[itri_ini].v[inoel_c_ini] == ipoin );
			unsigned int itri0= itri_ini;
			unsigned int inoel_c0 = inoel_c_ini;
			unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
			for(;;){
				assert( itri0 < aTri.size() );
				assert( inoel_c0 < 3 );
				assert( aTri[itri0].v[inoel_c0] == ipoin );
				{
					aIndexTriAround.push_back( std::make_pair(itri0,inoel_c0) );
				}
				if( aTri[itri0].g2[inoel_b0] == -2 ){
					unsigned int itri1 = aTri[itri0].s2[inoel_b0];
                    const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
                    unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
                    unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
					assert( itri1 < aTri.size() );
                    assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );
					assert( aTri[itri1].v[inoel_c1] == ipoin );
					if( itri1 == itri_ini ) break;
					itri0 = itri1;
					inoel_c0 = inoel_c1;
					inoel_b0 = inoel_b1;
				}
				else{	// この点は境界上の点だから動かしてはならない。
					is_bound_flg = true;
					break;
				}
			}
		}
		if( is_bound_flg ) continue;		

		/////////////////////////////////

		CVector2D vec_delta(0.0,0.0);
		{
			for(unsigned int iitri=0;iitri<aIndexTriAround.size();iitri++){
				unsigned int itri0 = aIndexTriAround[iitri].first;
                unsigned int inotri0 = noelTriEdge[ aIndexTriAround[iitri].second ][0];
				unsigned int ipo0 = aTri[itri0].v[inotri0];
				{
					double x_delta = aPo[ipo0].p.x - aPo[ipoin].p.x;
					double y_delta = aPo[ipo0].p.y - aPo[ipoin].p.y;
					const double dist = sqrt( x_delta*x_delta+y_delta*y_delta );
					x_delta /= dist; y_delta /= dist;
					const double dist_ratio = dist / elen;
					const double dtmp1 = dist_ratio*dist_ratio*dist_ratio*dist_ratio;
					const double dtmp2 = (1-dtmp1)*exp(-dtmp1);
					vec_delta.x += 0.2*dtmp2*x_delta*elen;
					vec_delta.y += 0.2*dtmp2*y_delta*elen;
				}
			}
		}

		bool is_moved = false;
		for(int i=0;i<4;i++)
		{
			CVector2D vec_new;
			vec_new.x = aPo[ipoin].p.x + vec_delta.x;
			vec_new.y = aPo[ipoin].p.y + vec_delta.y;

			bool is_ok_flg = true;
			for(unsigned int iitri=0;iitri<aIndexTriAround.size();iitri++){
				unsigned int itri0 = aIndexTriAround[iitri].first;
        unsigned int inotri0 = noelTriEdge[ aIndexTriAround[iitri].second ][0];
				unsigned int ipo0 = aTri[itri0].v[inotri0];
        unsigned int inotri1 = noelTriEdge[ aIndexTriAround[iitri].second ][1];
				unsigned int ipo1 = aTri[itri0].v[inotri1];
				{
					const double area = Com::TriArea(vec_new, aPo[ipo0].p, aPo[ipo1].p);
					if( area < 1.0e-10 ){ is_ok_flg = false; break; }
				}
			}
			if( is_ok_flg ){
				aPo[ipoin].p = vec_new;		// 点の移動
				DelaunayAroundPoint(ipoin,aPo,aTri);	// 周囲をDelaunay分割する
				is_moved = true;
				break;
			}
			else{
				vec_delta.x *= 0.5;
				vec_delta.y *= 0.5;
			}
		}
		if( !is_moved ) continue;

        double node_extent=0, min_edge_len=0;
		unsigned int itri_min, iedtri_min;
		{
			unsigned int ntri_around = 0;
			////////////////
			for(unsigned int iitri=0;iitri<aIndexTriAround.size();iitri++){
				unsigned int itri0 = aIndexTriAround[iitri].first;
				const double area0 = Com::TriArea(aPo[aTri[itri0].v[0]].p,aPo[aTri[itri0].v[1]].p,aPo[aTri[itri0].v[2]].p);
				node_extent += area0;
				double edge_len;
				{
					unsigned int inotri0 = noelTriEdge[ aIndexTriAround[iitri].second ][0];
					unsigned int ipo0 = aTri[itri0].v[inotri0];
					const double len = sqrt( (aPo[ipo0].p.x-aPo[ipoin].p.x)*(aPo[ipo0].p.x-aPo[ipoin].p.x)+(aPo[ipo0].p.y-aPo[ipoin].p.y)*(aPo[ipo0].p.y-aPo[ipoin].p.y) );
					edge_len = len / elen;
				}
				if( iitri==0 || edge_len < min_edge_len ){
					min_edge_len = edge_len;
					itri_min = itri0;
                    iedtri_min = noelTriEdge[ aIndexTriAround[iitri].second ][1];
				}
				ntri_around++;
			}
			node_extent = 0.8*2.3 * node_extent / ( elen*elen*sqrt(ntri_around*(ntri_around-2.0) ) );
		}

//		std::cout << ipoin << " " << node_extent << " " << min_edge_len << std::endl;

		if( min_edge_len < 0.1 || node_extent < 1.0 ){
			// Delete This Node
			std::vector<unsigned int> aIndexPoAround;
			{
				for(unsigned int iitri=0;iitri<aIndexTriAround.size();iitri++){
					unsigned int itri0 = aIndexTriAround[iitri].first;
                    unsigned int inotri0 = noelTriEdge[ aIndexTriAround[iitri].second ][0];
					unsigned int ipo0 = aTri[itri0].v[inotri0];
					aIndexPoAround.push_back(ipo0);
				}
			}
			DeletePointFromMesh(ipoin,aPo,aTri);
			for(unsigned int iipo=0;iipo<aIndexPoAround.size();iipo++){
				unsigned int ipoin0 = aIndexPoAround[iipo];
				if( ipoin0 >= aPo.size() ) continue;
				DelaunayAroundPoint(ipoin0,aPo,aTri);
			}
		}
		/*
		else if( max_edge_extent > 1.0 ){
			// Splid This Edge
			CVector2D po_ins;
			{
				unsigned int ipo1 = aTri[itri_max].v[ noelTriEdge[iedtri_max][0] ];
				assert( aTri[itri_max].v[ noelTriEdge[iedtri_max][1] ] == ipoin );
				po_ins.x = aPo[ipoin].p.x*0.5+aPo[ipo1].p.x*0.5;
				po_ins.y = aPo[ipoin].p.y*0.5+aPo[ipo1].p.y*0.5;
			}
			unsigned int ipo_ins = aPo.size();
			aPo.resize( aPo.size()+1 );
			aPo[ipo_ins].p = po_ins;
			aPo[ipo_ins].e = -1;
			InsertPoint_ElemEdge(ipo_ins, 
				itri_max, iedtri_max,
				aPo,aTri);
			DelaunayAroundPoint(ipoin,aPo,aTri);
		}
		*/
	}	// end loop ipoin
}


void Msh::ColorCodeBarAry( const std::vector<SBar>& aBar,  const std::vector<CVector2D>& aVec, 
						  std::vector< std::vector<unsigned int> >& aIndBarAry )
{
	const unsigned int nbar = aBar.size();

	std::vector<int> elsuel;
	{
		const unsigned int nvec = aVec.size();
		unsigned int* elsup_ind = new unsigned int [nvec+1];
		unsigned int nelsup;
		unsigned int* elsup;
		Msh::MakePointSurBar(aBar,nvec,elsup_ind,nelsup,elsup);
		elsuel.resize(nbar*2,-1);
		for(unsigned int ivec=0;ivec<nvec;ivec++){
			if( elsup_ind[ivec+1]-elsup_ind[ivec] != 2 ) continue;
			const unsigned int ielsup0 = elsup_ind[ivec];
			assert( ielsup0 < nelsup );
			const int ielem0 = elsup[ielsup0  ];	assert( ielem0>=0 && ielem0<(int)nbar );
			const int ielem1 = elsup[ielsup0+1];	assert( ielem1>=0 && ielem1<(int)nbar );
			unsigned int inobar0 = ( aBar[ielem0].v[0]==ivec ) ? 0 : 1;
			unsigned int inobar1 = ( aBar[ielem1].v[0]==ivec ) ? 0 : 1;
			assert( aBar[ielem0].v[inobar0] == ivec );
			assert( aBar[ielem1].v[inobar1] == ivec );
			elsuel[ielem0*2+(1-inobar0)] = ielem1;
			elsuel[ielem1*2+(1-inobar1)] = ielem0;
		}
		delete[] elsup_ind;
		delete[] elsup;
	}

	std::vector<int> aColor;
	aColor.clear();
	aColor.resize( nbar, -1 );

	for(;;){
		int ibar_ker = -1;
		for(unsigned int ibar=0;ibar<nbar;ibar++)	// まだ色付けされていない線要素ibar_kerを見つける
		{	
			if( aColor[ibar] == -1 ){
				ibar_ker = ibar;
				break;
			}
		}
		if( ibar_ker == -1 ) break;
		aColor[ibar_ker] = aIndBarAry.size();
		std::deque<unsigned int> aBarColor;
		aBarColor.push_back(ibar_ker);
		for(unsigned int inobar_ker=0;inobar_ker<2;inobar_ker++)	// ibar_kerから左右に探索する
		{
			unsigned int inobar0 = inobar_ker;
			unsigned int ibar0 = ibar_ker;
			for(;;){
				if( elsuel[ibar0*2+inobar0] == -1 ) break;
				unsigned int ibar1 = elsuel[ibar0*2+inobar0];
				if( aColor[ibar1] != -1 ) break;
				const unsigned int iv0 = aBar[ibar0].v[  inobar0];
				const unsigned int iv1 = aBar[ibar0].v[1-inobar0];
				unsigned int inobar1 = (aBar[ibar1].v[0]==iv1) ? 1 : 0;
				assert( aBar[ibar1].v[1-inobar1] == iv1 );
				assert( elsuel[ibar1*2+inobar1] == (int)ibar0 );
				const unsigned int iv2 = aBar[ibar1].v[inobar1];
				bool iflag = true;
				{
					CVector2D vec01( aVec[iv1].x-aVec[iv0].x,aVec[iv1].y-aVec[iv0].y);
					CVector2D vec12( aVec[iv2].x-aVec[iv1].x,aVec[iv2].y-aVec[iv1].y);
					{
						double len01 = sqrt( SquareLength(vec01) );
						double len12 = sqrt( SquareLength(vec12) );
						vec01.x /= len01; vec01.y /= len01;
						vec12.x /= len12; vec12.y /= len12;
					}
					if( Dot(vec01,vec12) > 0.8 ){ iflag = true; }
					else{ iflag = false; }
				}
				if( !iflag ) break;
				aColor[ibar1] = aIndBarAry.size();
				if( inobar_ker == 0 ){ aBarColor.push_back( ibar1); }
				else{                  aBarColor.push_front(ibar1); }
				ibar0 = ibar1;
				inobar0 = 1-inobar1;
			}
		}
		{
			std::vector<unsigned int> aBar_add;
			std::deque<unsigned int>::iterator itr = aBarColor.begin();
			for(;itr!=aBarColor.end();itr++){
				unsigned int ibar0 = *itr;
				aBar_add.push_back(ibar0);
			}
			aIndBarAry.push_back( aBar_add );
		}
	}
}

