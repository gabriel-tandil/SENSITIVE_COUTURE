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
// eqn_stokes.cpp : ストークス流体方程式の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <math.h>

#include "delfem/field_world.h"

#include "delfem/femls/linearsystem_field.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"

#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_hex.h"
#include "delfem/femeqn/eqn_stokes.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

void MakeMat_Stokes2D_Static_P1P1(
    double alpha, double g_x, double g_y, 
    double coords[][2],
    double emat_uu[][3][2][2], double emat_up[][3][2], double emat_pu[][3][2], double emat_pp[][3], 
    double eres_u[][2] )
{
    const unsigned int nno = 3;
    const unsigned int ndim = 2;

	// 要素剛性行列、残差を０で初期化
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&emat_up[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&emat_pu[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno;          i++){ *(&emat_pp[0][0]+i)       = 0.0; }

	// 面積を求める
	const double area = TriArea(coords[0],coords[1],coords[2]);

	// 形状関数のｘｙ微分を求める
	double dldx[nno][ndim];	// 形状関数のxy微分
	double const_term[nno];	// 形状関数の定数項
	TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		emat_uu[ino][jno][0][0] = dtmp1;
		emat_uu[ino][jno][1][1] = dtmp1;
	}
	}

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		emat_up[ino][jno][0] += area*dldx[ino][0]*0.333333333333333333333333;
		emat_up[ino][jno][1] += area*dldx[ino][1]*0.333333333333333333333333;
	}
	}

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		emat_pu[ino][jno][0] += area*dldx[jno][0]*0.333333333333333333333333;
		emat_pu[ino][jno][1] += area*dldx[jno][1]*0.333333333333333333333333;
	}
	}

	double tau;
	{
		const double h = sqrt( area / 3.14 )*2;
		tau = -h*h/alpha*0.1;
//		tau = 0.0;
	}

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		emat_pp[ino][jno] = dtmp1;
	}
	}

	// 外力ベクトルを求める
	for(unsigned int ino=0;ino<nno;ino++){
		eres_u[ino][0] = area*g_x*0.3333333333333333333;
		eres_u[ino][1] = area*g_y*0.3333333333333333333;
	}
}

static bool AddLinSys_Stokes2D_Static_P1P1(
		double alpha, double g_x, double g_y, 
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Stokes2D Static Tri P1P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_pres = world.GetField(id_field_press);

	const CElemAry::CElemSeg& es_velo_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_velo_c_va = field_velo.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_pres_c_va = field_pres.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double emat_uu[nno][nno][ndim][ndim];
	double emat_pp[nno][nno];
	double emat_pu[nno][nno][ndim];
	double emat_up[nno][nno][ndim];
	double eres_u[nno][ndim];
	double eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 != ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER,world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo, CORNER, id_field_press,CORNER, world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo, CORNER, world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,world);

	const CNodeAry::CNodeSeg& ns_co    = field_velo.GetNodeSeg( CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo  = field_velo.GetNodeSeg( CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_pres = field_pres.GetNodeSeg(CORNER,true, world,VELOCITY);//na_press.GetSeg(id_ns_press);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no_v[nno];	// 要素節点の全体節点番号		
		// 要素の節点番号を取ってくる
		es_velo_c_co.GetNodes(ielem,no_v);
		// 節点の座標、値を取ってくる
		double coords[nno][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(no_v[ino],coords[ino]);
		}

		// 要素の節点番号を取ってくる
		es_velo_c_va.GetNodes(ielem,no_v);
		// 節点の座標、値を取ってくる
		double velo[nno][ndim];	// 要素節点の値
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(no_v[ino],velo[ino]);
		}

		unsigned int no_p[nno];	// 要素節点の全体節点番号		
		// 要素の節点番号を取ってくる
		es_pres_c_va.GetNodes(ielem,no_p);
		double press[nno];
		for(unsigned int ino=0;ino<nno;ino++){
			ns_pres.GetValue(no_p[ino],&press[ino]);
		}

        MakeMat_Stokes2D_Static_P1P1(alpha, g_x,g_y, 
            coords,
            emat_uu,emat_up,emat_pu,emat_pp, eres_u);

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1];
				eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
				eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
			}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]
					+emat_pu[ino][jno][1]*velo[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= emat_pp[ino][jno]*press[jno];
			}
        }

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,no_v,nno,no_v,	4,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,no_v,nno,no_p,	2,&emat_up[0][0][0]);
		mat_pu.Mearge(nno,no_p,nno,no_v,	2,&emat_pu[0][0][0]);
		mat_pp.Mearge(nno,no_p,nno,no_p,	1,&emat_pp[0][0]);
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue(no_v[ino],0,eres_u[ino][0]);
			res_u.AddValue(no_v[ino],1,eres_u[ino][1]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue(no_p[ino],0,eres_p[ino]);
		}
	}
	return true;
}

static bool AddLinSys_Stokes2D_Static_P1P1_Combined(
		double alpha, double g_x, double g_y, 
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Stokes2D Static Tri P1P1_Combined" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg& es_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double emat_uu[nno][nno][ndim][ndim];
	double emat_pp[nno][nno];
	double emat_pu[nno][nno][ndim];
	double emat_up[nno][nno][ndim];
	double eres_u[nno][ndim];
	double eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		==  ls.FindIndexArray_Seg(id_field_press,CORNER,world) );
	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) != -1);
	assert( ls.FindIndexArray_Seg(id_field_press,CORNER,world) != -1);
	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(  id_field_velo,CORNER,world);
	CVector_Blk&    res_u  = ls.GetResidual(id_field_velo,CORNER,world);

	const CNodeAry::CNodeSeg& ns_co    = field_velo.GetNodeSeg( CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo  = field_velo.GetNodeSeg( CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_press = field_press.GetNodeSeg(CORNER,true, world,VELOCITY);//na_press.GetSeg(id_ns_press);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int noes[nno];	// 要素節点の全体節点番号		
		// 要素の節点番号を取ってくる
		es_c_co.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		double coords[nno][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
		}

		// 要素の節点番号を取ってくる
		es_c_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		double velo[nno][ndim];	// 要素節点の値
		double press[nno];
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
		}

        MakeMat_Stokes2D_Static_P1P1(alpha, g_x,g_y, 
            coords,
            emat_uu,emat_up,emat_pu,emat_pp, eres_u);

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][0] -= emat_uu[ino][jno][0][0]*velo[jno][0]+emat_uu[ino][jno][0][1]*velo[jno][1];
				eres_u[ino][1] -= emat_uu[ino][jno][1][0]*velo[jno][0]+emat_uu[ino][jno][1][1]*velo[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][0] -= emat_up[ino][jno][0]*press[jno];
				eres_u[ino][1] -= emat_up[ino][jno][1]*press[jno];
			}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= emat_pu[ino][jno][0]*velo[jno][0]
					+emat_pu[ino][jno][1]*velo[jno][1];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= emat_pp[ino][jno]*press[jno];
			}
		}

		double emat[nno][nno][3][3];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
			emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
			emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
			emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
			emat[ino][jno][0][2] = emat_up[ino][jno][0];
			emat[ino][jno][1][2] = emat_up[ino][jno][1];
			emat[ino][jno][2][0] = emat_pu[ino][jno][0];
			emat[ino][jno][2][1] = emat_pu[ino][jno][1];
			emat[ino][jno][2][2] = emat_pp[ino][jno];
		}
		}
		mat_uu.Mearge(nno,noes,nno,noes,	9,&emat[0][0][0][0]);
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
			res_u.AddValue( noes[ino],2,eres_p[ino]);
		}
	}

	return true;
}



static bool AddLinSys_Stokes2D_Static_P1bP1(
		double alpha, double g_x, double g_y, 
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea)
{
//	std::cout << "Stokes2D Static Tri P1bP1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const unsigned int nno_c = 3;
//	const unsigned int nno_b = 1;
	const unsigned int ndim = 2;

	CMatDia_BlkCrs& mat_cucu = ls.GetMatrix(id_field_velo, CORNER,  world);
	CMatDia_BlkCrs& mat_bubu = ls.GetMatrix(id_field_velo, BUBBLE,  world);
	CMatDia_BlkCrs& mat_pp   = ls.GetMatrix(id_field_press,CORNER,  world);

	CMat_BlkCrs& mat_cubu = ls.GetMatrix(id_field_velo, CORNER,  id_field_velo, BUBBLE,  world);
	CMat_BlkCrs& mat_bucu = ls.GetMatrix(id_field_velo, BUBBLE,  id_field_velo, CORNER,  world);
	CMat_BlkCrs& mat_cup  = ls.GetMatrix(id_field_velo, CORNER,  id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_bup  = ls.GetMatrix(id_field_velo, BUBBLE,  id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_pcu  = ls.GetMatrix(id_field_press,CORNER,  id_field_velo, CORNER,  world);
	CMat_BlkCrs& mat_pbu  = ls.GetMatrix(id_field_press,CORNER,  id_field_velo, BUBBLE,  world);

	CVector_Blk& res_cu = ls.GetResidual(id_field_velo, CORNER,  world);
	CVector_Blk& res_bu = ls.GetResidual(id_field_velo, BUBBLE,  world);
	CVector_Blk& res_p  = ls.GetResidual(id_field_press,CORNER,  world); 

	const CElemAry::CElemSeg& es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_b_va = field_velo.GetElemSeg(id_ea,BUBBLE,true, world);
	const CElemAry::CElemSeg& es_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);

	const CNodeAry::CNodeSeg& ns_co      = field_velo.GetNodeSeg( CORNER,false,world,VALUE   );
	const CNodeAry::CNodeSeg& ns_press   = field_press.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_velo  = field_velo.GetNodeSeg( CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_b_velo  = field_velo.GetNodeSeg( BUBBLE,true, world,VELOCITY);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
	    unsigned int noes_c[nno_c];	// 要素節点の全体節点番号
		es_c_va.GetNodes(ielem,noes_c);
		// 節点の座標、値を取ってくる
	    double velo_c[nno_c][ndim];	// 要素節点の値
	    double press[nno_c];
		for(unsigned int ino=0;ino<nno_c;ino++){
			ns_c_velo.GetValue(noes_c[ino],velo_c[ino]);
			ns_press.GetValue(noes_c[ino],&press[ino]);
		}
		// 要素の節点番号を取ってくる
		es_c_co.GetNodes(ielem,noes_c);
	    double coords[nno_c][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno_c;ino++){
			ns_co.GetValue(noes_c[ino],coords[ino]);
		}
	    unsigned int noes_b;	// 要素節点の全体節点番号
		es_b_va.GetNodes(ielem,&noes_b);
	    double velo_b[ndim];	// 要素節点の値
		ns_b_velo.GetValue(noes_b,velo_b);

	    double emat_cucu[nno_c][nno_c][ndim][ndim];
	    double emat_cubu[nno_c][ndim][ndim];
	    double emat_cup[nno_c][nno_c][ndim];

	    double emat_bubu[ndim][ndim];
	    double emat_bucu[nno_c][ndim][ndim];
	    double emat_bup[nno_c][ndim];

	    double emat_pcu[nno_c][nno_c][ndim];
	    double emat_pbu[nno_c][ndim];
	    double emat_pp[nno_c][nno_c];

	    double eres_cu[nno_c][ndim], eres_bu[ndim], eres_p[nno_c];

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno_c*nno_c*ndim*ndim; i++){ *(&emat_cucu[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno_c*      ndim*ndim; i++){ *(&emat_cubu[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno_c*      ndim*ndim; i++){ *(&emat_bucu[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<            ndim*ndim; i++){ *(&emat_bubu[0][0]+i)       = 0.0; }
		for(unsigned int i=0;i<nno_c*nno_c*ndim; i++){ *(&emat_cup[0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno_c*      ndim; i++){ *(&emat_bup[0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno_c*nno_c*ndim; i++){ *(&emat_pcu[0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno_c*      ndim; i++){ *(&emat_pbu[0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno_c*nno_c;      i++){ *(&emat_pp[0][0]+i)     = 0.0; }

		// 面積を求める
		const double area = TriArea(coords[0],coords[1],coords[2]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno_c][ndim];	// 形状関数のxy微分
		double const_term[nno_c];	// 形状関数の定数項
		TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

		double vc_b[4];
		vc_b[0] = 1.0/3.0; vc_b[1] = 1.0/3.0; vc_b[2] = 1.0/3.0; vc_b[3] = 27.0;
		{
			const double tmp_val1 = area*vc_b[3]*vc_b[3]/180.0*( 
				dldx[0][0]*dldx[0][0]+dldx[0][1]*dldx[0][1]+
				dldx[1][0]*dldx[1][0]+dldx[1][1]*dldx[1][1]+
				dldx[2][0]*dldx[2][0]+dldx[2][1]*dldx[2][1] );
			for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				const double tmp1 = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1])
					 + alpha*vc_b[ino]*vc_b[jno]*tmp_val1;
				emat_cucu[ino][jno][0][0] = tmp1;
				emat_cucu[ino][jno][1][1] = tmp1;
			}
			}
			for(unsigned int ino=0;ino<nno_c;ino++){
				const double tmp1 = -alpha*vc_b[ino]*tmp_val1;
				emat_cubu[ino][0][0] = tmp1;
				emat_cubu[ino][1][1] = tmp1;
				emat_bucu[ino][0][0] = tmp1;
				emat_bucu[ino][1][1] = tmp1;
			}
			emat_bubu[0][0] = alpha*tmp_val1;
			emat_bubu[1][1] = alpha*tmp_val1;
		}

		{	// Grad_P1B
			const double tmp_val1 = area/3.0;
			const double tmp_val2 = area*vc_b[3]/60.0;
			for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				emat_cup[ino][jno][0] += dldx[ino][0]*tmp_val1 + dldx[jno][0]*vc_b[ino]*tmp_val2;
				emat_cup[ino][jno][1] += dldx[ino][1]*tmp_val1 + dldx[jno][1]*vc_b[ino]*tmp_val2;
			}
			}
			for(unsigned int jno=0;jno<nno_c;jno++){
				emat_bup[jno][0] += -dldx[jno][0]*tmp_val2;
				emat_bup[jno][1] += -dldx[jno][1]*tmp_val2;
			}
		}
		{	// Div_P1B
			const double tmp_val1 = area/3.0;
			const double tmp_val2 = area*vc_b[3]/60.0;
			for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				emat_pcu[ino][jno][0] += dldx[jno][0]*tmp_val1 + dldx[ino][0]*vc_b[jno]*tmp_val2;
				emat_pcu[ino][jno][1] += dldx[jno][1]*tmp_val1 + dldx[ino][1]*vc_b[jno]*tmp_val2;
			}
			}
            for(unsigned int ino=0;ino<nno_c;ino++){
				emat_pbu[ino][0] += -dldx[ino][0]*tmp_val2;
				emat_pbu[ino][1] += -dldx[ino][1]*tmp_val2;
			}
		}

		////////////////
		for(unsigned int ino=0;ino<nno_c;ino++)
		{
			eres_cu[ino][0] = g_x*area*11.0/60.0;
			eres_cu[ino][1] = g_y*area*11.0/60.0;
		}
		eres_bu[0] = g_x*area*27.0/60.0;
		eres_bu[1] = g_y*area*27.0/60.0;
		for(unsigned int ino=0;ino<nno_c;ino++){ eres_p[ino] = 0.0; }
		////////////////
		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				eres_cu[ino][0] -= emat_cucu[ino][jno][0][0]*velo_c[jno][0]+emat_cucu[ino][jno][0][1]*velo_c[jno][1];
				eres_cu[ino][1] -= emat_cucu[ino][jno][1][0]*velo_c[jno][0]+emat_cucu[ino][jno][1][1]*velo_c[jno][1];
			}
			eres_cu[ino][0] -= emat_cubu[ino][0][0]*velo_b[0]+emat_cubu[ino][0][1]*velo_b[1];
			eres_cu[ino][1] -= emat_cubu[ino][1][0]*velo_b[0]+emat_cubu[ino][1][1]*velo_b[1];
			for(unsigned int jno=0;jno<nno_c;jno++){
				eres_cu[ino][0] -= emat_cup[ino][jno][0]*press[jno];
				eres_cu[ino][1] -= emat_cup[ino][jno][1]*press[jno];
			}
		}
		////////
		for(unsigned int jno=0;jno<nno_c;jno++){
			eres_bu[0] -= emat_bucu[jno][0][0]*velo_c[jno][0];
			eres_bu[1] -= emat_bucu[jno][1][1]*velo_c[jno][1];
		}
		eres_bu[0] -= emat_bubu[0][0]*velo_b[0];
		eres_bu[1] -= emat_bubu[1][1]*velo_b[1];
		for(unsigned int jno=0;jno<nno_c;jno++){
			eres_bu[0] -= emat_bup[jno][0]*press[jno];
			eres_bu[1] -= emat_bup[jno][1]*press[jno];
		}
		////////
		for(unsigned int ino=0;ino<nno_c;ino++){
			for(unsigned int jno=0;jno<nno_c;jno++){
				eres_p[ino] -= emat_pcu[ino][jno][0]*velo_c[jno][0]
					         + emat_pcu[ino][jno][1]*velo_c[jno][1];
			}
			eres_p[ino] -= emat_pbu[ino][0]*velo_b[0]+
				           emat_pbu[ino][1]*velo_b[1];
			for(unsigned int jno=0;jno<nno_c;jno++){
				eres_p[ino] -= emat_pp[ino][jno]*press[jno];
			}
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_cucu.Mearge(nno_c, noes_c,  nno_c, noes_c,	 4,&emat_cucu[0][0][0][0]);
		mat_cubu.Mearge(nno_c, noes_c,  1,    &noes_b, 	 4,&emat_cubu[0][0][0]);
		mat_bucu.Mearge(1,    &noes_b,  nno_c, noes_c,	 4,&emat_bucu[0][0][0]);
		mat_bubu.Mearge(1,    &noes_b,  1,    &noes_b,	 4,&emat_bubu[0][0]);
		mat_cup .Mearge(nno_c, noes_c,  nno_c, noes_c,	 2,&emat_cup[0][0][0]);
		mat_pcu .Mearge(nno_c, noes_c,  nno_c, noes_c,	 2,&emat_pcu[0][0][0]);
		mat_bup .Mearge(1,    &noes_b,  nno_c, noes_c,	 2,&emat_bup[0][0]);
		mat_pbu .Mearge(nno_c, noes_c,  1,    &noes_b,	 2,&emat_pbu[0][0]);
		mat_pp  .Mearge(nno_c, noes_c,  nno_c, noes_c,	 1,&emat_pp[0][0]);

		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno_c;ino++){
			res_cu.AddValue( noes_c[ino],0,eres_cu[ino][0]);
			res_cu.AddValue( noes_c[ino],1,eres_cu[ino][1]);
		}
		res_bu.AddValue( noes_b,0,eres_bu[0]);
		res_bu.AddValue( noes_b,1,eres_bu[1]);
		for(unsigned int ino=0;ino<nno_c;ino++){
			res_p.AddValue( noes_c[ino],0,eres_p[ino]);
		}
	}
	return true;
}



bool Fem::Eqn::AddLinSys_Stokes2D_Static(
		double alpha, 
		double g_x, double g_y,
		CLinearSystem_Field& ls,
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea ) 
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;

	const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++)
	{
		const unsigned int id_ea = aIdEA[iiea];
		if(      field_velo.GetInterpolationType(id_ea,world) == TRI11   ){
            if(    ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_press,CORNER,world) ){
                AddLinSys_Stokes2D_Static_P1P1_Combined(
				alpha,   g_x,g_y,   
				ls,   id_field_velo,id_field_press,
				world,id_ea);
            }
            else{
                AddLinSys_Stokes2D_Static_P1P1(
				alpha,   g_x,g_y,   
				ls,   id_field_velo,id_field_press,
				world,id_ea);
            }
		}
		else if( field_velo.GetInterpolationType(id_ea,world) == TRI1011 ){
			AddLinSys_Stokes2D_Static_P1bP1(
				alpha,   g_x,g_y,
				ls,   id_field_velo,id_field_press,
				world,id_ea);
		}
		else{
			std::cout << "Interpolation not defined for stokes " << field_velo.GetInterpolationType(id_ea,world) << std::endl;
			assert(0);
		}
	}
	return true;
}

void AddMat_Stokes2D_NonStatic_Newmark_P1P1( 
    double alpha, double rho, double g_x, double g_y, 
    double coords[3][2],
    double eCmat_uu[][3][2][2], double eCmat_up[][3][2], double eCmat_pu[][3][2], double eCmat_pp[][3],
    double eMmat_uu[][3][2][2],
    double eqf_out_u[3][2])
{
//	std::cout << "AddMat_Stokes2D_NonStatic_Newmark_P1P1" << std::endl;

    const unsigned int nno = 3;
    const unsigned int ndim = 2;

	// 要素剛性行列、残差を０で初期化
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eCmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_up[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_pu[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno;          i++){ *(&eCmat_pp[0][0]+i)       = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*ndim;         i++){ *(&eqf_out_u[0][0]+i)      = 0.0; }

	// 面積を求める
	const double area = TriArea(coords[0],coords[1],coords[2]);

	// 形状関数のｘｙ微分を求める
	double dldx[nno][ndim];	// 形状関数のxy微分
	double const_term[nno];	// 形状関数の定数項
	TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = area*alpha*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		eCmat_uu[ino][jno][0][0] = dtmp1;
		eCmat_uu[ino][jno][1][1] = dtmp1;
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eCmat_up[ino][jno][0] += area*dldx[ino][0]*0.333333333333333333333;
		eCmat_up[ino][jno][1] += area*dldx[ino][1]*0.333333333333333333333;
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eCmat_pu[ino][jno][0] += area*dldx[jno][0]*0.333333333333333333333;
		eCmat_pu[ino][jno][1] += area*dldx[jno][1]*0.333333333333333333333;
	}
    }
	
	double tau;
	{
		const double h = sqrt( area / 3.14 )*2;
		tau = -h*h/alpha*0.1;
	}

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = area*tau*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		eCmat_pp[ino][jno] = dtmp1;
	}
	}

    ////////////////
	{
		const double dtmp1 = area*rho*0.0833333333333333333333333333;
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eMmat_uu[ino][jno][0][0] = dtmp1;
				eMmat_uu[ino][jno][1][1] = dtmp1;
			}
			eMmat_uu[ino][ino][0][0] += dtmp1;
			eMmat_uu[ino][ino][1][1] += dtmp1;
		}
    }

	// 外力ベクトルを求める
	for(unsigned int ino=0;ino<nno;ino++){
		eqf_out_u[ino][0] = area*g_x*0.33333333333333333333;
		eqf_out_u[ino][1] = area*g_y*0.33333333333333333333;
	}
}

static bool AddLinearSystem_Stokes2D_NonStatic_Newmark_P1P1(
		double rho, double alpha, 
		double g_x, double g_y,
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea )
{
	std::cout << "Storkes2D NonStatic_Newmark Tri P1P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg& es_velo_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_velo_c_va = field_velo.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_pres_c_va = field_press.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim];
	double emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eqf_out_u[nno][ndim], eres_u[nno][ndim];
	double eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 != ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo,CORNER,world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo,CORNER, id_field_press,CORNER, world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo,CORNER, world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no_v[nno];	// 要素節点の全体節点番号
		// 要素の節点番号を取ってくる
		es_velo_c_co.GetNodes(ielem,no_v);	
		double coords[nno][ndim];	// 要素節点の座標
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(no_v[ino],coords[ino]);
		}

		// 要素の節点番号を取ってくる
		es_velo_c_va.GetNodes(ielem,no_v);
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(no_v[ino],velo[ino]);
			ns_acc.GetValue(no_v[ino],acc[ino]);
		}
		
		unsigned int no_p[nno];	// 要素節点の全体節点番号
		// 要素の節点番号を取ってくる
		es_pres_c_va.GetNodes(ielem,no_p);
		double press[nno];
		double apress[nno];
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_press.GetValue(no_p[ino],&press[ino]);
			ns_apress.GetValue(no_p[ino],&apress[ino]);
		}

        AddMat_Stokes2D_NonStatic_Newmark_P1P1(alpha,rho, g_x,g_y,
            coords,eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp,
            eMmat_uu,eqf_out_u);

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres_u[ino][idim] = eqf_out_u[ino][idim];
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt*acc[jno][jdim])
					+ eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
			}
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt*apress[jno]);
			}
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt*acc[jno][jdim]);
			}
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
			}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					emat_uu[ino][jno][idim][jdim] 
						= eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
				}
				}
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
					emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
					emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
				}
			}
			}		
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
			}
			}
		}
        {
			// 要素剛性行列の全体剛性行列へのマージ
			mat_uu.Mearge(nno,no_v,nno,no_v,	4,&emat_uu[0][0][0][0]);
			mat_up.Mearge(nno,no_v,nno,no_p,	2,&emat_up[0][0][0]);
			mat_pu.Mearge(nno,no_p,nno,no_v,	2,&emat_pu[0][0][0]);
			mat_pp.Mearge(nno,no_p,nno,no_p,	1,&emat_pp[0][0]);
			// 残差ベクトルのマージ
			for(unsigned int ino=0;ino<nno;ino++){
				res_u.AddValue( no_v[ino],0,eres_u[ino][0]);
				res_u.AddValue( no_v[ino],1,eres_u[ino][1]);
			}
			for(unsigned int ino=0;ino<nno;ino++){
				res_p.AddValue( no_p[ino],0,eres_p[ino]);
			}
		}
	}
	return true;
}


static bool AddLinearSystem_Stokes2D_NonStatic_Newmark_P1P1_Combined(
		double rho, double alpha, 
		double g_x, double g_y,
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Storkes2D NonStatic_Newmark Tri P1P1 Combined" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg& es_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim];
	double emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eqf_out_u[nno][ndim], eres_u[nno][ndim];
	double eres_p[nno];
	
	assert( field_velo.GetIdElemSeg(id_ea,CORNER,true,world) 
		 == field_press.GetIdElemSeg(id_ea,CORNER,true,world) );

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix( id_field_velo, CORNER,world);
	CVector_Blk& res_u = ls.GetResidual( id_field_velo, CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int noes[nno];	// 要素節点の全体節点番号			
		es_c_co.GetNodes(ielem,noes);	// 要素の節点番号を取ってくる
		double coords[nno][ndim];	// 要素節点の座標
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
		}

		es_c_va.GetNodes(ielem,noes);// 要素の節点番号を取ってくる
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		double press[nno];
		double apress[nno];
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_acc.GetValue(noes[ino],acc[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
			ns_apress.GetValue(noes[ino],&apress[ino]);
		}

        AddMat_Stokes2D_NonStatic_Newmark_P1P1(alpha,rho, g_x,g_y,
            coords,eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp,
            eMmat_uu,eqf_out_u);

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres_u[ino][idim] = eqf_out_u[ino][idim];
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_u[ino][idim] -= eCmat_uu[ino][jno][idim][jdim]*(velo[jno][jdim]+dt*acc[jno][jdim])
					+ eMmat_uu[ino][jno][idim][jdim]*acc[jno][jdim];
			}
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_u[ino][idim] -= eCmat_up[ino][jno][idim]*(press[jno]+dt*apress[jno]);
			}
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_p[ino] -= eCmat_pu[ino][jno][jdim]*(velo[jno][jdim]+dt*acc[jno][jdim]);
			}
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
			}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					emat_uu[ino][jno][idim][jdim] 
						= eMmat_uu[ino][jno][idim][jdim]+dtmp1*eCmat_uu[ino][jno][idim][jdim];
				}
				}
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
					emat_up[ino][jno][idim] = dtmp1*eCmat_up[ino][jno][idim];
					emat_pu[ino][jno][idim] = dtmp1*eCmat_pu[ino][jno][idim];
				}
			}
			}		
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				emat_pp[ino][jno] = dtmp1*eCmat_pp[ino][jno];
			}
			}
		}

		double emat[nno][nno][3][3];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno][0][0] = emat_uu[ino][jno][0][0];
			emat[ino][jno][0][1] = emat_uu[ino][jno][0][1];
			emat[ino][jno][1][0] = emat_uu[ino][jno][1][0];
			emat[ino][jno][1][1] = emat_uu[ino][jno][1][1];
			emat[ino][jno][0][2] = emat_up[ino][jno][0];
			emat[ino][jno][1][2] = emat_up[ino][jno][1];
			emat[ino][jno][2][0] = emat_pu[ino][jno][0];
			emat[ino][jno][2][1] = emat_pu[ino][jno][1];
			emat[ino][jno][2][2] = emat_pp[ino][jno];
		}
		}
		mat_uu.Mearge(nno,noes,nno,noes,	9,&emat[0][0][0][0]);
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
			res_u.AddValue( noes[ino],2,eres_p[ino]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_Stokes2D_NonStatic_Newmark(
		double rho, double alpha, 
		double g_x, double g_y,
		double gamma, double dt,
		CLinearSystem_Field& ls,
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;

	if( id_ea != 0 ){   // 特定のＥＡにマージ
		bool res;   
		if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_velo,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_press,CORNER,world) ){
			    res = AddLinearSystem_Stokes2D_NonStatic_Newmark_P1P1_Combined(
				    rho,alpha,g_x,g_y,
				    gamma,dt,
				    ls,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
            else{
			    res = AddLinearSystem_Stokes2D_NonStatic_Newmark_P1P1(
				    rho,alpha,g_x,g_y,
				    gamma,dt,
				    ls,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
		}
		else{
			res = false;
			assert(0);
		}
		return res;
	}
	else{   // fieldに属する全てのＥＡにマージ
		const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res;
			res = AddLinSys_Stokes2D_NonStatic_Newmark(
				rho,alpha,g_x,g_y,
				gamma,dt,
				ls,
				id_field_velo,id_field_press,world,
				id_ea);
			return res;
		}
	}
	return true;
}









static bool AddLinearSystem_Stokes3D_Static_P1P1(
		double alpha, double g_x, double g_y, double g_z,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Stokes3D Tetrahedra 4-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg& es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 4;
	const unsigned int ndim = 3;

	unsigned int noes[nno];	// 要素節点の全体節点番号

	double velo[nno][ndim];	// 要素節点の値
	double press[nno];
	double coords[nno][ndim];	// 要素節点の座標
				
	double dldx[nno][ndim];	// 形状関数のxy微分
	double const_term[nno];	// 形状関数の定数項

	double emat_uu[nno][nno][ndim][ndim];
	double emat_pp[nno][nno];
	double emat_pu[nno][nno][ndim];
	double emat_up[nno][nno][ndim];
	double eqf_out_u[nno][ndim], eres_u[nno][ndim];
	double eres_p[nno];

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER,  world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo,CORNER,    id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER,   id_field_velo, CORNER,  world);	
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,  world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,  world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);

	assert( ns_co.Length() == ndim );
	assert( ns_velo.Length() == ndim );
	assert( ns_press.Length() == 1 );

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_c_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
		}

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&emat_up[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&emat_pu[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno*nno;          i++){ *(&emat_pp[0][0]+i)       = 0.0; }
		for(unsigned int i=0;i<nno*ndim;         i++){ *(&eqf_out_u[0][0]+i)     = 0.0; }

		// 面積を求める
		double vol = TetVolume(coords[0],coords[1],coords[2],coords[3]);

		// 形状関数のｘｙ微分を求める
		TetDlDx(dldx, const_term,   coords[0],coords[1],coords[2],coords[3]);

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			const double dtmp1 = vol*alpha*(
				 dldx[jno][0]*dldx[ino][0]
				+dldx[jno][1]*dldx[ino][1]
				+dldx[jno][2]*dldx[ino][2]);
			emat_uu[ino][jno][0][0] = dtmp1;
			emat_uu[ino][jno][1][1] = dtmp1;
			emat_uu[ino][jno][2][2] = dtmp1;
		}
		}

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat_up[ino][jno][0] += vol*dldx[ino][0]/4.0;
			emat_up[ino][jno][1] += vol*dldx[ino][1]/4.0;
			emat_up[ino][jno][2] += vol*dldx[ino][2]/4.0;
		}
		}

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat_pu[ino][jno][0] += vol*dldx[jno][0]/4.0;
			emat_pu[ino][jno][1] += vol*dldx[jno][1]/4.0;
			emat_pu[ino][jno][2] += vol*dldx[jno][2]/4.0;
		}
		}

		double tau;
		{
//			const double h = pow( vol / 3.14, 1/3 )*2;
//			tau = -h*h/alpha*0.1;
			tau = 0.001;
		}

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			const double dtmp1 = vol*tau*(
				 dldx[jno][0]*dldx[ino][0]
				+dldx[jno][1]*dldx[ino][1]
				+dldx[jno][2]*dldx[ino][2]);
			emat_pp[ino][jno] = -dtmp1;
		}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out_u[ino][0] = vol*g_x/4.0;
			eqf_out_u[ino][1] = vol*g_y/4.0;
			eqf_out_u[ino][2] = vol*g_z/4.0;
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_u[ino][0] = eqf_out_u[ino][0];
			eres_u[ino][1] = eqf_out_u[ino][1];
			eres_u[ino][2] = eqf_out_u[ino][2];
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres_u[ino][idim] -= 
					 emat_uu[ino][jno][idim][0]*velo[jno][0]
					+emat_uu[ino][jno][idim][1]*velo[jno][1]
					+emat_uu[ino][jno][idim][2]*velo[jno][2];
			}
			}
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres_u[ino][idim] -= emat_up[ino][jno][idim]*press[jno];
			}
			}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= 
					 emat_pu[ino][jno][0]*velo[jno][0]
					+emat_pu[ino][jno][1]*velo[jno][1]
					+emat_pu[ino][jno][2]*velo[jno][2];
			}
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p[ino] -= emat_pp[ino][jno]*press[jno];
			}
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,noes,nno,noes, ndim*ndim,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,noes,nno,noes,	   ndim,&emat_up[0][0][0]);
		mat_pu.Mearge(nno,noes,nno,noes,	   ndim,&emat_pu[0][0][0]);
		mat_pp.Mearge(nno,noes,nno,noes,	      1,&emat_pp[0][0]);
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
			res_u.AddValue( noes[ino],2,eres_u[ino][2]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue( noes[ino],0,eres_p[ino]);
		}
	}
	return true;
}


static bool AddLinearSystem_Stokes3D_Static_Hex81(
		double alpha, 
		double rho, double g_x, double g_y, double g_z,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea)
{
//	std::cout << "Poisson3D Hexahedra 8-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == HEX );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

//	const Fem::Field::CField::CElemInterpolation& ei_v = field_velo.GetElemInterpAry()[iei];
//	const Fem::Field::CField::CElemInterpolation& ei_p = field_press.GetElemInterpAry()[iei];
//	assert( ei_v.id_ea == ei_p.id_ea );
//	const unsigned int id_ea = ei_v.id_ea;

	const CElemAry::CElemSeg& es_c_v = field_velo.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg& es_b_p = field_press.GetElemSeg(id_ea,BUBBLE,true,world);
	assert( es_c_v.Length() == 8 );
	assert( es_b_p.Length() == 1 );

	unsigned int num_integral = 1;
	const unsigned int nInt = NIntLineGauss[num_integral];
	const double (*Gauss)[2] = LineGauss[num_integral];
	double detjac, detwei;

	const unsigned int nno = 8;
	const unsigned int ndim = 3;

	unsigned int noes[nno];	// 要素節点の全体節点番号
	unsigned int noes_b;

	double velo[nno][ndim];	// 要素節点の値
	double press;
	double coords[nno][ndim];	// 要素節点の座標
				
	double dndx[nno][ndim];	// 形状関数のxy微分
	double an[nno];	// 形状関数の定数項

	double emat_uu[nno][nno][ndim][ndim];
	double emat_pp;
	double emat_pu[nno][ndim];
	double emat_up[nno][ndim];
	double eqf_out_u[nno][ndim], eres_u[nno][ndim];
	double eres_p;

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER,  world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,BUBBLE,  world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo, CORNER,  id_field_press,BUBBLE,  world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,BUBBLE,  id_field_velo, CORNER,  world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,  world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,BUBBLE,  world);

	const CNodeAry::CNodeSeg& ns_co    = field_velo.GetNodeSeg(CORNER,false,world);
	const CNodeAry::CNodeSeg& ns_velo  = field_velo.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_press = field_press.GetNodeSeg(BUBBLE,true,world,VELOCITY);

	assert( ns_co.Length()    == ndim );
	assert( ns_velo.Length()  == ndim );
	assert( ns_press.Length() == 1 );

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_c_v.GetNodes(ielem,noes);
		es_b_p.GetNodes(ielem,&noes_b);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(   noes[ino], coords[ino]);
			ns_velo.GetValue( noes[ino], velo[ino]);
			ns_press.GetValue(noes_b,    &press);
		}

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&emat_uu[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*ndim;         i++){ *(&emat_up[0][0]+i)       = 0.0; }
		for(unsigned int i=0;i<nno*ndim;         i++){ *(&emat_pu[0][0]+i)       = 0.0; }
		emat_pp = 0;
		for(unsigned int i=0;i<nno*ndim;         i++){ *(&eqf_out_u[0][0]+i)     = 0.0; }

		double vol = 0.0;
		for(unsigned int ir1=0;ir1<nInt;ir1++){
		for(unsigned int ir2=0;ir2<nInt;ir2++){
		for(unsigned int ir3=0;ir3<nInt;ir3++){
			const double r1 = Gauss[ir1][0];
			const double r2 = Gauss[ir2][0];
			const double r3 = Gauss[ir3][0];
			ShapeFunc_Hex8(r1,r2,r3,coords,detjac,dndx,an);
			detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
			vol += detwei;
            for(unsigned int ino=0;ino<nno;ino++){
            for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
					dtmp1 += dndx[ino][idim]*dndx[jno][idim];
				}
				for(unsigned int idim=0;idim<ndim;idim++){
					emat_uu[ino][jno][idim][idim] += detwei*alpha*dtmp1;
				}
/*				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					emat_uu[ino][jno][idim][jdim] += detwei*alpha*dndx[ino][jdim]*dndx[jno][idim];
				}
				}*/
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
				emat_up[ino][0] += detwei*dndx[ino][0];
				emat_up[ino][1] += detwei*dndx[ino][1];
				emat_up[ino][2] += detwei*dndx[ino][2];
			}
			for(unsigned int ino=0;ino<nno;ino++){
				emat_pu[ino][0] += detwei*dndx[ino][0];
				emat_pu[ino][1] += detwei*dndx[ino][1];
				emat_pu[ino][2] += detwei*dndx[ino][2];
			}
			emat_pp = -0.001;
			// 要素節点等価外力ベクトルを積n分毎に足し合わせる
			for(unsigned int ino=0;ino<nno;ino++){
				eqf_out_u[ino][0] += detwei*rho*g_x*an[ino];
				eqf_out_u[ino][1] += detwei*rho*g_y*an[ino];
				eqf_out_u[ino][2] += detwei*rho*g_z*an[ino];
			}
		}
		}
		}

		////////////////////////////////////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_u[ino][0] = eqf_out_u[ino][0];
			eres_u[ino][1] = eqf_out_u[ino][1];
			eres_u[ino][2] = eqf_out_u[ino][2];
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eres_u[ino][idim] -= 
					 emat_uu[ino][jno][idim][0]*velo[jno][0]
					+emat_uu[ino][jno][idim][1]*velo[jno][1]
					+emat_uu[ino][jno][idim][2]*velo[jno][2];
			}
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				eres_u[ino][idim] -= emat_up[ino][idim]*press;
			}
		}

		// 要素内残差ベクトルを求める
		{
			eres_p = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_p -= 
					 emat_pu[jno][0]*velo[jno][0]
					+emat_pu[jno][1]*velo[jno][1]
					+emat_pu[jno][2]*velo[jno][2];
			}
			eres_p -= emat_pp*press;
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,noes,     nno,noes,     ndim*ndim, &emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,noes,     1,  &noes_b,  ndim,      &emat_up[0][0]);
		mat_pu.Mearge(1,  &noes_b,  nno,noes,	  ndim,      &emat_pu[0][0]);
		mat_pp.Mearge(1,  &noes_b,  1,  &noes_b,  1,         &emat_pp);

		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
			res_u.AddValue( noes[ino],2,eres_u[ino][2]);
		}
		res_p.AddValue( noes_b,0,eres_p);
	}
	return true;
}

bool Fem::Eqn::AddLinSys_Stokes3D_Static(
		double alpha, 
		double rho, double g_x, double g_y, double g_z,
		CLinearSystem_Field& ls,
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType() != VECTOR3 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;
		
	const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if(    field_velo.GetInterpolationType( id_ea,world) == TET11
			&& field_press.GetInterpolationType(id_ea,world) == TET11 ){
			AddLinearSystem_Stokes3D_Static_P1P1(
				alpha,
				g_x,g_y,g_z,
				ls,
				id_field_velo,id_field_press,world,
                id_ea);
		}
		else if( field_velo.GetInterpolationType(id_ea,world) == HEX11 
			&&  field_press.GetInterpolationType(id_ea,world) == HEX1001 ){
			AddLinearSystem_Stokes3D_Static_Hex81(
				alpha,
				rho, g_x,g_y,g_z,
				ls,
				id_field_velo,id_field_press,world,
                id_ea);
		}
		else{ assert(0); }
	}
	return true;
}

