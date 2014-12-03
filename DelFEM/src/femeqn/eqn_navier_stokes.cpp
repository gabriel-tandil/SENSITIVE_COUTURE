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
// eqn_navier_stokes.cpp : ナビア・ストークス流体方程式の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <math.h>

#include "delfem/field_world.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femeqn/eqn_navier_stokes.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

void MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1(
    double dt, double gamma, 
    double rho, double myu, double g_x, double g_y,
    const double coords[][2],
    const double velo[][2], const double acc[][2], const double press[], const double apress[],
    double eres_u[3][2], double eres_p[3],
    double eCmat_uu[][3][2][2], double eCmat_up[][3][2], double eCmat_pu[][3][2], double eCmat_pp[][3],
    double eMmat_uu[][3][2][2], double eMmat_pu[][3][2])   
{
    const unsigned int nno = 3;
    const unsigned int ndim = 2;

	// 要素剛性行列、残差を０で初期化
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eCmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_up[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_pu[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno;          i++){ *(&eCmat_pp[0][0]+i)       = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eMmat_pu[0][0][0]+i)    = 0.0; }

	// 面積を求める
	const double area = TriArea(coords[0],coords[1],coords[2]);

	// 形状関数のｘｙ微分を求める
    double dldx[3][2];
    double const_term[3];
	TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

	// 粘性項
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eCmat_uu[ino][jno][0][0] += area*myu*dldx[ino][0]*dldx[jno][0];
		eCmat_uu[ino][jno][0][1] += area*myu*dldx[ino][1]*dldx[jno][0];
		eCmat_uu[ino][jno][1][0] += area*myu*dldx[ino][0]*dldx[jno][1];
		eCmat_uu[ino][jno][1][1] += area*myu*dldx[ino][1]*dldx[jno][1];
		const double dtmp1 = area*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		eCmat_uu[ino][jno][0][0] += dtmp1;
		eCmat_uu[ino][jno][1][1] += dtmp1;
	}
	}

	{	// 移流項を追加
		const double dtmp0[2] = { 
			velo[0][0]+velo[1][0]+velo[2][0], 
			velo[0][1]+velo[1][1]+velo[2][1] };
		for(unsigned int jno=0;jno<nno;jno++){
			const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]);
			for(unsigned int ino=0;ino<nno;ino++){
				double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]);
				dtmp2 *= area*rho*0.083333333333333;
				eCmat_uu[ino][jno][0][0] += dtmp2;
				eCmat_uu[ino][jno][1][1] += dtmp2;
			}
		}
	}

	{   // 圧力勾配・非圧縮項を追加
		const double dtmp1 = area*0.33333333333333333;
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_up[ino][jno][0] += -dtmp1*dldx[ino][0];
			eCmat_up[ino][jno][1] += -dtmp1*dldx[ino][1];
			eCmat_pu[ino][jno][0] +=  dtmp1*dldx[jno][0];
			eCmat_pu[ino][jno][1] +=  dtmp1*dldx[jno][1];
		}
		}
	}

	////////////////

	{	// 慣性行列を作る
		const double dtmp1 = area*rho*0.0833333333333333;
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eMmat_uu[ino][jno][0][0] += dtmp1;
				eMmat_uu[ino][jno][1][1] += dtmp1;
			}
			eMmat_uu[ino][ino][0][0] += dtmp1;
			eMmat_uu[ino][ino][1][1] += dtmp1;
		}
	}

    // 外力ベクトルを求める
	for(unsigned int ino=0;ino<nno;ino++){
		eres_u[ino][0] = area*rho*g_x*0.3333333333333333333;
		eres_u[ino][1] = area*rho*g_y*0.3333333333333333333;
	}
	for(unsigned int ino=0;ino<nno;ino++){
		eres_p[ino] = 0.0;
	}

	////////////////
	// Calc Stabilization Parameter
	double tau;
	{
		const double velo_ave[2] = { 
			(velo[0][0]+velo[1][0]+velo[2][0])/3.0, 
			(velo[0][1]+velo[1][1]+velo[2][1])/3.0 };
		const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
		const double h = sqrt( area / 3.14 )*2;
		const double tau_c = h*0.5/norm_v;
		const double cou_c = norm_v*dt/h;
		if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
			const double dtmp1 = 1/(cou_c*cou_c)+1;
			tau = tau_c / sqrt(dtmp1);
		}
		else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
			tau = h*h*rho*0.5/myu;
		}
		else{
			const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
			const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
			tau = tau_c / sqrt(dtmp1);
		}
	}


	// 慣性項に対するSUPGを求める
	for(unsigned int jno=0;jno<nno;jno++){
		double tmp_vec[ndim] = { 0.0, 0.0 };
		for(unsigned int kno=0;kno<nno;kno++){
			tmp_vec[0] += velo[kno][0];
			tmp_vec[1] += velo[kno][1];
		}
		tmp_vec[0] += velo[jno][0];
		tmp_vec[1] += velo[jno][1];
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1])*rho*tau*area*0.083333333333333;
			eMmat_uu[ino][jno][0][0] += dtmp1;
			eMmat_uu[ino][jno][1][1] += dtmp1;
		}
	}

    {	// SUPGの移流項
		double tmp_mat[ndim][ndim] = { {0,0}, {0,0} };
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
				tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
				tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
				tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
			}
			tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
			tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			double dtmp1 = 0.0;
			dtmp1 += dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
			        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
			        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
			        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1];
			dtmp1 *= tau*rho*area*0.083333333333333333333;
			eCmat_uu[ino][jno][0][0] += dtmp1;
			eCmat_uu[ino][jno][1][1] += dtmp1;
		}
		}
	}

	double ave_velo[2];
	{
		ave_velo[0] = 0.0;
		ave_velo[1] = 0.0;
		for(unsigned int ino=0;ino<nno;ino++){
			ave_velo[0] += velo[ino][0];
			ave_velo[1] += velo[ino][1];
		}
		ave_velo[0] *= 0.3333333333333333333;
		ave_velo[1] *= 0.3333333333333333333;
	}

    // SUPGの圧力勾配項への適応
	for(unsigned int ino=0;ino<nno;ino++){
		double dtmp1 = (dldx[ino][0]*ave_velo[0]+dldx[ino][1]*ave_velo[1])*tau*area;
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
			eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
		}
	}
/*		
		// SUPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
			eres_u[ino][0] += dtmp1*g_x;
			eres_u[ino][1] += dtmp1*g_y;
		}
*/
	// PSPGの慣性項への適応を代入
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eMmat_pu[ino][jno][0] -= tau*area*0.33333333333333333*dldx[ino][0];
		eMmat_pu[ino][jno][1] -= tau*area*0.33333333333333333*dldx[ino][1];
	}
	}

	// PSPGの移流項への適応
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = (dldx[jno][0]*ave_velo[0]+dldx[jno][1]*ave_velo[1])*tau*area;
		for(unsigned int ino=0;ino<nno;ino++){
			eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
			eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
		}
	}

	// PSPGの圧力項
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eCmat_pp[ino][jno] += area*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
	}
	}
/*
		// PSPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
		}
*/
}


static bool AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1(
		double gamma, double dt, CLinearSystem_Field& ls, 
		double rho, double myu, double g_x, double g_y, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea )
{
//    std::cout << "NavierStorkes2D_NonStatic_Newmark Triangle 3-point 1st order " << gamma << " " << dt << " " << rho << " " << myu << " " << id_ea << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg es_velo_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg es_velo_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg es_pres_c_va = field_press.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
        !=  ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER,  world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo,CORNER, id_field_press,CORNER,  world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo,CORNER,  world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,  world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,  world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no_v[nno];	// 要素節点の全体節点番号				
		es_velo_c_co.GetNodes(ielem,no_v);	// 要素の節点番号を取ってくる
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
		es_pres_c_va.GetNodes(ielem,no_p);	// 要素の節点番号を取ってくる
		double press[nno];
		double apress[nno];		
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_press.GetValue(no_p[ino],&press[ino]);
			ns_apress.GetValue(no_p[ino],&apress[ino]);
		}

        ////////////////////////////////

        MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1(dt,gamma,
            rho,myu,g_x,g_y,
            coords,   velo,acc,   press,apress,
            eres_u,eres_p, 
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, eMmat_uu,eMmat_pu);

        ////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
			}
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,no_v,nno,no_v,	4,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,no_v,nno,no_p,	2,&emat_up[0][0][0]   );
		mat_pu.Mearge(nno,no_p,nno,no_v,	2,&emat_pu[0][0][0]   );
		mat_pp.Mearge(nno,no_p,nno,no_p,	1,&emat_pp[0][0]      );
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( no_v[ino],0,eres_u[ino][0]);
			res_u.AddValue( no_v[ino],1,eres_u[ino][1]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue( no_p[ino],0,eres_p[ino]);
		}
	}
	return true;
}


static bool AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1_Combined(
		double gamma, double dt, CLinearSystem_Field& ls, 
		double rho, double myu, double g_x, double g_y, 
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world, 
		const unsigned int id_ea )
{
//    std::cout << "NavierStorkes2D_NonStatic_Newmark Triangle 3-point 1st order Combined" << gamma << " " << dt << " " << rho << " " << myu << " " << id_ea << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg es_c_co = field_velo.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(  id_field_velo, CORNER,world);
	CVector_Blk&    res_u  = ls.GetResidual(id_field_velo, CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int noes[nno];	// 要素節点の全体節点番号
		
		// 要素の節点番号を取ってくる
		es_c_co.GetNodes(ielem,noes);
		double coords[nno][ndim];	// 要素節点の座標
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
		}
        // 要素の節点番号を取ってくる(id_esは流速と圧力で同じはず)
		es_c_va.GetNodes(ielem,noes);
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


        MakeMat_NavierStokes2D_NonStatic_Newmark_P1P1(dt,gamma,
            rho,myu,g_x,g_y,
            coords,   velo,acc,   press,apress,
            eres_u,eres_p, 
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, eMmat_uu,eMmat_pu);


		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
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


bool Fem::Eqn::AddLinSys_NavierStokes2D_NonStatic_Newmark(
		double dt, double gamma, CLinearSystem_Field& ls,
		double rho, double alpha, double g_x, double g_y,
		const unsigned int id_field_velo, unsigned int id_field_press, const CFieldWorld& world,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;

	if( id_ea != 0 ){
		bool res;
		if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_velo,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_press,CORNER,world) )
            {
			    res = AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1_Combined(
				    gamma,dt,  ls,
				    rho,alpha, g_x,g_y,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
            else{
			    res = AddLinSys_NavierStokes2D_NonStatic_Newmark_P1P1(
				    gamma,dt,  ls,
				    rho,alpha, g_x,g_y,
				    id_field_velo,id_field_press,world,
                    id_ea);
            }
		}
		else{
			assert(0); 
			res = false;
		}
		return res;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰文
			bool res = AddLinSys_NavierStokes2D_NonStatic_Newmark(
				dt, gamma, ls,
				rho, alpha, g_x, g_y,
				id_field_velo, id_field_press, world,
				id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}










static bool AddLinearSystem_NavierStokesALE2D_NonStatic_Newmark_P1P1(
		double rho, double myu, double g_x, double g_y, 
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_msh_velo,
		const CFieldWorld& world, 
		unsigned int id_ea)
{
//	std::cout << "NavierStorkes2D_NonStatic_Newmark Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_msh_velo) ) return false;
	const CField& field_msh_velo = world.GetField(id_field_msh_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	const CElemAry::CElemSeg es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	unsigned int noes[nno];	// 要素節点の全体節点番号

	double velo[nno][ndim];	// 要素節点の値
	double velo_msh[nno][ndim];	// 要素節点の値
	double velo_r[nno][ndim];	// relative velocity
	double acc[nno][ndim];	// 要素節点の値
	double press[nno];
	double apress[nno];
	double coords[nno][ndim];	// 要素節点の座標
				
	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo, CORNER, world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER, world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo, CORNER, id_field_press,CORNER, world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo, CORNER, world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER, world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER, world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_msh_velo = field_msh_velo.GetNodeSeg(CORNER,true, world,VELOCITY);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_c_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_acc.GetValue(noes[ino],acc[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
			ns_apress.GetValue(noes[ino],&apress[ino]);
			ns_msh_velo.GetValue(noes[ino],velo_msh[ino]);
		}
		for(unsigned int i=0;i<nno*ndim;i++){ 
			(&velo_r[0][0])[i] = (&velo[0][0])[i]-(&velo_msh[0][0])[i];
		}

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eCmat_uu[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_up[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_pu[0][0][0]+i)    = 0.0; }
		for(unsigned int i=0;i<nno*nno;          i++){ *(&eCmat_pp[0][0]+i)       = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_uu[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eMmat_pu[0][0][0]+i)    = 0.0; }

		// 面積を求める
		const double area = TriArea(coords[0],coords[1],coords[2]);
		// 形状関数のｘｙ微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

		// 粘性項
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_uu[ino][jno][0][0] += area*myu*dldx[ino][0]*dldx[jno][0];
			eCmat_uu[ino][jno][0][1] += area*myu*dldx[ino][1]*dldx[jno][0];
			eCmat_uu[ino][jno][1][0] += area*myu*dldx[ino][0]*dldx[jno][1];
			eCmat_uu[ino][jno][1][1] += area*myu*dldx[ino][1]*dldx[jno][1];
			const double dtmp1 = area*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
			eCmat_uu[ino][jno][0][0] += dtmp1;
			eCmat_uu[ino][jno][1][1] += dtmp1;
		}
		}

		{	// 移流項を追加
			const double dtmp0[2] = { 
				velo_r[0][0]+velo_r[1][0]+velo_r[2][0], 
				velo_r[0][1]+velo_r[1][1]+velo_r[2][1] };
			for(unsigned int jno=0;jno<nno;jno++){
				const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]);
				for(unsigned int ino=0;ino<nno;ino++){
					double dtmp2 = dtmp1 + (dldx[jno][0]*velo_r[ino][0]+dldx[jno][1]*velo_r[ino][1]);
					dtmp2 *= area*rho*0.083333333333333;
					eCmat_uu[ino][jno][0][0] += dtmp2;
					eCmat_uu[ino][jno][1][1] += dtmp2;
				}
			}
		}

		{
			const double dtmp1 = area*0.33333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eCmat_up[ino][jno][0] += -dtmp1*dldx[ino][0];
				eCmat_up[ino][jno][1] += -dtmp1*dldx[ino][1];
				eCmat_pu[ino][jno][0] +=  dtmp1*dldx[jno][0];
				eCmat_pu[ino][jno][1] +=  dtmp1*dldx[jno][1];
			}
			}
		}

		////////////////

		{	// 慣性行列を作る
			const double dtmp1 = area*rho*0.0833333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat_uu[ino][jno][0][0] += dtmp1;
					eMmat_uu[ino][jno][1][1] += dtmp1;
				}
				eMmat_uu[ino][ino][0][0] += dtmp1;
				eMmat_uu[ino][ino][1][1] += dtmp1;
			}
		}

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_u[ino][0] = area*rho*g_x*0.3333333333333333333;
			eres_u[ino][1] = area*rho*g_y*0.3333333333333333333;
		}
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] = 0.0;
		}

		////////////////

		// Calc Stabilization Parameter
		double tau;
		{
			const double velo_ave[2] = { 
				(velo_r[0][0]+velo_r[1][0]+velo_r[2][0])/3.0, 
				(velo_r[0][1]+velo_r[1][1]+velo_r[2][1])/3.0 };
			const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
			const double h = sqrt( area / 3.14 )*2;
			const double tau_c = h*0.5/norm_v;
			const double cou_c = norm_v*dt/h;
			if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
				const double dtmp1 = 1/(cou_c*cou_c)+1;
				tau = tau_c / sqrt(dtmp1);
			}
			else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
				tau = h*h*rho*0.5/myu;
			}
			else{
				const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
				const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
				tau = tau_c / sqrt(dtmp1);
			}
		}


		// 慣性項に対するSUPGを求める
		for(unsigned int jno=0;jno<nno;jno++){
			double tmp_vec[ndim] = { 0.0, 0.0 };
			for(unsigned int kno=0;kno<nno;kno++){
				tmp_vec[0] += velo_r[kno][0];
				tmp_vec[1] += velo_r[kno][1];
			}
			tmp_vec[0] += velo_r[jno][0];
			tmp_vec[1] += velo_r[jno][1];
			for(unsigned int ino=0;ino<nno;ino++){
				const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1])*rho*tau*area*0.083333333333333;
				eMmat_uu[ino][jno][0][0] += dtmp1;
				eMmat_uu[ino][jno][1][1] += dtmp1;
			}
		}

		{	// SUPGの移流項
			double tmp_mat[ndim][ndim] = { {0,0}, {0,0} };
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					tmp_mat[0][0] += velo_r[ino][0]*velo_r[jno][0];
					tmp_mat[0][1] += velo_r[ino][0]*velo_r[jno][1];
					tmp_mat[1][0] += velo_r[ino][1]*velo_r[jno][0];
					tmp_mat[1][1] += velo_r[ino][1]*velo_r[jno][1];
				}
				tmp_mat[0][0] += velo_r[ino][0]*velo_r[ino][0];
				tmp_mat[1][1] += velo_r[ino][1]*velo_r[ino][1];
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				dtmp1 += dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
				        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
				        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
				        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1];
				dtmp1 *= tau*rho*area*0.083333333333333333333;
				eCmat_uu[ino][jno][0][0] += dtmp1;
				eCmat_uu[ino][jno][1][1] += dtmp1;
			}
			}
		}

		double ave_velo_r[2];
		{
			ave_velo_r[0] = 0.0;
			ave_velo_r[1] = 0.0;
			for(unsigned int ino=0;ino<nno;ino++){
				ave_velo_r[0] += velo_r[ino][0];
				ave_velo_r[1] += velo_r[ino][1];
			}
			ave_velo_r[0] *= 0.3333333333333333333;
			ave_velo_r[1] *= 0.3333333333333333333;
		}

		// SUPGの圧力勾配項への適応
		for(unsigned int ino=0;ino<nno;ino++){
			double dtmp1 = (dldx[ino][0]*ave_velo_r[0]+dldx[ino][1]*ave_velo_r[1])*tau*area;
			for(unsigned int jno=0;jno<nno;jno++){
				eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
				eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
			}
		}
		/*
		// SUPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = area*tau*rho*(ave_velo_r[0]*dldx[ino][0]+ave_velo_r[1]*dldx[ino][1]);
			eres_u[ino][0] += dtmp1*g_x;
			eres_u[ino][1] += dtmp1*g_y;
		}
		*/

		// PSPGの慣性項への適応を代入
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eMmat_pu[ino][jno][0] -= tau*area*0.33333333333333333*dldx[ino][0];
			eMmat_pu[ino][jno][1] -= tau*area*0.33333333333333333*dldx[ino][1];
		}
		}

		// PSPGの移流項への適応
		for(unsigned int jno=0;jno<nno;jno++){
			const double dtmp1 = (dldx[jno][0]*ave_velo_r[0]+dldx[jno][1]*ave_velo_r[1])*tau*area;
			for(unsigned int ino=0;ino<nno;ino++){
				eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
				eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
			}
		}

		// PSPGの圧力項
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_pp[ino][jno] += area*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		}
		}

		/*
		// PSPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
		}
		*/

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
			}
		}

		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,noes,nno,noes,	4,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,noes,nno,noes,	2,&emat_up[0][0][0]   );
		mat_pu.Mearge(nno,noes,nno,noes,	2,&emat_pu[0][0][0]   );
		mat_pp.Mearge(nno,noes,nno,noes,	1,&emat_pp[0][0]      );

		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue( noes[ino],0,eres_p[ino]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_NavierStokesALE2D_NonStatic_Newmark(
		double dt, double gamma, CLinearSystem_Field& ls,
		double rho, double alpha, double g_x, double g_y,
		unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_msh_velo,
		const CFieldWorld& world )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( !world.IsIdField(id_field_msh_velo) ) return false;
	const CField& field_msh_velo = world.GetField(id_field_msh_velo);

	if( field_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_msh_velo.GetFieldType() != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR ) return false;

	const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinearSystem_NavierStokesALE2D_NonStatic_Newmark_P1P1(
				rho,alpha, g_x,g_y,
				gamma,dt,
				ls,
				id_field_velo,id_field_press,id_field_msh_velo,
				world,id_ea);
		}
	}
	return true;
}




void MakeMat_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1(
    double dt,  // 安定化のために必要
    double rho, double myu, double g_x, double g_y,
    const double coords[][2],
    const double velo[][2],
    const double temp[],
    double eCmat_uu[][3][2][2], double eCmat_up[][3][2], double eCmat_pu[][3][2], double eCmat_pp[][3],
    double eMmat_uu[][3][2][2], double eMmat_pu[][3][2], 
    double eres_u[][2], double eres_p[])
{

    const unsigned int nno = 3;
    const unsigned int ndim = 2;
    
    // 要素剛性行列、残差を０で初期化
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eCmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_up[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eCmat_pu[0][0][0]+i)    = 0.0; }
	for(unsigned int i=0;i<nno*nno;          i++){ *(&eCmat_pp[0][0]+i)       = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_uu[0][0][0][0]+i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*ndim;     i++){ *(&eMmat_pu[0][0][0]+i)    = 0.0; }

	// 面積を求める
	const double area = TriArea(coords[0],coords[1],coords[2]);

	// 形状関数のｘｙ微分を求める
	double dldx[nno][ndim];	// 形状関数のxy微分
	double const_term[nno];	// 形状関数の定数項
	TriDlDx(dldx, const_term,   coords[0], coords[1], coords[2]);

	// 粘性項
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = area*myu*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
		eCmat_uu[ino][jno][0][0] += dtmp1;
		eCmat_uu[ino][jno][1][1] += dtmp1;
	}
	}

	{	// 移流項を追加
		const double dtmp0[2] = { 
			velo[0][0]+velo[1][0]+velo[2][0], 
			velo[0][1]+velo[1][1]+velo[2][1] };
		for(unsigned int jno=0;jno<nno;jno++){
			const double dtmp1 = (dldx[jno][0]*dtmp0[0]+dldx[jno][1]*dtmp0[1]);
			for(unsigned int ino=0;ino<nno;ino++){
				double dtmp2 = dtmp1 + (dldx[jno][0]*velo[ino][0]+dldx[jno][1]*velo[ino][1]);
				dtmp2 *= area*rho*0.083333333333333;
				eCmat_uu[ino][jno][0][0] += dtmp2;
				eCmat_uu[ino][jno][1][1] += dtmp2;
			}
		}
	}

	{
		const double dtmp1 = area*0.33333333333333333;
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_up[ino][jno][0] += -dtmp1*dldx[ino][0];
			eCmat_up[ino][jno][1] += -dtmp1*dldx[ino][1];
			eCmat_pu[ino][jno][0] +=  dtmp1*dldx[jno][0];
			eCmat_pu[ino][jno][1] +=  dtmp1*dldx[jno][1];
		}
		}
	}

	////////////////

	{	// 慣性行列を作る
		const double dtmp1 = area*rho*0.0833333333333333;
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eMmat_uu[ino][jno][0][0] += dtmp1;
				eMmat_uu[ino][jno][1][1] += dtmp1;
			}
			eMmat_uu[ino][ino][0][0] += dtmp1;
			eMmat_uu[ino][ino][1][1] += dtmp1;
		}
	}

	// 外力ベクトルを求める
	for(unsigned int ino=0;ino<nno;ino++){
		eres_u[ino][0] = area*rho*g_x*0.3333333333333333333;
		eres_u[ino][1] = area*rho*g_y*0.3333333333333333333;
	}
	for(unsigned int ino=0;ino<nno;ino++){
		eres_p[ino] = 0.0;
	}
	for(unsigned int ino=0;ino<nno;ino++){
		eres_u[ino][1] += area*(temp[ino]+temp[0]+temp[1]+temp[2])*0.083333333333333;
	}

	////////////////

	// Calc Stabilization Parameter
	double tau;
	{
		const double velo_ave[2] = { 
			(velo[0][0]+velo[1][0]+velo[2][0])/3.0, 
			(velo[0][1]+velo[1][1]+velo[2][1])/3.0 };
		const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
		const double h = sqrt( area / 3.14 )*2;
		const double tau_c = h*0.5/norm_v;
		const double cou_c = norm_v*dt/h;
		if( norm_v*h*rho*1.0e-30 > myu ){ // Re = \infty
			const double dtmp1 = 1/(cou_c*cou_c)+1;
			tau = tau_c / sqrt(dtmp1);
		}
		else if( norm_v*h*rho < myu*1.0e-30 ){ // Re = 0
			tau = h*h*rho*0.5/myu;
		}
		else{
			const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
			const double dtmp1 = 1/(cou_c*cou_c)+1+1/(re_c*re_c);
			tau = tau_c / sqrt(dtmp1);
		}
	}


	// 慣性項に対するSUPGを求める
	for(unsigned int jno=0;jno<nno;jno++){
		double tmp_vec[ndim] = { 0.0, 0.0 };
		for(unsigned int kno=0;kno<nno;kno++){
			tmp_vec[0] += velo[kno][0];
			tmp_vec[1] += velo[kno][1];
		}
		tmp_vec[0] += velo[jno][0];
		tmp_vec[1] += velo[jno][1];
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = (dldx[ino][0]*tmp_vec[0]+dldx[ino][1]*tmp_vec[1])*rho*tau*area*0.083333333333333;
			eMmat_uu[ino][jno][0][0] += dtmp1;
			eMmat_uu[ino][jno][1][1] += dtmp1;
		}
	}

    {	// SUPGの移流項
		double tmp_mat[ndim][ndim] = { {0,0}, {0,0} };
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				tmp_mat[0][0] += velo[ino][0]*velo[jno][0];
				tmp_mat[0][1] += velo[ino][0]*velo[jno][1];
				tmp_mat[1][0] += velo[ino][1]*velo[jno][0];
				tmp_mat[1][1] += velo[ino][1]*velo[jno][1];
			}
			tmp_mat[0][0] += velo[ino][0]*velo[ino][0];
			tmp_mat[1][1] += velo[ino][1]*velo[ino][1];
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			double dtmp1 = 0.0;
			dtmp1 += dldx[ino][0]*dldx[jno][0]*tmp_mat[0][0]
			        +dldx[ino][0]*dldx[jno][1]*tmp_mat[0][1]
			        +dldx[ino][1]*dldx[jno][0]*tmp_mat[1][0]
			        +dldx[ino][1]*dldx[jno][1]*tmp_mat[1][1];
			dtmp1 *= tau*rho*area*0.083333333333333333333;
			eCmat_uu[ino][jno][0][0] += dtmp1;
			eCmat_uu[ino][jno][1][1] += dtmp1;
		}
		}
	}

	double ave_velo[2];
	{
		ave_velo[0] = 0.0;
		ave_velo[1] = 0.0;
		for(unsigned int ino=0;ino<nno;ino++){
			ave_velo[0] += velo[ino][0];
			ave_velo[1] += velo[ino][1];
		}
		ave_velo[0] *= 0.3333333333333333333;
		ave_velo[1] *= 0.3333333333333333333;
	}

    // SUPGの圧力勾配項への適応
	for(unsigned int ino=0;ino<nno;ino++){
		double dtmp1 = (dldx[ino][0]*ave_velo[0]+dldx[ino][1]*ave_velo[1])*tau*area;
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat_up[ino][jno][0] += dtmp1*dldx[jno][0];
			eCmat_up[ino][jno][1] += dtmp1*dldx[jno][1];
		}
	}
/*		
		// SUPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			const double dtmp1 = area*tau*rho*(ave_velo[0]*dldx[ino][0]+ave_velo[1]*dldx[ino][1]);
			eres_u[ino][0] += dtmp1*g_x;
			eres_u[ino][1] += dtmp1*g_y;
		}
*/
	// PSPGの慣性項への適応を代入
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eMmat_pu[ino][jno][0] -= tau*area*0.33333333333333333*dldx[ino][0];
		eMmat_pu[ino][jno][1] -= tau*area*0.33333333333333333*dldx[ino][1];
	}
	}

	// PSPGの移流項への適応
	for(unsigned int jno=0;jno<nno;jno++){
		const double dtmp1 = (dldx[jno][0]*ave_velo[0]+dldx[jno][1]*ave_velo[1])*tau*area;
		for(unsigned int ino=0;ino<nno;ino++){
			eCmat_pu[ino][jno][0] += dtmp1*dldx[ino][0];
			eCmat_pu[ino][jno][1] += dtmp1*dldx[ino][1];
		}
	}

	// PSPGの圧力項
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		eCmat_pp[ino][jno] += area*tau/rho*(dldx[jno][0]*dldx[ino][0]+dldx[jno][1]*dldx[ino][1]);
	}
	}
/*
		// PSPGの外力項
		for(unsigned int ino=0;ino<nno;ino++){
			eres_p[ino] += area*tau*(dldx[ino][0]*g_x+dldx[ino][1]*g_y);
		}
*/



}











static bool AddLinearSystem_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1(
		double rho, double myu, double g_x, double g_y, 
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_temp,
		const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "NavierStorkesThermalBuoy2D_NonStatic_Newmark Tri11" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( !world.IsIdField(id_field_temp) ) return false;
	const CField& field_temp = world.GetField(id_field_temp);

	const CElemAry::CElemSeg es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 != ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

	CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo,CORNER,world);
	CMatDia_BlkCrs& mat_pp = ls.GetMatrix(id_field_press,CORNER,world);
	CMat_BlkCrs& mat_up = ls.GetMatrix(id_field_velo,CORNER, id_field_press,CORNER, world);
	CMat_BlkCrs& mat_pu = ls.GetMatrix(id_field_press,CORNER, id_field_velo,CORNER, world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,world);
	CVector_Blk& res_p = ls.GetResidual(id_field_press,CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);
	const CNodeAry::CNodeSeg& ns_temp = field_temp.GetNodeSeg(CORNER,true, world,VALUE);//na_velo.GetSeg(id_ns_velo);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		unsigned int noes[nno];	// 要素節点の全体節点番号
		es_c_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		double press[nno];
		double apress[nno];
		double temp[nno];
		double coords[nno][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_acc.GetValue(noes[ino],acc[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
			ns_apress.GetValue(noes[ino],&apress[ino]);
			ns_temp.GetValue(noes[ino],&temp[ino]);
		}

        ////////////////////////////////

        MakeMat_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1(dt,rho,myu,g_x,g_y,
            coords,velo,temp,
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, 
            eMmat_uu,eMmat_pu,
            eres_u,eres_p);

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
			}
		}
    
		// 要素剛性行列の全体剛性行列へのマージ
		mat_uu.Mearge(nno,noes,nno,noes,	4,&emat_uu[0][0][0][0]);
		mat_up.Mearge(nno,noes,nno,noes,	2,&emat_up[0][0][0]   );
		mat_pu.Mearge(nno,noes,nno,noes,	2,&emat_pu[0][0][0]   );
		mat_pp.Mearge(nno,noes,nno,noes,	1,&emat_pp[0][0]      );
		// 残差ベクトルのマージ
		for(unsigned int ino=0;ino<nno;ino++){
			res_u.AddValue( noes[ino],0,eres_u[ino][0]);
			res_u.AddValue( noes[ino],1,eres_u[ino][1]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_p.AddValue( noes[ino],0,eres_p[ino]);
		}
	}
	return true;
}

static bool AddLinearSystem_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1_Combined(
		double rho, double myu, double g_x, double g_y, 
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_temp,
		const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "NavierStorkesThermalBuoy2D_NonStatic_Newmark Tri11 Combined" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( !world.IsIdField(id_field_temp) ) return false;
	const CField& field_temp = world.GetField(id_field_temp);

	const CElemAry::CElemSeg es_c_va = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	double eCmat_uu[nno][nno][ndim][ndim], eCmat_pp[nno][nno], eCmat_pu[nno][nno][ndim], eCmat_up[nno][nno][ndim];
	double eMmat_uu[nno][nno][ndim][ndim], eMmat_pu[nno][nno][ndim];
	double  emat_uu[nno][nno][ndim][ndim],  emat_pp[nno][nno],  emat_pu[nno][nno][ndim],  emat_up[nno][nno][ndim];
	double eres_u[nno][ndim], eres_p[nno];

	assert( ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_press,CORNER,world) );

    CMatDia_BlkCrs& mat_uu = ls.GetMatrix(id_field_velo,CORNER,world);
	CVector_Blk& res_u = ls.GetResidual(id_field_velo, CORNER,world);

	const CNodeAry::CNodeSeg& ns_co   = field_velo.GetNodeSeg(CORNER,false,world,VALUE);//na_co.GetSeg(id_ns_co);
	const CNodeAry::CNodeSeg& ns_velo = field_velo.GetNodeSeg(CORNER,true, world,VELOCITY);//na_velo.GetSeg(id_ns_velo);
	const CNodeAry::CNodeSeg& ns_acc  = field_velo.GetNodeSeg(CORNER,true, world,ACCELERATION);//na_velo.GetSeg(id_ns_acc);
	const CNodeAry::CNodeSeg& ns_press  = field_press.GetNodeSeg(CORNER,true,world,VELOCITY);//na_press.GetSeg(id_ns_press);
	const CNodeAry::CNodeSeg& ns_apress = field_press.GetNodeSeg(CORNER,true,world,ACCELERATION);//na_press.GetSeg(id_ns_apress);
	const CNodeAry::CNodeSeg& ns_temp = field_temp.GetNodeSeg(CORNER,true, world,VALUE);//na_velo.GetSeg(id_ns_velo);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		unsigned int noes[nno];	// 要素節点の全体節点番号
		es_c_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		double velo[nno][ndim];	// 要素節点の値
		double acc[nno][ndim];	// 要素節点の値
		double press[nno];
		double apress[nno];
		double temp[nno];
		double coords[nno][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_co.GetValue(noes[ino],coords[ino]);
			ns_velo.GetValue(noes[ino],velo[ino]);
			ns_acc.GetValue(noes[ino],acc[ino]);
			ns_press.GetValue(noes[ino],&press[ino]);
			ns_apress.GetValue(noes[ino],&apress[ino]);
			ns_temp.GetValue(noes[ino],&temp[ino]);
		}

        ////////////////////////////////

        MakeMat_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1(dt,rho,myu,g_x,g_y,
            coords,velo,temp,
            eCmat_uu,eCmat_up,eCmat_pu,eCmat_pp, 
            eMmat_uu,eMmat_pu,
            eres_u,eres_p);

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_u[ino][0] -= eCmat_uu[ino][jno][0][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][0][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][0][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][0][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][0]*(press[jno]+dt*apress[jno]);
			eres_u[ino][1] -= eCmat_uu[ino][jno][1][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_uu[ino][jno][1][0]*acc[jno][0]
			                 +eCmat_uu[ino][jno][1][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_uu[ino][jno][1][1]*acc[jno][1]
			                 +eCmat_up[ino][jno][1]*(press[jno]+dt*apress[jno]);
		}
		}

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_p[ino] -= eCmat_pu[ino][jno][0]*(velo[jno][0]+dt*acc[jno][0]) + eMmat_pu[ino][jno][0]*acc[jno][0]
			              +eCmat_pu[ino][jno][1]*(velo[jno][1]+dt*acc[jno][1]) + eMmat_pu[ino][jno][1]*acc[jno][1]
			              +eCmat_pp[ino][jno]*(press[jno]+dt*apress[jno]);
		}
		}

		////////////////////////////////

		{
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_uu[0][0][0][0])[i] = (&eMmat_uu[0][0][0][0])[i]+dtmp1*(&eCmat_uu[0][0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_up[0][0][0])[i] = dtmp1*(&eCmat_up[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno*ndim;i++){
				(&emat_pu[0][0][0])[i] = (&eMmat_pu[0][0][0])[i]+dtmp1*(&eCmat_pu[0][0][0])[i];
			}
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat_pp[0][0])[i] = dtmp1*(&eCmat_pp[0][0])[i];
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


bool Fem::Eqn::AddLinSys_NavierStokes2DThermalBuoy_NonStatic_Newmark(
	double rho, double alpha, double g_x, double g_y,
	double gamma, double dt,
	Fem::Ls::CLinearSystem_Field& ls, 
	const unsigned int id_field_velo, unsigned int id_field_press, unsigned int id_field_temp,
	const Fem::Field::CFieldWorld& world )
{
	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( !world.IsIdField(id_field_press) ) return false;
	const CField& field_press = world.GetField(id_field_press);

	if( field_velo.GetFieldType()  != VECTOR2 ) return false;
	if( field_press.GetFieldType() != SCALAR  ) return false;

	const std::vector<unsigned int>& aIdEA = field_velo.GetAryIdEA();
	if(    ls.FindIndexArray_Seg(id_field_velo, CORNER,world) 
		== ls.FindIndexArray_Seg(id_field_press,CORNER,world) ){
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
				AddLinearSystem_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1_Combined(
					rho,alpha, g_x,g_y,
					gamma,dt,
					ls, id_field_velo,id_field_press,id_field_temp,
					world,id_ea);
			}
		}
		return true;
	}
	////////////////
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		const unsigned int id_ea = aIdEA[iiea];
		if( field_velo.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinearSystem_NavierStokesThermalBuoy2D_NonStatic_Newmark_P1P1(
				rho,alpha, g_x,g_y,
				gamma,dt,
				ls, id_field_velo,id_field_press,id_field_temp,
				world,id_ea);
		}
	}
	return true;

}


