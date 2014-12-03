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
// eqn_diffusion.cpp : 拡散方程式の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#include "delfem/field_world.h"

#include "delfem/femeqn/eqn_diffusion.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/bcflag_blk.h"

#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

void Set_RhoTri_CB_Scalar(double mat_cc[][3],  double mat_cb[], double mat_bc[], double& mat_bb,
				  const double area, const double dldx[][2], const double vc_b[],
				  const double rho){
	const double tmp1 = rho * area;
	const double tmp2 = 1.0 / 12.0;
	const double tmp3 = vc_b[3] / 180.0;
	const double tmp4 = vc_b[3] * vc_b[3] / 2520.0;

	int inoel,jnoel;
	double tmp_val1;

	for(inoel=0;inoel<3;inoel++){
		for(jnoel=0;jnoel<3;jnoel++){
			tmp_val1 = tmp1*( tmp2 - (vc_b[inoel]+vc_b[jnoel])*tmp3 + vc_b[inoel]*vc_b[jnoel]*tmp4 );
			mat_cc[inoel][jnoel] = tmp_val1;
		}
		mat_cc[inoel][inoel] = tmp1*tmp2;
	}
	for(inoel=0;inoel<3;inoel++){
		tmp_val1 = tmp1*(tmp3-vc_b[inoel]*tmp4);
		mat_bc[inoel] = tmp_val1;
		mat_cb[inoel] = tmp_val1;
	}
	tmp_val1 = tmp1*tmp4;
	mat_bb = tmp_val1;
	/*
	const double tmp1 = rho*area;
	const double tmp2 = tmp1*-19.0/280.0;
	const double tmp3 = tmp1/12.0;
	const double tmp4 = tmp1*15.0/280.0;
	const double tmp5 = tmp1*81.0/280.0;

	int inoel,jnoel;
	for(inoel=0;inoel<3;inoel++){
		for(jnoel=0;jnoel<3;jnoel++){
			mat_cc[inoel][jnoel] = tmp3+tmp2;
		}
		mat_cc[inoel][inoel] += tmp3;
	}
	for(inoel=0;inoel<3;inoel++){
		mat_cb[inoel] = tmp4;
		mat_bc[inoel] = tmp4;
	}
	mat_bb = tmp5;*/
}


static bool AddLinearSystem_Diffusion_P1b(
		double rho, double alpha, double source,
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea )
{	
//	std::cout << "Diffusion2D Tri P1b" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);
	const CElemAry::CElemSeg& es_b = field_val.GetElemSeg(id_ea,BUBBLE,true,world);

	const unsigned int nno_c = 3;
	const unsigned int nno_b = 1;
	const unsigned int ndim = 2;

	unsigned int no_c[nno_c];
	unsigned int no_b;

	double val_c[nno_c], val_b;
	double vval_c[nno_c], vval_b;
	double coord_c[nno_c][ndim];
				
	double dldx[nno_c][ndim];
	double const_term[nno_c];

	double eCmat_cc[nno_c][nno_c], eCmat_cb[nno_c], eCmat_bc[nno_c], eCmat_bb;
	double eMmat_cc[nno_c][nno_c], eMmat_cb[nno_c], eMmat_bc[nno_c], eMmat_bb;
	double eqf_out_c[nno_c], eqf_out_b;
	double eqf_in_c[nno_c], eqf_in_b;

	double emat_cc[nno_c][nno_c], emat_cb[nno_c], emat_bc[nno_c], emat_bb;
	double eres_c[nno_c], eres_b;	// 要素節点等価内力、外力、残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(id_field_val,CORNER, world);
	CMatDia_BlkCrs& mat_bb = ls.GetMatrix(id_field_val,BUBBLE, world);
	CMat_BlkCrs&    mat_cb = ls.GetMatrix(id_field_val,CORNER, id_field_val, BUBBLE, world);
	CMat_BlkCrs&    mat_bc = ls.GetMatrix(id_field_val,BUBBLE, id_field_val, CORNER, world);
	////////////////
	CVector_Blk& res_c = ls.GetResidual(id_field_val,CORNER, world);
	CVector_Blk& res_b = ls.GetResidual(id_field_val,BUBBLE, world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);//na_c_val.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_vval = field_val.GetNodeSeg(CORNER,true,world,VELOCITY);//na_c_val.GetSeg(id_ns_c_vval);
	const CNodeAry::CNodeSeg& ns_b_val = field_val.GetNodeSeg(BUBBLE,true,world,VALUE);//na_b_val.GetSeg(id_ns_b_val);
	const CNodeAry::CNodeSeg& ns_b_vval = field_val.GetNodeSeg(BUBBLE,true,world,VELOCITY);//na_b_val.GetSeg(id_ns_b_vval);
	const CNodeAry::CNodeSeg& ns_c_co = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_val.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から節点セグメントの節点番号を取り出す
		es_c.GetNodes(ielem,no_c);
		es_b.GetNodes(ielem,&no_b);
		// 節点の値を取ってくる
		for(unsigned int inoes=0;inoes<nno_c;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
			ns_c_val.GetValue(no_c[inoes],&val_c[inoes]);
			ns_c_vval.GetValue(no_c[inoes],&vval_c[inoes]);
		}
		ns_b_val.GetValue(no_b,&val_b);
		ns_b_vval.GetValue(no_b,&vval_b);

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);
		{	// 要素剛性行列を作る
			double vc_b[4];
			vc_b[0] = 1.0/3.0; vc_b[1] = 1.0/3.0; vc_b[2] = 1.0/3.0; vc_b[3] = 27.0;
			const double tmp_val1 = vc_b[3]*vc_b[3]*alpha*area/180.0*( 
				dldx[0][0]*dldx[0][0]+dldx[0][1]*dldx[0][1]+
				dldx[1][0]*dldx[1][0]+dldx[1][1]*dldx[1][1]+
				dldx[2][0]*dldx[2][0]+dldx[2][1]*dldx[2][1] );
            for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
            for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
				eCmat_cc[ino_c][jno_c] 
					= alpha*area*(dldx[ino_c][0]*dldx[jno_c][0]+dldx[ino_c][1]*dldx[jno_c][1])
					+vc_b[ino_c]*vc_b[jno_c]*tmp_val1;
			}
			}
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				const double tmp1 = -1.0*vc_b[ino_c]*tmp_val1;
				eCmat_cb[ino_c] = tmp1;
				eCmat_bc[ino_c] = tmp1;
			}
			eCmat_bb = tmp_val1;
			Set_RhoTri_CB_Scalar(eMmat_cc,eMmat_cb,eMmat_bc,eMmat_bb,  area, dldx,vc_b, rho);
		}
		// 要素外力ベクトルを求める
        for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
			eqf_out_c[ino_c] = source*area*11.0/60.0;
		}
		eqf_out_b = source*area*27.0/60.0;

		////////////////////////////////////////////////////////////////

		// 要素内力ベクトルを求める
		for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
			eqf_in_c[ino_c] = 0.0;
			for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
				eqf_in_c[ino_c] += eCmat_cc[ino_c][jno_c]*(val_c[jno_c]+dt*vval_c[jno_c])
					 + eMmat_cc[ino_c][jno_c]*vval_c[jno_c];
			}
			eqf_in_c[ino_c] += eCmat_cb[ino_c]*(val_b+dt*vval_b) + eMmat_cb[ino_c]*vval_b;
		}
		eqf_in_b = 0.0;
		for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
			eqf_in_b += eCmat_bc[jno_c]*(val_c[jno_c]+dt*vval_c[jno_c]) + eMmat_bc[jno_c]*vval_c[jno_c];
		}
		eqf_in_b += eCmat_bb*(val_b+dt*vval_b) + eMmat_bb*vval_b;

		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno_c;i++){
				for(unsigned int j=0;j<nno_c;j++){
					emat_cc[i][j] = eMmat_cc[i][j]+dtmp1*eCmat_cc[i][j];
				}
				emat_cb[i] = eMmat_cb[i]+dtmp1*eCmat_cb[i];
				emat_bc[i] = eMmat_bc[i]+dtmp1*eCmat_bc[i];
			}
			emat_bb = eMmat_bb+dtmp1*eCmat_bb;
		}
		////////////////////////////////////////////////////////////////

		// 要素残差ベクトルを求める
        for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
			eres_c[ino_c] = eqf_out_c[ino_c] - eqf_in_c[ino_c];
		}
		eres_b = eqf_out_b - eqf_in_b;
		// 要素剛性行列のマージ
		mat_cc.Mearge(nno_c,no_c,nno_c,no_c,	1,&emat_cc[0][0]);
		mat_cb.Mearge(nno_c,no_c,nno_b,&no_b,	1,&emat_cb[0]   );
		mat_bc.Mearge(nno_b,&no_b,nno_c,no_c,	1,&emat_bc[0]   );
		mat_bb.Mearge(nno_b,&no_b,nno_b,&no_b,	1,&emat_bb      );
		// 要素残差ベクトルのマージ
		for(unsigned int inoes=0;inoes<nno_c;inoes++){
			res_c.AddValue( no_c[inoes],0,eres_c[inoes]);
		}
		res_b.AddValue( no_b,0,eres_b );
	}
	return true;
}


static bool AddLinearSystem_Diffusion2D_P1(
		double rho, double alpha, double source,
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea)
{
//	std::cout << "Diffusion2D Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	double val_c[nno];		// 要素節点の値
	double vval_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
				
	double emat[nno][nno];
	double eCmat[nno][nno];	// 要素剛性行列
	double eMmat[nno][nno];	// 要素剛性行列
	double eres_c[nno];	// 残差ベクトル
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);//na_c_val.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_vval = field_val.GetNodeSeg(CORNER,true,world,VELOCITY);//na_c_vval.GetSeg(id_ns_c_vval);
	const CNodeAry::CNodeSeg& ns_c_co = field_val.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_co.GetNodes(ielem,no_c);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
		}
		es_c_va.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no_c[inoes],&val_c[inoes]);
			ns_c_vval.GetValue(no_c[inoes],&vval_c[inoes]);
		}

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = rho*area*0.08333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno] = dtmp1;
				}
				eMmat[ino][ino] += dtmp1;
			}
		}
		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = source*area*0.333333333333333333;
		}

		////////////////////////////////////////////////////////////////

		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno;i++){ 
				(&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i]; 
			}
		}
		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_c[ino]	-= eCmat[ino][jno]*(val_c[jno]+dt*vval_c[jno])
				         + eMmat[ino][jno]*vval_c[jno];
		}
		}
		// 要素剛性行列にマージする
		mat_cc.Mearge(nno,no_c,nno,no_c,1,&emat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno;inoes++){
			res_c.AddValue( no_c[inoes],0,eres_c[inoes]);
		}
	}
	return true;
}


bool Fem::Eqn::AddLinSys_Diffusion(
		double dt, double gamma, 
        Fem::Ls::CLinearSystem_Field& ls,
		double rho, double alpha, double source,
		const CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != SCALAR ) return false;
	if( id_ea != 0 ){
		INTERPOLATION_TYPE intp_type = field_val.GetInterpolationType(id_ea,world);
		if( intp_type == TRI11 ){
			AddLinearSystem_Diffusion2D_P1(
				rho,alpha,source,
				gamma, dt,
				ls, id_field_val, world,
				id_ea);
		}
		else if( intp_type == TRI1011 ){
			AddLinearSystem_Diffusion_P1b(
				rho,alpha,source,
				gamma, dt,
				ls, id_field_val, world,
				id_ea);
		}
		else{
			std::cout << "Error!-->NotImplimented" << std::endl;
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_Diffusion(
					dt, gamma, 
					ls,
					rho, alpha, source,
					world, id_field_val,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

static bool AddLinearSystem_Diffusion2D_AxSym_P1(
		double rho, double alpha, double source,
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea)
{
//	std::cout << "Diffusion2D Axial Symmetry Tri P1" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_vval = field_val.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_co = field_val.GetNodeSeg(CORNER,false,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no[nno];	// 要素節点の全体節点番号
		es_c_co.GetNodes(ielem,no);
		// 座標を取り出す
		double coord[nno][ndim];	// 要素節点の座標
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(no[ino],coord[ino]);
		}
		es_c_va.GetNodes(ielem,no);
		// 節点の値を取って来る
		double val_c[nno];		// 要素節点の値
		double vval_c[nno];		// 要素節点の値
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no[inoes],&val_c[inoes]);
			ns_c_vval.GetValue(no[inoes],&vval_c[inoes]);
		}

		const double rad[3] = {
			fabs( coord[0][0] ),
			fabs( coord[1][0] ),
			fabs( coord[2][0] )
		};
		const double ave_rad = (rad[0]+rad[1]+rad[2])*0.33333333333333333333;

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double area = TriArea(coord[0],coord[1],coord[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord[0],coord[1],coord[2]);
		// 要素剛性行列を作る
		double eCmat[nno][nno];	// 要素剛性行列
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat[ino][jno] = alpha*area*ave_rad*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		double eMmat[nno][nno];	// 要素剛性行列
		{
			const double dtmp1 = rho*area/60.0;
			eMmat[0][0] = dtmp1*(6*rad[0] + 2*rad[1] + 2*rad[2]);
			eMmat[1][1] = dtmp1*(2*rad[0] + 6*rad[1] + 2*rad[2]);
			eMmat[2][2] = dtmp1*(2*rad[0] + 2*rad[1] + 6*rad[2]);
			
			eMmat[0][1] = dtmp1*(2*rad[0] + 2*rad[1] + 1*rad[2]);
			eMmat[1][0] = eMmat[0][1];
			eMmat[0][2] = dtmp1*(2*rad[0] + 1*rad[1] + 2*rad[2]);
			eMmat[2][0] = eMmat[0][2];
			eMmat[1][2] = dtmp1*(1*rad[0] + 2*rad[1] + 2*rad[2]);
			eMmat[2][1] = eMmat[1][2];
		}
		double eres_c[nno];	// 残差ベクトル
		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = source*area*0.333333333333333333;
		}

		////////////////////////////////////////////////////////////////

		double emat[nno][nno];
		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno;i++){ 
				(&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i]; 
			}
		}
		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_c[ino]	-= eCmat[ino][jno]*(val_c[jno]+dt*vval_c[jno])
				         + eMmat[ino][jno]*vval_c[jno];
		}
		}
		// 要素剛性行列にマージする
		mat_cc.Mearge(nno,no,nno,no,1,&emat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			res_c.AddValue( no[ino],0,eres_c[ino]);
		}
	}
	return true;
}



bool Fem::Eqn::AddLinSys_Diffusion_AxSym(
		double dt, double gamma, 
		CLinearSystem_Field& ls,
		double rho, double alpha, double source,
		const CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != SCALAR ) return false;
	if( id_ea != 0 ){
		INTERPOLATION_TYPE intp_type = field_val.GetInterpolationType(id_ea,world);
		if( intp_type == TRI11 ){
			AddLinearSystem_Diffusion2D_AxSym_P1(
				rho,alpha,source,
				gamma, dt,
				ls, id_field_val, world,
				id_ea);
		}
		else{
			std::cout << "Error!-->NotImplimented" << std::endl;
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_Diffusion(
					dt, gamma, 
					ls,
					rho, alpha, source,
					world, id_field_val,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


static bool AddLinearSystem_Diffusion2D_P1(
		double rho, double alpha, double source,
		CLinearSystem_SaveDiaM_Newmark& ls, 
		unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea)
{
//	std::cout << "Diffusion2D Tri P1 savemat " << gamma << " " << dt << std::endl;

	const double gamma = ls.GetGamma();
	const double dt = ls.GetDt();

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_val = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co  = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	double coord_c[nno][ndim];	// 要素節点の座標
				
	double emat[nno][nno];
	double eCmat[nno][nno];	// 要素剛性行列
	double eMmat[nno];	// 要素剛性行列
	double eqf_out_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc   = ls.GetMatrix(id_field_val,CORNER,world);
	CVector_Blk&    force_c  = ls.GetForce( id_field_val,CORNER,world);
	
	CMat_BlkCrs& mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,  id_field_val,CORNER,  world);
	const CNodeAry::CNodeSeg& ns_c_co   = field_val.GetNodeSeg(CORNER,false,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_co.GetNodes(ielem,no_c);
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(no_c[ino],coord_c[ino]);
		}
		es_c_val.GetNodes(ielem,no_c);

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat[ino][jno] = alpha*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = rho*area/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eMmat[ino] = dtmp1;
			}
		}
		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out_c[ino] = source*area*0.33333333333333;
		}

		////////////////////////////////////////////////////////////////

		{	// 要素係数行列を求める
			const double dtmp1 = gamma*dt;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					emat[ino][jno] = dtmp1*eCmat[ino][jno];
				}
				emat[ino][ino] += eMmat[ino];
			}
		}

		// 剛性行列にマージする
		mat_cc      .Mearge(nno,no_c,nno,no_c,  1, &emat[0][0]);
		mat_cc_bound.Mearge(nno,no_c,nno,no_c,  1,&eCmat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			force_c.AddValue( no_c[ino],0,eqf_out_c[ino]);
		}
	}
	return true;
}

static bool AddLinearSystem_Diffusion3D_P1(
		double rho, double alpha, double source,
		CLinearSystem_SaveDiaM_Newmark& ls, 
		unsigned int id_field_val, const CFieldWorld& world,
		unsigned int id_ea )
{
//	std::cout << "Diffusion3D Tet savemat" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const double gamma = ls.GetGamma();
	const double dt = ls.GetDt();

	const CElemAry::CElemSeg& es_c = field_val.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nno = 4;
	const unsigned int ndim = 3;

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	double val_c[nno];		// 要素節点の値
	double vval_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
				
	double emat[nno][nno];
	double eCmat[nno][nno];	// 要素剛性行列
	double eMmat[nno][nno];	// 要素剛性行列
	double eqf_out_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc  = ls.GetMatrix(id_field_val,CORNER,world);
	CVector_Blk&    force_c = ls.GetForce( id_field_val,CORNER,world);
	
	CMat_BlkCrs& mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,  id_field_val,CORNER,  world);
	CDiaMat_Blk& mat_mass     = ls.GetDiaMassMatrix(  id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world,VALUE);	
	const CNodeAry::CNodeSeg& ns_c_vval = field_val.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_co = field_val.GetNodeSeg(CORNER,false,world,VALUE);	

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(no_c[ino],coord_c[ino]);
			ns_c_val.GetValue(no_c[ino],&val_c[ino]);
			ns_c_vval.GetValue(no_c[ino],&vval_c[ino]);
		}

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double vol = TetVolume(coord_c[0],coord_c[1],coord_c[2],coord_c[3]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TetDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2],coord_c[3]);
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat[ino][jno] = alpha*vol*( dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]+dldx[ino][2]*dldx[jno][2]);
		}
		}
		{
			const double dtmp1 = rho*vol*0.05;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno] = dtmp1;
				}
				eMmat[ino][ino] += dtmp1;
			}
		}
		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out_c[ino] = source*vol*0.25;
		}

		////////////////////////////////////////////////////////////////

		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
			}
		}

		// 剛性行列にマージする
		mat_cc      .Mearge(nno,no_c,nno,no_c,  1, &emat[0][0]);
		mat_cc_bound.Mearge(nno,no_c,nno,no_c,  1,&eCmat[0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			mat_mass.Mearge( no_c[ino], 1, &eMmat[0][0] );
        }

		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			force_c.AddValue( no_c[ino],0,eqf_out_c[ino]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_Diffusion(
		CLinearSystem_SaveDiaM_Newmark& ls,
		double rho, double alpha, double source,
		const CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);
	if( field_val.GetFieldType() != SCALAR ) return false;

	if( id_ea != 0 ){
		INTERPOLATION_TYPE intp_type = field_val.GetInterpolationType(id_ea,world);
		if( intp_type == TRI11 ){
			AddLinearSystem_Diffusion2D_P1(
				rho,alpha,source,
				ls, id_field_val, world,
				id_ea);
		}
		else if( intp_type == TET11 ){
			AddLinearSystem_Diffusion3D_P1(
				rho,alpha,source,
				ls, id_field_val, world,
				id_ea);
		}
		else{
			std::cout << "NotImplimented" << std::endl;
			getchar();
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res =  Fem::Eqn::AddLinSys_Diffusion(
					ls,
					rho, alpha, source,
					world,
					id_field_val,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}





////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/*

#include "../Mat/Preconditioner.h"
#include "../Mat/Solver_CG.h"


static bool AddLinearSystem_ReactionDiffusion2D_P1(
		double diff1, double diff2, 
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		unsigned int id_val1, unsigned int id_val2, const CFieldWorld& world, 
		int* tmp_buffer, const CField::CElemInterpolation& ei, bool IsInitial )
{
	std::cout << "ReactionDiffusion2D Triangle 3-point 1st order" << std::endl;

	if( !world.IsIdField(id_val1) ) return false;
	const CField& field1 = world.GetField(id_val1);

	if( !world.IsIdField(id_val2) ) return false;
	const CField& field2 = world.GetField(id_val2);

	// 角節点の節点配列ID
	unsigned int id_na_c_co=0,id_ns_c_co=0,  id_na_c_val1=0, id_ns_c_val1=0,id_ns_c_velo1=0;
	{
		unsigned int id_ns_c_acc=0;
		field1.GetID_NodeSeg_CoordValue(CORNER, 
			id_na_c_co,id_ns_c_co, 
			id_na_c_val1, 
			id_ns_c_val1,id_ns_c_velo1,id_ns_c_acc );
		assert( id_na_c_co!=0 && id_ns_c_co!=0 );
		assert( id_na_c_val1!=0 && id_ns_c_val1!=0 && id_ns_c_velo1!=0 );
	}
//	std::cout << id_val1 << " " << id_na_c_val1 << " " << id_ns_c_val1 << " " << id_ns_c_velo1 << std::endl;
	// 角節点の節点配列ID
	unsigned int id_na_c_val2=0, id_ns_c_val2=0,id_ns_c_velo2=0;
	{
		unsigned int id_na_c_co2=0,id_ns_c_co2=0;
		unsigned int id_ns_c_acc=0;
		field2.GetID_NodeSeg_CoordValue(CORNER, 
			id_na_c_co2,id_ns_c_co2, 
			id_na_c_val2, 
			id_ns_c_val2,id_ns_c_velo2,id_ns_c_acc );
		assert( id_na_c_co2==id_na_c_co && id_ns_c_co2==id_ns_c_co );
		assert( id_na_c_val2!=0 && id_ns_c_val2!=0 && id_ns_c_velo2!=0 );
	}
//	std::cout << id_val2 << " " << id_na_c_val2 << " " << id_ns_c_val2 << " " << id_ns_c_velo2 << std::endl;

	const unsigned int id_ea = ei.id_ea;
	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	unsigned int id_es_c = ei.id_es_c_va;	// 角節点の要素セグメントID

	const unsigned int nno_c = 3;
	const unsigned int ndim = 2;

	assert( id_es_c != 0 );

	unsigned int no_c[nno_c];	// 要素節点の全体節点番号
	int buf_no_c[nno_c];

	double val_c1[nno_c], velo_c1[nno_c];		// 要素節点の値
	double val_c2[nno_c], velo_c2[nno_c];		// 要素節点の値
	double coord_c[nno_c][ndim];	// 要素節点の座標
				
	double dldx[nno_c][ndim];	// 形状関数のxy微分
	double const_term[nno_c];	// 形状関数の定数項

	double emat_c1c1[nno_c][nno_c];
	double emat_c2c2[nno_c][nno_c];
	double eCmat[nno_c][nno_c];	// 要素剛性行列
	double eMmat[nno_c][nno_c];	// 要素剛性行列
	double eqf_out_c1[nno_c], eqf_in_c1[nno_c], eres_c1[nno_c];	// 要素節点等価内力、外力、残差ベクトル
	double eqf_out_c2[nno_c], eqf_in_c2[nno_c], eres_c2[nno_c];	// 要素節点等価内力、外力、残差ベクトル

	CMatDia_BlkCrs* mat_c1c1 = ls.GetMatrix(  id_val1,CORNER,world); assert( mat_c1c1!=0 ); // 要素剛性行列(コーナ-コーナー)				
	CMatDia_BlkCrs* mat_c2c2 = ls.GetMatrix(  id_val2,CORNER,world); assert( mat_c2c2!=0 ); // 要素剛性行列(コーナ-コーナー)
	CVector_Blk*    res_c1   = ls.GetResidual(id_val1,CORNER,world); assert( res_c1!=0  ); // 要素残差ベクトル(コーナー)
	CVector_Blk*    res_c2   = ls.GetResidual(id_val2,CORNER,world); assert( res_c2!=0  ); // 要素残差ベクトル(コーナー)

	const CNodeAry& na_c_val1 = world.GetNA(id_na_c_val1);
	const CNodeAry& na_c_val2 = world.GetNA(id_na_c_val2);
	const CNodeAry& na_c_co  = world.GetNA(id_na_c_co );
	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		{	// 要素配列から要素セグメントの節点番号を取り出す
			unsigned int len;
			for(unsigned int ilen=0;ilen<len;ilen++){ 
				no_c[ilen] = (unsigned int)buf_no_c[ilen];
			}
		}
		{	// 節点の値を取って来る
			double val;
			for(unsigned int inoes=0;inoes<nno_c;inoes++){
				for(unsigned int idim=0;idim<ndim;idim++){
					na_c_co.GetValueFromNode(no_c[inoes],id_ns_c_co,idim,val);
					coord_c[inoes][idim] = val;
				}
			}
			for(unsigned int inoes=0;inoes<nno_c;inoes++){
				na_c_val1.GetValueFromNode(no_c[inoes],id_ns_c_val1 ,0,val);
				val_c1[ inoes] = val;
				na_c_val1.GetValueFromNode(no_c[inoes],id_ns_c_velo1,0,val);
				velo_c1[inoes] = val;
			}
			for(unsigned int inoes=0;inoes<nno_c;inoes++){
				na_c_val2.GetValueFromNode(no_c[inoes],id_ns_c_val2 ,0,val);
				val_c2[ inoes] = val;
				na_c_val2.GetValueFromNode(no_c[inoes],id_ns_c_velo2,0,val);
				velo_c2[inoes] = val;
			}
		}

//		std::cout << id_ns_c_val1 << " " << val_c1[0] << " " << val_c1[1] << " " << val_c1[2] << std::endl;
//		std::cout << id_ns_c_val2 << " " << val_c2[0] << " " << val_c2[1] << " " << val_c2[2] << std::endl;

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);
		// 要素剛性行列を作る
		for(int ino_c=0;ino_c<nno_c;ino_c++){
		for(int jno_c=0;jno_c<nno_c;jno_c++){
			eCmat[ino_c][jno_c] 
				= area*(dldx[ino_c][0]*dldx[jno_c][0]+dldx[ino_c][1]*dldx[jno_c][1]);
		}
		}
		{
			const double dtmp1 = area/12.0;
			for(int ino_c=0;ino_c<nno_c;ino_c++){
			for(int jno_c=0;jno_c<nno_c;jno_c++){
				eMmat[ino_c][jno_c] = dtmp1;
			}
			eMmat[ino_c][ino_c] += dtmp1;
			}
		}

		{
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_out_c1[ino_c] = 100.0*( val_c1[ino_c]*(1.0-val_c1[ino_c])*(val_c1[ino_c]-0.2)-val_c2[ino_c] )*area/3.0;
			}
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_out_c2[ino_c] = 25.0/8.0*(8.0/25.0*(val_c1[ino_c]-0.4)+6.0/125.0 - val_c2[ino_c])*area/3.0;
			}
		}

		////////////////////////////////////////////////////////////////

		if( IsInitial ){	// 要素節点等価内力ベクトルを求める
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_in_c1[ino_c] = 0.0;
				for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
					eqf_in_c1[ino_c]
						+= diff1*eCmat[ino_c][jno_c]*(val_c1[jno_c]+dt*velo_c1[jno_c])
						 + eMmat[ino_c][jno_c]*velo_c1[jno_c];
				}
			}
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_in_c2[ino_c] = 0.0;
				for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
					eqf_in_c2[ino_c]
						+= diff2*eCmat[ino_c][jno_c]*(val_c2[jno_c]+dt*velo_c2[jno_c])
						 + eMmat[ino_c][jno_c]*velo_c2[jno_c];
				}
			}
		}
		else{
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_in_c1[ino_c] = 0.0;
				for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
					eqf_in_c1[ino_c]
						+= diff1*eCmat[ino_c][jno_c]*val_c1[jno_c] + eMmat[ino_c][jno_c]*velo_c1[jno_c];
				}
			}
			for(unsigned int ino_c=0;ino_c<nno_c;ino_c++){
				eqf_in_c2[ino_c] = 0.0;
				for(unsigned int jno_c=0;jno_c<nno_c;jno_c++){
					eqf_in_c2[ino_c]
						+= diff2*eCmat[ino_c][jno_c]*val_c2[jno_c] + eMmat[ino_c][jno_c]*velo_c2[jno_c];
				}
			}
		}

		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno_c;i++){
			for(unsigned int j=0;j<nno_c;j++){
				emat_c1c1[i][j] = eMmat[i][j]+dtmp1*diff1*eCmat[i][j];
			}
			}
			for(unsigned int i=0;i<nno_c;i++){	
			for(unsigned int j=0;j<nno_c;j++){
				emat_c2c2[i][j] = eMmat[i][j]+dtmp1*diff2*eCmat[i][j];
			}
			}
		}
		// 要素節点等価残差ベクトルを求める
		for(int ino_c=0;ino_c<nno_c;ino_c++){
			eres_c1[ino_c] = eqf_out_c1[ino_c] - eqf_in_c1[ino_c];
		}
		for(int ino_c=0;ino_c<nno_c;ino_c++){
			eres_c2[ino_c] = eqf_out_c2[ino_c] - eqf_in_c2[ino_c];
		}
		// 要素剛性行列にマージする
		mat_c1c1->Mearge(nno_c,no_c,nno_c,no_c,1,&emat_c1c1[0][0],tmp_buffer);
		mat_c2c2->Mearge(nno_c,no_c,nno_c,no_c,1,&emat_c2c2[0][0],tmp_buffer);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno_c;inoes++){
			res_c1->AddValue( no_c[inoes],0,eres_c1[inoes]);
		}
		for(unsigned int inoes=0;inoes<nno_c;inoes++){
			res_c2->AddValue( no_c[inoes],0,eres_c2[inoes]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinearSystem_ReactionDiffusion(
		double m_diffuse_a, double m_diffuse_b,
		double gamma, double dt,
		CLinearSystem_Field& ls,
		unsigned int id_val_a, unsigned int id_val_b, const CFieldWorld& world, bool IsInitial )
{
	if( !world.IsIdField(id_val_a) ) return false;
	const CField& field_a = world.GetField(id_val_a);
	if( field_a.GetFieldType() != SCALAR ) return false;

	const std::vector<CField::CElemInterpolation>& aEI = field_a.GetElemInterpAry();

	int* tmp_buffer;
	{
		const unsigned int ntmp = ls.GetTmpBufferSize();
		tmp_buffer= new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){
			tmp_buffer[itmp] = -1;
		}
	}
	for(unsigned int iei=0;iei<aEI.size();iei++){
		INTERPOLATION_TYPE intp_type = field_a.GetInterpolationType(iei,world);
		if( intp_type == TRI_100_100 ){
			AddLinearSystem_ReactionDiffusion2D_P1(
				m_diffuse_a,m_diffuse_b,
				gamma, dt,
				ls, 
				id_val_a,id_val_b, world,
				tmp_buffer,aEI[iei], IsInitial );
		}
		else{
			assert(0);
		}
	}
	delete[] tmp_buffer;
	return true;
}

CEqn_ReactionDiffusion::CEqn_ReactionDiffusion() 
: m_diffuse_a(0.25), m_diffuse_b(1.0),
	m_gamma_newmark(1.0), m_dt(0.01)
{
}

bool CEqn_ReactionDiffusion::SetDomain(Fem::Field::CFieldWorld& world){
	m_id_a = world.MakeField_AllRegion(SCALAR,VELOCITY|VALUE,CORNER);
	m_id_b = world.MakeField_AllRegion(SCALAR,VELOCITY|VALUE,CORNER);
	this->ClearLinearSystem();
	return true;
}

void CEqn_ReactionDiffusion::ClearLinearSystem()
{
	if( pLS != 0 ){ delete pLS; pLS=0; }
	if( pPrec != 0 ){ delete pPrec; pPrec=0; }
}

bool CEqn_ReactionDiffusion::InitializeLinearSystem(Fem::Field::CFieldWorld& world)
{	
	if( !world.IsIdField(m_id_a) ) return false;
	if( !world.IsIdField(m_id_b) ) return false;
	// 連立一次方程式クラスの設定
	pLS = new CLinearSystem;
	pLS->AddPattern(m_id_a,world);	// val_fieldからできる全体剛性行列を追加する
	pLS->AddPattern(m_id_b,world);	// val_fieldからできる全体剛性行列を追加する
	// 前処理クラスの作成
	pPrec = new CPreconditioner_ILU;
	pPrec->SetLinearSystem(*pLS);
	return true;
}

double CEqn_ReactionDiffusion::MakeLinearSystem(Fem::Field::CFieldWorld& world, bool IsInitial )
{	
	if( pLS==0 || pPrec==0 ){	
		if( !this->InitializeLinearSystem(world) ) return -1.0;
	}	
	// 連立一次方程式を作る
	pLS->InitializeMarge();	// 連立一次方程式を初期化する
	Fem::Eqn::AddLinearSystem_ReactionDiffusion(	// 全体剛性行列に拡散方程式を足し合わせる
		m_diffuse_a, m_diffuse_b,
		m_gamma_newmark,m_dt,
		*pLS,
		m_id_a,m_id_b,world,
		IsInitial );
	const double norm_res = pLS->FinalizeMarge();
	// 前処理行列を作る
	pPrec->SetValue(*pLS);
	return norm_res;
}


bool CEqn_ReactionDiffusion::Solve(Fem::Field::CFieldWorld& world)
{
	for(unsigned int i=0;i<1;i++)
	{
		double norm_res;
		if( i == 0 ){ norm_res = MakeLinearSystem(world,true); }
		else{ norm_res = MakeLinearSystem(world,false); }
		if( norm_res < 0.0 ) return false;
		std::cout << " Norm Res : " << norm_res << std::endl;

		// Solve Matrix
		{
			double conv_ratio = 1.0e-6;
			unsigned int max_iter = 1000;
			// Solve with Preconditioned Conjugate Gradient
			Fem::Sol::Solve_PCG(conv_ratio,max_iter,*pLS,*pPrec);	
			// Solve with Conjugate Gradient
		//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
			std::cout << max_iter << " " << conv_ratio << std::endl;
		}

		// Update Solution
		if( i== 0 ){
			pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_a,world,VELOCITY,true);
			pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_b,world,VELOCITY,true);
		}
		else{
			pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_a,world,VELOCITY,false);
			pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_b,world,VELOCITY,false);
		}
	}
	return true;
}
*/
