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
// Eqn_AdvectionDiffusion.cpp : 移流拡散方程式の要素剛性作成関数の実装
////////////////////////////////////////////////////////////////

#include <math.h>

#include "delfem/field_world.h"

#include "delfem/femeqn/eqn_advection_diffusion.h"

#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/bcflag_blk.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

static bool AddLinSys_AdvectionDiffusion_Static_P1P1(
		double myu, double source, 
		CLinearSystem_Field& ls, 
		const unsigned int id_field_val, const unsigned int id_field_velo, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Advection Diffusion Static 2D Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	// 角節点の節点配列ID
	unsigned int id_na_c_co = val_field.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_ns_c_co = val_field.GetNodeSegInNodeAry(CORNER).id_ns_co;
	unsigned int id_na_c_val = val_field.GetNodeSegInNodeAry(CORNER).id_na_va;
	unsigned int id_ns_c_val = val_field.GetNodeSegInNodeAry(CORNER).id_ns_va;

	unsigned int id_na_c_velo = field_velo.GetNodeSegInNodeAry(CORNER).id_na_va;
	unsigned int id_ns_c_velo = field_velo.GetNodeSegInNodeAry(CORNER).id_ns_ve;
	assert( id_na_c_co != 0 && id_ns_c_co != 0 );
	assert( id_na_c_val != 0 && id_ns_c_val != 0 );
	assert( id_na_c_velo != 0 );
	assert( id_ns_c_velo != 0 );

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	const CElemAry::CElemSeg& es_c_val = field_velo.GetElemSeg(id_ea,CORNER,true,world);

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	double value_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
	double velo_c[nno][ndim];
	
	double emat[nno][nno];	// 要素剛性行列
	double eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry& na_c_val = world.GetNA(id_na_c_val);
	const CNodeAry::CNodeSeg& ns_c_val = na_c_val.GetSeg(id_ns_c_val);
	const CNodeAry& na_c_velo = world.GetNA(id_na_c_velo);
	const CNodeAry::CNodeSeg& ns_c_velo = na_c_velo.GetSeg(id_ns_c_velo);
	const CNodeAry& na_c_co = world.GetNA(id_na_c_co);
	const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_val.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
			ns_c_val.GetValue(no_c[inoes],&value_c[inoes]);
			ns_c_velo.GetValue(no_c[inoes],velo_c[inoes]);
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
			emat[ino][jno] = myu*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = area*0.08333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				const double dtmp_0 = dtmp1*(velo_c[0][0]+velo_c[1][0]+velo_c[2][0]+velo_c[ino][0]);
				const double dtmp_1 = dtmp1*(velo_c[0][1]+velo_c[1][1]+velo_c[2][1]+velo_c[ino][1]);
				for(unsigned int jno=0;jno<nno;jno++){
					emat[ino][jno] += dldx[jno][0]*dtmp_0+dldx[jno][1]*dtmp_1;
				}
			}
		}

		// Calc Stabilization Parameter
		double tau;
		{
			const double velo_ave[2] = { 
				(velo_c[0][0]+velo_c[1][0]+velo_c[2][0])/3.0, 
				(velo_c[0][1]+velo_c[1][1]+velo_c[2][1])/3.0 };
			const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
			const double velo_dir[2] = { velo_ave[0]/norm_v, velo_ave[1]/norm_v };
			// calc element length along the direction of velocity
			double h;
			{
				double dtmp1 = 0;
                for(unsigned int inode=0;inode<nno;inode++){
					dtmp1 += fabs(velo_dir[0]*dldx[inode][0]+velo_dir[1]*dldx[inode][1]);
				}
				h = 2.0/dtmp1;
			}
			// calc stabilization parameter
			if( myu > 1.0e-20 ){
				const double re_c = 0.5*norm_v*h/myu;	// 0.5*norm_v*h*rho/myu;
				if(  re_c < 3.0 ){ tau = h * 0.5 / norm_v * re_c / 3.0; }
				else{ tau = h * 0.5 / norm_v; }
			}
			else{ tau = h * 0.5 / norm_v; }
		}

		{
			double tmp_mat[ndim][ndim];
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				double dtmp1 = 0.0;
				for(unsigned int ino=0;ino<nno;ino++){
					for(unsigned int jno=0;jno<nno;jno++){
						dtmp1 += velo_c[ino][idim]*velo_c[jno][jdim];
					}
					dtmp1 += velo_c[ino][idim]*velo_c[ino][jdim];
				}
				tmp_mat[idim][jdim] = area*tau*dtmp1/12.0;
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dtmp1 += dldx[ino][idim]*dldx[jno][jdim]*tmp_mat[idim][jdim];
				}
				}
				emat[ino][jno] += dtmp1;
			}
			}
		}


		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = source*area*0.3333333333333333333;
		}

		////////////////////////////////////////////////////////////////

		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){			
			for(unsigned int jno=0;jno<nno;jno++){
				eres_c[ino] -= emat[ino][jno]*value_c[jno];
			}	
		}
		// 要素剛性行列にマージする
		mat_cc.Mearge(nno,no_c,nno,no_c,1,&emat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			res_c.AddValue( no_c[ino],0,eres_c[ino]);
		}
	}
	return true;
}


bool Fem::Eqn::AddLinSys_AdvectionDiffusion_Static(
        Fem::Ls::CLinearSystem_Field& ls,
		double myu, double source,
		const CFieldWorld& world,
		const unsigned int id_field_val, unsigned int id_field_velo, 
		unsigned int id_ea )
{
	
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( field_val.GetFieldType() != Fem::Field::SCALAR ) return false;
	if( field_velo.GetFieldType() != VECTOR2 ) return false;


	if( id_ea != 0 ){
		Fem::Field::INTERPOLATION_TYPE intp_type_val  = field_val.GetInterpolationType(id_ea,world);
		Fem::Field::INTERPOLATION_TYPE intp_type_velo = field_velo.GetInterpolationType(id_ea,world);
		if( intp_type_val == TRI11 && intp_type_velo == TRI11 ){
			AddLinSys_AdvectionDiffusion_Static_P1P1(
				myu,source,
				ls,
				id_field_val,id_field_velo,world,
				id_ea);
		}
		else{
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_AdvectionDiffusion_Static(
				ls,
				myu, source,
				world,
				id_field_val, id_field_velo, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}






static bool AddLinSys_AdvectionDiffusion_Static_P1P1(
		double myu, double source, 
		CLinearSystem_Save& ls, 
		unsigned int id_field_val, unsigned int id_field_velo, const CFieldWorld& world, 
		unsigned int id_ea)
{
//	std::cout << "Advection Diffusion Static 2D Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	// 角節点の節点配列ID
	unsigned int id_na_c_co = val_field.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_ns_c_co = val_field.GetNodeSegInNodeAry(CORNER).id_ns_co;
	unsigned int id_na_c_val = val_field.GetNodeSegInNodeAry(CORNER).id_na_va;
	unsigned int id_ns_c_val = val_field.GetNodeSegInNodeAry(CORNER).id_ns_va;

	unsigned int id_na_c_velo = field_velo.GetNodeSegInNodeAry(CORNER).id_na_va;
	unsigned int id_ns_c_velo = field_velo.GetNodeSegInNodeAry(CORNER).id_ns_ve;
	assert( id_na_c_co != 0 && id_ns_c_co != 0 );
	assert( id_na_c_val != 0 && id_ns_c_val != 0 );
	assert( id_na_c_velo != 0 );
	assert( id_ns_c_velo != 0 );

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	const CElemAry::CElemSeg& es_c_val = val_field.GetElemSeg(id_ea,CORNER,true,world);

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	double value_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
	double velo_c[nno][ndim];
				
	double emat[nno][nno];	// 要素剛性行列
	double eqf_out_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);
	
	CMat_BlkCrs& mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,id_field_val,CORNER,world);

	const CNodeAry& na_c_val = world.GetNA(id_na_c_val);
	const CNodeAry::CNodeSeg& ns_c_val = na_c_val.GetSeg(id_ns_c_val);
	const CNodeAry& na_c_velo = world.GetNA(id_na_c_velo);
	const CNodeAry::CNodeSeg& ns_c_velo = na_c_velo.GetSeg(id_ns_c_velo);
	const CNodeAry& na_c_co = world.GetNA(id_na_c_co);
	const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_val.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
			ns_c_val.GetValue(no_c[inoes],&value_c[inoes]);
			ns_c_velo.GetValue(no_c[inoes],velo_c[inoes]);
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
			emat[ino][jno] = myu*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = area*0.083333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				const double dtmp_0 = dtmp1*(velo_c[0][0]+velo_c[1][0]+velo_c[2][0]+velo_c[ino][0]);
				const double dtmp_1 = dtmp1*(velo_c[0][1]+velo_c[1][1]+velo_c[2][1]+velo_c[ino][1]);
				for(unsigned int jno=0;jno<nno;jno++){
					emat[ino][jno] += dldx[jno][0]*dtmp_0+dldx[jno][1]*dtmp_1;
				}
			}
		}

		// Calc Stabilization Parameter
		double tau;
		{
			const double velo_ave[2] = { 
				(velo_c[0][0]+velo_c[1][0]+velo_c[2][0])/3.0, 
				(velo_c[0][1]+velo_c[1][1]+velo_c[2][1])/3.0 };
			const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
			const double velo_dir[2] = { velo_ave[0]/norm_v, velo_ave[1]/norm_v };
			// calc element length along the direction of velocity
			double h;
			{
				double dtmp1 = 0;
                for(unsigned int inode=0;inode<nno;inode++){
					dtmp1 += fabs(velo_dir[0]*dldx[inode][0]+velo_dir[1]*dldx[inode][1]);
				}
				h = 2.0/dtmp1;
			}
			// calc stabilization parameter
			if( myu > 1.0e-20 ){
				const double re_c = 0.5*norm_v*h/myu;	// 0.5*norm_v*h*rho/myu;
				if(  re_c < 3.0 ){ tau = h * 0.5 / norm_v * re_c / 3.0; }
				else{ tau = h * 0.5 / norm_v; }
			}
			else{ tau = h * 0.5 / norm_v; }
		}

		{
			double tmp_mat[ndim][ndim];
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				double dtmp1 = 0.0;
				for(unsigned int ino=0;ino<nno;ino++){
					for(unsigned int jno=0;jno<nno;jno++){
						dtmp1 += velo_c[ino][idim]*velo_c[jno][jdim];
					}
					dtmp1 += velo_c[ino][idim]*velo_c[ino][jdim];
				}
				tmp_mat[idim][jdim] = area*tau*dtmp1/12.0;
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dtmp1 += dldx[ino][idim]*dldx[jno][jdim]*tmp_mat[idim][jdim];
				}
				}
				emat[ino][jno] += dtmp1;
			}
			}
		}
		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out_c[ino] = source*area*0.333333333333333333333;
		}

		////////////////////////////////////////////////////////////////

		// 要素剛性行列にマージする
		mat_cc      .Mearge(nno,no_c,nno,no_c,1,&emat[0][0]);
		mat_cc_bound.Mearge(nno,no_c,nno,no_c,1,&emat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			res_c.AddValue( no_c[ino],0,eqf_out_c[ino]);
		}
	}
	return true;
}




bool Fem::Eqn::AddLinSys_AdvectionDiffusion_Static(
		CLinearSystem_Save& ls,
		double myu, double source,
		const CFieldWorld& world,
		const unsigned int id_field_val, unsigned int id_field_velo, 
		unsigned int id_ea )
{	
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( field_val.GetFieldType() != SCALAR ) return false;
	if( field_velo.GetFieldType() != VECTOR2 ) return false;

	if( id_ea != 0 )
	{
		Fem::Field::INTERPOLATION_TYPE intp_type_val  = field_val.GetInterpolationType(id_ea,world);
		Fem::Field::INTERPOLATION_TYPE intp_type_velo = field_velo.GetInterpolationType(id_ea,world);
		if( intp_type_val == TRI11 && intp_type_velo == TRI11 ){
			return AddLinSys_AdvectionDiffusion_Static_P1P1(
				myu,source,
				ls,
				id_field_val,id_field_velo,world,
				id_ea);
		}
		else{
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰呼び出し
			bool res =  AddLinSys_AdvectionDiffusion_Static(
				ls,
				myu, source,
				world, 
				id_field_val, id_field_velo, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}

static bool AddLinSys_AdvectionDiffusion_NonStatic_Newmark_P1P1(
		double rho, double myu, double source, 
		double gamma, double dt,
		CLinearSystem_Field& ls, 
		const unsigned int id_field_val, const unsigned int id_field_velo, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Advection Diffusion NonStatic 2D Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	const CElemAry::CElemSeg& es_c_val = val_field.GetElemSeg(id_ea,CORNER,true,world);

	double val_c[nno];		// 要素節点の値
	double vval_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
	double velo_c[nno][ndim];
	double eCmat[nno][nno];
	double eMmat[nno][nno];
	double emat[nno][nno];	// 要素剛性行列
	double eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_val,CORNER,world);
	CVector_Blk&    res_c  = ls.GetResidual(id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = val_field.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_vval = val_field.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_velo = field_velo.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_co = field_velo.GetNodeSeg(CORNER,false,world,VALUE);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no_c[nno];	// 要素節点の全体節点番号
		es_c_val.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
			ns_c_val.GetValue(no_c[inoes],&val_c[inoes]);
			ns_c_vval.GetValue(no_c[inoes],&vval_c[inoes]);
			ns_c_velo.GetValue(no_c[inoes],velo_c[inoes]);
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
			eCmat[ino][jno] = myu*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = rho*area*0.0833333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				const double dtmp_0 = dtmp1*(velo_c[0][0]+velo_c[1][0]+velo_c[2][0]+velo_c[ino][0]);
				const double dtmp_1 = dtmp1*(velo_c[0][1]+velo_c[1][1]+velo_c[2][1]+velo_c[ino][1]);
				for(unsigned int jno=0;jno<nno;jno++){
					eCmat[ino][jno] += dldx[jno][0]*dtmp_0+dldx[jno][1]*dtmp_1;
				}
			}
		}

		// Calc Stabilization Parameter
		double tau;
		{
			const double velo_ave[2] = { 
				(velo_c[0][0]+velo_c[1][0]+velo_c[2][0])*0.3333333333333333, 
				(velo_c[0][1]+velo_c[1][1]+velo_c[2][1])*0.3333333333333333 };
			const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
			if( norm_v < 1.0e-10 ){ tau = 0.0; }
			else{
				const double velo_dir[2] = { velo_ave[0]/norm_v, velo_ave[1]/norm_v };
				// calc element length along the direction of velocity
				double h;
				{
					double dtmp1 = 0;
					for(int inode=0;inode<3;inode++){
						dtmp1 += fabs(velo_dir[0]*dldx[inode][0]+velo_dir[1]*dldx[inode][1]);
					}
					h = 2.0/dtmp1;
				}
				// calc stabilization parameter
				if( norm_v*h*rho < 6.0*myu ){
					const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
					tau = h * 0.5 / norm_v * re_c / 3.0;
				}
				else{ tau = h * 0.5 / norm_v; }
				tau *= 0.5;
			}
		}
		{
			double tmp_mat[ndim][ndim];
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				double dtmp1 = 0.0;
				for(unsigned int ino=0;ino<nno;ino++){
					for(unsigned int jno=0;jno<nno;jno++){
						dtmp1 += velo_c[ino][idim]*velo_c[jno][jdim];
					}
					dtmp1 += velo_c[ino][idim]*velo_c[ino][jdim];
				}
				tmp_mat[idim][jdim] = area*tau*dtmp1*0.0833333333333333;
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dtmp1 += dldx[ino][idim]*dldx[jno][jdim]*tmp_mat[idim][jdim];
				}
				}
				eCmat[ino][jno] += dtmp1*rho;
			}
			}
		}

		{
			const double dtmp1 = rho*area*0.083333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno] = dtmp1;
				}
				eMmat[ino][ino] += dtmp1;
			}
		}

		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = source*area*0.333333333333333;
		}

		////////////////////////////////////////////////////////////////

		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eres_c[ino] -= eCmat[ino][jno]*(val_c[jno]+dt*vval_c[jno])
				         + eMmat[ino][jno]*vval_c[jno];
		}
		}
		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
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

bool Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
	double dt, double gamma, 
	Fem::Ls::CLinearSystem_Field& ls,
	double rho, double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( field_val.GetFieldType() != SCALAR ) return false;
	if( field_velo.GetFieldType() != VECTOR2 ) return false;

	if( id_ea != 0 ){
		Fem::Field::INTERPOLATION_TYPE intp_type_val  = field_val.GetInterpolationType(id_ea,world);
//		Fem::Field::INTERPOLATION_TYPE intp_type_velo = field_velo.GetInterpolationType(id_ea,world);
		if( intp_type_val == TRI11 ){
			AddLinSys_AdvectionDiffusion_NonStatic_Newmark_P1P1(
				rho,myu,source,
				gamma,dt,
				ls,
				id_field_val,id_field_velo,world,
				id_ea);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			// 再帰呼び出し
			bool res = Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
				dt, gamma,
				ls,
				rho, myu, source,
				world,
				id_field_val, id_field_velo, 
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


static bool AddLinSys_AdvectionDiffusion_NonStatic_Newmark_P1P1(
		double rho, double myu, double source, 
		CLinearSystem_SaveDiaM_Newmark& ls, 
		const unsigned int id_field_val, const unsigned int id_field_velo, const CFieldWorld& world, 
		unsigned int id_ea )
{
//	std::cout << "Poisson2D Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	const double gamma = ls.GetGamma();
	const double dt = ls.GetDt();

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	const CElemAry::CElemSeg& es_c_val = val_field.GetElemSeg(id_ea,CORNER,true,world);

	double emat[nno][nno];	// 要素剛性行列
	double eqf_out_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CMatDia_BlkCrs& mat_cc  = ls.GetMatrix(id_field_val,CORNER,world);
	CVector_Blk&    force_c = ls.GetForce( id_field_val,CORNER,world);

	CMat_BlkCrs& mat_cc_bound = ls.GetMatrix_Boundary(id_field_val,CORNER,id_field_val,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_val = val_field.GetNodeSeg(CORNER,true,world,VALUE);		//na_c_val.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_vval = val_field.GetNodeSeg(CORNER,true,world,VELOCITY);	//na_c_val.GetSeg(id_ns_c_vval);
	const CNodeAry::CNodeSeg& ns_c_velo = field_velo.GetNodeSeg(CORNER,true,world,VELOCITY);//na_c_velo.GetSeg(id_ns_c_velo);
	const CNodeAry::CNodeSeg& ns_c_co = field_velo.GetNodeSeg(CORNER,false,world,VALUE);	//na_c_co.GetSeg(id_ns_c_co);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no_c[nno];	// 要素節点の全体節点番号
		es_c_val.GetNodes(ielem,no_c);
		// 節点の値を取って来る
		double val_c[nno], vval_c[nno];		// 要素節点の値
		double coord_c[nno][ndim];	// 要素節点の座標
		double velo_c[nno][ndim];	// advection velocity
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_val.GetValue(no_c[ino],&val_c[ino]);
			ns_c_vval.GetValue(no_c[ino],&vval_c[ino]);
			ns_c_velo.GetValue(no_c[ino],velo_c[ino]);
			ns_c_co.GetValue(no_c[ino],coord_c[ino]);
		}

		////////////////////////////////////////////////////////////////

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);

		double eCmat[nno][nno];
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			eCmat[ino][jno] = myu*area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			const double dtmp1 = area/12.0;
			for(unsigned int ino=0;ino<nno;ino++){
				const double dtmp_0 = dtmp1*(velo_c[0][0]+velo_c[1][0]+velo_c[2][0]+velo_c[ino][0]);
				const double dtmp_1 = dtmp1*(velo_c[0][1]+velo_c[1][1]+velo_c[2][1]+velo_c[ino][1]);
                for(unsigned int jno=0;jno<nno;jno++){
					eCmat[ino][jno] += dldx[jno][0]*dtmp_0+dldx[jno][1]*dtmp_1;
				}
			}
		}

		// Calc Stabilization Parameter
		double tau;
		{
			const double velo_ave[2] = { 
				(velo_c[0][0]+velo_c[1][0]+velo_c[2][0])/3.0, 
				(velo_c[0][1]+velo_c[1][1]+velo_c[2][1])/3.0 };
			const double norm_v = sqrt(velo_ave[0]*velo_ave[0]+velo_ave[1]*velo_ave[1]);
			const double velo_dir[2] = { velo_ave[0]/norm_v, velo_ave[1]/norm_v };			
			double h;
			{	// calc element length along the direction of velocity
				double dtmp1 = 0;
				for(int inode=0;inode<3;inode++){
					dtmp1 += fabs(velo_dir[0]*dldx[inode][0]+velo_dir[1]*dldx[inode][1]);
				}
				h = 2.0/dtmp1;
			}
			// calc stabilization parameter
			if( norm_v*h*rho < 6.0*myu ){
				const double re_c = 0.5*norm_v*h*rho/myu;	// 0.5*norm_v*h*rho/myu;
				tau = h * 0.5 / norm_v * re_c / 3.0;
			}
			else{ tau = h * 0.5 / norm_v; }
//			tau *= 0.1;
		}

		{
			double tmp_mat[ndim][ndim];
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				double dtmp1 = 0.0;
				for(unsigned int ino=0;ino<nno;ino++){
					for(unsigned int jno=0;jno<nno;jno++){
						dtmp1 += velo_c[ino][idim]*velo_c[jno][jdim];
					}
					dtmp1 += velo_c[ino][idim]*velo_c[ino][jdim];
				}
				tmp_mat[idim][jdim] = area*tau*dtmp1/12.0;
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				double dtmp1 = 0.0;
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dtmp1 += dldx[ino][idim]*dldx[jno][jdim]*tmp_mat[idim][jdim];
				}
				}
				eCmat[ino][jno] += dtmp1;
			}
			}
		}

		double eMmat[nno][nno];
		{
			const double dtmp1 = rho*area*0.08333333333333333333333333;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat[ino][jno] = dtmp1;
				}
				eMmat[ino][ino] += dtmp1;
			}
		}

		// 要素節点等価外力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eqf_out_c[ino] = source*area*0.333333333333333333333;
		}

		////////////////////////////////////////////////////////////////

		{	// 要素係数行列を求める
			double dtmp1 = gamma*dt;
			for(unsigned int i=0;i<nno*nno;i++){
				(&emat[0][0])[i] = (&eMmat[0][0])[i]+dtmp1*(&eCmat[0][0])[i];
			}
		}
		// 要素剛性行列にマージする
		mat_cc      .Mearge(nno,no_c,nno,no_c,1, &emat[0][0]);
		mat_cc_bound.Mearge(nno,no_c,nno,no_c,1,&eCmat[0][0]);
		// 残差ベクトルにマージする
		for(unsigned int ino=0;ino<nno;ino++){
			force_c.AddValue( no_c[ino],0,eqf_out_c[ino]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
	Fem::Ls::CLinearSystem_SaveDiaM_Newmark& ls,
	double rho, double myu, double source,
	const Fem::Field::CFieldWorld& world,
	unsigned int id_field_val, unsigned int id_field_velo, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	if( !world.IsIdField(id_field_velo) ) return false;
	const CField& field_velo = world.GetField(id_field_velo);

	if( field_val.GetFieldType() != SCALAR ) return false;
	if( field_velo.GetFieldType() != VECTOR2 ) return false;

	if( id_ea != 0 ){
		Fem::Field::INTERPOLATION_TYPE intp_type_val  = field_val.GetInterpolationType(id_ea,world);
		Fem::Field::INTERPOLATION_TYPE intp_type_velo = field_velo.GetInterpolationType(id_ea,world);
		if( intp_type_val == TRI11 && intp_type_velo == TRI11 ){
			AddLinSys_AdvectionDiffusion_NonStatic_Newmark_P1P1(
				rho,myu,source,
				ls,
				id_field_val,id_field_velo,world,
				id_ea);
		}
		else{
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_val.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
				ls,
				rho, myu, source,
				world,
				id_field_val, id_field_velo, 
				id_ea );
			if( !res ) return false;
		}
		return true;
	}



	return true;
}
