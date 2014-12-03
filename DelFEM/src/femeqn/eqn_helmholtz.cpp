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
// eqn_helmholtz.cpp : ヘルムホルツ方程式の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif


#include <math.h>

#include "delfem/field_world.h"

#include "delfem/matvec/zmatdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/zvector_blk.h"
#include "delfem/matvec/bcflag_blk.h"

#include "delfem/femeqn/eqn_helmholtz.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_hex.h"

#include "delfem/femls/zlinearsystem.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

static bool AddLinearSystem_Helmholtz2D_P1(
		CZLinearSystem& ls, 
		double wave_length,
		const unsigned int id_field_val, const CFieldWorld& world, int* tmp_buffer, 
		const unsigned int id_ea )
{

	std::cout << "Helmholtz2D Triangle 3-point 1st order" << std::endl;

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

	Com::Complex value_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
				
	Com::Complex emat[nno][nno];	// 要素剛性行列
	Com::Complex eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CZMatDia_BlkCrs* mat_cc = ls.GetMatrixPtr(id_field_val,CORNER,world);	assert( mat_cc!=0 );	// 要素剛性行列(コーナ-コーナー)
	CZVector_Blk* res_c = ls.GetResidualPtr(id_field_val,CORNER,world); assert( res_c!=0 );		// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_co.GetNodes(ielem,no_c);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
		}
		// 節点の値を取って来る
		es_c_va.GetNodes(ielem,no_c);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no_c[inoes],&value_c[inoes]);
		}
//		std::cout << "VAL: " << value_c[0] << " " << value_c[1] << " " << value_c[2] << std::endl;

		// 面積を求める
		const double area = TriArea(coord_c[0],coord_c[1],coord_c[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord_c[0],coord_c[1],coord_c[2]);
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno] = area*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{
			double k = 2*3.1416/wave_length;
			double tmp_val = k*k*area/12.0;
			for(unsigned int ino=0;ino<nno;ino++){
				emat[ino][ino] -= tmp_val;
				for(unsigned int jno=0;jno<nno;jno++){
					emat[ino][jno] -= tmp_val;
				}
			}
		}
		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_c[ino] -= emat[ino][jno]*value_c[jno];
			}	
		}
		// 要素剛性行列にマージする
		mat_cc->Mearge(nno,no_c,nno,no_c,1,&emat[0][0],tmp_buffer);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno;inoes++){
			res_c->AddValue( no_c[inoes],0,eres_c[inoes]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_Helmholtz(
		CZLinearSystem& ls,
		double wave_length,
		const CFieldWorld& world,
		const unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( val_field.GetFieldType() != ZSCALAR ) return false;

	if( id_ea != 0 ){
		const unsigned int ntmp = ls.GetTmpBufferSize();
		int* tmp_buffer = new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){ tmp_buffer[itmp] = -1; }
		if( val_field.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinearSystem_Helmholtz2D_P1(ls,
				wave_length,
				id_field_val,world,tmp_buffer,id_ea);
		}
		delete[] tmp_buffer;
	}
	else{
		const std::vector<unsigned int> aIdEA = val_field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_Helmholtz(
					ls,
					wave_length,
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
////////////////////////////////////////////////////////////////

static bool AddLinearSystem_Helmholtz2D_AxalSym_P1(
		CZLinearSystem& ls, 
		double wave_length,
		const unsigned int id_field_val, const CFieldWorld& world, int* tmp_buffer, 
		const unsigned int id_ea )
{

	std::cout << "Helmholtz2D Axal_Symmetry Triangle(1st)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;


	Com::Complex emat[nno][nno];	// 要素剛性行列
	Com::Complex eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CZMatDia_BlkCrs* mat_cc = ls.GetMatrixPtr(id_field_val,CORNER,world);	assert( mat_cc!=0 );	// 要素剛性行列(コーナ-コーナー)
	CZVector_Blk* res_c = ls.GetResidualPtr(id_field_val,CORNER,world); assert( res_c!=0 );		// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	const double k = 2*3.1416/wave_length;
	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no_c[nno];	// 要素節点の全体節点番号
		es_c_co.GetNodes(ielem,no_c);
		double coord[nno][ndim];	// 要素節点の座標
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord[inoes]);
		}
		// 節点の値を取って来る
		es_c_va.GetNodes(ielem,no_c);
		Com::Complex value_c[nno];		// 要素節点の値	
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no_c[inoes],&value_c[inoes]);
		}
//		std::cout << "VAL: " << value_c[0] << " " << value_c[1] << " " << value_c[2] << std::endl;

		const double rad[3] = {
			fabs( coord[0][0] ),
			fabs( coord[1][0] ),
			fabs( coord[2][0] )
		};
		const double ave_rad = (rad[0]+rad[1]+rad[2])*0.33333333333333333333;

		// 面積を求める
		const double area = TriArea(coord[0],coord[1],coord[2]);
		// 形状関数の微分を求める
		double dldx[nno][ndim];	// 形状関数のxy微分
		double const_term[nno];	// 形状関数の定数項
		TriDlDx(dldx,const_term,coord[0],coord[1],coord[2]);
		// 要素剛性行列を作る
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			emat[ino][jno] = area*ave_rad*(dldx[ino][0]*dldx[jno][0]+dldx[ino][1]*dldx[jno][1]);
		}
		}
		{	
			const double dtmp1 = k*k*area/60.0;
			emat[0][0] -= dtmp1*(6*rad[0] + 2*rad[1] + 2*rad[2]);
			emat[1][1] -= dtmp1*(2*rad[0] + 6*rad[1] + 2*rad[2]);
			emat[2][2] -= dtmp1*(2*rad[0] + 2*rad[1] + 6*rad[2]);
			const double dtmp2 = dtmp1*(2*rad[0] + 2*rad[1] + 1*rad[2]);
			emat[0][1] -= dtmp2; 
			emat[1][0] -= dtmp2;
			const double dtmp3 = dtmp1*(2*rad[0] + 1*rad[1] + 2*rad[2]);
			emat[0][2] -= dtmp3;
			emat[2][0] -= dtmp3;
			const double dtmp4 = dtmp1*(1*rad[0] + 2*rad[1] + 2*rad[2]);
			emat[1][2] -= dtmp4;
			emat[2][1] -= dtmp4;
		}
		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_c[ino] -= emat[ino][jno]*value_c[jno];
			}	
		}
		// 要素剛性行列にマージする
		mat_cc->Mearge(nno,no_c,nno,no_c,1,&emat[0][0],tmp_buffer);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno;inoes++){
			res_c->AddValue( no_c[inoes],0,eres_c[inoes]);
		}
	}
	return true;
}

bool Fem::Eqn::AddLinSys_Helmholtz_AxalSym(
		CZLinearSystem& ls,
		double wave_length,
		const CFieldWorld& world,
		const unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( val_field.GetFieldType() != ZSCALAR ) return false;

	if( id_ea != 0 ){
		const unsigned int ntmp = ls.GetTmpBufferSize();
		int* tmp_buffer = new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){ tmp_buffer[itmp] = -1; }
		if( val_field.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinearSystem_Helmholtz2D_AxalSym_P1(ls,
				wave_length,
				id_field_val,world,tmp_buffer,id_ea);
		}
		delete[] tmp_buffer;
	}
	else{
		const std::vector<unsigned int> aIdEA = val_field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_Helmholtz_AxalSym(
					ls,
					wave_length,
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
////////////////////////////////////////////////////////////////

static bool AddLinearSystem_MassMatrixEigen_AxalSym_P1(
		CZLinearSystem_GeneralEigen& ls, 
		const unsigned int id_field_val, const CFieldWorld& world,
		const unsigned int id_ea )
{

	std::cout << "Helmholtz2D Axal_Symmetry Triangle(1st)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

//	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
    const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;
				
	CDiaMat_Blk* mat_cc = ls.GetDiaMassMatrixPtr(id_field_val,CORNER,world);	assert( mat_cc!=0 );	// 要素剛性行列(コーナ-コーナー)

//	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no[nno];	// 要素節点の全体節点番号
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];	// 要素節点の座標
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no[inoes],coord[inoes]);
		}
		const double rad[3] = {
			fabs( coord[0][0] ),
			fabs( coord[1][0] ),
			fabs( coord[2][0] )
		};
		const double ave_rad = (rad[0]+rad[1]+rad[2])*0.33333333333333333333;

		// 面積を求める
		const double area = TriArea(coord[0],coord[1],coord[2]);
		double emat[nno];	// 要素剛性行列
		{	
			const double dtmp1 = area/3.0*ave_rad;
			emat[0] = dtmp1;
			emat[1] = dtmp1;
			emat[2] = dtmp1;
		}
		// 要素剛性行列にマージする
		mat_cc->Mearge(no[0],1,&emat[0]);
		mat_cc->Mearge(no[1],1,&emat[1]);
		mat_cc->Mearge(no[2],1,&emat[2]);
	}
	return true;
}

bool Fem::Eqn::AddLinSys_MassMatrixEigen_AxalSym(
		CZLinearSystem_GeneralEigen& ls,
		const CFieldWorld& world,
		const unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( val_field.GetFieldType() != ZSCALAR ) return false;

	if( id_ea != 0 ){
		const unsigned int ntmp = ls.GetTmpBufferSize();
		int* tmp_buffer = new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){ tmp_buffer[itmp] = -1; }
		if( val_field.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinearSystem_MassMatrixEigen_AxalSym_P1(
				ls,
				id_field_val,world,
				id_ea);
		}
		delete[] tmp_buffer;
	}
	else{
		const std::vector<unsigned int> aIdEA = val_field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_MassMatrixEigen_AxalSym(
					ls,
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
////////////////////////////////////////////////////////////////


static bool AddLinearSystem_SommerfeltRadiationBC2D_B1(
		CZLinearSystem& ls, 
		double wave_length,
		const unsigned int id_field_val, const CFieldWorld& world, int* tmp_buffer, 
		const unsigned int id_ea )
{

	std::cout << "SommerfeltRadiationBC2D B1 2-point 2nd order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == LINE );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 2;
	const unsigned int ndim = 2;

	unsigned int no_c[nno];	// 要素節点の全体節点番号

	Com::Complex value_c[nno];		// 要素節点の値
	double coord_c[nno][ndim];	// 要素節点の座標
				
	Com::Complex emat[nno][nno];	// 要素剛性行列
	Com::Complex eres_c[nno];	// 要素節点等価内力、外力、残差ベクトル
				
	CZMatDia_BlkCrs* mat_cc = ls.GetMatrixPtr(id_field_val,CORNER,world);	assert( mat_cc!=0 );	// 要素剛性行列(コーナ-コーナー)
	CZVector_Blk* res_c = ls.GetResidualPtr(id_field_val,CORNER,world); assert( res_c!=0 );		// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		es_c_co.GetNodes(ielem,no_c);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no_c[inoes],coord_c[inoes]);
		}
		// 節点の値を取って来る
		es_c_va.GetNodes(ielem,no_c);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no_c[inoes],&value_c[inoes]);
		}
//		std::cout << "VAL: " << value_c[0] << " " << value_c[1] << " " << value_c[2] << std::endl;
		const double elen = sqrt( (coord_c[0][0]-coord_c[1][0])*(coord_c[0][0]-coord_c[1][0]) + (coord_c[0][1]-coord_c[1][1])*(coord_c[0][1]-coord_c[1][1]) );

		{
			const double k = 2*3.1416/wave_length;
			Com::Complex tmp_val1 = (k/6.0*elen)*Com::Complex(0,1);
			Com::Complex tmp_val2 = -1/(2.0*elen*k)*Com::Complex(0,1);
//			Com::Complex tmp_val2 = 0.0;
			emat[0][0] = tmp_val1*2+tmp_val2;
			emat[0][1] = tmp_val1  -tmp_val2;
			emat[1][0] = tmp_val1  -tmp_val2;
			emat[1][1] = tmp_val1*2+tmp_val2;
		}
		// 要素節点等価内力ベクトルを求める
		for(unsigned int ino=0;ino<nno;ino++){
			eres_c[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres_c[ino] -= emat[ino][jno]*value_c[jno];
			}	
		}
		// 要素剛性行列にマージする
		mat_cc->Mearge(nno,no_c,nno,no_c,1,&emat[0][0],tmp_buffer);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno;inoes++){
			res_c->AddValue( no_c[inoes],0,eres_c[inoes]);
		}
	}
	return true;
}

// Helmholtz方程式
bool Fem::Eqn::AddLinSys_SommerfeltRadiationBC(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( val_field.GetFieldType() != ZSCALAR ) return false;

	if( id_ea != 0 ){
		const unsigned int ntmp = ls.GetTmpBufferSize();
		int* tmp_buffer = new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){ tmp_buffer[itmp] = -1; }
		if( val_field.GetInterpolationType(id_ea,world) == LINE11 ){
			AddLinearSystem_SommerfeltRadiationBC2D_B1(ls,
				wave_length,
				id_field_val,world,tmp_buffer,id_ea);
		}
		delete[] tmp_buffer;
	}
	else{
		const std::vector<unsigned int> aIdEA = val_field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_SommerfeltRadiationBC(
					ls,
					wave_length,
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
////////////////////////////////////////////////////////////////


static bool AddLinearSystem_SommerfeltRadiationBC2D_AxalSym_B1(
		CZLinearSystem& ls, 
		double wave_length,
		const unsigned int id_field_val, const CFieldWorld& world, int* tmp_buffer, 
		const unsigned int id_ea )
{

	std::cout << "SommerfeltRadiationBC2D B1 2-point 2nd order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == LINE );

	if( !world.IsIdField(id_field_val) ) return false;
	const CField& field_val = world.GetField(id_field_val);

	const CElemAry::CElemSeg& es_c_va = field_val.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_val.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 2;
	const unsigned int ndim = 2;

	CZMatDia_BlkCrs* mat_cc = ls.GetMatrixPtr(id_field_val,CORNER,world);	assert( mat_cc!=0 );	// 要素剛性行列(コーナ-コーナー)
	CZVector_Blk* res_c = ls.GetResidualPtr(id_field_val,CORNER,world); assert( res_c!=0 );		// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_val.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_val.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素配列から要素セグメントの節点番号を取り出す
		unsigned int no[nno];	// 要素節点の全体節点番号
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];	// 要素節点の座標
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no[inoes],coord[inoes]);
		}
		// 節点の値を取って来る
		es_c_va.GetNodes(ielem,no);
		Com::Complex value[nno];		// 要素節点の値	
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_val.GetValue(no[inoes],&value[inoes]);
		}
//		std::cout << "VAL: " << value_c[0] << " " << value_c[1] << " " << value_c[2] << std::endl;
		const double elen = sqrt( (coord[0][0]-coord[1][0])*(coord[0][0]-coord[1][0]) + (coord[0][1]-coord[1][1])*(coord[0][1]-coord[1][1]) );

		const double ave_len = fabs( (coord[0][0]+coord[1][0])*0.5 );

		Com::Complex emat[nno][nno];	// 要素剛性行列
		{
			const double k = 2*3.1416/wave_length;
			Com::Complex tmp_val1 = (k/6.0*elen)*Com::Complex(0,1)*ave_len;
			Com::Complex tmp_val2 = -1/(2.0*elen*k)*Com::Complex(0,1)*ave_len;
			emat[0][0] = tmp_val1*2+tmp_val2;
			emat[0][1] = tmp_val1  -tmp_val2;
			emat[1][0] = tmp_val1  -tmp_val2;
			emat[1][1] = tmp_val1*2+tmp_val2;
		}

		// 要素節点等価内力ベクトルを求める
		Com::Complex eres[nno];	// 要素節点等価内力、外力、残差ベクトル
		for(unsigned int ino=0;ino<nno;ino++){
			eres[ino] = 0.0;
			for(unsigned int jno=0;jno<nno;jno++){
				eres[ino] -= emat[ino][jno]*value[jno];
			}	
		}
		// 要素剛性行列にマージする
		mat_cc->Mearge(nno,no,nno,no,1,&emat[0][0],tmp_buffer);
		// 残差ベクトルにマージする
		for(unsigned int inoes=0;inoes<nno;inoes++){
			res_c->AddValue( no[inoes],0,eres[inoes]);
		}
	}
	return true;
}

// Helmholtz方程式
bool Fem::Eqn::AddLinSys_SommerfeltRadiationBC_AxalSym(
		Fem::Ls::CZLinearSystem& ls,
		double wave_length,
		const Fem::Field::CFieldWorld& world,
		unsigned int id_field_val,
		unsigned int id_ea )
{
	if( !world.IsIdField(id_field_val) ) return false;
	const CField& val_field = world.GetField(id_field_val);

	if( val_field.GetFieldType() != ZSCALAR ) return false;

	if( id_ea != 0 ){
		const unsigned int ntmp = ls.GetTmpBufferSize();
		int* tmp_buffer = new int [ntmp];
		for(unsigned int itmp=0;itmp<ntmp;itmp++){ tmp_buffer[itmp] = -1; }
		if( val_field.GetInterpolationType(id_ea,world) == LINE11 ){
			AddLinearSystem_SommerfeltRadiationBC2D_AxalSym_B1(ls,
				wave_length,
				id_field_val,world,tmp_buffer,id_ea);
		}
		delete[] tmp_buffer;
	}
	else{
		const std::vector<unsigned int> aIdEA = val_field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_SommerfeltRadiationBC_AxalSym(
					ls,
					wave_length,
					world,
					id_field_val,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}
