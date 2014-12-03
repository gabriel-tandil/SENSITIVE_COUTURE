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
// eqn_st_venant.cpp : St.Venant-Kirchhoff体の要素剛性作成部の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"

#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_hex.h"
#include "delfem/femeqn/eqn_st_venant.h"

#include "delfem/field_world.h"
#include "delfem/field.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;



void AddElemMatFin_StVenant2D( 
		double detwei, double myu, double lambda, 
		const unsigned int nnoes, const double dudx[][2], const double dndx[][2],
		double eKMat[][2][2], double eForce_in[][2] )
{
	const unsigned int ndim = 2;

	double strain[ndim][ndim];
	for(unsigned int idim=0;idim<ndim;idim++){
	for(unsigned int jdim=0;jdim<ndim;jdim++){
		strain[idim][jdim] = 0.5*( dudx[idim][jdim] + dudx[jdim][idim] );
		for(unsigned int kdim=0;kdim<ndim;kdim++){
			strain[idim][jdim] += 0.5*dudx[kdim][idim]*dudx[kdim][jdim];
		}
	}
	}
	assert( nnoes <= 16 );
	double disp2strain[16][ndim][ndim][ndim];
	{
		double z_mat[ndim][ndim];
		for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				z_mat[idim][jdim] = dudx[idim][jdim];
			}
			z_mat[idim][idim] += 1.0;
		}
		for(unsigned int knoes=0;knoes<nnoes;knoes++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
				for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int kdim=0;kdim<ndim;kdim++){
				disp2strain[knoes][kdim][idim][jdim] = dndx[knoes][jdim]*z_mat[kdim][idim];
			}
			}
		}
		}
	}
	double stress[ndim][ndim];
	{	// 応力を求める
		double dtmp1 = 0.0;
		for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				stress[idim][jdim] = 2.0*myu*strain[idim][jdim];
			}
			dtmp1 += strain[idim][idim];
		}
		for(unsigned int idim=0;idim<ndim;idim++){
			stress[idim][idim] += lambda*dtmp1;
		}
	}
	// 内力を求める
	for(unsigned int knoes=0;knoes<nnoes;knoes++){
	for(unsigned int kdim=0;kdim<ndim;kdim++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eForce_in[knoes][kdim]
				+= detwei*stress[idim][jdim]*disp2strain[knoes][kdim][idim][jdim];
		}
		}
	}
	}
	// 剛性行列を求める
	for(unsigned int ino=0;ino<nnoes;ino++){
	for(unsigned int jno=0;jno<nnoes;jno++){
		// 初期剛性行列を求める
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			{
				double dtmp1 = 0.0;
				for(unsigned int gdim=0;gdim<ndim;gdim++){
				for(unsigned int hdim=0;hdim<ndim;hdim++){
					dtmp1 += disp2strain[ino][idim][gdim][hdim]*disp2strain[jno][jdim][gdim][hdim]
						+disp2strain[ino][idim][gdim][hdim]*disp2strain[jno][jdim][hdim][gdim];
				}
				}
				eKMat[ino*nnoes+jno][idim][jdim] += detwei*myu*dtmp1;
			}
			{
				double dtmp1=0.0, dtmp2=0.0;
				for(unsigned int gdim=0;gdim<ndim;gdim++){
					dtmp1+=disp2strain[ino][idim][gdim][gdim];
					dtmp2+=disp2strain[jno][jdim][gdim][gdim];
				}
				eKMat[ino*nnoes+jno][idim][jdim] += detwei*lambda*dtmp1*dtmp2;
			}
		}
		}
		{	// 初期応力行列を求める
			double dtmp1 = 0.0;
			for(unsigned int gdim=0;gdim<ndim;gdim++){
			for(unsigned int hdim=0;hdim<ndim;hdim++){
				dtmp1 += stress[gdim][hdim]*dndx[ino][hdim]*dndx[jno][gdim];
			}
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				eKMat[ino*nnoes+jno][idim][idim] += detwei*dtmp1;
			}
		}
	}
	}
}

////////////////////////////////
////////////////////////////////


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


bool AddLinSys_StVenant2D_Static_P1(
        Fem::Eqn::ILinearSystem_Eqn& ls,
		double lambda, double myu,
		double  rho, double g_x, double g_y,
		const unsigned int id_field_disp, const CFieldWorld& world, 
		const unsigned int id_ea)
{
	
//	std::cout << "St.Venant2D Static Tri3point" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );
	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);

	const unsigned int nnoes = 3;
	const unsigned int ndim = 2;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eforce_ex[nnoes][ndim];		// 要素内外力ベクトル
	double eforce_in[nnoes][ndim];		// 要素内内力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double ecoords[nnoes][ndim];		// 要素節点座標
	double edisp[  nnoes][ndim];		// 要素節点変位

	double dldx[nnoes][ndim];		// 形状関数の空間微分
	double zero_order_term[nnoes];	// 形状関数の定数項
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);	// 要素剛性行列(コーナ-コーナー)
	CVector_Blk&     res_c = ls.GetResidual(id_field_disp,CORNER,world);	// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_disp.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);
/*
	unsigned int num_integral = 0;
	const unsigned int nInt = NIntTriGauss[num_integral];
	const double (*Gauss)[3] = TriGauss[num_integral];
*/
	double g[2] = { g_x, g_y };

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		es_co.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(noes[ino],ecoords[ino]);
		}
		// 要素の節点番号を取ってくる
		es_va.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_val.GetValue(noes[ino],edisp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_ex[0][0] +i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_in[0][0] +i) = 0.0; }

		// 面積を求める
		const double area = TriArea(ecoords[0],ecoords[1],ecoords[2]);

		// 形状関数のｘｙ微分を求める
		TriDlDx(dldx, zero_order_term,   ecoords[0], ecoords[1], ecoords[2]);

        double dudx[ndim][ndim] = { { 0.0, 0.0}, {0.0, 0.0} };
		for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				dudx[idim][jdim] += edisp[ino][idim]*dldx[ino][jdim];
			}
			}
		}

		AddElemMatFin_StVenant2D( area, myu,lambda, nnoes,dudx,dldx,  emat[0],eforce_in );

		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eforce_ex[ino][idim] += area*rho*g[idim]/nnoes;
		}
		}

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = eforce_ex[ino][idim] - eforce_in[ino][idim];
		}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
		for(unsigned int idim=0;idim<ndim;idim++){
			res_c.AddValue(noes[ino],idim,eres[ino][idim]);
		}
		}
	}

	return true;
}


bool Fem::Eqn::AddLinSys_StVenant2D_Static
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y,
 const CFieldWorld& world,
 const unsigned int id_field_disp,
 unsigned int id_ea )
{
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR2 ) return false;
	
	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TRI11 ){
			return AddLinSys_StVenant2D_Static_P1
			(ls,
			 lambda, myu,
			 rho, g_x, g_y,
			 id_field_disp,world,
			 id_ea);
		}
        assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_StVenant2D_Static
			(ls,
			 lambda, myu, rho,  g_x, g_y,
			 world, id_field_disp,
			 id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

///////////////////////////////////////////


bool AddLinSys_StVenant2D_NonStatic_NewmarkBeta_P1
(double gamma, double beta, double dt,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double rho, double g_x, double g_y,
 const unsigned int id_field_disp, const CFieldWorld& world, 
 bool is_initial,
 const unsigned int id_ea )
{
//	std::cout << "St.Venant2D Triangle 3-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_disp.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 3;
	const unsigned int ndim = 2;

	unsigned int noes[nnoes];	// elment node nuber 2 grobal node number
	
	double eKmat[nnoes][nnoes][ndim][ndim];	// element stiffness matrix
	double eMmat[nnoes][nnoes][ndim][ndim];	// element mass matrix		
	double eforce_in[nnoes][ndim];		// element residual vector		
	double eforce_ex[nnoes][ndim];		// element external force vector
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(id_field_disp,CORNER,world);	// matrix
	CVector_Blk&     res_c = ls.GetResidual(id_field_disp,CORNER,world);// residual vector
	
	const CNodeAry::CNodeSeg& ns_c_val  = field_disp.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc  = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);
	assert( ns_c_val.Length() == ndim );
	assert( ns_c_co.Length()  == ndim );

	double g[2] = { g_x, g_y };

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		es_co.GetNodes(ielem,noes);	// get global node number of coordinate
		double ecoords[nnoes][ndim];	// element node coordinate
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino], ecoords[ino]);
		}
		
		es_va.GetNodes(ielem,noes);	// get global node nubmer of value
		double edisp[  nnoes][ndim];		// displacement 
		double evelo[  nnoes][ndim];		// velocity
		double eacc[  nnoes][ndim];			// acceleration
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_val.GetValue( noes[ino], edisp[ino]);
			ns_c_velo.GetValue(noes[ino], evelo[ino]);
			ns_c_acc.GetValue( noes[ino], eacc[ ino]);
		}

		////////////////////////////////
		
		// set 0 to elemnet stiffness matrix, element residual
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_ex[0][0]  +i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_in[0][0]  +i) = 0.0; }

		const double area = TriArea(ecoords[0],ecoords[1],ecoords[2]);

		double dldx[nnoes][ndim];		// spacial derivative of shape function
		double zero_order_term[nnoes];	// constant term of shape function
		TriDlDx(dldx, zero_order_term,   ecoords[0], ecoords[1], ecoords[2]);
		
		{	// calc stiffness matrix
            double dudx[ndim][ndim] = { { 0.0, 0.0}, {0.0, 0.0} };
			for(unsigned int ino=0;ino<nnoes;ino++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dudx[idim][jdim] += edisp[ino][idim]*dldx[ino][jdim];
				}
				}
			}
			AddElemMatFin_StVenant2D( area, myu,lambda, nnoes,dudx,dldx,  eKmat[0],eforce_in );
		}
		{	// calc mass matrix
			const double tmp1 = rho*area/12.0;
			for(unsigned int ino=0;ino<nnoes;ino++){
				for(unsigned int jno=0;jno<nnoes;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
					eMmat[ino][jno][idim][idim] += tmp1;
				}
				}
				for(unsigned int idim=0;idim<ndim;idim++){
					eMmat[ino][ino][idim][idim] += tmp1;
				}
			}
		}
		// call external force
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eforce_ex[ino][idim] += area*rho*g[idim]/(double)nnoes;
		}
		}

		////////////////////////////////
		
		double emat[nnoes][nnoes][ndim][ndim];	// element coeff matrix	
		double eres[nnoes][ndim];		// element residual vector		
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){
			(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i] + beta*dt*dt*(&eKmat[0][0][0][0])[i];
		}
		// get element residual vector
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = eforce_ex[ino][idim] - eforce_in[ino][idim];
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres[ino][idim] -= eMmat[ino][jno][idim][jdim]*eacc[jno][jdim];
			}
			}
		}
		}
		if( is_initial ){
			for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jno=0;jno<nnoes;jno++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres[ino][idim] -= dt*eKmat[ino][jno][idim][jdim]*evelo[jno][jdim]
						+ 0.5*dt*dt*eKmat[ino][jno][idim][jdim]*eacc[jno][jdim];
				}
				}
			}
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, ndim*ndim,&emat[0][0][0][0] );	// marge element stiffness
		for(unsigned int ino=0;ino<nnoes;ino++){	// marge element residual	
		for(unsigned int idim=0;idim<ndim;idim++){
			res_c.AddValue(noes[ino],idim,eres[ino][idim]);
		}
		}
	}
	return true;
}


// Dynamic
bool Fem::Eqn::AddLinSys_StVenant2D_NonStatic_NewmarkBeta
(double dt, double gamma, double beta,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y,
 const Fem::Field::CFieldWorld& world,
 unsigned int id_field_disp,
 bool is_initial, 
 unsigned int id_ea )
{
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR2 ) return false;
	
	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_StVenant2D_NonStatic_NewmarkBeta_P1
			(gamma,beta,dt,ls,
			 lambda, myu,
			 rho, g_x, g_y,
			 id_field_disp,world,
			 is_initial,
			 id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_StVenant2D_NonStatic_NewmarkBeta
			(dt, gamma, beta,
			 ls,
			 lambda, myu, rho, g_x, g_y,
			 world, id_field_disp, 
			 is_initial, 
			 id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}

////////////////////////////////////////////////

bool AddLinSys_StVenant2D_NonStatic_BackwardEular_P1
(double dt,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double rho, double g_x, double g_y,
 const unsigned int id_field_disp, const CFieldWorld& world, 
 const MatVec::CVector_Blk& velo_pre,
 bool is_initial,
 const unsigned int id_ea )
{
	//	std::cout << "St.Venant2D Triangle 3-point 1st order" << std::endl;
	
	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );
	
	const CField& field_disp = world.GetField(id_field_disp);
	
	const CElemAry::CElemSeg& es_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);
	const CElemAry::CElemSeg& es_va = field_disp.GetElemSeg(id_ea,CORNER,true,world);
	
	const unsigned int nno = 3;
	const unsigned int ndim = 2;
	
	unsigned int no[nno];	// elment node nuber 2 grobal node number
	
	double eKmat[nno][nno][ndim][ndim];	// element stiffness matrix
	double eMmat[nno][nno][ndim][ndim];	// element mass matrix		
	double eforce_in[nno][ndim];		// element residual vector		
	double eforce_ex[nno][ndim];		// element external force vector
	
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(id_field_disp,CORNER,world);	// matrix
	CVector_Blk&     res_c = ls.GetResidual(id_field_disp,CORNER,world);// residual vector
	
	const CNodeAry::CNodeSeg& ns_c_val  = field_disp.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_co   = field_disp.GetNodeSeg(CORNER,false,world);
	assert( ns_c_val.Length() == ndim );
	assert( ns_c_co.Length()  == ndim );
	
	double g[2] = { g_x, g_y };
	
	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		es_co.GetNodes(ielem,no);	// get global node number of coordinate
		double ecoords[nno][ndim];	// element node coordinate
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_co.GetValue(  no[ino], ecoords[ino]);
		}
		
		es_va.GetNodes(ielem,no);	// get global node nubmer of value
		double edisp[nno][ndim];	// displacement 
		double evelo[nno][ndim];	// velocity
		for(unsigned int ino=0;ino<nno;ino++){
			ns_c_val.GetValue( no[ino], edisp[ino]);
			ns_c_velo.GetValue(no[ino], evelo[ino]);
		}
		
		////////////////////////////////
		
		// set 0 to elemnet stiffness matrix, element residual
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eKmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eforce_ex[0][0]  +i) = 0.0; }
		for(unsigned int i=0;i<         nno*ndim;i++){ *(&eforce_in[0][0]  +i) = 0.0; }
		
		const double area = TriArea(ecoords[0],ecoords[1],ecoords[2]);
		
		double dldx[nno][ndim];		// spacial derivative of shape function
		double zero_order_term[nno];	// constant term of shape function
		TriDlDx(dldx, zero_order_term,   ecoords[0], ecoords[1], ecoords[2]);
		
		{	// calc stiffness matrix
            double dudx[ndim][ndim] = { { 0.0, 0.0}, {0.0, 0.0} };
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int idim=0;idim<ndim;idim++){
					for(unsigned int jdim=0;jdim<ndim;jdim++){
						dudx[idim][jdim] += edisp[ino][idim]*dldx[ino][jdim];
					}
				}
			}
			AddElemMatFin_StVenant2D( area, myu,lambda, nno,dudx,dldx,  eKmat[0],eforce_in );
		}
		{	// calc mass matrix
			const double tmp1 = rho*area/12.0;
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					for(unsigned int idim=0;idim<ndim;idim++){
						eMmat[ino][jno][idim][idim] += tmp1;
					}
				}
				for(unsigned int idim=0;idim<ndim;idim++){
					eMmat[ino][ino][idim][idim] += tmp1;
				}
			}
		}
		// call external force
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eforce_ex[ino][idim] += area*rho*g[idim]*0.33333333333333333333;
		}
		}
		
		////////////////////////////////
		
		double emat[nno][nno][ndim][ndim];	// element coeff matrix	
		double eres[nno][ndim];		// element residual vector		
		for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
			(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i] + dt*dt*(&eKmat[0][0][0][0])[i];
		}
		// get element residual vector
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = (eforce_ex[ino][idim] - eforce_in[ino][idim])*dt;
		}
		}
		if( is_initial ){
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres[ino][idim] -= eKmat[ino][jno][idim][jdim]*evelo[jno][jdim]*dt*dt;
				}
				}
			}
			}
		}
		else{
			double velo0[nno][ndim];
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				velo0[ino][idim] = velo_pre.GetValue(no[ino],idim);
			}
			}
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres[ino][idim] -= eMmat[ino][jno][idim][jdim]*(evelo[jno][jdim]-velo0[jno][jdim]);
				}
				}
			}
			}
		}
		
		////////////////////////////////
		
		mat_cc.Mearge(nno,no, nno,no, ndim*ndim,&emat[0][0][0][0] );	// marge element stiffness
		for(unsigned int ino=0;ino<nno;ino++){	// marge element residual	
			for(unsigned int idim=0;idim<ndim;idim++){
				res_c.AddValue(no[ino],idim,eres[ino][idim]);
			}
		}
	}
	return true;
}




// Dynamic
bool Fem::Eqn::AddLinSys_StVenant2D_NonStatic_BackwardEular
(double dt, 
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double rho, double g_x, double g_y,
 const Fem::Field::CFieldWorld& world,
 unsigned int id_field_disp, 
 const MatVec::CVector_Blk& velo_pre,
 bool is_initial,
 unsigned int id_ea )
{
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR2 ) return false;
	
	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TRI11 ){
			AddLinSys_StVenant2D_NonStatic_BackwardEular_P1
			(dt,
			 ls,
			 lambda, myu,
			 rho, g_x, g_y,
			 id_field_disp,world,
			 velo_pre,
			 is_initial,
			 id_ea);
		}
		else{ assert(0); }
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_StVenant2D_NonStatic_BackwardEular
			(dt,
			 ls,
			 lambda, myu, rho, g_x, g_y,
			 world, id_field_disp, 
			 velo_pre,
			 is_initial, 
			 id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}






////////////////////////////////////////////////////////////////

void SetElemMatFin_StVenant3D
(double detwei, double myu, double lambda, 
 const unsigned int nno, const double dudx[][3], const double dndx[][3],
 double eKMat[][3][3], double eForce_in[][3] )
{
	const unsigned int ndim = 3;

	double stress2[6];  // { s_00, s_11, s_22, s_01, s_12, s_20 }
    {
	    double strain2[6];  // { e_00, e_11, e_22, e_01, e_12, e_20 }
		strain2[0] = dudx[0][0] +              0.5*( dudx[0][0]*dudx[0][0] + dudx[1][0]*dudx[1][0] + dudx[2][0]*dudx[2][0] );
		strain2[1] = dudx[1][1] +              0.5*( dudx[0][1]*dudx[0][1] + dudx[1][1]*dudx[1][1] + dudx[2][1]*dudx[2][1] );
		strain2[2] = dudx[2][2] +              0.5*( dudx[0][2]*dudx[0][2] + dudx[1][2]*dudx[1][2] + dudx[2][2]*dudx[2][2] );
		strain2[3] = 0.5*( dudx[0][1] + dudx[1][0] + dudx[0][0]*dudx[0][1] + dudx[1][0]*dudx[1][1] + dudx[2][0]*dudx[2][1] );
		strain2[4] = 0.5*( dudx[1][2] + dudx[2][1] + dudx[0][1]*dudx[0][2] + dudx[1][1]*dudx[1][2] + dudx[2][1]*dudx[2][2] );
		strain2[5] = 0.5*( dudx[2][0] + dudx[0][2] + dudx[0][2]*dudx[0][0] + dudx[1][2]*dudx[1][0] + dudx[2][2]*dudx[2][0] );
	    const double dtmp1 = strain2[0] + strain2[1] + strain2[2];
        stress2[0] = 2.0*myu*strain2[0] + lambda*dtmp1;
        stress2[1] = 2.0*myu*strain2[1] + lambda*dtmp1;
        stress2[2] = 2.0*myu*strain2[2] + lambda*dtmp1;
        stress2[3] = 2.0*myu*strain2[3];
        stress2[4] = 2.0*myu*strain2[4];
        stress2[5] = 2.0*myu*strain2[5];
    }
	assert( nno <= 32 );
	double disp2strain2[32][ndim][6];     // { e_00, e_11, e_22, 2e_01, 2e_12, 2e_20 }
	{
        const double z_mat[ndim][ndim] = {
            { dudx[0][0]+1, dudx[0][1],     dudx[0][2] },
            { dudx[1][0],   dudx[1][1]+1,   dudx[1][2] },
            { dudx[2][0],   dudx[2][1],     dudx[2][2]+1 } };
		for(unsigned int kno=0;kno<nno;kno++){
		for(unsigned int kdim=0;kdim<ndim;kdim++){
			disp2strain2[kno][kdim][0] =  dndx[kno][0]*z_mat[kdim][0];
			disp2strain2[kno][kdim][1] =  dndx[kno][1]*z_mat[kdim][1];
			disp2strain2[kno][kdim][2] =  dndx[kno][2]*z_mat[kdim][2];
			disp2strain2[kno][kdim][3] = (dndx[kno][0]*z_mat[kdim][1] + dndx[kno][1]*z_mat[kdim][0] );
			disp2strain2[kno][kdim][4] = (dndx[kno][1]*z_mat[kdim][2] + dndx[kno][2]*z_mat[kdim][1] );
			disp2strain2[kno][kdim][5] = (dndx[kno][2]*z_mat[kdim][0] + dndx[kno][0]*z_mat[kdim][2] );
	    }
		}
	}
    ////////////////
	// 内力を求める
	for(unsigned int kno=0;kno<nno;kno++){
	for(unsigned int kdim=0;kdim<ndim;kdim++){
        double dtmp1 = 0;
    	dtmp1 += stress2[0]*disp2strain2[kno][kdim][0];
    	dtmp1 += stress2[1]*disp2strain2[kno][kdim][1];
    	dtmp1 += stress2[2]*disp2strain2[kno][kdim][2];
    	dtmp1 += stress2[3]*disp2strain2[kno][kdim][3];
    	dtmp1 += stress2[4]*disp2strain2[kno][kdim][4];
    	dtmp1 += stress2[5]*disp2strain2[kno][kdim][5];
        eForce_in[kno][kdim] = dtmp1*detwei;
	}
	}
	// 初期変位項を求める
    // 非対角ブロック成分
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<ino;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			double dtmp1 = 0.0;
			dtmp1 += disp2strain2[ino][idim][0]*disp2strain2[jno][jdim][0];
			dtmp1 += disp2strain2[ino][idim][1]*disp2strain2[jno][jdim][1];
			dtmp1 += disp2strain2[ino][idim][2]*disp2strain2[jno][jdim][2];
			dtmp1 += disp2strain2[ino][idim][3]*disp2strain2[jno][jdim][3]*0.5;
			dtmp1 += disp2strain2[ino][idim][4]*disp2strain2[jno][jdim][4]*0.5;
			dtmp1 += disp2strain2[ino][idim][5]*disp2strain2[jno][jdim][5]*0.5;
			const double dtmp2 = disp2strain2[ino][idim][0] + disp2strain2[ino][idim][1] + disp2strain2[ino][idim][2];
            const double dtmp3 = disp2strain2[jno][jdim][0] + disp2strain2[jno][jdim][1] + disp2strain2[jno][jdim][2];
            const double dtmp4 = detwei*(myu*2*dtmp1 + lambda*dtmp2*dtmp3);
		    eKMat[ino*nno+jno][idim][jdim] = dtmp4;
		    eKMat[jno*nno+ino][jdim][idim] = dtmp4;
		}
		}
    }
    }
    // 対角ブロック成分
	for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			double dtmp1 = 0.0;
			dtmp1 += disp2strain2[ino][idim][0]*disp2strain2[ino][jdim][0];
			dtmp1 += disp2strain2[ino][idim][1]*disp2strain2[ino][jdim][1];
			dtmp1 += disp2strain2[ino][idim][2]*disp2strain2[ino][jdim][2];
			dtmp1 += disp2strain2[ino][idim][3]*disp2strain2[ino][jdim][3]*0.5;
			dtmp1 += disp2strain2[ino][idim][4]*disp2strain2[ino][jdim][4]*0.5;
			dtmp1 += disp2strain2[ino][idim][5]*disp2strain2[ino][jdim][5]*0.5;
			const double dtmp2 = disp2strain2[ino][idim][0] + disp2strain2[ino][idim][1] + disp2strain2[ino][idim][2];
            const double dtmp3 = disp2strain2[ino][jdim][0] + disp2strain2[ino][jdim][1] + disp2strain2[ino][jdim][2];
		    eKMat[ino*nno+ino][idim][jdim] = detwei*(myu*2*dtmp1 + lambda*dtmp2*dtmp3);
		}
		}
    }
    // 初期応力項を求める
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		double dtmp1 = 0.0;
		dtmp1 += stress2[0]* dndx[ino][0]*dndx[jno][0];
		dtmp1 += stress2[1]* dndx[ino][1]*dndx[jno][1];
		dtmp1 += stress2[2]* dndx[ino][2]*dndx[jno][2];
		dtmp1 += stress2[3]*(dndx[ino][0]*dndx[jno][1] + dndx[ino][1]*dndx[jno][0] );
		dtmp1 += stress2[4]*(dndx[ino][1]*dndx[jno][2] + dndx[ino][2]*dndx[jno][1] );
		dtmp1 += stress2[5]*(dndx[ino][2]*dndx[jno][0] + dndx[ino][0]*dndx[jno][2] );
		eKMat[ino*nno+jno][0][0] += detwei*dtmp1;
		eKMat[ino*nno+jno][1][1] += detwei*dtmp1;
		eKMat[ino*nno+jno][2][2] += detwei*dtmp1;
	}
	}
}

////////////////

bool AddLinSys_StVenant3D_Static_P1
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y, double g_z,
 const unsigned int id_field_disp, const CFieldWorld& world, 
 unsigned int id_ea)
{
//	std::cout << "StVenant3D Tet 4-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_c = field_disp.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 4;	assert( nnoes == es_c.Length() );
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eforce_ex[nnoes][ndim];		// 要素内外力ベクトル
	double eforce_in[nnoes][ndim];		// 要素内内力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double ecoords[nnoes][ndim];		// 要素節点座標
	double edisp[  nnoes][ndim];		// 要素節点変位

	double dldx[nnoes][ndim];		// 形状関数の空間微分
	double zero_order_term[nnoes];	// 形状関数の定数項
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);// 要素剛性行列(コーナ-コーナー)
	CVector_Blk&    res_c  = ls.GetResidual(id_field_disp,CORNER,world);// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_disp.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

	double g[ndim] = { g_x, g_y, g_z };

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(noes[ino],ecoords[ino]);
			ns_c_val.GetValue(noes[ino],edisp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_ex[0][0] +i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_in[0][0] +i) = 0.0; }

		// 面積を求める
		const double vol = TetVolume(ecoords[0],ecoords[1],ecoords[2],ecoords[3]);

		// 形状関数のｘｙ微分を求める
		TetDlDx(dldx, zero_order_term,   ecoords[0],ecoords[1],ecoords[2],ecoords[3]);

		{	// 要素剛性行列，要素内力ベクトルを作る
			double dudx[ndim][ndim] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
			for(unsigned int ino=0;ino<nnoes;ino++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					dudx[idim][jdim] += edisp[ino][idim]*dldx[ino][jdim];
				}
				}
			}	
			SetElemMatFin_StVenant3D( vol,  myu,lambda,  nnoes,dudx,dldx,  emat[0],eforce_in );
		}
		// 外力ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				eforce_ex[ino][idim] += vol*rho*g[idim]/nnoes;
			}
		}

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = eforce_ex[ino][idim] - eforce_in[ino][idim];
		}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim, &emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			for(unsigned int idim=0;idim<ndim;idim++){
				res_c.AddValue(noes[ino],idim,eres[ino][idim]);
			}
		}
	}
	return true;
}

bool AddLinSys_StVenant3D_Static_Q1
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y, double g_z,
 const unsigned int id_field_disp, const CFieldWorld& world, 
 const unsigned int	id_ea)
{
//	std::cout << "StVenant3D Hex 8-point 1st order" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == HEX );

	assert( world.IsIdField(id_field_disp) );
	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_c = field_disp.GetElemSeg(id_ea,CORNER,true,world);

	unsigned int num_integral = 1;
	const unsigned int nInt = NIntLineGauss[num_integral];
	const double (*Gauss)[2] = LineGauss[num_integral];
	double detjac, detwei;

	const unsigned int nnoes = 8;	assert( nnoes == es_c.Length() );
	const unsigned int ndim = 3;

	unsigned int noes[nnoes];	// 要素内の節点の節点番号

	double emat[nnoes][nnoes][ndim][ndim];	// 要素剛性行列
	double eforce_ex[nnoes][ndim];		// 要素内外力ベクトル
	double eforce_in[nnoes][ndim];		// 要素内内力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	double ecoords[nnoes][ndim];		// 要素節点座標
	double edisp[  nnoes][ndim];		// 要素節点変位

	double dndx[nnoes][ndim];		// 形状関数の空間微分
	double an[nnoes];
				
	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);// 要素剛性行列(コーナ-コーナー)
	CVector_Blk&     res_c = ls.GetResidual(id_field_disp,CORNER,world);// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val = field_disp.GetNodeSeg(CORNER,true,world,VALUE);//.GetSeg(id_ns_c_val);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world,VALUE);//na_c_co.GetSeg(id_ns_c_co);

//	double g[ndim] = { g_x, g_y, g_z };

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		// 要素の節点番号を取ってくる
		es_c.GetNodes(ielem,noes);
		// 節点の座標、値を取ってくる
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(noes[ino],ecoords[ino]);
			ns_c_val.GetValue(noes[ino],edisp[ino]);
		}

		////////////////////////////////

		// 要素剛性行列、残差を０で初期化
		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&emat[0][0][0][0]+i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_ex[0][0] +i) = 0.0; }
		for(unsigned int i=0;i<           nnoes*ndim;i++){ *(&eforce_in[0][0] +i) = 0.0; }

		double vol = 0.0;
		for(unsigned int ir1=0;ir1<nInt;ir1++){
		for(unsigned int ir2=0;ir2<nInt;ir2++){
		for(unsigned int ir3=0;ir3<nInt;ir3++){
			const double r1 = Gauss[ir1][0];
			const double r2 = Gauss[ir2][0];
			const double r3 = Gauss[ir3][0];
			ShapeFunc_Hex8(r1,r2,r3,ecoords,detjac,dndx,an);
			detwei = detjac*Gauss[ir1][1]*Gauss[ir2][1]*Gauss[ir3][1];
			vol += detwei;
			{	// 要素剛性行列，要素内力ベクトルを作る
				double dudx[ndim][ndim] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
				for(unsigned int ino=0;ino<nnoes;ino++){
					for(unsigned int idim=0;idim<ndim;idim++){
					for(unsigned int jdim=0;jdim<ndim;jdim++){
						dudx[idim][jdim] += edisp[ino][idim]*dndx[ino][jdim];
					}
					}
				}	
				SetElemMatFin_StVenant3D( vol,  myu,lambda,  nnoes,dudx,dndx,  emat[0],eforce_in );
			}
			// 要素節点等価外力ベクトルを積n分毎に足し合わせる
			for(unsigned int ino=0;ino<nnoes;ino++){
				eforce_ex[ino][0] += detwei*rho*g_x;
				eforce_ex[ino][1] += detwei*rho*g_y;
				eforce_ex[ino][2] += detwei*rho*g_z;
			}
		}
		}
		}

		////////////////////////////////

		// 要素内残差ベクトルを求める
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = eforce_ex[ino][idim] - eforce_in[ino][idim];
		}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, // 全体剛性行列に要素剛性行列をマージ
			ndim*ndim, &emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){	// 要素内残差をマージ
			for(unsigned int idim=0;idim<ndim;idim++){
				res_c.AddValue(noes[ino],idim,eres[ino][idim]);
			}
		}
	}

	return true;
}

bool Fem::Eqn::AddLinSys_StVenant3D_Static
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y, double g_z,
 const CFieldWorld& world,
 unsigned int id_field_disp, 
 unsigned int id_ea )
{
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR3 ){ assert(0); return false; }
	
	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TET11 ){
			return AddLinSys_StVenant3D_Static_P1
			(ls,
			 lambda, myu,
			 rho, g_x, g_y, g_z,
			 id_field_disp,world,
			 id_ea);
		}
		else if( field_disp.GetInterpolationType(id_ea,world) == HEX11 ){
			return AddLinSys_StVenant3D_Static_Q1
			(ls,
			 lambda, myu,
			 rho, g_x, g_y, g_z,
			 id_field_disp,world,id_ea);
		}
		assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_StVenant3D_Static
			(ls,
			 lambda, myu, rho,  g_x, g_y, g_z,
			 world, id_field_disp, 
			 id_ea );
			if( !res ) return false;
		}
		return true;
	}
	return true;
}


/////////////////////////////////////////////////


bool AddLinSys_StVenant3D_NonStatic_NewmarkBeta_P1
(double gamma, double beta, double dt,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double rho, double g_x, double g_y, double g_z,
 const unsigned int id_field_disp, const CFieldWorld& world, 
 bool is_initial,
 const unsigned int id_ea)
{
//	std::cout << "St.Venant3D Tet1st" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TET );

	const CField& field_disp = world.GetField(id_field_disp);

	const CElemAry::CElemSeg& es_c = field_disp.GetElemSeg(id_ea,CORNER,true,world);

	const unsigned int nnoes = 4;
	const unsigned int ndim = 3;

	double emat[nnoes][nnoes][ndim][ndim];	// coefficient element matrix
	double eMmat[nnoes][nnoes][ndim][ndim];	// mass element matrix
	double eKmat[nnoes][nnoes][ndim][ndim];	// stiffness element matrix

	double eforce_in[nnoes][ndim];		// 要素内内力ベクトル
	double eres[nnoes][ndim];		// 要素内残差ベクトル

	CMatDia_BlkCrs& mat_cc = ls.GetMatrix(  id_field_disp,CORNER,world);// 要素剛性行列(コーナ-コーナー)
	CVector_Blk&     res_c = ls.GetResidual(id_field_disp,CORNER,world);// 要素残差ベクトル(コーナー)

	const CNodeAry::CNodeSeg& ns_c_val  = field_disp.GetNodeSeg(CORNER,true,world,VALUE);
	const CNodeAry::CNodeSeg& ns_c_velo = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_acc  = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);
	assert( ns_c_val.Length() == ndim );
	assert( ns_c_co.Length()  == ndim );

	const double g[ndim] = { g_x, g_y, g_z };

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
    {
	    unsigned int noes[nnoes];
		es_c.GetNodes(ielem,noes);
	    double ecoords[nnoes][ndim];		// coordinate 
	    double edisp[  nnoes][ndim];		// displacement 
	    double evelo[  nnoes][ndim];		// velocity
	    double eacc[  nnoes][ndim];			// acceleration
		for(unsigned int ino=0;ino<nnoes;ino++){
			ns_c_co.GetValue(  noes[ino], ecoords[ino]);
			ns_c_val.GetValue( noes[ino], edisp[ino]);
			ns_c_velo.GetValue(noes[ino], evelo[ino]);
			ns_c_acc.GetValue( noes[ino], eacc[ ino]);
		}

		////////////////////////////////

		const double vol = TetVolume(ecoords[0],ecoords[1],ecoords[2],ecoords[3]);

	    double dldx[nnoes][ndim];		// derivative of shape function
	    double zero_order_term[nnoes];	// constant term of shape function
		TetDlDx(dldx, zero_order_term,  ecoords[0],ecoords[1],ecoords[2],ecoords[3]);
		
		{	// calc stiffness matrix
			double dudx[ndim][ndim] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
			for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				dudx[idim][jdim] += edisp[ino][idim]*dldx[ino][jdim];
			}
			}
			}
			SetElemMatFin_StVenant3D( vol, myu,lambda, nnoes,dudx,dldx,  eKmat[0],eforce_in );
		}

		{	// calc mass matrix
		    for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){ *(&eMmat[0][0][0][0]+i) = 0.0; }
			const double tmp1 = rho*vol*0.05;
			for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int jno=0;jno<nnoes;jno++){
				eMmat[ino][jno][0][0] += tmp1;
				eMmat[ino][jno][1][1] += tmp1;
				eMmat[ino][jno][2][2] += tmp1;
            }
            }
			for(unsigned int ino=0;ino<nnoes;ino++){
				eMmat[ino][ino][0][0] += tmp1;
				eMmat[ino][ino][1][1] += tmp1;
				eMmat[ino][ino][2][2] += tmp1;
			}
		}

		////////////////////////////////

		for(unsigned int i=0;i<nnoes*nnoes*ndim*ndim;i++){
			(&emat[0][0][0][0])[i] = (&eMmat[0][0][0][0])[i] + beta*dt*dt*(&eKmat[0][0][0][0])[i];
		}
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres[ino][idim] = vol*rho*g[idim]*0.25 - eforce_in[ino][idim];
			for(unsigned int jno=0;jno<nnoes;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres[ino][idim] -= eMmat[ino][jno][idim][jdim]*eacc[jno][jdim];
			}
			}
		}
		}
		if( is_initial ){
			for(unsigned int ino=0;ino<nnoes;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jno=0;jno<nnoes;jno++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres[ino][idim] -= dt*eKmat[ino][jno][idim][jdim]*evelo[jno][jdim]
						+ 0.5*dt*dt*eKmat[ino][jno][idim][jdim]*eacc[jno][jdim];
				}
				}
			}
			}
		}

		////////////////////////////////
		
		mat_cc.Mearge(nnoes,noes, nnoes,noes, ndim*ndim,&emat[0][0][0][0] );
		for(unsigned int ino=0;ino<nnoes;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			res_c.AddValue(noes[ino],idim,eres[ino][idim]);
		}
		}
	}
	return true;
}


// 非定常
bool Fem::Eqn::AddLinSys_StVenant3D_NonStatic_NewmarkBeta
(double dt, double gamma, double beta,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double lambda, double myu,
 double  rho, double g_x, double g_y, double g_z,
 const Fem::Field::CFieldWorld& world,
 unsigned int id_field_disp,
 bool is_initial, 
 unsigned int id_ea )
{
	const CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != VECTOR3 ) return false;
	
	if( id_ea != 0 ){
		if( field_disp.GetInterpolationType(id_ea,world) == TET11 ){
			return AddLinSys_StVenant3D_NonStatic_NewmarkBeta_P1
			(gamma,beta,dt,ls,
			 lambda, myu, rho, g_x, g_y, g_z,
			 id_field_disp,world,
			 is_initial,
			 id_ea);
		}
		assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int>& aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinSys_StVenant3D_NonStatic_NewmarkBeta
			(dt, gamma, beta,
			 ls,
			 lambda, myu, rho,  g_x, g_y, g_z,
			 world, id_field_disp,
			 is_initial, 
			 id_ea );
			if( !res ) return false;
		}
		return true;
	}
	
	return true;
}




