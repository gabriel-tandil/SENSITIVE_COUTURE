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
	#pragma warning( disable : 4786 )
#endif

#include <math.h>
#include <time.h>

#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/eqn_dkt.h"
#include "delfem/eqnsys_shell.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// 3D equation

Fem::Eqn::CEqnSystem_DKT::CEqnSystem_DKT(unsigned int id_base, Fem::Field::CFieldWorld& world) 
: m_IsGeomNonlin(false), m_IsStationary(true), m_IsSaveStiffMat(false)
{
	isnt_direct = false;
	m_young = 1.0;
	m_poisson = 0.0;
	m_thick = 0.1;
	m_arealrho = 1.0;
	m_press= 0.0;
	
	m_g_x = 0.0;
	m_g_y = 0.0;
	m_g_z = 0.0;

	this->SetDomain_Field(id_base, world);
}

CEqnSystem_DKT::CEqnSystem_DKT()
: m_IsGeomNonlin(false), m_IsStationary(true), m_IsSaveStiffMat(false)
{
	isnt_direct = false;
//	isnt_direct = true;
	m_young = 1.0;
	m_poisson = 0.0;
	m_arealrho = 1.0;
	m_thick = 0.1;
	m_press= 0.0;
	
	m_g_x = 0.0;
	m_g_y = 0.0;
	m_g_z = 0.0;
}

double CEqnSystem_DKT::MakeLinearSystem(const Fem::Field::CFieldWorld& world, bool is_initial)
{	
	if( pLS==0 || pPrec==0 ) this->InitializeLinearSystem(world);

	clock_t t1;
    clock_t t2;

	t1 = clock();

	// 連立一次方程式を作る
	pLS->InitializeMarge();	// 連立一次方程式を初期化する
	if( this->m_IsGeomNonlin ){
		if( this->m_IsStationary ){
			Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_Static(
				*pLS,
				m_young, m_poisson, m_thick, m_arealrho,
				m_g_x,  m_g_y,  m_g_z, m_press,
				world,m_id_disp,m_id_rot);
		}
		else{ // 動的St.Vennat-Kirchhoff体
			Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_NonStatic(
				m_dt, m_gamma_newmark, m_beta_newmark, is_initial,
				*pLS,
				m_young, m_poisson, m_thick, m_arealrho,
				m_g_x,  m_g_y,  m_g_z, m_press,
				world,m_id_disp,m_id_rot);
		}
	}
	else{
		if( this->m_IsSaveStiffMat ){
			if( this->m_IsStationary ){
				Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static_Save(
					*(CLinearSystem_Save*)pLS, 
					m_young, m_poisson, m_thick, m_arealrho,
					m_g_x,  m_g_y,  m_g_z, m_press,
					world,m_id_disp,m_id_rot);
			}
			else{
				(*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS).SetNewmarkParameter(this->m_beta_newmark,this->m_gamma_newmark,this->m_dt);
				Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic_Save(
					*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS,
					m_young, m_poisson, m_thick, m_arealrho,
					m_g_x,  m_g_y,  m_g_z, m_press,
					world,m_id_disp,m_id_rot);
			}
		}
		else{
			if( this->m_IsStationary ){
				Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static(
					*pLS, 
					m_young, m_poisson, m_thick, m_arealrho,
					m_g_x,  m_g_y,  m_g_z, m_press,
					world,m_id_disp,m_id_rot);
			}
			else{ // 動的線形弾性体
				Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic(
					m_dt, m_gamma_newmark, m_beta_newmark,
					*pLS,
					m_young, m_poisson, m_thick, m_arealrho,
					m_g_x,  m_g_y,  m_g_z, m_press,
					world,m_id_disp,m_id_rot);
			}
		}
	}
	double norm_res;
	if( this->m_IsSaveStiffMat ){
		if( this->m_IsStationary ){
			norm_res = (*(CLinearSystem_SaveDiaM_Newmark*)pLS).FinalizeMarge();
		}
		else{
			norm_res = (*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS).FinalizeMarge();
		}
	}
	else{
		norm_res = pLS->FinalizeMarge();
	}
	this->m_is_cleared_value_ls = false;

    t2 = clock();
//    printf("time = %10.3f\n", (double)(t2 - t1)/CLOCKS_PER_SEC );

	// 前処理行列を作る
	if( !isnt_direct ){
		pPrec->SetValue((*pLS).m_ls);
	}

	return norm_res;
}

bool CEqnSystem_DKT::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{
	if( pLS!=0 || pPrec!=0 ) ClearLinearSystemPreconditioner();

	// 連立一次方程式クラスの作成
	assert( pLS == 0 );
	if( this->m_IsSaveStiffMat ){
		if( this->m_IsStationary ){
			pLS = new CLinearSystem_Save;
		}
		else{
			pLS = new CLinearSystem_SaveDiaM_NewmarkBeta;
		}
	}
	else{
		pLS = new CLinearSystem_Field;
	}

	if(!isnt_direct){
		if( !m_IsSaveStiffMat ){
			pLS->AddPattern_CombinedField(m_id_disp,m_id_rot,world);
		}
		else{
			if( this->m_IsStationary ){
				(*(CLinearSystem_Save*)pLS).AddPattern_CombinedField(m_id_disp,m_id_rot,world);
			}
			else{
				(*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS).AddPattern_CombinedField(m_id_disp,m_id_rot,world);
			}
		}
		for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
			const unsigned int id_field = m_aIdFixField[idf].first; 
			const int idof     = m_aIdFixField[idf].second;
			if( idof == -1 ){
				pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_fieldを固定境界条件に設定	
				pLS->SetFixedBoundaryCondition_Field( id_field, 1, world ); // bc0_fieldを固定境界条件に設定	
				pLS->SetFixedBoundaryCondition_Field( id_field, 2, world ); // bc0_fieldを固定境界条件に設定
			}
			else{
				pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_fieldを固定境界条件に設定	
			}
		}
    pPrec = new LsSol::CPreconditioner_ILU;
//    ((LsSol::CPreconditioner_ILU*)pPrec)->SetOrdering(true);
    ((LsSol::CPreconditioner_ILU*)pPrec)->SetFillInLevel(2);
		pPrec->SetLinearSystem((*pLS).m_ls);
	}
	else{
		// 連立一次方程式クラスの設定
		pLS->AddPattern_Field(m_id_disp,world);	// val_fieldからできる全体剛性行列を追加する
		pLS->AddPattern_Field(m_id_rot,m_id_disp,world);	// val_fieldからできる全体剛性行列を追加する
		for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
			const unsigned int id_field = m_aIdFixField[idf].first; 
			const int idof     = m_aIdFixField[idf].second;
			if( idof == -1 ){
				pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_fieldを固定境界条件に設定	
				pLS->SetFixedBoundaryCondition_Field( id_field, 1, world ); // bc0_fieldを固定境界条件に設定	
				pLS->SetFixedBoundaryCondition_Field( id_field, 2, world ); // bc0_fieldを固定境界条件に設定
			}
			else{
				pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_fieldを固定境界条件に設定	
			}
		}
		// 前処理クラスの作成
		assert( pPrec == 0 );
//		pPrec = new CPreconditioner_ILU(*pLS);
	}
		

	return true;
}

bool CEqnSystem_DKT::Solve(Fem::Field::CFieldWorld& world)
{
	if( this->m_IsGeomNonlin ){	// nonlinear
		assert( !this->m_IsSaveStiffMat );
		if( pLS == 0 || pPrec == 0 ){ this->InitializeLinearSystem(world); }
		double ini_norm_res;	
		for(unsigned int iitr=0;iitr<1;iitr++){
//			std::cout << iitr << std::endl;
			double norm_res = this->MakeLinearSystem(world,iitr==0);
			if( iitr==0 ){
				if( norm_res < 1.0e-20 ) break;
				ini_norm_res = norm_res;
			}
			if( norm_res < ini_norm_res * 1.0e-6 ) break;
//			std::cout << iitr << " " << norm_res << " " << ini_norm_res << " " << norm_res / ini_norm_res << std::endl;
			{	// solve matrix
				double conv_ratio = 1.0e-5;
				unsigned int max_iter = 4000;
				// Solve with Preconditioned Conjugate Gradient
				if( isnt_direct ){
					LsSol::Solve_CG(conv_ratio,max_iter,*pLS);
				}
				else{
          LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
					LsSol::Solve_PCG(conv_ratio,max_iter,lsp);
				}
				// Solve with Conjugate Gradient
//				std::cout << max_iter << " " << conv_ratio << std::endl;
			}
			if( this->m_IsStationary ){
				pLS->UpdateValueOfField(m_id_disp,world,VALUE);
//				pLS->UpdateValueOfField_Rotate(m_id_rot,world,VALUE);
				pLS->UpdateValueOfField(m_id_rot,world,VALUE);
			}
			else{
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
                    this->m_id_disp,world, iitr==0 );
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
					m_id_rot,world, iitr==0 );
			}
		}
	}
	else{  
		if( pLS == 0 ){
			this->InitializeLinearSystem(world);
			double res = this->MakeLinearSystem(world,true);
			std::cout << "Residual0 : " << res << std::endl;
		}	
		else{
			if( !this->m_IsSaveStiffMat ){ 
				double res = this->MakeLinearSystem(world,true);
				std::cout << "Residual1 : " << res << std::endl;
			}
		}
		if( this->m_IsSaveStiffMat ){ 
			if( this->m_is_cleared_value_ls ){
				this->MakeLinearSystem(world,true);
			}
			double res = 0;
			if( this->m_IsStationary ){
				res = ((CLinearSystem_Save*)pLS)->MakeResidual(world); 
			}
			else{
				res = ((CLinearSystem_SaveDiaM_NewmarkBeta*)pLS)->MakeResidual(world); 
			}
			std::cout << "Residual2 : " << res << std::endl;
		}
		assert( this->m_is_cleared_value_ls   == false );

		{	// 行列を解く
			double conv_ratio = 1.0e-5;
			unsigned int max_iter = 4000;
			if( isnt_direct ){
				LsSol::Solve_CG(conv_ratio,max_iter,*pLS);
			}
			else{
                LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
				LsSol::Solve_PCG(conv_ratio,max_iter,lsp);
			}
			std::cout << max_iter << " " << conv_ratio << std::endl;
		}

		// 解を更新する
		if( this->m_IsStationary ){
			if( !this->m_IsSaveStiffMat ){
				pLS->UpdateValueOfField(m_id_disp,world,VALUE);
				pLS->UpdateValueOfField(m_id_rot, world,VALUE);
			}
			else{
				*(CLinearSystem_Save*)pLS->UpdateValueOfField(m_id_disp,world,VALUE);
				*(CLinearSystem_Save*)pLS->UpdateValueOfField(m_id_rot, world,VALUE);
			}
//			pLS->UpdateValueOfField_Rotate(m_id_rot,world,VALUE);	
		}
		else{
			if( !this->m_IsSaveStiffMat ){
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
					m_id_disp,world,true);
				pLS->UpdateValueOfField_NewmarkBeta(m_gamma_newmark,m_beta_newmark,m_dt,
					m_id_rot,world,true);
			}
			else{
				(*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS).UpdateValueOfField(m_id_disp,world,ACCELERATION);
				(*(CLinearSystem_SaveDiaM_NewmarkBeta*)pLS).UpdateValueOfField(m_id_rot, world,ACCELERATION);
			}
		}
	}
	return true;
}

bool CEqnSystem_DKT::SetDomain_FieldElemAry(unsigned int id_base, unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
  {   // 入力フィールドの座標節点セグメントのdofが３かどうかチェックする
    assert( world.IsIdField(id_base) );
    Fem::Field::CField& field = world.GetField(id_base);
    if( field.GetNDimCoord() != 3 ){
      this->ClearLinearSystemPreconditioner();
      m_id_disp = 0;
      m_id_rot = 0;
      return false;
    }
  }
	m_id_disp = world.MakeField_FieldElemAry(id_base,id_ea,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
	m_id_rot  = world.MakeField_FieldElemAry(id_base,id_ea,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
	this->ClearLinearSystemPreconditioner();
	return true;
}

bool CEqnSystem_DKT::SetDomain_Field(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
  {   // 入力フィールドの座標節点セグメントのdofが３かどうかチェックする
    //        const unsigned int id_base = world.GetFieldBaseID();
    assert( world.IsIdField(id_base) );
    Fem::Field::CField& field = world.GetField(id_base);
    if( field.GetNDimCoord() != 3 ){
      this->ClearLinearSystemPreconditioner();
      m_id_disp = 0;
      m_id_rot = 0;
      return false;
    }
  }
	m_id_disp = world.MakeField_FieldElemDim(id_base,2,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
	m_id_rot  = world.MakeField_FieldElemDim(id_base,2,VECTOR3,VALUE|VELOCITY|ACCELERATION,CORNER);
	this->ClearLinearSystemPreconditioner();
	return true;
}


bool CEqnSystem_DKT::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return true;
}

unsigned int CEqnSystem_DKT::AddFixElemAry( 
		unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

bool CEqnSystem_DKT::ClearFixElemAry(
		unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdEA( id_ea ) ) return false;
	for(unsigned int ifix=0;ifix<m_aIdFixField.size();ifix++){
		const unsigned int id_field_fix = m_aIdFixField[ifix].first;
		const Fem::Field::CField& field = world.GetField(id_field_fix);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		if( aIdEA.size() != 1 ){
			std::cout << "Error!-->Not Implimented" << std::endl;
			assert(0);
		}
		if( aIdEA[0] == id_ea ){
			m_aIdFixField.erase( m_aIdFixField.begin()+ifix );
			this->ClearLinearSystem();
		}
	}
	return true;
}


void CEqnSystem_DKT::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

unsigned int CEqnSystem_DKT::AddFixElemAry( 
		const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_id_disp, aIdEA );
	if( id_field == 0 ) return 0;
	assert( world.IsIdField(id_field) );
	{
		CField& field = world.GetField(id_field);
		unsigned int nlen_val = field.GetNLenValue();
		for(unsigned int ilen=0;ilen<nlen_val;ilen++){
      Fem::Field::SetFieldValue_Constant(id_field,ilen,Fem::Field::VALUE,       world,0);
      Fem::Field::SetFieldValue_Constant(id_field,ilen,Fem::Field::VELOCITY,    world,0);
      Fem::Field::SetFieldValue_Constant(id_field,ilen,Fem::Field::ACCELERATION,world,0);
		}
	}
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return id_field;
}
