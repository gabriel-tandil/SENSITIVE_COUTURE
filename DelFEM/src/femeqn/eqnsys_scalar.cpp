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
// eqnsys_scalar.cpp : スカラー型の連立方程式クラス
// (Fem::Eqn::CEqnSystem_Scalar2D,Fem::Eqn::CEqn_Scalar2D,Fem::Eqn::CEqn_Scalar3D)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )   // C4786なんて表示すんな( ﾟДﾟ)ｺﾞﾙｧ
#endif
#define for if(0); else for

#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/solver_ls_iter.h"
#include "delfem/ls/preconditioner.h"

#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"

#include "delfem/femeqn/eqn_diffusion.h"
#include "delfem/femeqn/eqn_poisson.h"
#include "delfem/femeqn/eqn_advection_diffusion.h"

#include "delfem/eqnsys_scalar.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;

bool CEqn_Scalar2D::AddLinSys( Fem::Ls::CLinearSystem_Field& ls, const CFieldWorld& world )
{
	if( this->m_IdFieldAdvec != 0  ){
		return Fem::Eqn::AddLinSys_AdvectionDiffusion_Static(
			ls,
			m_alpha, m_source,
			world,
			m_IdFieldVal, m_IdFieldAdvec,
			m_id_ea );
	}
	else{
		return Fem::Eqn::AddLinSys_Poisson(
			ls,
			m_alpha, m_source,
			world,
			m_IdFieldVal,
			m_id_ea );
	}
	return false;
}

bool CEqn_Scalar2D::AddLinSys_Newmark( 
		double dt, double gamma, CLinearSystem_Field& ls, 
		bool is_ax_sym,
		const CFieldWorld& world )
{
	if( this->m_IdFieldAdvec != 0 ){
		if( is_ax_sym ){
			std::cout << "Error!-->Not implimented" << std::endl;
			assert(0);
			abort();
		}
		return Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
			dt, gamma, 
			ls,
			m_capa, m_alpha, m_source,
			world,
			m_IdFieldVal, m_IdFieldAdvec,
			m_id_ea );
	}
	else{
		if( is_ax_sym ){
			return Fem::Eqn::AddLinSys_Diffusion_AxSym(
				dt, gamma,
				ls,
				m_capa, m_alpha, m_source,
				world,
				m_IdFieldVal,
				m_id_ea );
		}
		else{		
			return Fem::Eqn::AddLinSys_Diffusion(
				dt, gamma,
				ls,
				m_capa, m_alpha, m_source,
				world,
				m_IdFieldVal,
				m_id_ea );
		}
	}
	return false;
}

bool CEqn_Scalar2D::AddLinSys_Save( CLinearSystem_Save& ls, const CFieldWorld& world )
{
	if( this->m_IdFieldAdvec != 0 ){
		return Fem::Eqn::AddLinSys_AdvectionDiffusion_Static(
			ls,
			m_alpha, m_source,
			world, m_IdFieldVal,m_IdFieldAdvec, 
			m_id_ea );
	}
	else{
		return Fem::Eqn::AddLinSys_Poisson(
			ls,
			m_alpha, m_source,
			world, m_IdFieldVal,
			m_id_ea );
	}
	return false;
}

bool CEqn_Scalar2D::AddLinSys_SaveKDiaC( CLinearSystem_SaveDiaM_Newmark& ls, const CFieldWorld& world )
{
	if( this->m_IdFieldAdvec != 0 ){
		return Fem::Eqn::AddLinSys_AdvectionDiffusion_NonStatic_Newmark(
			ls,
			m_capa, m_alpha, m_source,
			world,
			m_IdFieldVal,m_IdFieldAdvec,
			m_id_ea );
	}
	else{
		return Fem::Eqn::AddLinSys_Diffusion(	// 全体剛性行列に拡散方程式を足し合わせる
			ls,
			m_capa, m_alpha, m_source,
			world,
			m_IdFieldVal,
			m_id_ea );
	}
	return true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


CEqnSystem_Scalar2D::CEqnSystem_Scalar2D()
{
	m_IsSaveStiffMat = false;
	m_IsAxialSymmetry = false;
}

bool CEqnSystem_Scalar2D::SetDomain_Field(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
	m_IdFieldVal = world.MakeField_FieldElemDim(id_base,2,SCALAR,VELOCITY|VALUE,CORNER);
	{
		m_aEqn.clear();
		const CField& field = world.GetField(m_IdFieldVal);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			m_aEqn.push_back( CEqn_Scalar2D(id_ea,m_IdFieldVal) );
		}
	}
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}

bool CEqnSystem_Scalar2D::SetDomain_FieldElemAry(unsigned int id_base, unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
	m_IdFieldVal = world.MakeField_FieldElemAry(id_base,id_ea,SCALAR,VELOCITY|VALUE,CORNER);
	{
		m_aEqn.clear();
		const CField& field = world.GetField(m_IdFieldVal);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			m_aEqn.push_back( CEqn_Scalar2D(id_ea,m_IdFieldVal) );
		}
	}
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}

bool CEqnSystem_Scalar2D::MatrixProperty(bool& is_c, bool& is_m, bool& is_asym)
{
	if( m_aEqn.size() == 0 ) return false;
	is_c = false;
	is_m = false;
	is_asym = false;
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if(  m_aEqn[ieqn].IsAdvection()  ) is_asym = true;
	}
	return true;
}

bool CEqnSystem_Scalar2D::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{	
	////////////////////////////////
	// 連立一次方程式クラスの設定
	if( this->pLS != 0 ){ delete pLS; pLS = 0; }
	if( this->m_IsSaveStiffMat ){
		if( !m_IsStationary ){ pLS = new CLinearSystem_SaveDiaM_Newmark; }
		else{ pLS = new CLinearSystem_Save; }
	}
	else{ pLS = new CLinearSystem_Field; }
	pLS->AddPattern_Field(m_IdFieldVal,world);	// val_fieldからできる全体剛性行列を追加する
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		unsigned int id_field = m_aIdFixField[idf].first; 
		pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_fieldを固定境界条件に設定		
	}

	////////////////////////////////
	// 前処理クラスの作成
	if( this->pPrec != 0 ){ delete pPrec; pPrec = 0; }
    pPrec = new LsSol::CPreconditioner_ILU;
	pPrec->SetLinearSystem( (*pLS).m_ls );
	return true;
}

double CEqnSystem_Scalar2D::MakeLinearSystem( const Fem::Field::CFieldWorld& world)
{	
	if( pLS==0 || pPrec==0 ){ this->InitializeLinearSystem(world); }
	////////////////////////////////
	// 連立一次方程式を初期化する
	pLS->InitializeMarge();	
	////////////////////////////////
	// 連立一次方程式に要素剛性をマージする
	if( this->m_IsSaveStiffMat ){
		if( !m_IsStationary ){
			(*(Fem::Ls::CLinearSystem_SaveDiaM_Newmark*)pLS).SetNewmarkParameter(m_gamma_newmark, m_dt);
			for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
				m_aEqn[ieqn].AddLinSys_SaveKDiaC( (*(Fem::Ls::CLinearSystem_SaveDiaM_Newmark*)pLS), world);
			}
		}
		else{
			for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
				m_aEqn[ieqn].AddLinSys_Save( (*(Fem::Ls::CLinearSystem_Save*)pLS), world);
			}
		}
	}
	else{
		if( !m_IsStationary ){
			for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
				m_aEqn[ieqn].AddLinSys_Newmark( m_dt, m_gamma_newmark, *pLS, 
					this->m_IsAxialSymmetry,
					world);
			}
		}
		else{
			for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
				m_aEqn[ieqn].AddLinSys( *pLS, world);
			}
		}
			
	}
	////////////////
	// マージの終了
	if( this->m_IsSaveStiffMat ){
		if( this->m_IsStationary ){
			(*(Fem::Ls::CLinearSystem_Save*)pLS).FinalizeMarge();
		}
		else{
			(*(Fem::Ls::CLinearSystem_SaveDiaM_Newmark*)pLS).FinalizeMarge();
		}
	}
	else{
		pLS->FinalizeMarge();
	}
	this->m_is_cleared_value_ls   = false;

	// 前処理行列を作る
	pPrec->SetValue( (*pLS).m_ls );
	this->m_is_cleared_value_prec = false;

	return 0.0;	// TODO : 本当は残差の値を返さなくてはならない
}

bool CEqnSystem_Scalar2D::Solve(Fem::Field::CFieldWorld& world)
{
	////////////////////////////////
	// 行列を作る
	if( pLS == 0 || pPrec == 0 ){
		this->InitializeLinearSystem(world);
		this->MakeLinearSystem(world);
	}
	else{
		if( !this->m_IsSaveStiffMat ){ this->MakeLinearSystem(world); }
	}
	if( this->m_IsSaveStiffMat ){ 
		if( this->m_is_cleared_value_ls || this->m_is_cleared_value_prec ){
			this->MakeLinearSystem(world);  
		}
		double res = pLS->MakeResidual(world); 
//		std::cout << "Residual : " << res << std::endl;
	}
	assert( this->m_is_cleared_value_ls   == false );
	assert( this->m_is_cleared_value_prec == false );

	////////////////////////////////
	// 行列を解く
	{
		double conv_ratio = 1.0e-6;
		unsigned int max_iter = 1000;
		bool is_c, is_m, is_asym;
		this->MatrixProperty(is_c,is_m,is_asym);
		if( is_asym ){
            LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
			LsSol::Solve_PBiCGSTAB(conv_ratio,max_iter,lsp);	
		//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
		}
		else{
            LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
			LsSol::Solve_PCG(conv_ratio,max_iter,lsp);	
		//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
		}
		this->m_aItrNormRes.clear();
		this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
//		std::cout << max_iter << " " << conv_ratio << std::endl;
	}

	////////////////////////////////
	// 解を更新する

	if( this->m_IsSaveStiffMat ){
		if( m_IsStationary ){
			pLS->UpdateValueOfField(m_IdFieldVal,world,VALUE); 
		}
		else{
			pLS->UpdateValueOfField(m_IdFieldVal,world,VELOCITY); 
		}
	}
	else{
		if( m_IsStationary ){
			pLS->UpdateValueOfField(m_IdFieldVal,world,VALUE); 
		}
		else{ 
			pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_IdFieldVal,world,VELOCITY,true); 
		}
	}

	return true;
}


bool CEqnSystem_Scalar2D::SetEquation( const CEqn_Scalar2D& eqn )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == eqn.GetIdEA() ){
			m_aEqn[ieqn] = eqn;
			this->m_is_cleared_value_ls   = true;
			this->m_is_cleared_value_prec = true;
			return true;
		}
	}
	return false;
}

CEqn_Scalar2D CEqnSystem_Scalar2D::GetEquation(unsigned int id_ea) const
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == id_ea ){ return m_aEqn[ieqn]; }
	}
	return CEqn_Scalar2D(0,0);
}

void CEqnSystem_Scalar2D::SetAlpha(double alpha)
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetAlpha(alpha);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Scalar2D::SetCapacity(double capa)
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetCapacity(capa);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Scalar2D::SetSource(double source)
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetSource(source);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Scalar2D::SetAdvection(unsigned int id_field_advec){
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		m_aEqn[ieqn].SetAdvection(id_field_advec);
	}
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}

void CEqnSystem_Scalar2D::SetStationary(bool is_stat){
	this->m_IsStationary = is_stat;
	this->ClearLinearSystemPreconditioner();
}


void CEqnSystem_Scalar2D::SetAxialSymmetry(bool is_ax_sym){
	this->m_IsAxialSymmetry = is_ax_sym;
	this->m_is_cleared_value_ls   = true;
	this->m_is_cleared_value_prec = true;
}


bool CEqnSystem_Scalar2D::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return true;
}

unsigned int CEqnSystem_Scalar2D::AddFixElemAry( 
		unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

bool CEqnSystem_Scalar2D::ClearFixElemAry(
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


void CEqnSystem_Scalar2D::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

unsigned int CEqnSystem_Scalar2D::AddFixElemAry( 
		const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_IdFieldVal, aIdEA );
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



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



CEqn_Scalar3D::CEqn_Scalar3D() 
: m_alpha(1.0), m_source(0.0)
{
}

bool CEqn_Scalar3D::SetDomain(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
	m_IdFieldVal = world.MakeField_FieldElemDim(id_base,3,SCALAR,VALUE,CORNER);
	this->ClearLinearSystemPreconditioner();
	this->m_aIdFixField.clear();
	return true;
}


bool CEqn_Scalar3D::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{	
	// 連立一次方程式クラスの設定
	pLS = new CLinearSystem_Field;
	pLS->AddPattern_Field(m_IdFieldVal,world);	// val_fieldからできる全体剛性行列を追加する
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		unsigned int id_field = m_aIdFixField[idf].first; 
		pLS->SetFixedBoundaryCondition_Field( id_field, 0, world ); // bc0_fieldを固定境界条件に設定		
	}
	// 前処理クラスの作成
    pPrec = new LsSol::CPreconditioner_ILU;
	pPrec->SetLinearSystem( (*pLS).m_ls );
	return true;
}

double CEqn_Scalar3D::MakeLinearSystem(const Fem::Field::CFieldWorld& world)
{	
	if( pLS==0 || pPrec==0 ) this->InitializeLinearSystem(world);
	// 連立一次方程式を作る
	pLS->InitializeMarge();	// 連立一次方程式を初期化する
	Fem::Eqn::AddLinSys_Poisson(
		*pLS,
		m_alpha, m_source,
		world,
		m_IdFieldVal);
	pLS->FinalizeMarge();
	this->m_is_cleared_value_ls = false;

	// 前処理行列を作る
	pPrec->SetValue( (*pLS).m_ls );
	this->m_is_cleared_value_prec = false;

	return true;
}

bool CEqn_Scalar3D::Solve(Fem::Field::CFieldWorld& world)
{
	////////////////////////////////
	// 行列を作る
	if( pLS == 0 || pPrec == 0 ){
		this->InitializeLinearSystem(world);
	}
	this->MakeLinearSystem(world);
	assert( this->m_is_cleared_value_ls   == false );
	assert( this->m_is_cleared_value_prec == false );

	////////////////////////////////
	// 行列を解く
	{
		double conv_ratio = 1.0e-6;
		unsigned int max_iter = 1000;
		// Solve with Preconditioned Conjugate Gradient
        LsSol::CLinearSystemPreconditioner lsp( (*pLS).m_ls, *pPrec );
		LsSol::Solve_PCG(conv_ratio,max_iter,lsp);
		// Solve with Conjugate Gradient
		//	Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);
		this->m_aItrNormRes.clear();
		this->m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
//		std::cout << max_iter << " " << conv_ratio << std::endl;
	}

	////////////////////////////////
	// 解を更新する
	pLS->UpdateValueOfField(m_IdFieldVal,world,VALUE); 
	return true;
}

bool CEqn_Scalar3D::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	this->ClearLinearSystem();
	return true;
}

unsigned int CEqn_Scalar3D::AddFixElemAry( 
		unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

bool CEqn_Scalar3D::ClearFixElemAry(
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


void CEqn_Scalar3D::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

unsigned int CEqn_Scalar3D::AddFixElemAry( 
		const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_IdFieldVal, aIdEA );
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

