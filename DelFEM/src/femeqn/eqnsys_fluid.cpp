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
// eqnsys_fluid.cpp : 
// 流体の連立方程式クラス(CEqnSystem_Fluid2D,CEqn_Fluid2D,CEqn_Fluid3D)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/linearsystem_interface_solver.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/ls/solver_ls_iter.h"

#include "delfem/femls/linearsystem_field.h"

#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_tet.h"

#include "delfem/femeqn/eqn_stokes.h"
#include "delfem/femeqn/eqn_navier_stokes.h"

#include "delfem/eqnsys_fluid.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;

using namespace Fem::Eqn;

bool CEqn_Fluid2D::AddLinSys(Fem::Ls::CLinearSystem_Field& ls,
                             const Fem::Field::CFieldWorld& world )
{
//	std::cout << "CEqn_Fluid2D::AddLinSys" << std::endl;
	if( !world.IsIdEA(m_id_ea) ) return false;
	return Fem::Eqn::AddLinSys_Stokes2D_Static(
		m_myu, m_g_x, m_g_y,
		ls,
		this->m_IdFieldVelo, this->m_IdFieldPress,world,
		m_id_ea );
//	return false;
}

// 連立一次方程式マージメソッド
bool CEqn_Fluid2D::AddLinSys_NewmarkBetaAPrime( double dt,
    double gamma, double beta, bool is_initial,
	Fem::Ls::CLinearSystem_Field& ls, const Fem::Field::CFieldWorld& world )
{
//    std::cout << "CEqn_Fluid2D::AddLinSys_NewmakrBetaAPrime " << m_id_ea << std::endl;
	if( !world.IsIdEA(m_id_ea) ) return false;
	if( m_IsNavierStokes ){
		assert( !this->m_IsNavierStokesALE );
        assert( !this->m_IsStokes );
//        std::cout << "CEqn_Fluid2D::AddLinSys_NewmarkBetaAPrime  ns " << m_id_ea << std::endl;
		return Fem::Eqn::AddLinSys_NavierStokes2D_NonStatic_Newmark(
			dt, gamma, ls,
			m_rho,m_myu, m_g_x,m_g_y,
			this->m_IdFieldVelo,this->m_IdFieldPress,world,
			m_id_ea );
	}
	else if( m_IsNavierStokesALE ){
		assert( !this->m_IsStokes );
		assert( !this->m_IsNavierStokes );
		return Fem::Eqn::AddLinSys_NavierStokesALE2D_NonStatic_Newmark(
			dt, gamma, ls,
			m_rho,m_myu, m_g_x,m_g_y,
			this->m_IdFieldVelo,this->m_IdFieldPress,this->m_IdFieldMshVelo, world );
	}
	else if( this->m_IsStokes ){
		assert( !this->m_IsNavierStokes );
		assert( !this->m_IsNavierStokesALE );
//        std::cout << "CEqn_Fluid2D::AddLinSys_NewmarkBetaAPrime  stokes " << m_id_ea << std::endl;
		return Fem::Eqn::AddLinSys_Stokes2D_NonStatic_Newmark(
			m_rho,m_myu,m_g_x,m_g_y,
			gamma,dt,
			ls,
			this->m_IdFieldVelo,this->m_IdFieldPress,world,
			m_id_ea);
	}
	else{
		assert(0);
	}
	return false;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



CEqnSystem_Fluid2D::CEqnSystem_Fluid2D(unsigned int id_field, Fem::Field::CFieldWorld& world)
{
	this->Clear();
	this->UpdateDomain_Field(id_field,world);
}

CEqnSystem_Fluid2D::CEqnSystem_Fluid2D()
{
	this->Clear();
}

bool CEqnSystem_Fluid2D::SetEquation( const CEqn_Fluid2D& eqn )
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == eqn.GetIdEA() ){
			m_aEqn[ieqn] = eqn;
			this->ClearLinearSystemPreconditioner();
			return true;
		}
	}
	m_aEqn.push_back( eqn );
	this->ClearLinearSystemPreconditioner();
	return true;
}

CEqn_Fluid2D CEqnSystem_Fluid2D::GetEquation(unsigned int id_ea) const
{
	for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
		if( m_aEqn[ieqn].GetIdEA() == id_ea ){
			return m_aEqn[ieqn];
		}
	}
	return CEqn_Fluid2D(0,0,0);
}

double CEqnSystem_Fluid2D::MakeLinearSystem(const Fem::Field::CFieldWorld& world)
{	
	if( pLS==0 || pPrec==0 ) this->InitializeLinearSystem(world);
	// 連立一次方程式を作る
	pLS->InitializeMarge();	// 連立一次方程式を初期化する(0クリア)

	if( m_IsStationary ){
		for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
			m_aEqn[ieqn].AddLinSys(*pLS,world);
		}
	}
	else{
		for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
//            std::cout << "CEqnSystem_Fluid2D::MakeLinearSystem " << ieqn << " " << m_aEqn[ieqn].GetIdEA() << std::endl;
			m_aEqn[ieqn].AddLinSys_NewmarkBetaAPrime(m_dt,m_gamma_newmark,m_beta_newmark,true, 
				*pLS, world);
		}
	}
    if( world.IsIdField( this->m_id_force ) )
    {
        assert( world.IsIdField(this->m_id_force) );
        const Fem::Field::CField& ff = world.GetField(this->m_id_force);
        const Fem::Field::CNodeAry::CNodeSeg& nsf_v = ff.GetNodeSeg(Fem::Field::CORNER,true, world,Fem::Field::VELOCITY);
        const unsigned int nno = nsf_v.Size();
        assert( nsf_v.Length() == 2 );
        ////////////////
        assert( world.IsIdField(this->m_id_velo) );
        const Fem::Field::CField& fv = world.GetField(this->m_id_velo);
        const Fem::Field::CNodeAry::CNodeSeg& nsv_v = fv.GetNodeSeg(Fem::Field::CORNER,true, world,Fem::Field::VELOCITY);
        assert( nsv_v.Size() == nno );
        assert( nsv_v.Length() == 2 );
        ////////////////
        MatVec::CVector_Blk& vec_res = pLS->GetResidual(this->m_id_velo,Fem::Field::CORNER,world);
        assert( vec_res.NBlk() == nno );
        assert( vec_res.Len() >= 2 );    // Combineのときは2以上
        double nres_add = 0;
        for(unsigned int ino=0;ino<nno;ino++)
        {
            double force[2];
            nsf_v.GetValue(ino,force);
            vec_res.AddValue(ino,0,force[0]);
            vec_res.AddValue(ino,1,force[1]);
            nres_add += force[0]*force[0]+force[1]*force[1];
        }
//        std::cout << res << " " << nres_add << std::endl;
    }
	const double norm_res = pLS->FinalizeMarge();
	pPrec->SetValue((*pLS).m_ls);
	return norm_res;
}

bool CEqnSystem_Fluid2D::UpdateDomain_Field(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
	m_id_press = world.MakeField_FieldElemDim(id_base,2,
        Fem::Field::SCALAR, VELOCITY|ACCELERATION,CORNER);
    assert( world.IsIdField(m_id_press) );
//    std::cout << "press " << m_id_press << std::endl;

	if( m_IsntInterpolationBubble){
		std::cout << "corner intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemDim(id_base,2,
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER);
        m_IsntCombine = false;
	}
	else{
		std::cout << "bubble intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemDim(id_base,2,
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER|BUBBLE);
        m_IsntCombine = true;
	}
//    std::cout << "velo : " << m_id_velo << std::endl;
    assert( world.IsIdField(m_id_velo) );

	{	// 同じ要素配列IDを持つ方程式があったら，それを使う．なければ新規に追加
		std::vector<CEqn_Fluid2D> aEqn_old = m_aEqn;
		m_aEqn.clear();
		const CField& field = world.GetField(m_id_velo);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			unsigned int ieqn0=0;
			for(;ieqn0<aEqn_old.size();ieqn0++){
				if( aEqn_old[ieqn0].GetIdEA() == id_ea ){ 
					const unsigned int ieqn1 = m_aEqn.size();
					m_aEqn.push_back( aEqn_old[ieqn0] );
					m_aEqn[ieqn1].SetIdFieldVelocity(m_id_velo);
					m_aEqn[ieqn1].SetIdFieldPressure(m_id_press);
					break; 
				}
			}
			if( ieqn0 != aEqn_old.size() ){ continue; }
            CEqn_Fluid2D eqn1(id_ea,m_id_velo,m_id_press);
            eqn1.SetRho(m_rho_back);
            eqn1.SetMyu(m_myu_back);
            if( this->m_is_stokes_back ){	eqn1.SetStokes(); }
            else{							eqn1.SetNavierStokes(); }
            m_aEqn.push_back( eqn1 );
		}
	}
/*
    std::cout << "Size Eqn : " << m_aEqn.size() << std::endl;
    for(unsigned int ieqn=0;ieqn<m_aEqn.size();ieqn++){
        const CEqn_Fluid2D& eqn = m_aEqn[ieqn];
        std::cout << ieqn << " " << eqn.GetMyu() << " " << eqn.GetRho() << " " << eqn.IsNavierStokes() << std::endl;
    }
*/
	if( !world.IsIdField(m_id_velo) ) return false;
	this->ClearLinearSystemPreconditioner();
	return true;
}


bool CEqnSystem_Fluid2D::UpdateDomain_FieldVeloPress(
		unsigned int id_base_field_v, unsigned int id_base_field_p, Fem::Field::CFieldWorld& world)
{
	m_id_press = world.MakeField_FieldElemDim(id_base_field_p,2,
        Fem::Field::SCALAR, VELOCITY|ACCELERATION,CORNER);
    assert( world.IsIdField(m_id_press) );
//    std::cout << "press " << m_id_press << std::endl;

    m_IsntCombine = true;
	if( m_IsntInterpolationBubble){
		std::cout << "corner intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemDim(id_base_field_v,2,
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER);
	}
	else{
		std::cout << "bubble intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemDim(id_base_field_v,2,
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER|BUBBLE);
	}
//    std::cout << "velo : " << m_id_velo << std::endl;
    assert( world.IsIdField(m_id_velo) );

	{	// 同じ要素配列IDを持つ方程式があったら，それを使う．なければ新規に追加
		std::vector<CEqn_Fluid2D> aEqn_old = m_aEqn;
		m_aEqn.clear();
		const CField& field = world.GetField(m_id_velo);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			unsigned int ieqn0=0;
			for(;ieqn0<aEqn_old.size();ieqn0++){
				if( aEqn_old[ieqn0].GetIdEA() == id_ea ){ 
					const unsigned int ieqn1 = m_aEqn.size();
					m_aEqn.push_back( aEqn_old[ieqn0] );
					m_aEqn[ieqn1].SetIdFieldVelocity(m_id_velo);
					m_aEqn[ieqn1].SetIdFieldPressure(m_id_press);
					break; 
				}
			}
			if( ieqn0 != aEqn_old.size() ){ continue; }
            CEqn_Fluid2D eqn1(id_ea,m_id_velo,m_id_press);
            eqn1.SetRho(m_rho_back);
            eqn1.SetMyu(m_myu_back);
            if( this->m_is_stokes_back ){
//                std::cout << "Stokes" << std::endl;
                eqn1.SetStokes();
            }
            else{
//                std::cout << "Navier-Stokes" << std::endl;
                eqn1.SetNavierStokes();
            }
            m_aEqn.push_back( eqn1 );
		}
	}
	if( !world.IsIdField(m_id_velo) ) return false;
	this->ClearLinearSystemPreconditioner();
	return true;
}

bool CEqnSystem_Fluid2D::UpdateDomain_FieldElemAry(unsigned int id_base,unsigned int id_ea, Fem::Field::CFieldWorld& world)
{
	m_id_press = world.MakeField_FieldElemAry(id_base,id_ea, 
        Fem::Field::SCALAR, VELOCITY|ACCELERATION,CORNER);

	if( m_IsntInterpolationBubble){
		std::cout << "not bubble intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemAry(id_base,id_ea, 
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER);
        m_IsntCombine = false;
	}
	else{
		std::cout << "bubble intp" << std::endl;
		m_id_velo  = world.MakeField_FieldElemAry(id_base,id_ea, 
            Fem::Field::VECTOR2,VELOCITY|ACCELERATION,CORNER|BUBBLE);
        m_IsntCombine = true;
	}
	
	{	// 同じ要素配列IDを持つ方程式があったら，それを使う．なければ新規に追加
		std::vector<CEqn_Fluid2D> aEqn_old = m_aEqn;
		m_aEqn.clear();
		const CField& field = world.GetField(m_id_velo);
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			unsigned int ieqn0=0;
			for(;ieqn0<aEqn_old.size();ieqn0++){
				if( aEqn_old[ieqn0].GetIdEA() == id_ea ){ 
					const unsigned int ieqn1 = m_aEqn.size();
					m_aEqn.push_back( aEqn_old[ieqn0] );
					m_aEqn[ieqn1].SetIdFieldVelocity(m_id_velo);
					m_aEqn[ieqn1].SetIdFieldPressure(m_id_press);
					break; 
				}
			}
			if( ieqn0 != aEqn_old.size() ){ continue; }
            CEqn_Fluid2D eqn1(id_ea,m_id_velo,m_id_press);
            eqn1.SetRho(m_rho_back);
            eqn1.SetMyu(m_myu_back);
            if( this->m_is_stokes_back ){ eqn1.SetStokes(); }
            else{ eqn1.SetNavierStokes(); }
            m_aEqn.push_back( eqn1 );
		}
	}

	this->ClearLinearSystemPreconditioner();
	return true;
}

bool CEqnSystem_Fluid2D::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,-1) );
	return true;
}


unsigned int CEqnSystem_Fluid2D::AddFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world );
}

unsigned int CEqnSystem_Fluid2D::AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
	const unsigned int id_field = world.GetPartialField(m_id_velo, aIdEA );
	if( id_field == 0 ) return 0;
	assert( world.IsIdField(id_field) );
	CField& field = world.GetField(id_field);
  Fem::Field::SetFieldValue_Constant(id_field,0,Fem::Field::VELOCITY,world,0);
  Fem::Field::SetFieldValue_Constant(id_field,1,Fem::Field::VELOCITY,world,0);
	m_aIdFixField.push_back( std::make_pair(id_field,-1) );
	return id_field;
}


bool CEqnSystem_Fluid2D::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{
//	std::cout << "InitializeLinearSystem" << std::endl;
	if( pLS!=0 || pPrec!=0 ) ClearLinearSystemPreconditioner();

	// 連立一次方程式クラスの設定
	if( m_IsntCombine ){
		pLS = new CLinearSystem_Field;
		pLS->AddPattern_Field(m_id_velo,world);	// val_fieldからできる全体剛性行列を追加する
		pLS->AddPattern_Field(m_id_press,m_id_velo,world);
	}
    else{   // 流速-圧力結合自由度で行列を作る
        assert( this->m_IsntInterpolationBubble );
		pLS = new CLinearSystem_Field;
		pLS->AddPattern_CombinedField(m_id_velo,m_id_press,world);
	}
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		const unsigned int id_field = m_aIdFixField[idf].first; 
		const int idof = m_aIdFixField[idf].second; 
		if( idof == -1 ){
			pLS->SetFixedBoundaryCondition_Field( id_field, world ); // bc0_fieldを固定境界条件に設定	
		}
		else{
			pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_fieldを固定境界条件に設定	
		}
	}
	
	// 前処理クラスの作成
	assert( pPrec == 0 );
    if( m_IsntCombine ){    // 流速-圧力を分離する
        pPrec = new LsSol::CPreconditioner_ILU( (*pLS).m_ls, 1 );
	}
    else{   // 流速-圧力を結合させる
        assert( this->m_IsntInterpolationBubble );
        pPrec = new LsSol::CPreconditioner_ILU;
//        (*(LsSol::CPreconditioner_ILU*)pPrec).SetOrdering(false);
        (*(LsSol::CPreconditioner_ILU*)pPrec).SetFillInLevel(0);
//		(*(CPreconditioner_ILU*)pPrec).SetOrdering(false);
//		(*(CPreconditioner_ILU*)pPrec).SetFillInLevel(0);
		pPrec->SetLinearSystem((*pLS).m_ls);
	}

	return true;
}

void CEqnSystem_Fluid2D::SetIsStationary(bool is_stat)
{
	if( this->m_IsStationary == is_stat ) return;
	this->m_IsStationary = is_stat;
	this->m_is_cleared_value_ls = true;
	this->m_is_cleared_value_prec = true;
}

bool CEqnSystem_Fluid2D::Solve(Fem::Field::CFieldWorld& world)
{
	if( pLS == 0 || pPrec == 0 ) this->InitializeLinearSystem(world);
    double res = this->MakeLinearSystem(world);
	{
		double conv_ratio = 1.0e-5;
		unsigned int max_iter = 100;
		bool is_asym;
		this->EqnationProperty(is_asym);
		if( is_asym ){
			assert( !this->m_IsStationary );
            LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
			LsSol::Solve_PBiCGSTAB(conv_ratio,max_iter,lsp);
		}
		else{
            LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
			LsSol::Solve_PCG(conv_ratio,max_iter,lsp);	// Solve with Preconditioned Conjugate Gradient
		}
		this->m_aItrNormRes.clear();
		m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
//        std::cout << max_iter << " " << conv_ratio << std::endl;
	}

	if( this->m_IsStationary ){
		pLS->UpdateValueOfField(m_id_velo, world,VELOCITY);
		pLS->UpdateValueOfField(m_id_press,world,VELOCITY);
	}
	else{
		pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_velo, world,ACCELERATION,true);
		pLS->UpdateValueOfField_Newmark(m_gamma_newmark,m_dt,m_id_press,world,ACCELERATION,true);
	}
	return true;
}

bool CEqnSystem_Fluid2D::ClearFixElemAry(
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


void CEqnSystem_Fluid2D::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}






////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////









CEqn_Fluid3D::CEqn_Fluid3D(unsigned int id_base, Fem::Field::CFieldWorld& world)
: m_myu(1.0), m_rho(1.0), m_g_x(0.0), m_g_y(0.0), m_g_z(0.0)
//m_dt(0.5), m_gamma_newmark(0.5)
{
	this->SetDomain(id_base,world);
}

CEqn_Fluid3D::CEqn_Fluid3D()
: m_myu(1.0), m_rho(1.0), m_g_x(0.0), m_g_y(0.0), m_g_z(0.0)
//m_dt(0.5), m_gamma_newmark(0.5)
{}

double CEqn_Fluid3D::MakeLinearSystem(const Fem::Field::CFieldWorld& world)
{	
	if( pLS==0 || pPrec==0 ) this->InitializeLinearSystem(world);
	// 連立一次方程式を作る
	pLS->InitializeMarge();	// 連立一次方程式を初期化する(0クリア)
	Fem::Eqn::AddLinSys_Stokes3D_Static(
		m_myu, 
		m_rho, m_g_x, m_g_y, m_g_z,
		*pLS,
		m_id_velo,m_id_press,world);
	const double norm_res = pLS->FinalizeMarge();
	pPrec->SetValue((*pLS).m_ls);
	return norm_res;
}

bool CEqn_Fluid3D::SetDomain(unsigned int id_base, Fem::Field::CFieldWorld& world)
{
	bool is_mixed = false;
	Field::ELEM_TYPE elem_type = (Field::ELEM_TYPE)0;
	{
        assert( world.IsIdField(id_base) );
		const CField& field_base = world.GetField(id_base);
		const std::vector<unsigned int>& aIdEA = field_base.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea = aIdEA[iiea];
			CElemAry& ea = world.GetEA(id_ea);
			if( ea.ElemType() == TET ){
				if( elem_type == 0 ){ elem_type = TET; }
				else if( elem_type == HEX ){ is_mixed = true; break; }
			}
			else if( ea.ElemType() == HEX ){
				if( elem_type == 0 ){ elem_type = HEX; }
				else if( elem_type == TET ){ is_mixed = true; break; }
			}
		}
	}
	if( is_mixed ){
		assert(0);
		std::cout << "Not Implimented" << std::endl;
		m_id_press = 0;
		m_id_velo = 0;
		return false;
	}
	if( elem_type == TET ){
		m_id_press = world.MakeField_FieldElemDim(id_base,3,
            Fem::Field::SCALAR, VELOCITY|ACCELERATION,CORNER);
	}
	else if( elem_type == HEX ){
		m_id_press = world.MakeField_FieldElemDim(id_base,3,
            Fem::Field::SCALAR, VELOCITY|ACCELERATION,BUBBLE);
	}
	m_id_velo  = world.MakeField_FieldElemDim(id_base,3,
        Fem::Field::VECTOR3,VELOCITY|ACCELERATION,CORNER);
	this->ClearLinearSystemPreconditioner();
	return true;
}

bool CEqn_Fluid3D::AddFixField(const unsigned int id_field, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdField( id_field ) ) return false;
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	return true;
}

unsigned int CEqn_Fluid3D::AddFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world, int idof)
{
	if( !world.IsIdEA( id_ea ) ) return 0;
	std::vector<unsigned int> aIdEA;
	aIdEA.push_back(id_ea);
	return this->AddFixElemAry( aIdEA, world, idof );
}

unsigned int CEqn_Fluid3D::AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof)
{
	for(unsigned int iid_ea=0;iid_ea<aIdEA.size();iid_ea++){
		if( !world.IsIdEA( aIdEA[iid_ea] ) ) return 0;
	}
    if( aIdEA.size() == 0 ) return 0;
	const unsigned int id_field = world.GetPartialField(m_id_velo, aIdEA );
	if( id_field == 0 ) return 0;
	assert( world.IsIdField(id_field) );
//	CField& field = world.GetField(id_field);
  Fem::Field::SetFieldValue_Constant(id_field,0,Fem::Field::VELOCITY,world,0);
  Fem::Field::SetFieldValue_Constant(id_field,1,Fem::Field::VELOCITY,world,0);
  Fem::Field::SetFieldValue_Constant(id_field,2,Fem::Field::VELOCITY,world,0);
	m_aIdFixField.push_back( std::make_pair(id_field,idof) );
	return id_field;
}


bool CEqn_Fluid3D::InitializeLinearSystem(const Fem::Field::CFieldWorld& world)
{
	if( pLS!=0 || pPrec!=0 ) ClearLinearSystemPreconditioner();
	// 連立一次方程式クラスの設定
	pLS = new CLinearSystem_Field;
	pLS->AddPattern_Field(m_id_velo,world);	// val_fieldからできる全体剛性行列を追加する
	pLS->AddPattern_Field(m_id_press,m_id_velo,world);
	for(unsigned int idf=0;idf<m_aIdFixField.size();idf++){
		const unsigned int id_field = m_aIdFixField[idf].first; 
		const int idof = m_aIdFixField[idf].second;
		if( idof == -1 ){
			pLS->SetFixedBoundaryCondition_Field( id_field, world ); // bc0_fieldを固定境界条件に設定
		}
		else{
			pLS->SetFixedBoundaryCondition_Field( id_field, idof, world ); // bc0_fieldを固定境界条件に設定	
		}
	}
	
	// 前処理クラスの作成
	assert( pPrec == 0 );
    pPrec = new LsSol::CPreconditioner_ILU((*pLS).m_ls);

	return true;
}
			
bool CEqn_Fluid3D::Solve(Fem::Field::CFieldWorld& world)
{
	if( pLS == 0 || pPrec == 0 ) this->InitializeLinearSystem(world);
	this->MakeLinearSystem(world);
	{
		double conv_ratio = 1.0e-6;
		unsigned int max_iter = 1000;
        LsSol::CLinearSystemPreconditioner lsp((*pLS).m_ls,*pPrec);
		LsSol::Solve_PCG(conv_ratio,max_iter,lsp);	// Solve with Preconditioned Conjugate Gradient
//		Fem::Sol::Solve_CG(conv_ratio,max_iter,ls);		// Solve with Conjugate Gradient
//		std::cout << max_iter << " " << conv_ratio << std::endl;
		m_aItrNormRes.clear();
		m_aItrNormRes.push_back( std::make_pair(max_iter,conv_ratio) );
	}

	pLS->UpdateValueOfField(m_id_velo, world,VELOCITY);
	pLS->UpdateValueOfField(m_id_press,world,VELOCITY);

	return true;
}


bool CEqn_Fluid3D::ClearFixElemAry(
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


void CEqn_Fluid3D::ClearFixElemAry()
{
	m_aIdFixField.clear();
	this->ClearLinearSystem();
}

