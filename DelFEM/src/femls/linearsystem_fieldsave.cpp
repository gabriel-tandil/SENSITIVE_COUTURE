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
// LinearSystem.cpp : 連立一次方程式クラス(CLinearSystem.h)の実装
////////////////////////////////////////////////////////////////


#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#ifndef for 
#define for if(0); else for
#endif

#include "math.h"

#include "delfem/field_world.h"

#include "delfem/femls/linearsystem_fieldsave.h"

#include "delfem/indexed_array.h"
#include "delfem/vector3d.h"
#include "delfem/quaternion.h"

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/bcflag_blk.h"

using namespace MatVec;
using namespace Fem::Ls;
using namespace Fem::Field;

CLinearSystem_Save::~CLinearSystem_Save(){
	unsigned int nlss = m_Force.size();
	for(unsigned int ilss=0;ilss<nlss;ilss++){
	for(unsigned int jlss=0;jlss<nlss;jlss++){
		if( m_Matrix_Boundary[ilss][jlss] != 0 ) delete m_Matrix_Boundary[ilss][jlss];
	}
	}
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		if( m_Force[ilss] != 0 ) delete m_Force[ilss];
	}
}


void CLinearSystem_Save::AddBoundaryMatrixForce(unsigned int id_field, 
												const Fem::Field::CFieldWorld& world, 
												unsigned int nlss0, unsigned int nlss1)
{
	std::cout << "CLinearSystem_Save::Add Boundary Matrix Force " << nlss0 << " " << nlss1 << std::endl;
	m_Matrix_Boundary.resize(nlss1);
	for(unsigned int ilss=0;ilss<nlss1;ilss++){
		m_Matrix_Boundary[ilss].resize(nlss1,0);
	}

	m_Force.resize(nlss1);
	for(unsigned int ilss=nlss0;ilss<nlss1;ilss++){
		const unsigned int nblk = m_ls.m_Residual[ilss]->NBlk();
		const unsigned int len = m_ls.m_Residual[ilss]->Len();
		m_Force[ilss] = new CVector_Blk(nblk,len);
	}

	for(unsigned int ilss=0;ilss<nlss1;ilss++){
	for(unsigned int jlss=0;jlss<nlss1;jlss++){
		if( ilss < nlss0 && jlss < nlss0 ) continue; 
		if( ilss == jlss ){
			if( m_ls.m_Matrix_Dia[ilss] == 0 ) continue;
			const unsigned int nblk = m_ls.m_Matrix_Dia[ilss]->NBlkMatCol();
			const unsigned int len  = m_ls.m_Matrix_Dia[ilss]->LenBlkCol();
            this->m_Matrix_Boundary[ilss][jlss] = new MatVec::CMat_BlkCrs(nblk,len, nblk,len);
		}
		else{
			if( m_ls.m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
			const unsigned int nblk_col = m_ls.m_Matrix_NonDia[ilss][jlss]->NBlkMatCol();
			const unsigned int nblk_row = m_ls.m_Matrix_NonDia[ilss][jlss]->NBlkMatRow();
			const unsigned int len_col  = m_ls.m_Matrix_NonDia[ilss][jlss]->LenBlkCol();
			const unsigned int len_row  = m_ls.m_Matrix_NonDia[ilss][jlss]->LenBlkRow();
			this->m_Matrix_Boundary[ilss][jlss] = new CMat_BlkCrs(nblk_col,len_col, nblk_row,len_row);
		}
	}
	}
}

// fieldで初期化する、fieldの中の非ゼロパターンを作る
bool CLinearSystem_Save::AddPattern_Field(const unsigned int id_field, 
										  const Fem::Field::CFieldWorld& world)
{
	const unsigned int nlss0 = this->GetNLynSysSeg();
	assert( m_Force.size() == nlss0 );
	CLinearSystem_Field::AddPattern_Field(id_field,world);
	const unsigned int nlss1 = this->GetNLynSysSeg();
	assert( nlss1 >= nlss0 );
	this->AddBoundaryMatrixForce(id_field,world,nlss0,nlss1);
	return true;
}

// fieldとfield2がパターンが同じだとして，ブロックが結合された一つの行列を作る
bool CLinearSystem_Save::AddPattern_CombinedField(unsigned id_field, unsigned int id_field2, 
												  const Fem::Field::CFieldWorld& world)
{
	const unsigned int nlss0 = this->GetNLynSysSeg();
	assert( m_Force.size() == nlss0 );
	CLinearSystem_Field::AddPattern_CombinedField(id_field,id_field2,world);
	const unsigned int nlss1 = this->GetNLynSysSeg();
	assert( nlss1 >= nlss0 );
	this->AddBoundaryMatrixForce(id_field,world,nlss0,nlss1);
	return true;
}

// fieldで初期化する、fieldとfield-field2の中の非ゼロパターンを作る
bool CLinearSystem_Save::AddPattern_Field(unsigned int id_field, unsigned int id_field2, 
										  const Fem::Field::CFieldWorld& world)
{
	const unsigned int nlss0 = this->GetNLynSysSeg();
	assert( m_Force.size() == nlss0 );
	CLinearSystem_Field::AddPattern_Field(id_field,id_field2,world);
	const unsigned int nlss1 = this->GetNLynSysSeg();
	assert( nlss1 >= nlss0 );
	this->AddBoundaryMatrixForce(id_field,world,nlss0,nlss1);
	return true;
}



// 境界行列をゲットする
CMat_BlkCrs& CLinearSystem_Save::GetMatrix_Boundary(
		unsigned int id_field_col, Fem::Field::ELSEG_TYPE elseg_type_col,
		unsigned int id_field_row, Fem::Field::ELSEG_TYPE elseg_type_row,
		const Fem::Field::CFieldWorld& world)
{
	int ils_col = FindIndexArray_Seg(id_field_col,elseg_type_col,world);
    if( ils_col < 0 || ils_col >= (int)m_aSegField.size() ){ assert(0); throw 0; }
	int ils_row = FindIndexArray_Seg(id_field_row,elseg_type_row,world);
    if( ils_row < 0 || ils_row >= (int)m_aSegField.size() ){ assert(0); throw 0; }
	assert( (unsigned int)ils_col < m_Matrix_Boundary.size() );
	assert( (unsigned int)ils_row < m_Matrix_Boundary[ils_col].size() );
	return *m_Matrix_Boundary[ils_col][ils_row];
}

CVector_Blk& CLinearSystem_Save::GetForce(
	unsigned int id_field, Fem::Field::ELSEG_TYPE elseg_type,
	const Fem::Field::CFieldWorld& world)
{
	int ilss = FindIndexArray_Seg(id_field,elseg_type,world);
    if( ilss < 0 || ilss >= (int)m_aSegField.size() ){ assert(0); throw 0; }
	assert( (unsigned int)ilss < m_Force.size() );
	return *m_Force[ilss];
}

// マージ前の初期化（基底クラスの隠蔽）
void CLinearSystem_Save::InitializeMarge()
{
	CLinearSystem_Field::InitializeMarge();	// 基底クラスを呼ぶ
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	this->MakeBoundaryPattern();
	assert( m_Matrix_Boundary.size() == nlss );

	for(unsigned int ilss=0;ilss<nlss;ilss++){
		m_Force[ilss]->SetVectorZero();
	}
	
	for(unsigned int ilss=0;ilss<nlss;ilss++){
	for(unsigned int jlss=0;jlss<nlss;jlss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
		m_Matrix_Boundary[ilss][jlss]->SetZero();
	}
	}
}

// マージ後の処理（基底クラスの隠蔽）
double CLinearSystem_Save::FinalizeMarge(){
	std::cout << "CLinearSystem_Save::FinalizeMearge()" << std::endl;
	CLinearSystem_Field::FinalizeMarge();	// 基底クラスを呼ぶ
	// ブロックサイズが１より大きい場合は以下が必要。
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	if( this->m_Matrix_Boundary.size() != nlss ) return true;
	for(unsigned int jlss=0;jlss<nlss;jlss++){
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
		m_Matrix_Boundary[ilss][jlss]->SetBoundaryConditionInverse_Colum(m_ls.GetBCFlag(jlss));//*m_ls.m_BCFlag[jlss]);
//		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Colum(*m_BCFlag[jlss]);
		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Row(m_ls.GetBCFlag(ilss));//*m_ls.m_BCFlag[ilss]);
	}
	}
	return 0.0;
}

double CLinearSystem_Save::MakeResidual(const Fem::Field::CFieldWorld& world)
{
	const unsigned int nseg = this->m_aSegField.size();
	if( nseg == 0 )	return 0.0;

	////////////////////////////////
	// updateに値をセットする
	for(unsigned  int ilss=0;ilss<nseg;ilss++){
		const CLinSysSeg_Field& ls = this->m_aSegField[ilss];
		CVector_Blk* update = this->m_ls.m_Update[ilss];
		unsigned int ilen1 = 0;
		{
			const Fem::Field::CField& field = world.GetField(ls.id_field);
			ilen1 = field.GetNLenValue();
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			unsigned int id_ns_val = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			na.GetValueFromNodeSegment(id_ns_val,*update);
		}
		if( world.IsIdField(ls.id_field2) ){
			const Fem::Field::CField& field = world.GetField(ls.id_field2);
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			unsigned int id_ns_val = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			na.GetValueFromNodeSegment(id_ns_val,*update,ilen1);
		}
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*m_ls.m_Residual[iseg]) = (*m_Force[iseg]);
	}

	for(unsigned int iseg=0;iseg<nseg;iseg++){
	for(unsigned int jseg=0;jseg<nseg;jseg++){
		if( m_Matrix_Boundary[iseg][jseg] == 0 ) continue;
		m_Matrix_Boundary[iseg][jseg]->MatVec( -1.0, (*m_ls.m_Update[jseg]), 1.0, (*m_ls.m_Residual[iseg]), true );
	}
	}
	return sqrt(this->DOT(-1,-1));
}

bool CLinearSystem_Save::UpdateValueOfField( 
		unsigned int id_field, Fem::Field::CFieldWorld& world,FIELD_DERIVATION_TYPE fdt)
{
	if( !world.IsIdField(id_field) ){ assert(0); return false; }
	const CField& field = world.GetField(id_field);

	{	// コーナー節点について値を更新
		unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		unsigned int id_ns_val = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
		if( id_na_val != 0 ){
			const int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
            assert( ilss0 >= 0 && ilss0 < (int)this->GetNLynSysSeg() );
			const CLinSysSeg_Field& lss0 = this->m_aSegField[ilss0];
			CVector_Blk& upd = this->m_ls.GetVector(-2,ilss0);  // 固定境界条件をセットするので，constにできない
			const CBCFlag& bc_flag = m_ls.GetBCFlag(ilss0);//this->m_ls.m_BCFlag[ilss0];
//			const Fem::Field::CField& field = world.GetField(id_field);
			CNodeAry& na = world.GetNA(id_na_val);
			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
			if( lss0.id_field == id_field ){
				{	// Updateに固定境界条件の値をセット
					const unsigned int nblk = na.Size();
					const unsigned int nlen = ns.Length();
					assert( nblk == bc_flag.NBlk() );
					double* val = new double [nlen];
					for(unsigned int iblk=0;iblk<nblk;iblk++){
						ns.GetValue(iblk,val);
						for(unsigned int ilen=0;ilen<nlen;ilen++){
							if( bc_flag.GetBCFlag(iblk,ilen) == 0 ) continue;
							upd.SetValue(iblk,ilen,val[ilen]);
						}
					}
					delete[] val;
				}
				na.SetValueToNodeSegment(id_ns_val,upd);
			}
			else{
				assert( lss0.id_field2 == id_field );
				const Fem::Field::CField& field1 = world.GetField(lss0.id_field);
				const unsigned int ilen1 = field1.GetNLenValue();
				{	// Updateに固定境界条件の値をセット
					const unsigned int nblk = na.Size();
					const unsigned int nlen = ns.Length();
					assert( nblk == bc_flag.NBlk() );
                    assert( (int)(nlen+ilen1) == bc_flag.LenBlk() );
					double* val = new double [nlen];
					for(unsigned int iblk=0;iblk<nblk;iblk++){
						ns.GetValue(iblk,val);
						for(unsigned int ilen=0;ilen<nlen;ilen++){
							if( bc_flag.GetBCFlag(iblk,ilen+ilen1) == 0 ) continue;
							upd.SetValue(iblk,ilen+ilen1,val[ilen]);
						}
					}
					delete[] val;
				}
				na.SetValueToNodeSegment(id_ns_val,upd,ilen1);
			}
		}
	}
	{	// 辺節点について値を更新
		unsigned int id_na_val = field.GetNodeSegInNodeAry(EDGE).id_na_va;
		unsigned int id_ns_val = field.GetNodeSegInNodeAry(EDGE).id_ns_va;
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,EDGE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na = world.GetNA(id_na_val);
			na.AddValueToNodeSegment(id_ns_val,upd,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
//			ns.is_updated = true;
		}
	}
	{	// バブル節点について値を更新
		unsigned int id_na_val = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		unsigned int id_ns_val = field.GetNodeSegInNodeAry(BUBBLE).id_ns_va;
		if( id_na_val != 0 ){
            const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,EDGE,world);
			const CVector_Blk& upd = m_ls.GetVector(-2,ilss0);
			CNodeAry& na = world.GetNA(id_na_val);
			na.AddValueToNodeSegment(id_ns_val,upd,1.0);
//			CNodeAry::CNodeSeg& ns = na.GetSeg(id_ns_val);
//			ns.is_updated = true;
		}
	}
	return true;
}

void CLinearSystem_Save::MakeBoundaryPattern()
{
	std::cout << "MakeBoundaryPattern " << std::endl;
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	assert( this->m_Force.size() == nlss );
//	assert( this->m_ls.m_BCFlag.size() == nlss );
	////////////////
	if( m_Matrix_Boundary.size() != nlss ){ assert(0); }
	////////////////
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		assert( m_Matrix_Boundary[ilss].size() == nlss );
		for(unsigned int jlss=0;jlss<nlss;jlss++){
			if( ilss == jlss ) continue;
			if( m_ls.m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
			m_Matrix_Boundary[ilss][jlss]->SetPatternBoundary(
                *m_ls.m_Matrix_NonDia[ilss][jlss],
                m_ls.GetBCFlag(ilss),//*m_ls.m_BCFlag[ilss],
                m_ls.GetBCFlag(jlss));//*m_ls.m_BCFlag[jlss]);
		} 
		if( m_ls.m_Matrix_Dia[ilss] == 0 ) continue;
		{	// 対角サブ行列の境界部分のパターンを作る
			m_Matrix_Boundary[ilss][ilss]->SetPatternBoundary(
                *m_ls.m_Matrix_Dia[ilss],
                m_ls.GetBCFlag(ilss),//*m_ls.m_BCFlag[ilss],
                m_ls.GetBCFlag(ilss));//*m_ls.m_BCFlag[ilss]);
		}
	}
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

bool CLinearSystem_SaveDiaM_Newmark::AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	// オーバーロードしている親クラスの関数を呼ぶ
	CLinearSystem_Field::AddPattern_Field(id_field,world);

	const int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
    assert( ilss0 >= 0 && ilss0 < (int)this->GetNLynSysSeg() );
    assert( ilss0 == (int)this->m_DiaMassMatrix.size() );
	const unsigned int len  = this->m_ls.m_Update[ilss0]->Len();
	const unsigned int nblk = this->m_ls.m_Update[ilss0]->NBlk();
	////////////////
	m_Matrix_Boundary.resize(ilss0+1);
    for(unsigned int ilss=0;ilss<=(unsigned int)ilss0;ilss++){
		m_Matrix_Boundary[ilss].resize(ilss0+1,0);
	}
    for(unsigned int ilss=0;ilss<(unsigned int)ilss0;ilss++){
		if( m_ls.m_Matrix_NonDia[ilss][ilss0] != 0 ){
			this->m_Matrix_Boundary[ilss][ilss0] = new CMat_BlkCrs(*m_ls.m_Matrix_NonDia[ilss][ilss0],false,true);
		}
		if( m_ls.m_Matrix_NonDia[ilss0][ilss] != 0 ){
			this->m_Matrix_Boundary[ilss0][ilss] = new CMat_BlkCrs(*m_ls.m_Matrix_NonDia[ilss0][ilss],false,true);
		}
	}
	if( m_ls.m_Matrix_Dia[ilss0] != 0 ){
		this->m_Matrix_Boundary[ilss0][ilss0] = new CMat_BlkCrs(nblk,len,nblk,len);
		this->m_Matrix_Boundary[ilss0][ilss0]->SetPatternDia( *m_ls.m_Matrix_Dia[ilss0] );
	}
	////////////////
//	m_DiaMassMatrix.resize(ilss0+1,0);
//	m_DiaMassMatrix[ilss0] = new CDiaMat_Blk(nblk,len);
	////////////////
	m_Force.resize(ilss0+1,0);
	m_Force[ilss0] = new CVector_Blk(nblk,len);
	return true;
}

// マージ前の初期化（基底クラスの隠蔽）
void CLinearSystem_SaveDiaM_Newmark::InitializeMarge()
{	
	CLinearSystem_Field::InitializeMarge();	// 基底クラスを呼ぶ
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();

	assert( m_Force.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		m_Force[ilss]->SetVectorZero();
	}
/*
	assert( m_DiaMassMatrix.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		m_DiaMassMatrix[ilss]->SetZero();
	}
*/
	assert( m_Matrix_Boundary.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
	for(unsigned int jlss=0;jlss<nlss;jlss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
		m_Matrix_Boundary[ilss][jlss]->SetZero();
	}
	}
}

// マージ後の処理（基底クラスの隠蔽）
double CLinearSystem_SaveDiaM_Newmark::FinalizeMarge(){
	std::cout << "CLinearSystem_SaveDiaM::FinalizeMearge()" << std::endl;
	CLinearSystem_Field::FinalizeMarge();	// 基底クラスを呼ぶ
	// ブロックサイズが１より大きい場合は以下が必要。
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	if( this->m_Matrix_Boundary.size() != nlss ) return true;
	for(unsigned int jlss=0;jlss<nlss;jlss++){
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
//		m_Matrix_Boundary[ilss][jlss]->SetBoundaryConditionInverse_Colum(*m_BCFlag[jlss]);
//		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Colum(*m_BCFlag[jlss]);
		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Row(m_ls.GetBCFlag(ilss));//*m_ls.m_BCFlag[ilss]);
	}
	}
	return 0.0;
}

double CLinearSystem_SaveDiaM_Newmark::MakeResidual(const Fem::Field::CFieldWorld& world)
{
	const unsigned int nseg = this->m_aSegField.size();
	if( nseg == 0 )	return 0.0;

	// 外力項をセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*m_ls.m_Residual[iseg]) = (*m_Force[iseg]);
	}

	////////////////////////////////
	// updateに値をセットする
	for(unsigned int ilss=0;ilss<nseg;ilss++){
		const CLinSysSeg_Field& lss = m_aSegField[ilss];
		CVector_Blk* update = this->m_ls.m_Update[ilss];
		{
			const CField& field = world.GetField(lss.id_field);
			unsigned int id_ns_u   = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			unsigned int id_ns_v   = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			na.GetValueFromNodeSegment(             id_ns_u,*update);
			na.AddValueFromNodeSegment(dt*(1-gamma),id_ns_v,*update);
		}
	}
//	this->GCMV(-1, -2, 1, -1);

	// 行列ベクトル積により残差を計算する
	for(unsigned int iseg=0;iseg<nseg;iseg++){
	for(unsigned int jseg=0;jseg<nseg;jseg++){
		if( m_Matrix_Boundary[iseg][jseg] == 0 ) continue;
		m_Matrix_Boundary[iseg][jseg]->MatVec( -1.0, *m_ls.m_Update[jseg], 1.0, *m_ls.m_Residual[iseg], true );
	}
	}
	
/*
	////////////////////////////////
	// updateにMを掛けるため加速度値をセットする
	for(unsigned int ilss=0;ilss<nseg;ilss++){
		const CLinSysSeg_Field& ls = this->m_aSegField[ilss];
		CVector_Blk* update = this->m_Update[ilss];
		unsigned int ilen1 = 0;
		{
			const Fem::Field::CField& field = world.GetField(ls.id_field);
			ilen1 = field.GetNLenValue();
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_ve = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			na.GetValueFromNodeSegment(id_ns_ve,*update);
		}
		if( world.IsIdField(ls.id_field2) ){
			const Fem::Field::CField& field = world.GetField(ls.id_field2);
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_ve = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			na.GetValueFromNodeSegment(id_ns_ve,*update,ilen1);
		}
	}
	// 行列ベクトル積により残差を計算する
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		m_DiaMassMatrix[iseg]->MatVec(-1.0, *m_Update[iseg], 1.0, *m_Residual[iseg] );
	}
	*/
	
	// 残差への境界条件のセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
//		m_ls.m_BCFlag[iseg]->SetZeroToBCDof(*m_ls.m_Residual[iseg]);
		m_ls.GetBCFlag(iseg).SetZeroToBCDof(*m_ls.m_Residual[iseg]);
	}

	return sqrt(this->DOT(-1,-1));
}

bool CLinearSystem_SaveDiaM_Newmark::UpdateValueOfField(
		unsigned int id_field, Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
	if( !world.IsIdField(id_field) ) return false;
	CField& field = world.GetField(id_field);

	{
		unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		if( id_na_val != 0 ){
			const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
			CVector_Blk* pUpdate = this->m_ls.m_Update[ilss0];
			assert( pUpdate != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			const unsigned int nblk = na.Size();
			CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true,world,VALUE);
			CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true,world,VELOCITY);
			assert( ns_v.Length() == ns_u.Length() );
			const unsigned int nlen = ns_u.Length();
			const CBCFlag& bc_flag = m_ls.GetBCFlag(ilss0);//this->m_ls.m_BCFlag[ilss0];
			assert( nblk == bc_flag.NBlk() );
            assert( (int)nlen == bc_flag.LenBlk() );
			double* velo0 = new double [nlen];
			double* velo1 = new double [nlen];
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				ns_v.GetValue(iblk,velo0);
				ns_u.GetValue(iblk,velo1);
				for(unsigned int ilen=0;ilen<nlen;ilen++){
					if( bc_flag.GetBCFlag(iblk,ilen) == 0 ){
						const double velo1 = pUpdate->GetValue(iblk,ilen);
						// 値を更新
						ns_u.AddValue(iblk,ilen,dt*(1-gamma)*velo0[ilen]+dt*gamma*velo1);
					}
					else{
						pUpdate->SetValue(iblk,ilen,velo0[ilen]);
//						pUpdate->SetValue(iblk,ilen,velo1[ilen]);
					}
				}
			}
			delete[] velo0;
			const unsigned int id_ns_v = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			na.SetValueToNodeSegment(id_ns_v,*pUpdate);
		}
	}
	return true;
}

CDiaMat_Blk& CLinearSystem_SaveDiaM_Newmark::GetDiaMassMatrix(
	unsigned int id_field, 
	Fem::Field::ELSEG_TYPE elseg_type, 
	const Fem::Field::CFieldWorld& world)
{
	int ilss = this->FindIndexArray_Seg(id_field,elseg_type,world);
    if( ilss == -1 ){ assert(0); throw 0; }
	return *m_DiaMassMatrix[ilss];
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CLinearSystem_SaveDiaM_NewmarkBeta::~CLinearSystem_SaveDiaM_NewmarkBeta(){
	unsigned int nlss = this->m_DiaMassMatrix.size();
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		if( m_DiaMassMatrix[ilss]!= 0 ){ delete m_DiaMassMatrix[ilss]; }
	}
}

void CLinearSystem_SaveDiaM_NewmarkBeta::AddBoundaryMatrixForce(unsigned int id_field, 
												const Fem::Field::CFieldWorld& world, 
												unsigned int nlss0, unsigned int nlss1)
{
	std::cout << "Add Boundary Matrix Force " << nlss0 << " " << nlss1 << std::endl;
	assert( nlss0 == this->m_DiaMassMatrix.size() );

	m_Matrix_Boundary.resize(nlss1);
	for(unsigned int ilss=0;ilss<nlss1;ilss++){
		m_Matrix_Boundary[ilss].resize(nlss1,0);
	}

	m_Force.resize(nlss1);
	for(unsigned int ilss=nlss0;ilss<nlss1;ilss++){
		const unsigned int nblk = m_ls.m_Residual[ilss]->NBlk();
		const unsigned int len = m_ls.m_Residual[ilss]->Len();
		m_Force[ilss] = new CVector_Blk(nblk,len);
	}

	m_DiaMassMatrix.resize(nlss1);
	for(unsigned int ilss=nlss0;ilss<nlss1;ilss++){
		const unsigned int nblk = m_ls.m_Residual[ilss]->NBlk();
		const unsigned int len = m_ls.m_Residual[ilss]->Len();
		m_DiaMassMatrix[ilss] = new CDiaMat_Blk(nblk,len);
	}

	for(unsigned int ilss=0;ilss<nlss1;ilss++){
	for(unsigned int jlss=0;jlss<nlss1;jlss++){
		if( ilss < nlss0 && jlss < nlss0 ) continue; 
		if( ilss == jlss ){
			if( m_ls.m_Matrix_Dia[ilss] == 0 ) continue;
			const unsigned int nblk = m_ls.m_Matrix_Dia[ilss]->NBlkMatCol();
			const unsigned int len  = m_ls.m_Matrix_Dia[ilss]->LenBlkCol();
			this->m_Matrix_Boundary[ilss][ilss] = new CMat_BlkCrs(nblk,len,nblk,len);
			this->m_Matrix_Boundary[ilss][ilss]->SetPatternDia( *m_ls.m_Matrix_Dia[ilss] );
		}
		else{
			if( m_ls.m_Matrix_NonDia[ilss][jlss] == 0 ) continue;
//			const unsigned int nblk_col = m_ls.m_Matrix_NonDia[ilss][jlss]->NBlkMatCol();
//			const unsigned int nblk_row = m_ls.m_Matrix_NonDia[ilss][jlss]->NBlkMatRow();
//			const unsigned int len_col  = m_ls.m_Matrix_NonDia[ilss][jlss]->LenBlkCol();
//			const unsigned int len_row  = m_ls.m_Matrix_NonDia[ilss][jlss]->LenBlkRow();
			this->m_Matrix_Boundary[ilss][jlss] = new CMat_BlkCrs(*m_ls.m_Matrix_NonDia[ilss][jlss],false,true);
		}
	}
	}
}

bool CLinearSystem_SaveDiaM_NewmarkBeta::AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	const unsigned int nlss0 = this->GetNLynSysSeg();
	assert( m_Force.size() == nlss0 );
	// オーバーロードしている親クラスの関数を呼ぶ
	CLinearSystem_Field::AddPattern_Field(id_field,world);
	const unsigned int nlss1 = this->GetNLynSysSeg();
	assert( nlss1 >= nlss0 );
	this->AddBoundaryMatrixForce(id_field,world,nlss0,nlss1);
	return true;
}

bool CLinearSystem_SaveDiaM_NewmarkBeta::AddPattern_CombinedField(
		const unsigned int id_field, const unsigned int id_field1, 
		const Fem::Field::CFieldWorld& world)
{
	const unsigned int nlss0 = this->GetNLynSysSeg();
	assert( m_Force.size() == nlss0 );
	// オーバーロードしている親クラスの関数を呼ぶ
	CLinearSystem_Field::AddPattern_CombinedField(id_field,id_field1,world);
	const unsigned int nlss1 = this->GetNLynSysSeg();
	assert( nlss1 >= nlss0 );
	this->AddBoundaryMatrixForce(id_field,world,nlss0,nlss1);
	return true;
}

// マージ前の初期化（基底クラスの隠蔽）
void CLinearSystem_SaveDiaM_NewmarkBeta::InitializeMarge()
{	
	std::cout << "CLinearSystem_SaveDiaM_NewmarkBeta::InitializeMearge()" << std::endl;
	CLinearSystem_Field::InitializeMarge();	// 基底クラスを呼ぶ
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	////////////////
	assert( m_Force.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		m_Force[ilss]->SetVectorZero();
	}
	////////////////
	assert( m_Matrix_Boundary.size() == nlss );
	for(unsigned int ilss=0;ilss<nlss;ilss++){
	for(unsigned int jlss=0;jlss<nlss;jlss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
		m_Matrix_Boundary[ilss][jlss]->SetZero();
	}
	}
	////////////////
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		m_DiaMassMatrix[ilss]->SetZero();
	}
}

// マージ後の処理（基底クラスの隠蔽）
double CLinearSystem_SaveDiaM_NewmarkBeta::FinalizeMarge(){
	std::cout << "CLinearSystem_SaveDiaM_NewmarkBeta::FinalizeMearge()" << std::endl;
	CLinearSystem_Field::FinalizeMarge();	// 基底クラスを呼ぶ
	// ブロックサイズが１より大きい場合は以下が必要。
	unsigned int nlss = CLinearSystem_Field::GetNLynSysSeg();
	if( this->m_Matrix_Boundary.size() != nlss ) return true;
	for(unsigned int jlss=0;jlss<nlss;jlss++){
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		if( m_Matrix_Boundary[ilss][jlss] == 0 ) continue;
//		m_Matrix_Boundary[ilss][jlss]->SetBoundaryConditionInverse_Colum(*m_BCFlag[jlss]);
//		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Colum(*m_BCFlag[jlss]);
		m_Matrix_Boundary[ilss][jlss]->SetBoundaryCondition_Row(m_ls.GetBCFlag(ilss));//*m_ls.m_BCFlag[ilss]);
	}
	}
	return 0.0;
}


double CLinearSystem_SaveDiaM_NewmarkBeta::MakeResidual(const Fem::Field::CFieldWorld& world)
{
	const unsigned int nseg = this->m_aSegField.size();
	if( nseg == 0 )	return 0.0;

	// 外力項をセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		(*m_ls.m_Residual[iseg]) = (*m_Force[iseg]);
	}
	
	////////////////////////////////
	// updateにKを掛けるための値をセットする
	for(unsigned int ilss=0;ilss<nseg;ilss++){
		const CLinSysSeg_Field& ls = this->m_aSegField[ilss];
		CVector_Blk* update = this->m_ls.m_Update[ilss];
		unsigned int ilen1 = 0;
		{
			const Fem::Field::CField& field = world.GetField(ls.id_field);
			ilen1 = field.GetNLenValue();
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_u = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			unsigned int id_ns_v = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			unsigned int id_ns_a = field.GetNodeSegInNodeAry(CORNER).id_ns_ac;
			na.GetValueFromNodeSegment(                 id_ns_u,*update);
			na.AddValueFromNodeSegment(dt,              id_ns_v,*update);
			na.AddValueFromNodeSegment(dt*dt*(0.5-beta),id_ns_a,*update);
		}
		if( world.IsIdField(ls.id_field2) ){
			const Fem::Field::CField& field = world.GetField(ls.id_field2);
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_u = field.GetNodeSegInNodeAry(CORNER).id_ns_va;
			unsigned int id_ns_v = field.GetNodeSegInNodeAry(CORNER).id_ns_ve;
			unsigned int id_ns_a = field.GetNodeSegInNodeAry(CORNER).id_ns_ac;
			na.GetValueFromNodeSegment(                 id_ns_u,*update,ilen1);
			na.AddValueFromNodeSegment(dt,              id_ns_v,*update,ilen1);
			na.AddValueFromNodeSegment(dt*dt*(0.5-beta),id_ns_a,*update,ilen1);
		}
	}

	// 行列ベクトル積により残差を計算する
	for(unsigned int iseg=0;iseg<nseg;iseg++){
	for(unsigned int jseg=0;jseg<nseg;jseg++){
		if( m_Matrix_Boundary[iseg][jseg] == 0 ) continue;
		m_Matrix_Boundary[iseg][jseg]->MatVec( -1.0, *m_ls.m_Update[jseg], 1.0, *m_ls.m_Residual[iseg], true );
	}
	}
/*
	////////////////////////////////
	// updateにMを掛けるため加速度値をセットする
	for(unsigned int ilss=0;ilss<nseg;ilss++){
		const CLinSysSeg_Field& ls = this->m_aSegField[ilss];
		CVector_Blk* update = this->m_Update[ilss];
		unsigned int ilen1 = 0;
		{
			const Fem::Field::CField& field = world.GetField(ls.id_field);
			ilen1 = field.GetNLenValue();
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_a = field.GetNodeSegInNodeAry(CORNER).id_ns_ac;
			na.GetValueFromNodeSegment(id_ns_a,*update);
		}
		if( world.IsIdField(ls.id_field2) ){
			const Fem::Field::CField& field = world.GetField(ls.id_field2);
			unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
			assert( world.IsIdNA(id_na_val) );
			const CNodeAry& na = world.GetNA(id_na_val);
			unsigned int id_ns_a = field.GetNodeSegInNodeAry(CORNER).id_ns_ac;
			na.GetValueFromNodeSegment(id_ns_a,*update,ilen1);
		}
	}

	// 行列ベクトル積により残差を計算する
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		m_DiaMassMatrix[iseg]->MatVec( -1.0, *m_Update[iseg], 1.0, *m_Residual[iseg] );
	}
*/
	// 残差への境界条件のセット
	for(unsigned int iseg=0;iseg<nseg;iseg++){
//		m_ls.m_BCFlag[iseg]->SetZeroToBCDof(*m_ls.m_Residual[iseg]);
		m_ls.GetBCFlag(iseg).SetZeroToBCDof(*m_ls.m_Residual[iseg]);
	}

	return sqrt(this->DOT(-1,-1));
}

bool CLinearSystem_SaveDiaM_NewmarkBeta::UpdateValueOfField(
		unsigned int id_field, Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
//	std::cout << "CLinearSystem_SaveDiaM_NewmarkBeta::UpdateValueOfField" << std::endl;
	if( !world.IsIdField(id_field) ) return false;
	CField& field = world.GetField(id_field);

	assert( fdt == ACCELERATION );
	{
		unsigned int id_na_val = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		unsigned int id_ns_a = field.GetNodeSegInNodeAry(CORNER).id_ns_ac;
		if( id_na_val != 0 ){
			const unsigned int ilss0 = this->FindIndexArray_Seg(id_field,CORNER,world);
			const CLinSysSeg_Field& lss = this->m_aSegField[ilss0];
			CVector_Blk* pUpdate = this->m_ls.m_Update[ilss0];
			assert( pUpdate != 0 );
			CNodeAry& na = world.GetNA(id_na_val);
			const unsigned int nblk = na.Size();
			CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true,world,VALUE);
			CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true,world,VELOCITY);
			CNodeAry::CNodeSeg& ns_a = field.GetNodeSeg(CORNER,true,world,ACCELERATION);
			assert( ns_v.Length() == ns_u.Length() );
			const unsigned int nlen = ns_u.Length();
			const CBCFlag& bc_flag = m_ls.GetBCFlag(ilss0);//this->m_ls.m_BCFlag[ilss0];
			assert( nblk == bc_flag.NBlk() );
			double* velo0 = new double [nlen];
			double* acc0  = new double [nlen];
			if( lss.id_field == id_field ){
				for(unsigned int iblk=0;iblk<nblk;iblk++){
					ns_v.GetValue(iblk,velo0);
					ns_a.GetValue(iblk, acc0);
					for(unsigned int ilen=0;ilen<nlen;ilen++){
						if( bc_flag.GetBCFlag(iblk,ilen) == 0 ){	// 固定境界条件が入ってなかったら
							const double acc1 = pUpdate->GetValue(iblk,ilen);	// 更新後の速度
							// 値を更新
							ns_u.AddValue(iblk,ilen,dt*velo0[ilen]+dt*dt*(0.5-beta)*acc0[ilen]+dt*dt*beta*acc1);
							ns_v.AddValue(iblk,ilen,dt*(1-gamma)*acc0[ilen]+gamma*dt*acc1);
							ns_a.AddValue(iblk,ilen,acc1);
						}
						else{	// 固定境界条件が入ってなかったら
							// Updateに固定境界条件つき速度を入れる
							pUpdate->SetValue(iblk,ilen,acc0[ilen]);
						}
					}
				}
				na.SetValueToNodeSegment(id_ns_a,*pUpdate);
			}
			else{
				assert( lss.id_field2 == id_field );
				unsigned int nlen0;
				{
					const unsigned int id_field0 = lss.id_field;
					const CField& field = world.GetField(id_field0);
					nlen0 = field.GetNLenValue();
				}
				for(unsigned int iblk=0;iblk<nblk;iblk++){
					ns_v.GetValue(iblk,velo0);
					ns_a.GetValue(iblk, acc0);
					for(unsigned int ilen=0;ilen<nlen;ilen++){
						if( bc_flag.GetBCFlag(iblk,ilen+nlen0) == 0 ){	// 固定境界条件が入ってなかったら
							const double acc1 = pUpdate->GetValue(iblk,ilen+nlen0);	// 更新後の速度
							// 値を更新
							ns_u.AddValue(iblk,ilen,dt*velo0[ilen]+dt*dt*(0.5-beta)*acc0[ilen]+dt*dt*beta*acc1);
							ns_v.AddValue(iblk,ilen,dt*(1-gamma)*acc0[ilen]+gamma*dt*acc1);
							ns_a.AddValue(iblk,ilen,acc1);
						}
						else{	// 固定境界条件が入ってなかったら
							// Updateに固定境界条件つき速度を入れる
							pUpdate->SetValue(iblk,ilen+nlen0,acc0[ilen]);
						}
					}
				}
				na.SetValueToNodeSegment(id_ns_a,*pUpdate,nlen0);
			}
			delete[] velo0;
			delete[] acc0;
			// 加速度の更新
		}
	}
	return true;
}


CDiaMat_Blk& CLinearSystem_SaveDiaM_NewmarkBeta::GetDiaMassMatrix(
	unsigned int id_field, 
	Fem::Field::ELSEG_TYPE elseg_type, 
	const Fem::Field::CFieldWorld& world)
{
	int ilss = this->FindIndexArray_Seg(id_field,elseg_type,world);
    if( ilss == -1 ){ assert(0); throw 0; }
	return *m_DiaMassMatrix[ilss];
}





////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CLinearSystem_Eigen::Clear(){
	CLinearSystem_Field::Clear();
	const unsigned int nlss = this->GetNLynSysSeg();
	for(unsigned int ilss=0;ilss<nlss;ilss++){
		delete m_DiaMassMatrix[ilss];
	}
	m_DiaMassMatrix.clear();
}

bool CLinearSystem_Eigen::AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	// オーバーロードしている親クラスの関数を呼ぶ
	CLinearSystem_Field::AddPattern_Field(id_field,world);

	int ilss = this->FindIndexArray_Seg(id_field,CORNER,world);
    assert( ilss >= 0 && ilss < (int)this->GetNLynSysSeg() );
    assert( (unsigned int)ilss == this->m_DiaMassMatrix.size() );
	const unsigned int len = this->m_ls.m_Update[ilss]->Len();
	const unsigned int nblk = this->m_ls.m_Update[ilss]->NBlk();
	m_DiaMassMatrix.resize(ilss+1,0);
	m_DiaMassMatrix[ilss] = new CDiaMat_Blk(nblk,len);
	return true;
}


void CLinearSystem_Eigen::InitializeMarge()
{
	// オーバーロードしている親クラスの関数を呼ぶ
	CLinearSystem_Field::InitializeMarge();
	assert( this->GetNLynSysSeg() == 1 );
	m_DiaMassMatrix[0]->SetZero();
}

CDiaMat_Blk& CLinearSystem_Eigen::GetDiaMassMatrix(
	unsigned int id_field, Fem::Field::ELSEG_TYPE elseg_type, const Fem::Field::CFieldWorld& world)
{
	int ilss = this->FindIndexArray_Seg(id_field,elseg_type,world);
    if( ilss == -1 ){ assert(0); throw 0; }
	return *m_DiaMassMatrix[ilss];
}


bool CLinearSystem_Eigen::MultVecMassDecomp(int ivec){
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	CVector_Blk& vec = m_ls.GetVector(ivec,0);
/*	if( ivec >= 0 ){
		assert( ivec < this->GetTmpVectorArySize() );
		pUpdate = this->m_ls.m_TmpVectorArray[ivec][0];
	}
	else if( ivec == -2 ){
		pUpdate = this->m_ls.m_Update[0];
	}
	else{ assert(0); }*/
	CDiaMat_Blk* dmat = this->m_DiaMassMatrix[0];
	const unsigned int nblk = dmat->NBlk();
	const unsigned int nlen = dmat->LenBlk();
	assert( nlen == 3 );
	if( nlen == 3 ){
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			const double* pValMi = dmat->GetPtrValDia(iblk);
			double val0 = vec.GetValue(iblk,0);
			double val1 = vec.GetValue(iblk,1);
			double val2 = vec.GetValue(iblk,2);
			const double val3 = val0/pValMi[0];
			const double val4 = val1/pValMi[4];
			const double val5 = val2/pValMi[8];
			vec.SetValue(iblk,0,val3);
			vec.SetValue(iblk,1,val4);
			vec.SetValue(iblk,2,val5);
		}
	}
	return true;
}

bool CLinearSystem_Eigen::MultUpdateInvMassDecomp()
{
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	{
		CVector_Blk* pUpdate = this->m_ls.m_Update[0];
		CDiaMat_Blk* dmat = this->m_DiaMassMatrix[0];
		const unsigned int nblk = dmat->NBlk();
		const unsigned int nlen = dmat->LenBlk();
		if( nlen == 1 ){
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				const double* pValMi = dmat->GetPtrValDia(iblk);
				double val = pUpdate->GetValue(iblk,0);
				val *= (*pValMi);
				pUpdate->SetValue(iblk,0,val);
			}
		}
		else if( nlen == 2 ){
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				const double* pValMi = dmat->GetPtrValDia(iblk);
				double val0 = pUpdate->GetValue(iblk,0);
				double val1 = pUpdate->GetValue(iblk,1);
				const double val2 = pValMi[0]*val0+pValMi[1]*val1;
				const double val3 = pValMi[2]*val0+pValMi[3]*val1;
				pUpdate->SetValue(iblk,0,val2);
				pUpdate->SetValue(iblk,1,val3);
			}
		}
		else if( nlen == 3 ){
			for(unsigned int iblk=0;iblk<nblk;iblk++){
				const double* pValMi = dmat->GetPtrValDia(iblk);
				double val0 = pUpdate->GetValue(iblk,0);
				double val1 = pUpdate->GetValue(iblk,1);
				double val2 = pUpdate->GetValue(iblk,2);
				const double val3 = pValMi[0]*val0+pValMi[1]*val1+pValMi[2]*val2;
				const double val4 = pValMi[3]*val0+pValMi[4]*val1+pValMi[5]*val2;
				const double val5 = pValMi[6]*val0+pValMi[7]*val1+pValMi[8]*val2;
				pUpdate->SetValue(iblk,0,val3);
				pUpdate->SetValue(iblk,1,val4);
				pUpdate->SetValue(iblk,2,val5);
			}
		}
	}
	return true;
}
/*
void CLinearSystem_Eigen::RemoveConstant(int iv1)
{
	const unsigned int nseg = this->m_aSeg.size();
	if( nseg == 0 )	return;

	std::vector< CVector_Blk* >* p_vec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = &m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) p_vec1 = &this->m_Residual;
	else if( iv1 == -2 ) p_vec1 = &this->m_Update;
	else{ assert(0); }

	assert( p_vec1->size() == nseg );

	unsigned int ntotdf = 0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		CVector_Blk& vec = (*(*p_vec1)[iseg]);
		ntotdf += vec.BlkLen()*vec.BlkVecLen();
	}
	double sum = 0;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		CVector_Blk& vec = (*(*p_vec1)[iseg]);
		const unsigned int nblk = vec.BlkVecLen();
		const unsigned int nlen = vec.BlkLen();
		for(unsigned int iblk=0;iblk<nblk;iblk++){
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			sum += vec.GetValue(iblk,ilen);
		}
		}
	}
	const double ave = sum / (double)ntotdf;
	for(unsigned int iseg=0;iseg<nseg;iseg++){
		CVector_Blk& vec = (*(*p_vec1)[iseg]);
		const unsigned int nblk = vec.BlkVecLen();
		const unsigned int nlen = vec.BlkLen();
		for(unsigned int iblk=0;iblk<nblk;iblk++){
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			const double val = vec.GetValue(iblk,ilen);
			vec.SetValue(iblk,ilen,val-ave);
		}
		}
	}
	return;
}
*/
void CLinearSystem_Eigen::OffsetDiagonal(double lambda){
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	CMatDia_BlkCrs* dmat = this->m_ls.m_Matrix_Dia[0];
	const unsigned int nblk = dmat->NBlkMatCol();
	const unsigned int nlen = dmat->LenBlkCol();
	for(unsigned int iblk=0;iblk<nblk;iblk++){
		double* pVal = dmat->GetPtrValDia(iblk);
		for(unsigned int ilen=0;ilen<nlen;ilen++){
			pVal[ilen*nlen+ilen] -= lambda;
		}
	}
}

bool CLinearSystem_Eigen::DecompMultMassMatrix()
{
	const unsigned int nlss = this->GetNLynSysSeg();
	assert( nlss == 1 );
	m_DiaMassMatrix[0]->CholeskyDecomp();
	{
		CMatDia_BlkCrs* dmat = this->m_ls.m_Matrix_Dia[0];
		const unsigned int nblk = dmat->NBlkMatCol();
		const unsigned int nlen = dmat->LenBlkCol();
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			const double* pVal_Mi = m_DiaMassMatrix[0]->GetPtrValDia(iblk);
			unsigned int npsup;
			const unsigned int* psup = dmat->GetPtrIndPSuP(iblk,npsup);
			for(unsigned int ipsup=0;ipsup<npsup;ipsup++){
				const unsigned int jblk = psup[ipsup];
				assert( iblk != jblk );
				unsigned int npsup0;
				double* pVal_Cij = dmat->GetPtrValPSuP(iblk,npsup0);
				assert( npsup0 == npsup );
				const double* pVal_Mj = m_DiaMassMatrix[0]->GetPtrValDia(jblk);
				if( nlen == 1 ){
					pVal_Cij[ipsup] = (*pVal_Mi)*pVal_Cij[ipsup]*(*pVal_Mj);
				}
				else if( nlen == 2 ){
					const double C_old[4] = { pVal_Cij[ipsup*4  ], pVal_Cij[ipsup*4+1], 
						                      pVal_Cij[ipsup*4+2], pVal_Cij[ipsup*4+3] };
					for(unsigned int ilen=0;ilen<nlen;ilen++){
					for(unsigned int jlen=0;jlen<nlen;jlen++){
						pVal_Cij[ipsup*nlen*nlen+nlen*ilen+jlen] = 0;
						for(unsigned int klen=0;klen<nlen;klen++){
						for(unsigned int llen=0;llen<nlen;llen++){
							pVal_Cij[ipsup*nlen*nlen+nlen*ilen+jlen]
								+= pVal_Mi[nlen*ilen+klen]*C_old[nlen*klen+llen]*pVal_Mj[nlen*llen+jlen];
						}
						}
					}
					}
				}
				else if( nlen == 3 ){
					const double C_old[9] = 
					{ pVal_Cij[ipsup*9  ], pVal_Cij[ipsup*9+1], pVal_Cij[ipsup*9+2], 
					  pVal_Cij[ipsup*9+3], pVal_Cij[ipsup*9+4], pVal_Cij[ipsup*9+5], 
					  pVal_Cij[ipsup*9+6], pVal_Cij[ipsup*9+7], pVal_Cij[ipsup*9+8] };
					for(unsigned int ilen=0;ilen<nlen;ilen++){
					for(unsigned int jlen=0;jlen<nlen;jlen++){
						pVal_Cij[ipsup*9+3*ilen+jlen] = 0;
						pVal_Cij[ipsup*9+3*ilen+jlen]
							+= pVal_Mi[3*ilen+ilen]*C_old[3*ilen+jlen]*pVal_Mj[3*jlen+jlen];
					}
					}
				}
			}
			{
				double* pVal_C = dmat->GetPtrValDia(iblk);
				if( nlen == 1 ){
					(*pVal_C) = (*pVal_Mi)*(*pVal_C)*(*pVal_Mi);
				}
				else if( nlen == 2 ){
					const double C_old[4] = { pVal_C[0], pVal_C[1], pVal_C[2], pVal_C[3] };
					for(unsigned int ilen=0;ilen<nlen;ilen++){
					for(unsigned int jlen=0;jlen<nlen;jlen++){
						pVal_C[nlen*ilen+jlen] = 0;
						for(unsigned int klen=0;klen<nlen;klen++){
						for(unsigned int llen=0;llen<nlen;llen++){
							pVal_C[nlen*ilen+jlen] 
								+= pVal_Mi[nlen*ilen+klen]*C_old[nlen*klen+llen]*pVal_Mi[nlen*llen+jlen];
						}
						}
					}
					}
				}
				else if( nlen == 3 ){
					const double C_old[9] = { pVal_C[0], pVal_C[1], pVal_C[2], 
						                      pVal_C[3], pVal_C[4], pVal_C[5], 
											  pVal_C[6], pVal_C[7], pVal_C[8], };
					for(unsigned int ilen=0;ilen<nlen;ilen++){
					for(unsigned int jlen=0;jlen<nlen;jlen++){
						pVal_C[nlen*ilen+jlen] = 0;
						pVal_C[nlen*ilen+jlen] 
							+= pVal_Mi[nlen*ilen+ilen]*C_old[nlen*ilen+jlen]*pVal_Mi[nlen*jlen+jlen];
/*						for(unsigned int klen=0;klen<nlen;klen++){
						for(unsigned int llen=0;llen<nlen;llen++){
							pVal_C[nlen*ilen+jlen] 
								+= pVal_Mi[nlen*ilen+klen]*C_old[nlen*klen+llen]*pVal_Mi[nlen*llen+jlen];
						}
						}*/
					}
					}
				}
			}
		}
	}
	return true;
}


bool CLinearSystem_Eigen::SetVector_fromField(int iv1,
	unsigned int id_field, const Fem::Field::CFieldWorld& world, Fem::Field::FIELD_DERIVATION_TYPE fdt )
{
	if( !world.IsIdField(id_field) ){ assert(0); return false; }
//	const CField& field = world.GetField(id_field);
	const unsigned int nseg = this->m_aSegField.size();
	if( nseg == 0 )	return true;

    assert( nseg == 1 );

    MatVec::CVector_Blk& vec = m_ls.GetVector(iv1,0);
/*	std::vector< CVector_Blk* >* paVec1 = 0;
	if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) paVec1 = &m_ls.m_TmpVectorArray[iv1];
	else if( iv1 == -1 ) paVec1 = &this->m_ls.m_Residual;
	else if( iv1 == -2 ) paVec1 = &this->m_ls.m_Update;
	else assert(0);*/


	// コーナー節点について値を更新
	{
		const CField& field = world.GetField(id_field);
		unsigned int id_na = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		const CNodeAry& na = world.GetNA(id_na);
		const CNodeAry::CNodeSeg& ns = field.GetNodeSeg(CORNER,true,world,fdt);
		const unsigned int nblk = na.Size();
		const unsigned int nlen = ns.Length();
		double* val = new double [nlen];
		for(unsigned int iblk=0;iblk<nblk;iblk++){
			ns.GetValue(iblk,val);
			for(unsigned int ilen=0;ilen<nlen;ilen++){
				vec.SetValue(iblk,ilen,val[ilen]);
			}
		}
		delete[] val;
	}
	return true;
}
