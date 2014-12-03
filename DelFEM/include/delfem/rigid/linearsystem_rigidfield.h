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

/*! @file
@brief 剛体クラス(Com::CRigidBody3D)の実装
@author Nobuyuki Umetani
*/

#if !defined(LINEAR_SYSTEM_RIGID_FIELD_H)
#define LINEAR_SYSTEM_RIGID_FIELD_H

#include <vector>
#include <cassert>
#include <math.h>

#include "delfem/vector3d.h"
#include "delfem/ls/preconditioner.h"
#include "delfem/rigid/linearsystem_rigid.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/linearsystem_interface_eqnsys.h"


////////////////////////////////////////////////////////////////

namespace Rigid{
    class CRigidBody3D;
    class CConstraint;
}

namespace Ls
{
/*
class CPreconditioner_RigidField;
class CLinearSystem_RigidField 
: public CLinearSystem_RigidBody, 
  public Sol::CLinearSystem_SolInterface
{
public:
    CLinearSystem_RigidField(){}
    virtual ~CLinearSystem_RigidField(){ this->Clear(); }
    void Clear(){ 
        ls_rigid.Clear();
        ls_field.Clear();
    }
    void SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB, 
        const std::vector<Rigid::CConstraint*>& aConst)
    {
        ls_rigid.SetRigidSystem(aRB,aConst);
    }
    void AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world)
    {
        ls_field.AddPattern_Field(id_field,world);
    }
    void AddPattern_Field(const unsigned int id_field, const unsigned int id_field0, 
        const Fem::Field::CFieldWorld& world )
    {
        ls_field.AddPattern_Field(id_field,id_field0,world);
    }
    unsigned int GetSizeRigidBody() const { return ls_rigid.GetSizeRigidBody(); }
    void SetFixedBoundaryCondition_Field(unsigned int id_disp_fix0, Fem::Field::CFieldWorld& world){
        ls_field.SetFixedBoundaryCondition_Field(id_disp_fix0,world);
    }
    void UpdateValueOfField_NewmarkBeta(double newmark_gamma, double newmark_beta, double dt,
        unsigned int id_disp  , Fem::Field::CFieldWorld& world, bool is_first){
        ls_field.UpdateValueOfField_NewmarkBeta(newmark_gamma, newmark_beta, dt,
            id_disp, world, is_first );
    }

    ////////////////////////////////////////////////////////////////
    void InitializeMarge(){ 
        ls_rigid.InitializeMarge(); 
        ls_field.InitializeMarge();
    }
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
        const Com::CVector3D& vres, double d ){
        ls_rigid.AddResidual(ind,is_rb,offset,vres,d);
    }
    virtual void SubResidual(const unsigned int ind, bool is_rb, const double* res){
        ls_rigid.SubResidual(ind,is_rb,res);
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans ){
        ls_rigid.AddMatrix(indr,is_rb_r,offsetr,   indl,is_rb_l,offsetl,   m,d,isnt_trans);
    }
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_column){
        ls_rigid.AddMatrix_Vector(indr,is_rb_r,offsetr,   indl,is_rb_l,offsetl,   vec,d,is_column);
    }
    double FinalizeMarge(){
        const double res_r = ls_rigid.FinalizeMarge();
        const double res_f = ls_field.FinalizeMarge();
        return sqrt( res_r*res_r + res_f*res_f );
    }
    ////////////////////////////////////////////////////////////////
    virtual unsigned int GetTmpVectorArySize() const{ 
        assert( ls_rigid.GetTmpVectorArySize() == ls_field.GetTmpVectorArySize() );
        return ls_field.GetTmpVectorArySize();
    }
    virtual bool ReSizeTmpVecSolver(unsigned int isize){ 
        ls_rigid.ReSizeTmpVecSolver(isize);
        ls_field.ReSizeTmpVecSolver(isize);
        return true; 
    }
    virtual double DOT(int iv0,int iv1){
        double d0 = ls_field.DOT(iv0,iv1);
        d0       += ls_rigid.DOT(iv0,iv1);
        return d0;
    }
    virtual bool COPY(int iv0,int iv1){
        ls_field.COPY(iv0,iv1);
        ls_rigid.COPY(iv0,iv1);
        return true;
    }
    virtual bool SCAL(double d,int iv){ 
        ls_field.SCAL(d,iv);
        ls_rigid.SCAL(d,iv);
        return true;
    }
    virtual bool AXPY(double d,int iv0,int iv1){ 
        ls_field.AXPY(d,iv0,iv1);
        ls_rigid.AXPY(d,iv0,iv1);
        return true; 
    }
    virtual bool MATVEC(double a,int iv0,double b,int iv1){ 
        ls_field.MATVEC(a,iv0,b,iv1);
        ls_rigid.MATVEC(a,iv0,b,iv1);
        return true; 
    }

    ////////////////////////////////////////////////////////////////

    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta,         bool is_first) const
    {
        return ls_rigid.UpdateValueOfRigidSystem(aRB,aConst,  dt,newmark_gamma,newmark_beta, is_first);
    }

    void Solve(Ls::CPreconditioner_RigidField& prec);

    ////////////////////////////////////////////////////////////////

    CLinearSystem_RigidBody_CRS2& GetLinearSystemRigid(){ return ls_rigid; }
    const CLinearSystem_RigidBody_CRS2& GetLinearSystemRigid() const { return ls_rigid; }

    Fem::Ls::CLinearSystem_Field& GetLinearSystemField(){ return ls_field; }
    const Fem::Ls::CLinearSystem_Field& GetLinearSystemField() const { return ls_field; }
private:
    CLinearSystem_RigidBody_CRS2 ls_rigid;
    Fem::Ls::CLinearSystem_Field ls_field;
};


class CPreconditioner_RigidField
{
public:
    CPreconditioner_RigidField(){
        prec_field.SetFillInLevel(0);
    }
    virtual ~CPreconditioner_RigidField(){
        this->prec_field.Clear();
    }
    ////////////////
    void SetLinearSystem(const CLinearSystem_RigidField& ls)
    {
        std::cout << "prec set ls rigid" << std::endl;
        prec_rigid.SetLinearSystem( ls.GetLinearSystemRigid() );
        std::cout << "prec set ls rigid end" << std::endl;
        prec_field.SetLinearSystem( ls.GetLinearSystemField().m_ls );
    }
    void SetValue(const CLinearSystem_RigidField& ls)
    {
        std::cout << "prec set val rigid" << std::endl;
        prec_rigid.SetValue(ls.GetLinearSystemRigid());
        std::cout << "prec set val rigid end" << std::endl;
        prec_field.SetValue(ls.GetLinearSystemField().m_ls );
    }
    void SolvePrecond(CLinearSystem_RigidField& ls, int iv){
        prec_field.SolvePrecond( ls.GetLinearSystemField().m_ls, iv );
        prec_rigid.Solve( ls.GetLinearSystemRigid().GetVector(iv) );
    }
    const CPreconditioner_RigidBody_CRS2& GetPreconditionerRigid() const { return prec_rigid; }
    CPreconditioner_RigidBody_CRS2& GetPreconditionerRigid(){ return prec_rigid; }
    const LsSol::CPreconditioner_ILU&  GetPreconditionerField() const { return prec_field; }
private:
    CPreconditioner_RigidBody_CRS2 prec_rigid;
    LsSol::CPreconditioner_ILU prec_field;
};

//! 前処理行列クラスの抽象クラス
class CLinearSystemPreconditioner_RigidField: public Sol::CLinearSystemPreconditioner_SolInterface
{
public:
    CLinearSystemPreconditioner_RigidField(CLinearSystem_RigidField& ls, 
        CPreconditioner_RigidField& prec) : ls(ls), prec(prec){}
public:
	//! ソルバに必要な作業ベクトルの数を得る
    virtual unsigned int GetTmpVectorArySize() const{ return ls.GetTmpVectorArySize(); }
	//! ソルバに必要な作業ベクトルの数を設定
    virtual bool ReSizeTmpVecSolver(unsigned int size_new){ return ls.ReSizeTmpVecSolver(size_new); }

    virtual double DOT(int iv1, int iv2){ return ls.DOT(iv1,iv2); }
    virtual bool COPY(int iv1, int iv2){ return ls.COPY(iv1,iv2); }
    virtual bool SCAL(double alpha, int iv1){ return ls.SCAL(alpha,iv1); }
    virtual bool AXPY(double alpha, int iv1, int iv2){ return ls.AXPY(alpha,iv1,iv2); }
    virtual bool MATVEC(double alpha, int iv1, double beta, int iv2){ return ls.MATVEC(alpha,iv1,beta,iv2); }

    virtual bool SolvePrecond(int iv){
        prec.SolvePrecond(ls,iv);
        return true;
    }
private:
    CPreconditioner_RigidField& prec;
    CLinearSystem_RigidField& ls;
};
*/

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CLinearSystem_RigidField2 
: public Fem::Eqn::ILinearSystem_Eqn, 
  public CLinearSystem_RigidBody, 
  public LsSol::ILinearSystem_Sol
{
public:
    CLinearSystem_RigidField2(){ ilss_rigid = -1; }
    virtual ~CLinearSystem_RigidField2(){ this->Clear(); }
    void Clear(){ m_ls.Clear(); this->m_aSegRF.clear(); }

    ////////////////////////////////////////////////////////////////
    // 場へのインターフェース
    int FindIndexArray_Seg( unsigned int id_field, Fem::Field::ELSEG_TYPE type, const Fem::Field::CFieldWorld& world );
    int GetIndexSegRigid(){ return ilss_rigid; }

    bool AddPattern_Field(const unsigned int id_field, const Fem::Field::CFieldWorld& world);
    bool AddPattern_Field(const unsigned int id_field, const unsigned int id_field0, 
        const Fem::Field::CFieldWorld& world );
    bool SetFixedBoundaryCondition_Field(unsigned int id_disp_fix0, const Fem::Field::CFieldWorld& world);
    bool UpdateValueOfField_NewmarkBeta(double newmark_gamma, double newmark_beta, double dt,
        unsigned int id_disp  , Fem::Field::CFieldWorld& world, bool is_first);

    virtual MatVec::CVector_Blk& GetResidual( 
        unsigned int id_field, 
        Fem::Field::ELSEG_TYPE type, 
        const Fem::Field::CFieldWorld& world ){
        int ils = this->FindIndexArray_Seg(id_field,type,world);
        assert( ils >= 0 && ils < m_ls.GetNLinSysSeg() );
        return m_ls.GetVector(-1,ils);
    }
    virtual MatVec::CMatDia_BlkCrs& GetMatrix( 
        unsigned int id_field, 
        Fem::Field::ELSEG_TYPE type, 
        const Fem::Field::CFieldWorld& world ){
        int ils = this->FindIndexArray_Seg(id_field,type,world);
        assert( ils >= 0 && ils < m_ls.GetNLinSysSeg() );
        return m_ls.GetMatrix(ils);
    }
    virtual MatVec::CMat_BlkCrs& GetMatrix( 
        unsigned int id_field1, Fem::Field::ELSEG_TYPE type1, 
        unsigned int id_field2, Fem::Field::ELSEG_TYPE type2, 
        const Fem::Field::CFieldWorld& world )
    {
        int ils1 = this->FindIndexArray_Seg(id_field1,type1,world);
        int ils2 = this->FindIndexArray_Seg(id_field2,type2,world);
        assert( ils1 >= 0 && ils1 < m_ls.GetNLinSysSeg() );
        assert( ils2 >= 0 && ils2 < m_ls.GetNLinSysSeg() );
        return m_ls.GetMatrix(ils1,ils2);
    }

    ////////////////////////////////////////////////////////////////
    // 剛体へのインターフェース
    void SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB, 
        const std::vector<Rigid::CConstraint*>& aConst);
    unsigned int GetSizeRigidBody() const {
        assert( ilss_rigid < m_aSegRF.size() );
        assert( m_aSegRF[ilss_rigid].is_rigid );
        return m_aSegRF[ilss_rigid].nRB;
    }
    unsigned int GetSizeConstraint() const {
        assert( ilss_rigid < m_aSegRF.size() );
        assert( m_aSegRF[ilss_rigid].is_rigid );
        return m_aSegRF[ilss_rigid].nConst;
    }
    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta,         bool is_first) const;

    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
        const Com::CVector3D& vres, double d );
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
        const double* eres, double d );
    virtual void SubResidual(const unsigned int ind, bool is_rb, const double* res);
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans );
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_column);
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
                           const double* emat, double d);

    ////////////////////////////////
    // 剛体弾性体連成インターフェース
    bool AddPattern_RigidField(
        const unsigned int id_field, const unsigned int id_field0, const Fem::Field::CFieldWorld& world, 
        unsigned int irb, std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst);

    ////////////////////////////////

    virtual MatVec::CMatDia_BlkCrs& GetMatrix( unsigned int ils ){
        assert( ils < m_ls.GetNLinSysSeg() );
        return m_ls.GetMatrix(ils);
    }
    virtual MatVec::CMat_BlkCrs& GetMatrix( unsigned int ils, unsigned int jls ){
        assert( ils < m_ls.GetNLinSysSeg() );
        assert( jls < m_ls.GetNLinSysSeg() );
        return m_ls.GetMatrix(ils,jls);
    }
    virtual MatVec::CVector_Blk& GetResidual( unsigned int ils ){
        assert( ils < m_ls.GetNLinSysSeg() );
        return m_ls.GetVector(-1,ils);
    }
    unsigned int GetNLinSysSeg() const { 
        assert( m_aSegRF.size() == m_ls.GetNLinSysSeg() );
        return m_aSegRF.size(); 
    }
    void InitializeMarge(){ m_ls.InitializeMarge(); }
    double FinalizeMarge(){ return m_ls.FinalizeMarge(); }

    ////////////////////////////////////////////////////////////////
    // Solverへのインターフェース
    virtual unsigned int GetTmpVectorArySize() const{ return m_ls.GetTmpVectorArySize(); }
    virtual bool ReSizeTmpVecSolver(unsigned int isize){ return m_ls.ReSizeTmpVecSolver(isize); }
    virtual double DOT(int iv0,int iv1){ return m_ls.DOT(iv0,iv1); }
    virtual bool COPY(int iv0,int iv1){ return m_ls.COPY(iv0,iv1); }
    virtual bool SCAL(double d,int iv){ return m_ls.SCAL(d,iv); }
    virtual bool AXPY(double d,int iv0,int iv1){ return m_ls.AXPY(d,iv0,iv1); }
    virtual bool MATVEC(double a,int iv0,double b,int iv1){ return m_ls.MATVEC(a,iv0,b,iv1); }

    ////////////////////////////////////////////////////////////////

private:
	class CLinSysSegRF{
    public:
        CLinSysSegRF(unsigned int id_f, Fem::Field::ELSEG_TYPE iconf){
            this->is_rigid = false;
            id_field = id_f;
            node_config = iconf;
        }
        CLinSysSegRF(unsigned int nRB, unsigned int nConst){
            this->is_rigid = true;
            this->nRB = nRB;
            this->nConst = nConst;
        }
	public:
        bool is_rigid;
		unsigned int id_field;	// parent_fieldでなければならない
        Fem::Field::ELSEG_TYPE node_config;
        ////////////////
        // nRB+nConstが行や列のブロックの数
        unsigned int nRB;
        unsigned int nConst;
	};
    unsigned int ilss_rigid;
public:
	std::vector< CLinSysSegRF > m_aSegRF;
    LsSol::CLinearSystem m_ls;
};


class CPreconditioner_RigidField2
{
public:
    CPreconditioner_RigidField2(){
        prec.SetFillInLevel(0);
    }
    virtual ~CPreconditioner_RigidField2(){ this->Clear(); }
    ////////////////
	void SetFillInLevel(int lev, int ilss0 = -1){ 
        prec.SetFillInLevel(lev,ilss0);
    }
    void Clear(){ prec.Clear(); }
    void SetLinearSystem(const CLinearSystem_RigidField2& ls)
    {
        prec.SetLinearSystem( ls.m_ls );
    }
    void SetValue(const CLinearSystem_RigidField2& ls)
    {
        prec.SetValue( ls.m_ls );
    }
    void SolvePrecond(CLinearSystem_RigidField2& ls, int iv){
        prec.SolvePrecond( ls.m_ls,iv );
    }
private:
    LsSol::CPreconditioner_ILU prec;
};


//! 前処理行列クラスの抽象クラス
class CLinearSystemPreconditioner_RigidField2: public LsSol::ILinearSystemPreconditioner_Sol
{
public:
    CLinearSystemPreconditioner_RigidField2(CLinearSystem_RigidField2& ls, 
        CPreconditioner_RigidField2& prec) : ls(ls), prec(prec){}
public:
	//! ソルバに必要な作業ベクトルの数を得る
    virtual unsigned int GetTmpVectorArySize() const{ return ls.GetTmpVectorArySize(); }
	//! ソルバに必要な作業ベクトルの数を設定
    virtual bool ReSizeTmpVecSolver(unsigned int size_new){ return ls.ReSizeTmpVecSolver(size_new); }

    virtual double DOT(int iv1, int iv2){ return ls.DOT(iv1,iv2); }
    virtual bool COPY(int iv1, int iv2){ return ls.COPY(iv1,iv2); }
    virtual bool SCAL(double alpha, int iv1){ return ls.SCAL(alpha,iv1); }
    virtual bool AXPY(double alpha, int iv1, int iv2){ return ls.AXPY(alpha,iv1,iv2); }
    virtual bool MATVEC(double alpha, int iv1, double beta, int iv2){ return ls.MATVEC(alpha,iv1,beta,iv2); }

    virtual bool SolvePrecond(int iv){
        prec.SolvePrecond(ls,iv);
        return true;
    }
private:
    CPreconditioner_RigidField2& prec;
    CLinearSystem_RigidField2& ls;
};






} // end namespace LS

#endif