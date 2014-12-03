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

#if !defined(LINEAR_SYSTEM_RIGID_H)
#define LINEAR_SYSTEM_RIGID_H

#include <vector>
#include <cassert>
#include <math.h>

#include "delfem/vector3d.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/ls/linearsystem_interface_solver.h"
#include "delfem/matvec/matdiafrac_blkcrs.h"


////////////////////////////////////////////////////////////////

namespace Rigid{
    class CRigidBody3D;
    class CConstraint;
}

namespace Ls
{

class CLinearSystem_RigidBody
{
public:
    virtual unsigned int GetSizeRigidBody() const = 0;
    
    ////////////////
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
        const Com::CVector3D& vres, double d ) = 0;
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
        const double* eres, double d ) = 0;
    virtual void SubResidual(const unsigned int ind, bool is_rb, 
        const double* res) = 0;
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans) = 0;
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                                  unsigned int indl, bool is_rb_l, unsigned int offsetl,
                                  const Com::CVector3D& vec, double d, bool is_colum) = 0;
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
                           const double* emat, double d) = 0;
    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta, 
        bool is_first) const = 0;
};
/*
class CLinearSystem_RigidBody_Full : public CLinearSystem_RigidBody
{
public:
    CLinearSystem_RigidBody_Full(const std::vector<Rigid::CRigidBody3D>& aRB,
        const std::vector<Rigid::CConstraint*>& aConst);
    ~CLinearSystem_RigidBody_Full(){
        delete[] mat;
        delete[] residual;
        delete[] update;
    }
    unsigned int GetSizeRigidBody() const { return nRB; }

    ////////////////////////////////////////////////////////////////
    
    void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, const Com::CVector3D& vres, double d ){
        unsigned int idof = (is_rb) ? aIndDof[ind]+offset : aIndDof[ind+nRB]+offset;
        assert( idof+3 <= this->ndof );
        residual[idof+0] += vres.x*d;
        residual[idof+1] += vres.y*d;
        residual[idof+2] += vres.z*d;
    }
    void SubResidual(const unsigned int ind, bool is_rb, const double* res){
        const unsigned int ind0 = (is_rb) ? ind : ind + nRB;
        const unsigned int idof0 = aIndDof[ind0  ];
        const unsigned int n     = aIndDof[ind0+1] - idof0;
        for(unsigned int i=0;i<n;i++){
            residual[i+idof0] -= res[i];
        }
    }
    void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                   unsigned int indl, bool is_rb_l, unsigned int offsetl,
                   const Com::CMatrix3& m, double d, bool isnt_trans)
    {
        unsigned int idofr = (is_rb_r) ? aIndDof[indr]+offsetr : aIndDof[indr+nRB]+offsetr;
        unsigned int idofl = (is_rb_l) ? aIndDof[indl]+offsetl : aIndDof[indl+nRB]+offsetl;
        assert( idofr+3 <= this->ndof );
        assert( idofl+3 <= this->ndof );
        if( isnt_trans ){
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                mat[(idofr+i)*ndof+idofl+j] += m.mat[i*3+j]*d;
            }
            }
        }
        else{
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                mat[(idofr+i)*ndof+idofl+j] += m.mat[j*3+i]*d;
            }
            }
        }
    }
    void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                         unsigned int indl, bool is_rb_l, unsigned int offsetl,
                         const Com::CVector3D& vec, double d, bool is_column)
    {
        const unsigned int idofr = (is_rb_r) ? aIndDof[indr]+offsetr : aIndDof[indr+nRB]+offsetr;
        const unsigned int idofl = (is_rb_l) ? aIndDof[indl]+offsetl : aIndDof[indl+nRB]+offsetl;
        if( is_column ){
            mat[(idofr+0)*ndof+idofl] += vec.x*d;
            mat[(idofr+1)*ndof+idofl] += vec.y*d;
            mat[(idofr+2)*ndof+idofl] += vec.z*d;
        }
        else{
            mat[idofr*ndof+idofl+0] += vec.x*d;
            mat[idofr*ndof+idofl+1] += vec.y*d;
            mat[idofr*ndof+idofl+2] += vec.z*d;
        }
    }

    ////////////////////////////////////////////////////////////////
    void InitializeMarge(){
        unsigned int i;
        for(i=0;i<ndof*ndof;i++){ mat[i] = 0; }
        for(i=0;i<ndof     ;i++){ update[i] = 0; }
        for(i=0;i<ndof     ;i++){ residual[i] = 0; }
    }
    virtual double FinalizeMarge(){
        double norm_res = 0;
        for(unsigned int i=0;i<ndof;i++){
            norm_res += residual[i]*residual[i];
        }
        return sqrt( norm_res );
    }
    ////////////////////////////////////////////////////////////////
    void Solve();
    
    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta, 
        bool is_first) const;
private:
    const double* GetUpdate(unsigned int ind_obj, bool is_rb ) const{
        if( is_rb ){ return &update[ aIndDof[ind_obj    ] ]; }
        else{        return &update[ aIndDof[ind_obj+nRB] ]; }
    }
private:
    unsigned int nRB;
    unsigned int nConst;
    std::vector<unsigned int> aIndDof;
    unsigned int ndof;
    double* mat;
    double* update;
    double* residual;
};

class CMat_CrsBlkFlex{
public:
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



class CLinearSystem_RigidBody_CRS : public CLinearSystem_RigidBody, public Sol::CLinearSystem_SolInterface
{
    friend class CPreconditioner_RigidBody_CRS;
public:
    CLinearSystem_RigidBody_CRS(){
        ndof = 0;
        ColInd = 0; RowPtr = 0; ValPtr = 0; ValList = 0;
        Update = 0; Residual = 0;
    }
    CLinearSystem_RigidBody_CRS(const std::vector<Rigid::CRigidBody3D>& aRB,
        const std::vector<Rigid::CConstraint*>& aConst)
    {
        ndof = 0;
        ColInd = 0; RowPtr = 0; ValPtr = 0; ValList = 0;
        Update = 0; Residual = 0;
        this->SetRigidSystem(aRB,aConst);
    }
    virtual ~CLinearSystem_RigidBody_CRS(){
        this->Clear();
    }
    void Clear(){   // データを全て削除
        for(unsigned int ivec=0;ivec<aTmpVec.size();ivec++){
            delete[] aTmpVec[ivec];
        }
        aTmpVec.clear();
        if( ColInd   != 0 ){ delete[] ColInd;   ColInd   = 0; }
        if( RowPtr   != 0 ){ delete[] RowPtr;   RowPtr   = 0; }
        if( ValPtr   != 0 ){ delete[] ValPtr;   ValPtr   = 0; }
        if( ValList  != 0 ){ delete[] ValList;  ValList  = 0; }
        if( Update   != 0 ){ delete[] Update;   Update   = 0; }
        if( Residual != 0 ){ delete[] Residual; Residual = 0; }
    }
    void SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB,
        const std::vector<Rigid::CConstraint*>& aConst);
    unsigned int GetSizeRigidBody() const { return nRB; }

    ////////////////////////////////////////////////////////////////
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
        const Com::CVector3D& vres, double d )
    {
        unsigned int idof = (is_rb) ? aIndDof[ind]+offset : aIndDof[ind+nRB]+offset;
        assert( idof+3 <= this->ndof );
        Residual[idof+0] += vres.x*d;
        Residual[idof+1] += vres.y*d;
        Residual[idof+2] += vres.z*d;
    }
    virtual void SubResidual(const unsigned int ind, bool is_rb, const double* res)
    {
        const unsigned int ind0 = (is_rb) ? ind : ind + nRB;
        const unsigned int idof0 = aIndDof[ind0  ];
        const unsigned int n     = aIndDof[ind0+1] - idof0;
        for(unsigned int i=0;i<n;i++){
            Residual[i+idof0] -= res[i];
        }
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans )
    {
        const unsigned int ind0 = (is_rb_r) ? indr : indr + nRB;
        const unsigned int ind1 = (is_rb_l) ? indl : indl + nRB;
        unsigned int icrs0 = ncrs;
        for(unsigned int icrs=ColInd[ind0];icrs<ColInd[ind0+1];icrs++){
            if( RowPtr[icrs] == ind1 ){
                icrs0 = icrs;
                break;
            }
        }
        assert( icrs0 < ncrs );
        const unsigned int ind_val = ValPtr[icrs0];
        const unsigned int lenrow = aIndDof[ind1+1] - aIndDof[ind1];
        if( isnt_trans ){
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                ValList[ind_val + (i+offsetr)*lenrow + (j+offsetl)] += m.mat[i*3+j]*d;
            }
            }
        }
        else{
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                ValList[ind_val + (i+offsetr)*lenrow + (j+offsetl)] += m.mat[j*3+i]*d;
            }
            }
        }
    }
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                         unsigned int indl, bool is_rb_l, unsigned int offsetl,
                         const Com::CVector3D& vec, double d, bool is_column)
    {
        const unsigned int ind0 = (is_rb_r) ? indr : indr + nRB;
        const unsigned int ind1 = (is_rb_l) ? indl : indl + nRB;
        unsigned int icrs0 = ncrs;
        for(unsigned int icrs=ColInd[ind0];icrs<ColInd[ind0+1];icrs++){
            if( RowPtr[icrs] == ind1 ){
                icrs0 = icrs;
                break;
            }
        }
        assert( icrs0 < ncrs );
        const unsigned int ind_val = ValPtr[icrs0];
        const unsigned int lenrow = aIndDof[ind1+1] - aIndDof[ind1];
        if( is_column ){
            ValList[ind_val + (0+offsetr)*lenrow + (offsetl)] += vec.x*d;
            ValList[ind_val + (1+offsetr)*lenrow + (offsetl)] += vec.y*d;
            ValList[ind_val + (2+offsetr)*lenrow + (offsetl)] += vec.z*d;
        }
        else{
            ValList[ind_val + (offsetr)*lenrow + (0+offsetl)] += vec.x*d;
            ValList[ind_val + (offsetr)*lenrow + (1+offsetl)] += vec.y*d;
            ValList[ind_val + (offsetr)*lenrow + (2+offsetl)] += vec.z*d;
        }
    }

    virtual void InitializeMarge(){
        if( ValPtr != 0 ){
            const unsigned int ntotval = ValPtr[ncrs];
            for(unsigned int i=0;i<ntotval;i++){ ValList[i] = 0; }
        }
        if( ndof != 0 ){
            int i;
            for(i=0;i<ndof     ;i++){ Update[i] = 0; }
            for(i=0;i<ndof     ;i++){ Residual[i] = 0; }
        }
    }
    virtual double FinalizeMarge(){
        double norm_res = 0;
        for(unsigned int i=0;i<ndof;i++){
            norm_res += Residual[i]*Residual[i];
        }
        return sqrt( norm_res );
    }
    ////////////////////////////////////////////////////////////////
    virtual unsigned int GetTmpVectorArySize() const{ 
        return aTmpVec.size();
    }
    virtual bool ReSizeTmpVecSolver(unsigned int isize){
        if( aTmpVec.size() > isize ){
            const unsigned int isize0 = aTmpVec.size();
            for(unsigned int i=isize;i<isize0;i++){
                delete aTmpVec[i];
            }
            aTmpVec.resize(isize);
        }
        else if( aTmpVec.size() < isize ){
            const unsigned int isize0 = aTmpVec.size();
            aTmpVec.resize(isize,0);
            for(unsigned int i=isize0;i<isize;i++){
                aTmpVec[i] = new double [ndof];
            }
        }
        return true; 
    }
    virtual double DOT(int iv1,int iv2); // {v2}*{v1}
    virtual bool COPY(int iv1,int iv2);  // {v2} := {v1}
    virtual bool SCAL(double d,int iv);  // {v1} := alpha * {v1}
    virtual bool AXPY(double d,int iv1,int iv2); // {v2} := alpha*{v1} + {v2}
    virtual bool MATVEC(double a,int iv1,double b,int iv2);  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}

    double* GetVector(unsigned int iv){
    	double* p_vec = 0;
	    if( iv >= 0 && iv < (int)this->GetTmpVectorArySize() ) p_vec = aTmpVec[iv];
	    else if( iv == -1 ) p_vec = this->Residual;
	    else if( iv == -2 ) p_vec = this->Update;
	    else assert(0);
        return p_vec;
    }

    virtual void Solve(const Ls::CPreconditioner_RigidBody_CRS& prec);
    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta, 
        bool is_first) const;
    virtual const double* GetUpdate(unsigned int ind_obj, bool is_rb ) const{ 
        if( is_rb ){ return &Update[ aIndDof[ind_obj    ] ]; }
        else{        return &Update[ aIndDof[ind_obj+nRB] ]; }
    }
    ////////////////////////////////////////////////////////////////
private:
    unsigned int nRB;
    unsigned int nConst;
    // nRB+nConstが行や列のブロックの数
    std::vector<unsigned int> aIndDof;
    unsigned int ndof;
    unsigned int ncrs;
    unsigned int* ColInd;   // ブロックCRSのColumnIndex
    unsigned int* RowPtr;   // ブロックCRSのRowPointer
    unsigned int* ValPtr;   // CRS番号からブロック行列の値の先頭番地を返す
    std::vector<double*> aTmpVec;
    double* ValList;
    double* Update;
    double* Residual;
};


class CPreconditioner_RigidBody_CRS
{
public:
    CPreconditioner_RigidBody_CRS(){}
    virtual ~CPreconditioner_RigidBody_CRS(){
        delete[] ColInd;
        delete[] RowPtr;
        delete[] ValPtr;
        delete[] ValList;
    }
    void SetLinearSystem(const CLinearSystem_RigidBody_CRS& ls);
    void SetValue(const CLinearSystem_RigidBody_CRS& ls);
//
//    void Factorization();
    bool Solve(double* vec) const;
private:
    unsigned int nRB;
    unsigned int nConst;
    // nRB+nConstが行や列のブロックの数
    std::vector<unsigned int> aIndDof;
    unsigned int ndof;

    unsigned int ncrs;
    unsigned int* ColInd;   // ブロックCRSのColumnIndex
    unsigned int* RowPtr;   // ブロックCRSのRowPointer
    unsigned int* DiaInd;   // 対角ブロックへのCRS番号
    unsigned int* ValPtr;   // CRS番号からブロック行列の値の先頭番地を返す
    double* ValList;
};


//! 前処理行列クラスの抽象クラス
class CLinearSystemPreconditioner_RigidBody_CRS: public Sol::CLinearSystemPreconditioner_SolInterface
{
public:
    CLinearSystemPreconditioner_RigidBody_CRS(CLinearSystem_RigidBody_CRS& ls, 
        CPreconditioner_RigidBody_CRS& prec) : ls(ls), prec(prec){}
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
        return prec.Solve( ls.GetVector(iv) );
    }
private:
    CPreconditioner_RigidBody_CRS& prec;
    CLinearSystem_RigidBody_CRS& ls;
};

*/

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



class CLinearSystem_RigidBody_CRS2 : public CLinearSystem_RigidBody, public LsSol::ILinearSystem_Sol
{
    friend class CPreconditioner_RigidBody_CRS;
public:
    CLinearSystem_RigidBody_CRS2(){}
    CLinearSystem_RigidBody_CRS2(const std::vector<Rigid::CRigidBody3D>& aRB,
        const std::vector<Rigid::CConstraint*>& aConst){
        this->SetRigidSystem(aRB,aConst);
    }
    virtual ~CLinearSystem_RigidBody_CRS2(){
        this->Clear();
    }
    void Clear(){   // データを全て削除
    }
    void SetRigidSystem(const std::vector<Rigid::CRigidBody3D>& aRB,
        const std::vector<Rigid::CConstraint*>& aConst);

    unsigned int GetSizeRigidBody() const { return nRB; }

    ////////////////////////////////////////////////////////////////
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, 
        const Com::CVector3D& vres, double d )
    {
        unsigned int iblk = (is_rb) ? ind : ind+nRB;
        assert( offset+3 <= m_residual.Len(iblk) );
        m_residual.AddValue(iblk,offset+0,vres.x*d);
        m_residual.AddValue(iblk,offset+1,vres.y*d);
        m_residual.AddValue(iblk,offset+2,vres.z*d);
    }
    virtual void AddResidual(unsigned int ind, bool is_rb, unsigned int offset, unsigned int size,
        const double* eres, double d )
    {
        unsigned int iblk = (is_rb) ? ind : ind+nRB;
        assert( offset+size <= m_residual.Len(iblk) );
        for(unsigned int i=0;i<size;i++){
            m_residual.AddValue(iblk,offset+i,eres[i]*d);
        }
    }
    virtual void SubResidual(const unsigned int ind, bool is_rb, const double* res)
    {
        unsigned int iblk = (is_rb) ? ind : ind+nRB;
        unsigned int len = m_residual.Len(iblk);
        for(unsigned int i=0;i<len;i++){
            m_residual.AddValue(iblk,i,res[i]*-1);
        }
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl,
                           const Com::CMatrix3& m, double d, bool isnt_trans )
    {
        const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
        const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
        const unsigned int lencol = m_mat.LenBlkCol(iblk0);
        const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
        assert( lencol*lenrow <= 36 );
        double tmp[36];
        for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
        if( isnt_trans ){
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[i*3+j]*d;
            }
            }
        }
        else{
            for(unsigned int i=0;i<3;i++){
            for(unsigned int j=0;j<3;j++){
                tmp[(i+offsetr)*lenrow + (j+offsetl)] += m.mat[j*3+i]*d;
            }
            }
        }
        m_mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
    }
    virtual void AddMatrix(unsigned int indr, bool is_rb_r, unsigned int offsetr, unsigned int sizer,
                           unsigned int indl, bool is_rb_l, unsigned int offsetl, unsigned int sizel,
                           const double* emat, double d){
        const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
        const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
        const unsigned int lencol = m_mat.LenBlkCol(iblk0);
        const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
        assert( lencol*lenrow <= 36 );
        double tmp[36];
        for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
        for(unsigned int i=0;i<sizer;i++){
        for(unsigned int j=0;j<sizel;j++){
            tmp[(i+offsetr)*lenrow + (j+offsetl)] += emat[i*sizel+j]*d;
        }
        }
        m_mat.Mearge(1,&iblk0, 1,&iblk1, lencol*lenrow, tmp );
    }
    virtual void AddMatrix_Vector(unsigned int indr, bool is_rb_r, unsigned int offsetr,
                         unsigned int indl, bool is_rb_l, unsigned int offsetl,
                         const Com::CVector3D& vec, double d, bool is_column)
    {
        const unsigned int iblk0 = (is_rb_r) ? indr : indr + nRB;
        const unsigned int iblk1 = (is_rb_l) ? indl : indl + nRB;
        const unsigned int lencol = m_mat.LenBlkCol(iblk0);
        const unsigned int lenrow = m_mat.LenBlkRow(iblk1);
        assert( lencol*lenrow <= 36 );
        double tmp[36];
        for(unsigned int i=0;i<lencol*lenrow;i++){ tmp[i] = 0; }
        if( is_column ){
            tmp[ (0+offsetr)*lenrow + (offsetl)] += vec.x*d;
            tmp[ (1+offsetr)*lenrow + (offsetl)] += vec.y*d;
            tmp[ (2+offsetr)*lenrow + (offsetl)] += vec.z*d;
        }
        else{
            tmp[ (offsetr)*lenrow + (0+offsetl)] += vec.x*d;
            tmp[ (offsetr)*lenrow + (1+offsetl)] += vec.y*d;
            tmp[ (offsetr)*lenrow + (2+offsetl)] += vec.z*d;
        }
        m_mat.Mearge(1,&iblk0, 1, &iblk1, lencol*lenrow, tmp );
    }
    virtual void InitializeMarge(){
        m_mat.SetZero();
        m_update.SetVectorZero();
        m_residual.SetVectorZero();
    }
    virtual double FinalizeMarge(){
        double norm_res = m_residual.GetSquaredVectorNorm();
        return sqrt( norm_res );
    }
    ////////////////////////////////////////////////////////////////
    virtual unsigned int GetTmpVectorArySize() const{ return m_aTmpVec.size(); }
    virtual bool ReSizeTmpVecSolver(unsigned int isize){
        if( m_aTmpVec.size() > isize ){ m_aTmpVec.resize(isize); }
        else if( m_aTmpVec.size() < isize ){
            const unsigned int isize0 = m_aTmpVec.size();
            m_aTmpVec.resize(isize);
            for(unsigned int i=isize0;i<isize;i++){
                m_aTmpVec[i].Initialize(nRB+nConst,m_aBlkSize);
            }
        }
        return true; 
    }
    virtual double DOT(int iv1,int iv2); // {v2}*{v1}
    virtual bool COPY(int iv1,int iv2);  // {v2} := {v1}
    virtual bool SCAL(double d,int iv);  // {v1} := alpha * {v1}
    virtual bool AXPY(double d,int iv1,int iv2); // {v2} := alpha*{v1} + {v2}
    virtual bool MATVEC(double a,int iv1,double b,int iv2);  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}

    MatVec::CVector_Blk& GetVector(unsigned int iv){
	    if( iv >= 0 && iv < (int)this->GetTmpVectorArySize() ) return m_aTmpVec[iv];
	    else if( iv == -1 ) return this->m_residual;
        assert( iv == -2 );
        return this->m_update;
    }
    const MatVec::CMatDia_BlkCrs& GetMatrix() const { return m_mat; }
        
    virtual bool UpdateValueOfRigidSystem(
        std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
        double dt, double newmark_gamma, double newmark_beta, 
        bool is_first) const;
    ////////////////////////////////////////////////////////////////
private:

private:
    unsigned int nRB;
    unsigned int nConst;
    // nRB+nConstが行や列のブロックの数

    std::vector<unsigned int> m_aBlkSize;
    MatVec::CMatDia_BlkCrs m_mat;
    MatVec::CVector_Blk m_residual;
    MatVec::CVector_Blk m_update;
    std::vector< MatVec::CVector_Blk > m_aTmpVec;
};


class CPreconditioner_RigidBody_CRS2
{
public:
    CPreconditioner_RigidBody_CRS2(){}
    virtual ~CPreconditioner_RigidBody_CRS2(){}
    void SetLinearSystem(const CLinearSystem_RigidBody_CRS2& ls, int ilev = 0){
        m_mat.MakePattern_Initialize( ls.GetMatrix() );
        m_mat.AddFracPtn(-1);
        m_mat.MakePatternFinalize();
    }
    void SetValue(const CLinearSystem_RigidBody_CRS2& ls){
        m_mat.SetValue( ls.GetMatrix() );
    }
    bool Solve( MatVec::CVector_Blk& vec ) const{
        return m_mat.Solve( vec );
    }
private:
    MatVec::CMatDiaFrac_BlkCrs m_mat;
};

//! 前処理行列クラスの抽象クラス
class CLinearSystemPreconditioner_RigidBody_CRS2: public LsSol::ILinearSystemPreconditioner_Sol
{
public:
    CLinearSystemPreconditioner_RigidBody_CRS2(CLinearSystem_RigidBody_CRS2& ls, 
        CPreconditioner_RigidBody_CRS2& prec) : ls(ls), prec(prec){}
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
        assert( iv == -1 || iv == -2 || (iv>=0 && iv<ls.GetTmpVectorArySize()) );
        MatVec::CVector_Blk& vec = ls.GetVector(iv);
        return prec.Solve( vec );
        return true;
    }
private:
    CPreconditioner_RigidBody_CRS2& prec;
    CLinearSystem_RigidBody_CRS2& ls;
};


}

#endif