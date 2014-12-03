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
@brief çÑëÃÉNÉâÉX(Com::CRigidBody3D)ÇÃé¿ëï
@author Nobuyuki Umetani
*/

#if !defined(RIGID_BODY_H)
#define RIGID_BODY_H

#include <vector>
#include <cassert>
#include <math.h>

#include "delfem/vector3d.h"

namespace Ls{
    class CLinearSystem_RigidBody;
}

namespace Rigid
{

class CRigidBody3D
{
public:
    CRigidBody3D(){
		// èâä˙ê›íË
        this->ini_pos_cg.SetVector(0,0,0);
		this->mineatia[0] = 1.0; 
		this->mineatia[1] = 1.0;
		this->mineatia[2] = 1.0;
        this->mass = 1.0;

		// ílÇÃèâä˙âª
        this->disp_cg.SetVector(0,0,0);
        this->velo_cg.SetVector(0,0,0);
        this->acc_cg.SetVector(0,0,0);
        crv[0] = 0;   crv[1] = 0;  crv[2] = 0;
        Omega.SetVector(0,0,0);
        dOmega.SetVector(0,0,0);
    }
    unsigned int GetDOF() const { return 6; }
    void GetInvRotMatrix44(double* rot) const;
    void GetRotMatrix33(double* rot) const;
    Com::CMatrix3 GetRotMatrix() const;
//    void Draw() const;
    void Clear(){
        crv[0] = 0; crv[1] = 0; crv[2] = 0;
        Omega.SetVector(0,0,0);
        dOmega.SetVector(0,0,0);
        velo_cg.SetVector(0,0,0);
        disp_cg.SetVector(0,0,0);
        acc_cg.SetVector(0,0,0);
    }
    Com::CVector3D GetPositionFromInital(const Com::CVector3D& vec) const
    {
        const Com::CMatrix3& rot = this->GetRotMatrix();
        return rot.MatVec(vec-ini_pos_cg)+ini_pos_cg+disp_cg;
    }
    void AddRotation( const double* rot );
    void UpdateSolution(
        const double* upd,
        const double dt, const double newmark_gamma, const double newmark_beta, 
        bool is_first_iter);
    void AddLinearSystem(
        Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const Com::CVector3D& gravity, 
        bool is_first);
public:
    double crv[3];
    Com::CVector3D Omega;
    Com::CVector3D dOmega;
    Com::CVector3D velo_cg;
    Com::CVector3D disp_cg;
    Com::CVector3D acc_cg;

    // èâä˙ê›íË
    Com::CVector3D ini_pos_cg;
    double mass;
    double mineatia[3];
};

////////////////////////////////////////////////////////////////

class CConstraint
{
public:
    virtual unsigned int GetDOF() const = 0;
    virtual void Clear() = 0;
    virtual void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta ) = 0;
    virtual void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB, 
        bool is_initial ) const = 0;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const = 0;
    const std::vector<unsigned int>& GetAry_IndexRB(){ return aIndRB; }
public:
    std::vector<unsigned int> aIndRB;
};



class CFix_Spherical : public CConstraint
{
public:
    CFix_Spherical(unsigned int irb){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
        ini_pos_fix.SetVector(0,0,0);
        aIndRB.push_back(irb);
    }
    virtual unsigned int GetDOF() const { return 3; }
    virtual void Clear(){ 
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
    }
    void SetIniPosFix(double x, double y, double z){ ini_pos_fix.SetVector(x,y,z); }
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB, bool is_initial = false ) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[3];
    Com::CVector3D ini_pos_fix;
private:
};

class CFix_Hinge : public CConstraint
{
public:
    CFix_Hinge(unsigned int irb){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
        axis.SetVector(1,0,0);
        loc_coord[0].SetVector(0,0,0);
        loc_coord[1].SetVector(0,0,0);
        ini_pos_fix.SetVector(0,0,0);
        aIndRB.push_back(irb);
    }
    virtual unsigned int GetDOF() const { return 5; }
    virtual void Clear(){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;  lambda[3] = 0;  lambda[4] = 0;
    }
    void SetIniPosFix(double x, double y, double z){
        this->ini_pos_fix.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB, bool is_initial = false ) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[5];
    Com::CVector3D ini_pos_fix;
    Com::CVector3D loc_coord[2];
private:
    Com::CVector3D axis;
};

class CFix_HingeRange : public CConstraint
{
public:
    CFix_HingeRange(unsigned int irb){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
        axis.SetVector(1,0,0);
        loc_coord[0].SetVector(0,0,0);
        loc_coord[1].SetVector(0,0,0);
        ini_pos_fix.SetVector(0,0,0);
        aIndRB.push_back(irb);
        min_t = -180;
        max_t =  180;
    }
    virtual unsigned int GetDOF() const { return 6; }
    virtual void Clear(){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;  lambda[3] = 0;  lambda[4] = 0; lambda[5] = 0;
    }
    void SetIniPosFix(double x, double y, double z){
        this->ini_pos_fix.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void SetRange(double min_t, double max_t){
        this->min_t = min_t;
        this->max_t = max_t;
    }
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB, 
        bool is_initial = false ) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[6];
    Com::CVector3D ini_pos_fix;
    Com::CVector3D loc_coord[2];
    double min_t;
    double max_t;
private:
    Com::CVector3D axis;
};

class CJoint_Spherical : public CConstraint
{
public:
    CJoint_Spherical(unsigned int irb0, unsigned int irb1){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
        ini_pos_joint.SetVector(0,0,0);
        aIndRB.push_back(irb0);
        aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 3; }
    virtual void Clear(){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
        this->ini_pos_joint.SetVector(x,y,z);
    }
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB, 
        bool is_initial = false) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[3];
    Com::CVector3D ini_pos_joint;
private:
};

class CJoint_Hinge : public CConstraint
{
public:
    CJoint_Hinge(unsigned int irb0, unsigned int irb1){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
        ini_pos_joint.SetVector(0,0,0);
        aIndRB.push_back(irb0);
        aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 5; }
    virtual void Clear(){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
        this->ini_pos_joint.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB,
        bool is_initial = false ) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[5];
    Com::CVector3D ini_pos_joint;
    Com::CVector3D loc_coord[2];
private:
    Com::CVector3D axis;
};

class CJoint_HingeRange : public CConstraint
{
public:
    CJoint_HingeRange(unsigned int irb0, unsigned int irb1){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
        ini_pos_joint.SetVector(0,0,0);
        aIndRB.push_back(irb0);
        aIndRB.push_back(irb1);
    }
    virtual unsigned int GetDOF() const { return 6; }
    virtual void Clear(){
        lambda[0] = 0;  lambda[1] = 0;  lambda[2] = 0; lambda[3] = 0; lambda[4] = 0; lambda[5] = 0;
    }
    void SetIniPosJoint(double x, double y, double z){
        this->ini_pos_joint.SetVector(x,y,z);
    }
    void SetAxis(double ax, double ay, double az);
    void SetRange(double min_t, double max_t){
        this->min_t = min_t;
        this->max_t = max_t;
    }
    void UpdateSolution(const double* upd,
        double dt, double newmark_gamma, double newmark_beta);
    void AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
        const double dt, const double newmark_gamma, const double newmark_beta,
        const std::vector<CRigidBody3D>& aRB,
        bool is_initial = false ) const;
//    virtual void Draw(const std::vector<CRigidBody3D>& aRB) const;
public:
    double lambda[6];
    Com::CVector3D ini_pos_joint;
    Com::CVector3D loc_coord[2];
    double min_t;
    double max_t;
private:
    Com::CVector3D axis;
};



}

#endif