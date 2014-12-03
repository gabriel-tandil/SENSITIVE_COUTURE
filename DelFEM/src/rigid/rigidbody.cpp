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
#pragma warning ( disable : 4786 ) 
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif  /* _WIN32 */
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

//#include <GL/glut.h>



#include <iostream>
#include "delfem/rigid/rigidbody.h"
#include "delfem/rigid/linearsystem_rigid.h"

/*
bool GetVertical2Vector(const Com::CVector3D& vec0, Com::CVector3D& vec1, Com::CVector3D& vec2)
{
    Com::CVector3D vec_a = vec0;
    const double len = vec_a.Length();
    if( len < 1.0e-20 ){ 
        assert(0);
        return false; 
    }
    vec_a *= 1.0/len;
    Com::CVector3D vec_tmp0(0,1,0);
    vec1 = Com::Cross(vec_tmp0,vec_a);
    if( vec1.Length() < 1.0e-20 ){
        // vec_a ?y??????
        Com::CVector3D vec_tmp0(1,0,0);
        vec1 = Com::Cross(vec_tmp0,vec_a);  // z????
        vec2 = Com::Cross(vec_a,vec1);  // x????
    }
    else{
        vec1.Normalize();
        vec2 = Com::Cross(vec_a,vec1);
    }
    return true;
}
*/

void MakeRotMatrix33_CartesianRotationVector(double* rot, const Com::CVector3D& vec)
{
    const double theta = vec.Length();
	if( theta < 1.0e-30 ){
		rot[0] = 1;		rot[1] = 0;		rot[2] = 0;
		rot[3] = 0;		rot[4] = 1;		rot[5] = 0;
		rot[6] = 0;		rot[7] = 0;		rot[8] = 1;
		return;
	}
    const double inv_theta = 1.0/theta;
    const double uvec[3] = { vec.x*inv_theta, vec.y*inv_theta, vec.z*inv_theta };

    const double st = sin(theta);
    const double ct = 1-cos(theta);

    const double w[9] = { 0,-uvec[2],+uvec[1], +uvec[2],0,-uvec[0],  -uvec[1],+uvec[0], 0 };

    for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            rot[i*3+j]  = st*w[i*3+j] + ct*(w[i*3+0]*w[0*3+j]+w[i*3+1]*w[1*3+j]+w[i*3+2]*w[2*3+j]);
        }
        rot[i*3+i] += 1;
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/*
void Rigid::CRigidBody3D::Draw() const
{
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        ::glPushMatrix();
        ::glTranslated( 
            ini_pos_cg.x + disp_cg.x,
            ini_pos_cg.y + disp_cg.y,
            ini_pos_cg.z + disp_cg.z );
        {
            double rot0[16];
		    GetInvRotMatrix44(rot0);
            glMultMatrixd(rot0);
        }
        ::glColor3d(0,0,0);
	    ::glutSolidCube(0.2);
    //    glutSolidTeapot(1.0);
    //    glutWireTeapot(0.2);
    //    glutWireOctahederon();
    //    glutWireDodecahedron();
        ::glLineWidth(1);
        ::glBegin(GL_LINES);
        ::glColor3d(1,0,0);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0.4,0,0);
        ::glColor3d(0,1,0);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0,0.4,0);
        ::glColor3d(0,0,1);
        ::glVertex3d(0,0,0);
        ::glVertex3d(0,0,0.4);
        ::glEnd();
	    ::glPopMatrix();
    }
}
*/

void Rigid::CRigidBody3D::GetInvRotMatrix44(double* rot) const 
{
    const double c0 = 0.125*( 16.0 - crv[0]*crv[0] - crv[1]*crv[1] - crv[2]*crv[2] );
    const double tmp = 1.0/( (4.0-c0)*(4.0-c0) );
    rot[0*4+0] = tmp*( (c0*c0+8*c0-16) + 2*crv[0]*crv[0] );
    rot[1*4+0] = tmp*(                   2*crv[0]*crv[1] - 2*c0*crv[2] );
    rot[2*4+0] = tmp*(                   2*crv[0]*crv[2] + 2*c0*crv[1] );
    rot[3*4+0] = 0;
    rot[0*4+1] = tmp*(                   2*crv[1]*crv[0] + 2*c0*crv[2] );
    rot[1*4+1] = tmp*( (c0*c0+8*c0-16) + 2*crv[1]*crv[1] );
    rot[2*4+1] = tmp*(                   2*crv[1]*crv[2] - 2*c0*crv[0] );
    rot[3*4+1] = 0;
    rot[0*4+2] = tmp*(                   2*crv[2]*crv[0] - 2*c0*crv[1] );
    rot[1*4+2] = tmp*(                   2*crv[2]*crv[1] + 2*c0*crv[0] );
    rot[2*4+2] = tmp*( (c0*c0+8*c0-16) + 2*crv[2]*crv[2] );
    rot[3*4+2] = 0;
    rot[0*4+3] = 0;
    rot[1*4+3] = 0;
    rot[2*4+3] = 0;
    rot[3*4+3] = 1;
}

// CRV‚Ì‰ñ“]s—ñ‚ð“¾‚é
void Rigid::CRigidBody3D::GetRotMatrix33(double* rot) const
{
    const double c0 = 0.125*( 16.0 - crv[0]*crv[0] - crv[1]*crv[1] - crv[2]*crv[2] );
    const double tmp = 1.0/( (4.0-c0)*(4.0-c0) );
    rot[0*3+0] = tmp*( (c0*c0+8*c0-16) + 2*crv[0]*crv[0] );
    rot[0*3+1] = tmp*(                   2*crv[0]*crv[1] - 2*c0*crv[2] );
    rot[0*3+2] = tmp*(                   2*crv[0]*crv[2] + 2*c0*crv[1] );
    rot[1*3+0] = tmp*(                   2*crv[1]*crv[0] + 2*c0*crv[2] );
    rot[1*3+1] = tmp*( (c0*c0+8*c0-16) + 2*crv[1]*crv[1] );
    rot[1*3+2] = tmp*(                   2*crv[1]*crv[2] - 2*c0*crv[0] );
    rot[2*3+0] = tmp*(                   2*crv[2]*crv[0] - 2*c0*crv[1] );
    rot[2*3+1] = tmp*(                   2*crv[2]*crv[1] + 2*c0*crv[0] );
    rot[2*3+2] = tmp*( (c0*c0+8*c0-16) + 2*crv[2]*crv[2] );
}

Com::CMatrix3 Rigid::CRigidBody3D::GetRotMatrix() const
{
    Com::CMatrix3 m;
    const double c0 = 0.125*( 16.0 - crv[0]*crv[0] - crv[1]*crv[1] - crv[2]*crv[2] );
    const double tmp = 1.0/( (4.0-c0)*(4.0-c0) );
    assert( fabs(4.0-c0) > 1.0e-30 );
    m.mat[0*3+0] = tmp*( (c0*c0+8*c0-16) + 2*crv[0]*crv[0] );
    m.mat[0*3+1] = tmp*(                   2*crv[0]*crv[1] - 2*c0*crv[2] );
    m.mat[0*3+2] = tmp*(                   2*crv[0]*crv[2] + 2*c0*crv[1] );
    m.mat[1*3+0] = tmp*(                   2*crv[1]*crv[0] + 2*c0*crv[2] );
    m.mat[1*3+1] = tmp*( (c0*c0+8*c0-16) + 2*crv[1]*crv[1] );
    m.mat[1*3+2] = tmp*(                   2*crv[1]*crv[2] - 2*c0*crv[0] );
    m.mat[2*3+0] = tmp*(                   2*crv[2]*crv[0] - 2*c0*crv[1] );
    m.mat[2*3+1] = tmp*(                   2*crv[2]*crv[1] + 2*c0*crv[0] );
    m.mat[2*3+2] = tmp*( (c0*c0+8*c0-16) + 2*crv[2]*crv[2] );
    return m;
}

void Rigid::CRigidBody3D::AddRotation( const double* rot_a )
{
    double rot0[9];
    this->GetRotMatrix33(rot0);
    double rot1[9];
    for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
        rot1[i*3+j] = rot0[i*3+0]*rot_a[0*3+j] + rot0[i*3+1]*rot_a[1*3+j] + rot0[i*3+2]*rot_a[2*3+j];
    }
    }

    const double smat[16] = {
        1+rot1[0*3+0]+rot1[1*3+1]+rot1[2*3+2],  
        rot1[2*3+1]-rot1[1*3+2],
        rot1[0*3+2]-rot1[2*3+0],
        rot1[1*3+0]-rot1[0*3+1],
        rot1[2*3+1]-rot1[1*3+2],
        1+rot1[0*3+0]-rot1[1*3+1]-rot1[2*3+2],
        rot1[0*3+1]+rot1[1*3+0],
        rot1[0*3+2]+rot1[2*3+0],
        rot1[0*3+2]-rot1[2*3+0],
        rot1[1*3+0]+rot1[0*3+1],
        1-rot1[0*3+0]+rot1[1*3+1]-rot1[2*3+2],
        rot1[1*3+2]+rot1[2*3+1],
        rot1[1*3+0]-rot1[0*3+1],
        rot1[0*3+2]+rot1[2*3+0],
        rot1[1*3+2]+rot1[2*3+1],
        1-rot1[0*3+0]-rot1[1*3+1]+rot1[2*3+2],
    };

    unsigned int imax;
    imax = ( smat[0   *4+   0] > smat[1*4+1] ) ? 0    : 1;
    imax = ( smat[imax*4+imax] > smat[2*4+2] ) ? imax : 2;
    imax = ( smat[imax*4+imax] > smat[3*4+3] ) ? imax : 3;

    double eparam2[4];
    eparam2[imax] = 0.5*sqrt(smat[imax*4+imax]);
    for(unsigned int k=0;k<4;k++){
        if( k==imax ) continue;
        eparam2[k] = smat[imax*4+k]*0.25/eparam2[imax];
    }
    /*
	double eparam[4];
	{
		const double tr = rot1[0]+rot1[4]+rot1[8];
		double vect_rot1[3] = {
			0.5*(rot1[7]-rot1[5]), 
			0.5*(rot1[2]-rot1[6]), 
			0.5*(rot1[3]-rot1[1]) 
        };
        if( 1+tr <= 0 ){ 
            std::cout << 1+tr << std::endl;
            eparam[0] = 0; 
        }
        else{            eparam[0] = 0.5*sqrt(1+tr); }
        for(unsigned int i=0;i<3;i++){
            double sign;
            if( fabs(vect_rot1[i]) < 1.0e-30 ){ sign = 1.0; }
            else{ sign = vect_rot1[i]/fabs(vect_rot1[i]); }
             if( 1+2*rot1[i*4]-tr <= 0 ){ eparam[i+1] = 0; }
            else{ eparam[i+1] = 0.5*sign*sqrt(1+2*rot1[i*4]-tr); }
        }
	}
    crv[0] = 4*eparam[1]/(1+eparam[0]);
    crv[1] = 4*eparam[2]/(1+eparam[0]);
    crv[2] = 4*eparam[3]/(1+eparam[0]);
    */
    crv[0] = 4*eparam2[1]/(1+eparam2[0]);
    crv[1] = 4*eparam2[2]/(1+eparam2[0]);
    crv[2] = 4*eparam2[3]/(1+eparam2[0]);
}

void Rigid::CRigidBody3D::UpdateSolution(const double* upd, 
                                double dt, double newmark_gamma, double newmark_beta, 
                                bool is_first_iter )
{
    if( is_first_iter ){
        const unsigned int imode = 1;
        if( imode == 0 ){
            disp_cg += velo_cg*dt + acc_cg*0.5*dt*dt;
            velo_cg += acc_cg*dt;
            ////////////////
            {
                Com::CVector3D rot_update = dt*Omega + 0.5*dt*dt*dOmega; 
		        double derot[9];
		        MakeRotMatrix33_CartesianRotationVector( derot, rot_update );
		        AddRotation( derot );
            }
            Omega += dOmega*dt;
        }
        else{
            disp_cg += velo_cg*dt + acc_cg*(0.5-newmark_beta)*dt*dt;
            velo_cg += acc_cg*(1-newmark_gamma)*dt;
            acc_cg.SetZero();
            ////////////////
            {
                Com::CVector3D rot_update = dt*Omega + (0.5-newmark_beta)*dt*dt*dOmega; 
		        double derot[9];
		        MakeRotMatrix33_CartesianRotationVector( derot, rot_update );
		        AddRotation( derot );
            }
            Omega += dOmega*(1-newmark_gamma)*dt;
            dOmega.SetZero();
        }
        return;
    }
    const Com::CVector3D acc_disp_update( upd[0], upd[1], upd[2] );
    const Com::CVector3D acc_rot_update(  upd[3], upd[4], upd[5] );
    const double dtmp1 = dt*dt*newmark_beta;
    const double dtmp2 = dt*newmark_gamma;

	// •ÏˆÊ‚ÌƒAƒbƒvƒf[ƒg
    disp_cg += dtmp1*acc_disp_update;
    // ‘¬“x‚ÌƒAƒbƒvƒf[ƒg
    velo_cg += dtmp2*acc_disp_update;
    // ‰Á‘¬“x‚ÌƒAƒbƒvƒf[ƒg
    acc_cg += acc_disp_update;

    ////////////////////////////////////////////////
	{	// ‰ñ“]‚ÌƒAƒbƒvƒf[ƒg
        Com::CVector3D rot_update = dtmp1*acc_rot_update;
		double derot[9];
		MakeRotMatrix33_CartesianRotationVector( derot, rot_update );
		AddRotation( derot );
	}
    const Com::CVector3D oo = Omega;
    const Com::CMatrix3 woo(oo);
    {   // Šp‘¬“x‚ÌƒAƒbƒvƒf[ƒg
        Omega += dtmp2*acc_rot_update;
        Omega += dtmp1*woo.MatVec(acc_rot_update);
    }
    {   // Šp‰Á‘¬“x‚ÌƒAƒbƒvƒf[ƒg
        const Com::CMatrix3 wdoo(dOmega);
        const Com::CMatrix3 woowoo = woo.MatMat(woo);
        dOmega += acc_rot_update;
        dOmega += dtmp1*( 0.5*wdoo.MatVec(acc_rot_update)-1.0/6.0*woowoo.MatVec(acc_rot_update) );
    }
}
                               

void Rigid::CRigidBody3D::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int irb,
    const double dt, const double newmark_gamma, const double newmark_beta,
    const Com::CVector3D& gravity, 
    bool is_initial)
{
	// •Ài‚ÌŽc·
    ls.AddResidual( irb,true,0,  gravity - acc_cg, mass );
	{   // Šµ«—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·
        double res[6] = {0,0,0, 0,0,0};
		res[3] = mineatia[0]*dOmega.x;
        res[4] = mineatia[1]*dOmega.y;
        res[5] = mineatia[2]*dOmega.z;
		const double wOmg[9] = { 0,-Omega.z,+Omega.y, +Omega.z,0,-Omega.x,  -Omega.y,+Omega.x, 0 };
        const double jo[3] = { mineatia[0]*Omega.x, mineatia[1]*Omega.y, mineatia[2]*Omega.z };
		for(unsigned int i=0;i<3;i++){
		for(unsigned int j=0;j<3;j++){
			res[i+3] += wOmg[i*3+j]*jo[j];
		}
    	}
        ls.SubResidual(irb,true, res);
    }
	{   // „«s—ñ‚ðì‚é
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity(mass);
        ls.AddMatrix(irb,true,0,  irb,true,0,  mtmp,1, true);
        mtmp.mat[0] = mineatia[0]; mtmp.mat[1] = 0;           mtmp.mat[2] = 0;
        mtmp.mat[3] = 0;           mtmp.mat[4] = mineatia[1]; mtmp.mat[5] = 0;
        mtmp.mat[6] = 0;           mtmp.mat[7] = 0;           mtmp.mat[8] = mineatia[2];
        ls.AddMatrix(irb,true,3,  irb,true,3,  mtmp,1, true);
	}
}

////////////////////////////////////////////////////////////////

void Rigid::CFix_Spherical::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                     const double dt, const double newmark_gamma, const double newmark_beta, 
                     const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
    const unsigned int irb = aIndRB[0];
    assert( irb < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb = aRB[irb];

    ////////////////////////////////////////////////
    const Com::CMatrix3& mrot = rb.GetRotMatrix();
    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
    const Com::CVector3D& Xdistfix = rb.ini_pos_cg - ini_pos_fix;

    // S‘©—ÍŽc·
	// •Ài‚ÌŽc·
    ls.AddResidual(irb,true,0, vlambda, 1);
	{   // S‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·
        Com::CVector3D tmp_v = mrot.MatVecTrans(vlambda);
        Com::CMatrix3 wX(Xdistfix);
        const Com::CVector3D& vec_q = wX.MatVec(tmp_v);
        ls.AddResidual(irb,true,3, vec_q, -1);
	}
	{   // S‘©ðŒ‚Ì•Ï•ª
        const Com::CVector3D& rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
        ls.AddResidual(icst,false,0,  rb.ini_pos_cg+rb.disp_cg-rot_pos_cg, 1 );
	}
    ////////////////////////////////
    Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
    Com::CMatrix3 wXwRtL;
	{
        const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix);
        wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity(-dtmp1);
        ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   1,     true );
        ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   1,     true );
        ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL, dtmp1, true );
        ls.AddMatrix(irb, true, 3, icst,false,0, RotwX, -dtmp1, false);
        ls.AddMatrix(icst,false,0, irb, true, 3, RotwX, -dtmp1, true );
	}
}

void Rigid::CFix_Spherical::UpdateSolution(const double* upd,
    double dt, double newmark_gamma, double newmark_beta)
{
    const double tmp1 = dt*dt*newmark_beta;
//    const double* upd = ls.GetUpdate(icst,false);
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
}

/*
void Rigid::CFix_Spherical::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 2;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        ::glColor3d(1,1,0);
        ::glBegin(GL_POINTS);
        ::glVertex3d(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
        ::glEnd();
    }
    else if( imode == 2 ){
        const unsigned int irb = aIndRB[0];
        const CRigidBody3D& rb = aRB[irb];
        const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
        const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.disp_cg;

        ::glPushMatrix();
        ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();

        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
        ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
        ::glEnd();
    }
}
*/

////////////////////////////////////////////////////////////////

void Rigid::CFix_Hinge::SetAxis(double ax, double ay, double az){
    axis = Com::CVector3D(ax,ay,az);
    if( axis.Length() < 1.0e-20 ){ assert(0); return; }
    axis.Normalize();
    GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

void Rigid::CFix_Hinge::UpdateSolution(const double* upd,
                    double dt, double newmark_gamma, double newmark_beta)
{
//    const double* upd = ls.GetUpdate(icst,false);
    const double tmp1 = dt*dt*newmark_beta;
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
    lambda[3] += tmp1*upd[3];
    lambda[4] += tmp1*upd[4];
}

void Rigid::CFix_Hinge::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
    const double dt, const double newmark_gamma, const double newmark_beta,
    const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
    const unsigned int irb = aIndRB[0];
    assert( irb < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb = aRB[irb];
    
    const Com::CMatrix3& mrot = rb.GetRotMatrix();
    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
    const Com::CVector3D& Xdistfix = rb.ini_pos_cg - ini_pos_fix;

    ////////////////////////////////////////////////
    // S‘©—ÍŽc·
	// •Ài‚ÌŽc·
    ls.AddResidual(irb,true,0, vlambda, 1);
	{   // ˆÊ’uS‘©—Í‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& tmp_vec = mrot.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix);
        ls.AddResidual(irb,true,3,  wX.MatVec(tmp_vec), -1 );
	}
	{   // S‘©ŠðŒ‚Ì•Ï•ª
        double res[5] = {0,0,0,0,0};
        // •ÏˆÊS‘©ðŒ‚Ì•Ï•ª
        Com::CVector3D rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
        Com::CVector3D res_disp = - rb.ini_pos_cg - rb.disp_cg + rot_pos_cg;
        res[0] = res_disp.x; 
        res[1] = res_disp.y; 
        res[2] = res_disp.z;
        // ‰ñ“]S‘©ðŒ‚Ì•Ï•ª
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        res[3] = Com::Dot(Rta,loc_coord[0]);
        res[4] = Com::Dot(Rta,loc_coord[1]);
        ls.SubResidual(icst,false,res);
	}
    {   // ‰ñ“]‚‘¬‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        const Com::CMatrix3 wTmp(lambda[3]*loc_coord[0] + lambda[4]*loc_coord[1]);
        ls.AddResidual(irb,true,3, wTmp.MatVec(Rta),-1 );
    }
    ////////////////////////////////
    Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
    Com::CMatrix3 wXwRtL;
	{
        const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix);
        wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    Com::CVector3D wm0Rta, wm1Rta;
    Com::CMatrix3 lwLCwRta;
    {
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        const Com::CMatrix3 wm0(loc_coord[0]);
        const Com::CMatrix3 wm1(loc_coord[1]);
        wm0Rta = wm0.MatVec(Rta);
        wm1Rta = wm1.MatVec(Rta);
        Com::CMatrix3 wTmp(lambda[3]*loc_coord[0]+lambda[4]*loc_coord[1]);
        lwLCwRta = wTmp.MatMat( Com::CMatrix3(Rta) );
    }
	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity();
        ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   -dtmp1, true);
        ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   -dtmp1, true);

        ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL,  dtmp1, true);
        ls.AddMatrix(irb, true, 3, irb, true, 3, lwLCwRta,dtmp1, true);

        ls.AddMatrix(irb, true, 3, icst,false,0, RotwX,  -dtmp1, false);
        ls.AddMatrix(icst,false,0, irb, true, 3, RotwX,  -dtmp1, true );

        ls.AddMatrix_Vector(irb, true, 3, icst,false,3, wm0Rta,dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3, irb, true, 3, wm0Rta,dtmp1, false);
        ls.AddMatrix_Vector(irb, true, 3, icst,false,4, wm1Rta,dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4, irb, true, 3, wm1Rta,dtmp1, false);
	}
}
/*
void Rigid::CFix_Hinge::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        const unsigned int irb = aIndRB[0];
        const CRigidBody3D& rb = aRB[irb];
        const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
        const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.disp_cg;

        ::glPushMatrix();
        ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();

        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
        ::glEnd();

        const Com::CVector3D& lcb0 = this->loc_coord[0];
        const Com::CVector3D& lcb1 = this->loc_coord[1];
        unsigned int ndiv = 16;
        const double dtheta = 2*3.1416/ndiv;
        const double radius = 0.4;
        ::glColor3d(0,1,1);
        ::glBegin(GL_TRIANGLE_FAN);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
            const Com::CVector3D v0 = vec_j + sin(idiv*dtheta  )*lcb0*radius + cos(idiv*dtheta  )*lcb1*radius;
            ::glVertex3d(v0.x,  v0.y,  v0.z);
        }
        ::glEnd();
    }
}
*/

////////////////////////////////////////////////////////////////

void Rigid::CFix_HingeRange::SetAxis(double ax, double ay, double az){
    axis = Com::CVector3D(ax,ay,az);
    if( axis.Length() < 1.0e-20 ){ assert(0); return; }
    axis.Normalize();
    GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

void Rigid::CFix_HingeRange::UpdateSolution(const double* upd,
                    double dt, double newmark_gamma, double newmark_beta)
{
//    const double* upd = ls.GetUpdate(icst,false);
    const double tmp1 = dt*dt*newmark_beta;
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
    lambda[3] += tmp1*upd[3];
    lambda[4] += tmp1*upd[4];
//    lambda[5] += tmp1*upd[5];
    const double tmp2 = dt*newmark_gamma;
    lambda[5] += tmp2*upd[5];
}

void Rigid::CFix_HingeRange::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
    const double dt, double newmark_gamma, double newmark_beta,
    const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
    const unsigned int irb = aIndRB[0];
    assert( irb < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb = aRB[irb];
    
    const Com::CMatrix3& mrot = rb.GetRotMatrix();
    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
    const Com::CVector3D& Xdistfix = rb.ini_pos_cg - ini_pos_fix;

    ////////////////////////////////////////////////
    // S‘©—ÍŽc·
	// •Ài‚ÌŽc·
    ls.AddResidual(irb,true,0, vlambda, 1);
	{   // ˆÊ’uS‘©—Í‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& tmp_vec = mrot.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix);
        ls.AddResidual(irb,true,3,  wX.MatVec(tmp_vec), -1 );
	}
	{   // S‘©ŠðŒ‚Ì•Ï•ª
        double res[6] = {0,0,0,0,0,0};
        // •ÏˆÊS‘©ðŒ‚Ì•Ï•ª
        Com::CVector3D rot_pos_cg = mrot.MatVec(Xdistfix) + ini_pos_fix;
        Com::CVector3D res_disp = - rb.ini_pos_cg - rb.disp_cg + rot_pos_cg;
        res[0] = res_disp.x; 
        res[1] = res_disp.y; 
        res[2] = res_disp.z;
        // ‰ñ“]S‘©ðŒ‚Ì•Ï•ª
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        res[3] = Com::Dot(Rta,loc_coord[0]);
        res[4] = Com::Dot(Rta,loc_coord[1]);
        res[5] = 0;
        ls.SubResidual(icst,false,res);
	}
    {   // ‰ñ“]‚‘¬‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        const Com::CMatrix3 wTmp(lambda[3]*loc_coord[0] + lambda[4]*loc_coord[1]);
        ls.AddResidual(irb,true,3, wTmp.MatVec(Rta),-1 );
    }
    ////////////////////////////////
    Com::CMatrix3 RotwX = mrot.MatMat( Com::CMatrix3(Xdistfix) );
    Com::CMatrix3 wXwRtL;
	{
        const Com::CVector3D RtL = mrot.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix);
        wXwRtL = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    Com::CVector3D wm0Rta, wm1Rta;
    Com::CMatrix3 lwLCwRta;
    {
        const Com::CVector3D& Rta = mrot.MatVecTrans(axis);
        const Com::CMatrix3 wm0(loc_coord[0]);
        const Com::CMatrix3 wm1(loc_coord[1]);
        wm0Rta = wm0.MatVec(Rta);
        wm1Rta = wm1.MatVec(Rta);
        Com::CMatrix3 wTmp(lambda[3]*loc_coord[0]+lambda[4]*loc_coord[1]);
        lwLCwRta = wTmp.MatMat( Com::CMatrix3(Rta) );
    }
    {   // Šp“x‚É‚Â‚¢‚Ä‚ÌKKTðŒ‚ð“ü‚ê‚é
        double m0Rm0 = Com::Dot(loc_coord[0],mrot.MatVec(loc_coord[0]));
        double m1Rm0 = Com::Dot(loc_coord[1],mrot.MatVec(loc_coord[0]));
        const double PI = 3.14159265358979323846;
        double u_min[2] = { cos(PI*min_t/180.0), sin(PI*min_t/180.0) };
        double u_max[2] = { cos(PI*max_t/180.0), sin(PI*max_t/180.0) };
        double v_min[2] = { -u_min[1],  u_min[0] };
        double v_max[2] = {  u_max[1], -u_max[0] };
        bool flg = false;
        double v[2];
        double f_value;
        const double delta_f = 1.0e-5;
        {
            double dmin = v_min[0]*m0Rm0 + v_min[1]*m1Rm0;
            double dmax = v_max[0]*m0Rm0 + v_max[1]*m1Rm0;
//            std::cout << dmin << " " << dmax << " " << lambda[5] << std::endl;
            if( dmin < -delta_f || dmax < -delta_f ) flg = true;
            const double emin = u_min[0]*m0Rm0 + u_min[1]*m1Rm0;
            const double emax = u_max[0]*m0Rm0 + u_max[1]*m1Rm0;
            if( emin > emax ){ 
                v[0] = v_min[0]; v[1] = v_min[1];
                f_value = dmin;
            }
            else{                              
                v[0] = v_max[0]; v[1] = v_max[1];
                f_value = dmax;
            }
        }
/*        if( flg==true || lambda[5]<-1.0e-5 ){
            Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0Rtm = wm0.MatVec( Rtm );
            double c = 10000;
            const double df = Com::Dot(wm0Rtm,rb.Omega);
            std::cout << "vm0 : " << df << std::endl;
            std::cout << "Omega : " << rb.Omega.x << " " << rb.Omega.y << " " << rb.Omega.z << std::endl;
            Com::CVector3D res_t(c*wm0Rtm.x*df, c*wm0Rtm.y*df, c*wm0Rtm.z*df);
            ls.AddResidual(irb, true, 3,   res_t,-1 );
            std::cout << "wm0Rtm : " <<  wm0Rtm.x << " " << wm0Rtm.y << " " << wm0Rtm.z << std::endl;
            const double mat[9] = {
                wm0Rtm.x*wm0Rtm.x, wm0Rtm.x*wm0Rtm.y, wm0Rtm.x*wm0Rtm.z,
                wm0Rtm.y*wm0Rtm.x, wm0Rtm.y*wm0Rtm.y, wm0Rtm.y*wm0Rtm.z,
                wm0Rtm.z*wm0Rtm.x, wm0Rtm.z*wm0Rtm.y, wm0Rtm.z*wm0Rtm.z,
            };
            ls.AddMatrix(irb,true,3,3,  irb,true,3,3,  mat, dt*newmark_gamma*c);
//          Com::CVector3D res_t = c*rb.Omega;
        }
		const double dtmp1 = dt*dt*newmark_beta;
        else{
            std::cout << "not contact" << std::endl;
            const double val0 = 1;
            const double val1 = lambda[5];
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp2);
            ls.AddResidual(icst,false,5,1, &val1,-1 );
        }*/
		const double dtmp1 = dt*dt*newmark_beta;
        const double dtmp2 = dt*newmark_gamma;
        if( flg == true && is_initial ){
            std::cout << "consider contact" << std::endl;
            double df;
            Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0Rtm = wm0.MatVec( Rtm );
            Com::CVector3D res_t = lambda[5]*wm0Rtm;
            df = Com::Dot(wm0Rtm,rb.Omega);
//            std::cout << " velo f : " << df << std::endl;
            ls.AddResidual(irb, true, 3,   res_t,-1 );
            ls.AddResidual(icst,false,5,1, &df,-1 );
            ls.AddMatrix_Vector(irb, true, 3, icst,false,5, wm0Rtm,dtmp2, true );
            ls.AddMatrix_Vector(icst,false,5, irb, true, 3, wm0Rtm,dtmp2, false);
        }
        else if( !is_initial && lambda[5]<-1.0e-5 ){
//            std::cout << "after impact" << std::endl;
            Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0Rtm = wm0.MatVec( Rtm );
            Com::CVector3D res_t = lambda[5]*wm0Rtm;
            ls.AddResidual(irb, true, 3,   res_t,-1 );
            ////////////////
            double val0 = 0.0;
            ls.AddResidual(icst,false,5,1, &val0,-1 );
            double val1 = 1.0;
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1, &val1, dtmp2);
        }
        else{
//            std::cout << "not contact" << std::endl;
            const double val0 = 1;
            const double val1 = lambda[5];
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp2);
            ls.AddResidual(icst,false,5,1, &val1,-1 );
        }
/*		const double dtmp1 = dt*dt*newmark_beta;
        if( flg == true || lambda[5] < 0  )
        {
            std::cout << "consider contact" << std::endl;
            Com::CVector3D Rtm = mrot.MatVecTrans( v[0]*loc_coord[0] + v[1]*loc_coord[1] );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0Rtm = wm0.MatVec(Rtm);
            Com::CVector3D res_t = lambda[5]*wm0Rtm;
            ls.AddResidual(irb, true, 3,   res_t,-1 );
            ls.AddResidual(icst,false,5,1, &f_value,-1 );
            ls.AddMatrix_Vector(irb, true, 3, icst,false,5, wm0Rtm,dtmp1, true );
            ls.AddMatrix_Vector(icst,false,5, irb, true, 3, wm0Rtm,dtmp1, false);
            Com::CMatrix3 wvRtm( lambda[5]*Rtm );
            Com::CMatrix3 ktt = wm0.MatMat( wvRtm );
            ls.AddMatrix(irb,true,3,  irb,true,3, ktt, dtmp1, true);
        }
        else{
            std::cout << "not contact" << std::endl;
            const double val0 = 1;
            const double val1 = lambda[5];
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp1);
            ls.AddResidual(icst,false,5,1, &val1,-1 );
        }*/
    }
	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity();
        ls.AddMatrix(irb, true, 0, icst,false,0, mtmp,   -dtmp1, true);
        ls.AddMatrix(icst,false,0, irb, true, 0, mtmp,   -dtmp1, true);

        ls.AddMatrix(irb, true, 3, irb, true, 3, wXwRtL,  dtmp1, true);
        ls.AddMatrix(irb, true, 3, irb, true, 3, lwLCwRta,dtmp1, true);

        ls.AddMatrix(irb, true, 3, icst,false,0, RotwX,  -dtmp1, false);
        ls.AddMatrix(icst,false,0, irb, true, 3, RotwX,  -dtmp1, true );

        ls.AddMatrix_Vector(irb, true, 3, icst,false,3, wm0Rta,dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3, irb, true, 3, wm0Rta,dtmp1, false);
        ls.AddMatrix_Vector(irb, true, 3, icst,false,4, wm1Rta,dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4, irb, true, 3, wm1Rta,dtmp1, false);

//        ls.AddMatrix_Vector(icst,false,5, irb, true, 5, 1,dtmp1, false);
	}
}
/*
void Rigid::CFix_HingeRange::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        const unsigned int irb = aIndRB[0];
        const CRigidBody3D& rb = aRB[irb];
        const Com::CVector3D& vec_j = rb.GetPositionFromInital(ini_pos_fix);
        const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.disp_cg;

        ::glPushMatrix();
        ::glTranslated(ini_pos_fix.x, ini_pos_fix.y, ini_pos_fix.z);
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();

        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        ::glVertex3d(vec_cg.x,vec_cg.y,vec_cg.z);
        ::glEnd();

        const Com::CVector3D& lcb0 = this->loc_coord[0];
        const Com::CVector3D& lcb1 = this->loc_coord[1];
        unsigned int ndiv_t = 32;
        unsigned int ndiv0 = ndiv_t*(max_t-min_t      )/360.0 + 1;
        const double dtheta0 = 2*3.1416/ndiv0*(max_t-min_t)/360.0;
        const double radius = 1;
        ::glColor3d(0,1,1);
        ::glBegin(GL_TRIANGLE_FAN);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        for(unsigned int idiv=0;idiv<ndiv0+1;idiv++){
            const Com::CVector3D v0 = vec_j + 
                sin(idiv*dtheta0 - max_t*3.14/180)*lcb0*radius + 
                cos(idiv*dtheta0 - max_t*3.14/180)*lcb1*radius;
            ::glVertex3d(v0.x,  v0.y,  v0.z);
        }
        ::glEnd();
    }
}
*/
////////////////////////////////////////////////////////////////

void Rigid::CJoint_Spherical::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
                     const double dt, const double newmark_gamma, const double newmark_beta, 
                     const std::vector<CRigidBody3D>& aRB, bool is_initial ) const 
{
    const unsigned int irb0 = aIndRB[0];
    assert( irb0 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb0 = aRB[irb0];
    const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
    const Com::CVector3D& Xdistfix0 = rb0.ini_pos_cg-ini_pos_joint;

    const unsigned int irb1 = aIndRB[1];
    assert( irb1 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb1 = aRB[irb1];
    const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
    const Com::CVector3D& Xdistfix1 = rb1.ini_pos_cg-ini_pos_joint;

    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);

    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Oj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // S‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Oj
        Com::CVector3D tmp_vec = mrot0.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        ls.AddResidual(irb0,true,3, wX.MatVec(tmp_vec),-1);
	}
	{   // S‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Oj
        const Com::CVector3D& rot_pos_cg = mrot0.MatVec(Xdistfix0) + ini_pos_joint;
        ls.AddResidual( icst,false,0,   rb0.ini_pos_cg+rb0.disp_cg-rot_pos_cg,1 );
	}
    Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
    Com::CMatrix3 wXwRtL0;
	{
        const Com::CVector3D& RtL = mrot0.MatVec(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Pj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // S‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Pj
        const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        ls.AddResidual(irb1,true,3, wX.MatVec(tmp_vec),1 ); 
	}
	{   // S‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Pj
        const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1) + ini_pos_joint;
        ls.AddResidual(icst,false,0,  rb1.ini_pos_cg+rb1.disp_cg-rot_pos_cg,-1 );
	}
    Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
    Com::CMatrix3 wXwRtL1;
	{
        const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity();
        ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
        ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
        ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
        ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);

        ls.AddMatrix(irb0,true,3, irb0,true,3, wXwRtL0,  dtmp1, true);
        ls.AddMatrix(irb1,true,3, irb1,true,3, wXwRtL1, -dtmp1, true);

        ls.AddMatrix(irb0,true, 3,  icst,false,0,  RotwX0, -dtmp1, false);
        ls.AddMatrix(icst,false,0,  irb0,true, 3,  RotwX0, -dtmp1, true );
        ls.AddMatrix(irb1,true, 3,  icst,false,0,  RotwX1, +dtmp1, false);
        ls.AddMatrix(icst,false,0,  irb1,true, 3,  RotwX1, +dtmp1, true );
	}
}

void Rigid::CJoint_Spherical::UpdateSolution(const double* upd,
    double dt, double newmark_gamma, double newmark_beta)
{
    const double tmp1 = dt*dt*newmark_beta;
//    const double* upd = ls.GetUpdate(icst,false);
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
}
/*
void Rigid::CJoint_Spherical::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        const unsigned int irb0 = aIndRB[0];
        const CRigidBody3D& rb0 = aRB[irb0];
        const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.disp_cg;

        const unsigned int irb1 = aIndRB[1];
        const CRigidBody3D& rb1 = aRB[irb1];
        const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.disp_cg;
    
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);

        ::glPushMatrix();
        ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();
        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
        ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
        ::glEnd();
    }
}
*/
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////

void Rigid::CJoint_Hinge::SetAxis(double ax, double ay, double az){
    axis = Com::CVector3D(ax,ay,az);
    if( axis.Length() < 1.0e-20 ){ assert(0); return; }
    axis.Normalize();
    GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}

/*
void Rigid::CJoint_Hinge::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 2;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        const unsigned int irb0 = aIndRB[0];
        const CRigidBody3D& rb0 = aRB[irb0];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
        ::glColor3d(1,1,0);
        ::glBegin(GL_POINTS);
        ::glVertex3d(vec_j.x, vec_j.y, vec_j.z );
        ::glEnd();
    }
    else if( imode == 2 ){
        const unsigned int irb0 = aIndRB[0];
        const unsigned int irb1 = aIndRB[1];
        const CRigidBody3D& rb0 = aRB[irb0];
        const CRigidBody3D& rb1 = aRB[irb1];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
        const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.disp_cg;
        const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.disp_cg;

        ::glPushMatrix();
        ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();
        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
        ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
        ::glEnd();
    
        const Com::CMatrix3& mrot = rb0.GetRotMatrix();
        const Com::CVector3D& lcb0 = mrot.MatVec(this->loc_coord[0]);
        const Com::CVector3D& lcb1 = mrot.MatVec(this->loc_coord[1]);
        unsigned int ndiv = 16;
        const double dtheta = 2*3.1416/ndiv;
        const double radius = 0.4;
        ::glColor3d(0,1,1);
        ::glBegin(GL_TRIANGLE_FAN);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        for(unsigned int idiv=0;idiv<ndiv+1;idiv++){
            const Com::CVector3D v0 = vec_j + sin(idiv*dtheta  )*lcb0*radius + cos(idiv*dtheta  )*lcb1*radius;
            ::glVertex3d(v0.x,  v0.y,  v0.z);
        }
        ::glEnd();
    }
}
*/

void Rigid::CJoint_Hinge::UpdateSolution(const double* upd,
                    double dt, double newmark_gamma, double newmark_beta)
{
    const double tmp1 = dt*dt*newmark_beta;
//    const double* upd = ls.GetUpdate(icst,false);
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
    lambda[3] += tmp1*upd[3];
    lambda[4] += tmp1*upd[4];
}

void Rigid::CJoint_Hinge::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
    const double dt, const double newmark_gamma, const double newmark_beta,
    const std::vector<CRigidBody3D>& aRB, bool is_initial ) const
{
    const unsigned int irb0 = aIndRB[0];
    assert( irb0 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb0 = aRB[irb0];

    const unsigned int irb1 = aIndRB[1];
    assert( irb1 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb1 = aRB[irb1];

    const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
    const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
    const Com::CVector3D Xdistfix0 = rb0.ini_pos_cg - ini_pos_joint;
    const Com::CVector3D Xdistfix1 = rb1.ini_pos_cg - ini_pos_joint;

    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Oj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // •ÏˆÊS‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Oj
        const Com::CVector3D& tmp_vec = mrot0.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        ls.AddResidual(irb0,true,3,  wX.MatVec(tmp_vec),-1);
	}
	{   // •ÏˆÊS‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Oj
        Com::CVector3D rot_pos_cg = mrot0.MatVec(Xdistfix0) + ini_pos_joint;
        ls.AddResidual(icst,false,0,   rb0.ini_pos_cg+rb0.disp_cg-rot_pos_cg,1 );
	}
    Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
    Com::CMatrix3 wXwRtL0;
	{
        const Com::CVector3D& RtL = mrot0.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Pj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // S‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Pj
        const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        ls.AddResidual(irb1,true,3,  wX.MatVec(tmp_vec),1 );
	}
	{   // S‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Pj
        const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1) + ini_pos_joint;
        ls.AddResidual(icst,false,0,  rb1.ini_pos_cg+rb1.disp_cg-rot_pos_cg,-1);
	}
    Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
    Com::CMatrix3 wXwRtL1;
	{
        const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    
    ////////////////////////////////

    Com::CVector3D Rtwma[2][2];
    {   // ‰ñ“]S‘©‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& R1a = mrot1.MatVec(axis);
        const Com::CVector3D& R0m0 = mrot0.MatVec(loc_coord[0]);
        const Com::CVector3D& R0m1 = mrot0.MatVec(loc_coord[1]);
        const Com::CMatrix3 wR0m0(R0m0);
        const Com::CMatrix3 wR0m1(R0m1);
        const Com::CVector3D& wR0m0R1a = wR0m0.MatVec(R1a);
        const Com::CVector3D& wR0m1R1a = wR0m1.MatVec(R1a);
        Rtwma[0][0] = mrot0.MatVecTrans(wR0m0R1a);
        Rtwma[0][1] = mrot0.MatVecTrans(wR0m1R1a);
        Rtwma[1][0] = mrot1.MatVecTrans(wR0m0R1a);
        Rtwma[1][1] = mrot1.MatVecTrans(wR0m1R1a);
        ls.AddResidual(irb0,true,3, lambda[3]*Rtwma[0][0] + lambda[4]*Rtwma[0][1],-1 );
        ls.AddResidual(irb1,true,3, lambda[3]*Rtwma[1][0] + lambda[4]*Rtwma[1][1], 1 );
        double res[5] = {0,0,0,0,0};
        res[3] = Com::Dot(R1a,R0m0);
        res[4] = Com::Dot(R1a,R0m1);
        ls.SubResidual(icst,false,res);
    }
    ////////////////////////////////////////////////////////////////

	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity();
        ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
        ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
        ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
        ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);

        ls.AddMatrix(irb0,true, 3, irb0,true, 3, wXwRtL0, dtmp1, true);
        ls.AddMatrix(irb1,true, 3, irb1,true, 3, wXwRtL1,-dtmp1, true);

        ls.AddMatrix(irb0,true, 3, icst,false,0,  RotwX0,-dtmp1, false);
        ls.AddMatrix(icst,false,0, irb0,true, 3,  RotwX0,-dtmp1, true );
        ls.AddMatrix(irb1,true, 3, icst,false,0,  RotwX1, dtmp1, false);
        ls.AddMatrix(icst,false,0, irb1,true, 3,  RotwX1, dtmp1, true );

        ls.AddMatrix_Vector(irb0,true, 3,  icst,false,3,  Rtwma[0][0], dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3,  irb0,true, 3,  Rtwma[0][0], dtmp1, false);
        ls.AddMatrix_Vector(irb0,true, 3,  icst,false,4,  Rtwma[0][1], dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4,  irb0,true, 3,  Rtwma[0][1], dtmp1, false);
        ls.AddMatrix_Vector(irb1,true, 3,  icst,false,3,  Rtwma[1][0],-dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3,  irb1,true, 3,  Rtwma[1][0],-dtmp1, false);
        ls.AddMatrix_Vector(irb1,true, 3,  icst,false,4,  Rtwma[1][1],-dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4,  irb1,true, 3,  Rtwma[1][1],-dtmp1, false);
	}
}






////////////////////////////////////////////////////////////////

void Rigid::CJoint_HingeRange::SetAxis(double ax, double ay, double az){
    axis = Com::CVector3D(ax,ay,az);
    if( axis.Length() < 1.0e-20 ){ assert(0); return; }
    axis.Normalize();
    GetVertical2Vector(axis,loc_coord[0],loc_coord[1]);
}
/*
void Rigid::CJoint_HingeRange::Draw(const std::vector<CRigidBody3D>& aRB) const
{
    const unsigned int imode = 2;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        const unsigned int irb0 = aIndRB[0];
        const CRigidBody3D& rb0 = aRB[irb0];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
        ::glColor3d(1,1,0);
        ::glBegin(GL_POINTS);
        ::glVertex3d(vec_j.x, vec_j.y, vec_j.z );
        ::glEnd();
    }
    else if( imode == 2 ){
        const unsigned int irb0 = aIndRB[0];
        const unsigned int irb1 = aIndRB[1];
        const CRigidBody3D& rb0 = aRB[irb0];
        const CRigidBody3D& rb1 = aRB[irb1];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(ini_pos_joint);
        const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.disp_cg;
        const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.disp_cg;

        ::glPushMatrix();
        ::glTranslated( vec_j.x, vec_j.y, vec_j.z );
        ::glColor3d(1,1,0);
        ::glutSolidSphere(0.1,10,5);
	    ::glPopMatrix();
        ::glColor3d(1,1,1);
        ::glLineWidth(2);
        ::glBegin(GL_LINES);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z  );
        ::glVertex3d(vec_cg0.x,vec_cg0.y,vec_cg0.z);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        ::glVertex3d(vec_cg1.x,vec_cg1.y,vec_cg1.z);
        ::glEnd();
    
        const Com::CMatrix3& mrot = rb0.GetRotMatrix();
        const Com::CVector3D& lcb0 = mrot.MatVec(this->loc_coord[0]);
        const Com::CVector3D& lcb1 = mrot.MatVec(this->loc_coord[1]);
        unsigned int ndiv_t = 32;
        unsigned int ndiv0 = ndiv_t*(max_t-min_t)/360.0 + 1;
        const double dtheta = 2*3.1416/ndiv0*(max_t-min_t)/360.0;
        const double radius = 1;
        ::glColor3d(0,1,1);
        ::glBegin(GL_TRIANGLE_FAN);
        ::glVertex3d(vec_j.x,  vec_j.y,  vec_j.z);
        for(unsigned int idiv=0;idiv<ndiv0+1;idiv++){
            const Com::CVector3D v0 = vec_j 
                + sin(idiv*dtheta - max_t*3.1416/180.0 )*lcb0*radius 
                + cos(idiv*dtheta - max_t*3.1416/180.0 )*lcb1*radius;
            ::glVertex3d(v0.x,  v0.y,  v0.z);
        }
        ::glEnd();
    }
}
*/

void Rigid::CJoint_HingeRange::UpdateSolution(const double* upd,
                    double dt, double newmark_gamma, double newmark_beta)
{
    const double tmp1 = dt*dt*newmark_beta;
//    const double* upd = ls.GetUpdate(icst,false);
    lambda[0] += tmp1*upd[0];
    lambda[1] += tmp1*upd[1];
    lambda[2] += tmp1*upd[2];
    lambda[3] += tmp1*upd[3];
    lambda[4] += tmp1*upd[4];
//    lambda[5] += tmp1*upd[5];
    const double tmp2 = dt*newmark_gamma;
    lambda[5] += tmp2*upd[5];
}

void Rigid::CJoint_HingeRange::AddLinearSystem(Ls::CLinearSystem_RigidBody& ls, unsigned int icst,
    const double dt, const double newmark_gamma, const double newmark_beta,
    const std::vector<CRigidBody3D>& aRB, bool is_initial ) const
{
    const unsigned int irb0 = aIndRB[0];
    assert( irb0 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb0 = aRB[irb0];

    const unsigned int irb1 = aIndRB[1];
    assert( irb1 < ls.GetSizeRigidBody() );
    const CRigidBody3D& rb1 = aRB[irb1];

    const Com::CMatrix3& mrot0 = rb0.GetRotMatrix();
    const Com::CMatrix3& mrot1 = rb1.GetRotMatrix();
    const Com::CVector3D vlambda(lambda[0],lambda[1],lambda[2]);
    const Com::CVector3D Xdistfix0 = rb0.ini_pos_cg - ini_pos_joint;
    const Com::CVector3D Xdistfix1 = rb1.ini_pos_cg - ini_pos_joint;

    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Oj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb0,true,0,  vlambda,1);
	{   // •ÏˆÊS‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Oj
        const Com::CVector3D& tmp_vec = mrot0.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        ls.AddResidual(irb0,true,3,  wX.MatVec(tmp_vec),-1);
	}
	{   // •ÏˆÊS‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Oj
        Com::CVector3D rot_pos_cg = mrot0.MatVec(Xdistfix0) + ini_pos_joint;
        ls.AddResidual(icst,false,0,   rb0.ini_pos_cg+rb0.disp_cg-rot_pos_cg,1 );
	}
    Com::CMatrix3 RotwX0 = mrot0.MatMat( Com::CMatrix3(Xdistfix0) );
    Com::CMatrix3 wXwRtL0;
	{
        const Com::CVector3D& RtL = mrot0.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix0);
        wXwRtL0 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    ////////////////////////////////////////////////
    // S‘©—ÍŽc·i•¨‘Ì‚Pj
	// •Ài‚ÌŽc·
    ls.AddResidual(irb1,true,0,  vlambda,-1);
	{   // S‘©—Í‚©‚ç—ˆ‚é‰ñ“]‚ÌŽc·i•¨‘Ì‚Pj
        const Com::CVector3D& tmp_vec = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        ls.AddResidual(irb1,true,3,  wX.MatVec(tmp_vec),1 );
	}
	{   // S‘©ðŒ‚Ì•Ï•ªi•¨‘Ì‚Pj
        const Com::CVector3D& rot_pos_cg = mrot1.MatVec(Xdistfix1) + ini_pos_joint;
        ls.AddResidual(icst,false,0,  rb1.ini_pos_cg+rb1.disp_cg-rot_pos_cg,-1);
	}
    Com::CMatrix3 RotwX1 = mrot1.MatMat( Com::CMatrix3(Xdistfix1) );
    Com::CMatrix3 wXwRtL1;
	{
        const Com::CVector3D& RtL = mrot1.MatVecTrans(vlambda);
        const Com::CMatrix3 wX(Xdistfix1);
        wXwRtL1 = wX.MatMat( Com::CMatrix3(RtL) );
	}   
    
    ////////////////////////////////

    Com::CVector3D Rtwma[2][2];
    {   // ‰ñ“]S‘©‚©‚ç—ˆ‚é•¨‘Ì‰ñ“]‚ÌŽc·
        const Com::CVector3D& R1a = mrot1.MatVec(axis);
        const Com::CVector3D& R0m0 = mrot0.MatVec(loc_coord[0]);
        const Com::CVector3D& R0m1 = mrot0.MatVec(loc_coord[1]);
        const Com::CMatrix3 wR0m0(R0m0);
        const Com::CMatrix3 wR0m1(R0m1);
        const Com::CVector3D& wR0m0R1a = wR0m0.MatVec(R1a);
        const Com::CVector3D& wR0m1R1a = wR0m1.MatVec(R1a);
        Rtwma[0][0] = mrot0.MatVecTrans(wR0m0R1a);
        Rtwma[0][1] = mrot0.MatVecTrans(wR0m1R1a);
        Rtwma[1][0] = mrot1.MatVecTrans(wR0m0R1a);
        Rtwma[1][1] = mrot1.MatVecTrans(wR0m1R1a);
        ls.AddResidual(irb0,true,3, lambda[3]*Rtwma[0][0] + lambda[4]*Rtwma[0][1],-1 );
        ls.AddResidual(irb1,true,3, lambda[3]*Rtwma[1][0] + lambda[4]*Rtwma[1][1], 1 );
        double res[6] = {0,0,0,0,0,0};
        res[3] = Com::Dot(R1a,R0m0);
        res[4] = Com::Dot(R1a,R0m1);
        ls.SubResidual(icst,false,res);
    }
    {   // Šp“x‚É‚Â‚¢‚Ä‚ÌKKTðŒ‚ð“ü‚ê‚é
        bool flg = false;
        double v[2];
        double f_value;
        const double delta_f = 1.0e-8;
        {
            const double mint = -max_t;
            const double maxt = -min_t;
            double R0m0R1m0 = Com::Dot(mrot0.MatVec(loc_coord[0]), mrot1.MatVec(loc_coord[0]));
            double R0m0R1m1 = Com::Dot(mrot0.MatVec(loc_coord[0]), mrot1.MatVec(loc_coord[1]));
            const double PI = 3.14159265358979323846;
            double u_min[2] = { cos(PI*mint/180.0), sin(PI*mint/180.0) };
            double u_max[2] = { cos(PI*maxt/180.0), sin(PI*maxt/180.0) };
            double v_min[2] = { -u_min[1],  u_min[0] };
            double v_max[2] = {  u_max[1], -u_max[0] };
            double dmin = v_min[0]*R0m0R1m0 + v_min[1]*R0m0R1m1;
            double dmax = v_max[0]*R0m0R1m0 + v_max[1]*R0m0R1m1;
//            std::cout << dmin << " " << dmax << " " << lambda[5] << std::endl;
            if( dmin < -delta_f || dmax < -delta_f ) flg = true;
            const double emin = u_min[0]*R0m0R1m0 + u_min[1]*R0m0R1m1;
            const double emax = u_max[0]*R0m0R1m0 + u_max[1]*R0m0R1m1;
            if( emin > emax ){ 
                v[0] = v_min[0]; v[1] = v_min[1];
                f_value = dmin;
            }
            else{                              
                v[0] = v_max[0]; v[1] = v_max[1];
                f_value = dmax;
            }
        }
        const double dtmp2 = dt*newmark_gamma;
        if( flg == true && is_initial )
        {
//            std::cout << "consider contact" << std::endl;
            Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
            ////////////////
            Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
            Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m;
            ls.AddResidual(irb0, true, 3,   res_t0,-1 );
            ls.AddMatrix_Vector(irb0,true, 3, icst,false,5, wm0R0tR1m,dtmp2, true );
            ls.AddMatrix_Vector(icst,false,5, irb0,true, 3, wm0R0tR1m,dtmp2, false);
            ////////////////
            Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
            Com::CMatrix3 wm( m );
            Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
            Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0;
            ls.AddResidual(irb1, true, 3,   res_t1,-1 );
            ls.AddMatrix_Vector(irb1,true, 3, icst,false,5, wmR1tR0m0,dtmp2, true );
            ls.AddMatrix_Vector(icst,false,5, irb1,true, 3, wmR1tR0m0,dtmp2, false);
            ////////////////
            const double df = Com::Dot(wm0R0tR1m,rb0.Omega) + Com::Dot(wmR1tR0m0,rb1.Omega);
            ls.AddResidual(icst,false,5,1, &df,-1 );
        }
        else if( !is_initial && lambda[5]<-1.0e-5 )
        {
//            std::cout << "after impact" << std::endl;
            Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
            ////////////////
            Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
            Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m;
            ls.AddResidual(irb0, true, 3,   res_t0,-1 );
            ////////////////
            Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
            Com::CMatrix3 wm( m );
            Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
            Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0;
            ls.AddResidual(irb1, true, 3,   res_t1,-1 );
            ////////////////
            double val0 = 0.0;
            ls.AddResidual(icst,false,5,1, &val0,-1 );
            double val1 = 1.0;
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1, &val1, dtmp2);
        }
        else{
//            std::cout << "not contact" << std::endl;
            const double val0 = 1;
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp2);
            const double val1 = lambda[5];
            ls.AddResidual(icst,false,5,1, &val1,-1 );
        }
/*		const double dtmp1 = dt*dt*newmark_beta;
        if( flg == true && is_initial )
        {
            std::cout << "consider contact" << std::endl;
            Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
            ////////////////
            Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
            Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m;
            ls.AddResidual(irb0, true, 3,   res_t0,-1 );
            ls.AddMatrix_Vector(irb0,true, 3, icst,false,5, wm0R0tR1m,dtmp1, true );
            ls.AddMatrix_Vector(icst,false,5, irb0,true, 3, wm0R0tR1m,dtmp1, false);
            ////////////////
            Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
            Com::CMatrix3 wm( m );
            Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
            Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0;
            ls.AddResidual(irb1, true, 3,   res_t1,-1 );
            ls.AddMatrix_Vector(irb1,true, 3, icst,false,5, wmR1tR0m0,dtmp1, true );
            ls.AddMatrix_Vector(icst,false,5, irb1,true, 3, wmR1tR0m0,dtmp1, false);
            ////////////////
            ls.AddResidual(icst,false,5,1, &f_value,-1 );
            ////////////////
            {
                Com::CMatrix3 wm0wR0tR1m = wm0.MatMat( Com::CMatrix3(R0tR1m)  );
                ls.AddMatrix(irb0,true,3,  irb0,true,3,  wm0wR0tR1m, dt*dt*newmark_beta, true); 
                Com::CMatrix3 wmwR1tR0m0 = wm.MatMat(  Com::CMatrix3(R1tR0m0) );
                ls.AddMatrix(irb1,true,3,  irb1,true,3,  wmwR1tR0m0, dt*dt*newmark_beta, true); 
//                Com::CMatrix3 wm0R0tR1wm = wm0.MatMat( mrot0.MatMatTrans( mrot1.MatMat( wm ) ) );
//                ls.AddMatrix(irb0,true,3,  irb1,true,3,  wm0R0tR1wm, -dt*dt*newmark_beta, true); 
//                Com::CMatrix3 wmR1tR0wm0 = wm.MatMat( mrot1.MatMatTrans( mrot0.MatMat( wm0 ) ) );
//                ls.AddMatrix(irb1,true,3,  irb0,true,3,  wmR1tR0wm0, -dt*dt*newmark_beta, true); 
            }
        }
        else if( !is_initial && lambda[5]<-1.0e-5 )
        {
            Com::CVector3D m = v[0]*loc_coord[0] + v[1]*loc_coord[1];
            Com::CVector3D R0tR1m = mrot0.MatVecTrans( mrot1.MatVec( m ) );
            Com::CMatrix3 wm0( loc_coord[0] );
            Com::CVector3D wm0R0tR1m = wm0.MatVec( R0tR1m );
            Com::CVector3D res_t0 = lambda[5]*wm0R0tR1m*newmark_gamma;
            ls.AddResidual(irb0, true, 3,   res_t0,-1 );
            Com::CVector3D R1tR0m0 = mrot1.MatVecTrans( mrot0.MatVec( loc_coord[0] ) );
            Com::CMatrix3 wm( m );
            Com::CVector3D wmR1tR0m0 = wm.MatVec( R1tR0m0 );
            Com::CVector3D res_t1 = lambda[5]*wmR1tR0m0*newmark_gamma;
            ls.AddResidual(irb1, true, 3,   res_t1,-1 );
            ////////////////
            double val0 = 0.0;
            ls.AddResidual(icst,false,5,1, &val0,-1 );
            double val1 = 1.0;
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1, &val1, dtmp1);
            ////////////////
            {
                Com::CMatrix3 wm0wR0tR1m = wm0.MatMat( Com::CMatrix3(R0tR1m)  );
                ls.AddMatrix(irb0,true,3,  irb0,true,3,  wm0wR0tR1m, dt*dt*newmark_beta, true); 
                Com::CMatrix3 wmwR1tR0m0 = wm.MatMat(  Com::CMatrix3(R1tR0m0) );
                ls.AddMatrix(irb1,true,3,  irb1,true,3,  wmwR1tR0m0, dt*dt*newmark_beta, true); 
//                Com::CMatrix3 wm0R0tR1wm = wm0.MatMat( mrot0.MatMatTrans( mrot1.MatMat( wm ) ) );
//                ls.AddMatrix(irb0,true,3,  irb1,true,3,  wm0R0tR1wm, -dt*dt*newmark_beta, true); 
//                Com::CMatrix3 wmR1tR0wm0 = wm.MatMat( mrot1.MatMatTrans( mrot0.MatMat( wm0 ) ) );
//                ls.AddMatrix(irb1,true,3,  irb0,true,3,  wmR1tR0wm0, -dt*dt*newmark_beta, true); 
            }
            {
                const double c = 100;
                const double df = Com::Dot(wm0R0tR1m,rb0.Omega) + Com::Dot(wmR1tR0m0,rb1.Omega);
                Com::CVector3D res_t0 = c*df*wm0R0tR1m;
                ls.AddResidual(irb0, true, 3,   res_t0,-1 );
                Com::CVector3D res_t1 = c*df*wmR1tR0m0;
                ls.AddResidual(irb1, true, 3,   res_t1,-1 );
                {
                    const double mat[9] = {
                        wm0R0tR1m.x*wm0R0tR1m.x, wm0R0tR1m.x*wm0R0tR1m.y, wm0R0tR1m.x*wm0R0tR1m.z,
                        wm0R0tR1m.y*wm0R0tR1m.x, wm0R0tR1m.y*wm0R0tR1m.y, wm0R0tR1m.y*wm0R0tR1m.z,
                        wm0R0tR1m.z*wm0R0tR1m.x, wm0R0tR1m.z*wm0R0tR1m.y, wm0R0tR1m.z*wm0R0tR1m.z,
                    };
                    ls.AddMatrix(irb0,true,3,3,  irb0,true,3,3,  mat, dt*newmark_gamma*c);
                }
                {
                    const double mat[9] = {
                        wmR1tR0m0.x*wmR1tR0m0.x, wmR1tR0m0.x*wmR1tR0m0.y, wmR1tR0m0.x*wmR1tR0m0.z,
                        wmR1tR0m0.y*wmR1tR0m0.x, wmR1tR0m0.y*wmR1tR0m0.y, wmR1tR0m0.y*wmR1tR0m0.z,
                        wmR1tR0m0.z*wmR1tR0m0.x, wmR1tR0m0.z*wmR1tR0m0.y, wmR1tR0m0.z*wmR1tR0m0.z,
                    };
                    ls.AddMatrix(irb1,true,3,3,  irb1,true,3,3,  mat, dt*newmark_gamma*c);
                }
            }
        }
        else{
            std::cout << "not contact" << std::endl;
            const double val0 = 1;
            ls.AddMatrix(icst,false,5,1,  icst,false,5,1,  &val0, dtmp1);
            const double val1 = lambda[5];
            ls.AddResidual(icst,false,5,1, &val1,-1 );
        }*/
    }
    ////////////////////////////////////////////////////////////////

	{   // „«s—ñ‚ðì‚é
		const double dtmp1 = dt*dt*newmark_beta;
        Com::CMatrix3 mtmp;
        mtmp.SetIdentity();
        ls.AddMatrix(irb0,true, 0, icst,false,0, mtmp,   -dtmp1, true);
        ls.AddMatrix(icst,false,0, irb0,true, 0, mtmp,   -dtmp1, true);
        ls.AddMatrix(irb1,true, 0, icst,false,0, mtmp,    dtmp1, true);
        ls.AddMatrix(icst,false,0, irb1,true, 0, mtmp,    dtmp1, true);

        ls.AddMatrix(irb0,true, 3, irb0,true, 3, wXwRtL0, dtmp1, true);
        ls.AddMatrix(irb1,true, 3, irb1,true, 3, wXwRtL1,-dtmp1, true);

        ls.AddMatrix(irb0,true, 3, icst,false,0,  RotwX0,-dtmp1, false);
        ls.AddMatrix(icst,false,0, irb0,true, 3,  RotwX0,-dtmp1, true );
        ls.AddMatrix(irb1,true, 3, icst,false,0,  RotwX1, dtmp1, false);
        ls.AddMatrix(icst,false,0, irb1,true, 3,  RotwX1, dtmp1, true );

        ls.AddMatrix_Vector(irb0,true, 3,  icst,false,3,  Rtwma[0][0], dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3,  irb0,true, 3,  Rtwma[0][0], dtmp1, false);
        ls.AddMatrix_Vector(irb0,true, 3,  icst,false,4,  Rtwma[0][1], dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4,  irb0,true, 3,  Rtwma[0][1], dtmp1, false);
        ls.AddMatrix_Vector(irb1,true, 3,  icst,false,3,  Rtwma[1][0],-dtmp1, true );
        ls.AddMatrix_Vector(icst,false,3,  irb1,true, 3,  Rtwma[1][0],-dtmp1, false);
        ls.AddMatrix_Vector(irb1,true, 3,  icst,false,4,  Rtwma[1][1],-dtmp1, true );
        ls.AddMatrix_Vector(icst,false,4,  irb1,true, 3,  Rtwma[1][1],-dtmp1, false);
	}
}


