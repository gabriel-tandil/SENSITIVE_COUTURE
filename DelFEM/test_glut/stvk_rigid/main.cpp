
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif
#define for if(0);else for

#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <deque>

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesh3d.h"

#include "delfem/ls/solver_ls_iter.h"
#include "delfem/matvec/mat_blkcrs.h"
#include "delfem/matvec/matdia_blkcrs.h"

#include "delfem/rigid/rigidbody.h"
#include "delfem/rigid/linearsystem_rigid.h"
#include "delfem/rigid/linearsystem_rigidfield.h"

#include "delfem/femeqn/eqn_linear_solid3d.h"
#include "delfem/femeqn/eqn_st_venant.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"
#include "delfem/drawer_field_face.h"

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
int imodifier;



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using namespace Fem::Field;
using namespace MatVec;

Fem::Field::View::CDrawerArrayField drawer_ary;

bool AddLinSys_AbovePlane_NonStatic_NewmarkBeta(
	double gamma, double beta, double dt,  bool is_inital,
	Fem::Eqn::ILinearSystem_Eqn& ls,
	unsigned int id_field_disp, 
	unsigned int id_field_lambda, 
	Fem::Field::CFieldWorld& world )
{
	if( !world.IsIdField(id_field_disp) ) return false;
	const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != Fem::Field::VECTOR3 ) return false;
	////////////////
	if( !world.IsIdField(id_field_lambda) ) return false;
	const Fem::Field::CField& field_lambda = world.GetField(id_field_lambda);
	if( field_lambda.GetFieldType() != Fem::Field::VECTOR3 ) return false;

	unsigned int id_na_co = field_lambda.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_na_lambda = field_lambda.GetNodeSegInNodeAry(CORNER).id_na_va;

	const unsigned int ndim = 3;

	CMat_BlkCrs& pmat_dl = ls.GetMatrix(id_field_disp,  CORNER, id_field_lambda,CORNER,world);
	assert( pmat_dl.LenBlkCol() == ndim );
	assert( pmat_dl.LenBlkRow() == ndim );
	CMat_BlkCrs& pmat_ld = ls.GetMatrix(id_field_lambda,CORNER, id_field_disp  ,CORNER,world);
	assert( pmat_ld.LenBlkCol() == ndim );
	assert( pmat_ld.LenBlkRow() == ndim );
	CMatDia_BlkCrs& pmat_ll = ls.GetMatrix(id_field_lambda,CORNER,world);
	assert( pmat_ll.LenBlkCol() == ndim );
	assert( pmat_ll.LenBlkRow() == ndim );
	CVector_Blk& res_l = ls.GetResidual(id_field_lambda,CORNER,world); 
	CVector_Blk& res_d = ls.GetResidual(id_field_disp,CORNER,world); 

	double coord[ndim];
	double udisp[ndim], vdisp[ndim], adisp[ndim];
	double ulambda[ndim], vlambda[ndim], alambda[ndim];
	const CNodeAry& na_l = world.GetNA(id_na_lambda);
	const CNodeAry::CNodeSeg& ns_co      = field_disp.GetNodeSeg(  CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_udisp   = field_disp.GetNodeSeg(  CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vdisp   = field_disp.GetNodeSeg(  CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_adisp   = field_disp.GetNodeSeg(  CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_ulambda = field_lambda.GetNodeSeg(CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vlambda = field_lambda.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_alambda = field_lambda.GetNodeSeg(CORNER,true, world,ACCELERATION);

	for(unsigned int inode=0;inode<na_l.Size();inode++)
	{
		const unsigned int inode_co = field_lambda.GetMapVal2Co(inode);
		ns_co.GetValue(inode_co,coord);
		ns_udisp.GetValue(inode_co,udisp);
		ns_vdisp.GetValue(inode_co,vdisp);
		ns_adisp.GetValue(inode_co,adisp);

		ns_ulambda.GetValue(inode,ulambda);
		ns_vlambda.GetValue(inode,vlambda);
		ns_alambda.GetValue(inode,alambda);

		double eKmat_dl[ndim][ndim],eKmat_ld[ndim][ndim], eKmat_ll[ndim][ndim];
		double eres_d[ndim], eres_l[ndim];

        const double normal[3] = { 0, 0, -1 };
        const double cur_pos[3] = { udisp[0]+coord[0], udisp[1]+coord[1], udisp[2]+coord[2]+0.3 };
        const double f_value = cur_pos[0]*normal[0]+cur_pos[1]*normal[1]+cur_pos[2]*normal[2];

        if( ulambda[0]  > 0.000001 || f_value > 0.0000001 ){
			eKmat_dl[0][0] = normal[0];   eKmat_dl[0][1] = 0;           eKmat_dl[0][2] = 0;
			eKmat_dl[1][0] = normal[1];   eKmat_dl[1][1] = 0;           eKmat_dl[1][2] = 0;
			eKmat_dl[2][0] = normal[2];   eKmat_dl[2][1] = 0;           eKmat_dl[2][2] = 0;
			eKmat_ld[0][0] = normal[0];   eKmat_ld[0][1] = normal[1];   eKmat_ld[0][2] = normal[2];
			eKmat_ld[1][0] = 0;           eKmat_ld[1][1] = 0;           eKmat_ld[1][2] = 0;
			eKmat_ld[2][0] = 0;           eKmat_ld[2][1] = 0;           eKmat_ld[2][2] = 0;
			eKmat_ll[0][0] = 0;           eKmat_ll[0][1] = 0;           eKmat_ll[0][2] = 0;
			eKmat_ll[1][0] = 0;           eKmat_ll[1][1] = 1;           eKmat_ll[1][2] = 0;
			eKmat_ll[2][0] = 0;           eKmat_ll[2][1] = 0;           eKmat_ll[2][2] = 1;
			eres_d[0] = -normal[0]*ulambda[0];   eres_d[1] = -normal[1]*ulambda[0];   eres_d[2] = -normal[2]*ulambda[0];
			eres_l[0] = -f_value;             eres_l[1] = 0;                    eres_l[2] = 0;
		}
        else{
			eKmat_dl[0][0] = 0;   eKmat_dl[0][1] = 0;   eKmat_dl[0][2] = 0;
			eKmat_dl[1][0] = 0;   eKmat_dl[1][1] = 0;   eKmat_dl[1][2] = 0;
			eKmat_dl[2][0] = 0;   eKmat_dl[2][1] = 0;   eKmat_dl[2][2] = 0;
			eKmat_ld[0][0] = 0;   eKmat_ld[0][1] = 0;   eKmat_ld[0][2] = 0;
			eKmat_ld[1][0] = 0;   eKmat_ld[1][1] = 0;   eKmat_ld[1][2] = 0;
			eKmat_ld[2][0] = 0;   eKmat_ld[2][1] = 0;   eKmat_ld[2][2] = 0;
			eKmat_ll[0][0] = 1;   eKmat_ll[0][1] = 0;   eKmat_ll[0][2] = 0;
			eKmat_ll[1][0] = 0;   eKmat_ll[1][1] = 1;   eKmat_ll[1][2] = 0;
			eKmat_ll[2][0] = 0;   eKmat_ll[2][1] = 0;   eKmat_ll[2][2] = 1;
			eres_d[0] = 0;           eres_d[1] = 0;   eres_d[2] = 0;
			eres_l[0] = -ulambda[0]; eres_l[1] = 0;   eres_l[2] = 0;
        }
		////////////////
		double emat_dl[3][3], emat_ld[3][3], emat_ll[3][3];
		{
			const double dtmp = dt*dt*beta;
            for(unsigned int i=0;i<9;i++){
                (&emat_dl[0][0])[i] = dtmp*(&eKmat_dl[0][0])[i];
                (&emat_ld[0][0])[i] = dtmp*(&eKmat_ld[0][0])[i];
                (&emat_ll[0][0])[i] = dtmp*(&eKmat_ll[0][0])[i];
            }
		}
		if( is_inital ){
            for(unsigned int i=0;i<ndim;i++){
            for(unsigned int j=0;j<ndim;j++){
			    eres_d[i] -= dt*eKmat_dl[i][j]*vlambda[j];
                eres_l[i] -= dt*(eKmat_ld[i][j]*vdisp[j]+eKmat_ll[i][j]*vlambda[j]);
            }
            }
            for(unsigned int i=0;i<ndim;i++){
            for(unsigned int j=0;j<ndim;j++){
                eres_d[i] -= 0.5*dt*dt*eKmat_dl[i][j]*alambda[j];
                eres_l[i] -= 0.5*dt*dt*(eKmat_ld[i][j]*adisp[j]+eKmat_ll[i][j]*alambda[j]);
            }
            }
		}
		pmat_dl.Mearge(1,&inode_co,  1,&inode,     9,  &emat_dl[0][0]);
		pmat_ld.Mearge(1,&inode   ,  1,&inode_co,  9,  &emat_ld[0][0]);
		pmat_ll.Mearge(1,&inode,     1,&inode,     9,  &emat_ll[0][0]);
		res_d.AddValue(inode_co,0,eres_d[0]);
		res_d.AddValue(inode_co,1,eres_d[1]);
		res_d.AddValue(inode_co,2,eres_d[2]);
		res_l.AddValue(inode,   0,eres_l[0]);
		res_l.AddValue(inode,   1,eres_l[1]);
		res_l.AddValue(inode,   2,eres_l[2]);
	}

	return true;
}

bool AddLinSys_ConnectRigid_NonStatic_NewmarkBeta(
	double gamma, double beta, double dt,  bool is_initial,
	Ls::CLinearSystem_RigidField2& ls,
	unsigned int id_field_disp, unsigned int id_field_lambda, Fem::Field::CFieldWorld& world, 
    unsigned int irb, const std::vector<Rigid::CRigidBody3D>& aRB, const std::vector<Rigid::CConstraint*>& aConst )
{
	if( !world.IsIdField(id_field_disp) ) return false;
	const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != Fem::Field::VECTOR3 ) return false;
	////////////////
	if( !world.IsIdField(id_field_lambda) ) return false;
	const Fem::Field::CField& field_lambda = world.GetField(id_field_lambda);
	if( field_lambda.GetFieldType() != Fem::Field::VECTOR3 ) return false;

	unsigned int id_na_co = field_lambda.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_na_lambda = field_lambda.GetNodeSegInNodeAry(CORNER).id_na_va;

	const unsigned int ndim = 3;
    const int ilss_d = ls.FindIndexArray_Seg(id_field_disp,  CORNER,world);
    const int ilss_l = ls.FindIndexArray_Seg(id_field_lambda,CORNER,world);
    const int ilss_r = ls.GetIndexSegRigid();

	CMat_BlkCrs& pmat_dl = ls.GetMatrix(id_field_disp,  CORNER, id_field_lambda,CORNER,world);
	assert( pmat_dl.LenBlkCol() == ndim );
	assert( pmat_dl.LenBlkRow() == ndim );
	CMat_BlkCrs& pmat_ld = ls.GetMatrix(id_field_lambda,CORNER, id_field_disp  ,CORNER,world);
	assert( pmat_ld.LenBlkCol() == ndim );
	assert( pmat_ld.LenBlkRow() == ndim );
	CMatDia_BlkCrs& pmat_ll = ls.GetMatrix(id_field_lambda,CORNER,world);
	assert( pmat_ll.LenBlkCol() == ndim );
	assert( pmat_ll.LenBlkRow() == ndim );
	CMat_BlkCrs& pmat_lr = ls.GetMatrix(ilss_l,ilss_r);
	CMat_BlkCrs& pmat_rl = ls.GetMatrix(ilss_r,ilss_l);
	CMatDia_BlkCrs& pmat_rr = ls.GetMatrix(ilss_r);
	CVector_Blk& res_l = ls.GetResidual(id_field_lambda,CORNER,world); 
	CVector_Blk& res_d = ls.GetResidual(id_field_disp,CORNER,world);
	CVector_Blk& res_r = ls.GetResidual(ilss_r);

	const CNodeAry& na_l = world.GetNA(id_na_lambda);
	const CNodeAry::CNodeSeg& ns_co      = field_disp.GetNodeSeg(  CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_udisp   = field_disp.GetNodeSeg(  CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vdisp   = field_disp.GetNodeSeg(  CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_adisp   = field_disp.GetNodeSeg(  CORNER,true, world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_ulambda = field_lambda.GetNodeSeg(CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vlambda = field_lambda.GetNodeSeg(CORNER,true, world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_alambda = field_lambda.GetNodeSeg(CORNER,true, world,ACCELERATION);

    const Rigid::CRigidBody3D& rigid = aRB[irb];

	for(unsigned int inode=0;inode<na_l.Size();inode++)
	{
		const unsigned int inode_co = field_lambda.GetMapVal2Co(inode);

        Com::CVector3D Co0, Co1, L1;
        {
	        double coord[ndim];
	        double udisp[ndim], vdisp[ndim], adisp[ndim];
	        double ulambda[ndim], vlambda[ndim], alambda[ndim];
            ////////////////
		    ns_co.GetValue(inode_co,coord);
		    ns_udisp.GetValue(inode_co,udisp);
            Co0 = Com::CVector3D(coord[0],coord[1],coord[2]);
            Co1 = Com::CVector3D(coord[0]+udisp[0],coord[1]+udisp[1],coord[2]+udisp[2]);
            if( is_initial ){
		        ns_vdisp.GetValue(inode_co,vdisp);
		        ns_adisp.GetValue(inode_co,adisp);
                Co1.x += dt*vdisp[0] + 0.5*dt*dt*adisp[0];
                Co1.y += dt*vdisp[1] + 0.5*dt*dt*adisp[1];
                Co1.z += dt*vdisp[2] + 0.5*dt*dt*adisp[2];
            }
            ////////////////
		    ns_ulambda.GetValue(inode,ulambda);
            L1 = Com::CVector3D(ulambda[0],ulambda[1],ulambda[2]);
            if( is_initial ){
		        ns_vlambda.GetValue(inode,vlambda);
		        ns_alambda.GetValue(inode,alambda);
                L1.x += dt*vlambda[0] + 0.5*dt*dt*alambda[0];
                L1.y += dt*vlambda[1] + 0.5*dt*dt*alambda[1];
                L1.z += dt*vlambda[2] + 0.5*dt*dt*alambda[2];
            }
        }

        const Com::CVector3D D1 = Co0-rigid.ini_pos_cg;
        const Com::CMatrix3 rot = rigid.GetRotMatrix();
        Com::CVector3D f_value;
        {
            const Com::CVector3D pos1_X = rigid.ini_pos_cg + rigid.disp_cg;
            const Com::CVector3D RD1 = rot.MatVec(D1);
            f_value = Co1 - pos1_X - RD1;
        }
        const Com::CMatrix3 wD1(D1);
        const Com::CVector3D Rtl = rot.MatVecTrans( L1 );
        const Com::CVector3D wD1Rtl = wD1.MatVec(Rtl);
        const Com::CMatrix3 RwD1 = rot.MatMat( wD1 );
        const Com::CMatrix3 wD1wRtl = wD1.MatMat( Com::CMatrix3(Rtl) );

		double eKmat_dl[ndim][ndim],eKmat_ld[ndim][ndim];

        {
			eKmat_dl[0][0] = 1; eKmat_dl[0][1] = 0; eKmat_dl[0][2] = 0;
			eKmat_dl[1][0] = 0; eKmat_dl[1][1] = 1; eKmat_dl[1][2] = 0;
			eKmat_dl[2][0] = 0; eKmat_dl[2][1] = 0; eKmat_dl[2][2] = 1;
        }
        {
			eKmat_ld[0][0] = 1; eKmat_ld[0][1] = 0; eKmat_ld[0][2] = 0;
			eKmat_ld[1][0] = 0; eKmat_ld[1][1] = 1; eKmat_ld[1][2] = 0;
			eKmat_ld[2][0] = 0; eKmat_ld[2][1] = 0; eKmat_ld[2][2] = 1;
        }
        double eKmat_lr[ndim][6], eKmat_rl[6][ndim];
        {
			eKmat_lr[0][0] =-1;             eKmat_lr[0][1] = 0;             eKmat_lr[0][2] = 0;
			eKmat_lr[1][0] = 0;             eKmat_lr[1][1] =-1;             eKmat_lr[1][2] = 0;
			eKmat_lr[2][0] = 0;             eKmat_lr[2][1] = 0;             eKmat_lr[2][2] =-1;
			eKmat_lr[0][3] = RwD1.mat[0];   eKmat_lr[0][4] = RwD1.mat[1];   eKmat_lr[0][5] = RwD1.mat[2];
			eKmat_lr[1][3] = RwD1.mat[3];   eKmat_lr[1][4] = RwD1.mat[4];   eKmat_lr[1][5] = RwD1.mat[5];
			eKmat_lr[2][3] = RwD1.mat[6];   eKmat_lr[2][4] = RwD1.mat[7];   eKmat_lr[2][5] = RwD1.mat[8];
        }
        {
			eKmat_rl[0][0] =-1;             eKmat_rl[0][1] = 0;             eKmat_rl[0][2] = 0;
			eKmat_rl[1][0] = 0;             eKmat_rl[1][1] =-1;             eKmat_rl[1][2] = 0;
			eKmat_rl[2][0] = 0;             eKmat_rl[2][1] = 0;             eKmat_rl[2][2] =-1;
			eKmat_rl[3][0] = RwD1.mat[0];   eKmat_rl[3][1] = RwD1.mat[3];   eKmat_rl[3][2] = RwD1.mat[6];
			eKmat_rl[4][0] = RwD1.mat[1];   eKmat_rl[4][1] = RwD1.mat[4];   eKmat_rl[4][2] = RwD1.mat[7];
			eKmat_rl[5][0] = RwD1.mat[2];   eKmat_rl[5][1] = RwD1.mat[5];   eKmat_rl[5][2] = RwD1.mat[8];
        }
        double eKmat_rr[6][6];
        for(unsigned int i=0;i<36;i++){ (&eKmat_rr[0][0])[i] = 0.0; }
        {
            eKmat_rr[3][3] = -wD1wRtl.mat[0];   eKmat_rr[3][4] = -wD1wRtl.mat[1];   eKmat_rr[3][5] = -wD1wRtl.mat[2];
            eKmat_rr[4][3] = -wD1wRtl.mat[3];   eKmat_rr[4][4] = -wD1wRtl.mat[4];   eKmat_rr[4][5] = -wD1wRtl.mat[5];
            eKmat_rr[5][3] = -wD1wRtl.mat[6];   eKmat_rr[5][4] = -wD1wRtl.mat[7];   eKmat_rr[5][5] = -wD1wRtl.mat[8];
        }
		double eres_d[ndim], eres_l[ndim], eres_r[6];
        {
			eres_d[0] = -L1.x;      eres_d[1] = -L1.y;      eres_d[2] = -L1.z;
            eres_r[0] =  L1.x;      eres_r[1] =  L1.y;      eres_r[2] =  L1.z;
            eres_r[3] =  wD1Rtl.x;  eres_r[4] =  wD1Rtl.y;  eres_r[5] =  wD1Rtl.z;
			eres_l[0] = -f_value.x; eres_l[1] = -f_value.y; eres_l[2] = -f_value.z;
		}
		////////////////
		double emat_dl[3][3], emat_ld[3][3];
		double emat_lr[3][6], emat_rl[6][3];
        double emat_rr[6][6];
		{
			const double dtmp = dt*dt*beta;
            for(unsigned int i=0;i<9;i++){
                (&emat_dl[0][0])[i] = dtmp*(&eKmat_dl[0][0])[i];
                (&emat_ld[0][0])[i] = dtmp*(&eKmat_ld[0][0])[i];
            }
            for(unsigned int i=0;i<18;i++){
                (&emat_rl[0][0])[i] = dtmp*(&eKmat_rl[0][0])[i];
                (&emat_lr[0][0])[i] = dtmp*(&eKmat_lr[0][0])[i];
            }
            for(unsigned int i=0;i<36;i++){
                (&emat_rr[0][0])[i] = dtmp*(&eKmat_rr[0][0])[i];
            }
		}
		pmat_dl.Mearge(1,&inode_co,  1,&inode,     9,  &emat_dl[0][0]);
		pmat_ld.Mearge(1,&inode   ,  1,&inode_co,  9,  &emat_ld[0][0]);
		pmat_lr.Mearge(1,&inode   ,  1,&irb,      18,  &emat_lr[0][0]);
		pmat_rl.Mearge(1,&irb   ,    1,&inode,    18,  &emat_rl[0][0]);
		pmat_rr.Mearge(1,&irb   ,    1,&irb,      36,  &emat_rr[0][0]);
		res_d.AddValue(inode_co,0,eres_d[0]);
		res_d.AddValue(inode_co,1,eres_d[1]);
		res_d.AddValue(inode_co,2,eres_d[2]);
		res_l.AddValue(inode,   0,eres_l[0]);
		res_l.AddValue(inode,   1,eres_l[1]);
		res_l.AddValue(inode,   2,eres_l[2]);
        for(unsigned int i=0;i<6;i++){
		    res_r.AddValue(irb,i,eres_r[i]);
        }
	}

	return true;
}

bool is_animation = true;

const double dt = 0.03;
const double newmark_gamma = 0.6;
const double newmark_beta = 0.25*(0.5+newmark_gamma)*(0.5+newmark_gamma);

const Com::CVector3D gravity( 0, 0.0, -9.8 );

std::vector<Rigid::CRigidBody3D> aRB;
std::vector<Rigid::CConstraint*> apFix;
/*
void StepTime2()
{
    ls_rf.InitializeMarge();
    ls_rf.UpdateValueOfRigidSystem(aRB,apFix,   dt,newmark_gamma,newmark_beta, true);
    double norm_res0;
	for(unsigned int itr=0;itr<10;itr++){
        ////////////////
		ls_rf.InitializeMarge();
        for(unsigned int irb=0;irb<aRB.size();irb++){
            aRB[irb].AddLinearSystem(ls_rf,irb,      
                dt,newmark_gamma,newmark_beta,      gravity );
        }
        for(unsigned int ifix=0;ifix<apFix.size();ifix++){
            apFix[ifix]->AddLinearSystem(ls_rf,ifix,         
                dt,newmark_gamma,newmark_beta,       aRB);
        }
        {
			Fem::Eqn::AddLinSys_StVenant3D_NonStatic_NewmarkBeta(
				dt, newmark_gamma, newmark_beta,
				ls_rf.GetLinearSystemField(),
				0.00, 70.0,
				0.1, gravity.x, gravity.y, gravity.z,
				world, id_disp, 
				itr==0);
        }
        {
			AddLinSys_AbovePlane_NonStatic_NewmarkBeta(
				newmark_gamma, newmark_beta, dt, itr==0,
				ls_rf.GetLinearSystemField(), id_disp, id_lambda,
				world );
        }
	    const double res = ls_rf.FinalizeMarge();
        ////////////////
        if( res < 1.0e-30 ) break;
        if( itr==0 ){ norm_res0 = res; }
        std::cout << "itr : " << itr << "     Residual : " << res << std::endl;
        ////////////////
//        {
//            ls_rf.GetLinearSystemRigid().COPY(-1,-2);
//            prec_rf.SetValue(ls_rf);
//            ls_rf.Solve(prec_rf);
//        }
        {
    		double conv_ratio = 1.0e-6;
		    unsigned int max_iter = 1000;
            prec_rf.SetValue(ls_rf);
            Ls::CLinearSystemPreconditioner_RigidField lsp(ls_rf,prec_rf);
            Sol::Solve_PBiCGSTAB(conv_ratio,max_iter,lsp);
//            Sol::Solve_BiCGSTAB(conv_ratio,max_iter,ls_rf);
            std::cout << "       solver itr : " << max_iter << "  conv : " << conv_ratio << std::endl;
        }
        ////////////////////////////////
        ls_rf.GetLinearSystemRigid().UpdateValueOfRigidSystem(aRB,apFix,dt,newmark_gamma, newmark_beta, false);
		ls_rf.UpdateValueOfField_NewmarkBeta(newmark_gamma,newmark_beta,dt,
			id_disp  ,world,itr==0);
		ls_rf.UpdateValueOfField_NewmarkBeta(newmark_gamma,newmark_beta,dt,
			id_lambda,world,itr==0);
        drawer_ary.Update(world);
        if( res < norm_res0*1.0e-8 ) return;
	}
}
*/
Ls::CLinearSystem_RigidField2 ls_rf2;
Ls::CPreconditioner_RigidField2 prec_rf2;
Fem::Field::CFieldWorld world;
double cur_time = 0;

class RigidElasticConnection{
public:
  RigidElasticConnection(unsigned int irigid, unsigned int id_lambda, unsigned int id_disp){
    this->irigid = irigid;
    this->id_lambda = id_lambda;
    this->id_disp = id_disp;
  }
  unsigned int id_disp;
  unsigned int id_lambda;
  unsigned int irigid;  
};
std::vector<unsigned int> aIdDisp;
std::vector< RigidElasticConnection > aRigidElastic;
//unsigned int id_lambda;


void StepTime4()
{
//    Ls::CLinearSystem_RigidField2 ls_rf2;
//    ls_rf2.SetRigidSystem(aRB,apFix);
    ////////////////
    ls_rf2.InitializeMarge();
    ls_rf2.UpdateValueOfRigidSystem(aRB,apFix,dt,newmark_gamma,newmark_beta,     true);
    double norm_res0;
	for(unsigned int itr=0;itr<10;itr++){
        ls_rf2.InitializeMarge();
        for(unsigned int irb=0;irb<aRB.size();irb++){
            aRB[irb].AddLinearSystem(ls_rf2,irb,      dt,newmark_gamma,newmark_beta,      gravity, itr==0 );
        }
        for(unsigned int ifix=0;ifix<apFix.size();ifix++){
            apFix[ifix]->AddLinearSystem(ls_rf2,ifix,         dt,newmark_gamma,newmark_beta,       aRB, itr==0 );
        }
        for(unsigned int idisp=0;idisp<aIdDisp.size();idisp++){
            unsigned int id_disp = aIdDisp[idisp];
			Fem::Eqn::AddLinSys_StVenant3D_NonStatic_NewmarkBeta(
				dt, newmark_gamma, newmark_beta,
				ls_rf2,
				0.00, 2000.0,
				0.3, gravity.x, gravity.y, gravity.z,
				world, id_disp, 
				itr==0);
        }
/*        {
			AddLinSys_AbovePlane_NonStatic_NewmarkBeta(
				newmark_gamma, newmark_beta, dt, itr==0,
				ls_rf2, id_disp, id_lambda,
				world );
        }*/
        for(unsigned int ire=0;ire<aRigidElastic.size();ire++){
            const unsigned int irb = aRigidElastic[ire].irigid;
            const unsigned int id_lambda = aRigidElastic[ire].id_lambda;
            const unsigned int id_disp = aRigidElastic[ire].id_disp;
			AddLinSys_ConnectRigid_NonStatic_NewmarkBeta(
				newmark_gamma, newmark_beta, dt, itr==0,
				ls_rf2, 
                id_disp, id_lambda, world, 
                irb, aRB, apFix );
        }
        ////////////////////////////////
        const double normres = ls_rf2.FinalizeMarge();
        if( normres < 1.0e-30 ) break;
        if( itr==0 ){ norm_res0 = normres; }
        std::cout << "itr : " << itr << "     Residual : " << normres << std::endl;
        /*
        ls_rf2.ReSizeTmpVecSolver(1);
        ls_rf2.COPY(-1,0);
        ls_rf2.COPY(-1,-2);
        {
//            Ls::CPreconditioner_RigidField2 prec_rf2;
//            prec_rf2.SetLinearSystem(ls_rf2);
            prec_rf2.SetValue(ls_rf2);
            ls_rf2.COPY(-1,-2);
            prec_rf2.SolvePrecond( ls_rf2, -2 );
        }
        ls_rf2.MATVEC(-1.0, -2, 1, 0);
        std::cout << "dot : " << ls_rf2.DOT(0,0) << std::endl;
        */
        {
//            Ls::CPreconditioner_RigidField2 prec_rf2;
//            prec_rf2.SetLinearSystem(ls_rf2);
    		double conv_ratio = 1.0e-6;
		    unsigned int max_iter = 100;
            prec_rf2.SetValue(ls_rf2);
            Ls::CLinearSystemPreconditioner_RigidField2 lsp2(ls_rf2,prec_rf2);
            LsSol::Solve_PBiCGSTAB(conv_ratio,max_iter,lsp2);
//            Sol::Solve_BiCGSTAB(conv_ratio,max_iter,ls_rf);
            std::cout << "       solver itr : " << max_iter << "  conv : " << conv_ratio << std::endl;
        }
        ////////////////////////////////
        for(unsigned int idisp=0;idisp<aIdDisp.size();idisp++){
            unsigned int id_disp = aIdDisp[idisp];
		    ls_rf2.UpdateValueOfField_NewmarkBeta(newmark_gamma,newmark_beta,dt,
    			id_disp  ,world,itr==0);
        }
        for(unsigned int ire=0;ire<aRigidElastic.size();ire++){
            const unsigned int id_lambda = aRigidElastic[ire].id_lambda;
		    ls_rf2.UpdateValueOfField_NewmarkBeta(newmark_gamma,newmark_beta,dt,
			    id_lambda,world,itr==0);
        }
        ls_rf2.UpdateValueOfRigidSystem(aRB,apFix,   dt,newmark_gamma,newmark_beta,     false);
        drawer_ary.Update(world);
        if( normres < norm_res0*1.0e-8 ) return;
	}
}



void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

////////////////////////////////

void DrawRigidBody(const Rigid::CRigidBody3D& r){
    const unsigned int imode = 1;
    if( imode == 0 ){
    }
    else if( imode == 1 ){
        ::glPushMatrix();
        ::glTranslated( 
            r.ini_pos_cg.x + r.disp_cg.x,
            r.ini_pos_cg.y + r.disp_cg.y,
            r.ini_pos_cg.z + r.disp_cg.z );
        {
            double rot0[16];
			r.GetInvRotMatrix44(rot0);
            ::glMultMatrixd(rot0);
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

void DrawConstraint(const Rigid::CConstraint* c, const std::vector<Rigid::CRigidBody3D>& aRB)
{	
	if(      const Rigid::CFix_Hinge* cfh = dynamic_cast<const Rigid::CFix_Hinge*>(c) ){
//		std::cout << "CFix_Hinge" << std::endl;
        const unsigned int irb = cfh->aIndRB[0];
        const Rigid::CRigidBody3D& rb = aRB[irb];
        const Com::CVector3D& vec_j = rb.GetPositionFromInital(cfh->ini_pos_fix);
        const Com::CVector3D& vec_cg = rb.ini_pos_cg + rb.disp_cg;
		const Com::CVector3D& ini_pos_fix = cfh->ini_pos_fix;

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

        const Com::CVector3D& lcb0 = cfh->loc_coord[0];
        const Com::CVector3D& lcb1 = cfh->loc_coord[1];
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
		return;
	}
	else if( const Rigid::CFix_HingeRange* cfhr = dynamic_cast<const Rigid::CFix_HingeRange*>(c) ){
//		std::cout << "CFix_HingeRange" << std::endl;
        const unsigned int irb = cfhr->aIndRB[0];
        const Rigid::CRigidBody3D& rb = aRB[irb];
		const Com::CVector3D& ini_pos_fix = cfhr->ini_pos_fix;
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

		const double max_t=cfhr->max_t;
		const double min_t=cfhr->min_t;
        const Com::CVector3D& lcb0 = cfhr->loc_coord[0];
        const Com::CVector3D& lcb1 = cfhr->loc_coord[1];
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
		return;
	}
	else if( const Rigid::CFix_Spherical* cfs = dynamic_cast<const Rigid::CFix_Spherical*>(c) ){
//		std::cout << "CFix_HingeSpherical" << std::endl;
        const unsigned int irb = cfs->aIndRB[0];
        const Rigid::CRigidBody3D& rb = aRB[irb];
		const Com::CVector3D& ini_pos_fix = cfs->ini_pos_fix;
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
		return;
	}
	else if( const Rigid::CJoint_Hinge* cjh = dynamic_cast<const Rigid::CJoint_Hinge*>(c) ){
//		std::cout << "CJoint_Hinge" << std::endl;
        const unsigned int irb0 = cjh->aIndRB[0];
        const unsigned int irb1 = cjh->aIndRB[1];
        const Rigid::CRigidBody3D& rb0 = aRB[irb0];
        const Rigid::CRigidBody3D& rb1 = aRB[irb1];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(cjh->ini_pos_joint);
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
        const Com::CVector3D& lcb0 = mrot.MatVec(cjh->loc_coord[0]);
        const Com::CVector3D& lcb1 = mrot.MatVec(cjh->loc_coord[1]);
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
		return;
	}
	else if( const Rigid::CJoint_HingeRange* cjhr = dynamic_cast<const Rigid::CJoint_HingeRange*>(c) ){
//		std::cout << "CJoint_HingeRange" << std::endl;
        const unsigned int irb0 = cjhr->aIndRB[0];
        const unsigned int irb1 = cjhr->aIndRB[1];
        const Rigid::CRigidBody3D& rb0 = aRB[irb0];
        const Rigid::CRigidBody3D& rb1 = aRB[irb1];
        const Com::CVector3D& vec_j = rb0.GetPositionFromInital(cjhr->ini_pos_joint);
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
        const Com::CVector3D& lcb0 = mrot.MatVec(cjhr->loc_coord[0]);
        const Com::CVector3D& lcb1 = mrot.MatVec(cjhr->loc_coord[1]);
        unsigned int ndiv_t = 32;
		const double max_t=cjhr->max_t;
		const double min_t=cjhr->min_t;
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
		return;
	}
	else if( const Rigid::CJoint_Spherical* cjs = dynamic_cast<const Rigid::CJoint_Spherical*>(c) ){
//		std::cout << "CJoint_Spherical" << std::endl;			
		const unsigned int irb0 = cjs->aIndRB[0];
		const Rigid::CRigidBody3D& rb0 = aRB[irb0];
		const Com::CVector3D& vec_cg0 = rb0.ini_pos_cg + rb0.disp_cg;

		const unsigned int irb1 = cjs->aIndRB[1];
		const Rigid::CRigidBody3D& rb1 = aRB[irb1];
		const Com::CVector3D& vec_cg1 = rb1.ini_pos_cg + rb1.disp_cg;
    
		const Com::CVector3D& vec_j = rb0.GetPositionFromInital(cjs->ini_pos_joint);

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
		return;
	}
}

void myGlutDisplay(void)
{
	::glClearColor(0.2f, .7f, 0.7f, 1.0f);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);
  
  if( is_animation ){
    // 剛体の運動を解く
    //        StepTime();
    //        StepTime2();
    StepTime4();
    cur_time += dt;
    std::cout << cur_time << std::endl;
  }
  
  ShowBackGround();
  
  ::glColor3d(0.7,0.7,0.7);
  ::glBegin(GL_QUADS);
  ::glVertex3d(-2,-2,0);
  ::glVertex3d( 2,-2,0);
  ::glVertex3d( 2, 2,0);
  ::glVertex3d(-2, 2,0);
  ::glEnd();
  {
    ::glLineWidth(1);
    ::glBegin(GL_LINES);
    ::glColor3d(1,0,0);
    ::glVertex3d(0,0,0);
    ::glVertex3d(1,0,0);
    ::glColor3d(0,1,0);
    ::glVertex3d(0,0,0);
    ::glVertex3d(0,1,0);
    ::glColor3d(0,0,1);
    ::glVertex3d(0,0,0);
    ::glVertex3d(0,0,1);
    ::glEnd();
  }
  
  for(unsigned int irb=0;irb<aRB.size();irb++){
    //	    aRB[irb].Draw();
		DrawRigidBody(aRB[irb]);
  }
  for(unsigned int ifix=0;ifix<apFix.size();ifix++){
		DrawConstraint(apFix[ifix],aRB);
    //	    apFix[ifix]->Draw(aRB);
  }
  
  drawer_ary.Draw();
  
  ShowFPS();
  
	glutSwapBuffers();  
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	if( imodifier == GLUT_ACTIVE_CTRL ){
		camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	else if( imodifier == GLUT_ACTIVE_SHIFT ){
		camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
  imodifier = ::glutGetModifiers();
}





void SetProblem()
{
	const unsigned int nprob = 3;	// 問題数
	static int iprob = 0;

    for(unsigned int icst=0;icst<apFix.size();icst++){ delete apFix[icst]; }
    apFix.clear(); 
    aRB.clear();
    cur_time = 0;

/*    if( iprob == 0 ){
        double tot_len = 3.0;
        const unsigned int nRB = 6; 
        aRB.resize(nRB);
        const double div_len = tot_len / nRB;
        for(unsigned int irb=0;irb<nRB;irb++){
            aRB[irb].ini_pos_cg.SetVector(div_len*(irb+1),0,0);
            if( irb == 0 ){
                Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(irb);
                pFix->SetIniPosFix(0,0,0);
                apFix.push_back( pFix );
            }
            else{
                Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(irb-1,irb);
                pFix->SetIniPosJoint(div_len*(irb+0.5),0,0);
                apFix.push_back( pFix );
            }
        }
    }*/
	if( iprob == 0 )
	{
        aRB.resize(1);
        aRB[0].ini_pos_cg.SetVector(-1.0,0,0);      aRB[0].mass = 0.5;
//        aRB[1].ini_pos_cg.SetVector( 3.0,0.0,0);    aRB[1].mass = 0.5;
//        aRB[2].ini_pos_cg.SetVector( 2.0,0.0,0);    aRB[2].mass = 0.5;
//        aRB[3].ini_pos_cg.SetVector( 1.0,0.0,0);    aRB[3].mass = 0.5;
        {
            Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
            pFix->SetIniPosFix(-2,0,0);
            apFix.push_back( pFix );
        }
/*        {
            Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(1);
            pFix->SetIniPosFix(1,1,0);
            apFix.push_back( pFix );
        }*/
/*        {
            Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(1,2);
            pFix->SetIniPosJoint(2.5,0,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Hinge* pFix = new Rigid::CJoint_Hinge(2,3);
            pFix->SetIniPosJoint(1.5,0,0);
            pFix->SetAxis(0,1,1);
            apFix.push_back( pFix );
        }*/

        ////////////////

		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
            std::vector<Com::CVector2D> vec_ary;
            vec_ary.push_back( Com::CVector2D(-1,-0.0) );
            vec_ary.push_back( Com::CVector2D( 2,-0.0) );
            vec_ary.push_back( Com::CVector2D( 2, 0.5) );
            vec_ary.push_back( Com::CVector2D(-1, 0.5) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
        Msh::CMesher2D mesh2d = Msh::CMesher2D(cad_2d,0.2);
        Msh::CMesh3D_Extrude mesh3d;
        mesh3d.Extrude(mesh2d,0.5,0.2);
		const unsigned int id_base = world.AddMesh( mesh3d );
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		////////////////		
		const unsigned int id_disp = world.MakeField_FieldElemDim(id_base,3,VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
        aIdDisp.push_back(id_disp);
/*        unsigned int id_disp_fix0;
        {
		    std::vector<unsigned int> aIdEAFix;
            aIdEAFix.push_back( conv.GetIdEA_fromCad(2,1,2) );
            id_disp_fix0 = world.GetPartialField(id_disp,aIdEAFix);
        }*/
		const unsigned int id_lambda1 = world.MakeField_FieldElemAry(id_base, conv.GetIdEA_fromCad(4,Cad::EDGE,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
//		const unsigned int id_lambda2 = world.MakeField_ElemAry( conv.GetIdEA_fromCad(2,1,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);

		////////////////
		ls_rf2.Clear();
		ls_rf2.AddPattern_Field(id_disp,world);
        ls_rf2.SetRigidSystem(aRB,apFix);
		ls_rf2.AddPattern_RigidField(id_lambda1,id_disp,world,0,aRB,apFix);
//		ls_rf2.AddPattern_RigidField(id_lambda2,id_disp,world,3,aRB,apFix);
//		ls_rf2.AddPattern_Field(id_lambda,id_disp,world);

        aRigidElastic.clear();
        aRigidElastic.push_back( RigidElasticConnection(0,id_lambda1,id_disp) );
        aIdDisp.clear();
        aIdDisp.push_back(id_disp);

//		ls_rf2.SetFixedBoundaryCondition_Field(id_disp_fix0,world);
        prec_rf2.Clear();
        prec_rf2.SetFillInLevel(1);
        prec_rf2.SetLinearSystem(ls_rf2);
        
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerEdge(id_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_disp,false,world) );
	}
	else if( iprob == 1 )
	{
        aRB.resize(2);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        aRB[0].mass = 0.1;
        aRB[1].ini_pos_cg.SetVector(1.5,1.0,0);
        {
            Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
            pFix->SetIniPosFix(0,0,0);
            apFix.push_back( pFix );
        }
        {
            Rigid::CJoint_Hinge* pFix = new Rigid::CJoint_Hinge(0,1);
//            Rigid::CJoint_Spherical* pFix = new Rigid::CJoint_Spherical(0,1);
            pFix->SetIniPosJoint(0,1,0);
            pFix->SetAxis(0,0,1);
            apFix.push_back( pFix );
        }

        ////////////////

		Cad::CCadObj2D cad_2d;
 		{	// 形を作る
            std::vector<Com::CVector2D> vec_ary;
            vec_ary.push_back( Com::CVector2D(1.5,-0.5) );
            vec_ary.push_back( Com::CVector2D(5,  -0.5) );
            vec_ary.push_back( Com::CVector2D(5,   1.0) );
            vec_ary.push_back( Com::CVector2D(1.5, 1.0) );
			cad_2d.AddPolygon( vec_ary );
		}
		world.Clear();
        Msh::CMesher2D mesh2d = Msh::CMesher2D(cad_2d,0.3);
        Msh::CMesh3D_Extrude mesh3d;
        mesh3d.Extrude(mesh2d,0.3,0.3);
		const unsigned int id_base = world.AddMesh( mesh3d );
        const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);

		////////////////		
		const unsigned int id_disp = world.MakeField_FieldElemDim(id_base,3,VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
        aIdDisp.push_back(id_disp);
        unsigned int id_disp_fix0;
        {
		    std::vector<unsigned int> aIdEAFix;
            aIdEAFix.push_back( conv.GetIdEA_fromCad(2,Cad::EDGE,2) );
            id_disp_fix0 = world.GetPartialField(id_disp,aIdEAFix);
        }
		const unsigned int id_lambda = world.MakeField_FieldElemAry(id_base,conv.GetIdEA_fromCad(4,Cad::EDGE,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);

		////////////////
		ls_rf2.Clear();
		ls_rf2.AddPattern_Field(id_disp,world);
        ls_rf2.SetRigidSystem(aRB,apFix);
		ls_rf2.AddPattern_RigidField(id_lambda,id_disp,world,1,aRB,apFix);
//		ls_rf2.AddPattern_Field(id_lambda,id_disp,world);
        
        aRigidElastic.clear();            
        aRigidElastic.push_back( RigidElasticConnection(1,id_lambda,id_disp) );
        aIdDisp.clear();
        aIdDisp.push_back(id_disp);

//		ls_rf2.SetFixedBoundaryCondition_Field(id_disp_fix0,world);
        prec_rf2.Clear();
        prec_rf2.SetFillInLevel(1);
        prec_rf2.SetLinearSystem(ls_rf2);
        
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerEdge(id_disp,false,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_disp,false,world) );
	}
	else if( iprob == 2 )
	{
        aRB.resize(3);
        aRB[0].ini_pos_cg.SetVector(1.0,0,0);
        aRB[0].mass = 0.1;
        aRB[1].ini_pos_cg.SetVector(3.0,0,0);
        aRB[1].mass = 0.1;
        aRB[2].ini_pos_cg.SetVector(4.0,0,0);
        aRB[2].mass = 0.1;
        {
//            Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
            Rigid::CFix_HingeRange* pFix = new Rigid::CFix_HingeRange(0);
            pFix->SetIniPosFix(0,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-30,30);
            apFix.push_back( pFix );
        }
        {
//            Rigid::CFix_Spherical* pFix = new Rigid::CFix_Spherical(0);
            Rigid::CJoint_HingeRange* pFix = new Rigid::CJoint_HingeRange(1,2);
            pFix->SetIniPosJoint(3.5,0,0);
            pFix->SetAxis(0,1,0);
            pFix->SetRange(-45,45);
            apFix.push_back( pFix );
        }

        ////////////////

        unsigned int id_base1, id_base2;
		world.Clear();
        {
		    Cad::CCadObj2D cad_2d;
 		    {	// 形を作る
                std::vector<Com::CVector2D> vec_ary;
                vec_ary.push_back( Com::CVector2D(1.0,-0.5) );
                vec_ary.push_back( Com::CVector2D(3.0,-0.5) );
                vec_ary.push_back( Com::CVector2D(3.0, 0.5) );
                vec_ary.push_back( Com::CVector2D(1.0, 0.5) );
			    cad_2d.AddPolygon( vec_ary );
		    }
            Msh::CMesher2D mesh2d = Msh::CMesher2D(cad_2d,0.3);
            Msh::CMesh3D_Extrude mesh3d;
            mesh3d.Extrude(mesh2d,0.3,0.3);
		    id_base1 = world.AddMesh( mesh3d );
        }
        {
		    Cad::CCadObj2D cad_2d;
 		    {	// 形を作る
                std::vector<Com::CVector2D> vec_ary;
                vec_ary.push_back( Com::CVector2D(4.0,-0.5) );
                vec_ary.push_back( Com::CVector2D(6.0,-0.5) );
                vec_ary.push_back( Com::CVector2D(6.0, 0.5) );
                vec_ary.push_back( Com::CVector2D(4.0, 0.5) );
			    cad_2d.AddPolygon( vec_ary );
		    }
            Msh::CMesher2D mesh2d = Msh::CMesher2D(cad_2d,0.3);
            Msh::CMesh3D_Extrude mesh3d;
            mesh3d.Extrude(mesh2d,0.3,0.3);
		    id_base2 = world.AddMesh( mesh3d );
        }

		////////////////		
		const unsigned int id_disp1 = world.MakeField_FieldElemDim(id_base1,3, VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
        const Fem::Field::CIDConvEAMshCad& conv1 = world.GetIDConverter(id_base1);
        unsigned int id_disp1_fix0;
        {
		    std::vector<unsigned int> aIdEAFix;
            aIdEAFix.push_back( conv1.GetIdEA_fromCad(2,Cad::EDGE,2) );
            id_disp1_fix0 = world.GetPartialField(id_disp1,aIdEAFix);
        }
		const unsigned int id_lambda1 = world.MakeField_FieldElemAry(id_base1,conv1.GetIdEA_fromCad(4,Cad::EDGE,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
		const unsigned int id_lambda2 = world.MakeField_FieldElemAry(id_base1,conv1.GetIdEA_fromCad(2,Cad::EDGE,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);

		////////////////		
		const unsigned int id_disp2 = world.MakeField_FieldElemDim(id_base2,3, VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);
        const Fem::Field::CIDConvEAMshCad& conv2 = world.GetIDConverter(id_base2);
        unsigned int id_disp2_fix0;
        {
		    std::vector<unsigned int> aIdEAFix;
            aIdEAFix.push_back( conv2.GetIdEA_fromCad(2,Cad::EDGE,2) );
            id_disp2_fix0 = world.GetPartialField(id_disp2,aIdEAFix);
        }
		const unsigned int id_lambda3 = world.MakeField_FieldElemAry(id_base2,conv2.GetIdEA_fromCad(4,Cad::EDGE,2), VECTOR3, VALUE|VELOCITY|ACCELERATION, CORNER);

		////////////////
		ls_rf2.Clear();
		ls_rf2.AddPattern_Field(id_disp1,world);
		ls_rf2.AddPattern_Field(id_disp2,world);
        ls_rf2.SetRigidSystem(aRB,apFix);
		ls_rf2.AddPattern_RigidField(id_lambda1,id_disp1,world,0,aRB,apFix);
		ls_rf2.AddPattern_RigidField(id_lambda2,id_disp1,world,1,aRB,apFix);
		ls_rf2.AddPattern_RigidField(id_lambda3,id_disp2,world,2,aRB,apFix);
        
        aRigidElastic.clear();
        aRigidElastic.push_back( RigidElasticConnection(0,id_lambda1,id_disp1) );
        aRigidElastic.push_back( RigidElasticConnection(1,id_lambda2,id_disp1) );
        aRigidElastic.push_back( RigidElasticConnection(2,id_lambda3,id_disp2) );
        aIdDisp.clear();
        aIdDisp.push_back(id_disp1);
        aIdDisp.push_back(id_disp2);

        prec_rf2.Clear();
        prec_rf2.SetFillInLevel(1);
        prec_rf2.SetLinearSystem(ls_rf2);
        
		// 描画オブジェクトの登録
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerEdge(id_disp1,false,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_disp1,false,world) );

		drawer_ary.PushBack( new View::CDrawerEdge(id_disp2,false,world) );
		drawer_ary.PushBack( new View::CDrawerFace(id_disp2,false,world) );
	}

	iprob++;
	if( iprob == nprob ){
		iprob = 0;
	}

  camera.SetRotationMode(Com::View::ROT_3D);
}


void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' は ESC の ASCII コード */
	  exit(0);
	  break;
  case 'a':
      is_animation = !is_animation;
	  break;
  case 's':
      StepTime4();
      break;
  case 'c':
      for(unsigned int irb=0;irb<aRB.size();irb++){
          aRB[irb].Clear();
      }
      for(unsigned int icst=0;icst<apFix.size();icst++){
          apFix[icst]->Clear();
      }
      break;
  case ' ':
      SetProblem();
  default:
    break;
  }
}

void myGlutSpecialFunc(int key, int x, int y){
  switch(key){
    case GLUT_KEY_PAGE_UP :     
      camera.SetScale( camera.GetScale()*0.9    );
      break;
    case GLUT_KEY_PAGE_DOWN :
      camera.SetScale( camera.GetScale()*1.1111 );
      break;
  }
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
}

int main(int argc,char* argv[])
{
	// Initailze GLUT
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 300);
	::glutInit(&argc, argv);	
    ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("Cad View");

	// Set callback function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
  ::glutSpecialFunc(myGlutSpecialFunc);
	::glutIdleFunc(myGlutIdle);

  camera.SetRotationMode(Com::View::ROT_3D);
//    camera.SetIsPers(true);
  {
    Com::CBoundingBox3D bb(-2,2,-2,2,-2,2);
    camera.SetObjectBoundingBox(bb);
    camera.Fit(bb);      
  }
  SetProblem();

	// Enter main loop
	::glutMainLoop();
	return 0;
}
