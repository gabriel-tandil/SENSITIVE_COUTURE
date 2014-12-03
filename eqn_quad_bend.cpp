/*
 *  eqn_quadric_bend.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 12/1/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */


#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/linearsystem_interface_eqnsys.h"
#include "delfem/vector3d.h"

#include "emat_cst.h"
#include "emat_quad_bend.h"
#include "eqn_quad_bend.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace MatVec;

bool AddLinSys_QuadBend_CST_BackWardEular
(double dt, double damp_coeff,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_bend, double lambda, double myu,
 double rho, double g_x, double g_y, double g_z, 
 const unsigned int id_field_disp, const unsigned int id_field_hinge,
 const Fem::Field::CFieldWorld& world, 
 double& kinetic_energy,
 double& strain_energy,
 double& potential_energy)
{
  kinetic_energy   = 0;
  strain_energy    = 0;
  potential_energy = 0;
	const Fem::Field::CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
  MatVec::CMatDia_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
  MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);
	if( stiff_bend > 1.0e-30 ){
		const CField& field_hinge = world.GetField(id_field_hinge);
		std::vector<unsigned int> aIdEA = field_hinge.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry::CElemSeg& es = field_hinge.GetElemSeg(id_ea,CORNER,true,world);
			for(unsigned int ielem=0;ielem<es.Size();ielem++){
				unsigned int no[4];	es.GetNodes(ielem,no);
        //                std::cout << iiea << " " << ielem << "  -->  " << no[0] << " " << no[1] << " " << no[2] << " " << no[3] << std::endl;
				double C[4][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);	ns_c.GetValue(no[2],C[2]);	ns_c.GetValue(no[3],C[3]);
				double u[4][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);	ns_u.GetValue(no[2],u[2]);	ns_u.GetValue(no[3],u[3]);
				double v[4][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);	ns_v.GetValue(no[2],v[2]);	ns_v.GetValue(no[3],v[3]);
				double Kmat[4][4][3][3], Res[4][3];
				GetMatRes_QuadBend_BackwardEular(Kmat,Res, C,u,v, stiff_bend, dt, 
                                         strain_energy);
				mat.Mearge(4,no, 4,no, 9, &Kmat[0][0][0][0]);
				for(unsigned int ino=0;ino<4;ino++){
					res.AddValue(no[ino],0,-Res[ino][0]);
					res.AddValue(no[ino],1,-Res[ino][1]);
					res.AddValue(no[ino],2,-Res[ino][2]);
				}
			}
		}
	}
	////////////////
	{
		const CField& field_disp = world.GetField(id_field_disp);
		std::vector<unsigned int> aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry::CElemSeg& es = field_disp.GetElemSeg(id_ea,CORNER,true,world);
			for(unsigned int ielem=0;ielem<es.Size();ielem++){
				unsigned int no[3];	es.GetNodes(ielem,no);
				double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);	ns_c.GetValue(no[2],C[2]);
				double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);	ns_u.GetValue(no[2],u[2]);
				double v[3][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);	ns_v.GetValue(no[2],v[2]);
				double Kmat[3][3][3][3], Res[3][3];
				GetMatRes_MassCST_BackwardEular(Kmat,Res, C,u,v, lambda,myu,
                                        rho,g_x,g_y,g_z, 
                                        dt, damp_coeff,
                                        kinetic_energy,strain_energy,potential_energy);
				mat.Mearge(3,no, 3,no, 9, &Kmat[0][0][0][0]);
				for(unsigned int ino=0;ino<3;ino++){
					res.AddValue(no[ino],0,-Res[ino][0]);
					res.AddValue(no[ino],1,-Res[ino][1]);
					res.AddValue(no[ino],2,-Res[ino][2]);
				}
			}
		}
	}
	return true;
}
 
static void GetdFdC_Tri
(double dFdC[3][3][3][3],
 const double C[3][3],
 const double rho, double gx, double gy, double gz)
{  
  double Gd[3][3] = {
		{ C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] },
		{ C[2][0]-C[0][0], C[2][1]-C[0][1], C[2][2]-C[0][2] }, { 0,0,0 } };
  double Area;
  Com::UnitNormalAreaTri3D(Gd[2], Area, C[0], C[1], C[2]);
  
  double DAdC[3][3];
  {
    const double D[3][3] = {
      { C[2][0]-C[1][0], C[2][1]-C[1][1], C[2][2]-C[1][2] },
      { C[0][0]-C[2][0], C[0][1]-C[2][1], C[0][2]-C[2][2] },
      { C[1][0]-C[0][0], C[1][1]-C[0][1], C[1][2]-C[0][2] } };
    Com::Cross3D(DAdC[0], Gd[2], D[0]);
    Com::Cross3D(DAdC[1], Gd[2], D[1]);
    Com::Cross3D(DAdC[2], Gd[2], D[2]);
  }  
  for(unsigned int jno=0;jno<3;jno++){
    for(unsigned int jdim=0;jdim<3;jdim++){
      const double tmp_r = DAdC[jno][jdim]*0.5/3.0*rho;
      dFdC[0][jno][0][jdim] = gx*tmp_r;
      dFdC[0][jno][1][jdim] = gy*tmp_r;
      dFdC[0][jno][2][jdim] = gz*tmp_r;
      dFdC[1][jno][0][jdim] = gx*tmp_r;
      dFdC[1][jno][1][jdim] = gy*tmp_r;
      dFdC[1][jno][2][jdim] = gz*tmp_r;
      dFdC[2][jno][0][jdim] = gx*tmp_r;
      dFdC[2][jno][1][jdim] = gy*tmp_r;
      dFdC[2][jno][2][jdim] = gz*tmp_r;    
    }
  }
}

bool AddLinSys_QuadBend_CST_Sensitivity_FictitousBending
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_bend, double lambda, double myu,
 double rho, double g_x, double g_y, double g_z,
 const unsigned int id_field_disp, const unsigned int id_field_hinge,
 const unsigned int id_field_senseX, const unsigned int id_field_senseY,
 const unsigned int id_field_lamX, const unsigned int id_field_lamY,
 Fem::Field::CFieldWorld& world)
{
	const Fem::Field::CField& field_lx = world.GetField(id_field_lamX);
	const Fem::Field::CNodeAry::CNodeSeg& ns_lx = field_lx.GetNodeSeg(CORNER,true, world,VALUE);
  
	const Fem::Field::CField& field_ly = world.GetField(id_field_lamY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_ly = field_ly.GetNodeSeg(CORNER,true, world,VALUE);  
  
	Fem::Field::CField& field_sx = world.GetField(id_field_senseX);
	Fem::Field::CNodeAry::CNodeSeg& ns_sx = field_sx.GetNodeSeg(CORNER,true, world,VALUE);  
  ns_sx.SetZero();
  Fem::Field::CField& field_sy = world.GetField(id_field_senseY);
	Fem::Field::CNodeAry::CNodeSeg& ns_sy = field_sy.GetNodeSeg(CORNER,true, world,VALUE);
  ns_sy.SetZero();  
  
	const Fem::Field::CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
  MatVec::CMatDia_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
  MatVec::CVector_Blk&    res = ls.GetResidual(id_field_disp,CORNER,world);
	if( stiff_bend > 1.0e-30 ){
		const CField& field_hinge = world.GetField(id_field_hinge);
		std::vector<unsigned int> aIdEA = field_hinge.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry::CElemSeg& es = field_hinge.GetElemSeg(id_ea,CORNER,true,world);
			for(unsigned int ielem=0;ielem<es.Size();ielem++){
				unsigned int no[4];	es.GetNodes(ielem,no);
				double C[4][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);	ns_c.GetValue(no[2],C[2]);	ns_c.GetValue(no[3],C[3]);
				double u[4][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);	ns_u.GetValue(no[2],u[2]);	ns_u.GetValue(no[3],u[3]);
        double c[4][3] = {
          { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
          { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
          { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] },          
          { C[3][0]+u[3][0], C[3][1]+u[3][1], C[3][2]+u[3][2] } };
				double Kmat[4][4][3][3], Res[4][3], dRdC[4][4][3][3];
				GetKmatResdRdC_QuadBend(Kmat,Res,dRdC, c,c,stiff_bend);
				mat.Mearge(4,no, 4,no, 9, &Kmat[0][0][0][0]);
        for(unsigned int ino=0;ino<4;ino++){
          res.AddValue(no[ino],0,Res[ino][0]);
          res.AddValue(no[ino],1,Res[ino][1]);
          res.AddValue(no[ino],2,Res[ino][2]);          
        }
			}
		}
	}
	////////////////
	{
		const CField& field_disp = world.GetField(id_field_disp);
		std::vector<unsigned int> aIdEA = field_disp.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry::CElemSeg& es = field_disp.GetElemSeg(id_ea,CORNER,true,world);
			for(unsigned int ielem=0;ielem<es.Size();ielem++){
				unsigned int no[3];	es.GetNodes(ielem,no);
				double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);	ns_c.GetValue(no[2],C[2]);
				double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);	ns_u.GetValue(no[2],u[2]);
        double c[3][3] = {
          { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
          { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },
          { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };
        double dRdC[3][3][3][3];        
        {
          double Kmat[3][3][3][3],Res[3][3];
          GetKmatResdRdC_CST(Kmat,Res,dRdC, C,c, lambda,myu);        
          mat.Mearge(3,no, 3,no, 9, &Kmat[0][0][0][0]);          
          for(unsigned int ino=0;ino<3;ino++){
            res.AddValue(no[ino],0,Res[ino][0]);
            res.AddValue(no[ino],1,Res[ino][1]);
            res.AddValue(no[ino],2,Res[ino][2]);          
          }                  
        }
        {
          double dFdC[3][3][3][3];
          GetdFdC_Tri(dFdC,C, rho, g_x, g_y, g_z);
          for(unsigned int i=0;i<81;i++){ (&dRdC[0][0][0][0])[i] -= (&dFdC[0][0][0][0])[i]; }
        }          
        ////
        double lx[3][2];	ns_lx.GetValue(no[0],lx[0]);	ns_lx.GetValue(no[1],lx[1]);	ns_lx.GetValue(no[2],lx[2]);
        for(unsigned int ino=0;ino<3;ino++){
          double rx[3] = {0,0,0};
          for(unsigned int jno=0;jno<3;jno++){
            rx[0] -= (dRdC[ino][jno][0][0]*lx[jno][0] + dRdC[ino][jno][0][1]*lx[jno][1]);
            rx[1] -= (dRdC[ino][jno][1][0]*lx[jno][0] + dRdC[ino][jno][1][1]*lx[jno][1]);
            rx[2] -= (dRdC[ino][jno][2][0]*lx[jno][0] + dRdC[ino][jno][2][1]*lx[jno][1]);
          }
          ns_sx.AddValue(no[ino],0,rx[0]);
          ns_sx.AddValue(no[ino],1,rx[1]);
          ns_sx.AddValue(no[ino],2,rx[2]);          
        }            
        double ly[3][2];	ns_ly.GetValue(no[0],ly[0]);	ns_ly.GetValue(no[1],ly[1]);	ns_ly.GetValue(no[2],ly[2]);        
        for(unsigned int ino=0;ino<3;ino++){           
          double ry[3] = {0,0,0};
          for(unsigned int jno=0;jno<3;jno++){
            ry[0] -= (dRdC[ino][jno][0][0]*ly[jno][0] + dRdC[ino][jno][0][1]*ly[jno][1]);
            ry[1] -= (dRdC[ino][jno][1][0]*ly[jno][0] + dRdC[ino][jno][1][1]*ly[jno][1]);
            ry[2] -= (dRdC[ino][jno][2][0]*ly[jno][0] + dRdC[ino][jno][2][1]*ly[jno][1]);
          }
          ns_sy.AddValue(no[ino],0,ry[0]);
          ns_sy.AddValue(no[ino],1,ry[1]);
          ns_sy.AddValue(no[ino],2,ry[2]);          
        }                            
			}
		}
	}
	return true;  
}


