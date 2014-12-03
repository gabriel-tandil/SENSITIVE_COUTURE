/*
 *  eqn_contact3d.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 7/21/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <math.h>


#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/vector3d.h"

#include "eqn_contact3d.h"
#include "contact_target.h"


using namespace Fem::Ls;
using namespace Fem::Field;

void Update_FrictionalContact
(const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp,
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint )
{	
	assert( world.IsIdField(id_field_disp) );
	const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
	assert( field_disp.GetFieldType() == Fem::Field::VECTOR3 );
	////////////////
	
	const unsigned int ndim = 3;
	const CNodeAry::CNodeSeg& ns_co      = field_disp.GetNodeSeg(  CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_udisp   = field_disp.GetNodeSeg(  CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vdisp   = field_disp.GetNodeSeg(  CORNER,true, world,VELOCITY);
	
	assert( aFrictionPoint.size() == ns_co.Size() );
	for(unsigned int inode=0;inode<ns_co.Size();inode++)
	{
		double Co[ndim];	ns_co.GetValue(   inode,Co);
		double ud[ndim];	ns_udisp.GetValue(inode,ud);
		double uv[ndim];	ns_vdisp.GetValue(inode,uv);
		
		double co[3] = { Co[0]+ud[0], Co[1]+ud[1], Co[2]+ud[2] };
		
		double n1[3];
		const double pd1 = ct.Projection(co[0],co[1],co[2],n1)+offset;
		
		CFrictionPoint& fp = aFrictionPoint[inode];
		const double pd0 = fp.pd;

    if( fp.is_pin ) continue;
		
		if( pd1 < 0 ){} // not contact
		else if( pd0 < 0 ){	// h1>0 && h0<0	: contact in next step
			//			std::cout << "contact next time step" << std::endl;
			// put anchor at projected point
			for(unsigned int i=0;i<3;i++){ fp.aloc[i] = co[i]+pd1*n1[i]; }
		}
		else{	// h1>0 && h0>0
			// update the anchor that spring give force equal to dynamic friction direction to the velocity
			if( fp.itype_contact == 2 ){	
				// update anchor
				double v_t[3];
				{	// tangent direction
					const double t = Com::Dot3D(n1,uv);
					for(unsigned int i=0;i<3;i++){ v_t[i] = uv[i] - t*n1[i]; }
				}
				const double len_vt = Com::Length3D(v_t);
				const double invlen_vt = 1.0/len_vt;
				for(unsigned int i=0;i<3;i++){ v_t[i] *= invlen_vt; }
				const double dist = pd1*stiff_n*myu_k/stiff_f;
				for(unsigned int i=0;i<3;i++){ fp.aloc[i] = co[i]-v_t[i]*dist; }
			}
		}
	}
}

bool AddLinSys_FrictionalContact_Penalty_NonStatic_Sensitivity
(Fem::Ls::CLinearSystem_Field& ls, 
 const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp, 
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint )
{
	if( !world.IsIdField(id_field_disp) ) return false;
	const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != Fem::Field::VECTOR3 ) return false;
	////////////////
	
	MatVec::CMatDia_BlkCrs& pmat_dd = ls.GetMatrix(id_field_disp,  CORNER,world);
	MatVec::CVector_Blk& res_d = ls.GetResidual(id_field_disp,  CORNER,world);   
	
	const unsigned int ndim = 3;
	const CNodeAry::CNodeSeg& ns_co      = field_disp.GetNodeSeg(  CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_udisp   = field_disp.GetNodeSeg(  CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vdisp   = field_disp.GetNodeSeg(  CORNER,true, world,VELOCITY);
	
	assert( aFrictionPoint.size() == ns_co.Size() );
	for(unsigned int inode=0;inode<ns_co.Size();inode++)
	{
		double Co[ndim];	ns_co.GetValue(   inode,Co);
		double ud[ndim];	ns_udisp.GetValue(inode,ud);
		double co[3] = { Co[0]+ud[0], Co[1]+ud[1], Co[2]+ud[2] };
    CFrictionPoint& fp = aFrictionPoint[inode];    
    double n0[3];
    const double pd = ct.Projection(co[0],co[1],co[2], n0)+offset;
    fp.pd = pd;
    
    if( pd < 0 ){
      fp.itype_contact = 0;
      continue;
    }
    double eKmat[3][3];
    for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){
        eKmat[i][j] = stiff_n*n0[i]*n0[j];
      }
    }
    double eres_d[3];
    eres_d[0] = stiff_n*n0[0]*pd;
    eres_d[1] = stiff_n*n0[1]*pd;
    eres_d[2] = stiff_n*n0[2]*pd;
        
    for(unsigned int i=0;i<3;i++){
      for(unsigned int j=0;j<3;j++){
        eKmat[i][j] += -n0[i]*n0[j]*stiff_f;
      }
      eKmat[i][i] += stiff_f;      
    }
    double ap_t[3] = { co[0]-fp.aloc[0], co[1]-fp.aloc[1], co[2]-fp.aloc[2] };    
    for(unsigned int i=0;i<3;i++){ eres_d[i] += -stiff_f*ap_t[i]; }
    pmat_dd.Mearge(1,&inode,  1,&inode,     9,  &eKmat[0][0]);
    res_d.AddValue(inode,0,eres_d[0]);
    res_d.AddValue(inode,1,eres_d[1]);
    res_d.AddValue(inode,2,eres_d[2]);    
  }  
  return true;
}
    


bool AddLinSys_FrictionalContact_Penalty_NonStatic_BackwardEular
(double dt,
 Fem::Ls::CLinearSystem_Field& ls, 
 const CContactTarget3D& ct, double stiff_n, double stiff_f, double myu_s, double myu_k,
 double offset,
 unsigned int id_field_disp, 
 Fem::Field::CFieldWorld& world,
 std::vector<CFrictionPoint>& aFrictionPoint )
{
	if( !world.IsIdField(id_field_disp) ) return false;
	const Fem::Field::CField& field_disp = world.GetField(id_field_disp);
	if( field_disp.GetFieldType() != Fem::Field::VECTOR3 ) return false;
	////////////////
	
	MatVec::CMatDia_BlkCrs& pmat_dd = ls.GetMatrix(id_field_disp,  CORNER,world);
	MatVec::CVector_Blk& res_d = ls.GetResidual(id_field_disp,  CORNER,world); 
	
	const unsigned int ndim = 3;
	const CNodeAry::CNodeSeg& ns_co      = field_disp.GetNodeSeg(  CORNER,false,world,VALUE);
	const CNodeAry::CNodeSeg& ns_udisp   = field_disp.GetNodeSeg(  CORNER,true, world,VALUE);
	const CNodeAry::CNodeSeg& ns_vdisp   = field_disp.GetNodeSeg(  CORNER,true, world,VELOCITY);
	
	assert( aFrictionPoint.size() == ns_co.Size() );
	for(unsigned int inode=0;inode<ns_co.Size();inode++)
	{
		double Co[ndim];	ns_co.GetValue(   inode,Co);
		double ud[ndim];	ns_udisp.GetValue(inode,ud);
		double uv[ndim];	ns_vdisp.GetValue(inode,uv);		
		double co[3] = { Co[0]+ud[0], Co[1]+ud[1], Co[2]+ud[2] };
    CFrictionPoint& fp = aFrictionPoint[inode];
    if( fp.is_pin )
    {
      double n[3] = { fp.aloc[0]-co[0], fp.aloc[1]-co[1], fp.aloc[2]-co[2] };
      double eKmat[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
      for(unsigned int i=0;i<3;i++){ eKmat[i][i] = stiff_n; }
      double eres_d[3];
      eres_d[0] = stiff_n*n[0]*dt;
      eres_d[1] = stiff_n*n[1]*dt;
      eres_d[2] = stiff_n*n[2]*dt;

      double emat_dd[3][3];
      for(unsigned int i=0;i<9;i++){ (&emat_dd[0][0])[i] = dt*dt*(&eKmat[0][0])[i]; }
      {
        eres_d[0] -= (eKmat[0][0]*uv[0]+eKmat[0][1]*uv[1]+eKmat[0][2]*uv[2])*dt*dt;
        eres_d[1] -= (eKmat[1][0]*uv[0]+eKmat[1][1]*uv[1]+eKmat[1][2]*uv[2])*dt*dt;
        eres_d[2] -= (eKmat[2][0]*uv[0]+eKmat[2][1]*uv[1]+eKmat[2][2]*uv[2])*dt*dt;
      }
      pmat_dd.Mearge(1,&inode,  1,&inode,     9,  &emat_dd[0][0]);
      res_d.AddValue(inode,0,eres_d[0]);
      res_d.AddValue(inode,1,eres_d[1]);
      res_d.AddValue(inode,2,eres_d[2]);
      continue;
    }

    double n0[3];
    const double pd = ct.Projection(co[0],co[1],co[2], n0)+offset;
    fp.pd = pd;

    if( pd < 0 ){
      fp.itype_contact = 0;
      continue;
    }
    double eKmat[3][3];
    for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
      eKmat[i][j] = stiff_n*n0[i]*n0[j];
    }
    }
    double eres_d[3];
    eres_d[0] = stiff_n*n0[0]*pd*dt;
    eres_d[1] = stiff_n*n0[1]*pd*dt;
    eres_d[2] = stiff_n*n0[2]*pd*dt;

    // friction handling
    double ap_t[3] = { co[0]-fp.aloc[0], co[1]-fp.aloc[1], co[2]-fp.aloc[2] };
    {	// tangent vector from anchor to point
      const double t = Com::Dot3D(n0,ap_t);
      for(unsigned int i=0;i<3;i++){ ap_t[i] -= t*n0[i]; }
    }
    double velo_t[3] = { uv[0],uv[1],uv[2] };
    {	// tangent velocity
      const double t = Com::Dot3D(n0,velo_t);
      for(unsigned int i=0;i<3;i++){ velo_t[i] -= t*n0[i]; }
    }
    const double len_ap_t = Com::Length3D(ap_t);
    const double len_velo_t = Com::Length3D(velo_t);
    const double force_f = len_ap_t*stiff_f;
    const double force_n = pd*stiff_n;
    if( force_f < force_n*myu_s && len_velo_t < 1.0e-1 ){
      fp.itype_contact = 1;
      for(unsigned int i=0;i<3;i++){ eres_d[i] += -dt*stiff_f*ap_t[i]; }
      for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
          eKmat[i][j] += -n0[i]*n0[j]*stiff_f;
        }
        eKmat[i][i] += stiff_f;
      }
    }
    else{
      //			std::cout << "dynamic friction" << std::endl;
      fp.itype_contact = 2;
      if( len_velo_t > 1.0e-10 ){
        const double invlen = 1.0/len_velo_t;
        for(unsigned int i=0;i<3;i++){ velo_t[i] *= invlen; }
        for(unsigned int i=0;i<3;i++){ eres_d[i] += -dt*velo_t[i]*force_n*myu_k; }
        for(unsigned int i=0;i<3;i++){
          for(unsigned int j=0;j<3;j++){
            eKmat[i][j] += -velo_t[i]*velo_t[j]*force_n*myu_k*invlen;
          }
          eKmat[i][i] += force_n*myu_k*invlen;
        }
      }
    }

    ////////////////
    double emat_dd[3][3];
    for(unsigned int i=0;i<9;i++){ (&emat_dd[0][0])[i] = dt*dt*(&eKmat[0][0])[i]; }
    {
      eres_d[0] -= (eKmat[0][0]*uv[0]+eKmat[0][1]*uv[1]+eKmat[0][2]*uv[2])*dt*dt;
      eres_d[1] -= (eKmat[1][0]*uv[0]+eKmat[1][1]*uv[1]+eKmat[1][2]*uv[2])*dt*dt;
      eres_d[2] -= (eKmat[2][0]*uv[0]+eKmat[2][1]*uv[1]+eKmat[2][2]*uv[2])*dt*dt;
    }
    pmat_dd.Mearge(1,&inode,  1,&inode,     9,  &emat_dd[0][0]);
    res_d.AddValue(inode,0,eres_d[0]);
    res_d.AddValue(inode,1,eres_d[1]);
    res_d.AddValue(inode,2,eres_d[2]);
	}
	return true;
}


