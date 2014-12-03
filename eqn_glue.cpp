/*
 *  eqn_glue.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 7/27/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <math.h>

#include "delfem/field_world.h"
#include "delfem/field.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/linearsystem_interface_eqnsys.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/vector3d.h"

#include "eqn_glue.h"
#include "emat_glue.h"

using namespace Fem::Ls;
using namespace Fem::Field;



unsigned int MakeField_GlueEdge_Lambda
(CFieldWorld& world,
 unsigned int id_field_base, 
 unsigned int id_ea1, unsigned int id_ea2, 
 const int derivative_type )
{
	assert( world.IsIdField(id_field_base) );
	const CField& field_base = world.GetField(id_field_base);
  
	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
  
	unsigned int nno_v = world.GetEA(id_ea1).Size();
	unsigned int id_na_d = field_base.GetNodeSegInNodeAry(CORNER).id_na_va;
	unsigned int id_na_v = world.AddNodeAry(nno_v);
	unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;
	unsigned int id_ns_co = field_base.GetNodeSegInNodeAry(CORNER).id_ns_co;
	Fem::Field::CField::CNodeSegInNodeAry na_c, na_b;
	{
		na_c.id_na_va = id_na_v;
		na_c.id_na_co = 0;
		na_c.id_ns_co = 0;
		na_c.is_part_va = false;
		na_c.is_part_co = true;	
	}
	{
		na_b.id_na_co = id_na_co;
		na_b.id_ns_co = id_ns_co;
	}
  
	unsigned int id_ea_add = world.AddElemAry(nno_v,POINT);
	unsigned int id_es_c, id_es_v;
	{
		CElemAry& ea = world.GetEA(id_ea_add);
		std::vector< std::pair<unsigned int,CElemAry::CElemSeg> > aEs;
		aEs.push_back( std::make_pair(0,CElemAry::CElemSeg(id_na_d,EDGE  )) );
		aEs.push_back( std::make_pair(0,CElemAry::CElemSeg(id_na_v,CORNER)) );
		std::vector<int> lnods;
		lnods.resize(nno_v*3);
		{	
			CElemAry& ea1 = world.GetEA(id_ea1);		
			unsigned int id_es1 = 0;
			{
				std::vector<unsigned int> aIdES = ea1.GetAry_SegID();
				for(unsigned int iies=0;iies<aIdES.size();iies++){
					id_es1 = aIdES[iies];
					const CElemAry::CElemSeg& es = ea1.GetSeg(id_es1);
					if( es.GetIdNA() == id_na_co ) break;
				}
			}
			CElemAry& ea2 = world.GetEA(id_ea2);
			unsigned int id_es2 = 0;
			{
				std::vector<unsigned int> aIdES = ea2.GetAry_SegID();
				for(unsigned int iies=0;iies<aIdES.size();iies++){
					id_es2 = aIdES[iies];
					const CElemAry::CElemSeg& es = ea2.GetSeg(id_es2);
					if( es.GetIdNA() == id_na_co ) break;
				}
			}
      //			std::cout << id_es1 << " " << id_es2 << std::endl;
			const CElemAry::CElemSeg& es1 = ea1.GetSeg(id_es1);
			const CElemAry::CElemSeg& es2 = ea2.GetSeg(id_es2);
			assert( es1.Size() == es2.Size() );
			for(unsigned int ino=0;ino<nno_v;ino++){
				unsigned int no1[2]; es1.GetNodes(ino,        no1);
				unsigned int no2[2]; es2.GetNodes(nno_v-ino-1,no2);
				unsigned int no[2]  = { no1[0], no2[1] };
        //				std::cout << no[0] << " " << no[1] << std::endl;
				lnods[ino*3+0] = no[0];
				lnods[ino*3+1] = no[1];
				lnods[ino*3+2] = ino;
			}
		}
		std::vector<int> res = ea.AddSegment(aEs,lnods);
		id_es_c = res[0];
		id_es_v = res[1];
	}
	aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea_add, id_es_v,0, 0,0, 0,id_es_c) );
  
	unsigned int id_field = world.AddField(0,aElemIntp,na_c,na_b,0);
  assert( id_field != 0 );
  CField& field = world.GetField(id_field);
  field.SetValueType(field_base.GetFieldType(),field_base.GetFieldDerivativeType(),world);
  assert( field.AssertValid(world) );
	return id_field;
  
}


////////


unsigned int MakePartialField_GlueEdge_Penalty
(Fem::Field::CFieldWorld& world,
 unsigned int id_field_base, 
 unsigned int id_ea1, unsigned int id_ea2, bool is_same_dir,
 bool& is_reordering_needed )
{
//  std::cout << "MakePartialField_GluEdge_Penalty " << id_field_base << std::endl;
	assert( world.IsIdField(id_field_base) );
	const CField& field_base = world.GetField(id_field_base);
	const unsigned int id_na_co = field_base.GetNodeSegInNodeAry(CORNER).id_na_co;  
  
  assert( world.IsIdEA(id_ea1) );
  CElemAry& ea1 = world.GetEA(id_ea1);
  unsigned int id_es1 = 0;
  {
    std::vector<unsigned int> aIdES = ea1.GetAry_SegID();
    for(unsigned int iies=0;iies<aIdES.size();iies++){
      id_es1 = aIdES[iies];
      const CElemAry::CElemSeg& es = ea1.GetSeg(id_es1);
      if( es.GetIdNA() == id_na_co ) break;
    }
  }
  CElemAry& ea2 = world.GetEA(id_ea2);
  unsigned int id_es2 = 0;
  {
    std::vector<unsigned int> aIdES = ea2.GetAry_SegID();
    for(unsigned int iies=0;iies<aIdES.size();iies++){
      id_es2 = aIdES[iies];
      const CElemAry::CElemSeg& es = ea2.GetSeg(id_es2);
      if( es.GetIdNA() == id_na_co ) break;
    }
  }
  const CElemAry::CElemSeg& es1 = ea1.GetSeg(id_es1);
  const CElemAry::CElemSeg& es2 = ea2.GetSeg(id_es2);
  unsigned int ne1 = es1.Size();
  unsigned int ne2 = es2.Size();
  bool node_order_larger1;
  {
    unsigned int icnt1 = 0;
    for(unsigned int i=0;i<ne1;i++){      
        unsigned int no[2]; es1.GetNodes(i,no);
      icnt1 = (no[0]>icnt1) ? no[0]:icnt1;
      icnt1 = (no[1]>icnt1) ? no[1]:icnt1;      
    }        
    unsigned int icnt2 = 0;
    for(unsigned int i=0;i<ne2;i++){      
      unsigned int no[2]; es2.GetNodes(i,no);
      icnt2 = (no[0]>icnt2) ? no[0]:icnt2;
      icnt2 = (no[1]>icnt2) ? no[1]:icnt2;
    }            
    node_order_larger1 = icnt1 > icnt2;
    if( (ne1!=ne2) && (ne1 < ne2) == node_order_larger1 ){
      is_reordering_needed = true;
      std::cout << "this matrix is difficult" << std::endl;
    }
    else{
      is_reordering_needed = false;
    }
  }
  std::vector<int> lnods;
  std::vector<double> ratios;
  if( is_same_dir ){
    if( ne1 < ne2 ){
//      std::cout << "type1" << std::endl;
      for(unsigned int i1=0;i1<ne1+1;i1++){      
        unsigned int in1;
        if( i1 == ne1 ){  
          unsigned int no1[2]; es1.GetNodes(ne1-1,no1);
          in1 = no1[1];
        }
        else{
          unsigned int no1[2]; es1.GetNodes(i1,no1);
          in1 = no1[0];
        }
        unsigned int ie2 = floor((double)ne2/(double)ne1*i1);
        if( i1 == 0 ){ ie2 = 0; }
        if( i1 == ne1 ){ ie2 = ne2-1; }
        double r2 = (double)ne2/(double)ne1*i1-ie2;
        unsigned int no2[2]; es2.GetNodes(ie2,no2);
        if( in1 == no2[0] || in1 == no2[1] ) continue;
        //      std::cout << in1 << " " << no2[0] << " " << no2[1] << " " << r2 << std::endl;
        lnods.push_back(in1);
        lnods.push_back(no2[0]);
        lnods.push_back(no2[1]);      
        ratios.push_back(r2);      
      }          
    }
    else{      
//      std::cout << "type2" << std::endl;      
      for(unsigned int i2=0;i2<ne2+1;i2++){
        unsigned int in2;
        if( i2 == ne2 ){
          unsigned int no2[2]; es2.GetNodes(ne2-1,no2);
          in2 = no2[1];
        }
        else{
          unsigned int no2[2]; es2.GetNodes(i2,no2);
          in2 = no2[0];
        }
        unsigned int ie1 = floor((double)ne1/(double)ne2*i2);
        if( i2 == 0 ){ ie1 = 0; }
        if( i2 == ne2 ){ ie1 = ne1-1; }
        double r1 = (double)ne1/(double)ne2*i2-ie1;
        unsigned int no1[2]; es1.GetNodes(ie1,no1);
        if( in2 == no1[0] || in2 == no1[1] ) continue;
        //      std::cout << in2 << " " << no1[0] << " " << no1[1] << " " << r1 << std::endl;
        lnods.push_back(in2);
        lnods.push_back(no1[0]);
        lnods.push_back(no1[1]);      
        ratios.push_back(r1);              
      }       
    }
  }
  else{
    if( ne1 < ne2 ){
//      std::cout << "type3" << std::endl;
//      std::cout << "different dir ne1<ne2 :" << ne1 << " " << ne2 << std::endl;
      for(unsigned int i1=0;i1<ne1+1;i1++){      
        unsigned int in1;
        if( i1 == ne1 ){  
          unsigned int no1[2]; es1.GetNodes(ne1-1,no1);
          in1 = no1[1];
        }
        else{
          unsigned int no1[2]; es1.GetNodes(i1,no1);
          in1 = no1[0];
        }
        unsigned int ie2 = floor(ne2-(double)ne2/(double)ne1*i1);
        if( i1 == 0 ){ ie2 = ne2-1; }
        if( i1 == ne1 ){ ie2 = 0; }
        double r2 = 1.0-(ne2-(double)ne2/(double)ne1*i1-ie2);
        unsigned int no2[2]; es2.GetNodes(ie2,no2);
//        std::cout << "ino e r" << in1 << " " << ie2 << " " << r2 << std::endl;
        if( in1 == no2[0] || in1 == no2[1] ) continue;
//        std::cout << " no retio : " << in1 << " " << no2[1] << " " << no2[0] << " " << r2 << std::endl;
        lnods.push_back(in1);
        lnods.push_back(no2[1]);
        lnods.push_back(no2[0]);      
        ratios.push_back(r2);
      }      
    }
    else{
//      std::cout << "type4 " << id_ea1 << " " << id_ea2 << std::endl;
//      std::cout << "different dir ne1>ne2 :" << ne1 << " " << ne2 << std::endl;      
      for(unsigned int i2=0;i2<ne2+1;i2++){      
        unsigned int in2;
        if( i2 == ne2 ){  
          unsigned int no2[2]; es2.GetNodes(ne2-1,no2);
          in2 = no2[1];
        }
        else{
          unsigned int no2[2]; es2.GetNodes(i2,no2);
          in2 = no2[0];
        }
        unsigned int ie1 = floor(ne1-(double)ne1/(double)ne2*i2);
        if( i2 == 0 ){ ie1 = ne1-1; }
        if( i2 == ne2 ){ ie1 = 0; }
        double r1 = 1.0-(ne1-(double)ne1/(double)ne2*i2-ie1);
        unsigned int no1[2]; es1.GetNodes(ie1,no1);      
//        std::cout << "ino e r" << in2 << " " << ie1 << " " << r1 << std::endl;
        if( in2 == no1[0] || in2 == no1[1] ) continue;
//        std::cout << " no retio : " << in2 << " " << no1[1] << " " << no1[0] << " " << r1 << std::endl;
        lnods.push_back(in2);
        lnods.push_back(no1[1]);
        lnods.push_back(no1[0]);      
        ratios.push_back(r1);      
//        ratios.push_back(1);              
      }
    }
  }
  
  Fem::Field::CField::CNodeSegInNodeAry na_c = field_base.GetNodeSegInNodeAry(CORNER);
  Fem::Field::CField::CNodeSegInNodeAry na_b;  
  {
    assert(ratios.size()==lnods.size()/3);
    unsigned int id_na_r = world.AddNodeAry(ratios.size());
    CNodeAry& na_r = world.GetNA(id_na_r);
    CNodeAry::CNodeSeg ns_r(1,"");
    std::vector< std::pair<unsigned int,CNodeAry::CNodeSeg> > aIdSeg;
    aIdSeg.push_back( std::make_pair(0,ns_r) );
    const std::vector<int>& res = na_r.AddSegment(aIdSeg,ratios);
    const unsigned int id_ns_r = res[0];
    assert( na_r.IsSegID(id_ns_r) );
    na_b.id_na_co = id_na_r;
    na_b.id_ns_co = id_ns_r;
  }
  
	std::vector<Fem::Field::CField::CElemInterpolation> aElemIntp;
  {
    const unsigned int nelem_add = lnods.size()/3;
    unsigned int id_ea_add = world.AddElemAry(nelem_add,TRI); 
    unsigned int id_es_c;
    {
      CElemAry& ea_add = world.GetEA(id_ea_add);
      id_es_c = ea_add.AddSegment(0,CElemAry::CElemSeg(id_na_co,CORNER),lnods);
    }
    unsigned int id_es_b;
    {
      std::vector<int> lnods_b;
      for(unsigned int ino=0;ino<nelem_add;ino++){ lnods_b.push_back(ino); }
      CElemAry& ea_add = world.GetEA(id_ea_add);
      id_es_b = ea_add.AddSegment(0,CElemAry::CElemSeg(na_b.id_na_co,BUBBLE),lnods_b);
    }    
    aElemIntp.push_back( Fem::Field::CField::CElemInterpolation(id_ea_add, id_es_c,id_es_c, 0,0, 0,id_es_b) );
  }
  
	unsigned int id_field = world.AddField(id_field_base,aElemIntp,na_c,na_b,0);
  assert( id_field != 0 );  
  assert( world.IsIdField(id_field) );  
  const CField& field = world.GetField(id_field);  
  assert( field.AssertValid(world) );  
	return id_field;
}

bool AddLinSys_Glue_Lagrange_NewmarkBeta
(double dt, double gamma_newmark, double beta_newmark,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 unsigned int id_field_disp, unsigned int id_field_lambda,
 const Fem::Field::CFieldWorld& world,
 bool is_initial)
{
	const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_a = field.GetNodeSeg(CORNER,true, world,ACCELERATION);
//	MatVec::CMatDia_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
	MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);
	const CField& field_l = world.GetField(id_field_lambda);
	MatVec::CMatDia_BlkCrs& mat_ll = ls.GetMatrix(id_field_lambda,CORNER, world);
	MatVec::CMat_BlkCrs&    mat_ld = ls.GetMatrix(id_field_lambda,CORNER, id_field_disp,CORNER, world);
	MatVec::CMat_BlkCrs&    mat_dl = ls.GetMatrix(id_field_disp,CORNER, id_field_lambda,CORNER, world);
	MatVec::CVector_Blk& res_l = ls.GetResidual(id_field_lambda,CORNER,world);
	const Fem::Field::CNodeAry::CNodeSeg& ns_ul = field_l.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_vl = field_l.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_al = field_l.GetNodeSeg(CORNER,true, world,ACCELERATION);
	unsigned int id_ea = field_l.GetAryIdEA()[0];
	const Fem::Field::CElemAry::CElemSeg& es = field_l.GetElemSeg(id_ea,BUBBLE,false,world);
	for(unsigned int inol=0;inol<ns_ul.Size();inol++){
		unsigned int no[2];	es.GetNodes(inol,no);
		double C[2][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);
		double u[2][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);
		double v[2][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);
		double a[2][3];	ns_a.GetValue(no[0],a[0]);	ns_a.GetValue(no[1],a[1]);
		double ul[3]; ns_ul.GetValue(inol,ul);
		double vl[3]; ns_vl.GetValue(inol,vl);
		double al[3]; ns_al.GetValue(inol,al);
		////////////////
		double Kmat[2][2][3][3], Kmat_ll[3][3], Kmat_ld[2][3][3], Kmat_dl[2][3][3];
		double Res[2][3], Resl[3];
		GetMatRes_Glue_Lagrange_NewmarkBeta
		(Kmat,Kmat_ll,Kmat_ld,Kmat_dl, Res,Resl,
		 C, u,v,a, ul,vl,al,
		 dt,gamma_newmark,beta_newmark,
		 is_initial);
		////////////////
		mat_ll.Mearge(1,&inol, 1,&inol, 9,&Kmat_ll[0][0]);
		mat_ld.Mearge(1,&inol, 2,no,    9,&Kmat_ld[0][0][0]);
		mat_dl.Mearge(2,no,    1,&inol, 9,&Kmat_dl[0][0][0]);
		for(unsigned int ino=0;ino<2;ino++){
			res.AddValue(no[ino],0,-Res[ino][0]);
			res.AddValue(no[ino],1,-Res[ino][1]);
			res.AddValue(no[ino],2,-Res[ino][2]);
		}
		{
			res_l.AddValue(inol,0,-Resl[0]);
			res_l.AddValue(inol,1,-Resl[1]);
			res_l.AddValue(inol,2,-Resl[2]);
		}
	}
	return true;
}



void GetMatRes_Glue_Lagrange_BackwardEular
(double mat[2][2][3][3], double mat_ll[3][3],  
 double mat_ld[2][3][3], double mat_dl[2][3][3],
 double Res[2][3],  double Resl[3],
 double C[2][3], double u[2][3], double v[2][3],
 double ul[3], double vl[3], 
 double dt)
{
	const double c[2][3] = {
		{ C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
		{ C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] } };
	double eps = -0.0e-0;
	double Kmat_ll0[3][3] = { {eps,0,0}, {0,eps,0}, {0,0,eps} };
	const double k = 1.0e+0;
	double Kmat_dl0[2][3][3] = { { {+k,0,0}, {0,+k,0}, {0,0,+k} }, { {-k,0,0}, {0,-k,0}, {0,0,-k} } };
	double Kmat_ld0[2][3][3] = { { {+k,0,0}, {0,+k,0}, {0,0,+k} }, { {-k,0,0}, {0,-k,0}, {0,0,-k} } };
	const double d = 0.0;
	double Cmat_dl0[2][3][3] = { { {+d,0,0}, {0,+d,0}, {0,0,+d} }, { {-d,0,0}, {0,-d,0}, {0,0,-d} } };
	double Cmat_ld0[2][3][3] = { { {+d,0,0}, {0,+d,0}, {0,0,+d} }, { {-d,0,0}, {0,-d,0}, {0,0,-d} } };
	const double m = 0.0e-0;
	double Mmat_dl0[2][3][3] = { { {+m,0,0}, {0,+m,0}, {0,0,+m} }, { {-m,0,0}, {0,-m,0}, {0,0,-m} } };
	double Mmat_ld0[2][3][3] = { { {+m,0,0}, {0,+m,0}, {0,0,+m} }, { {-m,0,0}, {0,-m,0}, {0,0,-m} } };
	Res[0][0] = ( + ul[0]*k )*dt;
	Res[0][1] = ( + ul[1]*k )*dt;
	Res[0][2] = ( + ul[2]*k )*dt;
	Res[1][0] = -Res[0][0];					 
	Res[1][1] = -Res[0][1];
	Res[1][2] = -Res[0][2];
	Resl[0] = ( (c[0][0]-c[1][0])*k + (v[0][0]-v[1][0])*d )*dt;
	Resl[1] = ( (c[0][1]-c[1][1])*k + (v[0][1]-v[1][1])*d )*dt;
	Resl[2] = ( (c[0][2]-c[1][2])*k + (v[0][2]-v[1][2])*d )*dt;
	
	for(unsigned int i=0;i<2*2*3*3;i++){ (&mat[0][0][0][0])[i] = 0; }
	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_ld[0][0][0])[i] = dt*dt*(&Kmat_ld0[0][0][0])[i] + dt*(&Cmat_ld0[0][0][0])[i] + (&Mmat_ld0[0][0][0])[i]; }
	for(unsigned int i=0;i<2*3*3;  i++){ (&mat_dl[0][0][0])[i] = dt*dt*(&Kmat_dl0[0][0][0])[i] + dt*(&Cmat_dl0[0][0][0])[i] + (&Mmat_dl0[0][0][0])[i]; }
	for(unsigned int i=0;i<3*3;    i++){ (&mat_ll[0][0]   )[i] = dt*dt*(&Kmat_ll0[0][0]   )[i]; }
	
	{
		for(unsigned int ino=0;ino<2;ino++){
			for(unsigned int idim=0;idim<3;idim++){
				for(unsigned int jdim=0;jdim<3;jdim++){
					Res[ino][idim] += Kmat_dl0[ino][idim][jdim]*vl[jdim]*dt*dt;
				}
			}
		}
		for(unsigned int idim=0;idim<3;idim++){
			for(unsigned int jno=0;jno<2;jno++){
				for(unsigned int jdim=0;jdim<3;jdim++){
					Resl[idim] += Kmat_ld0[jno][idim][jdim]*v[jno][jdim]*dt*dt;
				}
			}
			for(unsigned int jdim=0;jdim<3;jdim++){
				Resl[idim] += Kmat_ll0[idim][jdim]*vl[jdim]*dt*dt;
			}
		}
	}
}



bool AddLinSys_Glue_Lagrange_BackwardEular
(double dt, 
 Fem::Eqn::ILinearSystem_Eqn& ls,
 unsigned int id_field_disp, unsigned int id_field_lambda,
 const Fem::Field::CFieldWorld& world,
 bool is_initial)
{	
	const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
	MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);
	const CField& field_l = world.GetField(id_field_lambda);
	MatVec::CMatDia_BlkCrs& mat_ll = ls.GetMatrix(id_field_lambda,CORNER, world);
	MatVec::CMat_BlkCrs&    mat_ld = ls.GetMatrix(id_field_lambda,CORNER, id_field_disp,CORNER, world);
	MatVec::CMat_BlkCrs&    mat_dl = ls.GetMatrix(id_field_disp,CORNER, id_field_lambda,CORNER, world);
	MatVec::CVector_Blk& res_l = ls.GetResidual(id_field_lambda,CORNER,world);
	const Fem::Field::CNodeAry::CNodeSeg& ns_ul = field_l.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_vl = field_l.GetNodeSeg(CORNER,true, world,VELOCITY);
	unsigned int id_ea = field_l.GetAryIdEA()[0];
	const Fem::Field::CElemAry::CElemSeg& es = field_l.GetElemSeg(id_ea,BUBBLE,false,world);
	for(unsigned int inol=0;inol<ns_ul.Size();inol++){
		unsigned int no[2];	es.GetNodes(inol,no);
		double C[2][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);
		double u[2][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);
		double v[2][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);
		double ul[3]; ns_ul.GetValue(inol,ul);
		double vl[3]; ns_vl.GetValue(inol,vl);
		////////////////
		double Kmat[2][2][3][3], Kmat_ll[3][3], Kmat_ld[2][3][3], Kmat_dl[2][3][3];
		double Res[2][3], Resl[3];
		GetMatRes_Glue_Lagrange_BackwardEular
		(Kmat,Kmat_ll,Kmat_ld,Kmat_dl, Res,Resl,
		 C, u,v, ul,vl,
		 dt);
		////////////////
		mat_ll.Mearge(1,&inol, 1,&inol, 9,&Kmat_ll[0][0]);
		mat_ld.Mearge(1,&inol, 2,no,    9,&Kmat_ld[0][0][0]);
		mat_dl.Mearge(2,no,    1,&inol, 9,&Kmat_dl[0][0][0]);
		for(unsigned int ino=0;ino<2;ino++){
			res.AddValue(no[ino],0,-Res[ino][0]);
			res.AddValue(no[ino],1,-Res[ino][1]);
			res.AddValue(no[ino],2,-Res[ino][2]);
		}
		{
			res_l.AddValue(inol,0,-Resl[0]);
			res_l.AddValue(inol,1,-Resl[1]);
			res_l.AddValue(inol,2,-Resl[2]);
		}
	}	
  return true;
}


bool AddLinSys_Glut_Penalty_NewmarkBeta
(double dt, double gamma_newmark, double beta_newmark,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world,
 bool is_initial)
{
	const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_a = field.GetNodeSeg(CORNER,true, world,ACCELERATION);
	MatVec::CMat_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
	MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);				
	const CField& field_dart = world.GetField(id_field_dart);
	const unsigned int id_ea = field_dart.GetAryIdEA()[0];
	const Fem::Field::CElemAry::CElemSeg& es = field_dart.GetElemSeg(id_ea,CORNER,true,world);
	const unsigned int nelem = es.Size();
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		unsigned int no[2];	es.GetNodes(ielem,no);
		double C[2][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);
		double u[2][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);
		double v[2][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);
		double a[2][3];	ns_a.GetValue(no[0],a[0]);	ns_a.GetValue(no[1],a[1]);
		double Kmat[2][2][3][3], Res[2][3];
		GetMatRes_Glue_Penalty_NewmarkBeta(Kmat,Res, C,u,v,a, stiff_dart,dt,gamma_newmark,beta_newmark,is_initial);
		mat.Mearge(2,no, 2,no, 9, &Kmat[0][0][0][0]);
		for(unsigned int ino=0;ino<2;ino++){
			res.AddValue(no[ino],0,-Res[ino][0]);
			res.AddValue(no[ino],1,-Res[ino][1]);
			res.AddValue(no[ino],2,-Res[ino][2]);
		}
	}
  return true;
}


////////////////////////////////////

bool AddLinSys_Glut_Penalty_Sensitivity
(Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world)
{
	const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
  const CField& field_dart = world.GetField(id_field_dart);
  const Fem::Field::CNodeAry::CNodeSeg& ns_r = field_dart.GetNodeSeg(BUBBLE,false,world,VALUE);
	MatVec::CMat_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
  
  double se = 0;
	const unsigned int id_ea = field_dart.GetAryIdEA()[0];
  const CElemAry& ea = world.GetEA(id_ea);
	const unsigned int nelem = ea.Size();  
	const Fem::Field::CElemAry::CElemSeg& es_c = field_dart.GetElemSeg(id_ea,CORNER,true,world);
	const Fem::Field::CElemAry::CElemSeg& es_b = field_dart.GetElemSeg(id_ea,BUBBLE,false,world);  
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		unsigned int no[3];	es_c.GetNodes(ielem,no);
		double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);  ns_c.GetValue(no[2],C[2]);
		double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);  ns_u.GetValue(no[2],u[2]);
    const double c[3][3] = {
      { C[0][0]+u[0][0], C[0][1]+u[0][1], C[0][2]+u[0][2] },
      { C[1][0]+u[1][0], C[1][1]+u[1][1], C[1][2]+u[1][2] },    
      { C[2][0]+u[2][0], C[2][1]+u[2][1], C[2][2]+u[2][2] } };      
    unsigned int ino_b; es_b.GetNodes(ielem,&ino_b);
    double r; ns_r.GetValue(ino_b,&r);
    double Kmat[3][3][3][3], Res[3][3];
    GetKmatRes_GluePointLine_Penalty(Kmat,Res, C, c, r, stiff_dart, se);
		mat.Mearge(3,no, 3,no, 9, &Kmat[0][0][0][0]);
	}    
  return true;
}

bool AddLinSys_Glut_Penalty_BackwardEular
(double dt, double damp_coeff,
 Fem::Eqn::ILinearSystem_Eqn& ls,
 double stiff_dart,
 unsigned int id_field_disp, unsigned int id_field_dart,
 const Fem::Field::CFieldWorld& world,
 double& strain_energy)
{
	const CField& field = world.GetField(id_field_disp);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = field.GetNodeSeg(CORNER,false,world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_u = field.GetNodeSeg(CORNER,true, world,VALUE);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = field.GetNodeSeg(CORNER,true, world,VELOCITY);
  const CField& field_dart = world.GetField(id_field_dart);
  const Fem::Field::CNodeAry::CNodeSeg& ns_r = field_dart.GetNodeSeg(BUBBLE,false,world,VALUE);
	MatVec::CMat_BlkCrs& mat = ls.GetMatrix(  id_field_disp,CORNER,world);
	MatVec::CVector_Blk& res = ls.GetResidual(id_field_disp,CORNER,world);				

	const unsigned int id_ea = field_dart.GetAryIdEA()[0];
  const CElemAry& ea = world.GetEA(id_ea);
	const unsigned int nelem = ea.Size();  
	const Fem::Field::CElemAry::CElemSeg& es_c = field_dart.GetElemSeg(id_ea,CORNER,true,world);
	const Fem::Field::CElemAry::CElemSeg& es_b = field_dart.GetElemSeg(id_ea,BUBBLE,false,world);  
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		unsigned int no[3];	es_c.GetNodes(ielem,no);
		double C[3][3];	ns_c.GetValue(no[0],C[0]);	ns_c.GetValue(no[1],C[1]);  ns_c.GetValue(no[2],C[2]);
		double u[3][3];	ns_u.GetValue(no[0],u[0]);	ns_u.GetValue(no[1],u[1]);  ns_u.GetValue(no[2],u[2]);
		double v[3][3];	ns_v.GetValue(no[0],v[0]);	ns_v.GetValue(no[1],v[1]);  ns_v.GetValue(no[2],v[2]);
    unsigned int ino_b; es_b.GetNodes(ielem,&ino_b);
    double r; ns_r.GetValue(ino_b,&r);
    double Kmat[3][3][3][3], Res[3][3];
	  GetMatRes_GluePointLine_Penalty_BackwardEular(Kmat,Res, C,u,v,r, stiff_dart, dt,damp_coeff,strain_energy);
		mat.Mearge(3,no, 3,no, 9, &Kmat[0][0][0][0]);
		for(unsigned int ino=0;ino<3;ino++){
			res.AddValue(no[ino],0,-Res[ino][0]);
			res.AddValue(no[ino],1,-Res[ino][1]);
			res.AddValue(no[ino],2,-Res[ino][2]);
		}
	}
  return true;
}



