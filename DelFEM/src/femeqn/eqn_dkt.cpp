#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif

#include <math.h>
#include <fstream>

#include "delfem/field_world.h"

#include "delfem/vector3d.h"
#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/diamat_blk.h"
#include "delfem/matvec/vector_blk.h"
#include "delfem/matvec/bcflag_blk.h"
#include "delfem/femls/linearsystem_field.h"
#include "delfem/femls/linearsystem_fieldsave.h"
#include "delfem/femeqn/eqn_dkt.h"
#include "delfem/femeqn/ker_emat_tri.h"
#include "delfem/femeqn/ker_emat_tet.h"
#include "delfem/femeqn/ker_emat_quad.h"
#include "delfem/femeqn/ker_emat_hex.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;
using namespace MatVec;

static void MakeCurveture(
	double B1[][3], double B2[][3][2],
	const double coord0[], const double coord1[], const double coord2[],
	const double l1, const double l2 )
{
	const double l0 = 1-l1-l2;
	const double vec0[2] = { coord2[0]-coord1[0], coord2[1]-coord1[1] };
	const double vec1[2] = { coord0[0]-coord2[0], coord0[1]-coord2[1] };
	const double vec2[2] = { coord1[0]-coord0[0], coord1[1]-coord0[1] };
	const double invsqlen0 = 1.0/(vec0[0]*vec0[0]+vec0[1]*vec0[1]);
	const double invsqlen1 = 1.0/(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
	const double invsqlen2 = 1.0/(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
	double p0=-6*vec0[0]*invsqlen0, q0=3*vec0[0]*vec0[1]*invsqlen0, r0=3*vec0[1]*vec0[1]*invsqlen0, t0=-6*vec0[1]*invsqlen0;
	double p1=-6*vec1[0]*invsqlen1, q1=3*vec1[0]*vec1[1]*invsqlen1, r1=3*vec1[1]*vec1[1]*invsqlen1, t1=-6*vec1[1]*invsqlen1;
	double p2=-6*vec2[0]*invsqlen2, q2=3*vec2[0]*vec2[1]*invsqlen2, r2=3*vec2[1]*vec2[1]*invsqlen2, t2=-6*vec2[1]*invsqlen2;

	double H1[4][3];
	H1[0][0]=-(l1-l0)*t2+l2*t1;	H1[0][1]=+(l1-l0)*t2+l2*t0;	H1[0][2]=-l2*(t0+t1);
	H1[1][0]= (l2-l0)*t1-l1*t2;	H1[1][1]=+l1*(t0+t2);		H1[1][2]=-l1*t0-(l2-l0)*t1;
	H1[2][0]= (l1-l0)*p2-l2*p1;	H1[2][1]=-(l1-l0)*p2-l2*p0;	H1[2][2]=+l2*(p0+p1);
	H1[3][0]=-(l2-l0)*p1+l1*p2;	H1[3][1]=-l1*(p0+p2);		H1[3][2]=+l1*p0+(l2-l0)*p1;

	double H2[4][3][2];
	H2[0][0][0]=-1+(l1-l0)*r2+l2*r1;		H2[0][0][1]=-(l1-l0)*q2-l2*q1; 
	H2[0][1][0]= 1+(l1-l0)*r2-l2*r0;		H2[0][1][1]=-(l1-l0)*q2+l2*q0;
	H2[0][2][0]=-l2*(r0-r1);				H2[0][2][1]= l2*(q0-q1);

	H2[1][0][0]=-1+l1*r2+(l2-l0)*r1;		H2[1][0][1]=-l1*q2-(l2-l0)*q1;
	H2[1][1][0]=-l1*(r0-r2);				H2[1][1][1]= l1*(q0-q2);
	H2[1][2][0]= 1-l1*r0+(l2-l0)*r1;		H2[1][2][1]= l1*q0-(l2-l0)*q1;

	H2[2][0][0]=-(l1-l0)*q2-l2*q1;			H2[2][0][1]= 2-6*l0-(l1-l0)*r2-l2*r1;	
	H2[2][1][0]=-(l1-l0)*q2+l2*q0;			H2[2][1][1]=-2+6*l1-(l1-l0)*r2+l2*r0;
	H2[2][2][0]= l2*(q0-q1);				H2[2][2][1]= l2*(r0-r1);

	H2[3][0][0]=-l1*q2-(l2-l0)*q1;			H2[3][0][1]= 2-6*l0-l1*r2-(l2-l0)*r1;
	H2[3][1][0]= l1*(q0-q2);				H2[3][1][1]= l1*(r0-r2);
	H2[3][2][0]= l1*q0-(l2-l0)*q1;			H2[3][2][1]=-2+6*l2+l1*r0-(l2-l0)*r1;

	double dldx[3][2];
	double const_term[3];
	TriDlDx(dldx,const_term,coord0,coord1,coord2);

	for(unsigned int i=0;i<3;i++){
		B1[0][i] =  dldx[1][0]*H1[2][i]+dldx[2][0]*H1[3][i];
		B1[1][i] = -dldx[1][1]*H1[0][i]-dldx[2][1]*H1[1][i];
		B1[2][i] =  dldx[1][1]*H1[2][i]+dldx[2][1]*H1[3][i] - dldx[1][0]*H1[0][i]-dldx[2][0]*H1[1][i];
	}
	for(unsigned int i=0;i<3;i++){
		B2[0][i][0] =  dldx[1][0]*H2[2][i][0]+dldx[2][0]*H2[3][i][0];
		B2[0][i][1] =  dldx[1][0]*H2[2][i][1]+dldx[2][0]*H2[3][i][1];
		B2[1][i][0] = -dldx[1][1]*H2[0][i][0]-dldx[2][1]*H2[1][i][0];
		B2[1][i][1] = -dldx[1][1]*H2[0][i][1]-dldx[2][1]*H2[1][i][1];
		B2[2][i][0] =  dldx[1][1]*H2[2][i][0]+dldx[2][1]*H2[3][i][0] - dldx[1][0]*H2[0][i][0]-dldx[2][0]*H2[1][i][0];
		B2[2][i][1] =  dldx[1][1]*H2[2][i][1]+dldx[2][1]*H2[3][i][1] - dldx[1][0]*H2[0][i][1]-dldx[2][0]*H2[1][i][1];
	}
}

static void MakeStiffMat_Membrane(double eKmat[][3][2][2], double eres[][2],
		const double young, const double poisson, const double thickness,
		double coords[][2], double disp[][2] )
{
	const unsigned int ndim = 2;
	const unsigned int nno = 3;

	const double lambda = young*poisson/( (1+poisson)*(1-2*poisson) );
	const double myu = young/( 1-poisson*poisson );
	
	// –ÊÏ‚ð‹‚ß‚é
	double area = TriArea(coords[0],coords[1],coords[2]);

	// Œ`óŠÖ”‚Ì‚˜‚™”÷•ª‚ð‹‚ß‚é
	double dldx[nno][ndim];		// Œ`óŠÖ”‚Ì‹óŠÔ”÷•ª
	double zero_order_term[nno];	// Œ`óŠÖ”‚Ì’è”€
	TriDlDx(dldx, zero_order_term,   coords[0], coords[1], coords[2]);

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
	   double dtmp1 = 0.0;
	   for(unsigned int idim=0;idim<ndim;idim++){
		  for(unsigned int jdim=0;jdim<ndim;jdim++){
			 eKmat[ino][jno][idim][jdim] = area*thickness*( lambda*dldx[ino][idim]*dldx[jno][jdim]
				                                              +myu*dldx[jno][idim]*dldx[ino][jdim] );
		  }
		  dtmp1 += dldx[ino][idim]*dldx[jno][idim];
	   }
	   for(unsigned int idim=0;idim<ndim;idim++){
		  eKmat[ino][jno][idim][idim] += area*thickness*myu*dtmp1;
	   }
	}
	}

	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int idim=0;idim<2;idim++){
		eres[ino][idim] = 0.0;
		for(unsigned int jno=0;jno<nno;jno++){
			eres[ino][idim] -= eKmat[ino][jno][idim][0]*disp[jno][0]
				                   + eKmat[ino][jno][idim][1]*disp[jno][1];
		}
	}
	}
}

static void MakeStiffMat_DKT_PlateBending(
		double emat_ww[][3], double emat_wr[][3][2], double emat_rw[][3][2], double emat_rr[][3][2][2], 
		double eres_w[], double eres_r[3][2],
		const double young, const double poisson, const double thickness, 
		const double coord[][2], const double w[], const double rot[][2])
{	
	const unsigned int ndim = 2;
	const unsigned int nno = 3;

	for(unsigned int i=0;i<nno*nno;    i++){ *(&emat_ww[0][0]      +i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_wr[0][0][0]   +i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*2;  i++){ *(&emat_rw[0][0][0]   +i) = 0.0; }
	for(unsigned int i=0;i<nno*nno*2*2;i++){ *(&emat_rr[0][0][0][0]+i) = 0.0; }

	double dmat[3][3];
	{
		const double dtmp1 = young*thickness*thickness*thickness/(12.0*(1.0-poisson*poisson));
		dmat[0][0] = dtmp1;			dmat[0][1] = dtmp1*poisson;	dmat[0][2] = 0;
		dmat[1][0] = dtmp1*poisson;	dmat[1][1] = dtmp1;			dmat[1][2] = 0;
		dmat[2][0] = 0.0;			dmat[2][1] = 0.0;			dmat[2][2] = dtmp1*(1-poisson)*0.5;
	}
	double B1[3][nno];
	double B2[3][nno][ndim];
	const double area = TriArea(coord[0],coord[1],coord[2]);
	const double dtmp1 = area/3.0;
	for(unsigned int iw=0;iw<3;iw++){
		if(      iw == 0 ){ MakeCurveture(B1,B2,coord[0],coord[1],coord[2],0.5,0.5); }
		else if( iw == 1 ){ MakeCurveture(B1,B2,coord[0],coord[1],coord[2],0.0,0.5); }
		else if( iw == 2 ){ MakeCurveture(B1,B2,coord[0],coord[1],coord[2],0.5,0.0); }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int k=0;k<3;k++){
		for(unsigned int l=0;l<3;l++){
			emat_ww[ino][jno] += dtmp1*B1[k][ino]*dmat[k][l]*B1[l][jno];
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_rr[ino][jno][idim][jdim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B2[l][jno][jdim];
			}
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				emat_rw[ino][jno][idim] += dtmp1*B2[k][ino][idim]*dmat[k][l]*B1[l][jno];
				emat_wr[ino][jno][idim] += dtmp1*B1[k][ino]*dmat[k][l]*B2[l][jno][idim];
			}
		}
		}
		}
		}
	}
	
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int idim=0;idim<2;idim++){
		eres_r[ino][idim] = 0.0;
		for(unsigned int jno=0;jno<nno;jno++){
			eres_r[ino][idim] -= emat_rr[ino][jno][idim][0]*rot[jno][0]
				               + emat_rr[ino][jno][idim][1]*rot[jno][1]
							   + emat_rw[ino][jno][idim]*w[jno];
		}
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
		eres_w[ino] = 0.0;
		for(unsigned int jno=0;jno<nno;jno++){
			eres_w[ino] -= emat_ww[ino][jno]*w[jno]
				         + emat_wr[ino][jno][0]*rot[jno][0]
						 + emat_wr[ino][jno][1]*rot[jno][1];
		}
	}
}


static inline void MakeRotation(
	const double coord0[], const double coord1[], const double coord2[],
	const double l1, const double l2,
	const double w[], const double rot[][2],
	double rot_ins[] )
{
	const double l0 = 1-l1-l2;
	const double vec0[2] = { coord2[0]-coord1[0], coord2[1]-coord1[1] };
	const double vec1[2] = { coord0[0]-coord2[0], coord0[1]-coord2[1] };
	const double vec2[2] = { coord1[0]-coord0[0], coord1[1]-coord0[1] };
	const double invsqlen0 = 1.0/(vec0[0]*vec0[0]+vec0[1]*vec0[1]);
	const double invsqlen1 = 1.0/(vec1[0]*vec1[0]+vec1[1]*vec1[1]);
	const double invsqlen2 = 1.0/(vec2[0]*vec2[0]+vec2[1]*vec2[1]);

	double rot_e[3][2];
	{
		double a0 = -vec0[0]*invsqlen0;
		double a1 = -vec1[0]*invsqlen1;
		double a2 = -vec2[0]*invsqlen2;
		double b0 =  0.75*vec0[0]*vec0[1]*invsqlen0;
		double b1 =  0.75*vec1[0]*vec1[1]*invsqlen1;
		double b2 =  0.75*vec2[0]*vec2[1]*invsqlen2;
		double c0 =  0.25-0.75*vec0[1]*vec0[1]*invsqlen0;
		double c1 =  0.25-0.75*vec1[1]*vec1[1]*invsqlen1;
		double c2 =  0.25-0.75*vec2[1]*vec2[1]*invsqlen2;
		double d0 = -vec0[1]*invsqlen0;
		double d1 = -vec1[1]*invsqlen1;
		double d2 = -vec2[1]*invsqlen2;
		double e0 =  0.75*vec0[1]*vec0[1]*invsqlen0-0.5;
		double e1 =  0.75*vec1[1]*vec1[1]*invsqlen1-0.5;
		double e2 =  0.75*vec2[1]*vec2[1]*invsqlen2-0.5;
		rot_e[0][0] =  1.5*d0*w[1]-1.5*d0*w[2]-e0*(rot[1][0]+rot[2][0])+b0*(rot[1][1]+rot[2][1]);
		rot_e[1][0] =  1.5*d1*w[2]-1.5*d1*w[0]-e1*(rot[2][0]+rot[0][0])+b1*(rot[2][1]+rot[0][1]);
		rot_e[2][0] =  1.5*d2*w[0]-1.5*d2*w[1]-e2*(rot[0][0]+rot[1][0])+b2*(rot[0][1]+rot[1][1]);
		rot_e[0][1] = -1.5*a0*w[1]+1.5*a0*w[2]+b0*(rot[1][0]+rot[2][0])-c0*(rot[1][1]+rot[2][1]);
		rot_e[1][1] = -1.5*a1*w[2]+1.5*a1*w[0]+b1*(rot[2][0]+rot[0][0])-c1*(rot[2][1]+rot[0][1]);
		rot_e[2][1] = -1.5*a2*w[0]+1.5*a2*w[1]+b2*(rot[0][0]+rot[1][0])-c2*(rot[0][1]+rot[1][1]);
	}

	rot_ins[0] = l0*(2*l0-1)*rot[0][0] + l1*(2*l1-1)*rot[1][0] + l2*(2*l2-1)*rot[2][0]
		+ 4*l1*l2*rot_e[0][0] + 4*l2*l0*rot_e[1][0] + 4*l0*l1*rot_e[2][0];
	rot_ins[1] = l0*(2*l0-1)*rot[0][1] + l1*(2*l1-1)*rot[1][1] + l2*(2*l2-1)*rot[2][1]
		+ 4*l1*l2*rot_e[0][1] + 4*l2*l0*rot_e[1][1] + 4*l0*l1*rot_e[2][1];
}

static bool AddLinearSystem_DKT2D_P1(
		CLinearSystem_Field& ls, 
		const unsigned int id_field_deflect, const unsigned int id_field_rot,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT2D Triangle 3point" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);

	if( !world.IsIdField(id_field_rot) ) return false;
	const CField& field_rot = world.GetField(id_field_rot);

	const CElemAry::CElemSeg& es_c_va = field_deflect.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_deflect.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 2;

	unsigned int no[nno];

	double w[nno];
	double rot[nno][ndim];
	double coord[nno][ndim];
				
	double emat_ww[nno][nno];
	double emat_rr[nno][nno][ndim][ndim];
	double emat_wr[nno][nno][ndim];
	double emat_rw[nno][nno][ndim];
				
	CMatDia_BlkCrs& mat_ww = ls.GetMatrix(id_field_deflect,CORNER,world);
	CMatDia_BlkCrs& mat_rr = ls.GetMatrix(id_field_rot,    CORNER,world);
	CMat_BlkCrs& mat_wr = ls.GetMatrix(id_field_deflect,CORNER, id_field_rot,    CORNER, world);
	CMat_BlkCrs& mat_rw = ls.GetMatrix(id_field_rot,    CORNER, id_field_deflect,CORNER, world);
	CVector_Blk& res_w = ls.GetResidual(id_field_deflect,CORNER,world);
	CVector_Blk& res_r = ls.GetResidual(id_field_rot,    CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_w = field_deflect.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_deflect.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++){
		es_c_co.GetNodes(ielem,no);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_co.GetValue(no[inoes],coord[inoes]);
		}
		es_c_va.GetNodes(ielem,no);
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_w.GetValue(no[inoes],&w[inoes]);
		}
		for(unsigned int inoes=0;inoes<nno;inoes++){
			ns_c_r.GetValue(no[inoes],rot[inoes]);
		}
		////////////////
		double eres_r[nno][2], eres_w[nno];
		MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr, 
			eres_w,eres_r,  
			1.0,0.0,1,   
			coord,w,rot);
		////////////////
		for(unsigned int i=0;i<nno;              i++){ *(&eres_w[0]+i)           = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				eres_w[ino] -= emat_ww[ino][jno]*w[jno];
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres_w[ino] -= emat_wr[ino][jno][jdim]*rot[jno][jdim];
				}
			}
		}
		for(unsigned int i=0;i<nno*ndim;         i++){ *(&eres_r[0][0]+i)        = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jno=0;jno<nno;jno++){
				eres_r[ino][idim] -= emat_rw[ino][jno][idim]*w[jno];
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eres_r[ino][idim] -= emat_rr[ino][jno][idim][jdim]*rot[jno][jdim];
				}
			}
		}
		}
		mat_ww.Mearge(nno,no,nno,no,1,&emat_ww[0][0]);
		mat_rr.Mearge(nno,no,nno,no,ndim*ndim,&emat_rr[0][0][0][0]);
		mat_rw.Mearge(nno,no,nno,no,ndim,&emat_rw[0][0][0]);
		mat_wr.Mearge(nno,no,nno,no,ndim,&emat_wr[0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_w.AddValue(no[ino],0,eres_w[ino]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_r.AddValue(no[ino],0,eres_r[ino][0]);
			res_r.AddValue(no[ino],1,eres_r[ino][1]);
		}
	}
	return true;
}


bool Fem::Eqn::AddLinearSystem_DKT2D_Static(
	Fem::Ls::CLinearSystem_Field& ls,
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != SCALAR ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
			return AddLinearSystem_DKT2D_P1(ls,id_field_deflect,id_field_rot,world,id_ea);
		}
        assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT2D_Static(
					ls,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}

static void MakeLocalCoordBase(double loc_base[][3], double coord2[][2], const double coord[][3])
{
	Com::CVector3D vec01(coord[1][0]-coord[0][0],coord[1][1]-coord[0][1],coord[1][2]-coord[0][2]);
	Com::CVector3D vec02(coord[2][0]-coord[0][0],coord[2][1]-coord[0][1],coord[2][2]-coord[0][2]);
	Com::CVector3D e0  = vec01 / vec01.Length();
	Com::CVector3D tmp = vec02 / vec02.Length();
	Com::CVector3D e2 = Com::Cross(e0,tmp);
	e2.Normalize();
	Com::CVector3D e1 = Com::Cross(e2,e0);
	loc_base[0][0] = e0.x; loc_base[0][1] = e0.y; loc_base[0][2] = e0.z;
	loc_base[1][0] = e1.x; loc_base[1][1] = e1.y; loc_base[1][2] = e1.z;
	loc_base[2][0] = e2.x; loc_base[2][1] = e2.y; loc_base[2][2] = e2.z;

	coord2[0][0] = 0; 
	coord2[0][1] = 0;
	coord2[1][0] = loc_base[0][0]*(coord[1][0]-coord[0][0]) 
	 			 + loc_base[0][1]*(coord[1][1]-coord[0][1]) 
				 + loc_base[0][2]*(coord[1][2]-coord[0][2]);
	coord2[1][1] = 0;
	coord2[2][0] = loc_base[0][0]*(coord[2][0]-coord[0][0]) 
			     + loc_base[0][1]*(coord[2][1]-coord[0][1]) 
			     + loc_base[0][2]*(coord[2][2]-coord[0][2]);
	coord2[2][1] = loc_base[1][0]*(coord[2][0]-coord[0][0]) 
		         + loc_base[1][1]*(coord[2][1]-coord[0][1]) 
			     + loc_base[1][2]*(coord[2][2]-coord[0][2]);
}


static bool AddLinearSystem_DKT3D_Linear_P1(
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT 3dim Linear" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

    CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,world);
	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,world);
    CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp, CORNER, id_field_theta,CORNER, world);
	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER, id_field_disp, CORNER, world);
    CVector_Blk& res_d = ls.GetResidual(id_field_disp, CORNER,world);
	CVector_Blk& res_t = ls.GetResidual(id_field_theta,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		////////////////

		// local coordinate base ( loc_coord_base ) ‚ðì‚é
		double loc_base[nno][ndim], coord2[nno][2];
		MakeLocalCoordBase(loc_base,coord2,coord);
		const double area = TriArea(coord2[0],coord2[1],coord2[2]);

		double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
		{
			double tmp_w[3], tmp_r[3][2];
			double w[3] = { 0,0,0 };
			double r[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
				tmp_w, tmp_r,
				young,poisson,thickness,
				coord2,w,r);
		}
		const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

		double emat_uu[3][3][2][2];
		{
			double eres_u[3][2];
			double disp[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_dd[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
					+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
				emat_dt[ino][jno][idim][jdim] 
					= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
					+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
				emat_td[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
					+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
				emat_tt[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
			}
			}	
			}
			for(unsigned int idim=0;idim<ndim;idim++){
					emat_tt[ino][ino][idim][idim] += loc_base[2][idim]*loc_base[2][idim]*torsion_stiff;
			}
		}

		double eres_d[nno][ndim];
		double eres_t[nno][ndim];
		////////////////		
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmp1*g_x;
				eres_d[ino][1] = dtmp1*g_y;
				eres_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmpx;
				eres_d[ino][1] = dtmpy;
				eres_d[ino][2] = dtmpz;
			}
		}

		for(unsigned int i=0;i<nno*ndim; i++){ *(&eres_t[0][0]+i) = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_d[ino][idim] -= emat_dd[ino][jno][idim][jdim]*disp[jno][jdim] 
					               + emat_dt[ino][jno][idim][jdim]*rot[ jno][jdim];
				eres_t[ino][idim] -= emat_td[ino][jno][idim][jdim]*disp[jno][jdim] 
					               + emat_tt[ino][jno][idim][jdim]*rot[ jno][jdim];
			}
			}
		}
		}

		mat_dd.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_t.AddValue(no[ino],0,eres_t[ino][0]);
			res_t.AddValue(no[ino],1,eres_t[ino][1]);
			res_t.AddValue(no[ino],2,eres_t[ino][2]);
		}
	}
	return true;
}


static bool AddLinearSystem_DKT3D_Linear_P1_Combined(
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D Linear (Comblined)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(  id_field_disp,CORNER,world);
	CVector_Blk&    res_d  = ls.GetResidual(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		////////////////

		// local coordinate base ( loc_coord_base ) ‚ðì‚é
		double loc_base[nno][ndim], coord2[nno][2];
		MakeLocalCoordBase(loc_base,coord2,coord);
		const double area = TriArea(coord2[0],coord2[1],coord2[2]);

		double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
		{
			double tmp_w[3], tmp_r[3][2];
			double w[3] = { 0,0,0 };
			double r[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
				tmp_w, tmp_r,
				young,poisson,thickness,
				coord2,w,r);
		}
		const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

		double emat_uu[3][3][2][2];
		{
			double eres_u[3][2];
			double disp[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_dd[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
					+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
				emat_dt[ino][jno][idim][jdim] 
					= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
					+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
				emat_td[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
					+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
				emat_tt[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
			}
			}	
			}
			for(unsigned int idim=0;idim<ndim;idim++){
					emat_tt[ino][ino][idim][idim] += loc_base[2][idim]*loc_base[2][idim]*torsion_stiff;
			}
		}

		double eres_d[nno][ndim];
		double eres_t[nno][ndim];
		////////////////		
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmp1*g_x;
				eres_d[ino][1] = dtmp1*g_y;
				eres_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmpx;
				eres_d[ino][1] = dtmpy;
				eres_d[ino][2] = dtmpz;
			}
		}

		for(unsigned int i=0;i<nno*ndim; i++){ *(&eres_t[0][0]+i) = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				eres_d[ino][idim] -= emat_dd[ino][jno][idim][jdim]*disp[jno][jdim] 
					               + emat_dt[ino][jno][idim][jdim]*rot[ jno][jdim];
				eres_t[ino][idim] -= emat_td[ino][jno][idim][jdim]*disp[jno][jdim] 
					               + emat_tt[ino][jno][idim][jdim]*rot[ jno][jdim];
			}
			}
		}
		}

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = emat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = emat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = emat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = emat_tt[ino][jno][idim][jdim];
		}
		}
		}
		}
		mat_dd.Mearge(nno,no,nno,no,6*6,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
			res_d.AddValue(no[ino],3,eres_t[ino][0]);
			res_d.AddValue(no[ino],4,eres_t[ino][1]);
			res_d.AddValue(no[ino],5,eres_t[ino][2]);
		}
	}
	return true;
}




bool Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static(
	Fem::Ls::CLinearSystem_Field& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press, 
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != VECTOR3 ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if(    ls.FindIndexArray_Seg(id_field_deflect,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_rot,    CORNER,world) ){
			    return AddLinearSystem_DKT3D_Linear_P1_Combined(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_Linear_P1(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
		}
		else{
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static(
					ls,
					young, poisson, thickness, arearho,
					g_x, g_y, g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

static bool AddLinearSystem_DKT3D_Linear_P1_Save(
		CLinearSystem_Save& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT 3dim Linear Save" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,world);
	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,world);
	CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp, CORNER,  id_field_theta,CORNER,  world);
	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER,  id_field_disp, CORNER,  world);
	CMat_BlkCrs& mat_dd_bound = ls.GetMatrix_Boundary(id_field_disp, CORNER,  id_field_disp, CORNER,  world);
	CMat_BlkCrs& mat_tt_bound = ls.GetMatrix_Boundary(id_field_theta,CORNER,  id_field_theta,CORNER,  world);
	CMat_BlkCrs& mat_dt_bound = ls.GetMatrix_Boundary(id_field_disp, CORNER,  id_field_theta,CORNER,  world);
	CMat_BlkCrs& mat_td_bound = ls.GetMatrix_Boundary(id_field_theta,CORNER,  id_field_disp, CORNER,  world);
	CVector_Blk& force_d = ls.GetForce(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		////////////////

		// local coordinate base ( loc_coord_base ) ‚ðì‚é
		double loc_base[nno][ndim], coord2[nno][2];
		MakeLocalCoordBase(loc_base,coord2,coord);
		const double area = TriArea(coord2[0],coord2[1],coord2[2]);

		double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
		{
			double tmp_w[3], tmp_r[3][2];
			double w[3] = { 0,0,0 };
			double r[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
				tmp_w, tmp_r,
				young,poisson,thickness,
				coord2,w,r);
		}
		const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

		double emat_uu[3][3][2][2];
		{
			double eres_u[3][2];
			double disp[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_dd[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
					+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
				emat_dt[ino][jno][idim][jdim] 
					= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
					+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
				emat_td[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
					+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
				emat_tt[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
			}
			}	
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				emat_tt[ino][ino][idim][idim] += loc_base[2][idim]*loc_base[2][idim]*torsion_stiff;
			}
		}

		double eforce_d[nno][ndim];
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmp1*g_x;
				eforce_d[ino][1] = dtmp1*g_y;
				eforce_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmpx;
				eforce_d[ino][1] = dtmpy;
				eforce_d[ino][2] = dtmpz;
			}
		}

		mat_dd.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);

        mat_dd_bound.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt_bound.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td_bound.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt_bound.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			force_d.AddValue(no[ino],0,eforce_d[ino][0]);
			force_d.AddValue(no[ino],1,eforce_d[ino][1]);
			force_d.AddValue(no[ino],2,eforce_d[ino][2]);
		}
	}
	return true;
}




static bool AddLinearSystem_DKT3D_Linear_P1_Save_Combined(
		CLinearSystem_Save& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT 3dim Linear Save (Combined)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert( ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd    = ls.GetMatrix(         id_field_disp,CORNER,world);
	CMat_BlkCrs& mat_dd_bound = ls.GetMatrix_Boundary(id_field_disp,CORNER,id_field_disp,CORNER, world);
	CVector_Blk& force_d      = ls.GetForce(          id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		////////////////

		// local coordinate base ( loc_coord_base ) ‚ðì‚é
		double loc_base[nno][ndim], coord2[nno][2];
		MakeLocalCoordBase(loc_base,coord2,coord);
		const double area = TriArea(coord2[0],coord2[1],coord2[2]);

		double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
		{
			double tmp_w[3], tmp_r[3][2];
			double w[3] = { 0,0,0 };
			double r[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
				tmp_w, tmp_r,
				young,poisson,thickness,
				coord2,w,r);
		}
		const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

		double emat_uu[3][3][2][2];
		{
			double eres_u[3][2];
			double disp[3][2] = { {0,0}, {0,0}, {0,0} };
			MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_dd[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
					+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
				emat_dt[ino][jno][idim][jdim] 
					= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
					+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
				emat_td[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
					+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
				emat_tt[ino][jno][idim][jdim] 
					= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
					+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
					+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
					+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
			}
			}	
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				emat_tt[ino][ino][idim][idim] += loc_base[2][idim]*loc_base[2][idim]*torsion_stiff;
			}
		}

		double eforce_d[nno][ndim];
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmp1*g_x;
				eforce_d[ino][1] = dtmp1*g_y;
				eforce_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmpx;
				eforce_d[ino][1] = dtmpy;
				eforce_d[ino][2] = dtmpz;
			}
		}

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = emat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = emat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = emat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = emat_tt[ino][jno][idim][jdim];
		}
		}
		}
		}
		mat_dd.Mearge(      nno,no, nno,no, 6*6, &emat[0][0][0][0]);
		mat_dd_bound.Mearge(nno,no, nno,no, 6*6, &emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			force_d.AddValue(no[ino],0,eforce_d[ino][0]);
			force_d.AddValue(no[ino],1,eforce_d[ino][1]);
			force_d.AddValue(no[ino],2,eforce_d[ino][2]);
		}
	}
	return true;
}




bool Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static_Save(
	Fem::Ls::CLinearSystem_Save& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press, 
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != VECTOR3 ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_deflect,CORNER,world)
                == ls.FindIndexArray_Seg(id_field_rot,CORNER,world) )
            {
			    return AddLinearSystem_DKT3D_Linear_P1_Save_Combined(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_Linear_P1_Save(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
		}
		assert(0);
		return false;
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_Linear_Static_Save(
					ls,
					young, poisson, thickness, arearho,
					g_x, g_y, g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


static bool AddLinearSystem_DKT3D_LinearNonStatic_P1(
		double dt, double gamma_newmark, double beta_newmark,
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D LinearNonstatic" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

/*	const bool isnt_combined = (
		     ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
          != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );*/

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,  world);
	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,  world);
	CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp, CORNER, id_field_theta,CORNER,  world);
	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER, id_field_disp, CORNER,  world);
	CVector_Blk& res_d = ls.GetResidual(id_field_disp, CORNER,  world); 
	CVector_Blk& res_t = ls.GetResidual(id_field_theta,CORNER,  world);

	const CNodeAry::CNodeSeg& ns_c_d  = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r  = field_rot.GetNodeSeg( CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double disp_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_da.GetValue( no[ino],disp_acc[ ino]); }
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_dv.GetValue( no[ino],disp_velo[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		double rot_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_rv.GetValue( no[ino],rot_velo[ ino]); }
		double rot_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_ra.GetValue( no[ino],rot_acc[ ino]); }
		////////////////

		double eKmat_dd[nno][nno][ndim][ndim];
		double eKmat_tt[nno][nno][ndim][ndim];
		double eKmat_dt[nno][nno][ndim][ndim];
		double eKmat_td[nno][nno][ndim][ndim];
		double area;
		double loc_base[nno][ndim];
		{
			double coord2[nno][2];
			// local coordinate base ( loc_coord_base ) ‚ðì‚é
			MakeLocalCoordBase(loc_base,coord2,coord);
			area = TriArea(coord2[0],coord2[1],coord2[2]);

			double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
			{
				double tmp_w[3], tmp_r[3][2];
				double w[3] = { 0,0,0 };
				double r[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
					tmp_w, tmp_r,
					young,poisson,thickness,
					coord2,w,r);
			}
			const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

			double emat_uu[3][3][2][2];
			{
				double eres_u[3][2];
				double disp[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
			}

			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eKmat_dd[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
						+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
					eKmat_dt[ino][jno][idim][jdim] 
						= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
						+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
					eKmat_td[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
						+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
					eKmat_tt[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
				}
				}	
				}
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
						eKmat_tt[ino][ino][idim][jdim] += loc_base[2][idim]*loc_base[2][jdim]*torsion_stiff;
				}
				}
			}
		}
		double eMmat_dd[nno][nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area/12.0;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat_dd[ino][jno][0][0] += dtmp0;
					eMmat_dd[ino][jno][1][1] += dtmp0;
					eMmat_dd[ino][jno][2][2] += dtmp0;
				}
				eMmat_dd[ino][ino][0][0] += dtmp0;
				eMmat_dd[ino][ino][1][1] += dtmp0;
				eMmat_dd[ino][ino][2][2] += dtmp0;
			}
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		double eres_d[nno][ndim];
		double eres_t[nno][ndim];
		{
			double dtmp1 = beta_newmark*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_dd[0][0][0][0])[i] = (&eMmat_dd[0][0][0][0])[i]+dtmp1*(&eKmat_dd[0][0][0][0])[i];
				(&emat_dt[0][0][0][0])[i] = dtmp1*(&eKmat_dt[0][0][0][0])[i];
				(&emat_td[0][0][0][0])[i] = dtmp1*(&eKmat_td[0][0][0][0])[i];
				(&emat_tt[0][0][0][0])[i] = dtmp1*(&eKmat_tt[0][0][0][0])[i];
			}
		}

		////////////////		
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmp1*g_x;
				eres_d[ino][1] = dtmp1*g_y;
				eres_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmpx;
				eres_d[ino][1] = dtmpy;
				eres_d[ino][2] = dtmpz;
			}
		}

		for(unsigned int i=0;i<nno*ndim; i++){ *(&eres_t[0][0]+i) = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= eKmat_dd[ino][jno][idim][jdim]*disp[jno][jdim] 
				               + eKmat_dt[ino][jno][idim][jdim]*rot[ jno][jdim];
			eres_d[ino][idim] -= eMmat_dd[ino][jno][idim][jdim]*disp_acc[jno][jdim]; 
			eres_t[ino][idim] -= eKmat_td[ino][jno][idim][jdim]*disp[jno][jdim] 
				               + eKmat_tt[ino][jno][idim][jdim]*rot[ jno][jdim];	
		}
		}
		}
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= dt*(       eKmat_dd[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_dt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_d[ino][idim] -= 0.5*dt*dt*(eKmat_dd[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_dt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
			////////////////////////////////
			eres_t[ino][idim] -= dt*(       eKmat_td[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_tt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_t[ino][idim] -= 0.5*dt*dt*(eKmat_td[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_tt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
				                         
		}
		}
		}
		}

		mat_dd.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_t.AddValue(no[ino],0,eres_t[ino][0]);
			res_t.AddValue(no[ino],1,eres_t[ino][1]);
			res_t.AddValue(no[ino],2,eres_t[ino][2]);
		}
	}
	return true;
}





static bool AddLinearSystem_DKT3D_LinearNonStatic_P1_Combined(
		double dt, double gamma_newmark, double beta_newmark,
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D LinearNonstatic (Combined)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert( ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		 == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(  id_field_disp,CORNER,world);
	CVector_Blk&    res_d  = ls.GetResidual(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d  = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r  = field_rot.GetNodeSeg( CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double disp_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_da.GetValue( no[ino],disp_acc[ ino]); }
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_dv.GetValue( no[ino],disp_velo[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		double rot_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_rv.GetValue( no[ino],rot_velo[ ino]); }
		double rot_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_ra.GetValue( no[ino],rot_acc[ ino]); }
		////////////////

		double eKmat_dd[nno][nno][ndim][ndim];
		double eKmat_tt[nno][nno][ndim][ndim];
		double eKmat_dt[nno][nno][ndim][ndim];
		double eKmat_td[nno][nno][ndim][ndim];
		double area;
		double loc_base[nno][ndim];
		{
			double coord2[nno][2];
			// local coordinate base ( loc_coord_base ) ‚ðì‚é
			MakeLocalCoordBase(loc_base,coord2,coord);
			area = TriArea(coord2[0],coord2[1],coord2[2]);

			double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
			{
				double tmp_w[3], tmp_r[3][2];
				double w[3] = { 0,0,0 };
				double r[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
					tmp_w, tmp_r,
					young,poisson,thickness,
					coord2,w,r);
			}
			const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

			double emat_uu[3][3][2][2];
			{
				double eres_u[3][2];
				double disp[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
			}

			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eKmat_dd[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
						+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
					eKmat_dt[ino][jno][idim][jdim] 
						= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
						+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
					eKmat_td[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
						+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
					eKmat_tt[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
				}
				}	
				}
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
						eKmat_tt[ino][ino][idim][jdim] += loc_base[2][idim]*loc_base[2][jdim]*torsion_stiff;
				}
				}
			}
		}
		double eMmat_dd[nno][nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area/12.0;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat_dd[ino][jno][0][0] += dtmp0;
					eMmat_dd[ino][jno][1][1] += dtmp0;
					eMmat_dd[ino][jno][2][2] += dtmp0;
				}
				eMmat_dd[ino][ino][0][0] += dtmp0;
				eMmat_dd[ino][ino][1][1] += dtmp0;
				eMmat_dd[ino][ino][2][2] += dtmp0;
			}
		}

		double emat_dd[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		double eres_d[nno][ndim];
		double eres_t[nno][ndim];
		{
			double dtmp1 = beta_newmark*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_dd[0][0][0][0])[i] = (&eMmat_dd[0][0][0][0])[i]+dtmp1*(&eKmat_dd[0][0][0][0])[i];
				(&emat_dt[0][0][0][0])[i] = dtmp1*(&eKmat_dt[0][0][0][0])[i];
				(&emat_td[0][0][0][0])[i] = dtmp1*(&eKmat_td[0][0][0][0])[i];
				(&emat_tt[0][0][0][0])[i] = dtmp1*(&eKmat_tt[0][0][0][0])[i];
			}
		}

		////////////////		
		{
			double dtmp1 = area/3.0*arearho;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmp1*g_x;
				eres_d[ino][1] = dtmp1*g_y;
				eres_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eres_d[ino][0] = dtmpx;
				eres_d[ino][1] = dtmpy;
				eres_d[ino][2] = dtmpz;
			}
		}

		for(unsigned int i=0;i<nno*ndim; i++){ *(&eres_t[0][0]+i) = 0.0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= eKmat_dd[ino][jno][idim][jdim]*disp[jno][jdim] 
				               + eKmat_dt[ino][jno][idim][jdim]*rot[ jno][jdim];
			eres_d[ino][idim] -= eMmat_dd[ino][jno][idim][jdim]*disp_acc[jno][jdim]; 
			eres_t[ino][idim] -= eKmat_td[ino][jno][idim][jdim]*disp[jno][jdim] 
				               + eKmat_tt[ino][jno][idim][jdim]*rot[ jno][jdim];	
		}
		}
		}
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= dt*(       eKmat_dd[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_dt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_d[ino][idim] -= 0.5*dt*dt*(eKmat_dd[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_dt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
			////////////////////////////////
			eres_t[ino][idim] -= dt*(       eKmat_td[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_tt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_t[ino][idim] -= 0.5*dt*dt*(eKmat_td[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_tt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
				                         
		}
		}
		}
		}

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = emat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = emat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = emat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = emat_tt[ino][jno][idim][jdim];
		}
		}
    	}
		}
		mat_dd.Mearge(nno,no,nno,no,6*6,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
			res_d.AddValue(no[ino],3,eres_t[ino][0]);
			res_d.AddValue(no[ino],4,eres_t[ino][1]);
			res_d.AddValue(no[ino],5,eres_t[ino][2]);
		}
	}
	return true;
}




bool Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic(
	double dt, double gamma_newmark, double beta_newmark,
	Fem::Ls::CLinearSystem_Field& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press,
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != VECTOR3 ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_deflect,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_rot,CORNER,world) )
            {
			    return AddLinearSystem_DKT3D_LinearNonStatic_P1_Combined(
				    dt, gamma_newmark, beta_newmark, ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_LinearNonStatic_P1(
				    dt, gamma_newmark, beta_newmark, ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
		}
        assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic(
					dt, gamma_newmark, beta_newmark,
					ls,
					young, poisson, thickness, arearho,
					g_x, g_y, g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


static bool AddLinearSystem_DKT3D_LinearNonStatic_P1_Save(
		CLinearSystem_SaveDiaM_NewmarkBeta& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D LinearNonStatic P1 Save" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

//	const double dt = ls.GetDt();
//	const double gamma_newmark = ls.GetGamma();
//	const double beta_newmark = ls.GetBeta();

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

//	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,world);
//	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,world);
//	CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp, CORNER, id_field_theta,CORNER, world);
//	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER, id_field_disp, CORNER, world);
//	CMat_BlkCrs& mat_dd_bound = ls.GetMatrix_Boundary(id_field_disp,CORNER,  id_field_disp,CORNER,  world);
//	CMat_BlkCrs& mat_tt_bound = ls.GetMatrix_Boundary(id_field_theta,CORNER, id_field_theta,CORNER, world);
//	CMat_BlkCrs& mat_dt_bound = ls.GetMatrix_Boundary(id_field_disp,CORNER,  id_field_theta,CORNER, world);
//	CMat_BlkCrs& mat_td_bound = ls.GetMatrix_Boundary(id_field_theta,CORNER, id_field_disp,CORNER,  world);
//	CVector_Blk& force_d = ls.GetForce(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d  = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r  = field_rot.GetNodeSeg( CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double disp_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_da.GetValue( no[ino],disp_acc[ ino]); }
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_dv.GetValue( no[ino],disp_velo[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		double rot_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_rv.GetValue( no[ino],rot_velo[ ino]); }
		double rot_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_ra.GetValue( no[ino],rot_acc[ ino]); }
		////////////////

		double eKmat_dd[nno][nno][ndim][ndim];
		double eKmat_tt[nno][nno][ndim][ndim];
		double eKmat_dt[nno][nno][ndim][ndim];
		double eKmat_td[nno][nno][ndim][ndim];
		double area;
		double loc_base[nno][ndim];
		{
			double coord2[nno][2];
			// local coordinate base ( loc_coord_base ) ‚ðì‚é
			MakeLocalCoordBase(loc_base,coord2,coord);
			area = TriArea(coord2[0],coord2[1],coord2[2]);

			double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
			{
				double tmp_w[3], tmp_r[3][2];
				double w[3] = { 0,0,0 };
				double r[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
					tmp_w, tmp_r,
					young,poisson,thickness,
					coord2,w,r);
			}
			const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

			double emat_uu[3][3][2][2];
			{
				double eres_u[3][2];
				double disp[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
			}

			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eKmat_dd[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
						+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
					eKmat_dt[ino][jno][idim][jdim] 
						= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
						+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
					eKmat_td[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
						+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
					eKmat_tt[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
				}
				}	
				}
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
						eKmat_tt[ino][ino][idim][jdim] += loc_base[2][idim]*loc_base[2][jdim]*torsion_stiff;
				}
				}
			}
		}
		double eMmat_dd[nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area/3.0;
			for(unsigned int i=0;i<nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				eMmat_dd[ino][0][0] = dtmp0;
				eMmat_dd[ino][1][1] = dtmp0;
				eMmat_dd[ino][2][2] = dtmp0;
			}
		}

		////////////////		
		double eforce_d[nno][ndim];
		{
			double dtmp1 = area*arearho/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmp1*g_x;
				eforce_d[ino][1] = dtmp1*g_y;
				eforce_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmpx;
				eforce_d[ino][1] = dtmpy;
				eforce_d[ino][2] = dtmpz;
			}
		}

        // –¢Š®¬

		assert(0);
/*		mat_dd->Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0],tmp_buffer);
		mat_dt->Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0],tmp_buffer);
		mat_td->Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0],tmp_buffer);
		mat_tt->Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0],tmp_buffer);
		for(unsigned int ino=0;ino<nno;ino++){
			force_d->AddValue(no[ino],0,eforce_d[ino][0]);
			force_d->AddValue(no[ino],1,eforce_d[ino][1]);
			force_d->AddValue(no[ino],2,eforce_d[ino][2]);
		}*/
	}
	return true;
}



static bool AddLinearSystem_DKT3D_LinearNonStatic_P1_Save_Combined(
		CLinearSystem_SaveDiaM_NewmarkBeta& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press, 
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D LinearNonStatic P1 Save (Combined)" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const double dt = ls.GetDt();
//	const double gamma_newmark = ls.GetGamma();
	const double beta_newmark = ls.GetBeta();

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd    = ls.GetMatrix(id_field_disp,CORNER,world);
	CMat_BlkCrs& mat_dd_bound = ls.GetMatrix_Boundary(id_field_disp,CORNER,  id_field_disp,CORNER,  world);
	CDiaMat_Blk& mat_mass     = ls.GetDiaMassMatrix(id_field_disp,CORNER,world);
	CVector_Blk& force_d      = ls.GetForce(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d  = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r  = field_rot.GetNodeSeg( CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);
		double coord[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],coord[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue( no[ino],disp[ ino]); }
		double disp_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_da.GetValue( no[ino],disp_acc[ ino]); }
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_dv.GetValue( no[ino],disp_velo[ ino]); }
		double rot[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot[  ino]); }
		double rot_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_rv.GetValue( no[ino],rot_velo[ ino]); }
		double rot_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_ra.GetValue( no[ino],rot_acc[ ino]); }
		////////////////

		double eKmat_dd[nno][nno][ndim][ndim];
		double eKmat_tt[nno][nno][ndim][ndim];
		double eKmat_dt[nno][nno][ndim][ndim];
		double eKmat_td[nno][nno][ndim][ndim];
		double area;
		double loc_base[nno][ndim];
		{
			double coord2[nno][2];
			// local coordinate base ( loc_coord_base ) ‚ðì‚é
			MakeLocalCoordBase(loc_base,coord2,coord);
			area = TriArea(coord2[0],coord2[1],coord2[2]);

			double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
			{
				double tmp_w[3], tmp_r[3][2];
				double w[3] = { 0,0,0 };
				double r[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,  
					tmp_w, tmp_r,
					young,poisson,thickness,
					coord2,w,r);
			}
			const double torsion_stiff = young*thickness*thickness*thickness*area*1.0e-6;

			double emat_uu[3][3][2][2];
			{
				double eres_u[3][2];
				double disp[3][2] = { {0,0}, {0,0}, {0,0} };
				MakeStiffMat_Membrane(emat_uu,eres_u,   1.0,0.0,1.0,   coord2,disp);
			}

			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					eKmat_dd[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_uu[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_uu[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_uu[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_uu[ino][jno][1][1]
						+ loc_base[2][idim]*loc_base[2][jdim]*emat_ww[ino][jno];
					eKmat_dt[ino][jno][idim][jdim] 
						= loc_base[2][idim]*loc_base[0][jdim]*emat_wr[ino][jno][0]
						+ loc_base[2][idim]*loc_base[1][jdim]*emat_wr[ino][jno][1];
					eKmat_td[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[2][jdim]*emat_rw[ino][jno][0]
						+ loc_base[1][idim]*loc_base[2][jdim]*emat_rw[ino][jno][1];
					eKmat_tt[ino][jno][idim][jdim] 
						= loc_base[0][idim]*loc_base[0][jdim]*emat_rr[ino][jno][0][0]
						+ loc_base[0][idim]*loc_base[1][jdim]*emat_rr[ino][jno][0][1]
						+ loc_base[1][idim]*loc_base[0][jdim]*emat_rr[ino][jno][1][0]
						+ loc_base[1][idim]*loc_base[1][jdim]*emat_rr[ino][jno][1][1];
				}
				}	
				}
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
						eKmat_tt[ino][ino][idim][jdim] += loc_base[2][idim]*loc_base[2][jdim]*torsion_stiff;
				}
				}
			}
		}
		double eMmat_dd[nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area/3.0;
			for(unsigned int i=0;i<nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				eMmat_dd[ino][0][0] = dtmp0;
				eMmat_dd[ino][1][1] = dtmp0;
				eMmat_dd[ino][2][2] = dtmp0;
			}
		}

		////////////////		
		double eforce_d[nno][ndim];
		{
			double dtmp1 = area*arearho/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmp1*g_x;
				eforce_d[ino][1] = dtmp1*g_y;
				eforce_d[ino][2] = dtmp1*g_z;
			}
			const double dtmpx = loc_base[2][0]*area*press/3.0;
			const double dtmpy = loc_base[2][1]*area*press/3.0;
			const double dtmpz = loc_base[2][2]*area*press/3.0;
			for(unsigned int ino=0;ino<nno;ino++){
				eforce_d[ino][0] = dtmpx;
				eforce_d[ino][1] = dtmpy;
				eforce_d[ino][2] = dtmpz;
			}
		}

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = eKmat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = eKmat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = eKmat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = eKmat_tt[ino][jno][idim][jdim];
		}
		}
		}
		}
		mat_dd_bound.Mearge(nno,no, nno,no, 6*6, &emat[0][0][0][0]);

		////////////////////////////////

		double eMmat[nno][6][6];
		for(unsigned int i=0;i<nno*6*6;i++){ *((&eMmat[0][0][0])+i) = 0; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eMmat[ino][idim][jdim] = eMmat_dd[ino][idim][jdim];
		}
		}
		}
		for(unsigned int ino=0;ino<nno;ino++){
			mat_mass.Mearge(no[ino],6*6,&eMmat[ino][0][0]);
		}
		double dtmp1 = beta_newmark*dt*dt;
		for(unsigned int i=0;i<nno*nno*6*6;i++){ *((&emat[0][0][0][0])+i) *= dtmp1; }
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][ino][idim][jdim] += eMmat_dd[ino][idim][jdim];
		}
		}
		}
		mat_dd.Mearge(      nno,no,nno,no,6*6,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			force_d.AddValue(no[ino],0,eforce_d[ino][0]);
			force_d.AddValue(no[ino],1,eforce_d[ino][1]);
			force_d.AddValue(no[ino],2,eforce_d[ino][2]);
		}
	}
	return true;
}




bool Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic_Save(
	Fem::Ls::CLinearSystem_SaveDiaM_NewmarkBeta& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press, 
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
    if( !world.IsIdField(id_field_deflect) ){ assert(0); return false; }
	const CField& field_deflect = world.GetField(id_field_deflect);
    if( field_deflect.GetFieldType() != VECTOR3 ){ assert(0); return false; }

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_deflect,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_rot,CORNER,world) ){
			    return AddLinearSystem_DKT3D_LinearNonStatic_P1_Save_Combined(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_LinearNonStatic_P1_Save(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,id_ea);
            }
		}
		else{
			assert(0);
        }
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_Linear_NonStatic_Save(
					ls,
					young, poisson, thickness, arearho,
					g_x, g_y, g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////




static inline void MakeDerivativeRotateVector(double drdt[][3], const Com::CVector3D& arot, const Com::CVector3D& vec0){
	const double t = arot.Length();
	Com::CVector3D wrot;
	const double tant = tan(t*0.5);
	if( t < 1.0e-30 ){ wrot = arot; }
	else{ wrot = arot * (2*tant/t); }
	double drdw[3][3];
	{
		const double dtmp1 = 1.0 / ( 1.0 + 0.25*wrot.DLength() );
		const double dtmp2 = Com::Dot(vec0,wrot);
		drdw[0][0] = dtmp1*(       -0.5*wrot.x*vec0.x+0.5*dtmp2);
		drdw[0][1] = dtmp1*(+vec0.z-0.5*wrot.x*vec0.y          );
		drdw[0][2] = dtmp1*(-vec0.y-0.5*wrot.x*vec0.z          );
		drdw[1][0] = dtmp1*(-vec0.z-0.5*wrot.y*vec0.x          );
		drdw[1][1] = dtmp1*(       -0.5*wrot.y*vec0.y+0.5*dtmp2);
		drdw[1][2] = dtmp1*(+vec0.x-0.5*wrot.y*vec0.z          );
		drdw[2][0] = dtmp1*(+vec0.y-0.5*wrot.z*vec0.x          );
		drdw[2][1] = dtmp1*(-vec0.x-0.5*wrot.z*vec0.y          );
		drdw[2][2] = dtmp1*(       -0.5*wrot.z*vec0.z+0.5*dtmp2);
	}
	double dwdt[3][3];
	{
		double dtmp1;
		if( t > 1.0e-30 ){ 	dtmp1 = (1+tant*tant)/(t*t) - 2*tant/(t*t*t); }
		else{ dtmp1 = 0.0; }
		double dtmp2;
		if( t > 1.0e-30 ){ dtmp2 = 2*tant/t; }
		else{ dtmp2 = 1.0; }
		dwdt[0][0] = dtmp1 * arot.x * arot.x + dtmp2;
		dwdt[0][1] = dtmp1 * arot.x * arot.y;
		dwdt[0][2] = dtmp1 * arot.x * arot.z;
		dwdt[1][0] = dtmp1 * arot.y * arot.x;
		dwdt[1][1] = dtmp1 * arot.y * arot.y + dtmp2;
		dwdt[1][2] = dtmp1 * arot.y * arot.z;
		dwdt[2][0] = dtmp1 * arot.z * arot.x;
		dwdt[2][1] = dtmp1 * arot.z * arot.y;
		dwdt[2][2] = dtmp1 * arot.z * arot.z + dtmp2;
	}
	for(unsigned int i=0;i<3;i++){
	for(unsigned int j=0;j<3;j++){
		drdt[i][j] = 0;
		for(unsigned int k=0;k<3;k++){
			drdt[i][j] += drdw[i][k]*dwdt[k][j];
		}
	}
	}
}


static void MakeLocalCoordBaseDerivative(double lbd[][3][3][3], const double coord[][3])
{
	Com::CVector3D vec01(coord[1][0]-coord[0][0],coord[1][1]-coord[0][1],coord[1][2]-coord[0][2]);
	Com::CVector3D vec02(coord[2][0]-coord[0][0],coord[2][1]-coord[0][1],coord[2][2]-coord[0][2]);
	const double len01 = vec01.Length();
	const double len02 = vec02.Length();
	const Com::CVector3D e0 = vec01/len01;
	const Com::CVector3D et = vec02*(1/len02);
	const Com::CVector3D eu = Com::Cross(e0,et);
	const double lenu = eu.Length();
	const Com::CVector3D e2 = eu*(1/lenu);
	const Com::CVector3D e1 = Com::Cross(e2,e0);

	const double dtmp0 = 1.0/len01;
	const double dtmp1 = 1.0/(len01*len01*len01);
	const double dtmp2 = 1.0/len02;
	const double dtmp3 = 1.0/(len02*len02*len02);
	const double dtmp4 = 1.0/lenu;
	const double dtmp5 = 1.0/(lenu*lenu*lenu);

	{
		Com::CVector3D e0c0x = Com::CVector3D(-dtmp0,0,0)+vec01*(+vec01.x*dtmp1);
		Com::CVector3D etc0x = Com::CVector3D(-dtmp2,0,0)+vec02*(+vec02.x*dtmp3);
		Com::CVector3D euc0x = Com::Cross(e0c0x,et) + Com::Cross(e0,etc0x);
		Com::CVector3D e2c0x = euc0x*dtmp4+eu*(-Com::Dot(euc0x,eu)*dtmp5);
		Com::CVector3D e1c0x = Com::Cross(e2c0x,e0) + Com::Cross(e2,e0c0x);
		lbd[0][0][0][0]=e0c0x.x; lbd[0][1][0][0]=e0c0x.y; lbd[0][2][0][0]=e0c0x.z; 
		lbd[1][0][0][0]=e1c0x.x; lbd[1][1][0][0]=e1c0x.y; lbd[1][2][0][0]=e1c0x.z; 
		lbd[2][0][0][0]=e2c0x.x; lbd[2][1][0][0]=e2c0x.y; lbd[2][2][0][0]=e2c0x.z; 
	}
	{
		Com::CVector3D e0c0y = Com::CVector3D(0,-dtmp0,0)+vec01*(+vec01.y*dtmp1);
		Com::CVector3D etc0y = Com::CVector3D(0,-dtmp2,0)+vec02*(+vec02.y*dtmp3);
		Com::CVector3D euc0y = Com::Cross(e0c0y,et) + Com::Cross(e0,etc0y);
		Com::CVector3D e2c0y = euc0y*dtmp4+eu*(-Com::Dot(euc0y,eu)*dtmp5);
		Com::CVector3D e1c0y = Com::Cross(e2c0y,e0) + Com::Cross(e2,e0c0y);
		lbd[0][0][0][1]=e0c0y.x; lbd[0][1][0][1]=e0c0y.y; lbd[0][2][0][1]=e0c0y.z;
		lbd[1][0][0][1]=e1c0y.x; lbd[1][1][0][1]=e1c0y.y; lbd[1][2][0][1]=e1c0y.z;
		lbd[2][0][0][1]=e2c0y.x; lbd[2][1][0][1]=e2c0y.y; lbd[2][2][0][1]=e2c0y.z;
	}
	{
		Com::CVector3D e0c0z = Com::CVector3D(0,0,-dtmp0)+vec01*(+vec01.z*dtmp1);
		Com::CVector3D etc0z = Com::CVector3D(0,0,-dtmp2)+vec02*(+vec02.z*dtmp3);
		Com::CVector3D euc0z = Com::Cross(e0c0z,et) + Com::Cross(e0,etc0z);
		Com::CVector3D e2c0z = euc0z*dtmp4+eu*(-Com::Dot(euc0z,eu)*dtmp5);
		Com::CVector3D e1c0z = Com::Cross(e2c0z,e0) + Com::Cross(e2,e0c0z);
		lbd[0][0][0][2]=e0c0z.x; lbd[0][1][0][2]=e0c0z.y; lbd[0][2][0][2]=e0c0z.z;
		lbd[1][0][0][2]=e1c0z.x; lbd[1][1][0][2]=e1c0z.y; lbd[1][2][0][2]=e1c0z.z;
		lbd[2][0][0][2]=e2c0z.x; lbd[2][1][0][2]=e2c0z.y; lbd[2][2][0][2]=e2c0z.z;
	}
	////////////////
	{
		Com::CVector3D e0c1x = Com::CVector3D( dtmp0,0,0)+vec01*(-vec01.x*dtmp1);
		Com::CVector3D etc1x(0,0,0);
		Com::CVector3D euc1x = Com::Cross(e0c1x,et) + Com::Cross(e0,etc1x);
		Com::CVector3D e2c1x = euc1x*dtmp4+eu*(-Com::Dot(euc1x,eu)*dtmp5);
		Com::CVector3D e1c1x = Com::Cross(e2c1x,e0) + Com::Cross(e2,e0c1x);
		lbd[0][0][1][0]=e0c1x.x; lbd[0][1][1][0]=e0c1x.y; lbd[0][2][1][0]=e0c1x.z; 
		lbd[1][0][1][0]=e1c1x.x; lbd[1][1][1][0]=e1c1x.y; lbd[1][2][1][0]=e1c1x.z; 
		lbd[2][0][1][0]=e2c1x.x; lbd[2][1][1][0]=e2c1x.y; lbd[2][2][1][0]=e2c1x.z; 
	}
	{
		Com::CVector3D e0c1y = Com::CVector3D(0, dtmp0,0)+vec01*(-vec01.y*dtmp1);
		Com::CVector3D etc1y(0,0,0);
		Com::CVector3D euc1y = Com::Cross(e0c1y,et) + Com::Cross(e0,etc1y);
		Com::CVector3D e2c1y = euc1y*dtmp4+eu*(-Com::Dot(euc1y,eu)*dtmp5);
		Com::CVector3D e1c1y = Com::Cross(e2c1y,e0) + Com::Cross(e2,e0c1y);
		lbd[0][0][1][1]=e0c1y.x; lbd[0][1][1][1]=e0c1y.y; lbd[0][2][1][1]=e0c1y.z;
		lbd[1][0][1][1]=e1c1y.x; lbd[1][1][1][1]=e1c1y.y; lbd[1][2][1][1]=e1c1y.z;
		lbd[2][0][1][1]=e2c1y.x; lbd[2][1][1][1]=e2c1y.y; lbd[2][2][1][1]=e2c1y.z;
	}
	{
		Com::CVector3D e0c1z = Com::CVector3D(0,0, dtmp0)+vec01*(-vec01.z*dtmp1);
		Com::CVector3D etc1z(0,0,0);
		Com::CVector3D euc1z = Com::Cross(e0c1z,et) + Com::Cross(e0,etc1z);
		Com::CVector3D e2c1z = euc1z*dtmp4+eu*(-Com::Dot(euc1z,eu)*dtmp5);
		Com::CVector3D e1c1z = Com::Cross(e2c1z,e0) + Com::Cross(e2,e0c1z);
		lbd[0][0][1][2]=e0c1z.x; lbd[0][1][1][2]=e0c1z.y; lbd[0][2][1][2]=e0c1z.z;
		lbd[1][0][1][2]=e1c1z.x; lbd[1][1][1][2]=e1c1z.y; lbd[1][2][1][2]=e1c1z.z;
		lbd[2][0][1][2]=e2c1z.x; lbd[2][1][1][2]=e2c1z.y; lbd[2][2][1][2]=e2c1z.z;
	}
	////////////////
	{
		Com::CVector3D e0c2x(0,0,0);
		Com::CVector3D etc2x = Com::CVector3D( dtmp2,0,0)+vec02*(-vec02.x*dtmp3);
		Com::CVector3D euc2x = Com::Cross(e0c2x,et) + Com::Cross(e0,etc2x);
		Com::CVector3D e2c2x = euc2x*dtmp4+eu*(-Com::Dot(euc2x,eu)*dtmp5);
		Com::CVector3D e1c2x = Com::Cross(e2c2x,e0) + Com::Cross(e2,e0c2x);
		lbd[0][0][2][0]=e0c2x.x; lbd[0][1][2][0]=e0c2x.y; lbd[0][2][2][0]=e0c2x.z; 
		lbd[1][0][2][0]=e1c2x.x; lbd[1][1][2][0]=e1c2x.y; lbd[1][2][2][0]=e1c2x.z; 
		lbd[2][0][2][0]=e2c2x.x; lbd[2][1][2][0]=e2c2x.y; lbd[2][2][2][0]=e2c2x.z; 
	}
	{
		Com::CVector3D e0c2y(0,0,0);
		Com::CVector3D etc2y = Com::CVector3D(0, dtmp2,0)+vec02*(-vec02.y*dtmp3);
		Com::CVector3D euc2y = Com::Cross(e0c2y,et) + Com::Cross(e0,etc2y);
		Com::CVector3D e2c2y = euc2y*dtmp4+eu*(-Com::Dot(euc2y,eu)*dtmp5);
		Com::CVector3D e1c2y = Com::Cross(e2c2y,e0) + Com::Cross(e2,e0c2y);
		lbd[0][0][2][1]=e0c2y.x; lbd[0][1][2][1]=e0c2y.y; lbd[0][2][2][1]=e0c2y.z; 
		lbd[1][0][2][1]=e1c2y.x; lbd[1][1][2][1]=e1c2y.y; lbd[1][2][2][1]=e1c2y.z; 
		lbd[2][0][2][1]=e2c2y.x; lbd[2][1][2][1]=e2c2y.y; lbd[2][2][2][1]=e2c2y.z; 
	}
	{
		Com::CVector3D e0c2z(0,0,0);
		Com::CVector3D etc2z = Com::CVector3D(0,0, dtmp2)+vec02*(-vec02.z*dtmp3);
		Com::CVector3D euc2z = Com::Cross(e0c2z,et) + Com::Cross(e0,etc2z);
		Com::CVector3D e2c2z = euc2z*dtmp4+eu*(-Com::Dot(euc2z,eu)*dtmp5);
		Com::CVector3D e1c2z = Com::Cross(e2c2z,e0) + Com::Cross(e2,e0c2z);
		lbd[0][0][2][2]=e0c2z.x; lbd[0][1][2][2]=e0c2z.y; lbd[0][2][2][2]=e0c2z.z; 
		lbd[1][0][2][2]=e1c2z.x; lbd[1][1][2][2]=e1c2z.y; lbd[1][2][2][2]=e1c2z.z; 
		lbd[2][0][2][2]=e2c2z.x; lbd[2][1][2][2]=e2c2z.y; lbd[2][2][2][2]=e2c2z.z; 
	}
}

static void MakeLocalCoordBaseSecondDerivative(double lcbd2[][3][3][3][3][3], const double coord[][3])
{
	const unsigned int ndim = 3;
	const unsigned int nno = 3;
	Com::CVector3D v01( coord[1][0]-coord[0][0], coord[1][1]-coord[0][1], coord[1][2]-coord[0][2] );
	const double len01 = v01.Length();
	const double inv3len01 = 1.0/(len01*len01*len01);
	const double inv5len01 = 1.0/(len01*len01*len01*len01*len01);
	Com::CVector3D e0 = v01 / len01;
	double dotvec01[ndim][nno][ndim];
	{
		for(unsigned int i=0;i<ndim*3*ndim;i++){ *(&dotvec01[0][0][0]+i)=0.0; }
		for(unsigned int idim=0;idim<ndim;idim++){
			dotvec01[idim][0][idim] = -1;
			dotvec01[idim][1][idim] =  1;
		}
	}

	Com::CVector3D v0t( coord[2][0]-coord[0][0], coord[2][1]-coord[0][1], coord[2][2]-coord[0][2] );
	Com::CVector3D v2t = Com::Cross(e0,v0t);
	const double len2t = v2t.Length();
	const double inv3len2t = 1.0/(len2t*len2t*len2t);
	const double inv5len2t = 1.0/(len2t*len2t*len2t*len2t*len2t);
	Com::CVector3D e2 = v2t / len2t;
	double dotvec0t[ndim][nno][ndim];
	{
		for(unsigned int i=0;i<ndim*3*ndim;i++){ *(&dotvec0t[0][0][0]+i)=0.0; }
		for(unsigned int idim=0;idim<ndim;idim++){
			dotvec0t[idim][0][idim] = -1;
			dotvec0t[idim][2][idim] = +1;
		}
	}
	double lcbd[ndim][ndim][nno][ndim];
	MakeLocalCoordBaseDerivative(lcbd,coord);
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			if( jno*ndim+jdim < ino*ndim+idim ) continue;	// compute only lower triangle matrix
			Com::CVector3D dd_e0;
			{
				Com::CVector3D omega( dotvec01[0][ino][idim], dotvec01[1][ino][idim], dotvec01[2][ino][idim] );
				Com::CVector3D epsil( dotvec01[0][jno][jdim], dotvec01[1][jno][jdim], dotvec01[2][jno][jdim] );
				const double dtmp1 = Com::Dot(omega,v01);
				const double dtmp2 = Com::Dot(epsil,v01);
				const double dtmp3 = Com::Dot(omega,epsil);
				dd_e0 = (3*dtmp1*dtmp2*inv5len01-dtmp3*inv3len01)*v01 - dtmp2*inv3len01*omega - dtmp1*inv3len01*epsil;
			}
			lcbd2[0][0][ino][jno][idim][jdim] = dd_e0.x;
			lcbd2[0][1][ino][jno][idim][jdim] = dd_e0.y;
			lcbd2[0][2][ino][jno][idim][jdim] = dd_e0.z;
			////////////////
			Com::CVector3D d0_e0( lcbd[0][0][ino][idim], lcbd[0][1][ino][idim], lcbd[0][2][ino][idim] );
			Com::CVector3D d1_e0( lcbd[0][0][jno][jdim], lcbd[0][1][jno][jdim], lcbd[0][2][jno][jdim] );
			Com::CVector3D dd_e2;
			{
				Com::CVector3D d0_v0t( dotvec0t[0][ino][idim], dotvec0t[1][ino][idim], dotvec0t[2][ino][idim] );
				Com::CVector3D d1_v0t( dotvec0t[0][jno][jdim], dotvec0t[1][jno][jdim], dotvec0t[2][jno][jdim] );
				Com::CVector3D d0_v2t = Com::Cross(d0_e0,v0t) + Com::Cross(   e0,d0_v0t);
				Com::CVector3D d1_v2t = Com::Cross(d1_e0,v0t) + Com::Cross(   e0,d1_v0t);
				Com::CVector3D dd_v2t = Com::Cross(dd_e0,v0t) + Com::Cross(d0_e0,d1_v0t) + Com::Cross(d1_e0,d0_v0t);
				const double dtmp0 = 1.0/len2t;
				const double dtmp1 = -Com::Dot(d1_v2t,v2t)*inv3len2t;
				const double dtmp2 = -Com::Dot(d0_v2t,v2t)*inv3len2t;
				const double dtmp3 =  Com::Dot(d0_v2t,v2t)*Com::Dot(d1_v2t,v2t)*3*inv5len2t - ( Com::Dot(dd_v2t,v2t) + Com::Dot(d0_v2t,d1_v2t) )*inv3len2t;
				dd_e2 = dtmp0*dd_v2t + dtmp1*d0_v2t + dtmp2*d1_v2t + dtmp3*v2t;
				lcbd2[2][0][ino][jno][idim][jdim] = dd_e2.x;
				lcbd2[2][1][ino][jno][idim][jdim] = dd_e2.y;
				lcbd2[2][2][ino][jno][idim][jdim] = dd_e2.z;
			}
			////////////////
			Com::CVector3D d0_e2( lcbd[2][0][ino][idim], lcbd[2][1][ino][idim], lcbd[2][2][ino][idim] );
			Com::CVector3D d1_e2( lcbd[2][0][jno][jdim], lcbd[2][1][jno][jdim], lcbd[2][2][jno][jdim] );
			Com::CVector3D dd_e1;
			{
				dd_e1 = Com::Cross(dd_e2,e0) + Com::Cross(d1_e2,d0_e0) + Com::Cross(d0_e2,d1_e0) + Com::Cross(e2,dd_e0);
				lcbd2[1][0][ino][jno][idim][jdim] = dd_e1.x;
				lcbd2[1][1][ino][jno][idim][jdim] = dd_e1.y;
				lcbd2[1][2][ino][jno][idim][jdim] = dd_e1.z;
			}
		}
		}
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			if( jno*ndim+jdim <= ino*ndim+idim ) continue;
			lcbd2[0][0][jno][ino][jdim][idim] = lcbd2[0][0][ino][jno][idim][jdim];
			lcbd2[0][1][jno][ino][jdim][idim] = lcbd2[0][1][ino][jno][idim][jdim];
			lcbd2[0][2][jno][ino][jdim][idim] = lcbd2[0][2][ino][jno][idim][jdim];
			lcbd2[1][0][jno][ino][jdim][idim] = lcbd2[1][0][ino][jno][idim][jdim];
			lcbd2[1][1][jno][ino][jdim][idim] = lcbd2[1][1][ino][jno][idim][jdim];
			lcbd2[1][2][jno][ino][jdim][idim] = lcbd2[1][2][ino][jno][idim][jdim];
			lcbd2[2][0][jno][ino][jdim][idim] = lcbd2[2][0][ino][jno][idim][jdim];
			lcbd2[2][1][jno][ino][jdim][idim] = lcbd2[2][1][ino][jno][idim][jdim];
			lcbd2[2][2][jno][ino][jdim][idim] = lcbd2[2][2][ino][jno][idim][jdim];
		}
		}
	}
	}
/*	double delta = 1.0e-10;
	double cord_delta[3] = { 1.3*delta, -0.7*delta, 0.8*delta };
	for(unsigned int inode=0;inode<nno;inode++){
		double cord3d_tmp[3][3];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			cord3d_tmp[ino][idim] = coord[ino][idim];
		}
		}
		for(unsigned int idim=0;idim<ndim;idim++){
			cord3d_tmp[inode][idim] += cord_delta[idim];
		}
		double lcbd_tmp[ndim][ndim][nno][ndim];
		MakeLocalCoordBaseDerivative(lcbd_tmp,cord3d_tmp);
		double lcbd[ndim][ndim][nno][ndim];
		MakeLocalCoordBaseDerivative(lcbd,coord);
		const unsigned int jno0 = 2;
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			const double dtmp0 = lcbd_tmp[jno0][0][ino][idim] - lcbd[jno0][0][ino][idim];
			const double dtmp1 = lcbd_tmp[jno0][1][ino][idim] - lcbd[jno0][1][ino][idim];
			const double dtmp2 = lcbd_tmp[jno0][2][ino][idim] - lcbd[jno0][2][ino][idim];
			double dtmp3=0,dtmp4=0,dtmp5=0;
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				dtmp3 += lcbd2[jno0][0][ino][inode][idim][jdim]*cord_delta[jdim];
				dtmp4 += lcbd2[jno0][1][ino][inode][idim][jdim]*cord_delta[jdim];
				dtmp5 += lcbd2[jno0][2][ino][inode][idim][jdim]*cord_delta[jdim];
			}
			const double dtmp6 = sqrt( (dtmp0-dtmp3)*(dtmp0-dtmp3) + (dtmp1-dtmp4)*(dtmp1-dtmp4) + (dtmp2-dtmp5)*(dtmp2-dtmp5) )/delta;
			std::cout << inode << " " << ino << " " << idim << " --> " << dtmp6 << std::endl;
//			std::cout << dtmp0 << " " << dtmp3 << std::endl;
//			std::cout << dtmp1 << " " << dtmp4 << std::endl;
//			std::cout << dtmp2 << " " << dtmp5 << std::endl;
		}
		}
	}*/
}


void MakeElemStiffMatDKT3D_Nonlinear(
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press,
	const double cord3d0[][3], const double disp3d[][3], const double rot3d[][3],
	double& area0,
	double emat_dd[][3][3][3], double emat_dt[][3][3][3], double emat_td[][3][3][3], double emat_tt[][3][3][3], 
	double eres_d[][3], double eres_t[][3] )
{
	

	const unsigned int nno = 3;
	const unsigned int ndim = 3;
	
	double cord3d1[nno][ndim];
	for(unsigned int i=0;i<nno*ndim;i++){ *(&cord3d1[0][0]+i) = *(&cord3d0[0][0]+i) + *(&disp3d[0][0]+i); }

	double loc_base0[nno][ndim], cord2d0[nno][2];
	MakeLocalCoordBase(loc_base0,cord2d0,cord3d0);
	area0 = TriArea(cord2d0[0],cord2d0[1],cord2d0[2]);

	double loc_base1[nno][ndim], cord2d1[nno][2];
	MakeLocalCoordBase(loc_base1,cord2d1,cord3d1);

	double rot2d1[nno][2], w2d1[nno], disp2d1[nno][2], arot2d1[nno];
	double loc_base2[nno][2][ndim];
	for(unsigned int ino=0;ino<nno;ino++){
		const Com::CVector3D dirc3d0(loc_base0[2][0],loc_base0[2][1],loc_base0[2][2]);
		Com::CVector3D rot3di(rot3d[ino][0],rot3d[ino][1],rot3d[ino][2]);
		Com::CVector3D dirc3d1 = Com::RotateVector(dirc3d0,rot3di);
		disp2d1[ino][0] =  cord2d1[ino][0]-cord2d0[ino][0];
		disp2d1[ino][1] =  cord2d1[ino][1]-cord2d0[ino][1];
		w2d1[   ino]    =  0;
		rot2d1[ ino][0] = -loc_base1[1][0]*dirc3d1.x-loc_base1[1][1]*dirc3d1.y-loc_base1[1][2]*dirc3d1.z;
		rot2d1[ ino][1] =  loc_base1[0][0]*dirc3d1.x+loc_base1[0][1]*dirc3d1.y+loc_base1[0][2]*dirc3d1.z;
		arot2d1[ino]    =  Com::Dot(dirc3d0,rot3di); // Ž²•ûŒü‰ñ“]
		{	// loc_base2(”÷¬‰ñ“]‚µ‚½Œã‚É‰ñ“])‚ðì‚é
			double dtdo[2][2];
			const Com::CVector3D ax3d0(loc_base0[0][0],loc_base0[0][1],loc_base0[0][2]);
			const Com::CVector3D ay3d0(loc_base0[1][0],loc_base0[1][1],loc_base0[1][2]);
			const Com::CVector3D ax3d1 = Com::RotateVector(ax3d0,rot3di);
			const Com::CVector3D ay3d1 = Com::RotateVector(ay3d0,rot3di);
			dtdo[0][0] =  loc_base1[1][0]*ay3d1.x+loc_base1[1][1]*ay3d1.y+loc_base1[1][2]*ay3d1.z;
			dtdo[0][1] = -loc_base1[1][0]*ax3d1.x-loc_base1[1][1]*ax3d1.y-loc_base1[1][2]*ax3d1.z;
			dtdo[1][0] = -loc_base1[0][0]*ay3d1.x-loc_base1[0][1]*ay3d1.y-loc_base1[0][2]*ay3d1.z;
			dtdo[1][1] =  loc_base1[0][0]*ax3d1.x+loc_base1[0][1]*ax3d1.y+loc_base1[0][2]*ax3d1.z;
			loc_base2[ino][0][0] = dtdo[0][0]*loc_base0[0][0] + dtdo[0][1]*loc_base0[1][0];
			loc_base2[ino][0][1] = dtdo[0][0]*loc_base0[0][1] + dtdo[0][1]*loc_base0[1][1];
			loc_base2[ino][0][2] = dtdo[0][0]*loc_base0[0][2] + dtdo[0][1]*loc_base0[1][2];
			loc_base2[ino][1][0] = dtdo[1][0]*loc_base0[0][0] + dtdo[1][1]*loc_base0[1][0];
			loc_base2[ino][1][1] = dtdo[1][0]*loc_base0[0][1] + dtdo[1][1]*loc_base0[1][1];
			loc_base2[ino][1][2] = dtdo[1][0]*loc_base0[0][2] + dtdo[1][1]*loc_base0[1][2];
		}
	}
	
	double lcbd[3][3][3][3];
	MakeLocalCoordBaseDerivative(lcbd,cord3d1);

	double d2ud3u[3][nno][ndim];
	{
		for(unsigned int i=0;i<3;i++){ 
			d2ud3u[0][0][i]=-loc_base1[0][i]; d2ud3u[0][1][i]=loc_base1[0][i]; d2ud3u[0][2][i]=0; 
			d2ud3u[1][0][i]=-loc_base1[0][i]; d2ud3u[1][1][i]=0;               d2ud3u[1][2][i]=loc_base1[0][i];
			d2ud3u[2][0][i]=-loc_base1[1][i]; d2ud3u[2][1][i]=0;               d2ud3u[2][2][i]=loc_base1[1][i];
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			double dtmp0 = 0, dtmp1 = 0, dtmp2 = 0;
			for(unsigned int i=0;i<ndim;i++){
				dtmp0 += (cord3d1[1][i]-cord3d1[0][i])*lcbd[0][i][ino][idim]; 
				dtmp1 += (cord3d1[2][i]-cord3d1[0][i])*lcbd[0][i][ino][idim]; 
				dtmp2 += (cord3d1[2][i]-cord3d1[0][i])*lcbd[1][i][ino][idim]; 
			}
			d2ud3u[0][ino][idim] += dtmp0;
			d2ud3u[1][ino][idim] += dtmp1;
			d2ud3u[2][ino][idim] += dtmp2;
		}
		}
	}
	double lcbd2[ndim][ndim][nno][nno][ndim][ndim];
	MakeLocalCoordBaseSecondDerivative(lcbd2,cord3d1);

	double dd2ud3ud3u[3][nno][nno][ndim][ndim];
	{
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			dd2ud3ud3u[0][0][0][idim][jdim] = -lcbd[0][idim][0][jdim] - lcbd[0][jdim][0][idim];
			dd2ud3ud3u[0][0][1][idim][jdim] = -lcbd[0][idim][1][jdim] + lcbd[0][jdim][0][idim];
			dd2ud3ud3u[0][0][2][idim][jdim] = -lcbd[0][idim][2][jdim];
			dd2ud3ud3u[0][1][0][idim][jdim] =  lcbd[0][idim][0][jdim] - lcbd[0][jdim][1][idim];
			dd2ud3ud3u[0][1][1][idim][jdim] =  lcbd[0][idim][1][jdim] + lcbd[0][jdim][1][idim];
			dd2ud3ud3u[0][1][2][idim][jdim] =  lcbd[0][idim][2][jdim];
			dd2ud3ud3u[0][2][0][idim][jdim] =                         - lcbd[0][jdim][2][idim];
			dd2ud3ud3u[0][2][1][idim][jdim] =                         + lcbd[0][jdim][2][idim];
			dd2ud3ud3u[0][2][2][idim][jdim] =  0;
			////////////////
			dd2ud3ud3u[1][0][0][idim][jdim] = -lcbd[0][idim][0][jdim] - lcbd[0][jdim][0][idim];
			dd2ud3ud3u[1][0][1][idim][jdim] = -lcbd[0][idim][1][jdim];
			dd2ud3ud3u[1][0][2][idim][jdim] = -lcbd[0][idim][2][jdim] + lcbd[0][jdim][0][idim];
			dd2ud3ud3u[1][1][0][idim][jdim] =                         - lcbd[0][jdim][1][idim];
			dd2ud3ud3u[1][1][1][idim][jdim] =  0;
			dd2ud3ud3u[1][1][2][idim][jdim] =                         + lcbd[0][jdim][1][idim]; 
			dd2ud3ud3u[1][2][0][idim][jdim] =  lcbd[0][idim][0][jdim] - lcbd[0][jdim][2][idim]; 
			dd2ud3ud3u[1][2][1][idim][jdim] =  lcbd[0][idim][1][jdim]; 
			dd2ud3ud3u[1][2][2][idim][jdim] =  lcbd[0][jdim][2][idim] + lcbd[0][jdim][2][idim];
			////////////////
			dd2ud3ud3u[2][0][0][idim][jdim] = -lcbd[1][idim][0][jdim] - lcbd[1][jdim][0][idim];
			dd2ud3ud3u[2][0][1][idim][jdim] = -lcbd[1][idim][1][jdim];
			dd2ud3ud3u[2][0][2][idim][jdim] = -lcbd[1][idim][2][jdim] + lcbd[1][jdim][0][idim];
			dd2ud3ud3u[2][1][0][idim][jdim] =                         - lcbd[1][jdim][1][idim];
			dd2ud3ud3u[2][1][1][idim][jdim] =  0;
			dd2ud3ud3u[2][1][2][idim][jdim] =                         + lcbd[1][jdim][1][idim]; 
			dd2ud3ud3u[2][2][0][idim][jdim] =  lcbd[1][idim][0][jdim] - lcbd[1][jdim][2][idim]; 
			dd2ud3ud3u[2][2][1][idim][jdim] =  lcbd[1][idim][1][jdim]; 
			dd2ud3ud3u[2][2][2][idim][jdim] =  lcbd[1][jdim][2][idim] + lcbd[1][jdim][2][idim];
		}
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				double dtmp0 = 0, dtmp1 = 0, dtmp2 = 0;
				for(unsigned int i=0;i<3;i++){
					dtmp0 += (cord3d1[1][i]-cord3d1[0][i])*lcbd2[0][i][ino][jno][idim][jdim];
					dtmp1 += (cord3d1[2][i]-cord3d1[0][i])*lcbd2[0][i][ino][jno][idim][jdim];
					dtmp2 += (cord3d1[2][i]-cord3d1[0][i])*lcbd2[1][i][ino][jno][idim][jdim];
				}
				dd2ud3ud3u[0][ino][jno][idim][jdim] += dtmp0;
				dd2ud3ud3u[1][ino][jno][idim][jdim] += dtmp1;
				dd2ud3ud3u[2][ino][jno][idim][jdim] += dtmp2;
			}
			}
		}
		}
	}

	////////////////
	//DKT”Â‹È‚°—v‘f‚ð“ü‚ê‚é
	double emat_ww[3][3], emat_wr[3][3][2], emat_rw[3][3][2], emat_rr[3][3][2][2];
	double res_rot2d[nno][2], res_w2d[nno];
	MakeStiffMat_DKT_PlateBending(emat_ww,emat_wr,emat_rw,emat_rr,   res_w2d,res_rot2d, 
		young,poisson,thickness,
		cord2d0, w2d1, rot2d1 );
	////////////////
	//‚Ë‚¶‚è„«‚ð“ü‚ê‚é
	const double torsion_stiff = young*thickness*thickness*thickness*area0*1.0e-6;
	double res_arot2d[nno];
	for(unsigned int ino=0;ino<nno;ino++){
		res_arot2d[ino] = -torsion_stiff * arot2d1[ino];
	}
	////////////////
	//CST–Œ—v‘f‚ðì‚é
	double emat_uu[nno][nno][2][2];
	double res_disp2d[3][2];
	MakeStiffMat_Membrane(emat_uu,res_disp2d,   
		young,poisson,thickness,   
		cord2d0,disp2d1);

	////////////////////////////////////////////////////////////////

	double dodu[nno][2][nno][ndim];
	double ddodudt[nno][2][nno][ndim][ndim]; // [ a, n, ino, idim, jdim ]
	Com::CVector3D zrx(0,0,0),zry(0,0,0);
	{	
		for(unsigned int jno=0;jno<nno;jno++){
			const Com::CVector3D ax3d0(loc_base0[0][0],loc_base0[0][1],loc_base0[0][2]);
			const Com::CVector3D ay3d0(loc_base0[1][0],loc_base0[1][1],loc_base0[1][2]);
			const Com::CVector3D az3d0(loc_base0[2][0],loc_base0[2][1],loc_base0[2][2]);
			const Com::CVector3D rot3di(rot3d[jno][0],rot3d[jno][1],rot3d[jno][2]);
			const Com::CVector3D ax3d1 = Com::RotateVector(ax3d0,rot3di);
			const Com::CVector3D ay3d1 = Com::RotateVector(ay3d0,rot3di);
			const Com::CVector3D az3d1 = Com::RotateVector(az3d0,rot3di);
			zrx += az3d1 * res_rot2d[jno][0];
			zry += az3d1 * res_rot2d[jno][1];
			double ddodudr[2][2];
			for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int idim=0;idim<ndim;idim++){
				dodu[jno][0][ino][idim] = -lcbd[1][0][ino][idim]*az3d1.x - lcbd[1][1][ino][idim]*az3d1.y - lcbd[1][2][ino][idim]*az3d1.z;
				dodu[jno][1][ino][idim] =  lcbd[0][0][ino][idim]*az3d1.x + lcbd[0][1][ino][idim]*az3d1.y + lcbd[0][2][ino][idim]*az3d1.z;
				ddodudr[0][0] =  lcbd[1][0][ino][idim]*ay3d1.x + lcbd[1][1][ino][idim]*ay3d1.y + lcbd[1][2][ino][idim]*ay3d1.z;
				ddodudr[0][1] = -lcbd[1][0][ino][idim]*ax3d1.x - lcbd[1][1][ino][idim]*ax3d1.y - lcbd[1][2][ino][idim]*ax3d1.z;
				ddodudr[1][0] = -lcbd[0][0][ino][idim]*ay3d1.x - lcbd[0][1][ino][idim]*ay3d1.y - lcbd[0][2][ino][idim]*ay3d1.z;
				ddodudr[1][1] =  lcbd[0][0][ino][idim]*ax3d1.x + lcbd[0][1][ino][idim]*ax3d1.y + lcbd[0][2][ino][idim]*ax3d1.z;
				ddodudt[jno][0][ino][idim][0] = ddodudr[0][0]*loc_base0[0][0] + ddodudr[0][1]*loc_base0[1][0];
				ddodudt[jno][0][ino][idim][1] = ddodudr[0][0]*loc_base0[0][1] + ddodudr[0][1]*loc_base0[1][1];
				ddodudt[jno][0][ino][idim][2] = ddodudr[0][0]*loc_base0[0][2] + ddodudr[0][1]*loc_base0[1][2];
				ddodudt[jno][1][ino][idim][0] = ddodudr[1][0]*loc_base0[0][0] + ddodudr[1][1]*loc_base0[1][0];
				ddodudt[jno][1][ino][idim][1] = ddodudr[1][0]*loc_base0[0][1] + ddodudr[1][1]*loc_base0[1][1];
				ddodudt[jno][1][ino][idim][2] = ddodudr[1][0]*loc_base0[0][2] + ddodudr[1][1]*loc_base0[1][2];
			}
			}
		}
	}
	/*	
		double emat_dd[nno][nno][ndim][ndim];
		double emat_td[nno][nno][ndim][ndim];
		double emat_tt[nno][nno][ndim][ndim];
		double emat_dt[nno][nno][ndim][ndim];
	*/
/*
		for(unsigned int ino=0;ino<nno;ino++){
			for(unsigned int jno=0;jno<nno;jno++){
				for(unsigned int idim=0;idim<ndim;idim++){
				for(unsigned int jdim=0;jdim<ndim;jdim++){
					emat_dd[ino][jno][idim][jdim] 
						= d2ud3u[0][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[1][1][0][0]
						+ d2ud3u[0][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[1][2][0][0]
						+ d2ud3u[0][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[1][2][0][1]
						+ d2ud3u[1][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[2][1][0][0]
						+ d2ud3u[1][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[2][2][0][0]
						+ d2ud3u[1][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[2][2][0][1]
						+ d2ud3u[2][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[2][1][1][0]
						+ d2ud3u[2][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[2][2][1][0]
						+ d2ud3u[2][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[2][2][1][1]
						+ loc_base1[2][idim]*loc_base1[2][jdim]*emat_ww[ino][jno];
					emat_dt[ino][jno][idim][jdim] 
						= loc_base1[2][idim]*loc_base2[jno][0][jdim]*emat_wr[ino][jno][0]
						+ loc_base1[2][idim]*loc_base2[jno][1][jdim]*emat_wr[ino][jno][1];
					emat_td[ino][jno][idim][jdim] 
						= loc_base2[ino][0][idim]*loc_base1[2][jdim]*emat_rw[ino][jno][0]
						+ loc_base2[ino][1][idim]*loc_base1[2][jdim]*emat_rw[ino][jno][1];
					emat_tt[ino][jno][idim][jdim] 
						= loc_base2[ino][0][idim]*loc_base2[jno][0][jdim]*emat_rr[ino][jno][0][0]
						+ loc_base2[ino][0][idim]*loc_base2[jno][1][jdim]*emat_rr[ino][jno][0][1]
						+ loc_base2[ino][1][idim]*loc_base2[jno][0][jdim]*emat_rr[ino][jno][1][0]
						+ loc_base2[ino][1][idim]*loc_base2[jno][1][jdim]*emat_rr[ino][jno][1][1];
				}
				}	
			}
			for(unsigned int idim=0;idim<ndim;idim++){
				emat_tt[ino][ino][idim][idim] += loc_base0[2][idim]*loc_base0[2][idim]*torsion_stiff;
			}
		}
*/
	////////////////////////////////
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
	for(unsigned int idim=0;idim<ndim;idim++){
	for(unsigned int jdim=0;jdim<ndim;jdim++){
		if( ino*ndim+idim < jno*ndim+jdim ) continue;
		emat_dd[ino][jno][idim][jdim] 
			= d2ud3u[0][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[1][1][0][0]
			+ d2ud3u[0][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[1][2][0][0]
			+ d2ud3u[0][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[1][2][0][1]
			+ d2ud3u[1][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[2][1][0][0]
			+ d2ud3u[1][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[2][2][0][0]
			+ d2ud3u[1][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[2][2][0][1]
			+ d2ud3u[2][ino][idim]*d2ud3u[0][jno][jdim]*emat_uu[2][1][1][0]
			+ d2ud3u[2][ino][idim]*d2ud3u[1][jno][jdim]*emat_uu[2][2][1][0]
			+ d2ud3u[2][ino][idim]*d2ud3u[2][jno][jdim]*emat_uu[2][2][1][1];
		for(unsigned int kno=0;kno<nno;kno++){
		for(unsigned int lno=0;lno<nno;lno++){
			emat_dd[ino][jno][idim][jdim] 
				+= dodu[kno][0][ino][idim]*dodu[lno][0][jno][jdim]*emat_rr[kno][lno][0][0]
				 + dodu[kno][0][ino][idim]*dodu[lno][1][jno][jdim]*emat_rr[kno][lno][0][1]
				 + dodu[kno][1][ino][idim]*dodu[lno][0][jno][jdim]*emat_rr[kno][lno][1][0]
				 + dodu[kno][1][ino][idim]*dodu[lno][1][jno][jdim]*emat_rr[kno][lno][1][1];
		}
		}
		{
			const double dtmp1 = -lcbd2[1][0][ino][jno][idim][jdim]*zrx.x
							     -lcbd2[1][1][ino][jno][idim][jdim]*zrx.y
							     -lcbd2[1][2][ino][jno][idim][jdim]*zrx.z;
			const double dtmp2 =  lcbd2[0][0][ino][jno][idim][jdim]*zry.x
							     +lcbd2[0][1][ino][jno][idim][jdim]*zry.y
							     +lcbd2[0][2][ino][jno][idim][jdim]*zry.z;
			emat_dd[ino][jno][idim][jdim] -= dtmp1+dtmp2;
		}
		if( ino*ndim+idim != jno*ndim+jdim ){
			emat_dd[jno][ino][jdim][idim] = emat_dd[ino][jno][idim][jdim];
		}
	}
	}
	}
	}

	////////////////////////////////
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
	for(unsigned int idim=0;idim<ndim;idim++){
	for(unsigned int jdim=0;jdim<ndim;jdim++){
		for(unsigned int kno=0;kno<nno;kno++){
			emat_dt[ino][jno][idim][jdim] 
				= dodu[kno][0][ino][idim]*loc_base2[jno][0][jdim]*emat_rr[kno][jno][0][0]
				+ dodu[kno][0][ino][idim]*loc_base2[jno][1][jdim]*emat_rr[kno][jno][0][1]
				+ dodu[kno][1][ino][idim]*loc_base2[jno][0][jdim]*emat_rr[kno][jno][1][0]
				+ dodu[kno][1][ino][idim]*loc_base2[jno][1][jdim]*emat_rr[kno][jno][1][1];
		}
		{
			const double dtmp1 = ddodudt[jno][0][ino][idim][jdim]*res_rot2d[jno][0];
			const double dtmp2 = ddodudt[jno][1][ino][idim][jdim]*res_rot2d[jno][1];
			emat_dt[ino][jno][idim][jdim] -= dtmp1+dtmp2;
		}
	}
	}
	}
	}
	
	////////////////////////////////
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
	for(unsigned int idim=0;idim<ndim;idim++){
	for(unsigned int jdim=0;jdim<ndim;jdim++){
		emat_td[ino][jno][idim][jdim] = emat_dt[jno][ino][jdim][idim];
	}
	}
	}
	}
	
	////////////////////////////////
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int jno=0;jno<nno;jno++){
	for(unsigned int idim=0;idim<ndim;idim++){
	for(unsigned int jdim=0;jdim<ndim;jdim++){
		emat_tt[ino][jno][idim][jdim] 
			= loc_base2[ino][0][idim]*loc_base2[jno][0][jdim]*emat_rr[ino][jno][0][0]
			+ loc_base2[ino][0][idim]*loc_base2[jno][1][jdim]*emat_rr[ino][jno][0][1]
			+ loc_base2[ino][1][idim]*loc_base2[jno][0][jdim]*emat_rr[ino][jno][1][0]
			+ loc_base2[ino][1][idim]*loc_base2[jno][1][jdim]*emat_rr[ino][jno][1][1];
	}	
	}	
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int idim=0;idim<ndim;idim++){
		emat_tt[ino][ino][idim][idim] += loc_base0[2][idim]*loc_base0[2][idim]*torsion_stiff;
	}
	}
	for(unsigned int ino=0;ino<nno;ino++){
//		const Com::CVector3D ax3d0(loc_base0[0][0],loc_base0[0][1],loc_base0[0][2]);
//		const Com::CVector3D ay3d0(loc_base0[1][0],loc_base0[1][1],loc_base0[1][2]);
		const Com::CVector3D az3d0(loc_base0[2][0],loc_base0[2][1],loc_base0[2][2]);
		const Com::CVector3D rot3di(rot3d[ino][0],rot3d[ino][1],rot3d[ino][2]);
//		const Com::CVector3D ax3d1 = Com::RotateVector(ax3d0,rot3di);
//		const Com::CVector3D ay3d1 = Com::RotateVector(ay3d0,rot3di);
		const Com::CVector3D az3d1 = Com::RotateVector(az3d0,rot3di);
		const Com::CVector3D rb1x(loc_base1[0][0],loc_base1[0][1],loc_base1[0][2]);
		const Com::CVector3D rb1y(loc_base1[1][0],loc_base1[1][1],loc_base1[1][2]);
		Com::CVector3D yrxxry = (-res_rot2d[ino][0])*rb1y+(res_rot2d[ino][1])*rb1x;
//		const double dtmp0 = Com::Dot(yrxxry,ax3d1);
//		const double dtmp1 = Com::Dot(yrxxry,ay3d1);
		const double dtmp2 = Com::Dot(yrxxry,az3d1);
		if( dtmp2 > 0 ){
			for(unsigned int idim=0;idim<ndim;idim++){
				const double dti0 = loc_base0[0][idim]; 
				const double dti1 = loc_base0[1][idim];
				const double val = -dtmp2*(dti0*dti0+dti1*dti1);
				emat_tt[ino][ino][idim][idim] -= val;
			}
		}
	}

	{
		double rtmp0=0, rtmp1=0, rtmp2=0;
		rtmp0 = ( res_disp2d[1][0] < 0 ) ? res_disp2d[1][0] : 0;
		rtmp1 = ( res_disp2d[2][0] < 0 ) ? res_disp2d[2][0] : 0;
//		rtmp1 =   res_disp2d[2][0];
		rtmp2 = ( res_disp2d[2][1] < 0 ) ? res_disp2d[2][1] : 0;
//		rtmp0 = res_disp2d[1][0];
//		rtmp1 = res_disp2d[2][0];
//		rtmp2 = res_disp2d[2][1];
/*			
		const double th0 = 0.000000000001;
		if( rtmp0 > 0 ){ rtmp0 = 0; }
		else if( rtmp0 > -th0 ){ rtmp0 = -1/(th0*th0) * rtmp0*rtmp0*(rtmp0+2*th0); }
		if( rtmp2 > 0 ){ rtmp2 = 0; }
		else if( rtmp2 > -th0 ){ rtmp2 = -1/(th0*th0) * rtmp2*rtmp2*(rtmp2+2*th0); }*/

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
			for(unsigned int idim=0;idim<ndim;idim++){
			for(unsigned int jdim=0;jdim<ndim;jdim++){
				emat_dd[ino][jno][idim][jdim] -= 
					 rtmp0*dd2ud3ud3u[0][ino][jno][idim][jdim]
					+rtmp1*dd2ud3ud3u[1][ino][jno][idim][jdim]
					+rtmp2*dd2ud3ud3u[2][ino][jno][idim][jdim];
			}
			}
		}
		}
	}

//	std::cout << res_disp2d[1][0] << " " << res_disp2d[2][0] << " " << res_disp2d[2][1] << std::endl;

	////////////////
//	double eres_d[nno][ndim];
	{
		double dtmp1 = area0/3.0*arearho;
		for(unsigned int ino=0;ino<nno;ino++){
			eres_d[ino][0] = dtmp1*g_x;
			eres_d[ino][1] = dtmp1*g_y;
			eres_d[ino][2] = dtmp1*g_z;
		}
		const double dtmpx = loc_base1[2][0]*area0*press/3.0;
		const double dtmpy = loc_base1[2][1]*area0*press/3.0;
		const double dtmpz = loc_base1[2][2]*area0*press/3.0;
		for(unsigned int ino=0;ino<nno;ino++){
			eres_d[ino][0] += dtmpx;
			eres_d[ino][1] += dtmpy;
			eres_d[ino][2] += dtmpz;
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int idim=0;idim<ndim;idim++){
			eres_d[ino][idim] += d2ud3u[0][ino][idim]*res_disp2d[1][0]
							   + d2ud3u[1][ino][idim]*res_disp2d[2][0]
							   + d2ud3u[2][ino][idim]*res_disp2d[2][1];
			double dtmp0=0, dtmp1=0;
			for(unsigned int jno=0;jno<nno;jno++){
				dtmp0 += dodu[jno][0][ino][idim]*res_rot2d[jno][0];
				dtmp1 += dodu[jno][1][ino][idim]*res_rot2d[jno][1];
			}
			eres_d[ino][idim] += dtmp0 + dtmp1;
		}
		}
	}

//	double eres_t[nno][ndim];
	for(unsigned int ino=0;ino<nno;ino++){
	for(unsigned int idim=0;idim<ndim;idim++){
		eres_t[ino][idim] = loc_base2[ino][0][idim]*res_rot2d[ ino][0]
			              + loc_base2[ino][1][idim]*res_rot2d[ ino][1]
						  + loc_base0[2][idim]*res_arot2d[ino];
	}
	}
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



static bool AddLinearSystem_DKT3D_NonLinear_Static_P1(
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta, const CFieldWorld& world,
		const unsigned int id_ea )
{
//	std::cout << "DKT 3dim NonLinear" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,world);
	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,world);
	CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp,CORNER, id_field_theta,CORNER, world);
	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER, id_field_disp,CORNER, world);
	CVector_Blk& res_d = ls.GetResidual(id_field_disp, CORNER,world);
	CVector_Blk& res_t = ls.GetResidual(id_field_theta,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);	
		double cord3d0[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],cord3d0[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue(no[ino],disp3d[ino]); }
		double rot3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue(no[ino],rot3d[ino]); }


		double emat_dd[3][3][3][3], emat_dt[3][3][3][3], emat_td[3][3][3][3], emat_tt[3][3][3][3];
		double eres_d[3][3], eres_t[3][3];

		double area0;
		MakeElemStiffMatDKT3D_Nonlinear(
			young, poisson, thickness, arearho, 
			g_x, g_y, g_z, press,
			cord3d0, disp3d, rot3d,
			area0,
			emat_dd, emat_dt, emat_td, emat_tt, 
			eres_d, eres_t );

		mat_dd.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
    	}
		for(unsigned int ino=0;ino<nno;ino++){
			res_t.AddValue(no[ino],0,eres_t[ino][0]);
			res_t.AddValue(no[ino],1,eres_t[ino][1]);
			res_t.AddValue(no[ino],2,eres_t[ino][2]);
		}
	}
	return true;
}


static bool AddLinearSystem_DKT3D_NonLinear_Static_P1_Combined(
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta, const CFieldWorld& world,
		const unsigned int id_ea )
{
//	std::cout << "DKT 3dim NonLinear" << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp,CORNER,world);
	CVector_Blk& res_d = ls.GetResidual(id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);	
		double cord3d0[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],cord3d0[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue(no[ino],disp3d[ino]); }
		double rot3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue(no[ino],rot3d[ino]); }


		double emat_dd[3][3][3][3], emat_dt[3][3][3][3], emat_td[3][3][3][3], emat_tt[3][3][3][3];
		double eres_d[3][3], eres_t[3][3];

		double area0;
		MakeElemStiffMatDKT3D_Nonlinear(
			young, poisson, thickness, arearho, 
			g_x, g_y, g_z, press,
			cord3d0, disp3d, rot3d,
			area0,
			emat_dd, emat_dt, emat_td, emat_tt, 
			eres_d, eres_t );

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = emat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = emat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = emat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = emat_tt[ino][jno][idim][jdim];
		}
		}
		}
		}
		mat_dd.Mearge(nno,no,nno,no,6*6,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
			res_d.AddValue(no[ino],3,eres_t[ino][0]);
			res_d.AddValue(no[ino],4,eres_t[ino][1]);
			res_d.AddValue(no[ino],5,eres_t[ino][2]);
		}
	}

	return true;
}


bool Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_Static(
	Fem::Ls::CLinearSystem_Field& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z,  double press, 
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != VECTOR3 ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_deflect,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_rot,CORNER,world) ){
			    return AddLinearSystem_DKT3D_NonLinear_Static_P1_Combined(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,
				    id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_NonLinear_Static_P1(
				    ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,
				    id_ea);
            }
		}
		assert(0);
        return false;
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_Static(
					ls,
					young, poisson, thickness, arearho,
					g_x, g_y, g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////






static bool AddLinearSystem_DKT3D_NonLinear_NonStatic_P1(
		double dt, double gamma_newmark, double beta_newmark, bool is_first_itr,
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
	std::cout << "DKT3D NonLinear NonStatic : " << dt << " " << gamma_newmark << " " << beta_newmark << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  != ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp, CORNER,world);
	CMatDia_BlkCrs& mat_tt = ls.GetMatrix(id_field_theta,CORNER,world);
	CMat_BlkCrs& mat_dt = ls.GetMatrix(id_field_disp,CORNER, id_field_theta,CORNER, world);
	CMat_BlkCrs& mat_td = ls.GetMatrix(id_field_theta,CORNER, id_field_disp,CORNER, world);
	CVector_Blk& res_d = ls.GetResidual(id_field_disp, CORNER,world); 
	CVector_Blk& res_t = ls.GetResidual(id_field_theta,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);	
		double cord3d0[nno][ndim];
        for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],cord3d0[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp3d[nno][ndim];
		double disp_acc[nno][ndim];
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ 
            ns_c_d .GetValue(no[ino],disp3d[   ino]); 
		    ns_c_da.GetValue(no[ino],disp_acc[ ino]); 
		    ns_c_dv.GetValue(no[ino],disp_velo[ino]); 
        }
		double rot3d[   nno][ndim];
		double rot_velo[nno][ndim];
		double rot_acc[ nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ 
		    ns_c_r .GetValue(no[ino],rot3d[   ino]); 
		    ns_c_rv.GetValue(no[ino],rot_velo[ino]); 
		    ns_c_ra.GetValue(no[ino],rot_acc[ ino]); 
        }

		double eKmat_dd[3][3][3][3], eKmat_dt[3][3][3][3], eKmat_td[3][3][3][3], eKmat_tt[3][3][3][3];
		double eres_d[3][3], eres_t[3][3];

		double area0;
		MakeElemStiffMatDKT3D_Nonlinear(
			young, poisson, thickness, arearho,
			g_x, g_y, g_z, press,
			cord3d0, disp3d, rot3d,
			area0,
			eKmat_dd, eKmat_dt, eKmat_td, eKmat_tt, 
			eres_d, eres_t );

		double eMmat_dd[nno][nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area0/12.0;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat_dd[ino][jno][0][0] += dtmp0;
					eMmat_dd[ino][jno][1][1] += dtmp0;
					eMmat_dd[ino][jno][2][2] += dtmp0;
				}
				eMmat_dd[ino][ino][0][0] += dtmp0;
				eMmat_dd[ino][ino][1][1] += dtmp0;
				eMmat_dd[ino][ino][2][2] += dtmp0;
			}
		}

		double emat_dd[3][3][3][3], emat_dt[3][3][3][3], emat_td[3][3][3][3], emat_tt[3][3][3][3];
		{
			double dtmp1 = beta_newmark*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_dd[0][0][0][0])[i] = (&eMmat_dd[0][0][0][0])[i]+dtmp1*(&eKmat_dd[0][0][0][0])[i];
				(&emat_dt[0][0][0][0])[i] = dtmp1*(&eKmat_dt[0][0][0][0])[i];
				(&emat_td[0][0][0][0])[i] = dtmp1*(&eKmat_td[0][0][0][0])[i];
				(&emat_tt[0][0][0][0])[i] = dtmp1*(&eKmat_tt[0][0][0][0])[i];
			}
		}

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= eMmat_dd[ino][jno][idim][jdim]*disp_acc[jno][jdim]; 
		}
		}
		}
		}
        if( is_first_itr ){
		    for(unsigned int ino=0;ino<nno;ino++){
		    for(unsigned int jno=0;jno<nno;jno++){
		    for(unsigned int idim=0;idim<ndim;idim++){
		    for(unsigned int jdim=0;jdim<ndim;jdim++){
			    eres_d[ino][idim] -= dt*(       eKmat_dd[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							                   +eKmat_dt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			    eres_d[ino][idim] -= 0.5*dt*dt*(eKmat_dd[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                                   +eKmat_dt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
			    ////////////////////////////////
			    eres_t[ino][idim] -= dt*(       eKmat_td[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							                   +eKmat_tt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			    eres_t[ino][idim] -= 0.5*dt*dt*(eKmat_td[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                                   +eKmat_tt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
				                             
		    }
		    }
		    }
		    }
        }

		mat_dd.Mearge(nno,no,nno,no,ndim*ndim,&emat_dd[0][0][0][0]);
		mat_dt.Mearge(nno,no,nno,no,ndim*ndim,&emat_dt[0][0][0][0]);
		mat_td.Mearge(nno,no,nno,no,ndim*ndim,&emat_td[0][0][0][0]);
		mat_tt.Mearge(nno,no,nno,no,ndim*ndim,&emat_tt[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
		}
		for(unsigned int ino=0;ino<nno;ino++){
			res_t.AddValue(no[ino],0,eres_t[ino][0]);
			res_t.AddValue(no[ino],1,eres_t[ino][1]);
			res_t.AddValue(no[ino],2,eres_t[ino][2]);
		}
	}
	return true;
}



static bool AddLinearSystem_DKT3D_NonLinear_NonStatic_P1_Combined(
		double dt, double gamma_newmark, double beta_newmark, bool is_first_itr, 
		CLinearSystem_Field& ls, 
		double young, double poisson, double thickness, double arearho,
		double g_x, double g_y, double g_z, double press,
		const unsigned int id_field_disp, const unsigned int id_field_theta,
		const CFieldWorld& world,
		const unsigned int id_ea )
{
    std::cout << "DKT3D NonLinear NonStatic(Combined)  dt:" << dt << std::endl;

	assert( world.IsIdEA(id_ea) );
	const CElemAry& ea = world.GetEA(id_ea);
	assert( ea.ElemType() == TRI );

	if( !world.IsIdField(id_field_disp) ) return false;
	const CField& field_disp = world.GetField(id_field_disp);

	if( !world.IsIdField(id_field_theta) ) return false;
	const CField& field_rot = world.GetField(id_field_theta);

	const CElemAry::CElemSeg& es_c_va = field_disp.GetElemSeg(id_ea,CORNER,true, world);
	const CElemAry::CElemSeg& es_c_co = field_disp.GetElemSeg(id_ea,CORNER,false,world);

	const unsigned int nno = 3;
	const unsigned int ndim = 3;

	assert(  ls.FindIndexArray_Seg(id_field_disp, CORNER,world) 
		  == ls.FindIndexArray_Seg(id_field_theta,CORNER,world) );

	CMatDia_BlkCrs& mat_dd = ls.GetMatrix(id_field_disp,CORNER,world);
	CVector_Blk& res_d = ls.GetResidual(  id_field_disp,CORNER,world);

	const CNodeAry::CNodeSeg& ns_c_d = field_disp.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_dv = field_disp.GetNodeSeg(CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_da = field_disp.GetNodeSeg(CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_r = field_rot.GetNodeSeg(CORNER,true,world);
	const CNodeAry::CNodeSeg& ns_c_rv = field_rot.GetNodeSeg( CORNER,true,world,VELOCITY);
	const CNodeAry::CNodeSeg& ns_c_ra = field_rot.GetNodeSeg( CORNER,true,world,ACCELERATION);
	const CNodeAry::CNodeSeg& ns_c_co  = field_disp.GetNodeSeg(CORNER,false,world);

	for(unsigned int ielem=0;ielem<ea.Size();ielem++)
	{
		unsigned int no[nno];
		es_c_co.GetNodes(ielem,no);	
		double cord3d0[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_co.GetValue(no[ino],cord3d0[ino]); }
		es_c_va.GetNodes(ielem,no);
		double disp3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_d.GetValue(no[ino],disp3d[ino]); }
		double disp_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_da.GetValue( no[ino],disp_acc[ ino]); }
		double disp_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_dv.GetValue( no[ino],disp_velo[ ino]); }
		double rot3d[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_r.GetValue( no[ino],rot3d[  ino]); }
		double rot_velo[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_rv.GetValue( no[ino],rot_velo[ ino]); }
		double rot_acc[nno][ndim];
		for(unsigned int ino=0;ino<nno;ino++){ ns_c_ra.GetValue( no[ino],rot_acc[ ino]); }

		double eKmat_dd[3][3][3][3], eKmat_dt[3][3][3][3], eKmat_td[3][3][3][3], eKmat_tt[3][3][3][3];
		double eres_d[3][3], eres_t[3][3];

		double area0;
		MakeElemStiffMatDKT3D_Nonlinear(
			young, poisson, thickness, arearho,
			g_x, g_y, g_z, press,
			cord3d0, disp3d, rot3d,
			area0,
			eKmat_dd, eKmat_dt, eKmat_td, eKmat_tt, 
			eres_d, eres_t );

		double eMmat_dd[nno][nno][ndim][ndim];
		{
			const double dtmp0 = arearho*area0/12.0;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){ *(&eMmat_dd[0][0][0][0]+i) = 0; }
			for(unsigned int ino=0;ino<nno;ino++){
				for(unsigned int jno=0;jno<nno;jno++){
					eMmat_dd[ino][jno][0][0] += dtmp0;
					eMmat_dd[ino][jno][1][1] += dtmp0;
					eMmat_dd[ino][jno][2][2] += dtmp0;
				}
				eMmat_dd[ino][ino][0][0] += dtmp0;
				eMmat_dd[ino][ino][1][1] += dtmp0;
				eMmat_dd[ino][ino][2][2] += dtmp0;
			}
		}

		double emat_dd[3][3][3][3], emat_dt[3][3][3][3], emat_td[3][3][3][3], emat_tt[3][3][3][3];
		{
			double dtmp1 = beta_newmark*dt*dt;
			for(unsigned int i=0;i<nno*nno*ndim*ndim;i++){
				(&emat_dd[0][0][0][0])[i] = (&eMmat_dd[0][0][0][0])[i]+dtmp1*(&eKmat_dd[0][0][0][0])[i];
				(&emat_dt[0][0][0][0])[i] = dtmp1*(&eKmat_dt[0][0][0][0])[i];
				(&emat_td[0][0][0][0])[i] = dtmp1*(&eKmat_td[0][0][0][0])[i];
				(&emat_tt[0][0][0][0])[i] = dtmp1*(&eKmat_tt[0][0][0][0])[i];
			}
		}

		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= eMmat_dd[ino][jno][idim][jdim]*disp_acc[jno][jdim]; 
		}
		}
		}
		}
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			eres_d[ino][idim] -= dt*(       eKmat_dd[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_dt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_d[ino][idim] -= 0.5*dt*dt*(eKmat_dd[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_dt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
			////////////////////////////////
			eres_t[ino][idim] -= dt*(       eKmat_td[ino][jno][idim][jdim]*disp_velo[jno][jdim]
							               +eKmat_tt[ino][jno][idim][jdim]*rot_velo[ jno][jdim]);
			eres_t[ino][idim] -= 0.5*dt*dt*(eKmat_td[ino][jno][idim][jdim]*disp_acc[ jno][jdim]
			                               +eKmat_tt[ino][jno][idim][jdim]*rot_acc[  jno][jdim]);
				                         
		}
		}
		}
		}

		double emat[nno][nno][6][6];
		for(unsigned int ino=0;ino<nno;ino++){
		for(unsigned int jno=0;jno<nno;jno++){
		for(unsigned int idim=0;idim<ndim;idim++){
		for(unsigned int jdim=0;jdim<ndim;jdim++){
			emat[ino][jno][idim     ][jdim     ] = emat_dd[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim     ] = emat_td[ino][jno][idim][jdim];
			emat[ino][jno][idim     ][jdim+ndim] = emat_dt[ino][jno][idim][jdim];
			emat[ino][jno][idim+ndim][jdim+ndim] = emat_tt[ino][jno][idim][jdim];
		}
		}
		}
		}
		mat_dd.Mearge(nno,no,nno,no,6*6,&emat[0][0][0][0]);
		for(unsigned int ino=0;ino<nno;ino++){
			res_d.AddValue(no[ino],0,eres_d[ino][0]);
			res_d.AddValue(no[ino],1,eres_d[ino][1]);
			res_d.AddValue(no[ino],2,eres_d[ino][2]);
			res_d.AddValue(no[ino],3,eres_t[ino][0]);
			res_d.AddValue(no[ino],4,eres_t[ino][1]);
			res_d.AddValue(no[ino],5,eres_t[ino][2]);
		}
	}
	return true;
}



bool Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_NonStatic(
	double dt, double gamma_newmark, double beta_newmark, bool is_first_itr,
	Fem::Ls::CLinearSystem_Field& ls,
	double young, double poisson, double thickness, double arearho,
	double g_x, double g_y, double g_z, double press,
	const Fem::Field::CFieldWorld& world,
	const unsigned int id_field_deflect, unsigned int id_field_rot, 
	unsigned int id_ea )
{
	if( !world.IsIdField(id_field_deflect) ) return false;
	const CField& field_deflect = world.GetField(id_field_deflect);
	if( field_deflect.GetFieldType() != VECTOR3 ) return false;

	if( id_ea != 0 ){
		if( field_deflect.GetInterpolationType(id_ea,world) == TRI11 ){
            if( ls.FindIndexArray_Seg(id_field_deflect,CORNER,world) 
                == ls.FindIndexArray_Seg(id_field_rot,CORNER,world) ){
			    return AddLinearSystem_DKT3D_NonLinear_NonStatic_P1_Combined(
				    dt, gamma_newmark, beta_newmark, is_first_itr, ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,
				    id_ea);
            }
            else{
			    return AddLinearSystem_DKT3D_NonLinear_NonStatic_P1(
				    dt, gamma_newmark, beta_newmark, is_first_itr, ls,
				    young, poisson, thickness, arearho,
				    g_x, g_y, g_z, press,
				    id_field_deflect,id_field_rot,world,
				    id_ea);
            }
		}
		else{
			assert(0);
		}
	}
	else{
		const std::vector<unsigned int> aIdEA = field_deflect.GetAryIdEA();
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			const unsigned int id_ea = aIdEA[iiea];
			bool res = Fem::Eqn::AddLinearSystem_DKT3D_NonLinear_NonStatic(
					dt, gamma_newmark, beta_newmark, is_first_itr, ls,
					young, poisson, thickness, arearho,   g_x,g_y,g_z, press,
					world,
					id_field_deflect, id_field_rot,
					id_ea );
			if( !res ) return false;
		}
		return true;
	}

	return true;
}
