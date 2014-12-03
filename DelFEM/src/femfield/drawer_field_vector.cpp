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
// drawer_field_vector.cpp : implementation of class (CDrawer_Vector)
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
    #pragma warning ( disable : 4786 )
#endif

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#include <assert.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <memory>

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif


#include "delfem/drawer_field_vector.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

CDrawerVector::CDrawerVector(){
	this->pData = 0;
	this->npo = 0;
}

CDrawerVector::CDrawerVector(unsigned int id_field, const Fem::Field::CFieldWorld& world){
	this->pData = 0;
	this->npo = 0;
	this->Set(id_field, world);
}

CDrawerVector::~CDrawerVector(){
	if( pData != 0 ){ delete pData; }
}

Com::CBoundingBox3D CDrawerVector::GetBoundingBox( double rot[] ) const{
//	if( this->pData == 0 ) return Com::CBoundingBox();
//	return m_paVer->GetBoundingBox(rot);
	return Com::CBoundingBox3D();
}

// update value of vector
bool CDrawerVector::Update_VECTOR(const Fem::Field::CFieldWorld& world)
{
	assert( world.IsIdField(id_field) );
	if( !world.IsIdField(id_field) ) return false;
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "Error!->Not Implemented" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	assert( field.GetFieldType() == VECTOR2 || field.GetFieldType() == VECTOR3 );
	npo = 0;
	{
		// add number of node in corner
		unsigned int id_na_val_c = field.GetNodeSegInNodeAry(CORNER).id_na_va;
		if( id_na_val_c !=  0 ){
			assert( world.IsIdNA(id_na_val_c) );
			const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_c);
			npo += na_val.Size();
		}
		// add number of node in bubble
		unsigned int id_na_val_b = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		if( id_na_val_b != 0 ){
			assert( world.IsIdNA(id_na_val_b) );
			const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_b);
			npo += na_val.Size();
		}
	}

	int ilayer_min, ilayer_max;
	{
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		if( aIdEA.size() > 0 ){
			ilayer_min = field.GetLayer(aIdEA[0]);
			ilayer_max = ilayer_min;
		}
		else{ ilayer_min=0; ilayer_max=0; }
		for(unsigned int iiea=1;iiea<aIdEA.size();iiea++){
			int ilayer = field.GetLayer(aIdEA[iiea]);
			ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
			ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
		}
	}

//	std::cout << ilayer_min << " " << ilayer_max << std::endl;

	unsigned int ndim_co0 = field.GetNDimCoord();
	if( pData == 0 ){
		if( ilayer_min == ilayer_max ){
			ndim_co = ndim_co0;
			ndim_va = ndim_co;
		}
		else{
			assert( ndim_co0 == 2 );
			ndim_co = 3;
			ndim_va = 2;
		}
		pData = new double [(ndim_co+ndim_va)*npo];
	}
	else{
		if( ilayer_min == ilayer_max ){ 
			assert( ndim_co == ndim_co0 ); 
			assert( ndim_va == ndim_co );
		}
		else{ 
			assert( ndim_co == 3 && ndim_va == 2 ); 			
		}
	}
  
//  std::cout << ndim_co0 << " " << ndim_co << std::endl;

	unsigned int icoun = 0;
	// Make Data of Coord
	if( field.GetNodeSegInNodeAry(CORNER).id_na_va != 0 ){
		const CField::CNodeSegInNodeAry& nsna_c = field.GetNodeSegInNodeAry(CORNER);
		assert( world.IsIdNA(nsna_c.id_na_va) );
		const Fem::Field::CNodeAry& na_c_val = world.GetNA(nsna_c.id_na_va);
		const unsigned int npoin_va = na_c_val.Size();
		unsigned int id_ns_c_v;
		{
			if(      nsna_c.id_ns_va != 0 ) id_ns_c_v = nsna_c.id_ns_va;
			else if( nsna_c.id_ns_ve != 0 ) id_ns_c_v = nsna_c.id_ns_ve;
			else if( nsna_c.id_ns_ac != 0 ) id_ns_c_v = nsna_c.id_ns_ac;
			else{ assert(0); }
		}
		assert( na_c_val.IsSegID(id_ns_c_v) );
		const Fem::Field::CNodeAry::CNodeSeg& ns_c_val = na_c_val.GetSeg(id_ns_c_v);
		assert( world.IsIdNA(nsna_c.id_na_co) );
		const Fem::Field::CNodeAry& na_c_co = world.GetNA(nsna_c.id_na_co);
		const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(nsna_c.id_ns_co);
		double coord[3],value[3];
		for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
			unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
			ns_c_co.GetValue(ipoin_co,coord);
			ns_c_val.GetValue(ipoin,value);
			for(unsigned int idim=0;idim<ndim_co0;idim++){	// 
				pData[ipoin*(ndim_co+ndim_va)+idim] = coord[idim];
			}
			for(unsigned int idim=0;idim<ndim_va;idim++){
				pData[ipoin*(ndim_co+ndim_va)+ndim_co+idim] = value[idim];
			}
			if( ilayer_min != ilayer_max ){
				assert( ndim_co0 == 2 );
				assert( ndim_va  == 2 );
				assert( ndim_co  == 3 );
				pData[ipoin*5+2] = 0.01;
			}
		}
		if( ilayer_min != ilayer_max ){ // consider layer
			const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
			for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
				const unsigned int id_ea = aIdEA[iiea];
				const int ilayer = field.GetLayer(id_ea);
				const double height = (ilayer+0.01-ilayer_min)/(1+ilayer_max-ilayer_min);
				const CElemAry::CElemSeg& es = field.GetElemSeg(id_ea,CORNER,true,world);
				const unsigned int nnoes = es.Length();
				assert( nnoes < 16 );
				unsigned int noes[16];
				for(unsigned int ielem=0;ielem<es.Size();ielem++){
					es.GetNodes(ielem,noes);
					for(unsigned int inoes=0;inoes<nnoes;inoes++){
						const unsigned int ipo0 = noes[inoes];
						pData[ipo0*5+2] = height;
					}
				}
			}
		}
		icoun = npoin_va;
	}
	if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 ){
		const CField::CNodeSegInNodeAry& nsna_b = field.GetNodeSegInNodeAry(BUBBLE);
		assert( world.IsIdNA(nsna_b.id_na_va) );
		const Fem::Field::CNodeAry& na_val = world.GetNA(nsna_b.id_na_va);
		const unsigned int npoin_va = na_val.Size();
		unsigned int id_ns_v;
		{
			if(      nsna_b.id_ns_va != 0 ) id_ns_v = nsna_b.id_ns_va;
			else if( nsna_b.id_ns_ve != 0 ) id_ns_v = nsna_b.id_ns_ve;
			else if( nsna_b.id_ns_ac != 0 ) id_ns_v = nsna_b.id_ns_ac;
			else{ assert(0); }
		}
		assert( na_val.IsSegID(id_ns_v) );
		const Fem::Field::CNodeAry::CNodeSeg& ns_val = na_val.GetSeg(id_ns_v);
		assert( ndim_co == ns_val.Length() );
		if( world.IsIdNA(nsna_b.id_na_co) ){
//			std::cout << "bubble with coord" << std::endl;
			const Fem::Field::CNodeAry& na_co = world.GetNA(nsna_b.id_na_co);
			const CNodeAry::CNodeSeg& ns_co = na_co.GetSeg(nsna_b.id_ns_co);
			double coord[3],value[3];
			for(unsigned int ipoin=0;ipoin<npoin_va;ipoin++){
				unsigned int ipoin_co = field.GetMapVal2Co(ipoin);
				ns_co.GetValue(ipoin_co,coord);
				ns_val.GetValue(ipoin,value);
				for(unsigned int idim=0;idim<ndim_co0;idim++){	// coord
					pData[(icoun+ipoin)*(ndim_co+ndim_va)+idim] = coord[idim];
				}
				for(unsigned int idim=0;idim<ndim_va;idim++){	// value
					pData[(icoun+ipoin)*(ndim_co+ndim_va)+ndim_va+idim] = value[idim];
				}
			}
		}
		else{
//			std::cout << "bubble with gravity center" << std::endl;
			unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
			unsigned int id_ns_c_co = field.GetNodeSegInNodeAry(CORNER).id_ns_co;
			assert( world.IsIdNA(id_na_c_co) );
			const CNodeAry& na_c_co = world.GetNA(id_na_c_co);
			assert( na_c_co.IsSegID(id_ns_c_co) );
			const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);
			assert( world.IsIdNA(nsna_b.id_na_va) );
			const CNodeAry& na_va = world.GetNA(nsna_b.id_na_va);
			assert( na_va.IsSegID(id_ns_v) );
			unsigned int noes[64];
			const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
			for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
				unsigned int id_ea = aIdEA[iiea];
				const CElemAry& ea = world.GetEA(id_ea);
				const CElemAry::CElemSeg& es_c_co = field.GetElemSeg(id_ea,CORNER,false,world);
				assert( es_c_co.GetIdNA() == id_na_c_co );
				const CElemAry::CElemSeg& es_b_va = field.GetElemSeg(id_ea,BUBBLE,true,world);
				assert( es_b_va.GetIdNA() == nsna_b.id_na_va );
        const unsigned int nnoes = es_c_co.Length();        
        const double invnnoes = 1.0/(double)nnoes;
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					double coord_cnt[3];
					{	// gravity center
						for(unsigned int idim=0;idim<ndim_co0;idim++){ coord_cnt[idim] = 0.0; }
						es_c_co.GetNodes(ielem,noes);
						double coord[3];
						for(unsigned int inoes=0;inoes<nnoes;inoes++){
							unsigned int ipoi0 = noes[inoes]; assert( ipoi0 < na_c_co.Size() );
							ns_c_co.GetValue(ipoi0,coord);
							for(unsigned int idim=0;idim<ndim_co0;idim++){ coord_cnt[idim] += coord[idim]; }
						}
						for(unsigned int idim=0;idim<ndim_co0;idim++){ coord_cnt[idim] *= invnnoes; }
//            std::cout << coord_cnt[0] << " " << coord_cnt[1] << " " << coord_cnt[2] << std::endl;
					}
					double value[3];
					{	// get value
						es_b_va.GetNodes(ielem,noes);
						unsigned int ipoi0 = noes[0];
						assert( ipoi0 < na_va.Size() );
						ns_val.GetValue(ipoi0,value);
					}
					for(unsigned int idim=0;idim<ndim_co0;idim++){	// coord
						pData[icoun*(ndim_co+ndim_va)+idim] = coord_cnt[idim];
					}
					for(unsigned int idim=0;idim<ndim_va;idim++){	// value
						pData[icoun*(ndim_co+ndim_va)+ndim_co+idim] = value[idim];
					}
					icoun++;
				} // end ielem
			} // end iei
		} // end if
	}
	return true;
}

void GetPrincipleVector_STSR2(const double sstr[3], 
							  double pvecs[2],double& ls, 
							  double pvecl[2],double& ll )
{
	{
		const double tmp1 = sqrt( (sstr[0]-sstr[1])*(sstr[0]-sstr[1])+4*sstr[2]*sstr[2] );
		const double tmp2 = sstr[0]+sstr[1];
		const double l1 = 0.5*(tmp2-tmp1);
		const double l2 = 0.5*(tmp2+tmp1);
		if( fabs(l1) > fabs(l2) ){ ll=l1; ls=l2; }
		else{                      ll=l2; ls=l1; }
	}
	{
		const double candl1[2] = { -sstr[2], sstr[0]-ls };
		const double candl2[2] = { sstr[1]-ls, -sstr[2] };
		const double sqlen1 = candl1[0]*candl1[0]+candl1[1]*candl1[1];
		const double sqlen2 = candl2[0]*candl2[0]+candl2[1]*candl2[1];
		if( sqlen1 > sqlen2 ){	pvecs[0] = candl1[0];	pvecs[1] = candl1[1];	}
		else{					pvecs[0] = candl2[0];	pvecs[1] = candl2[1];	}
		const double sqlen = sqrt(pvecs[0]*pvecs[0]+pvecs[1]*pvecs[1]);
		if( sqlen < 1.0e-10 ){
			pvecs[0] = 0;	pvecs[1] = 0;
		}
		else{
			const double normalizer = ls/sqlen;
			pvecs[0] *= normalizer;
			pvecs[1] *= normalizer;
		}
	}
	{
		const double candl1[2] = { -sstr[2], sstr[0]-ll };
		const double candl2[2] = { sstr[1]-ll, -sstr[2] };
		const double sqlen1 = candl1[0]*candl1[0]+candl1[1]*candl1[1];
		const double sqlen2 = candl2[0]*candl2[0]+candl2[1]*candl2[1];
		if( sqlen1 > sqlen2 ){	pvecl[0] = candl1[0];	pvecl[1] = candl1[1];	}
		else{					pvecl[0] = candl2[0];	pvecl[1] = candl2[1];	}
		const double sqlen = sqrt(pvecl[0]*pvecl[0]+pvecl[1]*pvecl[1]);
		if( sqlen < 1.0e-10 ){
			pvecl[0] = 0;	pvecl[1] = 0;
		}
		else{
			const double normalizer = ll/sqlen;
			pvecl[0] *= normalizer;
			pvecl[1] *= normalizer;
		}
	}
}

// visualize symetric tensor with principle vector
bool CDrawerVector::Update_SSTR2(const Fem::Field::CFieldWorld& world)
{
//	std::cout << "Update_SSTR2" << std::endl;
	assert( world.IsIdField(id_field) );
	if( !world.IsIdField(id_field) ) return false;
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "Error!-->Not Implemented" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	assert( field.GetFieldType() == STSR2 );
	if( field.GetNodeSegInNodeAry(BUBBLE).id_na_va == 0 ) return false;
	npo= 0;
	{	// get number of line (BUBBLE)
		unsigned int id_na_val_b = field.GetNodeSegInNodeAry(BUBBLE).id_na_va;
		assert( world.IsIdNA(id_na_val_b) );
		const Fem::Field::CNodeAry& na_val = world.GetNA(id_na_val_b);
		npo = na_val.Size();
	}
	int ilayer_min, ilayer_max;
	{
		const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
		if( aIdEA.size() > 0 ){
			ilayer_min = field.GetLayer(aIdEA[0]);
			ilayer_max = ilayer_min;
		}
		else{ ilayer_min=0; ilayer_max=0; }
		for(unsigned int iiea=1;iiea<aIdEA.size();iiea++){
			int ilayer = field.GetLayer(aIdEA[iiea]);
			ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
			ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
		}
	}
	assert( ilayer_min == ilayer_max );

	ndim_co = field.GetNDimCoord();
	ndim_va = 6;
	if( pData == 0 ){
		pData = new double [npo*(ndim_co+ndim_va)]; 
	}

	unsigned int icoun = 0;

	assert( field.GetNodeSegInNodeAry(BUBBLE).id_na_va != 0 );
	{
		const CField::CNodeSegInNodeAry& nsna_b = field.GetNodeSegInNodeAry(BUBBLE);
		assert( world.IsIdNA(nsna_b.id_na_va) );
		const Fem::Field::CNodeAry& na_val = world.GetNA(nsna_b.id_na_va);
		unsigned int id_ns_v;
		{
			if(      nsna_b.id_ns_va != 0 ) id_ns_v = nsna_b.id_ns_va;
			else if( nsna_b.id_ns_ve != 0 ) id_ns_v = nsna_b.id_ns_ve;
			else if( nsna_b.id_ns_ac != 0 ) id_ns_v = nsna_b.id_ns_ac;
			else{ assert(0); }
		}
		assert( na_val.IsSegID(id_ns_v) );
		const Fem::Field::CNodeAry::CNodeSeg& ns_val = na_val.GetSeg(id_ns_v);
		assert( ns_val.Length() == ndim_co*(ndim_co+1)/2 );
		{
			unsigned int id_na_c_co = field.GetNodeSegInNodeAry(CORNER).id_na_co;
			unsigned int id_ns_c_co = field.GetNodeSegInNodeAry(CORNER).id_ns_co;
			assert( world.IsIdNA(id_na_c_co) );
			const CNodeAry& na_c_co = world.GetNA(id_na_c_co);
			assert( na_c_co.IsSegID(id_ns_c_co) );
			const CNodeAry::CNodeSeg& ns_c_co = na_c_co.GetSeg(id_ns_c_co);
			assert( world.IsIdNA(nsna_b.id_na_va) );
			const CNodeAry& na_va = world.GetNA(nsna_b.id_na_va);
			assert( na_va.IsSegID(id_ns_v) );
			unsigned int noes[64];
			const std::vector<unsigned int>& aIdEA = field.GetAryIdEA();
			for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
				unsigned int id_ea = aIdEA[iiea];
				const CElemAry& ea = world.GetEA(id_ea);
				const CElemAry::CElemSeg& es_c_co = field.GetElemSeg(id_ea,CORNER,false,world);
				assert( es_c_co.GetIdNA() == id_na_c_co );
				const CElemAry::CElemSeg& es_b_va = field.GetElemSeg(id_ea,BUBBLE,true,world);
				assert( es_b_va.GetIdNA() == nsna_b.id_na_va );
				for(unsigned int ielem=0;ielem<ea.Size();ielem++){
					double coord_cnt[3];
					{	// get gravity center
						for(unsigned int idim=0;idim<ndim_co;idim++){ coord_cnt[idim] = 0.0; }
						const unsigned int nnoes = es_c_co.Length();
						es_c_co.GetNodes(ielem,noes);
						double coord[3];
						for(unsigned int inoes=0;inoes<nnoes;inoes++){
							unsigned int ipoi0 = noes[inoes];
							assert( ipoi0 < na_c_co.Size() );
							ns_c_co.GetValue(ipoi0,coord);
							for(unsigned int idim=0;idim<ndim_co;idim++){ coord_cnt[idim] += coord[idim]; }
						}
						for(unsigned int idim=0;idim<ndim_co;idim++){ coord_cnt[idim] /= nnoes; }
					}
					double value[6];	assert( ns_val.Length() <= 6 );
					{	// get tensor value
						es_b_va.GetNodes(ielem,noes);
						unsigned int ipoi0 = noes[0];
						assert( ipoi0 < na_va.Size() );
						ns_val.GetValue(ipoi0,value);
					}
					if(      ndim_co == 2 ){
						double pvs[2],ls, pvl[2],ll;
						GetPrincipleVector_STSR2(value,pvs,ls, pvl,ll);
						double* ptr = pData + icoun*(ndim_co+ndim_va);
						ptr[0] = coord_cnt[0];
						ptr[1] = coord_cnt[1];
						ptr[2] = pvs[0];
						ptr[3] = pvs[1];
						ptr[4] = ls;
						ptr[5] = pvl[0];
						ptr[6] = pvl[1];
						ptr[7] = ll;
					}
					else if( ndim_co == 3 ){
					}
					icoun++;
				} // end ielem
			} // end iei
		} // end if
	}
	
	assert( icoun == npo );
	
	return true;
}

bool CDrawerVector::Update(const Fem::Field::CFieldWorld& world)
{
	const Fem::Field::CField& field = world.GetField(id_field);
	if( field.IsPartial() ){
		std::cout  << "Error!-->Not Implemented" << std::endl;
		getchar();
		assert(0);
	}
	////////////////
	if( field.GetFieldType() == VECTOR2 || field.GetFieldType() == VECTOR3 ){ 
		itype = 0;
		return this->Update_VECTOR(world); 
	}
	else if( field.GetFieldType() == STSR2 ){ 
		itype = 1;
		return this->Update_SSTR2(world);  
	}
	return true;
}

bool CDrawerVector::Set(unsigned int id_field, const Fem::Field::CFieldWorld& world)
{
	if( !world.IsIdField(id_field) ) return false;
	this->id_field = id_field;
	if( pData != 0 ){ delete pData; pData=0; }
	this->Update(world);
	{
		const CField& field = world.GetField(id_field);
		unsigned int ndim = field.GetNDimCoord();
		if(      ndim == 2 ){ this->sutable_rot_mode = 1; }
		else if( ndim == 3 ){ this->sutable_rot_mode = 3; }
	}
	return true;
}

void CDrawerVector::Draw() const{
	if( npo == 0 ) return;
	assert( pData != 0 );
   ::glEnable(GL_DEPTH_TEST);
	::glDisable(GL_TEXTURE_2D);
	::glColor3d(0.0,0.0,0.0);
	::glLineWidth(2);
	::glBegin(GL_LINES);
	if (ndim_co == 2) {
		if(      itype == 0 ){	// 2d vector
			assert( ndim_va == 2 );
			for(unsigned int ipo=0;ipo<npo;ipo++){
				double* pCo = pData + ipo*(ndim_co+ndim_va);
				double* pVa = pData + ipo*(ndim_co+ndim_va) + ndim_co;
				::glVertex2dv(pCo);
				assert( ndim_va == 2 );
				::glVertex2d(pCo[0]+pVa[0],pCo[1]+pVa[1]);
			}
		}
		else if( itype == 1 ){	// 2d symmetric tensor
			assert( ndim_va == 6 );
			for(unsigned int ipo=0;ipo<npo;ipo++){
				double* pCo = pData + ipo*(ndim_co+ndim_va);
				double* pVa = pData + ipo*(ndim_co+ndim_va) + ndim_co;
				if(pVa[2]>0){	::glColor3d(0,0,1); }
				else{			::glColor3d(1,0,0); }
				::glVertex2dv(pCo);
				::glVertex2d(pCo[0]+pVa[0],pCo[1]+pVa[1]);
				::glVertex2dv(pCo);
				::glVertex2d(pCo[0]-pVa[0],pCo[1]-pVa[1]);
				////
				if(pVa[5]>0){	::glColor3d(0,0,1); }
				else{			::glColor3d(1,0,0); }				
				::glVertex2dv(pCo);
				::glVertex2d(pCo[0]+pVa[3],pCo[1]+pVa[4]);
				::glVertex2dv(pCo);
				::glVertex2d(pCo[0]-pVa[3],pCo[1]-pVa[4]);
			}		
		}
	}
	else if ( ndim_co == 3 ){
		if( itype == 0 ){
			if(      ndim_va == 2 ){	// 2d vector with 3d coord
				for(unsigned int ipo=0;ipo<npo;ipo++){
					double* pCo = pData + ipo*(ndim_co+ndim_va);
					double* pVa = pData + ipo*(ndim_co+ndim_va) + ndim_co;
					::glVertex3dv(pCo);
					::glVertex3d(pCo[0]+pVa[0],pCo[1]+pVa[1],pCo[2]);
				}
			}
			else if( ndim_va == 3 ){	// 3d vector
				for(unsigned int ipo=0;ipo<npo;ipo++){
					double* pCo = pData + ipo*(ndim_co+ndim_va);
					double* pVa = pData + ipo*(ndim_co+ndim_va) + ndim_co;
					::glVertex3dv(pCo);
					::glVertex3d(pCo[0]+pVa[0],pCo[1]+pVa[1],pCo[2]+pVa[2]);
				}
			}		
		}
	}
	::glEnd();
}
