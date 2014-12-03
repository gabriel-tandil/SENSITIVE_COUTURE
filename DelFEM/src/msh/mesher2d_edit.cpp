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
#pragma warning(disable: 4786)
#pragma warning(disable: 4996)
#endif
#define for if(0);else for

#include <stdio.h>
#include <set>
#include <vector>
#include <stack>
#include <cassert>
#include <math.h>

#include <iostream>
#include <fstream>

#include "delfem/msh/meshkernel2d.h"
#include "delfem/mesher2d_edit.h"

using namespace Msh;
using namespace Com;

void CMesher2D_Edit::SmoothingMesh_Laplace(unsigned int num_iter, double elen)
{
	std::vector<int> flg_vec;
	flg_vec.resize(aVec2D.size(),0);
	for(unsigned int ivertex=0;ivertex<m_aVertex.size();ivertex++){
		const unsigned int iv = m_aVertex[ivertex].v;
		flg_vec[iv] = 2;
	}
	for(unsigned int iiter=0;iiter<num_iter;iiter++){
		if( iiter != 0 ){
			for(unsigned int ivec=0;ivec<aVec2D.size();ivec++){	
				if( flg_vec[ivec] == 1  ){ flg_vec[ivec] = 0; }
			}
		}
		for(unsigned int itriary=0;itriary<m_aTriAry.size();itriary++){
			const std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;
			for(unsigned int itri=0;itri<aTri.size();itri++){	// search point around point
			for(unsigned int inotri=0;inotri<3;inotri++){
				const unsigned int ipoin = aTri[itri].v[inotri];
				assert( ipoin < aVec2D.size() );
				if( flg_vec[ipoin] != 0 ) continue;
				flg_vec[ipoin] = 1;

				////////////////////////////////
				// kernel triangle and direction
				const unsigned int itri_ini = itri;
				const unsigned int inoel_c_ini = inotri;
				// current triangle and direction
				unsigned int itri0= itri_ini;
				unsigned int inoel_c0 = inoel_c_ini;
				unsigned int inoel_b0 = noelTriEdge[inoel_c0][0];
				unsigned int inoel_f0 = noelTriEdge[inoel_c0][1];
				////////////////
				CVector2D vec_delta(0, 0);
				unsigned int ntri_around = 1;
				bool is_bound_flg = false;
				for(;;){
					assert( itri0 < aTri.size() );
					assert( inoel_c0 < 3 );
					assert( aTri[itri0].v[inoel_c0] == ipoin );
					{
						unsigned int ipo0 = aTri[itri0].v[inoel_b0];
						double x_delta = aVec2D[ipo0].x - aVec2D[ipoin].x;
						double y_delta = aVec2D[ipo0].y - aVec2D[ipoin].y;
						double dist = sqrt( x_delta*x_delta+y_delta*y_delta );
						if( dist < elen*1.0e-3 ) dist = 1.0e-3*elen;
						x_delta /= dist; y_delta /= dist;
						const double ratio = dist / elen;
						double dtmp2;
						{
							if( ratio >= 0.0 && ratio < 0.3 ){
								dtmp2 = 0;
							}
							else if( ratio >= 0.3 && ratio < 0.5 ){
								dtmp2 = 0.5*cos( 3.14*5*(ratio-0.5) )+0.5;
							}
							else if( ratio >= 0.5 && ratio < 1.0 ){ 
								dtmp2 = 1.0; 
							}
							else if( ratio >= 1.0 && ratio < 2.0 ){
								dtmp2 = 1.4+0.4*cos( 3.14*(ratio-2.0) );
							}
							else{ 
								dtmp2 = 1.8;
							}
//							dtmp2 = 1.0;
						}
						vec_delta.x += dtmp2*dist*x_delta;
						vec_delta.y += dtmp2*dist*y_delta;
						ntri_around++;
					}
					if( aTri[itri0].g2[inoel_b0] == -2 ){
						unsigned int itri1 = aTri[itri0].s2[inoel_b0];
            const unsigned int rel01 = aTri[itri0].r2[inoel_b0];
						unsigned int inoel_c1 = relTriTri[rel01][inoel_c0];
						unsigned int inoel_b1 = relTriTri[rel01][ noelTriEdge[inoel_c0][1] ];
						unsigned int inoel_f1 = relTriTri[rel01][ noelTriEdge[inoel_c0][0] ];
						assert( itri1 < aTri.size() );
            assert( relTriTri[rel01][inoel_b0] < 4 );
            if( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] != itri0 ){
              std::cout << itri0 << " " << itri1 << " " << aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] << "  " << relTriTri[rel01][inoel_b0] << " " << rel01 << " " << inoel_b0 << std::endl;              
            }
						assert( aTri[itri1].s2[ relTriTri[rel01][inoel_b0] ] == itri0 );  // this line is sometime problem
						assert( aTri[itri1].v[inoel_c1] == ipoin ); // this line is sometimes problem also
						if( itri1 == itri_ini ) break;
						itri0 = itri1;
						inoel_c0 = inoel_c1;
						inoel_b0 = inoel_b1;
						inoel_f0 = inoel_f1;
						}
						else{	// this point is on the boundary : please dont move
						is_bound_flg = true;
						flg_vec[ipoin] = 2;
						break;
					}
				}
				if( is_bound_flg ) continue;
					aVec2D[ipoin].x += vec_delta.x / ntri_around;
					aVec2D[ipoin].y += vec_delta.y / ntri_around;
			}	// inotri
			}	// itri
		}	// itriary
	}	// iiter
}


// relocate 2 two end point of edge to the destination 
bool CMesher2D_Edit::MovePointsOnEdge
(const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad, 
 std::vector<int>& aflg_ismov,
 const CVector2D& dist_mov_s, const CVector2D& dist_mov_e )
{
	assert( cad_2d.IsElemID(Cad::EDGE,id_e_cad) );

	unsigned int id_v_cad_s = cad_2d.GetIdVertex_Edge(id_e_cad,true );
  unsigned int id_v_cad_e = cad_2d.GetIdVertex_Edge(id_e_cad,false);
	
	unsigned int ibarary;
	{
		unsigned int itype;
    if( !this->FindElemLocType_CadIDType(ibarary,itype,id_e_cad,Cad::EDGE) ){ return false; }
    assert( itype == 1 );
	}

	const std::vector<SBar>& aBar = m_aBarAry[ibarary].m_aBar;
	unsigned int iv_s, iv_e;
	{
		unsigned int itype;
		this->FindElemLocType_CadIDType(iv_s,itype,id_v_cad_s,Cad::VERTEX); assert( itype == 0 );
		this->FindElemLocType_CadIDType(iv_e,itype,id_v_cad_e,Cad::VERTEX); assert( itype == 0 );
	}

	{
		const unsigned int nbar = aBar.size();
		assert( aBar[0].v[0] == iv_s );
		for(unsigned int ibar=0;ibar<nbar-1;ibar++){
			assert( aBar[ibar].v[1] == aBar[ibar+1].v[0] );
		}
		assert( aBar[nbar-1].v[1] == iv_e );
	}
  
  std::vector<Com::CVector2D> aCo;
  cad_2d.GetCurveAsPolyline(id_e_cad,aCo,aBar.size(),dist_mov_s,dist_mov_e);
  assert( aCo.size() == aBar.size()-1 );
  
  {
    const unsigned int nbar = aBar.size();
    if( aflg_ismov[iv_s] == 0 ){
      aVec2D[iv_s] = dist_mov_s;
      aflg_ismov[iv_s] = 1;
    }
    for(unsigned int ibar=1;ibar<nbar;ibar++){
      unsigned int iv0 = aBar[ibar].v[0];
      if( aflg_ismov[iv0] != 0 ) continue;
      aflg_ismov[iv0] = 1;
      aVec2D[iv0] = aCo[ibar-1];
    }
    if( aflg_ismov[iv_e] == 0 ){
      aVec2D[iv_e] = dist_mov_e;
      aflg_ismov[iv_e] = 1;
    }
  }
  return true;  

}

// move one point (id_v_cad_mov) of edge
bool CMesher2D_Edit::MovePointsOnEdge
(const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad, 
 std::vector<int>& aflg_ismov,
 unsigned int id_v_cad_mov, const CVector2D& dist_mov )
{
	assert( cad_2d.IsElemID(Cad::EDGE,id_e_cad) );

	unsigned int id_v_cad_s = cad_2d.GetIdVertex_Edge(id_e_cad,true );
  unsigned int id_v_cad_e = cad_2d.GetIdVertex_Edge(id_e_cad,false);
	
	if( id_v_cad_mov == id_v_cad_s ){
		return this->MovePointsOnEdge(cad_2d,id_e_cad,
			aflg_ismov, 
			dist_mov, 
			cad_2d.GetVertexCoord(id_v_cad_e) );
	}
	else{
		assert( id_v_cad_mov == id_v_cad_e );
		return this->MovePointsOnEdge(cad_2d,id_e_cad,
			aflg_ismov, 
			cad_2d.GetVertexCoord(id_v_cad_s),
			dist_mov );
	}
	return true;
}

// move poins on edge following cad
bool CMesher2D_Edit::MovePointsOnEdge
(const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad, 
 std::vector<int>& aflg_ismov )
{
	assert( cad_2d.IsElemID(Cad::EDGE,id_e_cad) );

	unsigned int id_v_cad_s = cad_2d.GetIdVertex_Edge(id_e_cad,true );
  unsigned int id_v_cad_e = cad_2d.GetIdVertex_Edge(id_e_cad,false);	
	const CVector2D& vec_s = cad_2d.GetVertexCoord(id_v_cad_s);
	const CVector2D& vec_e = cad_2d.GetVertexCoord(id_v_cad_e);
	return this->MovePointsOnEdge(cad_2d,id_e_cad,
                                aflg_ismov, 
                                vec_s, vec_e);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

bool CMesher2D_Edit::FitMeshToCad_Edge
(const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad_mov,  
 unsigned int& itype_operation )
{
	unsigned int id_v_cad_s = cad_2d.GetIdVertex_Edge(id_e_cad_mov,true );
  unsigned int id_v_cad_e = cad_2d.GetIdVertex_Edge(id_e_cad_mov,false);
	
  std::set<unsigned int> setIdL;
  std::vector<unsigned int> aIdCadE_s;
  for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v_cad_s);!itrv.IsEnd();itrv++){
    unsigned int id_e0;   bool is_same_dir0;
    itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
    if( id_e0 != id_e_cad_mov ) aIdCadE_s.push_back(id_e0);
    setIdL.insert(itrv.GetIdLoop());
  }
  std::vector<unsigned int> aIdCadE_e;
  for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v_cad_e);!itrv.IsEnd();itrv++){
    unsigned int id_e0;   bool is_same_dir0;
    itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
    if( id_e0 != id_e_cad_mov ) aIdCadE_e.push_back(id_e0);
    setIdL.insert(itrv.GetIdLoop());
  }  
  
  const double ave_edge_len = this->GetAverageEdgeLength(cad_2d,setIdL);

	bool is_precomp = false;
  std::vector< Com::CVector2D > aVec_tmp = aVec2D;	// copy in case of mesh inversion
	{
		CVector2D vec_dist_s = cad_2d.GetVertexCoord(id_v_cad_s);
		CVector2D vec_dist_e = cad_2d.GetVertexCoord(id_v_cad_e);
		std::vector<int> aflg_ismov;	// flag if the coord is updated
		aflg_ismov.resize(aVec2D.size(),0);
		Com::CVector2D delta_s, delta_e;
		{
			unsigned int ind_v_mov_s, itype_s;
			if( !this->FindElemLocType_CadIDType(ind_v_mov_s,itype_s,id_v_cad_s,Cad::VERTEX) ){
				return false;   // this point is not exist in mesh
			}
			assert( itype_s == 0 );
			assert( ind_v_mov_s < this->m_aVertex.size() );
			////////////////
			unsigned int ind_v_mov_e, itype_e;
			if( !this->FindElemLocType_CadIDType(ind_v_mov_e,itype_e,id_v_cad_e,Cad::VERTEX) ){
				return false;   // this point is not exist in mesh
			}
			assert( itype_e == 0 );
			assert( ind_v_mov_e < this->m_aVertex.size() );
			////////////////
			unsigned iv_s = m_aVertex[ind_v_mov_s].v;
			unsigned iv_e = m_aVertex[ind_v_mov_e].v;
			delta_s = vec_dist_s - aVec2D[iv_s];
			delta_e = vec_dist_e - aVec2D[iv_e];
		}
		if(  SquareLength(delta_s-delta_e)<ave_edge_len*ave_edge_len*1.0e-10
			&& SquareLength(delta_s)>ave_edge_len*ave_edge_len*1.0e-10
			&& aLamTnsr.size() == aVec2D.size()*4 ){	// Translate
			for(unsigned int iv=0;iv<aVec2D.size();iv++){
        const double l[4] =  { aLamTnsr[iv*4+0], aLamTnsr[iv*4+1], aLamTnsr[iv*4+2], aLamTnsr[iv*4+3] };        
				this->aVec2D[iv].x += l[0]*delta_s.x + l[1]*delta_s.y;
				this->aVec2D[iv].y += l[2]*delta_s.x + l[3]*delta_s.y;
			}
			is_precomp = true;
		}
		else{
			for(unsigned int iid_e=0;iid_e<aIdCadE_s.size();iid_e++){
				const unsigned int id_e0 = aIdCadE_s[iid_e];
				this->MovePointsOnEdge(cad_2d,id_e0,aflg_ismov,  id_v_cad_s,vec_dist_s);
			}
			for(unsigned int iid_e=0;iid_e<aIdCadE_e.size();iid_e++){
				const unsigned int id_e0 = aIdCadE_e[iid_e];
				this->MovePointsOnEdge(cad_2d,id_e0,aflg_ismov,  id_v_cad_e,vec_dist_e);
			}
			this->MovePointsOnEdge(cad_2d,id_e_cad_mov,aflg_ismov);
			this->SmoothingMesh_Laplace(2,ave_edge_len);
		}
	}

  ////////////////

	bool is_inverted;
	double max_aspect;
	this->CheckMeshQuality(is_inverted,max_aspect,ave_edge_len);
	itype_operation = 0;
	if( is_inverted ){ aVec2D = aVec_tmp; }
	else{
		itype_operation += 1;
		if( max_aspect < 15.0 ){ return true; }
	}

	unsigned int num_reconnect;
	this->SmoothingMesh_Delaunay(num_reconnect);
	if( num_reconnect != 0 && is_precomp ){
		this->SmoothingMesh_Laplace(4,ave_edge_len);
		this->Precomp_FitMeshToCad(cad_2d,aLamTnsr_Vtx);
	}
//	std::cout << "Num Reconnect " << num_reconnect << std::endl;
	if( num_reconnect != 0 ){ itype_operation += 2; }
	if( !is_inverted ) return true;
	return false;
}

bool CMesher2D_Edit::FitMeshToCad_Loop
(const Cad::CCadObj2D& cad_2d, unsigned int id_l_cad_mov, 
 unsigned int& itype_operation )
{
	double ave_edge_len;
	{
    std::set<unsigned int> setIdL;
		for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l_cad_mov);!itrl.IsEnd();itrl++){
			unsigned int id_e_cad;   bool is_same_dir;
			itrl.GetIdEdge(id_e_cad,is_same_dir);
			unsigned int id_l_l, id_l_r;
			cad_2d.GetIdLoop_Edge(id_l_l,id_l_r,id_e_cad);
			setIdL.insert(id_l_l);
			setIdL.insert(id_l_r);
		}
    ave_edge_len = this->GetAverageEdgeLength(cad_2d,setIdL);
	}
  bool is_precomp = false;
	std::vector< CVector2D > aVec_tmp = aVec2D;
//  std::cout << aLamTnsr.size() << " " << aVec2D.size()*4 << std::endl;
  if( aLamTnsr.size() == aVec2D.size()*4 ){
    Com::CVector2D delta(0,0);
    { // calc average of the vertex movement
      unsigned int icnt=0;
      for(std::auto_ptr<Cad::IItrLoop> itrl=cad_2d.GetPtrItrLoop(id_l_cad_mov);!itrl->IsEndChild();itrl->ShiftChildLoop()){
        for(;!itrl->IsEnd();(*itrl)++){
          unsigned int id_v_cad = itrl->GetIdVertex();		
          CVector2D vec_dist = cad_2d.GetVertexCoord(id_v_cad);
          unsigned int ind_v_mov, itype;
          if( !this->FindElemLocType_CadIDType(ind_v_mov,itype,id_v_cad,Cad::VERTEX) ){ continue; }
          assert( itype == 0 );
          assert( ind_v_mov < this->m_aVertex.size() );
          unsigned iv = m_aVertex[ind_v_mov].v;          
          delta.x += vec_dist.x - aVec2D[iv].x;
          delta.y += vec_dist.y - aVec2D[iv].y;          
          icnt++;
        }
      } 
      delta *= (1.0/icnt);
    }
    for(unsigned int iv=0;iv<aVec2D.size();iv++){
			const double l[4] =  { aLamTnsr[iv*4+0], aLamTnsr[iv*4+1], aLamTnsr[iv*4+2], aLamTnsr[iv*4+3] };      
//      std::cout << iv << " " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
      this->aVec2D[iv].x += l[0]*delta.x + l[1]*delta.y;
      this->aVec2D[iv].y += l[2]*delta.x + l[3]*delta.y;
    }
    is_precomp = true;    
  }
	else{
		std::vector<int> aflg_ismov;
		aflg_ismov.resize(aVec2D.size(),0);
		for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l_cad_mov);!itrl.IsEndChild();itrl.ShiftChildLoop()){
			for(;!itrl.IsEnd();itrl++){
				unsigned int id_e_cad;   bool is_same_dir;
				itrl.GetIdEdge(id_e_cad,is_same_dir);
				unsigned int id_v_cad = itrl.GetIdVertex();		
				this->MovePointsOnEdge(cad_2d,id_e_cad,aflg_ismov);		// move this edge
				// points around this edge
				CVector2D vec_dist = cad_2d.GetVertexCoord(id_v_cad);
				const unsigned int idv_ahead  = itrl.GetIdVertex_Ahead();
				const unsigned int idv_behind = itrl.GetIdVertex_Behind();
				for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v_cad);!itrv.IsEnd();itrv++){
					unsigned int id_e0;   bool is_same_dir0;
					itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
					// check whether this edge belong to the loop
					unsigned int id_vs = cad_2d.GetIdVertex_Edge(id_e0,true );
          unsigned int id_ve = cad_2d.GetIdVertex_Edge(id_e0,false);
					if( id_vs == id_v_cad && (id_ve == idv_ahead || id_ve == idv_behind) ) continue;
					if( id_ve == id_v_cad && (id_vs == idv_ahead || id_vs == idv_behind) ) continue;
					this->MovePointsOnEdge(cad_2d,id_e0,aflg_ismov,  id_v_cad,vec_dist);
				}
			}
		}
    this->SmoothingMesh_Laplace(2,ave_edge_len);    
	}  

	bool is_inverted;
	double max_aspect;
	this->CheckMeshQuality(is_inverted,max_aspect,ave_edge_len);
	////////////////
//////
//	itype_operation = 1;
//  return true;
//////
//	std::cout << max_edge_len_ratio << " " << max_aspect << std::endl;
	if( is_inverted ){ aVec2D = aVec_tmp; }
	else{
		itype_operation += 1;
		if( max_aspect < 15.0 ){ return true; }
	}

	unsigned int num_reconnect;
	this->SmoothingMesh_Delaunay(num_reconnect);
  if( num_reconnect != 0 && is_precomp ){
		this->SmoothingMesh_Laplace(4,ave_edge_len);
		this->Precomp_FitMeshToCad(cad_2d,aLamTnsr_Vtx);
	}  
	if( num_reconnect != 0 ){ itype_operation += 2; }
	if( !is_inverted ) return true;
	return false;
}

bool CMesher2D_Edit::FitMeshToCad_Slider_UsingPrecomp
(const Cad::CCadObj2D& cad_2d, unsigned int& itype_operation,
 double delta)
{    
	double ave_edge_len;
	{
		std::set<unsigned int> setIdL;
		const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);
		for(unsigned int iidl=0;iidl<aIdL.size();iidl++){ setIdL.insert(aIdL[iidl]); }
		ave_edge_len = this->GetAverageEdgeLength(cad_2d,setIdL);
	}

  ////
	std::vector< CVector2D > aVec_tmp = aVec2D;
  //  std::cout << aLamTnsr.size() << " " << aVec2D.size()*4 << std::endl;
  if( aLamTnsr.size() != aVec2D.size()*4 ){ return false; }    
  for(unsigned int iv=0;iv<aVec2D.size();iv++){
    const double l[4] =  { aLamTnsr[iv*4+0], aLamTnsr[iv*4+1], aLamTnsr[iv*4+2], aLamTnsr[iv*4+3] };      
    //      std::cout << iv << " " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    this->aVec2D[iv].x += l[0]*delta;
    this->aVec2D[iv].y += l[2]*delta;
  }
  
	bool is_inverted;
	double max_aspect;
	this->CheckMeshQuality(is_inverted,max_aspect,ave_edge_len);
	////////////////
  
  //	itype_operation = 1;
  //  return true;
  //	std::cout << max_edge_len_ratio << " " << max_aspect << std::endl;
	if( is_inverted ){ aVec2D = aVec_tmp; }
	else{
		itype_operation += 1;
		if( max_aspect < 15.0 ){ return true; }
	}
  
	unsigned int num_reconnect;
	this->SmoothingMesh_Delaunay(num_reconnect);
  if( num_reconnect != 0 ){
		this->SmoothingMesh_Laplace(4,ave_edge_len);
		this->Precomp_FitMeshToCad(cad_2d,aLamTnsr_Vtx);
	}  
	if( num_reconnect != 0 ){ itype_operation += 2; }
	if( !is_inverted ) return true;
	return false;  
}

bool CMesher2D_Edit::FitMeshToCad_All(const Cad::CCadObj2D& cad_2d, unsigned int& itype_operation )
{
	double ave_edge_len;
	{
		std::set<unsigned int> setIdL;
		const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);
		for(unsigned int iidl=0;iidl<aIdL.size();iidl++){ setIdL.insert(aIdL[iidl]); }
		ave_edge_len = this->GetAverageEdgeLength(cad_2d,setIdL);
	}

	std::vector< CVector2D > aVec_tmp = aVec2D;
	{
		std::vector<int> aflg_ismov;
		aflg_ismov.resize(aVec2D.size(),0);
		std::vector<unsigned int> aIdEdge = cad_2d.GetAryElemID(Cad::EDGE);
		for(unsigned int iid_e=0;iid_e<aIdEdge.size();iid_e++){
			unsigned int id_e_cad = aIdEdge[iid_e];
			this->MovePointsOnEdge(cad_2d,id_e_cad,aflg_ismov);
		}
	}	
	
	this->SmoothingMesh_Laplace(2,ave_edge_len);
	
	bool is_inverted;
	double max_aspect;
	this->CheckMeshQuality(is_inverted,max_aspect,ave_edge_len);
	////////////////
	itype_operation = 0;
	if( is_inverted ){ aVec2D = aVec_tmp; }
	else{
		itype_operation += 1;
		if( max_aspect < 15.0 ){ return true; }
	}
	unsigned int num_reconnect;
	this->SmoothingMesh_Delaunay(num_reconnect);
	itype_operation += 2;
	if( !is_inverted ) return true;
	return false;
}

bool CMesher2D_Edit::FitMeshToCad_Vertex
(const Cad::CCadObj2D& cad_2d, unsigned int id_v_cad_mov, 
 unsigned int& itype_operation )
{
  if( !cad_2d.IsElemID(Cad::VERTEX,id_v_cad_mov) ) return false;
	
  {
    unsigned int iv_mov, itype;
    if( !this->FindElemLocType_CadIDType(iv_mov,itype,id_v_cad_mov,Cad::VERTEX) ){
      return false;   // this point is not exist in mesh
    }
    assert( itype == 0 );
    assert( iv_mov < this->m_aVertex.size() );
  }
  
  std::vector<unsigned int> aEdgeID;
	double ave_edge_len;
	{
		std::set<unsigned int> setIdL;
		for(Cad::CBRepSurface::CItrVertex itrv=cad_2d.GetItrVertex(id_v_cad_mov);!itrv.IsEnd();itrv++){
      setIdL.insert(itrv.GetIdLoop());
			unsigned int id_e0;   bool is_same_dir0;
			itrv.GetIdEdge_Behind(id_e0,is_same_dir0);
      if( cad_2d.IsElemID(Cad::EDGE,id_e0) ){ aEdgeID.push_back(id_e0); }
		}
		ave_edge_len = this->GetAverageEdgeLength(cad_2d,setIdL);	// calc averge mesh length
	}

	bool is_followed = true;	// is the mesh can follow cad movement
	CVector2D dist_mov0 = cad_2d.GetVertexCoord(id_v_cad_mov);	// put dist_mov0 the cad vertex position
	CVector2D dist_mov1 = dist_mov0;	// the distination to the vertex relocation
  
	bool is_precomp = false;
	std::vector< CVector2D > aVec_tmp = aVec2D;	// store all vertices before movement
//  std::cout << aLamTnsr.size() << " " << aVec2D.size() << std::endl;
	if( aLamTnsr.size() == aVec2D.size()*4 ){
		Com::CVector2D delta;
		{
			unsigned int ind_v_mov, itype;
			if( !this->FindElemLocType_CadIDType(ind_v_mov,itype,id_v_cad_mov,Cad::VERTEX) ){
				return false;   // this point is not exist in mesh
			}
			assert( itype == 0 );
			assert( ind_v_mov < this->m_aVertex.size() );
			unsigned iv = m_aVertex[ind_v_mov].v;
			delta = dist_mov0 - aVec2D[iv];
		}
		for(unsigned int iv=0;iv<aVec2D.size();iv++){
			const double l[4] =  { aLamTnsr[iv*4+0], aLamTnsr[iv*4+1], aLamTnsr[iv*4+2], aLamTnsr[iv*4+3] };
//      std::cout << iv << " " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
			this->aVec2D[iv].x += l[0]*delta.x + l[1]*delta.y;
			this->aVec2D[iv].y += l[2]*delta.x + l[3]*delta.y;
		}
		is_precomp = true;
	}
	else {	// moveedge
		std::vector<int> aflg_ismov;
		aflg_ismov.resize(aVec2D.size(),0);
		for(unsigned int iid_e=0;iid_e<aEdgeID.size();iid_e++){
			const unsigned int id_e = aEdgeID[iid_e];	// edge around vertices
			this->MovePointsOnEdge(cad_2d,id_e,aflg_ismov,  id_v_cad_mov,dist_mov1);	// move mesh point on edge
		}
		unsigned int iv_mov, itype;
		this->FindElemLocType_CadIDType(iv_mov,itype,id_v_cad_mov,Cad::VERTEX);	// mesh point on vertex
		assert( itype == 0 );
		this->aVec2D[iv_mov].x = dist_mov1.x;
		this->aVec2D[iv_mov].y = dist_mov1.y;
		this->SmoothingMesh_Laplace(2,ave_edge_len);
	}

	bool is_inverted;
	double max_aspect;
	this->CheckMeshQuality(is_inverted,max_aspect,ave_edge_len);
//  std::cout << "mesh qualitity " << is_inverted << " " << max_aspect << std::endl;

	itype_operation = 0;
	if( is_inverted ){ aVec2D = aVec_tmp; } // if negative area, reset mesh poisition
	else{
		itype_operation += 1;
		if( max_aspect < 15.0 ){ // if mesh is destorted, reconnect mesh even if mesh can follow cad
			return is_followed; // no negative or distorted triangle. great! :)
		}
	}
	// Reconnect mesh
	unsigned int num_reconnect;
	this->SmoothingMesh_Delaunay(num_reconnect);
	if( num_reconnect != 0 && is_precomp ){
		this->SmoothingMesh_Laplace(4,ave_edge_len);
		this->Precomp_FitMeshToCad(cad_2d,aLamTnsr_Vtx);
	}
	if( num_reconnect != 0 ){ itype_operation += 2; }
	if( !is_inverted ) return is_followed;
	return false;
}

////////////////

void CalcWeight_MVC(std::vector<double>& aLamTnsr,
                    unsigned int iv0,
                    ////
                    unsigned int nloop,                       
                    const std::vector<unsigned int>& aIndLoop, // nloop+1
                    const std::vector<CVector2D>& aVtx, 
                    const std::vector<double>& aValVtx, // nvtx
                    const CVector2D& v0)
{
  const unsigned int nvtx = aVtx.size();  
  std::vector<double> buff;
  buff.resize(nvtx);
  assert( aIndLoop.size() == nloop+1 );
  assert( aValVtx.size() == nvtx*4 );
  assert( buff.size() == nvtx );
  for(unsigned int iloop=0;iloop<nloop;iloop++){
    unsigned int ivtx0 = aIndLoop[iloop];
    unsigned int ivtx1 = aIndLoop[iloop+1];
    for(unsigned int ivtx=ivtx0;ivtx<ivtx1;ivtx++){
      const CVector2D& vb = (ivtx==ivtx0) ? aVtx[ivtx1-1] : aVtx[ivtx-1];
      const CVector2D& vm = aVtx[ivtx];
      const CVector2D& vf = (ivtx==ivtx1-1) ? aVtx[ivtx0] : aVtx[ivtx+1];
      double b0 = (vb-v0).Length();
      double m0 = (vm-v0).Length();
      double f0 = (vf-v0).Length();
      double bm = (vm-vb).Length();
      double fm = (vm-vf).Length();
      const double sin_b0m = TriArea(vb, v0, vm);
      const double sin_m0f = TriArea(vm, v0, vf);      
      const double cos_b0m = (b0*b0+m0*m0-bm*bm)/(2*b0*m0);
      const double cos_m0f = (f0*f0+m0*m0-fm*fm)/(2*f0*m0);      
      double tanh_b0m = sqrt( (1-cos_b0m)/(1+cos_b0m) );
      double tanh_m0f = sqrt( (1-cos_m0f)/(1+cos_m0f) );
      tanh_b0m = ( sin_b0m > 0 ) ? tanh_b0m : -tanh_b0m;
      tanh_m0f = ( sin_m0f > 0 ) ? tanh_m0f : -tanh_m0f;      
      buff[ivtx] = (tanh_b0m+tanh_m0f)/m0;      
    }
  }
  double sum = 0;
  for(unsigned int ivtx=0;ivtx<nvtx;ivtx++){ sum += buff[ivtx]; }
  for(unsigned int ivtx=0;ivtx<nvtx;ivtx++){ buff[ivtx] /= sum; }  
  for(unsigned int ivtx=0;ivtx<nvtx;ivtx++){ 
    aLamTnsr[iv0*4+0] += aValVtx[ivtx*4+0]*buff[ivtx]; 
    aLamTnsr[iv0*4+1] += aValVtx[ivtx*4+1]*buff[ivtx]; 
    aLamTnsr[iv0*4+2] += aValVtx[ivtx*4+2]*buff[ivtx]; 
    aLamTnsr[iv0*4+3] += aValVtx[ivtx*4+3]*buff[ivtx];     
  }  
}

class CVisiblePoint
{
public:
  CVisiblePoint(){
    is_visible = true;
    sqdist = -1;
  }
public:
  bool is_visible;
  unsigned int ivtx0;
  unsigned int ivtx1;
  unsigned int ivtx_nxt;
  double ratio;
  double sqdist;
};



void CalcWeight_PMVC
(std::vector<double>& aLamTnsr,
 unsigned int iv0,
 /////
 unsigned int nloop, 
 const std::vector<unsigned int>& aIndLoop, // nloop+1
 const std::vector<CVector2D>& aVtx, 
 const std::vector<double>& aValVtx, // nvtx
 const CVector2D& v0)
{
  const unsigned int nvtx = aVtx.size();
  assert( aIndLoop.size() == nloop+1 );
  assert( aValVtx.size() == nvtx*4 );
  std::vector< CVisiblePoint > aVis;
  aVis.resize(nvtx);
  for(unsigned int iloop=0;iloop<nloop;iloop++){
    const unsigned int ivtxs = aIndLoop[iloop];
    const unsigned int ivtxe = aIndLoop[iloop+1];
    for(unsigned int ivtx=ivtxs;ivtx<ivtxe;ivtx++){
      const unsigned int ivtxf = (ivtx==ivtxe-1) ? ivtxs : ivtx+1;
      aVis[ivtx].ivtx_nxt = ivtxf;
//      std::cout << "coord : " << ivtx << " " << aVtx[ivtx].x << " " << aVtx[ivtx].y << std::endl;
    }
  }
  for(unsigned int ivtx=0;ivtx<nvtx;ivtx++){
    for(unsigned int jloop=0;jloop<nloop;jloop++){
      const unsigned int jvtxs = aIndLoop[jloop];
      const unsigned int jvtxe = aIndLoop[jloop+1];
      for(unsigned int jvtx=jvtxs;jvtx<jvtxe;jvtx++){
        const unsigned int fvtx = (jvtx==jvtxe-1) ? jvtxs : jvtx+1;
        if( jvtx == ivtx || fvtx == ivtx ) continue;
        const double a0jf = TriArea(v0,aVtx[jvtx],aVtx[fvtx]); 
        if( a0jf < 0 ) continue;
        const double a0ij = TriArea(v0,aVtx[ivtx],aVtx[jvtx]);
        const double a0if = TriArea(v0,aVtx[ivtx],aVtx[fvtx]); 
        if( (a0ij > 0) == (a0if > 0) ) continue;
        const double aijf = TriArea(aVtx[ivtx],aVtx[jvtx],aVtx[fvtx]);
        if( aijf < 0 ){ // blocked 
          aVis[ivtx].is_visible = false;
          break;
        }
        else if( a0jf > aijf ){ // back
          const double ratio = a0if/(a0if-a0ij);
          const CVector2D& vm = ratio*aVtx[jvtx] + (1-ratio)*aVtx[fvtx];
          double sqdist = SquareLength(vm,v0);
          if( aVis[ivtx].sqdist < 0 || sqdist < aVis[ivtx].sqdist ){
            aVis[ivtx].sqdist = sqdist;
            aVis[ivtx].ivtx0 = jvtx;
            aVis[ivtx].ivtx1 = fvtx;
            aVis[ivtx].ratio = ratio;
          }
        }
      }
      if( !aVis[ivtx].is_visible ) break;
    }        
  }
  std::map<double,unsigned int> tmap;
  for(unsigned int ivtx=0;ivtx<nvtx;ivtx++){
    if( !aVis[ivtx].is_visible ) continue;
    double dx = aVtx[ivtx].x - v0.x;
    double dy = aVtx[ivtx].y - v0.y;
    double theta = atan2(dy, dx);
    tmap.insert( std::make_pair(theta,ivtx) );
  }  
  ////
  std::vector<unsigned int> aivp;  // array index visible point
  aivp.reserve(nvtx);
  for(std::map<double,unsigned int>::iterator itr=tmap.begin();itr!=tmap.end();itr++){
    unsigned int ivtx0 = itr->second;
    aivp.push_back(ivtx0);
  }
  aLamTnsr[iv0*4+0] = 0;
  aLamTnsr[iv0*4+1] = 0;
  aLamTnsr[iv0*4+2] = 0;
  aLamTnsr[iv0*4+3] = 0;
  double sum = 0;
  for(unsigned int iivp=0;iivp<aivp.size();iivp++){
    CVector2D vs, ve;
    double vals[4],vale[4]; 
    unsigned int ivtx0 = aivp[iivp];
    unsigned int ivtx1 = ( iivp == aivp.size()-1 ) ? aivp[0] : aivp[iivp+1];
    if( aVis[ivtx0].ivtx_nxt == ivtx1 ){
      vs = aVtx[ivtx0];
      ve = aVtx[ivtx1];        
      vals[0] = aValVtx[ivtx0*4+0];
      vals[1] = aValVtx[ivtx0*4+1];
      vals[2] = aValVtx[ivtx0*4+2];
      vals[3] = aValVtx[ivtx0*4+3];      
      vale[0] = aValVtx[ivtx1*4+0];
      vale[1] = aValVtx[ivtx1*4+1];
      vale[2] = aValVtx[ivtx1*4+2];
      vale[3] = aValVtx[ivtx1*4+3];
    }
    else{
//      assert( aVis[ivtx0].sqdist > 0 || aVis[ivtx1].sqdist > 0 );
      if( aVis[ivtx0].sqdist > 0 ){
        unsigned int jvtx0s = aVis[ivtx0].ivtx0; // visible
        unsigned int jvtx0e = aVis[ivtx0].ivtx1;
        double ratio0 = aVis[ivtx0].ratio;      
        vs   = ratio0*aVtx[   jvtx0s] + (1-ratio0)*aVtx[   jvtx0e];
        vals[0] = ratio0*aValVtx[jvtx0s*4+0] + (1-ratio0)*aValVtx[jvtx0e*4+0]; 
        vals[1] = ratio0*aValVtx[jvtx0s*4+1] + (1-ratio0)*aValVtx[jvtx0e*4+1]; 
        vals[2] = ratio0*aValVtx[jvtx0s*4+2] + (1-ratio0)*aValVtx[jvtx0e*4+2]; 
        vals[3] = ratio0*aValVtx[jvtx0s*4+3] + (1-ratio0)*aValVtx[jvtx0e*4+3]; 
      }
      else{
        vs = aVtx[ivtx0];
        vals[0] = aValVtx[ivtx0*4+0];          
        vals[1] = aValVtx[ivtx0*4+1];
        vals[2] = aValVtx[ivtx0*4+2];
        vals[3] = aValVtx[ivtx0*4+3];
      }
      if( aVis[ivtx1].sqdist > 0 ){
        unsigned int jvtx1s = aVis[ivtx1].ivtx0; // visible
        unsigned int jvtx1e = aVis[ivtx1].ivtx1;
        double ratio1 = aVis[ivtx1].ratio;      
        ve   = ratio1*aVtx[   jvtx1s] + (1-ratio1)*aVtx[   jvtx1e];          
        vale[0] = ratio1*aValVtx[jvtx1s*4+0] + (1-ratio1)*aValVtx[jvtx1e*4+0];
        vale[1] = ratio1*aValVtx[jvtx1s*4+1] + (1-ratio1)*aValVtx[jvtx1e*4+1];
        vale[2] = ratio1*aValVtx[jvtx1s*4+2] + (1-ratio1)*aValVtx[jvtx1e*4+2];
        vale[3] = ratio1*aValVtx[jvtx1s*4+3] + (1-ratio1)*aValVtx[jvtx1e*4+3];
      }
      else{
        ve = aVtx[ivtx1]; 
        vale[0] = aValVtx[ivtx1*4+0]; 
        vale[1] = aValVtx[ivtx1*4+1]; 
        vale[2] = aValVtx[ivtx1*4+2]; 
        vale[3] = aValVtx[ivtx1*4+3]; 
      }
    }  
    double s0 = (vs-v0).Length();
    double e0 = (ve-v0).Length();
    double se = (ve-vs).Length();
    const double cos_s0e = (s0*s0+e0*e0-se*se)/(2*s0*e0);
    const double tanh_s0e = sqrt( (1-cos_s0e)/(1+cos_s0e) );
    aLamTnsr[iv0*4+0] += tanh_s0e*(vals[0]/s0+vale[0]/e0);
    aLamTnsr[iv0*4+1] += tanh_s0e*(vals[1]/s0+vale[1]/e0);
    aLamTnsr[iv0*4+2] += tanh_s0e*(vals[2]/s0+vale[2]/e0);
    aLamTnsr[iv0*4+3] += tanh_s0e*(vals[3]/s0+vale[3]/e0);
    sum += tanh_s0e*(1.0/s0+1.0/e0);
  }
  aLamTnsr[iv0*4+0] /= sum;
  aLamTnsr[iv0*4+1] /= sum;
  aLamTnsr[iv0*4+2] /= sum;
  aLamTnsr[iv0*4+3] /= sum;
}


bool CMesher2D_Edit::LambdaLoop_MVC
(const Cad::CCadObj2D& cad_2d, unsigned int id_l, 
 const std::vector<double>& aValVtx,
 std::vector<int>& aflg_ismov)
{
	for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl.ShiftChildLoop()){
    for(itrl.Begin();!itrl.IsEnd();itrl++){
      unsigned int id_e0;   bool is_same_dir0;
      itrl.GetIdEdge(id_e0,is_same_dir0);
      unsigned int ibarary0, itype0;
      if( !this->FindElemLocType_CadIDType(ibarary0,itype0,id_e0,Cad::EDGE) ){ return false; }
      assert( itype0 == 1 );
      assert( ibarary0 < m_aBarAry.size() );
      const std::vector<SBar>& aBar = m_aBarAry[ibarary0].m_aBar;
      for(unsigned int ibar=0;ibar<aBar.size();ibar++){
        const unsigned int iv0 = aBar[ibar].v[0];
        const unsigned int iv1 = aBar[ibar].v[1];
        aflg_ismov[iv0] = 1;
        aflg_ismov[iv1] = 1;
      }      
    }
	}
	std::vector< Com::CVector2D > aVtxLoop;
	std::vector<double> aValVtxLoop;  
  unsigned int nloop = 0;  
  std::vector< unsigned int > aIndLoop;    
  aIndLoop.push_back(0);
	// get coordinate of verticies around loop
	for(Cad::CBRepSurface::CItrLoop itrl=cad_2d.GetItrLoop(id_l);!itrl.IsEnd();itrl.ShiftChildLoop()){
    for(itrl.Begin();!itrl.IsEnd();itrl++){            
      unsigned int id_v = itrl.GetIdVertex();
      assert( cad_2d.IsElemID(Cad::VERTEX,id_v) );
      aVtxLoop.push_back( cad_2d.GetVertexCoord(id_v) );
      aValVtxLoop.push_back( aValVtx[id_v*4+0] );
      aValVtxLoop.push_back( aValVtx[id_v*4+1] );
      aValVtxLoop.push_back( aValVtx[id_v*4+2] );
      aValVtxLoop.push_back( aValVtx[id_v*4+3] );      
      unsigned int id_e; bool is_same_dir;
      if( !itrl.GetIdEdge(id_e,is_same_dir) ) continue;
      if( cad_2d.GetEdgeCurveType(id_e) == 0 ){ continue; } // line
      else if( cad_2d.GetEdgeCurveType(id_e) == 1 ){ continue; assert(0); } // arc
      else if( cad_2d.GetEdgeCurveType(id_e) == 2 ){ // polyline
        unsigned int id_vs = cad_2d.GetIdVertex_Edge(id_e,true );
        unsigned int id_ve = cad_2d.GetIdVertex_Edge(id_e,false);
        Com::CVector2D po_s = cad_2d.GetVertexCoord(id_vs);
        Com::CVector2D po_e = cad_2d.GetVertexCoord(id_ve);        
        std::vector<double> aRelCo;
        cad_2d.GetCurve_Polyline(id_e,aRelCo);
        unsigned int nno = aRelCo.size()/2;
        for(unsigned int ino=0;ino<nno;ino++){
          double hr,vr;
          if( is_same_dir ){
            hr = aRelCo[ino*2+0];
            vr = aRelCo[ino*2+1];
          }
          else{
            hr = aRelCo[(nno-1-ino)*2+0];
            vr = aRelCo[(nno-1-ino)*2+1];
          }
          double d1 = aValVtx[id_vs*4+0]*(1-hr) + aValVtx[id_ve*4+0]*hr + (aValVtx[id_vs*4+2]-aValVtx[id_ve*4+2])*vr;
          double d2 = aValVtx[id_vs*4+1]*(1-hr) + aValVtx[id_ve*4+1]*hr + (aValVtx[id_vs*4+3]-aValVtx[id_ve*4+3])*vr;
          double d3 = aValVtx[id_vs*4+2]*(1-hr) + aValVtx[id_ve*4+2]*hr - (aValVtx[id_vs*4+0]-aValVtx[id_ve*4+0])*vr;
          double d4 = aValVtx[id_vs*4+3]*(1-hr) + aValVtx[id_ve*4+3]*hr - (aValVtx[id_vs*4+1]-aValVtx[id_ve*4+1])*vr;
          aValVtxLoop.push_back( d1 );
          aValVtxLoop.push_back( d2 );
          aValVtxLoop.push_back( d3 );
          aValVtxLoop.push_back( d4 );                
          const Com::CVector2D& vh = po_e - po_s;
          const Com::CVector2D vv(-vh.y, vh.x);
          Com::CVector2D v0 = po_s + vh*hr + vv*vr;
          aVtxLoop.push_back(v0);
        }
      }
    }      
    nloop++;
    aIndLoop.push_back( aVtxLoop.size() );
	}
	unsigned int itriary;
	{
		unsigned int itype;
		if( !this->FindElemLocType_CadIDType(itriary,itype,id_l,Cad::LOOP) ){ return false; }
		assert( itype == 2 );
	}
	std::vector<double> buff;
	buff.resize(aVtxLoop.size(),0);  
	assert( itriary < m_aTriAry.size() );
	const std::vector<STri2D>& aTri = m_aTriAry[itriary].m_aTri;
	for(unsigned int itri=0;itri<aTri.size();itri++){
	for(unsigned int inotri=0;inotri<3;inotri++){
		unsigned int iv0 = aTri[itri].v[inotri];
		if( aflg_ismov[iv0] == 1 ) continue;
		assert( iv0 < aVec2D.size() );
		const CVector2D& v0 = aVec2D[iv0];
		CalcWeight_PMVC(aLamTnsr,iv0,  nloop,aIndLoop, aVtxLoop,aValVtxLoop, v0);
//    CalcWeight_MVC(aLamTnsr,iv0,  nloop, aIndLoop, aVtxLoop, aValVtxLoop, v0);
		aflg_ismov[iv0] = 1;
	}
	}
	return true;
}


void SetBlendedValue
(std::vector<double>& au,
 unsigned int iu,
 double ratio, double hr,
 const std::vector<double>& av, 
 unsigned int iv1, unsigned int iv2)
{
  au[iu*4+0] = av[iv1*4+0]*(1-ratio) + av[iv2*4+0]*ratio + (av[iv1*4+2]-av[iv2*4+2])*hr;
  au[iu*4+1] = av[iv1*4+1]*(1-ratio) + av[iv2*4+1]*ratio + (av[iv1*4+3]-av[iv2*4+3])*hr;
  au[iu*4+2] = av[iv1*4+2]*(1-ratio) + av[iv2*4+2]*ratio - (av[iv1*4+0]-av[iv2*4+0])*hr;
  au[iu*4+3] = av[iv1*4+3]*(1-ratio) + av[iv2*4+3]*ratio - (av[iv1*4+1]-av[iv2*4+1])*hr;
}

// if( is_same_dir == true ) then the start point of id_e_cad is 1 and end point is 0.
bool CMesher2D_Edit::LambdaEdge
(const Cad::CCadObj2D& cad_2d, unsigned int id_e_cad, 
 const std::vector<double>& aValVtx,
 std::vector<int>& aflg_ismov)
{
	assert( cad_2d.IsElemID(Cad::EDGE,id_e_cad) );

	unsigned int id_v_cad_s = cad_2d.GetIdVertex_Edge(id_e_cad,true );
  unsigned int id_v_cad_e = cad_2d.GetIdVertex_Edge(id_e_cad,false);    
	unsigned int ibarary;
	{
		unsigned int itype;
    if( !this->FindElemLocType_CadIDType(ibarary,itype,id_e_cad,Cad::EDGE) ){ return false; }
    assert( itype == 1 );
	}

	const std::vector<SBar>& aBar = m_aBarAry[ibarary].m_aBar;
	unsigned int iv_s, iv_e;
	{
		unsigned int itype;
		this->FindElemLocType_CadIDType(iv_s,itype,id_v_cad_s,Cad::VERTEX); assert( itype == 0 );
		this->FindElemLocType_CadIDType(iv_e,itype,id_v_cad_e,Cad::VERTEX); assert( itype == 0 );
	}

	{ // assert the adjacent number
		const unsigned int nbar = aBar.size();
		assert( aBar[0].v[0] == iv_s );
		for(unsigned int ibar=0;ibar<nbar-1;ibar++){
			assert( aBar[ibar].v[1] == aBar[ibar+1].v[0] );
		}
		assert( aBar[nbar-1].v[1] == iv_e );
	}
    
  if( cad_2d.GetEdgeCurveType(id_e_cad) == 0 ){  // this edge is line
    const Com::CVector2D& vs = cad_2d.GetVertexCoord(id_v_cad_s);
    const Com::CVector2D& ve = cad_2d.GetVertexCoord(id_v_cad_e);    
    const unsigned int nbar = aBar.size();
    if( aflg_ismov[iv_s] == 0 ){ 
			aflg_ismov[iv_s] = 1;
      SetBlendedValue(aLamTnsr,iv_s, 0,0, aValVtx,id_v_cad_s,id_v_cad_e);
    }
    for(unsigned int ibar=1;ibar<nbar;ibar++){
      unsigned int iv0 = aBar[ibar].v[0];
      if( aflg_ismov[iv0] != 0 ) continue;
			const double len_s = Com::Distance(vs,aVec2D[iv0]);
			const double len_e = Com::Distance(ve,aVec2D[iv0]);
			const double rs = len_s/(len_s+len_e);
      aflg_ismov[iv0] = 1;
      SetBlendedValue(aLamTnsr,iv0 ,rs,0,  aValVtx,id_v_cad_s,id_v_cad_e);      
    }
    if( aflg_ismov[iv_e] == 0 ){
      aflg_ismov[iv_e] = 1;
      SetBlendedValue(aLamTnsr,iv_e, 1,0,  aValVtx,id_v_cad_s,id_v_cad_e);      
    }
  }
  else{ // curved edge
    const Com::CVector2D& vs = cad_2d.GetVertexCoord(id_v_cad_s);
    const Com::CVector2D& ve = cad_2d.GetVertexCoord(id_v_cad_e);    
    const Com::CVector2D& vse = ve-vs;
    const double invlen_se = 1.0/vse.Length();
    const unsigned int nbar = aBar.size();
    if( aflg_ismov[iv_s] == 0 ){ 
			aflg_ismov[iv_s] = 1;
      SetBlendedValue(aLamTnsr,iv_s, 0,0,  aValVtx,id_v_cad_s,id_v_cad_e);
    }
    for(unsigned int ibar=1;ibar<nbar;ibar++){
      unsigned int iv0 = aBar[ibar].v[0];
      if( aflg_ismov[iv0] != 0 ) continue;
      double area = Com::TriArea(vs,ve,aVec2D[iv0]);
      double hr = area*invlen_se*2;
			const double len_s = Com::Dot(vse,ve-aVec2D[iv0]);
			const double len_e = Com::Dot(vse,aVec2D[iv0]-vs);
			const double rs = len_e/(len_s+len_e);
      aflg_ismov[iv0] = 1;
      SetBlendedValue(aLamTnsr,iv0 ,rs,hr,   aValVtx,id_v_cad_s,id_v_cad_e);      
    }
    if( aflg_ismov[iv_e] == 0 ){
      aflg_ismov[iv_e] = 1;
      SetBlendedValue(aLamTnsr,iv_e, 1,0,    aValVtx,id_v_cad_s,id_v_cad_e);      
    }    
  }
  return true;  
}

void CMesher2D_Edit::Precomp_FitMeshToCad(const Cad::CCadObj2D& cad_2d,
//                                          Cad::CAD_ELEM_TYPE itype_elem, unsigned int id_elem, 
                                          const std::vector<double>& aLamTnsr_Vtx)
{
  this->aLamTnsr_Vtx = aLamTnsr_Vtx;
//  this->move_cad_elem_type = itype_elem;
//  this->move_cad_elem_id = id_elem;
  const unsigned int nvec = this->aVec2D.size();
  std::vector<int> aflg_ismov;  
  aflg_ismov.resize(nvec,0);
	aLamTnsr.clear(); 
  this->aLamTnsr.resize(nvec*4,0);  
  {
    const std::vector<unsigned int>& aIdE = cad_2d.GetAryElemID(Cad::EDGE);
    for(unsigned int iie=0;iie<aIdE.size();iie++){
      const unsigned int id_e = aIdE[iie];
      if( !cad_2d.IsElemID(Cad::EDGE,id_e) ) continue;
			this->LambdaEdge(cad_2d,id_e,aLamTnsr_Vtx,aflg_ismov);
    }
  }
  {
    const std::vector<unsigned int>& aIdL = cad_2d.GetAryElemID(Cad::LOOP);
    for(unsigned int iil=0;iil<aIdL.size();iil++){
      const unsigned int id_l = aIdL[iil];
			this->LambdaLoop_MVC(cad_2d,id_l,aLamTnsr_Vtx,aflg_ismov);
    }
  }  
}