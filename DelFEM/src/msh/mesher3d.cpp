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
#include <queue>
#include <cassert>
#include <math.h>
#include <cstdlib> //(rand, RAND_MAX)
#include <string>
#include <iostream>
#include <fstream>

#include "delfem/msh/meshkernel2d.h"
#include "delfem/msh/meshkernel3d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"

using namespace Com;

bool Msh::CMesher3D::ReadFile_PLY( const std::string& file_name )
{
	FILE* fp;
	const unsigned int buff_size = 512;
	char stmp1[buff_size];

	if( (fp = fopen(file_name.c_str(),"r")) == NULL ){
		std::cout << "Error!-->Cannot Open File" << std::endl;
		assert(0);
		return false;
	}

	std::cout << "Ply File Open : " << file_name.c_str() << std::endl;

	fgets(stmp1,buff_size,fp);
	assert( strcmp(stmp1,"ply\n") == 0 );

	unsigned int nface;
	unsigned int nvector;
	for(;;){
		fgets(stmp1,buff_size,fp);
		if( strcmp(stmp1,"end_header\n") == 0 ){ 
			break;
		}
		else if( strncmp(stmp1,"element vertex",14) == 0 ){
			sscanf(stmp1,"%*s%*s%d",&nvector);
		}
		else if( strncmp(stmp1,"element face",12) == 0 ){
			sscanf(stmp1,"%*s%*s%d",&nface);
		}
	}
//	std::cout << nvector << " " << nface << std::endl;
	{	// Read Vertex
		this->aVec.resize(nvector);
		double vx,vy,vz;
		for(unsigned int ivector=0;ivector<nvector;ivector++){
			if( fgets(stmp1,buff_size,fp) == NULL ){ assert(0); }
			sscanf(stmp1,"%lf%lf%lf",&vx,&vy,&vz);
//			std::cout << vx << " " << vy << " " << vz << std::endl;
			aVec[ivector].x = vx;
			aVec[ivector].y = vy;
			aVec[ivector].z = vz;
		}
	}
	{	// Read Face
		unsigned int i1,i2,i3,i4;
		unsigned int ntri=0, nquad=0;
		unsigned int* tmp_buffer = new unsigned int [nface*5];
		for(unsigned int iface=0;iface<nface;iface++){
			if( fgets(stmp1,buff_size,fp) == NULL ){ assert(0); }
			if( stmp1[0] == '3' ){
				sscanf(stmp1,"%*d%d%d%d",&i1,&i2,&i3);
				tmp_buffer[iface*5  ] = 3;
				tmp_buffer[iface*5+1] = i1;
				tmp_buffer[iface*5+2] = i2;
				tmp_buffer[iface*5+3] = i3;
				ntri++;
			}
			else if( stmp1[0] == '4' ){
				sscanf(stmp1,"%*d%d%d%d%d",&i1,&i2,&i3,&i4);
				tmp_buffer[iface*5  ] = 4;
				tmp_buffer[iface*5+1] = i1;
				tmp_buffer[iface*5+2] = i2;
				tmp_buffer[iface*5+3] = i3;
				tmp_buffer[iface*5+4] = i4;
				nquad++;
			}
			else{ assert(0); }
		}
//		std::cout << nquad << " " << ntri << std::endl;
		{
			// TODO : Quadにも対応(段階を踏むのがいいかも．まずQUADを２分割してTRIで対応する）
			this->m_aTriAry.resize(1);
			m_aTriAry[0].id = 1;
			m_aTriAry[0].m_aTri.resize(ntri+nquad*2);
			unsigned int itri0 = 0;
			for(unsigned int iface=0;iface<nface;iface++){
				if( tmp_buffer[iface*5] == 3 ){
					m_aTriAry[0].m_aTri[itri0].v[0] = tmp_buffer[iface*5+1];
					m_aTriAry[0].m_aTri[itri0].v[1] = tmp_buffer[iface*5+2];
					m_aTriAry[0].m_aTri[itri0].v[2] = tmp_buffer[iface*5+3];
					itri0++;
				}
				else if( tmp_buffer[iface*5] == 4 ){
					m_aTriAry[0].m_aTri[itri0].v[0] = tmp_buffer[iface*5+1];
					m_aTriAry[0].m_aTri[itri0].v[1] = tmp_buffer[iface*5+2];
					m_aTriAry[0].m_aTri[itri0].v[2] = tmp_buffer[iface*5+4];
					itri0++;
					m_aTriAry[0].m_aTri[itri0].v[0] = tmp_buffer[iface*5+2];
					m_aTriAry[0].m_aTri[itri0].v[1] = tmp_buffer[iface*5+3];
					m_aTriAry[0].m_aTri[itri0].v[2] = tmp_buffer[iface*5+4];
					itri0++;
				}
			}
		}
		delete[] tmp_buffer;
	}
	fclose(fp);
	this->MakeElemLocationType();
	this->EraseSameLocationPointOnSurface();
	return true;
}


bool Msh::CMesher3D::HomogenizeSurface(double elen)
{
	assert( elen > 0.0 );
	assert( m_aTriAry.size() == 1 );
	std::vector<STri3D>& aTri = m_aTriAry[0].m_aTri;
	/*
	COctTree oct_tree;
	oct_tree.SetBoundingBox(aVec);
	for(unsigned int ivec=0;ivec<aVec.size();ivec++){
		const int iret = oct_tree.InsertPoint( ivec,aVec );
		assert( iret == -1 );
	}
	*/
	aVec.reserve( aVec.size()+aTri.size()*10 );
	
	for(unsigned int iloop_out=0;iloop_out<10;iloop_out++){
	for(unsigned int itri=0;itri<aTri.size();itri++){
		const double area = Com::TriArea( aTri[itri].v[0], aTri[itri].v[1], aTri[itri].v[2], aVec );
		if( area < elen*elen*0.5*2.5 ){
//			std::cout << "A : " << itri << " " << area << " " << elen*elen*0.5*3.0 << std::endl;
			continue;
		}
		unsigned int itried;	// 長くて特徴線で三角形が扁平な辺を探す
		for(itried=0;itried<Msh::nEdTri;itried++){
			const unsigned int ipo0 = aTri[itri].v[ itried ];
			const unsigned int ipo1 = aTri[itri].v[ noelTriEdge[itried][0] ];
			const unsigned int ipo2 = aTri[itri].v[ noelTriEdge[itried][1] ];
			const double len0 = Com::Distance( aVec[ipo1], aVec[ipo2] );
			if( len0 < elen*1.4 ){
//				std::cout << "B : " << len0 << " " << elen << std::endl;
				continue;
			}
/*			const double len1 = Length( aVec[ipo0], aVec[ipo2] );
			const double len2 = Length( aVec[ipo0], aVec[ipo1] );
			if( len0 < (len1+len2)*0.8 ){
//				std::cout << "C : " << len0 << " " << len1 << " " << len2 << std::endl;
				continue;
			}*/
			unsigned int ipo_d;
			{
				const unsigned int itri_d = aTri[itri].s2[itried];
				const unsigned int* rel_d = relTriTri[ aTri[itri].r2[itried] ];
				const unsigned int itried_d = rel_d[itried];
				ipo_d = aTri[itri_d].v[itried_d];
			}
			{
				Com::CVector3D norm1, norm2;
				Com::UnitNormal(norm1,aVec[ipo0],aVec[ipo1],aVec[ipo2]);
				Com::UnitNormal(norm2,aVec[ipo_d],aVec[ipo2],aVec[ipo1]);
				const double dot = Com::Dot(norm1,norm2);
				if( dot > 1.0-1.0e-5 ){
//					std::cout << "D : " << dot << std::endl;
					continue;
				}
			}
			break;
		}
		if( itried != Msh::nEdTri ){	// 辺に点を加える
			const unsigned int ipo1 = aTri[itri].v[ noelTriEdge[itried][0] ];
			const unsigned int ipo2 = aTri[itri].v[ noelTriEdge[itried][1] ];
			const double ins_x = aVec[ipo1].x*0.5+aVec[ipo2].x*0.5;
			const double ins_y = aVec[ipo1].y*0.5+aVec[ipo2].y*0.5;
			const double ins_z = aVec[ipo1].z*0.5+aVec[ipo2].z*0.5;
			CVector3D po_ins(ins_x,ins_y,ins_z);
			const unsigned int ivec_add = aVec.size();
			aVec.push_back( po_ins );
			Msh::AddPointTri_Edge(ivec_add,itri,itried,aTri);
//			assert( Msh::CheckTri(aTri) );
			assert( aTri[itri].v[0] == ivec_add );
			Msh::DelaunayAroundPointTri(itri,0,aVec,aTri);
//			assert( Msh::CheckTri(aTri) );
		}
		else{	// 面の中心に点を加える
			const unsigned int ipo0 = aTri[itri].v[0];
			const unsigned int ipo1 = aTri[itri].v[1];
			const unsigned int ipo2 = aTri[itri].v[2];
			const double ins_x = aVec[ipo0].x/3.0+aVec[ipo1].x/3.0+aVec[ipo2].x/3.0;
			const double ins_y = aVec[ipo0].y/3.0+aVec[ipo1].y/3.0+aVec[ipo2].y/3.0;
			const double ins_z = aVec[ipo0].z/3.0+aVec[ipo1].z/3.0+aVec[ipo2].z/3.0;
			CVector3D po_ins(ins_x,ins_y,ins_z);
			const unsigned int ivec_add = aVec.size();
			aVec.push_back( po_ins );
			Msh::AddPointTri_Face(ivec_add,itri,aTri);
//			assert( Msh::CheckTri(aTri) );
			assert( aTri[itri].v[0] == ivec_add );
			Msh::DelaunayAroundPointTri(itri,0,aVec,aTri);
//			assert( Msh::CheckTri(aTri) );
		}
	}
	}

	////////////////////////////////////////////////////////////////

	assert( Msh::CheckTri(aTri) );

	std::vector<int> vec_flg;
	vec_flg.reserve( aVec.size()*2 );
	vec_flg.resize( aVec.size(), 0 );
	for(unsigned int iloop=1;iloop<2;iloop++){
		for(unsigned int itri=0;itri<aTri.size();itri++){
		for(unsigned int inotri=0;inotri<Msh::nNoTri;inotri++){
			const unsigned int ivec0 = aTri[itri].v[inotri];
			assert( ivec0<aVec.size() );
			if( vec_flg[ivec0] == (int)iloop || vec_flg[ivec0] == -1 ) continue;
			unsigned int ntri_around;
			bool is_same_norm = true;
			{	// このivec0の周りの要素の法線が同じかどうか調べる
				CVector3D norm_i;
				Msh::UnitNormal(norm_i,itri,aTri,aVec);
				unsigned int jtri0 = itri;
				unsigned int jnotri0 = inotri;
				ntri_around = 1;
				for(;;){
					unsigned int jtri_nex = aTri[jtri0].s2[ noelTriEdge[jnotri0][0] ];
					const unsigned int* rel = relTriTri[ aTri[jtri0].r2[ noelTriEdge[jnotri0][0] ] ];
					unsigned int jnotri_nex = rel[ jnotri0 ];
					assert( aTri[jtri_nex].v[jnotri_nex] == ivec0 );
					if( jtri_nex == itri ) break;
					jtri0 = jtri_nex;
					jnotri0 = jnotri_nex;
					{	// ここからitriを含まないivec0を囲む要素が順にjtri0に入るはず
						ntri_around++;
						CVector3D norm_j;
						Msh::UnitNormal(norm_j,jtri0,aTri,aVec);
						const double dot = Dot(norm_i,norm_j);
						if( dot < 1.0-1.0e-10 ){
							is_same_norm = false;
							break;
						}
					}
				}
			}
			if( !is_same_norm ){
				vec_flg[ivec0] = -1;
				continue;
			}

			CVector3D vec_delta;
			{	// Bossen Heckbertの方法による平滑化
				unsigned int jtri0 = itri;
				unsigned int jnotri0 = inotri;
				for(;;){
					{	// ここからjtri0はivec0を囲む全ての三角形になるはず
						const unsigned int ivec1 = aTri[jtri0].v[ noelTriEdge[jnotri0][0] ];
						double coeff;
//						if( iloop < 10 ){	
							coeff = 1.0 / ntri_around;
/*						}
						else{
							const double dist = Length( aVec[ivec0], aVec[ivec1] );
							const double ratio = dist / elen;
							const double dtmp1 = ratio*ratio*ratio*ratio;
							coeff = 0.2/dist*(1.0-dtmp1)*exp(-dtmp1);
						}*/
						vec_delta += (aVec[ivec1] + -1.0*aVec[ivec0])*coeff;
					}
					unsigned int jtri_nex = aTri[jtri0].s2[ noelTriEdge[jnotri0][0] ];
					const unsigned int* rel = relTriTri[ aTri[jtri0].r2[ noelTriEdge[jnotri0][0] ] ];
					unsigned int jnotri_nex = rel[ jnotri0 ];
					assert( aTri[jtri_nex].v[jnotri_nex] == ivec0 );
					if( jtri_nex == itri ) break;
					jtri0 = jtri_nex;
					jnotri0 = jnotri_nex;
				}
			}

			bool is_tri_ok_flg = true;
			{
				CVector3D vec_new = aVec[ivec0] + vec_delta;
				CVector3D norm_i;
				Com::Normal(norm_i, vec_new, aVec[aTri[itri].v[1]], aVec[aTri[itri].v[2]]);
				unsigned int jtri0 = itri;
				unsigned int jnotri0 = inotri;
				for(;;){
					{	// ここからjtri0はivec0を囲む全ての三角形になるはず
						const unsigned int ivec1 = aTri[jtri0].v[ noelTriEdge[jnotri0][0] ];
						const unsigned int ivec2 = aTri[jtri0].v[ noelTriEdge[jnotri0][1] ];
//						const double dist = Length( aVec[ivec0], aVec[ivec1] );
						CVector3D norm1;
						Com::Normal(norm1, vec_new,aVec[ivec1],aVec[ivec2]);
						double dot = Dot(norm1,norm_i);
						if( dot < 0.0 ){
							is_tri_ok_flg = false;
							break;
						}
					}
					unsigned int jtri_nex = aTri[jtri0].s2[ noelTriEdge[jnotri0][0] ];
					const unsigned int* rel = relTriTri[ aTri[jtri0].r2[ noelTriEdge[jnotri0][0] ] ];
					unsigned int jnotri_nex = rel[ jnotri0 ];
					assert( aTri[jtri_nex].v[jnotri_nex] == ivec0 );
					if( jtri_nex == itri ) break;
					jtri0 = jtri_nex;
					jnotri0 = jnotri_nex;
				}
			}

			vec_flg[ivec0] = iloop;

			if( is_tri_ok_flg ){
				// 節点の値を更新
				aVec[ivec0] += vec_delta;
				Msh::DelaunayAroundPointTri(itri,inotri,aVec,aTri);
//				assert( Msh::CheckTri(aTri) );
				{
					unsigned int jnotri;
					for(jnotri=0;jnotri<3;jnotri++){
						if( aTri[itri].v[jnotri] == ivec0 ){ break; }
					}
					assert( jnotri != 3 );
					inotri = jnotri;
				}
			}
			else{
				std::cout << "not ok " << ivec0 << std::endl;
			}

			ntri_around = 0;
			double node_extent = 0.0;
			double max_edge_extent = 0.0;
            unsigned int itri_max=0, inotri_max=0;
			{
				CVector3D vec_new = aVec[ivec0] + vec_delta;
				CVector3D norm_i;
				Com::Normal(norm_i, vec_new, aVec[aTri[itri].v[1]], aVec[aTri[itri].v[2]]);
				unsigned int jtri0 = itri;
				unsigned int jnotri0 = inotri;
				for(;;){
					unsigned int jtri_nex = aTri[jtri0].s2[ noelTriEdge[jnotri0][0] ];
					const unsigned int* rel = relTriTri[ aTri[jtri0].r2[ noelTriEdge[jnotri0][0] ] ];
					unsigned int jnotri_nex = rel[ jnotri0 ];
					assert( aTri[jtri_nex].v[jnotri_nex] == ivec0 );
					{	// ここからjtri0はivec0を囲む全ての三角形になるはず
						const double area_j = Com::TriArea( aVec[aTri[jtri0].v[0]], aVec[aTri[jtri0].v[1]], aVec[aTri[jtri0].v[2]] );
						node_extent += area_j;
						ntri_around += 1;
						const double area_nex = Com::TriArea( aVec[aTri[jtri_nex].v[0]], aVec[aTri[jtri_nex].v[1]], aVec[aTri[jtri_nex].v[2]] );
						const double edge_extent = 0.8*(area_j+area_nex)/(elen*elen);
						if( max_edge_extent < edge_extent ){
							itri_max = jtri0;
							inotri_max = noelTriEdge[jnotri0][0];
							max_edge_extent = edge_extent;
						}
					}
					if( jtri_nex == itri ) break;
					jtri0 = jtri_nex;
					jnotri0 = jnotri_nex;
				}
				node_extent = 2.3 * node_extent / ( elen*elen*sqrt(ntri_around*(ntri_around-2.0)) );
			}

			std::cout << node_extent << " " << max_edge_extent << std::endl;
			if( node_extent < 1.0 ){
				std::cout << "****************" << std::endl;
			}

//			continue;

			if( node_extent < 1.0 ){	// Delete this node

			}
			else if( max_edge_extent > 1.0 ){	// Add Node to Edgde
				unsigned int ivec_ins;
				{
					ivec_ins = aVec.size();
					vec_flg.resize( aVec.size()+1,iloop);
					aVec.resize( aVec.size()+1 );
					const unsigned int ivec1 = aTri[itri_max].v[ noelTriEdge[inotri_max][0] ];
					assert( ivec0 == aTri[itri_max].v[ noelTriEdge[inotri_max][1] ] );
					aVec[ivec_ins].x = aVec[ivec0].x*0.5 + aVec[ivec1].x*0.5;
					aVec[ivec_ins].y = aVec[ivec0].y*0.5 + aVec[ivec1].y*0.5;
					aVec[ivec_ins].z = aVec[ivec0].z*0.5 + aVec[ivec1].z*0.5;
				}
				Msh::AddPointTri_Edge( ivec_ins, itri_max,inotri_max, aTri );
//				assert( Msh::CheckTri(aTri) );
				Msh::DelaunayAroundPointTri(itri_max,0,aVec,aTri);
//				assert( Msh::CheckTri(aTri) );
			}
		}
		}
	}
	return true;
}




bool Msh::CMesher3D::EraseSameLocationPointOnSurface()
{	
	std::cout << "EraseSameLocationPointOnSurface" << std::endl;
	const unsigned int nvec_old = aVec.size();

	std::vector<int> vec_old2new;
	std::vector<int> vec_old2old;
	vec_old2new.resize( nvec_old );
	vec_old2old.resize( nvec_old );

	unsigned int nvec_new = 0;
	COctTree oct_tree;
	{
		double x_min,x_max,  y_min,y_max,  z_min,z_max;
		{
			const unsigned int nvec = aVec.size();	assert( nvec > 0 );
			x_min = aVec[0].x; x_max = x_min;
			y_min = aVec[0].y; y_max = y_min;
			z_min = aVec[0].z; z_max = z_min;
			for(unsigned int ivec=1;ivec<nvec;ivec++){
				x_min = ( aVec[ivec].x < x_min ) ? aVec[ivec].x : x_min;
				x_max = ( aVec[ivec].x > x_max ) ? aVec[ivec].x : x_max;
				y_min = ( aVec[ivec].y < y_min ) ? aVec[ivec].y : y_min;
				y_max = ( aVec[ivec].y > y_max ) ? aVec[ivec].y : y_max;
				z_min = ( aVec[ivec].z < z_min ) ? aVec[ivec].z : z_min;
				z_max = ( aVec[ivec].z > z_max ) ? aVec[ivec].z : z_max;
			}
			double x_off = (x_max - x_min)*0.01;
			x_min -= x_off;	x_max += x_off;
			double y_off = (y_max - y_min)*0.01;
			y_min -= y_off;	y_max += y_off;
			double z_off = (z_max - z_min)*0.01;
			z_min -= z_off;	z_max += z_off;
		}
		oct_tree.SetBoundingBox( CBoundingBox3D(x_min,x_max, y_min,y_max, z_min,z_max) );
	}
	for(unsigned int ivec=0;ivec<aVec.size();ivec++){
//		std::cout << ivec << " " << aVec.size() << std::endl;
		const int iret = oct_tree.InsertPoint( ivec, aVec[ivec] );
		assert( iret != -2 );
		if( iret >= 0 ){ 
			assert( (unsigned int)iret < aVec.size() );
			vec_old2old[ivec] = iret; 
			vec_old2new[ivec] = -1;
		}
		else{ 
			assert( iret == -1 );
			vec_old2old[ivec] = ivec; 
			vec_old2new[ivec] = nvec_new;
			nvec_new++;
		}
	}

//	std::cout << nvec_new << std::endl;

	std::vector<unsigned int> vec_new2old;
	vec_new2old.resize( nvec_new );
	for(unsigned int ivec=0;ivec<nvec_old;ivec++){
		if( vec_old2new[ivec] == -1 ) continue;
		const unsigned int ivec_new = vec_old2new[ivec];
		assert( ivec_new < nvec_new );
		vec_new2old[ivec_new] = ivec;
	}

	assert( m_aTriAry.size() == 1 );
	for(unsigned int itri=0;itri<m_aTriAry[0].m_aTri.size();itri++){
	for(unsigned int inotri=0;inotri<3;inotri++){
		unsigned int ivec_old0 = m_aTriAry[0].m_aTri[itri].v[inotri];
		unsigned int ivec_old1 = vec_old2old[ ivec_old0 ];
		unsigned int ivec_new = vec_old2new[ ivec_old1 ];
		m_aTriAry[0].m_aTri[itri].v[inotri] = ivec_new;
	}
	}

	for(unsigned int ivec=0;ivec<nvec_new;ivec++){
		unsigned int ivec_old = vec_new2old[ivec];
		assert( ivec <= ivec_old );
		aVec[ivec] = aVec[ivec_old];
	}
	aVec.resize(nvec_new);

	Msh::MakeTriSurTri(m_aTriAry[0].m_aTri);
	assert( Msh::CheckTri(m_aTriAry[0].m_aTri) );

	return true;
}




bool Msh::CMesher3D::CutMesh(double elen)
{
	assert( elen > 0.0 );
	//
	// ここで表面を細かく切るルーティンを呼ぶ
	//
	std::vector<CPoint3D> aPo3D;
	std::vector<int> vec2po;
	{
		// 要素分割する領域の節点　aPo2Dを作成
		// 面に属する節点の全体番号から、要素分割する領域のローカル番号への対応(vec2po)を作成
		////////////////////////////////
		vec2po.resize( aVec.size(), -1 );
		{	// vec2poを作る、aPo2Dを確保する
//			int ipo=0;
			for(unsigned int ivec=0;ivec<aVec.size();ivec++){
				vec2po[ivec] = ivec+8;
			}
			aPo3D.resize( aVec.size()+8 );
		}
		for(unsigned int ivec=0;ivec<aVec.size();ivec++){
			const unsigned int ipoin = vec2po[ivec];
			assert( ipoin < aPo3D.size() );
			aPo3D[ipoin].p.x = aVec[ivec].x;
			aPo3D[ipoin].p.y = aVec[ivec].y;
			aPo3D[ipoin].p.z = aVec[ivec].z;
			aPo3D[ipoin].e = -1;
			aPo3D[ipoin].poel = 0;
		}
	}

	COctTree oct_tree;
	{
		double x_min,x_max,  y_min,y_max,  z_min,z_max;
		{
			const unsigned int nvec = aVec.size();	assert( nvec > 0 );
			x_min = aPo3D[8].p.x; x_max = x_min;
			y_min = aPo3D[8].p.y; y_max = y_min;
			z_min = aPo3D[8].p.z; z_max = z_min;
			for(unsigned int ivec=9;ivec<nvec;ivec++){
				x_min = ( aPo3D[ivec].p.x < x_min ) ? aPo3D[ivec].p.x : x_min;
				x_max = ( aPo3D[ivec].p.x > x_max ) ? aPo3D[ivec].p.x : x_max;
				y_min = ( aPo3D[ivec].p.y < y_min ) ? aPo3D[ivec].p.y : y_min;
				y_max = ( aPo3D[ivec].p.y > y_max ) ? aPo3D[ivec].p.y : y_max;
				z_min = ( aPo3D[ivec].p.z < z_min ) ? aPo3D[ivec].p.z : z_min;
				z_max = ( aPo3D[ivec].p.z > z_max ) ? aPo3D[ivec].p.z : z_max;
			}
			double x_off = (x_max - x_min)*0.01;
			x_min -= x_off;	x_max += x_off;
			double y_off = (y_max - y_min)*0.01;
			y_min -= y_off;	y_max += y_off;
			double z_off = (z_max - z_min)*0.01;
			z_min -= z_off;	z_max += z_off;
		}
//		std::cout << x_min << " " << x_max << "  " << y_min << " " << y_max << "  " << z_min << " " << z_max << std::endl;
		oct_tree.SetBoundingBox( CBoundingBox3D(x_min,x_max, y_min,y_max, z_min,z_max) );
		for(unsigned int ipoin=8;ipoin<aPo3D.size();ipoin++){
			oct_tree.InsertPoint(ipoin,aPo3D[ipoin].p);
		}
	}

	std::vector<STet> aTet;
	{	// 与えられた点群を含む大きな四面体を作る
		double x_min,x_max,  y_min,y_max,  z_min,z_max;
		{
			const unsigned int nvec = aVec.size();	assert( nvec > 0 );
			x_min = aVec[0].x; x_max = x_min;
			y_min = aVec[0].y; y_max = y_min;
			z_min = aVec[0].z; z_max = z_min;
			for(unsigned int ivec=1;ivec<nvec;ivec++){
				x_min = ( aVec[ivec].x < x_min ) ? aVec[ivec].x : x_min;
				x_max = ( aVec[ivec].x > x_max ) ? aVec[ivec].x : x_max;
				y_min = ( aVec[ivec].y < y_min ) ? aVec[ivec].y : y_min;
				y_max = ( aVec[ivec].y > y_max ) ? aVec[ivec].y : y_max;
				z_min = ( aVec[ivec].z < z_min ) ? aVec[ivec].z : z_min;
				z_max = ( aVec[ivec].z > z_max ) ? aVec[ivec].z : z_max;
			}
			double x_off = (x_max - x_min)*0.5;
			x_min -= x_off;	x_max += x_off;
			double y_off = (y_max - y_min)*0.5;
			y_min -= y_off;	y_max += y_off;
			double z_off = (z_max - z_min)*0.5;
			z_min -= z_off;	z_max += z_off;
		}

		aPo3D[0].p.x = x_min;	aPo3D[0].p.y = y_min;	aPo3D[0].p.z = z_min;	aPo3D[0].e = 0;	aPo3D[0].poel = 1;
		aPo3D[1].p.x = x_min;	aPo3D[1].p.y = y_min;	aPo3D[1].p.z = z_max;	aPo3D[1].e = 0;	aPo3D[1].poel = 0;
		aPo3D[2].p.x = x_min;	aPo3D[2].p.y = y_max;	aPo3D[2].p.z = z_min;	aPo3D[2].e = 5;	aPo3D[2].poel = 0;
		aPo3D[3].p.x = x_min;	aPo3D[3].p.y = y_max;	aPo3D[3].p.z = z_max;	aPo3D[3].e = 3;	aPo3D[3].poel = 3;
		aPo3D[4].p.x = x_max;	aPo3D[4].p.y = y_min;	aPo3D[4].p.z = z_min;	aPo3D[4].e = 0;	aPo3D[4].poel = 3;
		aPo3D[5].p.x = x_max;	aPo3D[5].p.y = y_min;	aPo3D[5].p.z = z_max;	aPo3D[5].e = 2;	aPo3D[5].poel = 0;
		aPo3D[6].p.x = x_max;	aPo3D[6].p.y = y_max;	aPo3D[6].p.z = z_min;	aPo3D[6].e = 3;	aPo3D[6].poel = 0;
		aPo3D[7].p.x = x_max;	aPo3D[7].p.y = y_max;	aPo3D[7].p.z = z_max;	aPo3D[7].e = 0;	aPo3D[7].poel = 2;

		aTet.resize(6);
		aTet[0].v[0] = 1;	aTet[0].v[1] = 0;	aTet[0].v[2] = 7;	aTet[0].v[3] = 4;
		aTet[0].g[0] = -2;	aTet[0].g[1] = -2;	aTet[0].g[2] = -1;	aTet[0].g[3] = -2;
		aTet[0].s[0] =  1;	aTet[0].s[1] =  2;	aTet[0].s[2] =  0;	aTet[0].s[3] =  4;
		aTet[0].f[0] =  2;	aTet[0].f[1] =  3;	aTet[0].f[2] =  0;	aTet[0].f[3] =  2;
		aTet[1].v[0] = 6;	aTet[1].v[1] = 7;	aTet[1].v[2] = 0;	aTet[1].v[3] = 4;
		aTet[1].g[0] = -2;	aTet[1].g[1] = -1;	aTet[1].g[2] = -1;	aTet[1].g[3] = -2;
		aTet[1].s[0] =  0;	aTet[1].s[1] =  0;	aTet[1].s[2] =  0;	aTet[1].s[3] =  3;
		aTet[1].f[0] =  2;	aTet[1].f[1] =  0;	aTet[1].f[2] =  0;	aTet[1].f[3] =  2;
		aTet[2].v[0] = 5;	aTet[2].v[1] = 1;	aTet[2].v[2] = 7;	aTet[2].v[3] = 4;
		aTet[2].g[0] = -2;	aTet[2].g[1] = -1;	aTet[2].g[2] = -1;	aTet[2].g[3] = -1;
		aTet[2].s[0] =  0;	aTet[2].s[1] =  0;	aTet[2].s[2] =  0;	aTet[2].s[3] =  0;
		aTet[2].f[0] =  3;	aTet[2].f[1] =  0;	aTet[2].f[2] =  0;	aTet[2].f[3] =  0;
		aTet[3].v[0] = 6;	aTet[3].v[1] = 0;	aTet[3].v[2] = 7;	aTet[3].v[3] = 3;
		aTet[3].g[0] = -2;	aTet[3].g[1] = -1;	aTet[3].g[2] = -2;	aTet[3].g[3] = -2;
		aTet[3].s[0] =  4;	aTet[3].s[1] =  0;	aTet[3].s[2] =  5;	aTet[3].s[3] =  1;
		aTet[3].f[0] =  2;	aTet[3].f[1] =  0;	aTet[3].f[2] =  8;	aTet[3].f[3] =  2;
		aTet[4].v[0] = 1;	aTet[4].v[1] = 7;	aTet[4].v[2] = 0;	aTet[4].v[3] = 3;
		aTet[4].g[0] = -2;	aTet[4].g[1] = -1;	aTet[4].g[2] = -1;	aTet[4].g[3] = -2;
		aTet[4].s[0] =  3;	aTet[4].s[1] =  0;	aTet[4].s[2] =  0;	aTet[4].s[3] =  0;
		aTet[4].f[0] =  2;	aTet[4].f[1] =  0;	aTet[4].f[2] =  0;	aTet[4].f[3] =  2;
		aTet[5].v[0] = 2;	aTet[5].v[1] = 0;	aTet[5].v[2] = 6;	aTet[5].v[3] = 3;
		aTet[5].g[0] = -2;	aTet[5].g[1] = -1;	aTet[5].g[2] = -1;	aTet[5].g[3] = -1;
		aTet[5].s[0] =  3;	aTet[5].s[1] =  0;	aTet[5].s[2] =  0;	aTet[5].s[3] =  0;
		aTet[5].f[0] =  8;	aTet[5].f[1] =  0;	aTet[5].f[2] =  0;	aTet[5].f[3] =  0;

		assert( Msh::CheckTet(aTet,aPo3D) );
	}

	for(unsigned int ipoin=0;ipoin<aPo3D.size();ipoin++){
		if( aPo3D[ipoin].e >= 0 ) continue; 
		for(unsigned int itet=0;itet<aTet.size();itet++){
			
			unsigned int inum_loc = 0;
			{	
				const double vol = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], aTet[itet].v[2], aTet[itet].v[3], aPo3D);
				const double tol0 = 1.0e-20;
				const double ratio0 = Msh::TetVolume( ipoin, aTet[itet].v[1], aTet[itet].v[2], aTet[itet].v[3], aPo3D)/vol;
				if( ratio0 < -tol0 ) continue;
				const double ratio1 = Msh::TetVolume( aTet[itet].v[0], ipoin, aTet[itet].v[2], aTet[itet].v[3], aPo3D)/vol;
				if( ratio1 < -tol0 ) continue;
				const double ratio2 = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], ipoin, aTet[itet].v[3], aPo3D)/vol;
				if( ratio2 < -tol0 ) continue;
				const double ratio3 = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], aTet[itet].v[2], ipoin, aPo3D)/vol;
				if( ratio3 < -tol0 ) continue;

				const double tol1 = (1.0e-10/vol < 1.0e-5) ? 1.0e-10/vol : 1.0e-5;
//				std::cout << tol1 << std::endl;
				if( ratio0 < tol1 ) inum_loc += 1;
				if( ratio1 < tol1 ) inum_loc += 2;
				if( ratio2 < tol1 ) inum_loc += 4;
				if( ratio3 < tol1 ) inum_loc += 8;
//				std::cout << ipoin << " " << itet << " " << ratio0 << " " << ratio1 << " " << ratio2 << " " << ratio3 << std::endl;
			}

			if( inum_loc == 0 ){
				std::cout << "InsertElem " << std::endl;
				Msh::AddPointTet_Elem(itet,ipoin,aPo3D,aTet);
			}
			else if( inum_loc == 1 || inum_loc == 2 || inum_loc == 4 || inum_loc == 8 ){
				unsigned int iface = 0;
				if(      inum_loc == 1 ){ iface = 0; }
				else if( inum_loc == 2 ){ iface = 1; }
				else if( inum_loc == 4 ){ iface = 2; }
				else if( inum_loc == 8 ){ iface = 3; }
				{
					const unsigned int itet_a = aTet[itet].s[iface];
					const unsigned * rel = tetRel[ aTet[itet].f[iface] ];
					const unsigned int ipo_a = aTet[itet_a].v[ rel[iface] ];
					assert( aTet[itet_a].s[ rel[iface] ] == itet );
					double vol0 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][0]], aTet[itet].v[noelTetFace[iface][1]],aPo3D);
					double vol1 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][1]], aTet[itet].v[noelTetFace[iface][2]],aPo3D);
					double vol2 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][2]], aTet[itet].v[noelTetFace[iface][0]],aPo3D);
					std::cout << vol0 << " " << vol1 << " " << vol2 << std::endl;
					if( vol0 < 0.0 || vol1 < 0.0 || vol2 < 0.0 ) continue;
				}
				std::cout << "InsertFace" << iface << std::endl;
				Msh::AddPointTet_Face(itet,iface,ipoin,aPo3D,aTet);
			}
			else if( inum_loc == 12 || inum_loc == 10 || inum_loc == 6 || inum_loc == 9 || inum_loc == 5 || inum_loc == 3 ){
				unsigned int idedge = 0;
				if(      inum_loc == 12 ){ idedge = 0; }
				else if( inum_loc == 10 ){ idedge = 1; }
				else if( inum_loc == 6  ){ idedge = 2; }
				else if( inum_loc == 9  ){ idedge = 4; }
				else if( inum_loc == 5  ){ idedge = 5; }
				else if( inum_loc == 3  ){ idedge = 8; }
				std::cout << "InsertEdge " << idedge << std::endl;
				ElemAroundEdge elared;
				Msh::MakeElemAroundEdge( elared, itet, idedge, aTet );
				Msh::AddPointTet_Edge( elared,ipoin,aPo3D,aTet);
			}
			else{ 
				std::vector<CPoint3D>::iterator itr = aPo3D.begin()+ipoin;
				aPo3D.erase(itr);
				std::cout << "Duplicated Point " << ipoin << " " << aPo3D[ipoin].p.x << " " << aPo3D[ipoin].p.y << " " << aPo3D[ipoin].p.z << std::endl;
//				getchar();
				break;
			}
			assert( aPo3D[ipoin].e != -1 );
//			assert( Msh::CheckTet(aTet,aPo3D) );
			Msh::DelaunayAroundPointTet(ipoin,aPo3D,aTet);
//			assert( Msh::CheckTet(aTet,aPo3D) );
			break;
		}
		if( aPo3D[ipoin].e == -1 ){
			std::cout << " Cannot Insert : " << ipoin << std::endl;
//			getchar();
		}
	}

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );


	{
		double min_len = 10.0;
		for(unsigned int itet=0;itet<aTet.size();itet++){
			for(unsigned int isedge=0;isedge<nSEdgeTet;isedge++){
				unsigned int ipo0 = aTet[itet].v[ sEdge2Noel[isedge][0] ];
				unsigned int ipo1 = aTet[itet].v[ sEdge2Noel[isedge][1] ];
				const double len = Distance( aPo3D[ipo0].p, aPo3D[ipo1].p );
				if( len < min_len ){
					min_len = len;
					std::cout << min_len << std::endl;
				}
			}
		}
	}


	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////



	for(int i=0;i<5000;i++){
		double x = -0.95 + rand()*1.9/(1.0+RAND_MAX);
		double y = -0.95 + rand()*1.9/(1.0+RAND_MAX);
		double z = -0.95 + rand()*1.9/(1.0+RAND_MAX);
		if( x*x+y*y > 0.95*0.95 ) continue;
		
		CVector3D vec_ins(x,y,z);
		if( oct_tree.IsPointInSphere(0.1,vec_ins) ) continue;

		assert( oct_tree.Check() );

		std::cout << x << " " << y << " " << z << std::endl;

		unsigned int ipo = aPo3D.size();
		
		aPo3D.resize( aPo3D.size()+1 );
		aPo3D[ipo].p = vec_ins;
		aPo3D[ipo].e = -1;
		aPo3D[ipo].poel = 0;

		const int iret = oct_tree.InsertPoint( ipo, vec_ins );
		assert( iret == -1 );
	}



	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////

	{
		double min_len = 10.0;
		for(unsigned int itet=0;itet<aTet.size();itet++){
			for(unsigned int isedge=0;isedge<nSEdgeTet;isedge++){
				unsigned int ipo0 = aTet[itet].v[ sEdge2Noel[isedge][0] ];
				unsigned int ipo1 = aTet[itet].v[ sEdge2Noel[isedge][1] ];
				const double len = Distance( aPo3D[ipo0].p, aPo3D[ipo1].p );
				if( len < min_len ){
					min_len = len;
					std::cout << min_len << std::endl;
				}
			}
		}
	}



	for(unsigned int ipoin=0;ipoin<aPo3D.size();ipoin++){
		if( aPo3D[ipoin].e >= 0 ) continue; 
		for(unsigned int itet=0;itet<aTet.size();itet++){
			
			unsigned int inum_loc = 0;
			{	
				const double vol = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], aTet[itet].v[2], aTet[itet].v[3], aPo3D);
				const double tol0 = 1.0e-20;
				const double ratio0 = Msh::TetVolume( ipoin, aTet[itet].v[1], aTet[itet].v[2], aTet[itet].v[3], aPo3D)/vol;
				if( ratio0 < -tol0 ) continue;
				const double ratio1 = Msh::TetVolume( aTet[itet].v[0], ipoin, aTet[itet].v[2], aTet[itet].v[3], aPo3D)/vol;
				if( ratio1 < -tol0 ) continue;
				const double ratio2 = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], ipoin, aTet[itet].v[3], aPo3D)/vol;
				if( ratio2 < -tol0 ) continue;
				const double ratio3 = Msh::TetVolume( aTet[itet].v[0], aTet[itet].v[1], aTet[itet].v[2], ipoin, aPo3D)/vol;
				if( ratio3 < -tol0 ) continue;

				const double tol1 = (1.0e-10/vol < 1.0e-5) ? 1.0e-10/vol : 1.0e-5;
				std::cout << tol1 << std::endl;
				if( ratio0 < tol1 ) inum_loc += 1;
				if( ratio1 < tol1 ) inum_loc += 2;
				if( ratio2 < tol1 ) inum_loc += 4;
				if( ratio3 < tol1 ) inum_loc += 8;
				std::cout << ipoin << " " << itet << " " << ratio0 << " " << ratio1 << " " << ratio2 << " " << ratio3 << std::endl;
			}

			if( inum_loc == 0 ){
				std::cout << "InsertElem " << std::endl;
				Msh::AddPointTet_Elem(itet,ipoin,aPo3D,aTet);
			}
			else if( inum_loc == 1 || inum_loc == 2 || inum_loc == 4 || inum_loc == 8 ){
				unsigned int iface = 0;
				if(      inum_loc == 1 ){ iface = 0; }
				else if( inum_loc == 2 ){ iface = 1; }
				else if( inum_loc == 4 ){ iface = 2; }
				else if( inum_loc == 8 ){ iface = 3; }
				{
					const unsigned int itet_a = aTet[itet].s[iface];
					const unsigned * rel = tetRel[ aTet[itet].f[iface] ];
					const unsigned int ipo_a = aTet[itet_a].v[ rel[iface] ];
					assert( aTet[itet_a].s[ rel[iface] ] == itet );
					double vol0 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][0]], aTet[itet].v[noelTetFace[iface][1]],aPo3D);
					double vol1 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][1]], aTet[itet].v[noelTetFace[iface][2]],aPo3D);
					double vol2 = TetVolume(ipoin,ipo_a,aTet[itet].v[noelTetFace[iface][2]], aTet[itet].v[noelTetFace[iface][0]],aPo3D);
					std::cout << vol0 << " " << vol1 << " " << vol2 << std::endl;
					if( vol0 < 0.0 || vol1 < 0.0 || vol2 < 0.0 ) continue;
				}
				std::cout << "InsertFace" << iface << std::endl;
				Msh::AddPointTet_Face(itet,iface,ipoin,aPo3D,aTet);
			}
			else if( inum_loc == 12 || inum_loc == 10 || inum_loc == 6 || inum_loc == 9 || inum_loc == 5 || inum_loc == 3 ){
				unsigned int idedge = 0;
				if(      inum_loc == 12 ){ idedge = 0; }
				else if( inum_loc == 10 ){ idedge = 1; }
				else if( inum_loc == 6  ){ idedge = 2; }
				else if( inum_loc == 9  ){ idedge = 4; }
				else if( inum_loc == 5  ){ idedge = 5; }
				else if( inum_loc == 3  ){ idedge = 8; }
				std::cout << "InsertEdge " << idedge << std::endl;
				ElemAroundEdge elared;
				Msh::MakeElemAroundEdge( elared, itet, idedge, aTet );
				Msh::AddPointTet_Edge( elared,ipoin,aPo3D,aTet);
			}
			else{ 
				std::vector<CPoint3D>::iterator itr = aPo3D.begin()+ipoin;
				aPo3D.erase(itr);
				std::cout << "Duplicated Point " << ipoin << " " << aPo3D[ipoin].p.x << " " << aPo3D[ipoin].p.y << " " << aPo3D[ipoin].p.z << std::endl;
//				getchar();
				break;
			}
			assert( aPo3D[ipoin].e != -1 );
//			assert( Msh::CheckTet(aTet,aPo3D) );
			Msh::DelaunayAroundPointTet(ipoin,aPo3D,aTet);
//			assert( Msh::CheckTet(aTet,aPo3D) );
			break;
		}
		if( aPo3D[ipoin].e == -1 ){
			std::cout << " Cannot Insert : " << ipoin << std::endl;
//			getchar();
		}
	}

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );

	Reconnect(aTet,aPo3D);
	assert( Msh::CheckTet(aTet,aPo3D) );

	{
		double min_len = 10.0;
		for(unsigned int itet=0;itet<aTet.size();itet++){
			for(unsigned int isedge=0;isedge<nSEdgeTet;isedge++){
				unsigned int ipo0 = aTet[itet].v[ sEdge2Noel[isedge][0] ];
				unsigned int ipo1 = aTet[itet].v[ sEdge2Noel[isedge][1] ];
				const double len = Distance( aPo3D[ipo0].p, aPo3D[ipo1].p );
				if( len < min_len ){
					min_len = len;
					std::cout << min_len << std::endl;
				}
			}
		}
	}


	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////


	std::vector<STet> aTet_in;
	{
		aTet_in.reserve( aTet.size() );
//		unsigned int icnt_tet = 0;
		for(unsigned int itet=0;itet<aTet.size();itet++){
			bool is_out_flg = false;
			for(unsigned int inotet=0;inotet<8;inotet++){
				unsigned int ipo0 = aTet[itet].v[inotet];
				if( ipo0 < 8 ){ 
					is_out_flg = true;
					break;
				}
			}
			if( is_out_flg ) continue;
			aTet_in.push_back( aTet[itet] );
		}
	}

	std::vector<int> po2vec;
	{	// aVecの最後に追加された節点を追加
		po2vec.resize(aPo3D.size(),-1);
		assert( vec2po.size() == aVec.size() );
		for(unsigned int ivec=0;ivec<aVec.size();ivec++){
			const unsigned int ipo0 = vec2po[ivec];
			po2vec[ipo0] = ivec;
		}

		unsigned int nvec_new = 0;
		for(unsigned int ipoin=0;ipoin<aPo3D.size();ipoin++){
			if( po2vec[ipoin] == -1 ){
				nvec_new++;
			}
		}
		const unsigned int nvec_old = aVec.size();
		aVec.resize( aVec.size()+nvec_new );
		unsigned int icnt = 0;
		for(unsigned int ipoin=0;ipoin<aPo3D.size();ipoin++){
			if( po2vec[ipoin] == -1 ){
				po2vec[ipoin] = nvec_old+icnt;
				aVec[nvec_old+icnt].x = aPo3D[ipoin].p.x;
				aVec[nvec_old+icnt].y = aPo3D[ipoin].p.y;
				aVec[nvec_old+icnt].z = aPo3D[ipoin].p.z;
				icnt++;
			}
		}
		assert( nvec_new == icnt );	
	}

	m_aTetAry.resize( 1 );
	m_aTetAry[0].m_aTet.resize( aTet_in.size() );
	for(unsigned int itet=0;itet<aTet_in.size();itet++){
	for(unsigned int inotet=0;inotet<4;inotet++){
		unsigned int ipo = aTet_in[itet].v[inotet];
		assert( ipo < aPo3D.size() );
		unsigned int ivec = po2vec[ ipo ];
		assert( ivec < aVec.size() );
		m_aTetAry[0].m_aTet[itet].v[inotet] = ivec;
	}
	}

	return true;
}
