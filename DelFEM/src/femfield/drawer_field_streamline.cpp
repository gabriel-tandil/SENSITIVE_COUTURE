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
// DrawerField.cpp : 場可視化クラス(DrawerField)の実装
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

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include "delfem/drawer_field_streamline.h"
#include "delfem/elem_ary.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer.h"
#include "delfem/vector3d.h"

using namespace Fem::Field::View;
using namespace Fem::Field;

/*
// 点を含む（辺上でもよい）三角形の番号を返す
// 見つからなければ-1を返す
int FindVelocityAtPoint(const double co[], double velo[], int ielem_stat,
						const CElSuP& elsup, 
						unsigned int id_field_velo, unsigned int id_ea_field, const Fem::Field::CFieldWorld& world)
{
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea_field,CORNER,true,world);
	////////////////
	const unsigned int nnoes = 3;
	int ielem_inc = -1;
	if( ielem_stat != -1 ){
		unsigned int noes[nnoes];
		es_v.GetNodes(ielem_stat,noes);
		for(unsigned int inoes=0;inoes<nnoes;inoes++){
			const unsigned int ino0 = noes[inoes];
			assert ( ino0 < elsup.nno );
			for(unsigned int ielsup=elsup.elsup_ind[ino0];ielsup<elsup.elsup_ind[ino0+1];ielsup++){
				const unsigned int ielem0 = elsup.elsup[ielsup];
				unsigned int noes0[nnoes];
				es_v.GetNodes(ielem0,noes0);
				////////////////
				double ev[nnoes][2];
				double ec[nnoes][2];
				for(unsigned int inoes=0;inoes<nnoes;inoes++){
					const unsigned int ino = noes0[inoes];
					ns_v.GetValue(ino,ev[inoes]);
					ns_c.GetValue(ino,ec[inoes]);
				}
				const double at = TriArea2D(ec[0],ec[1],ec[2]);
				const double a0 = TriArea2D(co,ec[1],ec[2]);
				const double a1 = TriArea2D(co,ec[2],ec[0]);
				const double a2 = TriArea2D(co,ec[0],ec[1]);
				if( a0 > -at*1.0e-3 && a1 > -at*1.0e-3 && a2 > -at*1.0e-3 ){ 
					velo[0] = (a0*ev[0][0] + a1*ev[1][0] + a2*ev[2][0])/at;
					velo[1] = (a0*ev[0][1] + a1*ev[1][1] + a2*ev[2][1])/at;
					return ielem0; 
				}
			}
		}
	}
	for(unsigned int ielem=0;ielem<es_v.GetSizeElem();ielem++)
	{	
		unsigned int noes[nnoes];
		es_v.GetNodes(ielem,noes);
		////////////////
		double ev[nnoes][2];
		double ec[nnoes][2];
		for(unsigned int inoes=0;inoes<nnoes;inoes++){
			const unsigned int ino = noes[inoes];
			ns_v.GetValue(ino,ev[inoes]);
			ns_c.GetValue(ino,ec[inoes]);
		}
		const double at = TriArea2D(ec[0],ec[1],ec[2]);
		const double a0 = TriArea2D(co,ec[1],ec[2]);
		const double a1 = TriArea2D(co,ec[2],ec[0]);
		const double a2 = TriArea2D(co,ec[0],ec[1]);
		if( a0 > -at*1.0e-3 && a1 > -at*1.0e-3 && a2 > -at*1.0e-3 ){ 
			velo[0] = (a0*ev[0][0] + a1*ev[1][0] + a2*ev[2][0])/at;
			velo[1] = (a0*ev[0][1] + a1*ev[1][1] + a2*ev[2][1])/at;
			ielem_inc = ielem; 
			return ielem_inc;
		}
	}
	return -1;
}
*/

// 三角形の面積を求める関数
static double TriArea2D(const double p0[], const double p1[], const double p2[]){
	return 0.5*( (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]) );
}


bool SeekPosition(bool iflg_forward, double dt, 
				  const double co_cur[2], const double ve_cur[2], 
				  unsigned int& id_ea_stat, unsigned int& ielem_stat, 
				  double co_nxt[2], double ve_nxt[2],
				  unsigned int id_field_velo, const Fem::Field::CFieldWorld& world )
{
	double fb = ( iflg_forward ) ? 1 : -1;
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
//	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
//	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	double r1,r2;
    const unsigned int imode = 1;
	if( imode == 0 ){
		co_nxt[0] = co_cur[0]+fb*ve_cur[0]*dt;
		co_nxt[1] = co_cur[1]+fb*ve_cur[1]*dt;
		return  fv.FindVelocityAtPoint(ve_nxt,  id_ea_stat,ielem_stat,r1,r2,  co_nxt,world);
	}
	else if( imode == 1 ){
		double co0[2] = { co_cur[0]+fb*ve_cur[0]*dt, co_cur[1]+fb*ve_cur[1]*dt };
		double ve0[2];
		if( !fv.FindVelocityAtPoint(ve0,  id_ea_stat,ielem_stat,r1,r2,  co0,world) ){ return false; }

		co_nxt[0] = co_cur[0]+fb*(ve_cur[0]+ve0[0])*0.5*dt;
		co_nxt[1] = co_cur[1]+fb*(ve_cur[1]+ve0[1])*0.5*dt;
		return fv.FindVelocityAtPoint(ve_nxt,  id_ea_stat,ielem_stat,r1,r2,  co_nxt, world);
	}
	else if( imode == 2 ){
		double co0[2] = { co_cur[0]+fb*ve_cur[0]*dt*0.5, co_cur[1]+fb*ve_cur[1]*dt*0.5 };
		double ve0[2];
		if( !fv.FindVelocityAtPoint(ve0,  id_ea_stat,ielem_stat,r1,r2,  co0,world) ){ return false; }

		double co1[2] = { co_cur[0]+fb*ve0[0]*dt*0.5, co_cur[1]+fb*ve0[1]*dt*0.5 };
		double ve1[2];
		if( !fv.FindVelocityAtPoint(ve1,  id_ea_stat,ielem_stat,r1,r2,  co1,world) ){ return false; }

		double co2[2] = { co_cur[0]+fb*ve1[0]*dt, co_cur[1]+fb*ve1[1]*dt };
		double ve2[2];
		if( !fv.FindVelocityAtPoint(ve2,  id_ea_stat,ielem_stat,r1,r2,  co2,world) ){ return false; }

		co_nxt[0] = co_cur[0]+fb*(ve_cur[0]+2*ve0[0]+2*ve1[0]+ve2[0])/6.0*dt;
		co_nxt[1] = co_cur[1]+fb*(ve_cur[1]+2*ve0[1]+2*ve1[1]+ve2[1])/6.0*dt;
		return fv.FindVelocityAtPoint(ve_nxt,  id_ea_stat,ielem_stat,r1,r2,  co_nxt,world);
	}
	return false;
}

void AddStreamLine(unsigned int id_ea_stat, unsigned int ielem_stat, double r1, double r2,
				   unsigned int id_field_velo, const Fem::Field::CFieldWorld& world,
				   std::vector< std::vector<unsigned int> >& aMaskElem,
				   std::vector<double>& aXYStreamline, 
				   const std::vector<double>& idea2height)
{
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);

	// 種の位置
	double co_stat[2];
	{	
		const Fem::Field::CElemAry::CElemSeg& es_c = fv.GetElemSeg(id_ea_stat,CORNER,false,world);
		unsigned int no[3];
		es_c.GetNodes(ielem_stat,no);
		double co[3][2];
		ns_c.GetValue(no[0],co[0]);
		ns_c.GetValue(no[1],co[1]);
		ns_c.GetValue(no[2],co[2]);
		co_stat[0] = (1-r1-r2)*co[0][0]+r1*co[1][0]+r2*co[2][0];
		co_stat[1] = (1-r1-r2)*co[0][1]+r1*co[1][1]+r2*co[2][1];
	}
	// 種の速度
	double  ve_stat[2];
	{
		const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea_stat,CORNER,true,world);
		unsigned int no[3];
		es_v.GetNodes(ielem_stat,no);
		double ve[3][2];
		ns_v.GetValue(no[0],ve[0]);
		ns_v.GetValue(no[1],ve[1]);
		ns_v.GetValue(no[2],ve[2]);
		ve_stat[0] = (1-r1-r2)*ve[0][0]+r1*ve[1][0]+r2*ve[2][0];
		ve_stat[1] = (1-r1-r2)*ve[0][1]+r1*ve[1][1]+r2*ve[2][1];
	}
	const double norm_ve_stat = sqrt(ve_stat[0]*ve_stat[0]+ve_stat[1]*ve_stat[1]);
	if( norm_ve_stat < 1.0e-20 ) return;

	// dtをいくらにとるか
	double dt;
	{
		const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea_stat,CORNER,true,world);
		unsigned int noes[3];
		es_v.GetNodes(ielem_stat,noes);
		double co[3][2];
		ns_c.GetValue(noes[0],co[0]);
		ns_c.GetValue(noes[1],co[1]);
		ns_c.GetValue(noes[2],co[2]);
		double area = TriArea2D(co[0],co[1],co[2]);
        dt = 0.5*sqrt(area)/norm_ve_stat;
	}

	// 流れに沿う流線を作る
	std::vector<double> st_f;
	{
		double co_cur[2] = { co_stat[0], co_stat[1] };
		double ve_cur[2] = { ve_stat[0], ve_stat[1] };
		unsigned int id_ea_cur = id_ea_stat;
        unsigned int ielem_cur = ielem_stat;
		int icnt = 0;
		for(;;){
			st_f.push_back( co_cur[0] );
			st_f.push_back( co_cur[1] );
			st_f.push_back( idea2height[id_ea_cur] );
            assert( id_ea_cur < aMaskElem.size() );
            if( ielem_cur >= aMaskElem[id_ea_cur].size() ){
                std::cout << ielem_cur << " " << aMaskElem[id_ea_cur].size() << std::endl;
            }
            assert( ielem_cur < aMaskElem[id_ea_cur].size() );
			aMaskElem[id_ea_cur][ielem_cur] = 1;
			double co_nxt[2], ve_nxt[2];
			unsigned int id_ea_nxt=id_ea_cur;
			unsigned int ielem_nxt=ielem_cur;
            if( !SeekPosition(true,dt, co_cur,ve_cur,  id_ea_nxt,ielem_nxt,  co_nxt,ve_nxt, 
					id_field_velo,world) ){ break; }
            assert( id_ea_nxt < aMaskElem.size() );
            assert( ielem_nxt < aMaskElem[id_ea_nxt].size() );
			if( aMaskElem[id_ea_nxt][ielem_nxt] != 0 && ielem_nxt != ielem_cur ) break;
            if( ielem_nxt == ielem_cur ){ icnt++; if( icnt > 5 ) break; }
			else{ icnt = 0; }
			const double norm_ve = sqrt(ve_nxt[0]*ve_nxt[0]+ve_nxt[1]*ve_nxt[1]);
            if( fabs(norm_ve-norm_ve_stat) > norm_ve_stat*0.5 ) break;
            id_ea_cur = id_ea_nxt;
			ielem_cur = ielem_nxt;            
			co_cur[0]=co_nxt[0];  co_cur[1]=co_nxt[1];
			ve_cur[0]=ve_nxt[0];  ve_cur[1]=ve_nxt[1];
		}
	}
	// 流れに逆らう方向の流線を作る
	std::vector<double> st_b;
	{
		double co_cur[2] = { co_stat[0], co_stat[1] };
		double ve_cur[2] = { ve_stat[0], ve_stat[1] };
		unsigned int id_ea_cur = id_ea_stat;
		unsigned int ielem_cur = ielem_stat;
		int icnt = 0;
		for(;;){
			double co_nxt[2], ve_nxt[2];
			unsigned int id_ea_nxt=id_ea_cur;
			unsigned int ielem_nxt=ielem_cur;
			if( !SeekPosition(false,dt, co_cur,ve_cur,   id_ea_nxt,ielem_nxt,   co_nxt,ve_nxt,  
					id_field_velo,world) ){ break; }
            assert( id_ea_nxt < aMaskElem.size() );
            assert( ielem_nxt < aMaskElem[id_ea_nxt].size() );
			if( aMaskElem[id_ea_nxt][ielem_nxt] != 0 && ielem_nxt != ielem_cur ) break;
            if( ielem_nxt == ielem_cur ){ icnt++; if( icnt > 5 ) break; }
			else{ icnt = 0; }
			const double norm_ve = sqrt(ve_nxt[0]*ve_nxt[0]+ve_nxt[1]*ve_nxt[1]);
            if( fabs(norm_ve-norm_ve_stat) > norm_ve_stat*0.5 ) break;
			id_ea_cur = id_ea_nxt;
			ielem_cur = ielem_nxt;
			co_cur[0]=co_nxt[0];  co_cur[1]=co_nxt[1];
			ve_cur[0]=ve_nxt[0];  ve_cur[1]=ve_nxt[1];
			st_b.push_back( co_cur[0] );
			st_b.push_back( co_cur[1] );
			st_b.push_back( idea2height[id_ea_cur] );
            assert( id_ea_cur < aMaskElem.size() );
            assert( ielem_cur < aMaskElem[id_ea_cur].size() );
			aMaskElem[id_ea_cur][ielem_cur] = 1;
		}
	}
	////////////////
	aXYStreamline.resize( st_b.size() + st_f.size() );
	unsigned int nno_b = st_b.size()/3;
	unsigned int nno_f = st_f.size()/3;
	for(unsigned int ino=0;ino<nno_b;ino++){
		aXYStreamline[ino*3  ] = st_b[(nno_b-1-ino)*3  ];
		aXYStreamline[ino*3+1] = st_b[(nno_b-1-ino)*3+1];
		aXYStreamline[ino*3+2] = st_b[(nno_b-1-ino)*3+2];
	}
	for(unsigned int ino=0;ino<nno_f;ino++){
		aXYStreamline[nno_b*3+ino*3  ] = st_f[ino*3  ];
		aXYStreamline[nno_b*3+ino*3+1] = st_f[ino*3+1];
		aXYStreamline[nno_b*3+ino*3+2] = st_f[ino*3+2];
	}	
}


void Fem::Field::View::MakeStreamLine(
					unsigned int id_field_velo, const Fem::Field::CFieldWorld& world,
					std::vector< std::vector<double> >& aStLine)
{
//	std::cout << "MakeStreamline" << std::endl;
	aStLine.clear();
	const Fem::Field::CField& fv = world.GetField(id_field_velo);
	assert( fv.GetNLenValue() == 2 );
	assert( fv.GetFieldType() == VECTOR2 );
	const Fem::Field::CNodeAry::CNodeSeg& ns_v = fv.GetNodeSeg(CORNER,true, world,VELOCITY);
	const Fem::Field::CNodeAry::CNodeSeg& ns_c = fv.GetNodeSeg(CORNER,false,world,VELOCITY);
	////////////////
	double sq_max_tot = 0;
	unsigned int ino_max = 0;
	{
		for(unsigned int ino=0;ino<ns_v.Size();ino++){
			double velo[2];
			ns_v.GetValue(ino,velo);
			if( velo[0]*velo[0]+velo[1]*velo[1] > sq_max_tot ){
				ino_max = ino;
				sq_max_tot = velo[0]*velo[0]+velo[1]*velo[1];
			}
		}
	}
	if( ino_max == 0 ) return;
	////////////////
	std::vector<double> idea2height;
	{
		int ilayer_min, ilayer_max;
		const std::vector<unsigned int>& aIdEA = fv.GetAryIdEA();
		if( aIdEA.size() > 0 ){
			ilayer_min = fv.GetLayer(aIdEA[0]);
			ilayer_max = ilayer_min;
		}
		else{ ilayer_min=0; ilayer_max=0; }
		for(unsigned int iiea=1;iiea<aIdEA.size();iiea++){
			int ilayer = fv.GetLayer(aIdEA[iiea]);
			ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
			ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
		}
		unsigned int max_id = 0;
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea = aIdEA[iiea];
			max_id = ( max_id > id_ea ) ? max_id : id_ea;
		}
		idea2height.resize(max_id+1);
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea = aIdEA[iiea];
			int ilayer = fv.GetLayer(id_ea);
			double height = (ilayer-ilayer_min+0.01)/(ilayer_max-ilayer_min+1);
			idea2height[id_ea] = height;
		}
	}
	////////////////////////////////////////////////
	// マスクを作る
	std::vector< std::vector<unsigned int> > aMaskElem;
	aMaskElem.clear();
	const std::vector<unsigned int> aIdEA = fv.GetAryIdEA();
	{
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea = aIdEA[iiea];
			aMaskElem.resize( id_ea+1 );
		}
		for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
			unsigned int id_ea = aIdEA[iiea];
			const Fem::Field::CElemAry& ea = world.GetEA(id_ea);
			aMaskElem[id_ea].resize(ea.Size(),0);
		}
	}

	////////////////////////////////
	// ３頂点ともratio以内の速さだったらマスクする
	for(unsigned int iiea=0;iiea<aIdEA.size();iiea++){
		unsigned int id_ea = aIdEA[iiea];
		const Fem::Field::CElemAry::CElemSeg& es_v = fv.GetElemSeg(id_ea,CORNER,true,world);
    assert( aMaskElem[id_ea].size() == es_v.Size() );
		for(unsigned int ielem=0;ielem<es_v.Size();ielem++){
			unsigned int noes[3];
			es_v.GetNodes(ielem,noes);
			const double ratio = 0.003;
			double ve[2];
			ns_v.GetValue(noes[0],ve); if( ve[0]*ve[0]+ve[1]*ve[1] > sq_max_tot*ratio ) continue;
			ns_v.GetValue(noes[1],ve); if( ve[0]*ve[0]+ve[1]*ve[1] > sq_max_tot*ratio ) continue;
			ns_v.GetValue(noes[2],ve); if( ve[0]*ve[0]+ve[1]*ve[1] > sq_max_tot*ratio ) continue;
			aMaskElem[id_ea][ielem] = 1;
		}
	}
	////////////////////////////////
	{
		double co_max[2];
		double ve_max[2];
		ns_c.GetValue(ino_max,co_max);
		unsigned int id_ea_max, ielem_max;
		double r1,r2;
        if( fv.FindVelocityAtPoint(ve_max,  id_ea_max,ielem_max,r1,r2,  co_max,world) ){
            assert( id_ea_max < aMaskElem.size() );
            assert( ielem_max < aMaskElem[id_ea_max].size() );
			std::vector<double> aXYHStreamline;
            AddStreamLine(id_ea_max,ielem_max,r1,r2,  id_field_velo,world,  aMaskElem,aXYHStreamline, idea2height);
			if( aXYHStreamline.size() >= 10 ){
				aStLine.push_back( aXYHStreamline );
			}
        }
	}
	////////////////////////////////
	std::vector<unsigned int> aMaskNode;
	aMaskNode.resize( ns_v.Size(), 0 );
	for(;;){
		double sq_max = 0;
		int ino_max = -1;
		for(unsigned int ino=0;ino<ns_v.Size();ino++){
			if( aMaskNode[ino] != 0 ) continue;
			double ve[2], co[2];
			ns_v.GetValue(ino,ve);
			ns_c.GetValue(ino,co);
			if( ve[0]*ve[0]+ve[1]*ve[1] < sq_max ){ continue; }
			if( ve[0]*ve[0]+ve[1]*ve[1] < sq_max_tot*0.005 ){
				aMaskNode[ino] = 1;
				continue;
			}
			bool iflg = false;
			for(unsigned int ist=0;ist<aStLine.size();ist++){
				const unsigned int nno = aStLine[ist].size()/2;
				for(unsigned int ino=0;ino<nno;ino++){
					double x0 = aStLine[ist][ino*2  ];
					double y0 = aStLine[ist][ino*2+1];
					const double sq_dist = (x0-co[0])*(x0-co[0]) + (y0-co[1])*(y0-co[1]);
                    if( sq_dist < 0.01 ){ iflg = true; break; }
				}
				if( iflg ) break;
			}
			if( iflg ){ aMaskNode[ino] = 1; continue; }
			ino_max = ino;
			sq_max = ve[0]*ve[0]+ve[1]*ve[1];
		}
//		std::cout << "Max : " << sq_max << " " << ino_max << std::endl;
		if( ino_max == -1 ) break;
		aMaskNode[ino_max] = 1;
		{
			double co_max[2];
			double ve_max[2];
			double r1,r2;
			unsigned int id_ea_max, ielem_max;
			ns_c.GetValue(ino_max,co_max);
			fv.FindVelocityAtPoint(ve_max,  id_ea_max,ielem_max,r1,r2,  co_max,world);
			assert( id_ea_max < aMaskElem.size() );
			assert( ielem_max < aMaskElem[id_ea_max].size() );
			std::vector<double> aXYHStreamline;
			AddStreamLine(id_ea_max,ielem_max,r1,r2,  id_field_velo,world, aMaskElem, aXYHStreamline, idea2height);
            if( aXYHStreamline.size() >= 50 ){
				aStLine.push_back( aXYHStreamline );
			}
		}
	}
//	std::cout << "MakeStreamline end" << std::endl;
}


Com::CBoundingBox3D Fem::Field::View::CDrawerStreamline::GetBoundingBox( double rot[] ) const
{
	if( aStLine.empty() ){ return Com::CBoundingBox3D(); }
	if( aStLine[0].empty() ){ return Com::CBoundingBox3D(); }
	double x_min = aStLine[0][0];
	double y_min = aStLine[0][1];
	double x_max = x_min;
	double y_max = y_min;
	for(unsigned int ist=0;ist<aStLine.size();ist++){
		const unsigned int nno = aStLine[ist].size()/2;
		for(unsigned int ino=0;ino<nno;ino++){
			double x1 = aStLine[ist][ino*2  ];
			double y1 = aStLine[ist][ino*2+1];
			x_min = ( x1 < x_min ) ? x1 : x_min;
			x_max = ( x1 > x_max ) ? x1 : x_max;
			y_min = ( y1 < y_min ) ? y1 : y_min;
			y_max = ( y1 > y_max ) ? y1 : y_max;
		}
	}
//		std::cout << x_min << " " << x_max << "  " << y_min << " " << y_max << std::endl;
	return Com::CBoundingBox3D(x_min,x_max, y_min,y_max,-1,1);
}

void Fem::Field::View::CDrawerStreamline::Draw() const 
{
	::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable(GL_BLEND);
    ::glColor3d(0,0,0);
	////////////////
    for(unsigned int ist=0;ist<aStLine.size();ist++){
        const unsigned int nno = aStLine[ist].size()/3;
        const double x0 = aStLine[ist][(nno-1)*3  ];
        const double y0 = aStLine[ist][(nno-1)*3+1];
        const double h0 = aStLine[ist][(nno-1)*3+2];
        double xe2=0, ye2=0;
        for(unsigned int i=2;i<nno;i++){
            const double x1 = aStLine[ist][(nno-i)*3  ];
            const double y1 = aStLine[ist][(nno-i)*3+1];
            const double len0 = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1) );
            xe2 = (x1-x0)/len0*0.025;
            ye2 = (y1-y0)/len0*0.025;
            if( len0 > 0.025 ) break;
        }
        ::glBegin(GL_TRIANGLES);
        ::glVertex3d( x0-xe2*0.4, y0-ye2*0.4, h0 );
        ::glVertex3d( x0+xe2+ye2*0.4, y0+ye2-xe2*0.4, h0 );
        ::glVertex3d( x0+xe2-ye2*0.4, y0+ye2+xe2*0.4, h0 );
        ::glEnd();
    }
    ::glLineWidth(5);
    if( this->m_is_anti_aliasing ){ // アンチエリアシングの導入
        ::glEnable(GL_LINE_SMOOTH);
        ::glEnable(GL_BLEND);
        ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        ::glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    }
	for(unsigned int ist=0;ist<aStLine.size();ist++){
		const unsigned int nno = aStLine[ist].size()/3;
		::glBegin(GL_LINE_STRIP);
		for(unsigned int ino=0;ino<nno;ino++){
			::glVertex3d( aStLine[ist][ino*3  ], aStLine[ist][ino*3+1], aStLine[ist][ino*3+2] );
		}
		::glEnd();
    }
}
