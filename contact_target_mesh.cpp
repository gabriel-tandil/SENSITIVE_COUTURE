/*
 *  contact_target_mesh.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 7/26/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */

#include <fstream>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <vector>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/vector3d.h"
#include "contact_target.h"

CContactTarget3D_Mesh::CContactTarget3D_Mesh(){
	nnode_ = 0;	pXYZs_ = 0;
	ntri_ = 0;	aTri_ = 0;
	pBoxel_ = 0;
	is_hole = false;
}

CContactTarget3D_Mesh::~CContactTarget3D_Mesh(){
	if( pXYZs_  != 0 ){ delete pXYZs_; }
	if( aTri_   != 0 ){ delete aTri_; }
	if( pBoxel_ != 0 ){ delete pBoxel_; }
}

void CContactTarget3D_Mesh::GetCenterWidth(double& cx, double& cy, double& cz,  
                                           double& wx, double& wy, double& wz)
{  
	double x_min = pXYZs_[0], x_max = pXYZs_[0];
  double y_min = pXYZs_[1], y_max = pXYZs_[1];
  double z_min = pXYZs_[2], z_max = pXYZs_[2];    
	for(unsigned int ino=1;ino<nnode_;ino++){
		x_min = ( x_min < pXYZs_[ino*3+0] ) ? x_min : pXYZs_[ino*3+0];
		x_max = ( x_max > pXYZs_[ino*3+0] ) ? x_max : pXYZs_[ino*3+0];			
		y_min = ( y_min < pXYZs_[ino*3+1] ) ? y_min : pXYZs_[ino*3+1];
		y_max = ( y_max > pXYZs_[ino*3+1] ) ? y_max : pXYZs_[ino*3+1];			
		z_min = ( z_min < pXYZs_[ino*3+2] ) ? z_min : pXYZs_[ino*3+2];
		z_max = ( z_max > pXYZs_[ino*3+2] ) ? z_max : pXYZs_[ino*3+2];
	}
  cx = (x_min+x_max)*0.5;
  cy = (y_min+y_max)*0.5;
  cz = (z_min+z_max)*0.5;  
  wx = x_max-x_min;
  wy = y_max-y_min;
  wz = z_max-z_min;    
}

/*
void CContactTarget3D_Mesh::Load_Off(const std::string& fname)
{
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
		std::cout << "Fail Read Fail" << std::endl;
		return;
	}		
	std::string str;
	fin >> str;
	fin >> nnode_ >> ntri_ >> str;
	std::cout << "Load Off Nnode: ntri :" << nnode_ << " " << ntri_ << std::endl;
	pXYZs_ = new double [nnode_*3];
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x,y,z;
		fin >> x >> y >> z;
//		std::cout << ino << " " << x << " " << y << " " << z << std::endl;
		pXYZs_[ino*3+0] = x;
		pXYZs_[ino*3+1] = y;
		pXYZs_[ino*3+2] = z;
	}		
	for(unsigned int itri=0;itri<ntri_;itri++){			
		int itmp, i1, i2, i3;
		fin >> itmp >> i1 >> i2 >> i3;
		aTri_[itri*3+0] = i1;
		aTri_[itri*3+1] = i2;
		aTri_[itri*3+2] = i3;
//		std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
}

void CContactTarget3D_Mesh::Load_Gmv(const std::string& fname)
{
	std::cout << "File load " << fname << std::endl;
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
		std::cout << "Fail Read Fail" << std::endl;
		return;
	}		
	std::string str;
	fin >> str;	
	fin >> str; 
	fin >> str; 
	fin >> str;	
	fin >> str; 
	fin >> str; 
	fin >> str;	
	fin >> str; 
	fin >> str; 
	fin >> nnode_;
	std::cout << "Nnode " << nnode_ << std::endl;
	pXYZs_ = new double [nnode_*3];
  
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x;
		fin >> x;
		pXYZs_[ino*3+0] = x;
	}		
	for(unsigned int ino=0;ino<nnode_;ino++){
		double y;
		fin >> y;
		pXYZs_[ino*3+1] = y;
	}		
	for(unsigned int ino=0;ino<nnode_;ino++){
		double z;
		fin >> z;
		pXYZs_[ino*3+2] = z;
	}		
	
	fin >> str;
	fin >> ntri_;
	std::cout << "Ntri " << ntri_ << std::endl;
	aTri_ = new unsigned int [ntri_*3];
	for(unsigned int itri=0;itri<ntri_;itri++){			
		int itmp, i1, i2, i3;
		fin >> str >> itmp >> i1 >> i2 >> i3;
		aTri_[itri*3+0] = i1-1;
		aTri_[itri*3+1] = i2-1;
		aTri_[itri*3+2] = i3-1;
		//			std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
}	

void CContactTarget3D_Mesh::Load_Ply(const std::string& fname)
{  
	std::cout << "File load " << fname << std::endl;
	std::ifstream fin;
	fin.open(fname.c_str());
	if( fin.fail() ){
		std::cout << "Fail Read Fail" << std::endl;
		return;
	}		
  const unsigned int nbuff = 256;
  char buff[nbuff], buff1[nbuff], buff2[nbuff];
	std::string str1,str2;
  fin.getline(buff,nbuff);  // ply
  fin.getline(buff,nbuff);  // format asi 1.0
  for(;;){
    fin.getline(buff,nbuff);
    if( strncmp(buff, "comment ", 8) != 0 ){ break; }
  }
  /////
  sscanf(buff,"%s %s %d",buff1,buff2,&nnode_);
	std::cout << "Nnode " << nnode_ << std::endl;  
  ////
  for(;;){
    fin.getline(buff,nbuff);
    if( strncmp(buff, "property ", 9) != 0 ){ break; }
  }
  sscanf(buff,"%s %s %d",buff1,buff2,&ntri_);
  std::cout << "NTri " << ntri_ << std::endl;  
  /////  
  fin.getline(buff,nbuff);  // property list int int vertex_indices
  fin.getline(buff,nbuff);  // end header
  ////
	pXYZs_ = new double [nnode_*3];  
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x,y,z;
		fin >> x >> y >> z;
    //		std::cout << ino << " " << x << " " << y << " " << z << std::endl;
		pXYZs_[ino*3+0] = x;
		pXYZs_[ino*3+1] = y;
		pXYZs_[ino*3+2] = z;
	}		
	aTri_ = new unsigned int [ntri_*3];  
	for(unsigned int itri=0;itri<ntri_;itri++){			
		int itmp, i1, i2, i3;
		fin >> itmp >> i1 >> i2 >> i3;
		aTri_[itri*3+0] = i1;
		aTri_[itri*3+1] = i2;
		aTri_[itri*3+2] = i3;
    //		std::cout << itri << " " << itmp << " " << i1 << " " << i2 << " " << i3 << std::endl;
	}
}
*/

void CContactTarget3D_Mesh::SetMesh(const std::vector<unsigned int>& aTri, const std::vector<double>& aXYZ)
{  
  ntri_ = aTri.size()/3;
	aTri_ = new unsigned int [ntri_*3];  
  for(unsigned int i=0;i<ntri_*3;i++){ this->aTri_[i] = aTri[i]; }
  nnode_ = aXYZ.size()/3;
	pXYZs_ = new double [nnode_*3];  
  for(unsigned int i=0;i<nnode_*3;i++){ this->pXYZs_[i] = aXYZ[i]; }
}

void CContactTarget3D_Mesh::Draw() const
{

	::glColor3d(1,0,0);
	::glBegin(GL_LINES);
	for(unsigned int itri=0;itri<ntri_;itri++){
		unsigned int i1 = aTri_[itri*3+0];
		unsigned int i2 = aTri_[itri*3+1];
		unsigned int i3 = aTri_[itri*3+2];
		::glVertex3dv(pXYZs_+i1*3);
		::glVertex3dv(pXYZs_+i2*3);
		
		::glVertex3dv(pXYZs_+i2*3);
		::glVertex3dv(pXYZs_+i3*3);
		
		::glVertex3dv(pXYZs_+i3*3);
		::glVertex3dv(pXYZs_+i1*3);			
	}
	::glEnd();
}

// return penetration depth (inside is positive)
double CContactTarget3D_Mesh::Projection
(double px, double py, double pz,
 double n[3]) const // normal outward
{
	unsigned int inout;
	double dist;	
	if( pBoxel_ != 0 ){
		inout = this->FindInOut_Boxel(px,py,pz);
		dist = this->Distance_Mesh_Boxel(px,py,pz, n);
	}
	else {
		inout = this->FindInOut(px,py,pz);
		dist = this->Distance_Mesh(px,py,pz, n);
	}
	if(      inout == 0 ){ return  dist; }
	else if( inout == 1 ){ return -dist; }
	return -dist;	// if not sure assume out
}

void CContactTarget3D_Mesh::GetMesh
(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ, double elen) const
{
  aTri.resize(ntri_*3);
  for(unsigned int i=0;i<ntri_*3;i++){ aTri[i] = this->aTri_[i]; }
  aXYZ.resize(nnode_*3);
  for(unsigned int i=0;i<nnode_*3;i++){ aXYZ[i] = this->pXYZs_[i]; }
}


// 0:in 1:out 2:not sure
double CContactTarget3D_Mesh::Distance_Mesh
(double px, double py, double pz,
 double n[3]) const
{
	double p0[3] = {px,py,pz};
	n[0] = pXYZs_[0]-px;
	n[1] = pXYZs_[1]-py;
	n[2] = pXYZs_[2]-pz;
	double dist = Com::Length3D(n);
	for(unsigned int ino=1;ino<nnode_;ino++){		
		const double d0 = Com::Distance3D(p0, pXYZs_+ino*3);
		if( d0 >= dist ){ continue; }
		dist = d0;
		n[0] = pXYZs_[ino*3+0]-px;
		n[1] = pXYZs_[ino*3+1]-py;
		n[2] = pXYZs_[ino*3+2]-pz;
	}
	for(unsigned int itri=0;itri<ntri_;itri++){
		const unsigned int i1 = aTri_[itri*3+0];
		const unsigned int i2 = aTri_[itri*3+1];
		const unsigned int i3 = aTri_[itri*3+2];
		const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		const double sign = ( v0 > 0 ) ? 1 : -1;
		double normal[3];
    Com::NormalTri3D(normal,pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		{
			const double invlen = 1.0/Com::Length3D(normal);
			normal[0] *= invlen;
			normal[1] *= invlen;
			normal[2] *= invlen;					
		}
		const double p1[3] = { px+normal[0], py+normal[1], pz+normal[2] };
		const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
		const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
		const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
		if( fabs(v1+v2+v3) < 1.0e-10 ){	continue; } // p0 and p1 is on the triangle
		double inv_v4 = 1.0/(v1+v2+v3);
		const double r1 = v1*inv_v4;
		const double r2 = v2*inv_v4;
		const double r3 = v3*inv_v4;
		const double tol = 1.0e-4;		
		if( r1 > -tol && r2 > -tol && r3 > -tol )
		{
			const double dir2[3] = {
				pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
				pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
				pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
			//			const double dotdir = Dot3D(normal,dir2);
			const double d1 = Com::Length3D(dir2);
			if( d1 >= dist ) continue;
			dist = d1;
			n[0] = dir2[0];
			n[1] = dir2[1];
			n[2] = dir2[2];
		}
	}
	{
		const double invlen = 1.0/Com::Length3D(n);
		n[0] *= invlen;
		n[1] *= invlen;
		n[2] *= invlen;		
	}
	return dist;
}


// 0:in 1:out 2:not sure
double CContactTarget3D_Mesh::Distance_Mesh_Boxel
(double px, double py, double pz,
 double n[3]) const
{
	assert( this->pBoxel_ != 0 );
	double p0[3] = {px,py,pz};
	double dist = pBoxel_->GetWidth()*3;
	
	pBoxel_->Find_NearestTriCand(p0,aIndTriCand);
	//	std::cout << aIndTriCand.size() << std::endl;
	
	for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
		unsigned int itri = aIndTriCand[iitri];
		assert( itri < ntri_ );
		const unsigned int i1 = aTri_[itri*3+0];
		const unsigned int i2 = aTri_[itri*3+1];
		const unsigned int i3 = aTri_[itri*3+2];
		{	
			const double d1 = Com::Distance3D(p0, pXYZs_+i1*3);
			if( d1 < dist ){
				dist = d1;
				n[0] = pXYZs_[i1*3+0]-px;
				n[1] = pXYZs_[i1*3+1]-py;
				n[2] = pXYZs_[i1*3+2]-pz;			
			}			
			const double d2 = Com::Distance3D(p0, pXYZs_+i2*3);
			if( d2 < dist ){ 
				dist = d2;
				n[0] = pXYZs_[i2*3+0]-px;
				n[1] = pXYZs_[i2*3+1]-py;
				n[2] = pXYZs_[i2*3+2]-pz;			
			}
			const double d3 = Com::Distance3D(p0, pXYZs_+i3*3);
			if( d3 < dist ){ 
				dist = d3;
				n[0] = pXYZs_[i3*3+0]-px;
				n[1] = pXYZs_[i3*3+1]-py;
				n[2] = pXYZs_[i3*3+2]-pz;			
			}			
		}
		const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		const double sign = ( v0 > 0 ) ? 1 : -1;
		double normal[3];
    Com::NormalTri3D(normal,pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		{
			const double invlen = 1.0/Com::Length3D(normal);
			normal[0] *= invlen;
			normal[1] *= invlen;
			normal[2] *= invlen;					
		}
		const double p1[3] = { px+normal[0], py+normal[1], pz+normal[2] };
		const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
		const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
		const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
		if( fabs(v1+v2+v3) < 1.0e-10 ){	continue; } // p0 and p1 is on the triangle
		double inv_v4 = 1.0/(v1+v2+v3);
		const double r1 = v1*inv_v4;
		const double r2 = v2*inv_v4;
		const double r3 = v3*inv_v4;
		const double tol = 1.0e-4;		
		if( r1 > -tol && r2 > -tol && r3 > -tol )
		{
			const double dir2[3] = {
				pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
				pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
				pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
			//			const double dotdir = Dot3D(normal,dir2);
			const double d1 = Com::Length3D(dir2);
			if( d1 >= dist ) continue;
			dist = d1;
			n[0] = dir2[0];
			n[1] = dir2[1];
			n[2] = dir2[2];
		}
	}
	{
		const double invlen = 1.0/Com::Length3D(n);
		n[0] *= invlen;
		n[1] *= invlen;
		n[2] *= invlen;		
	}
	return dist;
}


// 0:in 1:out 2:not sure
unsigned int CContactTarget3D_Mesh::FindInOut_IntersectionRay
(double px, double py, double pz,
 const double dir[3]) const
{
	double p0[3] = {px,py,pz};
	double p1[3] = {px+dir[0],py+dir[1],pz+dir[2]};
	unsigned int icnt = 0;
	for(unsigned int itri=0;itri<ntri_;itri++){
		unsigned int i1 = aTri_[itri*3+0];
		unsigned int i2 = aTri_[itri*3+1];
		unsigned int i3 = aTri_[itri*3+2];
		const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		const double sign = ( v0 > 0 ) ? 1 : -1;
		const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
		const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
		const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
		if( fabs(v1+v2+v3) < 1.0e-10 ) return 2;	// p0 and p1 is on the triangle
		double inv_v4 = 1.0/fabs(v1+v2+v3);
		const double r1 = v1*inv_v4;
		const double r2 = v2*inv_v4;
		const double r3 = v3*inv_v4;
		const double tol = 1.0e-2;
		if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
		if( r1 < tol || r2 < tol || r3 < tol ) return 2;	// on the edge
		double dir2[3] = {
			pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
			pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
			pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
		double dotdir = Com::Dot3D(dir,dir2);
		if( dotdir > 0 ) icnt++;
	}
	if( icnt % 2 == 0 ) return 1;
	return 0;
}


// 0:in 1:out 2:not sure
unsigned int CContactTarget3D_Mesh::FindInOut_IntersectionRay_Boxel
(double px, double py, double pz,
 const double dir[3]) const
{
	assert( pBoxel_ != 0 );
	double p0[3] = {px,py,pz};
	pBoxel_->Find_IntersecTriCand(p0,dir,aIndTriCand);	
	//	std::cout << aIndTriCand_.size() << std::endl;
	double p1[3] = {px+dir[0],py+dir[1],pz+dir[2]};	  
	unsigned int icnt = 0;
	for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
		unsigned int itri = aIndTriCand[iitri];
		if( aFlgTriUsed[itri] == 1 ) continue;
		aFlgTriUsed[itri] = 1;
		unsigned int i1 = aTri_[itri*3+0];
		unsigned int i2 = aTri_[itri*3+1];
		unsigned int i3 = aTri_[itri*3+2];
		const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
		const double sign = ( v0 > 0 ) ? 1 : -1;
		const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
		const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
		const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
		if( fabs(v1+v2+v3) < 1.0e-10 ) goto AMBIGUOUS;	// p0 and p1 is on the triangle
		double inv_v4 = 1.0/fabs(v1+v2+v3);
		const double r1 = v1*inv_v4;
		const double r2 = v2*inv_v4;
		const double r3 = v3*inv_v4;
		const double tol = 1.0e-2;
		if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
		if( r1 < tol || r2 < tol || r3 < tol ) goto AMBIGUOUS;	// on the edge
		double dir2[3] = {
			pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3 - px,
			pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3 - py,
			pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 - pz};
		double dotdir = Com::Dot3D(dir,dir2);
		if( dotdir > 0 ) icnt++;
	}
	for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
		unsigned int itri = aIndTriCand[iitri];
		aFlgTriUsed[itri] = 0;
	}
	//	std::cout << "Cunt" << icnt << std::endl;
	if( icnt % 2 == 0 ) return 1;
	return 0;
AMBIGUOUS:
	for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
		unsigned int itri = aIndTriCand[iitri];
		aFlgTriUsed[itri] = 0;
	}	
	return 2;
}




void CContactTarget3D_Mesh::BuildBoxel()
{
	//	double c[3] = {0.0542,-0.04374532,0.06234};
	double c[3] = {0,0,0};	
  double w[3] = {0,0,0};
  this->GetCenterWidth(c[0],c[1],c[2], w[0],w[1],w[2]);
  double width = w[0];
  width = ( w[1] > width ) ? w[1] : width;
  width = ( w[2] > width ) ? w[2] : width;  
	if( pBoxel_ != 0 ){ delete pBoxel_; }
	pBoxel_ = new CSpatialHash_Grid3D(32,c,width*0.5*1.13454325);
	for(unsigned int itri=0;itri<ntri_;itri++){
		unsigned int i0 = aTri_[itri*3+0];
		unsigned int i1 = aTri_[itri*3+1];
		unsigned int i2 = aTri_[itri*3+2];		
		pBoxel_->AddTri(itri, pXYZs_+i0*3, pXYZs_+i1*3, pXYZs_+i2*3);
	}
	if( !is_hole ){ pBoxel_->BuildOutFlg(); }
	aFlgTriUsed.clear();
	aFlgTriUsed.resize(ntri_,0);
	aIndTriCand.reserve(2048);
}




unsigned int CContactTarget3D_Mesh::FindInOut(double px, double py, double pz) const
{
	unsigned int icnt_in  = 0;
	unsigned int icnt_out = 0;	
	for(unsigned int i=0;i<10;i++){
		const double theta = i*6.28/10+0.1431432154;
		//		const double dir[3] = { sin(theta)*cos(theta*2), sin(theta)*sin(theta*2), cos(theta) };
		const double dir[3] = { 1, 0, 0 };		
		unsigned int ires = FindInOut_IntersectionRay(px,py,pz, dir);
		if( ires != 2 && !is_hole ){ return ires; }
		if( ires == 0 ) icnt_in++;
		if( ires == 1 ) icnt_out++;
	}	
	if( icnt_in > 5 )  return 0;
	if( icnt_out > 5 ) return 1;
	/*	for(unsigned int i=1;i<20;i++){
	 const double theta = i*6.28*543260;
	 const double dir[3] = { sin(theta)*cos(theta*3), cos(theta), sin(theta)*sin(theta*3) };
	 unsigned int ires = FindInOut_IntersectionRay(px,py,pz, dir);
	 if( ires == 0 ) icnt_in++;
	 if( ires == 1 ) icnt_out++;
	 if( icnt_in  - icnt_out > 5 ) return 0;
	 if( icnt_out - icnt_in  > 5 ) return 1;
	 }*/
	return 2;
}	

unsigned int CContactTarget3D_Mesh::FindInOut_Boxel
(double px, double py, double pz) const
{
	assert( pBoxel_ != 0 );
	double p[3] = { px, py, pz };
	if( pBoxel_->IsOut(p) && !is_hole ){ return 1; }
	unsigned int icnt_in  = 0;
	unsigned int icnt_out = 0;	
	for(unsigned int i=0;i<10;i++){
		const double theta = i*6.28/10+0.15432452356673;
		const double dir[3] = { sin(theta)*cos(theta*2), sin(theta)*sin(theta*2), cos(theta) };
		//		const double dir[3] = { -0.35, 0.1342, 0.3 };
		unsigned int ires  = FindInOut_IntersectionRay_Boxel(px,py,pz, dir);
		//		unsigned int ires1 = FindInOut_IntersectionRay(px,py,pz, dir);
		/*		if( ires != ires1 ){
		 std::cout << "hoge " << px << " " << py << " " << pz << std::endl;			
		 }*/
		if( ires != 2 && !is_hole ){ return ires; }
		if( ires == 0 ) icnt_in++;
		if( ires == 1 ) icnt_out++;
	}	
	if( icnt_in > 5 )  return 0;
	if( icnt_out > 5 ) return 1;
	return 2;
}

void CContactTarget3D_Mesh::Translate(double x, double y, double z)
{
	for(unsigned int ino=0;ino<nnode_;ino++){
		pXYZs_[ino*3+0] += x;
		pXYZs_[ino*3+1] += y;
		pXYZs_[ino*3+2] += z;		
	}
}

bool CContactTarget3D_Mesh::IntersectionPoint
(double p[3], 
 const double org[3], const double dir[3]) const
{
  if( pBoxel_ != 0 ){
    std::vector<unsigned int> aIndTriCand;
    pBoxel_->Find_IntersecTriCand(org,dir, aIndTriCand);
    if( aIndTriCand.empty() ) return false;    
    bool iflg = false;
    double min_dist;
    const double p0[3] = {org[0],org[1],org[2]};
    const double p1[3] = {org[0]+dir[0],org[1]+dir[1],org[2]+dir[2]};	  
    for(unsigned int iitri=0;iitri<aIndTriCand.size();iitri++){
      unsigned int itri = aIndTriCand[iitri];
      unsigned int i1 = aTri_[itri*3+0];
      unsigned int i2 = aTri_[itri*3+1];
      unsigned int i3 = aTri_[itri*3+2];
      const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
      const double sign = ( v0 > 0 ) ? 1 : -1;
      const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
      const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
      const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
      double inv_v4 = 1.0/fabs(v1+v2+v3);
      const double r1 = v1*inv_v4;
      const double r2 = v2*inv_v4;
      const double r3 = v3*inv_v4;
      const double tol = 1.0e-2;
      if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
      double end[3] = {
        pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3,
        pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3,
        pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 };
      double dir2[3] = { end[0]-org[0], end[1]-org[1], end[2]-org[2] };
      double dotdir = Com::Dot3D(dir,dir2);
      if( dotdir < 0 ) continue;
      if( !iflg ){ 
        min_dist = dotdir; 
        p[0]=end[0];  p[1]=end[1]; p[2]=end[2];        
      }
      else if( dotdir < min_dist ){
        min_dist = dotdir;
        p[0]=end[0];  p[1]=end[1]; p[2]=end[2];
      }
      iflg = true;      
    }
    if( iflg ) return true;
  } 
  //////////////////////////////////////
  else {    
    bool iflg = false;
    double min_dist;
    const double p0[3] = {org[0],org[1],org[2]};
    const double p1[3] = {org[0]+dir[0],org[1]+dir[1],org[2]+dir[2]};	  
    for(unsigned int itri=0;itri<ntri_;itri++){
      unsigned int i1 = aTri_[itri*3+0];
      unsigned int i2 = aTri_[itri*3+1];
      unsigned int i3 = aTri_[itri*3+2];
      const double v0 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);
      const double sign = ( v0 > 0 ) ? 1 : -1;
      const double v1 = Com::TetVolume3D(p0, pXYZs_+i2*3, pXYZs_+i3*3, p1)*sign;
      const double v2 = Com::TetVolume3D(p0, pXYZs_+i3*3, pXYZs_+i1*3, p1)*sign;
      const double v3 = Com::TetVolume3D(p0, pXYZs_+i1*3, pXYZs_+i2*3, p1)*sign;
      double inv_v4 = 1.0/fabs(v1+v2+v3);
      const double r1 = v1*inv_v4;
      const double r2 = v2*inv_v4;
      const double r3 = v3*inv_v4;
      const double tol = 1.0e-2;
      if( r1 < -tol || r2 < -tol || r3 < -tol ) continue;	// need tol  ( compare with fabs(v1+v2+v3)? )
      double end[3] = {
        pXYZs_[i1*3+0]*r1 + pXYZs_[i2*3+0]*r2 + pXYZs_[i3*3+0]*r3,
        pXYZs_[i1*3+1]*r1 + pXYZs_[i2*3+1]*r2 + pXYZs_[i3*3+1]*r3,
        pXYZs_[i1*3+2]*r1 + pXYZs_[i2*3+2]*r2 + pXYZs_[i3*3+2]*r3 };
      double dir2[3] = { end[0]-org[0], end[1]-org[1], end[2]-org[2] };
      double dotdir = Com::Dot3D(dir,dir2);
      if( dotdir < 0 ) continue;
      if( !iflg ){ 
        min_dist = dotdir; 
        p[0]=end[0];  p[1]=end[1]; p[2]=end[2];        
      }
      else if( dotdir < min_dist ){
        min_dist = dotdir;
        p[0]=end[0];  p[1]=end[1]; p[2]=end[2];
      }
      iflg = true;                
    }
    if( iflg ) return true;
  }

  return false;
}

////////////////

