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

#include <fstream>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <vector>
#include <cstring>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/msh/surface_mesh_reader.h"




void CSurfaceMeshReader::Load_Off(const std::string& fname)
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
  /*
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
  //    std::cout << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << std::endl;
	double max_width = (x_max-x_min) > (y_max-y_min) ? (x_max-x_min) : (y_max-y_min);
	max_width = max_width > (z_max-z_min) ? max_width : (z_max-z_min);
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x = pXYZs_[ino*3+0];
		pXYZs_[ino*3+0] = (x-(x_min+x_max)*0.5)*2.0/max_width;
		double y = pXYZs_[ino*3+1];
		pXYZs_[ino*3+1] = (y-(y_min+y_max)*0.5)*2.0/max_width;
		double z = pXYZs_[ino*3+2];
		pXYZs_[ino*3+2] = (z-(z_min+z_max)*0.5)*2.0/max_width;		
	}
   */
  if( is_norm_ ){ this->MakeNormal(); }
}

void CSurfaceMeshReader::Load_Gmv(const std::string& fname)
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
	/*		for(unsigned int ino=0;ino<nnode_;ino++){
	 double x,y,z;
	 fin >> x >> y >> z;
	 pXYZs_[ino*3+0] = x;
	 pXYZs_[ino*3+1] = y;
	 pXYZs_[ino*3+2] = z;
	 }*/
  
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
  /*
	double x_min = pXYZs_[0], x_max = pXYZs_[0];
  double y_min = pXYZs_[1], y_max = pXYZs_[1];
  double z_min = pXYZs_[2], z_max = pXYZs_[2];    
	for(unsigned int ino=0;ino<nnode_;ino++){
		x_min = ( x_min < pXYZs_[ino*3+0] ) ? x_min : pXYZs_[ino*3+0];
		x_max = ( x_max > pXYZs_[ino*3+0] ) ? x_max : pXYZs_[ino*3+0];			
		y_min = ( y_min < pXYZs_[ino*3+1] ) ? y_min : pXYZs_[ino*3+1];
		y_max = ( y_max > pXYZs_[ino*3+1] ) ? y_max : pXYZs_[ino*3+1];			
		z_min = ( z_min < pXYZs_[ino*3+2] ) ? z_min : pXYZs_[ino*3+2];
		z_max = ( z_max > pXYZs_[ino*3+2] ) ? z_max : pXYZs_[ino*3+2];
	}
	double max_width = (x_max-x_min) > (y_max-y_min) ? (x_max-x_min) : (y_max-y_min);
	max_width = max_width > (z_max-z_min) ? max_width : (z_max-z_min);
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x = pXYZs_[ino*3+0];
		pXYZs_[ino*3+0] = (x-(x_min+x_max)*0.5)*2.0/max_width;
		double y = pXYZs_[ino*3+1];
		pXYZs_[ino*3+1] = (y-(y_min+y_max)*0.5)*2.0/max_width;
		double z = pXYZs_[ino*3+2];
		pXYZs_[ino*3+2] = (z-(z_min+z_max)*0.5)*2.0/max_width;		
	}
   */
  if( is_norm_ ){ this->MakeNormal(); }  
}	


void CSurfaceMeshReader::Load_Ply(const std::string& fname)
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
  if( is_norm_ ){ this->MakeNormal(); }    
}

void CSurfaceMeshReader::GetCenterWidth(double& cx, double& cy, double& cz, 
                                        double& wx, double& wy, double& wz) const
{  
  if( pXYZs_ == 0 ){
    cx=cy=cz=0;
    wx=wy=wz=1;
    return;
  }
	double x_min = pXYZs_[0], x_max = pXYZs_[0];
  double y_min = pXYZs_[1], y_max = pXYZs_[1];
  double z_min = pXYZs_[2], z_max = pXYZs_[2];    
	for(unsigned int ino=0;ino<nnode_;ino++){
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
/*  
  std::cout << x_max-x_min << " " << y_max-y_min << " " << z_max-z_min << std::endl;
	double max_width = (x_max-x_min) > (y_max-y_min) ? (x_max-x_min) : (y_max-y_min);
	max_width = max_width > (z_max-z_min) ? max_width : (z_max-z_min);
  std::cout << max_width << std::endl;
	for(unsigned int ino=0;ino<nnode_;ino++){
		double x = pXYZs_[ino*3+0];
		pXYZs_[ino*3+0] = (x-(x_min+x_max)*0.5)*2.0/max_width;
		double y = pXYZs_[ino*3+1];
		pXYZs_[ino*3+1] = (y-(y_min+y_max)*0.5)*2.0/max_width;
		double z = pXYZs_[ino*3+2];
		pXYZs_[ino*3+2] = (z-(z_min+z_max)*0.5)*2.0/max_width;		
	}    
 */
}

static void MatVec3(const double m[9], const double x[3], double y[3]){
  y[0] = m[0]*x[0] + m[1]*x[1] + m[2]*x[2];
  y[1] = m[3]*x[0] + m[4]*x[1] + m[5]*x[2];
  y[2] = m[6]*x[0] + m[7]*x[1] + m[8]*x[2];
}

static void VecMat3(const double x[3], const double m[9],  double y[3]){
  y[0] = m[0]*x[0] + m[3]*x[1] + m[6]*x[2];
  y[1] = m[1]*x[0] + m[4]*x[1] + m[7]*x[2];
  y[2] = m[2]*x[0] + m[5]*x[1] + m[8]*x[2];
}

static inline double Length3D(const double v[3]){
  return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

inline double Length3D(const double p0[3], const double p1[3]){
	return sqrt( (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]) );
}

inline void  UnitNormalAreaTri3D(double n[3], double& a, const double v1[3], const double v2[3], const double v3[3]){
	n[0] = ( v2[1] - v1[1] )*( v3[2] - v1[2] ) - ( v3[1] - v1[1] )*( v2[2] - v1[2] );
	n[1] = ( v2[2] - v1[2] )*( v3[0] - v1[0] ) - ( v3[2] - v1[2] )*( v2[0] - v1[0] );
	n[2] = ( v2[0] - v1[0] )*( v3[1] - v1[1] ) - ( v3[0] - v1[0] )*( v2[1] - v1[1] );
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])*0.5;
	const double invlen = 0.5/a;
	n[0]*=invlen;	n[1]*=invlen;	n[2]*=invlen;
}

void CSurfaceMeshReader::MakeNormal()
{
  if( aNorm_ != 0 ){ delete[] aNorm_; }
  aNorm_ = new double [nnode_*3];
  for(unsigned int i=0;i<nnode_*3;i++){ aNorm_[i] = 0; }
  for(unsigned int itri=0;itri<ntri_;itri++){
    unsigned int i1 = aTri_[itri*3+0];
    unsigned int i2 = aTri_[itri*3+1];
    unsigned int i3 = aTri_[itri*3+2];
    double un[3], area;    
    UnitNormalAreaTri3D(un,area, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);    
    aNorm_[i1*3+0] += un[0];  aNorm_[i1*3+1] += un[1];  aNorm_[i1*3+2] += un[2];
    aNorm_[i2*3+0] += un[0];  aNorm_[i2*3+1] += un[1];  aNorm_[i2*3+2] += un[2];    
    aNorm_[i3*3+0] += un[0];  aNorm_[i3*3+1] += un[1];  aNorm_[i3*3+2] += un[2];    
  }
  for(unsigned int ino=0;ino<nnode_;ino++){
    double invlen = 1.0/Length3D(aNorm_+ino*3);
    aNorm_[ino*3+0] *= invlen;
    aNorm_[ino*3+1] *= invlen;
    aNorm_[ino*3+2] *= invlen;    
  }  
}


void CSurfaceMeshReader::Draw() const
{  
  bool is_lighting = ::glIsEnabled(GL_LIGHTING);
  ::glEnable(GL_LIGHTING);
  bool is_texture  = ::glIsEnabled(GL_TEXTURE_2D);
  ::glDisable(GL_TEXTURE_2D);    
  {
    float gray[4] = {0.3,0.3,0.3,1};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, gray);
    float shine[4] = {0,0,0,0};
    ::glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shine);
    ::glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 127.0);
    //    ::glColor3d(1,1,1);
  }
  if( is_norm_ ){
    ::glBegin(GL_TRIANGLES);
    for(unsigned int itri=0;itri<ntri_;itri++){
      const unsigned int i1 = aTri_[itri*3+0];
      const unsigned int i2 = aTri_[itri*3+1];
      const unsigned int i3 = aTri_[itri*3+2];
      ::glNormal3dv(aNorm_+i1*3);   ::glVertex3dv(pXYZs_+i1*3);
      ::glNormal3dv(aNorm_+i2*3);   ::glVertex3dv(pXYZs_+i2*3);
      ::glNormal3dv(aNorm_+i3*3);   ::glVertex3dv(pXYZs_+i3*3);
    }
    ::glEnd();    
  }
  else{
    ::glDisable(GL_LIGHTING);
    ::glBegin(GL_TRIANGLES);
    for(unsigned int itri=0;itri<ntri_;itri++){
      const unsigned int i1 = aTri_[itri*3+0];
      const unsigned int i2 = aTri_[itri*3+1];
      const unsigned int i3 = aTri_[itri*3+2];
      ::glVertex3dv(pXYZs_+i1*3);
      ::glVertex3dv(pXYZs_+i2*3);
      ::glVertex3dv(pXYZs_+i3*3);
    }
    ::glEnd();        
  }
  if( is_lighting ){ ::glEnable( GL_LIGHTING); }
  else{              ::glDisable(GL_LIGHTING); }
  if(  is_texture  ){ ::glEnable(GL_TEXTURE_2D); } 
}

void CSurfaceMeshReader::SetIsNormal(bool is_norm)
{
  if( is_norm == is_norm_ ) return;
  is_norm_ = is_norm;
  if( is_norm ){    
    this->MakeNormal();
  }
  else{
    delete[] aNorm_;
    aNorm_ = 0;
  }  
}

void CSurfaceMeshReader::Translate(double tx, double ty, double tz)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] += tx;
    pXYZs_[ino*3+1] += ty;
    pXYZs_[ino*3+2] += tz;
  }
}

void CSurfaceMeshReader::Scale(double s)
{
  for(unsigned int ino=0;ino<nnode_;ino++){
    pXYZs_[ino*3+0] *= s;
    pXYZs_[ino*3+1] *= s;
    pXYZs_[ino*3+2] *= s;
  }  
}

void CSurfaceMeshReader::Rot_Bryant(double phi, double theta, double psi)
{
  phi   *= 3.1416/180.0;
  theta *= 3.1416/180.0;
  psi   *= 3.1416/180.0;  
  ////
  const double mat[9] = {
    cos(psi)*cos(theta),	cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi),
    sin(psi)*cos(theta),	sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi),
    -sin(theta),		    	cos(theta)*sin(phi),							              cos(theta)*cos(phi)};  
  ////  
  double res[3];  
  for(unsigned int ino=0;ino<nnode_;ino++){
    double* p = pXYZs_+ino*3;
    MatVec3(mat,p,res);
    p[0]=res[0]; p[1]=res[1]; p[2]=res[2];
  }
  if( is_norm_ && aNorm_!=0 ){
    for(unsigned int ino=0;ino<nnode_;ino++){  
      double* n = aNorm_+ino*3;
      MatVec3(mat,n,res);
      n[0]=res[0]; n[1]=res[1]; n[2]=res[2];
    }
  }
}


void CSurfaceMeshReader::GetMesh
(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ) const
{
  aTri.resize(ntri_*3);
  for(unsigned int i=0;i<ntri_*3;i++){ aTri[i] = this->aTri_[i]; }
  aXYZ.resize(nnode_*3);
  for(unsigned int i=0;i<nnode_*3;i++){ aXYZ[i] = this->pXYZs_[i]; }
}

bool CSurfaceMeshReader::WriteSTL(const std::string& fname,double scale) const
{
	FILE *fp;
	if( (fp = ::fopen(fname.c_str(),"w"))== NULL ){
		fclose(fp);
		assert(0);
		return false;
	}
  fprintf(fp, "solid cloth_pattern\n");
  for(unsigned int itri=0;itri<ntri_;itri++){
    unsigned int i1 = aTri_[itri*3+0];    
    unsigned int i2 = aTri_[itri*3+1];
    unsigned int i3 = aTri_[itri*3+2];
    double un[3], area;    
    UnitNormalAreaTri3D(un,area, pXYZs_+i1*3, pXYZs_+i2*3, pXYZs_+i3*3);    
    fprintf(fp, "facet normal %lf %lf %lf\n",un[0],un[1],un[2]);
    fprintf(fp, "outer loop\n");
    fprintf(fp, "vertex %lf %lf %lf\n",pXYZs_[i1*3+0]*scale,pXYZs_[i1*3+1]*scale,pXYZs_[i1*3+2]*scale); 
    fprintf(fp, "vertex %lf %lf %lf\n",pXYZs_[i2*3+0]*scale,pXYZs_[i2*3+1]*scale,pXYZs_[i2*3+2]*scale); 
    fprintf(fp, "vertex %lf %lf %lf\n",pXYZs_[i3*3+0]*scale,pXYZs_[i3*3+1]*scale,pXYZs_[i3*3+2]*scale);     
    fprintf(fp, "endloop\n");    
    fprintf(fp, "endfacet\n");        
  }    
  fprintf(fp,"endsolid\n");
  fclose(fp);
  return true;
}

