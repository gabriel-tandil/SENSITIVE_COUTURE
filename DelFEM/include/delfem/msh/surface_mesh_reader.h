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

#if !defined(SURFACE_MESH_READER_H)
#define SURFACE_MESH_READER_H

#include <string>

class CSurfaceMeshReader
{
public:  
  CSurfaceMeshReader(){
    nnode_ = 0;
    pXYZs_ = 0;
    ntri_  = 0;
    aTri_  = 0;
    aNorm_ = 0;
    is_norm_ = false;
  }
  void Clear(){
    if( pXYZs_ != 0 ){ delete[] pXYZs_; pXYZs_ = 0; }
    if( aTri_  != 0 ){ delete[] aTri_;  aTri_  = 0; }
    if( aNorm_ != 0 ){ delete[] aNorm_; aNorm_ = 0; }        
    nnode_ = 0;
    ntri_  = 0;    
  }
  /////
  void Load_Off(const std::string& fname);
  void Load_Gmv(const std::string& fname);
  void Load_Ply(const std::string& fname);    
  ////
  void Translate(double tx, double ty, double tz);
  void Scale(double s);
  void Rot_Bryant(double rx, double ry, double rz);
  ////
  void GetCenterWidth(double& cx, double& cy, double& cz, 
                      double& wx, double& wy, double& wz) const;  
  void SetIsNormal(bool is_norm);
  void Draw() const;
  void GetMesh(std::vector<unsigned int>& aTri, std::vector<double>& aXYZ) const;
  bool IsEmpty() const{
    if( nnode_ == 0 && ntri_ == 0 ) return true;
    return false;
  }
  bool WriteSTL(const std::string& fname,double scale) const;
private:
  void MakeNormal();
private:
  unsigned int nnode_;
  double* pXYZs_;  
  unsigned int ntri_;
	unsigned int* aTri_;  
  bool is_norm_;
  double* aNorm_;  
};

#endif