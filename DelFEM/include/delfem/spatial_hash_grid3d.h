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

#if !defined(SPATIAL_HASH_GRID_3D_H)
#define SPATIAL_HASH_GRID_3D_H

#include <math.h>
#include <vector>

class CSpatialHash_Grid3D
{
public:
  CSpatialHash_Grid3D(unsigned int ndiv, double center[3], double half_width);
  inline void GetIndex(const double p[3], int ip[3]) const
  {
    ip[0] = (int)floor((p[0]-org_[0])*invcellwidth_);
    ip[1] = (int)floor((p[1]-org_[1])*invcellwidth_);
    ip[2] = (int)floor((p[2]-org_[2])*invcellwidth_);
  }
  void AddTri(int itri, double p0[3], double p1[3], double p2[3]);
  void Find_NearestTriCand(double p[3], std::vector<unsigned int>& aIndTriCand);
  void Find_IntersecTriCand(const double p[3], const double d[3], std::vector<unsigned int>& aIndTriCand);
  void BuildOutFlg();
  double GetWidth() const { return width_; }
  // true:out  false:not obvious
  bool IsOut(double p[3]) const;
private:
  void AddData(unsigned int i, unsigned int j, unsigned int k, int idata);
  void SetData(int i, int j, int k,  std::vector<unsigned int>& aIndTriCand);
  void Find_Intersec_Layer_X
  (unsigned int i,  // x_num
   const double p0[2],  // yz
   const double p1[2],  // yz
   std::vector<unsigned int>& aIndTriCand);
  void Find_Intersec_Layer_Y
  (unsigned int j,  // y_num
   const double p0[2],  // zx
   const double p1[2],  // zx
   std::vector<unsigned int>& aIndTriCand);
  void Find_Intersec_Layer_Z
  (unsigned int k,  // z_num
   const double p0[2],  // xy
   const double p1[2],  // xy
   std::vector<unsigned int>& aIndTriCand);  
private:		
  unsigned int ndiv_;
  double org_[3];
  double width_;
  double invcellwidth_;
  // size 8 each: (0:positive dat size,negative -nxt ptr) (1-7:dat)		
  // i*ndiv*ndiv+j*ndiv+k
  std::vector<int> aIndTri_;
  std::vector<unsigned int> aFlgOut;
};

#endif