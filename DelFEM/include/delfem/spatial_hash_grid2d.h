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

#if !defined SPATIAL_HASH_GRID_2D_H
#define SPATIAL_HASH_GRID_2D_H


#include <vector>
#include <math.h>

class CSpatialHash_Grid2D
{
public:
  CSpatialHash_Grid2D();
  CSpatialHash_Grid2D(unsigned int ndiv, const double center[2], double half_width);
  inline void GetIndex(const double p[2], int ip[2]) const
  {
    ip[0] = (int)floor((p[0]-org_[0])*invcellwidth_);
    ip[1] = (int)floor((p[1]-org_[1])*invcellwidth_);
  }
  void AddTri(int itri, double p0[2], double p1[2], double p2[2]);
  void Find_IncludeTriCand(double p[2], std::vector<unsigned int>& aIndTriCand);
private:
  void AddData(unsigned int i, unsigned int j, int idata);
  void SetData(int i, int j, std::vector<unsigned int>& aIndTriCand);
private:		
  unsigned int ndiv_;
  double org_[2];
  double width_;
  double invcellwidth_;
  // size 8 each: (0:positive dat size,negative -nxt ptr) (1-7:dat)		
  // i*ndiv*ndiv+j*ndiv+k
  std::vector<int> aIndTri_;
};

#endif