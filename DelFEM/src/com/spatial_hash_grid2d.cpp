/*
 *  spatial_hash_grid2d.cpp
 *  cad
 *
 *  Created by Nobuyuki Umetani on 11/14/10.
 *  Copyright 2010 The University of Tokyo and Colubmia University. All rights reserved.
 *
 */

#include <assert.h>

#include "delfem/spatial_hash_grid2d.h"


CSpatialHash_Grid2D::CSpatialHash_Grid2D()
{
  ndiv_ = 0;
  org_[0] = 0;
  org_[1] = 0;
  width_ = 1;
  invcellwidth_ = 1;
}

CSpatialHash_Grid2D::CSpatialHash_Grid2D
(unsigned int ndiv, const double center[2], double half_width)
{
	this->ndiv_ = ndiv;
	org_[0] = center[0]-half_width;
	org_[1] = center[1]-half_width;
	width_ = half_width*2;
	invcellwidth_ = ndiv/width_;
	aIndTri_.reserve(8*ndiv*ndiv*ndiv*1.5);
	aIndTri_.resize(8*ndiv*ndiv*ndiv,0);
}

void CSpatialHash_Grid2D::AddTri(int itri, double p0[2], double p1[2], double p2[2])
{
  double min[2] = { p0[0], p0[1] };
  min[0] = ( min[0] < p1[0] ) ? min[0] : p1[0];
  min[0] = ( min[0] < p2[0] ) ? min[0] : p2[0];
  min[1] = ( min[1] < p1[1] ) ? min[1] : p1[1];
  min[1] = ( min[1] < p2[1] ) ? min[1] : p2[1];
  
  double max[2] = { p0[0], p0[1] };
  max[0] = ( max[0] > p1[0] ) ? max[0] : p1[0];
  max[0] = ( max[0] > p2[0] ) ? max[0] : p2[0];  
  max[1] = ( max[1] > p1[1] ) ? max[1] : p1[1];
  max[1] = ( max[1] > p2[1] ) ? max[1] : p2[1];  
  ////
  const double eps = width_/ndiv_*0.001;
  min[0] -= eps;
  min[1] -= eps;  
  max[0] += eps;
  max[1] += eps;  
	int imin[2]; GetIndex(min, imin);
	int imax[2]; GetIndex(max, imax);
	for(unsigned int i=imin[0];i<=imax[0];i++){
  for(unsigned int j=imin[1];j<=imax[1];j++){
    this->AddData(i,j, itri);
  }
	}  
}

void CSpatialHash_Grid2D::Find_IncludeTriCand
(double p[2], std::vector<unsigned int>& aIndTriCand)
{
  aIndTriCand.resize(0);
  int ip[2]; GetIndex(p, ip);  
  SetData(ip[0], ip[1], aIndTriCand);
}

void CSpatialHash_Grid2D::AddData(unsigned int i, unsigned int j, int idata)
{
	assert( i < ndiv_ );
	assert( j < ndiv_ );
	int ind = (i*ndiv_+j);
	for(;;){
		int isize = aIndTri_[ind*8];
		if (isize>=0 ){
			if( isize<7) {
				aIndTri_[ind*8+isize+1] = idata;
				aIndTri_[ind*8] += 1;
				return;
			}
			else if(isize==7){
				assert( aIndTri_.size() % 8 == 0 );
				unsigned int ind1 = aIndTri_.size()/8;
				assert(ind1!=0);
				aIndTri_[ind*8] = -ind1;
				aIndTri_.resize((ind1+1)*8,0);
				aIndTri_[ind1*8+0] = 1;
				aIndTri_[ind1*8+1] = idata;
				return;
			}
		}
		else {
			assert(isize<0);
			assert(aIndTri_.size()>-isize*8);
			ind = -isize;
		}		
	}  
}

void CSpatialHash_Grid2D::SetData(int i, int j, std::vector<unsigned int>& aIndTriCand)
{
  
	if( i < 0 || i >= ndiv_ ) return;
	if( j < 0 || j >= ndiv_ ) return;	
  //  std::cout << i << " " << j << " " << k << " " << aIndTri_[(i*ndiv_*ndiv_+j*ndiv_+k)*8] << std::endl;
	if( aIndTri_[(i*ndiv_+j)*8] == 0 ) return;
	// extract cell if the minimum distance between a certain point in ijk and in ip is smaller than idist+1
	int ind = (i*ndiv_+j);
	for(;;){
		int isize = aIndTri_[ind*8];
		assert( isize != 0 );
		if( isize > 0 ){
			assert( isize < 8 );	// the first is used for size
			for(unsigned int idat=0;idat<isize;idat++){
				aIndTriCand.push_back(aIndTri_[ind*8+idat+1]);
			}
			break;
		}
		else{
			for(unsigned int idat=0;idat<7;idat++){
				aIndTriCand.push_back(aIndTri_[ind*8+idat+1]);						
			}
			ind = -isize;
		}
	}	  
}
