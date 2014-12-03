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

/*! @file
@brief プリミティブ(６面体、円筒)に対するメッシュ
@author Nobuyuki Umetani
*/

#if !defined(MESH_PREMITIVE_H)
#define MESH_PREMITIVE_H

#pragma warning( disable : 4786 )

#include <vector>
#include <map>

#include "delfem/mesh_interface.h"

////////////////////////////////////////////////

namespace Msh{

////////////////////////////////////////////////
/*!
@ingroup Msh3D
*/
//! 直方体を分割したメッシュ
class CMesh_Primitive_Hexahedra : public IMesh
{
public:
  CMesh_Primitive_Hexahedra(double l_x, double l_y, double l_z,  
                            unsigned int div_x, unsigned int div_y, unsigned int div_z ){
    this->lx = l_x;  this->ly = l_y;  this->lz = l_z;
    this->nx = div_x; this->ny = div_y; this->nz = div_z;
    trans_x = 0;    trans_y = 0;    trans_z = 0;
  }
public:
  virtual unsigned int GetDimention() const{ return 3; }
	virtual void GetInfo(unsigned int id_msh,
                       unsigned int& id_cad, unsigned int& id_msh_before_ext, unsigned int& inum_ext, 
                       int& ilayer) const{
    id_cad = 0;
    id_msh_before_ext = 0;
    inum_ext = 0;
		ilayer = 0;
  }
  void Translate( double tx, double ty, double tz ){
    trans_x = tx;
    trans_y = ty;
    trans_z = tz;
  }
  virtual void GetCoord(std::vector<double>& coord) const{
    const double theta = 3.1416*2/360*30;
    coord.resize( (nx+1)*(ny+1)*(nz+1)*3 );
    const double dx = lx / nx;
    const double dy = ly / ny;
    const double dz = lz / nz;
    for(unsigned int iz=0;iz<nz+1;iz++){
      for(unsigned int iy=0;iy<ny+1;iy++){
        for(unsigned int ix=0;ix<nx+1;ix++){
          const unsigned int ip = (nx+1)*(ny+1)*iz + (nx+1)*iy + ix;
          const double x = ix*dx;
          const double y = iy*dy;
          const double z = iz*dz;
          coord[ip*3+0] = x + trans_x - lx*0.5;
          coord[ip*3+1] = y + trans_y - ly*0.5;
          coord[ip*3+2] = z + trans_z - lz*0.5;
        }
      }
    }
  }
  virtual MSH_TYPE GetConnectivity(unsigned int id_msh, std::vector<int>& lnods) const{
    if( id_msh == 1 ){  // 内部
      lnods.resize( nx*ny*nz*8 );
      for(unsigned int iz=0;iz<nz;iz++){
        for(unsigned int iy=0;iy<ny;iy++){
          for(unsigned int ix=0;ix<nx;ix++){
            const unsigned int ie = nx*ny*iz + ny*ix + iy;
            lnods[ie*8+0] = (nx+1)*(ny+1)*iz     + (nx+1)*iy     + ix;
            lnods[ie*8+1] = (nx+1)*(ny+1)*iz     + (nx+1)*iy     + (ix+1);
            lnods[ie*8+2] = (nx+1)*(ny+1)*iz     + (nx+1)*(iy+1) + (ix+1);
            lnods[ie*8+3] = (nx+1)*(ny+1)*iz     + (nx+1)*(iy+1) + ix;
            lnods[ie*8+4] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*iy     + ix;
            lnods[ie*8+5] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*iy     + (ix+1);
            lnods[ie*8+6] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*(iy+1) + (ix+1);
            lnods[ie*8+7] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*(iy+1) + ix;
          }
        }
      }
      return HEX;
    }
    else if( id_msh == 2 ){ // 底面
      lnods.resize( nx*ny*4 );
      for(unsigned int iy=0;iy<ny;iy++){
        for(unsigned int ix=0;ix<nx;ix++){
          const unsigned int ie = ny*ix + iy;
          lnods[ie*4+0] = (nx+1)*iy     + ix;
          lnods[ie*4+1] = (nx+1)*(iy+1) + ix;
          lnods[ie*4+2] = (nx+1)*(iy+1) + (ix+1);
          lnods[ie*4+3] = (nx+1)*iy     + (ix+1);
        }
      }
      return QUAD;
    }
    else if( id_msh == 3 ){ // 上面
      lnods.resize( nx*ny*4 );
      for(unsigned int iy=0;iy<ny;iy++){
        for(unsigned int ix=0;ix<nx;ix++){
          const unsigned int ie = ny*ix + iy;
          lnods[ie*4+0] = (nx+1)*(ny+1)*nz + (nx+1)*iy     + ix;
          lnods[ie*4+1] = (nx+1)*(ny+1)*nz + (nx+1)*iy     + (ix+1);
          lnods[ie*4+2] = (nx+1)*(ny+1)*nz + (nx+1)*(iy+1) + (ix+1);
          lnods[ie*4+3] = (nx+1)*(ny+1)*nz + (nx+1)*(iy+1) + ix;
        }
      }
      return QUAD;
    }
    else if( id_msh == 4 ){ // ｙが小さい面
      lnods.resize( nx*nz*4 );
      for(unsigned int iz=0;iz<nz;iz++){
        for(unsigned int ix=0;ix<nx;ix++){
          const unsigned int ie = nx*iz + ix;
          lnods[ie*4+0] = (nx+1)*(ny+1)*iz     + ix;
          lnods[ie*4+1] = (nx+1)*(ny+1)*iz     + (ix+1);
          lnods[ie*4+2] = (nx+1)*(ny+1)*(iz+1) + (ix+1);
          lnods[ie*4+3] = (nx+1)*(ny+1)*(iz+1) + ix;
        }
      }
      return QUAD;
    }
    else if( id_msh == 5 ){ // yが大きい面
      lnods.resize( nx*nz*4 );
      for(unsigned int iz=0;iz<nz;iz++){
        for(unsigned int ix=0;ix<nx;ix++){
          const unsigned int ie = nx*iz + ix;
          lnods[ie*4+0] = (nx+1)*(ny+1)*iz     + (nx+1)*ny + ix;
          lnods[ie*4+1] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*ny + ix;
          lnods[ie*4+2] = (nx+1)*(ny+1)*(iz+1) + (nx+1)*ny + (ix+1);
          lnods[ie*4+3] = (nx+1)*(ny+1)*iz     + (nx+1)*ny + (ix+1);
        }
      }
      return QUAD;
    }
    return (MSH_TYPE)0;
  }
  virtual std::vector<unsigned int> GetAry_ID() const{
    std::vector<unsigned int> ret;
    const unsigned int nid = 5; 
    for(unsigned int i=0;i<nid;i++){ ret.push_back(i+1); }
    return ret;
  }
  virtual std::vector<unsigned int> GetIncludeElemIDAry(unsigned int id_msh) const{
    std::vector<unsigned int> ret;
    if( id_msh == 1 ){
      ret.push_back(2);
      ret.push_back(3);
      ret.push_back(4);
      ret.push_back(5);
    }
    return ret;
  }
public:
  double lx,ly,lz;	//!< length in xyz direction
  unsigned int nx,ny,nz;	//!< number of divide in xyz direction
  double trans_x, trans_y, trans_z;	//!< translation in xyz direction
};    
  
/*!
 @ingroup Msh3D
 */
//! 円筒を分割したメッシュ
class CMesh_Primitive_ThickCylinder : public IMesh
{
public:
  CMesh_Primitive_ThickCylinder(double r_min, double r_max, double l_z,  
                                unsigned int div_r, unsigned int div_t, unsigned int div_z ){
    this->rmin = r_min;  this->rmax = r_max;  this->lz = l_z;
    this->nr = div_r; this->nt = div_t; this->nz = div_z;
  }
public:
	virtual void GetInfo(unsigned int id_msh,
                       unsigned int& id_cad, unsigned int& id_msh_before_ext, unsigned int& inum_ext, 
                       int& ilayer) const{
    id_cad = 0;
    id_msh_before_ext = 0;
    inum_ext = 0;
		ilayer = 0;
  }  
  virtual unsigned int GetDimention() const{ return 3; }
  virtual void GetInfo(unsigned int id_msh,
                       unsigned int& id_cad, unsigned int& id_msh_before_ext, unsigned int& inum_ext) const{
    id_cad = 0;
    id_msh_before_ext = 0;
    inum_ext = 0;
  }
  void Translate( double tx, double ty, double tz ){
    trans_x = tx;
    trans_y = ty;
    trans_z = tz;
  }
  
  virtual void GetCoord(std::vector<double>& coord) const{
    coord.resize( (nr+1)*nt*(nz+1)*3 );
    const double dr = (rmax-rmin) / nr;
    const double dt = 3.1416*2 / nt;
    const double dz = lz / nz;
    const double theta = 0.2;//3.14*0.25;
    for(unsigned int iz=0;iz<nz+1;iz++){
      for(unsigned int it=0;it<nt;  it++){
        for(unsigned int ir=0;ir<nr+1;ir++){
          const unsigned int ip = (nr+1)*nt*iz + (nr+1)*it + ir;
          const double x = (rmin+dr*ir)*cos(dt*it);
          const double y = (rmin+dr*ir)*sin(dt*it);
          const double z = iz*dz;
          coord[ip*3+0] = x + trans_x;
          coord[ip*3+1] = y + trans_y;
          coord[ip*3+2] = z + trans_z - lz*0.5;
        }
      }
    }
  }
  virtual MSH_TYPE GetConnectivity(unsigned int id_msh, std::vector<int>& lnods) const{
    if( id_msh == 1 ){  // interior
      lnods.resize( nr*nt*nz*8 );
      for(unsigned int iz=0;iz<nz;iz++){
        for(unsigned int it=0;it<nt-1;it++){
          for(unsigned int ir=0;ir<nr;ir++){
            const unsigned int ie = nt*nr*iz + nr*it + ir;
            lnods[ie*8+0] = nt*(nr+1)*iz     + (nr+1)*it     + ir;
            lnods[ie*8+1] = nt*(nr+1)*iz     + (nr+1)*it     + (ir+1);
            lnods[ie*8+2] = nt*(nr+1)*iz     + (nr+1)*(it+1) + (ir+1);
            lnods[ie*8+3] = nt*(nr+1)*iz     + (nr+1)*(it+1) + ir;
            lnods[ie*8+4] = nt*(nr+1)*(iz+1) + (nr+1)*it     + ir;
            lnods[ie*8+5] = nt*(nr+1)*(iz+1) + (nr+1)*it     + (ir+1);
            lnods[ie*8+6] = nt*(nr+1)*(iz+1) + (nr+1)*(it+1) + (ir+1);
            lnods[ie*8+7] = nt*(nr+1)*(iz+1) + (nr+1)*(it+1) + ir;
          }
        }
        for(unsigned int ir=0;ir<nr;ir++){
          const unsigned int ie = nt*nr*iz + nr*(nt-1) + ir;
          lnods[ie*8+0] = nt*(nr+1)*iz     + (nr+1)*(nt-1) + ir;
          lnods[ie*8+1] = nt*(nr+1)*iz     + (nr+1)*(nt-1) + (ir+1);
          lnods[ie*8+2] = nt*(nr+1)*iz     + (nr+1)*0      + (ir+1);
          lnods[ie*8+3] = nt*(nr+1)*iz     + (nr+1)*0      + ir;
          lnods[ie*8+4] = nt*(nr+1)*(iz+1) + (nr+1)*(nt-1) + ir;
          lnods[ie*8+5] = nt*(nr+1)*(iz+1) + (nr+1)*(nt-1) + (ir+1);
          lnods[ie*8+6] = nt*(nr+1)*(iz+1) + (nr+1)*0      + (ir+1);
          lnods[ie*8+7] = nt*(nr+1)*(iz+1) + (nr+1)*0      + ir;
        }
      }
      return HEX;
    }
    else if( id_msh == 2 ){ // bottom
      lnods.resize( nr*nt*4 );
      for(unsigned int it=0;it<nt-1;it++){
        for(unsigned int ir=0;ir<nr;ir++){
          const unsigned int ie = nr*it + ir;
          lnods[ie*4+0] = nt*(nr+1)*nz + (nr+1)*it     + ir;
          lnods[ie*4+1] = nt*(nr+1)*nz + (nr+1)*(it+1) + ir;
          lnods[ie*4+2] = nt*(nr+1)*nz + (nr+1)*(it+1) + (ir+1);
          lnods[ie*4+3] = nt*(nr+1)*nz + (nr+1)*it     + (ir+1);
        }
      }
      for(unsigned int ir=0;ir<nr;ir++){
        const unsigned int ie = nr*(nt-1) + ir;
        lnods[ie*4+0] = nt*(nr+1)*nz + (nr+1)*(nt-1) + ir;
        lnods[ie*4+1] = nt*(nr+1)*nz + (nr+1)*0      + ir;
        lnods[ie*4+2] = nt*(nr+1)*nz + (nr+1)*0      + (ir+1);
        lnods[ie*4+3] = nt*(nr+1)*nz + (nr+1)*(nt-1) + (ir+1);
      }
      return QUAD;
    }
    else if( id_msh == 3 ){ // top
      lnods.resize( nr*nt*4 );
      for(unsigned int it=0;it<nt-1;it++){
        for(unsigned int ir=0;ir<nr;ir++){
          const unsigned int ie = nr*it + ir;
          lnods[ie*4+0] = (nr+1)*it     + ir;
          lnods[ie*4+1] = (nr+1)*it     + (ir+1);
          lnods[ie*4+2] = (nr+1)*(it+1) + (ir+1);
          lnods[ie*4+3] = (nr+1)*(it+1) + ir;
        }
      }
      for(unsigned int ir=0;ir<nr;ir++){
        const unsigned int ie = nr*(nt-1) + ir;
        lnods[ie*4+0] = (nr+1)*(nt-1) + ir;
        lnods[ie*4+1] = (nr+1)*(nt-1) + (ir+1);
        lnods[ie*4+2] = (nr+1)*0      + (ir+1);
        lnods[ie*4+3] = (nr+1)*0      + ir;
      }
      return QUAD;
    }
    return (MSH_TYPE)0;
  }
  virtual std::vector<unsigned int> GetAry_ID() const{
    std::vector<unsigned int> ret;
    const unsigned int nid = 3;
    for(unsigned int i=0;i<nid;i++){ ret.push_back(i+1); }
    return ret;
  }
  virtual std::vector<unsigned int> GetIncludeElemIDAry(unsigned int id_msh) const{
    std::vector<unsigned int> ret;
    if( id_msh == 1 ){
      ret.push_back(2);
      ret.push_back(3);
    }
    return ret;
  }
public:
  double rmin, rmax, lz;
  unsigned int nr, nt, nz;
  double trans_x, trans_y, trans_z;
};  
  
// @}
};


  
#endif
