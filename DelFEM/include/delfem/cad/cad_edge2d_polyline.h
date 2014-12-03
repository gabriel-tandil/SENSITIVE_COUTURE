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
@brief ÇQéüå≥ÇbÇ`ÇcÉNÉâÉX(Cad::CCadObj2D)ÇÃÉCÉìÉ^Å[ÉtÉFÅ[ÉX
@author Nobuyuki Umetani
*/

#if !defined(CAD_EDGE_2D_POLYLINE_H)
#define CAD_EDGE_2D_POLYLINE_2D_H

#include <assert.h>
#include <vector>

#include "delfem/cad_obj2d.h"

class CTriDiaMat3;
// ÉgÉâÉXÉNÉâÉX
class CCadEdge2DPolyline
{
public:  
  CCadEdge2DPolyline();
  ~CCadEdge2DPolyline();  
  void SetCadEdge(const Cad::CCadObj2D& cad_2d, unsigned int id_e, const Com::CVector2D& pick_pos);
  void Drag(Cad::CCadObj2D& cad_2d, const Com::CVector2D& dist_pos);
private:                                                                                                  
  void SolveLinearStatic();
  void ClearMemory();
  void SetFixedBoundaryFlag(unsigned int ino, unsigned int idim){
    assert( ino < nno );
    assert( idim < 3 );
    assert( bc_flag != 0 );
    bc_flag[ino*3+idim] = 1;
  }  
  unsigned int GetSizeNode() const { return nno; }
  void ProjectPoint(double x_in, double y_in, int& idiv_min,
                    double& alpha, double& ndist, double& norm_x, double& norm_y);
  void SetDisp(unsigned int ino, unsigned int idim, double disp){
    if( ino >= nno ) return;
    assert( idim < 3 );
    assert( ut != 0 );
    ut[ino*3+idim] = disp;
  }
  void GetValueNode(unsigned int ino, double& x, double& y, double&t ) const
  {
    assert( ino < nno );
    x = ut[ino*3+0]+ini_x[ino*2+0];
    y = ut[ino*3+1]+ini_x[ino*2+1];
    t = ut[ino*3+2];
  }
private:
//  void SetInitial(const std::vector<double>& aXYs);
  inline double FindNearestPointParam_Line_Point(double pc[2], double ps[2], double pe[2])
  {
    const double es[2] = { pe[0]-ps[0], pe[1]-ps[1] };
    const double sc[2] = { ps[0]-pc[0], ps[1]-pc[1] };
    const double a = es[0]*es[0]+es[1]*es[1];
    const double b = es[0]*sc[0]+es[1]*sc[1];
    return - b/a;
  }  

private:
  Com::CVector2D pick_pos;
  int idiv_picked;
  unsigned int m_IdECad;
  ////
	unsigned int nno;	// êﬂì_êî
	double *ut;	// ïœà ÅCïœà ë¨ìx
	double *ini_x;	// èâä˙à íu
	int* bc_flag;	// ã´äEèåè
	////////////////
	// ï®ê´íl
	double EI;		// ã»Ç∞ÇÃçÑê´
	double ARho;	// éøó 
	double AE;		// êLèkÇÃçÑê´
private:
	// çÑê´çsóÒ
	double *dut, *Res;	// ïœà ÇÃëùï™ó ÅCécç∑
  CTriDiaMat3* m_mat;	// åWêîçsóÒ
};

#endif
