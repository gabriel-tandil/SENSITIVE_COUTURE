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

#include "delfem/drawer.h"
#include "delfem/camera.h"
#include "delfem/vector3d.h"
#include "delfem/quaternion.h"

using namespace Com::View;

Com::CBoundingBox3D CDrawerArray::GetBoundingBox( double rm[] ) const
{
	if( m_drawer_ary.empty() ){
		return CBoundingBox3D(-0.5,0.5, -0.5,0.5, -0.5,0.5);
	}	
    CBoundingBox3D bb = m_drawer_ary[0]->GetBoundingBox(rm);
	for(unsigned int idraw=1;idraw<m_drawer_ary.size();idraw++){
        bb += m_drawer_ary[idraw]->GetBoundingBox(rm);
	}
	return bb;
}

// rot is 3 by 3 matrix for rotation
// if rot is 0 this function don't perform rotation in measuring size
Com::CBoundingBox3D CVertexArray::GetBoundingBox( double rot[] ) const
{
	if( pVertexArray == 0 ){ return Com::CBoundingBox3D(); }
	if( rot == 0 )  // object axis alligned bounding box
  {
		if( ndim == 2 ){		
			Com::CBoundingBox3D bb;
			{
				const double x1 = pVertexArray[0];
				const double y1 = pVertexArray[1];
				const double z1 = 0.0;
				bb = Com::CBoundingBox3D(x1,x1, y1,y1, z1,z1);
			}
			for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
				const double x1 = pVertexArray[ipoin*2  ];
				const double y1 = pVertexArray[ipoin*2+1];
				const double z1 = 0.0;
				bb.x_max = ( x1 > bb.x_max ) ? x1 : bb.x_max;  bb.x_min = ( x1 < bb.x_min ) ? x1 : bb.x_min;			
				bb.y_max = ( y1 > bb.y_max ) ? y1 : bb.y_max;  bb.y_min = ( y1 < bb.y_min ) ? y1 : bb.y_min;			
				bb.z_max = ( z1 > bb.z_max ) ? z1 : bb.z_max;  bb.z_min = ( z1 < bb.z_min ) ? z1 : bb.z_min;
			}
			return bb;
		}
		if( ndim == 3 ){
			Com::CBoundingBox3D bb;
			{
				const double x1 = pVertexArray[0];
				const double y1 = pVertexArray[1];
				const double z1 = 0.0;
				bb = Com::CBoundingBox3D(x1,x1, y1,y1, z1,z1);
			}
			for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
				const double x1 = pVertexArray[ipoin*3  ];
				const double y1 = pVertexArray[ipoin*3+1];
				const double z1 = pVertexArray[ipoin*3+2];
				bb.x_max = ( x1 > bb.x_max ) ? x1 : bb.x_max;  bb.x_min = ( x1 < bb.x_min ) ? x1 : bb.x_min;
				bb.y_max = ( y1 > bb.y_max ) ? y1 : bb.y_max;  bb.y_min = ( y1 < bb.y_min ) ? y1 : bb.y_min;
				bb.z_max = ( z1 > bb.z_max ) ? z1 : bb.z_max;  bb.z_min = ( z1 < bb.z_min ) ? z1 : bb.z_min;
			}
			return bb;
		}		
	}
	if( ndim == 2 ) // view axis alligned bounding box
  {		
    double x_min,x_max, y_min,y_max, z_min,z_max;
		{
			const double x1 = pVertexArray[0];
			const double y1 = pVertexArray[1];
			const double z1 = 0.0;
			x_min = x_max = x1*rot[0]+y1*rot[1]+z1*rot[2];
			y_min = y_max = x1*rot[3]+y1*rot[4]+z1*rot[5];
			z_min = z_max = x1*rot[6]+y1*rot[7]+z1*rot[8];
		}
		for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
			const double x1 = pVertexArray[ipoin*2  ];
			const double y1 = pVertexArray[ipoin*2+1];
			const double z1 = 0.0;
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			x_max = ( x2 > x_max ) ? x2 : x_max;  x_min = ( x2 < x_min ) ? x2 : x_min;			
			y_max = ( y2 > y_max ) ? y2 : y_max;  y_min = ( y2 < y_min ) ? y2 : y_min;			
			z_max = ( z2 > z_max ) ? z2 : z_max;  z_min = ( z2 < z_min ) ? z2 : z_min;
		}    
    const double c1x = (x_min + x_max)*0.5;
    const double c1y = (y_min + y_max)*0.5;
    const double c1z = (z_min + z_max)*0.5;    
    const double c2x = c1x*rot[0]+c1y*rot[3]+c1z*rot[6];
    const double c2y = c1x*rot[1]+c1y*rot[4]+c1z*rot[7];
    const double c2z = c1x*rot[2]+c1y*rot[5]+c1z*rot[8];    
    const double hx = (x_max - x_min)*0.5;
    const double hy = (y_max - y_min)*0.5;
    const double hz = (z_max - z_min)*0.5;    
		return Com::CBoundingBox3D(c2x-hx,c2x+hx, c2y-hy,c2y+hy, c2z-hz,c2z+hz);
	}
	if( ndim == 3 ) // view axis alligned bounding box
  {
    double x_min,x_max, y_min,y_max, z_min,z_max;
    {
      const double x1 = pVertexArray[0];
      const double y1 = pVertexArray[1];
      const double z1 = pVertexArray[2];
			x_min = x_max = x1*rot[0]+y1*rot[1]+z1*rot[2];
			y_min = y_max = x1*rot[3]+y1*rot[4]+z1*rot[5];
			z_min = z_max = x1*rot[6]+y1*rot[7]+z1*rot[8];
		}
		for(unsigned int ipoin=1;ipoin<npoin;ipoin++){
			const double x1 = pVertexArray[ipoin*3  ];
			const double y1 = pVertexArray[ipoin*3+1];
			const double z1 = pVertexArray[ipoin*3+2];
			const double x2 = x1*rot[0]+y1*rot[1]+z1*rot[2];
			const double y2 = x1*rot[3]+y1*rot[4]+z1*rot[5];
			const double z2 = x1*rot[6]+y1*rot[7]+z1*rot[8];
			x_max = ( x2 > x_max ) ? x2 : x_max;  x_min = ( x2 < x_min ) ? x2 : x_min;
			y_max = ( y2 > y_max ) ? y2 : y_max;  y_min = ( y2 < y_min ) ? y2 : y_min;
			z_max = ( z2 > z_max ) ? z2 : z_max;  z_min = ( z2 < z_min ) ? z2 : z_min;
		}
    const double c1x = (x_min + x_max)*0.5;
    const double c1y = (y_min + y_max)*0.5;
    const double c1z = (z_min + z_max)*0.5;    
    const double c2x = c1x*rot[0]+c1y*rot[3]+c1z*rot[6];
    const double c2y = c1x*rot[1]+c1y*rot[4]+c1z*rot[7];
    const double c2z = c1x*rot[2]+c1y*rot[5]+c1z*rot[8];    
    const double hx = (x_max - x_min)*0.5;
    const double hy = (y_max - y_min)*0.5;
    const double hz = (z_max - z_min)*0.5;    
		return Com::CBoundingBox3D(c2x-hx,c2x+hx, c2y-hy,c2y+hy, c2z-hz,c2z+hz);
	}
	return Com::CBoundingBox3D();
}

void CDrawerArray::InitTrans(Com::View::CCamera& camera ){
	{	// get suitable rot mode
		unsigned int irot_mode = 0;
		for(unsigned int idraw=0;idraw<m_drawer_ary.size();idraw++){
			unsigned int irot_mode0 = m_drawer_ary[idraw]->GetSutableRotMode();
			irot_mode = (irot_mode0>irot_mode) ? irot_mode0 : irot_mode;
		}
		if(      irot_mode == 1 ){ camera.SetRotationMode(ROT_2D);  }
		else if( irot_mode == 2 ){ camera.SetRotationMode(ROT_2DH); }
		else if( irot_mode == 3 ){ camera.SetRotationMode(ROT_3D);  }
	}
	// set object size to the transformation
	double rot[9];	camera.RotMatrix33(rot);
	Com::CBoundingBox3D bb = this->GetBoundingBox( rot );
	camera.Fit(bb);
}
