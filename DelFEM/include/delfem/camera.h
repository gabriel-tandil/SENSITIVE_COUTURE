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

// AUTHOR
// Nobuyuki Umetani

// DESCRIPTION
// This file implements the class CCamera. 
// This class have a data represents model-view transformation.
// This class doesn't have actual transformation function so as to keep indipendency with OpenGL
// If you want to apply transformation use functions in drawer_gl_utility.h

#if !defined(CAMERA_H)
#define CAMERA_H

#include <assert.h>

#include "delfem/vector3d.h"
#include "delfem/quaternion.h"

#include <iostream>

namespace Com{
namespace View{

//! modes of rotation transformation
enum ROTATION_MODE{	
	ROT_2D,		//!< 2dim rotation
	ROT_2DH,	//!< z axis is allways pararell to the upright direction of screan
	ROT_3D		//!< track ball rotation
};

class CCamera
{
public:
	CCamera()
	{
		win_aspect = 1.0; 

		half_view_height = 1.0;
		inv_scale = 1.0;
		win_center_x = 0.0;  win_center_y = 0.0;

		obj_w = 1.0;  obj_h = 1.0;  obj_d = 1.0;
		obj_center.SetZero();
		
		rot_mode = ROT_2D;
		theta = 0;
    phi = -60.0*3.1415926/180.0;  
    rot_quat.SetUnit();
		
		is_pers = false;
		fov_y = 30.0 * 3.1415926 / 180.0; 
		dist = half_view_height / tan( fov_y*0.5 ) + obj_d * 0.5;
	}
	CCamera(ROTATION_MODE rotation_mode){
		win_aspect = 1.0; 

		half_view_height = 1.0;
		inv_scale = 1.0;
		win_center_x = 0.0;  win_center_y = 0.0;

		obj_w = 1.0;  obj_h = 1.0;  obj_d = 1.0;
		obj_center.SetZero();

		rot_mode = rotation_mode;		
		theta = 0.0;  
    phi = -60.0*3.1415926/180.0;  rot_quat.SetUnit();
		
		is_pers = false;
		fov_y = 30.0 * 3.1415926 / 180.0; 
		dist = half_view_height / tan( fov_y*0.5 ) + obj_d * 0.5;
	}

	void GetObjectCenter(double& x, double& y, double& z) const 
	{ x=obj_center.x; y=obj_center.y; z=obj_center.z; }
  CVector3D GetObjectCenter() const { return obj_center; }
	// Set the center of object witch is used for Fit
	void SetObjectCenter(const double& x, const double& y, const double& z)
	{ obj_center.x=x; obj_center.y=y; obj_center.z=z; }

	// Set the size of object in xyz direction
	void SetObjectSize(double w, double h, double d){ obj_w = w;  obj_h = h;  obj_d = d; }

	void SetObjectBoundingBox(const CBoundingBox3D& bb){
		obj_center.x = (bb.x_min + bb.x_max)*0.5;
		obj_center.y = (bb.y_min + bb.y_max)*0.5;
		obj_center.z = (bb.z_min + bb.z_max)*0.5;
		obj_w = bb.x_max - bb.x_min;
		obj_h = bb.y_max - bb.y_min;
		obj_d = bb.z_max - bb.z_min;
	}

	void Fit(){
		const double margin = 1.5;
		const double obj_aspect = obj_w / obj_h;
		inv_scale = 1.0;
		if( obj_aspect < win_aspect ){ half_view_height = obj_h*0.5*margin; }
		else{
			const double tmp_h = obj_w / win_aspect;
			half_view_height = tmp_h*0.5*margin;
		}		
		dist = half_view_height * inv_scale / tan( fov_y*0.5 ) + obj_d * 0.5;
		win_center_x = 0.0; win_center_y = 0.0;
	}
	void Fit(const Com::CBoundingBox3D& bb ){
		obj_center.x = (bb.x_min + bb.x_max)*0.5;
		obj_center.y = (bb.y_min + bb.y_max)*0.5;
		obj_center.z = (bb.z_min + bb.z_max)*0.5;
		obj_w = bb.x_max - bb.x_min;
		obj_h = bb.y_max - bb.y_min;
		obj_d = bb.z_max - bb.z_min;
		this->Fit();
	}
	
	//! get scale :  
	double GetScale() const { return 1.0 / inv_scale; }
	//! set scale
	void SetScale(const double& scale){
		if( scale < 0.01 ){ this->inv_scale = 100; }
		else if( scale > 100.0 ){ this->inv_scale = 0.01; }
		else{ this->inv_scale = 1.0 / scale; }
		dist = half_view_height * inv_scale / tan( fov_y*0.5 ) + obj_d * 0.5;
	}

	////////////////////////
	// Perspective
	
	bool IsPers() const { return is_pers; }
	void SetIsPers(bool is_pers){
		if( this->rot_mode == ROT_2D ) return;
		this->is_pers = is_pers;
		if( is_pers ){ dist = half_view_height * inv_scale / tan( fov_y*0.5 ) + obj_d * 0.5; }
	}	
	double GetFovY() const { return fov_y*180.0/3.1415926; };	//!< Get the view angle
	//! set the view angle
	void SetFovY(const double& fov_y)
	{
		if( !is_pers ) return;
		if( fov_y < 15.0 ){ this->fov_y = 15.0; }
		else if( fov_y > 90.0 ){ this->fov_y = 90.0; }
		else{ this->fov_y = fov_y; }
		this->fov_y = this->fov_y*3.1415926/180.0;
		dist = half_view_height * inv_scale / tan( this->fov_y*0.5 ) + obj_d * 0.5;
	}
	void GetPerspective(double& fov_y, double& aspect, double& clip_near, double& clip_far) const {
		fov_y = this->fov_y*180.0/3.1415926;	
		aspect = this->win_aspect;				
		clip_near = 0.001;
		clip_far = dist*2.0+obj_d*20.0+100.0;
	}
	void GetOrtho(double& hw, double& hh, double& hd) const {
		hw = half_view_height*inv_scale*win_aspect;
		hh = half_view_height*inv_scale;
    hd = ( obj_d*20.0 > 1.0e-4 ) ? obj_d*20.0 : 1.0e-4;
    hd = ( obj_w*20.0 > hd ) ? obj_w*20.0 : hd;
    hd = ( obj_h*20.0 > hd ) ? obj_h*20.0 : hd;
	}

	////////////////
	// 
	double GetHalfViewHeight() const { return half_view_height; }

	// get the 3D location where the center is located
	void GetCenterPosition(double& x, double& y, double& z) const {
		if( is_pers ){ x = 0.0;  y = 0.0;  z = -dist; }
		else{
			x = half_view_height*win_center_x*win_aspect;
			y = half_view_height*win_center_y;
			z = 0.0; 
		}
	}
  CVector3D GetCenterPosition(){
    double x,y,z; this->GetCenterPosition(x,y,z);
    return CVector3D(x,y,z);
  }
	void MousePan(const double& mov_begin_x, const double& mov_begin_y,
		const double& mov_end_x, const double& mov_end_y ){
			win_center_x += (mov_end_x - mov_begin_x)*inv_scale;
			win_center_y += (mov_end_y - mov_begin_y)*inv_scale;
	}

	void SetWindowAspect(double asp){ win_aspect = asp; }	// set window aspect ratio
	double GetWindowAspect() const { return win_aspect; }	// get window aspect ratio

	////////////////
	// rotation
	void RotMatrix33(double m[]) const {
		if( rot_mode == ROT_2D ){
			const double ct = cos(theta);
			const double st = sin(theta);
			m[0] =  ct; m[1] = -st; m[2] =  0.0;
			m[3] =  st; m[4] =  ct; m[5] =  0.0;
			m[6] = 0.0; m[7] = 0.0; m[8] =  1.0;
		}
		else if( rot_mode == ROT_2DH ){
			const double ct = cos(theta);
			const double st = sin(theta);
			const double cp = cos(phi);
			const double sp = sin(phi);
			m[0] =    ct; m[1] =   -st; m[2] = 0.0;
			m[3] = cp*st; m[4] = cp*ct; m[5] = -sp;
			m[6] = sp*st; m[7] = sp*ct; m[8] =  cp;
		}
		else if( rot_mode == ROT_3D ){  rot_quat.RotMatrix33(m);  }
	}
  CMatrix3 RotMatrix33() const{
    double r[9];  this->RotMatrix33(r);
    return CMatrix3(r);
  }
	void RotMatrix44Trans(double m[]) const {
		double m1[9];
		this->RotMatrix33(m1);
		m[0] = m1[0];  m[1] = m1[3];  m[2] = m1[6];  m[3] = 0.0;
		m[4] = m1[1];  m[5] = m1[4];  m[6] = m1[7];  m[7] = 0.0;
		m[8] = m1[2];  m[9] = m1[5];  m[10]= m1[8];  m[11]= 0.0;
		m[12]=   0.0;  m[13]=   0.0;  m[14]=   0.0;  m[15]= 1.0;
	}
	void SetRotationMode(ROTATION_MODE rot_mode){ this->rot_mode = rot_mode; }
	ROTATION_MODE GetRotationMode(){ return this->rot_mode; }
	void MouseRotation(const double& mov_begin_x, const double& mov_begin_y,
		const double& mov_end_x, const double& mov_end_y ){
		if( rot_mode == ROT_3D  ){
			if( fabs(mov_begin_x-mov_end_x) < 1.0e-8 && fabs(mov_begin_y-mov_end_y) < 1.0e-8 ){ return; }
			CVector3D p1;
			{
				const double norm2 = mov_begin_x*mov_begin_x+mov_begin_y*mov_begin_y;
				double z;
				if( norm2 < 0.5 ){ z = sqrt(1.0-norm2); }
				else{ z = 0.343443391 / norm2; } 
				p1.SetVector(mov_begin_x,mov_begin_y,z);
			}
			CVector3D p2;
			{
				const double norm2 = mov_end_x*mov_end_x+mov_end_y*mov_end_y;
				double z;
				if( norm2 < 0.5 ){ z = sqrt(1.0-norm2); }
				else{ z = 0.343443391 / norm2; } 
				p2.SetVector(mov_end_x,mov_end_y,z);
			}
			CVector3D norm = Cross(p1,p2);
			norm.Normalize();
			p2 -= p1;
			double t = p2.Length() / 2.0;
			if( t > 1.0 ){ t = 1.0; }
			else if( t < -1.0 ){ t = -1.0; }
			const double phi = 5.0*asin(t);
			CQuaternion tmp_quat = CQuaternion(norm*phi)*rot_quat;
			rot_quat = tmp_quat;
			rot_quat.Normalize();
		}
    else if( rot_mode == ROT_2DH ){
      const double tmp = 3.1415926/180.0;
      theta += (mov_end_x-mov_begin_x)*tmp*100.0;
      phi -= (mov_end_y-mov_begin_y)*tmp*100.0;
      if( phi > +0.0*tmp ){ phi = +0.0*tmp; }
      else if( phi < -180.0*tmp ){ phi = -180.0*tmp; }
		}
		else if( rot_mode == ROT_2D ){
			{
				const double diff_x = mov_end_x - mov_begin_x;
				const double diff_y = mov_end_y - mov_begin_y;
				const double q_len_diff = diff_x*diff_x + diff_y*diff_y;
				if( fabs(q_len_diff) < 1.0e-16 ) return;
			}
			const double cnt_x = win_center_x;
			const double cnt_y = win_center_y;
			const double q_len_end = (mov_end_x-cnt_x)*(mov_end_x-cnt_x) + (mov_end_y-cnt_y)*(mov_end_y-cnt_y);
			const double q_len_begin = (mov_begin_x-cnt_x)*(mov_begin_x-cnt_x) + (mov_begin_y-cnt_y)*(mov_begin_y-cnt_y);
			const double d_sin = ((mov_end_x-cnt_x)*(mov_begin_y-cnt_y)-(mov_end_y-cnt_y)*(mov_begin_x-cnt_x))/sqrt(q_len_end*q_len_begin);
			double d_theta = asin(d_sin);
			theta -= d_theta;
		}
	}
	
  Com::CVector3D ProjectionOnPlane(
		const double pos_x,   const double pos_y,
		double plane_x=0, double plane_y=0, double plane_z=0,
		double norm_x=0,  double norm_y=0,  double norm_z=1 ) const
	{
		assert( fabs(plane_x) < 1.0e-10 );
		assert( fabs(plane_y) < 1.0e-10 );
		assert( fabs(plane_z) < 1.0e-10 );
		
		assert( fabs(norm_x) < 1.0e-10 );
		assert( fabs(norm_y) < 1.0e-10 );
		assert( fabs(norm_z-1.0) < 1.0e-10 );

		assert( this->rot_mode == ROT_2D );

		const double hw = half_view_height*win_aspect;
		const double hh = half_view_height;

		const double obj_cnt_x=obj_center.x, obj_cnt_y=obj_center.y, obj_cnt_z=obj_center.z;

    const double ox =  hw*cos(theta)*(pos_x*inv_scale-win_center_x)+hh*sin(theta)*(pos_y*inv_scale-win_center_y)+obj_cnt_x;
    const double oy = -hw*sin(theta)*(pos_x*inv_scale-win_center_x)+hh*cos(theta)*(pos_y*inv_scale-win_center_y)+obj_cnt_y;
    const double oz = +obj_cnt_z;

    return Com::CVector3D(ox,oy,oz);
	}
  void SetRotation2DH(double theta, double phi){
    this->theta = theta;
    this->phi = phi;
  }
  void SetWindowCenter(double winx, double winy){
    this->win_center_x = winx; 
    this->win_center_y = winy;   
  }
  /*
  void SetBryantAngle(double phi, double theta, double psi){
    rot_quat = Com::CQuaternion();
    rot_quat *= CQuaternion( CVector3D(phi,0,0) );
    rot_quat *= CQuaternion( CVector3D(0,theta,0) );
    rot_quat *= CQuaternion( CVector3D(0,0,psi) );
  }
   */


private:
	double win_aspect;

	double win_center_x, win_center_y;
	double inv_scale;
	double half_view_height;

	double obj_w, obj_h, obj_d;
	CVector3D obj_center;

	ROTATION_MODE rot_mode;
	CQuaternion rot_quat;
	double theta;
	double phi;

	bool is_pers;
	double fov_y;
	double dist;
};

}	// end namespace View
}	// end namespace Com

#endif
