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

#define for if(0);else for

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif
#endif  /* _WIN32 */

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#include <math.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>	// strlen
#include <stdio.h>	// sprintf
#include <cstdlib>	// atoi


#include "delfem/camera.h"
#include "delfem/uglyfont.h"
#include "delfem/quaternion.h"
#include "delfem/drawer_gl_utility.h"

using namespace Com::View;

struct SPickedObject{
	unsigned int name_depth;
	int name[4];
	double min_depth;
	double max_depth;
};

void Com::View::PickPre(
		unsigned int size_buffer, unsigned int* select_buffer,
		unsigned int point_x, unsigned int point_y,
		unsigned int delX, unsigned int delY,
		const View::CCamera& mvp_trans)
//		int win_width, int win_height)
{
    // Initailze Selection
    glSelectBuffer(size_buffer, select_buffer);
    glRenderMode(GL_SELECT);

    // Get View Port
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    glInitNames();

    // Projection Transform From Here
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPickMatrix(point_x,viewport[3]-point_y,  delX,delY,  viewport);
	View::SetProjectionTransform(mvp_trans);

    // Model-View  Transform From Here
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	View::SetModelViewTransform(mvp_trans);

}

std::vector<SSelectedObject> Com::View::PickPost(
		unsigned int* const select_buffer,
		unsigned int point_x, unsigned int point_y,
		const View::CCamera& mvp_trans)
//		int win_width, int win_height)
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

	std::vector<SSelectedObject> aSelecObj;

    int nhits = glRenderMode(GL_RENDER);    // return value is number of hits
    if(nhits<=0){ return aSelecObj; }

    std::vector<SPickedObject> picked_object;
    {   // get picked_object name and its depth
        picked_object.resize(nhits);
        unsigned int* ptr = select_buffer;
		for(unsigned int i=0; i<picked_object.size(); i++){
			const unsigned int name_depth = static_cast<unsigned int>(*ptr);
			assert(name_depth<=4);
			picked_object[i].name_depth = name_depth;
			ptr++;
			picked_object[i].min_depth = (float) *ptr/0x7fffffff;
			ptr++;
			picked_object[i].max_depth = (float) *ptr/0x7fffffff;
			ptr++;
			for(unsigned int j=0; j<name_depth; j++){
				picked_object[i].name[j] = *ptr;
				ptr++;
			}
		}
		// sort picked object in the order of min depth
		for(unsigned int i=0;i<picked_object.size();i++){
			for(unsigned int j=i+1;j<picked_object.size();j++){
				if( picked_object[i].min_depth > picked_object[j].min_depth ){
					SPickedObject tmp = picked_object[i];
					picked_object[i] = picked_object[j];
					picked_object[j] = tmp;
				}
			}
		}
	}
/*
	std::cout << "Picked Object " << nhits << std::endl;
	for(unsigned int i=0; i<picked_object.size(); i++){
		std::cout << i << " ";
		std::cout << picked_object[i].name_depth << " ";
		std::cout << picked_object[i].min_depth << " ";
		std::cout << picked_object[i].max_depth << std::endl;
		std::cout << "   ";
		for(unsigned int j=0; j<picked_object[i].name_depth; j++){
			printf("%d ",picked_object[i].name[j]);
		}
		printf("\n");
	}
*/
	aSelecObj.clear();
	for(unsigned int i=0; i<picked_object.size(); i++){
		aSelecObj.resize( aSelecObj.size()+1 );
		SSelectedObject& selec_obj = aSelecObj[ aSelecObj.size()-1 ];
		assert(picked_object[i].name_depth<=4);
		selec_obj.name_depth = 3;
		selec_obj.name[0] = picked_object[i].name[0];
		selec_obj.name[1] = picked_object[i].name[1];
		selec_obj.name[2] = picked_object[i].name[2];

		double ox,oy,oz;
		{   
			GLdouble mvMatrix[16],pjMatrix[16];
			GLint viewport[4]; 
			::glGetIntegerv(GL_VIEWPORT, viewport);
			::glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
			::glGetDoublev(GL_PROJECTION_MATRIX, pjMatrix);
			::gluUnProject(
				(double)point_x, (double)viewport[3]-point_y,
				picked_object[i].min_depth*0.5,
				mvMatrix, pjMatrix,	viewport,
				&ox, &oy, &oz);
		}
		selec_obj.picked_pos.x = ox;
		selec_obj.picked_pos.y = oy;
		selec_obj.picked_pos.z = oz;
	}
	return aSelecObj;
}



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void Com::View::SetProjectionTransform(const Com::View::CCamera& mvp_trans)
{
  //	::glMatrixMode(GL_PROJECTION);
	if( mvp_trans.IsPers() ){	// 透視投影変換
		double fov_y,aspect,clip_near,clip_far;
		mvp_trans.GetPerspective(fov_y,aspect,clip_near,clip_far);
		::gluPerspective(fov_y,aspect,clip_near,clip_far);
	}
	else{	// 正規投影変換
		const double inv_scale = 1.0/mvp_trans.GetScale();
		const double asp = mvp_trans.GetWindowAspect();
		const double h_h = mvp_trans.GetHalfViewHeight()*inv_scale;
		const double h_w = mvp_trans.GetHalfViewHeight()*inv_scale*asp;
		const double depth = 2.0*(h_h+h_w);
		::glOrtho(-h_w,h_w, -h_h, h_h, -depth, depth);
	}
}

void Com::View::SetModelViewTransform(const View::CCamera& mvp_trans)
{
	{	// pan the object
		double x,y,z;
		mvp_trans.GetCenterPosition(x,y,z);
		::glTranslated( x, y, z );
    //    std::cout << x << " " << y << " " << z << std::endl;
	}
  /*	{	// 物体の中心を原点にする
   double x,y,z;
   mvp_trans.GetObjectCenter(x,y,z);
   ::glTranslated( +x, +y, +z );
   std::cout << x << " " << y << " " << z << std::endl;    
   }	*/
	{	// rotate
		double rot[16];
		mvp_trans.RotMatrix44Trans(rot);
		::glMultMatrixd(rot);
	}
  {	// put the origin at the center of object
		double x,y,z;
		mvp_trans.GetObjectCenter(x,y,z);
		::glTranslated( -x, -y, -z );
    //    std::cout << x << " " << y << " " << z << std::endl;        
	}	  
}



bool Com::View::ReadPPM_SetTexture(const std::string& fname, 
                                   unsigned int& texName, 
                                   unsigned int& texWidth, unsigned int& texHeight)
{
  std::cout << "ReadPPM " << std::endl;
  FILE* fp = fopen(fname.c_str(),"r");
  if( fp == NULL ){ 
    std::cout << "Read PPM Fail" << std::endl;
    return false; 
  }
  
  int w, h;
  std::vector<char> aRGB;
  {
    const unsigned int buffSize = 256;
    char buff[buffSize];
    fgets(buff,buffSize,fp);
    fgets(buff,buffSize,fp);
    sscanf(buff,"%d%d",&w,&h);			
    fgets(buff,buffSize,fp);  // read 255
  }
  std::cout << "tex size : " << w << " " << h << std::endl;
  assert( w >= 0 && h >=0 );
  aRGB.resize(w*h*3);
  const unsigned int buffSize = (unsigned int)(4*3*w*1.2);	// ÇøÇÂÇ¡Ç∆ó]ï™ñ⁄Ç…Ç∆Ç¡ÇƒÇ®Ç≠
  char* buff = new char [buffSize];
  unsigned int icnt = 0;
  while (icnt<w*h*3) {
    fgets(buff,buffSize,fp);
    char* pCur = buff;
    char* pNxt;
    for(;;){
      if(      pCur[0] == ' ' ){ assert(0); }
      if(      pCur[1] == ' ' || pCur[1] == '\n' ){ pCur[1]='\0'; pNxt=pCur+2; }
      else if( pCur[2] == ' ' || pCur[2] == '\n' ){ pCur[2]='\0'; pNxt=pCur+3; }
      else if( pCur[3] == ' ' || pCur[3] == '\n' ){ pCur[3]='\0'; pNxt=pCur+4; }
      else{ assert(0); }
      unsigned int val = atoi(pCur);
      unsigned int ih = icnt/(w*3);
      unsigned int iw = icnt-ih*w*3;
      aRGB[(h-ih-1)*w*3+iw] = val;
      icnt++;
      if( pNxt[0] == '\n' || pNxt[0] == '\0') break;
      pCur = pNxt;
    }
  }
  delete[] buff;
  //	this->SetImage(w,h,aRGB);
  
  std::cout << "width height : " << w << " " << h << std::endl;
  
  ////////////////
  
  GLubyte* inputRGB = new GLubyte [w*h*3];
  for(unsigned int i=0;i<w*h*3;i++){ inputRGB[i] = aRGB[i]; }

  texWidth = 1;    
  for(;;){
    if( w <= texWidth ) break;
    texWidth *= 2;
  }    
  ////
  texHeight = 1;    
  for(;;){
    if( w <= texHeight ) break;
    texHeight *= 2;
  }          
  
  GLubyte* scaledRGB;
  if( w == texWidth && h == texHeight ){
    scaledRGB = inputRGB;
  }
  else{
    scaledRGB = new GLubyte [texWidth*texHeight*3];
    gluScaleImage( GL_RGB, w, h, GL_UNSIGNED_BYTE, inputRGB,
                  texWidth, texHeight, GL_UNSIGNED_BYTE, scaledRGB );
    delete [] inputRGB;
  }
  
  glEnable(GL_TEXTURE_2D);
  glGenTextures(1 , &texName);
  glBindTexture(GL_TEXTURE_2D , texName);
  glTexImage2D(GL_TEXTURE_2D , 0 , 3 , texWidth, texHeight,
               0 , GL_RGB , GL_UNSIGNED_BYTE , scaledRGB );
  delete[] scaledRGB;
  
//  std::cout << m_texName << std::endl;
  
  return true;
}


bool WritePPM_ScreenImage(const  std::string& fname)
{
  int viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport);
  void* image = malloc(3 * viewport[2] * viewport[3]);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, image);  
  unsigned int width = viewport[2];
  unsigned int height = viewport[3];
  std::ofstream fout;
  //  //  fname << "out";
  fout.open(fname.c_str(),std::ios::out);
  fout << "P3\n";
  fout << width << " " << height << "\n";
  fout << "255\n";
  //  fout << "255\n";
  //  fout << "255\n";
  char* img = (char*)image;  
  for(unsigned int ih=0;ih<height;ih++){    
    for(unsigned int iw=0;iw<width;iw++){    
      unsigned int i = (height-1-ih)*width+iw;
      int r = (unsigned char)img[i*3+0];
      int g = (unsigned char)img[i*3+1];
      int b = (unsigned char)img[i*3+2];
      fout << r << " " << g << " " << b << "\n";
      //    std::cout << i << " " << r << " "<< g << " "<< b << std::endl;
    }
  }
  fout.close();  
  return true;
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CDrawerCoord::CDrawerCoord(const CCamera& trans, unsigned int win_h)
{
	m_win_h = win_h;
	this->SetTrans(trans, win_h);
}

void CDrawerCoord::SetTrans(const CCamera& trans, int win_h)
{
	this->sutable_rot_mode = 1;
  
	double hh, hw, hd;
	trans.GetOrtho(hw,hh,hd);
	double ox,oy,oz;
	trans.GetCenterPosition(ox,oy,oz);
  
	if( win_h == 0 ){
		x_axis_coord.clear();
		x_axis_name.clear();
		y_axis_coord.clear();
		y_axis_name.clear();
		return ;
	}
  
	if( win_h != -1 ){
		m_win_h = win_h;
	}
	
	m_tex_scale = (double)hh/m_win_h;
  
	{
		coord_len = fabs(-ox+hw);
		coord_len = (coord_len>fabs(-ox-hw)) ? coord_len : fabs(-ox-hw);
		coord_len = (coord_len>fabs(-oy+hh)) ? coord_len : fabs(-oy+hh);
		coord_len = (coord_len>fabs(-oy-hh)) ? coord_len : fabs(-oy-hh);
		coord_len *= 1.5;
	}
  
	double div;
	{
		double len = (hh+hw)*0.5;
		if( len > 1.0e-10 ){
			double scale = log10(len);
			double scale_num = floor(scale);
			div = pow(10.0,scale_num);
			if( len / div > 8.0 ){ div = div * 5.0; }
			else if( len / div > 5.0 ){ div = div*2.0; }
			else if( len / div < 2.0 ){ div = div * 0.5; }
		}
		else{
			x_axis_coord.clear();
			x_axis_name.clear();
			y_axis_coord.clear();
			y_axis_name.clear();
			return ;
		}
	}

	{
		x_axis_coord.clear();
		x_axis_name.clear();
		char tmp_buff[256];
		double x_min = -ox - hw;
		double x_max = -ox + hw;
		double x0 = x_min;
		{
			double dix0 = ( x0 < 0.0 ) ? ceil(x0/div) : floor(x0/div);
			x0 = dix0*div;
		}
		for(;;){
			x_axis_coord.push_back(x0);
			sprintf(tmp_buff,"%f",x0);
			const unsigned int ilen = strlen(tmp_buff);
			for(int i=ilen-1;i>=0;i--){
				if( tmp_buff[i] != '0' ){
					tmp_buff[i+1] = '\0';
					break;
				}
			}
			x_axis_name.push_back(tmp_buff);
			x0 += div;
			if( x0 > x_max ) break;
		}
	}

	{
		y_axis_coord.clear();
		y_axis_name.clear();
		char tmp_buff[256];
		double y_min = -oy - hh;
		double y_max = -oy + hh;
		double y0 = y_min;
		{
			double diy0 = ( y0 < 0.0 ) ? ceil(y0/div) : floor(y0/div);
			y0 = diy0*div;
		}
		for(;;){
			y_axis_coord.push_back(y0);

			sprintf(tmp_buff,"%f",y0);
			const unsigned int ilen = strlen(tmp_buff);
			for(int i=ilen-1;i>=0;i--){
				if( tmp_buff[i] != '0' ){
					tmp_buff[i+1] = '\0';
					break;
				}
			}
			y_axis_name.push_back(tmp_buff);
			y0 += div;
			if( y0 > y_max ) break;
		}
	}
}
void CDrawerCoord::Draw() const
{
  if( !this->is_show ) return;

	::glDisable(GL_DEPTH_TEST);
	::glLineWidth(1);
	::glColor3d(0.5,0.5,0.5);

	for(unsigned int i=0;i<x_axis_coord.size();i++){
		double d = x_axis_coord[i];
		::glPushMatrix();
		::glTranslated(d,0.0,0.0);
		::glScaled(15.0*m_tex_scale,20.0*m_tex_scale,1.0);
		YsDrawUglyFont((char*)(x_axis_name[i].c_str()),0);
		::glPopMatrix();
	}

	for(unsigned int i=0;i<y_axis_coord.size();i++){
		double d = y_axis_coord[i];
		::glPushMatrix();
		::glTranslated(0.0,d,0.0);
		::glScaled(15.0*m_tex_scale,20.0*m_tex_scale,1.0);
		YsDrawUglyFont((char*)(y_axis_name[i].c_str()),0);
		::glPopMatrix();
	}

	::glBegin(GL_LINES);
	::glVertex3d(-coord_len, 0.0,-0.1);
	::glVertex3d( coord_len, 0.0,-0.1);
	::glVertex3d( 0.0,-coord_len,-0.1);
	::glVertex3d( 0.0, coord_len,-0.1);
	::glEnd();
	::glEnable(GL_DEPTH_TEST);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CDrawerRect::CDrawerRect(double pos_x, double pos_y, unsigned int imode){
	end_x = pos_x;
	end_y = pos_y;
	begin_x = pos_x;
	begin_y = pos_y;
	m_imode = imode;
}

void  CDrawerRect::SetPosition(double pos_x, double pos_y){
	end_x = pos_x;
	end_y = pos_y;
}

void CDrawerRect::Draw() const
{
	if( m_imode == 0 ){
		::glMatrixMode(GL_PROJECTION);
		::glPushMatrix();
		::glLoadIdentity();
		::glMatrixMode(GL_MODELVIEW);
		::glPushMatrix();
		::glLoadIdentity();
		{
			::glColor3d(0.0,0.0,0.0);
			::glEnable(GL_LINE_STIPPLE);
			::glLineStipple(1 , 0x3333);
			::glLineWidth(1);
			::glBegin(GL_LINES);
			::glVertex2d(begin_x,begin_y);
			::glVertex2d(begin_x,end_y);
			::glVertex2d(begin_x,end_y);
			::glVertex2d(end_x,end_y);
			::glVertex2d(end_x,end_y);
			::glVertex2d(end_x,begin_y);
			::glVertex2d(end_x,begin_y);
			::glVertex2d(begin_x,begin_y);
			::glEnd();
			::glDisable(GL_LINE_STIPPLE);
		}
		::glPopMatrix();
		::glMatrixMode(GL_PROJECTION);
		::glPopMatrix();
		::glMatrixMode(GL_MODELVIEW);
	}
	else if( m_imode == 1 ){
		::glColor3d(0.0,0.0,0.0);
    ::glLineWidth(1);
    double offset = 0.1;
    ::glBegin(GL_LINES);
    ::glVertex3d(begin_x,begin_y, offset);
    ::glVertex3d(begin_x,end_y,   offset);
    ::glVertex3d(begin_x,end_y,   offset);
    ::glVertex3d(end_x,  end_y,   offset);
    ::glVertex3d(end_x,  end_y,   offset);
    ::glVertex3d(end_x,  begin_y, offset);
    ::glVertex3d(end_x,  begin_y, offset);
    ::glVertex3d(begin_x,begin_y, offset);
		::glEnd();
	}
}


void CDrawerRect::GetCenterSize(double& cent_x, double& cent_y, double& size_x, double& size_y)
{
	cent_x = (begin_x+end_x)*0.5;
	cent_y = (begin_y+end_y)*0.5;
	size_x = fabs(begin_x-end_x);
	size_y = fabs(begin_y-end_y);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


bool CDrawerImageTexture::ReadPPM(const std::string& fname)
{
//	std::cout << "ReadPPM " << std::endl;
	FILE* fp = fopen(fname.c_str(),"r");
	if( fp == NULL ){ return false; }
	
	int w, h;
	std::vector<char> aRGB;
	{
		const unsigned int buffSize = 256;
		char buff[buffSize];
		fgets(buff,buffSize,fp);
		fgets(buff,buffSize,fp);
		sscanf(buff,"%d%d",&w,&h);			
		fgets(buff,buffSize,fp);
	}
	assert( w >= 0 && h >=0 );
	aRGB.resize(w*h*3);
	const unsigned int buffSize = (unsigned int)(4*3*w*1.2);	// ちょっと余分目にとっておく
	char* buff = new char [buffSize];
	for(unsigned int ih=0;ih<(unsigned int)h;ih++){
		fgets(buff,buffSize,fp);
		char* pCur = buff;
		char* pNxt;
		unsigned int icnt = 0;
		for(;;){
			if(      pCur[0] == ' ' ){ assert(0); }
			if(      pCur[1] == ' ' ){ pCur[1]='\0'; pNxt=pCur+2; }
			else if( pCur[2] == ' ' ){ pCur[2]='\0'; pNxt=pCur+3; }
			else if( pCur[3] == ' ' ){ pCur[3]='\0'; pNxt=pCur+4; }
			else{ assert(0); }
			unsigned int val = atoi(pCur);
			aRGB[(h-ih-1)*w*3+icnt] = val;
			icnt++;
			if( pNxt[0] == '\n' || pNxt[0] == '\0') break;
			pCur = pNxt;
		}
		assert( icnt == (unsigned int)w*3 );
	}
	delete[] buff;
	this->SetImage(w,h,aRGB);
	return true;
}

void CDrawerImageTexture::DeleteTexture()
{
    if( m_texName==0 ) return;
    glDeleteTextures(1 , &m_texName);
	m_texName = 0;
}

// ２のべき乗に合わせるために，対数の少数点切り捨てたものを返す
static unsigned int TextureCuttOffSize(unsigned int v)
{	
	if(      v <   2 ){ return    0; }
	else if( v <   4 ){ return    2; }
	else if( v <   8 ){ return    4; }
	else if( v <  16 ){ return    8; }
	else if( v <  32 ){ return   16; }
	else if( v <  64 ){ return   32; }
	else if( v < 128 ){ return   64; }
	else if( v < 256 ){ return  128; }
	else if( v < 512 ){ return  256; }
	else if( v <1024 ){ return  512; }
	else if( v <2048 ){ return 1024; }
	else if( v <4096 ){ return 2048; }
	else if( v <8192 ){ return 4096; }
	assert(0);
	return 0;
}

bool CDrawerImageTexture::SetImage(unsigned int w, unsigned int h, const std::vector<char>& aRGB)
{
//	std::cout << "SetImageTexture" << std::endl;
	this->DeleteTexture();
	assert( aRGB.size() == w*h*3 );

	GLubyte* inputRGB = new GLubyte [w*h*3];
	for(unsigned int i=0;i<w*h*3;i++){ inputRGB[i] = aRGB[i]; }

	m_texWidth = TextureCuttOffSize(w);
	m_texHight = TextureCuttOffSize(h);
//	std::cout << m_texWidth << " " << m_texHight << std::endl;

	GLubyte* scaledRGB;
	if( w == m_texWidth && h == m_texHight ){
		scaledRGB = inputRGB;
	}
	else{
		scaledRGB = new GLubyte [m_texWidth*m_texHight*3];
		gluScaleImage( GL_RGB, w, h, GL_UNSIGNED_BYTE, inputRGB,
			m_texWidth, m_texHight, GL_UNSIGNED_BYTE, scaledRGB );
		delete [] inputRGB;
	}

	glEnable(GL_TEXTURE_2D);
	glGenTextures(1 , &m_texName);
	glBindTexture(GL_TEXTURE_2D , m_texName);
	glTexImage2D(
		GL_TEXTURE_2D , 0 , 3 , m_texWidth, m_texHight,
		0 , GL_RGB , GL_UNSIGNED_BYTE , scaledRGB );
	delete[] scaledRGB;

	x_min = 0;
	x_max = (double)w/h;

	return true;
}


void CDrawerImageTexture::Draw() const
{
    if( m_texName == 0 ) return;
	// テクスチャマップの方法を設定
	::glColor3d(1.0, 1.0, 1.0);
	glEnable(GL_TEXTURE_2D);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glShadeModel(GL_FLAT);
	glBindTexture(GL_TEXTURE_2D , m_texName);
    glBegin(GL_QUADS);
        glTexCoord2d(0 , 0); glVertex3d(x_min, y_min, -0.5);
        glTexCoord2d(1 , 0); glVertex3d(x_max, y_min, -0.5);
        glTexCoord2d(1 , 1); glVertex3d(x_max, y_max, -0.5);
        glTexCoord2d(0 , 1); glVertex3d(x_min, y_max, -0.5);
	glEnd();
}
