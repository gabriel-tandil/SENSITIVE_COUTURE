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
@brief GLUTの便利な関数，クラス群
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_GL_UTILITY_H)
#define DRAWER_GL_UTILITY_H

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include <vector>
#include <assert.h>

#include "delfem/drawer.h"
#include "delfem/vector3d.h"

namespace Com{
namespace View{

//! 選択オブジェクト
struct SSelectedObject{
  unsigned int name_depth;	//!< 名前の深さ
  int name[4];	//!< 名前を格納している配列(4つまでしかないのは問題じゃね？)
  Com::CVector3D picked_pos;	//!< ピックされた点の３次元的な位置
};
  
  
class CCamera;
/*!
 @brief Projection変換をする
 @remark glLoadIdentity()はこの関数に含まれないので，この関数を呼ぶ前に行うこと(ピックで用いるため)
 */
void SetProjectionTransform(const CCamera& mvp_trans);	
/*!
 @brief ModelView変換をする
 @remark glLoadIdentity()はこの関数に含まれないので，この関数を呼ぶ前に行うこと(ピックで用いるため)
 */
void SetModelViewTransform(const CCamera& mvp_trans);
  
//! ピック前処理
void PickPre(unsigned int size_buffer, unsigned int* select_buffer,
             unsigned int point_x, unsigned int point_y,
             unsigned int delX, unsigned int delY,
             const View::CCamera& mvp_trans);
  
//! ピック後処理		
std::vector<SSelectedObject> PickPost(unsigned int* const select_buffer,
                                      unsigned int point_x, unsigned int point_y,
                                      const View::CCamera& mvp_trans);
  
bool ReadPPM_SetTexture(const std::string& fname, 
                        unsigned int& texName, 
                        unsigned int& texWidth, unsigned int& texHeight);
  
bool WritePPM_ScreenImage(const  std::string& fname);
  
//! Draw coordinate
class CDrawerCoord : public CDrawer{
public:
  CDrawerCoord(){}
  CDrawerCoord(const CCamera& trans, unsigned int win_h );
  virtual void Draw() const;
  
  // virutal functions which do nothing
  virtual void DrawSelection(unsigned int idraw) const {}
  virtual Com::CBoundingBox3D GetBoundingBox(double rot[]) const { return Com::CBoundingBox3D(); }
  virtual void AddSelected(const int selec_flag[]){}
  virtual void ClearSelected(){}
  
  // non-virtual functions
  void SetTrans(const CCamera& trans, int win_h = -1);
  void SetIsShown(bool is_show){ this->is_show = is_show; }
  bool GetIsShown() const { return is_show; }
private:
  bool is_show;
  
  std::vector<double> x_axis_coord;
  std::vector<std::string> x_axis_name;
  
  std::vector<double> y_axis_coord;
  std::vector<std::string> y_axis_name;
  
  unsigned int m_win_h;
  double m_tex_scale;
  double coord_len;
};
  
  
//! Draweing Rectangular Box for selection or specifying region
class CDrawerRect : public CDrawer
{
public:
  CDrawerRect(){
    begin_x = 0;    begin_y = 0;
    this->m_imode = 0;
  }
  CDrawerRect(double x, double y, unsigned int imode = 0);
  virtual void Draw() const;
  void SetInitialPositionMode(double x, double y, unsigned int imode){
    begin_x = x;    begin_y = y;
    this->m_imode = imode;
  }
  void SetPosition(double x, double y);
  void GetCenterSize(double& cent_x, double& cent_y, double& size_x, double& size_y);
  void GetPosition(double& x0, double& y0, double& x1, double& y1){
    x0 = begin_x;   y0 = begin_y;
    x1 = end_x;   y1 = end_y;
  }
  
  // vertual function
  virtual void DrawSelection(unsigned int idraw) const {}
  virtual Com::CBoundingBox3D GetBoundingBox(double rot[]) const { return Com::CBoundingBox3D(); }
  virtual void AddSelected(const int selec_flag[]){}
  virtual void ClearSelected(){}
  private :
  double begin_x,begin_y;
  double end_x, end_y;
  unsigned int m_imode;
};  
  
  
//! Draw texture in the background
class CDrawerImageTexture : public CDrawer
{
public:
  CDrawerImageTexture(){
    m_texName = 0;
    m_texWidth = 0;
    m_texHight = 0;
    x_min=0;	x_max=1;
    y_min=0;	y_max=1;
  }
  bool IsTexture(){ return m_texName!=0;}
  bool ReadPPM(const std::string& fname);
  void DeleteTexture();
  bool SetImage(unsigned int w, unsigned int h, const std::vector<char>& aRGB);
  virtual void Draw() const;
  
  // 以下のvirtual関数は実装されない
  virtual void DrawSelection(unsigned int idraw) const {}
  virtual Com::CBoundingBox3D GetBoundingBox(double rot[]) const { 
    return Com::CBoundingBox3D(x_min,x_max, y_min,y_max, -1,1); 
  }
  virtual void AddSelected(const int selec_flag[]){}
  virtual void ClearSelected(){}
  private :
  unsigned int m_texName;	// if 0 then no texture
  unsigned int m_texWidth;	// should be power of 2
  unsigned int m_texHight;	// should be power of 2
  double x_min,x_max,  y_min,y_max;
};  
  
} // end namespace View
} // end namespace Com


#endif
