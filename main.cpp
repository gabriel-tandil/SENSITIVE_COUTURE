/*
 *  main.cpp
 *  sensitive couture
 *
 *  Created by Nobuyuki Umetani on 9/14/10.
 *  Copyright 2010 The University of Tokyo and Columbia University. All rights reserved.
 *
 */



//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#pragma warning( disable : 4786 )

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "delfem/camera.h"

#include "analysis2d_cloth_static.h"
#include "designer2d_cloth.h"

double mov_begin_x, mov_begin_y;	// ÉJÅ[É\Éãà íuÇÃê≥ãKç¿ïW
int imodifier;

bool is_lighting = true;
bool is_animation = true;
CDesigner2D_Cloth gui_listner;
CAnalysis2D_Cloth_Static* pAnalysis = 0;

Com::View::CCamera camera_l;
Com::View::CCamera camera_r;
Com::View::CDrawerCoord drawer_coord;
unsigned int iclass_type = 0;
int iwin_des, iwin_sim;
int imenu_right_click;

// 0:drag
// 1:dart cut
// 2:curve_edit
// 3:smooth
// 4:slider 0
// 5:dart_cut not sew
unsigned int imode_ope = 0;
int islider_active = -1;
std::string slider_name;
bool is_active_slider = false;
double pick_org[3];
double pick_sns[3];
//std::vector<double> aXY_Stroke;
Com::CVector2D pos_pick, pos_relese;
bool is_rubber_band = false;
bool is_tex_mouse_pos = false;
unsigned int m_texWidth, m_texHeight;
unsigned int m_texName = 0;

bool is_display_rotation = false;
double display_rotation_theta = 0;

void RenderBitmapString(float x, float y, void *font,char *string)
{
  char *c;
  ::glRasterPos2f(x, y);
  for (c=string; *c != '\0'; c++) {
	  ::glutBitmapCharacter(font, *c);
  }
}

void ShowTextDesign()
{
  const bool is_lighting = glIsEnabled(GL_LIGHTING);  
	int* pFont=(int*)GLUT_BITMAP_8_BY_13;

	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];

	::glMatrixMode(GL_PROJECTION);
	::glPushMatrix();
	::glLoadIdentity();
	::gluOrtho2D(0, win_w, 0, win_h);
	::glMatrixMode(GL_MODELVIEW);
	::glPushMatrix();
	::glLoadIdentity();
	::glScalef(1, -1, 1);
	::glTranslatef(0, -win_h, 0);
  ::glDisable(GL_LIGHTING);
//	::glDisable(GL_DEPTH_TEST);
  // 0:drag
  // 1:dart cut
  // 2:curve_edit
  // 3:smooth
  // 4:slider 0  
	char s_tmp[256];  
  {
    strcpy(s_tmp,"Press \' \' button to try different 3D model");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,15, pFont, s_tmp);            
  }
  {
    strcpy(s_tmp,"Click right button of mouse to change tool");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,45, pFont, s_tmp);            
  } 
  { 
    if( pAnalysis != 0 ){
      if( pAnalysis->GetMode() == CLOTH_INITIAL_LOCATION ){
        strcpy(s_tmp,"Press \'b\' button to start simulation"); 
      }
      else{
        strcpy(s_tmp,"Press \'b\' button to initialize");         
      }
      ::glColor3d(0.0, 0.0, 1.0);
      RenderBitmapString(10,30, pFont, s_tmp);      
    }
  }  
  {
    strcpy(s_tmp,"Pan : Shift+Left drag,  Zoom : Page Up/Down");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,60, pFont, s_tmp);            
  } 
  if(      imode_ope == 0 ){ strcpy(s_tmp,"Tool : Drag"); }
  else if( imode_ope == 1 ){ strcpy(s_tmp, "Tool : Dart Cutter"); }
  else if( imode_ope == 2 ){ strcpy(s_tmp, "Tool : Curve Edit"); }  
  else if( imode_ope == 3 ){ strcpy(s_tmp, "Tool : Smooth"); }    
  else if( imode_ope == 4 ){ strcpy(s_tmp, "Tool : Slider"); }      
  else if( imode_ope == 5 ){ strcpy(s_tmp, "Tool : Hole"); }        
	::glColor3d(1.0, 0.0, 0.0);  
	RenderBitmapString(10,80, pFont, s_tmp);
  if( imode_ope == 4 ){
    s_tmp[32];
    strcpy(s_tmp,slider_name.c_str());
    ::glColor3d(1,0,0);
    RenderBitmapString(10,win_h*0.95-13, pFont, s_tmp);    
  }
  {
    strcpy(s_tmp,"Please wait a second at the beginning of mouse drag.");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,win_h*0.95-26, pFont, s_tmp);            
  }   
	if( is_lighting ){ ::glEnable(GL_LIGHTING); }    
	::glEnable(GL_DEPTH_TEST);
	::glPopMatrix();
	::glMatrixMode(GL_PROJECTION);
	::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
}

void ShowTextSimulation()
{
	int* pFont=(int*)GLUT_BITMAP_8_BY_13;
	static char s_fps[32];
	{
		static int frame, timebase;
		int time;
		frame++;
		time=glutGet(GLUT_ELAPSED_TIME);
		if (time - timebase > 500) {
      //	if (time - timebase > 5000) {
			sprintf(s_fps,"FPS:%4.2f",frame*1000.0/(time-timebase));
			timebase = time;
			frame = 0;
		}
	}
  
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
  
	::glMatrixMode(GL_PROJECTION);
	::glPushMatrix();
	::glLoadIdentity();
	::gluOrtho2D(0, win_w, 0, win_h);
	::glMatrixMode(GL_MODELVIEW);
	::glPushMatrix();
	::glLoadIdentity();
	::glScalef(1, -1, 1);
	::glTranslatef(0, -win_h, 0);
	const bool is_lighting = glIsEnabled(GL_LIGHTING);
  ::glDisable(GL_LIGHTING);
  char s_tmp[256];
  {
    strcpy(s_tmp,"Pan : Shift+Left drag,  Zoom : Page Up/Down");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,15, pFont, s_tmp);            
  }
  {
    strcpy(s_tmp,"Rotation : Ctrl+Left drag");
    ::glColor3d(0.0, 0.0, 1.0);
    RenderBitmapString(10,30, pFont, s_tmp);            
  }   
  ::glColor3d(1.0, 0.0, 0.0);
	RenderBitmapString(10,45, pFont, s_fps);
  if( is_lighting ){ ::glEnable(GL_LIGHTING); }  
	::glEnable(GL_DEPTH_TEST);
	::glPopMatrix();
	::glMatrixMode(GL_PROJECTION);
	::glPopMatrix();
	::glMatrixMode(GL_MODELVIEW);
}


void myGlutMotion( int x, int y )
{
	int viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);	
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	////////////////
	if( imodifier == GLUT_ACTIVE_SHIFT ){ // pan
		if( glutGetWindow() == iwin_sim ){	
			camera_r.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 		
			drawer_coord.SetTrans(camera_r);
		}
		else{
			camera_l.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
			drawer_coord.SetTrans(camera_l);
		}
	}
	else if( imodifier == GLUT_ACTIVE_CTRL ){ //  rotation
		if( glutGetWindow() == iwin_sim ){
			camera_r.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 		
			drawer_coord.SetTrans(camera_r);
		}
	}
	else{ // pick and move
    if( glutGetWindow() == iwin_des ){
      if( imode_ope == 0 ){
        Com::CVector3D oloc0, oloc1;
        oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);
        if( is_tex_mouse_pos ){ 
          gui_listner.SetTextureCenter_FaceCAD(-0.5+oloc0.x,-0.5+oloc0.y); 
          pAnalysis->SetTextureCenter(-0.5+oloc0.x,-0.5+oloc0.y); 
        }
        oloc1 = camera_l.ProjectionOnPlane(mov_end_x,mov_end_y);
        Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
        double tmp_x, tmp_y;
        gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);
        if( id_cad_part != 0 ){
          bool res = gui_listner.Cad_Move( itype_cad_part, id_cad_part, Com::CVector2D(oloc0.x,oloc0.y), Com::CVector2D(oloc1.x,oloc1.y) ); 		        
          if( res && pAnalysis->GetMode() == SENSITIVITY_DONE ){           
            pAnalysis->GuessSolutionMouseMove(oloc0.x,oloc0.y, oloc1.x,oloc1.y);
          }
        }
      }
      else if( imode_ope == 1 || imode_ope == 5 ){
        const Com::CVector3D& oloc1 = camera_l.ProjectionOnPlane(mov_end_x,mov_end_y);      
        pos_relese = Com::CVector2D(oloc1.x,oloc1.y);         
      }
      else if( imode_ope == 2 ){
        Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
        double tmp_x, tmp_y;
        gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);          
        if( itype_cad_part == Cad::EDGE && gui_listner.GetCad().IsElemID(itype_cad_part,id_cad_part) ){
          const Com::CVector3D& oloc = camera_l.ProjectionOnPlane(mov_end_x,mov_end_y);  
          if(      gui_listner.GetCad().GetEdgeCurveType(id_cad_part) == 1 ){
            bool res = gui_listner.Cad_DragArc(id_cad_part,Com::CVector2D(oloc.x,oloc.y));
          }
          else if( gui_listner.GetCad().GetEdgeCurveType(id_cad_part) == 2 ){
            bool res = gui_listner.Cad_DragPolyline(id_cad_part,Com::CVector2D(oloc.x,oloc.y));
          }
          if( pAnalysis->GetMode() != SIMULATION_STATIC && 
              pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){ pAnalysis->PerformStaticSolver(); }
        }              
      }
      else if( imode_ope == 3 ){
        Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
        double tmp_x, tmp_y;
        gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);          
        if( itype_cad_part == Cad::EDGE && gui_listner.GetCad().IsElemID(itype_cad_part,id_cad_part) ){
          const Com::CVector3D& oloc = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y); 
          gui_listner.Cad_SmoothingPolylineEdge(id_cad_part, 3,
                                                Com::CVector2D(oloc.x,oloc.y), 0.1); 
          if( pAnalysis->GetMode() != SIMULATION_STATIC && 
              pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){ pAnalysis->PerformStaticSolver(); }          
        }      
      }        
      else if( imode_ope == 4 ){
        if( islider_active >= 0 && is_active_slider ){
          double min, max;
          const double val0 = gui_listner.GetValueSlider(islider_active,min,max);
          const double val1 = min+(mov_end_x+1)*0.5*(max-min);
          gui_listner.SetValueSlider(islider_active,val1);        
          const double val2 = gui_listner.GetValueSlider(islider_active,min,max);
          if( pAnalysis->GetMode() == SENSITIVITY_DONE ){           
            pAnalysis->GuessSolutionSliderMove(val0,val2);
          }          
        }
      }
      if( is_tex_mouse_pos ){
        Com::CVector3D oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);
        pAnalysis->SetTextureCenter(-0.5+oloc0.x,-0.5+oloc0.y);
        gui_listner.SetTextureCenter_FaceCAD(-0.5+oloc0.x,-0.5+oloc0.y);         
      }
    }
    else if( glutGetWindow() == iwin_sim ){      
      if( imode_ope == 0 ){
        if( pAnalysis->GetMode() == CLOTH_INITIAL_LOCATION ){
          Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
          double tmp_x, tmp_y;
          gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);
          if( itype_cad_part == Cad::LOOP && gui_listner.GetCad().IsElemID(itype_cad_part,id_cad_part) ){
            double dir[3];
            {
              const double hvh = camera_r.GetHalfViewHeight()/camera_r.GetScale();
              const double asp = camera_r.GetWindowAspect();
              double d0[3] = { (mov_end_x-mov_begin_x)*hvh*asp, (mov_end_y-mov_begin_y)*hvh, 0 };     
              double rot[9];
              camera_r.RotMatrix33(rot);
              dir[0] = rot[0]*d0[0]+rot[3]*d0[1]+rot[6]*d0[2];
              dir[1] = rot[1]*d0[0]+rot[4]*d0[1]+rot[7]*d0[2];
              dir[2] = rot[2]*d0[0]+rot[5]*d0[1]+rot[8]*d0[2];
            } 
            pAnalysis->MoveClothLoopInitialPosition(id_cad_part,dir);          
          }          
        }
      }
      if( imode_ope == 4 ){
        if( is_active_slider ){
          double rot[9];
          camera_r.RotMatrix33(rot);
          const double sns1[3] = {
            rot[0]*pick_sns[0] + rot[1]*pick_sns[1] + rot[2]*pick_sns[2],
            rot[3]*pick_sns[0] + rot[4]*pick_sns[1] + rot[5]*pick_sns[2],
            rot[6]*pick_sns[0] + rot[7]*pick_sns[1] + rot[8]*pick_sns[2] };
          const double hvh = camera_r.GetHalfViewHeight()/camera_r.GetScale();
//          const double hvh = camera_r.GetHalfViewHeight();          
          const double asp = camera_r.GetWindowAspect();          
          const double mus1[2] = { (mov_end_x-mov_begin_x)*hvh*asp, (mov_end_y-mov_begin_y)*hvh };
          double dot = (mus1[0]*sns1[0] + mus1[1]*sns1[1])/(sns1[0]*sns1[0] + sns1[1]*sns1[1]);
          double min, max;
          const double val0 = gui_listner.GetValueSlider(islider_active,min,max);
          gui_listner.SetValueSlider(islider_active,val0+dot);             
          const double val2 = gui_listner.GetValueSlider(islider_active,min,max);
          if( pAnalysis->GetMode() == SENSITIVITY_DONE ){           
            pAnalysis->GuessSolutionSliderMove(val0,val2);
          }                    
        }
      }
    }
  }
	////////////////
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}


void myGlutPassive( int x, int y )
{
  if( glutGetWindow() == iwin_sim ){ return; } 
	int viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);	
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
  ////
  const unsigned int size_buffer = 128;
  unsigned int select_buffer[size_buffer];
  std::vector<Com::View::SSelectedObject> aSelecObj;
  Com::View::PickPre(size_buffer, select_buffer, x,y,  12,12,  camera_l);
  gui_listner.DrawSelection();
  aSelecObj = Com::View::PickPost(select_buffer, x,y, camera_l);
  gui_listner.SetSelection(aSelecObj);
  ////
  Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
  double tmp_x, tmp_y;
  gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);  
  if( is_tex_mouse_pos ){
    Com::CVector3D oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);    
    pAnalysis->SetTextureCenter(-0.5+oloc0.x,-0.5+oloc0.y);
    gui_listner.SetTextureCenter_FaceCAD(-0.5+oloc0.x,-0.5+oloc0.y);     
  }
//  if( itype_cad_part == Cad::EDGE || itype_cad_part == Cad::VERTEX ){
    pAnalysis->SetHilight(itype_cad_part,id_cad_part);
//  }
  ////////////////
  mov_begin_x = mov_end_x;
  mov_begin_y = mov_end_y;
  ::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y)
{
  /////
	imodifier = glutGetModifiers();
	int viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/(double)win_w;
	mov_begin_y = (win_h-2.0*y)/(double)win_h;
  ////
  if( glutGetWindow() == iwin_des ){  
    if( imode_ope == 1 ){ // cut line
      if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
        is_rubber_band = true;
        const Com::CVector3D& oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);    
        pos_pick = Com::CVector2D(oloc0.x,oloc0.y);
        pos_relese = pos_pick;
      }
      if( button == GLUT_LEFT_BUTTON && state == GLUT_UP ){
        if( is_rubber_band ){
          gui_listner.Cad_AddCutLine(pos_pick,pos_relese);
          is_rubber_band = false;
        }
      }      
      return;
    }
    if( imode_ope == 5 ){ // cut line
      if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ){
        const Com::CVector3D& oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);    
        pos_pick = Com::CVector2D(oloc0.x,oloc0.y);
        pos_relese = pos_pick;
      }
      if( button == GLUT_LEFT_BUTTON && state == GLUT_UP ){
        if( is_rubber_band ){
          gui_listner.Cad_AddCutLine(pos_pick,pos_relese);
          is_rubber_band = false;
        }
      }      
      return;
    }    
    if( imode_ope == 2 ){ // polygon pull
      if( state == GLUT_DOWN ){
        const unsigned int size_buffer = 128;
        unsigned int select_buffer[size_buffer];
        std::vector<Com::View::SSelectedObject> aSelecObj;
        Com::View::PickPre(size_buffer, select_buffer, x,y,  8,8,  camera_l);
        gui_listner.DrawSelection();
        aSelecObj = Com::View::PickPost(select_buffer, x,y, camera_l);
        gui_listner.SetSelection(aSelecObj);
        /////
        Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
        double tmp_x, tmp_y;
        gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);    
        if( itype_cad_part == Cad::EDGE && gui_listner.GetCad().IsElemID(itype_cad_part,id_cad_part) ){
          const Com::CVector3D& oloc = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);  
          if( gui_listner.GetCad().GetEdgeCurveType(id_cad_part) == 2 ){
            bool res = gui_listner.Cad_PreCompDragPolyline(id_cad_part,Com::CVector2D(oloc.x,oloc.y));
          }      
          return; 
        }            
      }
    }
    if( imode_ope == 4 ){
      if( state == GLUT_DOWN ){
        if( mov_begin_y < -0.9 ){
          is_active_slider = true;
          double min, max;
          const double val = gui_listner.GetValueSlider(islider_active,min,max);
          const double ratio = (val-min)/(max-min)*2-1;              
          if( fabs(ratio-mov_begin_x) < 0.05 ){
            gui_listner.Msh_PrecompSlider(islider_active);   
            if( pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){    
              std::vector<double> har;
              gui_listner.Msh_GetHarmonicPrecomp(har);  // Get mesh harmonic function
              if( har.size() == 0 ) return;      
              pAnalysis->SetSensitivity_Slider(islider_active,har,val);
            }            
          }
        }
        else{
          is_active_slider = false;
        }
      }
      else if( state ==GLUT_UP ){ 
        is_active_slider = false;
        if( pAnalysis->GetMode() == SENSITIVITY_DONE ){
          pAnalysis->PerformStaticSolver();
          return;
        }        
      }
    }
    if( imode_ope == 0 ){
      if( state == GLUT_UP ){
        if( pAnalysis->GetMode() == SENSITIVITY_DONE ){
          pAnalysis->PerformStaticSolver();
          {
            double max_aspect = 0;
            bool is_inverted = false;
            pAnalysis->MeshQualityInfo(max_aspect,is_inverted);
            std::cout << "ratio " << max_aspect << " " << is_inverted << std::endl;
            if( max_aspect > 10 || is_inverted ){
              gui_listner.Solve_fromCad_InterpValue();        
            }
          }
          return;
        }
      }  
      if( state == GLUT_DOWN ){ 
        if( !pAnalysis->IsBlowUp() ){ gui_listner.SaveTimeStamp(); }        
        /////
        const unsigned int size_buffer = 128;
        unsigned int select_buffer[size_buffer];
        std::vector<Com::View::SSelectedObject> aSelecObj;
        Com::View::PickPre(size_buffer, select_buffer, x,y,  8,8,  camera_l);
        gui_listner.DrawSelection();
        aSelecObj = Com::View::PickPost(select_buffer, x,y, camera_l);
        gui_listner.SetSelection(aSelecObj);
        ////////////////
        Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
        double tmp_x, tmp_y;
        gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);
        if( itype_cad_part == Cad::NOT_SET ) return;
        gui_listner.Msh_PrecompDrag(itype_cad_part, id_cad_part);    
        if( pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION && !gui_listner.IsntMeshDeformSensitivity() ){
          std::vector<double> har;
          gui_listner.Msh_GetHarmonicPrecomp(har);  // Get mesh harmonic function
          if( har.size() == 0 ) return;      
          const Com::CVector3D& oloc0 = camera_l.ProjectionOnPlane(mov_begin_x,mov_begin_y);    
          pAnalysis->SetSensitivity(itype_cad_part,id_cad_part,har,oloc0.x,oloc0.y);
        }            
      }
    }        
  }
  if( glutGetWindow() == iwin_sim ){
    if( imode_ope == 0 ){      
      if( state == GLUT_DOWN ){      
        double trans0[3]; camera_r.GetObjectCenter(trans0[0],trans0[1],trans0[2]);
        trans0[0] = -trans0[0]; trans0[1] = -trans0[1]; trans0[2] = -trans0[2];
        double rot[9];    camera_r.RotMatrix33(rot);
        double trans1[3]; camera_r.GetCenterPosition(trans1[0],trans1[1],trans1[2]);
        const double hvh = camera_r.GetHalfViewHeight()/camera_r.GetScale();
        const double asp = camera_r.GetWindowAspect();
        unsigned int no[3];
        double r[3];
        unsigned int id_l;
        bool res = pAnalysis->Pick(hvh*asp*mov_begin_x,hvh*mov_begin_y,trans0,rot,trans1,
                                   no,r,id_l);//, dir,org); 
        if( !res ){ 
          gui_listner.HilightCadTypeID(Cad::NOT_SET,0);
          gui_listner.Cad_SetPicked(Cad::NOT_SET,0,0,0);          
          return; 
        }
        gui_listner.HilightCadTypeID(Cad::LOOP,id_l);      
        gui_listner.Cad_SetPicked(Cad::LOOP,id_l,0,0);
      }
    }
    if( imode_ope == 4 ){
      if( state == GLUT_DOWN ){      
        double trans0[3]; camera_r.GetObjectCenter(trans0[0],trans0[1],trans0[2]);
        trans0[0] = -trans0[0]; trans0[1] = -trans0[1]; trans0[2] = -trans0[2];
        double rot[9];    camera_r.RotMatrix33(rot);
        double trans1[3]; camera_r.GetCenterPosition(trans1[0],trans1[1],trans1[2]);
        const double hvh = camera_r.GetHalfViewHeight()/camera_r.GetScale();
        const double asp = camera_r.GetWindowAspect();
        unsigned int no[3];
        double r[3];
        unsigned int id_l;
        bool res = pAnalysis->Pick(hvh*asp*mov_begin_x,hvh*mov_begin_y,trans0,rot,trans1,
                                   no,r,id_l);//, dir,org); 
        if( !res ){ 
          std::cout << "pick none" << std::endl;
          is_active_slider = false;
          return; 
        }
        is_active_slider = true;
        gui_listner.Msh_PrecompSlider(islider_active);   
        if( pAnalysis->GetMode() != CLOTH_INITIAL_LOCATION ){    
          std::vector<double> har;
          gui_listner.Msh_GetHarmonicPrecomp(har);  // Get mesh harmonic function
          if( har.size() == 0 ) return;      
          double min, max;
          const double val = gui_listner.GetValueSlider(islider_active,min,max);
          pAnalysis->SetSensitivity_Slider(islider_active,har,val);
          pAnalysis->GetSensitivityElem_Slider(no,r, pick_org,pick_sns);
          //      std::cout << "org : " << org[0] << " " << org[1] << " " << org[2] << std::endl;
          //      std::cout << "sns : " << sns[0] << " " << sns[1] << " " << sns[2] << std::endl;                
        }
      }    
      else if( state == GLUT_UP ){      
        if( pAnalysis->GetMode() == SENSITIVITY_DONE ){
          pAnalysis->PerformStaticSolver();
          return;
        }
      }
    }      
  }
}

void myGlutIdle(){	
	if( is_animation ){
		gui_listner.FollowMshToCad_ifNeeded();
		gui_listner.Solve_ifNeeded();
    if( pAnalysis->IsBlowUp() ){
      std::cout << "BlowUp" << std::endl;
      gui_listner.LoadTimeStamp();
    }
	}
	// redraw left window
	::glutSetWindow(iwin_des);
	::glutPostRedisplay();
  if( pAnalysis->GetMode() != SIMULATION_STATIC ){
    ::glutSetCursor(GLUT_CURSOR_INHERIT);
  }
	
	// redraw right window
	::glutSetWindow(iwin_sim);
	::glutPostRedisplay();
  if( pAnalysis->GetMode() != SIMULATION_STATIC ){
    ::glutSetCursor(GLUT_CURSOR_INHERIT);
  }  
}
 
void myGlutResize(int w, int h)
{		
	if (glutGetWindow() == iwin_des) {	
    ::glViewport(0,0,w,h);    
		camera_l.SetWindowAspect((double)w/h);
		drawer_coord.SetTrans(camera_l,h);
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera_l);
	}
	else{
    ::glViewport(0,0,w,h);        
		camera_r.SetWindowAspect((double)w/h);
		drawer_coord.SetTrans(camera_r,h);
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera_r);					
	}
	::glutPostRedisplay();
}


void myGlutMenu_Des(int value)
{
  Cad::CAD_ELEM_TYPE itype_cad_part; unsigned int id_cad_part;
  double tmp_x, tmp_y;
  gui_listner.Cad_GetPicked(itype_cad_part, id_cad_part, tmp_x,tmp_y);
  
  switch(value) {
    case 5:
      imode_ope = 0;  
      std::cout << "Mode Drag" << std::endl;
      break;      
    case 1:
      imode_ope = 2;  
      std::cout << "Edit Curve" << std::endl;
      break;
    case 6:
      imode_ope = 3;  
      std::cout << "Smooth Polyline" << std::endl;
      break;      
    case 7:
      imode_ope = 1;  
      std::cout << "DartCutter" << std::endl;
      break;            
    case 8:
      imode_ope = 5;  
      std::cout << "DartCutter" << std::endl;
      break;                  
    case 2:	// ChangeToLine
      if( itype_cad_part != Cad::EDGE ) break;
      std::cout << "ChangeToLine" << std::endl;
      gui_listner.Cad_SetCurveType(id_cad_part,0);
      break;
    case 3:	// ChangeToArc
      if( itype_cad_part != Cad::EDGE ) break;
      std::cout << "ChangeToArc" << std::endl;
			gui_listner.Cad_SetCurveType(id_cad_part,1);
      break;
    case 4:	// ChangeToArc
      if( itype_cad_part != Cad::EDGE ) break;
      std::cout << "ChangeToPolyline" << std::endl;
			gui_listner.Cad_SetCurveType(id_cad_part,2);      
      break;      
  }
  if( value >= 8 ){
    imode_ope = 4;
    islider_active = value-8;
    is_active_slider = false;
    slider_name = gui_listner.GetNameSlider(islider_active);
  }
}

void SetNewProblem()
{  
  const unsigned int indprob[4] = {6,9,11,10};
	static int iprob =  0;
	unsigned int nprob = 4;  
	if( pAnalysis != 0 ){ 
    ::glutSetWindow(iwin_des);
    ::glutDetachMenu(GLUT_RIGHT_BUTTON);  
    ::glutDestroyMenu(imenu_right_click);
//    unsigned int nsl = gui_listner.GetNumberOfSlider();
//    for(unsigned int isl=0;isl<nsl;isl++){ glutDestroyMenu(isl+8); }
    delete pAnalysis; 
    pAnalysis = 0; 
  }
  pAnalysis = new CAnalysis2D_Cloth_Static();  
  gui_listner.SetAnalysisInitialize(pAnalysis,indprob[iprob]);  
  if( is_tex_mouse_pos ){ 
    pAnalysis->SetColor_FaceFEM(0.8,0.8,0.8); 
    gui_listner.SetColor_CadFace(0.8, 0.8, 0.8);        
    pAnalysis->SetTextureScale_FaceFEM(1);    
    gui_listner.SetTextureScale_CadFace(1);    
  }
  else{
    pAnalysis->SetColor_FaceFEM(1.0,1.0,1.0); 
    gui_listner.SetColor_CadFace(1, 1, 1);    
    //    pAnalysis->SetTextureScale_FaceFEM(1);
    //    gui_listner.SetTextureScale_CadFace(1);
    pAnalysis->SetTextureScale_FaceFEM(5);
    gui_listner.SetTextureScale_CadFace(5);    
  }  
	////////////////
	{
		camera_r.SetRotationMode(Com::View::ROT_2DH);
		double rot[9];	camera_r.RotMatrix33(rot);
		Com::CBoundingBox3D bb = pAnalysis->GetBoundingBox(rot);
    bb.z_max*=1.5;
		camera_r.SetObjectBoundingBox(bb);
		camera_r.Fit();	
    camera_r.SetScale(1.7);    
		////
		::glutSetWindow(iwin_sim);
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera_r);
		::glutPostRedisplay();
//		drawer_coord.SetTrans(camera_r);
	}
	{
    camera_l.SetRotationMode(Com::View::ROT_2D);
		double rot[9];	camera_l.RotMatrix33(rot);
		const Com::CBoundingBox3D& bb = gui_listner.GetBoundingBox(rot);
		camera_l.SetObjectBoundingBox(bb);
		camera_l.Fit();	
    camera_l.SetScale(1.2);        
//   camera_l.SetScale(1.0);            
    ////
		::glutSetWindow(iwin_des);
		::glMatrixMode(GL_PROJECTION);
		::glLoadIdentity();
		Com::View::SetProjectionTransform(camera_l);
		::glutPostRedisplay();
		drawer_coord.SetTrans(camera_l);
	}
  {
    ::glutSetWindow(iwin_des);
    imenu_right_click = ::glutCreateMenu(myGlutMenu_Des);
    ::glutAddMenuEntry("Drag", 5);  
    ::glutAddMenuEntry("Dart Cutter", 7);    
    ::glutAddMenuEntry("EditCurve", 1);
    ::glutAddMenuEntry("ChangeToLine", 2);
//    ::glutAddMenuEntry("ChangeToArc", 3);
    ::glutAddMenuEntry("ChangeToPolyline", 4);
    ::glutAddMenuEntry("SmoothPolyline", 6);
    ////
    unsigned int nsl = gui_listner.GetNumberOfSlider();
    for(unsigned int isl=0;isl<nsl;isl++){ 
      const std::string& name = gui_listner.GetNameSlider(isl);
      char str_sl[256];
      sprintf(str_sl, "Slider [ %s ]",name.c_str());
      ::glutAddMenuEntry(str_sl, isl+8);  
    }
    ::glutAttachMenu(GLUT_RIGHT_BUTTON);    
  }
  imode_ope = 0;
  pAnalysis->SetIsLighting(is_lighting);
  pAnalysis->PerformStaticSolver();
	////////////////
	iprob++;
	if( iprob == nprob ){ iprob = 0; }
}

void myGlutSpecial(int Key, int x, int y)
{
	Com::View::CCamera& camera = (glutGetWindow()==iwin_sim) ? camera_r : camera_l;
	switch(Key)
	{
	case GLUT_KEY_PAGE_UP:
		{
			const double tmp_scale = camera.GetScale() * 1.111;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_PAGE_DOWN:
		{
			const double tmp_scale = camera.GetScale() * 0.9;
			camera.SetScale( tmp_scale );
		}
		break;
	case GLUT_KEY_HOME :
		{
			double rot[9];
			camera.RotMatrix33(rot);
			Com::CBoundingBox3D bb = gui_listner.GetBoundingBox(rot);
			camera.SetObjectBoundingBox(bb);
//			drawer_ary.InitTrans(camera);
			camera.Fit();
		}
		break;
	default:
		break;
	}
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}
 
void myGlutDisplay(void)
{
  ::glClearColor(1.f , 1.f, 1.f ,1.0f);      
  
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );
  
	if( glutGetWindow() == iwin_sim ){	
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		Com::View::SetModelViewTransform(camera_r);		    
		gui_listner.Draw(4);    
    ShowTextSimulation();    
	}
	else{
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		Com::View::SetModelViewTransform(camera_l);				
		gui_listner.Draw(1);
    ShowTextDesign();
    drawer_coord.Draw();    
    if( is_rubber_band ){   
      ::glColor3d(1,0,0);
      ::glBegin(GL_LINES);
      ::glVertex2d(pos_pick.x,pos_pick.y);
      ::glVertex2d(pos_relese.x,pos_relese.y);      
      ::glEnd();
    }    
    if( imode_ope == 4 ){
      double min, max;
      const double val = gui_listner.GetValueSlider(islider_active,min,max);
      const double ratio = (val-min)/(max-min)*2-1;      
      ::glMatrixMode(GL_PROJECTION);
      ::glPushMatrix();
      ::glLoadIdentity();      
      ::glMatrixMode(GL_MODELVIEW);
      ::glPushMatrix();
      ::glLoadIdentity();
      ::glColor3d(0,0,0);      
      ::glBegin(GL_LINES);      
      ::glVertex2d(-1,-0.9);
      ::glVertex2d(+1,-0.9);      
      ::glEnd();
      if( is_active_slider == -1 ){ ::glColor3d(0,0,0); }
      else{                         ::glColor3d(1,0,0); }
      ::glBegin(GL_QUADS);
      ::glVertex2d(ratio-0.05,-1.0);
      ::glVertex2d(ratio+0.05,-1.0);
      ::glVertex2d(ratio+0.05,-0.9);
      ::glVertex2d(ratio-0.05,-0.9);
      ::glEnd();
      ::glMatrixMode(GL_MODELVIEW);      
      ::glPopMatrix();
      ::glMatrixMode(GL_PROJECTION);      
      ::glPopMatrix();      
    }
  }

	::glutSwapBuffers();
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case '\033': 
      exit(0);
      break;
    case 'f':
      gui_listner.File_WriteDXF("cloth_pattern.dxf", 4);
      break;
    case 'g':
      if( pAnalysis != 0 ){
        pAnalysis->WriteObjMeshSTL("object.stl", 100);
      }      
      break;      
    case 'k':
      gui_listner.HilightCadTypeID(Cad::NOT_SET,0);
      if( imode_ope == 0 ){ imode_ope = 1; }
      else { imode_ope = 0; }
      break;
    case 'h':
//      gui_listner.Serialize( Com::CSerializer("hoge.txt",true ) );
      break;
    case 'r':
      gui_listner.Solve_fromCad_InitValue();        
      break;
    case 't':
      gui_listner.Solve_fromCad_InterpValue();        
      break;      
    case 'y':
      gui_listner.LoadTimeStamp();
      break;
    case 'u':
      gui_listner.SaveTimeStamp();
      break; 
    case ' ':
      SetNewProblem();
      break;
    case 'a':
      is_animation = !is_animation;
//      gui_listner.Enable_SolveCadChange(is_animation);
      break;
    case 'b':
      if( pAnalysis->GetMode()==CLOTH_INITIAL_LOCATION ){ 
        gui_listner.Solve_fromCad_InitValue();
        pAnalysis->PerformStaticSolver();
        ::glutSetWindow(iwin_des);  ::glutSetCursor(GLUT_CURSOR_WAIT);
        ::glutSetWindow(iwin_sim);  ::glutSetCursor(GLUT_CURSOR_WAIT);
      }
      else{ 
        pAnalysis->SetClothPiecePlacingMode();
      }
      break;
    case 'c':   
      break;
    case 'l':      
      is_lighting = !is_lighting;
      if( is_lighting ){ ::glEnable( GL_LIGHTING); }
      else{              ::glDisable(GL_LIGHTING); }
      pAnalysis->SetIsLighting(is_lighting);
    case 'e':
      pAnalysis->SetIsShowEdge( !pAnalysis->IsShowEdge() );
      break;      
    case 'd':   
      if( pAnalysis->IsDetail() ){
        pAnalysis->SetIsDetail(false,gui_listner.GetCad(),gui_listner.GetMesh());
      }
      else{
        pAnalysis->SetIsDetail(true, gui_listner.GetCad(),gui_listner.GetMesh());        
      }
      break;            
    case 'p':
      is_display_rotation = !is_display_rotation;
      if( is_display_rotation ){
        display_rotation_theta = -10;
        pAnalysis->SetIsDrawPatternBoundary(false);        
        {
          camera_r.SetRotationMode(Com::View::ROT_2DH);
          double rot[9];	camera_r.RotMatrix33(rot);
          Com::CBoundingBox3D bb = pAnalysis->GetBoundingBox(rot);
          camera_r.SetObjectBoundingBox(bb);
          camera_r.Fit();	
          camera_r.SetScale(1.8);
          ////
          ::glutSetWindow(iwin_sim);
          ::glMatrixMode(GL_PROJECTION);
          ::glLoadIdentity();
          Com::View::SetProjectionTransform(camera_r);
          ::glutPostRedisplay();
          //		drawer_coord.SetTrans(camera_r);
        }
        gui_listner.Enable_SolveCadChange(false);
      }
      else{
        pAnalysis->SetIsDrawPatternBoundary(true);        
        gui_listner.Enable_SolveCadChange(true);                
      }
      break;
    default:
      break;
  }
}



void Initialize_OpenGL_Texture(int ival)
{
  //    ReadPPM_SetTexture("cb.ppm");
  std::string tex_name;
  if(      ival == 0 ){ tex_name = "texture/red_point.ppm"; }
  else if( ival == 1 ){ tex_name = "texture/tex1.ppm"; }
  else if( ival == 2 ){ tex_name = "texture/tex2.ppm"; }
  else if( ival == 3 ){ tex_name = "texture/tex3.ppm"; }
  else if( ival == 4 ){ tex_name = "texture/tex4.ppm"; }
  else if( ival == 5 ){ tex_name = "texture/tex5.ppm"; }
  else if( ival == 6 ){ tex_name = "texture/tex6.ppm"; }
  else if( ival == 7 ){ tex_name = "texture/tex7.ppm"; }
  else if( ival == 8 ){ tex_name = "texture/tex8.ppm"; }
  else if( ival == 9 ){ tex_name = "texture/tex9.ppm"; }
  else if( ival == 10){ tex_name = "texture/tex10.ppm"; }  
  if( ival == 0 ){ is_tex_mouse_pos = true; }
  else{ is_tex_mouse_pos = false; }
  
  ::glutSetWindow(iwin_sim);
  Com::View::ReadPPM_SetTexture(tex_name.c_str(),
                                m_texName, m_texWidth,m_texHeight);    
  if( is_tex_mouse_pos ){
    float tex_border_Color[4] = { 1,1,1,1 };
    glTexParameterfv( GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, tex_border_Color ); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);        
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);          
  }
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  
  /////
  ::glutSetWindow(iwin_des);
  Com::View::ReadPPM_SetTexture(tex_name.c_str(),
                                m_texName, m_texWidth,m_texHeight);    
  if( is_tex_mouse_pos ){
    float tex_border_Color[4] = { 1,1,1,1 };
    glTexParameterfv( GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, tex_border_Color );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);        
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);          
  }
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);  
  
  if( pAnalysis == 0 ) return;
  
  
  if( is_tex_mouse_pos ){ 
    pAnalysis->SetColor_FaceFEM(0.8,0.8,0.8); 
    gui_listner.SetColor_CadFace(0.8, 0.8, 0.8);        
    pAnalysis->SetTextureScale_FaceFEM(1);    
    gui_listner.SetTextureScale_CadFace(1);    
  }
  else if( ival == 7 ){
    pAnalysis->SetColor_FaceFEM(1.0,1.0,1.0); 
    gui_listner.SetColor_CadFace(1, 1, 1);    
    pAnalysis->SetTextureScale_FaceFEM(2);
    gui_listner.SetTextureScale_CadFace(2);        
  }
  else if( ival == 10 || ival == 8 ){
    pAnalysis->SetColor_FaceFEM(1.0,1.0,1.0); 
    gui_listner.SetColor_CadFace(1, 1, 1);        
    pAnalysis->SetTextureScale_FaceFEM(4);
    gui_listner.SetTextureScale_CadFace(4);        
  }
  else{
    pAnalysis->SetColor_FaceFEM(1.0,1.0,1.0); 
    gui_listner.SetColor_CadFace(1, 1, 1);    
    pAnalysis->SetTextureScale_FaceFEM(5);
    gui_listner.SetTextureScale_CadFace(5);    
  }  
}


void myGlutMenu_Sim(int ivalue)
{
  if( ivalue < 11 ){
    Initialize_OpenGL_Texture(ivalue); 
  }
  else if( ivalue == 11 ){  // very thin
    std::cout << "very thin" << std::endl;        
    CClothParam cp = pAnalysis->GetParam_Cloth();
    cp.stiff_bend = 0;
    cp.rho = 0.02;    
    pAnalysis->SetParam_Cloth(cp);
  }
  else if( ivalue == 12 ){  // thin (default)
    std::cout << "thin" << std::endl;        
    CClothParam cp = pAnalysis->GetParam_Cloth();
    cp.stiff_bend = 1.0e-10;
    cp.rho = 0.02;    
    pAnalysis->SetParam_Cloth(cp);    
  }
  else if( ivalue == 13 ){  // medium
    std::cout << "mdium" << std::endl;    
    CClothParam cp = pAnalysis->GetParam_Cloth();
    cp.stiff_bend = 1.0e-5;
    cp.rho = 0.02;    
    pAnalysis->SetParam_Cloth(cp);        
  }
  else if( ivalue == 14 ){  // thick
    std::cout << "thick" << std::endl;
    CClothParam cp = pAnalysis->GetParam_Cloth();
    cp.stiff_bend = 1.0e-2;
    cp.stiff_myu = 20;
//    cp.rho = 0.001;
    // sent to takeo    
//    cp.stiff_bend = 5.0e-7;
//    cp.rho = 0.0001;
    pAnalysis->SetParam_Cloth(cp);        
  }
}

int main(int argc, char *argv[])
{
	
	::glutInit(&argc, argv);
  
	
	// left window
	::glutInitWindowPosition(50,50);
	::glutInitWindowSize(500,650);
  camera_l.SetWindowAspect(500.0/650);
	::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	iwin_des = glutCreateWindow("Design Window");
	::glutDisplayFunc(myGlutDisplay);
	::glutIdleFunc(myGlutIdle);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutSpecialFunc(myGlutSpecial);
	::glutMouseFunc(myGlutMouse);
  ::glutPassiveMotionFunc(myGlutPassive);  
	::glutMotionFunc(myGlutMotion);
  
	// right window
	::glutInitWindowPosition(600,50);
	::glutInitWindowSize(500,650);
  camera_r.SetWindowAspect(500.0/650);  
	::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	iwin_sim = glutCreateWindow("Simulation Window");
	::glutDisplayFunc(myGlutDisplay);
	::glutIdleFunc(myGlutIdle);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutSpecialFunc(myGlutSpecial);
	::glutMouseFunc(myGlutMouse);
	::glutMotionFunc(myGlutMotion);  

  {
    ::glutSetWindow(iwin_sim);
    ////
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_LIGHT0);            
    float light0pos[4] = {0,10,+20,0};
    ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);      
    float white[3] = {0.8,0.8,0.8};
    ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
    
    ::glEnable(GL_LIGHT1);    
    float light1pos[4] = {3,10,+10,0}; 
    ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);          
    float white2[3] = {0.4,0.4,0.4};
    ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white2);        
  }
  {
    ::glutSetWindow(iwin_des);        
    glEnable(GL_TEXTURE_2D);    
    ::glutSetWindow(iwin_sim);        
    glEnable(GL_TEXTURE_2D);    
  }
   
  {
    ::glutSetWindow(iwin_sim);
    int imenu_tex = ::glutCreateMenu(myGlutMenu_Sim);
    ::glutAddMenuEntry("Red Point", 0);  
    ::glutAddMenuEntry("Brown Check", 1);
    ::glutAddMenuEntry("Olive Check", 2);
    ::glutAddMenuEntry("Red Blue Line Check", 3);
    ::glutAddMenuEntry("Red Check", 4);
    ::glutAddMenuEntry("Blue Gingham Check", 5);
    ::glutAddMenuEntry("Pink Plait", 6);    
    ::glutAddMenuEntry("Floral Black/White", 7);    
    ::glutAddMenuEntry("Floral Color", 8);    
    ::glutAddMenuEntry("Blue Plait", 9);      
    ::glutAddMenuEntry("Hart", 10);          
    
    int imenu_bend = ::glutCreateMenu(myGlutMenu_Sim);
    ::glutAddMenuEntry("very thin",11);
    ::glutAddMenuEntry("thin (default)",12);
    ::glutAddMenuEntry("medium",13);
    ::glutAddMenuEntry("thick",14);
    
                     
    ::glutCreateMenu(myGlutMenu_Sim);
    ::glutAddSubMenu("texture",imenu_tex);
    ::glutAddSubMenu("cloth thickness",imenu_bend);    
    ::glutAttachMenu(GLUT_RIGHT_BUTTON);    
  }
    

  ::SetNewProblem();
  
  Initialize_OpenGL_Texture(4);

  std::cout << "enter main loop" << std::endl;
	::glutMainLoop();
	return 0;
}
