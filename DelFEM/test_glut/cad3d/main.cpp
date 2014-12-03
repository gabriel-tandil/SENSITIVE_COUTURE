
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdlib> //(exit)

#if defined(__APPLE__) && defined(__MACH__)
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "delfem/camera.h"
#include "delfem/drawer_gl_utility.h"
#include "delfem/glut_utility.h"

#include "delfem/serialize.h"
#include "delfem/cad_obj3d.h"
#include "delfem/drawer_cad3d.h"

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;
int imodifier;
Com::View::CDrawerArray drawer_ary;
Cad::CCadObj3D cad_3d;
unsigned int id_loop_selected = 0;
Com::CVector3D picked_pos, mouse_pos;
unsigned int imode = 0;

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

void myGlutDisplay(void)
{
//	::glClearColor(0.2, .7, 0.7, 1.0);
	::glClearColor(1.0, 1.0, 1.0, 1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

  ShowBackGround();  
  ShowFPS();
  

  {
    ::glDisable(GL_LIGHTING);
    ::glLineWidth(2);
    /*
    std::vector<unsigned int> aIdL = cad_3d.GetAryElemID(Cad::LOOP);
    for(unsigned int iil=0;iil<aIdL.size();iil++){
      unsigned int id_l = aIdL[iil];
      const Cad::CLoop3D& l = cad_3d.GetLoop(id_l);
      const Com::CVector3D& o = l.org; 
      const Com::CVector3D& n = l.normal;     
      double r = 0.1;
      ::glBegin(GL_LINES);
      ::glColor3d(1,0,0);
      ::glVertex3d(o.x,o.y,o.z);
      ::glVertex3d(o.x+n.x*r,o.y+n.y*r,o.z+n.z*r);
      ::glEnd();
    }
     */
    if( imode == 1 ){
      if( cad_3d.IsElemID(Cad::LOOP,id_loop_selected) ){
        const Cad::CLoop3D& l = cad_3d.GetLoop(id_loop_selected);
        const Com::CVector3D& o0 = l.org; 
        const Com::CVector3D& n0 = l.normal;
        Com::CVector3D x0 = l.dirx;
        Com::CVector3D y0 = Cross(n0,x0);
        x0 = Dot(mouse_pos-picked_pos,x0)*x0;
        y0 = Dot(mouse_pos-picked_pos,y0)*y0;
        const Com::CVector3D& p = picked_pos;
        const Com::CVector3D& px = picked_pos+x0;
        const Com::CVector3D& py = picked_pos+y0;
        const Com::CVector3D& pxy = picked_pos+x0+y0;        
        ::glLineWidth(1);
        ::glColor3d(0,0,0);
        ::glBegin(GL_LINES);
        ::glVertex3d(p.x,p.y,p.z);
        ::glVertex3d(px.x,px.y,px.z);
        ::glVertex3d(p.x,p.y,p.z);
        ::glVertex3d(py.x,py.y,py.z);        
        ::glVertex3d(pxy.x,pxy.y,pxy.z);
        ::glVertex3d(px.x,px.y,px.z);
        ::glVertex3d(pxy.x,pxy.y,pxy.z);
        ::glVertex3d(py.x,py.y,py.z);                
        ::glEnd();
      }
      /*
      ::glBegin(GL_POINTS);
      ::glColor3d(1,0,0);
      ::glVertex3d(picked_pos.x,picked_pos.y,picked_pos.z);
      ::glColor3d(1,0,1);      
      ::glVertex3d(mouse_pos.x,mouse_pos.y,mouse_pos.z);      
      ::glEnd();
       */
    }
    ::glEnable(GL_LIGHTING);        
  }

	drawer_ary.Draw();


	glutSwapBuffers();
}

void myGlutIdle(){
	::glutPostRedisplay();
}

void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	if(      imodifier == GLUT_ACTIVE_CTRL ){
		camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
	else if( imodifier == GLUT_ACTIVE_SHIFT ){
		camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
	}
  else if( cad_3d.IsElemID(Cad::LOOP,id_loop_selected) ){    
    if( imode == 0 ){
      if( cad_3d.IsElemID(Cad::LOOP,id_loop_selected) ){
        Com::CVector3D nv = cad_3d.GetLoop(id_loop_selected).normal;
        double n[3] = { nv.x, nv.y, nv.z };
        double r[9];  camera.RotMatrix33(r);
        double nr[3] = { r[0]*n[0]+r[1]*n[1]+r[2]*n[2], r[3]*n[0]+r[4]*n[1]+r[5]*n[2], r[6]*n[0]+r[7]*n[1]+r[8]*n[2] };
        double del[2] = { mov_end_x-mov_begin_x, mov_end_y-mov_begin_y };
        nr[0] /= camera.GetHalfViewHeight()*camera.GetWindowAspect();
        nr[1] /= camera.GetHalfViewHeight();
        if( nr[0]*nr[0] + nr[1]*nr[1] > 1.0e-5 ){
          double r = (nr[0]*del[0]+nr[1]*del[1])/(nr[0]*nr[0]+nr[1]*nr[1]);
          cad_3d.LiftLoop(id_loop_selected,nv*r);
          drawer_ary.Clear();
          drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );    
        }              
      }
    }
    if( imode == 1 ){      
      if( cad_3d.IsElemID(Cad::LOOP,id_loop_selected) ){
        const Cad::CLoop3D& l = cad_3d.GetLoop(id_loop_selected);
        const Com::CVector3D& o0 = l.org; 
        const Com::CVector3D& n0 = l.normal;
        Com::CMatrix3 r = camera.RotMatrix33();      
        Com::CVector3D d1 = r.MatVecTrans( Com::CVector3D(0,0,-1) );
        double hvh = camera.GetHalfViewHeight();
        double asp = camera.GetWindowAspect();
        Com::CVector3D cp = -camera.GetCenterPosition()+Com::CVector3D(hvh*asp*mov_end_x,hvh*mov_end_y,0);
        Com::CVector3D o1 = r.MatVecTrans( cp ) + camera.GetObjectCenter();
        double tmp0 = Dot(d1,n0);
        double tmp1 = Dot(o0-o1,n0);
        double tmp2 = tmp1/tmp0;
        mouse_pos = o1+d1*tmp2;
      }
    }
  }
	mov_begin_x = mov_end_x;
	mov_begin_y = mov_end_y;
	::glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	mov_begin_x = (2.0*x-win_w)/win_w;
	mov_begin_y = (win_h-2.0*y)/win_h;
  imodifier = ::glutGetModifiers();
	if( state == GLUT_DOWN && imodifier == 0 ){
		const unsigned int size_buffer = 2048;
		unsigned int selec_buffer[size_buffer];
		Com::View::PickPre(size_buffer,selec_buffer, x,y,5,5, camera);
		drawer_ary.DrawSelection();
		std::vector<Com::View::SSelectedObject> aSelecObj = Com::View::PickPost(selec_buffer, x,y, camera );
		drawer_ary.ClearSelected();
		if( aSelecObj.size() > 0 ){ 
      drawer_ary.AddSelected( aSelecObj[0].name ); 
      if( aSelecObj[0].name[1] == 3 ){ id_loop_selected = aSelecObj[0].name[2]; picked_pos = aSelecObj[0].picked_pos; }
      else{ id_loop_selected = 0; }      
      mouse_pos = picked_pos;
    }
	}
  if( state == GLUT_UP && imodifier == 0 ){
    if( imode == 1 ){
      if( cad_3d.IsElemID(Cad::LOOP,id_loop_selected) ){
        const Cad::CLoop3D& l = cad_3d.GetLoop(id_loop_selected);
        const Com::CVector2D v0 = l.Project(picked_pos);
        const Com::CVector2D v1 = l.Project(mouse_pos);      
        cad_3d.AddRectLoop(id_loop_selected,v0,v1);
        id_loop_selected = 0;
        drawer_ary.Clear();
        drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );          
      }      
    }
  }
}

bool SetNewProblem()
{
	const unsigned int nprob = 4;
	static unsigned int iprob = 0;
  
  std::cout << "SetNewProblem() " << iprob << std::endl;

  cad_3d.Clear();
	if( iprob == 0 ){
		unsigned int id_s0 = cad_3d.AddCuboid(1,1,1);        
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );    
		drawer_ary.InitTrans(camera);
  }
  if( iprob == 1 ){
		unsigned int id_s0 = cad_3d.AddCuboid(1,1,0.7);
    unsigned int id_l1 = 0;
    {
      std::vector<Com::CVector3D> aVec;
      aVec.push_back( Com::CVector3D(0.5,0.2,0.7) );
      aVec.push_back( Com::CVector3D(0.8,0.2,0.7) );
      aVec.push_back( Com::CVector3D(0.8,0.8,0.7) );
      aVec.push_back( Com::CVector3D(0.5,0.8,0.7) );
      id_l1 = cad_3d.AddPolygon(aVec,6);
    }
          
    unsigned int id_l2 = 0;
    {
      std::vector<Com::CVector3D> aVec;
      aVec.push_back( Com::CVector3D(0.2,1.0,0.3) );
      aVec.push_back( Com::CVector3D(0.2,1.0,0.5) );
      aVec.push_back( Com::CVector3D(0.8,1.0,0.5) );
      aVec.push_back( Com::CVector3D(0.8,1.0,0.3) );
      id_l2 = cad_3d.AddPolygon(aVec,3);
    }

    cad_3d.LiftLoop(id_l1,cad_3d.GetLoop(id_l1).normal* 0.2);
    cad_3d.LiftLoop(id_l2,cad_3d.GetLoop(id_l2).normal*-0.2);    

    {
      const unsigned int id_v1 = cad_3d.AddPoint( Cad::EDGE,1, Com::CVector3D(0,0.7,0  ) );
      const unsigned int id_v2 = cad_3d.AddPoint( Cad::EDGE,9, Com::CVector3D(0,0.7,0.7) );
      Cad::CBRepSurface::CResConnectVertex res = cad_3d.ConnectVertex(id_v1, id_v2);
      unsigned int id_l_lift = cad_3d.GetIdLoop_Edge(res.id_e_add,true);
      cad_3d.LiftLoop(id_l_lift,cad_3d.GetLoop(id_l_lift).normal*0.2);
    }
 
    {
      const unsigned int id_v1 = cad_3d.AddPoint( Cad::EDGE,9, Com::CVector3D(0,0.3,0.7 ) );
      const unsigned int id_v2 = cad_3d.AddPoint( Cad::EDGE,5, Com::CVector3D(0,0.0,0.2 ) );
      const unsigned int id_v3 = cad_3d.AddPoint( Cad::LOOP,17,Com::CVector3D(0,0.3,0.2 ) );
      cad_3d.ConnectVertex(id_v1, id_v3);          
      Cad::CBRepSurface::CResConnectVertex res = cad_3d.ConnectVertex(id_v2, id_v3);    
      unsigned int id_l_lift = res.id_l_add;
      cad_3d.LiftLoop(id_l_lift,cad_3d.GetLoop(id_l_lift).normal*-0.2);      
    }

    cad_3d.LiftLoop(6,Com::CVector3D(0,0,-0.1));
    
    {
      const unsigned int id_v1 = cad_3d.AddPoint( Cad::EDGE,3, Com::CVector3D(1.0,0.2,0.0 ) );
      const unsigned int id_v2 = cad_3d.AddPoint( Cad::LOOP,4, Com::CVector3D(1.0,0.2,0.2 ) );      
      const unsigned int id_v3 = cad_3d.AddPoint( Cad::LOOP,4, Com::CVector3D(1.0,0.6,0.2 ) );        
      const unsigned int id_v4 = cad_3d.AddPoint( Cad::EDGE,3, Com::CVector3D(1.0,0.6,0.0 ) );
      cad_3d.ConnectVertex(id_v1, id_v2);
      cad_3d.ConnectVertex(id_v2, id_v3);    
      Cad::CBRepSurface::CResConnectVertex res = cad_3d.ConnectVertex(id_v3, id_v4); 
      unsigned int id_l_lift = res.id_l_add;
      cad_3d.LiftLoop(id_l_lift,cad_3d.GetLoop(id_l_lift).normal*-0.1);
    }
    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );    
		drawer_ary.InitTrans(camera);
	}
  if( iprob == 2 ){
		unsigned int id_s0 = cad_3d.AddCuboid(1,1,1);        
    cad_3d.AddRectLoop(1, Com::CVector2D(0.2,0.2), Com::CVector2D(0.8,1.2) );
    cad_3d.LiftLoop(1,cad_3d.GetLoop(1).normal*0.1);        
//    cad_3d.LiftLoop(7,cad_3d.GetLoop(7).normal*0.1);
		drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );    
		drawer_ary.InitTrans(camera);
  }
  if( iprob == 3 ){
		unsigned int id_s0 = cad_3d.AddCuboid(1,1,1);       
    cad_3d.AddRectLoop(6, 
                       cad_3d.GetLoop(6).Project(Com::CVector3D(0.2,0.2,1)),
                       cad_3d.GetLoop(6).Project(Com::CVector3D(0.8,0.8,1)) );
    cad_3d.AddRectLoop(7, 
                       cad_3d.GetLoop(7).Project(Com::CVector3D(0.3,0.3,1)),
                       cad_3d.GetLoop(7).Project(Com::CVector3D(0.7,0.7,1)) );    
//    cad_3d.AddRectLoop(7, Com::CVector2D(0.3,0.3), Com::CVector2D(0.7,0.7) );    
    cad_3d.LiftLoop(7,cad_3d.GetLoop(7).normal*0.1);        
    drawer_ary.Clear();
		drawer_ary.PushBack( new Cad::View::CDrawer_Cad3D(cad_3d) );    
		drawer_ary.InitTrans(camera);    
  }
    
  
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);

	iprob++;
	if( iprob == nprob ) iprob=0;

	return true;
}

void myGlutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'q':
  case 'Q':
  case '\033':  /* '\033' ÇÕ ESC ÇÃ ASCII ÉRÅ[Éh */
	  exit(0);
	  break;
  case ' ':
	  SetNewProblem();
	  break;
  default:
    break;
  }
}


void myGlutSpecial(int Key, int x, int y)
{
	switch(Key)
	{
    case GLUT_KEY_PAGE_UP:
      if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
        if( camera.IsPers() ){
          const double tmp_fov_y = camera.GetFovY() + 10.0;
          camera.SetFovY( tmp_fov_y );
        }
      }
      else{
        const double tmp_scale = camera.GetScale() * 0.9;
        camera.SetScale( tmp_scale );
      }
      break;
    case GLUT_KEY_PAGE_DOWN:
      if( ::glutGetModifiers() && GLUT_ACTIVE_SHIFT ){
        if( camera.IsPers() ){
          const double tmp_fov_y = camera.GetFovY() - 10.0;
          camera.SetFovY( tmp_fov_y );
        }
      }
      else{
        const double tmp_scale = camera.GetScale() * 1.111;
        camera.SetScale( tmp_scale );
      }
      break;
    case GLUT_KEY_HOME :
      camera.Fit();
      break;
    case GLUT_KEY_END :
      if( camera.IsPers() ) camera.SetIsPers(false);
      else{ camera.SetIsPers(true); }
      break;
    default:
      break;
	}
	
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}


void myGlutMenu(int value)
{  
  switch(value) {
    case 0:
      imode = 0;  
      std::cout << "Drag Loop" << std::endl;
      break;      
    case 1:
      imode = 1;  
      std::cout << "Add Loop" << std::endl;
      break;
  }
}

int main(int argc,char* argv[])
{
	// Initailze GLUT
	::glutInitWindowPosition(200,200);
	::glutInitWindowSize(400, 300);
	::glutInit(&argc, argv);	
    ::glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	::glutCreateWindow("Cad 3D View");

	// Set callback function
	::glutMotionFunc(myGlutMotion);
	::glutMouseFunc(myGlutMouse);
	::glutDisplayFunc(myGlutDisplay);
	::glutReshapeFunc(myGlutResize);
	::glutKeyboardFunc(myGlutKeyboard);
	::glutIdleFunc(myGlutIdle);
  ::glutSpecialFunc(myGlutSpecial);

	SetNewProblem();
  
  {
    ::glEnable(GL_LIGHTING);
    ::glEnable(GL_LIGHT0);    
    float light0pos[4] = {0,0,+20,0};
    ::glLightfv(GL_LIGHT0, GL_POSITION, light0pos);      
    float white1[3] = {0.7,0.7,0.7};
    ::glLightfv(GL_LIGHT0, GL_DIFFUSE, white1);
    ////
    ::glEnable(GL_LIGHT1);    
    float light1pos[4] = {0,10,0,0};
    float white2[3] = {0.9,0.9,0.9};    
    ::glLightfv(GL_LIGHT1, GL_POSITION, light1pos);      
    ::glLightfv(GL_LIGHT1, GL_DIFFUSE, white2);        
  }  
  
  ::glutCreateMenu(myGlutMenu);
  ::glutAddMenuEntry("Drag Loop", 0);  
  ::glutAddMenuEntry("Add Loop", 1);   
  ::glutAttachMenu(GLUT_RIGHT_BUTTON);

	// Enter main loop
	::glutMainLoop();
	return 0;
}
