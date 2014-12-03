
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

#if defined(__VISUALC__)
#pragma warning( disable : 4786 ) 
#endif
#define for if(0);else for

#include <iostream>
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
#include "delfem/glut_utility.h"

#include "delfem/cad_obj2d.h"
#include "delfem/mesher2d.h"
#include "delfem/mesh3d.h"

#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/field_value_setter.h"
#include "delfem/drawer_field.h"
#include "delfem/drawer_field_face.h"
#include "delfem/drawer_field_edge.h"
#include "delfem/drawer_field_vector.h"

using namespace Fem::Field;

Com::View::CCamera camera;
double mov_begin_x, mov_begin_y;


void myGlutMotion( int x, int y ){
	GLint viewport[4];
	::glGetIntegerv(GL_VIEWPORT,viewport);
	const int win_w = viewport[2];
	const int win_h = viewport[3];
	const double mov_end_x = (2.0*x-win_w)/win_w;
	const double mov_end_y = (win_h-2.0*y)/win_h;
	camera.MouseRotation(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
//	camera.MousePan(mov_begin_x,mov_begin_y,mov_end_x,mov_end_y); 
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
}

void myGlutResize(int w, int h)
{
	camera.SetWindowAspect((double)w/h);
	glViewport(0, 0, w, h);
	::glMatrixMode(GL_PROJECTION);
	::glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	glutPostRedisplay();
}

View::CDrawerArrayField drawer_ary;

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
		drawer_ary.InitTrans(camera);
		camera.Fit();
		break;
	case GLUT_KEY_END :
		if( camera.IsPers() ) camera.SetIsPers(false);
		else{ camera.SetIsPers(true); }
		break;
	default:
		break;
	}
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	Com::View::SetProjectionTransform(camera);
	::glutPostRedisplay();
}

unsigned int id_field_val;
unsigned int id_base;
Fem::Field::CFieldWorld world;
std::vector<Fem::Field::CFieldValueSetter> field_value_setter_ary;

double cur_time = 0.0;
bool is_animation = true;

void myGlutDisplay(void)
{
	::glClearColor(1.0, 1.0, 1.0 ,1.0);
	::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	::glEnable(GL_DEPTH_TEST);

	::glEnable(GL_POLYGON_OFFSET_FILL );
	::glPolygonOffset( 1.1f, 4.0f );

	::glMatrixMode(GL_MODELVIEW);
	::glLoadIdentity();
	Com::View::SetModelViewTransform(camera);

	if( is_animation ){
		cur_time += 0.1;
    for(unsigned int iset=0;iset<field_value_setter_ary.size();iset++){
      field_value_setter_ary[iset].ExecuteValue(cur_time,world);
    }
		drawer_ary.Update(world);
	}
	drawer_ary.Draw();
	ShowFPS();

	glutSwapBuffers();
}

void myGlutIdle(){
	glutPostRedisplay();
}

bool SetNewProblem()
{
	const unsigned int nprob = 13;
	static unsigned int iprob = 0;

	static int id_val_bc0=0, id_val_bc1=0, id_val_bc2=0;
	
	if( iprob == 0 )
  {
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
    assert( world.IsIdField(id_field_val) );
    {
      field_value_setter_ary.clear();
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(10*sqrt(x^2+y^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
  }
	else if( iprob == 1 ){    
    {
      field_value_setter_ary.clear();
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(2*PI*x-t)*sin(2*PI*y-t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
	}
	else if( iprob == 2 )
	{
		Msh::CMesher2D msh_2d;
		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
//		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		world.Clear();
		id_base = world.AddMesh( msh_2d );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);        
    {
      field_value_setter_ary.clear();
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(0.5*sqrt((x+1)^2+y^2)-0.1*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 3 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		Msh::CMesher2D mesh_2d(cad_2d,0.07);
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d,1.0,0.07);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);  
    {
      field_value_setter_ary.clear();
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(10*sqrt(x^2+y^2+z^2)-PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
  else if( iprob == 4 )
	{
		Msh::CMesher2D mesh_2d;
		mesh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(mesh_2d,5.0,0.5);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR);
    {
      field_value_setter_ary.clear();      
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(1.0*sqrt(x^2+y^2+z^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 5 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
			cad_2d.AddVertex(Cad::EDGE, 1, Com::CVector2D(0.0, -0.5) );
			cad_2d.AddVertex(Cad::EDGE, 3, Com::CVector2D(0.0,  0.5) );
			cad_2d.ConnectVertex_Line(5,6);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.02) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR); assert( world.IsIdField(id_field_val) );
    unsigned int id_field0 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(1,Cad::LOOP)); assert( world.IsIdField(id_field0) ); 
    unsigned int id_field1 = world.GetPartialField(id_field_val,conv.GetIdEA_fromCad(2,Cad::LOOP)); assert( world.IsIdField(id_field1) ); 
    {
      field_value_setter_ary.clear();
      Fem::Field::CFieldValueSetter fvs0(id_field0,world);
      fvs0.SetMathExp("sin(10*sqrt((x+0.5)^2+y^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs0);
      Fem::Field::CFieldValueSetter fvs1(id_field1,world);
      fvs1.SetMathExp("sin(10*sqrt((x-0.5)^2+y^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs1); 
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 6 )
	{
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			cad_2d.AddPolygon( vec_ary );
			cad_2d.AddVertex(Cad::EDGE, 1, Com::CVector2D(0.0, -0.5) );
			cad_2d.AddVertex(Cad::EDGE, 3, Com::CVector2D(0.0,  0.5) );
			cad_2d.ConnectVertex_Line(5,6);
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.05) );
		const Fem::Field::CIDConvEAMshCad& conv = world.GetIDConverter(id_base);
    id_field_val = world.MakeField_FieldElemAry(id_base,conv.GetIdEA_fromCad(1,Cad::LOOP),Fem::Field::SCALAR);    
    {
      field_value_setter_ary.clear();      
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(10*sqrt((x+0.5)^2+y^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }        
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 7 ){
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,2,VECTOR2,VALUE,BUBBLE);
    {
      Fem::Field::CFieldValueSetter fvs(id_field_grad,world);
      fvs.SetMathExp("0.1*sin(t)", 0, Fem::Field::VALUE,world);
      fvs.SetMathExp("0.1*cos(t)", 1, Fem::Field::VALUE,world);      
      field_value_setter_ary.push_back(fvs);
    }            
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 8 ){
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_tet.msh");
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		id_field_val = world.MakeField_FieldElemDim(id_base,3,SCALAR);
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,3,VECTOR3,VALUE,BUBBLE);    
    field_value_setter_ary.clear();          
    {
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(t+0.5*x)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }                
    {
      Fem::Field::CFieldValueSetter fvs(id_field_grad,world);
      fvs.SetGradient(id_field_val, world);      
      field_value_setter_ary.push_back(fvs);
    }
		drawer_ary.Clear();
//		drawer_ary.PushBack( new View::CDrawerFaceContour(id_field_val,world,-1.0,1.0) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_grad,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 9 ){
		Msh::CMesher3D mesh_3d;
		mesh_3d.ReadFromFile_GiDMsh("../input_file/cylinder_hex.msh");
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
		id_field_val = world.MakeField_FieldElemDim(id_base,3,SCALAR);
		unsigned int id_field_grad = world.MakeField_FieldElemDim(id_field_val,3,VECTOR3,VALUE,BUBBLE);    
    field_value_setter_ary.clear();              
    {
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(t+0.5*x)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }
    {
      Fem::Field::CFieldValueSetter fvs(id_field_grad,world);
      fvs.SetGradient(id_field_val, world);      
      field_value_setter_ary.push_back(fvs);
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new Fem::Field::View::CDrawerVector(id_field_grad,world) );
		drawer_ary.PushBack( new Fem::Field::View::CDrawerEdge(id_field_grad,true,world) );
		drawer_ary.InitTrans( camera );
	}
	else if( iprob == 10 )
	{
		Msh::CMesher2D msh_2d;
//		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		world.Clear();
		id_base = world.AddMesh( msh_2d );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR,VALUE,BUBBLE);
    field_value_setter_ary.clear();                  
    {
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(x+y-0.1*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world,id_field_val,-1.0,1.0) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 11 )
	{
		Msh::CMesher2D msh_2d;
//		msh_2d.ReadFromFile_GiDMsh("../input_file/rect_quad.msh");
		msh_2d.ReadFromFile_GiDMsh("../input_file/hexa_tri.msh");
		Msh::CMesh3D_Extrude mesh_3d;
		mesh_3d.Extrude(msh_2d,5.0,1);
		world.Clear();
		id_base = world.AddMesh( mesh_3d );
    id_field_val = world.MakeField_FieldElemDim(id_base,3,Fem::Field::SCALAR,VALUE,BUBBLE);    
    field_value_setter_ary.clear();                      
    {
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(0.5*sqrt(x^2+y^2+z^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
//		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world) );
		drawer_ary.InitTrans( camera);
	}
	else if( iprob == 12 )
	{	// scalar field and vector field
		Cad::CCadObj2D cad_2d;
 		{
			std::vector<Com::CVector2D> vec_ary;
			vec_ary.push_back( Com::CVector2D(-0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5,-0.5) );
			vec_ary.push_back( Com::CVector2D( 0.5, 0.5) );
			vec_ary.push_back( Com::CVector2D(-0.5, 0.5) );
			const unsigned int id_l0 = cad_2d.AddPolygon( vec_ary ).id_l_add;
		}
		world.Clear();
		id_base = world.AddMesh( Msh::CMesher2D(cad_2d,0.1) );
    id_field_val = world.MakeField_FieldElemDim(id_base,2,Fem::Field::SCALAR,VALUE,CORNER|BUBBLE);  assert( world.IsIdField(id_field_val) );
		unsigned int id_field_vec = world.MakeField_FieldElemDim(id_base,2,Fem::Field::VECTOR2,VALUE,CORNER|BUBBLE); assert( world.IsIdField(id_field_vec) );    
    field_value_setter_ary.clear();                      
    {
      Fem::Field::CFieldValueSetter fvs(id_field_val,world);
      fvs.SetMathExp("sin(10*sqrt(x^2+y^2)-2*PI*t)", 0, Fem::Field::VALUE,world);
      field_value_setter_ary.push_back(fvs);
    }                   
    {
      Fem::Field::CFieldValueSetter fvs(id_field_vec,world);
      fvs.SetMathExp("0.05*sin(t)", 0, Fem::Field::VALUE,world);
      fvs.SetMathExp("0.05*cos(t)", 1, Fem::Field::VALUE,world);      
      field_value_setter_ary.push_back(fvs);
    }    
		drawer_ary.Clear();
		drawer_ary.PushBack( new View::CDrawerFace(id_field_val,true,world, id_field_val,-1.0,1.0) );
		drawer_ary.PushBack( new View::CDrawerVector(id_field_vec,world) );
		drawer_ary.InitTrans( camera );
	}

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
  case 'a' :
      is_animation = !is_animation;
      break;
  case ' ':
	  SetNewProblem();
	  ::glMatrixMode(GL_PROJECTION);
	  ::glLoadIdentity();
	  Com::View::SetProjectionTransform(camera);
  default:
    break;
  }
}

int main(int argc,char* argv[])
{
	SetNewProblem();

	glutInitWindowPosition(200,200);
	glutInitWindowSize(400, 400);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
	glutCreateWindow("Field View");

	glutKeyboardFunc(myGlutKeyboard);
	glutIdleFunc(myGlutIdle);
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutResize);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);
	glutSpecialFunc(myGlutSpecial);
	
	glutMainLoop();
	return 0;
}
