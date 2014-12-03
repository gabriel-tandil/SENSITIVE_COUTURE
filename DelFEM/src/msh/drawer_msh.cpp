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

////////////////////////////////////////////////////////////////
// DrawerMsh.cpp : メッシュ描画クラス(CDrawerMsh2D,CDrawerMsh3D)の実装
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif
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

#include <assert.h>
#include <iostream>

#include "delfem/drawer_msh.h"

using namespace Msh::View;
using namespace Com;

namespace Msh{
namespace View{

class CDrawPart
{
public : 
	CDrawPart(){ 
		is_selected = false; is_shown = true; 
        r = 0.8; g = 0.8; b = 0.8;
		height = 0;
		line_width = 1;
    }
	CDrawPart(const Msh::CTriAry2D& TriAry);
	CDrawPart(const Msh::CTriAry3D& TriAry);
	CDrawPart(const Msh::CBarAry& BarAry);
	CDrawPart(const Msh::CBarAry3D& BarAry);
	CDrawPart(const Msh::CQuadAry2D& QuadAry);
	CDrawPart(const Msh::CQuadAry3D& QuadAry);
	CDrawPart(const Msh::CTetAry& TetAry);
	CDrawPart(const Msh::CHexAry& HexAry);
	////////////////
	~CDrawPart(){
		if( pIA_Elem != 0 ){ delete[] pIA_Elem; }
		if( pIA_Edge != 0 ){ delete[] pIA_Edge; }
	}
	void DrawElements(){
		if(      type_elem == 1 ){ DrawElements_Bar();  }
		else if( type_elem == 2 ){ DrawElements_Tri();  }
		else if( type_elem == 3 ){ DrawElements_Quad(); }
		else if( type_elem == 4 ){ DrawElements_Tet();  }
		else if( type_elem == 5 ){ DrawElements_Hex();  }
	}
	void DrawElements_Select(){
		if(      type_elem == 1 ){ DrawElements_Select_Bar();  }
		else if( type_elem == 2 ){ DrawElements_Select_Tri();  }
		else if( type_elem == 3 ){ DrawElements_Select_Quad(); }
		else if( type_elem == 4 ){ DrawElements_Select_Tet();  }
		else if( type_elem == 5 ){ DrawElements_Select_Hex();  }
	}
	unsigned int GetElemDim() const {
		if(      type_elem == 1 ){ return 1; }
		else if( type_elem == 2 ){ return 2; }
		else if( type_elem == 3 ){ return 2; }
		else if( type_elem == 4 ){ return 3; }
		else if( type_elem == 5 ){ return 3; }
		return 0;
	}
	void SetHeight(double height){ this->height = height; }
	double GetHeight() const { return height; }
private:
	void DrawElements_Bar();
	void DrawElements_Tri();
	void DrawElements_Quad();
	void DrawElements_Tet();
	void DrawElements_Hex();
	////////////////
	void DrawElements_Select_Bar();
	void DrawElements_Select_Tri();
	void DrawElements_Select_Quad();
	void DrawElements_Select_Tet();
	void DrawElements_Select_Hex();
public:
	bool is_selected;
	bool is_shown;
	std::vector<unsigned int> selec_elem;
	unsigned int id_msh;
	unsigned int id_cad;
    ////////////////
    double r, g, b;
    unsigned int line_width;
	////////////////
	unsigned int nelem;
	unsigned int npoel;
	unsigned int* pIA_Elem;
	unsigned int nedge;
	unsigned int* pIA_Edge;
	////////////////
	double height;
private:
	unsigned int type_elem;  // bar:1，tri:2，quad:3, tet:4, hex:5
};

}
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

Msh::View::CDrawPart::CDrawPart(const Msh::CTetAry& TetAry)
{
	pIA_Elem = 0; 
	pIA_Edge = 0;
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	is_selected = false; is_shown = true; 
//	this->id_cad = TetAry.m_CadLoopID;
	this->id_msh = TetAry.id;			assert( id_msh != 0 );
	this->type_elem = 4;

	////////////////////////////////
	// 面のセット
	nelem = 0;
	for(unsigned int itet=0;itet<TetAry.m_aTet.size();itet++){
		for(unsigned int ielemtet=0;ielemtet<4;ielemtet++){
			if( TetAry.m_aTet[itet].g[ielemtet] == -2 ) continue;
			nelem++;
		}
	}
	pIA_Elem = new unsigned int [nelem*3];
	{
		unsigned int icnt = 0;
		for(unsigned int itet=0;itet<TetAry.m_aTet.size();itet++){
			for(unsigned int ifacetet=0;ifacetet<4;ifacetet++){
				if( TetAry.m_aTet[itet].g[ifacetet] == -2 ) continue;
                pIA_Elem[icnt*3+0] = TetAry.m_aTet[itet].v[ (int)Msh::noelTetFace[ifacetet][0] ];
                pIA_Elem[icnt*3+1] = TetAry.m_aTet[itet].v[ (int)Msh::noelTetFace[ifacetet][1] ];
                pIA_Elem[icnt*3+2] = TetAry.m_aTet[itet].v[ (int)Msh::noelTetFace[ifacetet][2] ];
				icnt++;
			}
		}
	}

	////////////////////////////////
	// 辺のセット
	nedge = nelem*3;
	pIA_Edge = new unsigned int [nedge*2];
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		pIA_Edge[(ielem*3  )*2+0] = pIA_Elem[ielem*3  ];
		pIA_Edge[(ielem*3  )*2+1] = pIA_Elem[ielem*3+1];
		////////////////
		pIA_Edge[(ielem*3+1)*2+0] = pIA_Elem[ielem*3+1];
		pIA_Edge[(ielem*3+1)*2+1] = pIA_Elem[ielem*3+2];
		////////////////
		pIA_Edge[(ielem*3+2)*2+0] = pIA_Elem[ielem*3+2];
		pIA_Edge[(ielem*3+2)*2+1] = pIA_Elem[ielem*3  ];
	}
}

CDrawPart::CDrawPart(const Msh::CHexAry& HexAry)
{
	pIA_Elem = 0; pIA_Edge = 0;
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
//	this->id_cad = TetAry.m_CadLoopID;
//	std::cout << HexAry.id << std::endl;
	this->id_msh = HexAry.id;			assert( id_msh != 0 );
	this->type_elem = 5;

	const unsigned int nfacehex = 6;
	
	////////////////////////////////
	// 面のセット
	nelem = 0;
	for(unsigned int ihex=0;ihex<HexAry.m_aHex.size();ihex++){
	for(unsigned int ifacehex=0;ifacehex<nfacehex;ifacehex++){
		if( HexAry.m_aHex[ihex].g[ifacehex] == -2 ) continue;
		nelem++;
	}
	}
	pIA_Elem = new unsigned int [nelem*4];
	{
		unsigned int icnt = 0;
		for(unsigned int ihex=0;ihex<HexAry.m_aHex.size();ihex++){
		for(unsigned int ifacehex=0;ifacehex<nfacehex;ifacehex++){
			if( HexAry.m_aHex[ihex].g[ifacehex] == -2 ) continue;
            pIA_Elem[icnt*4+0] = HexAry.m_aHex[ihex].v[ (int)Msh::noelHexFace[ifacehex][0] ];
            pIA_Elem[icnt*4+1] = HexAry.m_aHex[ihex].v[ (int)Msh::noelHexFace[ifacehex][1] ];
            pIA_Elem[icnt*4+2] = HexAry.m_aHex[ihex].v[ (int)Msh::noelHexFace[ifacehex][2] ];
            pIA_Elem[icnt*4+3] = HexAry.m_aHex[ihex].v[ (int)Msh::noelHexFace[ifacehex][3] ];
			icnt++;
		}
		}
		assert( icnt == nelem );
	}

	////////////////////////////////
	// 辺のセット
	nedge = nelem*4;
	pIA_Edge = new unsigned int [nedge*2];
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		pIA_Edge[(ielem*4  )*2+0] = pIA_Elem[ielem*4  ];
		pIA_Edge[(ielem*4  )*2+1] = pIA_Elem[ielem*4+1];
		////////////////
		pIA_Edge[(ielem*4+1)*2+0] = pIA_Elem[ielem*4+1];
		pIA_Edge[(ielem*4+1)*2+1] = pIA_Elem[ielem*4+2];
		////////////////
		pIA_Edge[(ielem*4+2)*2+0] = pIA_Elem[ielem*4+2];
		pIA_Edge[(ielem*4+2)*2+1] = pIA_Elem[ielem*4+3];
		////////////////
		pIA_Edge[(ielem*4+3)*2+0] = pIA_Elem[ielem*4+3];
		pIA_Edge[(ielem*4+3)*2+1] = pIA_Elem[ielem*4  ];
	}
}

CDrawPart::CDrawPart(const Msh::CQuadAry3D& QuadAry)
{
	pIA_Elem = 0; pIA_Edge = 0;
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
//	this->id_cad = QuadAry.m_CadLoopID;
	this->id_msh = QuadAry.id;			assert( id_msh != 0 );
	this->type_elem = 3;
	////////////////
	nelem = QuadAry.m_aQuad.size();
	{	// 面のセット
		pIA_Elem = new unsigned int [nelem*4];
		for(unsigned int iquad=0;iquad<nelem;iquad++){
			pIA_Elem[iquad*4+0] = QuadAry.m_aQuad[iquad].v[0];
			pIA_Elem[iquad*4+1] = QuadAry.m_aQuad[iquad].v[1];
			pIA_Elem[iquad*4+2] = QuadAry.m_aQuad[iquad].v[2];
			pIA_Elem[iquad*4+3] = QuadAry.m_aQuad[iquad].v[3];
		}
	}
	{	// 辺のセット
		nedge = nelem*4;
		pIA_Edge = new unsigned int [nedge*2];
		for(unsigned int iquad=0;iquad<nelem;iquad++){
			pIA_Edge[(iquad*4  )*2+0] = QuadAry.m_aQuad[iquad].v[0];
			pIA_Edge[(iquad*4  )*2+1] = QuadAry.m_aQuad[iquad].v[1];
			////////////////
			pIA_Edge[(iquad*4+1)*2+0] = QuadAry.m_aQuad[iquad].v[1];
			pIA_Edge[(iquad*4+1)*2+1] = QuadAry.m_aQuad[iquad].v[2];
			////////////////
			pIA_Edge[(iquad*4+2)*2+0] = QuadAry.m_aQuad[iquad].v[2];
			pIA_Edge[(iquad*4+2)*2+1] = QuadAry.m_aQuad[iquad].v[3];
			////////////////
			pIA_Edge[(iquad*4+3)*2+0] = QuadAry.m_aQuad[iquad].v[3];
			pIA_Edge[(iquad*4+3)*2+1] = QuadAry.m_aQuad[iquad].v[0];
		}
	}
}


CDrawPart::CDrawPart(const Msh::CQuadAry2D& QuadAry)
{
	pIA_Elem = 0; pIA_Edge = 0;
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	this->id_cad = QuadAry.id_l_cad;
	this->id_msh = QuadAry.id;			assert( id_msh != 0 );
	this->type_elem = 3;
	////////////////
	nelem = QuadAry.m_aQuad.size();
	{	// 面のセット
		pIA_Elem = new unsigned int [nelem*4];
		for(unsigned int iquad=0;iquad<nelem;iquad++){
			pIA_Elem[iquad*4+0] = QuadAry.m_aQuad[iquad].v[0];
			pIA_Elem[iquad*4+1] = QuadAry.m_aQuad[iquad].v[1];
			pIA_Elem[iquad*4+2] = QuadAry.m_aQuad[iquad].v[2];
			pIA_Elem[iquad*4+3] = QuadAry.m_aQuad[iquad].v[3];
		}
	}
	{	// 辺のセット
		nedge = nelem*4;
		pIA_Edge = new unsigned int [nedge*2];
		for(unsigned int iquad=0;iquad<nelem;iquad++){
			pIA_Edge[(iquad*4  )*2+0] = QuadAry.m_aQuad[iquad].v[0];
			pIA_Edge[(iquad*4  )*2+1] = QuadAry.m_aQuad[iquad].v[1];
			////////////////
			pIA_Edge[(iquad*4+1)*2+0] = QuadAry.m_aQuad[iquad].v[1];
			pIA_Edge[(iquad*4+1)*2+1] = QuadAry.m_aQuad[iquad].v[2];
			////////////////
			pIA_Edge[(iquad*4+2)*2+0] = QuadAry.m_aQuad[iquad].v[2];
			pIA_Edge[(iquad*4+2)*2+1] = QuadAry.m_aQuad[iquad].v[3];
			////////////////
			pIA_Edge[(iquad*4+3)*2+0] = QuadAry.m_aQuad[iquad].v[3];
			pIA_Edge[(iquad*4+3)*2+1] = QuadAry.m_aQuad[iquad].v[0];
		}
	}
}

CDrawPart::CDrawPart(const Msh::CTriAry2D& TriAry)
{
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	pIA_Elem = 0; pIA_Edge = 0;
	this->id_cad = TriAry.id_l_cad;
	this->id_msh = TriAry.id; assert( id_msh != 0 );
	this->type_elem = 2;
	////////////////
	nelem = TriAry.m_aTri.size();
	{	// 面のセット
		pIA_Elem = new unsigned int [nelem*3];
		for(unsigned int itri=0;itri<nelem;itri++){
			pIA_Elem[itri*3+0] = TriAry.m_aTri[itri].v[0];
			pIA_Elem[itri*3+1] = TriAry.m_aTri[itri].v[1];
			pIA_Elem[itri*3+2] = TriAry.m_aTri[itri].v[2];
		}
	}
	{	// 辺のセット
		nedge = nelem*3;
		pIA_Edge = new unsigned int [nedge*2];
		for(unsigned int itri=0;itri<nelem;itri++){
			pIA_Edge[(itri*3  )*2+0] = TriAry.m_aTri[itri].v[0];
			pIA_Edge[(itri*3  )*2+1] = TriAry.m_aTri[itri].v[1];
			pIA_Edge[(itri*3+1)*2+0] = TriAry.m_aTri[itri].v[1];
			pIA_Edge[(itri*3+1)*2+1] = TriAry.m_aTri[itri].v[2];
			pIA_Edge[(itri*3+2)*2+0] = TriAry.m_aTri[itri].v[2];
			pIA_Edge[(itri*3+2)*2+1] = TriAry.m_aTri[itri].v[0];
		}
	}
}

CDrawPart::CDrawPart(const Msh::CTriAry3D& TriAry)
{
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	pIA_Elem = 0; pIA_Edge = 0;
	this->id_msh = TriAry.id;			assert( id_msh != 0 );
	this->type_elem = 2;
	////////////////
	nelem = TriAry.m_aTri.size();
	{	// 面のセット
		pIA_Elem = new unsigned int [nelem*3];
		for(unsigned int itri=0;itri<nelem;itri++){
			pIA_Elem[itri*3+0] = TriAry.m_aTri[itri].v[0];
			pIA_Elem[itri*3+1] = TriAry.m_aTri[itri].v[1];
			pIA_Elem[itri*3+2] = TriAry.m_aTri[itri].v[2];
		}
	}
	{	// 辺のセット
		nedge = nelem*3;
		pIA_Edge = new unsigned int [nedge*2];
		for(unsigned int itri=0;itri<nelem;itri++){
			pIA_Edge[(itri*3  )*2+0] = TriAry.m_aTri[itri].v[0];
			pIA_Edge[(itri*3  )*2+1] = TriAry.m_aTri[itri].v[1];
			pIA_Edge[(itri*3+1)*2+0] = TriAry.m_aTri[itri].v[1];
			pIA_Edge[(itri*3+1)*2+1] = TriAry.m_aTri[itri].v[2];
			pIA_Edge[(itri*3+2)*2+0] = TriAry.m_aTri[itri].v[2];
			pIA_Edge[(itri*3+2)*2+1] = TriAry.m_aTri[itri].v[0];
		}
	}
}

CDrawPart::CDrawPart(const Msh::CBarAry& BarAry)
{
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	pIA_Elem = 0; pIA_Edge = 0;
	this->id_cad = BarAry.id_e_cad;
	this->id_msh = BarAry.id; assert( id_msh != 0 );
	this->type_elem = 1;

	////////////////
	nelem = BarAry.m_aBar.size();
	pIA_Elem = new unsigned int [nelem*2];
	for(unsigned int ibar=0;ibar<nelem;ibar++){
		pIA_Elem[ibar*2+0] = BarAry.m_aBar[ibar].v[0];
		pIA_Elem[ibar*2+1] = BarAry.m_aBar[ibar].v[1];
	}	
}

CDrawPart::CDrawPart(const Msh::CBarAry3D& BarAry)
{
	is_selected = false; is_shown = true; 
    r = 0.8; g = 0.8; b = 0.8;
	line_width = 1;
	height = 0;
	pIA_Elem = 0; pIA_Edge = 0;
	this->type_elem = 1;
	this->id_cad = BarAry.id_cad;
	this->id_msh = BarAry.id;			assert( id_msh != 0 );
	////////////////
	nelem = BarAry.m_aBar.size();
	pIA_Elem = new unsigned int [nelem*2];
	for(unsigned int ibar=0;ibar<nelem;ibar++){
		pIA_Elem[ibar*2+0] = BarAry.m_aBar[ibar].v[0];
		pIA_Elem[ibar*2+1] = BarAry.m_aBar[ibar].v[1];
	}
}

////////////////////////////////////////////////////////////////

void CDrawPart::DrawElements_Tet()
{
	if( !is_shown ) return;
	// 辺を描画
	::glLineWidth(1);
	if( this->is_selected ){ 
//		std::cout << "DrawElements_Tet selec" << std::endl;
		::glLineWidth(2); 
		::glColor3d(1.0,1.0,0.0); 
	}
	else{ 
//		std::cout << "DrawElements_Tet unselec" << std::endl;
		::glLineWidth(1); 
		::glColor3d(0.0,0.0,0.0);
	}
	::glDrawElements(GL_LINES,nedge*2,GL_UNSIGNED_INT,pIA_Edge); 

    // ピックされた面を描画
	::glColor3d(1.0,0.0,0.0);
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		::glBegin(GL_TRIANGLES);
		::glArrayElement(pIA_Elem[ielem0*3  ]);
		::glArrayElement(pIA_Elem[ielem0*3+1]);
		::glArrayElement(pIA_Elem[ielem0*3+2]);
		::glEnd();
	}

	// 面を描画
	::glColor3d(r,g,b);
	::glDrawElements(GL_TRIANGLES,nelem*3,GL_UNSIGNED_INT,pIA_Elem); 

}

void CDrawPart::DrawElements_Select_Tet(){
	if( !is_shown ) return;
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		::glPushName(ielem);
		::glBegin(GL_TRIANGLES);
		::glArrayElement(pIA_Elem[ielem*3  ]);
		::glArrayElement(pIA_Elem[ielem*3+1]);
		::glArrayElement(pIA_Elem[ielem*3+2]);
		::glEnd();
		::glPopName();
	}
}

////////////////////////////////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CDrawPart::DrawElements_Hex()
{
	if( !is_shown ) return;

    // ピックされた面を描画
	::glColor3d(1.0,0.0,0.0);
	::glBegin(GL_QUADS);
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		::glArrayElement(pIA_Elem[ielem0*4  ]);
		::glArrayElement(pIA_Elem[ielem0*4+1]);
		::glArrayElement(pIA_Elem[ielem0*4+2]);
		::glArrayElement(pIA_Elem[ielem0*4+3]);
	}
	::glEnd();

	// 面を描画
	::glColor3d(r,g,b);
	::glDrawElements(GL_QUADS,nelem*4,GL_UNSIGNED_INT,pIA_Elem); 

	// 辺を描画
	::glLineWidth(1);
	if( this->is_selected ){ 
		::glLineWidth(2); 
		::glColor3d(1.0,1.0,0.0); 
	}
	else{ 
		::glLineWidth(1); 
		::glColor3d(0.0,0.0,0.0);
	}
	::glDrawElements(GL_LINES,nedge*2,GL_UNSIGNED_INT,pIA_Edge); 

}

void CDrawPart::DrawElements_Select_Hex(){
	
	if( !is_shown ) return;

	// 面を描画
	for(unsigned int ielem=0;ielem<nelem;ielem++){
		::glPushName(ielem);
		::glBegin(GL_QUADS);
		::glArrayElement(pIA_Elem[ielem*4  ]);
		::glArrayElement(pIA_Elem[ielem*4+1]);
		::glArrayElement(pIA_Elem[ielem*4+2]);
		::glArrayElement(pIA_Elem[ielem*4+3]);
		::glEnd();
		::glPopName();
	}

	std::cout << "Picked Face " << std::endl;
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		std::cout << iielem << "  ";
		std::cout << pIA_Elem[ielem0*4  ] << " ";
		std::cout << pIA_Elem[ielem0*4+1] << " ";
		std::cout << pIA_Elem[ielem0*4+2] << " ";
		std::cout << pIA_Elem[ielem0*4+3] << std::endl;
	}
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void CDrawPart::DrawElements_Quad()
{	
	if( !is_shown ) return;

	// 辺を描画
	::glLineWidth(1);
	if( this->is_selected ){ 
		::glLineWidth(2); 
		::glColor3d(1.0,1.0,0.0); 
	}
	else{ 
		::glLineWidth(1); 
		::glColor3d(0.0,0.0,0.0);
	}
	::glDrawElements(GL_LINES,nedge*2,GL_UNSIGNED_INT,pIA_Edge); 

	::glColor3d(1.0,0.0,0.0);
	::glBegin(GL_QUADS);
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		::glArrayElement(pIA_Elem[ielem0*4  ]);
		::glArrayElement(pIA_Elem[ielem0*4+1]);
		::glArrayElement(pIA_Elem[ielem0*4+2]);
		::glArrayElement(pIA_Elem[ielem0*4+3]);
	}
	::glEnd();

	// 面を描画
	::glColor3d(r,g,b);
	::glDrawElements(GL_QUADS,nelem*4,GL_UNSIGNED_INT,pIA_Elem); 
}

void CDrawPart::DrawElements_Select_Quad()
{
	if( !is_shown ) return;
	// 面を描画
	for(unsigned int iquad=0;iquad<nelem;iquad++){
		::glPushName(iquad);
		::glBegin(GL_QUADS);
		::glArrayElement(pIA_Elem[iquad*4  ]);
		::glArrayElement(pIA_Elem[iquad*4+1]);
		::glArrayElement(pIA_Elem[iquad*4+2]);
		::glArrayElement(pIA_Elem[iquad*4+3]);
		::glEnd();
		::glPopName();
	}
}

////////////////////////////////
////////////////////////////////

void CDrawPart::DrawElements_Tri()
{
	if( !is_shown ) return;
	
//	std::cout << "DrawElements_Tri line_width" << line_width << std::endl;

	// 辺を描画
	if( this->is_selected ){ 
		::glLineWidth((int)line_width+1); 
		::glColor3d(1.0,1.0,0.0); 
	}
	else{ 
		::glLineWidth((int)line_width); 
		::glColor3d(0.0,0.0,0.0);
	}
	::glDrawElements(GL_LINES,    nedge*2,GL_UNSIGNED_INT,pIA_Edge);

	::glColor3d(1.0,0.0,0.0);
	::glBegin(GL_TRIANGLES);
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		::glArrayElement(pIA_Elem[ielem0*3  ]);
		::glArrayElement(pIA_Elem[ielem0*3+1]);
		::glArrayElement(pIA_Elem[ielem0*3+2]);
	}
	::glEnd();

	// 面を描画
	::glColor3d(0.8,0.8,0.8);
	::glDrawElements(GL_TRIANGLES,nelem*3, GL_UNSIGNED_INT,pIA_Elem );
}

void CDrawPart::DrawElements_Select_Tri()
{
	if( !is_shown ) return;
	// 面を描画
	for(unsigned int itri=0;itri<nelem;itri++){
		::glPushName(itri);
		::glBegin(GL_TRIANGLES);
		::glArrayElement(pIA_Elem[itri*3  ]);
		::glArrayElement(pIA_Elem[itri*3+1]);
		::glArrayElement(pIA_Elem[itri*3+2]);
		::glEnd();
		::glPopName();
	}
}

////////////////////////////////
////////////////////////////////



void CDrawPart::DrawElements_Bar()
{
	if( !is_shown ) return;
	::glLineWidth(2);

	::glColor3d(1.0,0.0,0.0);
	::glBegin(GL_LINES);
	for(unsigned int iielem=0;iielem<selec_elem.size();iielem++){
		unsigned int ielem0 = selec_elem[iielem];
		::glArrayElement(pIA_Elem[ielem0*2  ]);
		::glArrayElement(pIA_Elem[ielem0*2+1]);
	}
	::glEnd();

	if( this->is_selected ){ ::glColor3d(1.0,1.0,0.0); }
	else{ ::glColor3d(0.0,0.0,0.0); }
	::glDrawElements(GL_LINES,nelem*2,GL_UNSIGNED_INT,pIA_Elem); 
}

void CDrawPart::DrawElements_Select_Bar(){

	if( !is_shown ) return;

	for(unsigned int ibar=0;ibar<nelem;ibar++){
		::glPushName(ibar);
		::glBegin(GL_LINES);
		::glArrayElement(pIA_Elem[ibar*2  ]);
		::glArrayElement(pIA_Elem[ibar*2+1]);
		::glEnd();
		::glPopName();
	}
}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CDrawerMsh2D::~CDrawerMsh2D(){
//	std::cout << "~DrawerGL_MHS2D" << std::endl;
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ delete m_pIndexArySet[idp]; }
}

void CDrawerMsh2D::SetLineWidth(unsigned int iwidth){
//	std::cout << "SetLineWidth" << iwidth << std::endl;
//	getchar();
    for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){
        m_pIndexArySet[idp]->line_width = iwidth;
    }
}

bool CDrawerMsh2D::Set(const Msh::CMesher2D &mesh)
{
	this->sutable_rot_mode = 1;	// DrawMode 1 : 2D

    int ilayer_min=0, ilayer_max=0;
	{
		bool is_inited = false;
		const std::vector<Msh::CTriAry2D>& aTriAry = mesh.GetTriArySet();
		for(unsigned int itri=0;itri<aTriAry.size();itri++){
            int ilayer = aTriAry[itri].ilayer;
			if( is_inited ){
				ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
				ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
			}
			else{
				ilayer_min = ilayer;
				ilayer_max = ilayer;
				is_inited = true;
			}
		}
		const std::vector<Msh::CQuadAry2D>& aQuadAry = mesh.GetQuadArySet();
		for(unsigned int iquad=0;iquad<aQuadAry.size();iquad++){
            int ilayer = aQuadAry[iquad].ilayer;
			if( is_inited ){
				ilayer_min = ( ilayer < ilayer_min ) ? ilayer : ilayer_min;
				ilayer_max = ( ilayer > ilayer_max ) ? ilayer : ilayer_max;
			}
			else{
				ilayer_min = ilayer;
				ilayer_max = ilayer;
				is_inited = true;
			}
		}
	}
//	std::cout << ilayer_max << " " << ilayer_min << std::endl;
	double layer_height = 1.0/(ilayer_max-ilayer_min+1); 

	{	// 3角形要素をセット
		const std::vector<Msh::CTriAry2D>& aTriAry = mesh.GetTriArySet();
		for(unsigned int itri=0;itri<aTriAry.size();itri++){
			CDrawPart* dp = new CDrawPart( aTriAry[itri] );
			unsigned int ilayer = aTriAry[itri].ilayer;
			dp->SetHeight( (ilayer-ilayer_min)*layer_height );
			this->m_pIndexArySet.push_back( dp );
		}
	}

	{	// 4角形要素をセット
		const std::vector<Msh::CQuadAry2D>& aQuadAry = mesh.GetQuadArySet();
		for(unsigned int iquad=0;iquad<aQuadAry.size();iquad++){
			CDrawPart* dp = new CDrawPart( aQuadAry[iquad] );
			unsigned int ilayer = aQuadAry[iquad].ilayer;
			dp->SetHeight( (ilayer-ilayer_min)*layer_height );
			this->m_pIndexArySet.push_back( dp );
		}
	}

	{	// 線要素をセット
		const std::vector<Msh::CBarAry>& aBarAry = mesh.GetBarArySet();
		for(unsigned int ibar=0;ibar<aBarAry.size();ibar++){
			double height = 0;
			{
				const int ilayer = aBarAry[ibar].ilayer;
				height += (ilayer-ilayer_min)*layer_height;
				height += 0.01*layer_height;
			}
			CDrawPart* dp = new CDrawPart( aBarAry[ibar] );
			dp->SetHeight(height);
			this->m_pIndexArySet.push_back( dp );
		}
	}

	{	// 頂点をセット
		const std::vector<Msh::SVertex>& aVertex = mesh.GetVertexAry();
		for(unsigned int iver=0;iver<aVertex.size();iver++){
			double height = 0;
			{
				const int ilayer = aVertex[iver].ilayer;
				height += (ilayer-ilayer_min)*layer_height;
				height += 0.01*layer_height;
			}
			CDrawPart_MshVertex dpv;
			dpv.id_cad = aVertex[iver].id_v_cad;
			dpv.id_msh = aVertex[iver].id;
			dpv.id_v = aVertex[iver].v;
			dpv.height = height;
			dpv.is_selected = false;
			this->m_DrawPartCadVertex.push_back( dpv );
		}
	}
	
	{	// 座標をセット
		const std::vector<CVector2D>& aVec2D = mesh.GetVectorAry();
		const unsigned int ndim = 2;
		const unsigned int npoin = aVec2D.size();
		m_vertex_ary.SetSize(npoin,ndim);
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			m_vertex_ary.pVertexArray[ipoin*ndim  ] = aVec2D[ipoin].x;
			m_vertex_ary.pVertexArray[ipoin*ndim+1] = aVec2D[ipoin].y;
		}
	}

	return true;
}		


void CDrawerMsh2D::UpdateCoord( const Msh::CMesher2D& mesh ){
	{	// 座標をセット
		const std::vector<CVector2D>& aVec2D = mesh.GetVectorAry();
		const unsigned int ndim = 2;
		const unsigned int npoin = aVec2D.size();
		if( m_vertex_ary.NDim() != ndim ) return;
		if( m_vertex_ary.NPoin() != npoin ) return;
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			m_vertex_ary.pVertexArray[ipoin*ndim  ] = aVec2D[ipoin].x;
			m_vertex_ary.pVertexArray[ipoin*ndim+1] = aVec2D[ipoin].y;
		}
	}

}

void CDrawerMsh2D::Draw() const
{
	// ライティングの指定
	::glDisable(GL_LIGHTING);
	// 色の指定
	::glColor3d(0.8,0.8,0.8);
	 
	// 片面かどうかの指定
	::glEnable(GL_CULL_FACE); 
	::glCullFace(GL_BACK);
//	::glDisable(GL_CULL_FACE);
	::glEnable(GL_DEPTH_TEST);
//	::glDisable(GL_DEPTH_TEST);

	const unsigned int ndim = m_vertex_ary.NDim();

	// 頂点配列の設定
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);
	// 面と辺の描画
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ 
        if( !is_draw_face && (m_pIndexArySet[idp]->GetElemDim() == 2 ) ){ continue; }
		const double height = m_pIndexArySet[idp]->GetHeight();
		::glTranslated(0,0,+height);
	    m_pIndexArySet[idp]->DrawElements(); 
		::glTranslated(0,0,-height);
	}
	// 点の描画
	::glTranslated(0.0,0.0, 0.2);
	::glPointSize(5);
	::glBegin(GL_POINTS);
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		if( this->m_DrawPartCadVertex[iver].is_selected ){ ::glColor3d(1.0,1.0,0.0); }
		else{ ::glColor3d(0.0,0.0,0.0);	}
		unsigned int ipo0 = this->m_DrawPartCadVertex[iver].id_v;
		::glArrayElement(ipo0);
	}
	::glEnd();
	::glTranslated(0.0,0.0, -0.2);

	::glDisableClientState(GL_VERTEX_ARRAY);

	return;
}

void CDrawerMsh2D::DrawSelection(unsigned int idraw) const
{
	const unsigned int ndim = m_vertex_ary.NDim();

	// モデルの描画
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);

	::glPushName(idraw);
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ 
		::glPushName(idp);
		m_pIndexArySet[idp]->DrawElements_Select(); 
		::glPopName();
	}
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		unsigned int ipo0 = this->m_DrawPartCadVertex[iver].id_v;
		::glPushName(m_pIndexArySet.size()+iver);
		::glTranslated(0.0,0.0, 0.2);
		::glBegin(GL_POINTS);
		::glArrayElement(ipo0);
		::glEnd();
		::glTranslated(0.0,0.0,-0.2);
		::glPopName();
	}
	::glPopName();

	::glDisableClientState(GL_VERTEX_ARRAY);

	return;
}

void CDrawerMsh2D::AddSelected(const int selec_flag[]){
	const unsigned int idp0 = selec_flag[1];
	const unsigned int ielem0 = selec_flag[2];
	if( idp0 < m_pIndexArySet.size() ){
		m_pIndexArySet[idp0]->is_selected = true;
		std::vector<unsigned int>& selec_elem = m_pIndexArySet[idp0]->selec_elem;
		for(std::vector<unsigned int>::iterator itr=selec_elem.begin();itr!=selec_elem.end();itr++){
			if( *itr == ielem0 ){
				selec_elem.erase(itr);
				return;
			}
		}
		selec_elem.push_back(ielem0);
	}
	else{
		const unsigned int iver = idp0-m_pIndexArySet.size();
		this->m_DrawPartCadVertex[iver].is_selected = true;
	}
}

void CDrawerMsh2D::ClearSelected(){
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){
		m_pIndexArySet[idp]->is_selected = false;
		m_pIndexArySet[idp]->selec_elem.clear();
	}
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		this->m_DrawPartCadVertex[iver].is_selected = false;
	}
}


void CDrawerMsh2D::GetMshPartID(int selec_flag[], unsigned int& msh_part_id){
	const unsigned int idp0 = selec_flag[1];
//	const unsigned int ielem0 = selec_flag[2];
	if( idp0 < m_pIndexArySet.size() ){
		msh_part_id = m_pIndexArySet[idp0]->id_msh;
	}
	else{
		const unsigned int iver = idp0-m_pIndexArySet.size();
//	    std::cout << idp0 << " " << m_pIndexArySet.size() << " " << iver << " " << m_DrawPartCadVertex.size() << std::endl;
		assert( iver < m_DrawPartCadVertex.size() );
		msh_part_id = this->m_DrawPartCadVertex[iver].id_msh;
	}
}




////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CDrawerMsh3D::~CDrawerMsh3D(){
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ delete m_pIndexArySet[idp]; }
}

bool CDrawerMsh3D::Set(const Msh::CMesh3D &mesh)
{
	this->sutable_rot_mode = 3;	// DrawMode 3 : 3D

	{	// 線要素をセット
		const std::vector<Msh::CBarAry3D>& aBarAry = mesh.GetBarArySet();
		for(unsigned int ibar=0;ibar<aBarAry.size();ibar++){
			CDrawPart* dp = new CDrawPart( aBarAry[ibar] );
			this->m_pIndexArySet.push_back( dp );
		}
	}

	{	// 3角形要素をセット
		const std::vector<Msh::CTriAry3D>& aTriAry = mesh.GetTriArySet();
		for(unsigned int itri=0;itri<aTriAry.size();itri++){
			CDrawPart* dp = new CDrawPart( aTriAry[itri] );
			this->m_pIndexArySet.push_back( dp );
		}
	}

	{	// 4角形要素をセット
		const std::vector<Msh::CQuadAry3D>& aQuadAry = mesh.GetQuadArySet();
		for(unsigned int iaquad=0;iaquad<aQuadAry.size();iaquad++){
			CDrawPart* dp = new CDrawPart( aQuadAry[iaquad] );
			this->m_pIndexArySet.push_back( dp );
		}
	}


	////////////////////////////////////////////////////////////////
	// 面を含む要素を後に書くと選択時のハイライトで干渉しない．
/*
	{	// 四面体要素をセット
		const std::vector<Msh::CTetAry>& aTetAry = mesh.GetTetArySet();
		for(unsigned int itet=0;itet<aTetAry.size();itet++){
			CDrawPart_MshTet* dp = new CDrawPart_MshTet;
			dp->Set( aTetAry[itet] );
			this->m_pIndexArySet.push_back( dp );
		}
	}
*/
	{	// 六面体要素をセット
		const std::vector<Msh::CHexAry>& aHexAry = mesh.GetHexArySet();
		for(unsigned int ihex=0;ihex<aHexAry.size();ihex++){
			CDrawPart* dp = new CDrawPart( aHexAry[ihex] );
			this->m_pIndexArySet.push_back( dp );
		}
	}
	{	// 頂点をセット
		const std::vector<Msh::SVertex3D>& aVertex = mesh.GetVertexAry();
		for(unsigned int iver=0;iver<aVertex.size();iver++){
			CDrawPart_MshVertex dpv;
			dpv.id_msh = aVertex[iver].id;
			dpv.id_cad = aVertex[iver].id_cad;
			dpv.id_v = aVertex[iver].v;
			dpv.is_selected = false;
			this->m_DrawPartCadVertex.push_back( dpv );
		}
	}
	
	{	// 座標をセット
		const std::vector<CVector3D>& aVec3D = mesh.GetVectorAry();
		const unsigned int ndim = 3;
		const unsigned int npoin = aVec3D.size();
		m_vertex_ary.SetSize(npoin,ndim);
		for(unsigned int ipoin=0;ipoin<npoin;ipoin++){
			m_vertex_ary.pVertexArray[ipoin*ndim  ] = aVec3D[ipoin].x;
			m_vertex_ary.pVertexArray[ipoin*ndim+1] = aVec3D[ipoin].y;
			m_vertex_ary.pVertexArray[ipoin*ndim+2] = aVec3D[ipoin].z;
		}
	}

	return true;
}		
void CDrawerMsh3D::Draw() const
{
	::glEnable(GL_DEPTH_TEST);
	// ライティングの指定
	::glDisable(GL_LIGHTING);
	// 色の指定
	 
	// 片面かどうかの指定
//	::glDisable(GL_CULL_FACE);

	const unsigned int ndim = m_vertex_ary.NDim();

	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);

	// 辺、面の描画
	// TODO : これは辺を先に描画した方が選択したときに綺麗．

    ::glDisable(GL_CULL_FACE); 

	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ 
		const double height = m_pIndexArySet[idp]->height;
		::glTranslated(0,0, height);
		m_pIndexArySet[idp]->DrawElements(); 
		::glTranslated(0,0,-height);
	}

	// 点の描画
	::glPointSize(5);
	::glBegin(GL_POINTS);
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		if( this->m_DrawPartCadVertex[iver].is_selected ){ ::glColor3d(1.0,1.0,0.0); }
		else{ ::glColor3d(0.0,0.0,0.0);	}
		unsigned int ipo0 = this->m_DrawPartCadVertex[iver].id_v;
		::glArrayElement(ipo0);
	}
	::glEnd();

	::glDisableClientState(GL_VERTEX_ARRAY);

	return;
}

void CDrawerMsh3D::DrawSelection(unsigned int idraw) const
{
	const unsigned int ndim = m_vertex_ary.NDim();

	// モデルの描画
	::glEnableClientState(GL_VERTEX_ARRAY);
	::glVertexPointer(ndim,GL_DOUBLE,0,m_vertex_ary.pVertexArray);

	::glPushName(idraw);
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){ 
		::glPushName(idp);
		m_pIndexArySet[idp]->DrawElements_Select(); 
		::glPopName();
	}
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		::glPushName(m_pIndexArySet.size()+iver);
		const unsigned int ipo0 = this->m_DrawPartCadVertex[iver].id_v;
		::glBegin(GL_POINTS);
		::glArrayElement( ipo0 );
		::glEnd();
		::glPopName();
	}
	::glPopName();

	::glDisableClientState(GL_VERTEX_ARRAY);

	return;
}

void CDrawerMsh3D::AddSelected(const int selec_flag[]){
	const unsigned int idp = selec_flag[1];
	const unsigned int ielem = selec_flag[2];
	if( idp < m_pIndexArySet.size() ){
		m_pIndexArySet[idp]->is_selected = true;
		std::vector<unsigned int>& selec_elem = m_pIndexArySet[idp]->selec_elem;
		for(std::vector<unsigned int>::iterator itr=selec_elem.begin();itr!=selec_elem.end();itr++){
			if( *itr == ielem ){
				selec_elem.erase(itr);
				return;
			}
		}
		selec_elem.push_back(ielem);
	}
}

void CDrawerMsh3D::ClearSelected(){
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){
		m_pIndexArySet[idp]->is_selected = false;
		m_pIndexArySet[idp]->selec_elem.clear();
	}
	for(unsigned int iver=0;iver<this->m_DrawPartCadVertex.size();iver++){
		this->m_DrawPartCadVertex[iver].is_selected = false;
	}
}


void CDrawerMsh3D::Hide(unsigned int id_msh_hide){
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){
		if( m_pIndexArySet[idp]->id_msh == id_msh_hide ){
			m_pIndexArySet[idp]->is_shown = false;
		}
	}
}

void CDrawerMsh3D::SetColor(unsigned int id_msh, double r, double g, double b){
	for(unsigned int idp=0;idp<m_pIndexArySet.size();idp++){
		if( m_pIndexArySet[idp]->id_msh == id_msh ){
			m_pIndexArySet[idp]->r = r;
			m_pIndexArySet[idp]->g = g;
			m_pIndexArySet[idp]->b = b;
		}
	}
}
