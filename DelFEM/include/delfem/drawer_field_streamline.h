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
@brief 流線で速度場を可視化するクラス(Fem::Field::View::CDrawerFace, Fem::Field::View::CDrawerFaceContour)のインターフェース
@author Nobuyuki Umetani
*/

#if !defined(DRAWER_FIELD_STREAMLINE_H)
#define DRAWER_FIELD_STREAMLINE_H

#include <memory>

#include "delfem/drawer_field.h"
/*
class CElSuP{
public:
	CElSuP(){
		nno = 0;
		elsup_ind = 0;
		elsup = 0;
	}
	void Set(const Fem::Field::CElemAry::CElemSeg& es){
		nno = es.GetMaxNoes()+1;
		std::cout << "nno : " << nno << std::endl;
		elsup_ind = new unsigned int [nno+1];
		for(unsigned int ino=0;ino<nno+1;ino++){ elsup_ind[ino] = 0; }
		const unsigned int nnoes = es.GetSizeNoes();
		unsigned int* noes = new unsigned int [nnoes];
		for(unsigned int ielem=0;ielem<es.GetSizeElem();ielem++){
			es.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				elsup_ind[ noes[inoes]+1 ] += 1;
			}
		}
		for(unsigned int ino=0;ino<nno;ino++){ 
			elsup_ind[ino+1] += elsup_ind[ino];  
		}
		const unsigned int nelsup = elsup_ind[nno];
		elsup = new unsigned int [nelsup];
		for(unsigned int ielem=0;ielem<es.GetSizeElem();ielem++){
			es.GetNodes(ielem,noes);
			for(unsigned int inoes=0;inoes<nnoes;inoes++){
				const unsigned int ino0 = noes[inoes];
				const unsigned int ielsup0 = elsup_ind[ino0];
				elsup[ielsup0] = ielem;
				elsup_ind[ino0] += 1;
			}
		}
		delete[] noes;
		for(int ino=nno;ino>0;ino--){ 
			elsup_ind[ino] = elsup_ind[ino-1]; 
		}
		elsup_ind[0] = 0;
//		for(unsigned int ino=0;ino<nno;ino++){
//			std::cout << ino << " --> ";
//			for(unsigned int ielsup=elsup_ind[ino];ielsup<elsup_ind[ino+1];ielsup++){
//				std::cout << elsup[ielsup] << " ";
//			}
//			std::cout << std::endl;
//		}
	}
public:
	unsigned int nno;
	unsigned int* elsup_ind;
	unsigned int* elsup;
};
*/

namespace Fem{
namespace Field{
namespace View{
	
void MakeStreamLine(unsigned int id_field_velo, const Fem::Field::CFieldWorld& world,
					std::vector< std::vector<double> >& aStLine);

//! 流線を描画するクラス
class CDrawerStreamline : public CDrawerField
{
public:
	CDrawerStreamline();
	CDrawerStreamline(unsigned int id_field, const Fem::Field::CFieldWorld& world ){
		this->m_IdFieldVelo = id_field;
		this->Update(world);
	}
	virtual ~CDrawerStreamline(){}

	virtual bool Update(const Fem::Field::CFieldWorld& world){
		MakeStreamLine(m_IdFieldVelo, world, aStLine);
		return true;
	}
	Com::CBoundingBox3D GetBoundingBox( double rot[] ) const;
	virtual void AddSelected(const int selec_flag[]){}
	virtual void ClearSelected(){}
	virtual void DrawSelection(unsigned int idraw) const {}
	virtual void Draw() const;
private:
public:
private:
	unsigned int m_IdFieldVelo;	
	std::vector< std::vector<double> > aStLine;
	int ilayer_min, ilayer_max;
};


}
}
}

#endif
