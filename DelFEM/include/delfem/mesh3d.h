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
@brief ３次元メッシュクラス(Msh::CMesh3D)のインターフェース
@author Nobuyuki Umetani
*/
#if !defined(MESH_3D_H)
#define MESH_3D_H

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

#include <vector>

#include "delfem/vector2d.h"
#include "delfem/quaternion.h"

#include "delfem/mesh_interface.h"
#include "delfem/mesher2d.h"           //　ここをInterfaceクラスにできるように頑張る
#include "delfem/msh/meshkernel2d.h"
#include "delfem/msh/meshkernel3d.h"

////////////////////////////////////////////////

namespace Msh{

/*! 
@addtogroup Msh3D
*/
// @{

//! 頂点構造体
struct SVertex3D{
public:
	SVertex3D() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;   //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3) {０か奇数になるはず}
	unsigned int v;	//!< 点のID，今は一つだけだけど、そのうち配列にして複数扱えるようにしたいな。
};

//! 線要素配列
class CBarAry3D{
public:
	CBarAry3D() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;   //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	std::vector<SBar> m_aBar;	//!< ３次元線分要素配列
};

//! ４面体配列クラス	
class CTetAry{
public:
	CTetAry() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;   //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	std::vector<STet> m_aTet;	//!< ３次元四面体要素配列
};

//! 六面体配列クラス
class CHexAry{
public:
	CHexAry() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;  //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	std::vector<SHex> m_aHex;	//!< ３次元６面体要素配列
};

//! ３次元三角形配列クラス
class CTriAry3D{
public:
	CTriAry3D() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;  //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	std::vector<STri3D> m_aTri;	//!< ３次元３角形要素配列
};

//! ３次元４角形配列クラス
class CQuadAry3D{
public:
	CQuadAry3D() : id(0), id_cad(0), id_msh_before_extrude(0), inum_extrude(0){}	// IDを０で初期化しないと、GetMaxIDが失敗する
public:
	unsigned int id;	//!< ID
	unsigned int id_cad;	//!< CADの頂点ID（CADに関連されてなければ０）
	unsigned int id_msh_before_extrude;
    unsigned int inum_extrude;  //!< 突き出されてない(0), 底面(1), 側面(2), 上面(3)
	std::vector<SQuad3D> m_aQuad;	//!< ３次元４角形要素配列
};

////////////////////////////////

//! ３次元メッシュクラス
class CMesh3D  : public IMesh
{
public:
	CMesh3D(){}
	virtual ~CMesh3D(){ this->Clear(); }

	virtual unsigned int GetDimention() const{ return 3; }
	virtual void GetCoord(std::vector<double>& coord) const{
		unsigned int nnode = aVec.size();
		coord.resize( nnode*3 );
		for(unsigned int inode=0;inode<nnode;inode++){
			coord[inode*3  ] = aVec[inode].x;
			coord[inode*3+1] = aVec[inode].y;
			coord[inode*3+2] = aVec[inode].z;
		}
	}
	virtual void GetInfo(unsigned int id_msh, 
        unsigned int& id_cad, unsigned int& id_msh_be, unsigned int& inum_ext,
		int& ilayer) const
    { 
		const int itype = m_ElemType[id_msh];
		const int iloc = m_ElemLoc[id_msh];
		if(      itype == 0 ){ id_cad=m_aVertex[ iloc].id_cad;  id_msh_be=m_aVertex[ iloc].id_msh_before_extrude;  inum_ext=m_aVertex[ iloc].inum_extrude; }
		else if( itype == 1 ){ id_cad=m_aBarAry[ iloc].id_cad;  id_msh_be=m_aBarAry[ iloc].id_msh_before_extrude;  inum_ext=m_aBarAry[ iloc].inum_extrude; }
        else if( itype == 2 ){ id_cad=m_aTriAry[ iloc].id_cad;  id_msh_be=m_aTriAry[ iloc].id_msh_before_extrude;  inum_ext=m_aTriAry[ iloc].inum_extrude; }
        else if( itype == 3 ){ id_cad=m_aQuadAry[iloc].id_cad;  id_msh_be=m_aQuadAry[iloc].id_msh_before_extrude;  inum_ext=m_aQuadAry[iloc].inum_extrude; }
        else if( itype == 4 ){ id_cad=m_aTetAry[ iloc].id_cad;  id_msh_be=m_aTetAry[ iloc].id_msh_before_extrude;  inum_ext=m_aTetAry[ iloc].inum_extrude; }
        else if( itype == 5 ){ id_cad=m_aHexAry[ iloc].id_cad;  id_msh_be=m_aHexAry[ iloc].id_msh_before_extrude;  inum_ext=m_aHexAry[ iloc].inum_extrude; }
		else{ assert(0); }
		ilayer = 0;	// ３次元メッシュは層構造を持たないので
	}
	virtual MSH_TYPE GetConnectivity(unsigned int id_msh, std::vector<int>& lnods) const;
	virtual std::vector<unsigned int> GetIncludeElemIDAry(unsigned int id) const{
		{	// エラーの場合は空の配列を返す
			std::vector<unsigned int> id_ary;
			if( id >= m_ElemLoc.size() ) return id_ary;
			if( m_ElemLoc[id] == -1 ) return id_ary;
			if( id >= this->m_include_relation.size() ) return id_ary;
		}
		return m_include_relation[id];
	}
	virtual std::vector<unsigned int> GetAry_ID() const{
		std::vector<unsigned int> id_ary;
		for(unsigned int id=1;id<m_ElemLoc.size();id++){
			if( m_ElemLoc[id] == -1 ) continue;
			id_ary.push_back(id);
		}
		return id_ary;
	}

    void Scale(double r){
        for(unsigned int ivec=0;ivec<aVec.size();ivec++){
            aVec[ivec] *= r;
        }
    }
    void Translate(double x, double y, double z){
        for(unsigned int ivec=0;ivec<aVec.size();ivec++){
            aVec[ivec].x += x;
            aVec[ivec].y += y;
            aVec[ivec].z += z;
        }
    }
    void Rotate(double t, double x, double y, double z){
        Com::CVector3D n(x,y,z);
        n.Normalize();
        n *= t;
        Com::CQuaternion q(n);
        double m[9];
        q.RotMatrix33(m);
        for(unsigned int ivec=0;ivec<aVec.size();ivec++){
            Com::CVector3D vin = aVec[ivec];
            double x = m[0]*vin.x + m[1]*vin.y + m[2]*vin.z;
            double y = m[3]*vin.x + m[4]*vin.y + m[5]*vin.z;
            double z = m[6]*vin.x + m[7]*vin.y + m[8]*vin.z;
            aVec[ivec].x = x;
            aVec[ivec].y = y;
            aVec[ivec].z = z;
        }
    }


	////////////////////////////////

	void Clear();	//!< メッシュの情報を全てクリアして削除する

	// const関数
	bool IsID(unsigned int id) const;	//!< メッシュのIDを調べる

	// 要素配列、節点配列に関するGetメソッド
	const std::vector<CTetAry>& GetTetArySet()     const { return m_aTetAry;  }	//!< 四面体要素配列の配列を得る
	const std::vector<CHexAry>& GetHexArySet()     const { return m_aHexAry;  }	//!< 六面体要素配列の配列を得る
	const std::vector<CTriAry3D>& GetTriArySet()   const { return m_aTriAry;  }	//!< ３次元３角形要素配列の配列を得る
	const std::vector<CQuadAry3D>& GetQuadArySet() const { return m_aQuadAry; }	//!< ３次元３角形要素配列の配列を得る
	const std::vector<CBarAry3D>& GetBarArySet()    const { return m_aBarAry;  }	//!< 線分要素配列の配列を得る
	const std::vector<SVertex3D>& GetVertexAry()    const { return m_aVertex;  }	//!< 点要素の配列を得る
	const std::vector<Com::CVector3D>& GetVectorAry()  const { return aVec;  }	//!< 座標値の配列を得る

	////////////////////////////////////////////////////////////////
	// IOメソッド

	// 読み込み書き出し
	bool Serialize( Com::CSerializer& serialize );	
	// GiDメッシュの読み込み
	bool ReadFromFile_GiDMsh(const std::string& file_name);
    bool ReadFromFile_TetgenMsh(const std::string& file_name);
protected:	
	unsigned int FindMaxID() const;
	unsigned int GetFreeObjID();
	void MakeElemLocationType(); // 要素の場所と種類をハッシュ（IDが引数）している配列を初期化
protected:
	std::vector<int> m_ElemType;	// vertex(0) bar(1) tri(2) quad(3) tet(4) hex(5)	いつでも有効　
	std::vector<int> m_ElemLoc;		// index of elem_ary				いつでも有効　
	std::vector< std::vector<unsigned int> > m_include_relation;	// どの要素配列がどの要素配列を含んでるか

	std::vector<SVertex3D>	m_aVertex;	// type(0)
	std::vector<CBarAry3D>	m_aBarAry;	// type(1)
	std::vector<CTriAry3D>	m_aTriAry;	// type(2)
	std::vector<CQuadAry3D> m_aQuadAry;	// type(3)
	std::vector<CTetAry>	m_aTetAry;	// type(4)
	std::vector<CHexAry>	m_aHexAry;	// type(5)

	std::vector<Com::CVector3D> aVec;
};

////////////////////////////////

/*
// 本当は２Ｄと３Ｄを分けたくないが，描画がそれに対応していないので，とりあえず次元ごとにローダーを作る
class CMesh3D_FileLoader : public CMesh3D
{
public:
    CMesh3D_FileLoader(){
    }
    bool ReadUms(const std::string& fname)
private:
};
*/

////////////////////////////////

class CMesh3D_Extrude : public CMesh3D
{
public:
	bool UpdateMeshCoord(const CMesher2D& msh_2d);
	bool UpdateMeshConnectivity(const CMesher2D& msh_2d);

	// メッシュを突き出す
	// height : メッシュ高さ
	// elen : 高さ方向のメッシュ幅
	bool Extrude(const CMesher2D& msh_2d, const double height, const double elen);	// 本当は面を突き出すようにしたい。
private:
	bool MakeDirectionFlagForExtrusion(
		const CMesher2D& msh_2d,
		std::vector<std::vector<int> >& aTriFlgAry, 
		std::vector< std::vector<int> >& aBarFlgAry);
	bool MakeConnectivitryTriTetForExtrusion(
		const CMesher2D& msh_2d,
		const std::vector<std::vector<int> >& aTriFlgAry, 
		const std::vector< std::vector<int> >& aBarFlgAry,
		const unsigned int ndiv);
	bool Extrude_Tri( const CMesher2D& msh_2d, const double height, const double elen);
	bool Extrude_Quad(const CMesher2D& msh_2d, const double height, const double elen);
};

class CMesher3D : public CMesh3D
{
public:
	bool ReadFile_PLY( const std::string& file_name );
	bool CutMesh(double elen);	// メッシュを切る
	bool HomogenizeSurface(double elen);
private:
	bool EraseSameLocationPointOnSurface();
};

//! @}
};

#endif
