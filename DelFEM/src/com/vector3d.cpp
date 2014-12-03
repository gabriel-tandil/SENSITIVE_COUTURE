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
// Vector3D.cpp: CVector3D クラスのインプリメンテーション
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 )
#endif

#include <cassert>
#include <math.h>
#include <iostream>
#include <stack>

#include "delfem/vector3d.h"

using namespace Com;

////////////////////////////////////////////////////////////////////
// メンバ関数のフレンド関数
////////////////////////////////////////////////////////////////////

namespace Com{

bool operator == (const CVector3D& lhs, const CVector3D& rhs){
	if( fabs(lhs.x - rhs.x) < NEARLY_ZERO
		&& fabs(lhs.y - rhs.y) < NEARLY_ZERO
		&& fabs(lhs.z - rhs.z) < NEARLY_ZERO )
		return true;
	else return false;
}

bool operator != (const CVector3D& lhs, const CVector3D& rhs){
	if( lhs == rhs )	return false;
	else return true;
}


}

//////////////////////////////////////////////////////////////////////
//	メンバ関数の非フレンド関数
//////////////////////////////////////////////////////////////////////

namespace Com{

void CVector3D::Normalize()
{
	double mag;

	mag = Length();
	x /= mag;
	y /= mag;
	z /= mag;
}

void CVector3D::SetZero()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

}


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

COctTree::COctTree(){	
	m_aCell.reserve( 1024 );
	m_aCell.resize(1,0);
	m_aVec.reserve( 2048 );
}

void COctTree::SetBoundingBox( const CBoundingBox3D& bb )
{
	x_min = bb.x_min;  x_max = bb.x_max;
	y_min = bb.y_min;  y_max = bb.y_max;
	z_min = bb.z_min;  z_max = bb.z_max;
	return;
}

bool COctTree::Check() const
{
	CBoundingBox3D bb(0.0,0.0, 0.0,0.0, 0.0,0.0);
	std::vector<unsigned int> aIndVec;
	for(unsigned int icell=0;icell<m_aCell.size();icell++){
		this->GetBoundaryOfCell(icell, bb);
		aIndVec.clear();
		this->GetAllPointInCell(icell,aIndVec);
		for(unsigned int iivec=0;iivec<aIndVec.size();iivec++){
			unsigned int ivec = aIndVec[iivec];
			assert( bb.IsInside(m_aVec[ivec].second ) );
		}
	}
	return true;
}

void COctTree::GetAllPointInCell(unsigned int icell_in, std::vector<unsigned int>& aIndexVec ) const
{
	assert( icell_in < m_aCell.size() );
	if( m_aCell[icell_in].size == -1 ){
		for(unsigned int iloc=0;iloc<8;iloc++){
			unsigned int icell0 = m_aCell[icell_in].data[iloc];
			assert( icell0 < m_aCell.size() );
			assert( m_aCell[icell0].iparent == icell_in );
			this->GetAllPointInCell(icell0,aIndexVec);
		}
	}
	else{
		assert( m_aCell[icell_in].size >= 0 && m_aCell[icell_in].size <= 8 );
		unsigned int nloc = m_aCell[icell_in].size;
		for(unsigned int iloc=0;iloc<nloc;iloc++){
			aIndexVec.push_back( m_aCell[icell_in].data[iloc] );
		}
	}
}

int COctTree::GetIndexCell_IncludePoint( const CVector3D& VecIns ) const
{
	if( VecIns.x < x_min || VecIns.x > x_max || 
		VecIns.y < y_min || VecIns.y > y_max || 
		VecIns.z < z_min || VecIns.z > z_max ){ return -1; }

	unsigned int icell0 = 0;
	double x_min0=x_min,x_max0=x_max, y_min0=y_min,y_max0=y_max, z_min0=z_min,z_max0=z_max;
	for(;;){
		assert( VecIns.x >= x_min0 && VecIns.x <= x_max0 );
		assert( VecIns.y >= y_min0 && VecIns.y <= y_max0 );
		assert( VecIns.z >= z_min0 && VecIns.z <= z_max0 );
		if( m_aCell[icell0].size >= 0 ){	// このセルにはデータが入っている
//			std::cout << x_min0 << " " << x_max0 << "  " << y_min0 << " " << y_max0 << "  " << z_min0 << " " << z_max0 << std::endl;
			return icell0;
		}
		else{	// このセルは８つに分けられているので下のセルを探す
			const double x_cent=(x_max0+x_min0)*0.5, y_cent=(y_max0+y_min0)*0.5, z_cent=(z_max0+z_min0)*0.5;
			int icell1;
			if( VecIns.x <= x_cent ){
			if( VecIns.y <= y_cent ){
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[0];
				x_max0 = x_cent; y_max0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[1];
				x_max0 = x_cent; y_max0 = y_cent; z_min0 = z_cent;
			}
			}
			else{
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[2];
				x_max0 = x_cent; y_min0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[3];
				x_max0 = x_cent; y_min0 = y_cent; z_min0 = z_cent;
			}
			}
			}
			else{
			if( VecIns.y <= y_cent ){
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[4];
				x_min0 = x_cent; y_max0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[5];
				x_min0 = x_cent; y_max0 = y_cent; z_min0 = z_cent;
			}
			}
			else{
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[6];
				x_min0 = x_cent; y_min0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[7];
				x_min0 = x_cent; y_min0 = y_cent; z_min0 = z_cent;
			}
			}
			}
			icell0 = icell1;
		}
	}	
}

bool COctTree::IsPointInSphere( double radius, const CVector3D& vec ) const
{
	assert( radius > 0.0 );
	if( radius <= 0.0 ) return false;
	const int icell_in = this->GetIndexCell_IncludePoint(vec);
	if( icell_in != -1 ){
		assert( icell_in >= 0 && (unsigned int)icell_in < this->m_aCell.size() );
		assert( m_aCell[icell_in].size >= 0 );	// このセルにはデータが入っている
		CBoundingBox3D bb(0.0,0.0, 0.0,0.0, 0.0,0.0);
		this->GetBoundaryOfCell(icell_in,bb);
		assert( bb.IsInside(vec) );
		assert( m_aCell[icell_in].size >= 0 );
		// セルの中の点が球の中に入るか調べる
		for(unsigned int iloc=0;iloc<(unsigned int)m_aCell[icell_in].size;iloc++){
			unsigned int ivec0 = m_aCell[icell_in].data[iloc];
			assert( ivec0 < m_aVec.size() );
			const CVector3D& vec0 = m_aVec[ivec0].second;
			assert( bb.IsInside(vec0) );
			const double sq_dist = SquareDistance( vec, vec0 );
			if( sq_dist < radius*radius ) return true;
		}
		// 球がセルの中に収まっている場合
		if( radius < bb.x_max-vec.x && radius < vec.x-bb.x_min &&
			radius < bb.y_max-vec.y && radius < vec.y-bb.y_min && 
			radius < bb.z_max-vec.z && radius < vec.z-bb.z_min )
		{
			return false;
		}
	}

	// 干渉があるデータセルを全て検索する
	std::stack<unsigned int> stack_int_sect;
	{
		std::stack< std::pair<unsigned int,CBoundingBox3D> > stack_next;	// 干渉があるけど、データ型かポインタ型かわからないセル
		stack_next.push( std::make_pair(0,CBoundingBox3D(x_min,x_max,y_min,y_max,z_min,z_max)) );
		for(;!stack_next.empty();){
			const unsigned int icell_cur = stack_next.top().first;
			const CBoundingBox3D bb_cur = stack_next.top().second;
			stack_next.pop();
			assert( bb_cur.IsPossibilityIntersectSphere(vec,radius) );
			// 干渉を確認するルーティンがここに欲しい．
			if( m_aCell[icell_cur].size >= 0 ){	// データ型
				stack_int_sect.push( icell_cur );
				continue;
			}
			assert( m_aCell[icell_cur].size == -1 );
			const double x_cent=(bb_cur.x_min+bb_cur.x_max)*0.5;
			const double y_cent=(bb_cur.y_min+bb_cur.y_max)*0.5;
			const double z_cent=(bb_cur.z_min+bb_cur.z_max)*0.5;
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_max=x_cent;  bb.y_max=y_cent;  bb.z_max=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[0],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_max=x_cent;  bb.y_max=y_cent;  bb.z_min=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[1],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_max=x_cent;  bb.y_min=y_cent;  bb.z_max=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[2],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_max=x_cent;  bb.y_min=y_cent;  bb.z_min=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[3],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_min=x_cent;  bb.y_max=y_cent;  bb.z_max=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[4],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_min=x_cent;  bb.y_max=y_cent;  bb.z_min=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[5],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_min=x_cent;  bb.y_min=y_cent;  bb.z_max=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[6],bb) );
				}
			}
			{
				CBoundingBox3D bb = bb_cur;
				bb.x_min=x_cent;  bb.y_min=y_cent;  bb.z_min=z_cent;
				if( bb.IsPossibilityIntersectSphere(vec,radius) ){
					stack_next.push( std::make_pair(m_aCell[icell_cur].data[7],bb) );
				}
			}
		}
	}

	std::cout << "Hoge" << stack_int_sect.size() << std::endl;
	for(;!stack_int_sect.empty();){
		const unsigned int icell0 = stack_int_sect.top();
		stack_int_sect.pop();
		assert( icell0 >= 0 && (unsigned int)icell0 < this->m_aCell.size() );
		assert( m_aCell[icell0].size >= 0 );	// このセルにはデータが入っている
		CBoundingBox3D bb(0.0,0.0, 0.0,0.0, 0.0,0.0);
		this->GetBoundaryOfCell(icell0,bb);
		assert( m_aCell[icell0].size >= 0 );
		// セルの中の点が球の中に入るか調べる
		for(unsigned int iloc=0;iloc<(unsigned int)m_aCell[icell0].size;iloc++){
			unsigned int ivec0 = m_aCell[icell0].data[iloc];
			assert( ivec0 < m_aVec.size() );
			const CVector3D& vec0 = m_aVec[ivec0].second;
			assert( bb.IsInside(vec0) );
			const double sq_dist = SquareDistance( vec, vec0 );
			if( sq_dist < radius*radius ) return true;
		}
	}	
	return false;
}


void COctTree::GetBoundaryOfCell(unsigned int icell_in, CBoundingBox3D& bb ) const
{
	assert( icell_in < m_aCell.size() );
	std::stack<unsigned int> loc_stack;
	unsigned int icell0 = icell_in;
	for(;;){
		if( icell0 == 0 ) break;
		unsigned int icell_p = m_aCell[icell0].iparent;
		assert( icell_p < m_aCell.size() );
		assert( m_aCell[icell_p].size == -1 );
        unsigned int iloc0=0;
		for(;iloc0<8;iloc0++){
            if( m_aCell[icell_p].data[iloc0] == (int)icell0 ) break;
		}
		assert( iloc0 != 8 );
        loc_stack.push(iloc0);
		icell0 = icell_p;
	}
	assert( icell0 == 0 );
	bb.x_max=x_max; bb.x_min=x_min;   
	bb.y_max=y_max; bb.y_min=y_min;  
	bb.z_max=z_max; bb.z_min=z_min;
	for(;!loc_stack.empty();){
		unsigned int iloc0 = loc_stack.top();
		loc_stack.pop(); 
		assert( iloc0 < 8 );
		const double x_cent=(bb.x_min+bb.x_max)*0.5,  y_cent=(bb.y_min+bb.y_max)*0.5,  z_cent=(bb.z_min+bb.z_max)*0.5;
		if(      iloc0 == 0 ){ bb.x_max=x_cent;  bb.y_max=y_cent;  bb.z_max=z_cent; }
		else if( iloc0 == 1 ){ bb.x_max=x_cent;  bb.y_max=y_cent;  bb.z_min=z_cent; }
		else if( iloc0 == 2 ){ bb.x_max=x_cent;  bb.y_min=y_cent;  bb.z_max=z_cent; }
		else if( iloc0 == 3 ){ bb.x_max=x_cent;  bb.y_min=y_cent;  bb.z_min=z_cent; }
		else if( iloc0 == 4 ){ bb.x_min=x_cent;  bb.y_max=y_cent;  bb.z_max=z_cent; }
		else if( iloc0 == 5 ){ bb.x_min=x_cent;  bb.y_max=y_cent;  bb.z_min=z_cent; }
		else if( iloc0 == 6 ){ bb.x_min=x_cent;  bb.y_min=y_cent;  bb.z_max=z_cent; }
		else if( iloc0 == 7 ){ bb.x_min=x_cent;  bb.y_min=y_cent;  bb.z_min=z_cent; }
		else{ assert(0); }
		assert( m_aCell[icell0].size == -1 );
		unsigned int icell1 = m_aCell[icell0].data[iloc0];
		assert( m_aCell[icell1].iparent == icell0 );
		icell0 = icell1;
	}
	assert( icell0 == icell_in );
}


//  -1：成功
//  -2：範囲外
// 0～：ダブってる点の番号
int COctTree::InsertPoint( unsigned int ipo_ins, const CVector3D& VecIns )
{
//	std::cout << "InsertPoint : " << VecIns.x << " " << VecIns.y << " " << VecIns.z << std::endl;
	unsigned int icell0 = 0;
	double x_min0=x_min,x_max0=x_max, y_min0=y_min,y_max0=y_max, z_min0=z_min,z_max0=z_max;
	for(;;){
		assert( VecIns.x >= x_min0 && VecIns.x <= x_max0 );
		assert( VecIns.y >= y_min0 && VecIns.y <= y_max0 );
		assert( VecIns.z >= z_min0 && VecIns.z <= z_max0 );
		if( m_aCell[icell0].size >= 0 ){	// このセルにはデータが入っている
			int nloc = m_aCell[icell0].size;
			if( nloc > 0 ){
				unsigned int iloc_min = 0;
				double min_dist = Distance( VecIns, m_aVec[ m_aCell[icell0].data[0] ].second );
				for(unsigned int iloc=1;iloc<(unsigned int)nloc;iloc++){
					const double dist = Distance( VecIns, m_aVec[ m_aCell[icell0].data[iloc] ].second );
					if( dist < min_dist ){
						min_dist = dist;
						iloc_min = iloc;
					}
				}
				if( min_dist < 1.0e-20 ){ return m_aVec[ m_aCell[icell0].data[iloc_min] ].first; }
			}
			assert( nloc >= 0 );
			if( nloc < 8 ){
				m_aCell[icell0].data[nloc] = m_aVec.size();
				m_aVec.push_back( std::make_pair(ipo_ins,VecIns) );
				m_aCell[icell0].size += 1;	// データを追加して終了
				break; 
			}
			// このセルを分割する．
			const unsigned int ncell = m_aCell.size();
			m_aCell.resize(m_aCell.size()+8,icell0);
			double x_cent=(x_max0+x_min0)*0.5, y_cent=(y_max0+y_min0)*0.5, z_cent=(z_max0+z_min0)*0.5;
			for(int i=0;i<8;i++){
				const unsigned int ipo0 = m_aCell[icell0].data[i];
				if( m_aVec[ipo0].second.x <= x_cent ){
				assert( m_aVec[ipo0].second.x >= x_min0 );
				if( m_aVec[ipo0].second.y <= y_cent ){
				assert( m_aVec[ipo0].second.y >= y_min0 );
				if( m_aVec[ipo0].second.z <= z_cent ){
					assert( m_aVec[ipo0].second.z >= z_min0 );
					const unsigned int iloc1 = m_aCell[ncell+0].size;
					m_aCell[ncell+0].data[iloc1] = ipo0;
					m_aCell[ncell+0].size += 1;
				}
				else{
					assert( m_aVec[ipo0].second.z <= z_max0 );
					const unsigned int iloc1 = m_aCell[ncell+1].size;
					m_aCell[ncell+1].data[iloc1] = ipo0;
					m_aCell[ncell+1].size += 1;
				}
				}
				else{
				assert( m_aVec[ipo0].second.y <= y_max0 );
				if( m_aVec[ipo0].second.z <= z_cent ){
				assert( m_aVec[ipo0].second.z >= z_min0 );
					const unsigned int iloc1 = m_aCell[ncell+2].size;
					m_aCell[ncell+2].data[iloc1] = ipo0;
					m_aCell[ncell+2].size += 1;
				}
				else{
				assert( m_aVec[ipo0].second.z <= z_max0 );
					const unsigned int iloc1 = m_aCell[ncell+3].size;
					m_aCell[ncell+3].data[iloc1] = ipo0;
					m_aCell[ncell+3].size += 1;
				}
				}
				}
				else{
				assert( m_aVec[ipo0].second.x <= x_max0 );
				if( m_aVec[ipo0].second.y <= y_cent ){
				assert( m_aVec[ipo0].second.y >= y_min0 );
				if( m_aVec[ipo0].second.z <= z_cent ){
					assert( m_aVec[ipo0].second.z >= z_min0 );
					const unsigned int iloc1 = m_aCell[ncell+4].size;
					m_aCell[ncell+4].data[iloc1] = ipo0;
					m_aCell[ncell+4].size += 1;
				}
				else{
					assert( m_aVec[ipo0].second.z <= z_max0 );
					const unsigned int iloc1 = m_aCell[ncell+5].size;
					m_aCell[ncell+5].data[iloc1] = ipo0;
					m_aCell[ncell+5].size += 1;
				}
				}
				else{
				assert( m_aVec[ipo0].second.y <= y_max0 );
				if( m_aVec[ipo0].second.z <= z_cent ){
					assert( m_aVec[ipo0].second.z >= z_min0 );
					const unsigned int iloc1 = m_aCell[ncell+6].size;
					m_aCell[ncell+6].data[iloc1] = ipo0;
					m_aCell[ncell+6].size += 1;
				}
				else{
					assert( m_aVec[ipo0].second.z <= z_max0 );
					const unsigned int iloc1 = m_aCell[ncell+7].size;
					m_aCell[ncell+7].data[iloc1] = ipo0;
					m_aCell[ncell+7].size += 1;
				}
				}
				}
				m_aCell[icell0].data[i] = ncell+i;
				m_aCell[icell0].size = -1;
			}
		}
		else{	// このセルは８つに分けられているので下のセルを探す
			const double x_cent=(x_max0+x_min0)*0.5, y_cent=(y_max0+y_min0)*0.5, z_cent=(z_max0+z_min0)*0.5;
			int icell1;
			if( VecIns.x <= x_cent ){
			if( VecIns.y <= y_cent ){
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[0];
				x_max0 = x_cent; y_max0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[1];
				x_max0 = x_cent; y_max0 = y_cent; z_min0 = z_cent;
			}
			}
			else{
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[2];
				x_max0 = x_cent; y_min0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[3];
				x_max0 = x_cent; y_min0 = y_cent; z_min0 = z_cent;
			}
			}
			}
			else{
			if( VecIns.y <= y_cent ){
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[4];
				x_min0 = x_cent; y_max0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[5];
				x_min0 = x_cent; y_max0 = y_cent; z_min0 = z_cent;
			}
			}
			else{
			if( VecIns.z <= z_cent ){
				icell1 = m_aCell[icell0].data[6];
				x_min0 = x_cent; y_min0 = y_cent; z_max0 = z_cent;
			}
			else{
				icell1 = m_aCell[icell0].data[7];
				x_min0 = x_cent; y_min0 = y_cent; z_min0 = z_cent;
			}
			}
			}
			icell0 = icell1;
		}
	}	
//	assert( Check() );
	return -1;
}
