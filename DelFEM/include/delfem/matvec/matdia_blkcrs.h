/*
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
@brief 正方行列クラス(MatVec::CMatDia_BlkCrs)のインターフェイス
*/

#if !defined(MATDIA_CRS_H)
#define MATDIA_CRS_H

#include "delfem/matvec/mat_blkcrs.h"	// このクラスを継承している

namespace Com{
    class CIndexedArray;
}

namespace MatVec{

class CVector_Blk;
class CBCFlag;
class COrdering_Blk;

/*!
@brief square crs matrix class
@ingroup MatVec

データ構造:対角成分を持つブロックCRSクラス
*/
class CMatDia_BlkCrs  : public CMat_BlkCrs
{
public:
	/*!
	@brief サイズを引数とするコンストラクタ
	@param[in] nblk_colrow ブロック数
	@param[in] len_colrow ブロックのサイズ
	*/
	CMatDia_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow);
	//! デフォルトコンストラクタ
	CMatDia_BlkCrs();
	//! デストラクタ
	virtual ~CMatDia_BlkCrs();

    
    // Mat_BlkCrsクラスの隠蔽
    virtual bool Initialize(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
                            unsigned int nblk_row, const std::vector<unsigned int>& alen_row);
    bool Initialize(unsigned int nblk, const std::vector<unsigned int>& alen){
        return this->Initialize(nblk,alen, nblk,alen);
    }
    // Mat_BlkCrsクラスの隠蔽
    virtual bool Initialize(unsigned int nblk_col, unsigned int len_col,
                            unsigned int nblk_row, unsigned int len_row );
    bool Initialize(unsigned int nblk, unsigned int len){
        return this->Initialize(nblk,len, nblk,len);
    }

	//! パターンを全て消去　RowPtr,Valはメモリ解放
	virtual bool DeletePattern();
    bool AddPattern(const Com::CIndexedArray& crs);
	bool AddPattern(const CMatDia_BlkCrs& rhs, const bool isnt_trans);
	bool AddPattern(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order);
	bool AddPattern(const CMat_BlkCrs& m1, const CMatDia_BlkCrs& m2, const CMat_BlkCrs& m3);

	virtual bool SetValue(const CMatDia_BlkCrs& rhs, const bool isnt_trans);
	virtual bool SetValue(const CMat_BlkCrs& m1, const CMatDia_BlkCrs& m2, const CMat_BlkCrs& m3);	// := m1*m2*m3
	virtual bool SetValue(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order);

	//! 行列に０をセットするための関数(親クラスの隠蔽)
	virtual bool SetZero();
	//! マージするための関数(親クラスの隠蔽)
	virtual bool Mearge(
		unsigned int nblkel_col, const unsigned int* blkel_col,
		unsigned int nblkel_row, const unsigned int* blkel_row,
		unsigned int blksize, const double* emat);
	/*!
	@brief 行列ベクトル積(親クラスの隠蔽)
	{lhs} = beta*{lhs} + alpha*[A]{rhs}
	*/
	virtual bool MatVec(double alpha, const CVector_Blk& rhs, double beta, CVector_Blk& lhs) const;

	//! bc_flagが１の自由度の行と列を０に設定，但し対角成分は１を設定
	bool SetBoundaryCondition(const CBCFlag& bc_flag);

    ////////////////////////////////////////////////////////////////
    // 前処理用のクラス

	const double* GetPtrValDia(const unsigned int ipoin) const {
        assert( m_DiaValPtr == 0 );
		unsigned int blksize = this->LenBlkCol()*this->LenBlkRow();
		return &m_valDia_Blk[ipoin*blksize];
	}
	double* GetPtrValDia(const unsigned int ipoin){
        assert( m_DiaValPtr == 0 );
		unsigned int blksize = this->LenBlkCol()*this->LenBlkRow();
		return &m_valDia_Blk[ipoin*blksize];
	}

protected:
	double* m_valDia_Blk;	// 対角ブロック成分
    
  // Flex時に定義される値
  unsigned int* m_DiaValPtr;
};

}	// end namespace 'Ls'

#endif // MATDIA_CRS_H
