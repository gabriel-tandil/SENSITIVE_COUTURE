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
@brief Interface of class (MatVec::CMat_BlkCrs)
@author Nobuyuki Umetani
*/

#if !defined(MAT_CRS_H)
#define MAT_CRS_H

#include <assert.h>
#include <vector>

namespace Com{
    class CIndexedArray;
}
	
namespace MatVec{

class CVector_Blk;
class COrdering_Blk;
class CBCFlag;

/*! 
@brief ブロックCRS構造の行列クラス
@ingroup MatVec
*/
class CMat_BlkCrs  
{
public:
	/*!
	@brief サイズを引数としたコンストラクタ
	@param[in] nblk_col 列のブロック数
	@param[in] len_col  ブロックの列の長さ
	@param[in] nblk_row 行のブロック数
	@param[in] len_row  ブロックの行の長さ
	*/
	CMat_BlkCrs(unsigned int nblk_col, unsigned int len_col, 
		        unsigned int nblk_row, unsigned int len_row );
    
  CMat_BlkCrs(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
              unsigned int nblk_row, const std::vector<unsigned int>& alen_row );

	/*!
	@brief 行列をコピーするためのコンストラクタ
	@param[in] rhs 引数行列
	@param[in] is_value 値までコピーするかどうか
	@param[in] isnt_trans 転置としてコピーしないか
	*/
	CMat_BlkCrs(const CMat_BlkCrs& rhs, bool is_value, bool isnt_trans);

	CMat_BlkCrs(); 
	virtual ~CMat_BlkCrs();

	virtual bool Initialize(unsigned int nblk_col, unsigned int len_col, 
		                    unsigned int nblk_row, unsigned int len_row );
  // Diaで隠蔽するためにvirtualにしておく
  virtual bool Initialize(unsigned int nblk_col, const std::vector<unsigned int>& alen_col, 
                          unsigned int nblk_row, const std::vector<unsigned int>& alen_row );

	inline unsigned int NBlkMatCol() const { return m_nblk_MatCol; }	//!< 行に幾つブロックがあるか
	inline unsigned int NBlkMatRow() const { return m_nblk_MatRow; }	//!< 列に幾つブロックがあるか
	
	//! ブロックの行のサイズ(縦の長さ)を取得.もしもFixされて無い場合は-1を返す
	inline int LenBlkCol() const { return m_len_BlkCol; }
	//! ブロックの行のサイズ(縦の長さ)を取得.
	inline unsigned int LenBlkCol(unsigned int iblk) const { 
    assert( iblk < m_nblk_MatCol );
    if( m_len_BlkCol == -1 ){
      assert( m_DofPtrCol != 0 );
      return m_DofPtrCol[iblk+1]-m_DofPtrCol[iblk];
    }
    return m_len_BlkCol; 
  }
	//! ブロックの列のサイズ(横の長さ)を取得.もしもFixされて無い場合は-1を返す
	inline int LenBlkRow() const { return m_len_BlkRow; }
	//! ブロックの行のサイズ(横の長さ)を取得.
	inline unsigned int LenBlkRow(unsigned int iblk) const { 
    assert( iblk < m_nblk_MatRow );
    if( m_len_BlkRow == -1 ){
      assert( m_DofPtrRow != 0 );
      return m_DofPtrRow[iblk+1]-m_DofPtrRow[iblk];
    }
    return m_len_BlkRow; 
  }

	//! 値に０を設定．データ領域が確保されていなければ，確保する
	virtual bool SetZero();	
	//! パターンを全て消去　RowPtr,Valはメモリ解放
	virtual bool DeletePattern();	

	virtual void FillPattern();
  virtual bool AddPattern(const Com::CIndexedArray& crs);	//!< 非ゼロパターンを追加する
	virtual bool AddPattern(const CMat_BlkCrs& rhs, const bool isnt_trans);	//!< 非ゼロパターンを追加する
	virtual bool AddPattern(const CMat_BlkCrs& rhs, 
                          const COrdering_Blk& order_col, const COrdering_Blk& order_row);	//!< 非ゼロパターンを追加する

	// ここは派生クラスに任せるつもり
	virtual bool SetPatternBoundary(const CMat_BlkCrs& rhs, const CBCFlag& bc_flag_col, const CBCFlag& bc_flag_row);
	virtual bool SetPatternDia(const CMat_BlkCrs& rhs);
	
	virtual bool SetValue(const CMat_BlkCrs& rhs, const bool isnt_trans);
	virtual bool SetValue(const CMat_BlkCrs& rhs, 
                        const COrdering_Blk& order_col, const COrdering_Blk& order_row);
	
	//! 要素剛性行列をマージする
	virtual bool Mearge
  (unsigned int nblkel_col, const unsigned int* blkel_col,
   unsigned int nblkel_row, const unsigned int* blkel_row,
   unsigned int blksize, const double* emat);
  //! 作業用配列の領域を解放する
  void DeleteMargeTmpBuffer(){
    if( m_marge_tmp_buffer == 0 ) return;
    delete[]  m_marge_tmp_buffer;
    m_marge_tmp_buffer = 0;
  }

	//! 行列ベクトル積
	virtual bool MatVec(double alpha, const CVector_Blk& x, double beta, CVector_Blk& b, const bool isnt_trans) const;

	//! bc_flagが１の行の要素を０にする
	bool SetBoundaryCondition_Row(const CBCFlag& bc_flag);
	//! bc_flagが１の列の要素を０にする
	bool SetBoundaryCondition_Colum(const CBCFlag& bc_flag);
	//! bc_flagが０の要素を０にする（全体剛性行列保存による高速化法）
	bool SetBoundaryConditionInverse_Colum(const CBCFlag& bc_flag);
	
	////////////////////////////////////////////////
	// Crs Original Function ( Non-virtual )

	//! CRSのサイズ取得
	const unsigned int NCrs() const { return m_ncrs_Blk; }
	const unsigned int* GetPtrIndPSuP(const unsigned int ipoin, unsigned int& npsup) const
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_rowPtr_Blk[icrs];
	}
	const double* GetPtrValPSuP(const unsigned int ipoin, unsigned int& npsup) const
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_valCrs_Blk[icrs*m_len_BlkCol*m_len_BlkRow];
	}
	double* GetPtrValPSuP(const unsigned int ipoin, unsigned int& npsup)
	{
		if( m_rowPtr_Blk == 0 ){ npsup = 0; return 0; }
		assert( ipoin < m_nblk_MatCol );
		npsup = m_colInd_Blk[ipoin+1]-m_colInd_Blk[ipoin];
		const unsigned int icrs = m_colInd_Blk[ipoin];
		return &m_valCrs_Blk[icrs*m_len_BlkCol*m_len_BlkRow];
	}

private:
	unsigned int m_nblk_MatCol;	//!< 行に幾つブロックがあるか
	unsigned int m_nblk_MatRow;	//!< 列に幾つブロックがあるか

	int m_len_BlkCol;	//!< ブロックの列サイズ
	int m_len_BlkRow;	//!< ブロックの行サイズ

protected:
	unsigned int  m_ncrs_Blk;	//!< CRSのサイズ
	unsigned int* m_colInd_Blk;	//!< CRSのColum Index
	unsigned int* m_rowPtr_Blk;	//!< CRSのRow Pointer

  int* m_marge_tmp_buffer;	//!< マージの時に必要な作業バッファ

  // Flexの時に定義される値
  unsigned int* m_DofPtrCol;
  unsigned int* m_DofPtrRow;
  unsigned int* m_ValPtr;

	double* m_valCrs_Blk;	//!< 行列の値
};

}	// end namespace 'Ls'

#endif 
