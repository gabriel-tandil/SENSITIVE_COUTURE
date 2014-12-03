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
@brief interface of block vector class (MatVec::CVector_Blk)
@author Nobuyuki Umetani
*/

#if !defined(VECTOR_BLK_H)
#define VECTOR_BLK_H

#include <assert.h>
#include <vector>


namespace MatVec{

class CMat_BlkCrs;
class CMatDia_BlkCrs;

/*! 
@brief real value block vector class
@ingroup MatVec
*/
class CVector_Blk  
{
	friend class CMat_BlkCrs;
	friend class CMatDia_BlkCrs;
	friend class CMatDiaFrac_BlkCrs;
	friend class CMatDiaInv_BlkDia;
  friend double operator*(const CVector_Blk& lhs, const CVector_Blk& rhs);	//!< Dot Product
public:
	/*!
	@brief コンストラクタ
	@param[in] iblkveclen ブロックの数
	@param[in] iblklen １つのブロックのサイズ
	*/
  CVector_Blk(unsigned int nblk, unsigned int len) : m_nBlk(nblk), m_Len(len){
    m_DofPtr = 0;
		m_Value = new double [m_nBlk*m_Len];
	}
  CVector_Blk(unsigned int nblk, const std::vector<unsigned int>& aLen) : m_nBlk(nblk), m_Len(-1){
    m_DofPtr = 0;
    this->Initialize(nblk,aLen);
	}
  CVector_Blk(const CVector_Blk& vec){
    m_DofPtr = 0;
    m_nBlk = vec.NBlk();
    m_Len = vec.Len();
    m_Value = new double [m_nBlk*m_Len];
    for(unsigned int iblk=0;iblk<m_nBlk;iblk++){
    for(unsigned int ilen=0;ilen<m_Len;ilen++){      
      m_Value[iblk*m_Len+ilen] = vec.GetValue(iblk, ilen);
    }
    }
  }
	CVector_Blk() : m_nBlk(0), m_Len(0), m_Value(0), m_DofPtr(0){}	//!< default constructor
	virtual ~CVector_Blk(){ if( m_Value!=0) delete[] m_Value; }	//!< destructor

	bool Initialize(unsigned int nblk, unsigned int len ){
		m_nBlk = nblk;
		m_Len = len;
    if( m_DofPtr != 0 ){ delete[] m_DofPtr; m_DofPtr = 0; }
		if( m_Value == 0 ){ delete[] m_Value; }
		m_Value = new double [m_nBlk*m_Len];
		return true;
	}
  
  bool Initialize(unsigned int nblk, const std::vector<unsigned int>& alen ){
		m_nBlk = nblk;
		m_Len = -1;
    ////////////////
    if( m_DofPtr != 0 ){ delete[] m_DofPtr; m_DofPtr = 0; }
    m_DofPtr = new unsigned int [nblk+1];
    m_DofPtr[0] = 0;
    for(unsigned int iblk=0;iblk<nblk;iblk++){
      m_DofPtr[iblk+1] = m_DofPtr[iblk] + alen[iblk];
    }
    ////////////////
		if( m_Value == 0 ){ delete[] m_Value; }
    const unsigned int ni = m_DofPtr[nblk];
		m_Value = new double [ni];
		return true;
	}

	// @{
	CVector_Blk& operator*=(double d0);	//!< Scaler Product
	CVector_Blk& operator=(const CVector_Blk& rhs);	//!< Substitue Vector
	CVector_Blk& operator+=(const CVector_Blk& rhs); //!< Add 
	// @}


	////////////////////////////////
	// Menber function

	/*
	@brief ベクトルのスカラー倍の足し算
	{this} += alhpa*{rhs}
	*/
	CVector_Blk& AXPY(const double& alpha, const CVector_Blk& rhs);

	void SetVectorZero();	//!< Set 0 to Value
	double GetSquaredVectorNorm() const;	//!< ベクトルの２乗ノルムの２乗を計算する

	inline unsigned int NBlk() const { return m_nBlk; }	//!< Size of Bocks Vector have
	inline int Len() const { return m_Len; }	//!< Size of One Block(もしFlexサイズなら-1を返す)
  //! Size of One Block
	unsigned int Len(unsigned int iblk) const { 	
    assert( iblk < this->NBlk() );
    if( m_Len >= 0 ){ return m_Len; }
    return m_DofPtr[iblk+1]-m_DofPtr[iblk];
  }  
  
	/*!
	@brief 値を取得
	@param [in] iblk ブロックのインデックス
	@param [in] idofblk ブロックの中の自由度番号
	@retval 値
	*/
	inline const double& GetValue(const unsigned int iblk, const unsigned int idofblk) const{
		assert( iblk<m_nBlk );
    if( m_Len >= 0 ){
      assert( (int)idofblk<m_Len );
      return m_Value[iblk*m_Len+idofblk];
    }
    assert( m_Len == -1 );
    assert( idofblk < m_DofPtr[iblk+1]-m_DofPtr[iblk] );
    assert( m_DofPtr != 0 );
    return m_Value[ m_DofPtr[iblk]+idofblk ];    
	}
	/*!
	@brief 値をセット
	@param [in] iblk ブロックのインデックス
	@param [in] idofblk ブロックの中の自由度番号
	@param [in] val 値
	*/
	inline void SetValue(const unsigned int iblk, const unsigned int idofblk, const double val){
		assert( iblk<m_nBlk );
    if( m_Len >= 0 ){
      assert( (int)idofblk<m_Len );
      m_Value[iblk*m_Len+idofblk] = val;
      return;
    }
    assert( m_Len == -1 );
    assert( idofblk < m_DofPtr[iblk+1]-m_DofPtr[iblk] );
    assert( m_DofPtr != 0 );
    m_Value[ m_DofPtr[iblk]+idofblk ] = val;    
	}
    /*!
	@brief 特定の要素に値を足し合わせる
	@param [in] iblk ブロックのインデックス
	@param [in] idofblk ブロックの中の自由度番号
	@param [in] val 足し合わせる値
	*/
	inline void AddValue(const unsigned int iblk, const unsigned int idofblk, const double val){
		assert( iblk<m_nBlk );
    if( m_Len >= 0 ){
      assert( (int)idofblk<m_Len );
      m_Value[iblk*m_Len+idofblk] += val;
      return;
    }
    assert( m_Len == -1 );
    assert( idofblk < m_DofPtr[iblk+1]-m_DofPtr[iblk] );
    assert( m_DofPtr != 0 );
    m_Value[ m_DofPtr[iblk]+idofblk ] += val;    
	}
  double* GetValuePtr(unsigned int iblk){
		assert( iblk<m_nBlk );
    if( m_Len >= 0 ){
      return m_Value + iblk*m_Len;
    }
    assert( m_Len == -1 );
    assert( m_DofPtr != 0 );
    return m_Value + m_DofPtr[iblk];
  }
  const double* GetValuePtr(unsigned int iblk) const{
		assert( iblk<m_nBlk );
    if( m_Len >= 0 ){
      return m_Value + iblk*m_Len;
    }
    assert( m_Len == -1 );
    assert( m_DofPtr != 0 );
    return m_Value + m_DofPtr[iblk];
  }  
private:
  unsigned int GetTotalDofSize() const {
    if( m_Len == -1 ){
      assert( m_DofPtr != 0 );
      return m_DofPtr[m_nBlk];
    }
    assert( m_DofPtr == 0 );
    return m_nBlk*m_Len;      
  }
private:
	unsigned int m_nBlk;	//!< number of block
  int m_Len;  //!< degree of freedom par block(-1 if size is flex)
	double* m_Value;			//!< value array
  unsigned int* m_DofPtr; //!< 0 if blk size is fixed, nonzero if size is flex
};

}	// end namespace 'Ls'

#endif // VEC_H
