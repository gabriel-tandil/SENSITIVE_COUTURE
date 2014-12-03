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
@brief ILU前処理する正方行列クラス(MatVec::CMatDiaFrac_BlkCrs)のインターフェイス
@author Nobuyuki Umetani
*/

#if !defined(MAT_DIA_FRAC_BLK_CRS_H)
#define MAT_DIA_FRAC_BLK_CRS_H

#include <vector>

#include "delfem/matvec/matdia_blkcrs.h"
#include "delfem/matvec/matprecond_blk.h"
#include "delfem/matvec/ordering_blk.h"

namespace MatVec{

class CMatFrac_BlkCrs;
/*! 
@brief ILU前処理する正方行列クラス
@ingroup MatVec
*/
class CMatDiaFrac_BlkCrs : public CMatDia_BlkCrs, public CPrecond_Blk
{
	friend class CMatFrac_BlkCrs;
public:
    CMatDiaFrac_BlkCrs() : m_ConditionFlag(-1), m_DiaInd(0), m_pRowLev(0){}
	CMatDiaFrac_BlkCrs(const unsigned int nblk_colrow, const unsigned int len_colrow);
    //! ILU(0)によるパターンの初期化
	CMatDiaFrac_BlkCrs(const CMatDia_BlkCrs& rhs);
    //! ILU( lev_fill )　によるパターンの初期化
    CMatDiaFrac_BlkCrs(const int lev_fill, const CMatDia_BlkCrs& rhs);
    //! オーダーつきのコンストラクタ
	CMatDiaFrac_BlkCrs(
        const int lev_fill, const CMatDia_BlkCrs& rhs, const COrdering_Blk& order);
	virtual ~CMatDiaFrac_BlkCrs();

    ////////////////////////////////
    //! ILU分解実行
	virtual bool SetValue(const CMatDia_BlkCrs& rhs);
	virtual bool SetValue(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order);
    //! 連立一次方程式を解く関数
	virtual bool Solve(CVector_Blk& vec) const{
		this->ForwardSubstitution(vec);
		this->BackwardSubstitution(vec);
		return true;
	}
    /*!
    @brief 連立一次方程式を解く関数
    @remark CPrecond_Blkクラスから継承することによってPreconditionerとして抽象化
    */
	virtual bool SolvePrecond(const CMatDia_BlkCrs& mat, CVector_Blk& vec) const{
		return this->Solve(vec);
	}

	////////////////////////////////////////////////////////////////
	// 管理クラス(CPreconditionerみたいの)から呼ばれるであろうルーティン


    ////////////////////////////////
    // 非ゼロパターンをコピーして初期化
	bool MakePattern_Initialize(const CMatDia_BlkCrs& rhs);
	bool MakePattern_Initialize(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order);
	// 非ゼロパターンを追加
	bool AddFracPtn(const int lev_fill);
    bool AddFracPtn(const int lev_fill, const std::vector<unsigned int>& aBlk);
	bool AddFracPtnLowUp(const int lev_fill, const CMatFrac_BlkCrs& mat_low, const CMatFrac_BlkCrs& mat_up);
    //! 値配列を確保して非ゼロパターンの最終化
	bool MakePatternFinalize();
    ////////////////////////////////

    // 行列に値を代入する
	bool SetValue_Initialize(const CMatDia_BlkCrs& rhs){
        assert( m_ConditionFlag != -1 );
        if( m_ConditionFlag != 2 ){ this->MakePatternFinalize(); }
        assert( m_ConditionFlag == 2 );
        CMatDia_BlkCrs::SetZero();
		CMatDia_BlkCrs::SetValue(rhs,true);
/*
        {
          for(unsigned int inode=0;inode<this->NBlkMatCol();inode++){
                std::cout << inode << "-->" << std::endl;
                const unsigned int len = this->LenBlkCol();
                std::cout << "dia : ";
                for(unsigned int i=0;i<len*len;i++){
                    std::cout << m_valDia_Blk[inode*len*len+i] << "    ";
                }
                std::cout << std::endl;
		        for(unsigned int ijcrs=m_colInd_Blk[inode];ijcrs<m_colInd_Blk[inode+1];ijcrs++){
			        const int jblk0 = m_rowPtr_Blk[ijcrs];
			        std::cout << "   " << jblk0 << "  ";
                    for(unsigned int i=0;i<len*len;i++){
                        std::cout << m_valCrs_Blk[ijcrs*len*len+i] << " ";
                    }
                    std::cout << std::endl;
		        }
	        }
        }*/
		return true;
	}
	bool SetValue_Initialize(const CMatDia_BlkCrs& rhs, const COrdering_Blk& order){
        assert( m_ConditionFlag != -1 );
        if( m_ConditionFlag != 2 ){ this->MakePatternFinalize(); }
        assert( m_ConditionFlag == 2 );
        CMatDia_BlkCrs::SetZero();
		CMatDia_BlkCrs::SetValue(rhs,order);
		return true;
	}
	// ILU分解を実行する関数
	bool DoILUDecomp();
	bool DoILUDecompLowUp( 
		const CMatFrac_BlkCrs& mat_low, 
		const CMatFrac_BlkCrs& mat_up );

    ////////////////////////////////
    // 行列を解く
    //! 前進代入
	bool ForwardSubstitution( CVector_Blk& vec ) const;
    //! 後退代入
	bool BackwardSubstitution( CVector_Blk& vec ) const;

    ////////////////////////////////

    const double* GetValCrsPtr(unsigned int icrs) const{
        assert( icrs < NCrs() );
        assert( m_valCrs_Blk );
        if( m_ValPtr ){  // 対角の値インデックスがとりあえず判定基準
            assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
            unsigned int ipos = m_ValPtr[icrs];
            return &m_valCrs_Blk[ipos];
        }
        assert( LenBlkCol() >= 0 && LenBlkRow() >= 0 );
        return &m_valCrs_Blk[icrs*LenBlkCol()*LenBlkRow()];
    }
    const double* GetValDiaPtr(unsigned int iblk) const{
        assert( iblk < NBlkMatCol() );
        assert( m_valDia_Blk );
        if( m_DiaValPtr ){  // 対角の値インデックスがとりあえず判定基準
            assert( LenBlkCol() == -1 && LenBlkRow() == -1 );
            unsigned int ipos = m_DiaValPtr[iblk];
            return &m_valDia_Blk[ipos];
        }
        assert( LenBlkCol() >= 0 && LenBlkRow() >= 0 );
        return &m_valDia_Blk[iblk*LenBlkCol()*LenBlkRow()];
    }
	
private:
protected:
    //! パターンとそのレベルを表すときに使うクラス
	class CRowLev{
	public:
		CRowLev():row(0), lev(0){}
		CRowLev(unsigned int row, int lev) : row(row), lev(lev){}
		bool operator < (const CRowLev& rhs) const{
			return row < rhs.row;			
		}
	public:
		unsigned int row;
		int lev;
	};
    //! パターンを作る時に使うクラス
	struct SRowLevNext{
		unsigned int row;
		unsigned int lev;
		int next;
	};
    //! データは昇順に並ぶはず
    class CListRLN{
    public:
        CListRLN(unsigned int n) : max_size(n){ list.reserve(max_size); ipos_entry = -1; }
    public:
        void Initialize( const std::vector<CRowLev>& aRV, unsigned int is, unsigned int ie ){
            const unsigned int n = ie - is;
            assert( n < max_size );
            list.resize( n );
            for(unsigned int i=0;i<n;i++){
                list[i].row  = aRV[i+is].row;
                assert( aRV[i+is].row < max_size );
                list[i].lev  = aRV[i+is].lev;
                list[i].next = i+1;
            }
            if( n > 0 ){
                list[n-1].next = -1;
                ipos_entry = 0;
            }
            else{
                ipos_entry = -1;
            }
        }
        // 要素(jblk0,lev)を加える．加えた要素のindが返る．探索はjnz_curの次から始める．
        bool Insert(unsigned int jblk0, unsigned int lev, int& jnz_cur )
        {
            assert( jblk0 <= max_size );
            if( list.size() == 0 ){ // 何も無いところに付け加える場合
                assert( ipos_entry == -1 );
                list.resize(1);
                list[0].row = jblk0;
                list[0].lev = lev;
                list[0].next = -1;
                ipos_entry = 0;
                return true;
            }
            if( jblk0 < list[ipos_entry].row ){
			    const unsigned  int inz_last = list.size();
			    list.resize( list.size()+1 );
			    list[inz_last].row = jblk0;
			    list[inz_last].lev = lev;
			    list[inz_last].next = ipos_entry;
			    ipos_entry = inz_last;
                jnz_cur = ipos_entry;
                return true;    
            }
            if( jblk0 == list[ipos_entry].row ){ return false; }
            if( jnz_cur < 0 ){  // マイナスの値が入ったら全探索するということで
                jnz_cur = ipos_entry; 
                assert( jnz_cur >= 0 && jnz_cur < (int)list.size() );
            }
            assert( jnz_cur < (int)list.size() );
            if( jblk0 == list[jnz_cur].row ){ return false; }
            if( list[jnz_cur].row > jblk0 ){ jnz_cur = ipos_entry; }
            for(;;){
                int jnz_nex = list[jnz_cur].next;
                if( jnz_nex == -1 ) break;
                assert( jnz_nex >= 0 && jnz_nex < (int)list.size() );
                assert( list[jnz_cur].row < list[jnz_nex].row );
				if( list[jnz_nex].row  > jblk0 ){ break; }   // 飛び越えた
				if( list[jnz_nex].row == jblk0 ){ return false; }  // 既に要素が入っている
				assert( list[jnz_nex].row < jblk0 );
				jnz_cur = jnz_nex;
            }
            assert( jnz_cur < (int)list.size() );
			assert( list[jnz_cur].row < jblk0 );
			assert( list[jnz_cur].next == -1 || jblk0 < list[ list[jnz_cur].next ].row );
			const unsigned  int inz_last = list.size();
			list.resize( list.size()+1 );
			list[inz_last].row = jblk0;
			list[inz_last].lev = lev;
			list[inz_last].next = list[jnz_cur].next;
			list[jnz_cur].next = inz_last;
			jnz_cur = inz_last;
            return true;
        }
        void Finalize(unsigned int iblk, unsigned int* m_colInd_Blk,
            std::vector<CRowLev>& aRowLev, 
            unsigned int& ncrs, unsigned int* m_DiaInd)
        {            
            assert( list.size() <= max_size );
			if( m_colInd_Blk[iblk] + list.size() > aRowLev.max_size() ){
//				std::cout << "		overflow and memory reallocate in ilu frac" << std::endl;
			}
			aRowLev.resize( m_colInd_Blk[iblk] + list.size() );
			unsigned int icrs0 = m_colInd_Blk[iblk];
            if( list.size() > 0 ){
                unsigned int ipos_cur = ipos_entry;
                for(;;){
                    assert( ipos_cur < list.size() );
				    const unsigned int jblk = list[ipos_cur].row;
				    const unsigned int jlev = list[ipos_cur].lev;
				    assert( jblk<max_size );
				    assert( jblk != iblk );
				    aRowLev[icrs0].row = jblk;
				    aRowLev[icrs0].lev = jlev;
				    icrs0++;
                    if( list[ipos_cur].next == -1 ) break;
                    assert( list[ipos_cur].row < list[list[ipos_cur].next].row );   // 昇順に並んでるはず
                    ipos_cur = list[ipos_cur].next;
                }
            }
			m_colInd_Blk[iblk+1] = icrs0;
			ncrs = icrs0;
			m_DiaInd[iblk] = icrs0;
			for(unsigned  int ijcrs=m_colInd_Blk[iblk];ijcrs<m_colInd_Blk[iblk+1];ijcrs++){
				const unsigned int jblk0 = aRowLev[ijcrs].row;
				if( jblk0 > iblk ){
					m_DiaInd[iblk] = ijcrs;
					break;
				}
			}
        }
    public:
        std::vector<SRowLevNext> list;
        const unsigned int max_size;
        int ipos_entry;
    };    
protected:
    /*!
    -1 for no size
    0 for no pattern
    1 for during pattern
    2 for after pattern
    */
	int m_ConditionFlag;

	unsigned int* m_DiaInd;
    std::vector<CRowLev>* m_pRowLev; //! have value if condition flag=2
};

}	// end namespace 'Ls'

#endif // !defined(MAT_DIA_CRS_FRAC_H)
