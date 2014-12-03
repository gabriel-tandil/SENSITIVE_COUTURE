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
@brief 抽象方程式オブジェクトのインターフェース
@author Nobuyuki Umetani
*/

#if !defined(EQN_SYS_H)
#define EQN_SYS_H

#include <cassert>

#include <vector>
#include <map>

#if defined(__VISUALC__)
#pragma warning( disable : 4786 )
#endif

namespace LsSol{
	class CPreconditioner;	// 前処理行列
}
namespace Fem{
namespace Ls{
	class CLinearSystem_Field;					// 連立一次方程式
	class CLinearSystem_Save;					// 連立一次方程式(剛性行列保存)
	class CLinearSystem_SaveDiaM_NewmarkBeta;	// 連立一次方程式(NewmarkBeta法で剛性行列保存)
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn{

/*! 
@brief 抽象連成方程式クラス
@ingroup FemEqnSystem
*/
class CEqnSystem
{
public:		
	//! デフォルト・コンストラクタ
	CEqnSystem() : m_gamma_newmark(0.6), m_beta_newmark(0.3025), m_dt(0.1), pLS(0), pPrec(0){}
	//! デストラクタ
	virtual ~CEqnSystem(){ this->Clear(); }
	virtual void Clear();

	//! 方程式を解く
	virtual bool Solve(Fem::Field::CFieldWorld& world) = 0;

	const std::vector< std::pair<unsigned int, double> >& GetAry_ItrNormRes() const{ return m_aItrNormRes; }

	////////////////////////////////
	// 固定境界条件を追加する

	//! @{
	//! 場(id_field)の要素配列を全て固定境界にする
	virtual bool         AddFixField(                    unsigned int id_field, Fem::Field::CFieldWorld& world, int idof = -1) = 0;
	//! 場(m_id_val)の要素配列(id_ea)を固定境界にする
	virtual unsigned int AddFixElemAry(                  unsigned int id_ea,    Fem::Field::CFieldWorld& world, int idof = -1) = 0;
	//! 場(m_id_val)の要素配列の配列(aIdEA)を固定境界にする
	virtual unsigned int AddFixElemAry( const std::vector<unsigned int>& aIdEA, Fem::Field::CFieldWorld& world, int idof = -1) = 0;
	//! 固定境界条件を削除する
	virtual bool ClearFixElemAry( unsigned int id_ea, Fem::Field::CFieldWorld& world ) = 0;
	//! 全ての固定境界条件を削除する
	virtual void ClearFixElemAry() = 0;
	//! @}

	/*! 
	@brief 時間積分についてのパラメータを設定する
	@param [in] dt 時間刻み
	@param [in] gamma Newmark法のgamma(省略すれば0.6にセット)
	@param [in] beta Newmark法のbeta(gammaがあれば省略可)
	*/
	void SetTimeIntegrationParameter(double dt, double gamma=0.6, double beta=-1.0){
		m_dt = dt;
		m_gamma_newmark = gamma;
		if( beta < 0.0 ){ m_beta_newmark = 0.25*(m_gamma_newmark+0.5)*(m_gamma_newmark+0.5); }
		else{             m_beta_newmark = beta; }
	}

	////////////////////////////////
	// 連立一次方程式クラスや前処理クラスの再評価，再構築指定関数

	//! @{
	virtual void ClearValueLinearSystemPreconditioner(){ 
		this->ClearValueLinearSystem();
		this->ClearPreconditioner();
	}
	virtual void ClearValueLinearSystem(){   m_is_cleared_value_ls   = true; }	// フラグを立てるとSolveの時に値が再評価される
	virtual void ClearValuePreconditioner(){ m_is_cleared_value_prec = true; }	// フラグを立てるとSolveの時に値が再評価される

	virtual void ClearLinearSystemPreconditioner(){
		this->ClearLinearSystem(); 
		this->ClearPreconditioner();
	}
	virtual void ClearPreconditioner();
	virtual void ClearLinearSystem();
	//! @}
protected:
	std::vector< std::pair<unsigned int, double> > m_aItrNormRes;
	////////////////
	double m_gamma_newmark, m_beta_newmark, m_dt;
	Fem::Ls::CLinearSystem_Field* pLS;	// 連立一次方程式クラス
	LsSol::CPreconditioner* pPrec;	// 前処理クラス
	bool m_is_cleared_value_ls;
	bool m_is_cleared_value_prec;
};

}	// end namespace Eqn
}	// end namespace Fem

#endif
