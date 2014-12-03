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
@brief 数式評価クラス(Fem::Field::CEval)のインターフェイス
@author Nobuyuki Umetani
*/

#if !defined(EVAL_H)
#define EVAL_H

#include <iostream>
#include <string>
#include <vector>

namespace Fem{
namespace Field{

class CCmd;
/*!
@brief 数式評価クラス

文字列で表された数式を評価します．
*/
class CEval
{
public:
	class CKey
	{	// 数式の要素の値を保持するクラス
	public:
		CKey(const std::string& name, double val)
			: m_Name(name), m_Val(val){}
		std::string m_Name;
		std::vector<unsigned int> m_aiCmd;
		double m_Val;
	};
public:
	CEval(){ m_is_valid = false; }
	CEval(const std::string& exp){ this->SetExp(exp); }
	~CEval();

	////////////////
	bool SetExp(const std::string& key_name);
	void SetKey(const std::string& key_name, double val);
	bool IsKeyUsed(const std::string& key_name){
		for(unsigned int ikey=0;ikey<m_aKey.size();ikey++){
			if( m_aKey[ikey].m_Name == key_name ){
				if( !m_aKey[ikey].m_aiCmd.empty() ) return true;
				return false; 
			}
		}
		return false;
	}
	double Calc() const;
private:
	struct SExpCompo{
		std::string sOpe;
		int iOpeType;
		int iOpe;
	};
private:
	static bool MakeRPN(unsigned int icur, std::vector<SExpCompo>& exp_node_vec);
	static void RemoveExpressionSpaceNewline(std::string& exp);
	static void RemoveExpressionBracket(std::string& exp);
	static int GetLowestPriorityOperator(int& ibegin, int& iend, 
		int& iOpeType, int& iOpe, const std::string& exp);
	bool MakeCmdAry(std::vector<CCmd*>& cmd_vec, const std::vector<SExpCompo>& exp_node_vec);
private:
	bool m_is_valid;
	std::string m_sExp;
	std::vector<CCmd*> m_apCmd;	// 逆ポーランド記法された演算子や値の列
	std::vector<CKey> m_aKey;	// 文字列の名前と値が，どのIndexのコマンドに格納されているか
};

}
}

#endif	// !defind EVAL_H
