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
 @brief ２次元ベクトルクラス(Com::CVector2D)の実装
 @author Nobuyuki Umetani
 */

#if !defined(FIELD_VALUE_SETTER_H)
#define FIELD_VALUE_SETTER_H

#include "delfem/eval.h"
#include "delfem/field.h"
#include "delfem/field_world.h"

namespace Fem{
namespace Field{
  
//! set constant value to the field
bool SetFieldValue_Constant
  (unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
   Fem::Field::CFieldWorld& world,
   double val);
  
//! set mathematical expression to the field
bool SetFieldValue_MathExp
  (unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
   Fem::Field::CFieldWorld& world,                          
   std::string str_exp, double t=0);
  
//! set random field to the field
void SetFieldValue_Random
  (unsigned int id_field_to, unsigned int idofns, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
   Fem::Field::CFieldWorld& world,
   double ave, double range);
  
//! copy value to the field
void SetFieldValue_Copy
  (unsigned int id_field_to, Fem::Field::FIELD_DERIVATION_TYPE fdt, 
   Fem::Field::CFieldWorld& world,
   unsigned int id_field_from);
  
//! set gradient value to the field
bool SetFieldValue_Gradient
  (unsigned int id_field_to, Fem::Field::CFieldWorld& world,
   unsigned int id_field_from);

  
class CFieldValueSetter
{	//! set saved value to the entire field
public:
  CFieldValueSetter(){ id_field_ = 0; id_field_gradient_ = 0; }
  CFieldValueSetter(unsigned int id_field, Fem::Field::CFieldWorld& world);
  void Clear(){
    id_field_ = 0;
    id_field_gradient_ = 0;    
    aValueFieldDof_.clear();
  }
  void SetMathExp
  (const std::string& math_exp, unsigned int idof, Fem::Field::FIELD_DERIVATION_TYPE fdt,
   Fem::Field::CFieldWorld& world);
  void SetConstant
  (double val, unsigned int idof, Fem::Field::FIELD_DERIVATION_TYPE fdt,
   Fem::Field::CFieldWorld& world);  
  void SetGradient
  (unsigned int id_field_from, Fem::Field::CFieldWorld& world);
  
	bool ExecuteValue(double time, Fem::Field::CFieldWorld& world);
  
private:
	class CValueFieldDof{
	public:
		CValueFieldDof(double val){ this->SetValue(val); }
		CValueFieldDof(const std::string str){ 	this->SetValue(str); }
		CValueFieldDof(){ itype=0; }
		////////////////
		void SetValue(std::string str){
			itype = 2;
			math_exp = str;
		}
		void SetValue(double val){
			itype =1;
			this->val = val;
		}
		bool IsTimeDependent() const{
			if( itype != 2 ) return false;
			Fem::Field::CEval eval;
			eval.SetKey("t",0);
			eval.SetExp(math_exp);
			return eval.IsKeyUsed("t");
		}
    bool GetValue(double cur_t, double& value) const;
		const std::string GetString() const{
			if( itype == 1 ){
				char buff[16];
				sprintf(buff,"%lf",val);
				return std::string(buff);
			}
			else{
				return math_exp;
			}
		}
	public:
    // 0 : not set
    // 1 : const value 
    // 2 : math_exp
		int itype; 
		double val;
		std::string math_exp;
	};  
private:
  unsigned int id_field_;
  std::vector<CValueFieldDof> aValueFieldDof_;
  
  unsigned int id_field_gradient_;
};
  
}
}

#endif
