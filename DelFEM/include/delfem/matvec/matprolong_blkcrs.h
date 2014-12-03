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


// Prolongation_Crs.h: CProlongation_Crs クラスのインターフェイス
//
//////////////////////////////////////////////////////////////////////

#if !defined(PROLONGATION_CRS_H)
#define PROLONGATION_CRS_H

#include "delfem/matvec/mat_blkcrs.h"

namespace field{
	class CElemAry_SquarePow2;
}

namespace MatVec{

class CMatDia_BlkCrs;

class CMatProlong_BlkCrs : public CMat_BlkCrs  
{
private:
	struct SBucket{
		unsigned int* data;
		unsigned int size;
	};
public:
	CMatProlong_BlkCrs(const field::CElemAry_SquarePow2& ea_c, const field::CElemAry_SquarePow2& ea_f);
	CMatProlong_BlkCrs(const CMatDia_BlkCrs& mat,const double& theta);
	virtual ~CMatProlong_BlkCrs();

	bool MakeValueAMG(const CMatDia_BlkCrs& mat,const double& theta);
	bool ScatterCoaseValue(unsigned int* distination, const unsigned int* coarse_value);
	bool MaskCoarse(unsigned int* distination);
private:
	unsigned int* m_MapperCF;
};

}	// end of namespace ls

#endif // !defined(PROLONGATION_CRS_H)
