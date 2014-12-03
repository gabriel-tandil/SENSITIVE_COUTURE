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
// eqnsys.cpp : 抽象連立方程式オブジェクトのインターフェース
////////////////////////////////////////////////////////////////

#if defined(__VISUALC__)
	#pragma warning( disable : 4786 )
#endif

#include "delfem/eqnsys.h"

#include "delfem/femls/linearsystem_field.h"
#include "delfem/ls/preconditioner.h"

using namespace Fem::Eqn;
using namespace Fem::Field;
using namespace Fem::Ls;


void CEqnSystem::Clear(){
	m_gamma_newmark = 0.6;
	m_beta_newmark = 0.3025;
	m_dt = 0.1;
	this->m_aItrNormRes.clear();
	if( pLS != 0 ){ delete pLS; pLS=0; }
	if( pPrec != 0 ){ delete pPrec; pPrec=0; }
	m_is_cleared_value_ls = true;
	m_is_cleared_value_prec = true;
}

void CEqnSystem::ClearLinearSystem()
{
	if( pLS   != 0 ){ delete pLS;   pLS=0;   }
}

void CEqnSystem::ClearPreconditioner()
{
	if( pPrec != 0 ){ delete pPrec; pPrec=0; }
}
