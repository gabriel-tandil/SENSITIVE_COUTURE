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

#if defined(__VISUALC__)
#pragma warning ( disable : 4786 ) 
#pragma warning ( disable : 4996 )
#endif
#define for if(0);else for

#if defined(_WIN32)
#  include <windows.h>
#if defined(__VISUALC__)
#  pragma comment (lib, "winmm.lib")      /* link with Windows MultiMedia lib */
#  pragma comment (lib, "opengl32.lib")  /* link with Microsoft OpenGL lib */
#  pragma comment (lib, "glu32.lib")     /* link with Microsoft OpenGL Utility lib */
#endif  /* _WIN32 */
#endif

#if defined(__APPLE__) && defined(__MACH__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

//#include <GL/glut.h>



#include <iostream>
#include "delfem/rigid/linearsystem_rigid.h"
#include "delfem/rigid/rigidbody.h"
#include "delfem/indexed_array.h"

static void CalcInvMat(double* a, const unsigned int& n, int& info )
{
	double tmp1;

	info = 0;
	unsigned int i,j,k;
	for(i=0;i<n;i++){
		if( fabs(a[i+i*n]) < 1.0e-30 ){
			info = 1;
			return;
		}
		if( a[i+i*n] < 0.0 ){
			info--;
		}
		tmp1 = 1.0 / a[i+i*n];
		a[i+i*n] = 1.0;
		for(k=0;k<n;k++){
			a[i+k*n] *= tmp1;
		}
		for(j=0;j<n;j++){
			if( j!=i ){
				tmp1 = a[j+i*n];
				a[j+i*n] = 0.0;
				for(k=0;k<n;k++){
					a[j+k*n] -= tmp1*a[i+k*n];
				}
			}
		}
	}
}

/*
Ls::CLinearSystem_RigidBody_Full::CLinearSystem_RigidBody_Full(
    const std::vector<Rigid::CRigidBody3D>& aRB,
    const std::vector<Rigid::CConstraint*>& aConst)
{
    nRB = aRB.size();
    nConst = aConst.size();
    aIndDof.resize( nRB + nConst + 1 );
    aIndDof[0] = 0;
    for(unsigned int irb=0;irb<nRB;irb++){
        aIndDof[irb+1] = aIndDof[irb] + aRB[irb].GetDOF();
    }
    for(unsigned int icst=0;icst<nConst;icst++){
        aIndDof[nRB+icst+1] = aIndDof[nRB+icst] + aConst[icst]->GetDOF();
    }
    ndof = aIndDof[nRB+nConst];
    mat = new double [ndof*ndof];
    update = new double [ndof*ndof];
    residual = new double [ndof*ndof];
}


void Ls::CLinearSystem_RigidBody_Full::Solve()
{
	{
		int info;
		CalcInvMat(mat,ndof,info);
//        std::cout << info << std::endl;
		double norm_del = 0;
		for(unsigned int i=0;i<ndof;i++){
			update[i] = 0;
			for(unsigned int j=0;j<ndof;j++){
				update[i] += mat[i*ndof+j]*residual[j];
		    }
			norm_del += update[i]*update[i];
//			std::cout << "Delta : " << i << " " << update[i] << std::endl;
		}
		std::cout << "norm del : " << sqrt(norm_del) << std::endl;
    }
}


bool Ls::CLinearSystem_RigidBody_Full::UpdateValueOfRigidSystem(
    std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
    double dt, double newmark_gamma, double newmark_beta, 
    bool is_first) const
{
    for(unsigned int irb=0;irb<aRB.size();irb++){
        aRB[irb].UpdateSolution( this->GetUpdate(irb,true), 
            dt, newmark_gamma, newmark_beta, 
            is_first);
    }
    for(unsigned int icst=0;icst<aConst.size();icst++){
        aConst[icst]->UpdateSolution( this->GetUpdate(icst,false), 
            dt, newmark_gamma, newmark_beta);
    }
    return true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


void Ls::CLinearSystem_RigidBody_CRS::SetRigidSystem(
    const std::vector<Rigid::CRigidBody3D>& aRB,
    const std::vector<Rigid::CConstraint*>& aConst)
{
    this->Clear();
    nRB = aRB.size();
    nConst = aConst.size();
    aIndDof.resize( nRB + nConst + 1 );
    aIndDof[0] = 0;
    for(unsigned int irb=0;irb<nRB;irb++){
        aIndDof[irb+1] = aIndDof[irb] + aRB[irb].GetDOF();
    }
    for(unsigned int icst=0;icst<nConst;icst++){
        aIndDof[nRB+icst+1] = aIndDof[nRB+icst] + aConst[icst]->GetDOF();
    }
    ndof = aIndDof[nRB+nConst];

    ColInd = new unsigned int [nRB+nConst+1];
    {
        ColInd[0] = 0;
        for(unsigned int i=0;i<nRB;i++){ ColInd[i+1]=1; }
        for(unsigned int i=nRB;i<nRB+nConst;i++){ ColInd[i+1]=0; }
        for(unsigned int icst=0;icst<nConst;icst++){
            const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
            for(unsigned int i=0;i<aIndRB.size();i++){
                const unsigned int irb0 = aIndRB[i];
                ColInd[irb0+1] += 1;
            }
            ColInd[icst+nRB+1] += aIndRB.size();
        }
        for(unsigned int i=0;i<nRB+nConst;i++){ 
            ColInd[i+1] = ColInd[i+1] + ColInd[i];
        }
        ncrs = ColInd[nRB+nConst];
        RowPtr = new unsigned int [ncrs];
        for(unsigned int irb=0;irb<nRB;irb++){
            const unsigned int icrs = ColInd[irb];
            RowPtr[icrs] = irb;
            ColInd[irb] += 1;
        }
        for(unsigned int icst=0;icst<nConst;icst++){
            const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
            for(unsigned int i=0;i<aIndRB.size();i++){
                const unsigned int irb0 = aIndRB[i];
                const unsigned int icrs0 = ColInd[icst+nRB];
                RowPtr[icrs0] = irb0;
                ColInd[icst+nRB] += 1;
                const unsigned int icrs1 = ColInd[irb0];
                RowPtr[icrs1] = icst+nRB;
                ColInd[irb0] += 1;
            }
        }
        for(int i=nRB+nConst;i>0;i--){ 
            ColInd[i] = ColInd[i-1];
        }
        ColInd[0] = 0;
        assert( ColInd[nRB+nConst] == ncrs );
        ValPtr = new unsigned int [ncrs+1];
        unsigned int val_pos = 0;
        ////////////////
        for(unsigned int iblk=0;iblk<nRB+nConst;iblk++){ 
            const unsigned int lencol = aIndDof[iblk+1] - aIndDof[iblk];
            for(unsigned int icrs=ColInd[iblk];icrs<ColInd[iblk+1];icrs++){
                const unsigned int jblk0 = RowPtr[icrs];
                const unsigned int lenrow = aIndDof[jblk0+1] - aIndDof[jblk0];
                ValPtr[icrs] = val_pos;
                val_pos += lencol*lenrow;
            }
        }
        ValPtr[ncrs] = val_pos;
        const unsigned int ntotval = val_pos;
        ValList = new double [ntotval];

        ////////////////
    }
//  mat = new double [ndof*ndof];
    Update = new double [ndof*ndof];
    Residual = new double [ndof*ndof];
}


double Ls::CLinearSystem_RigidBody_CRS::DOT(int iv1,int iv2){

 	double* p_vec1 = 0;
    if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = aTmpVec[iv1];
    else if( iv1 == -1 ) p_vec1 = this->Residual;
    else if( iv1 == -2 ) p_vec1 = this->Update;
    else assert(0);

    double* p_vec2 = 0;
    if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = aTmpVec[iv2];
    else if( iv2 == -1 ) p_vec2 = this->Residual;
    else if( iv2 == -2 ) p_vec2 = this->Update;
    else assert(0);

    double d0 = 0;
    for(unsigned int i=0;i<ndof;i++){ d0 += p_vec2[i]*p_vec1[i]; }
    return d0;
}

bool Ls::CLinearSystem_RigidBody_CRS::COPY(int iv1,int iv2)  // {v2} := {v1}
{
    if( iv1 == iv2 ) return true;

   	double* p_vec1 = 0;
    if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = aTmpVec[iv1];
    else if( iv1 == -1 ) p_vec1 = this->Residual;
    else if( iv1 == -2 ) p_vec1 = this->Update;
    else assert(0);

    double* p_vec2 = 0;
    if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = aTmpVec[iv2];
    else if( iv2 == -1 ) p_vec2 = this->Residual;
    else if( iv2 == -2 ) p_vec2 = this->Update;
    else assert(0);

    for(unsigned int i=0;i<ndof;i++){ p_vec2[i] = p_vec1[i]; }
    return true;
}

bool Ls::CLinearSystem_RigidBody_CRS::SCAL(double d,int iv)  // {v1} := alpha * {v1}
{
    double* p_vec = 0;
	if( iv >= 0 && iv < (int)this->GetTmpVectorArySize() ) p_vec = aTmpVec[iv];
    else if( iv == -1 ) p_vec = this->Residual;
    else if( iv == -2 ) p_vec = this->Update;
    else assert(0);

    for(unsigned int i=0;i<ndof;i++){ p_vec[i] *= d; }
    return true;
}
   
bool Ls::CLinearSystem_RigidBody_CRS::AXPY(double d,int iv1,int iv2) // {v2} := alpha*{v1} + {v2}
{ 
   	double* p_vec1 = 0;
    if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = aTmpVec[iv1];
    else if( iv1 == -1 ) p_vec1 = this->Residual;
    else if( iv1 == -2 ) p_vec1 = this->Update;
    else assert(0);

    double* p_vec2 = 0;
    if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = aTmpVec[iv2];
    else if( iv2 == -1 ) p_vec2 = this->Residual;
    else if( iv2 == -2 ) p_vec2 = this->Update;
    else assert(0);

    for(unsigned int i=0;i<ndof;i++){ p_vec2[i] += d*p_vec1[i]; }
    return true;
}

bool Ls::CLinearSystem_RigidBody_CRS::MATVEC(double a,int iv1,double b,int iv2)  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
{
	double* p_vec1 = 0;
    if( iv1 >= 0 && iv1 < (int)this->GetTmpVectorArySize() ) p_vec1 = aTmpVec[iv1];
    else if( iv1 == -1 ) p_vec1 = this->Residual;
    else if( iv1 == -2 ) p_vec1 = this->Update;
    else assert(0);

    double* p_vec2 = 0;
    if( iv2 >= 0 && iv2 < (int)this->GetTmpVectorArySize() ) p_vec2 = aTmpVec[iv2];
    else if( iv2 == -1 ) p_vec2 = this->Residual;
    else if( iv2 == -2 ) p_vec2 = this->Update;
    else assert(0);

    for(unsigned int i=0;i<ndof;i++){ p_vec2[i] *= b; }
    const unsigned int nblk = nRB + nConst;
    for(unsigned int iblk=0;iblk<nblk;iblk++){
        const unsigned int idofcol = aIndDof[iblk];
        const unsigned int lencol = aIndDof[iblk+1]-aIndDof[iblk];
        for(unsigned int icrs=ColInd[iblk];icrs<ColInd[iblk+1];icrs++){
            const unsigned int jblk0 = RowPtr[icrs];
            const unsigned int idofrow = aIndDof[jblk0];
            const unsigned int lenrow = aIndDof[jblk0+1]-aIndDof[jblk0];
//                std::cout << iblk << " " << jblk0 << "   " << lencol << " " << lenrow << std::endl;
            const unsigned int idofstat = ValPtr[icrs];
            for(unsigned int i=0;i<lencol;i++){
            for(unsigned int j=0;j<lenrow;j++){
                p_vec2[idofcol+i] += a*ValList[idofstat+i*lenrow+j]*p_vec1[idofrow+j];
            }
            }
        }
    }
    return true; 
}

void Ls::CLinearSystem_RigidBody_CRS::Solve(const Ls::CPreconditioner_RigidBody_CRS& prec)
{
    ndof = this->ndof;
    for(unsigned int idof=0;idof<ndof;idof++){ 
        Update[idof] = Residual[idof];
    }
    prec.Solve(Update);
}


bool Ls::CLinearSystem_RigidBody_CRS::UpdateValueOfRigidSystem(
    std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
    double dt, double newmark_gamma, double newmark_beta, 
    bool is_first) const
{
    for(unsigned int irb=0;irb<aRB.size();irb++){
        aRB[irb].UpdateSolution( this->GetUpdate(irb,true), 
            dt, newmark_gamma, newmark_beta, 
            is_first);
    }
    for(unsigned int icst=0;icst<aConst.size();icst++){
        aConst[icst]->UpdateSolution( this->GetUpdate(icst,false), 
            dt, newmark_gamma, newmark_beta );
    }
    return true;
}






////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////






void Ls::CPreconditioner_RigidBody_CRS::SetLinearSystem(const Ls::CLinearSystem_RigidBody_CRS& ls)
{
    this->nRB = ls.nRB;
    this->nConst = ls.nConst;
    this->aIndDof = ls.aIndDof;
    this->ndof = ls.ndof;
    const unsigned int nblk = nRB + nConst;
    unsigned int nelem = 0;   // 総要素数(値のサイズを知る)
    std::vector<unsigned int> aPtn;
    {
        aPtn.resize(nblk*nblk,0);
        // パターンの初期化
        for(unsigned int iblk=0;iblk<nblk;iblk++){
            for(unsigned int icrs=ls.ColInd[iblk];icrs<ls.ColInd[iblk+1];icrs++){
                const unsigned int jblk0 = ls.RowPtr[icrs];
                aPtn[iblk*nblk+jblk0] = 1;
                nelem++;
            }
        }
        for(unsigned int iblk=0;iblk<nblk;iblk++){
            for(unsigned int jblk=0;jblk<nblk;jblk++){
                if( aPtn[iblk*nblk+jblk] != 0 ) continue;
                int kmax = (iblk<jblk) ? iblk : jblk;
                for(int kblk=0;kblk<kmax;kblk++){
                    const unsigned int levik = aPtn[iblk*nblk+kblk];
                    const unsigned int levkj = aPtn[kblk*nblk+jblk];
                    if( levik == 0 || levkj == 0 ) continue;
                    const unsigned int levij = ( levik < levkj ) ? levik+1 : levkj+1;
                    aPtn[iblk*nblk+jblk] = levij;
                    nelem++;
                }
            }
        }
    }
    RowPtr = new unsigned int [nelem];
    DiaInd = new unsigned int [nblk];
    ColInd = new unsigned int [nblk+1];
    ColInd[0] = 0;
    unsigned int icnt = 0;
    for(unsigned int iblk=0;iblk<nblk;iblk++){
        assert( aPtn[iblk*nblk+iblk] != 0 );
        bool iflg = false;
        for(unsigned int jblk=0;jblk<nblk;jblk++){
            if( aPtn[iblk*nblk+jblk] == 0 ) continue;
            RowPtr[icnt] = jblk;
            if( iblk == jblk ){
                DiaInd[iblk] = icnt;
                iflg = true;
            }
            icnt++;
        }
        assert( iflg );
        ColInd[iblk+1] = icnt;
    }
    ncrs = ColInd[nblk];
    ValPtr = new unsigned int [ncrs+1];
    unsigned int ntotdof = 0;
    for(unsigned int iblk=0;iblk<nRB+nConst;iblk++){ 
       const unsigned int lencol = aIndDof[iblk+1] - aIndDof[iblk];
        for(unsigned int icrs=ColInd[iblk];icrs<ColInd[iblk+1];icrs++){
            const unsigned int jblk0 = RowPtr[icrs];
            const unsigned int lenrow = aIndDof[jblk0+1] - aIndDof[jblk0];
            ValPtr[icrs] = ntotdof;
            ntotdof += lencol*lenrow;
        }
    }
    ValPtr[ncrs] = ntotdof;
    ValList = new double [ntotdof];
}

void Ls::CPreconditioner_RigidBody_CRS::SetValue(const CLinearSystem_RigidBody_CRS& ls)
{        
    const unsigned int nblk = nRB + nConst;
    const unsigned int ntotdof = ValPtr[ncrs];
    for(unsigned int i=0;i<ntotdof;i++){ ValList[i] = 0.0; }
    for(unsigned int iblk=0;iblk<nblk;iblk++){
        for(unsigned int icrs_in=ls.ColInd[iblk];icrs_in<ls.ColInd[iblk+1];icrs_in++){
            const unsigned int jblk0 = ls.RowPtr[icrs_in];
            unsigned icrs = ColInd[iblk];
            for(;icrs<ColInd[iblk+1];icrs++){
                if( RowPtr[icrs] == jblk0 ) break;
            }
            assert( icrs != ColInd[iblk+1] );
            double* val    = &ValList[       ValPtr[icrs   ] ];
            double* val_in = &ls.ValList[ ls.ValPtr[icrs_in] ];
            const unsigned int lenblk_i = aIndDof[iblk +1] - aIndDof[iblk ];
            const unsigned int lenblk_j = aIndDof[jblk0+1] - aIndDof[jblk0];
            for(unsigned int i=0;i<lenblk_i*lenblk_j;i++){ val[i] = val_in[i]; }
        }
    }
    ////////////////////////////////
    int* row2crs_f = new int [nblk];
	double pTmpBlk[6*6];
	for(unsigned int iblk=0;iblk<nblk;iblk++){
        const unsigned int lenblk_i = aIndDof[iblk+1] - aIndDof[iblk];
		for(unsigned int ijcrs=ColInd[iblk];ijcrs<ColInd[iblk+1];ijcrs++){
			assert( ijcrs<ncrs );
			const unsigned int jblk0 = RowPtr[ijcrs];
			assert( jblk0<nblk );
			row2crs_f[jblk0] = ijcrs;
		}
		for(unsigned int ikcrs=ColInd[iblk];ikcrs<DiaInd[iblk];ikcrs++){
			const unsigned int kblk = RowPtr[ikcrs];
            const unsigned int lenblk_k = aIndDof[kblk+1] - aIndDof[kblk];
			assert( kblk<nblk );
			const double* pVal_ik = &ValList[ ValPtr[ikcrs] ];
			for(unsigned int kjcrs=DiaInd[kblk]+1;kjcrs<ColInd[kblk+1];kjcrs++){
				const unsigned int jblk0 = RowPtr[kjcrs];
                const unsigned int lenblk_j = aIndDof[jblk0+1] - aIndDof[jblk0];
				assert( jblk0<nblk );
				const int ijcrs0 = row2crs_f[jblk0];
				if( ijcrs0 == -1 ) continue;
				double* pVal_kj = &ValList[ ValPtr[kjcrs ] ];
				double* pVal_ij = &ValList[ ValPtr[ijcrs0] ];
				assert( pVal_kj != 0 );
				assert( pVal_ij != 0 );
				for(unsigned int i=0;i<lenblk_i;i++){
				for(unsigned int j=0;j<lenblk_j;j++){
				for(unsigned int k=0;k<lenblk_k;k++){
					pVal_ij[i*lenblk_j+j] -= pVal_ik[i*lenblk_k+k]*pVal_kj[k*lenblk_j+j];
				}
				}
				}
//                std::cout << iblk << " " << kblk << " " << jblk0 << std::endl;
			}
		}
		{
			double* pVal_ii = &ValList[ ValPtr[ DiaInd[iblk] ] ];
//            for(unsigned int i=0;i<lenblk_i;i++){
//            for(unsigned int j=0;j<lenblk_i;j++){
//                std::cout << i << " " << j << " " << pVal_ii[i*lenblk_i+j] << std::endl;
//            }
//            }
			int info = 0;
			CalcInvMat(pVal_ii,lenblk_i,info);
			if( info==1 ){
				std::cout << "frac false" << iblk << std::endl;
				// TODO:?????????????????
//				getchar();
			}
		}
		// ????????????????????[U] = [1/D][U]
		for(unsigned int ijcrs=DiaInd[iblk]+1;ijcrs<ColInd[iblk+1];ijcrs++){
			assert( ijcrs<ncrs );
			double* pVal_ij = &ValList[ ValPtr[ijcrs] ];
			const double* pVal_ii = &ValList[ ValPtr[ DiaInd[iblk] ] ];
            const unsigned int jblk0 = RowPtr[ijcrs];
            const unsigned int lenblk_j = aIndDof[jblk0+1] - aIndDof[jblk0];
			for(unsigned int i=0;i<lenblk_i*lenblk_j;i++){ pTmpBlk[i] = pVal_ij[i]; }
			for(unsigned int i=0;i<lenblk_i;i++){
			for(unsigned int j=0;j<lenblk_j;j++){
				double dtmp0 = 0.0;
				for(unsigned int k=0;k<lenblk_i;k++){
					dtmp0 += pVal_ii[i*lenblk_i+k]*pTmpBlk[k*lenblk_j+j];
				}
				pVal_ij[i*lenblk_j+j] = dtmp0;
			}
			}
		}
		for(unsigned int ijcrs=ColInd[iblk];ijcrs<ColInd[iblk+1];ijcrs++){
			assert( ijcrs<ncrs );
			const unsigned int jblk0 = RowPtr[ijcrs];
			assert( jblk0<nblk );
			row2crs_f[jblk0] = -1;
		}
	}	// end iblk
	delete[] row2crs_f;
}

bool Ls::CPreconditioner_RigidBody_CRS::Solve( double* vec ) const
{
    const unsigned int nblk = nRB + nConst;
    double pTmpVec[6*6];
	for(unsigned int iblk=0;iblk<nblk;iblk++){
        const unsigned int ilenblk = aIndDof[iblk+1] - aIndDof[iblk];
        const unsigned int itotdof = aIndDof[iblk];
		for(unsigned int idof=0;idof<ilenblk;idof++){
			pTmpVec[idof] = vec[ itotdof+idof ];
		}
		for(unsigned int ijcrs=ColInd[iblk];ijcrs<DiaInd[iblk];ijcrs++){
			assert( ijcrs<ncrs );
			const unsigned int jblk0 = RowPtr[ijcrs];
            const unsigned int jlenblk = aIndDof[jblk0+1] - aIndDof[jblk0];
            const unsigned int jtotdof = aIndDof[jblk0];
			assert( jblk0<iblk );
			const double* pVal_ij = &ValList[ ValPtr[ijcrs] ];
			for(unsigned int idof=0;idof<ilenblk;idof++){
				for(unsigned int jdof=0;jdof<jlenblk;jdof++){
					pTmpVec[idof] -= pVal_ij[idof*jlenblk+jdof]*vec[ jtotdof+jdof ];
				}
			}
		}
		const double* pVal_ii = &ValList[ ValPtr[ DiaInd[iblk] ] ];
		for(unsigned int idof=0;idof<ilenblk;idof++){
			double dtmp1 = 0.0;
			for(unsigned int jdof=0;jdof<ilenblk;jdof++){
				dtmp1 += pVal_ii[idof*ilenblk+jdof]*pTmpVec[jdof];
			}
			vec[itotdof+idof] = dtmp1;
		}
	}
	for(int iblk=nblk-1;iblk>=0;iblk--){
        const unsigned int ilenblk = aIndDof[iblk+1] - aIndDof[iblk];
        const unsigned int itotdof = aIndDof[iblk];
		assert( (unsigned int)iblk < nblk );
		for(unsigned int idof=0;idof<ilenblk;idof++){
			pTmpVec[idof] = vec[itotdof+idof];
		}
		for(unsigned int ijcrs=DiaInd[iblk]+1;ijcrs<ColInd[iblk+1];ijcrs++){
			assert( ijcrs<ncrs );
			const unsigned int jblk0 = RowPtr[ijcrs];
            const unsigned int jlenblk = aIndDof[jblk0+1] - aIndDof[jblk0];
            const unsigned int jtotdof = aIndDof[jblk0];
			assert( jblk0>(unsigned int)iblk && jblk0<nblk );
			const double* pVal_ij = &ValList[ ValPtr[ijcrs] ];
			for(unsigned int idof=0;idof<ilenblk;idof++){
				for(unsigned int jdof=0;jdof<jlenblk;jdof++){
					pTmpVec[idof] -= pVal_ij[idof*jlenblk+jdof]*vec[jtotdof+jdof];
				}
			}
		}
		for(unsigned int idof=0;idof<ilenblk;idof++){
			vec[itotdof+idof] = pTmpVec[idof];
		}
	}
	return true;
}
*/

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void Ls::CLinearSystem_RigidBody_CRS2::SetRigidSystem(
    const std::vector<Rigid::CRigidBody3D>& aRB,
    const std::vector<Rigid::CConstraint*>& aConst)
{
  this->Clear();
  nRB = aRB.size();
  nConst = aConst.size();
  const unsigned int nblk = nRB+nConst;
  {
    m_aBlkSize.resize(nblk);
    for(unsigned int irb=0;irb<nRB;irb++){
      m_aBlkSize[irb] = aRB[irb].GetDOF();
    }
    for(unsigned int icst=0;icst<nConst;icst++){
      m_aBlkSize[nRB+icst] = aConst[icst]->GetDOF();
    }
    m_mat.Initialize(nblk,m_aBlkSize, nblk,m_aBlkSize);
    m_residual.Initialize(nblk,m_aBlkSize);
    m_update.Initialize(nblk,m_aBlkSize);
  }
  {
    Com::CIndexedArray crs;
    crs.InitializeSize(nblk);
    crs.index.resize(nblk+1);
    crs.index[0] = 0;
    for(unsigned int i=0;i<nRB;i++){ crs.index[i+1]=0; }
    for(unsigned int i=nRB;i<nRB+nConst;i++){ crs.index[i+1]=0; }
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        crs.index[irb0+1] += 1;
      }
      crs.index[icst+nRB+1] += aIndRB.size();
    }
    for(unsigned int i=0;i<nblk;i++){ 
      crs.index[i+1] = crs.index[i+1] + crs.index[i];
    }
    const unsigned int ncrs = crs.index[nblk];
    crs.array.resize(ncrs);
    for(unsigned int icst=0;icst<nConst;icst++){
      const std::vector<unsigned int>& aIndRB = aConst[icst]->GetAry_IndexRB();
      for(unsigned int i=0;i<aIndRB.size();i++){
        const unsigned int irb0 = aIndRB[i];
        const unsigned int icrs0 = crs.index[icst+nRB];
        crs.array[icrs0] = irb0;
        crs.index[icst+nRB] += 1;
        const unsigned int icrs1 = crs.index[irb0];
        crs.array[icrs1] = icst+nRB;
        crs.index[irb0] += 1;
      }
    }
    for(int i=nblk;i>0;i--){ 
      crs.index[i] = crs.index[i-1];
    }
    crs.index[0] = 0;
    crs.Sort();
    assert( crs.index[nRB+nConst] == ncrs );
    m_mat.AddPattern(crs);      
  }  
}


double Ls::CLinearSystem_RigidBody_CRS2::DOT(int iv1,int iv2)
{
    MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
    MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
    return vec1*vec2;
}

bool Ls::CLinearSystem_RigidBody_CRS2::COPY(int iv1,int iv2)  // {v2} := {v1}
{
    if( iv1 == iv2 ) return true;

    MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
    MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
    vec2 = vec1;
    return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::SCAL(double d,int iv)  // {v1} := alpha * {v1}
{
    MatVec::CVector_Blk& vec = this->GetVector(iv);
    vec *= d;
    return true;
}
   
bool Ls::CLinearSystem_RigidBody_CRS2::AXPY(double d,int iv1,int iv2) // {v2} := alpha*{v1} + {v2}
{ 
    MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
    MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
    vec2.AXPY(d,vec1);
    return true;
}

bool Ls::CLinearSystem_RigidBody_CRS2::MATVEC(double a,int iv1,double b,int iv2)  //  {v2} := alpha*[MATRIX]*{v1} + beta*{v2}
{
    MatVec::CVector_Blk& vec1 = this->GetVector(iv1);
    MatVec::CVector_Blk& vec2 = this->GetVector(iv2);
    m_mat.MatVec(a,vec1,b,vec2);
    return true; 
}



bool Ls::CLinearSystem_RigidBody_CRS2::UpdateValueOfRigidSystem(
    std::vector<Rigid::CRigidBody3D>& aRB, std::vector<Rigid::CConstraint*>& aConst, 
    double dt, double newmark_gamma, double newmark_beta, 
    bool is_first) const
{
    double tmp[6];
    for(unsigned int irb=0;irb<aRB.size();irb++){
        const unsigned int nlen = m_update.Len(irb);
        assert( nlen <= 6 );
        for(unsigned int ilen=0;ilen<nlen;ilen++){
            tmp[ilen] = m_update.GetValue(irb,ilen);
        }
        aRB[irb].UpdateSolution( tmp, 
            dt, newmark_gamma, newmark_beta, 
            is_first);
    }
    for(unsigned int icst=0;icst<aConst.size();icst++){
        const unsigned int nlen = m_update.Len(icst+nRB);
        assert( nlen <= 6 );
        for(unsigned int ilen=0;ilen<nlen;ilen++){
            tmp[ilen] = m_update.GetValue(icst+nRB,ilen);
        }
        aConst[icst]->UpdateSolution( tmp, 
            dt, newmark_gamma, newmark_beta );
    }
    return true;
}
