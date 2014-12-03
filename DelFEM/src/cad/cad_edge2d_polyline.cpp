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

#define for if(0);else for

#include <math.h>

#include "delfem/cad/cad_edge2d_polyline.h"

// ÉuÉçÉbÉN(ÇRÅ~ÇR)ÇRèdëŒäpçsóÒÇâÇ≠
class CTriDiaMat3
{
public:
	// initialize with block size n
	CTriDiaMat3(unsigned int n){
		this->n = n;
		v = new double [(n*3-2)*9];
	}
	~CTriDiaMat3(){ delete[] v; }
	// clear value
	void Clear(){ for(unsigned int i=0;i<(n*3-2)*9;i++){ v[i] = 0; } }
	// marge element stiffness matrix to the position (idiv,idiv+1)
	void Marge(unsigned int idiv, double eM[][2][3][3]){
		for(unsigned int i=0;i<36;i++){ v[idiv*27+i] += (&eM[0][0][0][0])[i]; }
	}
	// define fixed boudnary condition
	void FixBC(unsigned int ino, unsigned int idof){
		assert( idof < 3 && ino < n );
		if( ino != 0 ){
			double* pvu = v+ino*27-18;
			double* pvl = v+ino*27-9;
			pvu[0*3+idof] = 0;	pvu[1*3+idof] = 0;	pvu[2*3+idof] = 0;
			pvl[idof*3+0] = 0;	pvl[idof*3+1] = 0;	pvl[idof*3+2] = 0;
		}
		if( ino != n-1 ){
			double* pvu = v+ino*27+18;
			double* pvl = v+ino*27+9;
			pvu[0*3+idof] = 0;	pvu[1*3+idof] = 0;	pvu[2*3+idof] = 0;
			pvl[idof*3+0] = 0;	pvl[idof*3+1] = 0;	pvl[idof*3+2] = 0;
		}
		double* pvc = v+ino*27;
		pvc[0*3+idof] = 0;	pvc[1*3+idof] = 0;	pvc[2*3+idof] = 0;
		pvc[idof*3+0] = 0;	pvc[idof*3+1] = 0;	pvc[idof*3+2] = 0;
		pvc[idof*3+idof] = 1;
	}
	// execute ILU factorization
	void ILU_Frac()
	{
		double tmpBlk[9];
		for(unsigned int iblk=0;iblk<n;iblk++){
			if( iblk != 0 ){
				const double* pVal_ik = v+27*iblk-9;
				const double* pVal_kj = v+27*iblk-18;
				double* pVal_ij = v+27*iblk;
				for(unsigned int i=0;i<3;i++){
					pVal_ij[i*3+0] -= pVal_ik[i*3+0]*pVal_kj[0] + pVal_ik[i*3+1]*pVal_kj[3] + pVal_ik[i*3+2]*pVal_kj[6];
					pVal_ij[i*3+1] -= pVal_ik[i*3+0]*pVal_kj[1] + pVal_ik[i*3+1]*pVal_kj[4] + pVal_ik[i*3+2]*pVal_kj[7];
					pVal_ij[i*3+2] -= pVal_ik[i*3+0]*pVal_kj[2] + pVal_ik[i*3+1]*pVal_kj[5] + pVal_ik[i*3+2]*pVal_kj[8];
				}
			}
			{   // calc inverse of diagonal
				double* pVal_ii = v+27*iblk;
        CalcInvMat3(pVal_ii,tmpBlk);
			}
			// ëŒäpÇÃãtêîÇåvéZÇµÇƒè„éOäpçsóÒÇ…ä|ÇØÇÈÅB[U] = [1/D][U]
			if( iblk !=  n-1 ){
				double* pVal_ij = v+27*iblk+9;
				const double* pVal_ii = v+27*iblk;
				for(unsigned int i=0;i<9;i++){ tmpBlk[i] = pVal_ij[i]; }
				for(unsigned int i=0;i<3;i++){
					pVal_ij[i*3+0] = pVal_ii[i*3+0]*tmpBlk[0] + pVal_ii[i*3+1]*tmpBlk[3] + pVal_ii[i*3+2]*tmpBlk[6];
					pVal_ij[i*3+1] = pVal_ii[i*3+0]*tmpBlk[1] + pVal_ii[i*3+1]*tmpBlk[4] + pVal_ii[i*3+2]*tmpBlk[7];
					pVal_ij[i*3+2] = pVal_ii[i*3+0]*tmpBlk[2] + pVal_ii[i*3+1]*tmpBlk[5] + pVal_ii[i*3+2]*tmpBlk[8];
				}
			}
		}	// end iblk
	}
	// solve matrix
	void Solve(double* res){
		double pTmpVec[3];
		for(unsigned int iblk=0;iblk<n;iblk++){
			pTmpVec[0] = res[iblk*3+0];
			pTmpVec[1] = res[iblk*3+1];
			pTmpVec[2] = res[iblk*3+2];
			if( iblk != 0 ){
				const double* pVal_ij = v+iblk*27-9;
				const double valj0 = res[(iblk-1)*3+0];
				const double valj1 = res[(iblk-1)*3+1];
				const double valj2 = res[(iblk-1)*3+2];
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
				pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
				pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
			}
			const double* pVal_ii = v+27*iblk;
			res[iblk*3+0] = pVal_ii[0]*pTmpVec[0]+pVal_ii[1]*pTmpVec[1]+pVal_ii[2]*pTmpVec[2];
			res[iblk*3+1] = pVal_ii[3]*pTmpVec[0]+pVal_ii[4]*pTmpVec[1]+pVal_ii[5]*pTmpVec[2];
			res[iblk*3+2] = pVal_ii[6]*pTmpVec[0]+pVal_ii[7]*pTmpVec[1]+pVal_ii[8]*pTmpVec[2];
		}
		for(int iblk=n-1;iblk>=0;iblk--){
			pTmpVec[0] = res[iblk*3+0];
			pTmpVec[1] = res[iblk*3+1];
			pTmpVec[2] = res[iblk*3+2];
      if( iblk != (int)n-1 ){
				const double* pVal_ij = v+27*iblk+9;
				const double valj0 = res[(iblk+1)*3+0];
				const double valj1 = res[(iblk+1)*3+1];
				const double valj2 = res[(iblk+1)*3+2];
				pTmpVec[0] -= pVal_ij[0]*valj0+pVal_ij[1]*valj1+pVal_ij[2]*valj2;
				pTmpVec[1] -= pVal_ij[3]*valj0+pVal_ij[4]*valj1+pVal_ij[5]*valj2;
				pTmpVec[2] -= pVal_ij[6]*valj0+pVal_ij[7]*valj1+pVal_ij[8]*valj2;
			}
			res[iblk*3+0] = pTmpVec[0];
			res[iblk*3+1] = pTmpVec[1];
			res[iblk*3+2] = pTmpVec[2];
		}
	}
private:
	// ÇRÅ~ÇRÇÃãtçsóÒÇãÅÇﬂÇÈ
	static inline void CalcInvMat3(double a[], double t[] ){
		const double det = a[0]*a[4]*a[8] + a[3]*a[7]*a[2] + a[6]*a[1]*a[5]
						 - a[0]*a[7]*a[5] - a[6]*a[4]*a[2] - a[3]*a[1]*a[8];
		const double inv_det = 1.0/det;
		for(int i=0;i<9;i++){ t[i] = a[i]; }
		a[0] = inv_det*(t[4]*t[8]-t[5]*t[7]);
		a[1] = inv_det*(t[2]*t[7]-t[1]*t[8]);
		a[2] = inv_det*(t[1]*t[5]-t[2]*t[4]);
		a[3] = inv_det*(t[5]*t[6]-t[3]*t[8]);
		a[4] = inv_det*(t[0]*t[8]-t[2]*t[6]);
		a[5] = inv_det*(t[2]*t[3]-t[0]*t[5]);
		a[6] = inv_det*(t[3]*t[7]-t[4]*t[6]);
		a[7] = inv_det*(t[1]*t[6]-t[0]*t[7]);
		a[8] = inv_det*(t[0]*t[4]-t[1]*t[3]);
	}
private:
	unsigned int n;
	double* v;
};

CCadEdge2DPolyline::CCadEdge2DPolyline()
{
  this->EI = 30.0;
  this->ARho = 1;
  this->AE = 100000.0;  
  ////
	nno = 0;
  ////
	ut = 0;
  ini_x = 0;
	bc_flag = 0;
	dut = 0;
  Res = 0;
  m_mat = 0;
}

void CCadEdge2DPolyline::SetCadEdge(const Cad::CCadObj2D& cad_2d, const unsigned int id_e,
                                    const Com::CVector2D& pick_pos)
{
  this->pick_pos = pick_pos;
  this->m_IdECad = id_e;
  this->ClearMemory();
  std::vector<double> aRelCoPoly;
  cad_2d.GetCurve_Polyline(id_e,aRelCoPoly);
  std::vector<double> aXYs;
  const unsigned int ndiv = aRelCoPoly.size()/2+1;
  nno = ndiv+1;
  aXYs.reserve(nno*2);
  unsigned int id_vs = cad_2d.GetIdVertex_Edge(id_e,true );
  unsigned int id_ve = cad_2d.GetIdVertex_Edge(id_e,false);
  const Com::CVector2D& vs = cad_2d.GetVertexCoord(id_vs);
  const Com::CVector2D& ve = cad_2d.GetVertexCoord(id_ve);
  const Com::CVector2D hse = ve-vs;
  const Com::CVector2D vse = Com::CVector2D(-hse.y,hse.x);  
  aXYs.push_back(vs.x);
  aXYs.push_back(vs.y);    
  for(unsigned int i=0;i<ndiv-1;i++){
    Com::CVector2D v = vs+aRelCoPoly[i*2+0]*hse+aRelCoPoly[i*2+1]*vse;
    aXYs.push_back(v.x);
    aXYs.push_back(v.y);    
  }
  aXYs.push_back(ve.x);
  aXYs.push_back(ve.y);      
  //////////////////////////////////////
  this->m_mat = new CTriDiaMat3(nno);
  ini_x = new double [nno*2];
  ut = new double [nno*3];
  Res = new double [nno*3];  
  dut = new double [nno*3];
  bc_flag = new int [nno*3];  
  ////
  for(unsigned int i=0;i<nno*2;i++){ ini_x[i] = aXYs[i]; }  
  for(unsigned int i=0;i<nno*3;i++){ ut[i] = 0; }    
  for(unsigned int i=0;i<nno*3;i++){ bc_flag[i] = 0; }  
  //////////////////////////////////////
  this->SetFixedBoundaryFlag(0,0);
  this->SetFixedBoundaryFlag(0,1);
  this->SetFixedBoundaryFlag(nno-1,0);
  this->SetFixedBoundaryFlag(nno-1,1);    
  double alpha, ndist, norm_x, norm_y;
  ProjectPoint(pick_pos.x, pick_pos.y, idiv_picked,
               alpha, ndist, norm_x, norm_y);
  if( idiv_picked == -1 ) return;
  this->SetFixedBoundaryFlag(idiv_picked,0);
  this->SetFixedBoundaryFlag(idiv_picked,1); 
  this->SetFixedBoundaryFlag(idiv_picked+1,0);
  this->SetFixedBoundaryFlag(idiv_picked+1,1);     
}

void CCadEdge2DPolyline::Drag(Cad::CCadObj2D& cad_2d, const Com::CVector2D& dist_pos)
{
  const Com::CVector2D del = dist_pos-pick_pos;
  this->SetDisp(idiv_picked,   0, del.x);
  this->SetDisp(idiv_picked,   1, del.y);
  this->SetDisp(idiv_picked+1, 0, del.x);
  this->SetDisp(idiv_picked+1, 1, del.y);  
  
  this->SolveLinearStatic();
  std::vector<Com::CVector2D> aXYs;
  for(unsigned int i=1;i<nno-1;i++){
    double x,y,t;
    GetValueNode(i, x, y, t);
    aXYs.push_back( Com::CVector2D(x,y) );
  }
  cad_2d.SetCurve_Polyline(m_IdECad,aXYs);
}


void CCadEdge2DPolyline::ClearMemory()
{
  if( ini_x   != 0 ){ delete[] ini_x;   ini_x   = 0; }
  if( ut      != 0 ){ delete[] ut;      ut      = 0; }
  if( Res     != 0 ){ delete[] Res;     Res     = 0; }
  if( ut      != 0 ){ delete[] dut;     dut     = 0; }
  if( bc_flag != 0 ){ delete[] bc_flag; bc_flag = 0; }
  if( m_mat   != 0 ){ delete m_mat;     m_mat   = 0; }    
}


CCadEdge2DPolyline::~CCadEdge2DPolyline(){
  this->ClearMemory();
}

void CCadEdge2DPolyline::ProjectPoint(double x_in, double y_in, int& idiv_min,
    double& alpha, double& ndist, double& norm_x, double& norm_y)
{
  double x0[2] = { x_in, y_in };
  double dist = sqrt( (ini_x[0]+ut[0]-x0[0])*(ini_x[0]+ut[0]-x0[0])
                     + (ini_x[1]+ut[1]-x0[1])*(ini_x[1]+ut[1]-x0[1]) );
  idiv_min = -1;
  const unsigned int ndiv = nno-1;
  for(unsigned int idiv=0;idiv<ndiv;idiv++){
    double x1[2] = { ini_x[idiv*2+0]+ut[idiv*3+0], ini_x[idiv*2+1]+ut[idiv*3+1] };
    double x2[2] = { ini_x[idiv*2+2]+ut[idiv*3+3], ini_x[idiv*2+3]+ut[idiv*3+4] };    
    double t = FindNearestPointParam_Line_Point(x0,x1,x2);
    if( t < -0.001 || t > 1.001 ){ continue; }
    double x3[2] = { x1[0]*(1-t)+x2[0]*t, x1[1]*(1-t)+x2[1]*t };
    double d[2] = { x0[0]-x3[0], x0[1]-x3[1] };
    double elen = sqrt( (x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]) );
    double n[2] = { (x2[1]-x1[1])/elen, (x1[0]-x2[0])/elen };
    double dlen = sqrt( d[0]*d[0] + d[1]*d[1] );
    if( dlen < dist ){
      idiv_min = idiv;
      dist = dlen;
      alpha = t;
      ndist = (n[0]*d[0] + n[1]*d[1])/elen;
      norm_x = n[0];
      norm_y = n[1];
    }      
  }
  return;  
}

// îÒê¸å`ó¿óvëfÇÃóvëfçÑê´çsóÒ
void GetEMat_Linear(double EI, double AE,
                    const double x1[2], const double u1[2], double t1,
                    const double x2[2], const double u2[2], double t2,
                    double eQ[6], double eK[][6])
{
    const double eLen = sqrt( (x2[0]-x1[0])*(x2[0]-x1[0]) + (x2[1]-x1[1])*(x2[1]-x1[1]) );
    const double inv_eLen = 1.0/eLen;

    double eC[6][6];
    {
        const double tmp1 = EI/( eLen*eLen*eLen );
        const double tmp2 = AE/eLen;
        eC[0][0]= tmp2; eC[0][1]= 0;           eC[0][2]= 0;                eC[0][3]=-tmp2; eC[0][4]= 0;           eC[0][5]=0;
        eC[1][0]= 0;    eC[1][1]= tmp1*12;	   eC[1][2]= tmp1*eLen*6;	   eC[1][3]=0;     eC[1][4]=-tmp1*12;	  eC[1][5]= tmp1*eLen*6;
        eC[2][0]= 0;    eC[2][1]= tmp1*eLen*6; eC[2][2]= tmp1*eLen*eLen*4; eC[2][3]=0;     eC[2][4]=-tmp1*eLen*6; eC[2][5]= tmp1*eLen*eLen*2;
        eC[3][0]=-tmp2; eC[3][1]= 0;           eC[3][2]= 0;                eC[3][3]= tmp2; eC[3][4]= 0;           eC[3][5]=0;
        eC[4][0]= 0;    eC[4][1]=-tmp1*12;	   eC[4][2]=-tmp1*eLen*6;	   eC[4][3]=0;     eC[4][4]= tmp1*12;	  eC[4][5]=-tmp1*eLen*6;
        eC[5][0]= 0;    eC[5][1]= tmp1*eLen*6; eC[5][2]= tmp1*eLen*eLen*2; eC[5][3]=0;     eC[5][4]=-tmp1*eLen*6; eC[5][5]= tmp1*eLen*eLen*4;
    }
    const double T[2] = { (x2[0]-x1[0])*inv_eLen, (x2[1]-x1[1])*inv_eLen };
    double eR[6][6];
    {
        eR[0][0]= T[0]; eR[0][1]=-T[1]; eR[0][2]=0; eR[0][3]=    0; eR[0][4]=    0; eR[0][5]= 0;
        eR[1][0]= T[1]; eR[1][1]= T[0]; eR[1][2]=0; eR[1][3]=    0; eR[1][4]=    0; eR[1][5]= 0;
        eR[2][0]=    0; eR[2][1]=    0; eR[2][2]=1; eR[2][3]=    0; eR[2][4]=    0; eR[2][5]= 0;
        eR[3][0]=    0; eR[3][1]=    0; eR[3][2]=0; eR[3][3]= T[0]; eR[3][4]=-T[1]; eR[3][5]= 0;
        eR[4][0]=    0; eR[4][1]=    0; eR[4][2]=0; eR[4][3]= T[1]; eR[4][4]= T[0]; eR[4][5]= 0;
        eR[5][0]=    0; eR[5][1]=    0; eR[5][2]=0; eR[5][3]=    0; eR[5][4]=    0; eR[5][5]= 1;
    }
    for(unsigned int i=0;i<6;i++){
    for(unsigned int j=0;j<6;j++){
        eK[i][j] = 0;
        for(unsigned int k=0;k<6;k++){
        for(unsigned int l=0;l<6;l++){
            eK[i][j] += eR[i][k]*eC[k][l]*eR[j][l];
        }
        }
    }
    }
    for(unsigned int i=0;i<6;i++){
        eQ[i] = eK[i][0]*u1[0]+eK[i][1]*u1[1]+eK[i][2]*t1
               +eK[i][3]*u2[0]+eK[i][4]*u2[1]+eK[i][5]*t2;
    }
}

// îÒê¸å`ó¿óvëfÇÃå„ëﬁÉIÉCÉâÅ[éûä‘êœï™Ç…ÇÊÇÈóvëfåWêîçsóÒ
void GetCoeffMat_Linear(double EI, double AE,
                        double ARho, const double g[2],
                        const double x1[2], const double u1[2], double t1,
                        const double x2[2], const double u2[2], double t2,
                        double eRes[2][3], double eC[2][2][3][3])
{
    for(unsigned int i=0;i<6;i++){ (&eRes[0][0])[i] = 0; }
    for(unsigned int i=0;i<36;i++){ (&eC[0][0][0][0])[i] = 0; }
    ////////////////
    const double eLen = sqrt( (x2[0]-x1[0])*(x2[0]-x1[0]) + (x2[1]-x1[1])*(x2[1]-x1[1]) );
    const double erho = ARho*eLen*0.5;
    double dtmp_x = g[0]*erho;
    double dtmp_y = g[1]*erho;
    {	// ëÃêœóÕÇ∆äµê´óÕÇÃí«â¡
        eRes[0][0] += dtmp_x;
        eRes[0][1] += dtmp_y;
        eRes[1][0] += dtmp_x;
        eRes[1][1] += dtmp_y;
    }
    double eQ[6], eK[6][6];
    GetEMat_Linear(EI, AE, x1,u1,t1,  x2,u2,t2,  eQ,eK);
    // óvëfçÑê´çsóÒÇí«â¡
    for(unsigned int i=0;i<3;i++){
    for(unsigned int j=0;j<3;j++){
        eC[0][0][i][j] += eK[0+i][0+j];
        eC[0][1][i][j] += eK[0+i][3+j];
        eC[1][0][i][j] += eK[3+i][0+j];
        eC[1][1][i][j] += eK[3+i][3+j];
    }
    }
    // ì‡óÕÇécç∑Ç…í«â¡
    for(unsigned int i=0;i<3;i++){
        eRes[0][i] -= eQ[0+i];
        eRes[1][i] -= eQ[3+i];
    }
}


void CCadEdge2DPolyline::SolveLinearStatic()
{
  for(unsigned int i=0;i<nno*3;i++){ Res[i] = 0; }
  m_mat->Clear();
  ////////////////
  const unsigned int ndiv=nno-1;
  for(unsigned int idiv=0;idiv<ndiv;idiv++){
    const double* x1 = ini_x + idiv*2;
    const double* x2 = ini_x + idiv*2+2;
    const double* u1 = ut + idiv*3;
    const double* u2 = ut + idiv*3+3;
    const double t1 = ut[idiv*3+2];
    const double t2 = ut[idiv*3+5];
    double eC[2][2][3][3], eRes[2][3];
    const double g[2] = {0,0};
    GetCoeffMat_Linear(EI, AE, ARho, g,
                       x1, u1, t1,
                       x2, u2, t2,
                       eRes, eC);
    m_mat->Marge(idiv,eC);
    for(unsigned int i=0;i<3;i++){
      Res[idiv*3  +i] += eRes[0][i];
      Res[idiv*3+3+i] += eRes[1][i];
    }
  }
  // ã´äEèåèèàóù
  for(unsigned int ino=0;ino<nno;ino++){
    for(unsigned int idof=0;idof<3;idof++){
      if( bc_flag[ino*3+idof] == 0 ) continue;
      m_mat->FixBC(ino,idof);
      Res[ino*3+idof] = 0;
    }
  }
  {	// écç∑ÇÃÉmÉãÉÄÇãÅÇﬂÇÈ
    double nres = 0;
    for(unsigned int i=0;i<nno*3;i++){ nres += Res[i]*Res[i]; }
    //		std::cout << "itr : " << itr << "    Norm res : " << nres << std::endl;
  }
  m_mat->ILU_Frac();	// ILUï™âé¿çs
  {	// çsóÒÇâÇ¢Çƒïœà ÇÃëùï™ÇãÅÇﬂÇÈ
    for(unsigned int i=0;i<nno*3;i++){ dut[i] = Res[i]; }
    m_mat->Solve(dut);
  }
  for(unsigned int i=0;i<nno*3;i++){ ut[i] += dut[i]; }	// ïœà ÇçXêV
}

