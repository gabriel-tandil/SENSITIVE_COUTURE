
#if !defined(EIGEN_LANCZOS_H)
#define EIGEN_LANCZOS_H

#include <vector>

namespace LsSol{
class CLinearSystem;
class CPreconditioner;

//bool ArnoldiQR(Fem::Ls::CLinearSystem_Eigen& ls );
//bool EigenValue_Lanczos( unsigned int nlambda, std::vector<double>& aLambda, unsigned int num_iter, Fem::Ls::CLinearSystem_Eigen& ls);
double MinimumEigenValueVector_InvPower(
	LsSol::CLinearSystem& ls,
	LsSol::CPreconditioner& pls,
	const std::vector<unsigned int>& aIdVec, 
	unsigned int itr_invp,	// 逆べき乗法の最大反復回数
	unsigned int itr_lssol, // ICCG法の最大反復回数
	double conv_res_lssol,	// ICCG法の収束基準の相対残差
	int& iflag_conv );		// 0:正常に終了　1:ICCGが収束しなかった
}

#endif
