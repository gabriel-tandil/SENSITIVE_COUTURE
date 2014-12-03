
#if !defined(EQN_HYPER_H)
#define EQN_HYPER_H

#if defined(__VISUALC__)
    #pragma warning( disable : 4786 )
#endif

namespace Fem{
namespace Eqn{
	class ILinearSystem_Eqn;
}
namespace Field{
	class CField;
	class CFieldWorld;
}
namespace Eqn{
    // Static
	bool AddLinSys_Hyper2D_Static
	(Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y,
	 const unsigned int id_field_disp, const unsigned int id_field_lambda,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_ea = 0);

	// Static
	bool AddLinSys_Hyper3D_Static
	(Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,
	 const unsigned int id_field_disp, const unsigned int id_field_lambda,
	 const Fem::Field::CFieldWorld& world,
	 unsigned int id_ea = 0 );
		
	// Dynamic
	bool AddLinSys_Hyper3D_NonStatic_NewmarkBeta
	(double dt, double gamma, double beta,
	 Fem::Eqn::ILinearSystem_Eqn& ls,
	 double lambda, double myu,
	 double  rho, double g_x, double g_y, double g_z,    
	 const unsigned int id_field_disp, const unsigned int id_field_lambda,
	 const Fem::Field::CFieldWorld& world, 
	 bool is_initial,
	 unsigned int id_ea = 0 );
}
}

#endif
