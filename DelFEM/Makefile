CXX    = g++
CXXFLAGS = -O2 
LDFLAGS = 
INCLUDES = -Iinclude
ifeq ($(OS),Windows_NT)
	LIBS_GL = -lglu32 -lopengl32	# on Windows
else
	LIBS_GL = -lGL -lGLU	
endif

OBJS = drawer.o drawer_gl_utility.o quaternion.o uglyfont.o vector3d.o \
	spatial_hash_grid3d.o spatial_hash_grid2d.o tri_ary_topology.o \
	cad_obj2d.o cad_obj2d_move.o \
	cad_elem2d.o cad_edge2d_polyline.o \
	drawer_cad.o brep.o brep2d.o\
	drawer_msh.o mesh3d.o mesher2d.o mesher2d_edit.o \
	surface_mesh_reader.o mesher3d.o\
	meshkernel2d.o meshkernel3d.o mesh3d_extrude.o\
	drawer_field.o drawer_field_face.o drawer_field_edge.o drawer_field_vector.o elem_ary.o eval.o field.o field_world.o node_ary.o\
	mat_blkcrs.o matdia_blkcrs.o matdiafrac_blkcrs.o matdiainv_blkdia.o matfrac_blkcrs.o matprolong_blkcrs.o ordering_blk.o solver_mg.o solver_mat_iter.o vector_blk.o\
	zmat_blkcrs.o zmatdia_blkcrs.o zmatdiafrac_blkcrs.o zsolver_mat_iter.o zvector_blk.o\
	linearsystem.o preconditioner.o solver_ls_iter.o\
	linearsystem_field.o linearsystem_fieldsave.o zlinearsystem.o zsolver_ls_iter.o\
	rigidbody.o linearsystem_rigid.o linearsystem_rigidfield.o \
	eqn_advection_diffusion.o eqn_diffusion.o eqn_dkt.o eqn_helmholtz.o eqn_linear_solid2d.o eqn_linear_solid3d.o eqn_navier_stokes.o eqn_poisson.o eqn_stokes.o eqn_st_venant.o eqn_hyper.o\
	eqnsys.o eqnsys_fluid.o eqnsys_scalar.o eqnsys_shell.o eqnsys_solid.o ker_emat_tri.o

VPATH = src/com src/cad src/msh src/femfield\
	src/matvec src/femls src/femeqn src/femeqn src/ls src/rigid\
	lib\

all : lib/libdfm.a

lib/libdfm.a: $(patsubst %,lib/%,$(OBJS))
	ar r $@ $^
	ranlib lib/libdfm.a

lib/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -rf lib/*.o lib/*.a
