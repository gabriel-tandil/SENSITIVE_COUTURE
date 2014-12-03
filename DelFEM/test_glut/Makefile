MAKE = mingw32-make --no-print-directory

SUBDIRS = cad_view msh_view field_evaluate scalar2d scalar3d solid2d solid3d fluid2d dkt2d helmholtz2d hyper3d rigid fieldrigid

all :
	$(MAKE) -C cad_view
	$(MAKE) -C msh_view
	$(MAKE) -C field_evaluate
	$(MAKE) -C scalar2d
	$(MAKE) -C scalar3d
	$(MAKE) -C solid2d
	$(MAKE) -C solid3d
	$(MAKE) -C fluid2d
	$(MAKE) -C dkt2d
	$(MAKE) -C helmholtz2d
	$(MAKE) -C hyper3d
	$(MAKE) -C rigid
	$(MAKE) -C stvk_rigid

clean :
	$(MAKE) clean -C cad_view
	$(MAKE) clean -C msh_view
	$(MAKE) clean -C field_evaluate
	$(MAKE) clean -C scalar2d
	$(MAKE) clean -C scalar3d
	$(MAKE) clean -C solid2d
	$(MAKE) clean -C solid3d
	$(MAKE) clean -C fluid2d
	$(MAKE) clean -C dkt2d
	$(MAKE) clean -C helmholtz2d
	$(MAKE) clean -C hyper3d
	$(MAKE) clean -C rigid
	$(MAKE) clean -C stvk_rigid


