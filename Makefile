CXX    = g++
CFLAGS = -O2
LDFLAGS =
INCLUDES = -I./DelFEM/include
LIBS = -L./DelFEM/lib -ldfm

TARGET = main.out
ifeq ($(OS),Windows_NT) 
	LIBS_GL = -lfreeglut -lglu32 -lopengl32		#Windows
	TARGET = main.exe	
else ifeq ($(shell uname -s),Darwin)	
	LIBS_GL = -framework OpenGL -framework GLUT	#Mac
else	
	LIBS_GL = -lglut -lGL -lGLU			#Linux(defalut)
endif
OBJS = main.o \
	analysis2d_cloth_static.o \
	cloth_utility.o \
	contact_target_mesh.o 	contact_target_adf.o \
	cloth_handler.o \
	designer2d_cloth.o \
	emat_quad_bend.o eqn_quad_bend.o \
	eqn_glue.o \
	eqn_contact3d.o \
	emat_glue.o \
	emat_cst.o \
	stitch_array.o \

all: $(TARGET)
					
$(TARGET): $(OBJS) libdelfem.a
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(LIBS_GL)

clean:
	-rm -f $(OBJS) DelFEM/lib/libdelfem.a
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

libdelfem.a:
	make -C DelFEM
