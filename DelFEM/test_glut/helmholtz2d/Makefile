CXX    = g++
CFLAGS = -Wall -O2
LDFLAGS =
INCLUDES = -I../../include
LIBS = -L../../lib -ldfm

TARGET = main.out
ifeq ($(OS),Windows_NT) 
	LIBS_GL = -lfreeglut -lglu32 -lopengl32		#Windows
	TARGET = main.exe	
else ifeq ($(shell uname -s),Darwin)	
	LIBS_GL = -framework OpenGL -framework GLUT	#Mac
else	
	LIBS_GL = -lglut -lGL -lGLU			#Linux(defalut)
endif
OBJS = main.o

all: $(TARGET)
					
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(LIBS_GL)

clean:
	-rm -f $(OBJS)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<