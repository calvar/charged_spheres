BIN		:= chSphere

INCLUDES 	+= -I. 
LIBS		:= 

CXXFLAGS	:= -O3 -m64 #-g -fno-inline 
LDFLAGS		:= -lgsl -lgslcblas -lgomp

CXX		:= g++
LINKER		:= g++

C_SOURCES	:= main.cpp force.cpp particles.cpp
HEADERS		:= force.hpp particles.hpp 

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))


$(BIN): $(C_OBJS) $(HEADERS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CXX) -c $(C_SOURCES) $(CXXFLAGS) $(INCLUDES) #-fopenmp

clean: 
	rm $(BIN) *.o core*
