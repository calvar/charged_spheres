BIN		:= md_ch_sphere

CXXFLAGS	:= -O3 -m64 #-g -fno-inline
LDFLAGS		:= -lgsl -lgslcblas

CXX		:= g++
LINKER		:= g++

C_SOURCES	:= main.cpp particles.cpp initialize.cpp print_read.cpp force.cpp eq_of_motion.cpp timer.cpp
HEADERS		:= particles.hpp initialize.hpp print_read.hpp force.hpp eq_of_motion.hpp timer.hpp

C_OBJS		:= $(patsubst %.cpp, %.o, $(C_SOURCES))


$(BIN): $(C_OBJS) $(CU_OBJS) $(HEADERS)
	$(LINKER) -o $@ $(C_OBJS) $(LDFLAGS)

clean: 
	rm $(BIN) *.o core*
