# Makefile

CXX = g++
AR = ar

CXXFLAGS = -Iinclude -I${HOME}/.usr/include -I${HOME}/.usr/include/flint -g -std=gnu++11 -fPIC -g -fmax-errors=5
ARFLAGS = rcs
DELETE = rm -f

SRC = src/iint/*.cpp src/arb/*.cpp

OBJ = $(patsubst %.cpp,%.o,$(wildcard $(SRC)))
DEP = $(patsubst %.cpp,%.cpp.dep,$(wildcard $(SRC)))

DEPEND = .depend
OUT = libiint.a

all:	$(OUT)

dep:	$(DEPEND)
	$(DELETE) $(DEP)

clean:
	$(DELETE) $(OBJ)
	$(DELETE) $(OUT)
	$(DELETE) $(DEPEND)

$(OUT):	$(OBJ) 
	$(AR) $(ARFLAGS) $(OUT) $(OBJ)

$(DEPEND):	$(DEP)
	cat $(DEP) > $(DEPEND)

%.o:	%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

%.cpp.dep:	%.cpp
	$(CXX) -c -M -MT $*.o $< $(CXXFLAGS) > $@

ifneq ($(wildcard $(DEPEND)),)
include $(shell /bin/pwd)/$(DEPEND)
endif
