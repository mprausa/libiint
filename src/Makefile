# Makefile

CXX = g++

CXXFLAGS = -I../include -I${HOME}/.usr/include/flint -g -std=gnu++11
LDFLAGS=$(CXXFLAGS) -L${HOME}/.usr/lib -rdynamic -lflint -larb -lgmp -ldl -lyaml-cpp
DELETE = rm -f

SRC = *.cpp arb/*.cpp

OBJ = $(patsubst %.cpp,%.o,$(wildcard $(SRC)))
DEP = $(patsubst %.cpp,%.cpp.dep,$(wildcard $(SRC)))

DEPEND = .depend

OUT = iint

all:	$(OUT)

dep:	$(DEPEND)
	$(DELETE) $(DEP)

clean:
	$(DELETE) $(OBJ)
	$(DELETE) $(OUT)
	$(DELETE) $(DEPEND)

$(OUT):	$(OBJ) 
	$(CXX) $(OBJ) -o $(OUT) $(LDFLAGS) 

$(DEPEND):	$(DEP)
	cat $(DEP) > $(DEPEND)

%.o:	%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

%.cpp.dep:	%.cpp
	$(CXX) -c -M -MT $*.o $< $(CXXFLAGS) > $@

ifneq ($(wildcard $(DEPEND)),)
include $(shell /bin/pwd)/$(DEPEND)
endif
