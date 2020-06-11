MMADIR = /usr/local/Wolfram/Mathematica/11.3/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
MPREP = $(MMADIR)/mprep
CXX = g++
DELETE = rm -f

CXXFLAGS = -I../../libiint/include -I${HOME}/.usr/include -I${HOME}/.usr/include/flint -I$(MMADIR)
LDFLAGS = $(CXXFLAGS) -L../../libiint -L$(HOME)/.usr/lib -L$(MMADIR) -liint -lflint -larb -lgmp -lML64i4 -lrt -luuid -ldl

SRC = *.cpp
OBJ = $(patsubst %.cpp,%.o,$(wildcard $(SRC))) MathLink.o
OUT = IInt

all:	$(OUT)

clean:
	$(DELETE) $(OUT)
	$(DELETE) $(OBJ)
	$(DELETE) MathLink.cpp

$(OUT):	$(OBJ) ../../libiint/libiint.a
	$(CXX) $(OBJ) -o $(OUT) $(LDFLAGS) 

%.o:	%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) 

%.cpp:	%.tm
	$(MPREP) $< -o $@

