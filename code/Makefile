CXX 	= g++
CFLAGS  = -c -O2 -W -g
INCLUDE = -I./include -I/usr/include/suitesparse
LIBPATH = -L./lib
LDFLAGS = -lm -ldl -lfftw3 -lumfpack
OBJECTS	= MainCollect Data BasicOperators DualBox AuxInput FftwToolkit StableBulkPhase \
		  GenJacPoly ProjFourGJP CommonRotateProjBoxMat CommonFourGJP SysInitVal \
		  SparseOperators IterPrepare IterMethod DisplayResults MemFree

all:interface
./lib/lib%.so:./src/%.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(INCLUDE) $(LDFLAGS)

interface:$(addsuffix .so, $(addprefix ./lib/lib, $(OBJECTS))) interface.cpp
		$(CXX) -g $(INCLUDE) interface.cpp -o $@ $(LIBPATH) $(addprefix -l, $(OBJECTS)) $(LDFLAGS)
clean: 
	-rm ./lib/*.so
	-rm interface
#	-rm -r ./result
