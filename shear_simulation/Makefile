# Define targets
EXE = a.out

# -D DIM=xx: for simulation of different dimension.
# -D__triBox__: triclinic box.
# -D__orthBox__: orthogonal box.
# -D Harmonic/Hertzian: for different interaction.
# -D equalMass/equalDensity: all particles have the same mass/density.
CFLAGS = -D DIM=3 -D__triBox__ -D Hertzian -D equalDensity

# setting
debug ?= false

ENABLE_GSL = false
GSL_INSTALL = $(shell gsl-config --prefix)

ENABLE_BLAS = false
BLAS_INSTALL = 

ENABLE_SPQR = false
SPQR_INSTALL = 

ENABLE_VORO = false
VORO_INSTALL= 

# Source Files
cSRC := $(wildcard *.c)
cOBJ := $(patsubst %.c,build/%.o,$(cSRC))
cDep := $(patsubst %.c,build/%.dep,$(cSRC))

cppSRC := $(wildcard *.cpp)
cppOBJ := $(patsubst %.cpp,build/%.o,$(cppSRC))
cppDep := $(patsubst %.cpp,build/%.dep,$(cppSRC))

SRC := $(cSRC) $(cppSRC)
OBJ := $(cOBJ) $(cppOBJ)
DEP := $(cDep) $(cppDep) 

SourceFileName = $(SRC) $(wildcard *.h) $(wildcard *.hpp)
CompilePath = $(shell pwd)

#debug
ifeq ($(debug),true)
	CC = gcc -std=c11
    CFLAGS += -g -O0 -D_GNU_SOURCE -D __SourceFileName__="\"$(SourceFileName)\"" -D __CompilePath__="\"$(CompilePath)"\" -MMD -MP -MF build/$*.dep
    CXX = g++
    LINK = g++
	LIB_PATH = 
	LIB =
	ifeq ($(ENABLE_BLAS),true)
		ifeq ($(strip $(BLAS_INSTALL)),)
			$(error Please set BLAS_INSTALL directory!)
		endif
		LIB_PATH += -DOpenblas -I $(BLAS_INSTALL)/inclue
		LIB += -L $(BLAS_INSTALL)/lib -lopenblas -lm
	endif
else
# Define compiler, BLAS
ifneq ($(shell which icc 2>/dev/null),)
    # use the Intel icc compiler for C codes, and -qopenmp for OpenMP
    CC = icc -std=c11
    CFLAGS += -O3 -D_GNU_SOURCE -D __SourceFileName__="\"$(SourceFileName)\"" -D __CompilePath__="\"$(CompilePath)"\" -MMD -MP -MF build/$*.dep
    CXX = icpc
    LINK = icpc
	LIB_PATH =
	LIB =
	ifeq ($(ENABLE_BLAS),true)
		LIB_PATH += -DMkl
		LIB += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
	else
		BLAS_INSTALL =
	endif
else
    CC = gcc -std=c11
    CFLAGS += -O3 -D_GNU_SOURCE -D __SourceFileName__="\"$(SourceFileName)\"" -D __CompilePath__="\"$(CompilePath)"\" -MMD -MP -MF build/$*.dep
    CXX = g++
    LINK = g++
	LIB_PATH = 
	LIB =
	ifeq ($(ENABLE_BLAS),true)
		ifeq ($(strip $(BLAS_INSTALL)),)
			$(error Please set BLAS_INSTALL directory!)
		endif
		LIB_PATH += -DOpenblas -I $(BLAS_INSTALL)/inclue
		LIB += -L $(BLAS_INSTALL)/lib -lopenblas -lm
	endif
endif
endif

ifeq ($(shell uname), Linux)
	CFLAGS += -D__Linux__
	objSourceFile = build/sourceFile.obj
	objTool = objcopy
else
	CFLAGS += -D__MacOS__
	objSourceFile = 
	objTool = 
endif

#GSL
ifeq ($(ENABLE_GSL),true)
	LIB_PATH += -I $(GSL_INSTALL)/include
	LIB += -L $(GSL_INSTALL)/lib -lgsl -lgslcblas -lm
	CFLAGS += -D__gslMathLib__
endif

#SparseQR
ifeq ($(ENABLE_SPQR),true)
	ifeq ($(strip $(SPQR_INSTALL)),)
			$(error Please set SPQR_INSTALL directory!)
	endif
	LIB_PATH += -DSpqr -I $(SPQR_INSTALL)/include/suitesparse/
	LIB += -L $(SPQR_INSTALL)/lib64/ -lcholmod -lcamd -lccolamd -lspqr -lsuitesparseconfig -lamd -lcolamd
endif

#voro++
ifeq ($(ENABLE_VORO),true)
	ifeq ($(strip $(VORO_INSTALL)),)
			$(error Please set SPQR_INSTALL directory!)
	endif
	LIB_PATH += -DVoroPlusPlus -I $(VORO_INSTALL)/include
	LIB += -L $(VORO_INSTALL)/lib -lvoro++
endif

# Define linker and its settings
CXXFLAGS = $(CFLAGS)

#define size command
SIZE = size

# Path to src files 
vpath %.c .
vpath %.h .
vpath %.cpp .

# default target
all: build $(EXE)

# make directory
build:
	@mkdir -p build

$(EXE): $(OBJ) $(objSourceFile)
	$(LINK) -o $(EXE) $(OBJ) $(objSourceFile) $(LIB)
	$(SIZE) $(EXE)

# Compilation rules
build/%.o: %.c
	$(CC) $(CFLAGS) $(LIB_PATH) -c $< -o $@

build/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LIB_PATH) -c $< -o $@

$(objSourceFile): $(SourceFileName)
	echo "" > build/sourceFile
	echo "//===Starting of Source Files===" >> build/sourceFile
	for fname in $(SourceFileName); do echo "//===Source File:" $$fname "===" >> build/sourceFile; cat $$fname >> build/sourceFile; done
	$(objTool) -I binary -B i386 -O elf64-x86-64 build/sourceFile $(objSourceFile)
	echo "//===Ending of Source Files===" >> build/sourceFile
	echo "" >> build/sourceFile
	rm build/sourceFile

sinclude $(DEP)

# Clean command
clean:
	rm -rf $(OBJ)
	rm -rf $(DEP)
	rm -rf $(objSourceFile)
	rm -rf $(EXE)

# Phony Target
.PHONY: all clean build

