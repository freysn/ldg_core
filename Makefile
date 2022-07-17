### include common.mk

OSUPPER = $(shell uname -s 2>/dev/null | tr "[:lower:]" "[:upper:]")
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))
LIB =
#LIB += -LvolData 
ifneq ($(DARWIN),)
#ORTOOLS_DIR = /Users/freysn/dev/ext/or-tools.MacOsX64
ORTOOLS_DIR = /Users/freysn/dev/ext/ortools
#LIB += -lVolData_libc++
else
#ORTOOLS_DIR = /home/freysn/dev/ext/or-tools.Linux64
ORTOOLS_DIR = $(HOME)/dev/ext/ortools
#LIB += -lVolData
endif

HOST=$(shell hostname)

#LIB += -L/opt/local/lib


LIB += -lpthread
#-lpcl_kdtree

CXX_FLAGS = -std=c++11
ifeq ($(dbg),1)
CXX_FLAGS += -g
else
CXX_FLAGS += -O4
endif
#CXX_FLAGS += -DVERBOSE
CXX_FLAGS += -Wall -Wno-reorder -Wno-delete-non-virtual-dtor -Wno-sign-compare

# ifeq ($(DARWIN),)
# CXX_FLAGS += -fopenmp
# CXX_FLAGS += -DUSE_OMP_INTERVOL_INTERVOL
# endif

ifneq ($(no_ext_libs),1)
CXX_FLAGS += -Dcimg_use_png
ifeq ($(DARWIN),)
#LIB += -lboost_program_options
#-lboost_system -lboost_thread -lboost_program_options -lboost_serialization -lboost_chrono -lboost_filesystem 
else
#LIB += -lboost_program_options-mt
#LIB += -lboost_program_options
#-lboost_system-mt -lboost_thread-mt -lboost_program_options-mt -lboost_serialization-mt -lboost_chrono-mt -lboost_filesystem-mt 
endif
endif

#INCLUDES = -I/Developer/NVIDIA/CUDA-6.5/samples/common/inc -I/Developer/NVIDIA/CUDA-6.5/include -I/opt/local/include -I/usr/local/cuda/samples/common/inc -I/usr/local/cuda/include -I.. -I/opt/local/include/pcl-1.7/ -I/opt/local/include/eigen3 -I.

INCLUDES=-I.
#INCLUDES+=-I/opt/local/include

#LIB += -L/usr/local/opt/boost155/lib
#INCLUDES += -I/usr/local/opt/boost155/include

###

dbg ?= 0
ndbg ?= 0

PYTHON ?= python3.8

DARWIN_EXTRA_FLAGS=

LIB = -lbz2 -L/opt/local/lib -lpng -lcairo
MHOME=/home/freysn


ifneq ($(DARWIN),)
CXX=/opt/local/bin/clang-mp-9.0 -lc++
#CXX=/opt/local/bin/clang++-mp-10 -lc++
#ifeq "$(HOST)" "dupin.local"
CXX=/opt/local/bin/clang++-mp-11 -lc++
#CXX=clang -lc++
PYTHON = python3.9
#endif
DARWIN_EXTRA_FLAGS += -undefined dynamic_lookup -mlinker-version=450
#CXX=clang
INCLUDES += -I/opt/local/include
#INCLUDES += -I/Developer/NVIDIA/CUDA-10.2/include
else ifeq "$(HOST)" "zydeco"
CXX=g++-8
else ifeq "$(HOST)" "acappella.visus.uni-stuttgart.de"
PYTHON=python3
INCLUDES += -I$(MHOME)/dev/ext/pybind11/include
else
#ifeq "$(HOST)" "blackmetal01.visus.uni-stuttgart.de"
#MHOME=/tmp/freysn
#endif
CXX=g++
ifeq "$(HOST)" "piet"
PYTHON = python3
INCLUDES += -I$(MHOME)/.local/include/python3.8
else
INCLUDES += -I$(MHOME)/dev/ext/cimg
INCLUDES += -I$(MHOME)/dev/ext/pybind11/include
INCLUDES += -I$(MHOME)/dev/ext/libpng/include
INCLUDES += -I$(MHOME)/dev/ext/bzip2/include
INCLUDES += -I$(MHOME)/dev/ext/cairo/include
LIB += -L$(MHOME)/dev/ext/libpng/lib
LIB += -L$(MHOME)/dev/ext/bzip2/lib
LIB += -L$(MHOME)/dev/ext/cairo/lib
PYTHON = $(MHOME)/dev/ext/py3/bin/python3.8
endif
endif

CXX_FLAGS=-std=c++17 -Wall
CXX_FLAGS+=-fopenmp
#CXX_FLAGS+=-DNO_OMP
ifeq ($(dbg),1)
CXX_FLAGS += -g
else
ifneq ($(DARWIN),)
CXX_FLAGS += -O2
else
CXX_FLAGS += -O3
endif
endif

ifeq ($(ndbg),1)
CXX_FLAGS += -DNDEBUG
endif

# print dependencies instead of compiling or preprocessing. If you use it without other options, it includes all dependencies, including system headers
# lter by path, or use some options like -MF instead of plain -M if you want to find what subset of a large source tree is needed 
#CXX_FLAGS+=-M

#CXX_FLAGS+=-DUSE_STRAX_MCMC

EXEC_NAME=ldg_core

all:
	$(CXX) $(CXX_FLAGS) ${DARWIN_EXTRA_FLAGS} $(INCLUDES) supertiles_place.cpp -o ${EXEC_NAME} $(LIB)

pb:
	$(CXX) $(CXX_FLAGS) -shared -fPIC ${DARWIN_EXTRA_FLAGS} `${PYTHON} -m pybind11 --includes` supertiles_pybind11.cpp -o py/supertiles`${PYTHON}-config --extension-suffix` -I. ${INCLUDES} -I/opt/local/include -I$(MHOME)/dev/ext/pybind11/include $(LIB)
	# -Wall

clean:
	rm -f *~ ${EXEC_NAME}
