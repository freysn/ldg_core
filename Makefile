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
endif


INCLUDES=-I.

###

dbg ?= 0
ndbg ?= 0

PYTHON ?= python3.8

DARWIN_EXTRA_FLAGS=

LIB = -lbz2 -L/opt/local/lib -lpng -lcairo
MHOME=/home/freysn


ifneq ($(DARWIN),)

CXX=/opt/local/bin/clang++-mp-11 -lc++
PYTHON = python3.9
DARWIN_EXTRA_FLAGS += -undefined dynamic_lookup -mlinker-version=450
INCLUDES += -I/opt/local/include

else
PYTHON = python3
INCLUDES += -I$(MHOME)/.local/include/python3.8

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

#CXX_FLAGS+=-DUSE_STRAX_MCMC

PY_DIR=py/
EXEC_NAME=ldg_core

all:
	$(CXX) $(CXX_FLAGS) ${DARWIN_EXTRA_FLAGS} $(INCLUDES) supertiles_place.cpp -o ${EXEC_NAME} $(LIB)

pb:
	$(CXX) $(CXX_FLAGS) -shared -fPIC ${DARWIN_EXTRA_FLAGS} `${PYTHON} -m pybind11 --includes` ${PY_DIR}supertiles_pybind11.cpp -o ${PY_DIR}supertiles`${PYTHON}-config --extension-suffix` -I. ${INCLUDES} -I/opt/local/include -I$(MHOME)/dev/ext/pybind11/include $(LIB)

pb_test:
	$(CXX) $(CXX_FLAGS) ${DARWIN_EXTRA_FLAGS} ${PY_DIR}supertiles_pybind11.cpp -o ${PY_DIR}supertiles_pb_test -I. ${INCLUDES} -I/opt/local/include -I$(MHOME)/dev/ext/pybind11/include $(LIB) -DJUST_A_TEST

clean:
	rm -f *~ ${EXEC_NAME}
