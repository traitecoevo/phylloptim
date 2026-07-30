# Build and run the leaf test suite as plain C++ -- no R, no linking.
#
#   make          build and run
#   make build    build only
#   make clean
#
# Header locations are discovered from R if it is on the PATH (the normal case
# once BH and odelia are installed), and can be overridden:
#
#   make BH_INC=/opt/homebrew/include ODELIA_INC=../../../odelia/inst/include
#
# Falling back to a sibling odelia checkout means this works in the plant-family
# tree before anything is installed.

CXX      ?= c++
CXXSTD   ?= -std=c++20
OPT      ?= -O2
WARN     ?= -Wall -Wextra -Wno-unused-parameter

LEAF_INC   ?= ../../inst/include
SHIM_INC   ?= shim
BH_INC     ?= $(shell Rscript -e 'cat(system.file("include", package="BH"))' 2>/dev/null)
ODELIA_INC ?= $(shell Rscript -e 'cat(system.file("include", package="odelia"))' 2>/dev/null)

ifeq ($(strip $(BH_INC)),)
  BH_INC := /opt/homebrew/include
endif
ifeq ($(strip $(ODELIA_INC)),)
  ODELIA_INC := ../../../odelia/inst/include
endif

INCLUDES = -I$(LEAF_INC) -isystem $(SHIM_INC) -isystem $(ODELIA_INC) -isystem $(BH_INC)

.PHONY: all build clean
all: build
	./test_leaf

build: test_leaf

test_leaf: test_leaf.cpp $(wildcard $(LEAF_INC)/leaf/*.hpp) $(LEAF_INC)/leaf.hpp
	$(CXX) $(CXXSTD) $(OPT) $(WARN) $(INCLUDES) $< -o $@

clean:
	rm -f test_leaf
