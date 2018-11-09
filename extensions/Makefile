#===============================================================================
# Compiler and compiler flags.
#===============================================================================

CXXFLAGS += -c -std=c++11 -O3 -march=native -fopenmp
#CXXFLAGS += -O0 -g -Wall -pedantic -fbounds-check
LDDFLAGS += -fopenmp

#===============================================================================
# Add Python includes.
#===============================================================================

PY_CXXFLAGS = $(shell python3-config --includes)
PY_CXXFLAGS += $(shell python3 -c "import numpy; print('-I' + numpy.get_include())")

CXXFLAGS += $(PY_CXXFLAGS)

#===============================================================================
# Specify build directory and source files.
#===============================================================================

OBJ_DIR = build

SOURCES = backend.cpp
OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

#===============================================================================
# Primary rules.
#===============================================================================

all : my_extension.so

my_extension.so : $(OBJECTS)
	$(CXX) $(LDDFLAGS) -shared $(OBJECTS) -o my_extension.so

clean:
	@rm $(OBJ_DIR)/*.o *.so

.PHONY : all clean

#===============================================================================
# Other rules.
#===============================================================================

# Compiler rule.
$(OBJ_DIR)/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -fPIC -c $< -o $@

# These rules require that the build directory is existant before compiling.
$(OBJECTS) : | $(OBJ_DIR)

# Build directory rule.
$(OBJ_DIR) :
	mkdir $(OBJ_DIR)
