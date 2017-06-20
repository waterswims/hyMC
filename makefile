# Use intel compiler
CXX=icpc

# Specify subdirectories for source and objects
INC_PATH=include
OBJ_PATH=objects
LIB_PATH=lib

# C flags
override CXX+=--std=c++11 -W -Wall -pedantic -O3 -g -DMKL_ILP64 -I$(MKLROOT)/include
LDLIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

# Collect all the library source files from the library director
LIB_SOURCES=$(wildcard $(LIB_PATH)/*.cpp)

# Generate the names for each corresponding object file
OBJ_FILES=$(addprefix $(OBJ_PATH)/,$(notdir $(LIB_SOURCES:.cpp=.o)))

# Default target builds the static hymc library
default: hymc.a

hymc.a: $(OBJ_FILES)
	ar rcs $@ $^

# Build the object files
$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CXX)	$(CXXFLAGS) $(LDFLAGS) -c \
		-o $@ $<

clean:
	rm hymc.a
	rm $(OBJ_FILES)
