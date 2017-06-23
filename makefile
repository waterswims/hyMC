# Use intel compiler
CXX=icpc

# Specify subdirectories for source and objects
INC_PATH=include
OBJ_PATH=objects
LIB_PATH=lib
TEST_PATH=tests

GTEST_DIR=tests/googletest/googletest
GTEST_FLAGS=-isystem $(GTEST_DIR)/include

# C flags
CXXFLAGS=--std=c++11 -W -Wall -pedantic -O3 -g -DMKL_ILP64 -I$(MKLROOT)/include -use-intel-optimized-headers
LDLIBS=-use-intel-optimized-headers -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

# Collect all the library source files from the library director
LIB_SOURCES=$(wildcard $(LIB_PATH)/*.cpp)

# Generate the names for each corresponding object file
OBJ_FILES=$(addprefix $(OBJ_PATH)/,$(notdir $(LIB_SOURCES:.cpp=.o)))

# All of the test files that are needed
TEST_FILES=$(wildcard $(TEST_PATH)/*.hpp)

# Default target builds the static hymc library
default: hymc.a

hymc.a: $(OBJ_FILES)
	ar rcs 	$@ $^

tests: runtests

runtests: $(TEST_PATH)/test.cpp hymc.a $(TEST_PATH)/libs/gtest_main.a $(TEST_FILES)
	$(CXX) 	$(GTEST_FLAGS) $(CXXFLAGS) $< $(TEST_PATH)/libs/gtest_main.a \
		$(OBJ_FILES) \
		-o $@ \
		$(LDLIBS) hymc.a

# Build the object files
$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CXX)	$(CXXFLAGS) -c \
		-o $@ $<

# GTEST build
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
$(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

tests/libs/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) 	$(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ \
		$(GTEST_DIR)/src/gtest-all.cc

tests/libs/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) 	$(GTEST_FLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
		-o $@ \
		$(GTEST_DIR)/src/gtest_main.cc

tests/libs/gtest.a : tests/libs/gtest-all.o
	$(AR) 	$(ARFLAGS) $@ $^

tests/libs/gtest_main.a : tests/libs/gtest-all.o tests/libs/gtest_main.o
	$(AR) 	$(ARFLAGS) $@ $^

clean:
	rm -f hymc.a runtests
	rm -f $(OBJ_FILES)
	rm -f $(TEST_PATH)/libs/*
