# compiler options
#--------------------------------------------
COMPILER ?= g++
mode ?= dynamic  # Default to dynamic linking

# TODO: -g -ggdb3 -fsanitize=address -fno-omit-frame-pointer
CXXFLAGS = -std=c++17 -O3
WFLAGS += -Wno-unused-result -Wno-unused-command-line-argument -Wno-unknown-pragmas -Wno-undefined-inline # -Wall

INC = -Iexternal/CLI11/include/CLI \
			-Iexternal/parallel-hashmap \
			-Iexternal/boost/libs/math/include

# project files
#--------------------------------------------
PROGRAM = krepp
TEST_PROGRAM = build/omp_index_test
OBJECTS = build/common.o \
					build/MurmurHash3.o build/lshf.o \
					build/phytree.o	build/rqseq.o \
					build/index.o build/sketch.o \
					build/query.o build/seek.o \
					build/record.o build/table.o \
					build/krepp.o
# Everything but the CLI layer, which is what an embedding library links against.
TEST_OBJECTS = $(filter-out build/krepp.o,$(OBJECTS)) build/omp_index_test.o
# -MMD -MP writes these next to each object; without them a change to a header
# rebuilds nothing, and `make test` reports green on code it did not compile.
DEPS = $(OBJECTS:.o=.d) build/omp_index_test.d

# rules
#--------------------------------------------
.PHONY: all dynamic static clean test

all:
	$(MAKE) mode=dynamic $(PROGRAM)

dynamic:
	$(MAKE) mode=dynamic $(PROGRAM)

static:
	$(MAKE) mode=static $(PROGRAM)

test:
	$(MAKE) mode=dynamic $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

# Check for -lcurl
CURL_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -lcurl -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)

$(info ===== Build mode: $(mode) =====)
ifeq ($(mode),dynamic)
	LDLIBS = -lm -lz -lstdc++
else ifeq ($(mode),static)
	LDLIBS = --static -static-libgcc -static-libstdc++ -lm -lz
	CURL_SUPPORTED = no
else
	LDLIBS = -lm -lz -lstdc++
endif

OS := $(shell uname -s)
ifneq ($(OS),Darwin)
	LDLIBS += -lstdc++fs
	LDOMP += -lgomp
else
	OMPFLAGS = -Xclang
	LDOMP += -lomp
endif
OMPFLAGS += -fopenmp

# Check for -lgomp
GOMP_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) $(LDFLAGS) $(CXXFLAGS) $(OMPFLAGS) $(LDOMP) -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)

WLCURL = 0
WOPENMP = 0
ifneq ($(CURL_SUPPORTED),no)
  ifneq ($(mode),static)
	  LDLIBS += -lcurl
	  WLCURL = 1
  endif
endif
ifneq ($(GOMP_SUPPORTED),no)
	LDLIBS += $(LDOMP)
	CXXFLAGS += $(OMPFLAGS)
	WOPENMP = 1
endif
VARDEF= -D _WLCURL=$(WLCURL) -D _WOPENMP=$(WOPENMP)

ARCH := $(shell uname -m)
# Check for -mbmi2
BMI2_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -mbmi2 -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)
ifeq ($(filter $(ARCH),x86_64 i386),$(ARCH))
	ifneq ($(BMI2_SUPPORTED),no)
		CXXFLAGS += -mbmi2
	endif
endif

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) -MMD -MP $(VARDEF) $(INC) -c src/$*.cpp -o build/$*.o $(LDLIBS) 

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(VARDEF) $(LDFLAGS) $(INC) -o $@ $(LDLIBS) 

# Named explicitly rather than as a second build/%.o pattern rule: with two
# pattern rules for the same target, make silently takes the first that has a
# prerequisite, so a test sharing a name with a src file would quietly compile
# the wrong one.
build/omp_index_test.o: test/omp_index_test.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) -MMD -MP $(VARDEF) $(INC) -c $< -o $@ $(LDLIBS) 

-include $(DEPS)

$(TEST_PROGRAM): $(TEST_OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(VARDEF) $(LDFLAGS) $(INC) -o $@ $(LDLIBS) 

clean:
	rm -f $(PROGRAM) $(OBJECTS) $(TEST_PROGRAM) build/omp_index_test.o $(DEPS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
