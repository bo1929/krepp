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
UNIT_PROGRAM = build/krepp_tests
LIB_OBJECTS = $(filter-out build/krepp.o,$(OBJECTS))
TEST_OBJECTS = $(LIB_OBJECTS) build/omp_index_test.o
UNIT_SOURCES = $(wildcard test/unit/*.cpp)
UNIT_OBJECTS = $(patsubst test/unit/%.cpp,build/test/unit/%.o,$(UNIT_SOURCES))
DEPS = $(OBJECTS:.o=.d) build/omp_index_test.d $(UNIT_OBJECTS:.o=.d)

# rules
#--------------------------------------------
.PHONY: all dynamic static clean test test-unit test-regression test-threads coverage

all:
	$(MAKE) mode=dynamic $(PROGRAM)

dynamic:
	$(MAKE) mode=dynamic $(PROGRAM)

static:
	$(MAKE) mode=static $(PROGRAM)

test: test-unit test-regression test-threads

test-unit: $(UNIT_PROGRAM)
	./$(UNIT_PROGRAM)

test-regression: $(PROGRAM)
	bash test/regression/run_regression.sh ./$(PROGRAM)

test-threads: $(TEST_PROGRAM)
	./$(TEST_PROGRAM)

# Source-based coverage of the unit/integration suite (clang only).
# The report covers src/ only; the tests themselves are excluded.
COVERAGE_DIR = build/coverage
coverage:
	$(MAKE) clean
	$(MAKE) COVERAGE=yes $(UNIT_PROGRAM)
	rm -rf $(COVERAGE_DIR)
	mkdir -p $(COVERAGE_DIR)
	LLVM_PROFILE_FILE=$(COVERAGE_DIR)/krepp-%p.profraw ./$(UNIT_PROGRAM)
	xcrun llvm-profdata merge -sparse $(COVERAGE_DIR)/*.profraw -o $(COVERAGE_DIR)/krepp.profdata
	xcrun llvm-cov report ./$(UNIT_PROGRAM) -instr-profile=$(COVERAGE_DIR)/krepp.profdata \
		-sources src -ignore-filename-regex='(test|external)/'
	xcrun llvm-cov export ./$(UNIT_PROGRAM) -instr-profile=$(COVERAGE_DIR)/krepp.profdata \
		-sources src -ignore-filename-regex='(test|external)/' -format=lcov > $(COVERAGE_DIR)/lcov.info
	@echo "lcov report: $(COVERAGE_DIR)/lcov.info"

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
	OMPFLAGS = -fopenmp
	LDOMP = -lgomp
else
	OMPFLAGS = -Xclang -fopenmp
	OMP_PREFIX := $(firstword $(wildcard /opt/homebrew/opt/libomp /usr/local/opt/libomp))
	ifneq ($(OMP_PREFIX),)
		OMP_INC = -I$(OMP_PREFIX)/include
		OMP_LIB = -L$(OMP_PREFIX)/lib
	else ifneq ($(wildcard $(CONDA_PREFIX)/include/omp.h),)
		OMP_PREFIX = $(CONDA_PREFIX)
		OMP_INC = -I$(CONDA_PREFIX)/include
		OMP_LIB = -L$(CONDA_PREFIX)/lib
	endif
	LDOMP = $(OMP_LIB) -lomp
endif

# Check for -lgomp
GOMP_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) $(LDFLAGS) $(CXXFLAGS) -g0 $(OMPFLAGS) $(OMP_INC) $(LDOMP) -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)

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
	CXXFLAGS += $(OMPFLAGS) $(OMP_INC)
	WOPENMP = 1
endif
VARDEF= -D _WLCURL=$(WLCURL) -D _WOPENMP=$(WOPENMP)
$(info ===== OpenMP: $(WOPENMP)$(if $(OMP_PREFIX), (libomp at $(OMP_PREFIX)),) =====)

# Optional source-based coverage instrumentation for `make coverage`.
ifeq ($(COVERAGE),yes)
	CXXFLAGS += -fprofile-instr-generate -fcoverage-mapping
	LDFLAGS += -fprofile-instr-generate -fcoverage-mapping
endif

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

build/omp_index_test.o: test/omp_index_test.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) -MMD -MP $(VARDEF) $(INC) -c $< -o $@ $(LDLIBS) 

# doctest suite: links the whole library (everything but the CLI's main()).
build/test/unit/%.o: test/unit/%.cpp
	@mkdir -p build/test/unit
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) -MMD -MP $(VARDEF) $(INC) -Iexternal/doctest -Isrc -c $< -o $@ $(LDLIBS) 

-include $(DEPS)

$(UNIT_PROGRAM): $(UNIT_OBJECTS) $(LIB_OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(VARDEF) $(LDFLAGS) $(INC) -o $@ $(LDLIBS) 

$(TEST_PROGRAM): $(TEST_OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(VARDEF) $(LDFLAGS) $(INC) -o $@ $(LDLIBS) 

clean:
	rm -f $(PROGRAM) $(OBJECTS) $(TEST_PROGRAM) build/omp_index_test.o $(UNIT_PROGRAM) $(UNIT_OBJECTS) $(DEPS)
	@find build -type d -empty -delete 2>/dev/null || true
	@echo "Succesfully cleaned."
