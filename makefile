# compiler options
#--------------------------------------------
COMPILER ?= g++ # g++-14
mode ?= dynamic  # Default to dynamic linking

CXXFLAGS = -std=c++17 -O3 # -g
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument -Wno-unknown-pragmas -Wno-undefined-inline # -Wall

INC = -Iexternal/CLI11/include/CLI \
			-Iexternal/parallel-hashmap \
			-Iexternal/boost/libs/math/include

# project files
#--------------------------------------------
PROGRAM = krepp
OBJECTS = build/common.o \
					build/MurmurHash3.o build/lshf.o \
					build/phytree.o	build/rqseq.o \
					build/index.o build/sketch.o \
					build/query.o build/seek.o \
					build/record.o build/table.o \
					build/krepp.o

# rules
#--------------------------------------------
.PHONY: all dynamic static clean

all:
	$(MAKE) mode=dynamic $(PROGRAM)

dynamic:
	$(MAKE) mode=dynamic $(PROGRAM)

static:
	$(MAKE) mode=static $(PROGRAM)

# Check for -lcurl and -lgomp
CURL_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -lcurl -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)
GOMP_SUPPORTED := $(shell echo 'int main() { return 0; }' | $(COMPILER) -fopenmp -lgomp -x c++ -o /dev/null - 2>/dev/null && echo yes || echo no)

$(info ===== Build mode: $(mode) =====)
ifeq ($(mode),dynamic)
	LDLIBS = -lstdc++fs -lm -lz -lstdc++
else ifeq ($(mode),static)
	LDLIBS = --static -lstdc++fs -lm -lz -static-libgcc -static-libstdc++
	CURL_SUPPORTED = no
else
	LDLIBS = -lstdc++fs -lm -lz -lstdc++
endif

WLCURL = 0
WOPENMP = 0
ifneq ($(CURL_SUPPORTED),no)
	LDLIBS += -lcurl
	WLCURL = 1
endif
ifneq ($(GOMP_SUPPORTED),no)
	LDLIBS += -lgomp
	CXXFLAGS += -fopenmp
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
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(LDLIBS) $(VARDEF) $(INC) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(VARDEF) $(LDFLAGS) $(INC) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
