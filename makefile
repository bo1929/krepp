# compiler options
#--------------------------------------------
COMPILER ?= g++
mode ?= dynamic  # Default to dynamic linking

CXXFLAGS += -std=c++17 -O3 # TODO: remove -g
WFLAGS += -Wno-unused-result -Wno-unused-command-line-argument -Wno-unknown-pragmas -Wno-undefined-inline # -Wall

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
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(LDLIBS) $(VARDEF) $(INC) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(VARDEF) $(LDFLAGS) $(INC) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
