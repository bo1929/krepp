# compiler options
#--------------------------------------------
COMPILER = g++
WLCURL = 1
WOPENMP = 1
CSTATIC = 0

CXXFLAGS = -std=c++17 -O3

ifeq ($(CSTATIC), 0)
	LDLIBS = -lstdc++fs -lm -lz -lstdc++
else
	LDLIBS = --static -lstdc++fs -lm -lz -static-libgcc -static-libstdc++
endif

ifneq ($(WLCURL), 0)
	LDLIBS += -lcurl
endif

ifneq ($(WOPENMP), 0)
	LDLIBS += -lgomp
	CXXFLAGS += -fopenmp
endif

VARDEF= -D _WLCURL=$(WLCURL) -D _WOPENMP=$(WOPENMP)

INC = -Iexternal/CLI11/include/CLI \
			-Iexternal/parallel-hashmap \
			-Iexternal/boost/libs/math/include

WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument -Wno-unknown-pragmas

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
all: $(PROGRAM)

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
