# compiler options
#--------------------------------------------
WLCURL=1
COMPILER = g++
ifeq ($(WLCURL),0)
	LDLIBS= -lstdc++fs -lm -lz -lstdc++
	VARDEF= -D WLCURL=$(WLCURL)
else
	LDLIBS= -lstdc++fs -lm -lz -lstdc++ -lcurl
	VARDEF= 
endif
INC = -Iexternal/CLI11/include/CLI \
			-Iexternal/parallel-hashmap -Iexternal/boost/libs/math/include
CXXFLAGS = -std=c++17 -O3 -g -fopenmp
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument

# project files
#--------------------------------------------
PROGRAM = krepp
OBJECTS = build/common.o \
					build/MurmurHash3.o \
					build/lshf.o build/rqseq.o \
					build/index.o build/sketch.o \
					build/query.o build/compare.o \
					build/record.o build/phytree.o build/table.o \
					build/krepp.o

# rules
#--------------------------------------------
all: $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDLIBS) $(VARDEF) $(INC) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(VARDEF) $(CPPFLAGS) $(LDFLAGS) $(INC) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
