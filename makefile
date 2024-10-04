# compiler options
#--------------------------------------------
COMPILER = g++
LDLIBS = -lstdc++fs -lm -lz -lstdc++ -lcurl
INC = -Iexternal/CLI11/include/CLI \
			-Iexternal/parallel-hashmap -Iexternal/boost/libs/math/include
CXXFLAGS = -std=c++17 -O3 -g -fopenmp -funroll-loops
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument

# project files
#--------------------------------------------
PROGRAM = krepp
OBJECTS = build/common.o \
					build/MurmurHash3.o \
					build/lshf.o \
					build/library.o build/query.o build/rqseq.o \
					build/record.o build/phytree.o build/table.o \
					build/krepp.o

# rules
#--------------------------------------------
all: $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $(INC) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(CPPFLAGS) $(LDFLAGS) $(INC) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
