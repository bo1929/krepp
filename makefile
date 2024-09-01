# compiler options
#--------------------------------------------
COMPILER = g++
LDLIBS = -lstdc++fs -lm -lz -lstdc++ -lcurl
INC = -Iexternal/CLI11/include/CLI -Iexternal/eigen-3.4.0 \
			-Iexternal/l-bfgs-b/include -Iexternal/parallel-hashmap
CXXFLAGS = -std=c++17 -O3 -g -fopenmp
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument

# project files
#--------------------------------------------
PROGRAM = keremet
OBJECTS = build/common.o \
					build/MurmurHash3.o \
					build/lshf.o \
					build/library.o build/query.o build/rqseq.o \
					build/record.o build/phytree.o build/table.o \
					build/keremet.o

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
