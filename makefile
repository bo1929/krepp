# compiler options
#--------------------------------------------
COMPILER = g++
LDLIBS = -lm -lz -lstdc++ -lcurl
CXXFLAGS = -g -std=c++2a -O3 -fopenmp
WFLAGS = -Wno-unused-result -Wno-unused-command-line-argument

# project files
#--------------------------------------------
PROGRAM = keremet
OBJECTS = build/MurmurHash3.o build/common.o \
					build/lshf.o build/refseq.o build/subset.o build/tree.o build/table.o \
					build/keremet.o

# rules
#--------------------------------------------
all: $(PROGRAM)

# generic rule for compiling *.cpp -> *.o
build/%.o: src/%.cpp
	@mkdir -p build
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -c src/$*.cpp -o build/$*.o

$(PROGRAM): $(OBJECTS)
	$(COMPILER) $(WFLAGS) $(CXXFLAGS) $+ $(LDLIBS) $(CPPFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f $(PROGRAM) $(OBJECTS)
	@if [ -d build ]; then rmdir build; fi
	@echo "Succesfully cleaned."
