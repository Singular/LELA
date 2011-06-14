CXX = g++
CXXFLAGS = -g -O2
CPPFLAGS = -DDETAILED_PROFILE -DDEBUG
INCLUDES = -I..
LIBS = ../linbox/.libs/liblinbox.a -lgmp -lgivaro -lpng -lm4ri -Wl,-rpath,/usr/local/lib

TARGETS = test-solver test-matrix-domain-gf2 

all: $(TARGETS)

test-solver: test-solver.o
test-matrix-domain-gf2: test-matrix-domain-gf2.o

test-solver.o: test-solver.C solver.h
test-matrix-domain-gf2.o: test-matrix-domain-gf2.C

%.o : %.C
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

% : %.o
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f $(TARGETS) *.o *.out *~
