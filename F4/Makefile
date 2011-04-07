CXXFLAGS = -g -O2
CPPFLAGS = -DDETAILED_PROFILE -DDEBUG
INCLUDES = -I..
LIBS = ../linbox/.libs/liblinbox.a -lgmp -lgivaro -lpng -lm4ri -Wl,-rpath,/usr/local/lib

TARGETS = test-solver test-gauss-jordan test-matrix-domain-gf2 

all: $(TARGETS)

test-solver: test-solver.o
test-gauss-jordan: test-gauss-jordan.o
test-matrix-domain-gf2: test-matrix-domain-gf2.o

test-solver.o: test-solver.C gauss-jordan.h solver.h from-png.h from-png.inl
test-gauss-jordan.o: test-gauss-jordan.C gauss-jordan.h
test-matrix-domain-gf2.o: test-matrix-domain-gf2.C

%.o : %.C
	g++ -c $(CXXFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

% : %.o
	g++ $(CXXFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f $(TARGETS) *.o *.out *~
