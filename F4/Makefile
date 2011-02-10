CXXFLAGS = -g -O2
CPPFLAGS = -DDETAILED_PROFILE
INCLUDES = -I..
LIBS = ../linbox/.libs/liblinbox.a -lgmp -lgivaro -lpng -lm4ri -Wl,-rpath,/usr/local/lib

TARGETS = test-solver test-gauss-jordan test-hybrid-vector test-bit-subvector test-matrix-domain-gf2

all: $(TARGETS)

test-solver: test-solver.o
test-gauss-jordan: test-gauss-jordan.o
test-hybrid-vector: test-hybrid-vector.o
test-bit-subvector: test-bit-subvector.o
test-matrix-domain-gf2: test-matrix-domain-gf2.o

test-solver.o: test-solver.C gauss-jordan.h solver.h from-png.h from-png.inl
test-gauss-jordan.o: test-gauss-jordan.C gauss-jordan.h
test-hybrid-vector.o: test-hybrid-vector.C
test-bit-subvector.o: test-bit-subvector.C
test-matrix-domain-gf2.o: test-matrix-domain-gf2.C

%.o : %.C
	g++ -c $(CXXFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

% : %.o
	g++ $(CXXFLAGS) $^ $(LIBS) -o $@

clean:
	rm -f $(TARGETS) *.o
