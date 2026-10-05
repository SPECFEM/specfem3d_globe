# Run from the OpenMP build directory that 9i.test_gf_openmp.sh configures,
# with TESTDIR set to this directory: the flags, the module directory and
# the library are that build's, the sources this directory's.
include Makefile

default: test_gf_cache_omp

L := ./lib

# test_gf_cache.makefile's link, against the OpenMP build of the library.
# ${FCCOMPILE_CHECK} carries that build's FCFLAGS, OpenMP flag included.
test_gf_cache_omp:
	mkdir -p ./bin
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_cache \
		$(TESTDIR)/gf_manufactured.f90 $(TESTDIR)/test_gf_cache.f90 $L/libgf3d.a $(LDFLAGS)
