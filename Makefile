
include make.inc

all: lib testlib

lib: 
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) all )

testlib:
	( cd ./test; $(MAKE) run )

install:
	cp ./lib/libgaussian_spherical_harmonic.a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../gaussian_spherical_harmonic $(BIN_PATH)

clean: 
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )

.PHONY: all lib testlib install
