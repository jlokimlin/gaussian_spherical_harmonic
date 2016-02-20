
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

module_fast_fourier_transform.o: \
module_fast_fourier_transform.f90
	$(FC) $(FFLAGS) -c module_fast_fourier_transform.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GaussianSphericalHarmonic.o: \
module_fast_fourier_transform.o \
type_GaussianSphericalHarmonic.f90
	$(FC) $(FFLAGS) -c type_GaussianSphericalHarmonic.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
type_GaussianSphericalHarmonic.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	ar rc libgaussian_spherical_harmonic.a $(OBJ_DIRECTORY)/*.o
	mv libgaussian_spherical_harmonic.a $(LIB_DIRECTORY)

all:
	make clean
	make exec; $(EXECUTABLE)
	make lib

clean:
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)