
include ../make.inc

MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-L$(LIB_DIRECTORY) -lgaussian_spherical_harmonic \
-I$(MODULE_PATH)

EXECUTABLE = ./test.exe

test_spharmt.o: \
test_spharmt.f90
	$(FC) $(FFLAGS) -c test_spharmt.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

run:
	make clean
	make exec; $(EXECUTABLE)
	 
clean:
	rm -f $(OBJS)
	rm -f $(EXECUTABLE)
