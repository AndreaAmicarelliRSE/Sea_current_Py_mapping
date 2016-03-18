# Variables to be updated 
VERSION = 1_0
#    bin, debug
EXECUTION = bin
#    -O1 -I/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/include -L/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/lib -lnetcdf -lnetcdff  (ifort bin)
#    -O1 -I/cineca/prod/libraries/netcdf/4.3.2/gnu--4.9.2/include -L/cineca/prod/libraries/netcdf/4.3.2/gnu--4.9.2/lib -lnetcdf -lnetcdff  (gfortran bin)
#    -g -O0 -fbacktrace -C -I/cineca/prod/libraries/netcdf/4.3.2/gnu--4.9.2/include -L/cineca/prod/libraries/netcdf/4.3.2/gnu--4.9.2/lib -lnetcdf -lnetcdff  (gfortran debug)
#    -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all -I/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/include -L/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/lib -lnetcdf -lnetcdff
#    (ifort debug)
COMPILATION_FLAGS = -O1 -I/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/include -L/cineca/prod/libraries/netcdf/4.1.3/intel--cs-xe-2015--binary/lib -lnetcdf -lnetcdff
#    gfortran, ifort
COMPILER = ifort
#    -fopenmp, -openmp, (blank for debug)  
OMP_FLAG = -openmp
EXE_DIR = ../$(EXECUTION)/
SOURCES = \
$(wildcard ./*/*.f90)
# Other variables
MAIN_FILE_ROOT = ./Main/Sea_current_Py_mapping
CODE = $(EXE_DIR)Sea_current_Py_mapping_v_$(VERSION)_$(COMPILER)_$(EXECUTION)
OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))
# Targets
#    Correct sequence for compile+link:
#       make touch
#       make 
all: compile_main 
#    "compiling" cannot be launched together with "touching" as ".o" files
#    would not be available yet
compile_main:
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o
compile: $(SOURCES)
%.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS) -c 
link: $(OBJECTS)
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o
remove: $(OBJECTS)
	rm -f $^
	rm -f *.mod
touch:
	touch $(OBJECTS)
clean: 
	rm -f $(CODE).x 

