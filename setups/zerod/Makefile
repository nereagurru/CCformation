F90=gfortran
H5F=h5fc
FFLAGSCOV=-O2 -fdefault-real-8 -ffree-form -fopenmp -pg --coverage
FFLAGS=-O2 -fdefault-real-8 -ffree-form -fopenmp
#add preprocessors:
CFLAGS=
OPT_FILE ?= ifdef.opt
-include $(OPT_FILE)

SRC_DIR=../../src
OBJ_DIR=../../obj/zerod


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90 | $(OBJ_DIR)
	$(H5F) $(FFLAGS) $(CFLAGS) -J$(OBJ_DIR) -c $< -o $@

zerod: $(OBJ_DIR)/main_PB.o $(OBJ_DIR)/timestep.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/utils.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/ppplts.o
	$(H5F) $(FFLAGS) $(CFLAGS) $(LDFLAGS) $? -o $@

$(OBJ_DIR)/constants.o: $(SRC_DIR)/constants.F90
$(OBJ_DIR)/types.o: $(SRC_DIR)/types.F90
$(OBJ_DIR)/utils.o: $(SRC_DIR)/utils.F90
$(OBJ_DIR)/parameters.o: $(SRC_DIR)/parameters.F90 $(OBJ_DIR)/constants.o
$(OBJ_DIR)/discstruct_itp.o: $(SRC_DIR)/discstruct_itp.F90 $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/utils.o
$(OBJ_DIR)/advection.o: $(SRC_DIR)/advection.F90 $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/types.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/initproblem_itp.o: $(SRC_DIR)/initproblem_itp.F90 $(OBJ_DIR)/advection.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/types.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o
$(OBJ_DIR)/mrgrnk.o: $(SRC_DIR)/mrgrnk.F90 $(OBJ_DIR)/types.o
$(OBJ_DIR)/parallel_sort.o: $(SRC_DIR)/parallel_sort.F90 $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/types.o
$(OBJ_DIR)/grid.o: $(SRC_DIR)/grid.F90 $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o
$(OBJ_DIR)/collisions.o: $(SRC_DIR)/collisions.F90 $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/types.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/output.o: $(SRC_DIR)/output.F90 $(OBJ_DIR)/grid.o $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/timestep.o: $(SRC_DIR)/timestep.F90 $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/hdf5output.o: $(SRC_DIR)/hdf5output.F90 $(OBJ_DIR)/grid.o $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/ppplts.o: $(SRC_DIR)/ppplts.F90 $(OBJ_DIR)/parameters.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/initproblem_itp.o
$(OBJ_DIR)/main_PB.o: $(SRC_DIR)/main_PB.F90 $(OBJ_DIR)/timestep.o $(OBJ_DIR)/output.o $(OBJ_DIR)/hdf5output.o $(OBJ_DIR)/collisions.o $(OBJ_DIR)/advection.o $(OBJ_DIR)/grid.o $(OBJ_DIR)/parallel_sort.o $(OBJ_DIR)/mrgrnk.o $(OBJ_DIR)/initproblem_itp.o $(OBJ_DIR)/discstruct_itp.o $(OBJ_DIR)/parameters.o $(OBJ_DIR)/types.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/ppplts.o

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR)/ zerod *.out *.err *.o *.mod
cleancov:
	rm -rf *.gcno *.gcda *.gcov *.o *.mod *.dat *.out *.err zerod

