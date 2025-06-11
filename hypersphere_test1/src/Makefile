# Compiler
FC = gfortran
# Flags de compilation
FFLAGS = -O3 -Wall #for debugging -pedantic -fcheck=all -g

# Noms des exécutables et fichiers objets
EXEC = exe_main.out
MAIN = main
MODULE = equipartition_ND
OBJS = $(MAIN).o $(MODULE).o

# Règles de construction
all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(MAIN).o: $(MAIN).f90 $(MODULE).o
	$(FC) $(FFLAGS) -c $<

$(MODULE).o: $(MODULE).f90
	$(FC) $(FFLAGS) -c $<

# Nettoyage des fichiers objets et exécutables
clean:
	rm -f $(EXEC) *.o *.mod

