nat=2
ndim=3
nvec=10
init=3 #chose initialization 1= Golden spiral, 2=pure random, 3= Box-Muller
opt=1  #chose optimization 0= nothing, 1=forces, 2= dot_prod, 3= both
filepos=atomic_positions.xyz
fileout=output_equ_vect.xyz
filevect=input_known_vect.xyz

# Fill input variables file
echo "$nat"          >input_variables.dat   # Number of atoms
echo "$ndim"         >>input_variables.dat  # Number of dimensions
echo "$nvec"         >>input_variables.dat  # Number of vectors (not including fixed vectors)
echo "$init"         >>input_variables.dat  # Method for initialization
echo "$opt"          >>input_variables.dat  # optimization technique
echo "$fileout"      >>input_variables.dat  # Outfile where vectors are written 
echo "$filepos"      >>input_variables.dat  # File of atomic positions
echo "$filevect"     >>input_variables.dat  # File containing the known vectors 
echo "dummy"         >>input_variables.dat  # Unused
echo "dummy"         >>input_variables.dat  # Unused

# Fill atomic positions file
echo "2"             >$filepos
echo " "             >>$filepos
echo "1 0.0 0.0 0.0" >>$filepos
echo "2 1.0 0.0 0.0" >>$filepos

# Fill input known vectors (5 vectors here in nat*ndim dimensions)
echo "2"             >$filevect
echo ""              >>$filevect
echo "2.0 1.5 1.0"   >>$filevect
echo "0.0 0.0 0.0"   >>$filevect
echo "2"             >>$filevect
echo ""              >>$filevect
echo "1.0 0.0 0.0"   >>$filevect
echo "0.0 0.0 0.0"   >>$filevect
echo "2"             >>$filevect
echo ""              >>$filevect
echo "0.0 1.0 0.0"   >>$filevect
echo "0.0 0.0 0.0"   >>$filevect
echo "2"             >>$filevect
echo ""              >>$filevect
echo "0.0 0.0 0.0"   >>$filevect
echo "0.0 0.0 1.0"   >>$filevect
echo "2"             >>$filevect
echo ""              >>$filevect
echo "0.0 0.0 0.0"   >>$filevect
echo "1.0 3.0 7.0"   >>$filevect

#Run the code
./../../src/exe_main.out
