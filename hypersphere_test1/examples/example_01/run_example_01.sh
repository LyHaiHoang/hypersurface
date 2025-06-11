nat=1
ndim=7

for init in 1 2 3; do
  filepos="atomic_positions_${init}.xyz"
  fileout="output_equ_vect_${init}.xyz"
  
  #Fill the file of variables
  echo "$nat"      >input_variables.dat  # number of atoms
  echo "$ndim"    >>input_variables.dat  # number of dimensions
  echo "5000"     >>input_variables.dat  # number of vectors
  echo "$init"    >>input_variables.dat  # chose initialization 1= Golden spiral, 2=pure random, 3= Box-Muller
  echo "0"        >>input_variables.dat  # chose optimization 0= nothing, 1=forces, 2= dot_prod, 3= both
  echo "$fileout" >>input_variables.dat  # output file
  echo "$filepos" >>input_variables.dat  # File with atomic positions
  echo "dummy"    >>input_variables.dat  # Not needed there
  echo "dummy"    >>input_variables.dat  # Not needed there
  echo "dummy"    >>input_variables.dat  # Not needed there
  
  #Fill the file of atomic positions
  echo "$nat" >$filepos
  echo "">>$filepos
  for j in $(seq 1 $nat); do
    atom="${j}"
    for i in $(seq 1 $ndim); do atom="$atom "0.0; done
    echo $atom>>$filepos
  done
  
  #Run the code
  ./../../src/exe_main.out
  mv angles.dat "angles_${init}.dat"
done
