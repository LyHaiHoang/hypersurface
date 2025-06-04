filepos=atomic_positions.xyz
fileout=output_equ_vect.xyz
# Create the input parameters files
echo "1"        >input_variables.dat           # number of atoms
echo "8"        >>input_variables.dat          # number of dimensions
echo "120"      >>input_variables.dat          # number of vectors (not including fixed vectors)
echo "3"        >>input_variables.dat          # method for initialization
echo "$fileout" >>input_variables.dat
echo "$filepos" >>input_variables.dat
echo "dummy"    >>input_variables.dat
echo "dummy"    >>input_variables.dat
echo "dummy"    >>input_variables.dat

# Create atomic positions file
echo "1" >$filepos
echo " ">>$filepos
echo "2 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0">>$filepos
#./../../src/exe_main.out

for nvec in {2..50}; do
   echo "$nvec";
   # Modify the dimension in the parameter file 
   sed -i "3i $nvec" input_variables.dat
   sed -i '4d'      input_variables.dat 
   # Modify the number of coordinates in the atomic positions file
   #pos="2 ";
   #for i in $(seq 1 $dim ); do
   #    pos+=" 0.0"
   #done   
   #sed -i '3d' $filepos;
   #echo "$pos" >>$filepos
   # Do the equipartition and print
   dist=$(./../../src/exe_main.out)
   echo "$nvec $dist" >>results.dat
done
