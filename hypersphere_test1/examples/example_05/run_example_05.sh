dim=3
filepos=atomic_positions.xyz
filevec=input_known_vect.xyz
fileout=output_equ_vect.xyz
# Create the input variables file
echo "1"        > input_variables.dat    # number of atoms
echo "$dim"     >>input_variables.dat    # number of dimensions
echo "2"        >>input_variables.dat    # number of vectors (not including fixed vectors)
echo "3"        >>input_variables.dat    # method for initialization
echo "$fileout" >>input_variables.dat
echo "$filepos" >>input_variables.dat
echo "$filevec" >>input_variables.dat
echo "dummy"    >>input_variables.dat
echo "dummy"    >>input_variables.dat

# Create atomic positions file
echo "1" >$filepos
echo " ">>$filepos
echo "2 0.0 0.0 0.0">>$filepos
# Modify the number of coordinates in the atomic positions file
pos="2 ";
for i in $(seq 1 $dim ); do
    pos+=" 0.0"
done   
sed -i '3d' $filepos;
echo "$pos" >>$filepos

# Create fix vectors file
#echo "1" >$filevec
#echo "">>$filevec
#echo "1.0 0.0 0.0">>$filevec

for nv in {2..130}; do
   echo "$nv";
   # Modify the number of vectors in the parameter file 
   sed -i "3i $nv"  input_variables.dat
   sed -i '4d'      input_variables.dat 
   # Do the equipartition and print
   dist=$(./../../src/exe_main.out)
   echo "$nv $dist" >>results_$dim.dat
done
