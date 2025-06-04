fileperm=input_perm_matrix.dat
filesym=input_sym_matrix.dat
filepos=atomic_positions.xyz
filevect=input_known_vect.xyz

echo "1" >input_variables.dat        # number of atoms
echo "3" >>input_variables.dat       # number of dimensions
echo "10" >>input_variables.dat       # number of vectors (not including fixed vectors)
echo "2"        >>input_variables.dat  # Method for initialization
echo "output_equ_vect.xyz" >>input_variables.dat
echo "$filepos" >>input_variables.dat
echo "$filevect" >>input_variables.dat
echo "$filesym" >>input_variables.dat
echo "$fileperm" >>input_variables.dat

echo "1" >$filepos
echo " ">>$filepos
echo "2 0.0 0.0 0.0">>$filepos

echo "1" >$filesym
echo "">>$filesym
echo "1.000 0.000 0.000">>$filesym
echo "0.000 1.000 0.000">>$filesym
echo "0.000 0.000 -1.000">>$filesym

echo "1" >$fileperm
echo "">>$fileperm
echo "1">>$fileperm

echo "1" >$filevect
echo "">>$filevect
echo "0.0 0.0 1.0">>$filevect

./../../src/exe_main.out
