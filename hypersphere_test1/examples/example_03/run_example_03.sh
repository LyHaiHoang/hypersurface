#Edit name of the files
fileoutput=output_equ_vect.xyz
filepos=atomic_positions.xyz
filevect=input_known_vect.xyz
filesym=input_sym_matrix.dat
fileperm=input_perm_matrix.dat

#Create input file
echo "7"                  >input_variables.dat  # number of atoms
echo "3"                 >>input_variables.dat  # number of dimensions
echo "10"                >>input_variables.dat  # number of vectors (not including fixed vectors)
echo "3"                 >>input_variables.dat  # Method for initialization
echo "1"                 >>input_variables.dat  # Method for optimization 
echo "$fileoutput"       >>input_variables.dat  # File where are written the output vectors
echo "$filepos"          >>input_variables.dat  # File for the input atomic positions 
echo "$filevect"         >>input_variables.dat  # File for the input known fix vectors
echo "$filesym"          >>input_variables.dat  # File for the symetry matrixes
echo "$fileperm"         >>input_variables.dat  # File for the permutation matrices

#Fill the atomic positions file
echo "7"                  >$filepos
echo " "                 >>$filepos
echo "2  2.0  0.0  0.0"  >>$filepos
echo "2  1.0  2.0  0.0"  >>$filepos
echo "2 -1.0  2.0  0.0"  >>$filepos
echo "2 -2.0  0.0  0.0"  >>$filepos
echo "2 -1.0 -2.0  0.0"  >>$filepos
echo "2  1.0 -2.0  0.0"  >>$filepos
echo "2  0.0  0.0  0.0"  >>$filepos

#Fill the input known vectors (only one here)
echo "7"                  >$filevect
echo ""                  >>$filevect
echo "0.0 1.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect
echo "0.0 0.0 0.0"       >>$filevect

#Fill the 5 symetry matrixes (no identity) with size = the dimension of the system (here 3*3) 
echo "5"                    >$filesym
echo ""                    >>$filesym
echo " 0.500 -0.866 0.000" >>$filesym
echo " 0.866  0.500 0.000" >>$filesym
echo " 0.000  0.000 1.000" >>$filesym
echo "5"                   >>$filesym
echo ""                    >>$filesym
echo "-0.500 -0.866 0.000" >>$filesym
echo " 0.866 -0.500 0.000" >>$filesym
echo " 0.000  0.000 1.000" >>$filesym
echo "5"                   >>$filesym
echo ""                    >>$filesym
echo "-1.000  0.000 0.000" >>$filesym
echo " 0.000 -1.000 0.000" >>$filesym
echo " 0.000  0.000 1.000" >>$filesym
echo "5"                   >>$filesym
echo ""                    >>$filesym
echo "-0.500  0.866 0.000" >>$filesym
echo "-0.866 -0.500 0.000" >>$filesym
echo " 0.000  0.000 1.000" >>$filesym
echo "5"                   >>$filesym
echo ""                    >>$filesym
echo " 0.500 0.866 0.000"  >>$filesym
echo "-0.866 0.500 0.000"  >>$filesym
echo " 0.000 0.000 1.000"  >>$filesym

#Fill the 6 permutation matrixes with size = the number of atoms (here 7*7)
echo "5"                    >$fileperm
echo ""                    >>$fileperm
echo "0 0 0 0 0 1 0"       >>$fileperm
echo "1 0 0 0 0 0 0"       >>$fileperm
echo "0 1 0 0 0 0 0"       >>$fileperm
echo "0 0 1 0 0 0 0"       >>$fileperm
echo "0 0 0 1 0 0 0"       >>$fileperm
echo "0 0 0 0 1 0 0"       >>$fileperm
echo "0 0 0 0 0 0 1"       >>$fileperm
echo "5"                   >>$fileperm
echo ""                    >>$fileperm
echo "0 0 0 0 1 0 0"       >>$fileperm
echo "0 0 0 0 0 1 0"       >>$fileperm
echo "1 0 0 0 0 0 0"       >>$fileperm
echo "0 1 0 0 0 0 0"       >>$fileperm
echo "0 0 1 0 0 0 0"       >>$fileperm
echo "0 0 0 1 0 0 0"       >>$fileperm
echo "0 0 0 0 0 0 1"       >>$fileperm
echo "5"                   >>$fileperm
echo ""                    >>$fileperm
echo "0 0 0 1 0 0 0"       >>$fileperm
echo "0 0 0 0 1 0 0"       >>$fileperm
echo "0 0 0 0 0 1 0"       >>$fileperm
echo "1 0 0 0 0 0 0"       >>$fileperm
echo "0 1 0 0 0 0 0"       >>$fileperm
echo "0 0 1 0 0 0 0"       >>$fileperm
echo "0 0 0 0 0 0 1"       >>$fileperm
echo "5"                   >>$fileperm       
echo ""                    >>$fileperm
echo "0 0 1 0 0 0 0"       >>$fileperm
echo "0 0 0 1 0 0 0"       >>$fileperm
echo "0 0 0 0 1 0 0"       >>$fileperm
echo "0 0 0 0 0 1 0"       >>$fileperm
echo "1 0 0 0 0 0 0"       >>$fileperm
echo "0 1 0 0 0 0 0"       >>$fileperm
echo "0 0 0 0 0 0 1"       >>$fileperm
echo "5"                   >>$fileperm       
echo ""                    >>$fileperm
echo "0 1 0 0 0 0 0"       >>$fileperm
echo "0 0 1 0 0 0 0"       >>$fileperm
echo "0 0 0 1 0 0 0"       >>$fileperm
echo "0 0 0 0 1 0 0"       >>$fileperm
echo "0 0 0 0 0 1 0"       >>$fileperm
echo "1 0 0 0 0 0 0"       >>$fileperm
echo "0 0 0 0 0 0 1"       >>$fileperm

# execute
./../../src/exe_main.out
