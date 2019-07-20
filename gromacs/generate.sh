if [ -z "$1" ]
then

  echo "The script requires an argument. Example usage is:"
  echo "    generate.sh gro-filename-without-extension"
  
else

  read -p "press ENTER to process topology (step 1/3)"
  gmx pdb2gmx -f $1.gro -o $1_processed.gro -water spce
  
  read -p "press ENTER to write binary (step 2/3)"
  gmx grompp -f sample.mdp -c $1_processed.gro -p topol.top -o $1.tpr
  
  read -p "press ENTER to start simulation (step 3/3)"
  gmx mdrun -v -deffnm $1

fi
