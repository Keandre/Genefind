if [ -z "$1" ]
then

  echo "The script requires an argument. Example usage is:"
  echo "    simulate.sh gro-filename-without-extension"
  
else

  gmx mdrun -v -deffnm $1
  
fi