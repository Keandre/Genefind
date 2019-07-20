echo ''
echo 'This script will remove the following generated files:'

filesToRemove=()

for entry in $PWD/*
do
  filename=${entry##*/}
  fileext=${entry##*.}
  if [ $fileext != "sh" ] && [ $fileext != "gro" ] && [ $fileext != "mdp" ]
  then
    echo "    $filename"
    filesToRemove+=( $entry )
  else
    case $entry in
	*_processed.gro)
      echo "    $filename"
      filesToRemove+=( $entry )
	  ;;
	*mdout.mdp)
      echo "    $filename"
      filesToRemove+=( $entry )
	  ;;
    esac
  fi
done

echo ''
echo 'press ENTER to continue.'
read -p ''

for entry in ${filesToRemove[@]};
do
  filename=${entry##*/}
  echo "removing $filename"
  rm $entry
done

echo ''
echo 'Finished cleaning up.'