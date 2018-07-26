dirc1=debug_output

dirs1=output
dirs2=output_Ca48
dirs3=output_Ge76
dirs4=output_Se82
dirs5=output_to_javier
dirs6=output_to_mihai


CPclean(){
  local mydir=${1}
  rm -rf $mydir
  mkdir $mydir
  cd $mydir
  cp $IMASMS/zscripts/z*.sh .
  ./znewrecord.sh
  cd ..
}

CPsafe(){
  local mydir=${1}
  mkdir -p $mydir
  cd $mydir
  cp $IMASMS/zscripts/z*.sh .
  ./znewrecord.sh
  cd ..
}

inparam=${1}
inans='off'
if [ -z $inparam ]
then
  inparam='off'
  echo 'this script will do a znewrecord.sh (not recommended if runs are currently in progress) in each directory...'
  echo 'are you sure you want to do this? (Y/N)'
  read inans
fi
if [ $inans = 'y' ] || [ $inans = 'Y' ] || [ $inparam = 'true' ]
then
  CPclean $dirc1

  CPsafe $dirs1
  CPsafe $dirs2
  CPsafe $dirs3
  CPsafe $dirs4
  CPsafe $dirs5
  CPsafe $dirs6
else
  echo 'safely exiting...'
fi
