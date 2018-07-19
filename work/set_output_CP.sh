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


CPclean $dirc1

CPsafe $dirs1
CPsafe $dirs2
CPsafe $dirs3
CPsafe $dirs4
CPsafe $dirs5
CPsafe $dirs6
