dirc1=debug_output
dirs1=output
dirs2=output_to_javier

CPclean(){
  local mydir=${1}
  rm -rf $mydir
  mkdir $mydir
  cp $IMASMS/zscripts/z*.sh $mydir
}

CPsafe(){
  local mydir=${1}
  mkdir -p $mydir
  cp $IMASMS/zscripts/z*.sh $mydir
}

CPclean $dirc1
CPsafe $dirs1
CPsafe $dirs2
