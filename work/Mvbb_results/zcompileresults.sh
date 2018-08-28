## this bash script will compile all the Mvbb results into the file $myre
## NOTE: make sure that only the desired versions of the results exist in $IMAMYR
## by: Charlie Payne
## copyright (c): 2016-2018
gV=1.00 # change this accordingly!
gA=1.27 # " " "
precision=12


# function prototype: this function will get the lnum-th line of a file
getline(){
  local lnum=${1}
  local myfile=${2}
  local line=$(sed -n "${lnum}p" $myfile)
  echo $line
}

# function prototype: this function will get the place-th of a space separated string (aka: vector)
# make sure to put double qoutes around "vector", place counting starts from 0 (not 1)
getp(){
  local vector=${1}
  local place=${2}
  for ((p=0; p<$place; p++))
  do
    vector=$(echo $vector | sed 's/[^ ]* //') # remove everything before the current first space
  done
  vector=$(echo $vector | sed 's/\s.*$//') # remove everything after the current first space
  echo $vector
}

# function prototype: this function will undo a number from being in base e or E notation
# make sure that $precision is set somewhere above this prototype
bcE(){
  local value=${1}
  local check1=${value:0:1}
  local check2=${value:0:2}
  if [ $check1 = 'e' ] || [ $check1 = 'E' ]
  then
    value="1${value}"
  elif [ $check2 = '-e' ] || [ $check2 = '-E' ]
  then
    value=${value#-}
    value="-1${value}"
  fi
  local eval=${value#*[eE]}
  local epm=${eval:0:1}
  if [ $epm = '-' ]
  then
    value=$(echo ${value} | sed -e 's/[eE]-*/\/10\^/')
  elif [ $epm = '+' ]
  then
    value=$(echo ${value} | sed -e 's/[eE]+*/\*10\^/')
  elif [ "$eval" = "$value" ]
  then
    : # the value must not be in base e or E notation, derp...
  else
    echo "ERROR: value = $value has problems..."
  fi
  value=$(bc <<< "scale=${precision}; $value")
  echo $value
}

# function prototype: this function will take a number like .012348 and write it as 0.012348, because I ain't no tool!
append0(){
  local value=${1}
  if [ ${value:0:1} = '.' ]
  then
    value=$(echo $value | sed -e 's/./0./')
  elif [ ${value:0:2} = '-.' ]
  then
    value=$(echo $value | sed -e 's/-./-0./')
  fi
  echo $value
}

# function prototype: this function will take an absolute value of a number
absval(){
  local value=${1}
  local prefix=$(echo ${value:0:1})
  if [ $prefix = '-' ]
  then
    value=$(echo ${value#$prefix})
  fi
  echo $value
}


# set up the results file in $myre
basedir=$PWD
myre=zmyresults.txt
rm -f $myre # just in case it already exists
myline='================================================================================================================================================================================================'
echo "nuc.    GTbar/Fbar /Tbar     M0nu_evol_sp_int(NN)_int(3N)_emax_hw      |  M_GT              M_F                M_T              |  M0nu = M_GT - (g_V/g_A)^2 M_F + M_T,  g_V = ${gV}, g_A = ${gA}" >> $basedir/$myre
echo $myline >> $basedir/$myre

# loop through all the relevant directories with M0nu results in them
Zbar=zzzzz
Zre=z.zzzzzzzzzzzz
cd $basedir/M0nu
for nucpath in $basedir/M0nu/????
do
  [ -d $nucpath ] || continue # if not a directory, skip
  nuc=$(basename $nucpath)
  cd $nuc
  for mypath in $nucpath/M0nu_*
  do
    [ -d $mypath ] || continue # if not a directory, skip
    mydir=$(basename $mypath)
    GTbar=$Zbar
    GTre=$Zre
    GTphase='+' # nutbar gives random global phases, define GT as positive
    if [ -d $mypath/GT_* ]
    then
      cd $mypath/GT_*
      GTbar=$(basename $PWD)
      GTbar=$(echo ${GTbar:3:5}) # get the GT barcode
      GTre=$(getline '8' "nutbar_tensor*${GTbar}*.dat")
      GTre=$(getp "$GTre" '9')
      GTre=$(bcE $GTre) # convert the format
      if [ ${GTre:0:1} = '-' ]
      then
        GTre=$(absval $GTre)
        GTphase='-'
      fi
      GTre=$(append0 $GTre)
    fi
    Fbar=$Zbar
    Fre=$Zre
    if [ -d $mypath/F_* ]
    then
      cd $mypath/F_*
      Fbar=$(basename $PWD)
      Fbar=$(echo ${Fbar:2:5}) # get the F barcode
      Fre=$(getline '8' "nutbar_tensor*${Fbar}*.dat")
      Fre=$(getp "$Fre" '9')
      Fre=$(bcE $Fre) # convert the format
      if [ $GTphase = '-' ]
      then
        Fre=$(bc <<< "scale = ${precision}; -1.0*${Fre}")
      fi
      Fre=$(append0 $Fre)
    fi
    Tbar=$Zbar
    Tre=$Zre
    if [ -d $mypath/T_* ]
    then
      cd $mypath/T_*
      Tbar=$(basename $PWD)
      Tbar=$(echo ${Tbar:2:5}) # get the T barcode
      Tre=$(getline '8' "nutbar_tensor*${Tbar}*.dat")
      Tre=$(getp "$Tre" '9')
      Tre=$(bcE $Tre) # convert the format
      if [ $GTphase = '-' ]
      then
        Tre=$(bc <<< "scale = ${precision}; -1.0*${Tre}")
      fi
      Tre=$(append0 $Tre)
    fi
    TOTre=$Zre
    chMH=''
    if [ $GTre != $Zre ] && [ $Fre != $Zre ] && [ $Tre != $Zre ]
    then
      fact=0
      if [[ $mydir = *"_MH"* ]] # Mihai Hiroi (MH) specific parameters
      then
        fact=$(bc <<< "scale = ${precision}; (1.0*1.0)/(1.25*1.25)")
        chMH=', g_A = 1.25'
      else
        fact=$(bc <<< "scale = ${precision}; (${gV}*${gV})/(${gA}*${gA})")
      fi
      TOTre=$(bc <<< "scale = ${precision}; ${GTre} - (${fact}*${Fre}) + ${Tre}")
      TOTre=$(append0 $TOTre)
    fi
    echo "${nuc}    ${GTbar}/${Fbar}/${Tbar}    ${mydir}  |  ${GTre}    ${Fre}    ${Tre}  |  ${TOTre}${chMH}" >> $basedir/$myre
  done
  echo $myline >> $basedir/$myre
  cd $basedir/M0nu # this is necessary for next loop, derp
done

cat $myre

## FIN
