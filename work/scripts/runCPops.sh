#!/bin/bash

exe=../compiled/CPops
#exe="valgrind --track-origins=yes ../compiled/CPops"

emax=4
dirout=$IMAWRK/debug_output/
intfile="${dirout}" # this is hacky... but it works :/

A=76
hw=10.00
#A=48
#hw=16.00


## DO: M0nu_adpt

#Decay=GT
#Decay=F
Decay=T
Reduced=NR
#Reduced=R
Ec=7.72
SRC=none
#SRC=AV18
#SRC=CD-Bonn
#SRC=Miller-Spencer
#SRC=debug
Operators=M0nu_adpt_${dirout}_${Decay}_${Reduced}_${Ec}_${SRC}

## OR: M0nu_Print_Integrand

##Mode=0 # print first, integrate later
#Mode=1 # integrate first, print later
#Ec=7.72
#SRC=none
#nn=10
#ll=9
#np=10
#lp=9
#qa=0
#qb=1000
#Operators=M0nu_PrintIntegrand_${dirout}_${Mode}_${Ec}_${SRC}_${nn}_${ll}_${np}_${lp}_${qa}_${qb}
#whamit=1 # set this to one to automatically do wham.sh below


# the NumLine's below are for the M0nu TBME! (OLD)
#emax=1
#NumLine=20
#emax=2
#NumLine=272
#emax=3
#NumLine=2102
# the NumLine's below are for the scalar TBME! (OLD)
#emax=1
#NumLine=108
#emax=2
#NumLine=1596
#emax=3
#NumLine=12560
#emax=4
#NumLine=67156
#emax=5
#NumLine=274836
#emax=6
#NumLine=927680
#emax=7
#NumLine=2705004
#emax=8
#NumLine=7033580
#emax=9
#NumLine=16675208
#emax=10
#NumLine=36641676


valence_space="fp-shell" # has to be a valid space, even if not used, seems weird - should probs tell Ragnar...
reference="Ca$A" # same deal...

$exe emax=${emax} valence_space=${valence_space} hw=${hw} reference=${reference} Operators=${Operators} A=${A} intfile=${intfile}

if [ -z $whamit ]
then
  whamit=0
fi
if [ $whamit -eq 1 ]
then
  cd $dirout
  ./wham.sh ${nn} ${ll} ${np} ${lp}
  cd $IMARUN
fi
