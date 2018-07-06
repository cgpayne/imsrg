#!/bin/bash

exe=../compiled/CPops
#exe="valgrind --track-origins=yes ../compiled/CPops"

emax=4

#Operators=Engel
#intfile=$IMAJEA/new_mod_SRC_Argonne_F
#Operators=Mihai
#intfile=$IMAWRK/mihai_0vbb_data/vbbGT_N10_full_48Ca_E2

A=76
hw=10.00

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
#dirout=$IMAOUT/
dirout=$IMAWRK/debug_output/
intfile=${dirout}
Operators=M0nu_adpt_${dirout}_${Decay}_${Reduced}_${Ec}_${SRC}

#dirout=$IMAOUT/plotting_M0nu/zz_integrand/
##Mode=0 # print first, integrate later
#Mode=1 # integrate first, print later
#nn=10
#ll=9
#np=10
#lp=9
#qa=0
#qb=1000
#Operators=M0nu_PrintIntegrand_${dirout}_${Mode}_${Ec}_${SRC}_${nn}_${ll}_${np}_${lp}_${qa}_${qb}
#whamit=1

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

vnn=none
v3n=none
valence_space=sd-shell
#valence_space=Ca40
reference=Ca$A
scratch=SCRATCH

$exe 2bme=${vnn} 3bme=${v3n} emax=${emax} e3max=${e3max} valence_space=${valence_space} hw=${hw} smax=${smax} ${file3} omega_norm_max=${omega_norm_max} reference=${reference} Operators=${Operators} scratch=${scratch} A=${A} use_brueckner_bch=false intfile=${intfile}

if [ -z $whamit ]
then
  whamit=0
fi
if [ $whamit -eq 1 ]
then
  cd $IMAOUT/plotting_M0nu/zz_integrand
  ./wham.sh ${nn} ${ll} ${np} ${lp}
  cd $IMARUN
fi
