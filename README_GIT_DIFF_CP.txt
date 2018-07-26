LEGEND:
-> = main file(s) involving the described change
=> = secondary file(s) which were (slightly) altered in order to properly incorporate the described change

NOTE:
changes are listed in order from smallest to largest (measured wrt the amount of lines added/altered, or the "importance" of the alteration)

MODULES:
--for my .bashrc on cougar I load the modules (07/05/18):
module purge
module load gcc/5.2.0
module load gsl/1.16
module load boost/1.58.0
module load openblas/0.2.14
module load anaconda/2.4.0
module load hdf5/1.8.15
module load openmpi/gnu/1.10.0
module load git/2.17.0
module load openssl/1.0.2
--for my .bashrc on oak I load the modules (07/05/18):
module purge
module load StdEnv
module load gcc-7.3.0-gcc-4.8.5
module load openmpi/gcc_2.1.1
module load intel/17.0.4
module load intel-mkl-2018.2.199-gcc-4.8.5
module load gsl/2.4
module load anaconda2/2.4.0
module load boost/1.58.0

BONUS:
To go from the IMSRG evolved TBMEs of 0vBB to the NMEs (using NuShellX), see the fancy scripts I made in the repository:
https://github.com/cgpayne/SMscripts.git


My editor automatically trims trailing whitespace, so I recommend using 'git diff -w' in place of 'git diff' when checking the changes I made.
    -> many files...


Moved the git ignore file from src to the root repository and modified it accordingly.
    -> .gitignore


Symlinked Parameters.hh in work/compiled to src, pretty sure it's all outdated anyways...
    -> work/compiled/Parameters.hh


Added a very simple script to set up my preferred imsrg++ output directories, read it...
    -> work/set_output_CP.sh


Added this README (I love Goedelian loops) to described the changes I made from Ragnar's git repository (this is somewhat lacking/incomplete atm).
    -> README_GIT_DIFF_CP.txt


Added a directory where I store/compile duplicates of my 0vBB NME results.
    -> work/Mvbb_results/zcompileresults.sh


Added a version of Universal.cc called "CPops.cc", primarily to test/debug my (CP) operators (ops).
To compile it, first compile imsrg++ by going to src and doing:
make clean
make
make install
and then go back to work/scripts and enter:
make CPops
    -> work/compiled/CPops.cc
        => work/scripts/runCPops.sh


Added a bunch of output from CPops, which is very useful for calibrating different versions of the NDBD code.
In particular, if you compile imsrg++ with NDBD.hh set to the "JE" parameters (also change the M_NUCLEON in Constants.hh), running ./runCPops.sh will get the "CALIBRATOR" with the parameters:
emax=4
A=76
hw=10.00
Decay=T
Reduced=NR
Ec=7.72
SRC=none
dirout=<wherever_you_like>
If you can't reproduce the CALIBRATOR to within machine precision, something is wrong!
I spent countless hours benchmarking these TBMEs against those provided to me by Jon Engel (JE), who's an expert in the field (and a very nice dude).
Hence, I recommend being able to match the CALIBRATOR before submitting any large jobs. To be extra safe, also see if you can reproduce the TBMEs for the 'GT' and 'F' decays.
    -> work/output_runCPops/*
        => work/scripts/runCPops.sh


Added in a line to the Makefile in order to compile with MKL (in replace of OpenBLAS) on the oak cluster (see MODULES above).
For some reason, on oak there's some conflict between gcc/OpenMP/OpenBLAS which causes an infinite printing of an "OpenBLAS WARNING" which quickly consumes all the memory of a node.
To fix this, Ragnar invented this -DOPENBLAS_NOUSEOPENMP=1 flag, which turns off OpenMP (hence parallelization) for the IMSRG evolution step (which greatly depreciates runtime/speed).
I recommend using this OpenBLAS->MKL change because it allows full parallelization of the code (read: less runtime), although Ragnar doesn't prefer this because it's not open source (fair).
    -> src/Makefile


Added my own versions of Ragnar's python script in order to run imsrg++ from either the cougar cluster or the oak cluster.
As of July 26 2018, goOak.py is more up to date than goCougar.py.
It won't necessarily run on your user/machine, but with some tweaking, it should work!
    -> work/scripts/goCougar.sh
    -> work/scripts/goOak.sh
        => src/Makefile


I decided to keep the fundamental constants and alike (which were shared within imsrg++) in one place; i.e., making them consistent throughout the code (which you can see, they often aren't).
    -> src/Constants.hh
        => HartreeFock.*
        => IMSRG.hh
        => Makefile
        => Modelspace.hh
        => NDBD.*
        => Operator.*
        => ReadWrite.*
        => TwoBodyME.*
        => imsrg_util.*


Added in shell gap functionality, which actually doesn't make sense when comparing a VS approach to CC I think - so it's pretty much useless (for now).
    -> imsrg++.cc
        => src/Parameters.hh


Made "OperatorFromString" compatible with editing the opname string itself (was important for the NDBD operator name/string (see below), in particular).
Pretty confident Ragnar wasn't happy with this change - it is pretty hacky, and probably not advisable, but it was necessary.
    -> imsrg_util.cc
    -> imsrg++.cc


Neurtinoless Double-Beta Decay (NDBD)...
In order to calculate the TBMEs of this operator, it was necessary (accept this, trust me, #PandorasBox) to create a new class to interface with imsrg_util.* and alike.
I tried to note this code out as much as a I could, but if there's any confusing, maybe it's best just to email me (if not treat it as a blackbox): c.p0yne@gmail.com
    -> src/NDBD.*
        => imsrg_util.*
        => src/imsrg++.cc
        => src/ReadWrite.hh/cc

