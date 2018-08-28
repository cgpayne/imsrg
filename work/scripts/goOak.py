#!/usr/bin/env python

##########################################################################
##              adapted from Ragnar Stroberg's goCougar.py
##              things I should consistently change are:
##                  * wrkdir, MNU, smax, omega_norm_max
##                  file paramater list syle, ARGS['e3max']
##                  * (A,Z) in...
##                  valence_space, custom_valence_space
##                  * hw in [..., e in [..., MNU['Decay'] in [...
##                  2bme, 3bme, LECs
##  						-Charlie Payne
######################################################################

from os import path,environ,mkdir,remove,chdir,getcwd
from sys import argv
from subprocess import call,PIPE
from time import time,sleep
from datetime import datetime

### Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'SLURM'

### The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS=32
exe = '%s/bin/imsrg++'%(environ['HOME'])

### Flag to swith between submitting to the scheduler or running in the current shell
#batch_mode=False
batch_mode=True
if 'terminal' in argv[1:]: batch_mode=False

### Don't forget to change this. I don't want emails about your calculations...
mail_address = 'cgpayne@triumf.ca'

### This comes in handy if you want to loop over Z
ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

### MNU is a (string => string) dictionary of input variables that are specifically for the M0nu/M2nu operator (CP)
wrkdir = '/global/home/cgpayne/oakwork/imsrg_code/imsrg/work/'
#wrkdir = '/global/home/cgpayne/'
chdir(wrkdir)
print 'wrkdir    = ',wrkdir
MNU  =  {}
#MNU['outname'] = 'debug_output/'
MNU['outname'] = 'output/'
#MNU['outname'] = 'output_Ca48/'
#MNU['outname'] = 'output_Ge76/'
#MNU['outname'] = 'output_Se82/'
#MNU['outname'] = 'output_to_javier/'
#MNU['outname'] = 'output_to_mihai/'
MNU['dirout'] = wrkdir + MNU['outname']
MNU['Reduced'] = 'R'
MNU['Ec'] = '7.72' # Ca48
#MNU['Ec'] = '9.41' # Ge76,Se82
MNU['SRC'] = 'none'
#MNU['SRC'] = 'AV18'
#MNU['SRC'] = 'CD-Bonn'
#MNU['SRC'] = 'Miller-Spencer'
#MNU['int'] = 'BARE' # hw = 45*A^{-1/3} - 25*A^{-2/3}
MNU['int'] = 'magic' # hw = 16,24
#MNU['int'] = 'mag16' # hw = 16
#MNU['int'] = 'v3trans' # hw = 14
#MNU['BB'] = 'NN'
MNU['BB'] = '3N'
#MNU['BB'] = 'HF'
#MNU['BB'] = 'OS' # only for 'BARE' int
#MNU['gap'] = 'on' # this functionality might not actually be useful... talk to Ragnar...
MNU['gap'] = 'off'
print 'outname   = ',MNU['outname']
print 'dirout    = ',MNU['dirout']
print 'Reduced   = ',MNU['Reduced']
print 'Ec        = ',MNU['Ec']
print 'SRC       = ',MNU['SRC']
print 'BB        = ',MNU['BB']
print 'int       = ',MNU['int']
print 'gap?      = ',MNU['gap']

### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}

### Shell gap to match with Gaute (still in testing phase)
if MNU['gap'] == 'on':
  #ARGS['gapE'] = '-40'
  ARGS['gapE'] = '-30'
  #ARGS['gapE'] = '-20'
  #ARGS['gapE'] = '-10'
  #ARGS['gapE'] = '-5'
  #ARGS['gapE'] = '-2'
  ARGS['gapShell'] = '13'
  print 'running with an artificial gap of %s MeV for shells 0 through %s...'%(ARGS['gapE'],ARGS['gapShell'])
else:
  print 'running with NO artificial gap...'

### Maximum value of s, and maximum step size ds
if MNU['int'] == 'BARE' or MNU['BB'] == 'HF':
  ARGS['smax'] = '0'
else:
  #ARGS['smax'] = '10'
  ARGS['smax'] = '500'
print 'smax      = ',ARGS['smax']
ARGS['dsmax'] = '0.5'

### Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.25'
#ARGS['omega_norm_max'] = '500000' # for Jiangming
print 'o.norm.m. = ',ARGS['omega_norm_max']

### Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
#ARGS['scratch'] = 'SCRATCH'

### Generator for core decoupling, can be atan, white, imaginary-time.  (atan is default)
#ARGS['core_generator'] = 'white'
### Generator for valence deoupling, can be shell-model, shell-model-atan, shell-model-npnh, shell-model-imaginary-time (shell-model-atan is default)
#ARGS['valence_generator'] = 'shell-model-imaginary-time'
#ARGS['valence_generator'] = 'shell-model-atan-dNgt0'

### Solution method
ARGS['method'] = 'magnus'
#ARGS['method'] = 'NSmagnus'
#ARGS['method'] = 'brueckner'
#ARGS['method'] = 'flow'
#ARGS['method'] = 'HF' # not too sure what this is about...
#ARGS['method'] = 'MP3'
#ARGS['method'] = 'FCI'
print 'method    = ',ARGS['method']

### Tolerance for ODE solver if using flow solution method
#ARGS['ode_tolerance'] = '1e-5'

### Make the cluster run files (PBS or SLURM)
if BATCHSYS == 'PBS':
  FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q oak
#PBS -d %s
#PBS -l walltime=512:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=251gb
#PBS -m ae
#PBS -M %s
#PBS -j oe
#PBS -o imsrg_log/%s.o
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
%s
"""

elif BATCHSYS == 'SLURM':
  FILECONTENT = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --output=imsrg_log/%s.%%j
#SBATCH --time=%s
#SBATCH --mail-user=%s
#SBATCH --mail-type=END
cd $SLURM_SUBMIT_DIR
echo NTHREADS = %d
export OMP_NUM_THREADS=%d
time srun %s
"""

### Model space parameters used for reading Darmstadt-style or Johannes-style interaction files (this set-up is janky, and Ragnar knows it)
if MNU['int'] == 'BARE':
  print 'running BARE-style...'
  ARGS['file2e1max'] = '0 file2e2max=0 file2lmax=0' # I don't think...
  ARGS['file3e1max'] = '0 file3e2max=0 file3e3max=0' # ...this is necessary...
  ARGS['e3max'] = '0' # ...but, oh well
  if MNU['BB'] != 'OS':
    print 'for MNU[int] = BARE, choose MNU[BB] = OS'
    print 'exiting...'
    exit()
elif MNU['int'] == 'magic':
  print 'running Johannes-style int...'
  ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
  ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=14'
  ARGS['e3max'] = '14'
elif MNU['int'] == 'mag16':
  print 'running Johannes-style int, higher e3max...'
  ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
  ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
  ARGS['e3max'] = '16'
elif MNU['int'] == 'v3trans':
  print 'running GTquenching-style int...'
  ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=14'
  ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=16'
  ARGS['e3max'] = '16'
else:
  print 'this *int* has not been set up yet!'
  print 'exiting...'
  exit()
print 'e3max     = ',ARGS['e3max']

### Make a directory for the log files, if it doesn't already exist
if not path.exists('imsrg_log'): mkdir('imsrg_log')

### Loop over multiple jobs to submit
ARGS['BetaCM'] = '0'
ARGS['hwBetaCM'] = '16'
#for (A,Z) in [(40,20)]: # Ca40
for (A,Z) in [(48,20)]: # Ca48
#for (A,Z) in [(42,20),(44,20),(46,20),(48,20),(50,20),(52,20),(54,20),(56,20),(58,20),(60,20)]: # Ca42-60 CHAIN
#for (A,Z) in [(44,22),(46,22),(48,22),(50,22),(52,22),(54,22),(56,22),(58,22)]: # Ti42-60 CHAIN
#for (A,Z) in [(46,24),(48,24),(50,24),(52,24),(54,24),(56,24),(58,24),(60,24)]: # Cr42-60 CHAIN
#for (A,Z) in [(76,32)]: # Ge76
#for (A,Z) in [(82,34)]: # Se82
  print 'Z ========= ',Z
  print 'A ========= ',A
  ARGS['A'] = '%d'%A
  ### Set the reference
  reference = '%s%d'%(ELEM[Z],A)
  ARGS['reference'] = reference
  print 'ref       = ',ARGS['reference']
  ### Set the space (sp)
  print 'running for %s%d...'%(ELEM[Z],A)
  if A == 40 and Z == 20:
    ARGS['valence_space'] = 'fp-shell' # AKA: fp
    #ARGS['valence_space'] = 'Ca40' # for Jiangming
  elif A >= 42 and (Z == 20 or Z == 22 or Z == 24):
    ARGS['valence_space'] = 'fp-shell' # AKA: fp
    #ARGS['valence_space'] = 'Ca48' # for Jiangming
  elif A == 76 and Z == 32:
    ARGS['valence_space'] = 'pf5g9' # this is just a label when custom_valence_space is set
    ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9' # AKA: Ni56 core with jj44pn
  elif A == 82 and Z == 34:
    ARGS['valence_space'] = 'pf5g9' # this is just a label when custom_valence_space is set
    ARGS['custom_valence_space'] = 'Ni56,p0f5,n0f5,p1p3,n1p3,p1p1,n1p1,p0g9,n0g9' # AKA: Ni56 core with jj44pn
  else:
    print 'these *A and Z* have not been set up yet!'
    print 'exiting...'
    exit()
  print 'v.sp.     = ',ARGS['valence_space']
  #for hw in [10.49]:
  #for hw in [9.230]:
  #for hw in [9.033]:
  #for hw in [10]:
  #for hw in [8,12,16,20,24,28,32,36,40]:
  #for hw in [12,16,20,24,28]:
  #for hw in [12,20]:
  #for hw in [16,24]:
  #for hw in [14]:
  for hw in [16]:
    ARGS['hw'] = str(hw)
    print '======================='
    print 'hw        = ',hw
    print '======================='
    #for e in [6,8,10,12,14]:
    #for e in [7,9,11,13]
    #for e in [6,8,10,12]:
    #for e in [6,8,10]:
    #for e in [6,8]:
    for e in [13]:
      ARGS['emax'] = '%d'%e
      print '-----------------------'
      print 'e         = ',e
      print '-----------------------'
      #for MNU['Decay'] in ['GT','F','T','2','2c']:
      for MNU['Decay'] in ['GT','F','T']:
      #for MNU['Decay'] in ['GT','T']:
      #for MNU['Decay'] in ['GT','F']:
      #for MNU['Decay'] in ['GT']:
      #for MNU['Decay'] in ['F']:
      #for MNU['Decay'] in ['T']:
      #for MNU['Decay'] in ['2','2c']:
      #for MNU['Decay'] in ['2']:
      #for MNU['Decay'] in ['2c']:
        print 'Decay     = ',MNU['Decay']
        ### Set the M2nu or M0nu operator(s)
        if MNU['Decay'] == '2':
          ARGS['Operators'] = 'GamowTeller'
        elif MNU['Decay'] == '2c':
          ARGS['Operators'] = 'GamowTeller'
          ARGS['OperatorsFromFile'] = 'GTMEC^1_1_0_2^/global/scratch/cgpayne/interactions/misc/MEC_Park2003nl0_HebelerLamreg394_cD-0.316_N1max12_N12max24_hbo16.dat.gz' # note: hbo = \hbar\omega
        elif MNU['Decay'] == 'GT' or MNU['Decay'] == 'F' or MNU['Decay'] == 'T':
          ARGS['Operators'] = 'M0nu_adpt_%s_%s_%s_%s_%s'%(MNU['dirout'],MNU['Decay'],MNU['Reduced'],MNU['Ec'],MNU['SRC'])
        else:
          print 'this Decay has not been set up yet!'
          print 'exiting...'
          exit()
        ### Set the interaction (int)
        if MNU['int'] == 'BARE':
          print 'running BARE version (no IMSRG evolution)...' # interaction files don't matter for BARE operator
          ARGS['2bme'] = 'none'
          ARGS['3bme'] = 'none'
          ARGS['LECs'] = 'none'
          print 'running with OS basis (no IMSRG evolution)...'
          ARGS['basis'] = 'oscillator'
        else:
          if MNU['BB'] == 'HF':
            print 'running with HF basis (no IMSRG evolution)...'
          else:
            print 'running full IMSRG version...'
          if MNU['int'] == 'magic':
            ARGS['2bme'] = '/global/scratch/exch/me2j/fromJohannes/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
          if MNU['int'] == 'mag16':
            ARGS['2bme'] = '/global/scratch/exch/ME_share/vnn_hw%d.00_kvnn10_lambda1.80_mesh_kmax_7.0_100_pc_R15.00_N15.dat_to_me2j.gz'%(hw)
          elif MNU['int'] == 'v3trans':
            ARGS['2bme'] = '/global/scratch/cgpayne/interactions/misc/TBMEA2n4lo500-srg2.0_14_28.%d_TUD.int.gz'%(hw)
          if MNU['BB'] == '3N' or MNU['BB'] == 'HF':
            print 'running for 3N...'
            if MNU['int'] == 'magic':
              ARGS['3bme'] = '/global/scratch/exch/me3j/fromJohannes/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_14_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
              ARGS['LECs'] = 'EM1.8_2.0'
            if MNU['int'] == 'mag16':
              ARGS['3bme'] = '/global/scratch/exch/ME_share/jsTNF_Nmax_16_J12max_8_hbarOmega_%d.00_Fit_cutoff_2.00_nexp_4_c1_1.00_c3_1.00_c4_1.00_cD_1.00_cE_1.00_2pi_0.00_2pi1pi_0.00_2picont_0.00_rings_0.00_J3max_9_new_E3_16_e_14_ant_EM1.8_2.0.h5_to_me3j.gz'%(hw)
              ARGS['LECs'] = 'EM1.8_2.0'
            elif MNU['int'] == 'v3trans':
              ARGS['3bme'] = '/global/scratch/cgpayne/interactions/misc/v3trans_J3T3.int_NNn4lo5003Nloc650nonloc500cD045cEm003-srg2.0_from20_330_161615_%d_form.gz'%(hw)
              ARGS['LECs'] = 'N4LO_LNL'
          else:
            print 'running for 2N...'
            ARGS['3bme'] = 'none'
            if MNU['int'] == 'magic':
              ARGS['LECs'] = 'EM1.8_2.0'
            elif MNU['int'] == 'v3trans':
              ARGS['LECs'] = 'N4LO_LNL'
        print 'LECs      = ',ARGS['LECs']
        
        ### Make an estimate of how much time to request. Only used for slurm at the moment.
        time_request = '24:00:00'
        if   e <  5 : time_request = '00:10:00'
        elif e <  8 : time_request = '01:00:00'
        elif e < 10 : time_request = '04:00:00'
        elif e < 12 : time_request = '12:00:00'
        
        jobname  = '%s_%s_%s_%s_%s_%s_e%s_E%s_s%s_hw%s_A%s_'%(ARGS['valence_space'],MNU['int'],MNU['BB'],ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['e3max'],ARGS['smax'],ARGS['hw'],ARGS['A'])
        
        ### Some optional parameters that we probably want in the output name if we're using them
        if 'lmax3' in ARGS:  jobname  += 'l%d'%(ARGS['lmax3']) + '_'
        if 'eta_criterion' in ARGS: jobname += 'eta%s'%(ARGS['eta_criterion']) + '_'
        if 'core_generator' in ARGS: jobname += '' + ARGS['core_generator'] + '_'
        if 'gapE' in ARGS: jobname += 'gap' + ARGS['gapE'] + '_'
        if 'BetaCM' in ARGS: jobname += 'BCM' + ARGS['BetaCM'] + '_'
        if MNU['Decay'] == '2': jobname += 'forM2nu_'
        if MNU['Decay'] == '2c': jobname += 'forM2nuMEC_'
        ARGS['flowfile'] = MNU['outname'] + 'BCH_' + jobname + 'flow_' + datetime.fromtimestamp(time()).strftime('%y%m%d%H%M.dat')
        ARGS['intfile']  = MNU['outname'] + jobname
        
        logname = jobname + MNU['Decay'] + '_' + datetime.fromtimestamp(time()).strftime('%y%m%d%H%M.log')
        cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])
        
        ### Submit the job if we're running in batch mode, otherwise just run in the current shell
        if batch_mode==True:
          sfile = open(jobname+'.batch','w')
          if BATCHSYS == 'PBS':
            sfile.write(FILECONTENT%(jobname,getcwd(),NTHREADS,mail_address,logname,NTHREADS,cmd))
            sfile.close()
            call(['qsub', jobname+'.batch'])
          elif BATCHSYS == 'SLURM':
            sfile.write(FILECONTENT%(NTHREADS,jobname,time_request,mail_address,NTHREADS,NTHREADS,cmd))
            sfile.close()
            call(['sbatch', jobname])
          remove(jobname+'.batch') # delete the file
          sleep(0.1)
        else:
          call(cmd.split())  # Run in the terminal, rather than submitting
        
        print '~~~~~~~~~~~~~~~~~~~~~~~'
        sleep(0.5)

