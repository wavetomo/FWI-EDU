#!/bin/bash
###
 # @Author: QuYang Chen
 # @Date: 2025-9-23 18:38:04
 # @LastEditors: QuYang Chen
 # @LastEditTime: 2025-12-23 17:03:20
 # @FilePath: /FWILab2dCLangV1.2.5/demo/marm-ii/fwi.sh
 # @Description: 
 # 
 # Copyright (c) 2025 by WaveTomo, All Rights Reserved.
### 

#SBATCH -J fwi
# Time required (hh:mm:ss)
#SBATCH -t 72:00:00
# Resource requirements
#SBATCH -n 1
#SBATCH -p compute2
#SBATCH --cpus-per-task=56
# rename stdout and sterr
#SBATCH -o fwi.log
#SBATCH -e fwi.err

# Model parameters
nx=601
nz=221
nt=6001
dx=12.5
dz=12.5
dt=0.001

# Boundary conditions
freesurface=1
pmlThick=30

# Path and file
pwd=$(pwd)
files=$pwd/model_geometry
work_path=$pwd/fwi_free-surface
srcfile=$files/ricker_6001_1ms_10Hz_delay0.15ms.bin
vpfile=$files/vp_smooth_00221_00601_12.5m.bin
rhofile=
recordfile=$pwd/modeling_free-surface/shots_p_30_3001_18030_02.00ms.su

exefile=$pwd/../../bin/acoustic2d_fwi.exe

# Water layer 
Gardner=1
vpmax=4700
vpmin=1500
rhomax=2700
rhomin=1000
vp_water=1500.0
rho_water=1010.0
waterdepth=36
precondition=1

# Source parameters
Nshot=30

# Iteration parameters
freq=2,4,6,8,10
filter=0,0,0,0,2
iteration=10,5,5,5,5
rx=2,2,2,1,1
rz=2,2,2,1,1


# Output info
verbose=1

# output intervals for shot and time of wavefield
export shotInterval=5
export snapInterval=1000 # verbose=2

rm -rf $work_path
mkdir -p $work_path
cd $work_path

# Run executable
#valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-file=valgrind.log --show-leak-kinds=all \
$exefile \
    nx=$nx nz=$nz nt=$nt \
    dx=$dx dz=$dz dt=$dt \
    freesurface=$freesurface pmlThick=$pmlThick \
    srcfile=$srcfile vpfile=$vpfile rhofile=$rhofile recordfile=$recordfile \
    waterdepth=$waterdepth vp_water=$vp_water rho_water=$rho_water Gardner=$Gardner \
    Nshot=$Nshot \
    vpmin=$vpmin vpmax=$vpmax rhomax=$rhomax rhomin=$rhomin \
    freq=$freq filter=$filter iteration=$iteration rx=$rx rz=$rz \
    precondition=$precondition \
    verbose=$verbose
