#!/bin/bash
###
 # @Author: QuYang Chen
 # @Date: 2025-9-23 18:38:04
 # @LastEditors: Quyang Chen
 # @LastEditTime: 2025-12-29 14:31:32
 # @FilePath: /FWILab2dCLangV1.2.5/FWILab2dCLangV1.2.5/demo/marm-ii/modeling.sh
 # @Description: 
 # 
 # Copyright (c) 2025 by WaveTomo, All Rights Reserved.
### 

#SBATCH -J modeling
# Time required (hh:mm:ss)
#SBATCH -t 72:00:00
# Resource requirements
#SBATCH -n 1
#SBATCH -p compute2
#SBATCH --cpus-per-task=56
# rename stdout and sterr
#SBATCH -o modeling.log
#SBATCH -e modeling.err


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

# Model files
pwd=$(pwd)
files=$pwd/model_geometry
work_path=$pwd/modeling_free-surface
srcfile=$files/ricker_6001_1ms_10Hz_delay0.15ms.bin
vpfile=$files/vp_00221_00601_12.5m.bin
rhofile=

exefile=$pwd/../../bin/acoustic2d_modeling.exe \

# Water layer 
Gardner=1
vp_water=1500.0
rho_water=1010.0
waterdepth=36   #(unit: cell)

# Source parameters 
fsx=0          #(unit: cell)
ds=60          #(unit: cell)
sdepth=2       #(unit: cell)

# Receiver parameters
offsmin=0      #(unit: cell)
offsmax=600    #(unit: cell)
dtr=1          #(unit: cell)
gdepth=2       #(unit: cell)

Nshot=30     # independant parameter
sampleRate=2 # independant parameter -- reduce memory cost

# survey geometry in JSON file
jsonfile=$files/acquisition.json
# jsonfile=

# Output info
verbose=1

# output intervals for shot and time of wavefield
export shotInterval=5  # interval of outputing shot records
export snapInterval=1000 

export OMP_NUM_THREAD=1

rm -rf $work_path
mkdir -p $work_path
cd $work_path

# Run executable
# valgrind --tool=memcheck --leak-check=full --track-origins=yes --log-file=valgrind.log --show-leak-kinds=all \
$exefile \
    nx=$nx nz=$nz nt=$nt \
    dx=$dx dz=$dz dt=$dt \
    freesurface=$freesurface pmlThick=$pmlThick \
    srcfile=$srcfile vpfile=$vpfile rhofile=$rhofile \
    waterdepth=$waterdepth vp_water=$vp_water rho_water=$rho_water Gardner=$Gardner \
    Nshot=$Nshot fsx=$fsx ds=$ds sdepth=$sdepth \
    sampleRate=$sampleRate offsmin=$offsmin offsmax=$offsmax dtr=$dtr gdepth=$gdepth \
    jsonfile=$jsonfile \
    verbose=$verbose

