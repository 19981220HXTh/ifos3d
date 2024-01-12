#!/bin/bash

#### write to plot 
#### weiyw17@gmail.com
#### 2017.9.6

## set parameters
#your_job="/home/export/online1/wyw/workspace/GPIAG/swcbw/par/toy-small-40-46-40"
vpfile_path=$(cd "$(dirname "$0")"; pwd)
#vpfile_path=${your_job}/model
vpfile_true=toy.small.vp
vpfile_name=toy.vp_it60
#vpfile_name=$1
rsffile=toy.rsf
rsftrue=real.rsf
output_d2="/home/export/online1/wyw/workspace/GPIAG/swcbw/par/toy/model/d2.vpl"
output_d1="/home/export/online1/wyw/workspace/GPIAG/swcbw/par/toy/model/d1.vpl"
N1=184
N2=160
N3=160
y_pos=92
n2_pos=80
#N1=$2
#N2=$3
#N3=$4
#y_pos=$5
#n2_pos=$6

## get the right rsf file
#cd ${vpfile_path}
sfmath n1=${N1} n2=${N2} n3=${N3} output=0 > ${rsffile}
sfmath n1=${N1} n2=${N2} n3=${N3} output=0 > ${rsftrue}

sed "s/in=.*/in=${vpfile_name}/g" ${rsffile} > t.rsf
cat t.rsf > ${rsffile}
rm t.rsf

sed "s/in=.*/in=${vpfile_true}/g" ${rsftrue} > t.rsf
cat t.rsf > ${rsftrue}
rm t.rsf

## make graph

BIAS=`sfattr < ${rsffile} want='mean' | awk '{print $3}'`

# 2d
d2title="2D p-velocity in y=${y_pos}"
sfwindow f1=${y_pos} n1=1 < ${rsffile} | sfgrey bias=${BIAS} scalebar=y color=j \
label1="Depth" unit1=m label2="Distance" unit2=m \
title="${d2title}" \
| sfpen &
#> ${output_d2}

# 1d
d1title="1D p-velocity contrast in y=${y_pos} x=${n2_pos}"
sfwindow f1=${y_pos} n1=1 f2=${n2_pos} n2=1 < ${rsffile} > file1.rsf
sfwindow f1=${y_pos} n1=1 f2=${n2_pos} n2=1 < ${rsftrue} > file2.rsf
sfcat file1.rsf file2.rsf axis=2 \
| sfgraph color=j \
label2="Velocity" unit2=m/s label1="Depth" unit1=m \
title="${d1title}" \
dash=0,1  grid=y \
| sfpen
#> ${output_d1}

## clean 
sfrm file1.rsf file2.rsf
