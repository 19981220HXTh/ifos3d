## plot seismograms
## created by weiyw		contact: weiyw17@gmail.com
## 2017.9.10

#!/bin/bash
#cd ./su_obs
# namespaces
FILE_SU=$1
############################### modify ##############################
#FILE_SU=obs_toy_vy_it1.su.shot1
#####################################################################
FILE_RSF=${FILE_SU}.rsf
#sfsegyread tape=${FILE_SU} su=y endian=n format=5 > ${FILE_RSF}
sfsegyread tape=${FILE_SU} su=y endian=n format=5 datapath=/home/export/online1/mdt00/shisuan/sweq/eq/workspace/hmy/ifos3d/ifos3d_dma_test/par/toy-1/su_obs/ > ${FILE_RSF}
LINE_REC_NUM=`sfget < ${FILE_RSF} parform=n n2 `
#awk 'BEGIN{J2=sqrt(${LINE_REC_NUM});}'
J2=` echo "sqrt(${LINE_REC_NUM}) - 1" | bc `
echo ${J2}
sfwindow j2=${J2} < ${FILE_RSF} | sfgrey plicp=37 title="j" | sfpen &
sfwindow j2=21  < ${FILE_RSF} | sfgrey plicp=37 title="mj" | sfpen &
sfwindow n2=${J2} < ${FILE_RSF} | sfgrey plicp=37 title="n" | sfpen &
sfgrey < ${FILE_RSF} plicp=37 title="a" | sfpen
