#!/bin/bash
testresult=/home/export/online1/wanwb/lasefis_orig/1_ifos3d_test/par/toy/su_obs/obs_vx_it1.su.shot1.rsf 
original=/home/export/online1/wanwb/lasefis/ifos3d/par/toy/su_obs/obs_vx_it1.su.shot1.rsf
 sfattr < /home/export/online1/wanwb/lasefis_orig/1_ifos3d_test/par/toy/su_obs/obs_vx_it1.su.shot1.rsf 
 sfattr < /home/export/online1/wanwb/lasefis/ifos3d/par/toy/su_obs/obs_vx_it1.su.shot1.rsf
 sfin /home/export/online1/wanwb/lasefis_orig/1_ifos3d_test/par/toy/su_obs/obs_vx_it1.su.shot1.rsf 
 sfin < /home/export/online1/wanwb/lasefis/ifos3d/par/toy/su_obs/obs_vx_it1.su.shot1.rsf
 sfdiffcxx < ${testresult} match=${original} > diff.rsf
