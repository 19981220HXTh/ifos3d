#!/usr/bin/python
import os
import json
prefix = "toy"
param_fw = {
  'JSON_FILE' : 'fw.json',
  "MODELING PARAMETERS" : "comment",
  "PREFIX" : "%s" % prefix,
  "LOG_FILE" : "log/fw.log",

  "Note that y denotes the vertical direction !" : "comment",

  "Domain Decomposition" : "comment",
        "NPROCX" : "4",
        "NPROCY" : "4",
        "NPROCZ" : "4",

  "3-D Grid" : "comment",
        "NX" : "160",
        "NY" : "184",
        "NZ" : "160",
        "DX" : "0.8",
        "DY" : "0.8",
        "DZ" : "0.8",

  "FD order" : "comment",
        "FDORDER" : "4",
        "FDCOEFF" : "2",

  "Time Stepping" : "comment",
        "TIME" : "0.06",
        "DT" : "5.0e-05",

  "Source" : "comment",
        "SOURCE_SHAPE" : "4",
        "SOURCE_TYPE" : "4",
        "SRCREC" : "1",
        "SOURCE_FILE" : "input/sources.dat",
        "RUN_MULTIPLE_SHOTS" : "1",

  "Model" : "comment",
        "READMOD" : "0",
        "MFILE" : "model/%s" % prefix,

  "Q-approximation" : "comment",
        "L" : "0",

  "Free Surface" : "comment",
        "FREE_SURF" : "0",
        

  "Absorbing Boundary" : "comment",
        "ABS_TYPE" : "1",
        "FW" : "10",
        "VPPML" : "6200.0",
        "FPML" : "200.00",

  "Snapshots" : "comment",
        "SNAP" : "0",

  "Receiver" : "comment",
        "SEISMO" : "1",
        "SAVESU" : "1",
        "READREC" : "2",

  "Receiver array" : "comment",
        "REC_ARRAY" : "1",
        "REC_ARRAY_DEPTH" : "24.0",
        "REC_ARRAY_DIST" : "30.0",
        "DRX" : "10",
        "DRZ" : "10",

  "Seismograms" : "comment",
        "NDT" : "1",
        "SEIS_FORMAT" : "1",
        "SEIS_FILE" : "./su_obs/obs",


  "Method" : "comment",
        "METHOD" : "0",
        
        "MOD_OUT_FILE" : "./model/%s" % prefix
}

param_inv = {
  'JSON_FILE' : 'inv.json',
  "LOG_FILE" : "log/inv.log",
  "Seismograms" : "comment",
        "SEIS_FILE" : "./su/cal",
        "SAVESU" : "0",

  "Method" : "comment",
        "METHOD" : "1",

  "INVERSION PARAMETERS" : "comment",

  "In- and Output Files" : "comment",
        "GRADMO" : "0",
        "GRAD_FILE" : "./grad/grad",
        "MOD_OUT_FILE" : "./model/%s" % prefix,
        "SEIS_OBS_FILE" : "./su_obs/obs",
        "EXTOBS" : "0",
        "INV_FILE" : "input/workflow.dat",

  "General" : "comment",
        "ITMIN, ITMAX" : "1 , 2",
        "FILT" : "1",
        "NFMAX" : "5",
        "TAST" : "100",
        "VP0, VS0, RHO0" : "6200.0, 3600.0, 2800.0",
        "WEIGHT_VP,WEIGHT_VS,WEIGHT_RHO" : "1.0, 1.0, 0.0",

  "Steplength estimation" : "comment",
        "NSHOTS_STEP" : "4",
        "TESTSTEP" : "0.02",

  "Gradient preconditioning" : "comment",
        "DAMPTYPE" : "2",

  "Hessian" : "comment",
        "HESS" : "0",
   
  "L-BFGS" : "comment",
        "LBFGS" : "0"
}

with open(param_fw['JSON_FILE'], 'w') as json_file:
  json_file.write(json.dumps(param_fw, indent = 1))
param_fwt = param_fw.copy()
param_fwt.update(param_inv)  
with open(param_inv['JSON_FILE'], 'w') as json_file:
  json_file.write(json.dumps(param_fwt, indent = 1))
