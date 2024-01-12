#!/usr/bin/python
import os
import sys
import shutil
import json
import json_crt
FW_JSON = json_crt.param_fw['JSON_FILE']
INV_JSON = json_crt.param_inv['JSON_FILE']
#QUEUE_NAME = 'q_x86_share'
QUEUE_NAME = 'q_sw_expr'
JOB_NAME = 'hmy.ifos3d'

rm_dirs = ['grad', 'hess', 'su', 'su_obs']
mk_dirs = ['model', 'snap', 'log']
files = [FW_JSON, INV_JSON, 'bsub.out']

def init():
  for f in rm_dirs:
    if os.path.exists(f):
      shutil.rmtree(f)
      os.mkdir(f)
    else:
      os.mkdir(f)
  for f in mk_dirs:
    if not os.path.exists(f):
      os.mkdir(f)

def clean():
  for f in rm_dirs:
    if os.path.exists(f):
      shutil.rmtree(f)
  for f in files:
    if os.path.exists(f):
      os.system("rm -f %s" % f)

def read_json(json_file):
  with open(json_file) as f:
    data = json.load(f)
    data['JSON_FILE'] = json_file
    return data

def compile(path, para):# {{{
  cmd = "scons -C %s %s -j8" % (path, para)
  print cmd
  os.system(cmd)
# }}}  

def bsub_x86(queue, np, job, prog):# {{{
  cmd = "/usr/sw-mpp/bin/bsub -I -q %s -n %d -o ./log/`date +'%y%m%d%H%M'`.log -J %s %s" % (queue, np, job, prog)

  return cmd
# }}}
def bsub_sw(queue, np, job, prog):# {{{
  cmd="bsub -I -b -m 1 -p -q %s -host_stack 1024 -share_size 5000 -n %d -cgsp 64 -J %s %s" % (queue, np, job, prog)
  return cmd
# }}}
def run(submit_cmd, param):# {{{
  cmd = submit_cmd.rstrip() + '''\
    %s | tee  ./log/`date +%s`-%s
    ''' % (param['JSON_FILE'], "%y%m%d%H%M", param['LOG_FILE'])
  print cmd
  os.system(cmd)
# }}}

if '-c' in sys.argv:
  clean()
  exit()

exe="../../bin/ifos3d"
if 'fw' in sys.argv:
  init()
  compile('../../libcseife/', '')
  #compile('../../', '')
  compile('../../', 'model_scr=hh_toy_true.c')
  param = read_json(FW_JSON)

  if 'q' in QUEUE_NAME:
    run(bsub_sw(queue=QUEUE_NAME, np=int(param['NPROCX'])*int(param['NPROCY'])*int(param['NPROCZ']), job=JOB_NAME, prog=exe), param)
  elif 'x86' in QUEUE_NAME:
    run(bsub_x86(queue=QUEUE_NAME, np=int(param['NPROCX'])*int(param['NPROCY'])*int(param['NPROCZ']), job=JOB_NAME, prog=exe), param)

  shutil.copy('model/%s.vs_it0' % param['PREFIX'], 'model/%s.vs.true' % param['PREFIX'])
  shutil.copy('model/%s.vp_it0' % param['PREFIX'], 'model/%s.vp.true' % param['PREFIX'])
  shutil.copy('model/%s.rho_it0' % param['PREFIX'], 'model/%s.rho.true' % param['PREFIX'])

elif 'inv' in sys.argv:
  compile('../../libcseife', '')
  #compile('../../', '')
  compile('../../', 'model_scr=hh_toy_start.c')
  param = read_json(INV_JSON)

  if 'q' in QUEUE_NAME:
    run(bsub_sw(queue=QUEUE_NAME, np=int(param['NPROCX'])*int(param['NPROCY'])*int(param['NPROCZ']), job=JOB_NAME, prog=exe), param)
  elif 'x86' in QUEUE_NAME:
    run(bsub_x86(queue=QUEUE_NAME, np=int(param['NPROCX'])*int(param['NPROCY'])*int(param['NPROCZ']), job=JOB_NAME, prog=exe), param)

else:
  print 'Usage: ./run.py fw or ./run.py inv'
