from subprocess import Popen, PIPE
from datetime import datetime, timedelta
import os, sys


def get_snail_file_list():
  comp_web_path = os.environ["COMPMON_WEB"]
  run_mode = 'prex'
  if 'crex' in comp_web_path: run_mode = 'crex'

  snail_list = os.listdir(os.environ['COMPMON_SNAILS'])
  snail_list_trim = []
  for el in snail_list:
    if "~" not in el and '#' not in el: 
      snail_num = int(el.replace('.list', '').replace('snail', ''))
      if run_mode == 'prex' and (snail_num < 100 or snail_num == 500):
        snail_list_trim += [el]
      elif run_mode == 'crex' and (snail_num >= 100 and snail_num != 500):
        snail_list_trim += [el]

  snail_list = sorted(snail_list_trim, key=lambda fname: int(fname.replace("snail","").replace(".list", "")))
   
  return snail_list


def calc_end_time(run_num):
  time_file = open(os.environ['COMPMON_RUNPLOTS'] + '/Run' + str(run_num) + '_time.txt', 'r')
  lines = time_file.readlines()
  full_datetime = lines[0]
  run_events = int(lines[1])

  full_date = lines[0].split(' ')[0]
  full_time = lines[0].split(' ')[1]
  
  year = int(full_date.split('-')[0])
  month = int(full_date.split('-')[1])
  day = int(full_date.split('-')[2])

  hour = int(full_time.split(':')[0])
  minute = int(full_time.split(':')[1])
  second = int(full_time.split(':')[2])

  run_seconds = run_events*1.0/120
  if(run_num >= 4300 and run_num <= 4800):
    run_seconds = run_events*1.0/240

  start_time = datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
  add_time = timedelta(seconds=int(run_seconds))
  end_time = start_time + add_time

  return start_time, end_time


def check_run_ihwp(run_num):
  start_time, end_time = calc_end_time(run_num)

  fn_call = ['/adaqfs/home/apar/bin/myData', '-b', '' + str(start_time) + '', '-e', '' + str(end_time) + '', 'IGL1I00OD16_16']
  process = Popen(fn_call, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  output_lines = stdout.decode('utf-8').split('\n')
  error_lines = stderr.decode('utf-8').split('\n')
  if len(output_lines) > 3:
    print('IHWP Flip Detected for Run ' + str(run_num))
    print('  Run started at ' + str(start_time))
    print('  Run ended at ' + str(end_time))
    for i in range(2, len(output_lines) - 1):
      flip_info = filter(lambda x: x!= '', output_lines[i].split(' '))
      flip_date = flip_info[0].split('-')
      flip_time = flip_info[1].split(':')
      flip_datetime = datetime(year=int(flip_date[0]), month=int(flip_date[1]), day=int(flip_date[2]), hour=int(flip_time[0]), minute=int(flip_time[1]), second=int(flip_time[2]))
      print('  Flip occurred at: ' + str(flip_datetime))
      print('    ' + str((flip_datetime - start_time).seconds) + ' seconds after start of run')
      print('    ' + str((end_time - flip_datetime).seconds) + ' seconds before end of run')


def check_all_runs_ihwp():
  snail_file_list = get_snail_file_list()
  for fname in snail_file_list:
    snail_file = open(os.environ['COMPMON_SNAILS'] + '/' + fname, 'r')
    for line in snail_file.readlines():
      if line == '\n' or line == '': continue
      run_num = int(line.replace('\n', ''))
      check_run_ihwp(run_num)
      


check_all_runs_ihwp()

