from subprocess import Popen, PIPE
from datetime import datetime, timedelta
import os, sys

def write_position_file():
  print('Printing the BPM information...')
  infile_name = os.environ['COMPMON_RUNPLOTS'] + '/Run' + str(sys.argv[1]) + '_time.txt'
  infile = open(infile_name, 'r')

  lines = infile.readlines()
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
  cur_range = '110:160'
  if(int(sys.argv[1]) >= 4300 and int(sys.argv[1]) <= 4800):
    run_seconds = run_events*1.0/240
    cur_range = '65:88'

  start_time = datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
  add_time = timedelta(seconds=int(run_seconds))
  end_time = start_time + add_time

  fn_call = ['/adaqfs/home/apar/bin/myStats', '-b', '' + str(start_time) + '', '-e', '' + str(end_time) + '', '-c', 'IBC1H04CRCUR2', '-r', cur_range, '-l', 'IPM1P02A.XPOS,IPM1P02A.YPOS,IPM1P02B.XPOS,IPM1P02B.YPOS']
  process = Popen(fn_call, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  output_lines = stdout.decode('utf-8').split('\n')
  
  outfile_name = os.environ['COMPMON_RUNPLOTS'] + '/Run' + str(sys.argv[1]) + '_bpms.csv'
  outfile = open(outfile_name, 'w+')
  
  for i in range(2, 6):
    output_line = [j for j in output_lines[i].split(' ') if j is not '']
    outfile.write(output_line[2] + ',' + output_line[4] + '\n')

  infile.close()
  outfile.close()
  print('...Done!')


write_position_file()
