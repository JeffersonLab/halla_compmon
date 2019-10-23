import os, sys

spaces = '  '

preamble = """
plotsdir ./output
plotFormat pdf

guicolor mediumpurple

protorootfile $COMP_ROOTFILES/compmon_XXXXX.root
"""

def writeCFG(run_num):
  print('Writing config for run ' + run_num + '...')
  write_str = '\n\n'
  reg_plots = [('Essential Stats', 1), ('Snapshots', 2), ('Spectra', 3), ('Acc0', 4), ('Acc0 vs Time', 5), \
               ('PosHelAcc0', 6), ('NegHelAcc0', 7), ('Multiplet Diffs', 8), ('Multiplet Sums', 9), ('Multiplet Asyms', 10), \
               ('Asym and Polarization', 11), ('Asyms vs Time', 12), ('Acc4 Asym and Polarization', 13), \
               ('Acc4 Asyms vs Time', 14), ('Beam OFF Asym', 15)]
  #asym_plots = [('Asym and Polarization', 1), ('Asyms vs Time', 2), ('Acc4 Asym and Polarization', 4), \
  #              ('Acc4 Asyms vs Time', 5), ('Beam OFF Asym', 6)]

  for plot in reg_plots:
    write_str += 'newpage 1 1\n'
    write_str += 4*spaces + 'title ' + plot[0] + '\n'
    write_str += 4*spaces + 'macro $COMPMON_PANGUIN/macros/plotPanguin.C(' + run_num + ',' + str(plot[1]) + ');\n'

  #for plot in asym_plots:
  #  write_str += 'newpage 1 1\n'
  #  write_str += 4*spaces + 'title ' + plot[0] + '\n'
  #  write_str += 4*spaces + 'macro $COMPMON_PANGUIN/macros/asymOnlinePlot.C(' + run_num + ',' + str(plot[1]) + ');\n'

  path = os.environ["COMPMON_PANGUIN"]
  f = open(path + '/macros/prex_auto_runs.cfg', 'w+')
  f.write(preamble)
  f.write(write_str)
  f.close()

writeCFG(sys.argv[1])
