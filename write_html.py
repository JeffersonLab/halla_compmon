import os, sys

spaces = '  '

header_str = """
<!DOCTYPE html PUBLIC "-//w3c//dtd html 4.0 transitional//en">
<html>
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta name="GENERATOR" content="Mozilla/4.78 [en] (X11; U; SunOS 5.8 sun4u) [Netscape]">

"""
style_str = """
    <style type="text/css">
      body{background-color:white;}
      div{white-space:nowrap;}
      h1{font-family:"Arial"; font-size:24px; text-align:center; color:mediumslateblue;}
      h2{font-family:"Arial"; font-size:20px; color:mediumpurple;}
      h3{font-family:"Arial"; font-size:16px; color:mediumslateblue;}
      h4{font-family:"Arial"; font-size:14px; color:mediumslateblue;}
      p{font-family:"Arial"; font-size:16px;}table{border-collapse:collapse;}table,th, td{border: 1px solid black; font-family:"Arial";}td{padding: 10px ;}
      #header {background: #FFE4C4 ; border-bottom: solid 1px #aba; border-left: solid 1px #9a9;
        border-right: solid 1px #565;border-top: solid 1px #9a9; font: normal 80% \'Times New Roman\', Times, serif;
        letter-spacing: 0.0em; margin: 0; padding: 15px 10px 15px 60px;}
      a:link{color:Blue}
      a:visited{color:Blue}
    </style>

    <style>
      ul#menu {padding: 0;}
      ul#menu li {display: inline;}
      ul#menu li a {color: Mediumslateblue; padding: 10px 10px; font-family:"Arial"; font-size:16px; border-radius: 4px 4px 0 0;}
      ul#menu li a:hover {color: Peru;}
    </style>
  </head>

"""

main_body_str = """
  <body>
    <title>PREX-II Compton Analysis Plots </title>
    <h1>Compton </h1>
    <hr>
    <p>
      <h2>Compton Online Plots</h2>
    </p>

"""
end_str = """
    <hr>
    Webpage Last Updated: June 29 2019 15:21:22.<br>
  </body>
</html>
"""

dvcs_run_lo = 2063
dvcs_run_hi = 3108
test_run_lo = dvcs_run_lo + 1
test_run_hi = 4231
prex_run_lo = test_run_hi + 1
prex_run_hi = 9999

#Steps
#1. Open HTML file
#2. Look for current run number
#2a. If run number exists in file, ignore and leave HTML be
#2b. If not, proceed
#3. Create a <li> block for the run number's plots.
#4. Divide up HTML based on <!--SECTION--> tags
#5. Sections should be as follows:
#5a. 0 = Header (with meta tags)
#5b. 1 = Style blocks
#5c. 2 = Body opening block
#5d. 3 = Each <ul> block (inside a <div>)
#5e. 4 = HTML footer
#6. List run directories, add to block d each run that's there
#7. Stitch all the strings back together in order
#Run limits:
#2063 - 3108 = DVCS
#3109 - 4231 = Tests
#4232 - ???? = PREX

def create_plot_list(runs_to_write):
  list_elements = 3*spaces + '<ul>\n'
  list_files = [('snapshots.pdf', 'Snapshot Plots'), 
                ('sums.pdf', 'Triggered Sums'), 
                ('acc0.pdf', 'Acc0/NAcc0'), 
                ('acc0_time.pdf', 'Acc0/NAcc0 vs Time'),
                ('asym_hists.pdf', 'Quartet Histograms'),
                ('q_graphs.pdf', 'Quartet Vars vs Time')]
  for run in runs_to_write:
    list_elements += 4*spaces + '<li>Run ' + str(run) + ': &ensp;\n'
    for data in list_files:
      list_elements += 5*spaces + '<a href=\'runs/Run' + str(run) + '/' + data[0] + '\'>' + data[1] + '</a>&ensp;\n'
    list_elements += 4*spaces + '</li>\n'
  return list_elements + 3*spaces + '</ul>\n'

def create_prex_plot_list(runs_to_write):
  list_files = [('ess_stats.pdf', 'Essential Stats'),
                ('snapshots.pdf', 'Snapshot Plots'), 
                ('sums.pdf', 'Triggered Sums'), 
                ('acc0.pdf', 'Acc0/NAcc0'), 
                ('acc0_time.pdf', 'Acc0/NAcc0 vs Time'),
                ('asym_acc0_hists.pdf', 'Multiplet Acc0 Histograms'),
                ('q_acc0_graphs.pdf', 'Multiplet Acc0 Vars vs Time'),
				        ('asym_acc4_hists.pdf', 'Multiplet Acc4 Histograms'),
				        ('q_acc4_graphs.pdf', 'Multiplet Acc4 Vars vs Time'),
                ('cycle_qVars.pdf', 'Laser Cycles')]
  snail_list = os.listdir(os.environ['COMPMON_SNAILS']);
  snail_list_trim = []
  for el in snail_list:
    if "~" not in el and '#' not in el: snail_list_trim += [el]
  snail_list = sorted(snail_list_trim, key=lambda fname: int(fname.replace("snail","").replace(".list", ""))); snail_list.reverse()
  snails_and_runs = []
  snail_strs = {}
  for snail in snail_list:
   # print(snail)
    title = "Snail " + snail.replace("snail", "").replace(".list", "")
    snail_data = [title]
    snail_file = open(os.environ['COMPMON_SNAILS'] + '/' + snail)
    for line in snail_file.readlines():
      if line == '' or line == '\n' or line == ' ': continue
      snail_data += [int(line)]
    snail_data += [False]
    snails_and_runs += [snail_data]
    snail_strs[title] = 3*spaces + '<h4>' + title + '</h4>\n' + 3*spaces + '<ul>\n'
  
  for run in runs_to_write:
    for snail in snails_and_runs:
      if not snail[-1]:
        snail_strs[snail[0]] += 4*spaces + '<li>Snail Summary: &ensp;\n'
        #snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc0.pdf\'>Polarization (Acc0)</a>&ensp;\n'
        #snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/asymmetries_acc0.pdf\'>Asymmetries (Acc0)</a>&ensp;\n'
        #snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc4.pdf\'>Polarization (Acc4)</a>&ensp;\n'
        #snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/asymmetries_acc4.pdf\'>Asymmetries (Acc4)</a>&ensp;\n'
        snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc0_cycles.pdf\'>Polarization (Acc0)</a>&ensp;\n'
        snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc4_cycles.pdf\'>Polarization (Acc4)</a>&ensp;\n'
        snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc0.pdf\'>Polarization (Acc0, miniruns)</a>&ensp;\n'
        snail_strs[snail[0]] += 5*spaces + '<a href=\'snails/snail' + snail[0].split(' ')[1] + '/polarization_acc4.pdf\'>Polarization (Acc4, miniruns)</a>&ensp;\n'
        snail_strs[snail[0]] += 4*spaces + '</li>\n'
        snail[-1] = True
      if run in snail:
        snail_strs[snail[0]] += 4*spaces + '<li>Run ' + str(run) + ': &ensp;\n'
        for data in list_files:
          snail_strs[snail[0]] += 5*spaces + '<a href=\'runs/Run' + str(run) + '/' + data[0] + '\'>' + data[1] + '</a>&ensp;\n'
        snail_strs[snail[0]] += 4*spaces + '</li>\n'
  
  data_str = ''
  for snail in snails_and_runs:
    data_str += snail_strs[snail[0]] + 3*spaces + '</ul>\n'
  return data_str

def create_end_block(date, time):
  block_str = 2*spaces + '<hr>\n'
  block_str += 2*spaces + '<a href=\'https://hallaweb.jlab.org/parity/prex/onlinePlots/\'>Online Plots</a>&ensp;\n'
  block_str += 2*spaces + '<hr>\n'
  block_str += 2*spaces + 'Web page last edited: ' + date + ' at ' + time + '<br>\n'
  block_str += spaces + '</body>\n</html>'
  return block_str

def write_html():
  print('\nUpdating compton plots html...')
  split_str = '<!--SECTION-->\n'
  date = sys.argv[1]
  time = sys.argv[2]
  fname = sys.argv[3]
  comp_web_path = os.environ["COMPMON_WEB"]
  f = open(comp_web_path + "/" + fname)
  html = f.read()
  blocks = html.split(split_str)
  run_list = os.listdir(comp_web_path + "/runs/")
  dvcs_runs = []; test_runs = []; prex_runs = []
  for run in run_list:
    run_num = int(run.replace("Run", ""))
    if run_num < dvcs_run_lo: continue
    elif dvcs_run_lo <= run_num <= dvcs_run_hi: dvcs_runs += [run_num]
    elif test_run_lo <= run_num <= test_run_hi: test_runs += [run_num]
    elif prex_run_lo <= run_num <= prex_run_hi: prex_runs += [run_num]
    else: continue
  dvcs_runs.sort(); test_runs.sort(); prex_runs.sort()
  dvcs_runs.reverse(); test_runs.reverse(); prex_runs.reverse()
  dvcs_block = 2*spaces + '<div>\n' + 3*spaces + '<h3>DVCS Runs</h3>\n' + create_plot_list(dvcs_runs) + 2*spaces + '</div>\n'
  test_block = 2*spaces + '<div>\n' + 3*spaces + '<h3>Test Runs</h3>\n' + create_plot_list(test_runs) + 2*spaces + '</div>\n'
  prex_block = 2*spaces + '<div>\n' + 3*spaces + '<h3>PREX Runs</h3>\n' + create_prex_plot_list(prex_runs) + 2*spaces + '</div>\n'
  #list_block = prex_block + test_block + dvcs_block
  list_block = prex_block
  f.close();
  
  all_data = header_str + split_str + style_str + split_str + main_body_str + split_str
  all_data += list_block + split_str + create_end_block(date, time)

  fout = open(comp_web_path + "/" + fname, 'w+')
  fout.write(all_data + '\n')
  fout.close()
  print("...Done!")
  
write_html()
